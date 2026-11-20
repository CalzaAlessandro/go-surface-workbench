"""Resolve atomic C/O and optional group preferences to explicit group counts.

This selects a composition, not an equilibrium structure or measured group
distribution. Spatial placement remains the responsibility of build_sheet.
"""

from __future__ import annotations

import copy
import math
import random
from collections import Counter
from numbers import Real

import numpy as np

from .build import GROUPS, _bipartite_matching, _max_matching, build_sheet
from .model import adjacency

DEFAULT_WEIGHTS = {
    "OH": 1.0,
    "epoxide": 1.0,
    "carbonyl": 0.0,
    "COOH": 0.0,
    "COO": 0.0,
    "epoxide_eph": 0.0,
    "epoxide_ept": 0.0,
}


def _positive_number(value, label, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number")
    if value < 0 or (not allow_zero and value == 0):
        raise ValueError(f"{label} must be {'nonnegative' if allow_zero else 'positive'}")
    return float(value)


def _finite_real(value, label):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite number")
    return float(value)


def resolve_composition(
    recipe,
    c_over_o=None,
    preferences=None,
    *,
    o_over_c: float | list[float] | tuple[float, float] | None = None,
    priority: list[str] | tuple[str, ...] | None = None,
    mode: str = "balanced",
    enabled_groups: list[str] | tuple[str, ...] | None = None,
    seed: int | None = None,
    ph: float | None = None,
    pka: float = 4.5,
    deprotonation_fraction: float | None = None,
    _scaffold: dict | None = None,
    _jitter_seed: int | None = None,
    _weight_multipliers: dict[str, float] | None = None,
):
    """Return (exact-count recipe, composition report).

    A C/O (or O/C) interval is hard. A scalar allows half the largest enabled
    group's oxygen increment for rounding, measured against final carbon count.
    Preferences are nonnegative relative group-count weights: zero disables a
    group. Alternatively pass `priority` as an ordered list of groups in
    decreasing preference order (e.g. ["OH", "epoxide"]). `mode` may be
    "balanced" (default) or "random" (samples reproducible weights over enabled
    groups using `seed` or `recipe.get("seed", 0)`). Optionally specify `ph`
    (with `pka`, default 4.5) or `deprotonation_fraction` in [0, 1] to partition
    enabled carboxyl weight (`COOH` + `COO` + `carboxyl`) via Henderson-Hasselbalch.
    """
    try:
        from scipy.optimize import Bounds, LinearConstraint, milp
    except ImportError as error:
        raise ImportError("Composition resolution requires the 'composition' extra") from error
    if not isinstance(recipe, dict):
        raise ValueError("recipe must be an object")
    if recipe.get("groups"):
        raise ValueError(
            "Use C/O preferences or explicit groups; existing groups are not overwritten"
        )
    if c_over_o is not None and o_over_c is not None:
        raise ValueError("Specify either c_over_o or o_over_c, not both")
    if c_over_o is None and o_over_c is None:
        raise ValueError("Either c_over_o or o_over_c must be provided")
    if mode not in ("balanced", "random"):
        raise ValueError("mode must be 'balanced' or 'random'")
    if preferences is not None and priority is not None:
        raise ValueError("Specify either preferences or priority, not both")
    if ph is not None and deprotonation_fraction is not None:
        raise ValueError("Specify either ph or deprotonation_fraction, not both")

    effective_seed = recipe.get("seed", 0) if seed is None else seed
    if type(effective_seed) is not int or effective_seed < 0:
        raise ValueError("seed must be an integer >= 0")

    pka_val = _finite_real(pka, "pka")
    ph_val = None if ph is None else _finite_real(ph, "ph")
    target_deprot = None
    if deprotonation_fraction is not None:
        target_deprot = _finite_real(deprotonation_fraction, "deprotonation_fraction")
        if not 0.0 <= target_deprot <= 1.0:
            raise ValueError("deprotonation_fraction must be between 0 and 1")
    elif ph_val is not None:
        exponent = pka_val - ph_val
        if exponent > 300:
            target_deprot = 0.0
        elif exponent < -300:
            target_deprot = 1.0
        else:
            target_deprot = 1.0 / (1.0 + 10.0**exponent)

    requested_o_over_c = None
    if o_over_c is not None:
        interval = isinstance(o_over_c, (list, tuple))
        if interval:
            if len(o_over_c) != 2:
                raise ValueError("O/C interval must contain two positive values")
            o_low, o_high = (_positive_number(x, "O/C") for x in o_over_c)
            if o_low > o_high:
                raise ValueError("O/C interval lower bound exceeds upper bound")
            requested_o_over_c = [o_low, o_high]
            low, high = 1.0 / o_high, 1.0 / o_low
        else:
            o_ratio = _positive_number(o_over_c, "O/C")
            requested_o_over_c = o_ratio
            low = high = 1.0 / o_ratio
    else:
        interval = isinstance(c_over_o, (list, tuple))
        if interval:
            if len(c_over_o) != 2:
                raise ValueError("C/O interval must contain two positive values")
            low, high = (_positive_number(x, "C/O") for x in c_over_o)
            if low > high:
                raise ValueError("C/O interval lower bound exceeds upper bound")
        else:
            low = high = _positive_number(c_over_o, "C/O")

    allowed_prefs = set(GROUPS) | {"carboxyl"}
    if preferences is not None and (
        not isinstance(preferences, dict) or set(preferences) - allowed_prefs
    ):
        raise ValueError(
            "Unknown group preference; choose OH, epoxide, carbonyl, COOH, COO, carboxyl, epoxide_eph, or epoxide_ept"
        )
    if priority is not None:
        if (
            not isinstance(priority, (list, tuple))
            or not priority
            or any(not isinstance(g, str) or g not in allowed_prefs for g in priority)
        ):
            raise ValueError(
                "priority must be a non-empty list of valid group names"
            )
        if len(set(priority)) != len(priority):
            raise ValueError("priority must not contain duplicate group names")

    enabled_set = None
    if enabled_groups is not None:
        if (
            not isinstance(enabled_groups, (list, tuple))
            or not enabled_groups
            or any(not isinstance(g, str) or g not in allowed_prefs for g in enabled_groups)
        ):
            raise ValueError(
                "enabled_groups must be a non-empty list of valid group names"
            )
        enabled_set = set(enabled_groups)

    if enabled_set is not None and target_deprot is not None:
        allow_coo = ("COO" in enabled_set) or ("carboxyl" in enabled_set)
        allow_cooh = ("COOH" in enabled_set) or ("carboxyl" in enabled_set)
        if not allow_coo and (ph_val is not None or target_deprot > 0.0):
            raise ValueError(
                "enabled_groups excludes 'COO'; include 'COO' or 'carboxyl' to apply ph or deprotonation_fraction > 0"
            )
        if not allow_cooh and (ph_val is not None or target_deprot < 1.0):
            raise ValueError(
                "enabled_groups excludes 'COOH'; include 'COOH' or 'carboxyl' to apply ph or deprotonation_fraction < 1"
            )

    carboxyl_pref = 0.0
    if priority is not None:
        weights = dict.fromkeys(GROUPS, 0.0)
        n_prio = len(priority)
        for idx, group in enumerate(priority):
            w = float(2.0 ** (n_prio - 1 - idx))
            if group == "carboxyl":
                carboxyl_pref = w
            else:
                weights[group] = w
        if "carboxyl" in priority and target_deprot is None:
            raise ValueError("priority 'carboxyl' requires ph or deprotonation_fraction")
    elif preferences is not None:
        weights = dict(DEFAULT_WEIGHTS)
        for group, value in preferences.items():
            val = _positive_number(value, f"preference.{group}", allow_zero=True)
            if group == "carboxyl":
                carboxyl_pref = val
            else:
                weights[group] = val
        if "carboxyl" in preferences and target_deprot is None:
            raise ValueError("preference.carboxyl requires ph or deprotonation_fraction")
    elif enabled_set is not None:
        weights = dict.fromkeys(GROUPS, 0.0)
        for group in enabled_set:
            if group == "carboxyl":
                carboxyl_pref = 1.0
            else:
                weights[group] = 1.0
        if "carboxyl" in enabled_set and target_deprot is None:
            raise ValueError("enabled_groups 'carboxyl' requires ph or deprotonation_fraction")
    elif mode == "random":
        weights = dict.fromkeys(GROUPS, 0.0)
        for group in ("OH", "epoxide", "COOH", "COO"):
            weights[group] = 1.0
    elif target_deprot is not None:
        weights = dict.fromkeys(GROUPS, 0.0)
        weights["OH"] = 8.0
        weights["epoxide"] = 4.0
        weights["COOH"] = 2.0
        weights["COO"] = 1.0
    else:
        weights = dict(DEFAULT_WEIGHTS)

    if enabled_set is not None:
        if "carboxyl" in enabled_set and target_deprot is None:
            raise ValueError("enabled_groups 'carboxyl' requires ph or deprotonation_fraction")
        allowed_by_enabled = set(enabled_set)
        if "carboxyl" in enabled_set:
            allowed_by_enabled.update(("COOH", "COO"))
        has_carboxyl_enabled = bool(allowed_by_enabled & {"carboxyl", "COOH", "COO"})
        for group in GROUPS:
            if group not in allowed_by_enabled:
                weights[group] = 0.0
            elif (
                preferences is not None
                and group not in preferences
                and not (group in ("COOH", "COO") and "carboxyl" in preferences)
                and weights[group] == 0.0
            ):
                weights[group] = 1.0
        if not has_carboxyl_enabled:
            carboxyl_pref = 0.0

    if target_deprot is not None:
        carboxyl_total = weights["COOH"] + weights["COO"] + carboxyl_pref
        if carboxyl_total <= 0.0:
            raise ValueError(
                "ph/deprotonation_fraction requires positive COOH, COO, or carboxyl preference"
            )
        weights["COO"] = carboxyl_total * target_deprot
        weights["COOH"] = carboxyl_total * (1.0 - target_deprot)

    if not sum(weights.values()):
        raise ValueError("At least one group preference must be positive")

    if _scaffold is None:
        base = copy.deepcopy(recipe)
        base["groups"] = {}
        base.pop("composition", None)
        base.pop("c_over_o", None)
        base.pop("o_over_c", None)
        scaffold = build_sheet(base)
    else:
        scaffold = _scaffold

    carbons = {a["id"] for a in scaffold["atoms"] if a["element"] == "C"}
    degrees = Counter()
    for bond in scaffold["bonds"]:
        if bond["a"] in carbons and bond["b"] in carbons:
            degrees[bond["a"]] += 1
            degrees[bond["b"]] += 1
    n = len(carbons)
    basal_ids = {i for i in carbons if degrees[i] == 3}
    edge_ids = {i for i in carbons if degrees[i] == 2}
    basal = len(basal_ids)
    edge = len(edge_ids)
    ns = adjacency(scaffold)
    armchair_pairs = 0
    zigzag_edge_carbons = 0
    seen_edge: set[int] = set()
    for cid in sorted(edge_ids):
        if cid in seen_edge:
            continue
        comp = []
        stack = [cid]
        seen_edge.add(cid)
        while stack:
            curr = stack.pop()
            comp.append(curr)
            for nbr, _, _ in ns[curr]:
                if nbr in edge_ids and nbr not in seen_edge:
                    seen_edge.add(nbr)
                    stack.append(nbr)
        armchair_pairs += len(comp) // 2
        zigzag_edge_carbons += len(comp) % 2
    max_carboxyl = min(edge, armchair_pairs + zigzag_edge_carbons)

    det_rng = random.Random(0)
    max_epx = len(_max_matching(scaffold, basal_ids, det_rng))
    max_eph = len(_bipartite_matching(ns, sorted(basal_ids), edge_ids, det_rng))
    max_ept = len(_max_matching(scaffold, edge_ids, det_rng))
    capacities = {
        "OH": basal,
        "epoxide": max_epx,
        "carbonyl": edge,
        "COOH": max_carboxyl,
        "COO": max_carboxyl,
        "epoxide_eph": max_eph,
        "epoxide_ept": max_ept,
    }

    for g in GROUPS:
        if capacities[g] == 0:
            weights[g] = 0.0
    if not any(weights[g] > 0.0 for g in GROUPS):
        raise ValueError(
            "No integer composition satisfies the requested C/O, enabled groups and site capacities"
        )

    if mode == "random":
        comp_rng = random.Random(effective_seed)
        active = [g for g in GROUPS if weights[g] > 0.0 and capacities[g] > 0]
        has_user_weights = preferences is not None or priority is not None
        if target_deprot is not None and ("COOH" in active or "COO" in active):
            if has_user_weights:
                carb_w = (weights["COOH"] + weights["COO"]) * comp_rng.uniform(0.35, 2.2)
            else:
                carb_w = comp_rng.expovariate(1.0) + 0.05
            weights["COO"] = carb_w * target_deprot if "COO" in active else 0.0
            weights["COOH"] = carb_w * (1.0 - target_deprot) if "COOH" in active else 0.0
            for g in active:
                if g not in ("COOH", "COO"):
                    if has_user_weights:
                        weights[g] *= comp_rng.uniform(0.35, 2.2)
                    else:
                        weights[g] = comp_rng.expovariate(1.0) + 0.05
        else:
            for g in active:
                if has_user_weights:
                    weights[g] *= comp_rng.uniform(0.35, 2.2)
                else:
                    weights[g] = comp_rng.expovariate(1.0) + 0.05
        for g in GROUPS:
            if g not in active:
                weights[g] = 0.0

    if _jitter_seed is not None:
        j_rng = random.Random(_jitter_seed)
        if target_deprot is not None and (weights["COOH"] > 0.0 or weights["COO"] > 0.0):
            carb_w = (weights["COOH"] + weights["COO"]) * j_rng.uniform(0.65, 1.35)
            weights["COO"] = carb_w * target_deprot
            weights["COOH"] = carb_w * (1.0 - target_deprot)
            for g in GROUPS:
                if g not in ("COOH", "COO") and weights[g] > 0.0:
                    weights[g] *= j_rng.uniform(0.65, 1.35)
        else:
            for g in GROUPS:
                if weights[g] > 0.0:
                    weights[g] *= j_rng.uniform(0.65, 1.35)

    if _weight_multipliers is not None:
        for g in GROUPS:
            if g in _weight_multipliers and weights[g] > 0.0:
                weights[g] *= max(0.0, float(_weight_multipliers[g]))
        if not any(weights[g] > 0.0 for g in GROUPS):
            raise ValueError(
                "No integer composition satisfies the requested C/O, enabled groups and site capacities"
            )

    n_groups = len(GROUPS)
    oxygen = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0])
    extra_c = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
    fractions = np.array([weights[g] for g in GROUPS])
    fractions /= fractions.sum()
    rounding = float(np.max(oxygen[fractions > 0])) / 2
    middle = (low + high) / 2
    denominator = middle * float(oxygen @ fractions) - float(extra_c @ fractions)
    ideal = fractions * n / denominator if denominator > 0 else fractions * n
    # Integer counts (0..n_groups-1), then absolute deviations from preferred counts.
    rows = []
    lower = []
    upper = []

    def constrain(coeff, lo, hi):
        rows.append(np.r_[coeff, np.zeros(n_groups)])
        lower.append(lo)
        upper.append(hi)

    if interval:
        constrain(low * oxygen - extra_c, -np.inf, n)
        constrain(high * oxygen - extra_c, n, np.inf)
    else:
        constrain(low * oxygen - extra_c, n - low * rounding, n + low * rounding)
    constrain(np.array([1, 2, 0, 0, 0, 1, 0]), 0, basal)
    constrain(np.array([0, 0, 1, 1, 1, 1, 2]), 0, edge)
    constrain(np.array([0, 0, 0, 1, 1, 0, 0]), 0, max_carboxyl)
    constrain(oxygen, 1, np.inf)
    for j in range(n_groups):
        row = np.zeros(2 * n_groups)
        row[j] = 1
        row[n_groups + j] = -1
        rows.append(row)
        lower.append(-np.inf)
        upper.append(ideal[j])
        row = np.zeros(2 * n_groups)
        row[j] = -1
        row[n_groups + j] = -1
        rows.append(row)
        lower.append(-np.inf)
        upper.append(-ideal[j])
    maximum = np.array(
        [basal, max_epx, edge, max_carboxyl, max_carboxyl, max_eph, max_ept], dtype=float
    )
    maximum = np.where(np.array([weights[g] for g in GROUPS]) > 0, maximum, 0)
    solution = milp(
        c=np.r_[np.zeros(n_groups), 1 / (ideal + 1)],
        integrality=np.r_[np.ones(n_groups), np.zeros(n_groups)],
        bounds=Bounds(np.zeros(2 * n_groups), np.r_[maximum, np.full(n_groups, np.inf)]),
        constraints=LinearConstraint(np.array(rows), lower, upper),
        options={"time_limit": 10},
    )
    if not solution.success:
        if solution.status == 2:
            raise ValueError(
                "No integer composition satisfies the requested C/O, enabled groups and site capacities"
            )
        raise RuntimeError(f"Composition solver did not finish: {solution.message}")
    counts = np.rint(solution.x[:n_groups]).astype(int)
    final_c = n + int(extra_c @ counts)
    final_o = int(oxygen @ counts)
    actual = final_c / final_o
    actual_o_over_c = final_o / final_c
    if interval:
        assert low - 1e-8 <= actual <= high + 1e-8
    else:
        assert abs(final_o - final_c / low) <= rounding + 1e-8
    resolved = copy.deepcopy(recipe)
    resolved["groups"] = dict(zip(GROUPS, map(int, counts), strict=True))
    report = {
        "status": "INTEGER_FEASIBLE_NOT_PLACEMENT_VALIDATED",
        "requested_C_over_O": [low, high] if interval else low,
        "ratio_policy": "hard_interval" if interval else "half_group_oxygen_increment_rounding",
        "scalar_rounding_tolerance_oxygen_atoms": None if interval else rounding,
        "composition_mode": mode,
        "group_priority": list(priority) if priority else None,
        "group_preferences": weights,
        "preference_basis": "soft relative numbers of groups, not oxygen fractions",
        "default_preference_basis": "illustrative basal OH/epoxide 1:1; edge groups enabled only explicitly",
        "group_counts": resolved["groups"],
        "realized_group_fractions": dict(
            zip(GROUPS, map(float, counts / counts.sum()), strict=True)
        ),
        "scaffold_carbon_atoms": n,
        "basal_site_capacity": basal,
        "edge_site_capacity": edge,
        "max_disjoint_eph_capacity": max_eph,
        "max_disjoint_ept_capacity": max_ept,
        "predicted_carbon_atoms": final_c,
        "predicted_oxygen_atoms": final_o,
        "predicted_C_over_O": actual,
        "predicted_O_over_C": actual_o_over_c,
        "expected_formal_charge": -int(counts[4]),
        "experimental_group_assignment": False,
    }
    if requested_o_over_c is not None:
        report["requested_O_over_C"] = requested_o_over_c
    if mode == "random":
        report["composition_seed"] = effective_seed
    if target_deprot is not None:
        total_carboxyl = int(counts[3] + counts[4])
        report["protonation"] = {
            "mode": (
                "explicit_deprotonation_fraction"
                if deprotonation_fraction is not None
                else "henderson_hasselbalch"
            ),
            "ph": ph_val,
            "assumed_pka": None if deprotonation_fraction is not None else pka_val,
            "target_deprotonation_fraction": float(target_deprot),
            "realized_deprotonation_fraction": (
                float(counts[4] / total_carboxyl) if total_carboxyl > 0 else None
            ),
            "caveat": (
                "Single-pKa Henderson-Hasselbalch is an approximate model; "
                "graphene oxide surface carboxyl groups span an apparent pKa range of ~3-6."
            ),
        }
    return resolved, report
