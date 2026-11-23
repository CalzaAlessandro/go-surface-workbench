"""Graph-electronic (bipartite Kekule, Huckel pi-gap) and oxidation morphology analysis.

These diagnostics operate strictly on the constructed carbon/oxygen connectivity graph.
They screen closed-shell topological pi balance and cluster geometry; they do not
replace quantum-chemical calculations or dynamics.
"""

from __future__ import annotations

import math
from collections import deque
from typing import Any

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial import ConvexHull

from .model import _blossom_max_matching, _minimum_image_delta, adjacency


def _hex_cycles_in_graph(nodes: set[int], adj_map: dict[int, set[int]]) -> list[tuple[int, ...]]:
    """Find all distinct 6-membered cycles inside an induced degree-<=3 graph in O(N) time."""
    if len(nodes) < 6:
        return []
    cycles: dict[tuple[int, ...], tuple[int, ...]] = {}
    for u in sorted(nodes):
        nbrs_u = sorted(v for v in adj_map[u] if v in nodes and v > u)
        if len(nbrs_u) < 2:
            continue
        for i in range(len(nbrs_u)):
            a = nbrs_u[i]
            for j in range(i + 1, len(nbrs_u)):
                b = nbrs_u[j]
                for w1 in adj_map[a]:
                    if w1 not in nodes or w1 <= u or w1 == b:
                        continue
                    for w2 in adj_map[b]:
                        if w2 not in nodes or w2 <= u or w2 == a or w2 == w1:
                            continue
                        common = (adj_map[w1] & adj_map[w2]) & nodes
                        for x in common:
                            if x <= u or x in (a, b, w1, w2):
                                continue
                            # Exclude 6-cycles that also have a shorter chord inside the hex? Keep simple 6-cycles.
                            verts = (u, a, w1, x, w2, b)
                            if len(set(verts)) == 6:
                                key = tuple(sorted(verts))
                                cycles[key] = verts
    return list(cycles.values())


def _bipartite_coloring(
    nodes: list[int], adj_map: dict[int, set[int]], by_id: dict[int, dict[str, Any]]
) -> tuple[bool, dict[int, int]]:
    """Return (is_bipartite, {node_id: 0|1}), preferring stored graphene sublattice tags."""
    if nodes and all(by_id[n].get("metadata", {}).get("sublattice") in (0, 1) for n in nodes):
        color = {n: int(by_id[n]["metadata"]["sublattice"]) for n in nodes}
        node_set = set(nodes)
        for u in nodes:
            cu = color[u]
            for v in adj_map[u]:
                if v in node_set and color[v] == cu:
                    return False, color
        return True, color

    color = {}
    node_set = set(nodes)
    for start in nodes:
        if start in color:
            continue
        color[start] = 0
        dq = deque([start])
        while dq:
            u = dq.popleft()
            cu = color[u]
            for v in adj_map[u]:
                if v not in node_set:
                    continue
                if v not in color:
                    color[v] = 1 - cu
                    dq.append(v)
                elif color[v] == cu:
                    return False, color
    return True, color


def _identify_pi_system(
    model: dict[str, Any],
    by_id: dict[int, dict[str, Any]],
    bond_adj: dict[int, list[tuple[int, np.ndarray, float]]],
) -> tuple[set[int], list[int], list[tuple[int, int]]]:
    """Identify scaffold carbons, pure aromatic sp2 carbons, and carbonyl C=O pairs."""
    atoms = model["atoms"]
    scaffold_ids: set[int] = set()
    arom_ids: list[int] = []
    cab_pairs: list[tuple[int, int]] = []

    for atom in atoms:
        if atom["element"] != "C":
            continue
        aid = int(atom["id"])
        meta = atom.get("metadata", {})
        role = meta.get("role")
        group = meta.get("group")
        resname = atom.get("resname")
        nbrs = bond_adj[aid]
        deg = len(nbrs)

        is_sheet = role == "sheet_carbon"
        if role is None:
            c_nbrs = sum(1 for nbr, _, _ in nbrs if by_id[nbr]["element"] == "C")
            o_nbrs = sum(1 for nbr, _, _ in nbrs if by_id[nbr]["element"] == "O")
            is_sheet = c_nbrs >= 2 and not (o_nbrs >= 2 and deg == 3)

        if not is_sheet:
            continue
        scaffold_ids.add(aid)

        hybrid = meta.get("hybridization")
        is_sp3 = (hybrid == "sp3") or (hybrid is None and deg == 4)
        if is_sp3:
            continue

        o_nbrs = [int(nbr) for nbr, _, _ in nbrs if by_id[nbr]["element"] == "O"]
        is_carbonyl = (
            group == "carbonyl"
            or resname == "CAB"
            or (len(o_nbrs) == 1 and deg == 3)
        )
        if is_carbonyl and len(o_nbrs) == 1:
            cab_pairs.append((aid, o_nbrs[0]))
            continue

        if any(by_id[nbr]["element"] not in ("C", "H") for nbr, _, _ in nbrs):
            continue

        arom_ids.append(aid)

    return scaffold_ids, arom_ids, cab_pairs


def analyze_aromatic_electronic_structure(
    model: dict[str, Any],
    *,
    compute_huckel: bool = True,
    max_huckel_orbitals: int = 2500,
) -> dict[str, Any]:
    """Evaluate sp2 sublattice balance, Kekule matching, Huckel HOMO-LUMO gap, and intact hexagons."""
    by_id = {int(a["id"]): a for a in model["atoms"]}
    bond_adj = adjacency(model)
    c_adj: dict[int, set[int]] = {
        int(a["id"]): {
            int(nbr) for nbr, _, _ in bond_adj[int(a["id"])] if by_id[nbr]["element"] == "C"
        }
        for a in model["atoms"]
        if a["element"] == "C"
    }

    scaffold_ids, arom_ids, cab_pairs = _identify_pi_system(model, by_id, bond_adj)
    arom_set = set(arom_ids)
    arom_adj: dict[int, set[int]] = {
        cid: (c_adj[cid] & arom_set) for cid in arom_ids
    }

    bipartite, color = _bipartite_coloring(arom_ids, arom_adj, by_id)
    if bipartite:
        sub_a = [cid for cid in arom_ids if color[cid] == 0]
        sub_b = [cid for cid in arom_ids if color[cid] == 1]
        imbalance = abs(len(sub_a) - len(sub_b))
    else:
        sub_a = []
        sub_b = []
        imbalance = None

    # Exact maximum Kekule matching on pure aromatic sp2 network
    n_arom = len(arom_ids)
    if n_arom == 0:
        kekule_perfect = True
        kekule_matched = 0
    elif bipartite:
        if len(sub_a) != len(sub_b):
            kekule_perfect = False
            idx_a = {cid: i for i, cid in enumerate(sub_a)}
            idx_b = {cid: i for i, cid in enumerate(sub_b)}
            rows = [idx_a[cid] for cid in sub_a for nbr in arom_adj[cid] if nbr in idx_b]
            cols = [idx_b[nbr] for cid in sub_a for nbr in arom_adj[cid] if nbr in idx_b]
            mat = csr_matrix(([1] * len(rows), (rows, cols)), shape=(len(sub_a), len(sub_b)))
            perm = maximum_bipartite_matching(mat, perm_type="column")
            kekule_matched = int((perm >= 0).sum())
        else:
            idx_a = {cid: i for i, cid in enumerate(sub_a)}
            idx_b = {cid: i for i, cid in enumerate(sub_b)}
            rows = [idx_a[cid] for cid in sub_a for nbr in arom_adj[cid] if nbr in idx_b]
            cols = [idx_b[nbr] for cid in sub_a for nbr in arom_adj[cid] if nbr in idx_b]
            mat = csr_matrix(([1] * len(rows), (rows, cols)), shape=(len(sub_a), len(sub_b)))
            perm = maximum_bipartite_matching(mat, perm_type="column")
            kekule_matched = int((perm >= 0).sum())
            kekule_perfect = kekule_matched == len(sub_a)
    else:
        blossom = _blossom_max_matching(
            arom_ids, {cid: sorted(arom_adj[cid]) for cid in arom_ids}
        )
        kekule_matched = len(blossom) // 2
        kekule_perfect = len(blossom) == n_arom

    # Connected sp2 domains
    remaining = set(arom_ids)
    sp2_domains: list[int] = []
    while remaining:
        root = remaining.pop()
        stack = [root]
        sz = 0
        while stack:
            u = stack.pop()
            sz += 1
            for v in arom_adj[u]:
                if v in remaining:
                    remaining.remove(v)
                    stack.append(v)
        sp2_domains.append(sz)
    sp2_domains.sort(reverse=True)

    # Intact 6-membered aromatic rings and fused ring clusters (<= 6-ring Clar-like nuclei)
    total_hex_rings = _hex_cycles_in_graph(scaffold_ids, c_adj)
    intact_hex_rings = [cyc for cyc in total_hex_rings if all(v in arom_set for v in cyc)]
    intact_sets = [frozenset(cyc) for cyc in intact_hex_rings]
    bond_to_rings: dict[tuple[int, int], list[int]] = {}
    for ridx, cyc in enumerate(intact_hex_rings):
        for k in range(6):
            u, v = cyc[k], cyc[(k + 1) % 6]
            edge = (u, v) if u < v else (v, u)
            bond_to_rings.setdefault(edge, []).append(ridx)
    ring_adj: list[set[int]] = [set() for _ in range(len(intact_sets))]
    for ridxs in bond_to_rings.values():
        if len(ridxs) > 1:
            for a_idx in ridxs:
                for b_idx in ridxs:
                    if a_idx != b_idx:
                        ring_adj[a_idx].add(b_idx)
    seen_rings = [False] * len(intact_sets)
    ring_cluster_sizes: list[int] = []
    for ridx in range(len(intact_sets)):
        if seen_rings[ridx]:
            continue
        seen_rings[ridx] = True
        stack = [ridx]
        cnt = 0
        while stack:
            cur = stack.pop()
            cnt += 1
            for nxt in ring_adj[cur]:
                if not seen_rings[nxt]:
                    seen_rings[nxt] = True
                    stack.append(nxt)
        ring_cluster_sizes.append(cnt)
    ring_cluster_sizes.sort(reverse=True)

    # Per-component exact Huckel pi spectrum (pure aromatic sp2 + C=O carbonyls)
    huckel_gap_beta: float | None = None
    if compute_huckel:
        cab_c_ids = {c for c, _ in cab_pairs}
        cab_o_by_c = {c: o for c, o in cab_pairs}
        pi_c_set = arom_set | cab_c_ids
        pi_adj: dict[int, set[int]] = {cid: (c_adj[cid] & pi_c_set) for cid in pi_c_set}
        pi_remaining = set(pi_c_set)
        all_evals: list[float] = []
        too_large = False
        while pi_remaining:
            root = pi_remaining.pop()
            comp_c = [root]
            stack = [root]
            while stack:
                u = stack.pop()
                for v in pi_adj[u]:
                    if v in pi_remaining:
                        pi_remaining.remove(v)
                        comp_c.append(v)
                        stack.append(v)
            comp_cabs = [c for c in comp_c if c in cab_o_by_c]
            comp_nodes = comp_c + [cab_o_by_c[c] for c in comp_cabs]
            m = len(comp_nodes)
            if m > max_huckel_orbitals:
                too_large = True
                break
            if m == 1:
                # Single isolated sp2 radical site at alpha = 0
                c0 = comp_c[0]
                if c0 in cab_o_by_c:
                    all_evals.extend([-1.4, -1.0])
                else:
                    all_evals.append(0.0)
                continue
            idx = {nid: i for i, nid in enumerate(comp_nodes)}
            h_mat = np.zeros((m, m), dtype=float)
            for u in comp_c:
                iu = idx[u]
                for v in pi_adj[u]:
                    iv = idx[v]
                    h_mat[iu, iv] = -1.0
            for c_id in comp_cabs:
                ic, io = idx[c_id], idx[cab_o_by_c[c_id]]
                h_mat[ic, io] = -1.4
                h_mat[io, ic] = -1.4
                h_mat[io, io] = -1.0
            all_evals.extend(np.linalg.eigvalsh(h_mat).tolist())
        if not too_large:
            n_tot = len(all_evals)
            if n_tot >= 2 and n_tot % 2 == 0:
                all_evals.sort()
                huckel_gap_beta = float(max(0.0, all_evals[n_tot // 2] - all_evals[n_tot // 2 - 1]))
            elif n_tot >= 2:
                # Odd pi-electron count => SOMO half-filled => zero closed-shell HOMO-LUMO gap
                huckel_gap_beta = 0.0

    n_total_rings = len(total_hex_rings)
    return {
        "n_scaffold_carbons": len(scaffold_ids),
        "n_aromatic_sp2_carbons": n_arom,
        "n_carbonyl_pi_groups": len(cab_pairs),
        "bipartite": bipartite,
        "sublattice_a": len(sub_a) if bipartite else None,
        "sublattice_b": len(sub_b) if bipartite else None,
        "sublattice_imbalance": imbalance,
        "kekule_perfect": kekule_perfect,
        "kekule_matched_bonds": kekule_matched,
        "huckel_gap_beta": huckel_gap_beta,
        "largest_sp2_domain_carbons": sp2_domains[0] if sp2_domains else 0,
        "sp2_domain_count": len(sp2_domains),
        "total_hex_rings": n_total_rings,
        "intact_sp2_rings": len(intact_hex_rings),
        "intact_sp2_ring_fraction": (len(intact_hex_rings) / n_total_rings) if n_total_rings else 0.0,
        "largest_intact_sp2_ring_cluster": ring_cluster_sizes[0] if ring_cluster_sizes else 0,
        "intact_sp2_ring_clusters_ge6": sum(1 for s in ring_cluster_sizes if s >= 6),
        "intact_sp2_ring_cluster_sizes": ring_cluster_sizes,
    }


def compute_oxidation_morphology_metrics(
    model: dict[str, Any],
    *,
    compute_huckel: bool = True,
) -> dict[str, Any]:
    """Compute quantitative cluster, shape-irregularity, and sp2-domain metrics for GO models."""
    elec = analyze_aromatic_electronic_structure(model, compute_huckel=compute_huckel)
    by_id = {int(a["id"]): a for a in model["atoms"]}
    bond_adj = adjacency(model)
    c_adj: dict[int, set[int]] = {
        int(a["id"]): {
            int(nbr) for nbr, _, _ in bond_adj[int(a["id"])] if by_id[nbr]["element"] == "C"
        }
        for a in model["atoms"]
        if a["element"] == "C"
    }
    pbc = list(model.get("pbc", [False, False, False]))
    has_pbc = any(pbc)
    cell = np.asarray(model["cell"], dtype=float) if model.get("cell") is not None else None
    inv_cell = np.linalg.inv(cell) if (has_pbc and cell is not None) else None
    coords_by_id = {int(a["id"]): np.asarray(a["position"], dtype=float) for a in model["atoms"]}

    basal_groups = [
        g for g in model.get("metadata", {}).get("groups", []) if g.get("type") in ("OH", "epoxide")
    ]
    n_basal = len(basal_groups)
    if n_basal == 0:
        return {
            **elec,
            "n_basal_groups": 0,
            "mean_basal_group_nn_distance": None,
            "n_oxidized_clusters": 0,
            "singlet_basal_fraction": 0.0,
            "clustered_group_fraction_ge4": 0.0,
            "largest_oxidized_cluster_groups": 0,
            "mean_domain_asphericity": 0.0,
            "mean_domain_convex_circularity": 0.0,
            "mean_domain_hull_fill_fraction": 0.0,
            "oxidized_cluster_sizes": [],
        }

    # Midpoints for NN distance
    pts = []
    for g in basal_groups:
        sids = [int(s) for s in g["site_ids"]]
        if len(sids) == 1:
            pts.append(coords_by_id[sids[0]][:2])
        else:
            p0 = coords_by_id[sids[0]]
            d = coords_by_id[sids[1]] - p0
            if has_pbc and cell is not None and inv_cell is not None:
                d = _minimum_image_delta(d, pbc, cell, inv_cell)
            pts.append((p0 + 0.5 * d)[:2])
    pts_arr = np.asarray(pts, dtype=float)

    mean_nn: float | None = None
    if n_basal >= 2:
        nn_dists: list[float] = []
        for i in range(n_basal):
            delta = np.zeros((n_basal, 3), dtype=float)
            delta[:, :2] = pts_arr - pts_arr[i]
            if has_pbc and cell is not None and inv_cell is not None:
                delta = _minimum_image_delta(delta, pbc, cell, inv_cell)
            dists = np.linalg.norm(delta[:, :2], axis=1)
            dists[i] = float("inf")
            nn_dists.append(float(np.min(dists)))
        mean_nn = float(np.mean(nn_dists))

    # Build 2-hop C-scaffold adjacency among oxidized basal groups
    site_to_group: dict[int, int] = {}
    for gi, g in enumerate(basal_groups):
        for sid in g["site_ids"]:
            site_to_group[int(sid)] = gi

    group_adj: list[set[int]] = [set() for _ in range(n_basal)]
    for gi, g in enumerate(basal_groups):
        frontier_2hop: set[int] = set()
        for sid in g["site_ids"]:
            for n1 in c_adj.get(int(sid), ()):
                frontier_2hop.add(n1)
                frontier_2hop.update(c_adj.get(n1, ()))
        for c_id in frontier_2hop:
            gj = site_to_group.get(c_id)
            if gj is not None and gj != gi:
                group_adj[gi].add(gj)

    seen_g = [False] * n_basal
    clusters: list[list[int]] = []
    for gi in range(n_basal):
        if seen_g[gi]:
            continue
        seen_g[gi] = True
        comp = [gi]
        stack = [gi]
        while stack:
            cur = stack.pop()
            for nxt in group_adj[cur]:
                if not seen_g[nxt]:
                    seen_g[nxt] = True
                    comp.append(nxt)
                    stack.append(nxt)
        clusters.append(comp)
    clusters.sort(key=len, reverse=True)
    cluster_sizes = [len(c) for c in clusters]
    singlets = sum(sz for sz in cluster_sizes if sz == 1)
    ge4_groups = sum(sz for sz in cluster_sizes if sz >= 4)

    # Scaffold positions for hull fill fraction
    scaffold_ids = [
        int(a["id"])
        for a in model["atoms"]
        if a["element"] == "C" and a.get("metadata", {}).get("role", "sheet_carbon") == "sheet_carbon"
    ]
    scaffold_xy = np.asarray([coords_by_id[cid][:2] for cid in scaffold_ids], dtype=float)

    asphericities: list[float] = []
    circularities: list[float] = []
    fill_fractions: list[float] = []

    for comp in clusters:
        if len(comp) < 4:
            continue
        comp_sites = sorted({int(s) for gi in comp for s in basal_groups[gi]["site_ids"]})
        if has_pbc and cell is not None:
            # Unwrap via BFS over 2-hop carbon bridge around root carbon
            root = comp_sites[0]
            rel: dict[int, np.ndarray] = {root: coords_by_id[root][:2].copy()}
            comp_set = set(comp_sites)
            dq = deque([root])
            while len(rel) < len(comp_sites) and dq:
                u = dq.popleft()
                # Connect any unvisited cluster site within <= 2 hops first, then fallback nearest
                two_hop = {u}
                for n1 in c_adj.get(u, ()):
                    two_hop.add(n1)
                    two_hop.update(c_adj.get(n1, ()))
                for v in two_hop & comp_set:
                    if v not in rel:
                        d3 = _minimum_image_delta(
                            coords_by_id[v] - coords_by_id[u], pbc, cell, inv_cell
                        )
                        rel[v] = rel[u] + d3[:2]
                        dq.append(v)
                if not dq:
                    rem_sites = [s for s in comp_sites if s not in rel]
                    if rem_sites:
                        nxt = rem_sites[0]
                        d3 = _minimum_image_delta(
                            coords_by_id[nxt] - coords_by_id[root], pbc, cell, inv_cell
                        )
                        rel[nxt] = rel[root] + d3[:2]
                        dq.append(nxt)
            xy = np.asarray([rel[s] for s in comp_sites], dtype=float)
        else:
            xy = np.asarray([coords_by_id[s][:2] for s in comp_sites], dtype=float)

        centered = xy - np.mean(xy, axis=0)
        cov = np.cov(centered.T)
        evals = np.linalg.eigvalsh(cov)
        evals = np.maximum(evals, 0.0)
        s_eval = float(evals[0] + evals[1])
        asph = ((evals[1] - evals[0]) / s_eval) ** 2 if s_eval > 1e-12 else 0.0
        asphericities.append(float(asph))

        if len(xy) >= 3:
            try:
                hull = ConvexHull(xy)
                area = float(hull.volume)  # in 2D, scipy's volume is polygon area
                perim = float(hull.area)  # in 2D, scipy's area is polygon perimeter
                if perim > 1e-9:
                    circ = float((4.0 * math.pi * area) / (perim * perim))
                    circularities.append(min(1.0, circ))
                # Count scaffold carbons inside 2D convex hull to measure compact fill vs branched concavity
                hull_eq = hull.equations  # a*x + b*y + c <= 0 for interior
                if has_pbc and cell is not None:
                    # Check nearby scaffold atoms minimum-imaged around cluster centroid in raw frame
                    root_pos = coords_by_id[comp_sites[0]]
                    d_all = _minimum_image_delta(
                        np.asarray([coords_by_id[cid] for cid in scaffold_ids]) - root_pos,
                        pbc,
                        cell,
                        inv_cell,
                    )
                    s_xy = root_pos[:2] + d_all[:, :2]
                else:
                    s_xy = scaffold_xy
                inside = np.all(s_xy @ hull_eq[:, :-1].T + hull_eq[:, -1] <= 1e-7, axis=1)
                n_hull_c = int(np.sum(inside))
                if n_hull_c > 0:
                    fill_fractions.append(min(1.0, len(comp_sites) / float(n_hull_c)))
            except Exception:
                pass

    return {
        **elec,
        "n_basal_groups": n_basal,
        "mean_basal_group_nn_distance": mean_nn,
        "n_oxidized_clusters": len(clusters),
        "singlet_basal_fraction": singlets / float(n_basal),
        "clustered_group_fraction_ge4": ge4_groups / float(n_basal),
        "largest_oxidized_cluster_groups": cluster_sizes[0],
        "mean_domain_asphericity": float(np.mean(asphericities)) if asphericities else 0.0,
        "mean_domain_convex_circularity": float(np.mean(circularities)) if circularities else 0.0,
        "mean_domain_hull_fill_fraction": float(np.mean(fill_fractions)) if fill_fractions else 0.0,
        "oxidized_cluster_sizes": cluster_sizes,
    }


class IncrementalKekuleMatcher:
    """Maintain a maximum bipartite matching incrementally as carbon sites are oxidized."""

    def __init__(
        self,
        adj_map: dict[int, set[int]],
        color: dict[int, int],
        active_nodes: set[int] | list[int],
    ) -> None:
        self.adj: dict[int, set[int]] = {int(u): set(nbrs) for u, nbrs in adj_map.items()}
        self.color: dict[int, int] = dict(color)
        self.active: set[int] = set(active_nodes)
        self.match: dict[int, int] = {}

        sub_a = sorted(u for u in self.active if self.color[u] == 0)
        sub_b = sorted(u for u in self.active if self.color[u] == 1)
        if sub_a and sub_b:
            idx_a = {u: i for i, u in enumerate(sub_a)}
            idx_b = {v: j for j, v in enumerate(sub_b)}
            rows = [idx_a[u] for u in sub_a for v in self.adj[u] if v in idx_b]
            cols = [idx_b[v] for u in sub_a for v in self.adj[u] if v in idx_b]
            if rows:
                mat = csr_matrix(([1] * len(rows), (rows, cols)), shape=(len(sub_a), len(sub_b)))
                perm = maximum_bipartite_matching(mat, perm_type="column")
                for ai, bi in enumerate(perm):
                    if bi >= 0:
                        u = sub_a[int(ai)]
                        v = sub_b[int(bi)]
                        self.match[u] = v
                        self.match[v] = u

        self.free_a: set[int] = {u for u in sub_a if u not in self.match}
        self.free_b: set[int] = {v for v in sub_b if v not in self.match}

    def copy(self) -> IncrementalKekuleMatcher:
        clone = IncrementalKekuleMatcher.__new__(IncrementalKekuleMatcher)
        clone.adj = self.adj
        clone.color = self.color
        clone.active = set(self.active)
        clone.match = dict(self.match)
        clone.free_a = set(self.free_a)
        clone.free_b = set(self.free_b)
        return clone

    def _find_alt_path(
        self, start: int, targets: set[int], excluded: set[int]
    ) -> list[int] | None:
        if not targets:
            return None
        parent: dict[int, int | None] = {start: None}
        q = deque([start])
        while q:
            u = q.popleft()
            for y in self.adj[u]:
                if y in excluded or y not in self.active:
                    continue
                if y in targets:
                    path = [y, u]
                    cur = u
                    while parent[cur] is not None:
                        cur = parent[cur]
                        path.append(cur)
                    path.reverse()
                    return path
                mate = self.match.get(y)
                if (
                    mate is not None
                    and mate not in excluded
                    and mate in self.active
                    and mate not in parent
                ):
                    parent[y] = u
                    parent[mate] = y
                    q.append(mate)
        return None

    def _has_closing_partner(
        self,
        mate: int,
        excluded: set[int],
        eligible_closers: set[int] | None,
    ) -> bool:
        """Check whether removing an opposite-color neighbor later can re-close `mate`."""
        visited_same: set[int] = {mate}
        q = deque([mate])
        while q:
            u = q.popleft()
            for opp in self.adj[u]:
                if opp in excluded or opp not in self.active:
                    continue
                y = self.match.get(opp)
                if y is not None and y not in excluded and y in self.active:
                    if eligible_closers is None or y in eligible_closers:
                        return True
                    if y not in visited_same:
                        visited_same.add(y)
                        q.append(y)
        return False

    def _flip_path(self, path: list[int]) -> None:
        for i in range(0, len(path), 2):
            u, v = path[i], path[i + 1]
            self.match[u] = v
            self.match[v] = u

    def check_remove_pair(self, u: int, v: int) -> tuple[bool, Any]:
        """Check whether removing adjacent pair (u, v) keeps matching deficit non-increasing."""
        if u not in self.active or v not in self.active:
            return False, None
        cu, cv = self.color.get(u), self.color.get(v)
        if cu == cv:
            return False, None
        a, b = (u, v) if cu == 0 else (v, u)

        if self.match.get(a) == b:
            return True, ("matched", a, b)
        a_in = a in self.match
        b_in = b in self.match
        if not a_in and not b_in:
            return True, ("both_free", a, b)

        excl = {a, b}
        if a_in and not b_in:
            b_mate = self.match[a]
            path = self._find_alt_path(b_mate, self.free_a - {a}, excl)
            if path is not None:
                return True, ("half_a_aug", a, b, b_mate, path)
            return True, ("half_a_keep", a, b, b_mate)
        if b_in and not a_in:
            a_mate = self.match[b]
            path = self._find_alt_path(a_mate, self.free_b - {b}, excl)
            if path is not None:
                return True, ("half_b_aug", a, b, a_mate, path)
            return True, ("half_b_keep", a, b, a_mate)

        a_mate = self.match[b]
        b_mate = self.match[a]
        path = self._find_alt_path(a_mate, {b_mate}, excl)
        if path is not None:
            return True, ("bridge", a, b, a_mate, b_mate, path)
        if self.free_b and self.free_a:
            p1 = self._find_alt_path(a_mate, self.free_b, excl)
            if p1 is not None:
                p2 = self._find_alt_path(b_mate, self.free_a, excl | set(p1))
                if p2 is not None:
                    return True, ("double_aug", a, b, a_mate, b_mate, p1, p2)
        return False, None

    def commit_remove_pair(self, token: Any) -> None:
        kind = token[0]
        if kind == "matched":
            _, a, b = token
            self.active.discard(a)
            self.active.discard(b)
            self.match.pop(a, None)
            self.match.pop(b, None)
        elif kind == "both_free":
            _, a, b = token
            self.active.discard(a)
            self.active.discard(b)
            self.free_a.discard(a)
            self.free_b.discard(b)
        elif kind == "half_a_aug":
            _, a, b, b_mate, path = token
            self.active.discard(a)
            self.active.discard(b)
            self.free_b.discard(b)
            self.match.pop(a, None)
            self.match.pop(b_mate, None)
            self.free_a.discard(path[-1])
            self._flip_path(path)
        elif kind == "half_a_keep":
            _, a, b, b_mate = token
            self.active.discard(a)
            self.active.discard(b)
            self.free_b.discard(b)
            self.match.pop(a, None)
            self.match.pop(b_mate, None)
            self.free_b.add(b_mate)
        elif kind == "half_b_aug":
            _, a, b, a_mate, path = token
            self.active.discard(a)
            self.active.discard(b)
            self.free_a.discard(a)
            self.match.pop(b, None)
            self.match.pop(a_mate, None)
            self.free_b.discard(path[-1])
            self._flip_path(path)
        elif kind == "half_b_keep":
            _, a, b, a_mate = token
            self.active.discard(a)
            self.active.discard(b)
            self.free_a.discard(a)
            self.match.pop(b, None)
            self.match.pop(a_mate, None)
            self.free_a.add(a_mate)
        elif kind == "bridge":
            _, a, b, a_mate, b_mate, path = token
            self.active.discard(a)
            self.active.discard(b)
            self.match.pop(a, None)
            self.match.pop(b, None)
            self.match.pop(a_mate, None)
            self.match.pop(b_mate, None)
            self._flip_path(path)
        elif kind == "double_aug":
            _, a, b, a_mate, b_mate, p1, p2 = token
            self.active.discard(a)
            self.active.discard(b)
            self.match.pop(a, None)
            self.match.pop(b, None)
            self.match.pop(a_mate, None)
            self.match.pop(b_mate, None)
            self.free_b.discard(p1[-1])
            self.free_a.discard(p2[-1])
            self._flip_path(p1)
            self._flip_path(p2)

    def check_remove_single(
        self,
        x: int,
        *,
        allow_expose: bool = True,
        eligible_closers: set[int] | None = None,
    ) -> tuple[bool, Any]:
        """Check whether removing single site x preserves (or reversibly opens) Kekule match."""
        if x not in self.active:
            return False, None
        cx = self.color[x]
        if x not in self.match:
            return True, ("unmatched", x, cx)
        mate = self.match[x]
        self_free = self.free_a if cx == 0 else self.free_b
        if self_free:
            path = self._find_alt_path(mate, self_free, {x})
            if path is not None:
                return True, ("close", x, cx, mate, path)
        if allow_expose and self._has_closing_partner(mate, {x}, eligible_closers):
            return True, ("expose", x, cx, mate)
        return False, None

    def commit_remove_single(self, token: Any) -> None:
        kind = token[0]
        if kind == "unmatched":
            _, x, cx = token
            self.active.discard(x)
            (self.free_a if cx == 0 else self.free_b).discard(x)
        elif kind == "close":
            _, x, cx, mate, path = token
            self.active.discard(x)
            self.match.pop(x, None)
            self.match.pop(mate, None)
            (self.free_a if cx == 0 else self.free_b).discard(path[-1])
            self._flip_path(path)
        elif kind == "expose":
            _, x, cx, mate = token
            self.active.discard(x)
            self.match.pop(x, None)
            self.match.pop(mate, None)
            (self.free_b if cx == 0 else self.free_a).add(mate)
