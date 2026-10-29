"""Audited CPPTRAJ align→rms input migration; never rewrites trajectory data."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import tempfile
from pathlib import Path


def _mask(token):
    return token == "MASK" or token.startswith((":", "@", "!", "*", "("))


def rewrite_alignment(text, default_mask=None):
    """Keep the first-frame fit mask, but move all atoms using RMS.

    Blank legacy fits used the default all-atom mask. An explicit autoimage
    anchor supplies the requested GO mask for the new fit, without claiming
    those legacy trajectories suffered a masked-only transformation.
    """
    inferred = default_mask
    for line in text.splitlines():
        tokens = shlex.split(line, comments=True)
        if tokens and tokens[0].lower() == "autoimage":
            candidates = (
                tokens[2:] if len(tokens) > 1 and tokens[1].lower() == "anchor" else tokens[1:]
            )
            if candidates and _mask(candidates[0]):
                inferred = inferred or candidates[0]
    changes = []
    result = []
    for line in text.splitlines(keepends=True):
        match = re.match(r"^([ \t]*)align\s+first\b([^\r\n]*)(\r?\n)?$", line, re.I)
        if not match:
            result.append(line)
            continue
        indent, tail, ending = match.groups()
        body, separator, comment = tail.partition("#")
        tokens = shlex.split(body)
        fit = []
        mass = False
        move = None
        index = 0
        while index < len(tokens):
            token = tokens[index]
            if token.lower() == "mass":
                mass = True
            elif token.lower() == "move":
                if index + 1 >= len(tokens) or not _mask(tokens[index + 1]):
                    raise ValueError("align move requires an explicit atom mask")
                move = tokens[index + 1]
                index += 1
            elif _mask(token) and len(fit) < 2:
                fit.append(token)
            else:
                raise ValueError(f"Unsupported alignment argument; review manually: {token}")
            index += 1
        previous_fit = fit[0] if fit else "all_atoms_default"
        previous_move = move or previous_fit
        if not fit:
            if not inferred or not _mask(inferred):
                raise ValueError(
                    "Blank align first requires a known GO mask or explicit autoimage anchor"
                )
            fit = [inferred]
        risk = previous_move not in ("all_atoms_default", "*", ":*")
        replacement = indent + "rms first " + " ".join(fit) + (" mass" if mass else "")
        if separator:
            replacement += " #" + comment
        replacement += ending or ""
        changes.append(
            {
                "previous_command": line.strip(),
                "new_command": replacement.strip(),
                "previous_fit": previous_fit,
                "previous_move": previous_move,
                "fit_mask": fit[0],
                "relative_geometry_at_risk": risk,
            }
        )
        result.append(replacement)
    return {"text": "".join(result), "changes": changes}


def patch_inputs(paths, backup_dir, *, apply=False, default_mask=None):
    """Plan or atomically replace only listed inputs, with immutable backups.

    Read-only reference snapshots and historical archives should not be listed.
    No output trajectories/energies are regenerated or marked as recomputed.
    """
    backup_dir = Path(backup_dir).resolve()
    if apply and backup_dir.exists():
        raise FileExistsError(f"Choose a fresh backup directory: {backup_dir}")
    prepared = []
    seen = set()
    for value in paths:
        path = Path(value)
        if path.is_symlink():
            raise ValueError(f"Refusing a symlink input: {path}")
        path = path.resolve(strict=True)
        if path in seen:
            continue
        seen.add(path)
        old = path.read_bytes()
        rewritten = rewrite_alignment(old.decode(), default_mask)
        if not rewritten["changes"]:
            continue
        new = rewritten["text"].encode()
        digest = hashlib.sha256(old).hexdigest()
        record = {
            "path": str(path),
            "backup": f"backups/{digest}.in",
            "old_sha256": digest,
            "new_sha256": hashlib.sha256(new).hexdigest(),
            "changes": rewritten["changes"],
        }
        prepared.append((path, old, new, record))
    report = {
        "status": "PLANNED",
        "changed_files": len(prepared),
        "files": [entry[3] for entry in prepared],
        "trajectories_recomputed": False,
        "policy": "rms fits the GO mask and transforms the complete frame; original inputs backed up",
    }
    if not apply:
        return report
    backup_dir.mkdir(parents=True)
    (backup_dir / "backups").mkdir()
    for path, old, _, record in prepared:
        if path.read_bytes() != old:
            raise ValueError(f"Input changed after inventory: {path}")
        backup = backup_dir / record["backup"]
        if not backup.exists():
            backup.write_bytes(old)
    (backup_dir / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    for path, old, new, _record in prepared:
        if path.read_bytes() != old:
            raise ValueError(f"Input changed before patch: {path}")
        with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=".rms_alignment_", delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(new)
        os.chmod(temporary, path.stat().st_mode & 0o777)
        os.replace(temporary, path)
    report["status"] = "APPLIED"
    (backup_dir / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
