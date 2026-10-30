"""Explicit finite-model QM/RESP protocols; no implicit method substitution."""
from __future__ import annotations

import re

HF_METHOD = 'HF/6-31G*'
DFT_METHOD = 'wB97XD/6-311++G(d,p)'
RESP_METHODS = (HF_METHOD, DFT_METHOD)


def require_method(method):
    if method not in RESP_METHODS:
        raise ValueError(f'Unsupported QM/RESP method: {method!r}; choose one of {RESP_METHODS}')
    return method


def route_matches_method(route, method):
    """Match the complete method/basis token, including Gaussian print aliases.

    Whitespace wrapping inside a Gaussian route is insignificant. Require one
    method/basis token so a second, contradictory method cannot hide beside it.
    """
    require_method(method)
    names = (('HF/6-31G*', 'RHF/6-31G*', 'HF/6-31G(d)', 'RHF/6-31G(d)')
             if method == HF_METHOD else
             (DFT_METHOD, 'R'+DFT_METHOD, 'U'+DFT_METHOD))
    # Gaussian may wrap inside a keyword. Keep the boundary after the complete
    # method token so unsupported basis suffixes cannot be silently accepted.
    pattern = '|'.join(r'\s*'.join(re.escape(c) for c in name) for name in names)
    matched = re.search(r'(?<![A-Za-z0-9])(?:'+pattern+r')(?=\s|$)', route, re.I)
    without_iop = re.sub(r'\bIOp\s*\([^)]*\)', '', route, flags=re.I)
    return bool(matched) and without_iop.count('/') == 1


def energy_matches_method(label, method):
    require_method(method)
    pattern = r'(?:RHF|UHF|ROHF)' if method == HF_METHOD else r'(?:R|U)?wB97XD'
    return re.fullmatch(pattern, label, flags=re.I) is not None
