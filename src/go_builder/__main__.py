"""Module entrypoint (`python -m go_builder`) for installed package/wheel execution."""

from __future__ import annotations

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
