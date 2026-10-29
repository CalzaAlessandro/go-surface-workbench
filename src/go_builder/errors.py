"""Structured diagnostics when geometric search cannot achieve a requested loading."""

from __future__ import annotations


class PlacementError(ValueError):
    """A valid request could not be fulfilled; no partial model is returned.

    ``achieved`` describes a tried placement, not a proof of the maximum
    physically attainable loading. ``capacity`` contains only known bounds.
    """

    def __init__(
        self, message, *, requested, achieved=None, capacity=None, attempts=None, suggestions=None
    ):
        super().__init__(message)
        self.requested = dict(requested)
        self.achieved = dict(achieved) if achieved is not None else None
        self.capacity = dict(capacity) if capacity is not None else None
        self.attempts = attempts
        self.suggestions = (
            list(suggestions)
            if suggestions
            else [
                "Increase the sheet size or reduce the requested loading.",
                "Try a different seed or fragment conformation; geometric search is finite.",
            ]
        )

    def as_report(self):
        report = {"status": "INCOMPLETE", "message": str(self), "requested": self.requested}
        for key in ("achieved", "capacity", "attempts", "suggestions"):
            value = getattr(self, key)
            if value is not None:
                report[key] = value
        return report
