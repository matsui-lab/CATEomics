"""Silence a deprecation of a pinned dependency that the user cannot act on.

econml 0.16.0 passes ``force_all_finite=`` to scikit-learn's ``check_array``.
scikit-learn renamed that argument to ``ensure_all_finite`` in 1.6 and warns
that the old spelling goes away in 1.8. Both versions are pinned by
``R/basilisk.R``, so the warning is raised on every fit that goes through an
econml estimator -- eleven times in a two-method comparison -- and nothing a
user of this package can do will change it. The removal cannot reach them
either, because the pin holds scikit-learn at 1.6.1.

The filter matches that one message and that one category, so any other
warning the Python side raises still reaches the console. When econml renames
the argument the filter stops matching anything and can be deleted.
"""

import warnings

_APPLIED = False


def apply():
    """Register the filter once per Python session."""
    global _APPLIED
    if _APPLIED:
        return
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message=r".*force_all_finite.*",
    )
    _APPLIED = True
