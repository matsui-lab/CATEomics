"""Outcome model that sees the treatment interacted with the covariates.

econml's DRLearner fits its outcome nuisance on ``hstack([XW, T])`` and then
predicts it at ``T = 0`` and ``T = 1``. With a linear or penalised base
learner that design can only shift the intercept, so the two counterfactual
predictions differ by a constant and the nuisance carries no heterogeneity at
all. Expanding the design to ``[XW, T, XW * T]`` lets the treated and control
surfaces have different slopes, which is what the arm-wise nuisance of the
canonical DR-learner achieves. Without a penalty the two are exactly equal;
under one they differ only in whether the shrinkage is shared across arms.

CATEomics' S-learner fits the same expansion, so both go through this wrapper.

The wrapper takes the design econml hands it -- covariates first, the
treatment column(s) last -- and appends the products. ``predict_proba`` is
exposed only when the inner estimator has it, because econml tests for that
attribute to decide whether the outcome is discrete.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from sklearn.base import BaseEstimator, clone
from sklearn.utils.metaestimators import available_if


def _inner_has(attr: str):
    def check(self) -> bool:
        return hasattr(self.base_estimator, attr)
    return check


class InteractionAugmentedModel(BaseEstimator):
    """Fit ``base_estimator`` on ``[X, T, X * T]`` given ``[X, T]``.

    Parameters
    ----------
    base_estimator : sklearn estimator
        The learner that sees the expanded design.
    n_treatment_cols : int
        How many trailing columns of the input are the treatment. One for a
        binary treatment, which is what CATEomics supports.
    """

    def __init__(self, base_estimator: Any = None, n_treatment_cols: int = 1):
        self.base_estimator = base_estimator
        self.n_treatment_cols = n_treatment_cols

    def _expand(self, X) -> np.ndarray:
        X = np.asarray(X, dtype=np.float64)
        k = int(self.n_treatment_cols)
        if k < 1:
            raise ValueError("n_treatment_cols must be at least 1.")
        if X.shape[1] <= k:
            raise ValueError(
                f"Design has {X.shape[1]} column(s) but {k} of them were "
                "declared to be the treatment; there is nothing to interact."
            )
        covariates = X[:, :-k]
        treatment = X[:, -k:]
        products = np.hstack(
            [covariates * treatment[:, [j]] for j in range(k)]
        )
        return np.hstack([covariates, treatment, products])

    def fit(self, X, y, **kwargs):
        self.estimator_ = clone(self.base_estimator)
        self.estimator_.fit(self._expand(X), y, **kwargs)
        # Surface the attributes callers look for on a fitted estimator.
        if hasattr(self.estimator_, "classes_"):
            self.classes_ = self.estimator_.classes_
        return self

    def predict(self, X):
        return self.estimator_.predict(self._expand(X))

    @available_if(_inner_has("predict_proba"))
    def predict_proba(self, X):
        return self.estimator_.predict_proba(self._expand(X))

    @available_if(_inner_has("score"))
    def score(self, X, y, **kwargs):
        return self.estimator_.score(self._expand(X), y, **kwargs)


def make_interaction_augmented(base_estimator: Any,
                               n_treatment_cols: int = 1
                               ) -> InteractionAugmentedModel:
    """Convenience constructor used from the R side."""
    return InteractionAugmentedModel(base_estimator=base_estimator,
                                     n_treatment_cols=n_treatment_cols)
