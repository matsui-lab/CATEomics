"""
Learner factory: map cfomics base-learner names to sklearn estimators.

Used by Python-backed methods (AIPW / DoubleML, metalearners) so that the
`model_propensity`, `model_outcome`, and `model_final` arguments accept the
same string tokens as the R-side `.fit_base_learner()` / `.sl_lib_from_base_learner()`.

Tokens supported (mirrors R/base_learners.R):

    "linear"       -> LinearRegression (regression)
    "logistic"     -> LogisticRegression (classification; L2)
    "lasso"        -> LassoCV (regression) / LogisticRegressionCV penalty='l1'
    "ridge"        -> RidgeCV (regression) / LogisticRegressionCV penalty='l2'
    "elastic_net"  -> ElasticNetCV (regression) / LogisticRegressionCV penalty='elasticnet'
    "rf"           -> RandomForest{Regressor,Classifier}

Unknown tokens raise ValueError.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import (
    ElasticNetCV,
    LassoCV,
    LinearRegression,
    LogisticRegression,
    LogisticRegressionCV,
    RidgeCV,
)

_VALID = {"linear", "logistic", "lasso", "ridge", "elastic_net", "rf", "xgboost"}


def make_drlearner_final_linear(random_state: Optional[int] = None) -> Any:
    """Return econml's StatsModelsLinearRegression for DRLearner's final stage.

    When ``model_final = "linear"`` is requested at the cfomics level, the
    DR-Learner final stage uses econml's ``StatsModelsLinearRegression`` so
    that ``ate_interval()`` /
    ``effect_interval()`` return asymptotic Wald CIs based on the OLS
    sandwich variance (instead of failing silently as plain
    ``sklearn.linear_model.LinearRegression`` does).

    Parameters
    ----------
    random_state : int or None
        Forwarded for API symmetry; ``StatsModelsLinearRegression`` has no
        stochastic component but the argument is accepted for the factory
        signature.

    Returns
    -------
    StatsModelsLinearRegression
        Unfitted final-stage regressor with asymptotic-Wald inference.
    """
    # econml >= 0.14 ships StatsModelsLinearRegression under
    # ``econml.sklearn_extensions.linear_model``. Import lazily so the
    # factory module itself stays importable when econml is missing.
    from econml.sklearn_extensions.linear_model import (
        StatsModelsLinearRegression,
    )
    return StatsModelsLinearRegression()


def _check_token(name: str) -> None:
    if name not in _VALID:
        raise ValueError(
            f"Unknown base learner '{name}'. "
            f"Valid tokens: {sorted(_VALID)}"
        )


def make_regressor(name: str, random_state: Optional[int] = None) -> Any:
    """Return an unfitted sklearn regressor corresponding to ``name``.

    Mirrors the cfomics R-side base-learner vocabulary from
    ``R/base_learners.R`` so that the ``model_*`` arguments of Python-
    backed methods accept the same tokens.

    Parameters
    ----------
    name : str
        Base-learner token. One of ``"linear"`` (-> LinearRegression),
        ``"lasso"`` (-> LassoCV, cv=5), ``"ridge"`` (-> RidgeCV, no
        random_state), ``"elastic_net"`` (-> ElasticNetCV, cv=5,
        l1_ratio=0.5), ``"rf"`` (-> RandomForestRegressor,
        n_estimators=100). The ``"logistic"`` token falls back to
        LinearRegression since it is not meaningful for regression.
    random_state : int or None
        Forwarded to estimators that accept it. ``RidgeCV`` ignores it.

    Returns
    -------
    sklearn estimator
        An unfitted regressor. Caller is responsible for ``.fit()``.

    Raises
    ------
    ValueError
        If ``name`` is not in the supported token set.

    Notes
    -----
    Numerical outputs differ from the R-side ``.fit_base_learner()`` for
    the same token (different CV fold assignment, default regularization
    paths, and solver conventions).
    """
    _check_token(name)
    if name == "linear":
        return LinearRegression()
    if name == "lasso":
        return LassoCV(cv=5, random_state=random_state)
    if name == "ridge":
        return RidgeCV()  # RidgeCV has no random_state
    if name == "elastic_net":
        # max_iter = 5000 (raised from sklearn default 1000) to let
        # coordinate descent converge in n<<p settings (e.g. n=200, p=1001
        # after `[X | T | X*T]` expansion in cf_fit_slearner / tlearner /
        # xlearner). tol stays at sklearn default 1e-4 (precision preserved).
        # Mirrors the LogisticRegressionCV `max_iter=5000` below.
        return ElasticNetCV(
            cv=5, l1_ratio=0.5, max_iter=5000, random_state=random_state,
        )
    if name == "rf":
        return RandomForestRegressor(
            n_estimators=100,
            random_state=random_state,
        )
    if name == "xgboost":
        # Lazy import so factory module stays importable when xgboost
        # is missing (raises clear ImportError on construction).
        try:
            import xgboost as xgb
        except ImportError as e:
            raise ImportError(
                "Python xgboost is required for model_outcome='xgboost'. "
                "Install via `pip install xgboost` or "
                "`cf_install_python_env('unified')` (cfomics unified env "
                "includes xgboost)."
            ) from e
        return xgb.XGBRegressor(
            n_estimators=100,
            random_state=random_state,
            verbosity=0,
        )
    # "logistic" is not meaningful for regression -> fall back to linear.
    # Callers that pass "logistic" for an outcome model get a sensible default.
    return LinearRegression()


def make_classifier(name: str, random_state: Optional[int] = None) -> Any:
    """Return an unfitted sklearn classifier corresponding to ``name``.

    Mirrors the cfomics R-side base-learner vocabulary from
    ``R/base_learners.R`` so that the ``model_*`` arguments of Python-
    backed methods accept the same tokens.

    Parameters
    ----------
    name : str
        Base-learner token. One of ``"linear"`` / ``"logistic"`` (both
        alias to LogisticRegression, max_iter=1000), ``"lasso"`` (->
        LogisticRegressionCV, penalty='l1', solver='saga'), ``"ridge"``
        (-> LogisticRegressionCV, penalty='l2'), ``"elastic_net"`` (->
        LogisticRegressionCV, penalty='elasticnet', solver='saga',
        l1_ratios=[0.5]), ``"rf"`` (-> RandomForestClassifier,
        n_estimators=100).
    random_state : int or None
        Forwarded to all estimators.

    Returns
    -------
    sklearn estimator
        An unfitted binary classifier. Caller is responsible for ``.fit()``.

    Raises
    ------
    ValueError
        If ``name`` is not in the supported token set.

    Notes
    -----
    Numerical outputs differ from the R-side ``.fit_base_learner()`` for
    the same token (different CV fold assignment, default regularization
    paths, and solver conventions).
    """
    _check_token(name)
    if name in ("linear", "logistic"):
        return LogisticRegression(
            max_iter=1000, random_state=random_state,
        )
    if name == "lasso":
        return LogisticRegressionCV(
            penalty="l1", solver="saga", cv=5, max_iter=5000,
            random_state=random_state,
        )
    if name == "ridge":
        return LogisticRegressionCV(
            penalty="l2", cv=5, max_iter=5000, random_state=random_state,
        )
    if name == "elastic_net":
        return LogisticRegressionCV(
            penalty="elasticnet", solver="saga", l1_ratios=[0.5],
            cv=5, max_iter=5000, random_state=random_state,
        )
    if name == "rf":
        return RandomForestClassifier(
            n_estimators=100, random_state=random_state,
        )
    raise ValueError(f"Unknown base learner '{name}'")


class ProbaAsRegressor(BaseEstimator, RegressorMixin):
    """Wrap a binary classifier so it behaves like a regressor of ``P(Y=1|X)``.

    Binary-outcome support for estimators that only accept regressors (econml
    metalearners, DoubleML ``ml_g``). ``fit`` trains the wrapped
    classifier; ``predict`` returns the class-1 probability, so downstream code
    that expects ``E[Y|X]`` transparently receives the risk model and the CATE
    contract ``mu1 - mu0`` becomes a risk difference.
    """

    def __init__(self, classifier: Any):
        self.classifier = classifier

    def fit(self, X, y, **kwargs):
        self.classifier_ = clone(self.classifier)
        self.classifier_.fit(X, np.asarray(y).ravel(), **kwargs)
        self.classes_ = getattr(self.classifier_, "classes_", np.array([0.0, 1.0]))
        return self

    def predict(self, X):
        proba = np.asarray(self.classifier_.predict_proba(X))
        if proba.ndim == 2 and proba.shape[1] == 1:
            # Degenerate single-class fold: probability collapses to the label.
            return np.full(proba.shape[0], float(self.classes_[0]))
        # Locate the column for class == 1 (falls back to the last column).
        classes = np.asarray(self.classes_).ravel()
        idx = int(np.where(classes == 1)[0][0]) if np.any(classes == 1) else proba.shape[1] - 1
        return proba[:, idx].ravel()

    def get_params(self, deep: bool = True):
        return {"classifier": self.classifier}

    def set_params(self, **params):
        if "classifier" in params:
            self.classifier = params["classifier"]
        return self


def make_outcome_model(name: str, binary: bool = False,
                       random_state: Optional[int] = None) -> Any:
    """Outcome model factory.

    Returns a plain regressor, or ``ProbaAsRegressor(classifier)`` when
    ``binary`` is True so that binary outcomes are modelled as ``P(Y=1|X)``
    while keeping the regressor interface expected by econml / DoubleML.
    """
    if binary:
        return ProbaAsRegressor(make_classifier(name, random_state=random_state))
    return make_regressor(name, random_state=random_state)
