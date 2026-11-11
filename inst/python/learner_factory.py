"""
Learner factory: map CATEomics base-learner names to sklearn estimators.

Used by the econml-backed methods (the metalearners and the DR learner) so
that the `model_propensity`, `model_outcome`, and `model_final` arguments
accept the same string tokens as the R-side `.fit_base_learner()`.

Tokens supported (mirrors R/base_learners.R):

    "linear"       -> LinearRegression (regression)
    "logistic"     -> LogisticRegression (classification; L2)
    "lasso"        -> LassoCV (regression) / LogisticRegressionCV penalty='l1'
    "ridge"        -> RidgeCV (regression) / LogisticRegressionCV penalty='l2'
    "elastic_net"  -> ElasticNetCV (regression) / LogisticRegressionCV penalty='elasticnet'
    "rf"           -> RandomForest{Regressor,Classifier}
    "xgboost"      -> XGB{Regressor,Classifier}

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

# Hyperparameters the base-learner token itself fixes. Overriding one would
# turn e.g. "elastic_net" into a different learner than its name claims, so it
# is rejected rather than silently honoured.
_FROZEN = {
    "lasso": {"penalty", "l1_ratio", "l1_ratios"},
    "ridge": {"penalty", "l1_ratio", "l1_ratios"},
    "elastic_net": {"l1_ratio", "l1_ratios", "penalty"},
}

# Seed names, refused for every learner. The R side passes the one seed the
# whole run follows, and these are all in ``get_params()`` for the
# scikit-learn estimators, so without this they would pass the check below and
# take over from it for this one learner.
_RNG = {"random_state", "seed", "random_seed"}

# XGBoost carries two vocabularies for the same booster settings: the native
# names, which the booster reports in its own config, and the scikit-learn
# wrapper names, which are the only ones ``get_params()`` lists. The R side
# derives its accepted set from ``xgb.params()``, so it takes the native
# spelling, and the same key would be refused here. Upstream publishes no
# mapping, so the pairs are listed. ``seed`` is left out: it is refused as a
# seed name before this map is consulted.
_XGB_ALIAS = {
    # The R side takes `nrounds`, which is xgb.train's own argument for the
    # number of boosting rounds; the wrapper calls the same thing
    # `n_estimators`. Without this pair the number of rounds is the one
    # capacity setting that cannot be expressed on both routes.
    "nrounds": "n_estimators",
    "eta": "learning_rate",
    "min_split_loss": "gamma",
    "lambda": "reg_lambda",
    "alpha": "reg_alpha",
    "nthread": "n_jobs",
}


# Every penalty selection that forms a split is split the same number of ways,
# so one threshold describes all of them. A continuous target's ``"ridge"`` is
# a closed-form leave-one-out below and forms no split. The R side restates this as
# `.CF_CV_FOLDS` (it cannot be read across the process boundary) and a test
# reads this file to check the two have not drifted apart.
CV_FOLDS = 5


def _apply_params(est: Any, name: str, params: Optional[dict]) -> Any:
    """Apply user hyperparameters to a freshly built estimator.

    The accepted set is read off the estimator itself with ``get_params()``,
    so it is whatever this version of the upstream package exposes and no
    list is maintained here.

    The check is not left to ``set_params`` alone. Most scikit-learn
    estimators raise ValueError on an unknown keyword, but the XGBoost
    wrappers override ``set_params`` to forward unrecognised keywords to the
    booster, so a typo was accepted and only produced an xgboost-level
    "Parameters ... are not used" note on stderr at fit time -- the user's
    setting silently had no effect.
    """
    if not params:
        return est
    rng = sorted(set(params) & _RNG)
    if rng:
        raise ValueError(
            f"Hyperparameters cannot set {rng}; every draw follows the "
            f"random_state passed to cf_fit() / cf_compare()."
        )
    clash = sorted(set(params) & _FROZEN.get(name, set()))
    if clash:
        raise ValueError(
            f"Base learner '{name}' fixes {clash}; "
            f"choose a different learner instead of overriding it."
        )
    params = dict(params)
    if name == "xgboost":
        for native, sklearn_name in _XGB_ALIAS.items():
            if native in params:
                if sklearn_name in params:
                    raise ValueError(
                        f"'{native}' and '{sklearn_name}' are the same "
                        f"booster setting; pass one of them."
                    )
                params[sklearn_name] = params.pop(native)
    known = set(est.get_params())
    unknown = sorted(set(params) - known)
    if unknown:
        raise ValueError(
            f"Base learner '{name}' ({type(est).__name__}) has no parameter "
            f"{unknown}. It accepts: {sorted(known)}."
        )
    est.set_params(**params)
    return est


def make_drlearner_final_linear(random_state: Optional[int] = None) -> Any:
    """Return econml's StatsModelsLinearRegression for DRLearner's final stage.

    When ``model_final = "linear"`` is requested at the CATEomics level, the
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


def make_regressor(name: str, random_state: Optional[int] = None,
                   params: Optional[dict] = None) -> Any:
    """Build a regressor for ``name`` and apply ``params`` (see below)."""
    return _apply_params(_build_regressor(name, random_state), name, params)


def _build_regressor(name: str, random_state: Optional[int] = None) -> Any:
    """Return an unfitted sklearn regressor corresponding to ``name``.

    Mirrors the CATEomics R-side base-learner vocabulary from
    ``R/base_learners.R`` so that the ``model_*`` arguments of Python-
    backed methods accept the same tokens.

    Parameters
    ----------
    name : str
        Base-learner token. One of ``"linear"`` (-> LinearRegression),
        ``"lasso"`` (-> LassoCV, CV_FOLDS folds), ``"ridge"`` (-> RidgeCV,
        closed-form leave-one-out, no
        random_state), ``"elastic_net"`` (-> ElasticNetCV,
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
        return LassoCV(cv=CV_FOLDS, random_state=random_state)
    if name == "ridge":
        # Left at scikit-learn's default, which for this estimator alone is the
        # closed-form leave-one-out rather than a k-fold split: ridge has an
        # identity that gives the leave-one-out error from the hat matrix, and
        # the other two penalties do not. It forms no partition, so no number
        # of rows can leave a fold short, and R has no equivalent -- cv.glmnet
        # offers leave-one-out only as nfolds = n, by brute force. Overriding
        # it to k folds would discard a deliberate upstream choice to make one
        # number look uniform.
        return RidgeCV()  # RidgeCV has no random_state
    if name == "elastic_net":
        # max_iter = 5000 (raised from sklearn default 1000) to let
        # coordinate descent converge in n<<p settings (e.g. n=200, p=1001
        # after the `[X | T | X*T]` expansion that cf_fit_slearner and the
        # DR-learner's outcome nuisance use; the t- and x-learner fit econml's
        # per-arm models instead). tol stays at sklearn's default 1e-4.
        # Mirrors the LogisticRegressionCV `max_iter=5000` below.
        return ElasticNetCV(
            cv=CV_FOLDS, l1_ratio=0.5, max_iter=5000, random_state=random_state,
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
                "CATEomics provisions it in the environment basilisk builds on "
                "first use, so this means that environment is not the one "
                "running."
            ) from e
        return xgb.XGBRegressor(
            n_estimators=100,
            max_depth=4,
            learning_rate=0.1,
            min_child_weight=10,
            colsample_bytree=0.5,
            random_state=random_state,
            verbosity=0,
        )
    # "logistic" is not meaningful for regression -> fall back to linear.
    # Callers that pass "logistic" for an outcome model get a sensible default.
    return LinearRegression()


def make_classifier(name: str, random_state: Optional[int] = None,
                    params: Optional[dict] = None) -> Any:
    """Build a classifier for ``name`` and apply ``params`` (see below)."""
    return _apply_params(_build_classifier(name, random_state), name, params)


def _build_classifier(name: str, random_state: Optional[int] = None) -> Any:
    """Return an unfitted sklearn classifier corresponding to ``name``.

    Mirrors the CATEomics R-side base-learner vocabulary from
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
            penalty="l1", solver="saga", cv=CV_FOLDS, max_iter=5000,
            random_state=random_state,
        )
    if name == "ridge":
        return LogisticRegressionCV(
            penalty="l2", cv=CV_FOLDS, max_iter=5000, random_state=random_state,
        )
    if name == "elastic_net":
        return LogisticRegressionCV(
            penalty="elasticnet", solver="saga", l1_ratios=[0.5],
            cv=CV_FOLDS, max_iter=5000, random_state=random_state,
        )
    if name == "rf":
        return RandomForestClassifier(
            n_estimators=100, random_state=random_state,
        )
    if name == "xgboost":
        # Mirrors _build_regressor: the token is documented as one of the
        # seven base learners, and the R-side .fit_base_learner() supports it
        # for classification too ("binary:logistic"). Without this branch a
        # binary outcome with model_outcome = "xgboost" failed here, while
        # the same call on a continuous outcome worked.
        try:
            import xgboost as xgb
        except ImportError as e:
            raise ImportError(
                "Python xgboost is required for model_outcome='xgboost'. "
                "CATEomics provisions it in the environment basilisk builds on "
                "first use, so this means that environment is not the one "
                "running."
            ) from e
        return xgb.XGBClassifier(
            n_estimators=100,
            max_depth=4,
            learning_rate=0.1,
            min_child_weight=10,
            colsample_bytree=0.5,
            random_state=random_state,
            verbosity=0,
        )
    raise ValueError(f"Unknown base learner '{name}'")


class ProbaAsRegressor(BaseEstimator, RegressorMixin):
    """Wrap a binary classifier so it behaves like a regressor of ``P(Y=1|X)``.

    Binary-outcome support for estimators that only accept regressors (the
    econml metalearners). ``fit`` trains the wrapped
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
                       random_state: Optional[int] = None,
                       params: Optional[dict] = None) -> Any:
    """Outcome model factory.

    Returns a plain regressor, or ``ProbaAsRegressor(classifier)`` when
    ``binary`` is True so that binary outcomes are modelled as ``P(Y=1|X)``
    while keeping the regressor interface econml expects.
    """
    if binary:
        return ProbaAsRegressor(
            make_classifier(name, random_state=random_state, params=params))
    return make_regressor(name, random_state=random_state, params=params)
