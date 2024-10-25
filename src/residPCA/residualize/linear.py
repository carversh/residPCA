from typing import Literal, Union

import numpy as np
import scipy.linalg
from numpy import ndarray

from sklearn.preprocessing import StandardScaler
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline


def _test_onehot(C: ndarray) -> bool:
    """
    Test if the input matrix is one-hot encoded.

    Parameters
    ----------
    C : np.ndarray
        Covariates matrix.

    Returns
    -------
    bool
        Whether the input matrix is one-hot encoded.
    """
    return np.all(
        np.sum(C, axis=1) == 1
    ) and np.all(
        np.logical_or(C == 0, C == 1)
    )


def _linear_regression_effectsizes(
            X: ndarray,
            C: ndarray,
            method: Literal['onehot', 'lstsq', 'pinv'] = 'lstsq'
        ) -> Union[ndarray, ndarray]:
    """
    Perform linear regression and return the effect sizes.

    Parameters
    ----------
    X : np.ndarray
        Data matrix.
    C : np.ndarray
        Covariates matrix.
    method : {'onehot', 'lstsq', 'pinv'}
        Method to use for the regression. 'onehot' can be used when the
        covariates are one-hot encoded and is more efficient for large number
        of cells, 'lstsq' is used for the least squares solution (stable), and
        'pinv' is used for the pseudo-inverse solution (faster than 'lstsq' for
        large number of cells, but can be unstable for ill-conditioned
        matrices).

    Returns
    -------
    np.ndarray
        Effect sizes from the regression.
    """
    if method == 'onehot':
        if _test_onehot(C):
            beta = (C.T @ X) / np.sum(C, axis=0)[:, None]
        else:
            raise ValueError(
                'Data has to be one-hot encoded to use the onehot method.')
    elif method == 'lstsq':
        beta = scipy.linalg.lstsq(C, X)[0]
    elif method == 'pinv':
        beta = np.linalg.pinv(C) @ X
    else:
        raise ValueError(
            'Method not recognized. Use one of "onehot", "lstsq", or "pinv".')
    return beta


class _LinearRegressionResiduals(BaseEstimator, TransformerMixin):
    """Linear regression without scaling option"""

    def __init__(
            self,
            include_intercept: bool = False,
            method: Literal['onehot', 'lstsq', 'pinv'] = 'lstsq'
    ) -> None:
        self.include_intercept = include_intercept
        self.method = method

    def fit(self, X: ndarray, C: ndarray) -> '_LinearRegressionResiduals':
        if self.include_intercept:
            ones_column = np.ones((C.shape[0], 1))
            C = np.hstack((ones_column, C))
        self.beta = _linear_regression_effectsizes(X, C, self.method)
        return self

    def transform(self, X: ndarray, C: ndarray) -> ndarray:
        if self.include_intercept:
            ones_column = np.ones((X.shape[0], 1))
            C = np.hstack((ones_column, C))
        residual = X - C @ self.beta
        return residual

    def fit_transform(self, X: ndarray, C: ndarray) -> ndarray:
        return self.fit(X, C).transform(X, C)

    def get_effectsizes(self) -> ndarray:
        return self.beta


class LinearRegressionResiduals(Pipeline):
    """
    Get residuals from linear regression.

    Parameters
    ----------
    include_intercept : bool
        Whether to include an intercept in the regression. Not needed if the
        data matrix is already centered.
    standard_scale : bool
        Whether to standard scale the residuals before returning them.
    method : {'onehot', 'lstsq', 'pinv'}
        Method to use for the regression. 'onehot' can be used when the
        covariates are one-hot encoded and is more efficient for large number
        of cells, 'lstsq' is used for the least squares solution (stable), and
        'pinv' is used for the pseudo-inverse solution (faster than 'lstsq' for
        large number of cells, but can be unstable for ill-conditioned
        matrices).

    Attributes
    ----------
    beta : np.ndarray
        Effect sizes from the regression.
    scaler : sklearn.preprocessing.StandardScaler
        Scaler used to standardize the data.
    """

    def __init__(
            self,
            standard_scale: bool = True,
            standard_scale_residuals: bool = True,
            include_intercept: bool = False,
            method: Literal['onehot', 'lstsq', 'pinv'] = 'lstsq'
    ) -> None:
        """
        Initialize the transformer.

        Parameters
        ----------
        standard_scale : bool
            Whether to standard scale the residuals before returning them.
        include_intercept : bool
            Whether to include an intercept in the regression. Not needed if
            the data matrix is already centered.
        standard_scale_residuals : bool
            Whether to standard scale the residuals before returning them.
        method : {'onehot', 'lstsq', 'pinv'}
            Method to use for the regression. 'onehot' can be used when the
            covariates are one-hot encoded and is more efficient for large
            number of cells, 'lstsq' is used for the least squares solution
            (stable), and 'pinv' is used for the pseudo-inverse solution
            (faster than 'lstsq' for large number of cells, but can be unstable
            for ill-conditioned matrices).
        """
        steps = [
            ("linreg", _LinearRegressionResiduals(
                include_intercept=include_intercept, method=method)),
        ]
        if standard_scale:
            steps.insert(0, ["scaler", StandardScaler()])
        if standard_scale_residuals:
            steps.append(["residual scaler", StandardScaler()])
        super().__init__(steps)
        self.standard_scale = standard_scale
        self.standard_scale_residuals = standard_scale_residuals
        self.include_intercept = include_intercept
        self.method = method

    def get_effectsizes(self) -> ndarray:
        """
        Get the effect sizes from the regression.

        Returns
        -------
        np.ndarray
            Effect sizes from the regression.
        """
        return self.named_steps["linreg"].get_effectsizes()