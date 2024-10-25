from typing import Literal, Union

import numpy as np
import scipy.linalg
from numpy import ndarray

from sklearn.preprocessing import StandardScaler
from sklearn.base import TransformerMixin
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
            y: ndarray,
            method: Literal['onehot', 'lstsq', 'pinv'] = 'lstsq'
        ) -> Union[ndarray, ndarray]:
    """
    Perform linear regression and return the effect sizes.

    Parameters
    ----------
    X : np.ndarray
        Data matrix.
    y: np.ndarray
        Outcome matrix
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
        if _test_onehot(X):
            beta = (X.T @ y) / np.sum(X, axis=0)[:, None]
        else:
            raise ValueError(
                'Data has to be one-hot encoded to use the onehot method.')
    elif method == 'lstsq':
        beta = scipy.linalg.lstsq(X, y)[0]
    elif method == 'pinv':
        beta = np.linalg.pinv(X) @ y
    else:
        raise ValueError(
            'Method not recognized. Use one of "onehot", "lstsq", or "pinv".')
    return beta


class _LinearRegressionResiduals(TransformerMixin):
    """Linear regression without scaling option"""

    def __init__(
            self,
            include_intercept: bool = True,
            method: Literal['onehot', 'lstsq', 'pinv'] = 'lstsq'
    ) -> None:
        self.include_intercept = include_intercept
        self.method = method

    def fit(self, X: ndarray, y: ndarray) -> '_LinearRegressionResiduals':
        if self.include_intercept:
            if self.method == 'onehot':
                raise ValueError(
                    'One-hot encoding is not supported with an intercept. '
                    'Set intercept to False and standard scale outcome before.'
                )
            ones_column = np.ones((X.shape[0], 1))
            X = np.hstack((ones_column, X))
        self.beta = _linear_regression_effectsizes(X, y, self.method)
        return self

    def transform(self, X: ndarray, y: ndarray) -> ndarray:
        if self.include_intercept:
            if self.method == 'onehot':
                raise ValueError(
                    'One-hot encoding is not supported with an intercept.'
                    'Set intercept to False and standard_scale to True.'
                )
            ones_column = np.ones((X.shape[0], 1))
            X = np.hstack((ones_column, X))
        residual = y - X @ self.beta
        return residual

    def fit_transform(self, X: ndarray, y: ndarray) -> ndarray:
        return self.fit(X, y).transform(X, y)

    def get_effectsizes(self) -> ndarray:
        return self.beta


class LinearRegressionResiduals(Pipeline):
    """
    Get residuals from linear regression.

    Parameters
    ----------
    standard_scale_covars : bool, default=True
        Whether to standard scale the covariates matrix.
    include_intercept : bool, default=True
        Whether to include an intercept in the regression. Not needed if
        the data matrix is already centered and standard_scale_covars is
        False.
    method : {'onehot', 'lstsq', 'pinv'}
        Method to use for the regression.
            - 'lstsq' is used for the least squares solution. This method is
                stable with respect to ill-conditioned covariate matrices, but
                can be slow for large number of cells.
            - 'onehot' can be used when the covariates are one-hot encoded
                and is more efficient for large number of cells. The data
                matrix should be standardized before using this method, the
                options standard_scale_covars and include_intercept should be
                set to False.
            - 'pinv' is used for the pseudo-inverse solution. This method is
                faster than 'lstsq' for large number of cells, but can be
                unstable for ill-conditioned matrices.
    """

    def __init__(
            self,
            standard_scale_covars: bool = True,
            include_intercept: bool = True,
            method: Literal['onehot', 'lstsq', 'pinv'] = 'lstsq'
    ) -> None:
        """
        Initialize the transformer.

        Parameters
        ----------
        standard_scale_covars : bool, default=True
            Whether to standard scale the covariates matrix.
        include_intercept : bool, default=True
            Whether to include an intercept in the regression. Not needed if
            the data matrix is already centered and standard_scale_covars is
            False.
        method : {'onehot', 'lstsq', 'pinv'}
            Method to use for the regression.
             - 'lstsq' is used for the least squares solution. This method is
                stable with respect to ill-conditioned covariate matrices, but
                can be slow for large number of cells.
             - 'onehot' can be used when the covariates are one-hot encoded
                and is more efficient for large number of cells. The data
                matrix should be standardized before using this method, the
                options standard_scale_covars and include_intercept should be
                set to False.
              - 'pinv' is used for the pseudo-inverse solution. This method is
                faster than 'lstsq' for large number of cells, but can be
                unstable for ill-conditioned matrices.
        """
        steps = [
            ("linreg", _LinearRegressionResiduals(
                include_intercept=include_intercept, method=method)),
        ]
        if standard_scale_covars:
            if method == 'onehot':
                raise ValueError(
                    'One-hot encoding is not supported with standard scaling.'
                    'Set standard_scale_covars to False.'
                )
            steps.insert(0, ("covar scaler", StandardScaler()))
        super().__init__(steps)
        self.standard_scale_covars = standard_scale_covars
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