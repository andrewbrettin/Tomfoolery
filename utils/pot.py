"""
Peak-over-threshold method
==========================

This module implements the peak-over-threshold method for univariate
data.

Available functions
-------------------
Computation:
- select_extremes(data, q=0.95)
- fit_distribution(extremes, fix_loc=True)
- return_level(data, return_period, q=0.95, fitted_params=None)

Plotting:
- plot_threshold(data, q, ax=None)
- plot_extreme_value_distribution(
    extremes, fitted_params=None, bins=50, ax=None
  )
- plot_return_level(return_period, fitted_params=None, ax=None)

"""

__name__ = 'pot'
__author__ = '@andrewbrettin'

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import scipy as sc
from scipy.stats import genpareto

rng = np.random.default_rng()

def select_extremes(data, q=0.95):
    """
    Select an extreme subset of the given data.

    Parameters:
        data : array_like
            The data from which a subset should be selected.
        q : float
            Quantile to set as a threshold. Must be between 0 and 1.

    Returns:
        extremes : array-like
            Subset of extreme values.

    """

    threshold = np.quantile(data, q)
    extremes = np.array(data[data >= threshold])
    return extremes


def fit_distribution(extremes, fix_loc=True):
    """
    Fits distribution parameters to extreme values using the
    generalized pareto distribution.

    Parameters:
        extremes : array_like
            Extreme values.
        fix_loc : bool
            If True, set the location parameter to the threshold.
            This is done to avoid a poorly-fitted extreme value
            distribution to the data.
    Returns:
        fitted_params : tuple of floats
            Parameters corresponding to the fitted distribution.
            For generalized pareto distribution, the returned
            parameters are (shape, loc, scale).
    """

    # We may eventually include other extreme value distributions
    dist = getattr(sc.stats, 'genpareto')     
    fitted_params = dist.fit(extremes, floc=np.min(extremes))
    return fitted_params


def return_level(data, T, q=0.95, fitted_params=None):
    """
    Computes the return level for a given dataset.

    Parameters:
        data : array_like
            Data to compute return levels for.
        T : int
            Integer indicating the return period of the extreme event
            to compute. For instance, for daily data, to compute a
            one-year return level, set `return_period=365`.
        q : float
            Quantile threshold for extreme values.
        fitted_params : tuple of floats
    
    Returns : 
        z_T : float
            Return level.
    """
    
    extremes = select_extremes(data, q)
    
    if fitted_params is None:
        fitted_params = fit_distribution(extremes, fix_loc=True)
    shape, loc, scale = fitted_params
    
    p = 1/T
    z_T = loc + scale/shape * (((1-q)/p)**shape - 1)
    
    return z_T

def plot_threshold(data, q, ax=None, show_legend=True, **kwargs):
    """
    Plots a timeseries and quantile threshold.

    Parameters:
        data : array_like
            Timeseries to plot.
        q : float
            Quantile to set as a threshold. Must be between 0 and 1.
        ax : matplotlib.axes.Axes
            Plotting axis.
        kwargs : (optional)
            Additional keyword arguments to pass to matplotlib.
    Returns: 
        ax : matplotlib.axes.Axes
            Plotting axis.
    """

    threshold = np.quantile(data, q)

    if ax is None:
        fig, ax = plt.subplots()
    
    if isinstance(data, xr.DataArray):
        # Time is assumed to be the only dimension
        data.plot(ax=ax, label=data.name, **kwargs)
        if len(data.dims) != 1:
            raise ValueError('Input data has more than one dimension')
        
        ax.hlines(
            threshold, 0, len(data), ls='--',
            color='k', label='{}th percentile'.format(int(q*100))
        )
    else:
        ax.plot(data)
        N = len(data)
        ax.hlines(
            threshold, 0, N, ls='--', color='k',
            label='{}th percentile'.format(int(q*100))
        )
    if show_legend:
        ax.legend()
    
    return ax


def plot_extreme_value_distribution(extremes, fitted_params=None, bins=50,
        ax=None):
    """
    Plots the fitted extreme value distribution probability density
    function overlayed over a histogram of the empircal extreme
    values.
    
    Parameters:
        extremes : array_like
            Extreme values.
        fitted_params : tuple of floats
            parameters corresponding to the fitted distribution.
        bins : int
            Number of bins to use in histogram
        ax : matplotlib.axes.Axes
            Plotting axis.
        kwargs : (optional)
            Additional keyword arguments to pass to matplotlib.
    Returns: 
        ax : matplotlib.axes.Axes
            Plotting axis.
    """

    if ax is None:
        fig, ax = plt.subplots()

    if fitted_params is None:
        fitted_params = fit_distribution(extremes, fix_loc=True)
    
    # Plot histogram of extreme values
    ax.hist(extremes, bins=bins, density=True, label='Extreme values')

    # Plot MLE pdf
    quantiles = np.linspace(
        genpareto.ppf(0.01, *fitted_params),
        genpareto.ppf(0.99, *fitted_params),
        99
    )
    fitted_pdf = genpareto.pdf(quantiles, *fitted_params)
    ax.plot(quantiles, fitted_pdf, label='Fitted GPD pdf')

    ax.set(
        title='Best fit generalized Pareto distribution to extremes',
        xlabel='Extreme values [m]',
        ylabel='Density'
    )
    ax.legend(loc='upper right')

    return ax


def plot_return_level(data, T, fitted_params=None, q=0.95,
        ax=None, show_legend=True, **kwargs):
    """
    For a given timeseries dataset, computes a return level 
    for a given return period.
    
    Parameters:
        data : array_like
            The data from which a subset should be selected.
        T : int or array of ints
            Integers indicating the return period of the extreme event
            to plot. For instance, for daily data, to compute a
            one-year return level, set `return_period=365`.
        fitted_params : tuple of floats
            Parameters corresponding to the fitted distribution.
            For generalized pareto distribution, the returned
            parameters are (shape, loc, scale).
        q : float
            Quantile to set as a threshold. Must be between 0 and 1.
        ax : matplotlib.axes.Axes
            Plotting axis.
        kwargs : (optional)
            Additional keyword arguments to pass to matplotlib.
    Returns: 
        ax : matplotlib.axes.Axes
            Plotting axis.
    """
    
    if fitted_params is None:
        extremes = select_extremes(data, q)
        fitted_params = fit_distribution(extremes, fix_loc=True)
    
    if isinstance(T, int):
        return_levels = [
            return_level(data, T, q=0.95, fitted_params=fitted_params)
        ]
    else:
        return_levels = [
            return_level(data, level, q=0.95, fitted_params=fitted_params)
            for level in T
        ]
    
    if ax is None:
        fig, ax = plt.subplots()
    
    if isinstance(data, xr.DataArray):
        # Time is assumed to be the only dimension
        data.plot(ax=ax, label=data.name, **kwargs)
        if len(data.dims) != 1:
            raise ValueError('Input data has more than one dimension')
        
        ax.hlines(
            return_levels, 0, len(data), ls='--',
            color='k', label=''
        )
    else:
        ax.plot(data)
        N = len(data)
        ax.hlines(
            return_levels, 0, N, ls='--', color='k',
            label=''
        )
    if show_legend:
        ax.legend()
    
    return ax
    
    