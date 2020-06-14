# coding=utf-8
"""
Functions for plotting
"""
import warnings

import cmocean
import cv2
import matplotlib.colors as colors
import matplotlib.pylab as plt
import numpy as np
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import os

sys.path.append(os.path.dirname(__file__))
try:
    from model.load_model_output import *
except ModuleNotFoundError:
    print("No matched hydrographs here")

warnings.filterwarnings('ignore')

plt.style.use('ggplot')
import seaborn as sns
import pandas as pd

sns.set(font_scale=1.3)

sns.set_style("whitegrid", {'axes.linewidth': 1.0})
sns.set_context("notebook", font_scale=1.3,
                rc={"lines.linewidth": 1})
col_list = ["cool blue", "light grey",
            "viridian", "twilight blue",
            "dusty purple", "amber",
            "greyish", "faded green"]

col_list_palette = sns.xkcd_palette(col_list)

plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['grid.color'] = 'grey'
plt.rcParams['grid.alpha'] = 0.0
plt.rcParams['axes.linewidth'] = 0.5
plt.rc('axes', edgecolor='grey')

plt.rcParams['axes.spines.top'] = 0
plt.rcParams['axes.spines.right'] = 0
plt.rcParams['axes.spines.left'] = 1
plt.rcParams['axes.spines.bottom'] = 1
plt.rc('axes', edgecolor='grey')
plt.rcParams['image.cmap'] = 'Blues'

from itertools import cycle

lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)


def veg_points(veg, dx=1.0, veg_size=10, ax=None, alpha=0.75):
    """
    Parameters
    ----------
    veg
    dx
    veg_size
    ax
    alpha
    """
    xc, yc = create_grid(veg, dx)

    veg = (veg * veg_size).astype(float)

    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 6))

    ax.scatter(xc + dx / 2, yc + dx / 2.,
               s=veg.T,
               c="g", marker='o', alpha=alpha)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)


def veg_pcolor(veg, dx=1.0, ax=None):
    """

    Parameters
    ----------
    veg
    dx
    ax

    Returns
    -------

    """
    xc, yc = create_grid(veg, dx)

    veg = veg.astype(float)

    veg[veg == 0] = np.nan

    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 6))
    else:
        pass

    bounds = np.linspace(0, 2, 10)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

    ax.pcolormesh(xc + dx / 2, yc + dx / 2., veg,
                  norm=norm,
                  cmap="Greens")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)

    return ax


def plot_c(r, dx=1, ax='', linewidth=0.5):
    """

    Parameters
    ----------
    r
    dx
    ax
    linewidth

    Returns
    -------

    """

    thresh = r.astype(np.uint8)
    contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    if ax == '':
        fig, ax = plt.subplots()

    for n, contour in enumerate(contours):
        contour = np.squeeze(contour)

        contours[n] = contour
        if len(contour) > 2:
            ax.plot(contour[:, 1] * dx + dx, contour[:, 0] * dx + dx,
                    linewidth=linewidth, c='k')
    return ax


def plot_c_inv(r, dx=1, ax=''):
    """

    Parameters
    ----------
    r
    dx
    ax

    Returns
    -------
    """

    thresh = 1 - r.astype(np.uint8)
    contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    if ax == '':
        fig, ax = plt.subplots()

    for n, contour in enumerate(contours):
        contour = np.squeeze(contour)

        contours[n] = contour
        if len(contour) > 2:
            ax.plot(contour[:, 1] * dx + dx, contour[:, 0] * dx + dx,
                    linewidth=0.7, c='k')
    return ax


def arraycolor(array, ax=None,
             colorbar=True,
             clabel='',
             trim=0,
             cmin=None, cmax=None,
             cfontsize=16,
             cmap="Blues"):


    array = array[:, trim:]

    xc, yc = create_grid(array, dx= 1)


    if not ax:
        fig = plt.figure(figsize=(4, 6))
        ax = fig.add_subplot(111)
    else:
        fig = plt.gcf()

    if not cmin:
        cmin = np.nanmin(array.ravel())
    if not cmax:
        cmax = np.nanmax(array.ravel())

    bounds = np.linspace(cmin, cmax, 100)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

    if cmax == cmin:
        colorbar = False
    zinflplot = ax.pcolormesh(xc.T, yc.T, array.T,
                              norm=norm,
                              edgecolors="face",
                              cmap=cmap, alpha=1)

    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        cb = fig.colorbar(zinflplot, cax=cax, shrink=1)
        cb.set_label(clabel, fontsize=cfontsize)
        cb.ax.tick_params(labelsize=cfontsize)

        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()

    ax.set_ylim(yc.min(), yc.max())
    ax.set_xlim(xc.min(), xc.max())
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)

    return ax, zinflplot


def colormap(sim, array, ax=None,
             colorbar=True, veg_scale=False,
             clabel='',
             trim=0,
             cmin=None, cmax=None,
             veg_size=2, plot_veg=False,
             cfontsize=16,
             contour_type=None,
             cmap="Blues"):
    """

    Parameters
    ----------
    sim
    array
    ax
    colorbar
    veg_scale
    clabel
    trim
    cmin
    cmax
    veg_size
    plot_veg
    cfontsize
    contour_type
    cmap

    Returns
    -------

    """
    veg = sim['veg'][:, trim:].astype(float)
    array = array[:, trim:]
    dx = sim['dx']

    xc, yc = create_grid(array, dx)

    if veg.sum() == 0:
        veg_scale = False

    if not ax:
        fig = plt.figure(figsize=(4, 6))
        ax = fig.add_subplot(111)
    else:
        fig = plt.gcf()

    if veg_scale:
        scale_vals = array[veg > 0].ravel()
    else:
        scale_vals = array.ravel()

    if not cmin:
        cmin = np.nanmin(scale_vals)
    if not cmax:
        cmax = np.nanmax(scale_vals)

    bounds = np.linspace(cmin, cmax, 100)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

    if cmax == cmin:
        colorbar = False
    zinflplot = ax.pcolormesh(xc.T, yc.T, array.T,
                              norm=norm,
                              edgecolors="face",
                              cmap=cmap, alpha=1)

    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        cb = fig.colorbar(zinflplot, cax=cax, shrink=1)
        cb.set_label(clabel, fontsize=cfontsize)
        cb.ax.tick_params(labelsize=cfontsize)

        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()

    if plot_veg:
        ax.scatter(xc + sim.dx / 2., yc + sim.dx / 2.,
                   s=veg * veg_size, c='g',
                   marker='o', alpha=0.75)

    ax.set_ylim(yc.min(), yc.max())
    ax.set_xlim(xc.min(), xc.max())
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)

    if contour_type == 'veg':
        ax = plot_c(veg, dx=sim.dx, ax=ax)
    elif contour_type == 'inv':
        ax = plot_c_inv(veg, dx=sim.dx, ax=ax)

    return ax, zinflplot


"""
"""

def triptych(sim):
    """

    Parameters
    ----------
    sim

    Returns
    -------

    """

    
    if (sim.fV == 1) or (sim.fV==0):
        fig, axes = diptych(sim)
        return fig, axes

    fig, axes = plt.subplots(1, 3, figsize=(15, 3))
    plt.subplots_adjust(wspace=0.3)
    for i, label in enumerate(('Vegetation',
                               'Infiltration depth', 'Max velocity')):
        ax = axes[i]
        ax.set_title(label, fontsize=16)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)

    ax = axes[0]
    veg_pcolor(sim.veg, dx=sim.dx, ax=ax)

    ax = axes[1]
    infl_2d_cm = sim['infl_2d'].copy()*100
    colormap(sim, infl_2d_cm, ax=ax,
             clabel='$I$ (cm)',
             cmap=cmocean.cm.deep,
             colorbar=True, cmin=0)

    ax = axes[2]
    U_max = np.max(np.sqrt(sim.uc ** 2 + sim.vc ** 2), 0) * 100
    U_max[-1,-1] = 0    
    colormap(sim,U_max, ax=ax,
             clabel='velocity (cm/s)', cmin = 0,
             colorbar=True, cmap="Blues")
    return fig, axes

def diptych(sim):
    """

    Parameters
    ----------
    sim

    Returns
    -------

    """
    if sim.Ks < 0.001:
        fig, axes = plotUmax(sim)
        return fig, axes
        
    fig, axes = plt.subplots(1, 2, figsize=(11, 3))
    plt.subplots_adjust(wspace=0.4)
    
    for i, label in enumerate(('Infiltration depth', 'Max velocity')):
        ax = axes[i]
        ax.set_title(label, fontsize=16)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)


    ax = axes[0]
    infl_2d_cm = sim['infl_2d'].copy()*100
    colormap(sim, infl_2d_cm, ax=ax,
             clabel='$I$ (cm)',
             cmap=cmocean.cm.deep,
             colorbar=True, cmin=0)

    ax = axes[1]
    U_max = np.max(np.sqrt(sim.uc ** 2 + sim.vc ** 2), 0) * 100
    U_max[-1,-1] = 0    
    colormap(sim,U_max, ax=ax,
             clabel='velocity (cm/s)', cmin = 0,
             colorbar=True, cmap="Blues")
    return fig, axes

def plotUmax(sim):
    """

    Parameters
    ----------
    sim

    Returns
    -------
    """
    fig, axes = plt.subplots(1, figsize=(5, 3))
    plt.subplots_adjust(wspace=0.4)
    ax = axes    
    label = 'Max velocity'

    ax.set_title(label, fontsize=16)
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)

    U_max = np.max(np.sqrt(sim.uc ** 2 + sim.vc ** 2), 0) * 100
    U_max[-1,-1] = 0
    colormap(sim,U_max, ax=ax,
             clabel='velocity (cm/s)', cmin = 0,
             colorbar=True, cmap="Blues")
    return fig, axes


def triptych_micro(sim):
    """
    Plots microtopo in panel A, instead of vegetation
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 3))
    plt.subplots_adjust(wspace=0.3)

    for i, label in enumerate(('Microtopo',
                               'Infiltration depth', 'Max velocity')):
        axes[i].set_title(label, fontsize=16)

    ax = axes[0]
    #micro = sim['zc'] + sim['xc'].mean(0)*sim['So']
    micro = sim['zc'] - sim['zc'].mean(0)
    colormap(sim, micro, ax=ax,
             colorbar=True,
             cmap="gist_gray")

    ax = axes[1]
    colormap(sim, sim['infl_2d']*100, ax=ax,
             clabel='$I$ (cm)',
             colorbar=True)

    ax = axes[2]
    U_max = np.max(np.sqrt(sim.uc ** 2 + sim.vc ** 2), 0) * 100
    U_max[-1, -1] = 0

    colormap(sim, U_max, ax=ax,
             clabel='velocity (cm/s)',
             colorbar=True, cmap="Blues",
             )
    return fig, axes

    
def place_AB(axes, xloc=-0.1, yloc=1.15):
    """
    snippet to place As and Bs in figures 
    """
    axes = axes.ravel()
    labels = ['A)', 'B)', 'C)', 'D)', 'E)']
    for i, ax in enumerate(axes):
        ax.text(xloc, yloc, labels[i], transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top')
    return axes


def plot_hydrographs(core, nonzero=True, ax=None):
    """
    hydrographs (flux3) with negative values removed

    Parameters
    ----------
    ax
    nonzero
    core : pd.DataFrame
        SVE simulations

    """

    if not ax:
        fig, ax = plt.subplots(1)
    if type(core) == dict:
        core = pd.DataFrame(core).T

    for key in core.index:
        sim = core.loc[key]
        t_h = np.arange(len(sim['flux3']))
        flux3 = sim['flux3']
        if nonzero:
            flux3[flux3 < 0] = 0
        flux3_cm_hr = flux3 * 3.6e5 / sim["Ly"] / sim["Lx"]
        ax.plot(t_h / 60., flux3_cm_hr)

    ax.set_xlabel('minutes')
    ax.set_ylabel('cm/hr')
    return fig, ax


def plot_inflgraphs(core, nonzero=False, ax=None):
    """
    inflgraphs with negative values removed

    Parameters
    ----------
    ax
    nonzero
    core : pd.DataFrame
        SVE simulations

    """
    if not ax:
        ax = plt.gca()
    if type(core) == dict:
        core = pd.DataFrame(core).T

    for key in core.index:
        sim = core.loc[key]
        t_h = np.arange(len(sim['infl_1d']))
        infl_1d = sim['infl_1d']
        if nonzero:
            infl_1d[infl_1d < 0] = 0
        infl_1d_cm_hr = infl_1d * 3.6e5 / sim["Ly"] / sim["Lx"]
        ax.plot(t_h / 60., infl_1d_cm_hr, label = sim.Ks)
    
    ax.legend()
    ax.set_xlabel('minutes')
    ax.set_ylabel('cm/hr')
    return ax

    
def create_grid(array, dx):
    """
    Creates an x,y meshgrid, with size array and discretization dx

    Parameters
    ----------
    array
    dx

    Returns
    -------
    xc, yc
    """

    nrow = array.shape[0]
    ncol = array.shape[1]

    xc = np.arange(0, ncol * dx, dx) + dx / 2
    yc = np.arange(0, nrow * dx, dx) + dx / 2
    xc, yc = np.meshgrid(xc, yc)

    return xc, yc

"""
Gridded plots
"""
def plot_veg_grid(subset):
    """
    TODO: ignore cases where fV=1

    Parameters
    ----------
    subset

    Returns
    -------

    """
    
    num_sim = len(subset)

    smplkeys = subset.index
    num_col = min(4, len(subset))
    num_row = int(np.ceil(num_sim / 1. / num_col))

    figsize = (num_col * 4, num_row * 2.2)
    fig, axes = plt.subplots(num_row, num_col, figsize=figsize)
    plt.subplots_adjust(hspace=0.1)

    axes = axes.ravel()

    for i, ax in enumerate(axes[:num_sim]):
        key = smplkeys[i]
        sim = subset.loc[key]
        if (sim.fV == 1) or (sim.fV == 0):
            continue
        veg_pcolor(sim.veg, ax=ax)
    for i, ax in enumerate(axes[num_sim:]):
        ax.set_visible(False)
    return fig


def plot_infl_grid(subset):
    """
    Parameters
    ----------
    subset

    Returns
    -------
    """
    num_sim = len(subset)

    smplkeys = subset.index
    num_col = min(4, len(subset))
    num_row = int(np.ceil(num_sim / 1. / num_col))

    figsize = (num_col * 4, num_row * 2.2)
    fig, axes = plt.subplots(num_row, num_col, figsize=figsize)
    plt.subplots_adjust(hspace=0.1)

    axes = axes.ravel()

    cmax = subset.loc[smplkeys].infl_99th.max()
    zinflplot = []
    for i, ax in enumerate(axes[:num_sim]):
        key = smplkeys[i]
        sim = subset.loc[key]
        ax, zinflplot = colormap(sim, sim.infl_2d, ax=ax, cmax=cmax,
                                 colorbar=False)
    for i, ax in enumerate(axes[num_sim:]):
        ax.set_visible(False)

    cbaxes = fig.add_axes([.91, 0.125, 0.013, 0.755])
    cb = fig.colorbar(zinflplot, cax=cbaxes, shrink=1)
    cb.set_label('$I$ (m)', fontsize=18)
    cb.ax.tick_params(labelsize=18)

    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()

    return fig

"""
SVE mass balance checks
"""

def summarize_mass_balance(core, legend=False, ax = None):
    """

    Parameters
    ----------
    core
    legend

    Returns
    -------

    """
    if not ax:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = plt.gcf()

    for key in core.index:
        sim = core.loc[key]
        rain_vol = sim.rain_series*sim.Ly*sim.Lx
        name = sim.name[:-1].replace("/", ", ").replace(",", ", ")
        ax.plot(sim.t_h / 60., np.cumsum(rain_vol - sim.infl_1d +
                                              - sim.boundary_flux_1d - sim.dvol_1d), label=name)

    if legend:
        ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')

    ax.set_title("Cumulative mass balance error")
    ax.ticklabel_format(axis='xxx', style='sci', scilimits=(0, 0))
    ax.set_ylabel("m$^3$/s")
    ax.set_xlabel("Time (min)")

    return fig


def summarize_fluxes(sim, ):
    """
    Plot global mass balance components to check
    that they add to zero.

    Fluxes at each timestep:
        rain, infiltration

    scale: converts from /timestep to /hr

    """
    fig, ax = plt.subplots(figsize=(6, 4))

    # m3/s*3.6e5 (s cm)/(hr m) /m2
    scale = 3.6e5/sim.Lx/sim.Ly
    ax.plot(sim.t_h / 60.,  sim.infl_1d * scale, next(linecycler), 
        label="$I$: infiltration")
    ax.plot(sim.t_h / 60., sim.boundary_flux_1d * scale, next(linecycler), label="$q$: flux out")
    ax.plot(sim.t_h / 60., sim.dvol_1d * scale, next(linecycler), 
        label="$dV$: volume change ")
    ax.set_xlabel("Time (min)")

    ax.plot(sim.t_h / 60., sim.rain_series * 3.6e5, label="$p$: rain")
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')

    ax.set_ylabel("cm/hr")
    ax.set_title("Boundary fluxes")

    return fig


def plot_fluxes(sim):
    """
    Plot to check boundary fluxes
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    ax = axes[0]

    scale = 3.6e5/sim.Lx/sim.Ly
    ax.plot(sim.t_h / 60., sim.fluxin * scale,label="$1$")
    ax.plot(sim.t_h / 60., sim.flux2 * scale, label="$2$")
    ax.plot(sim.t_h / 60., sim.flux3 * scale, label="$3$")
    ax.plot(sim.t_h / 60., sim.flux4 * scale, label="$4$")
    ax.set_ylabel("cm/hr")
    ax.set_xlabel("Time (min)")
    ax.set_title("Fluxes 1-4")
    ax.legend()

    ax = axes[1]
    test_flux = - sim.flux1 + sim.flux2 + sim.flux3 - sim.flux4
    ax.plot(sim.t_h / 60., test_flux * scale, label="$q$(1-4)")
    ax.plot(sim.t_h / 60., sim.boundary_flux_1d  * scale, '--',
            label="$q$(total)")
    ax.plot(sim.t_h / 60., sim.flux3 * scale, '--',
            label="hydrograph")

    ax.legend()
    ax.set_xlabel("Time (min)")
    ax.set_title("Total flux out")


def check_cell_fluxes(sim, xi=-1, yi=1):
    """
    Plot all the fluxes in an out of cell (yi, xi)

    Notes:
    ------
    Sign convention is opposite to the direction of water flow:
        xflux0 = downslope flux into the cell
        xflux1 = downslope flux out of the cell

    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex='all', sharey='all')
    plt.subplots_adjust(wspace=0.4)

    flux_out = (sim.xflux1[:, yi, xi] - sim.xflux0[:, yi, xi]
                     + sim.yflux1[:, yi, xi] -  sim.yflux0[:, yi, xi]) # m3/s
    vol_this_cell = sim.hc[:, yi, xi] * sim.dx ** 2

    d_vol_this_cell = np.diff(vol_this_cell)
    d_vol_this_cell = np.insert(d_vol_this_cell, 0, 0)/sim.dt_print # m3/s

    rain_inputs = np.zeros(len(sim.t_print))
    rain_inputs[sim.t_print <= sim.t_rain] = sim.rain*sim.dx**2
    rain_inputs[sim.t_print == 0] = 0

    ax = axes[0]
    ax.plot(sim.t_print / 60., (rain_inputs - flux_out -
                            sim.infl_3d[:, yi, xi] - d_vol_this_cell) , '--')
    ax.set_title("Mass balance")
    ax.set_ylabel("m3/s")
    ax.set_xlabel("min")

    ax = axes[1]
    ax.plot(sim.t_print / 60.,  sim.infl_3d[:, yi, xi] , next(linecycler), label="$F$")
    ax.plot(sim.t_print / 60., flux_out , next(linecycler),
            label="$q$")
    ax.plot(sim.t_print / 60., rain_inputs, next(linecycler), label="$p$")
    ax.plot(sim.t_print / 60., d_vol_this_cell, next(linecycler),
            label=r"$\Delta V$")
    ax.set_title("Components ")
    ax.legend(loc="lower right")
    ax.set_xlabel("min")

    ax = axes[2]
    ax.plot(sim.t_print / 60.,  sim.xflux0[:, yi, xi] , next(linecycler),
            label=r"$q_{x}(in)$")
    ax.plot(sim.t_print / 60.,  sim.xflux1[:, yi, xi] , next(linecycler),
            label=r"$q_{x}(out)$")

    ax.set_xlabel("min")
    ax.legend(loc="best")
    ax.set_title("Downslope fluxes")

    return fig


def check_band_fluxes(sim, xi=0):
    """
    Plot the fluxes into and out of a given cross-slope band (`xi`).

    Left panel: downslope lateral flux (normalized to hillslope area)
    Right panel: mass balance (normalized to grid cell area)

    Notes:
    ------
    choose xi = -1 for hillslope outlet
    """

    flux_out = sim.xflux1[:, :, xi].mean(1) \
                     - sim.xflux0[:, :, xi].mean(1)

    vol_this_band =  sim.hc[:, :, xi].mean(1) * sim.dx ** 2
    d_vol_this_band = np.diff(vol_this_band)
    d_vol_this_band = np.insert(d_vol_this_band, 0, 0)/sim.dt_print

    rain_inputs = np.zeros(len(sim.t_print))
    rain_inputs[sim.t_print <= sim.t_rain] = sim.rain*sim.dx**2
    rain_inputs[sim.t_print == 0] = 0

    band_infl = sim.infl_3d[:, :, xi].mean(1)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharex='all', sharey='all')

    ax = axes[0]
    ax.plot(sim.t_print / 60., (rain_inputs - flux_out -
                             band_infl - d_vol_this_band) , '--')
    ax.set_title("Mass balance")
    ax.set_ylabel("m$^3$/s")
    ax.set_xlabel("min")


    ax = axes[1]
    ax.plot(sim.t_print / 60., band_infl, next(linecycler), label="$f$")
    ax.plot(sim.t_print / 60., flux_out , next(linecycler), label="$q$")
    ax.plot(sim.t_print / 60., rain_inputs, next(linecycler), label="$p$")
    ax.plot(sim.t_print / 60., d_vol_this_band , next(linecycler), label=r"$\Delta V$")

    ax.set_title("Components")
    ax.legend(loc='upper right')


    ax = axes[2]
    xflux1_band =  sim.xflux1[:, :, xi].mean(1)
    ax.plot(sim.t_print / 60., xflux1_band, next(linecycler),
            label=r"$q_{x}$(out)")

    xflux0_band = sim.xflux0[:, :, xi].mean(1)
    ax.plot(sim.t_print / 60., xflux0_band , next(linecycler),
            label=r"$q_{x}$(in)")
    ax.legend(loc = "upper right")

    ax.set_title("Downslope fluxes")
    ax.set_xlim(0, )
    ax.set_xlabel("min")

    return fig

def flux_box(sim):
    """
    Useful?
    """
    fig, axes = plt.subplots(2,2, figsize = (10,6), sharex= True)
    axes = axes.ravel()
    scale = 3.6e5/sim.Lx/sim.Ly
    axes[0].plot(sim.flux1*scale , label = "1")
    axes[1].plot(sim.flux2*scale,  label = "2")
    axes[2].plot(sim.flux3*scale,  label = "3")
    axes[3].plot(sim.flux4*scale,  label = "4")

    for i,ax in enumerate(axes):
        ax.set_title(i+1)
    axes[0].set_ylabel("cm/hr")
    axes[2].set_ylabel("cm/hr")

    return fig


"""
Inflow plots 
"""

def plot_inflow(sim,  t_f = None, N_profile=5):
    """
    Plot inflow profile at select times


    """
    fig, ax = plt.subplots(1, figsize=(7, 4.5))

    if t_f:
        inds = np.where(sim.t_print < t_f)[0]
    else:
        inds = np.where(np.diff(sim.hc.mean(1).mean(1)) > 1e-6)[0]

    freq = int(len(inds) / N_profile)

    for ind in inds[freq-1::freq]:
        ax.plot(sim.xc.mean(0), sim['hc'].mean(1)[ind], next(linecycler),
                label=np.round(sim.t_print[ind] / 60., 1))

    legend = ax.legend(title="time (min)")

    plt.setp(legend.get_title(),fontsize='medium')
    plt.xlabel('x')
    plt.ylabel("h (cm)")



def plot_all_inflow(sim, t_f = None, freq=3, label_axes=True, ax = None):
    """
    Plot inflow profile at select times

    Parameters:
    -----------
    sim : dict
        SVE simulation
    freq : int
        plot frequency
    """
    if not ax:
        fig, ax = plt.subplots(1, figsize=(7, 4.5))
    else:
        fig = plt.gcf()

    scale = 3.6e5/sim.Lx
    sim['bad'] = sim['hc'] * sim['uc']*scale
    sim['qc'] = sim['xflux1']/sim.dx*scale

    if t_f:
        inds = np.where(sim.t_print < t_f)[0]
    else:
        inds = np.where(np.diff(sim.hc.mean(1).mean(1)) > 1e-6)[0]

    for i in inds[::freq]:
        ax.plot(sim.xc.mean(0),  sim['qc'].mean(1)[i], 'b--')

    if label_axes:
        ax.set_xlabel('x')
        ax.set_ylabel("q (cm/hr)")
    else:
        ax.set_xticklabels("")
        ax.set_yticklabels("")
    ax.set_xlim(0, )
    return fig, ax


def plot_inflowgraphs(core, trim = 0, nonzero=False, ax=None):
    """
    flux1graphs with negative values removed

    Parameters
    ----------
    ax
    nonzero
    core : pd.DataFrame
        SVE simulations

    """
    if not ax:
        ax = plt.gca()
    if type(core) == dict:
        core = pd.DataFrame(core).T

    for key in core.index:
        sim = core.loc[key]
        t_h = sim.t_h
        flux1 = sim['flux1']
        if nonzero:
            flux1[flux1 < 0] = 0
        flux1_cm_hr = flux1 * 3.6e5 / sim["Ly"] / sim["Lx"]
        ax.plot(t_h[trim:] / 60., flux1_cm_hr[trim:])

    ax.set_xlabel('minutes')
    ax.set_ylabel('cm/hr')
    return ax


"""
Parameter sensitivy plots
"""
def plot_matched_hydrographs(core, match_var = "dt_sw", ):
    """
    Usage:
        fig1, fig2 = plot_matched_hydrographs(core)
    """
    match_vals = np.unique(core[match_var])
    match_to = match_vals[0]

    subset = core[core[match_var] == match_to]

    num_col = min(4, len(subset))
    num_row = int(np.ceil(len(subset) / 1. / num_col))

    figsize = (num_col * 4, num_row * 2.5)
    fig, axes = plt.subplots(num_row, num_col,
            figsize=figsize, sharey = True, sharex = True)
    plt.subplots_adjust(hspace=0.3)
    axes = axes.ravel()

    fig2, axes2 = plt.subplots(num_row, num_col,
        figsize=figsize, sharey = True, sharex = True)
    plt.subplots_adjust(hspace=0.2)
    axes2 = axes2.ravel()

    for i, ax in enumerate(axes[:len(subset)]):
        ax2 = axes2[i]
        key = subset.index[i]
        sim = subset.loc[key]
        ax.plot(sim.t_h/60., sim.flux3/sim.dx**0,
                label = "{0} = {1}".format(match_var, sim[match_var]))
        ax.set_title(sim.pretty)

        for other_val in match_vals[1:]:

            params = extract_match_params(sim)
            params[match_var] = other_val
            matched_sim = filter_core(core, params).iloc[0]
            ax.plot(matched_sim.t_h/60, matched_sim.flux3/matched_sim.dx**0,
                label = "{0} = {1}".format(match_var, matched_sim[match_var]))

            tf = min(len(sim.flux3), len(matched_sim.flux3))
            numer =  (sim.flux3[:tf]/sim.dx**0 - matched_sim.flux3[:tf]/matched_sim.dx**0)
            denom = np.max(sim.flux3[:tf]/sim.dx**0)
            ax2.plot(sim.t_h[:tf]/60,
                    numer/denom*100,
                     label = "{0} = {1}".format(match_var, matched_sim[match_var]))

    for i, ax in enumerate(axes[len(subset):]):
        ax.set_visible(False)

    for i, ax2 in enumerate(axes2[len(subset):]):
        ax2.set_visible(False)

    axes[0].set_ylabel("m3/s")
    axes[4].set_ylabel("m3/s") if len(subset) > 4 else 0
    axes2[0].set_ylabel("% difference")
    axes2[4].set_ylabel("% difference") if len(subset) > 4 else 0

    axes[0].legend()
    axes2[0].legend()
    return fig, fig2


def extract_match_params(sim):
    """
    Return sim parameters in a dictionary

    """
    keys = ['q1_m2hr', 'Ks', 'dt_sw', 'tr', 'p', 'tmax_scale', 'dt_print',
            'save_fluxes', 'save_sve', 'dx',
            'veg_type', 'fV', 'grad_fV', 'seed', 'sigma_scale', 'sigma',
            'stripe_count', 'downslope', 'spots',
            'topo', 'So', 'imodel', 's_scale', 'theta_r', 'theta_s',
            'theta_i', 'H_i', 'Ao', 'scheme', 'alpha', 'alphaB',
            'itype1', 'itype3', 'itype2', 'itype4',  'epsh']

    params = {}
    for key in keys:
        if key in sim.keys():
            params[key] = sim[key]
    return params
