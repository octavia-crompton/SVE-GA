# coding=utf-8
"""
TODO: add 3d plots with subplot index; default None
"""
from plot_functions import *
from matplotlib import cm

from mpl_toolkits.mplot3d import Axes3D


def fix_3D_axes(ax):
    """
    Set up 3D axes (and get rid of colored axis planes)
    """
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    ax.set_xticks([], [])
    ax.set_zticks([], [])
    ax.set_yticks([], [])
    ax.grid(False)

    return ax


"""
Plot veg patterns
"""

def plot_3D_veg_dots(sim):
    """
    Plot veg pattern as dots
    """
    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_subplot(111, projection='3d')

    ax = fix_3D_axes(ax)
    plt.Normalize()
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    ax.scatter(sim.xc[sim.veg == 1],
               sim.yc[sim.veg == 1],
               sim.zc[sim.veg == 1],
               c='g', marker='o', s=20, alpha=1)
    ax.set_ylim(0, sim.yc.max())
    ax.set_zlim(0, sim.zc.max())

    ax.view_init(25, 295)
    return fig, ax


def plot_3D_veg(sim, plot_infl=False):
    """
    Plot veg pattern as surface
    """
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(111, projection='3d')

    ax = fix_3D_axes(ax)

    norm = plt.Normalize()
    isveg = sim.veg.copy()
    isveg = isveg.astype(float)

    isveg[isveg == 0] = 0.1
    isveg[isveg == 1] = 0.9
    isveg[0, -1] = 0
    isveg[1, -1] = 1
    veg_colors = cm.Greens(norm(isveg))

    xc = sim.xc
    yc = sim.yc
    topo = sim.zc

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.plot_surface(xc, yc, topo,
                    facecolors=veg_colors,
                    rstride=1, cstride=1,
                    linewidth=0,
                    antialiased=True,
                    shade=False,
                    alpha=0.8)

    ax.view_init(25, 295)

    if plot_infl:
        norm = plt.Normalize(vmin=0)
        colors = cm.Blues(norm(sim.infl_2d))
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

        ax.plot_surface(xc, yc, topo, facecolors=colors,
                        rstride=1, cstride=1, linewidth=0,
                        antialiased=True, shade=False)
    return fig, ax


def plot_3D_infl(sim, vmax=None):
    """
    Plots infiltration `infl_2d` on incline

    Parameters
    ----------
    sim
    vmax
    """
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection='3d')
    ax = fix_3D_axes(ax)

    if not vmax:
        vmax = np.percentile(sim.infl_2d, 99)

    norm = plt.Normalize(vmin=0, vmax=vmax)
    colors = cm.Blues(norm(sim.infl_2d))

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    xc = sim.xc
    yc = sim.yc
    topo = sim.zc

    ax.plot_surface(xc, yc, topo, facecolors=colors,
                    rstride=1, cstride=1, linewidth=0,
                    antialiased=False,
                    shade=False)

    ax.view_init(25, 295)

    return fig


"""
Surface plots with a second plot on the horizontal plane
"""


def plot_3D_infl_veg(sim, infl_2d=None, trim_at_outlet=2):
    """
    Infiltration plot on the incline, and veg on the horizontal plane
    """
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(111, projection='3d')

    ax = fix_3D_axes(ax)
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    norm = plt.Normalize()
    isveg = sim.veg.copy()[:, trim_at_outlet:]
    isveg = isveg.astype(float)

    isveg[isveg == 0] = 0.1
    isveg[isveg == 1] = 0.8
    isveg[0, -1] = 0
    isveg[1, -1] = 1
    veg_colors = cm.Greens(norm(isveg))

    xc = sim.xc[:, trim_at_outlet:]
    yc = sim.yc[:, trim_at_outlet:]
    topo = sim.zc[:, trim_at_outlet:]

    if not infl_2d:
        infl_2d = sim.infl_2d

    norm = plt.Normalize(vmin=infl_2d.min(), vmax=infl_2d.max())
    colors = cm.Blues(norm(infl_2d[:, trim_at_outlet:]))
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    ax.plot_surface(xc, yc, topo, facecolors=colors,
                    rstride=1, cstride=1, linewidth=0,
                    antialiased=False, shade=False)

    ax.plot_surface(xc, yc, topo * 0, facecolors=veg_colors,
                    rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False,
                    alpha=0.8)

    ax.view_init(25, 295)

    return fig


def plot_h_surface(sim, scale=10):
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax = fix_3D_axes(ax)

    norm = plt.Normalize()

    veg = sim.veg.copy()
    veg[-1, 0] = -0.3
    veg[-1, -1] = 1.5
    veg_colors = cm.Greens(norm(veg))

    h_norm = colors.Normalize(vmin=10 * sim.hc.ravel().min() - .01,
                              vmax=10 * sim.hc.ravel().max())
    h_colors = cm.PuBu(h_norm(sim.hc[sim.i_tr] * 10.))

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.plot_surface(sim.xc, sim.yc + 1, sim.zc * 0, facecolors=veg_colors,
                    rstride=1, cstride=1, alpha=0.9,
                    linewidth=0, antialiased=True, shade=False)

    ax.plot_surface(sim.xc, sim.yc + 1, sim.zc + sim.hc[sim.i_tr] * scale,
                    facecolors=h_colors, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False, alpha=0.8)

    ax.view_init(25, 285)
    return fig, ax


def plot_U_surface(sim, scale=20):
    sim['umag'] = np.sqrt(sim.uc ** 2 + sim.vc ** 2)

    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax = fix_3D_axes(ax)

    norm = plt.Normalize()

    veg = sim.veg.copy()
    veg[-1, 0] = -0.3
    veg[-1, -1] = 1.5
    veg_colors = cm.Greens(norm(veg))

    h_norm = colors.Normalize(vmin=10 * sim.umag.ravel().min() - .01,
                              vmax=10 * sim.umag.ravel().max())
    h_colors = cm.Blues(h_norm(sim.umag[sim.i_tr] * 10.))

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.plot_surface(sim.xc, sim.yc + 1, sim.zc * 0, facecolors=veg_colors,
                    rstride=1, cstride=1, alpha=0.9,
                    linewidth=0, antialiased=True, shade=False)

    ax.plot_surface(sim.xc, sim.yc + 1, sim.zc + sim.hc[sim.i_tr] * scale,
                    facecolors=h_colors, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False, alpha=0.8)

    ax.view_init(25, 285)
    return fig, ax


def plot_3D_fhU(sim, h_scale=None, ind=None):
    """

    """
    if not h_scale:
        h_scale = int(sim.zc.max() / sim.hc.max() / 3.)
    if not ind:
        ind = sim.i_tr

    fig = plt.figure(figsize=(11, 6))

    ax = fig.add_subplot(111, projection='3d', )
    ax = fix_3D_axes(ax)

    veg = sim.veg.copy()
    veg[0, -1] = -0.2
    veg[1, -1] = 1.5

    norm = plt.Normalize()
    veg_colors = cm.Greens(norm(veg))

    f_norm = colors.Normalize(vmin=sim.infl_2d.ravel().min() - .01,
                              vmax=sim.infl_2d.ravel().max())

    f_colors = cmocean.cm.deep(f_norm(sim.infl_2d))

    U = np.sqrt(sim.uc ** 2 + sim.vc ** 2)
    U_norm = colors.Normalize(vmin=10 * U[ind].ravel().min() - .01,
                              vmax=10 * U[ind].ravel().max())

    U_colors = cm.Blues(U_norm(U[ind] * 10.))

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    ax.plot_surface(sim.xc, sim.yc, sim.zc, facecolors=veg_colors,
                    rstride=1, cstride=1, linewidth=0, antialiased=True,
                    shade=False, alpha=0.5)

    ax.plot_surface(sim.xc, sim.yc, sim.zc * 0, facecolors=f_colors,
                    rstride=1, cstride=1, linewidth=0, antialiased=True,
                    shade=False, alpha=0.8)

    ax.plot_surface(sim.xc, sim.yc, sim.zc + sim.hc[ind] * h_scale,
                    facecolors=U_colors, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False, alpha=0.8)

    ax.view_init(25, 295)
    return fig, ax
