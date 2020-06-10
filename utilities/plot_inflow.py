# coding=utf-8
from utilities.plot_functions import *


# def plot_inflow(sim,  t_f = None, N_profile=5):
#     """
#     Plot inflow profile at select times


#     """
#     fig, ax = plt.subplots(1, figsize=(7, 4.5))

    
#     if t_f:
#         inds = np.where(sim.t_print < t_f)[0]
#     else:
#         inds = np.where(np.diff(sim.hc.mean(1).mean(1)) > 1e-6)[0]

#     freq = int(len(inds) / N_profile)

#     for ind in inds[freq-1::freq]:
#         ax.plot(sim.xc.mean(0), sim['hc'].mean(1)[ind], next(linecycler),
#                 label=np.round(sim.t_print[ind] / 60., 1))

#     legend = ax.legend(title="time (min)")

#     plt.setp(legend.get_title(),fontsize='medium')
#     plt.xlabel('x')
#     plt.ylabel("h (cm)")



# def plot_all_inflow(sim, t_f = None, freq=3, label_axes=True):
#     """
#     Plot inflow profile at select times

#     Parameters:
#     -----------
#     sim : dict
#         SVE simulation
#     freq : int
#         plot frequency
#     """

#     fig, ax = plt.subplots(1, figsize=(7, 4.5))
#     plt.subplots_adjust(wspace=0.2)

#     scale = 3.6e5/sim.Lx
#     sim['bad'] = sim['hc'] * sim['uc']*scale
#     sim['qc'] = sim['xflux1']/sim.dx*scale

#     if t_f:
#         inds = np.where(sim.t_print < t_f)[0]
#     else:
#         inds = np.where(np.diff(sim.hc.mean(1).mean(1)) > 1e-6)[0]

#     for i in inds[::freq]:
#         ax.plot(sim.xc.mean(0),  sim['qc'].mean(1)[i], 'b--')
#         #ax.plot(sim.xc.mean(0),  sim['bad'].mean(1)[i], 'r--')

#     if label_axes:
#         ax.set_xlabel('x')
#         ax.set_ylabel("q (cm/hr)")
#     else:
#         ax.set_xticklabels("")
#         ax.set_yticklabels("")
#     ax.set_xlim(0, )
#     return fig, ax
