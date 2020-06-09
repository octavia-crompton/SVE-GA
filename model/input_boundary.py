import numpy as np
import scipy as sp
import os, sys


def write_boundary(path, params):
    """
    Writes `boundary.dat` for a rectangular domain 
    """
    nrow = params['nrow']
    ncol = params['ncol']
    itype3 = params['itype3']
    itype1 = params['itype1']
    itype2 = params['itype2']
    itype4 =  params['itype4']

    inum = np.zeros([nrow+1, ncol+1], dtype = int)
    inum[1:, 1] = 1
    inum[1:, -1]= 1
    inum[1, 1:] = 1
    inum[-1, 1:] = 1
    inum[1, 1] = 2
    inum[1, -1] = 2
    inum[-1, -1] = 2
    inum[-1, 1] = 2

    ipos = np.zeros( [nrow+1, ncol+1, 2], dtype = int)

    # bottom boundary
    ipos[2:-1, 1,0] = 1
    ipos[1, 1,1] = 1
    ipos[-1, 1,1] = 1

    # right boundary
    ipos[-1, 1:-1, 0] = 2
    ipos[-1, -1,1] = 2

    # left boundary
    ipos[1, 1:, 0] = 4

    # top boundary
    ipos[2:, -1,0] = 3
    ipos[1, -1,1] = 3

    itype = np.zeros([nrow+1, ncol+1, 2], dtype = int)
    # bottom boundary
    itype[2:-1, 1,0] = itype1
    itype[1, 1,1] = itype1
    itype[-1, 1,1] = itype1

    # right boundary
    itype[-1, 1:-1, 0] = itype2
    itype[-1, -1,1] = itype2

    # left boundary
    itype[1, 1:,0] = itype4

    # top boundary
    itype[2:, -1,0] = itype3
    itype[1, -1,1] = itype3

    nbcell = 2*nrow + 2*ncol - 4  # number of boundary cells

    fname = '{0}/input/boundary.dat'.format(path)
    f = open(fname, 'w')
    f.write('number of boundary cell \n')
    f.write('  {0} \n'.format(nbcell))
    f.write(' j    k          inum    itype             ipos \n')
    j = 1
    for k in range(1, ncol+1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1],
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    for j in range(2, nrow+1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1],
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    for k in range(ncol-1,0,-1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1],
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    for j in range(nrow-1,1,-1):
        if inum[j, k] == 2:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<8} {4:<9} {5:<8} {6:<6} \n'.format(
                        j, k, inum[j, k], itype[j, k, 0], itype[j, k, 1],
                         ipos[j, k, 0], ipos[j, k, 1]))
        else:
            f.write( '{0:<5} {1:<13} {2:<7} {3:<18} {4:<10}   \n'.format(
                         j, k, inum[j, k],  itype[j, k, 0],  ipos[j, k, 0], ))

    kbeg = np.ones(nrow+1, dtype = int)
    kend = np.ones(nrow+1, dtype = int)*ncol

    f.write('nrow\n')
    f.write("{0}\n".format(nrow))
    f.write('ncol\n')
    f.write("{0}\n".format(ncol))
    f.write('j     kbeg          kend \n')
    for j in range(1, nrow+1):
        f.write( '{0:>5}  {1:>5} {2:>13}   \n'.format(
                    j, kbeg[j],kend[k] ))

    f.close()
    boundary_fix(path, params)

    return inum, ipos, itype

def boundary_fix(path, params):
    """
    Writes subcritical fixed-flux boundary condition

    Parameters:
    -----------
    path : str
        Location to save boundary.dat

    Notes:
    -----
    q1 : subcritical boundary flux in m^2/s

    """
    fname = '{0}/input/boundary.dat'.format(path)
    f = open(fname, 'a')
    nrow = params['nrow']
    ncol = params['ncol']
    if params['itype1'] == 4:
        fixj = np.arange(1, nrow+1)
        fixk = np.ones(nrow, dtype = int)
        fixh = np.zeros(nrow, dtype = float )
        fixu = np.zeros(nrow, dtype = float)
        q1 = params['q1']
        fixv = np.ones(nrow, dtype = float)*q1
        ndir = len(fixj)
    else:
        ndir = 0

    f.write('number of fixed BC cells, ndir    \n ')
    f.write('{0}  \n '.format(ndir))
    f.write('j     k    fix h    fix v    fix u	\n')
    for i in range(ndir):
        f.write('{0}     {1}     {2}     {3}     {4}     \n'.format(
                fixj[i],fixk[i],fixh[i],fixu[i],fixv[i]))
    f.close()
