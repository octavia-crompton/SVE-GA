"""ravel_fxns_RF.py
Functions to create RF feature arrays!

The function `RF_patterns` is used by `base_ravel.py` to create a freature array `ravel`,  which is used by `fit_RF.py` to train the random forests.

`RF_patterns` is also used to apply a trained RF to a given vegetation
pattern. 
"""
import numpy as np
import pandas as pd
import scipy as sp
import scipy.ndimage
from scipy import ndimage
from scipy.ndimage.filters import gaussian_filter
import sys

np.seterr(divide='ignore', invalid='ignore')

def main(argv):
    pass 

def RF_patterns(isveg, ravel_params):        
    """
    Creates a dictionary containing feature arrays 
    
    Parameters
    ----------
    isveg: array_like
        binary vegetation array of shape (nrow, ncol)
    
    ravel_params: dict
    
    Returns 
    -------
    pattern_dict: dict
       dictionary of feature maps of shape (nrow, ncol)
    
    """

    isvegc = np.array(isveg, dtype = float) 
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    dx = ravel_params['dx']
    
    xc = np.arange(0, nrow*dx, dx)  + dx/2
    yc = np.arange(0, ncol*dx, dx)  + dx/2
    xc, yc = np.meshgrid(xc, yc)    
    xc = xc.T
    yc = yc.T
    
    edge = int(ravel_params['edge'])    
    saturate = int(ravel_params['saturate'])        
    gsigma = ravel_params['gsigma']
    
    bsigma_scale = ravel_params['bsigma_scale']
    upslopeLs = get_upslopeLs(ravel_params['window'])
    bsigma = get_bsigma(ravel_params['blur'])

    if gsigma == 0:
        gsigma = []
    elif type(gsigma) == int:
        gsigma = [gsigma]

    d2uB = func_d2uB(isvegc, edge,0)
    d2dB = func_d2dB(isvegc, edge, 0)     
    
    d2lB = func_d2lB(isvegc, edge)
    d2rB = func_d2rB(isvegc, edge)    
    d2xB = np.minimum(d2lB, d2rB)    
    d2xB[d2xB > saturate] = saturate
                  
    d2uV = func_d2uB(1-isvegc, edge,0) 
    d2dV = func_d2dB(1-isvegc, edge,0)
    d2lV = func_d2lV(isvegc, edge)
    d2rV = func_d2rV(isvegc, edge)   
    d2xV = np.minimum(d2lV, d2rV)
                  
    patchL,patchLB = get_patchL(isvegc, saturate) 
    bareL, bareLV = get_bareL(isvegc, saturate) 

    # assemble base pattern_dict
    pattern_dict = {'isvegc' : isvegc,                  
                    'd2uB' : d2uB, 
                    'd2dB' : d2dB, 
                    'd2xB' : d2xB,
                    'd2uV' : d2uV, 
                    'd2dV' : d2dV, 
                    'd2xV' : d2xV,
                    'patchLB' : patchLB,
                    'bareLV' : bareLV
                  }

    pos_sigma = [gs for gs in gsigma if gs > 0]
    
    for gs in pos_sigma:           
        for key in ['d2uV','d2dV','d2xV','bareLV']:  
            # smooth bare features
            pattern_dict[key + '_s{0}'.format(gs)] =   smoothB(pattern_dict[key], isvegc, gs)        
        
        for key in ['d2uB', 'd2dB', 'd2xB','patchLB']:
            # smooth vegetion features
            pattern_dict[key + '_s{0}'.format(gs)] = smoothV(pattern_dict[key], isvegc, gs)        
            
    if 'd6' in ravel_params['window']:
     upslopeLs = [6, 12, 24]         
    elif 'd4' in ravel_params['window']:
     upslopeLs = [4,8,12,16,20,24] 
    elif 'vary' in ravel_params['window']:
     upslopeLs = [2,4,6,8,10,12,14,16,20, 24,30]         
    elif  ravel_params['window'] == 'd2':
     upslopeLs = list(np.arange(2,30))    
    elif 'd2p' in ravel_params['window']:
     upslopeLs = list(np.arange(2,60))   
    elif 'd1' in ravel_params['window']:
     upslopeLs = list(np.arange(2,100))                       
    else: 
     upslopeLs = []           

    # compute upslope
    for L in upslopeLs:         
        upslopeL =  upslope_memory(isvegc,  min(ncol, int(L)))
        pattern_dict['upslope{0}'.format(L)] = upslopeL.copy()        
        pattern_dict['upslope{0}s'.format(L)] = gaussian_filter(upslopeL.astype(float), sigma=2)    

    for bs in bsigma: 
        blurred = gaussian_filter(isvegc.astype(float), sigma=(bs,bs))
        pattern_dict['blurred' + str(bs)] = blurred
        
        for bs_scale in bsigma_scale: 
            if bs_scale > 1:
                blurred_ani = gaussian_filter(isvegc.astype(float), 
                    sigma=(bs,bs*bs_scale))
                pattern_dict['blurred{0}_{1}'.format(bs, bs*bs_scale)] = \
                    blurred_ani 

                blurred_ani = gaussian_filter(isvegc.astype(float), 
                    sigma=(bs*bs_scale, bs))
                pattern_dict['blurred{0}_{1}'.format(bs*bs_scale,bs)] = \
                    blurred_ani 

            elif bs_scale == 1:
                pattern_dict['bdiff' + str(bs) + "_" + str(bs*bs_scale)] = blurred - isvegc.astype(float)

    # delete unsmoothed features if 0 not in gsigma
    if 0 not in gsigma:
        for key in ['d2uV','d2dV','d2xV','bareLV',  \
                    'd2uB', 'd2dB', 'd2xB','patchLB']: 
            del pattern_dict[key] 
        
        for L in upslopeLs: 
            del pattern_dict['upslope{0}'.format(L)] 
            
    pattern_dict['d2divide'] = ncol - yc/dx # in grid cells (not m)        
    pattern_dict['d2side'] = np.abs(nrow/2. - xc/dx )
    
    for key in pattern_dict.keys():
        dum = pattern_dict[key]
        dum[dum>saturate] = saturate
        pattern_dict[key] = dum
    
    return pattern_dict  

def smoothB(U, isvegc, gsigma):
    """
    Parameters
    ----------
    Returns 
    -------
    """
    U = U.astype(float)
    U[isvegc == 1] = np.nan
    V=U.copy()
    V[U!=U]=0
    VV=sp.ndimage.gaussian_filter(V,gsigma)

    W=0*U.copy()+1
    W[U!=U]=0
    WW=sp.ndimage.gaussian_filter(W,gsigma)

    Z=VV/WW
    Z = Z.astype(int)
    Z[isvegc ==1] = 0
    return Z

def smoothV(U, isvegc, gsigma):
    """
    Parameters
    ----------
    Returns 
    -------
    """
    U = U.astype(float)
    U[isvegc == 0] = np.nan
    V=U.copy()
    V[U!=U]=0
    VV=sp.ndimage.gaussian_filter(V,gsigma)

    W=0*U.copy()+1
    W[U!=U]=0
    WW=sp.ndimage.gaussian_filter(W,gsigma)

    Z=VV/WW
    Z = Z.astype(int)
    Z[isvegc ==0] = 0
    return Z
    
def func_d2uB(isvegc, edge, mask=0 ):
    """
    Computes the distance to nearest upslope bare cell    

    Parameters
    ----------
    Returns 
    -------
    
    Distane to nearest upslope bare cell
    =  0 for bare ground
    =  1 for veg cells with a neighboring bare cell upslope
    >  1 for veg cells with bare cells further upslope
    """

    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    arr = np.ones([nrow+ 2*edge, ncol + 2*edge]).T
    arr[edge:-edge, edge:-edge] = np.flipud(isvegc.T)
    
    if type(mask)!=int:
        newmask = np.ones([nrow+ 2*edge, ncol + 2*edge]).T
        newmask[edge:-edge, edge:-edge] = np.flipud(mask.T)  
        arr[newmask==0]= 1

    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = np.flipud(df1).T
    df1[isvegc == 0] = 0
    
    return df1


def func_d2dB(isvegc, edge, mask = 0):
    """
    Computes the distance to nearest downslope bare cell  

    Parameters
    ----------
    Returns 
    -------
    d2dB : arrayl_list
        array with distances to the nearest downslope bare cell
        ( 0 for bare ground,  1 for veg cells with a bare cell immediately downslope)
    """    
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    
    arr = np.ones([nrow+ 2*edge, ncol + 2*edge]).T
    arr[edge:-edge, edge:-edge] = isvegc.T

    if type(mask)!=int:
        newmask = np.ones([nrow+ 2*edge, ncol + 2*edge]).T
        newmask[edge:-edge, edge:-edge] = mask.T
        arr[newmask==0]= 1

    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = df1.T
    df1[isvegc == 0] = 0    
     
    return df1

def func_d2lB(isvegc, edge):
    """
    
    Parameters
    ----------
    Returns 
    -------
      d2lB : [nrow x ncol] array of distane to nearest left bare
        =  0 for bare cells
        =  1 for veg cells with a bare cell immediately left
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    arr = np.ones([nrow+ 2*edge, ncol + 2*edge])
    arr[edge:-edge, edge:-edge] = isvegc
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1[isvegc == 0] = 0
    
    return df1
        
def func_d2lV(isvegc, edge):
    """
    input: 
      isvegc : [nrow x ncol] array of vegetation field

    output : 
      d2lV : [nrow x ncol] array, distane to nearest veg cell to the left
        =  0 for veg cells
        =  1 for bare cells with a veg cell immediately left  
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    arr = np.ones([nrow+ 2*edge, ncol + 2*edge])
    arr[edge:-edge, edge:-edge] = 1 - isvegc
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1[isvegc == 1] = 0

    return df1

def func_d2rB(isvegc, edge):
    """
    input:
      isvegc : [nrow x ncol] array of vegetation field

    output :
      d2rB : [nrow x ncol] array;  distane to nearest bare cell to right
        =  0 for bare cells
        =  1 for veg cells with a veg cell immediately to right
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    arr = np.ones([nrow+ 2*edge, ncol + 2*edge])
    arr[edge:-edge, edge:-edge] = np.flipud(isvegc)
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = np.flipud(df1)
    df1[isvegc == 0] = 0
          
    return df1
    
def func_d2rV(isvegc, edge):
    """
    input: 
      isvegc : [nrow x ncol] array of vegetation field

    output : 
      d2rV : [nrow x ncol] array;  distane to nearest veg cell to right
        =  0 for veg cells
        =  1 for bare cells with a veg cell immediately to right  
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    arr = np.ones([nrow+ 2*edge, ncol + 2*edge])
    arr[edge:-edge, edge:-edge] = 1- np.flipud(isvegc)
    a = pd.DataFrame(arr) != 0

    df1 = a.cumsum()-a.cumsum().where(~a).ffill().fillna(0).astype(int)
    df1 = np.array(df1)[edge:-edge, edge:-edge]
    df1 = np.flipud(df1)
    df1[isvegc == 1] = 0

    return df1
 

def get_patchL(isvegc, saturate):
    """
    compute patch lengths

    Parameters:
    ----------
    isvegc : array_like
    saturate : int

    Returns:
    --------
    patchLv:  array_like
        vegetated patch length
    patchLB:  array_like 
        upslope interspace patch length (paired to veg patch)
    
    Notes:
    -----
    Ldict:  dict
        dictionary of veg patch lengths. 
        keys are downslope patch coordinates
    Bdict:  dict
        dictionary of paired upslope interspace lengths. 
  
    Usage:
    -----
    patchLv,patchLB,patchLc,Ldict,Bdict = get_patchL(isvegc)
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    patchLv = np.zeros(isvegc.shape, dtype = float)  # veg patch length
    patchLB = np.zeros(isvegc.shape, dtype = float)  # upslope interspace patch length (paired to veg patch)
    
    for i in range(nrow):  # loop over across-slope direction first
        count = 0           
        for j in range(ncol):    
            if isvegc[i, j] == 1:    #  if veg patch, add 1
                if j >= (ncol -1):  # if we're at the top of the hill                  
                  patchLv[i, j-count:] = count  # record veg patch length                  
                count += 1  
                                                        
            # if [i,j] is bare and the slope cell is vegetated, record.
            # each patch starts at [i,j-count] and ends at [i,j-1]
            elif isvegc[i,j] == 0 and isvegc[i, j-1] == 1:   
                if j > 0:
                  # veg patch starts at j-count and ends at j
                  patchLv[i, j-count:j] = count
                  try:
                      # find the nearest upslope veg cell
                      Lb = np.where(isvegc[i,j:] == 1)[0][0]                               
                      patchLB[i,j-count:j] = Lb
                  except IndexError:  # bare patch extends to top of hill
                      patchLB[i,j-count:j] = ncol - j
                  count = 0 
    patchLv[patchLv > saturate] = saturate
    patchLB[patchLB > saturate] = saturate
        
    return  patchLv, patchLB

def get_bareL(isvegc, saturate, skipflag = 0):
    """
    input : isvegc from get_source(df)  
    
    output : 
      bareL:  bare patch length
      bareLV:  upslope vegeted patch length (paired to bare patch)
      
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    bareLV = np.zeros(isvegc.shape, dtype = float)  # veg patch length
    bareL = np.zeros(isvegc.shape, dtype = float)  # upslope interspace patch length (paired to veg patch)
    
    for i in range(nrow):  # loop over across-slope direction first
        count = 0           
        for j in range(ncol):    
            if isvegc[i, j] == 0:    #  if bare, add 1
                if j >= (ncol -1):  # if we're at the top of the hill                  
                  bareL[i, j-count:] = count  # record bare length                  
                count += 1  
                                                        
            # if [i,j] is veg and the slope cell is bare, record.
            # each patch starts at [i,j-count] and ends at [i,j-1]
            elif isvegc[i,j] == 1 and isvegc[i, j-1] == 0:   
                if j > 0:
                  # veg patch starts at j-count and ends at j
                  bareL[i, j-count:j] = count
                  if skipflag == 1 and j-count ==0:
                      bareL[i, j-count:j] = 0
                  try:
                      # find the nearest upslope bare cell
                      Lb = np.where(isvegc[i,j:] == 0)[0][0]                               
                      bareLV[i,j-count:j] = Lb
                  except IndexError:  # bare patch extends to top of hill
                      bareLV[i,j-count:j] = ncol - j
                  count = 0 
    bareLV[bareLV > saturate] = saturate
    bareL[bareL > saturate] = saturate
        
    return  bareL, bareLV

def upslope_memory(isvegc,  memory = 3):
    """
    
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    
    dum = isvegc.copy()
    memory = int(memory)
    for k in range(int(ncol - memory)):
        dum[:, k] = isvegc[:, k:k+memory].sum(1)
    for k in range(1,memory+1):    
        dum[:, -k] = isvegc[:, -k:].sum(1)
    # dum[isvegc == 0] = 0
    return dum
    
                    
def mean_patch_len(array):
    """
    Given patchL or bareL, computes the mean along-slope length,
        but weights for multiple counts (i.e. only count each strip once.)
    """
    dum = array[ array > 0] 
    c = {}
    for n in np.unique( array[ array > 0]):
        c[n]= np.round((sum(dum == n)/n))
    L = 0
    for k in c.keys():
        L += k*c[k]
    return  L/np.sum(c.values())
     
def func_d2wB(isvegc, saturate, weight):
    """
    Distane to weighted nearest upslope bare cell
    =  0 for bare ground
    =  1 for veg cells with a neighboring bare cell upslope
    >  1 for veg cells with bare cells further upslope
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    
    res =  isvegc.copy()
    
    for i in range(ncol):
        d = isvegc.copy()
        d[:,:i+1] = 1
        res[:,i] = ndimage.distance_transform_edt(d, sampling = (weight, 1))[:, i]   
    res[isvegc ==0] = 0
    
    res[res>saturate] = saturate   
    return res

def func_d2wV(isvegc, saturate, weight):
    """
    Distane to weighted nearest upslope bare cell
    =  0 for bare ground
    =  1 for veg cells with a neighboring bare cell upslope
    >  1 for veg cells with bare cells further upslope
    """
    nrow = isvegc.shape[0]
    ncol = isvegc.shape[1]
    
    res =  isvegc.copy()
    
    for i in range(ncol):
        d = 1-isvegc.copy()
        d[:,:i+1] = 1
        res[:,i] = ndimage.distance_transform_edt(d, sampling = (weight, 1))[:, i]   
    res[isvegc ==1] = 0
    
    res[res>saturate] = saturate   
    return res 

def get_upslopeLs(window):
    """
    """
    if 'd6' in window:
        upslopeLs = [6, 12, 24]         
    elif 'd4' in window:
        upslopeLs = [4,8,12,16,20,24] 
    elif 'vary' in window:
        upslopeLs = [2,4,6,8,10,12,14,16,20, 24,30]         
    elif  window == 'd2':
        upslopeLs = list(np.arange(2,30))    
    elif 'd2p' in window:
        upslopeLs = list(np.arange(2,60)) 
    elif 'd21' in window:
        upslopeLs = list(np.arange(2,100,2))  
    elif 'd23' in window:
        upslopeLs = list(np.arange(2,100,3))   
    elif 'd123' in window:
        upslopeLs = [2,3,4,5] + list(np.arange(6,100,3))                 
    elif 'd1' in window:
        upslopeLs = list(np.arange(2,100))                       
    else: 
        upslopeLs = []

    return upslopeLs    

def get_bsigma(blur):
    """
    """
    if "1p" in blur:
        bsigma = [1,3,5,10,20,40]
    elif "3p" in blur:
        bsigma = [3,5,10,20,40]
    elif "12p" in blur:
        bsigma = [1,3,5,7,10,20,40]          
    elif "8p" in blur:
        bsigma = [1,3,5,7,10,15,20,40]            
    else:
        bsigma = []    
    return bsigma  
              

if __name__ == '__main__':
    main(sys.argv)
