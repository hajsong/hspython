from MITgcmutils import rdmds
import numpy as np

def loadgrid(gridpath, region=None, varname=None):
    """
    [Function]
    grd=loadgrid(gridname, region=None, varname=None):

    [Description]
    Read grid files and load them into the 'grd' object
    
    [Inputs]
    gridpath : A path to the grid files.
    region   : A list defining the boundary for the grid.
               [x0,x1,y0,y1] (default is [0,nx,0,ny])
    varname  : A list of strings for the grid data to load

    [Output]
    grd      : An object containing grid data.
               
    """

    if varname is None:
        varname = ['XC','YC','RAC','DXC','DYC','hFacC','hFacW','hFacS',\
                   'Depth','RC','RF','DRC','DRF','XG','YG','RAZ','DXG','DYG']

    class grd(object):
        for iv, vname in enumerate(varname):
            if region is None:
                exec('tmpvar = rdmds("'+gridpath+varname[iv]+'")')
                tmpvar = tmpvar.squeeze()
                exec(varname[iv]+' = tmpvar')
            else:
                if vname is 'RC' or vname is 'RF' or vname is 'DRC' or vname is 'DRF':
                    exec('tmpvar = rdmds("'+gridpath+varname[iv]+'")')
                else: 
                    exec('tmpvar = rdmds("'+gridpath+varname[iv]+\
                         '",region = '+str(region)+')')
                tmpvar = tmpvar.squeeze()
                exec(varname[iv]+' = tmpvar')
            if vname == 'hFacC':
                mskC = hFacC.copy()
                mskC[mskC==0] = np.nan
                mskC[np.isfinite(mskC)] = 1.
            if vname == 'hFacW':
                mskW = hFacW.copy()
                mskW[mskW==0] = np.nan
                mskW[np.isfinite(mskW)] = 1.
            if vname == 'hFacS':
                mskS = hFacS.copy()
                mskS[mskS==0] = np.nan
                mskS[np.isfinite(mskS)] = 1.
            del tmpvar
    del grd.iv,grd.vname

    return grd
