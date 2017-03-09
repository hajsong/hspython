import numpy as np

def calmld(grd,T,S,delrho=0.03,LK=50):
    """ 
    Computing MLD 
    NEED TO GENERALIZE THE CODE!!
    """
    import numpy as np
    from MITgcmutils import densjmd95
    [ly,lx]=grd.XC.shape

    # compute density
    rho=0*grd.hFacC;
    for k in xrange(LK):
      rho[k,:,:]=densjmd95(S[k,:,:],T[k,:,:],0);

    rho=rho*grd.mskC - 1000.0;

    # Compute MLD
    mld=np.zeros((800,2800))
    mld_dens = rho[0,:,:]+delrho;
    for ix in range(0,lx):
        for iy in range(0,ly):
            mdens = mld_dens[iy,ix];
            if np.isnan(mdens):
                mld[iy,ix]=float('nan');
            else:
                rhocol = rho[:,iy,ix];
                id=np.where(rhocol>=mdens)[0]
                if len(id)==0:
                    mld[iy,ix]=grd.Depth[iy,ix]
                elif len(id)==len(rhocol):
                    mld[iy,ix]=-grd.RC[0];
                else:
                    alpha=(mdens-rhocol[id[0]-1])/(rhocol[id[0]]-rhocol[id[0]-1]);
                    mld[iy,ix]=-(alpha*grd.RC[id[0]]+(1-alpha)*grd.RC[id[0]-1]);

    return mld 

def calc_conv(grd,fldU,fldV,xrng=None,yrng=None,ik=[0],list_factors=None):
    """
    object:    compute flow field convergent part (i.e. minus the divergence)
    inputs:    fldU and fldV are transport or velocity fields
    list_factors the list of factors that need to
                 be applied to fldU,fldV. By default it is empty (i.e. {}).
                 The most complete list wouldrbe ['dh','dz','hfac'].
    output:    fldDIV is the convergence
               (integrated, not averaged, over grid cell area)

    notes:     fldU,fldV  that may be either
               [A] a 3D vector field or
               [B] a 2D vector field

    in case [A], layer thicknesses = mygrid.DRF; in case [B] layer thickness = 1
    in any case, the global variable mygrid is supposed to be available
    """
    #
    #  grid size
    #
    [ly,lx]=grd.XC.shape;
    #
    #  prepare fldU/fldV:
    #
    fldU[np.isnan(fldU)]=0.; fldV[np.isnan(fldV)]=0.;
    #
    #  Apply grid factor if necessary
    #
    if list_factors!=None:
        facW=np.ones([Ly,Lx])
        facS=np.ones([Ly,Lx])
        for tmpstr in list_factors:
            if tmpstr=='dh': facW=facW*grd.DYG; facS=facS*grd.DXG;
            elif tmpstr=='dz': facW=facW*grd.DRF[ik]; facS=facS*grd.DRF[ik];
            elif tmpstr=='hfac':
                facW=facW*grd.hFacW[ik,...];
                facS=facS*grd.hFacS[ik,...];
        fldU=fldU*facW
        fldV=fldV*facS
    #
    #  if fldU is in 2D.
    #
    if len(fldU.shape)==2:
        FLDU=np.zeros((ly,lx+1));
        FLDU[:,:-1]=fldU;FLDU[:,-1]=fldU[:,-1]
        FLDV=np.zeros((ly+1,lx));
        FLDV[:-1,:]=fldV;FLDV[-1,:]=fldV[-1,:]
        fldDIV=np.zeros(fldU.shape);
        fldDIV=FLDU[:,0:-1]-FLDU[:,1:]+FLDV[0:-1,:]-FLDV[1:,:]
        if xrng!=None:
            fldDIV=fldDIV[yrng][:,xrng]    
    elif len(fldU.shape)==3:
        FLDU=np.zeros((len(ik),ly,lx+1));
        FLDU[:,:,:-1]=fldU;FLDU[:,:,-1]=fldU[:,:,-1];
        FLDV=np.zeros((len(ik),ly+1,lx));
        FLDV[:,:-1,:]=fldV;FLDV[:,-1,:]=fldV[:,-1,:];
        fldDIV=np.zeros(fldU.shape);
        fldDIV=FLDU[:,:,0:-1]-FLDU[:,:,1:]+FLDV[:,0:-1,:]-FLDV[:,1:,:]
        if xrng!=None:
            fldDIV=fldDIV[:,yrng][:,xrng]    
    
    return fldDIV
