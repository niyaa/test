# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:31:33 2016

@author: niya
"""


def CubVolMesh2D(CubOrd):
    [Ig.cR,Ig.cS,Ig.cW,Ig.cN] = Cubature2D(CubOrd);
     
    Ig.cVm=bb.InterpMatrix2D(cR,cS);
     
    [Ig.cDr,Ig.cDs]=bb.Dmatrices2D(gg.N,cR,cS,gg.V);
    
    [Ig.cRx, Ig.cSx, Ig.cRy, Ig.cSy, Ig.cJ]=bb.GeometricFactors2D(gg.x,gg.y,cDr,cDs);
    
    Ig.cChol= np.zeros((gg.Np,gg.Np,gg.K),dtype=np.int);
    
    for k in range (0,K):
        Ig.cChol[:,:,k]=cVm.T*np.diag((cJ[:,k]*cW))*cVm;
        
    cw=cW;  cW=np.matrix(cW)*np.ones((1,gg.K),dtype=np.int);
    ##check cW type
    cW= cW*cJ;
    cX=cVm*gg.x;cY=cVm*gg.y;
    
    return 0;

def Cubatrure2D(Corder):
    import integData;
    if(order<=28):
        temp=integData.Cub[order];
        cR= temp[:][0];
        cS= temp[:][1];
        cW= temp[:][2];
    else:
        cubNA= np.ceil((order+1)/2);
        [cA,cWA]=aa.JacobiGQ(0,0,cubNA-1);
        cubNB= np.ceil((order+1)/2);
        [cB,cWB]=aa.JacobiGQ(1,0,cubNB-1);
        
        cA=np.ones((np.size(cB),1),dtype=np.int)*np.matrix(cA);
        cB= np.matrix(cB).T*np.ones((1,cubNA),dtype=np.int);
        
        cA=np.array(cA);
        cB=np.array(cB);
        cR = 0.5*(1+cA)*(1-cB)-1;
        cS = cB;
        cW = 0.5*np.matrix(cWB).T*np.matrix(cWB);
        
        cR=np.ravel(cR,order='F')[:];
        cS=np.ravel(cS,order='F')[:];
        cW=np.ravel(cW,order='F')[:];
    Ncub = np.size(cubW);
    return cR, cS,cW,Ncub