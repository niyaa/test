# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:46:13 2016

@author: niyaa
"""

import Glob2 as gg;
import flux as ff;
import bb;
import numpy as np;

def Euler2D(Q,FinalTime,BC):
    #Initialize filter
    Filt =bb.CutOffFilter2D(N,0.95);
    #compute initial timestep
    gamma=1.4;
    dt=EulerDT2D(Q,gamma);
    time=0;tstep = 1;
    
    #storage for low storage RK time stepping
    rhsQ=0*Q;resQ=0*Q;
    
    #filter initial solution 
    for n in range (0,4):
        Q[:,:,n]=Filt*np.matrix(Q[:,:,n]);
    
    while(time<FinalTime):
        if(time+dt>FinalTime):
            dt=FinalTime - time;
            
        for INTRK in range (0,5):
            rhsQ=EulerRHS2D(Q,time,BC,tstep);
            for n in range (0,4):
                rhsQ[:,:,n]=Filt*np.matrix(rshQ[:,:,n]);
            resQ=gg.rk4a[INTRK]*np.matrix(resQ)+dt*np.matrix(rhsQ);
            Q=Q+gg.rk4b[INTRK]*np.matrix(resQ);
        
        time=time+dt;
        dt=EulerDT2D(Q,gamma);
        tstep=tstep+1;
        
    return Q;
    
def EulerDT2D(Q,gamma):
    rho=Q[:,:,0];
    rhou=Q[:,:,1];
    rhov=Q[:,:,2];
    Ener=Q[:,:,3];
    
    rho=rho[gg.vmapM];rhou=rhou[gg.vmapM];rhov=rhov[gg.vmapM];Ener=Ener[gg.vmapM];
    
    u=np.divide(rhou,rho);v=np.divide(rhov,rho);
    p=(gamma-1.0)*(Ener-rho*(u**2+v**2)/2);
    c=np.sqrt(np.abs(gamma*np.divide(p,rho)));
    
    dt=1/np.max(((gg.N+1)**2)*5*np.ravel(gg.Fscale,order='F')[:]*\
    (np.sqrt((np.ravel(u,order='F')[:])**2+(np.ravel(u,order='F')[:])**2)+\
    np.ravel(c,order='F')[:]));
    
    rhoprange=[np.min(rho),np.max(rho),np.min(p),np.max(p)];
    
    return dt;
    
def EulerRHS2D(Q,time,ExactSolutionBc,tstep):
    gg.vmapM=np.reshape(gg.vmapM,(gg.Nfp*gg.Nfaces,K),order='F');
    
    gg.vmapP=np.reshape(gg.vmapP,(gg.Nfp*gg.Nfaces,K),order='F');
    
    gamma=1.4;
    [F,G,rho,u,v,p]=flux.EulerFluxes2D(Q,gamma);
    
    for n in range (0,4):
        dFdr = gg.Drw*F[:,:,n];dFDs=gg.Dsw*F[:,:,n];
        dGdr = gg.Drw*G[:,:,n];dGds=gg.Dsw*G[:,:,n];
        rhsQ[:,:,n]=np.multiply(gg.rx,dFdr)+np.multiply(gg.sx,dFds)+\
        np.multiply(gg.ry,dGdr)+np.multiply(gg.sy,dGds);
        
    for n in range (0,4):
        Qn=Q[:,:,n];
        QM[:,:,n]=Qn[vmapM];QP[:,:,n]=Qn[vmapP];
        
    [fM,gM,rhoM,uM,vM,pM]=flux.EulerFluxes2D(QM,gamma);
    [fP,gP,rhoP,uP,vP,pP]=flux.EulerFluxes2D(QP,gamma);
    
    lamb = np.max(np.sqrt(uM**2+vM**2)+np.sqrt(np.abs(gamma*np.divide(pM,rhoM))),\
    np.sqrt(uP**2+vP**2)+np.sqrt(np.abs(gamma*np.divide(pP,rhoP))));
    
    lamb = np.reshape(lamb,(gg.Nfp,gg.Nfaces*gg.K));
    lamb = np.ones((gg.Nfp,1),dtype=np.float)*max
    lamb = np.reshape(lamb,(gg.Nfp*gg.Nfaces,gg.K));
    
    for n in range (0,4):
        nflux = np.multiply(nx,(fp[:,:,n]+fM[:,:,n]))+np.multiply(ny,(gp[:,:,n]+gM[:,:,n]))+\
        np.multiply(lamb,(QM[:,:,n]-QP[:,:,n]));
        rhsQ[:,:,n]=rhsQ[:,:,n]-gg.LIFT*(np.multiply(gg.Fscale,nflux/2));
        
        
    return 0;
    



def TimeStepEuler(Q,gamma):
    rho=gg.Q[:,:,0];rhou=Q[:,:,1];
    rhov=gg.Q[:,:,2];Ener=Q[:,:,3];
    rho=rho[gg.vmapM];rhou=rhou[gg.vmapM];rhov=rhov[gg.vmapM];Ener=Ener[gg.vmapM]
    u=np.divide(rhou,rho);
    v=np.divide(rhov,rho);
    
    p=(gamma-1.0)*(Ener-np.multiply(rho,(u**2+v**2))/2);
    c=np.sqrt(np.abs(gamma*np.divide(p,rho)));
    
    dt=1/max(((gg.N+1)**2)*5*np.ravel(gg.Fscale,order='F')*(np.sqrt(u**2+v**2)+c));
    
    rhoprange = [np.min(rho),np.max(rho),np.min(p),np.max(p)];
    
    return dt;
    
def EulerShock(N,mshFile):
    gg.N=3;
    [gg.Nv,gg.VX,gg.VY,gg.EToV,gg.K,gg.BCType]=bb.mshRead('SqWithCurv.msh'); 
    bb.Startup2D();
    Q=np.zeros((gg.Np,gg.K,1),dtype=np.float);
    
    [k,f]=np.where((gg.BCType.T==gg.Cyl));
    cylfaces=[k,f];
    curved=np.sort(np.unique(k));straight=np.setdiff1d(np.arange(0,gg.K),curved);
    ff.MakeCylinder2D(cylfaces,r,x0,y0);
    
        
    