# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 21:26:05 2016

@author: niya
"""

import numpy as np;
import glob1 as gg;
import aa;
def minmod(v):
    temp=np.shape(v);
    m=temp[0];
    mfunc=np.zeros((1,temp[1]),dtype=np.float);
    
    s=np.sum(np.sign(v),0)/m;
    ids= np.where(np.abs(s)==1);
    if(not all(ids)):
        mfunc[ids]=s[ids]*min(min(abs(v[:,ids])),1);
    
    return mfunc;


def EulerRHS1D(rho,rhou,Ener):
    gamma=1.4;
    temp=np.divide(1,rho);
    pres=(gamma-1.0)*(Ener-0.5*(rhou)**2*temp);
    cvel= np.sqrt(gamma*pres*temp);lm=abs(rhou*temp)+cvel;

    #copute fluxes
    rhof=rhou; 
    rhouf=rhou**2*temp+pres;
    Enerf=(Ener+pres)*rhou*temp;
    
    #Compute jumps at internal faces    
    drho=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.float);
    np.ravel(drho,order='F')[:]=rho(gg.vmapM) - rho(gg.vmapP);
    
    drhou=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.float);
    np.ravel(drhou,order='F')[:]=rhou(gg.vmapM) - rhou(gg.vmapP);

    
    dEner=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.float);
    np.ravel(dEner,order='F')[:]=Ener(gg.vmapM) - Ener(gg.vmapP);

    
    drhof=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.float);
    np.ravel(drhof,order='F')[:]=rhof(gg.vmapM) - rhof(gg.vmapP);

    
    
    drhouf=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.float);
    np.ravel(drhouf,order='F')[:]=rhouf(gg.vmapM) - rhouf(gg.vmapP);

    
    dEnerf=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.float);
    np.ravel(dEnerf,order='F')[:]=Enerf(gg.vmapM) - Enerf(gg.vmapP);

    
    LFc=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.float);
    np.ravel(LFc,order='F')[:]=max(lm(gg.vmapM) - lm(gg.vmapP));
    
    #Cmpute fluxes at interfaces
    


def Euler1D(rho,rhou,Ener,FinalTime):
    gamma=1.4;CFL=1.0;time=0;
    mindx=min(x[1,:],x[1,:]);
    return 


def SlopeLimitN(rho):
    #compute cell averages
    uh=gg.invV*u;
    uh[1:Np,:]=0;
    uavg = gg.V*uh;
    v=uavg[0,:];
    
    #Apply slope limiter as needed
    ulimit= u; eps0=1.0e-8;
    
    #find end values of each element
    ue1=u[0,:];ue2=u[-1,:];
    
    #find cell averages 
    vk=v; vkm1=np.array(v[0],v[0:gg.K-1]);
    vkp1= np.array(v[1:gg.K],v[gg.K-1]);
    
    ve1= vk - minmod(np.array([vk-ue1,vk-vkm1,vkp1-vk]));
    
    ve2= vk + minmod(np.array([ue2-vk,vk-vkm1,vkp1-vk]));
    
    ids = where(np.abs(ve1-ue1)>eps0|np.abs(ve2-ue2)>eps0);
    
    #Check to see if any elements require limiting
    
    if (not all(ids)):
        #create piecewise linear solution for limiting on specified elements
        uhl = gg.invV*u[:,ids];uhl[2:Np,:]=0;ul=V*uhl;
        ulimit[:,ids]=SlopeLimitLin(ul,x[:,ids],vkm1[ids],vk[ids],vkp[ids]);
   
    return ulimit;

def SlopeLimitLin(ul,xl,vm1,v0,vp1):
    #Apply slopelimited on linear function ul(Np,1) on x(Np,1)
    #(vm1,v0,vp1) are cell averages left, centre and right
    
    ulimit=ul;h=xl[gg.Np,:]-xl[0,:];
    x0=np.ones((Np,1),dtype=np.float)*(xl[0,:]+h/2);
    
    hN=np.ones((Np,1),dtype=np.int)*h;
    
    #limit function
    
    ux = np.divide(2,hN)*(Dr*ul);
    ulimit = np.ones((Np,1),dtype=np.int)*v0+(xl-x0)*(np.ones((Np,1),dtype=np.int)*minmod());
    return ulimit;    
    
    
    
def EulerDriver1D():
    gg.N=6;
    
    [gg.Nv,gg.VX,gg.K,gg.EToV]=aa.MeshGen1D(0.0,1.0,250);
    
    
   
    
  
    

    

    
    