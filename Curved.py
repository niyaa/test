# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:22:12 2016

@author: niya
"""
from numpy import linalg as la
import aa;
import bb;
import Glob2 as gg;
import numpy as np;
from math import atan2;
import flux;
def CurvedEulerRHS2D(Q,time,SolutionBC,fluxtype):
    cQ=np.zeros(gg.cubN,np.K,4);
    for n in range (0,4):
        cQ[:,:,n]=gg.cVm*Q[:,:,n];
    gamma=1.4;
    [F,G,rho,u,v,p]=flux.EulerFluxes2D(cQ,gamma);
    
    rhsQ=zeros(gg.Np,gg.K,4);
    for n in range (0,4):
        ddr = (gg.cDr.T)*np.multiply(gg.cW,(np.multiply(gg.cRx,F[:,:,n]))+\
        +np.multiply(gg.cRy,G[:,:,n]));
        dds = (gg.cDs.T)*np.multiply(gg.cW,(np.multiply(gg.cSx,F[:,:,n]))+\
        +np.multiply(gg.cSy,G[:,:,n]));
        rhsQ[:,:,n]=ddr+dds;
    
    gg.nx=