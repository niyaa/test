# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 18:04:50 2016

@author: niya
"""
import numpy as np;
class Ibc:
    
    a=0;
    def __init__(self):
        self.a=0;
        
    def IsentropicVortexIC2D(self,x,y,time):
        xo=5;yo=0;beta=5;gamma=1.4;
        rho=1;u=1;v=0;p=1;
        x=np.array(x);y=np.array(y);
        xmut=x-u*time;
        ymvt=y-v*time;
        
        r=np.sqrt((xmut-xo)**2+(ymvt-yo)**2);
        
        u=u-beta*np.exp(1-r**2)*(ymvt-yo)/(2*np.pi);
        v=v+beta*np.exp(1-r**2)*(xmut-xo)/(2*np.pi);
        rho1=(1-((gamma-1)*beta**2*np.exp(2*(1-r**2))/(16*gamma*np.pi*np.pi)))**(1/(gamma-1.0));
        p1=rho1**gamma;
        
        [t1,t2]=np.shape(rho1);
        Q=np.zeros((t1,t2,4),dtype=np.float);
        Q[:,:,0]=rho1;Q[:,:,1]=rho1*u;
        Q[:,:,2]=rho1*v;Q[:,:,3]=p1/(gamma-1)+0.5*rho1*(u**2+v**2);
        
        
  
        
        return Q;
        
        
    def IsentropicVortexBC2D(self,xin,yin,nxin,nyin,mapI,mapO,mapW,mapC,Q,time):
        Qbc=self.IsentropicVortexIC2D(xin,yin,time);
        mapB=np.hstack((mapI,mapO,mapW));
        for n in range (0,4):
            Qn=Q[:,:,n];
            Qbcn=Qbc[:,:,n];
            Qn[mapB]=Qbcn[mapB];
            Q[:,:,n]=Qn;
            
        return Q;
        

        
        