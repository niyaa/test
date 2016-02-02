# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 16:17:06 2015

@author: niyaa
"""
import numpy as np;
from numpy import linalg as la
import aa;
import bb;
import Glob2 as gg;
import numpy as np;
from math import atan2;

#MakeCylinder2D(cylfaces,r,x0,y0); r=0.2; xo=2;y=0.5;
gg.N=3;
[gg.Nv,gg.VX,gg.VY,gg.EToV,gg.K,gg.BCType]=bb.mshRead('cir.msh'); 
bb.Startup2D();
Q=np.zeros((gg.Np,gg.K,1),dtype=np.float);
xo=1;yo=0.5;ra=0.5;
[k1,f1]=np.where((gg.BCType.T==gg.Cyl));
cylfaces=[k1,f1];
curved=np.sort(np.unique(k1));straight=np.setdiff1d(np.arange(0,gg.K),curved);
faces=np.array(cylfaces);
NCurveFaces = np.size(faces,1);
vflag = np.zeros((len(gg.VX),1),dtype=np.int);
for n in range (0,NCurveFaces):
    
    k=faces[1,n];f=faces[0,n];
    v1=gg.EToV[k,f];v2=gg.EToV[k,np.mod((f+1),gg.Nfaces)];
 
    theta1 = atan2(gg.VY[v1]-yo,gg.VX[v1]-xo);
    theta2 = atan2(gg.VY[v2]-yo,gg.VX[v2]-xo);
    
    newx1=xo+ra*np.cos(theta1);newy1=yo+ra*np.sin(theta1);
    newx2=xo+ra*np.cos(theta2);newy2=yo+ra*np.sin(theta2);
    
    gg.VX[v1]=newx1;gg.VX[v2]=newx2;gg.VY[v1]=newy1;gg.VY[v2]=newy2;
    
    vflag[v1]=1; vflag[v2]=1;
vflag =  vflag[np.ravel(gg.EToV,order='F')[:]];
vflag= np.reshape(vflag,(gg.K,gg.Nfaces),order='F');
ks=np.where(np.sum(vflag,1)>0);
ks=np.array(ks).T;
va=gg.EToV[ks,0].T; vb=gg.EToV[ks,1].T;vc=gg.EToV[ks,2].T;
temp=0.5*(-np.matrix(gg.r+gg.s).T*np.matrix(gg.VX[va])+np.matrix(1+gg.r).T*np.matrix(gg.VX[vb])+np.matrix(1+gg.s).T*np.matrix(gg.VX[vc]));    
for i in range (0,len(ks)):
    gg.x[:,ks[i]]=temp[:,i]
temp=0.5*(-np.matrix(gg.r+gg.s).T*np.matrix(gg.VY[va])+np.matrix(1+gg.r).T*np.matrix(gg.VY[vb])+np.matrix(1+gg.s).T*np.matrix(gg.VY[vc]));                        
for i in range (0,len(ks)):
    gg.y[:,ks[i]]=temp[:,i]

for n in range (0,NCurveFaces):
    k=faces[1,n];f=faces[0,n];
    if(f==0): v1=gg.EToV[k,0];v2=gg.EToV[k,1];vr=gg.r;
    if(f==1): v1=gg.EToV[k,1];v2=gg.EToV[k,2];vr=gg.s;

    if(f==2): v1=gg.EToV[k,0];v2=gg.EToV[k,2];vr=gg.s;
    fr=vr[gg.Fmask[:,f]];
    x1=gg.VX[v1]; y1=gg.VY[v1];
    x2=gg.VX[v2]; y2=gg.VY[v2];
    
    theta1=atan2((y1-yo),(x1-xo));theta2=atan2((y2-yo),(x2-xo));
    if((theta2>0)&(theta1<0)):theta1 = theta1 + 2*np.pi;
    if((theta2>0)&(theta1<0)):theta2 = theta2 + 2*np.pi;
    
    theta= 0.5*theta1*(1-fr) + 0.5*theta2*(1+fr);
    
    fdx=  xo + (ra*np.cos(theta)).reshape(len(theta),1)-gg.x[gg.Fmask[:,f],k];
    fdy=  yo + (ra*np.sin(theta)).reshape(len(theta),1)-gg.y[gg.Fmask[:,f],k];
    
    
    Vface = aa.Vandermonde1D(gg.N,fr);
    Vvol  = aa.Vandermonde1D(gg.N,vr);
    
    vdx=Vvol*(la.solve(Vface,fdx));vdy= Vvol*(la.solve(Vface,fdy));
    
    ids= np.where(abs(1-vr)>1e-7);
    
    if(f==0): blend=-np.divide((gg.r[ids]+gg.s[ids]),(1-vr[ids]));
    if(f==1): blend= np.divide((gg.r[ids]+1),(1-vr[ids]));
    if(f==2): blend= -np.divide((gg.r[ids]+gg.s[ids]),(1-vr[ids]));
    
    gg.x[ids,k]=gg.x[ids,k]+blend*np.array(vdx[ids]).T;
    gg.y[ids,k]=gg.y[ids,k]+blend*np.array(vdy[ids]).T;
    
gg.Fx=gg.x[np.ravel(gg.Fmask,order='F'),:];gg.Fy=gg.y[np.ravel(gg.Fmask,order='F'),:];
    
[gg.rx,gg.sx,gg.ry,gg.sy,gg.J] = bb.GeometricFactors2D(gg.x,gg.y,gg.Dr,gg.Ds);
[gg.nx,gg.ny,gg.sJ]= bb.Normals2D();
gg.Fscale=np.divide(gg.sJ,gg.J[np.ravel(gg.Fmask,order='F'),:]);

