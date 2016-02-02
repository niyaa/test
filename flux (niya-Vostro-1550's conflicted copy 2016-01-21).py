# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 10:45:31 2016

@author: niya
"""
from numpy import linalg as la
import aa;
import bb;
import Glob2 as gg;
import numpy as np;
from math import atan2;
def MakeCylinder2D(faces,ra,xo,yo):
    NCurveFaces = np.size(faces,0);
    vflag = np.zeros((np.size(VX)),dtype=np.float);
    for n in range (0,NCurveFaces):
        k=faces(n,0);f=faces(n,1);
        v1=gg.EToV(k,f);v2=gg.EToV(k,np.mod(f,gg.Nfaces));
        theta1 = atan2(gg.VY(v1)-yo,VX(v1)-xo);
        theta1 = atan2(gg.VY(v2)-yo,VX(v2)-xo);
        
        newx1=xo+ra*cos(theta1);newy1=yo+ra*sin(theta1);
        newx2=xo+ra*cos(theta2);newy2=yo+ra*sin(theta2);
        
        gg.VX(v1)=newx1;gg.VX(v2)=newx2;gg.Vy(v1)=newy1;gg.VY(v2)=newy2;
        
        vflag(v1)=1; vflag(v2)=1;
    vflag =  vflag(gg.EToV);
    ks=np.whereis(np.sum(vflag,2)>0);
    
    va=gg.EToV(ks,0).T; vb=gg.EToV(ks,1).T;vc=gg.EToV(ks,2).T;
    gg.x[:,ks]=0.5*(-np.matrix(gg.r+gg.s).T*np.matrix(gg.VX[va])+np.matrix(1+gg.r).T*np.matrix(gg.VX[vb])+np.matrix(1+gg.s).T*np.matrix(gg.VX[vc]));    
    gg.y[:,ks]=0.5*(-np.matrix(gg.r+gg.s).T*np.matrix(gg.VY[va])+np.matrix(1+gg.r).T*np.matrix(gg.VY[vb])+np.matrix(1+gg.s).T*np.matrix(gg.VY[vc]));                        
    
    for n in range (0,NCurveFaces):
        k=faces(n,0);f=faces(n,1);
        if(f==1): v1=gg.EToV[k,0];v2=gg.EToV[k,1];vr=gg.r;
        if(f==2): v1=gg.EToV[k,1];v2=gg.EToV[k,2];vr=gg.s;
        if(f==3): v1=gg.EToV[k,0];v2=gg.EToV[k,2];vr=gg.s;
        fr=vr[gg.Fmask[:,f]];
        x1=gg.VX(v1); y1=gg.VY(v1);
        x2=gg.VX(v2); y2=gg.VY(v2);
        
        theta1=atan2(y1-y0,x1-xo);theta2=atan2(y2-y0,x2-x0);
        if((theta2>0)&(theta1<0)):theta1 = theta1 + 2*pi;
        if((theta2>0)&(theta1<0)):theta2 = theta2 + 2*pi;
        
        theta= 0.5*theta1*(1-fr) + 0.5*theta2*(1+fr);
        
        fdx=  xo + ra*np.cos(theta)-x[gg.Fmask[:,f],k];
        fdy=  yo + ra*np.sin(theta)-y[gg.Fmask[:,f],k];
        
        
        Vface = aa.Vandermonde1D(gg.N,fr);
        Vvol  = aa.Vandermonde1D(gg.N,vr);
        
        vdx=Vvol*(Vface*la.inv(fdx));vdy= Vvol*(Vface*la.inc(fdy));
        
        ids= np.where(abs(1-vr)>1e-7);
        
        if(f==1): blend=-np.divide((gg.r[ids]+gg.s[ids]),(1-vr[ids]));
        if(f==2): blend= np.divide((gg.r[ids]+1),(1-vr[ids]));
        if(f==3): blend= -np.divide((gg.r[ids]+gg.s[ids]),(1-vr[ids]));
        
        gg.x[ids,k]=gg.x[ids,k]+blend*vdx[ids];
        gg.y[ids,k]=gg.y[ids,k]+blend*vdy[ids];
        
    Fx=gg.x[Fmask[:],:];Fy=gg.y[Fmask[:],:];
        
    [rx,sx,ry,sy,J] = bb.GeometricFactors2D(gg.x,gg.y,gg.Dr,gg.Ds);
    [nx,ny,sJ]= bb.Normals2D();
    Fscale=np.divide(sj,J[Fmask,:]);
    
def EulerLimiter2D(Q,SolutionBC,time):
    gamma=1.4;
    AVE=np.sum(gg.MassMatrix)/2;
    dropAVE = np.identity(gg.Np)-np.ones((gg.Np,1),dtype=np.int);
    dx=dropAVE*x; dy=dropAVE*y;
    
    E1=gg.EToE[:,0].T; E2= gg.EToE[:,1].T;E3=gg.EToE[:,2].T;
    
    v1=gg.EToV[:,0]; xy1= gg.VX[v1]; yv1=gg.VY[v1];
    v2=gg.EToV[:,1]; xy2= gg.VX[v2]; yv2=gg.VY[v2];
    v3=gg.EToV[:,2]; xy3= gg.VX[v3]; yv3=gg.VY[v3];
    
    fnx=np.vstack(((yv2-yv1),(yv3-yv2),(yv1-yv3)));
    
    fny=np.vstack(((xv2-xv1),(xv3-xv2),(xv1-xv3)));
    
    fl=np.sqrt(fnx**2+fny**2); fnx=np.divide(fnx,fL);
    fny=np.divide(fny,fL);
    
    xc0 = AVE*gg.x; xc1=xc0[E1]; xc2=xc0[E2]; xc3=xc0[E3];
    yc0 = AVE*gg.y; yc1=yc0[E1]; yc2=yc0[E2]; yc3=yc0[E3];
    
    A0=AVE*gg.J*2/3; A1= A0+A0[E1]; A2=A0+A0[E2]; A3 = A0+A0[E3];
    
    id1= np.where(gg.BCType[:,0]); id2= np.where(gg.BCType[:,1]);
    id3= np.where(gg.BCType[:,2]);
    
    H1= np.divide(2*(A0[id1]),fL[0,id1]);
    xc2[id1]= xc1[id1]+2*fnx[0,id1]*H1;
    yc1[id1]= yc1[id1]+2*fny[0,id1]*H1;
    
    
    H2= np.divide(2*(A0[id1]),fL[1,id1]);
    xc2[id2]= xc2[id2]+2*fnx[1,id2]*H2;
    yc1[id2]= yc1[id2]+2*fny[1,id2]*H2;
    
    
    H1= np.divide(2*(A0[id3]),fL[2,id3]);
    xc2[id3]= xc1[id3]+2*fnx[2,id3]*H3;
    yc1[id3]= yc1[id3]+2*fny[2,id3]*H3;
    
    rho=Q[:,:,0]; rhou = Q[:,:,1]; 
    rhov = Q[:,:,2]; Ener = Q[:,:,3];
    
    rhoC = AVE*rho; rhouC = AVE*rhou; rhovC = AVE*rhov;
    EnerC = AVE*Ener;
    
    x=y;
    averhou = np.ones((gg.np,1),dtype=np.int)*rhouC;
    averhov = np.ones((gg.np,1),dtype=np.int)*rhovC;
    averho  = np.ones((gg.Np,1),dtype=np.int);
    averEner = np.ones((gg.Np,1),dtype= np.int)*EnerC;

    PC0[0,:,0]=rhoC; PC0[0,:,1]=np.divide(rhouC,rhoC);
    PC0[0,:,2]=np.divide(rhovC,rhoC);

    PC[:,:,3]=(gamma-1)*(EnerC - 0.5*np.divide((rhouC**2+rhovC**2),rhoC));    
    
    PC[:,:,0]=rhoC(np.ravel(gg.EToE,order='F')[:]);
    PC[:,:,1]=rhoC(np.ravel(gg.EToE,order='F')[:]);
    PC[:,:,2]=rhoC(np.ravel(gg.EToE,order='F')[:]);
    PC[:,:,3]=rhoC(np.ravel(gg.EToE,order='F')[:]);
    
    idW = np.where(BCType.T==Wall);idI = np.where(BCType.T==In);
    idO = np.where(BCType.T==Out); idC = np.where(BCType.T==Cyl);       
            
    
    ################
    
    ids = np.vstack((1,gg.Nfp,gg.Nfp+1,2*gg.Nfp,3*gg.Nfp,2*gg.Nfp+1))
    vmapP1 = np.reshape(gg.vmapP,(gg.Nfp*gg.Nfaces,gg.K),order='C');
    vmapP1 = gg.vmapP1[ids,:];
    vmapM1 = np.reshape(gg.vmapM,(gg.Nfp*gg.Nfaces,gg.K),order='C');
    vmapM1 = gg.vmapM1[ids,:];
    
    rhoA = (rho[vmapP1]+ rho[vmapM1])/2;
    EnerA = (Ener[vmapP1]+Ener[vmapM1])/2;
    rhouA = (rhou[vmapP1]+rhou[vmapM1])/2;
    rhovA = (rhov[vmapP1]+rhov[vmapM1])/2;

    uA = np.divide(rhouA,rhoA);
    vA = np.divide(rhovA,rhoA);
    pA = (gamma -1)*(EnerA - 0.5*rhoA*(uA**2+vA**2));
    PVA[:,:,0]= rhoA; PVA[:,:,1]=uA;
    PVA[:,:,2]= vA; PVA[:,:,3]= pA;
    
    aV = np.zeros((gg.Np,gg.K,4),dtype=np.float);
    dV = np.zeros((gg.Np,gg.K,4),dtype=np.float);
    
    for n in range (0,4):
        VC0 = PC0[0,:,n]; VC1=PC[0,:,n];
        VC2 = PC[1,:,n]; VC3 = PC[2,:,n];
        
       
        VA = PVA[:,:,n];
    
        dVdxE1 =  np.divide(0.5*( (VC1-VC0)*(yv2-yv1) + (VA[1,:]-VA[2,:])*(yc1 - yc0) ),A1);
        dVdyE1 = np.divide(-0.5*( (VC1-VC0)*(xv2-xv1) + (VA[1,:]-VA[2,:])*(xc1 - xc0) ),A1);
        dVdxE2 = np.divide( 0.5*( (VC2-VC0)*(yv3-yv2) + (VA[3,:]-VA[4,:])*(yc2 - yc0) ),A2); 
        dVdyE2 = np.divide(-0.5*( (VC2-VC0)*(xv3-xv2) + (VA[3,:]-VA[4,:])*(xc2 - xc0) ),A2);
        dVdxE3 = np.divide( 0.5*( (VC3-VC0)*(yv1-yv3) + (VA[5,:]-VA[6,:])*(yc3 - yc0) ),A3);
        dVdyE3 = np.divide(-0.5*( (VC3-VC0)*(xv1-xv3) + (VA[5,:]-VA[6,:])*(xc3 - xc0) ),A3);        
        
        dVdxC0 =  np.divide( (A1*dVdxE1 + A2*dVdxE2 + A3*dVdxE3),(A1+A2+A3));
        dVdyC0 = np.divide((A1*dVdyE1 + A2*dVdyE2 + A3*dVdyE3),(A1+A2+A3));
  
        dVdxC1 = dVdxC0[E1]; dVdxC2 = dVdxC0[E2]; dVdxC3 = dVdxC0[E3];
        dVdyC1 = dVdyC0[E1]; dVdyC2 = dVdyC0[E2]; dVdyC3 = dVdyC0[E3];
        
        g1 = (dVdxC1**2+dVdydC1**2); g2= (dVdxC2**2+dVdyC2**2);
        g3 = (dVdxC3**2+dVdyC3**2);
        
        epse = 1e-10; fac = g1**2 + g2**2 + g3**2;
        w1 = np.divide((g2**g3+epse),(fac+3*epse));
        w2 = np.divide((g1**g3+epse),(fac+3*epse));
        w3 = np.divide((g1**g2+epse),(fac+3*epse));

        LdVdxC0 = w1*dVdxC1 + w2*dVdxC2 + w3*dVdxC3;
        LdVdyC0 = w1*dVdyC1 + w2*dVdyC2 + w3*dVdyC3; 
        
        dV[:,:,n] = dx*(np.ones((gg.Np,1),dtype=np.int)*LdVdxC0) + dy*(np.ones((gg.Np,1),dtype=int)*LdVdyC0);
        aV[:,:,n] = ones((gg.Np,1),dtype=np.int)*VC0;
        
    avrerho = av[:,:,0]; aveu= aV[:,:,1]; avev = av[:,:,2];
    avep = av[:,:,3];
    
    drho = dV[:,:,0]; du = dV[:,:,1]; dv= dV[:,:,2];
    dp = dV[:,:,3];
    
    tol = 1e-2;
    Lrho = averho + drho; ids = np.where(np.min)
    
    Lrho = averho + drho; ids = np.where(np.min(Lrho,1))
    
    
    
    
    
        