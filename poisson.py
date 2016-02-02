# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 02:00:15 2015

@author: niya
"""
import Glob2 as gg,aa,bb;
import numpy as np;
from numpy import linalg as la;
def PoissonIPDG2D():
    massEdge=np.zeros((gg.Np,gg.Np,gg.Nfaces),dtype=np.float);
    Fm=gg.Fmask[:,0];faceR=gg.r[Fm];
    V1D=aa.Vandermonde1D(gg.N,faceR);
    temp=la.inv((V1D*V1D.T));
    for i in range (0,len(Fm)):
        for j in range (0,len(Fm)):
            massEdge[Fm[i],Fm[j],0]=temp[i,j]; 
    
    Fm=gg.Fmask[:,1];faceR=gg.r[Fm];
    V1D=aa.Vandermonde1D(gg.N,faceR);
    temp=la.inv((V1D*V1D.T));
    for i in range (0,len(Fm)):
        for j in range (0,len(Fm)):
            massEdge[Fm[i],Fm[j],1]=temp[i,j]; 
    
    Fm=gg.Fmask[:,2];faceS=gg.s[Fm];
    V1D=aa.Vandermonde1D(gg.N,faceR);
    temp=la.inv((V1D*V1D.T));
    for i in range (0,len(Fm)):
        for j in range (0,len(Fm)):
            massEdge[Fm[i],Fm[j],2]=temp[i,j]; 
            
            
    
    MM=np.zeros((gg.K*gg.Np*gg.Np,3));OP=np.zeros((gg.K*gg.Np*gg.Np*(1+gg.Nfaces),3));
    entries=np.arange(0,gg.Np*gg.Np);entriesMM=np.arange(0,gg.Np*gg.Np);
    
    for k1 in range(0,5):
        if not (np.mod(k1,1000)): k1=k1;
        rows1=np.arange((k1*gg.Np),(k1+1)*gg.Np).reshape(gg.Np,1)*np.ones((1,gg.Np),dtype=np.int);
        
        cols1=rows1.T;
        
        Dx=gg.rx[0,k1]*gg.Dr + gg.sx[0,k1]*gg.Ds;
        Dy=gg.ry[0,k1]*gg.Dr + gg.sy[0,k1]*gg.Ds;
        
        OP11 = gg.J[0,k1]*(Dx.T*gg.MassMatrix*Dx+Dy.T*gg.MassMatrix*Dy);
        
        for f1 in range (0,gg.Nfaces):
            k2=gg.EToE[k1,f1];f2=gg.EToF[k1,f1];
            rows2=np.arange((k2*gg.Np),(k2+1)*gg.Np).reshape(gg.Np,1)*np.ones((1,gg.Np),dtype=np.int);
            cols2=rows2.T;            
            fidM = (k1)*gg.Nfp*gg.Nfaces + (f1)*gg.Nfp + np.arange(0,gg.Nfp);
            gg.vidM = gg.vmapM[fidM];Fm1 = np.mod((gg.vidM),gg.Np);
            gg.vidP = gg.vmapP[fidM];Fm2 = np.mod((gg.vidP),gg.Np);
            
            id = 1+(f1)*gg.Nfp + (k1)*gg.Nfp*gg.Nfaces;
            lnx = np.ravel(gg.nx,order='F')[id]; lny =np.ravel(gg.ny,order='F')[id]; lsJ = np.ravel(gg.sJ,order='F')[id];
            hinv = max(np.ravel(gg.Fscale,order='F')[id],gg.Fscale[1+(f2)*gg.Nfp,k2]);
            
            Dx2 = gg.rx[0,k2]*gg.Dr +gg.sx[0,k2]*gg.Ds;
            Dy2 = gg.ry[0,k2]*gg.Dr + gg.sy[0,k2]*gg.Ds;
            
            
            Dn1 = lnx*Dx + lny*Dy;
            Dn2 = lnx*Dx2 + lny*Dy2;
            
            mmE = lsJ*massEdge[:,:,f1];
            
            gtau = 100*2*(gg.N+1)*(gg.N+1)*hinv;
            
            if (gg.BCType[k1,f1]==gg.Dirichlet):
                OP11 = OP11 + (gtau*mmE - mmE*Dn1 -Dn1.T*mmE);
            #elif(gg.BCType(k,f1)==Neuman):
                
            elif():
                OP11=OP11 + 0.5*(gtau*mmE - mmE*Dn1 - Dn1.T*mmE);
                OP12=np.zeros((gg.Np,gg.Np),dtype=np.float);
                OP12[:,Fm2]= -0.5*(gtau*mmE[:,Fm1]);
                for i in range (0,len(Fm1)):
                    OP12[Fm1[i],:]= OP12[Fm1[i],:]-0.5*(mmE[Fm1[i],Fm1])*Dn2[Fm2,:];
                OP12[:,Fm2]= OP12[:,Fm2]-0.5*(-Dn1.T*mmE[:,Fm1]);
              
                OP[entries,0]=np.ravel(rows1,order='F')[:];OP[entries,1]=np.ravel(cols2,order='F')[:];OP[entries,2]=np.ravel(OP12,order='F')[:];
                entries = entries + gg.Np*gg.Np;
                
            OP[entries,0]=np.ravel(rows1,order='F')[:];OP[entries,1]=np.ravel(cols1,order='F')[:];OP[entries,2]=np.ravel(OP11,order='F')[:];
            MM[entriesMM,0]=np.ravel(rows1,order='F')[:];MM[entriesMM,1]=np.ravel(cols1,order='F')[:];MM[entriesMM,2]=gg.J[0,k1]*np.ravel(gg.MassMatrix,order='F')[:];
            entries = entries+gg.Np*gg.Np; entriesMM = entriesMM + gg.Np*gg.Np;
            
    OP = OP[range(0,(np.max(entries)-gg.Np*gg.Np)),:];OP[np.abs(OP)<gg.eps]=0;
    MM = MM[range(0,(np.max(entries)-gg.Np*gg.Np)),:];MM[np.abs(OP)<gg.eps]=0;
    
    return OP, MM;
        
    

def PoissonIPDGbc2D(ubc,qbc):
    massEdge=np.zeros((gg.Np,gg.Np,gg.Nfaces),dtype=np.float);
    Fm=gg.Fmask[:,0];faceR=gg.r[Fm];
    V1D=aa.Vandermonde1D(gg.N,faceR);
    temp=la.inv((V1D*V1D.T));
    for i in range (0,len(Fm)):
        for j in range (0,len(Fm)):
            massEdge[Fm[i],Fm[j],0]=temp[i,j]; 
    
    Fm=gg.Fmask[:,1];faceR=gg.r[Fm];
    V1D=aa.Vandermonde1D(gg.N,faceR);
    temp=la.inv((V1D*V1D.T));
    for i in range (0,len(Fm)):
        for j in range (0,len(Fm)):
            massEdge[Fm[i],Fm[j],1]=temp[i,j]; 
    
    Fm=gg.Fmask[:,2];faceS=gg.s[Fm];
    V1D=aa.Vandermonde1D(gg.N,faceR);
    temp=la.inv((V1D*V1D.T));
    for i in range (0,len(Fm)):
        for j in range (0,len(Fm)):
            massEdge[Fm[i],Fm[j],2]=temp[i,j]; 
            
            
    bc=np.zeros((gg.Np,gg.K),dtype=int);
    for k1 in range (1,np.K+1):
        if not (np.mod(k1,1000)): 
            for f1 in range (0,gg.Nfaces):
               if(gg.BCType[k1,f1]):
                    Fm1= Fmask[:,f1];
                    fidM = (k1-1)*gg.Nfp*gg.Nfaces + (f1-1)*gg.Nfp + np.arange(1,gg.Nfp+1);
                    
                    id = 1+(f1-1)*gg.Nfp + (k1-1)*gg.Nfp*gg.Nfaces;
                    lnx = nx[id]; lny = ny[id]; lsJ = sJ[id]; hinv=gg.Fscale[id];
                    
                    Dx = gg.rx[0,k1]*gg.Dr + gg.sx[0,k1]*Ds;
                    Dy=  gg.ry[0,k1]*gg.Dr + gg.sy[0,k1]*Ds;
                    Dn1=lnx*gg.Dx + lny*gg.Dy;
                    
                    mmE = lsJ*massEdge[:,:,f1];
                    gtau = 100*2*(gg.N+1)*hinv;
                    
                    if (BCType[k1,f1]==Dirichlet):
                        bc[:,k1]=bc[:,k1]+(gtau*mmE[:,Fm1]-Dn1.T*mmE[:,Fm1])*ubc[fidM];
                    elif(BCType[k1,f1]==Neuman):
                        bc[:,k1]=bc[:,k1]+mmE[:,Fm1]*qbc[fidM];
    return bc;
                        
            
            
            
def PoissonDriver2D():
    gg.N =3;
    
    Startup2D();
    BuildBCMaps2D();
    [A, M]=PoissonIPDG2D();
    
    # Set up Dirichlet boundary conditions
    uD = np.zeros(gg.Nfp*gg.Nfaces,gg.K);
    uD[gg.mapD] = np.sin(np.pi*gg.Fx[mapD])*np.sin(np.pi*gg.Fy[mapD]);
    
    #Set up Neumann boundary coniditions
    qN = np.zeros((gg.Nfp*gg.Nfaces,K),dtype=np.float);
    qN[mapN] = nx[gg.mapN]*(np.pi*np.cos(np.pi*gg.Fx[gg.mapN])*sin(np.pi*gg.Fy[gg.mapN]))+ \
    ny[gg.mapN]*(np.pi*sin(np.pi*gg.Fx[gg.mapN]))*np.cos(np.pi*gg.Fy[gg.mapN]);
    
    Aqbc = PoissonIPDGbc2D (uD,qN);
    
    rhs = -2*(np.pi**2)*np.sin(np.pi*x)*np.sin(np.pi*y);
    
    rhs = -gg.MassMatrix*(np.multiply(gg.J,rhs)) + Aqbc;
    
    u= la.solve(A,np.ravel(rsh,order='F'));
    
    u= np.reshape(u,(Np,K),order='F');
    
    
    
    
    
    
def PoissonRHS2D(u):
    du=np.zeros(gg.nfp);

  
def Poisson2D():
    p1=gg.K*gg.Np;
    g=np.zeros((p1,1),dtype=np.float);
    A= csr_matrix((p1,p1));
    M= csr_matrix((p1,p1));
    
    for i in range (0,p1):
        g[i,0]=1.0;
        gmat=np.reshape(g,(Np,K),order='F');
        [Avec, Mvec]=PoissonRHS2D(gmat);
        ids=np.where(Avec);
        A[ids,i]=Avec[ids];
        ids=np.where(Mvec);
        M[ids,i]=Mvec(ids);
        g[i]=0.0;
        
    return A, M;




#def plot1qvtu(title, Ux):

def PoissonRHS2D(U):
    import Glob2 as gg;
    du=np.zeros(gg.Nfp*gg.Nfaces,gg.K);du()
    