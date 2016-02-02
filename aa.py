# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 17:51:18 2015

@author: niyaa
"""

import numpy as np;
from math import gamma
from scipy.special import gamma;
from numpy import linalg as la
eps=pow(10,-10);
import Glob1 as gg;

def JacobiP(x,alpha,beta,N):
    xp=np.array(x);
    dims=np.shape(xp);
    
    if (len(dims)>1):
        if (dims[1]>1): xp=xp.T;
    xp=np.reshape(xp,len(xp));
        
    from math import gamma
    #if dims[0] >1:
     #   xp=xp.T;
    PL=np.zeros((N+1,np.size(x)),dtype=np.float);
    gamma0=0.0;gamma1=0.0;
    gamma0 = 2**(alpha+beta+1.0)/(alpha+beta+1.0)*gamma(alpha+1.0)*gamma(beta+1.0)/gamma(alpha+beta+1.0);
     
    PL[0,:]= float(1.0/np.sqrt(gamma0));
    
   
    if N==0:
        P=PL.T
        P=np.reshape((P),(len(P),));
        
        return P;
    

    gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3.0)*gamma0;
       
    PL[1,:]= ((alpha+beta+2)*xp/2 + (alpha-beta)/2.0)/np.sqrt(gamma1);
    #print PL;
    if (N==1):
            P=np.transpose(PL[N,:]);
            
            return P;
 
    aold = 2.0/(2.0+alpha+beta)*np.sqrt((alpha+1)*(beta+1)/(alpha+beta+3.0));
    h1=0.0;anew=0.0;bnew=0.0;
    for i in range (1,N):
           h1 = 2*i+alpha+beta;
           anew = 2.0/(h1+2.0)*np.sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1.0)/(h1+3.0));
           bnew=-(alpha**2-beta**2)/float(h1)/(h1+2.0);       
           PL[i+1,:] = 1.0/(anew)*( -aold*PL[i-1,:] + (xp-bnew)*PL[i,:]);
           aold = anew;
           #print "h1,anew,bnew,aold,gamma0,gamma1"
           #print h1,anew,bnew,aold,gamma0,gamma1
           
    #print PL;       
        
    P=np.transpose(PL[N,:]);
    return P

           
        
       
def JacobiGQ(alpha,beta,N):
    x=np.zeros(N+1,dtype=np.float);w=np.zeros(N+1,dtype=np.float);
    if (N==0): x[0]=-(alpha-beta)/float(alpha+beta+2);w[0]=2;return (x,w);

    J=np.zeros(N+1,dtype=np.float);
    h1=2*np.arange(0,N+1)+alpha+beta;
    temp=-1/2.0*(pow(alpha,2)-pow(beta,2))/(h1+2.0)/(h1);
    temp1=np.diagflat(temp);
    temp2=2.0/(h1[0:N]+2.0)*np.sqrt(np.arange(1,N+1)*(np.arange(1,N+1)+alpha+beta)*(np.arange(1,N+1)+alpha)*(np.arange(1,N+1)+beta)/(h1[0:N]+1.0)/(h1[0:N]+3.0));
    temp3=np.diagflat(temp2,1);
    J=temp1+temp3;
    if (alpha+beta<eps):    J[0,0]=0.0;
    J= J + J.T;  
     
    D, V=la.eig(J);
    x1=D;
 
    w1=np.matrix((V[0,:]*V[0,:])).T*pow(2,(alpha+beta+1))/(alpha+beta+1.0)*gamma(alpha+1.0)*gamma(beta+1)/gamma(alpha+beta+1.0);
    if (N==1):  x=np.sort(x1);w=np.sort(w1);    return x, w    
    if (N%2 ==0):
        n=N/2+1;
    else:
        n=(N+1)/2;
    
    w[n:]=sorted(w1[n:],reverse=True);
    x[n:]=sorted(x1[n:]);
    for i in range(0,n):
        w[i]=w1[i];
        x[i]=x1[i];
       
    
    
    return x, w
  

def JacobiGL(alpha,beta,N):
    x=np.zeros((N+1,1),dtype=np.float);
    if (N==1):  x[0]=-1.0;x[1]=1.0; return x;
    [xint,w]=JacobiGQ(alpha+1,beta+1,N-2);
    x=np.concatenate(([-1], xint, [1]));
    return x

  
def Vandermonde1D(N,r):
    V1D=np.zeros((len(r),N+1),dtype=np.float);
    for i in range(0,N+1):
        V1D[:,i]=JacobiP(r,0,0,i);
    V1D=np.matrix(V1D);     
    return V1D    
    
    
def GradJacobiP(r,alpha,beta,N):
    dP=np.zeros((len(r),),dtype=np.float);
    if N==0:    dP[:]=0.0;    return dP;
    dP=np.sqrt(N*(N+alpha+beta+1))*JacobiP(r,alpha+1,beta+1,N-1);
    return dP;


def GradVandermonde1D(N,r):
    DVr=np.zeros((len(r),N+1),dtype=np.float);
    for i in range(0,N+1):
        y0=GradJacobiP(r,0,0,i);
        
        DVr[:,i]=y0;
    DVr=np.matrix(DVr);
    return DVr;
    
def Dmatrix1D(N,r,V):
    Vr=GradVandermonde1D(N,r);
    Dr=np.zeros((len(r),N+1),dtype=np.float);
    V=np.matrix(V);
    Vr=np.matrix(Vr);
    c=la.inv(V);
    #if (len(r)==N+1):
     #   Dr=la.solve(Vr,V);
      #  return Dr;
        
    Dr=Vr*c;
    return Dr; 
    

    
    
def Connect1D(EToV):
    gg.Nfaces = 2;
    from scipy.sparse import lil_matrix, csr_matrix;
    Nfaces = 2;
    
    
    K=np.size(EToV,0);
    Nv = K+1;
    TotalFaces = Nfaces*K;
   
    vn=np.array([0,1],dtype=np.int);
    
    SpFToV = csr_matrix((TotalFaces,Nv),dtype=np.int);
    SpFToV=lil_matrix(SpFToV);
   
    sk=0;
    for k in range (0,K):
        for face in range (0,Nfaces):
            t1=EToV[k,vn[face,:]];
        
            SpFToV[sk,t1] = 1;
            sk = sk+1;
    temp=np.identity(TotalFaces);
    temp1=csr_matrix(temp);
    SpFToV=SpFToV*SpFToV.T - temp1;
    [faces2, faces1]=np.where(SpFToV.toarray()==1);
    faces1=faces1+1;faces2=faces2+1;
    element1 = np.floor((faces1-1)/Nfaces)+1;
    element2 = np.floor((faces2-1)/Nfaces)+1;
    face1=np.remainder((faces1-1),Nfaces)+1;
    face2=np.remainder((faces2-1),Nfaces)+1;
    ind=np.zeros((len(face1),),dtype=np.int);
    for i in range (0,len(face1)):
        ind[i] = K*(face1[i]-1)+element1[i];
        
    EToE = np.arange(1,K+1).reshape(K,1)*np.ones((1,Nfaces),dtype=int);
    EToF = np.ones((K,1),dtype=int)*np.arange(1,Nfaces+1).reshape(1,Nfaces);    
    #[b,c]=np.unravel_index([ind-1],(16,3)) 
    EToE=EToE.flatten(order='F');
    EToE[ind-1]=element2;
    EToE=np.reshape(EToE,(K,Nfaces),order='F');
    EToF=EToF.flatten(order='F');
    EToF[ind-1]=face2;
    EToF=np.reshape(EToF,(K,Nfaces),order='F');
    return EToE, EToF;
    
    
def Filter1D(N,Nc,s):
    import Glob1 as gg;
    filterdiag = np.ones((N+1,1));
    alpha = - np.log(eps);
   
    for i in range (Nc-1,N):
        filterdiag[i]=np.exp(-alpha*((i-Nc)/(N - Nc))**s);
   
    F=gg.V*np.diag(filterdiag)*gg.invV;
    return F;      
      
      
      
    
def BuildMaps2D():
    import Glob1 as gg;
    nodeids=np.arange(0,gg.K*gg.Np).reshape(gg.Np,gg.K,order='F');
    gg.vmapM=np.zeros((gg.Nfp,gg.Nfaces,gg.K),dtype=np.int);    
    gg.vmapP=np.zeros((gg.Nfp,gg.Nfaces,gg.K),dtype=np.int);

  
    for k1 in range (0,gg.K):
        for f1 in range (0,gg.Nfaces):
            gg.vmapM[:,f1,k1]=nodeids[gg.Fmask[:,f1],k1];
               
   
    for k1 in range (0,gg.K):
        for f1 in range (0,gg.Nfaces):
            k2=gg.EToE[k1,f1];f2=gg.EToF[k1,f1];
            gg.vidM=gg.vmapM[:,f1,k1];gg.vidP=gg.vmapM[:,f2,k2];
            x1=np.ravel(gg.x,order='F')[gg.vidM];
            x2=np.ravel(gg.x,order='F')[gg.vidP];
            D=(x1-x2.T)**2;
            [idM,idP]=np.where(np.sqrt(np.abs(D))< gg.NODETOL*refd);
            #print idM, idP;
           #if(D<gg.eps):
                
    gg.vmapP=np.ravel(gg.vmapP,order='F');
    gg.vmapM=np.ravel(gg.vmapM,order='F');
    gg.mapB=np.where(gg.vmapP==gg.vmapM);
    gg.vmapB=gg.vmapM[gg.mapB];
    mapI=0;mapO=gg.K*gg.Nfaces-1;vmapI=0;vmapI=gg.K*gg.Np-1;
    return gg.vmapM,gg.vmapP,gg.vmapB,gg.mapB;
      

def MeshGen1D(xmin,xmax,K):
    Nv=K+1;
    VX=np.arange(0,Nv,dtype=np.float);
    for i in range(0,Nv):
        VX[i]=(xmax - xmin)*(i)/(Nv-1.0)+xmin;
        
    EToV=np.zeros((K,2),dtype=np.int);
    for k in range (0,K):
        EToV[k,0]=k;EToV[k,1]=k+1;
        
        
    return Nv, VX, K, EToV;
       
def Startup1D():
    gg.NODETOL=1e-10;
    gg.Nfp=1; gg.Np=gg.N+1;gg.Nfaces=2;
    gg.r = JacobiGL(0,0,gg.N);
    
    
    gg.V=Vandermonde1D(gg.N,gg.r);
    
    gg.invV = la.inv(gg.V);
    
    [gg.Dr]=Dmatrix1D(gg.N,gg.r,gg.V);    

    gg.LIFT=LIFT1D();
    
    va=gg.EToV[:,0].T;vb=gg.EToV[:,1].T;
    gg.x=np.ones((gg.N+1,1),dtype=np.float)*np.matrix(gg.VX[va]) + 0.5*np.matrix(gg.r+1).T*np.matrix(gg.VX[vb]-gg.VX[va]);
    fmask1 = np.where(np.abs(gg.r+1)<gg.NODETOL);
    fmask2 = np.where(np.abs(gg.r-1)<gg.NODETOL);
    
    gg.Fmask = np.concatenate([fmask1,fmask2]);
    gg.Fmask=gg.Fmask.T;
    
    Fx=x[np.ravel(gg.Fmask,order='F')[:],:];
    
    [gg.rx,gg.J]=GeometricFactors1D(gg.x,gg.Dr);
    gg.nx=Normals1D();
    
    gg.Fscale=np.divide(1,gg.J[gg.Fmask.flatten(order='F'),:])
    
    [gg.EToE, gg.EToF]=Connect1D(gg.EToV);
    gg.EToE-=1;gg.EToF-=1;
    
    BuildMaps2D();
   
    
def LIFT1D():
    import Glob1 as gg;
    Emat=np.zeros((gg.Np,gg.Nfaces),dtype=np.int);
    Emat[0,0]=1.0; Emat[gg.Np-1,1]=1.0;
    LIFT=gg.V*(gg.V.T*Emat);
    return LIFT;
    
def GeometricFactors1D(x,Dr):
    xr=gg.Dr*gg.x;gg.J=xr;gg.rx=np.divide(1,gg.J);
    return gg.rx, gg.J;

def Normals1D():
    gg.nx=np.zeros((gg.Nfp*gg.Nfaces,gg.K),dtype=np.int);
    gg.nx[0,:]=-1.0;gg.nx[1,:]=1.0;
    return gg.nx;
    
