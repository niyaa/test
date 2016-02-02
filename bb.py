# -*- oding: utf-8 -*-
"""
Created on Wed Dec  2 16:01:21 2015

@author: niyaa
"""

#global Np,Nfp,N,K,r,s,Dr,Ds,LIFT,Drw,Dsw;
#global MassMatrix,Fx,Fy,nx,ny,jac,Fscale,J;
#global vmapM, vmapP, vmapB, mapB,Fmask;
#global BCType, mapI, mapO, mapW, mapF, mapC, mapS, mapM, mapP;
#global mapD, mapN, vmapI, vmapO, vmapC,vmapS, vmapD, vmapN;
#global rx,ry,sx,sy,J,sJ;

import aa;
import Glob2 as gg;
import numpy as np;
import bb;
#import gmsh;
    
import sys
#obj=gmsh.Mesh();
#obj.read_msh('sample.msh');
   
gg.NODETOL = 1e-15;
from numpy import linalg as la



from math import cos, sin, pi
from numpy import linalg as la

gg.eps=1e-15;

def GradSimplex2DP(a,b,i,j):
    fa=aa.JacobiP(a,0,0,i); dfa=aa.GradJacobiP(a,0,0,i);
    gb=aa.JacobiP(b,2*i+1,0,j);dgb=aa.GradJacobiP(b,2*i+1,0,j);
    dmodedr = dfa*gb;
    if (i>0):
        dmodedr = dmodedr*((0.5*(1-b))**(i-1));
    dmodeds = dfa*(gb*(0.5*(1+a)));
    if (i>0):
        dmodeds = dmodeds*((0.5*(1-b))**(i-1));
    tmp=dgb*((0.5*(1-b))**i);
    if (i>0):
        tmp= tmp - 0.5*i*gb*((0.5*(1-b))**(i-1));
    dmodeds = dmodeds + fa*tmp;
    dmodedr = pow(2,(i+0.5))*dmodedr;dmodeds=pow(2,(i+0.5))*dmodeds;
    return dmodedr, dmodeds
    
def Simplex2DP(a,b,i,j):
    h1=aa.JacobiP(a,0,0,i);h2=aa.JacobiP(b,2*i+1,0,j);
    P=np.sqrt(2.0)*h1*h2*np.power((1-b),i);
    return P

            
            
def Nodes2D(N):
    #import Glob2 as gg;
    alpopt=np.array([0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, \
    1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258]);    
    if (N<16): 
        alpha = alpopt[N];
    else:
        alpha = 5.0/3.0;
        
    Np=(N+1)*(N+2)/2;
    
    L1=np.zeros((Np,1),dtype=np.float);L2=np.zeros((Np,1),dtype=np.float);L3=np.zeros((Np,1),dtype=np.float);
    sk=0;
    for n in range (1,N+2):
        for m in range (1,N+3-n):
            L1[sk]=(n-1.0)/N;L3[sk]=(m-1.0)/N;
            sk = sk+1;
            
    L2=1.0 - L1 - L3;
    x= -L2+L3;
    y= (-L2 -L3+2*L1)/np.sqrt(3.0);
    
    blend1 = 4*L2*L3;blend2 = 4*L1*L3; blend3 = 4*L1*L2;

    warpf1 = WarpFactor(N,L3-L2);warpf2 = WarpFactor(N,L1-L3);
    warpf3 = WarpFactor(N,L2-L1);

    warp1 = blend1*np.array(warpf1)*(1+(alpha*L1)**2);
    warp2 = blend2*np.array(warpf2)*(1+(alpha*L2)**2);
    warp3 = blend3*np.array(warpf3)*(1+(alpha*L3)**2); 

    
    x=x+1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
    y=y+0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
    x=np.reshape(x,(len(x),));
    y=np.reshape(y,(len(y),));
    return x,y
    
    
def WarpFactor(N, rout):
    LGLr = aa.JacobiGL(0,0,N);req = np.linspace(-1,1,N+1);
    Veq= aa.Vandermonde1D(N,req);
    Nr=len(rout);Pmat=np.zeros((N+1,Nr),dtype=np.float);
    for i in range (1,N+2):
        Pmat[i-1,:]=aa.JacobiP(rout,0,0,i-1);
    Pmat=np.matrix(Pmat);        
    Veq=np.matrix(Veq);
    Lmat = la.solve(Veq.T,Pmat);
    Lmat=Lmat.T;
    warp = Lmat*np.matrix(LGLr - req).T;
    zerof=(abs(rout)<1.0-1.0e-10);
    zerof=np.int8(zerof);
    sf=1.0 - (zerof*rout)**2;
    warp= np.divide(warp,sf)+np.array(warp)*(zerof-1);
    #return matrix
    return warp;
    
def rstoab(r,s):
    
    Np=len(r);
    a=np.zeros((Np),dtype=np.float);
    for n in range (0,Np):
        if(s[n]!=1):
            a[n] = 2*(1.0+r[n])/(1.0-s[n])-1.0;
        else:
            a[n] =-1;
    b=s;        
    return a, b;

def Normals2D():
    import Glob2 as gg;
    xr=np.matrix(gg.Dr)*np.matrix(gg.x); xs=np.matrix(gg.Ds)*np.matrix(gg.x);
    yr=np.matrix(gg.Dr)*np.matrix(gg.y); ys=np.matrix(gg.Ds)*np.matrix(gg.y);
    J=np.multiply(xr,ys)-np.multiply(xs,yr);
    fxr = xr[gg.Fmask.flatten(order='F'),:];fxs=xs[gg.Fmask.flatten(order='F'),:];
    fyr=yr[gg.Fmask.flatten(order='F'),:];fys=ys[gg.Fmask.flatten(order='F'),:];
    nx=np.zeros((3*gg.Nfp,gg.K),dtype=np.float);    
    ny=np.zeros((3*gg.Nfp,gg.K),dtype=np.float);    
    fid1=np.arange(0,gg.Nfp);fid2=np.arange(gg.Nfp,2*gg.Nfp);fid3=np.arange(2*gg.Nfp,3*gg.Nfp);
    #first face    
    nx[fid1,:]=fyr[fid1,:];ny[fid1,:]=-fxr[fid1,:];
    #face 2
    nx[fid2,:]= fys[fid2,:]-fyr[fid2,:];ny[fid2,:]=-fxs[fid2,:]+fxr[fid2,:];
    #face 3
    nx[fid3,:]= -fys[fid3,:]; ny[fid3,:]=fxs[fid3,:];
    
    sJ = np.sqrt(nx*nx+ny*ny); nx=np.divide(nx,sJ); ny=np.divide(ny,sJ);
    return nx, ny, sJ;    
    
        
       
def GeometricFactors2D(x, y, Dr, Ds):
   
    xr=np.matrix(Dr)*np.matrix(x); xs=np.matrix(Ds)*np.matrix(x);
    yr=np.matrix(Dr)*np.matrix(y); ys=np.matrix(Ds)*np.matrix(y);
    J=-np.multiply(xs,yr)+np.multiply(xr,ys);
    rx=np.divide(ys,J);sx=-np.divide(yr,J);
    ry=-np.divide(xs,J);sy=np.divide(xr,J);
    return rx,sx,ry,sy,J;
       
 
def Lift2D():
    Emat=np.zeros((gg.Np,gg.Nfaces*gg.Nfp));
    
    #face 1    
    faceR = gg.r[gg.Fmask[:,0]];
    V1D =aa.Vandermonde1D(gg.N,faceR);
    massEdge1 = la.inv(np.matrix(V1D)*np.matrix(V1D).T);
    Emat[gg.Fmask[:,0],0:gg.Nfp]=massEdge1;
    
    
    #face 2    
    faceR = gg.r[gg.Fmask[:,1]];
    V1D =aa.Vandermonde1D(gg.N,faceR);
    massEdge2 = la.inv(V1D*V1D.T);
    i=range(0,gg.Nfp);j=gg.Fmask[:,0];
    Emat[gg.Fmask[:,1],gg.Nfp:2*gg.Nfp]=massEdge2;
    
    
    #face 3    
    faceS = gg.s[gg.Fmask[:,2]];
    V1D =aa.Vandermonde1D(gg.N,faceS);
    massEdge3 = la.inv(V1D*V1D.T);
    Emat[gg.Fmask[:,2],2*gg.Nfp:3*gg.Nfp]=massEdge3;
    
    LIFT = gg.V*gg.V.T*Emat;
    return LIFT;
 
def InterpMatrix(rout, sout):
    Vout=aa.Vandermonde1D(gg.N, rout, sout);
    IM = Vout*gg.invV;
    return IM;
    
def Dmatrices2D(N,r,s,V):
    [Vr, Vs]=GradVandermonde2D(N,r,s);
    
    Dr=np.matrix(Vr)*np.matrix(la.inv(V));Ds=np.matrix(Vs)*np.matrix(la.inv(V));
    return Dr, Ds;
def Vandermonde2D(N,r,s):
    V2D = np.zeros((len(r),(N+1)*(N+2)/2),dtype=np.float);
    [a, b] = rstoab(r,s);
    sk=0;
    for i in range (0,N+1):
        for j in range (0,N-i+1):
            V2D[:,sk]=Simplex2DP(a,b,i,j);
            sk = sk+1;
    return V2D
            
            
    
    
def GradVandermonde2D(N,r,s):	
    V2Dr = np.zeros((len(r),(N+1)*(N+2)/2),dtype=np.float);
    V2Ds = np.zeros((len(r),(N+1)*(N+2)/2),dtype=np.float);
    [a, b]=rstoab(r,s);
    sk=0;
    for i in range(0,N+1):
        for j in range(0,N+1-i):
            [t1,t2]=GradSimplex2DP(a,b,i,j);
            V2Dr[:,sk]=t1;
            V2Ds[:,sk]=t2;
            sk=sk+1;
            
    #V2Dr=np.matrix(V2Dr);V2Ds=np.matrix(V2Ds);       
    return V2Dr,V2Ds;
    
def Div2D(u,v):
    ur=gg.Dr*u; us=gg.Ds*u; vr=gg.Dr*v; vs=gg.Ds*v;
    divu=np.multiply(gg.rx,ur)+np.multiply(gg.sx,us)+ \
    np.multiply(gg.rx,vr)+np.multiply(gg.ry,vs);
    return divu;
    
def EquiNodes2D(N):
    Np=(N+1)*(N+2)/2;
    L1= np.zeros((Np,1),dtype=np.float);
    L2= np.zeros((Np,1),dtype=np.float);    
    L3= np.zeros((Np,1),dtype=np.float);    
    sk=0;
    for n in range(1,N+2):
        for m in range(1,N+2-n+1):
            L1[sk]=(n-1.0)/N;L3[sk]=(m-1.0)/N;
            sk=sk+1;
    L2 = 1.0 - L1 - L3;
    x=-L2+L3; y=(-L2-L3+2*L1)/np.sqrt(3.0);
        
    return x, y;
    
def Filter2D(Norder,Nc,sp):
    filterdiag = np.ones(((Norder+1)*(Norder+2)/2));
    alpha = - np.log(aa.eps);
    
    sk=0;
    for i in range (0,Norder+1):
        for j in range (0,Norder-i+1):
            if (i+j>=Nc):
                filterdiag[sk]=np.exp(-alpha*((i+j-Nc)/(Norder - Nc))**sp);
                sk=sk+1
    F=gg.V*np.diag(filterdiag)*gg.invV;
    return F;       
    
def PhysDmatrices2D(x1,y1,interp):
    from scipy.sparse import dia_matrix
    IDr = interp*gg.Dr; IDs = interp*gg.Ds;
    [rx1,sx1,ry1,sy1,J1]=GeometricFactors2D(x1,y1,IDr,IDs);
    
    n=np.size(interp,0);
    Dx=dia_matrix((rx1,0),shape=(n,n))*IDr+dia_matrix((sx1,0),shape=(n,n))*IDs;
    Dy=dia_matrix((ry1,0),shape=(n,n))*IDr+dia_matrix((sy1,0),shape=(n,n))*IDs;
        
    
    
    return Dx,Dy;

def Grad2D(u):
    ur=gg.Dr*u; us=gg.Ds*u;
    ux=np.multiply(gg.rx,ur)+np.multiply(gg.sx,us)+np.multiply(gg.ry,ur)+ \
    np.multiply(gg.sy,us);
    
    
def xytors(x,y):
    L1=(np.sqrt(3.0)*y+1.0)/3.0;
    L2=(-3.0*x-np.sqrt(3.0)*y+2.0)/6.0;
    L3=(3.0*x-np.sqrt(3.0)*y+2.0)/6.0;
    r=-L2+L3-L1;s=-L2-L3+L1;
    
    return r,s;

def Curl2D(ux,uy,uz):
    uxr = gg.Dr*ux;uxs = gg.Ds*ux; uyr = gg.Dr*uy; uys=gg.Ds*uy;
    vz = np.multiply(gg.rx,uyr)+np.multiply(gg.sx,uys)-np.multiply(gg.ry,uxr)-\
    np.multiply(gg.sy,uxs);
    vx=np.matrix([]);vy=np.matrix([]);
    if uz:
        uzr = gg.Dr*uz;
        uzs = gg.Ds*uz;
        vx = np.multiply(gg.ry,uzr)+ np.multiply(gg.sy,uzs);
        vy=-np.multiply(gg.rx,uzr)-np.multiply(gg.sx,uzs);
    return vx,vy,vz;

def Connect2D(EToV):
    
    from scipy.sparse import lil_matrix, csr_matrix;
    Nfaces = 3;
    
   
    K=np.size(EToV,0);
    Nv = np.max(EToV)+1;
    TotalFaces = Nfaces*K;
   
    vn=np.array([[0,1],[1,2],[0,2]],dtype=np.int);
    
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
    SpFToV=SpFToV*SpFToV.T - 2*temp1;
    [faces2, faces1]=np.where(SpFToV.toarray()==2);
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
    
   
def tiConnect2D(EToV):
    Nfaces=3;
    K=np.size(EToV,0);
    Nnodes=np.amax(EToV);
    fnodes=np.concatenate((EToV[:,0:2],EToV[:,1:],EToV[:,[2,0]]),axis=0);
    fnodes=np.sort(fnodes,axis=-1)-1;
    EToE =np.arange(1,K+1).reshape(K,1)*np.ones((1,Nfaces),dtype=np.int);
    EToF = np.ones((K,1),dtype=int)*np.arange(1,Nfaces+1).reshape(1,Nfaces);  
    
    id = fnodes[:,0]*Nnodes + fnodes[:,1]+1;    
   
    spNodeToNode=np.zeros((K*Nfaces,4),dtype=np.int);    
    spNodeToNode[:,0]=id;
    spNodeToNode[:,1]=np.arange(1,Nfaces*K+1)
    spNodeToNode[:,2]=EToE.flatten(order='F');
    spNodeToNode[:,3]=EToF.flatten(order='F');
    
    t1=spNodeToNode.view(np.ndarray);
    np.lexsort((t1[:,0],));
    spsorted=t1[np.lexsort((t1[:,0],))];
    
  
    indices =np.where(spsorted[0:-1,0]==spsorted[1:,0]);
    indices = np.array(indices).reshape(np.size(indices),);
    
    matchL = np.concatenate((spsorted[indices],spsorted[indices+1]),axis=0);
    matchR = np.concatenate((spsorted[indices+1],spsorted[indices]),axis=0);              
    [b, c]=np.unravel_index([matchL[:,1]-1],(gg.K,gg.Nfaces),order='F');
    EToE[b,c]=matchR[:,2];
    EToF[b,c]=matchR[:,3];
    return EToE, EToF;
    
def FindLocalCoords2D(k,rout,yout):
    import Glob2 as gg;
    v1=gg.EToV[k,0];v2=gg.EToV[k,1];v3=gg.EToV[k,2];
    xy1 = np.array([VX[v1],VY[v1]]);
    xy2 = np.array([VX[v2],VY[v2]]);
    xy3 = np.array([VX[v3],VY[v3]]);    
    return K,rout,yout;
    
    

def readMesh():
    import gmsh as g;
    obj=g.Mesh();
    title='sample.msh';
    obj.read_msh(title);
    gg.VX=obj.Verts[:,0];
    gg.VY=obj.Verts[:,1];
    t1=obj.nElmts;
    gg.K=t1[2];
    gg.EToV=obj.Elmts[2][1];
    gg.EToV=gg.EToV ;
#    gg.BCType=
    
    
def BuildMaps2D():
    import Glob2 as gg;
    nodeids=np.arange(0,gg.K*gg.Np).reshape(gg.Np,gg.K,order='F');
    gg.vmapM=np.zeros((gg.Nfp,gg.Nfaces,gg.K),dtype=np.int);    
    gg.vmapP=np.zeros((gg.Nfp,gg.Nfaces,gg.K),dtype=np.int);
    gg.mapM =np.arange(gg.K*gg.Nfp*gg.Nfaces);
    b=gg.vmapM.copy();############### very nice concept***********
    gg.mapP =b.reshape((gg.Nfp,gg.Nfaces,gg.K),order='F');
    t1=np.zeros((gg.Nfaces*gg.K),dtype=np.int);
    t2=np.zeros((gg.Nfaces*gg.K),dtype=np.int);    
    for k1 in range (0,gg.K):
        for f1 in range (0,gg.Nfaces):
            gg.vmapM[:,f1,k1]=nodeids[gg.Fmask[:,f1],k1];
    sk=0;            
    one=np.ones((1,gg.Nfp),dtype=np.int);
    for k1 in range (0,gg.K):
        for f1 in range (0,gg.Nfaces):
            k2=gg.EToE[k1,f1];f2=gg.EToF[k1,f1]
            v1=gg.EToV[k1,f1];v2=gg.EToV[k1,np.mod((f1+1),gg.Nfaces)]
            t1[sk]=v1;t2[sk]=v2;sk+=1;   
            refd=np.sqrt( (gg.VX[v1]-gg.VX[v2])**2+(gg.VY[v1]-gg.VY[v2])**2);
         
            gg.vidM=gg.vmapM[:,f1,k1];gg.vidP=gg.vmapM[:,f2,k2];
            x1=np.ravel(gg.x,order='F')[gg.vidM];y1=np.ravel(gg.y,order='F')[gg.vidM];
            x2=np.ravel(gg.x,order='F')[gg.vidP];y2=np.ravel(gg.y,order='F')[gg.vidP];
            x1=np.reshape(x1,(len(x1),1))*one;y1=np.reshape(y1,(len(y1),1))*one;
            x2=np.reshape(x2,(len(x2),1))*one;y2=np.reshape(y2,(len(y2),1))*one;
            D=(x1-x2.T)**2+(y1-y2.T)**2;
            [idM,idP]=np.where(np.sqrt(np.abs(D))< gg.NODETOL*refd);
            #print idM, idP;
            gg.vmapP[idM,f1,k1]=gg.vidP[idP];
            gg.mapP[idM,f1,k1] = idP + (f2)*gg.Nfp+(k2)*gg.Nfaces*gg.Nfp;
            
    gg.vmapP=np.ravel(gg.vmapP,order='F');
    gg.vmapM=np.ravel(gg.vmapM,order='F');
    gg.mapP=np.ravel(gg.mapP,order='F');
    gg.mapB=np.where(gg.vmapP==gg.vmapM);
    gg.vmapB=gg.vmapM[gg.mapB];
    return gg.mapM,gg.mapP,gg.vmapM,gg.vmapP,gg.vmapB,gg.mapB;

def BuildBCMaps2D():
    import Glob2 as gg;
    bct =gg.BCType.T;
    bnodes = np.matrix(np.ones((gg.Nfp,1),dtype=np.int))*np.matrix(bct.ravel(order='F'));
    bnodes=np.ravel(bnodes,order='F');
    gg.mapI= np.where(bnodes==gg.In);gg.vmapI=gg.vmapM[gg.mapI]
    gg.mapO= np.where(bnodes==gg.Out);gg.vmapO=gg.vmapM[gg.mapO]
    gg.mapW= np.where(bnodes==gg.Wall);gg.vmapW=gg.vmapM[gg.mapW]
    gg.mapF= np.where(bnodes==gg.Far);gg.vmapF=gg.vmapM[gg.mapF]
    gg.mapC= np.where(bnodes==gg.Cyl); gg.vmapC=gg.vmapM[gg.mapC]
    gg.mapD= np.where(bnodes==gg.Dirichlet);gg.vmapD=gg.vmapM[gg.mapD]
    gg.mapN= np.where(bnodes==gg.Neuman);gg.vmapN=gg.vmapM[gg.mapN]
    gg.mapS= np.where(bnodes==gg.Slip);gg.vmapS=gg.vmapM[gg.mapS]
    
    

def Startup2D():
    
    gg.NODETOL=1e-12;
  
    gg.Nfp=gg.N+1; gg.Np=(gg.N+1)*(gg.N+2)/2;gg.Nfaces=3;
    
    [x,y]=Nodes2D(gg.N);
    [gg.r,gg.s]=xytors(x,y);
    gg.V=Vandermonde2D(gg.N,gg.r,gg.s);
    gg.V=np.matrix(gg.V);
    gg.invV = la.inv(gg.V);
    gg.MassMatrix = gg.invV.T*gg.invV;
    [gg.Dr,gg.Ds]=Dmatrices2D(gg.N,gg.r,gg.s,gg.V);    
    va=gg.EToV[:,0].T;vb=gg.EToV[:,1].T;vc=gg.EToV[:,2].T;
    gg.x=0.5*(-np.matrix(gg.r+gg.s).T*np.matrix(gg.VX[va])+np.matrix(1+gg.r).T*np.matrix(gg.VX[vb])+np.matrix(1+gg.s).T*np.matrix(gg.VX[vc]));    
    gg.y=0.5*(-np.matrix(gg.r+gg.s).T*np.matrix(gg.VY[va])+np.matrix(1+gg.r).T*np.matrix(gg.VY[vb])+np.matrix(1+gg.s).T*np.matrix(gg.VY[vc]));    
    
    fmask1 = np.where(np.abs(gg.s+1)<gg.NODETOL);
    fmask2 = np.where(np.abs(gg.r+gg.s)<gg.NODETOL);
    fmask3 = np.where(np.abs(gg.r+1)<gg.NODETOL);
    gg.Fmask = np.concatenate([fmask1,fmask2,fmask3]);
    gg.Fmask=gg.Fmask.T;
    gg.LIFT=Lift2D();
    [gg.rx,gg.sx,gg.ry,gg.sy,gg.J]=GeometricFactors2D(gg.x,gg.y,gg.Dr,gg.Ds);
    [gg.nx,gg.ny,gg.sJ]=Normals2D();
    
    gg.Fscale=np.divide(gg.sJ,gg.J[gg.Fmask.flatten(order='F'),:])
    
    [gg.EToE, gg.EToF]=tiConnect2D(gg.EToV);
    gg.EToE-=1;gg.EToF-=1;
    
    BuildMaps2D();
    [Vr, Vs]=GradVandermonde2D(gg.N,gg.r,gg.s);
    Vr=np.matrix(Vr);Vs=np.matrix(Vs);
    
    gg.Drw= (gg.V*Vr.T)*la.inv(gg.V*gg.V.T);
    
    gg.Dsw= (gg.V*Vs.T)*la.inv(gg.V*gg.V.T);
            
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

def CutOffFilter2D(Nc,frac):
    filterdiag = np.ones((gg.Np,1),dtype=np.float);
    
    sk=0;
    for i in range (0,N+1):
        for j in range(0,N+1-i):
            if(i+j>=Nc):
                filterdiag[sk]=frac;
    F=gg.V*np.matrix(np.diag(filtediag))*gg.invV;
    return F;
    
def GaussFaceMesh2D(NGauss):
    import aa;
    
    [gaussZ, gaussW]=aa.JacobiGQ(0,0,NGuass-1);
    face1r=gaussZ;face2r= -gaussZ;face3r = -np.ones((NGuass,1),dtype=np.float);
    face1s=-np.ones((NGuass,1),dtype=np.float);
    face2s=gaussZ;face3s=-guassZ;
    
    guassfinterp=np.zeros(NGauss,gg.Np,gg.Nfaces);
    V1=Vandermonde2D(N,face1r,face1s);
    gaussfinterp[:,:,0]=V1*gg.invV;
    
    V2=Vandermonde2D(N,face2r,face2s);
    gaussfinterp[:,:,1]=V2*gg.invV;
    
    V3=Vandermonde2D(N,face3r,face3s);
    gaussfinterp[:,:,2]=V1*gg.invV;
    
    gaussinterp=np.array((NGauss,Np,3),dtype=np.float);
    guassinterp[:,:,0]=guasfinterp[:,:,0];
    guassinterp[:,:,1]=guasfinterp[:,:,1];
    guassinterp[:,:,2]=guasfinterp[:,:,2];
        
        
              
              
    
            


    
def mshRead(mshfile):       
    
    
    
    elm_type = {}
    elm_type[1] = 2    # 2-node line
    elm_type[2] = 3    # 3-node triangle
    elm_type[3] = 4    # 4-node quadrangle
    elm_type[4] = 4    # 4-node tetrahedron
    elm_type[5] = 8    # 8-node hexahedron
    elm_type[6] = 6    # 6-node prism
    elm_type[7] = 5    # 5-node pyramid
    elm_type[8] = 3    # 3-node second order line
                            # (2 nodes at vertices and 1 with edge)
    elm_type[9] = 6    # 6-node second order triangle
                            # (3 nodes at vertices and 3 with edges)
    elm_type[10] = 9    # 9-node second order quadrangle
                            # (4 nodes at vertices,
                            #  4 with edges and 1 with face)
    elm_type[11] = 10   # 10-node second order tetrahedron
                            # (4 nodes at vertices and 6 with edges)
    elm_type[12] = 27   # 27-node second order hexahedron
                            # (8 nodes at vertices, 12 with edges,
                            #  6 with faces and 1 with volume)
    elm_type[13] = 18   # 18-node second order prism
                            # (6 nodes at vertices,
                            #  9 with edges and 3 with quadrangular faces)
    elm_type[14] = 14   # 14-node second order pyramid
                            # (5 nodes at vertices,
                            #  8 with edges and 1 with quadrangular face)
    elm_type[15] = 1    # 1-node point
    elm_type[16] = 8    # 8-node second order quadrangle
                            # (4 nodes at vertices and 4 with edges)
    elm_type[17] = 20   # 20-node second order hexahedron
                            # (8 nodes at vertices and 12 with edges)
    elm_type[18] = 15   # 15-node second order prism
                            # (6 nodes at vertices and 9 with edges)
    elm_type[19] = 13   # 13-node second order pyramid
                            # (5 nodes at vertices and 8 with edges)
    elm_type[20] = 9    # 9-node third order incomplete triangle
                            # (3 nodes at vertices, 6 with edges)
    elm_type[21] = 10   # 10-node third order triangle
                            # (3 nodes at vertices, 6 with edges, 1 with face)
    elm_type[22] = 12   # 12-node fourth order incomplete triangle
                            # (3 nodes at vertices, 9 with edges)
    elm_type[23] = 15   # 15-node fourth order triangle
                            # (3 nodes at vertices, 9 with edges, 3 with face)
    elm_type[24] = 15   # 15-node fifth order incomplete triangle
                            # (3 nodes at vertices, 12 with edges)
    elm_type[25] = 21   # 21-node fifth order complete triangle
                            # (3 nodes at vertices, 12 with edges, 6 with face)
    elm_type[26] = 4    # 4-node third order edge
                            # (2 nodes at vertices, 2 internal to edge)
    elm_type[27] = 5    # 5-node fourth order edge
                            # (2 nodes at vertices, 3 internal to edge)
    elm_type[28] = 6    # 6-node fifth order edge
                            # (2 nodes at vertices, 4 internal to edge)
    elm_type[29] = 20   # 20-node third order tetrahedron
                            # (4 nodes at vertices, 12 with edges,
                            #  4 with faces)
    elm_type[30] = 35   # 35-node fourth order tetrahedron
                            # (4 nodes at vertices, 18 with edges,
                            #  12 with faces, 1 in volume)
    elm_type[31] = 56   # 56-node fifth order tetrahedron
                            # (4 nodes at vertices, 24 with edges,
                            #  24 with faces, 4 in volume)
    Verts = [];Elmts = {};Phys = {};npts = 0;nElmts = {};nprops = 0
    
    
       # try:
    fid = open(mshfile, "r")
    #except IOError:
           # print "File '%s' not found." % (filename)
        #sys.exit()
    
    #line = 'start'
    #while line:
    line = fid.readline()
    
    if line.find('$MeshFormat') == 0:
        line = fid.readline()
    if line.split()[0][0] is not '2':
        print ('wrong gmsh version')
        sys.exit()
    line = fid.readline()
    if line.find('$EndMeshFormat') != 0:
        raise ValueError('expecting EndMeshFormat')
    line=fid.readline();
    if line.find('$PhysicalNames') == 0:
        line = fid.readline()
        nprops = int(line.split()[0])
        bctag=np.zeros((nprops),dtype=np.int);
        for i in range(0,  nprops):
            line = fid.readline()
            newkey = int(line.split()[1])
            bctag[i]= int(line.split()[1]);
            qstart = line.find('"')+1
            qend = line.find('"', -1, 0)-1
            Phys[i] = line[qstart:qend]
        line = fid.readline()
    if line.find('$EndPhysicalNames') != 0:
        raise ValueError('expecting EndPhysicalNames')
    line=fid.readline();
    if line.find('$Nodes') == 0:
        line = fid.readline()
        npts = int(line.split()[0])
        Verts =np.zeros(( npts, 3), dtype=float)
        for i in range(0,  npts):
            line = fid.readline()
            data = line.split()
            idx = int(data[0])-1  # fix gmsh 1-based indexing
            if i != idx:
                raise ValueError('problem with vertex ids')
            Verts[idx, :] = map(float, data[1:])
        line = fid.readline()
        if line.find('$EndNodes') != 0:
            raise ValueError('expecting EndNodes')
        line=fid.readline();
    if line.find('$Elements') == 0:
        line = fid.readline()
        nel = int(line.split()[0])
        bcidx=np.zeros((nel,),dtype=int);
        for i in range(0,  nel):
            line = fid.readline()
            data = line.split()
            idx = int(data[0])-1  # fix gmsh 1-based indexing
            if i != idx:
                raise ValueError('problem with elements ids')
            etype = int(data[1])           # element type
            nnodes =  elm_type[etype]   # lookup number of nodes
            ntags = int(data[2])           # number of tags following
            bcidx[i]=int(data[3]);                    
            k = 3
            if ntags > 0:                   # set physical id
                physid = int(data[k])
                if physid not in  Phys: #   Phys[physid] = 'Physical Entity %d' % physid
                    nprops += 1
                k += ntags
    
            verts = map(int, data[k:])
            verts = np.array(verts)-1  # fixe gmsh 1-based index
    
            if (etype not in  Elmts) or\
                    (len( Elmts[etype]) == 0):
                # initialize
                 Elmts[etype] = (physid, verts)
                 nElmts[etype] = 1
            else:
                # append
                 Elmts[etype] = \
                    (np.hstack(( Elmts[etype][0], physid)),
                     np.vstack(( Elmts[etype][1], verts)))
                 nElmts[etype] += 1
    
        line = fid.readline()
        if line.find('$EndElements') != 0:
            raise ValueError('expecting EndElements')
    fid.close();
    bcstr=['In','Out','Wall','Far','Cyl','Dirichlet','Neuman','Slip'];
    #bi=np.zeros((8,),dtype=np.int);
    VX=Verts[:,0];VY=Verts[:,1];t1=nElmts;K=t1[2];
    
       # EToV=np.zeros((K,3),dtype=np.int);
    EToV=Elmts[2][1];TwoSidT=Elmts[1][1];Nv=npts;
    BCType=np.zeros((K,3),dtype=np.int);
    t3=np.zeros((8,),dtype=np.int);
    
    
    
    sk=0;
    for i in range(0,8):
        if (sk in Phys):
            if (bcstr[i]==Phys[sk]):   t3[i]=bctag[sk]; sk+=1;
    for k in range (0,len(t3)):
        if (t3[k]>0):
            for i in range(0,len(TwoSidT)):
                for j in range (0,2):
                    [t1,t2]=np.where(EToV==TwoSidT[i][j]);
                    if (bcidx[i]==t3[k]): BCType[t1,t2]=k+1;                
    return Nv, VX, VY, EToV , K, BCType;
                    

                    

    
    
        
        
