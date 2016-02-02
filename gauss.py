# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 13:13:48 2016

@author: niya
"""

import Glob2 as gg;
import numpy as np;
import aa;
import bb;
def GaussFaceMesh2D(NGauss):
    Ng=NGauss;    
    [gZ, gW]= aa.JacobiGQ(0,0,Ng-1);
    face1r=gZ; face2r = -gZ; face3r=-np.ones((NGauss,1),dtype=int);
    face1s=-np.ones((Ng,1),dtype=np.int);
    face2s=gZ; face3s=-gZ;
    
    gfinterp=np.zeros((Ng,gg.Np,gg.Nfaces),dtype=np.float);
    
    V1=bb.Vandermonde2D(gg.N,face1r,face1s);gfinterp[:,:,0]=V1*gg.invV;
    V2=bb.Vandermonde2D(gg.N,face2r,face2s);gfinterp[:,:,1]=V2*gg.invV;

    V3=bb.Vandermonde2D(gg.N,face3r,face3s);gfinterp[:,:,2]=V3*gg.invV;

    ginterp = np.concatenate((gfinterp[:,:,0],gfinterp[:,:,1],gfinterp[:,:,2]),axis=1);
    
    gmapM = np.arange(0,NGauss*gg.Nfaces*gg.K).reshape(NGauss*gg.Nfaces,gg.K);
    gmapP = np.arange(0,NGauss*gg.Nfaces*gg.K).reshape(NGauss*gg.Nfaces,gg.K);
    
    zer=np.zeros(NGauss*gg.Nfaces,gg.K);
    gnx=zer; gny=zer; grx=zer; gry=zer; gsx=zer;gsy=zer; gJ =zer; gsJ=zer;

    for f1 in range (0,gg.Nfaces):
        VM=gfinterp[:,:,f1];
        dVMdr= VM*gg.Dr;dVMds=VM*gg.Ds;
        ids1= (f1-1)*NGauss+range(1,f1+1)*NGauss;
        for k1 in range (0,gg.K):
            [rx,sx,ry,sy,J]=bb.GeometricFactors2D(gg.x[:,k1],gg.y[:,k1],dVMdr,dVMds);
            
            if(f1==0): nx= -sx; ny= -sy; 
            if(f1==1): nx= rx+sx; ny=ry+sy; 
            if(f1==2): nx=-rx; ny= -ry;
            
            sJ = np.sqrt(nx**nx+ny**ny);
            nx = np.divide(nx,sJ);
            ny = np.divide(ny,sJ);
            
            sJ= np.multiply(sJ,J);
            
            gnx[ids1,k1]=nx; gny[ids1,k1]=ny; gsJ[ids1,k1]=sJ;
            grx[ids1,k1]=rx; gry[ids1,k1]=ry; gJ[ids1,k1]=J;
            gsx[ids1,k1]=sx; gsy[ids1,k1]=sy;
            
            k2=gg.EToE(k1,f1);f2=gg.EToF(gg.k1,gg.f1);
            ids2=np.arange(f2*gg.NGauss,(f2-1)*NGauss+1,-1);
            if(k1!=k2):
                gmapP[ids1,k1]=gmapM[ids2,k2];
            else:
                gmapP[ids,k1]=gmapM[ids1,k1];
                if(gg.BCType[k1,f1]==Wall):
                    gmapW=np.vstack((gmapW,gmapM[ids1,k1]));
                if(gg.BCType[k1,f1]==In):
                    gmapW=np.vstack((gmapI,gmapM[ids1,k1]));
                if(gg.BCType[k1,f1]==Out):
                    gmapW=np.vstack((gmapO,gmapM[ids1,k1]));
                if(gg.BCType[k1,f1]==Dirichlet): 
                    gmapW=np.vstack((gmapD,gmapM[ids1,k1]));
                if(gg.BCType[k1,f1]==Slip):
                    gmapW=np.vstack((gmapS,gmapM[ids1,k1]));    
                if(gg.BCType[k1,f1]==Neuman):
                    gmapW=np.vstack((gmapN,gmapM[ids1,k1]));
                if(gg.BCType[k1,f1]==Cyl):
                    gmapW=np.vstack((gmapC,gmapM[ids1,k1]));    
                    
    gx=ginterp*gg.x;
    gy=ginterp*gg.y;
    gW=np.vstack((gW,gW,gW))*np.matrix(np.ones((1,gg.K),dtype=np.float));
    gW=np.multiply(gW,gsJ);
    return 
                    
                    
                
                
            
                        
            
            
            
            