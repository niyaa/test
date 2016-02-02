# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 12:28:08 2016

@author: niya
"""
    mshfile='CyInCH.msh'
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