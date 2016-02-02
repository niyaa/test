def vtu_puv_write(title,a,b,c,d,e):
    import Glob2 as gg;
    import numpy as np;
    
    NumberOfPoints=gg.Np;
    NumberOfCells = gg.N*gg.N;
    with open (title,"W") as f:
        f.write('<?xml version="1.0"?>\n');
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
        f.write('  <UnstructuredGrid>\n');
        for ele in range (0,gg.K):
            f.write('<<Piece NumberOfPoints=%s NumberOfCells=%s>'%(str(NumberOfPoints),str(NumberOfCells)));
        
            f.write('\n');
            f.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii"> \n')
            for i in range (0,gg.Np):        
                f.write(str(gg.x[i,ele]),str(gg.y[i,ele]),str('\n'));
        
            f.write('</DataArray>\n');
            
            f.write('</Points>\n');
            f.write('<Cells>\n');
            f.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n');
        
        c=[];
        for i in range(0,gg.N):
            c.append(2*i-1);
            c.reverse();
            b=0;d=1;k=1;p=1;
            d1=[[0,0,0],[0,0,0]];
            for i in range (0,N):
                for j in range (0,c[i]):
                    if(np.remainder(j,2!=0)):
                        d1[0,0]=b;d[1,1]=b+1;d[1,2]=b+gg.N + p;
                        f.write(str(d[0,:]));
                        b=b+1;
                    else:
                        d[1,0]=d[0,2];d[1,1]=d[0,1];d[1,2]=d[1,0]+1;
                        f.write(str(d[1,:]));
                b=b+1;p=p-1;
            
            
            f.write('</DataArray>\n');
            f.write('<DataArray type="Int32" Name="offsets" format="ascii">\n');
            k=3;
            for i in range(0,NumberOfCells):
                f.write(str(k));
                k=k+3;
                
            f.write('</DataArray>\n');
            f.write('<DataArray type="Int32" Name="offsets" format="ascii">\n');
                k=5;
            
            for i in range (0,NumberOfCells):
                f.write(str(k));
            
             f.write('</DataArray>\n');
             f.write('</Cells>\n');
             f.write('<PointData>\n');
             f.write('<DataArray type="Float64" Name="Ux">\n');
             
             for i in range (0,gg.Np):
                 f.write(str(Ux[i,ele]));
             
             f.write('<PointData>\n');
             f.write('<DataArray type="Float64" Name="Uy">\n');
             
             for i in range (0,gg.Np):
                 f.write(str(Uy[i,ele]));
            
             f.write('<PointData>\n');
             f.write('<DataArray type="Float64" Name="Uy">\n');
             
             for i in range (0,gg.Np):
                 f.write(str(Pr[i,ele]))
            
             f.write('<PointData>\n');
             f.write('<DataArray type="Float64" Name="Uy">\n');
             
             for i in range (0,gg.Np):
                 f.write(str(vort[i,ele]));
             
             f.write('\n');
             f.write('</DataArray>\n');
             f.write('</PointData>\n');
             f.write('</Piece>\n');
             
             f.write('  </UnstructuredGrid>\n' );
             f.write('</VTKFile>\n' );
    f.close();
             

                        
            
    