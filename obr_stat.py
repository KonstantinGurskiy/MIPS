import numpy as numpy
import scipy as scipy
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import spatial
from numpy import sys
#np.genfromtxt("test", usecols=1, dtype=float)



def read_state(f):
    f.readline()                                                                        # ITEM: TIMESTEP
    f.readline()                                                                        # 0
    f.readline()                                                                        # ITEM: NUMBER OF ATOMS
    N=int(f.readline())
    f.readline()                                                                        # ITEM: BOX BOUNDS pp pp pp
    xsize = [float(item) for item in f.readline().split(' ')]
    ysize = [float(item) for item in f.readline().split(' ')]
    zsize = [float(item) for item in f.readline().split(' ')]
    f.readline()                                                                        # ITEM: ATOMS id ... 
    #print(f.readline().strip().split(' ')[2:-1])
    #d=[[float(item) for item in f.readline().strip().split(' ')[2:-1]] for i in range(N)]     # read data  [2:-1]
    d=[[float(item) for item in f.readline().strip().split(' ')[2:4]] for i in range(N)] 
    print(d)
    return [N,xsize,ysize,d]

def PolygonArea(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area



f=open('/run/media/softmatter/Новый том/MIPS/cluster-+-/longer/T_0.1Eps_0.5/T_0.1_rho_0.1_Dr_0.15_gamma_2.0_ts_0.0001_eps_0.5_Act_0.0.lammpstrj', "r")

    
for i_step in range(1):
    st=read_state(f)
    vor=scipy.spatial.Voronoi(st[-1])
    #----------------------------------
    # opredelenie pogranichnux regionov
    #----------------------------------
    xmin=st[1][0]
    xmax=st[1][1]
    ymin=st[2][0]
    ymax=st[2][1]
    ver_index=numpy.full(len(vor.vertices),True)
    for i in range(len(vor.vertices)):
        ver=vor.vertices[i]
        if ver[0]<xmin or ver[0]>xmax or ver[1]<ymin or ver[1]>ymax: 
            ver_index[i]=False
    #----------------------------------
    #----------------------------------
    #----------------------------------


    #----------------------------------
    # otsechenie pogranichnux regionov i poisk min max ploshadei
    #----------------------------------
    region_mask=numpy.full(len(vor.regions),True)
    for jj in range(len(vor.regions)):
        region=vor.regions[jj]
        if not -1 in region and not False in ver_index[region]:
            region_mask[jj]=True
            #polygon = [vor.vertices[i] for i in region]
            #pa=PolygonArea(polygon)
            #if pa>maxpa: maxpa=pa
            #if pa<minpa: minpa=pa
        else:
            region_mask[jj]=False
    #----------------------------------
    #----------------------------------
    #----------------------------------



    #----------------------------------
    # sostavlenie spiska blishaishix sosedei (neigh_id)
    #----------------------------------
    npoints=len(vor.points)
    neigh_id=[[] for i in range(npoints)]
    for rid in vor.ridge_points:
        neigh_id[rid[0]].append(rid[1])
        neigh_id[rid[1]].append(rid[0])
    #for it in neigh_id: 
    #    if len(it)-len(list(set(it)))!=0: print(1)     // proverka na dublikatu
    #----------------------------------
    #----------------------------------
    #----------------------------------



    #----------------------------------
    # vuchislenie pol9 R_field and S_field -- parametr regul9rnosti 9cheek
    #----------------------------------
    R_field=numpy.array([10. for i in range(npoints)])
    S_field=numpy.array([10. for i in range(npoints)])
    R_field_min=10.0
    R_field_max=-1.0
    S_field_min=(xmax-xmin)*(ymax-ymin)
    S_field_max=-1.0
    for jj in range(npoints):
        if region_mask[vor.point_region[jj]]:           # esli region ne otbroshen
            dr_t=vor.points[jj]-vor.points[neigh_id[jj]]
            r=numpy.array([numpy.sqrt(i.dot(i)) for i in dr_t])
            #tt=r.var()/(r.mean()**2)
            tt=r.var()
            R_field[jj]=tt
            if tt<R_field_min: R_field_min=tt
            if tt>R_field_max: R_field_max=tt
            polygon = [vor.vertices[i] for i in vor.regions[vor.point_region[jj]]]
            pa=PolygonArea(polygon)
            S_field[jj]=pa
            if pa<S_field_min: S_field_min=pa
            if pa>S_field_max: S_field_max=pa
        else:
            R_field[jj]=-1.0
            S_field[jj]=-1.0
    #----------------------------------
    #----------------------------------
    #----------------------------------



    #----------------------------------
    # rashet dispersii pol9 R_field
    #----------------------------------
    R_field_D=numpy.array([10. for i in range(npoints)])
    S_field_D=numpy.array([10. for i in range(npoints)])
    S_field_M=numpy.array([10. for i in range(npoints)])
    R_field_D_min=10.0
    R_field_D_max=-1.0
    S_field_D_min=10.0
    S_field_D_max=-1.0
    for jj in range(npoints):
        if region_mask[vor.point_region[jj]]:           # esli region ne otbroshen
            tt=numpy.array(R_field[neigh_id[jj]+[jj]])
            tt=tt[tt>-0.5]
            tt2=tt.var()
            R_field_D[jj]=tt2
            if tt2<R_field_D_min: R_field_D_min=tt2
            if tt2>R_field_D_max: R_field_D_max=tt2   
            tt=numpy.array(S_field[neigh_id[jj]+[jj]])
            tt=tt[tt>-0.5]
            tt2=tt.var()
            S_field_D[jj]=tt2
            S_field_M[jj]=tt.mean()
            if tt2<S_field_D_min: S_field_D_min=tt2
            if tt2>S_field_D_max: S_field_D_max=tt2 
        else:
            R_field_D[jj]=-1.0
            S_field_D[jj]=-1.0
            S_field_M[jj]=-1.0
    #----------------------------------
    #----------------------------------
    #----------------------------------    
            
    
    
    #----------------------------------
    # rashet pol9 Max-Min - Delta
    #----------------------------------        
    Delta_field=numpy.array([10. for i in range(npoints)])   
    Delta_field_max=-1.0
    for jj in range(npoints):
        if region_mask[vor.point_region[jj]]:  
            tt=numpy.array(R_field_D[neigh_id[jj]+[jj]])
            tt=tt[tt>-0.5]
            tt2=tt.max()-tt.min()
            Delta_field[jj]=tt2
            if tt2>Delta_field_max: Delta_field_max=tt2
        else: 
            Delta_field[jj]=-1.0
            
            
         
    #----------------------------------
    # dump to file
    #----------------------------------        
    #print("min R_field_D=",R_field_D_min)
    #print("max R_field_D=",R_field_D_max)
    #print("min S_field_D=",S_field_D_min)
    #print("max S_field_D=",S_field_D_max)
    
    #numpy.savetxt('Out/out_R_field_'+str(i_step)+'.txt',R_field)
    #numpy.savetxt('Out/out_S_field_'+str(i_step)+'.txt',S_field)
    #numpy.savetxt('Out/out_R_field_D_'+str(i_step)+'.txt',R_field_D)
    #numpy.savetxt('Out/out_S_field_D_'+str(i_step)+'.txt',S_field_D)
    #numpy.savetxt('Out/out_S_field_M_'+str(i_step)+'.txt',S_field_M)
    #numpy.savetxt('Out/out_Delta_field_'+str(i_step)+'.txt',Delta_field)
    
    n_RS=200;
    max_RS=10.
    min_RS=10.**(-7)
    log_h=(numpy.log(max_RS)-numpy.log(min_RS))/(n_RS-1)
    Data_RS=numpy.array([0. for i in range(n_RS)])
    Data_RS2=numpy.array([0. for i in range(n_RS)])
    Data_RS_n=numpy.array([0 for i in range(n_RS)])
    Setka_RS=numpy.array([numpy.exp(numpy.log(min_RS)+log_h*i) for i in range(n_RS)])
    for jj in range(len(S_field)):
        if R_field_D[jj]>0:
            index=(numpy.log(R_field_D[jj])-numpy.log(min_RS))//log_h
            if index>=0 and index<n_RS:
                index2=int(index)
                Data_RS_n[index2]+=1
                Data_RS[index2]+=S_field[jj]
    
    for i in range(len(Data_RS)):
        if Data_RS_n[i]>0:
            Data_RS2[i]=Data_RS[i]/Data_RS_n[i]
        else:
            Data_RS2[i]=-1
    

    #----------------------------------
    #----------------------------------
    #----------------------------------
    print(i_step)



for jj in range(len(vor.points)):
    region_id=vor.point_region[jj]
    if region_mask[region_id]:
        polygon = [vor.vertices[i] for i in vor.regions[region_id]]
        #cl=cm.get_cmap('RdYlBu')(R_field[jj]/R_field_max)          # R_field
        if ((R_field_D[neigh_id[jj]]<(0.00180547)).any() and (R_field_D[neigh_id[jj]]>(0.00180547)).any()):
#        if ((R_field_D[neigh_id[jj]]<(1.5*10**(-3.))).any() and (R_field_D[neigh_id[jj]]>(1.5*10**(-3))).any()):
            cl='orange'
        else:
            if R_field_D[jj]<0.00180547:
                cl='blue'
            else: 
                cl='red'
        #cl=cm.get_cmap('RdYlBu')(1.0-100*R_field_D[jj])
        #cl=cm.get_cmap('RdYlBu')(1.0-R_field_D[jj]/R_field_D_max)     # R_field_D
        #cl=cm.get_cmap('RdYlBu')(1.0-R_field_D[jj]/(4.*(10**(-4.))))     # R_field_D
        #cl=cm.get_cmap('RdYlBu')(1.0-100.0*Delta_field[jj]/Delta_field_max)     # Delta_field
        #cl=cm.get_cmap('RdYlBu')(1.0-S_field[jj]/S_field_max)     # S_field
        #cl=cm.get_cmap('RdYlBu')(1.0-S_field_D[jj]/S_field_D_max)     # S_field_D
        #cl=cm.get_cmap('RdYlBu')(1.0-R_field[jj]*S_field[jj]/(S_field_max*R_field_max))     # S_field*R_field
        #cl=cm.get_cmap('RdYlBu')(1.0-R_field_D[jj]*S_field_D[jj]/(S_field_D_max*R_field_D_max))     # S_field_D*R_field_D
        plt.fill(*zip(*polygon),facecolor=cl,edgecolor='black', linewidth=0.25)
plt.show()
