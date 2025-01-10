#%% to calculate the average simulation based on the restart file

import numpy as np
import os
import netCDF4 as nc
import struct
import params as pm

# forecast starts date
date=11
date_str = '{:02d}'.format(date)
# calculate average restart
data_sim=np.zeros((2,1320,1500,pm.ens_mem()))
data_da =np.zeros((2,1320,1500,pm.ens_mem()))

simpath = "/work/a06/yingying/camada/HydroDA/src/CaMa_in/restart/open/restart201901"
dapath   = "/work/a06/yingying/camada/HydroDA/src/CaMa_in/restart/assim/restart201901"
savepath = "/work/a06/yingying/camada/HydroDA/src/forecast/"+pm.runname(pm.mode())+"/"
def read_bin(filepath,ens,datatype):
    ens_str ="{:03d}".format(ens)
    with open(filepath+date_str+datatype+ens_str+".bin", 'rb') as file:
        binary_data = np.fromfile(file, dtype=np.float32)
    data_3d = binary_data.reshape((2, 1320, 1500))
    data = data_3d[0,:,:]
    fld  = data_3d[1,:,:]
    return data,fld

def write_average_bin():
    os.makedirs(savepath,exist_ok=True)
    for num in range(1,pm.ens_mem()+1):
        data_sim[0,:,:,num-1],data_sim[1,:,:,num-1] = read_bin(simpath,num,'C')
        data_da[0,:,:,num-1],data_da[1,:,:,num-1] = read_bin(dapath,num,'A')
    aver_sim = np.nanmean(data_sim,axis=3)
    da_sim   = np.nanmean(data_da,axis=3)
    return aver_sim,da_sim
aver_sim,da_sim = write_average_bin()

# # write .bin
def write_bin(var_each,filename):
    var_flat = var_each.flatten() 
    var_bin  = var_flat.astype('float32').tobytes()
    filesim  = savepath+filename
    with open(filesim,'wb') as file:
        file.write(var_bin)
write_bin(aver_sim,"open"+date_str+".bin")
write_bin(da_sim,"assim"+date_str+".bin")

def check_bin(filename):
    with open(savepath+filename, 'rb') as file:
        binary_data = np.fromfile(file, dtype=np.float32)
    data_3d = binary_data.reshape((2, 1320, 1500))
    data = data_3d[0,:,:]
    fld  = data_3d[1,:,:]
    print(data[583,1005])
    return data,fld 
check_bin("open"+date_str+".bin")
check_bin("assim"+date_str+".bin")
