# -*- coding: utf-8 -*-
"""
change .bin to .nc
@author: Lyy
"""

import numpy as np
import netCDF4 as nc

exp_name = 'wid'
interval = 12
time_len = int(20*24/interval)
interval_str = '{:02d}'.format(interval) 

#%% read obs
obs_path = '/work/a06/yingying/obs/2019h/'+'wlv'+'/'
def read_obs(day,hr):
    day_n    = '{:02d}'.format(day)
    hr_n     = '{:02d}'.format(hr)
    filename = '201910'+day_n+hr_n+'.bin'
    with open(obs_path+filename, 'rb') as f:
        day_array = np.fromfile(f, dtype=np.float32)
        day_array = day_array.reshape(1320,1500)
        day_array = np.where(day_array<-900,np.nan,day_array)
    f.close()
    return day_array
obs = np.full((time_len,1320,1500),np.nan)

#%%
# read .bin
filepath = "/work/a06/yingying/camada/HydroDA/src/CaMa_in_"+exp_name+interval_str+"/restart/open/restart201901"
dapath   = "/work/a06/yingying/camada/HydroDA/src/CaMa_in_"+exp_name+interval_str+"/restart/assim/restart201901"
wlvpath  = "/work/a06/yingying/camada/HydroDA/src/CaMa_out_"+exp_name+interval_str+"/201901"
onpath  = "/work/a06/yingying/camada/HydroDA/src/CaMa_out_man12/2019011300A002/"
def read_bin(day,hr,filepath,datatype):
    day_name="{:02d}".format(day) 
    hr_name ="{:02d}".format(hr) 
    with open(filepath+day_name+hr_name+datatype+"002.bin", 'rb') as file:
        binary_data = np.fromfile(file, dtype=np.float32)
    data_3d = binary_data.reshape((2, 1320, 1500))
    data = data_3d[0,:,:]
    fld  = data_3d[1,:,:]
    return data,fld
#%%
def read_wlv(day,hr,datatype):
    day_name="{:02d}".format(day)
    hr_name ="{:02d}".format(hr) 
    binf = wlvpath+day_name+hr_name+datatype+"002/rivdph2019.bin"
    with open(binf, 'rb') as file:
        binary_data = np.fromfile(file, dtype=np.float32)
    data_3d = binary_data.reshape((interval,1320, 1500))
    data = data_3d[-1,:,:] 
    return data
#%%
def read_on(filepath):
    with open(filepath+"rivdph2019.bin", 'rb') as file:
        binary_data = np.fromfile(file, dtype=np.float32)
    data = binary_data.reshape((interval,1320, 1500))
    return data

#%%
ind  = 0
for day in range(0+1,20+1):
    for hr in range(0,24,12):
        obs[ind,:,:]  = read_obs(day,24)
        restart,fld = read_bin(day,hr,filepath,"C")
        print('-----C-------')
        #print(restart[605,993])
        print(restart[603,1034])
        #print(fld[605,993])
        print(fld[603,1034])
        if day<20:
            wlv = read_wlv(day,hr,"C")
            print(wlv[603,1034])
            #print(wlv[605,993])
        print('-----A-------')
        restart,fld = read_bin(day,hr,dapath,"A")
        print(restart[603,1034])
        #print(restart[605,993])
        print(fld[603,1034])
        #print(fld[605,993])
        if day<20:
            wlv = read_wlv(day,hr,"A")
            #print(wlv[605,993])
            print(wlv[603,1034])
        print("restart file for day:",day,hr)
        ind = ind+1
    print('------------')
data_on = read_on(onpath)
#print("online da result:",data_on[605,993])
print("1013online da result:",data_on[:,603,1034])


