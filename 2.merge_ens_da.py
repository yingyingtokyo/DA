# merge ensemble DA results
import numpy as np
import netCDF4 as nc
import os
import re
import pandas as pd

ens_num = 20
da_path = '/work/a06/yingying/camada/HydroDA/src/assim_out/ens_xa/meanA/201901'
def read_bin(day):
    fday = '{:02d}'.format(day)
    with open(da_path+fday+'_xa.bin', 'rb') as f:
        day_array = np.fromfile(f, dtype=np.float32)
    day = day_array.reshape((1320,1500))
    day = np.where(day<-999,np.nan,day)
    return day

day = read_bin(12)
print(day[421-1,1095-1])
exit()
for day in range(1,20):
    da_data = np.full((1320,1500,ens_num),np.nan,dtype=np.float32)
    for ens in range(1,ens_num+1):
        da_data[:,:,int(ens-1)] = read_bin(day,ens)
    da_mean = np.nanmean(data_data,axis=2)
    write_bin(day,da_mean)

