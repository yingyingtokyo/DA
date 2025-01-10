#%% to calculate the average simulation based on the restart file

import numpy as np
import os
import netCDF4 as nc
import struct
import params as pm

def read_sim_data(numch):
    year=pm.spinup_end_year()
    month=pm.spinup_end_month()
    date=pm.spinup_end_date()
    yyyy="%04d"%(year)
    mm="%02d"%(month)
    dd="%02d"%(date)
    file = "/work/a06/yingying/camada/HydroDA/src/CaMa_out/"+yyyy+mm+dd+"C"+numch+"/2019-sp1/o_rivdph2019.nc"
    print(file)
    return file

def read_nc(file):
    # 20190301-20190501
    # average of 201904
    # check if the existence of data
    if os.path.exists(file):
        nf = nc.Dataset(file,'r')
        varname = nf.variables.keys()
        varname = list(varname)
        print(varname)
        #[time,lat,lon]
        var = np.array(nf.variables[varname[3]][31:-1,:,:])
        var = np.where(var>10000,np.nan,var)
        # average water surface elevation
        var_mean = np.nanmean(var,axis=0)
        print(var_mean[426,1016])
    return var_mean

file_path = "/work/a06/yingying/camada/HydroDA/src/simmean/"+pm.runname(pm.mode())+"/"
os.makedirs(file_path,exist_ok=True)
for num in range(1,pm.ens_mem()+1):
    numch='%03d'%num
    filename = read_sim_data(numch)
    var_mean = read_nc(filename)
    # # write .bin
    fileens ="C"+numch
    filebin  = file_path+fileens+'.bin'
    with open(filebin,'wb') as file:
        file.write(var_mean)

