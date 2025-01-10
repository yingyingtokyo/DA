import numpy as np



filepath = '/work/a06/yingying/camada/HydroDA/src/CaMa_out_hgt/20190112C001/'
filedate  = 'restart2019011300.bin'


def read_bin(filepath,filedate):
    with open(filepath+filedate, 'rb') as file:
        binary_data = np.fromfile(file, dtype=np.float32)
    data_3d = binary_data.reshape((2, 1320, 1500))
    data = data_3d[0,:,:]
    fld  = data_3d[1,:,:]
    return data,fld
data,fld = read_bin(filepath,filedate)
print(data[605,993])
print(fld[605,993])
