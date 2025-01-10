#!/opt/local/bin/python
# -*- coding: utf-8 -*-

#libralies
import os
import itertools
import numpy as np
import sys
import errno
from multiprocessing import Pool
from multiprocessing import Process
import datetime
import functools
import numpy.random as rd
import os.path
import datetime as dt
import glob
import shutil
import scipy.linalg as spla
from numpy import ma
import random
import re
import calendar
import math

#external python codes
import params as pm
# ################################################
#
# make folders
# Prepare input data sets [runoff]
# initial inflation parameter rho for assimilation
#
# ################################################
def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
###########################
def slink(src,dst):
    try:
        os.symlink(src,dst)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src,dst)
        else:
            raise
###########################
def initial(): #used
    # program for initialization

    # creating output folders
    dahour = pm.dahour()
    dahour_str = '{:02d}'.format(dahour)
    mkdir("CaMa_out_"+pm.runname(pm.mode())+dahour_str)
    mkdir("CaMa_in_"+pm.runname(pm.mode())+dahour_str)
    mkdir("CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart")
    mkdir("CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart/assim")
    mkdir("CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart/open")
    mkdir("CaMa_in_"+pm.runname(pm.mode())+dahour_str+"/restart/true")
    src_fold = '/work/a06/yingying/CaMa_v411/cmf_v411_pkg/inp'
    link_fold = '/work/a06/yingying/camada/HydroDA/src/CaMa_in_'+pm.runname(pm.mode())+dahour_str+'/inp'
    if os.path.islink(link_fold) or os.path.exists(link_fold):
            os.remove(link_fold)  # Remove the existing link
    os.symlink(src_fold,link_fold)

    mkdir("assim_"+pm.runname(pm.mode())+dahour_str)
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m")
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m/assim")
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m/open")
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/xa_m/true")

    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa")         # NEW v.1.1.0
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/assim")   # NEW v.1.1.0
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/open")    # NEW v.1.1.0
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/meanA")    # NEW v.1.1.0
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/ens_xa/meanC")    # NEW v.1.1.0

    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/nonassim")
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/nonassim/open")
    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/nonassim/assim")

    mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/rest_true")
    
    for var in ["rivdph","sfcelv","outflw","storge"]:
        mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/"+var)
        mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/"+var+"/assim")
        mkdir("assim_"+pm.runname(pm.mode())+dahour_str+"/"+var+"/open")

    # error output folder
    mkdir("err_out")

    # inflation parameter
    mkdir("inflation")

    os.system("touch assim_"+pm.runname(pm.mode())+dahour_str+"/__init__.py")
    mkdir("logout")

    return 0
###########################
def make_initial_infl():
    nx,ny,gsize=pm.map_dimension()
    parm_infl=np.ones([ny,nx],np.float32)*pm.initial_infl()
    start_year,start_month,start_date,start_hour=pm.starttime() # Start year month date
    yyyy='%04d' % (start_year)
    mm='%02d' % (start_month)
    dd='%02d' % (start_date)
    hh='%02d' % (start_hour)
    #hh='%02d' % (start_hour)
    parm_infl.tofile("./inflation/parm_infl"+yyyy+mm+dd+hh+".bin")
    return 0
###########################
def make_anomaly_data(mode=pm.mode()):
    # make directory for mean sfcelv
    mkdir("./assim_out/mean_sfcelv/")
    # copy the anomaly files
    if mode == 1:
        # for mean
        iname=pm.DA_dir()+"/dat/mean_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname="./assim_out/mean_sfcelv/mean_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
        # for std
        iname=pm.DA_dir()+"/dat/std_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname="./assim_out/mean_sfcelv/std_sfcelv.bin"
        os.system("cp "+iname+" "+oname)

    if mode == 2:
        # for mean
        iname=pm.DA_dir()+"/dat/mean_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname="./assim_out/mean_sfcelv/mean_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
        # for std
        iname=pm.DA_dir()+"/dat/std_sfcelv_E2O_"+pm.mapname()+"_1980-2014.bin"
        oname="./assim_out/mean_sfcelv/std_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
    
    if mode == 3:
        # for mean
        iname=pm.DA_dir()+"/dat/mean_sfcelv_VIC_BC_"+pm.mapname()+"_1979-2013.bin"
        oname="./assim_out/mean_sfcelv/mean_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
        # for std
        iname=pm.DA_dir()+"/dat/std_sfcelv_VIC_BC_"+pm.mapname()+"_1979-2013.bin"
        oname="./assim_out/mean_sfcelv/std_sfcelv.bin"
        os.system("cp "+iname+" "+oname)
    return 0
###########################
def save_statistic():
    # copy mean and std of simulated WSE
    # for anomaly and normalized assimilations
    mkdir("./assim_out/mean_sfcelv/")
    # mean
    inputlist=[]
    for ens in np.arange(1,pm.ens_mem(pm.mode())+1):
        ens_char="%03d"%(ens)
        iname = pm.DA_dir()+"/dat/mean_sfcelv_"+ens_char+".bin"
        oname = "./assim_out/mean_sfcelv/meansfcelvC"+ens_char+".bin"
        inputlist.append([iname,oname])

    # do parallel
    p=Pool(pm.para_nums())
    p.map(copy_stat,inputlist)
    p.terminate()
    
    #===========
    # std
    inputlist=[]
    for ens in np.arange(1,pm.ens_mem(pm.mode())+1):
        ens_char="%03d"%(ens)
        iname = pm.DA_dir()+"/dat/std_sfcelv_"+ens_char+".bin"
        oname = "./assim_out/mean_sfcelv/stdsfcelvC"+ens_char+".bin"
        inputlist.append([iname,oname])

    # do parallel
    p=Pool(pm.para_nums())#*cpu_nums())
    p.map(copy_stat,inputlist)
    p.terminate()
    return 0
###########################
def copy_stat(inputlist):
    iname = inputlist[0]
    oname = inputlist[1]
    print ("cp "+iname+" "+oname)
    os.system("cp "+iname+" "+oname)
    return 0
###########################
if __name__ == "__main__":
    initial()
    make_initial_infl()
    save_statistic()
