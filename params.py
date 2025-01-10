#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import re

########################
#
# parameters list
#
########################
# **************************************************************
# 0. HydroDA version
def version():
    return "HydroDA version 1.0.0"
    # version for WSE assimilation / observation localization
    # CaMa-Flood v396a used

# **************************************************************
# 1. experment type related definitions
def mode():
    return 3
    # parameter to change assimilation mode
    # runoff ensembles will change accordingly.
    # 1: ERA5, 2: tej, 3:ils, 4:topography, 5:wid, 6:hgt, 7:man 

def DAhour_setting():
    return "on"

# do 24h forecast
def DA_leading():
    return "on"

def dahour():
    # time interval for DA
    return 1 

def da_ehour():
    # length of simulation time
    if DA_leading()=='on':
        return 24
    # length of simulation time is the same as the DA interval
    else:
        return dahour()

def conflag():
    return 2
    # converstion flag for observation converstions 
    #  1 - Directly values 
    #  2 - Anomalies
    #  3 - Normalized values
    #  4 - Log converted values

def mapname():
    return "tej_01min"
    # return "amz_06min"
    # return "glb_15min"
    # realted CaMa-Flood map directory
    # [e.g. : glb_15min, glb_06min, Mkg_06min, etc.]
    # Check 

def map_dimension():
    fname=CaMa_dir()+"/map/"+mapname()+"/params.txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    #-------
    nx     = int(list(filter(None, re.split(" ",lines[0])))[0])
    ny     = int(list(filter(None, re.split(" ",lines[1])))[0])
    gsize  = float(list(filter(None, re.split(" ",lines[3])))[0])
    return nx,ny,gsize

# **************************************************************
# 2. Data assimilation related definitions
def assimS():
    # return -20
    return 48
    # data Assimilation's Region (South Edge at latitude)
    # *note: should be larger or equal to -80

def assimN():
    # return 5
    return 24
    # data Assimilation's Region (North Edge at latitude)
    # *note: should be smaller or equal to 80

def assimW():
    # return -80
    return 123
    #return -68.25 # use this for disabling west side of the Amazon basin's observation
    # data Assimilation's Region (West Edge at latitude)
    # *note: should be larger or equal to -170

def assimE():
    # return -45
    return 148
    # data Assimilation's Region (East Edge at latitude)
    # *note: should be smaller or equal to 170

def patch_size():
    # return 0
    return 100
    # the size of the local patch of LETKF(Local ** EnKF)
    # 0: only 1 pixel (the pixel itself) belongs to its local patch
    # 100: empirical local patch

def DA_dir():
    return "/work/a06/yingying/camada/HydroDA"
    # directory of HydroDA
    # where src, dat, sat, out exsits

def patch_dir():
    return "/work/a06/yingying/out_wlv2019/obs_patch/"

def patch_name():
    # return "amz_06min_S14FD_40/local_patch"
    # return "tej_01min_S14FD_60"
    return ""

def patch_id():
    # return "0.80"
    return "0.60"
    # return "0.40"
    # return "0.20"

def thersold():
    # return 0.80
    return 0.60
    # return 0.40
    # return 0.20
    # thersold to define the local patch

def initial_infl():
    return 1.08
    # initial inflation parameter

def rho():
    # return -1.0
    return 1.00
    # return 1.08
    # -1.0 : adaptive inflation will be used as in Myoshi et al (2011)
    # positive : fixed inflation parameter will be used
    # [E.g. 1.08, 1.10]

def sigma_b():
    return 0.0400000
    # bacground variance of inflation for adaptive inflation Myoshi et al (2011)

def ens_mem(mode=mode()):
    return 20
    # number of ensemble members

# **************************************************************
# 3. Experiment timimngs 
def timestep():
    return 3600 # outer timestep in seconds

def starttime(DAhour=DAhour_setting()):
    if DAhour == 'on':
        return (2019,1,10,0) # start date: [year,month,date]
    else:
        return (2019,1,1)

def endtime(DAhour=DAhour_setting()):
    if DAhour == 'on':
        return (2019,1,16,23) # end date: [year,month,date]
    else:
        return (2019,1,20)
                      # *note: this date is not included

# **************************************************************
# 4. Spinup options
def spinup_mode():
    return 0
    # 0: do spinup simulation for both (corrupted and true) simulation
    # 1: do spin up only at corrupted simulation
    # 2: do spin up only at true simulation
    # 3: no spinup simulation at all
    ### if initial restart file is ready, spinup simulation is no need

def spinup_end_year():
    return 2019

def spinup_end_month():
    return 12

def spinup_end_date():
    return 31

def calibrate():
    return "no"
    # return "yes"
    # return "corrupt"

# **************************************************************
# 5. Runoff forcing 
def runoff_dir():
    return "/work/a06/yingying/CaMa_v411/cmf_v411_pkg/inp"

def runname(num=mode()):
    if num == 1:
        return "ERA5"

    if num == 2:
        return "tej"

    if num == 3:
        return "ils"

    if num == 4:
        return "tpo"

    if num == 5:
        return "wid"

    if num == 6:
        return "hgt"

    if num == 7:
        return "man"

def input(num=mode()):
    if num==1:
        return "ERA5"

    if num==2:
        return "tej"
    # define the runoff data type.

    if num == 3:
        return "ils"

    if num == 4:
        return "tpo"
    
    if num == 5:
        return "wid"

    if num == 6:
        return "hgt"

    if num == 7:
        return "man"


def runoff_name(num=mode()):
    if num==1:
        return ""
    if num==2:
        return "TE_tej_SAT.nc"
    
def max_lat():
    return 80. # maximum latitude of assimilation
               # *note: SWOT ovservation is not available beyond 80 degs. this should be less or equal to 80
               ## modified 2018-06-05

def distopen(num):
    if num == 1:
        return 1.0

    if num == 2:
        return 1.0

    if num == 3:
        return 1.0
    #return 0.75 # not needed for ERA20CM
    # corrupted runoff's percentage
    # 0.75 for original Data Assimilation simulation (25% reduced)
    # 1.25 for 25% increased simulation
    # 1.00 for simulation using 1 year before runoff
    # *note: also editing and and re-compile of control_inp at CaMa-Flood is nessessary

def diststd(num):
    if num == 1:
        return 0.1

    if num == 2:
        return 0.1

    if num == 3:
        return 0.25

    #return 1.0 # not needed for ERA20CM
    # noise to make runoff input to scatter ensembles

def err_expansion():
    return 1.0
    # variance-covariance expansion
    # works well with 0.04

# # def rivman_error():
# #     return 0
# #     #define the experiment with or without rivman error
# #     # 0 : with out manning error
# #     # 1 : with manning error: Manning's n depend on river width
# #     # 2 : with manning error: Manning's n depend spatial covarience
# #     # 3 : with manning error: Manning's n randomly distributed subbasin
# #     # 4 : with manning error: Manning's n randomly distributed
# #     # 5 : with manning error: Manning's n depend on rivseq
# #     # 6 : with manning error: Manning's n depend on uparea
# #     # 7 : with manning error: Manning's n depend on

def run_flag():
    return 0
    # 0 run all simulations
    # 1 run only corrupted and assimilated simulations
    # 2 run only true and assimilated simulations
    # 3 run only assimilated simulation

def true_run(num):
    if num == 1:
        return 3 # ecmwf as true

    if num == 2:
        return 4 # ERA04 as true

    if num == 3: # get disturb runoff from ILS 
        return 3
# **************************************************************
# 6. CaMa-Flood settings
def CaMa_ver():
    # return "CaMa-Flood version 3.9.6"
    return "CaMa-Flood version 4.1.1"

def CaMa_dir():
    return "/work/a06/yingying/HydroDA/CaMa-Flood-DA_Japan20240520/CaMa_v411"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders

def CaMa_out(num=mode()):
    if num == 2:
        return "CaMa_out"
    if num == 3:
        return "CaMa_out_ils"
    if num == 4:
        return "CaMa_out_tpo"
    if num == 5:
        return "CaMa_out_wid"
    if num == 6:
        return "CaMa_out_hgt"
    if num == 7:
        return "CaMa_out_man"
    else:
        return "CaMa_out"

def CaMa_in(num=mode()):
    if num == 2:
        return "CaMa_in"
    if num == 3:
        return "CaMa_in_ils"
    if num == 4:
        return "CaMa_in_tpo"
    if num == 5:
        return "CaMa_in_wid"
    if num == 6:
        return "CaMa_in_hgt"
    if num == 7:
        return "CaMa_in_man"
    else:
        return "CaMa_in"

def MKLdir():
    return "/opt/intel/compilers_and_libraries_2016.3.170/mac/mkl"
    # directory of Intel MKL files
    # Intel MKL is needed for doing data assimilation
    # Please Download and Instal it to your System before running
    # for more information --> https://software.intel.com/en-us/qualify-for-free-software/academicresearcher

def output_er():
    return 0
    # setting for saving or deleting intermediate files
    # 0 for saving & 1 for deleting
    # those files may be more than 400GB, so erasing is recommended if not necessary

# **************************************************************
# 7. observations settings

def obs_name():
    # return "HydroWeb"
    return "MLIT"

def HydroWeb_dir():
    return "/cluster/data6/menaka/HydroWeb"

def obs_dir():
    return "/work/a06/yingying/obs/2019UST/wlv"

def make_log():
    return 1
    # setting for making log files
    # 1 is for making and 0 is for not making

def slack_notification():
    return 0
    # setting for validating slack notification
    # 1 for valid and 0 for invalid
    # 0 is a default if you are not familiar with slack
    # if you turn it to 1, you need to edit sendslack.py
    # for more information refer https://api.slack.com/incoming-webhooks

######### Not currently in use
# # def ens_at_non():
# #     return 1
# #     # * At Recent version, ensemble generating random number is constant for full simulation.
# #     # (For example, when ensemble 001 is corrupted with -0.1 at day 1, ensemble 001 will be always corrupted with 0.1 for full simulation period.)
# #     # Previously, ensemble mean was used as an assimilated value for non-observed location.
# #     # In this version, this treatment has changed and non-observed location is given with an ensemble value.
# #     # To enable this new feature, set the return of params.py method “ens_at_non()”, “1”(DEFAULT SETTING).
# #     # If you don’t want to use this, set it to “0”.

# # # functions for corrupting manning coeffcient ###################
# # # this is for corrupting manning coefficient at Corrupted Simulation
# # # manning coefficient will be corrupted with random numbers generated from following functions
# # # the random number is generated for each ensemble member
# # # random number is made by gaussian noise of average = corruptman_base(), stddev = corruptman_std()
# # def corruptman_base():
# #     return 0.03

# # def corruptman_std():
# #     return 0.015

# # def rivman_base():
# #     return 0.03

# # def rivman_min():
# #     return 0.025

# # def rivman_max():
# #     return 0.035

# # def corruptele_base():
# #     return 0.5 # not needed for ERA20CM

# # def corruptele_std():
# #     return 0.25 # not needed for ERA20CM

# # def non_hgt():
# #     return 7.0 # not needed for ERA20CM
# #     # nominal water height

# **************************************************************
# 8. parallel run settings
def para_nums():
    return 5 
    # setting number of parallels to run CaMa-Flood Model
    # defualt is 6, but may change depending on your system

def cpu_nums():
    return 20 
    # number of cpus used
