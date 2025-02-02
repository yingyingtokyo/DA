#!/opt/local/bin/python
# hourly on-line DA code
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
import netCDF4 as nc

#external python codes
import params as pm
#import src.letkf_lib as lb
#import src.calc_one_daybef as odb


########################
#
# main program
#
########################

#main code for LETKF
############
## main control function
############
def main_act():

    print (pm.runname(pm.mode()))

    # Set Time
    print ("Set Time")
    timestep=pm.timestep() #time step for assimilation
    start_year,start_month,start_date=pm.starttime() # Start year month date
    end_year,end_month,end_date=pm.endtime() # End year month date

    # Spin-up Simulation
    #print ("spin up simulation")
    #spin_up()

    # preapre the mean and std for anomaly/normalized assimilation
    print ("calculation of mean/std rivdph")
    calc_mean()

    # make initial restart
    print ("make intial restart")
    make_initial_restart()
    
    start_dt=datetime.date(start_year,start_month,start_date)
    end_dt=datetime.date(end_year,end_month,end_date)
    days_count=(end_dt-start_dt).days
    
    # run daily simulations
    for day in range(days_count): #Moiz: Day np.int64 doesn't work with timedelta
        running_dt=start_dt+datetime.timedelta(days=day)
        yyyy='%04d' % (running_dt.year)
        mm='%02d' % (running_dt.month)
        dd='%02d' % (running_dt.day)
        hh='%02d' % (running_dt.hour)

        one_day_loop(yyyy,mm,dd,day)
    
    # clean all intermediate files
    if pm.output_er()==1:
        os.system("rm -Rf ./CaMa_out/"+yyyy+"*")
        os.system("rm -Rf ./CaMa_in/"+pm.input()+"/Roff_CORR/Roff__"+yyyy+"*")
    
################
## single loop program
################
def one_day_loop(yyyy,mm,dd,day):
    print ("============== start loop of "+yyyy+" "+mm+" "+dd+" =======================")
    # Corrupted Simulation (Open Loop) ###################################################
    if pm.run_flag() == 0 or pm.run_flag() == 1:
        ODM_inputlist=[]
        #randlist=np.fromfile("./CaMa_out/randlist.bin",np.float32)

        # set for ensemble simulations
        ens_num=1
        for ens_num in np.arange(1,pm.ens_mem()+1):
            ODM_inputlist.append([yyyy,mm,dd,'%03d'%ens_num,"open"])

        # Run CaMa-Flood Model (ensemble simulations)
        p=Pool(pm.para_nums())
        p.map(one_day_sim,ODM_inputlist)
        p.terminate()

        # copy ensemble corrupted simulation forecasts from CaMa_out (xa)
        cprestart=[]
        for ens_num in np.arange(1,pm.ens_mem()+1):
            cprestart.append([yyyy,mm,dd,ens_num])
        p=Pool(pm.para_nums())
        p.map(copy_corrupted_sfcelv,cprestart)
        p.terminate()

        p=Pool(pm.para_nums())
        p.map(copy_corrupted_restart,cprestart)
        p.terminate()


    # Assimilated Simulation #############################################################
    ODM_inputlist=[]
    #randlist=np.fromfile("./CaMa_out/randlist.bin",np.float32)

    # set for ensemble simulations
    for ens_num in np.arange(1,pm.ens_mem()+1):
        ODM_inputlist.append([yyyy,mm,dd,'%03d'%ens_num,"assim"])


    # Run CaMa-Flood Model (ensemble simulations)
    p=Pool(pm.para_nums())
    p.map(one_day_sim,ODM_inputlist)
    p.terminate()

    # Calculate Ensemble Mean
    #os.system("./src/make_nonassim "+yyyy+mm+dd+" "+str(pm.ens_mem())+" "+"A")

    # make forecasted value for assimilated simulation
    # do assimilation (LETKF)
    data_assim(yyyy,mm,dd)

    # make restart MODIFIED v.1.1.0
    mkrestart=[]
    for ens_num in np.arange(1,pm.ens_mem()+1):
        mkrestart.append([yyyy,mm,dd,"assim",'%03d'%ens_num])

    # Modify the restart file
    p=Pool(pm.para_nums())
    p.map(make_restart,mkrestart)
    p.terminate()

    # store river variable files
    store_out(yyyy,mm,dd)

#    # make rivout @menaka
#    mkrivout=[]
#    for ens_num in np.arange(1,pm.ens_mem()+1):
#        mkrivout.append([yyyy,mm,dd,"assim",'%03d'%ens_num])
#    p=Pool(pm.para_nums())
#    p.map(make_rivout,mkrivout)
#    p.terminate()

    # clean files
    if pm.output_er()==1:
        bef_dt=datetime.date(int(yyyy),int(mm),int(dd))-datetime.timedelta(days=1)
        bef_yyyy='%04d' %bef_dt.year
        bef_mm='%02d' %bef_dt.month
        bef_dd='%02d' %bef_dt.day
        if bef_dt != datetime.date(bef_dt.year,12,31): # to keep restart file in last day of the year
            os.system("rm -Rf ./CaMa_out/"+bef_yyyy+bef_mm+bef_dd+"*")


#######################################################################################


############
## main program functions
############
def spin_up(): #used
    # run spin up simulation
    # 1 year spin up for calculating initial value
    # one simulation for true
    # ensmble simulation for open

    dir2=pm.CaMa_dir()
    cpunums = pm.cpu_nums()
    yyyy = "%04d"%(pm.spinup_end_year())
    print (pm.spinup_mode())
    if pm.spinup_mode()==3 or pm.spinup_mode()==4: # added flag 4 to copy restart
        return 0

    inputlist=[]

    if pm.spinup_mode()==0 or pm.spinup_mode()==1:
        for ens_num in np.arange(1,pm.ens_mem()+1):
            inputlist.append([yyyy,"open",'%03d'%ens_num])

    # Run spinup simulations
    p=Pool(pm.para_nums())
    p.map(spinup_loop,inputlist)
    p.terminate()

    print("======================= end spinup ==========================")

    return 0
###########################
def spinup_loop(inputlist):
    # Run spinup simulation
    yyyy=inputlist[0]
    loop=inputlist[1]
    ens_num=inputlist[2]
    dir2=pm.CaMa_dir()
    cpunums=pm.cpu_nums()
    mode=pm.mode()
    run_name=pm.runname(mode)
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    mapname=pm.mapname()
    cal=pm.calibrate()
    print  ("%s for %03d"%(loop,int(ens_num)))
    os.system("source "+pm.DA_dir()+"/src/spin_up.sh "+str(yyyy)+" "+str(loop)+" "+ens_num+" "
    +dir2+" "+str(cpunums)+" "+str(run_name)+" "+str(exp_dir)+" "+str(mapname)+" "+str(cal))
    return 0
###########################
def one_day_sim(inputlist):
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    ens_num=inputlist[3]
    looptype=inputlist[4]
    mode=pm.mode()
    run_name=pm.runname(mode)

    # program for running one day model

    bef_dt=datetime.date(int(yyyy),int(mm),int(dd))-datetime.timedelta(days=1)
    bef_yyyy='%04d' %bef_dt.year
    bef_mm='%02d' %bef_dt.month
    bef_dd='%02d' %bef_dt.day

    print ("oneday loop for",yyyy,mm,dd,ens_num,looptype)
    dir2=pm.CaMa_dir()

    print (yyyy+" "+mm+" "+dd+" "+ens_num+" "+dir2+" "+looptype)
    cpunums = pm.cpu_nums() # Moiz: OMP takes integers
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    mapname=pm.mapname()
    cal=pm.calibrate()
    DA_dir=pm.DA_dir()
    os.system("source "+pm.DA_dir()+"/src/oneday_sim.sh "+yyyy+" "+mm+" "+dd+" "+ens_num+" "+dir2
    +" "+looptype+" "+str(cpunums)+" "+str(run_name)+" "+str(exp_dir)+" "+str(mapname)+" "+str(cal)
    +" "+DA_dir) #+str(corrupt)+" "

    return 0
########################### # read nc file
def read_nc(file,nc_start):
    # 20190301-20190501
    # average of 201904
    # check if the existence of data
    if os.path.exists(file):
        nf = nc.Dataset(file,'r')
        varname = nf.variables.keys()
        varname = list(varname)
        print(varname)
        #[time,lat,lon]
        var = np.array(nf.variables[varname[3]][nc_start:-1,:,:])
        var = np.where(var>10000,np.nan,var)
        var = np.where(var<-100,np.nan,var)
        # average water surface elevation
        var_mean = np.nanmean(var,axis=0)
        print(var_mean[426,1016])
    return var_mean
########################### # modified to run paralle @Menaka 
def copy_corrupted_sfcelv(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1] 
    dd   = inputlist[2]
    num  = inputlist[3]
    numch='%03d'%num
    fname="./CaMa_out/"+yyyy+mm+dd+"C"+numch+"/sfcelv"+yyyy+".bin"
    os.system("cp "+fname+" ./assim_out/ens_xa/open/"+yyyy+mm+dd+"_"+numch+"_xa.bin")
    return 0
########################### # modified not to calculate restart again/ no chage in WSE in corrupted @Menaka
def copy_corrupted_restart(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    num  = inputlist[3]
    nxt_day = datetime.date(int(yyyy),int(mm),int(dd)) + datetime.timedelta(days=1)
    n_yyyy='%04d' % (nxt_day.year)
    n_mm='%02d' % (nxt_day.month)
    n_dd='%02d' % (nxt_day.day)
    numch='%03d'%num
    fname="./CaMa_out/"+yyyy+mm+dd+"C"+numch+"/restart"+n_yyyy+n_mm+n_dd+".bin"
    os.system("cp "+fname+" ./CaMa_in/restart/open/restart"+n_yyyy+n_mm+n_dd+"C"+numch+".bin") ## CaMa-Flood v4.1
    print ("copy restart",n_yyyy,n_mm,n_dd,"C"+numch)
    return 0
###########################
def copy_stoonly(iname,oname): # for CaMa_Flood v395b
    org=np.fromfile(iname,np.float32).reshape(6,-1)
    org[0:2].tofile(oname)
    return 0
###########################
def assim_at_fort(yyyy,mm,dd,day): #previous --> used
    dir1=pm.CaMa_dir()+"/"
    thisday=datetime.date(int(yyyy),int(mm),int(dd))
    nxt_day=thisday+datetime.timedelta(days=1)
    os.system(pm.DA_dir()+"/src/data_assim "+str(pm.assimN())+" "+str(pm.assimS())+" "+\
    str(pm.assimW())+" "+str(pm.assimE())+" "+yyyy+mm+dd+" "+str('%02d'%SWOT_day(yyyy,mm,dd))\
    +" "+str(pm.patch_size())+" "+str(pm.patch_name())+" "+str(pm.ens_mem())+" "+str(day)+" "\
    +str('%04d'%(nxt_day.year)+'%02d'%(nxt_day.month)+'%02d'%(nxt_day.day))+" "+str(pm.err_expansion())+" "+dir1)
    return 0
###########################
def data_assim(yyyy,mm,dd): # new data assimilation function (2020/05/18)
    parallels="%d"%(pm.para_nums()*pm.cpu_nums())
    os.environ['OMP_NUM_THREADS']=parallels
    exp_dir="./" 
    print ('****DA****')
    thisday=datetime.date(int(yyyy),int(mm),int(dd))
    nxt_day=thisday+datetime.timedelta(days=1)
    nyear=nxt_day.year
    nmon=nxt_day.month
    nday=nxt_day.day
    nxtyyyymmdd="%04d%02d%02d"%(nyear,nmon,nday)
    os.system(pm.DA_dir()+"/src/data_assim "+yyyy+mm+dd+" "\
    +str(pm.patch_size())+" "+str(pm.ens_mem(pm.mode()))+" "+nxtyyyymmdd+" "+pm.CaMa_dir()\
    +" "+str(pm.thersold())+" "+exp_dir+" "+pm.DA_dir()+" "+pm.patch_dir()+" "\
    +" "+pm.obs_dir()+" "+str(pm.rho())+" "+str(pm.sigma_b())\
    +" "+str(pm.conflag())+" "+(pm.calibrate())+" "+(pm.runname(pm.mode())))
    return 0
###########################
def make_init_storge():
    bef_yyyy='%04d' % (pm.spinup_end_year())
    bef_mm='%02d' % (pm.spinup_end_month())
    bef_dd='%02d' % (pm.spinup_end_date())
    os.system("./src/make_nonassim_init ./CaMa_out/spinup_open/storge"+str(pm.spinup_end_year())+".bin "+"./assim_out/nonassim/open/nonasmC"+bef_yyyy+bef_mm+bef_dd+".bin")
    os.system("./src/make_nonassim_init ./CaMa_out/spinup_open/storge"+str(pm.spinup_end_year())+".bin "+"./assim_out/nonassim/assim/nonasmA"+bef_yyyy+bef_mm+bef_dd+".bin")
    return 0
###########################
def make_initial_restart(): # updated the name
    start_year,start_month,start_date=pm.starttime()
    yyyy="%04d"%(start_year)
    mm_restart = "%02d"%(5)
    mm="%02d"%(start_month)
    dd="%02d"%(start_date)
    exp_dir="./" 
    if pm.spinup_mode()==4:
        return 0
    inputlist=[] # parallel 
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        spinup_open="%04d%2d%02dC%03d"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date(),num)
        infile_opn =exp_dir+"CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin"
        infile_asm =exp_dir+"CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin"
        outfile_opn=exp_dir+"CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin"
        outfile_asm=exp_dir+"CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin"
        inputlist.append([infile_opn,infile_asm,outfile_opn,outfile_asm])
        os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm_restart+dd+".bin  ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin")
        os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
    return 0
###########################
def make_initial_restart_one(): # updated the name
    # copy restartyyyymmddC001 as restart for all simulations 
    start_year,start_month,start_date=pm.starttime()
    yyyy="%04d"%(start_year)
    mm="%02d"%(start_month)
    dd="%02d"%(start_date)
    spinup_true="%04d%2d%02dT000"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date())
    copy_stoonly("./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    print("cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    spinup_open=spinup_true
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
###########################
def copy_out(loop): #old
    os.system("./src/copy_out "+loop)
    return 0
###########################
def reset_loop(yyyy): #old 
    os.system("rm -f ./assim_out/nonassim/*")
    os.system("rm -fR ./CaMa_out/global_15min_"+yyyy+"*")
###########################
def make_rand(std): #used
    # prepare random numbers for ensemable simulation
    randlist=rd.normal(0,std,pm.ens_mem())
    randlist=randlist.astype(np.float32)
    randlist.tofile("./CaMa_out/randlist.bin")
###########################
def initial(): #used
    # program for initialization

    # creating output folders
    mkdir("CaMa_out")
    mkdir("CaMa_in")
    mkdir("CaMa_in/restart")
    mkdir("CaMa_in/restart/assim")
    mkdir("CaMa_in/restart/open")
    mkdir("CaMa_in/restart/true")
    mkdir("assim_out")
    mkdir("assim_out/xa_m")
    mkdir("assim_out/xa_m/assim")
    mkdir("assim_out/xa_m/open")
    mkdir("assim_out/xa_m/true")

    mkdir("assim_out/ens_xa")         # NEW v.1.1.0
    mkdir("assim_out/ens_xa/assim")   # NEW v.1.1.0
    mkdir("assim_out/ens_xa/open")    # NEW v.1.1.0

    mkdir("assim_out/nonassim")
    mkdir("assim_out/nonassim/open")
    mkdir("assim_out/nonassim/assim")

    mkdir("assim_out/rest_true")
    mkdir("assim_out/rivout")
    mkdir("assim_out/rivout/open")
    mkdir("assim_out/rivout/assim")
    mkdir("assim_out/rivout/true")
    mkdir("assim_out/outflw")
    mkdir("assim_out/outflw/open")
    mkdir("assim_out/outflw/assim")
    mkdir("assim_out/outflw/true")
    mkdir("assim_out/fldout")
    mkdir("assim_out/fldout/open")
    mkdir("assim_out/fldout/assim")
    mkdir("assim_out/fldout/true")
    mkdir("assim_out/flddph")
    mkdir("assim_out/flddph/open")
    mkdir("assim_out/flddph/assim")
    mkdir("assim_out/flddph/true")
    mkdir("err_out")
    mkdir("assim_out/fldarea/")
    mkdir("assim_out/fldarea/open")
    mkdir("assim_out/fldarea/assim")
    mkdir("assim_out/fldarea/true")
    # inflation parameter
    mkdir("inflation")

    os.system("touch assim_out/__init__.py")
    mkdir("logout")
    return 0
###########################
def compile_func(): # old used
    # program for compiling
    # activate ifort
    os.system("ifort "+pm.DA_dir()+"/src/make_restart.f90 -o "+pm.DA_dir()+"/src/make_restart -O2 -assume byterecl -heap-arrays -nogen-interfaces -free -g -traceback  -lpthread -parallel")
    print ("compile data assimilation codes...")
    os.system("ifort "+pm.DA_dir()+"/src/data_assim.f90 -o "+pm.DA_dir()+"/src/data_assim -O3 -assume byterecl -heap-arrays -nogen-interfaces -free -mkl=parallel -g -traceback  -lpthread -parallel") #-openmp
    if pm.patch_size()==0:
         os.system("ifort "+pm.DA_dir()+"/src/data_assim_0.f90 -o "+pm.DA_dir()+"/src/data_assim -O3 -assume byterecl -heap-arrays -nogen-interfaces -free -mkl -g -traceback  -lpthread -parallel")
    return 0
###########################
def store_out(yyyy,mm,dd): # update on 2023/05/30
    # program for storing data #
    
    listCA = ["open","assim"]
    if pm.run_flag() == 3:
        listCA = ["assim"]
    #==========================#
    for looptype in listCA: # update on 2023/05/30
        if looptype == "open":
            CA = "C"
        else:
            CA = "A"

        for num in np.arange(1,pm.ens_mem()+1):
            numch = '%03d' % num 
            for var in ["rivdph","sfcelv","outflw","storge"]: 
                shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/"+var+yyyy+".bin","./assim_out/"+var+"/"+looptype+"/"+var+yyyy+mm+dd+"_"+numch+".bin")

#    looptype = "true"
#    # storing rivout
#    numch = "000" 
#    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/rivout"+yyyy+".bin","assim_out/rivout/"+looptype+"/rivout"+yyyy+mm+dd+".bin")

    # storing rivdph
    #shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/rivdph"+yyyy+".bin","assim_out/rivdph/"+looptype+"/rivdph"+yyyy+mm+dd+".bin")

#    # storing fldout
#    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/fldout"+yyyy+".bin","assim_out/fldout/"+looptype+"/fldout"+yyyy+mm+dd+".bin")
#
#    # storing flddph
#    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/flddph"+yyyy+".bin","assim_out/flddph/"+looptype+"/flddph"+yyyy+mm+dd+".bin")
#
#    # storing fldarea
#    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/fldare"+yyyy+".bin","assim_out/fldarea/"+looptype+"/fldarea"+yyyy+mm+dd+".bin")


#        if CA == "C": 
#        # storing rivout
#            for num in np.arange(1,pm.ens_mem()+1):
#                numch = '%03d' % num 
#                shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/rivdph"+yyyy+".bin","assim_out/rivdph/"+looptype+"/rivdph"+yyyy+mm+dd+"_"+numch+".bin")
        

#         # storing rivout
#         for num in np.arange(1,pm.ens_mem()+1):
#             numch = '%03d' % num 
#             shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/rivout"+yyyy+".bin","./assim_out/rivout/"+looptype+"/rivout"+yyyy+mm+dd+"_"+numch+".bin")

#         # storing outflw
#         for num in np.arange(1,pm.ens_mem()+1):
#             numch = '%03d' % num 
#             shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/outflw"+yyyy+".bin","./assim_out/outflw/"+looptype+"/outflw"+yyyy+mm+dd+"_"+numch+".bin")

# #        # storing fldout
# #        for num in np.arange(1,pm.ens_mem()+1):
# #            numch = '%03d' % num 
# #            shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/fldout"+yyyy+".bin","assim_out/fldout/"+looptype+"/fldout"+yyyy+mm+dd+"_"+numch+".bin")
# #
#         # storing flddph
#         for num in np.arange(1,pm.ens_mem()+1):
#             numch = '%03d' % num 
#             shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/flddph"+yyyy+".bin","assim_out/flddph/"+looptype+"/flddph"+yyyy+mm+dd+"_"+numch+".bin")

#         # storing fldarea
#         for num in np.arange(1,pm.ens_mem()+1):
#             numch = '%03d' % num 
#             shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/fldare"+yyyy+".bin","assim_out/fldarea/"+looptype+"/fldarea"+yyyy+mm+dd+"_"+numch+".bin")

    return 0
###########################    
def SWOT_day(yyyy,mm,dd):
    st_year,st_month,st_date=pm.starttime()
    start_time=datetime.date(st_year,st_month,st_date)
    this_time=datetime.date(int(yyyy),int(mm),int(dd))
    days=this_time-start_time
    days=days.days
    return days%21+1
###########################
def make_restart(inputlist):
    # new version
    # sfcelv >> restart
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    loop=inputlist[3]
    numch=inputlist[4]

    # make restart file
    # all in 4 bytes float (small endian) 1440*720
    # items are placed in the following order:
    #       river storage, floodplain storage, river outflow
    #       floodplain outflow, river depth, flood plain storage

    # built in hold
    print ("finish assimilating")
    print ("built in hold")
    print ("press enter")

    # get the date of one day before
    bef_y=calc_odb(yyyy,mm,dd,"year")
    bef_m=calc_odb(yyyy,mm,dd,"month")
    bef_d=calc_odb(yyyy,mm,dd,"date")

    yyyy_b="%04d"%bef_y
    mm_b="%02d"%bef_m
    dd_b="%02d"%bef_d

    # get the date of one day after
    nowdate=dt.datetime(int(yyyy),int(mm),int(dd))
    nextdate=nowdate+dt.timedelta(days=1)
    yyyy_n="%04d"%nextdate.year
    mm_n="%02d"%nextdate.month
    dd_n="%02d"%nextdate.day

    # calculate other variables from water storage
    exp_dir="./" #pm.DA_dir()+"/out/"+pm.experiment()
    os.system(pm.DA_dir()+"/src/make_restart "+yyyy+mm+dd+" "+yyyy_b+mm_b+dd_b+" "
    +yyyy_n+mm_n+dd_n+" "+loop+" "+pm.CaMa_dir()+" "+pm.mapname()+" "+str(pm.ens_mem())+" "
    +numch+" "+exp_dir)

    print ("finish restarting",numch)
###########################
def make_rivout(inputlist):
    # new version
    # sfcelv >> rivout
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    loop=inputlist[3]
    numch=inputlist[4]

    #print "calculate insantaneous discharge using Manning's equation"

    # calculate other variables from water storage
    os.system("./src/make_rivout "+yyyy+mm+dd+" "+loop+" "+pm.CaMa_dir()+" "+pm.mapname()+" "+str(pm.ens_mem())+" "+numch)
###########################
def copy_runoff(inputlist): #do it parallel
    iname=inputlist[0]
    oname=inputlist[1]
    distopen=float(inputlist[2])
    runoff=np.fromfile(iname,np.float32)*distopen
    runoff.tofile(oname)
    return 0
###########################
def make_corrupt_man_old():
    # prepare random numbers for ensemable simulation
    manrandlist=rd.normal(pm.corruptman_base(),pm.corruptman_std(),pm.ens_mem())
    manrandlist=manrandlist.astype(np.float32)

    f=open("./CaMa_out/manrandlist.txt","w")
    for i in np.arange(0,pm.ens_mem()):
        f.write(str(manrandlist[i]))
        f.write("\n")
    f.close()

    return 0
###########################
def make_corrpt_rivhgt():
    # prepare random numbers for ensemable simulation
    elerandlist=rd.normal(pm.corruptele_base(),pm.corruptele_std(),pm.ens_mem())
    elerandlist=elerandlist.astype(np.float32)

    f=open("./CaMa_out/elerandlist.txt","w")
    for i in np.arange(0,pm.ens_mem()):
        f.write(str(elerandlist[i]))
        f.write("\n")
    f.close()

    return 0
###########################
def intial_assim_rivhgt():
    print("initialize elevation")
    for ens in np.arange(1,pm.ens_mem()+1):
      #shutil.copy(pm.CaMa_dir()+"/map/global_15min/rivhgt_%03dC.bin"%(ens),pm.CaMa_dir()+"/map/global_15min/rivhgt19910101_%03dA.bin"%(ens))
      shutil.copy(pm.CaMa_dir()+"/map/global_15min/rivhgt_%03dC.bin"%(ens),pm.CaMa_dir()+"/map/global_15min/rivhgt_%03dA.bin"%(ens))

    return 0
###########################
def courrpt_rivhgt():
    print("make_corrupt_rivhgt.f90")
    os.system("./src/make_corrupt_rivhgt "+str('%4d'%(pm.spinup_end_year()))+" "+str(pm.non_hgt)+" "+str(pm.ens_mem())+" "+pm.CaMa_dir())
    return 0
###########################
def make_initial_infl():
    nx,ny,gsize=pm.map_dimension()
    parm_infl=np.ones([ny,nx],np.float32)*pm.initial_infl()
    start_year,start_month,start_date=pm.starttime() # Start year month date
    yyyy='%04d' % (start_year)
    mm='%02d' % (start_month)
    dd='%02d' % (start_date)
    parm_infl.tofile(pm.DA_dir()+"/out/"+pm.experiment()+"/inflation/parm_infl"+yyyy+mm+dd+".bin")
###########################
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
def multivariate_normal_sampler(mean,covariance,n_samples=1):
  # create multivarite samples 
  L=spla.cholesky(covariance)
  Z=np.random.normal(size=(n_samples,covariance.shape[0]))
  return Z.dot(L) + mean
###########################
def make_corrupt_man():
  # make multivariate normal distribution if corrupted manning
  rivnum = "./data/rivnum.bin"
  rivnum = np.fromfile(rivnum,np.int32).reshape(720,1440)
  #--
  corr_ens=np.ones([pm.manning_mem(),720,1440],np.float32)
  print(np.max(rivnum),np.shape(corr_ens))
  for num in np.arange(1,np.max(rivnum)):
    print(num)
    index=np.where(rivnum==num)
    l=np.shape(index)[1]
    #--
    fname="temp.txt"
    f = open(fname,"w")
    #--
    for i in np.arange(0,l):
      #print index[0][i], index[1][i]
      line="%04d    %04d\n"%(index[1][i]+1, index[0][i]+1)
      f.write(line)
    f.close()
    # find covariances
    print("/src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/")
    os.system("./src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/")
    # multivariate covariance
    print("multivariate covariance")
    mean=np.ones([l],np.float32)*pm.corruptman_base()
    covariance=np.fromfile("cov.bin",np.float32).reshape(l,l)
    corr_ens[:,index[0],index[1]]=multivariate_normal_sampler(mean,covariance,pm.manning_mem())
    #corr_ens[:,index[0],index[1]]=np.random.multivariate_normal(mean,covariance,pm.manning_mem())
   #--
  print("make ensembles")
  for ens in np.arange(1,pm.manning_mem()+1):
    fname="CaMa_out/corruptman%03d.bin"%(ens)
    corr_ens[ens-1].tofile(fname) 
    fname=pm.CaMa_dir()+"/map/glb_15min/rivmanC%03d.bin"%(ens)
    corr_ens[ens-1].tofile(fname)
  return 0
###########################
def make_corrlated_man():
  # make spatially correalted manning
  rivnum = "data/rivnum.bin"
  rivnum = np.fromfile(rivnum,np.int32).reshape(720,1440)
#--
  corr_ens=np.ones([720,1440],np.float32)
  print(np.max(rivnum),np.shape(corr_ens))
  for num in np.arange(1,np.max(rivnum)):
    print(num)
    index=np.where(rivnum==num)
    l=np.shape(index)[1]
    #--
    fname="temp.txt"
    f = open(fname,"w")
    #--
    for i in np.arange(0,l):
      #print index[0][i], index[1][i]
      line="%04d    %04d\n"%(index[1][i]+1, index[0][i]+1)
      f.write(line)
    f.close()
    # find covariances
    print("/src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/")
    os.system("./src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/")
    # calculate covariance
    print("calculate covariance")
    mean=np.ones([l],np.float32)*pm.corruptman_base()
    covariance=np.fromfile("cov.bin",np.float32).reshape(l,l)
    #--
    # Compute the Cholesky decomposition
    c=cholesky(covariance, lower=True)
    #--
    corr_ens[index[0],index[1]]=np.dot(c,mean)

  fname=pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin"
  corr_ens.tofile(fname)
  return 0
###########################
def make_corrupt_man_simple():
  # prepare random numbers for ensemable simulation
  manrandlist=rd.normal(pm.corruptman_base(),pm.corruptman_std(),pm.ens_mem())
  manrandlist=manrandlist.astype(np.float32)
  # ---
  for ens in np.arange(1,pm.manning_mem()+1):
     fname=pm.CaMa_dir()+"/map/glb_15min/rivmanC%03d.bin"%(ens)
     manrand=np.ones([720,1440],np.float32)*manrandlist[ens-1]
     manrand.tofile(fname)
  return 0
###########################
def cal_monthly_mean_ens(ens_num): #start_year,end_year,ens_num,months=24):
    # calc monthly mean value for two years
    start_year=pm.spinup_end_year()
    end_year,end_month,end_date=pm.starttime()
    #start_year=pm.start_year() #inputlist[0]
    #end_year=pm.end_year() #inputlist[1]
    #ens_num=inputlist[2]
    months=24 #inputlist[3]
    #threshold=inputlist[4]
    runname=pm.runname(pm.mode())
    if runname=="E2O":
        nx=1440
        ny=720
        threshold=0.1
    else:
        nx=360
        ny=180
        threshold=1.0e-8
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir("./CaMa_in/"+runname+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    #threshold=0.1

    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    ens_char="%03d"%(ens_num)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile(pm.DA_dir()+"/inp/"+runname+"/Roff/Roff__"+str(yyyy)+str(mm)+str(dd)+ens_char+".one",np.float32).reshape([ny,nx])
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        roff_mean.tofile("./CaMa_in/"+runname+"/mean_month/mean_"+ychar+mchar+ens_char+".bin")
###########################
def cal_monthly_total(start_year,end_year,months=24,threshold=0.1):
    # calc monthly mean value for two years
    runname=pm.runname(pm.mode())
    true_run=pm.true_run(pm.mode())
    if runname=="ELSE_KIM2009":
        nx,ny=360,180
        prefix="Roff____"
        suffix=".one"
    if runname=="E2O":
        nx,ny=1440,720
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    if runname=="ERA20CM":
        nx,ny=360,180
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    mkdir("./CaMa_in/"+runname+"/total_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    #os.system("rm -Rf ./CaMa_in/"+runname+"/total_month/*")
    #threshold=0.1
    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        #count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile("./CaMa_in/"+runname+"/Roff/"+prefix+str(yyyy)+str(mm)+str(dd)+suffix,np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            #count=count+(roff>threshold)
        roff_total=roff_mon #/(count+1e-20)
        roff_total=roff_total.astype(np.float32)
        roff_total=roff_total+threshold
        roff_total.tofile("./CaMa_in/"+runname+"/total_month/total_"+ychar+mchar+".bin")
###########################
def cal_monthly_mean(start_year,end_year,months=24):
    # calc monthly mean value for two years
    runname=pm.runname(pm.mode())
    true_run=pm.true_run(pm.mode())
    if runname=="ELSE_KIM2009":
        nx,ny=360,180
        prefix="Roff____"
        suffix=".one"
    if runname=="E2O":
        nx,ny=1440,720
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    if runname=="ERA20CM":
        nx,ny=360,180
        prefix="Roff__"
        suffix="%03d.one"%(true_run)
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir("./CaMa_in/"+runname+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    threshold=0.1

    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile("./CaMa_in/"+runname+"/Roff/"+prefix+str(yyyy)+str(mm)+str(dd)+suffix,np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        roff_mean.tofile("./CaMa_in/"+runname+"/mean_month/mean_"+ychar+mchar+".bin")
###########################
def direct_insert(yyyy,mm,dd,day):
    oday="%02d"%(SWOT_day(int(yyyy),int(mm),int(dd)))
    mesh_in=np.zeros([720,1440],np.float32)
    rivwth=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivwth.bin",np.float32).reshape(720,1440)
    mesh_in[40:680,:]=np.fromfile("data/mesh_day"+oday+".bin",np.float32).reshape(640,1440)
    mesh=(rivwth>50.0)*(mesh_in>=10)*(mesh_in<=60)*1.0
    mesh=mesh.astype(np.float32)
    #--
    rand=rd.normal(0,pm.ovs_err(),pm.ens_mem())
    rand=rand.astype(np.float32)
    sfcelv=np.fromfile("CaMa_out/"+yyyy+mm+dd+"T000/sfcelv"+yyyy+".bin",np.float32).reshape(720,1440)
    for ens_num in np.arange(1,pm.ens_mem()+1):
        ens_char="%03d"%(ens_num)
        forcst=np.fromfile("CaMa_out/"+yyyy+mm+dd+"A"+ens_char+"/sfcelv"+yyyy+".bin",np.float32).reshape(720,1440) 
        assim=forcst*(1.0-mesh)+((sfcelv*mesh))#+rand[ens_num-1])
        print(np.shape(assim))
        assim.tofile("assim_out/ens_xa/assim/"+yyyy+mm+dd+"_"+ens_char+"_xa.bin")
#################################################################################
def make_rivman():
    nx=1440
    ny=720
    if pm.rivman_error()==0:
        # copy rivman.bin for both rivmanTRUE.bin and rivmanCORR.bin
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin")
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanCORR.bin")
        # rivman
        rivman=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivman.bin",np.float32).reshape(ny,nx)
    elif pm.rivman_error()==1:
        # for rivwth depend calculation Pedonetti et al 2014
        base_man=pm.rivman_base()
        nmin=pm.rivman_min()#0.025
        nmax=pm.rivman_max()#0.035
        #---
        rivwth=pm.CaMa_dir()+"/map/glb_15min/rivwth_gwdlr.bin"
        rivwth=np.fromfile(rivwth,np.float32).reshape(ny,nx)

        # rivnum
        rivnum=pm.DA_dir()+"/dat/rivnum.bin"
        rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)

        # river mouth pixel
        rivmth={}
        fname=pm.DA_dir()+"/dat/rivmth.txt"
        f = open(fname,"r")
        lines = f.readlines()
        f.close()
        #---
        for line in lines[1::]:
            line    = list(filter(None, re.split(" ",line)))
            riverid = float(line[0])
            lon     = int(line[1])
            lat     = int(line[2])
            uparea2 = float(line[3])/1000000.0 # km^2
            rivmth[riverid]=[lon,lat,uparea2]

        # rivman
        rivman=np.ones([ny,nx],np.float32)*base_man
        #print np.amax(rivnum)
        for riv in np.arange(1,1000+1): #np.amax(rivnum)+1):
            #print riv
            num=float(riv)
            #river=ma.masked_where(rivnum!=riv,subbasin).filled(-9999.0)
            width=ma.masked_where(rivnum!=riv,rivwth).compressed()
            wmax=np.amax(width)
            wmin=np.amin(width)
            # get river mouth
            ix=rivmth[riv][0]-1
            iy=rivmth[riv][1]-1
            wmax=rivwth[iy,ix]
            #print wmin, wmax
            if wmin == wmax:
                continue
            index=np.where(rivnum==riv)
            for i in np.arange(len(index[0])):
                ix=index[1][i]
                iy=index[0][i]
                w=rivwth[iy,ix]
                rivman[iy,ix]=max(nmin+(nmax-nmin)*((wmax-w)/((wmax-wmin)+1.0e-20)),nmin)
        rivman.tofile(pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin")
        # copy rivman.bin as rivmanCORR.bin
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanCORR.bin")
    elif pm.rivman_error()==2:
        # rivman calculation considering the spatial covariance
        #nx=1440
        #ny=720
        base_man=pm.rivman_base()
        nmin=pm.rivman_min()#0.025
        nmax=pm.rivman_max()#0.035
        man=0.003
        method='cholesky'
        # subbasin
        subbasin="../dat/subbasin.bin"
        subbasin=np.rint(np.fromfile(subbasin,np.float32).reshape(ny,nx)*1.0e3)
        subbasin=subbasin.astype(int)
        # rivnum
        rivnum="../dat/rivnum.bin"
        rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)#*1e3
        # rivman
        rivman=np.ones([ny,nx],np.float32)*base_man
        rivman1=np.ones([ny,nx],np.float32)*base_man
        print(np.amax(rivnum))
        for riv in np.arange(1,np.amax(rivnum)+1):
            #print riv
            num=float(riv)
            river=ma.masked_where(rivnum!=riv,subbasin).filled(-9999.0)
            #print np.amax(river)#-int(round(num*1.0e3)) #int(round(np.amax(river)*1000.0)),int(round(num*1000.0)),int(round(np.amax(river)*1000.0))-int(round(num*1000.0))
            if np.amax(river) < 0:
                subs=1
            else:
                subs=np.amax(river)-int(round(num*1.0e3)) + 1
            print(subs)
            rdlist=np.random.normal(1.0,0.25,subs)
            #print rdlist
            for i in np.arange(0,subs):
                basin=int(round(num*1.0e3)+int(i))
                print(basin)
                index=np.where(subbasin==basin)
                rivman=ma.masked_where(subbasin==basin,rivman).filled(rdlist[i])
                rivman1=ma.masked_where(subbasin==basin,rivman1).filled(rdlist[i])
                if len(index[0]) != 0:
                    #print "zero pixel"
                    #continue
                    #print index[0],index[1],len(index[0])
                    #x=np.abs(norm.rvs(size=(len(index[0])))+1)
                    #x=np.random.normal(rdlist[i],0.1,len(index[0]))
                    x=np.ones([len(index[0])],np.float32)#*100.0#*rdlist[i]
                    print(len(x))
                    C=cov(index[0],index[1])#,0.01)
                    #print C
        #            try: #if method== 'cholesky':
        #                # cholesky decomposition
        #                c=cholesky(C,lower=True)
        #            except: #else:
        #                print "eigon values"
        #                # eigenvalues and eigenvectors 
        #                evals, evecs = eigh(C)
        #                c=np.dot(evecs, np.diag(np.sqrt(np.ma.masked_less(evals,0.0).filled(0.0))))
        #                c=np.nan_to_num(c)
                    #---
                    # eigenvalues and eigenvectors 
                    evals, evecs = eigh(C)
                    c=np.dot(evecs, np.diag(np.sqrt(np.ma.masked_less(evals,0.0).filled(0.0))))
                    c=np.nan_to_num(c)
                    print(rdlist[i], np.mean(np.abs(np.dot(c,x)*base_man*rdlist[i])))
                    rivman[index[0],index[1]]=np.abs(np.dot(c,x)*base_man)
        #--save rivman
        rivman.tofile(pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin")
        # copy rivman.bin as rivmanCORR.bin
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanCORR.bin")
        #rivman.tofile("rivman.bin")
        #rivman1
    elif pm.rivman_error()==3:
        # rivman calculation considering the subbasins simple random error
        #nx=1440
        #ny=720
        base_man=pm.rivman_base()
        nmin=pm.rivman_min()#0.025
        nmax=pm.rivman_max()#0.035
        man=0.003
        method='cholesky'
        # subbasin
        subbasin="../dat/subbasin.bin"
        subbasin=np.rint(np.fromfile(subbasin,np.float32).reshape(ny,nx)*1.0e3)
        subbasin=subbasin.astype(int)
        # rivnum
        rivnum="../dat/rivnum.bin"
        rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)#*1e3
        # rivman
        rivman=np.ones([ny,nx],np.float32)*base_man
        rivman1=np.ones([ny,nx],np.float32)*base_man
        print(np.amax(rivnum))
        for riv in np.arange(1,np.amax(rivnum)+1):
            #print riv
            num=float(riv)
            river=ma.masked_where(rivnum!=riv,subbasin).filled(-9999.0)
            #print np.amax(river)#-int(round(num*1.0e3)) #int(round(np.amax(river)*1000.0)),int(round(num*1000.0)),int(round(np.amax(river)*1000.0))-int(round(num*1000.0))
            if np.amax(river) < 0:
                subs=1
            else:
                subs=np.amax(river)-int(round(num*1.0e3)) + 1
            print(subs)
            rdlist=np.random.normal(1.0,0.25,subs)
            #print rdlist
            for i in np.arange(0,subs):
                basin=int(round(num*1.0e3)+int(i))
                print(basin)
                index=np.where(subbasin==basin)
                rivman=ma.masked_where(subbasin==basin,rivman).filled(rdlist[i])
                rivman1=ma.masked_where(subbasin==basin,rivman1).filled(rdlist[i])
                if len(index[0]) != 0:
                    #print "zero pixel"
                    #continue
                    #print index[0],index[1],len(index[0])
                    #x=np.abs(norm.rvs(size=(len(index[0])))+1)
                    #x=np.random.normal(rdlist[i],0.1,len(index[0]))
                    x=np.ones([len(index[0])],np.float32)#*100.0#*rdlist[i]
                    print(len(x))
                    C=cov(index[0],index[1])#,0.01)
                    #print C
        #            try: #if method== 'cholesky':
        #                # cholesky decomposition
        #                c=cholesky(C,lower=True)
        #            except: #else:
        #                print "eigon values"
        #                # eigenvalues and eigenvectors 
        #                evals, evecs = eigh(C)
        #                c=np.dot(evecs, np.diag(np.sqrt(np.ma.masked_less(evals,0.0).filled(0.0))))
        #                c=np.nan_to_num(c)
                    #---
                    # eigenvalues and eigenvectors 
                    evals, evecs = eigh(C)
                    c=np.dot(evecs, np.diag(np.sqrt(np.ma.masked_less(evals,0.0).filled(0.0))))
                    c=np.nan_to_num(c)
                    print(rdlist[i], np.mean(np.abs(np.dot(c,x)*base_man*rdlist[i])))
                    rivman[index[0],index[1]]=np.abs(np.dot(c,x)*base_man)
        #--save rivman
        rivman.tofile(pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin")
        # copy rivman.bin as rivmanCORR.bin
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanCORR.bin")
    elif pm.rivman_error()==4:
        #simple random normal distribution
        base_man=pm.rivman_base()
        nmin=pm.rivman_min()#0.025
        nmax=pm.rivman_max()#0.035
        # random values
        # considering the 3-sigma==99.7% of the range
        sigma=(nmax-base_man)/3.0
        rivman=np.random.normal(base_man,sigma,nx*ny).reshape(ny,nx)
        rivman=rivman.astype(np.float32)
        #--save rivman
        rivman.tofile(pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin")
        # copy rivman.bin as rivmanCORR.bin
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanCORR.bin")
    elif pm.rivman_error()==5:
        # rivman calcaulation according to rivseq
        base_man=pm.rivman_base()
        nmin=pm.rivman_min()#0.025
        nmax=pm.rivman_max()#0.035
        #---
        rivseq=pm.CaMa_dir()+"/map/glb_15min/rivseq.bin"
        rivseq=np.fromfile(rivseq,np.int32).reshape(ny,nx)

        # rivnum
        rivnum=pm.DA_dir()+"/dat/rivnum.bin"
        rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)

        # river mouth pixel
        rivmth={}
        fname=pm.DA_dir()+"/dat/rivmth.txt"
        f = open(fname,"r")
        lines = f.readlines()
        f.close()
        #---
        for line in lines[1::]:
            line    = list(filter(None, re.split(" ",line)))
            riverid = float(line[0])
            lon     = int(line[1])
            lat     = int(line[2])
            uparea2 = float(line[3])/1000000.0 # km^2
            rivmth[riverid]=[lon,lat,uparea2]

        # rivman
        rivman=np.ones([ny,nx],np.float32)*base_man
        #print np.amax(rivnum)
        for riv in np.arange(1,1000+1): #np.amax(rivnum)+1):
            #print riv
            num=float(riv)
            #river=ma.masked_where(rivnum!=riv,subbasin).filled(-9999.0)
            seq=ma.masked_where(rivnum!=riv,rivseq).compressed()
            smax=np.amax(seq)
            smin=np.amin(seq)
            # get river mouth
            ix=rivmth[riv][0]-1
            iy=rivmth[riv][1]-1
            smax=rivseq[iy,ix]
            #print wmin, wmax
            if smin == smax:
                continue
            index=np.where(rivnum==riv)
            for i in np.arange(len(index[0])):
                ix=index[1][i]
                iy=index[0][i]
                s=rivseq[iy,ix]
                rivman[iy,ix]=max(nmin+(nmax-nmin)*((smax-s)/((smax-smin)+1.0e-20)),nmin)
        rivman.tofile(pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin")
        # copy rivman.bin as rivmanCORR.bin
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanCORR.bin")
    elif pm.rivman_error()==6:
        # rivman calcaulation according to uparea
        base_man=pm.rivman_base()
        nmin=pm.rivman_min()#0.025
        nmax=pm.rivman_max()#0.035
        #---
        uparea=pm.CaMa_dir()+"/map/glb_15min/uparea.bin"
        uparea=np.fromfile(uparea,np.float32).reshape(ny,nx)

        # rivnum
        rivnum=pm.DA_dir()+"/dat/rivnum.bin"
        rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)

        # river mouth pixel
        rivmth={}
        fname=pm.DA_dir()+"/dat/rivmth.txt"
        f = open(fname,"r")
        lines = f.readlines()
        f.close()
        #---
        for line in lines[1::]:
            line    = list(filter(None, re.split(" ",line)))
            riverid = float(line[0])
            lon     = int(line[1])
            lat     = int(line[2])
            uparea2 = float(line[3])/1000000.0 # km^2
            rivmth[riverid]=[lon,lat,uparea2]

        # rivman
        rivman=np.ones([ny,nx],np.float32)*base_man
        #print np.amax(rivnum)
        for riv in np.arange(1,1000+1): #np.amax(rivnum)+1):
            #print riv
            num=float(riv)
            #river=ma.masked_where(rivnum!=riv,subbasin).filled(-9999.0)
            upa=ma.masked_where(rivnum!=riv,uparea).compressed()
            umax=np.amax(upa)
            umin=np.amin(upa)
            # get river mouth
            ix=rivmth[riv][0]-1
            iy=rivmth[riv][1]-1
            umax=uparea[iy,ix]
            #print wmin, wmax
            if umin == umax:
                continue
            index=np.where(rivnum==riv)
            for i in np.arange(len(index[0])):
                ix=index[1][i]
                iy=index[0][i]
                u=uparea[iy,ix]
                rivman[iy,ix]=max(nmin+(nmax-nmin)*((umax-u)/((umax-umin)+1.0e-20)),nmin)
        rivman.tofile(pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin")
        # copy rivman.bin as rivmanCORR.bin
        os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+pm.CaMa_dir()+"/map/glb_15min/rivmanCORR.bin")
    # rivman to assim_out/rivman
    exp_dir=pm.DA_dir()+"/out/"+pm.experiment()
    mkdir(exp_dir+"/assim_out/rivman")
    # rivman courrpted
    os.system("cp "+pm.CaMa_dir()+"/map/glb_15min/rivman.bin "+exp_dir+"/assim_out/rivman/rivmanCORR.bin")
    # rivman true
    rivman.tofile(exp_dir+"/assim_out/rivman/rivmanTRUE.bin")
###################################################################
def cov(ylist,xlist,sigma=1.0,T=1000.0):
    """covsriacne depend on the catersian distance"""
    C=np.zeros([len(xlist),len(xlist)],np.float32)
    for i in range(len(xlist)):
        ix=xlist[i]
        iy=ylist[i]
        for j in range(len(xlist)):
            iix=xlist[j]
            iiy=ylist[j]
            r=math.sqrt((ix-iix)**2+(iy-iiy)**2)
            C[i,j]=abs(sigma*math.exp((-r**2)/T))
    return C
###################################################################
def observation_error():
    """observation error of WSE depending on the L*W of each pixel
    used sigma*(1/l)*(1/w) l=k*L, w=q*W  Rodrigaz et al 2017:
    According to CaMa k=0.25, q=0.85"""
    k=1.00 # assume nearest part to the unit catchment
    q=1.00 # used 1.0 -> river width variability is 30%
    rivlen=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivlen.bin",np.float32).reshape(720,1440)
    rivwth=np.fromfile(pm.CaMa_dir()+"/map/glb_15min/rivwth_gwdlr.bin",np.float32).reshape(720,1440)
    nextx=(np.fromfile(pm.CaMa_dir()+"/map/glb_15min/nextxy.bin",np.int32).reshape(2,720,1440)[0]!=-9999)*1.0
    rivlen=1.0 #rivlen*1.0e-3 #used as one kilmeter
    rivwth=rivwth*1.0e-3
    area=(k*rivlen)*(q*rivwth)
    obs_err=pm.ovs_err()*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx
    #obs_err=pm.ovs_err()*(1/(area+1.0e-20))*nextx
    # if water area < 1.0 km2 -> 0.25
    obs_err=obs_err*(area>=1.0)*1.0+0.25*(1/(k*rivlen+1.0e-20))*(1/(q*rivwth+1.0e-20))*nextx*(area<1.0)*1.0
    obs_err=ma.masked_where(area<0.625,obs_err).filled(0.25) # 25cm for < 1km^2 water area
    obs_err=obs_err*nextx
    obs_err=obs_err.astype(np.float32)
    obs_err.tofile("data/obs_err.bin")
    return 0
#####################################################################
def err_rand():
    """make random values to add to true values"""
    obs_err=np.fromfile(pm.DA_dir()+"/sat/obs_err.bin",np.float32).reshape(720,1440)
    obs_err=obs_err*((obs_err<=0.25)*1.0) + 0.25*((obs_err>0.25)*1.0)
    #zeros=np.zeros([720,1440],np.float32)
    obs_err=obs_err.flatten()
    rand=[np.random.normal(0.0,obs_err[i],1) for i in np.arange(0,len(obs_err))]
    rand=np.array(rand,np.float32)
    rand.astype(np.float32)
    rand=rand.reshape(720,1440)
    #rand=np.random.normal(zeros,obs_err)
    rand.tofile("CaMa_out/errrand.bin")
    return 0
###########################
def prepare_input():
    # spinup start_year
    # simulation end_year
    start_year=pm.spinup_end_year()
    end_year,end_month,end_date=pm.endtime()
    start_dt=datetime.date(start_year,1,1)
    last_dt=datetime.date(end_year,end_month,end_date)
    start=0
    last=int((last_dt-start_dt).days)
    #---------
    # ERA5
    if pm.input()=="ERA5": #ECMWF ERA20CM
        #nXX,nYY=1440,720
        distopen=pm.distopen(2)
        diststd=pm.diststd(2)
        true_run=pm.true_run(2) # for true ensemble
        runname=pm.runname(2) # get runoff name
        # make courrpted runoff
        if(len(glob.glob("./CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff_CORR")
        #--
        # dist std for ensembles
        distopen_ranges={}
        #open_list=np.setdiff1d(np.arange(1,10+1),[true_run])
        open_list=np.arange(1,10+1)
        random_runs=random.sample(open_list,k=2)
        for runens in open_list:
            ens_char="%03d"%(runens)
            diststd_num=2
            distopen_range=rd.normal(1,diststd,diststd_num)
            distopen_range=np.sort(distopen_range)
            distopen_range=distopen_range.astype(np.float32)
            distopen_ranges[ens_char]=distopen_range#[0.25,1.00,1.25]
        #
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_num=1
            for runens in open_list: #np.arange(1,10+1):
                run_num="%03d"%(runens)
                iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff__"+yyyy+mm+dd+run_num+".one"
                distopens=distopen_ranges[run_num]
                for dist in distopens:
                    ens_char="C%03d"%(ens_num)
                    oname="./CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                    inputlist.append([iname,oname,str(dist)])
                    ens_num=ens_num+1
        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()
    #---------
    # tej
    if pm.input()=="tej": #tej
        nXX,nYY=1500,1320
        distopen=pm.distopen(2)
        diststd=pm.diststd(2)
        true_run=pm.true_run(2) # for true ensemble
        runname=pm.runname(2) # get runoff name
        # make courrpted runoff
        if(len(glob.glob("./CaMa_in/inp/"+runname+"/Runoff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        runoff_ori = read_nc('/work/a06/yingying/Runoff/'+runname+"/"+runoff_name,0)
        mean_runoff= np.nanmean(runoff_ori,axis=0)
        print(mean_runoff[426,1016])
        mkdir("./CaMa_in/inp/"+runname+"/Roff_random")
        #--
        # dist std for ensembles
        distopen_range=np.zeros([np.shape(runoff_ori)[0],nYY,nXX,pm.ens_mem()],np.float32)
        for iXX in range(nXX):
            for iYY in range(nYY):
                    #distopen_range[:,iYY,iXX]=np.sort(rd.normal(distopen,std_runoff[iYY,iXX],pm.ens_mem()))
                    #Log-normal model
                    #sk=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
                    sk=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
                    beta=0.0
                    E=std_runoff[iYY,iXX]/(mean_runoff[iYY,iXX]+1.0e-20)
                    #E=diststd
                    #distopen_range[mon,:,iYY,iXX]=((1+beta)/math.sqrt(E**2+1))*np.exp(math.sqrt(math.log(E**2+1))*sk)
                    #distopen_range[mon-1,:,iYY,iXX]=np.sort(rd.normal(distopen,E,pm.ens_mem()))
                    distopen_range[mon-1,:,iYY,iXX]=rd.normal(distopen,E,pm.ens_mem())
                    #distopen_range[:,iYY,iXX]=np.sort(rd.normal(distopen,diststd,pm.ens_mem()))
        #----------
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            mon=int(target_dt.month)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".one"
            #print iname
            roff=np.fromfile(iname,np.float32).reshape(nYY,nXX)
            #roff=np.ones([nYY,nXX],np.float32)*-9999.0
            fname="../../dat/std_runoff_E2O_1980-2000.bin"
            #print fname
            std_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            fname="../../dat/mean_runoff_E2O_1980-2000.bin"
            #print "L1481",fname
            mean_runoff=np.fromfile(fname,np.float32).reshape(nYY,nXX)
            #print mean_runoff
            std_runoff=ma.masked_where(std_runoff==-9999.0,std_runoff).filled(0.0)
            mean_runoff=ma.masked_where(mean_runoff==-9999.0,mean_runoff).filled(0.0)
            for ens_num in np.arange(pm.ens_mem()):
                ens_char="C%03d"%(ens_num+1)
                roffc=roff+distopen_range[mon-1,ens_num,:,:]*mean_runoff #*std_runoff #
                oname="./CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                roffc.tofile(oname)
    #--------------
    # -25% biased runoff experiment
    if pm.input()=="ELSE": # Kim 2009/E2O/ERA20CM
        #print "-25% biased runoff experiment", pm.true_run(3), pm.runname(3)
        distopen=pm.distopen(3)
        diststd=pm.diststd(3)
        true_run=pm.true_run(3) # for true ensemble
        runname=pm.runname(3) # get runoff name
        # copy for TRUE simulation
        if(len(glob.glob("./CaMa_in/"+runname+"/Roff_TRUE/Roff*"))!=0):
            pass
        # true input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff_TRUE")
        #print true_run
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_char="T000"
            true_char="%03d"%(true_run)
            if runname=="ELSE_KIM2009":
                nx,ny=360,180
                prefix="Roff____"
                suffix=".one"
            if runname=="E2O":
                nx,ny=1440,720
                prefix="Roff__"
                suffix="%03d.one"%(true_run)
            if runname=="ERA20CM":
                nx,ny=360,180
                prefix="Roff__"
                suffix="%03d.one"%(true_run)
            iname=pm.DA_dir()+"/inp/"+runname+"/Roff/"+prefix+yyyy+mm+dd+suffix
            oname="./CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
            inputlist.append([iname,oname,"1.00"])

        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()

        # calculate mothly mean
        cal_monthly_mean(start_year,end_year,24)

        # calculate monthly total
        cal_monthly_total(start_year,end_year,24)

        # make courrpted runoff
        if(len(glob.glob("./CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff_CORR")

        inputlist=[]
        std=rd.normal(0,diststd,pm.ens_mem())
        std=std.astype(np.float32)
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            # make random values
            #std=np.fromfile("./CaMa_out/randlist.bin",np.float32)
            #std=rd.normal(0,diststd,pm.ens_mem())
            #std=std.astype(np.float32)
            #print std
            ifile=np.fromfile("./CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+"T000.one",np.float32).reshape(ny,nx)
            roff_mean=np.fromfile("./CaMa_in/"+runname+"/mean_month/mean_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            roff_total=np.fromfile("./CaMa_in/"+runname+"/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            for ens in np.arange(1,pm.ens_mem()+1):
                ens_char="C%03d"%(ens)
                #print pm.distopen(),std[ens-1]
                ofile=ifile*distopen + roff_mean*std[ens-1]#*10.0
                #ofile=ifile*(distopen + std[ens-1])
                #ofile=ifile*distopen + roff_total*std[ens-1]
                #ofile=ifile*distopen + roff_total*std[ens-1]*0.75
                ofile.astype(np.float32)
                ofile.tofile("./CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")

        # do parallel
        #p=Pool(pm.para_nums())
        #p.map(copy_runoff,inputlist)
        #p.terminate()
    #--------------
    # blind runoff experiment
    if pm.input()=="ELSE_KIM2009": # differnt yeaer
        #distopen=1.0 #0.75 #pm.distopen()
        #diststd=0.5  #pm.diststd()
        distopen=pm.distopen(4)
        diststd=pm.diststd(4)
        true_run=pm.true_run(4) # for true ensemble
        runname=pm.runname(4) # get runoff name  
        # copy for TRUE simulation
        if(len(glob.glob("./CaMa_in/"+runname+"/Roff_TRUE/Roff*"))!=0):
            pass
        # true input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff_TRUE")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_char="T000"
            iname=pm.DA_dir()+"/inp/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".one"
            oname="./CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
            inputlist.append([iname,oname,"1.00"])

        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()

        # calculate mothly mean
        cal_monthly_mean(start_year,end_year,24)

        # calculate monthly total
        cal_monthly_total(start_year,end_year,24)

        # make courrpted runoff
        if(len(glob.glob("./CaMa_in/"+runname+"/Roff_CORR/Roff*"))!=0):
            pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/"+runname+"/Roff_CORR")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            # make random values
            #std=np.fromfile("./CaMa_out/randlist.bin",np.float32)
            std=rd.normal(0,diststd,pm.ens_mem())
            std=std.astype(np.float32)
            #print std
            ifile=np.fromfile("./CaMa_in/"+runname+"/Roff/Roff____"+yyyy+mm+dd+".one",np.float32).reshape(ny,nx)
            roff_mean=np.fromfile("./CaMa_in/"+runname+"/mean_month/mean_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            roff_total=np.fromfile("./CaMa_in/"+runname+"/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(ny,nx)
            for ens in np.arange(1,pm.ens_mem()+1):
                ens_char="C%03d"%(ens)
                #print pm.distopen(),std[ens-1]
                ofile=ifile*distopen + roff_mean*std[ens-1]*10.0
                #ofile=ifile*distopen + roff_total*std[ens-1]
                #ofile=ifile*distopen + roff_total*std[ens-1]*0.75
                ofile.astype(np.float32)
                ofile.tofile("./CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")

        # do parallel
        #p=Pool(pm.para_nums())
        #p.map(copy_runoff,inputlist)
        #p.terminate()
    return 0
###########################
def calc_odb(year,month,date,obj):
    year=int(year)
    month=int(month)
    date=int(date)
    #calc bef date
    if date-1<1:
        if month-1<1:
            bef_m=12
            bef_y=year-1
        else:
            bef_m=month-1
            bef_y=year

        #calc leap year
        leap_y=0 #0:not 1:leap year
        if bef_y%4==0:
            if bef_y%100==0 and bef_y%400!=0:
                leap_y=0
            else:
                leap_y=1

        #calc last day of the month
        if bef_m==1:
            last_d=31
        if bef_m==2:
            if leap_y==1:
                last_d=29
            else:
                last_d=28
        if bef_m==3:
            last_d=31
        if bef_m==4:
            last_d=30
        if bef_m==5:
            last_d=31
        if bef_m==6:
            last_d=30
        if bef_m==7:
            last_d=31
        if bef_m==8:
            last_d=31
        if bef_m==9:
            last_d=30
        if bef_m==10:
            last_d=31
        if bef_m==11:
            last_d=30
        if bef_m==12:
            last_d=31
        bef_d=last_d
    else:
        bef_y=year
        bef_m=month
        bef_d=date-1

    if obj=="year":
        return bef_y
    if obj=="month":
        return bef_m
    if  obj=="date":
        return bef_d
###########################
def sfcelv_mean(ens):
    # calculate yearly mean WSE
    year=pm.spinup_end_year()
    mon =pm.spinup_end_month()
    day =pm.spinup_end_date()
    ens =str(ens) # T000 or C0XX
    #--
    dz=days_year(year)
    nx,ny,gsize=pm.map_dimension()
    #create filename
    yyyy="%04d"%(year)
    mm="%2d"%(mon)
    dd="%2d"%(day)
    # sfcelv=pm.DA_dir()+"/out/"+pm.experiment()+"/CaMa_out/"+yyyy+mm+dd+ens+"/sfcelv"+yyyy+".bin"
    # sfcelv=np.fromfile(sfcelv,np.float32).reshape(dz,ny,nx)
    # sf_mean=np.mean(sfcelv,axis=0)
    # sf_std=np.std(sfcelv,axis=0)
    sfcelv=pm.DA_dir()+"/inp/"
    sfcelv=np.fromfile(sfcelv,np.float32).reshape(ny,nx)
    mkdir(pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/mean_sfcelv")
    fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/mean_sfcelv/meansfcelv"+ens+".bin"
    sf_mean.tofile(fname)
    fname=pm.DA_dir()+"/out/"+pm.experiment()+"/assim_out/mean_sfcelv/stdsfcelv"+ens+".bin"
    sf_std.tofile(fname)
    return 0
##########################
def days_year(year):
    year=int(year)
    days=365
    if calendar.isleap(year):
        days=366
    return days
##########################
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
def calc_mean():
    file_path = "/work/a06/yingying/camada/HydroDA/src/simmean/"+pm.runname(pm.mode())+"/"
    os.makedirs(file_path,exist_ok=True)
    for num in range(1,pm.ens_mem()+1):
        numch='%03d'%num
        filename = read_sim_data(numch)
        var_mean = read_nc(filename,31)
        # # write .bin
        fileens ="C"+numch
        filebin  = file_path+fileens+'.bin'
        with open(filebin,'wb') as file:
            file.write(var_mean)
##########################
def save_statistic():
    # copy mean and std of simulated WSE
    mkdir("./assim_out/mean_sfcelv/")
    if pm.input()=="E2O":
        os.system("cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1980-2014.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin")
        os.system("cp -r "+pm.DA_dir()+"/dat/std_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1980-2014.bin ./assim_out/mean_sfcelv/std_sfcelv.bin")
        print("cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1980-2014.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin")
    if pm.input()=="VIC_BC":
        os.system("cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin")
        os.system("cp -r "+pm.DA_dir()+"/dat/std_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin ./assim_out/mean_sfcelv/std_sfcelv.bin")
        print("cp -r "+pm.DA_dir()+"/dat/mean_sfcelv_"+pm.input()+"_"+pm.mapname()+"_1979-2013.bin ./assim_out/mean_sfcelv/mean_sfcelv.bin")
    return 0
