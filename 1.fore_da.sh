#!/bin/sh
#==========================================================
# CaMa-Flood sample go script (1) global 15min simulation
# -- Multi 1-year simulations (2000 spinup -> 2000 -> 2001)
# -- Daily runoff forcing (plain binary) at 1deg resolution
#
# (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   You may not use this file except in compliance with the License.
#   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is 
#  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and limitations under the License.
#==========================================================

#*** PBS setting when needed
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -m ea
#PBS -V
#========
#edited by Menaka@IIS for SWOTDA 2019/11/05
#================================================
# input settings
orgDIR=`pwd`

in_year=2019
in_month=10
in_date=13

ens_num=20
ens_num=$(printf "%03d" $ens_num)

CAMADIR="/work/a06/yingying/HydroDA/CaMa-Flood-DA_Japan20240520/CaMa_v411"

looptype="assim"

cpunums=20

runname='tej'

EXP_DIR='/work/a06/yingying/camada/HydroDA/src'

mapname='tej_01min'

cal='no'

DA_dir='/work/a06/yingying/camada/HydroDA'

# pwd
# echo ${CAMADIR}
echo $EXP_DIR
# years,months,dates in arranged digit
ar_year=$in_year
ar_month=$in_month
ar_date=$in_date
#----
if [ $looptype = "true" ];then
	ar_ens="T"$ens_num
elif [ $looptype = "open" ];then
	ar_ens="C"$ens_num
else
	ar_ens="A"$ens_num
fi
#================================================
# (0) Basic Setting (for workstation)

#*** 0a. Set CaMa-Flood base directory
BASE=$CAMADIR                         #   CaMa-Flood directory
PWDD=${BASE}/etc/reservoir_operation

OUTBASE="$EXP_DIR/forecast/$runname/$looptype"    		   #   base for output => CaMa_out
mkdir -p $OUTBASE
INBASE="/work/a06/yingying/camada/HydroDA/src/CaMa_in"                #   base for input => CaMa_in

# echo $BASE

#*** 0b. Set dynamic library if needed
export IFORTLIB="/opt/intel/lib:/opt/intel/mkl/lib"
export HDF5LIB="/opt/local/hdf5-1.10.5/lib"
export DYLD_LIBRARY_PATH="${IFORTLIB}:${DYLD_LIBRARY_PATH}"

#*** 0c. OpenMP thread number
export OMP_NUM_THREADS=$cpunums       # OpenMP cpu num

#================================================
# (1) Experiment setting
# -- some non-default options can be modified in NAMELIST section

#============================
#*** 1a. Experiment directory setting
EXP=$ar_year$ar_month$ar_date$ar_ens        # experiment name (output directory name)
RDIR=${OUTBASE}/${EXP}                      # directory to run CaMa-Flood
EXE="MAIN_cmf"                              # Execute file name
PROG=${BASE}/src/${EXE}                     # location of Fortran main program
NMLIST="./input_cmf.nam"                    # standard namelist
LOGOUT="./log_CaMa.txt"                     # standard log output


#============================
#*** 1b. Model physics option
DT=3600                                    # base DT (modified in physics loop by LADPSTP)
LADPSTP=".TRUE."                            # .TRUE. for adaptive time step

LFPLAIN=".TRUE."                            # .TRUE. to activate floodplain storage
LKINE=".FALSE."                             # .TRUE. to use kinematic wave equation
LFLDOUT=".TRUE."                            # .TRUE. to activate floodplain discharge
LPTHOUT=".FALSE."                            # .TRUE. to activate bifurcation flow, mainly for delta simulation
LPTHOUT=".TRUE."                       # .TRUE. to activate bifurcation flow, mainly for delta simulation
LDAMOUT=".FALSE."                           # .TRUE. to activate reservoir operation (under development)

#============================
#*** 1c. simulation time
# simulation for spinup is 1 day
SPINUP=1                                     # [0]: zero-storage start, [1]: from restart file
NSP=0                                        # spinup repeat time


#============================
#*** 1d. spinup setting

#* input restart file
LRESTART=".TRUE." # see (3) set each year   # TRUE. to use restart initial condition
if [ $looptype = "true" ];then
	CRESTSTO=$INBASE"/restart/true/restart"$ar_year$ar_month$ar_date"T000.bin" #restart file name
elif [ $looptype = "open" ];then
	#CRESTSTO=$INBASE"/restart/open/restart"$ar_year$ar_month$ar_date"C"$ens_num".bin" #restart file name
	CRESTSTO="/work/a06/yingying/camada/HydroDA/src/forecast/tej/mean/"$looptype$ar_date".bin" #restart file name
else
	CRESTSTO="/work/a06/yingying/camada/HydroDA/src/forecast/tej/mean/"$looptype$ar_date".bin" #restart file name
fi
echo $CRESTSTO
LSTOONLY=".TRUE."                          # .TRUE. for storage only restart (for assimilation)


#* output restart file
CRESTDIR="./"                               # output restart file directory
CVNREST="restart"                           # output restart file prefix
LRESTCDF=".FALSE."                          # .TRUE. to use netCDF restart file
LRESTDBL=".FALSE."                           # .TRUE. for binary restart double precision, This is new !!
IFRQ_RST="0"                                # output restat frequency.
                                            # [0]: only at last time, [1,2,3,...,24] hourly restart, [30]: monthly restart


#============================
#*** 1e. forcing setting
IFRQ_INP="1"                               # input forcing frequency: [1,2,3,...,24] hour
DROFUNIT="86400000"        # [mm/day->m/s]  # runoff unit conversion
if [ $runname = "E2O" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ECMWF000" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ECMWF050" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ERA20CM" ];then
     DROFUNIT="1000"       # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ELSE_KIM2009" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "VIC_BC" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "isimip3a" ];then
     DROFUNIT="1000"       # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ERA5" ];then
     DROFUNIT="86400000"      # [m/day->m/s]   # runoff unit conversion
elif [ $runname = "tej" ];then
     DROFUNIT="86400000"      # [kg/m2/s->m/s]   # runoff unit conversion
fi

#----- for plain binary runoff forcing
LINPCDF=".FALSE."                           # true for netCDF runoff
LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix

CROFDIR="${INBASE}/inp/$runname"    #   runoff directory
CROFPRE='Runoff_'
CROFSUF="$ens_num.bin" 
#CROFSUF="$ens_num.one" 
 
###** sub-surface runoff scheme (not available with plain binary runoff)
LROSPLIT=".FALSE."                          # .TRUE. for sub-surface runoff


#============================
#*** 1f. river map & topography
FMAP="${BASE}/map/${mapname}"                  # map directory
CDIMINFO="${FMAP}/diminfo_test-1deg.txt"      # dimention information file
CINPMAT="${FMAP}/inpmat_test-1deg.bin"        # runoff input matrix for interporlation
if [ $runname = "E2O" ];then
    CDIMINFO="${FMAP}/diminfo-15min.txt" # dimention information file
    CINPMAT="${FMAP}/inpmat-15min.bin"   # runoff input matrix for interporlation
elif [ ${runname} = "ERA5" ] ; then
#     CDIMINFO="${FMAP}/diminfo-06min.txt"  # dimention information file
#     CINPMAT="${FMAP}/inpmat-06min.bin"      # runoff input matrix for interporlation
     CDIMINFO="${FMAP}/diminfo_test-tej.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat_test-tej.bin"      # runoff input matrix for interporlation
elif [ ${runname} = "tej" ] ; then
     CDIMINFO="${FMAP}/diminfo_test-tej.txt"  # dimention information file
     CINPMAT=${FMAP}/inpmat_test-tej.bin      # runoff input matrix for interporlation
fi
#----- for plain binary map input
#** basic topography
LMAPCDF=".FALSE."                           # .TRUE. for netCDF map
CNEXTXY="${FMAP}/nextxy.bin"                # downstream xy (river network map)
CGRAREA="${FMAP}/ctmare.bin"                # unit-catchment area   [m2]
CELEVTN="${FMAP}/elevtn.bin"                # channel top elevation [m]
CNXTDST="${FMAP}/nxtdst.bin"                # downstream distance   [m]
CRIVLEN="${FMAP}/rivlen.bin"                # channel length        [m]
CFLDHGT="${FMAP}/fldhgt.bin"                # floodplain elevation profile (height above 'elevtn') [m]
# if [ $corrupt = 4 ] || [ $corrupt = 5 ];then
#      CFLDHGT="${FMAP}/fldhgt_corrupt.bin"   # floodplain elevation profile [m] (corrupted)
# fi
#** channel parameter
###CRIVWTH=${FMAP}/rivwth.bin"              # channel width [m] (empirical power-low)
CRIVWTH="${FMAP}/rivwth_gwdlr.bin"          # channel width [m] (GWD-LR + filled with empirical)
CRIVHGT="${FMAP}/rivhgt.bin"                # channel depth [m] (empirical power-low)
if [ $cal = "yes" ];then
  CRIVHGT="${FMAP}/rivhgt_Xudong.bin"       # channel depth [m] (Xudong et al,. 2022)
elif [ $cal = "corrupt" ];then
  CRIVHGT="${FMAP}/rivhgt_corrupt.bin"      # channel depth [m] (Corrupted rivhgt simple)
elif [ $cal = "no" ];then
  CRIVHGT="${FMAP}/rivhgt.bin"      # channel depth [m] (Corrupted rivhgt simple)
fi
CRIVMAN="${FMAP}/rivman.bin"                # manning coefficient river (The one in flood plain is a global parameter; set $PMANFLD below.)

#** bifurcation channel info
CPTHOUT="${FMAP}/bifprm.txt"                #   bifurcation channel list

###** groundwater delay (not available in plain binary runoff/map)
LGDWDLY=".FALSE."                           # .TRUE. to actuvate groundwater delay

###** mean sea level
LMEANSL=".FALSE."                           # .TRUE. to use mean sea level data

#============================
#*** 1g. Dynamic Boundary Sea Level (not default)
LSEALEV=".FALSE."                           # .TRUE. to activate dynamic sea level 

#============================
#*** 1h. Output Settings 
LOUTPUT=".TRUE."                            # .TRUE. to use CaMa-Flood standard output
IFRQ_OUT=1                                 # output frequency: [1,2,3,...,24] hour

LOUTCDF=".TRUE."                           # .TRUE. netCDF output, .FALSE. plain binary output
COUTDIR="./"                                # output directory 
CVARSOUT="rivdph,sfcelv,outflw,storge" #"rivout,rivsto,sfcelv,outflw,storge,fldout,fldsto,daminf,damsto,flddph,fldare" # list output variable (comma separated) # Moiz

COUTTAG=""  # see (3) set each year         #   output tag $(COUTDIR)/$(VARNAME)$(OUTTAG).bin

##### Model Parameters ################
PMANRIV="0.03D0"                            # manning coefficient river
PMANFLD="0.10D0"                            # manning coefficient floodplain
PCADP="0.7"                                 # satety coefficient for CFL condition
PDSTMTH="10000.D0"                          # downstream distance at river mouth [m]



#================================================
# (2) Initial setting

#*** 2a. create running dir 
mkdir -p ${RDIR}
cd ${RDIR}

#*** 2b. for new simulation, remove old files in running directory

if [ ${SPINUP} -eq 0 ]; then
  rm -rf ${RDIR}/????-sp*
  rm -rf ${RDIR}/*.bin
  rm -rf ${RDIR}/*.pth
  rm -rf ${RDIR}/*.vec
  rm -rf ${RDIR}/*.nc
  rm -rf ${RDIR}/*.log
  rm -rf ${RDIR}/*.txt
  rm -rf ${RDIR}/restart*
else
  NSP=0  # restart, no spinup
fi

#================================================
# (3) For each simulation year, modify setting
#--  loop 1-day simulation from $YSTART to $YEND
YSTA=$in_year
ISP=1           ## spinup count
IYR=${YSTA}     ## curent year
CYR=`printf %04d ${IYR}`


#*** 3b. update start-end year
SYEAR=$in_year
SMON=$in_month
SDAY=$in_date
SHOUR=0

EYEAR=`python ${DA_dir}/src/calc_end_fore.py $in_year $in_month $in_date "year"`
EMON=`python ${DA_dir}/src/calc_end_fore.py $in_year $in_month $in_date "month"`
EDAY=`python ${DA_dir}/src/calc_end_fore.py $in_year $in_month $in_date "date"`
EHOUR=0

ln -sf $PROG $EXE

#*** 3c. update input / output file data
CSYEAR=`printf %04d ${SYEAR}`
COUTTAG=${CSYEAR}                  # output file tag


#================================================
# (4) Create NAMELIST for simulation year
# it is OK to remove optional variables (set to default in CaMa-Flood)

rm -f ${NMLIST}

#*** 0. config
cat >> ${NMLIST} << EOF
&NRUNVER
LADPSTP  = ${LADPSTP}                  ! true: use adaptive time step
LFPLAIN  = ${LFPLAIN}                  ! true: consider floodplain (false: only river channel)
LKINE    = ${LKINE}                    ! true: use kinematic wave
LFLDOUT  = ${LFLDOUT}                  ! true: floodplain flow (high-water channel flow) active
LPTHOUT  = ${LPTHOUT}                  ! true: activate bifurcation scheme
LDAMOUT  = ${LDAMOUT}                  ! true: activate dam operation (under development)
LROSPLIT = ${LROSPLIT}                 ! true: input if surface (Qs) and sub-surface (Qsb) runoff
LGDWDLY  = ${LGDWDLY}                  ! true: Activate ground water reservoir and delay
LSLPMIX  = .FALSE.                     ! true: activate mixed kinematic and local inertia based on slope
LRESTART = ${LRESTART}                 ! true: initial condition from restart file
LSTOONLY = ${LSTOONLY}                 ! true: storage only restart (mainly for data assimilation)
LOUTPUT  = ${LOUTPUT}                  ! true: use standard output (to file)
LGRIDMAP = .TRUE.                      ! true: for standard XY gridded 2D map
LLEAPYR  = .TRUE.                      ! true: neglect leap year (Feb29 skipped)
LMAPEND  = .FALSE.                     ! true: for map data endian conversion
LBITSAFE = .FALSE.                     ! true: for Bit Identical simulation (avoid OSM ATOMIC)
/
&NDIMTIME
CDIMINFO = "${CDIMINFO}"               ! text file for dimention information
DT       = ${DT}                       ! time step length (sec)
IFRQ_INP = ${IFRQ_INP}                 ! input forcing update frequency (hour)
/
&NPARAM
PMANRIV  = ${PMANRIV}                  ! manning coefficient river
PMANFLD  = ${PMANFLD}                  ! manning coefficient floodplain
PGRV     = 9.8D0                       ! gravity accerelation
PDSTMTH  = ${PDSTMTH}                  ! downstream distance at river mouth [m]
PCADP    = ${PCADP}                    ! CFL coefficient
PMINSLP  = 1.D-5                       ! minimum slope (kinematic wave)
IMIS     = -9999                       ! missing value for integer
RMIS     = 1.E20                       ! missing value for real*4
DMIS     = 1.E20                       ! missing value for real*8
CSUFBIN  = '.bin'                      ! file suffix for plain binary 2D map
CSUFVEC  = '.vec'                      ! file suffix for plain binary 1D vector
CSUFPTH  = '.pth'                      ! file suffix for plain binary bifurcation channel
CSUFCDF  = '.nc'                       ! file suffix for netCDF
/
EOF

#*** 1. time
cat >> ${NMLIST} << EOF
&NSIMTIME
SYEAR   = ${SYEAR}                     ! start year
SMON    = ${SMON}                      !  month
SDAY    = ${SDAY}                      !  day
SHOUR   = ${SHOUR}                     !  hour
EYEAR   = ${EYEAR}                     ! end year
EMON    = ${EMON}                      !  month
EDAY    = ${EDAY}                      !  day
EHOUR   = ${EHOUR}                     !  hour
/
EOF

#*** 2. map
cat >> ${NMLIST} << EOF
&NMAP
LMAPCDF    = ${LMAPCDF}                ! * true for netCDF map input
CNEXTXY    = "${CNEXTXY}"              ! river network nextxy
CGRAREA    = "${CGRAREA}"              ! catchment area
CELEVTN    = "${CELEVTN}"              ! bank top elevation
CNXTDST    = "${CNXTDST}"              ! distance to next outlet
CRIVLEN    = "${CRIVLEN}"              ! river channel length
CFLDHGT    = "${CFLDHGT}"              ! floodplain elevation profile
CRIVWTH    = "${CRIVWTH}"              ! channel width
CRIVHGT    = "${CRIVHGT}"              ! channel depth
CRIVMAN    = "${CRIVMAN}"              ! river manning coefficient
CPTHOUT    = "${CPTHOUT}"              ! bifurcation channel table
CGDWDLY    = "${CGDWDLY}"              ! Groundwater Delay Parameter
CRIVCLINC  = "${CRIVCLINC}"            ! * river map netcdf
CRIVPARNC  = "${CRIVPARNC}"            ! * river parameter netcdf (width, height, manning, ground water delay)
/
EOF

#*** 3. restart
cat >> ${NMLIST} << EOF
&NRESTART
CRESTSTO = "${CRESTSTO}"               ! restart file
CRESTDIR = "${CRESTDIR}"               ! restart directory
CVNREST  = "${CVNREST}"                ! restart variable name
LRESTCDF = ${LRESTCDF}                 ! * true for netCDF restart file
LRESTDBL = ${LRESTDBL}                 ! * true for double precision binary restart  !! This is new
IFRQ_RST = ${IFRQ_RST}                 ! restart write frequency (1-24: hour, 0:end of run)
/
EOF

#*** 4. forcing
if [ ${LINPCDF} = ".FALSE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
LINPEND  = .FALSE.                     ! true for runoff endian conversion
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFDIR  = "${CROFDIR}"                ! runoff             input directory
CROFPRE  = "${CROFPRE}"                ! runoff             input prefix
CROFSUF  = "${CROFSUF}"                ! runoff             input suffix
CSUBDIR  = "${CSUBDIR}"                ! sub-surface runoff input directory
CSUBPRE  = "${CSUBPRE}"                ! sub-surface runoff input prefix
CSUBSUF  = "${CSUBSUF}"                ! sub-surface runoff input suffix
/
EOF

elif [ ${LINPCDF} = ".TRUE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
LINPEND  = .FALSE.                     ! true for runoff endian conversion
LITRPCDF = ${LINTERPCDF}               ! * true for netCDF input matrix
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFCDF  = "${CROFCDF}"                ! * netCDF input runoff file name
CVNROF   = "${CVNROF}"                 ! * netCDF input runoff variable name
CVNSUB   = "${CVNSUB}"                 ! * netCDF input runoff variable name
SYEARIN  = ${SYEARIN}                  ! * netCDF input start year
SMONIN   = ${SMONIN}                   ! * netCDF input start year
SDAYIN   = ${SDAYIN}                   ! * netCDF input start year
SHOURIN  = ${SHOURIN}                  ! * netCDF input start year
/
EOF

fi # (if LINPCDF)

#*** 5. outputs
cat >> ${NMLIST} << EOF
&NOUTPUT
COUTDIR  = "${COUTDIR}"                ! OUTPUT DIRECTORY
CVARSOUT = "${CVARSOUT}"               ! Comma-separated list of output variables to save 
COUTTAG  = "${COUTTAG}"                ! Output Tag Name for each experiment
LOUTVEC  = .FALSE                      ! TRUE FOR VECTORIAL OUTPUT, FALSE FOR NX,NY OUTPUT
LOUTCDF  = ${LOUTCDF}                  ! * true for netcdf outptu false for binary
NDLEVEL  = 0                           ! * NETCDF DEFLATION LEVEL 
IFRQ_OUT = ${IFRQ_OUT}                 ! output data write frequency (hour)
/
EOF

#================================================
# (5) Execute main program

echo "start: ${SYEAR}" `date`  >> log.txt
time ./${EXE}                  >> log.txt 
echo "end:   ${SYEAR}" `date`  >> log.txt

mv ${LOGOUT} log_CaMa-${CYR}.txt

# move restart file
RYR=`printf %04d ${EYEAR}`
RMO=`printf %02d ${EMON}`
RDA=`printf %02d ${EDAY}`
cp restart${RYR}${RMO}${RDA}00.bin restart${RYR}${RMO}${RDA}.bin

###################
# going back to original directory
cd $orgDIR

exit 0
