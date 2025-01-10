program data_assim
!**************************
! Data Assimilation using LETKF and empircal local patches [Revel et al,. (2019)]
! created by Ikeshima & Menaka
! Menaka@IIS 2020
! ifort data_assim.f90 -o data_assim -qmkl
!**************************
implicit none
integer,parameter               :: latpx=1320,lonpx=1500
character*128                   :: fname,buf,camadir,DAdir,expdir,patchdir,obsdir,cal,dataname,camaout_dir,camain_dir,interval_str
character*10                    :: loop
character*10                    :: yyyymmddhh,nxtyyyymmddhh,line,line1,line2
character*2                     :: dd
character(len=1)                :: loopchar
real                            :: lat,lon
real,allocatable                :: global_xa(:,:),global_xa_all(:,:,:),global_wlv_all(:,:,:),global_wlv(:,:),ens_xa(:,:),ens_xa_all(:,:,:)
integer*4                       :: ind_cent,lon_cent,lat_cent,patch_size,patch_side,i,j,k,countnum,patch_nums,countR,interval
integer,allocatable             :: local_obs(:),iwork(:),ifail(:),H(:,:)
real,allocatable                :: xf_m(:),xf(:,:),xa(:,:)
real*4,allocatable              :: global_sim(:,:,:,:),globalx(:,:,:),globaltrue(:,:),simmean(:,:,:)
integer                         :: i1,i2,ens_num,num,ios,ovs,info,info2,errflg,m
character*3                     :: numch
real,allocatable                :: xt(:),R(:,:),Rdiag(:),HETRHE(:,:)
real,allocatable                :: W(:,:),Pa(:,:),Pasqr(:,:),UNI(:,:),EfW(:,:)
real                            :: coef,VDVTmax
real,allocatable                :: work(:),la(:),U(:,:),Dinv(:,:),VDVT(:,:),Dsqr(:,:),yo(:),Hxf_m(:),Ef(:,:),la_p(:),U_p(:,:)
integer                         :: day
real,allocatable                :: K_(:,:)
!-for observation
integer,parameter               :: ind=29626  !wlv
real*4,dimension(lonpx,latpx)   :: ref
integer*4,dimension(ind)        :: lon_ind,lat_ind
integer                         :: conflag
!-for map
real,dimension(lonpx,latpx)     :: rivwth,rivlen,rivhgt,fldhgt,elevtn,nextdst,parm_infl,elvmean
integer,dimension(lonpx,latpx)  :: nextX,nextY,ocean,countp,targetp
integer,dimension(lonpx,latpx,2):: nextXY,countfile
real,dimension(lonpx,latpx)     :: errobs,errrand,errstd
real                            :: errfix
!-for inflation
real                            :: rho,rho_fixed ! covariance inflation parameter, 1.01
real,allocatable                :: dep(:),HEf(:,:),HEfR(:,:),HPH(:,:)
real,dimension(4)               :: parm
real                            :: sigma_o,gain
real,parameter                  :: rho_min=1.0d0
real                            :: thresold,sigma_b

real,allocatable                :: global_null(:,:)
real,allocatable                :: Wvec(:),lag(:),local_lag(:),wgt(:),local_wgt(:)
real                            :: wt,lag_dist
integer,allocatable             :: local_ocean(:)
real,allocatable                :: local_river(:)
integer*4                       :: i_m,j_m
integer,allocatable             :: xlist(:),ylist(:)
integer*4                       :: target_pixel,fn
character*8                     :: llon,llat

write(*,*) "--------------da---------------------"

call getarg(1,buf)
read(buf,*) yyyymmddhh

call getarg(2,buf)
read(buf,*) patch_size ! radius

call getarg(3,buf)
read(buf,*) ens_num ! number of ensemble

call getarg(4,buf)
read(buf,*) nxtyyyymmddhh

call getarg(5,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(6,buf)
read(buf,*) thresold

call getarg(7,buf)
read(buf,"(A)") expdir

call getarg(8,buf)
read(buf,"(A)") DAdir

call getarg(9,buf)
read(buf,"(A)") patchdir

call getarg(10,buf)
read(buf,"(A)") obsdir

call getarg(11,buf)
read(buf,*) rho_fixed

call getarg(12,buf)
read(buf,*) sigma_b

call getarg(13,buf)
read(buf,*) conflag

call getarg(14,buf)
read(buf,*) cal

call getarg(15,buf)
read(buf,*) dataname

call getarg(16,buf)
read(buf,*) camaout_dir

call getarg(17,buf)
read(buf,*) camain_dir

call getarg(18,buf)
read(buf,*) loop
! loop: "A" online DA
! loop: "C" offline DA

call getarg(19,buf)
read(buf,*) interval

call getarg(20,buf)
read(buf,*) interval_str

if(trim(adjustl(loop))=="open") loopchar="C"
if(trim(adjustl(loop))=="assim") loopchar="A"

write(*,*) "--------------",yyyymmddhh,'---------------------'
20 format(i4.4,2x,i4.4,2x,f8.4,2x,f8.4,2x,f8.4)
21 format(i4.4,2x,i4.4,2x,f12.7,2x,f12.7,2x,f12.7)
22 format(a4,2x,a4,2x,a8,2x,a8,2x,a8)
23 format(i4.4,2x,i4.4,2x,f10.7)

patch_side=patch_size*2+1
patch_nums=patch_side**2

! rho covariance inflation parameter
rho=1.0d0
!---
fname=trim(adjustl(expdir))//"/logout/errrand_"//yyyymmddhh//".log"
open(36,file=fname,status='replace')
errrand=-1
write(36,*) errrand
close(36)

fname=trim(adjustl(expdir))//"logout/assimLog_"//yyyymmddhh//".log"
open(78,file=fname,status='replace')

fname=trim(adjustl(expdir))//"logout/KLog_"//yyyymmddhh//".log"
open(84,file=fname,status='replace')

fname=trim(adjustl(expdir))//"logout/testLog"//yyyymmddhh//".log"
open(72,file=fname,status='replace')
write(72,22)"lon","lat","true","forcast","assim"

fname=trim(adjustl(expdir))//"logout/pixelLog_"//yyyymmddhh//".log"
open(79,file=fname,status='replace')
write(79,*) "lat","lon","valid pixels in emperical patch"

fname=trim(adjustl(expdir))//"logout/inflation_"//yyyymmddhh//".log"
open(73,file=fname,status='replace')

fname=trim(adjustl(expdir))//"logout/ensembles_"//yyyymmddhh//".log"
open(74,file=fname,status='replace')

fname=trim(adjustl(expdir))//"logout/error_"//yyyymmddhh//".log"
open(82,file=fname,status='replace')

! read river width
fname=trim(adjustl(camadir))//"/map/tej_01min/rivwth_gwdlr.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth"
end if
close(34)

! read river channel depth
if (trim(cal)=="yes") then
    fname=trim(adjustl(camadir))//"/map/tej_01min/rivhgt_Xudong.bin"
elseif (trim(cal)=="corrupt") then
    fname=trim(adjustl(camadir))//"/map/tej_01min/rivhgt_corrupt.bin"
else
    fname=trim(adjustl(camadir))//"/map/tej_01min/rivhgt.bin"
end if
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivhgt
    ! ocean is -9999
else
    write(*,*) "no file rivhgt"
end if
close(34)

! read flood plain height
fname=trim(adjustl(camadir))//"/map/tej_01min/fldhgt.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx*10,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) fldhgt
    ! ocean is -9999
else
    write(*,*) "no file fldhgt"
    write(82,*) "no file fldhgt at:",fname
    write(78,*) "no file fldhgt at:", fname
end if
close(34)

! read river elevation
fname=trim(adjustl(camadir))//"/map/tej_01min/elevtn.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) elevtn
    ! ocean is -9999
else
    write(*,*) "no file elevtn"
end if
close(34)

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"/map/tej_01min/nextxy.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) nextXY
else
    write(*,*) "no file nextXY at:",fname
end if
nextX=nextXY(:,:,1)
nextY=nextXY(:,:,2)
close(34)

! make ocean mask from nextx (1 is ocean; 0 is not ocean)
ocean=(nextX<=0) * (-1)

! read river length
fname=trim(adjustl(camadir))//"/map/tej_01min/rivlen.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivlen
    ! ocean is -9999
else
    write(*,*) "no file rivlen",fname
end if
close(34)

! read distance to next grid
fname=trim(adjustl(camadir))//"/map/tej_01min/nxtdst.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) nextdst
else
    write(*,*) "no file nextdst",fname
end if
close(34)

! read observations
! read reference data of observation
fname='/work/a06/yingying/obs/alloc/situ_elv.bin'
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
   read(34) ref
else
    write(*,*) "no obs"
end if
close(34)
! read observation data
allocate(globaltrue(lonpx,latpx))
!fname=trim(adjustl(obsdir))//"/201910"//dd//"24.bin"   daily online DA
fname=trim(adjustl(obsdir))//"/201910"//nxtyyyymmddhh(7:10)//".bin"
write(*,*) "obspath",fname
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
   read(34) globaltrue
else
    write(*,*) "no obs",fname
end if
close(34)
write(*,*) "globaltrue",globaltrue(994,606)
! change globaltrue into waterdepth
!globaltrue=globaltrue+ref-(elevtn-rivhgt)

! inflation parameter
fname=trim(adjustl(expdir))//"/inflation/parm_infl"//yyyymmddhh//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) parm_infl
else
    write(*,*) "no parm_infl"
end if
close(34)

! read WSE from all model
allocate(global_sim(lonpx,latpx,interval,ens_num))
allocate(globalx(lonpx,latpx,ens_num))
allocate(simmean(lonpx,latpx,ens_num))
globalx=0
global_sim=0
simmean=0
do num=1,ens_num
    write(numch,'(i3.3)') num
    fname="/work/a06/yingying/camada/HydroDA/src/"//trim(adjustl(camaout_dir))//"/"//yyyymmddhh//loopchar//trim(adjustl(numch))//"/rivdph"//yyyymmddhh(1:4)//".bin"
    open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
    if(ios==0)then
        read(34) global_sim(:,:,:,num)
    else
        write(*,*) "no x :", fname
        write(82,*) "no x at:", fname
        write(78,*) "no x at:", fname
    end if
    close(34)
    ! read simmean
    fname="/work/a06/yingying/camada/HydroDA/src/simmean/"//trim(adjustl(dataname))//"/C"//numch//".bin"
    open(35,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
    if(ios==0)then
       read(35) simmean(:,:,num)
    else
        write(*,*) "no simmean"
    end if
    close(35)
end do
globalx=global_sim(:,:,interval,:)
globalx = globalx-simmean
write(*,*) '---globalx---'

! change water level to water depth
fname=trim(adjustl(obsdir))//"/"//"elv_mean.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
   read(34) elvmean
else
    write(*,*) "no elv"
end if
close(34)
globaltrue = globaltrue-elvmean

! read countnum
fname=trim(adjustl(patchdir))//"countnum.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) countfile
else
    write(*,*) "no countp"
end if
countp=countfile(:,:,1)   !目标网格有多少个可用上下游观测数据
targetp=countfile(:,:,2)  !如果目标网格有站点观测，该站点观测在哪里
close(34)

! read obs with err
!fname=trim(adjustl(errdir))//"/obs_err/"//yyyymmddhh//".bin"
!open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
!if(ios==0)then
   !read(34,rec=1) errobs
!else
!    write(*,*) "no errobs"
!end if
!write(*,*) errobs(998,548)
!close(34)

! read err_std
fname="/work/a06/yingying/obs/2019h/obs_err.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
   read(34) errstd
else
    write(*,*) "no errstd"
end if
close(34)

! read site info
fname='/work/a06/yingying/Code/mcode/lon_dis.txt'
open(34,file=fname,form='formatted',iostat=ios)
j=0
if(ios==0)then
    do i=1,ind
        read(34,*) line
        read(line,*) lon_ind(j+1)
        j=j+1
    end do
else
    write(*,*) "no file lon_ind:",fname
end if
close(34)

fname='/work/a06/yingying/Code/mcode/lat_dis.txt'
open(34,file=fname,form='formatted',iostat=ios)
j=0
if(ios==0)then
    do i=1,ind
        read(34,*) line
        read(line,*) lat_ind(j+1)
        j=j+1
    end do
else
    write(*,*) "no file lat_ind:",fname
end if

allocate(global_xa(lonpx,latpx),global_null(lonpx,latpx),global_wlv_all(lonpx,latpx,ens_num),global_wlv(lonpx,latpx),global_xa_all(lonpx,latpx,ens_num))
global_xa = 0
global_xa_all = 0
global_wlv = 0
global_wlv_all = 0
global_null = 0.0

write(78,*) "====================================="
write(78,*) "Assimilation of each grid"

fname=trim(adjustl(expdir))//"assim_"//trim(adjustl(dataname))//trim(adjustl(interval_str))//"/ens_xa/mean"//loopchar//"/"//trim(adjustl(yyyymmddhh))//"_xa.bin"
open(89,file=fname,status="replace",access='stream')

fname=trim(adjustl(expdir))//"assim_"//trim(adjustl(dataname))//trim(adjustl(interval_str))//"/xa_log/"//trim(adjustl(yyyymmddhh))//"_xa.log"

do ind_cent = 1,ind
    lon_cent = lon_ind(ind_cent)
    lat_cent = lat_ind(ind_cent)
    !lon_cent = 522
    !lat_cent = 800
        if (ocean(lon_cent,lat_cent)==1) then
            cycle
        end if
        ! remove rivwth <= 0m
        if (rivwth(lon_cent,lat_cent) <=0.0) then
            cycle
        end if
        write(82,*) lon_cent, lat_cent

        ! countnum and target pixel
        ! countnum is the grids in the local moving window
        countnum=countp(lon_cent,lat_cent)
        target_pixel=targetp(lon_cent,lat_cent)
        ! open emperical local patch
        allocate(lag(patch_nums),xlist(countnum),ylist(countnum),wgt(countnum))
        write(llon,'(i4.4)') lon_cent
        write(llat,'(i4.4)') lat_cent
        write(*,*) "***************",lon_cent,lat_cent,"**************"
        fname=trim(adjustl(patchdir))//"patch"//trim(llon)//trim(llat)//".txt"
        open(34,file=fname,status='old',access='sequential',form='formatted',action='read')
        do i=1, countnum
          read(34,*) xlist(i),ylist(i),wgt(i) !wgt:gauss_weight
        end do
        close(34)

        if (countnum==0) then
            errflg=1
            stop
            write(82,*) lat_cent,lon_cent,"error",errflg
            deallocate(lag,xlist,ylist,wgt)
            write(*,*) "warning",errflg
            goto 9999
        else
            errflg=0
        end if

        allocate(xt(countnum),local_lag(countnum),local_wgt(countnum),local_obs(countnum))
        local_wgt=wgt(1:countnum)
        ! creating local ocean and river pixels
        xt=0
        local_obs=0
        ovs=0

        do i=1,countnum
            i_m=xlist(i)
            j_m=ylist(i)
            xt(i)=globaltrue(i_m,j_m) 
            if ((ocean(i_m,j_m)==0) .and. (globaltrue(i_m,j_m)>-30)) then
                local_obs(i)=1
                ovs=ovs+1
            end if
            !if (ocean(i_m,j_m) ==0) .and. (rivth(i_m,j_m)>50) then
        !    local_ocean(i)=ocean(i_m,j_m)
        !    local_river(i)=rivwth(i_m,j_m)
        end do

        if (ovs==0) then
            errflg=1
            write(82,*) lat_cent,lon_cent,"error",errflg
            deallocate(lag,xlist,ylist,wgt,xt,local_lag,local_wgt,local_obs)
            write(*,*) "warning",errflg
            goto 9999
        else
            ! observation available, make k ==========
            allocate(H(ovs,countnum))
            allocate(R(ovs,ovs),Rdiag(ovs))
            ! make xf =====================================
            allocate(xf(countnum,ens_num))
            R=0.
            countR=1
            H=0
            j=1
            xf=0
            errfix=0.01
            do i=1,countnum ! col number
                i_m=xlist(i)
                j_m=ylist(i)
                xf(i,:)=globalx(i_m,j_m,:)
                !errfix = 0.1*errstd(i_m,j_m)
                if(local_obs(i)==1)then
                    H(j,i)=1
                    wt=local_wgt(i)
                    j=j+1
                    ! errfix:观测值的标准差,即为高斯分布的标准差
                    Rdiag(countR)=errfix/wt
                    countR=countR+1
                end if
            end do
            do i=1,ovs
                R(i,i)=(Rdiag(i))**(-1.)
            end do
        end if

        if (target_pixel==0) then
            xf(countnum,:)=globalx(lon_cent,lat_cent,:)
        end if


        ! make xf_m ====================================
        !------------------------------------------------------------
        deallocate(lag,xlist,ylist,wgt)
        ! ========= reach here only when there is observation =======
        !------------------------------------------------------------
        write(78,*) "================================================"
        write(78,*) "******************",lon_cent,lat_cent,"*******************"
        write(78,*) "size",countnum
        write(78,*) "local obs",local_obs

        ! inflation parameter
        if (rho_fixed==-1.0) then
            rho=parm_infl(lon_cent,lat_cent)
        else
            rho=rho_fixed
        endif
        if (rho<rho_min) rho=rho_min

        ! make Ef ======================================
        allocate(Ef(countnum,ens_num),xf_m(countnum))   !有错
        Ef=0
        xf_m=0

        do i=1,countnum
            xf_m(i)=sum(xf(i,:))/(1.0e-20+real(ens_num))
        end do

        do k=1,ens_num
            Ef(:,k)=(xf(:,k)-xf_m(:))
        end do

        write(79,*) "================================================"
        write(79,*) lon_cent,lat_cent !lat,lon

        ! make W ====================================
        allocate(W(ens_num,ens_num),VDVT(ens_num,ens_num))
        allocate(Pa(ens_num,ens_num),Pasqr(ens_num,ens_num),UNI(ens_num,ens_num))
        allocate(la_p(ens_num),U_p(ens_num,ens_num),HETRHE(ens_num,ens_num))
        W=0
        VDVT=0
        Pa=0
        Pasqr=0
        UNI=0
        la_p=0
        U_p=0

        UNI=RESHAPE([(1,(0,i=1,ens_num),j=1,ens_num-1),1],[ens_num,ens_num])
        HETRHE=matmul(matmul(TRANSPOSE(matmul(H,Ef)),R),matmul(H,Ef)) 
        !VDVTmax=maxval(abs(HETRHE))

        allocate(work(1000),iwork(1000),ifail(1000))
        work=0
        iwork=0
        ifail=0
        info=0

        ! consider influation
        VDVT = real(ens_num-1.)*UNI/rho+HETRHE 

        call ssyevx("V","A","U",ens_num,VDVT,ens_num,1e-5,1e5,1,2,1.2e-38*2.,m,la_p,U_p,ens_num,work,1000,iwork,ifail,info)

        ! m:特征值个数
        if (m<ens_num)then
            errflg=2
            write(*,*) countnum,"Ef",Ef
            write(*,*) "Rdiag",Rdiag
            write(36,*) lat_cent,lon_cent,"error",errflg,"info",info
            write(*,*) "warning",errflg
            goto 9999 
        end if
        allocate(U(m,m),la(m))
        U=0
        la=0
        do i=1,m
            do j=1,m
                U(i,j)=U_p(i,j)
            end do
            la(i)=la_p(i)
        end do
        allocate(Dinv(m,m),Dsqr(m,m))
        Dinv=0
        Dsqr=0

        ! calc Dinv,Dsqr
        ! Dsqr: D-1
        if(info==0)then
            do i=1,m
                Dinv(i,i)=(la(i)+1.0e-20)**(-1)
            end do
            Dsqr=Dinv
            info2=0
            call spotrf("U",m,Dsqr,m,info2)
            if(info2/=0) then !error
                errflg=3
                write(*,*) "warning",errflg
                write(*,*) "====== ERROR cannot unpack Dsqr======",errflg
                write(82,*) lat_cent,lon_cent,"error",errflg
                goto 9999
            end if
            Pa   =matmul(matmul(U,Dinv),transpose(U))
            Pasqr=matmul(matmul(U,Dsqr),transpose(U))

            allocate(yo(ovs),Wvec(ens_num))
            yo=matmul(H,xt)

!Wvec: VD-1VT(HEf)T(R/w)-1(yo-Hxf) -> K(yo-Hxf)
            Hxf_m=matmul(H,xf_m)
            Wvec = matmul(matmul(matmul(Pa,TRANSPOSE(matmul(H,Ef))),R),yo-Hxf_m)
            do i = 1,ens_num
                W(:,i) = Wvec + sqrt(ens_num-1.)*Pasqr(:,i)
            end do

            deallocate(Wvec)
        else
            write(*,*) "NG INFO"
            W=0
            errflg=4
            write(*,*) 'warning',errflg
            write(82,*) lat_cent,lon_cent,"error",errflg,"info:",info
            goto 9999
        end if

        ! make xa ====================================
        allocate(xa(countnum,ens_num))
        xa=0
        allocate(EfW(countnum,ens_num))
        EfW=0

        EfW=matmul(Ef,W)

        do i=1,ens_num
            xa(:,i)=EfW(:,i)+xf_m(:)
        end do

        ! check center pixel ====================================
        if ((target_pixel==0) .or. (globaltrue(lon_cent,lat_cent)<-30)) then
            write(*,*) "no obs here"
            write(*,*) "forcast:",sum(globalx(lon_cent,lat_cent,:))/(ens_num+1.0e-20)
            write(*,*) "assimil:",sum(xa(countnum,:))/(ens_num+1.0e-20)
            write(*,*) "DA wlv:",sum(xa(countnum,:))/(ens_num+1.0e-20)+sum(simmean(lon_cent,lat_cent,:))/(ens_num+1.0e-20)
        else
            write(*,*) "number:",countnum
            write(*,*) "true   :",xt(target_pixel)
            write(*,*) "forcast:",sum(xf(target_pixel,:))/(ens_num+1.0e-20)
            write(*,*) "assimil:",sum(xa(target_pixel,:))/(ens_num+1.0e-20)
            write(*,*) "DA wlv:",sum(xa(target_pixel,:))/(ens_num+1.0e-20)+sum(simmean(lon_cent,lat_cent,:))/(ens_num+1.0e-20)
        end if

! xt: observation
! xf: ensemble forecast
! xa: assimilation result from each model
        write(78,*) "true   :",xt(target_pixel)!+sum(simmean(lon_cent,lat_cent,:))/(ens_num+1.0e-20))
        write(78,*) "forcast:",sum(xf(target_pixel,:))/(ens_num+1.0e-20)
        write(78,*) "assimil:",sum(xa(target_pixel,:))/(ens_num+1.0e-20)
        ! check K_ value (should be between 0-1) =======================================
        allocate(K_(countnum,ovs))
! K_: K kelman filter matrix
        K_ = matmul(Ef,matmul(matmul(Pa,TRANSPOSE(Ef)),R))
        write(78,*) "K:",sum(K_)
        write(72,*) lon_cent,lat_cent,xt(target_pixel),sum(xa(target_pixel,:))/(ens_num+1.0e-20)
        write(74,*) "+++++++++++++++++++++++++++++++++++++"
        write(74,*) lon_cent,lat_cent
        write(74,*) "true   :", xt(target_pixel)
        do i=1, ens_num
             write(74,*) "ensemble",i,"forcast:",xf(target_pixel,i),"assimil:",xa(target_pixel,i)
        end do
        !=====================================
        !-- cal inflation parameter estimation
        parm=0.0d0
        allocate(dep(ovs),HEf(ovs,ens_num),HEfR(ens_num,ovs),HPH(ovs,ovs))
        dep=yo-Hxf_m
        HEf=matmul(H,Ef)
        HEfR=matmul(TRANSPOSE(HEf),R)
        HPH=matmul(HEf,TRANSPOSE(HEf))
        do i=1,ovs
          parm(1)=parm(1)+dep(i)*dep(i)*Rdiag(i)**(-1.)
        end do
        !do num=1,ens_num
        do i=1,ovs
          parm(2)=parm(2)+HPH(i,i)*Rdiag(i)**(-1.)
        end do
        parm(2)=parm(2)/real(ens_num-1)
        parm(3)=sum(local_wgt)
        parm(4)=(parm(1)-parm(3))/(parm(2)+1.0e-20) - rho
        sigma_o=2.0d0/(parm(3)+1.0e-20) * ((rho*parm(2) + parm(3))/(parm(2)+1.0e-20))**2
        gain=sigma_b**2 / (sigma_o + sigma_b**2)
        write(73,*) "+++++++++++++++++++++++++++++++++++"
        write(73,*) gain, parm(4),parm(1),parm(2),parm(3)
        rho=rho+ gain* parm(4)
        if (rho<rho_min) rho=rho_min
        write(73,*) lon_cent,lat_cent,"rho",rho,"rho_min",rho_min
        !---
        parm_infl(lon_cent,lat_cent)=rho


! global_xa: store assimilation result from multi models
        if ((target_pixel==0) .or. (globaltrue(lon_cent,lat_cent)<-30)) then !no obs
            global_xa(lon_cent,lat_cent)=sum(xa(countnum,:))/(ens_num+1.0e-20)
            global_wlv(lon_cent,lat_cent)=sum(globalx(lon_cent,lat_cent,:))/(ens_num+1.0e-20)
            do num=1,ens_num !update wlv
            global_wlv_all(lon_cent,lat_cent,num)=globalx(lon_cent,lat_cent,num)+simmean(lon_cent,lat_cent,num)
            global_xa_all(lon_cent,lat_cent,num) = xa(countnum,num)+simmean(lon_cent,lat_cent,num)
            end do
        else
            global_xa(lon_cent,lat_cent) = sum(xa(target_pixel,:))/(ens_num+1.0e-20)
            global_wlv(lon_cent,lat_cent)= sum(xf(target_pixel,:))/(ens_num+1.0e-20)
            do num=1,ens_num !update wlv
            global_wlv_all(lon_cent,lat_cent,num)=globalx(lon_cent,lat_cent,num)+simmean(lon_cent,lat_cent,num)
            global_xa_all(lon_cent,lat_cent,num) = xa(target_pixel,num)+simmean(lon_cent,lat_cent,num)
            end do
        end if


!有值的地方存入global_null矩阵
        global_null(lon_cent,lat_cent) = 1.0
        if (sum(K_) > real(ens_num)) then 
            global_null(lon_cent,lat_cent) = 0.0
        end if

        if (sum(K_) < 0.0) then 
            global_null(lon_cent,lat_cent) = 0.0
        end if
        9999 continue

        ! clean memory ====================================
        if(errflg==0)then
            deallocate(Ef,R,Rdiag,W,VDVT,la,U,Dinv,Dsqr,Pa,Pasqr,UNI,work,EfW,H,iwork,ifail,U_p,la_p,xa,yo,dep,K_,HEf,HEfR,HPH,HETRHE,xf_m)
        elseif(errflg==2)then
            deallocate(Ef,R,Rdiag,W,VDVT,Pa,Pasqr,UNI,work,H,iwork,ifail,U_p,la_p,HETRHE,xf_m)
        elseif(errflg==3)then
            deallocate(Ef,R,Rdiag,W,VDVT,la,U,Dinv,Dsqr,Pa,Pasqr,UNI,work,H,iwork,ifail,U_p,la_p,HETRHE,xf_m)
        elseif(errflg==4)then
            deallocate(Ef,R,Rdiag,W,VDVT,la,U,Dinv,Dsqr,Pa,Pasqr,UNI,work,H,iwork,ifail,U_p,la_p,HETRHE,xf_m)
        end if
        if (errflg/=1) then
            deallocate(local_obs,xf,xt,local_lag,local_wgt)
        end if
end do

allocate(ens_xa(lonpx,latpx))
ens_xa = global_xa*global_null+global_wlv*(1-global_null)

! store all the output
allocate(ens_xa_all(lonpx,latpx,ens_num))
do num=1,ens_num
    ens_xa_all(:,:,num) = global_xa_all(:,:,num)*(global_null) + global_wlv_all(:,:,num)*(1-global_null)
end do
write(*,*) "global_xa_all",global_xa_all(1006,550:600,1)

! output ens_xa (assimilated + each ensemble result) NEW v.1.1.0
write(89) ens_xa
do num=1,ens_num
    write(numch,'(i3.3)') num
    fname=trim(adjustl(expdir))//"/assim_"//trim(adjustl(dataname))//trim(adjustl(interval_str))//"/ens_xa/assim/"//trim(adjustl(yyyymmddhh))//trim(adjustl(numch))//loopchar//"_xa.bin"
    open(35,file=fname,form="unformatted",access="stream",status="replace",iostat=ios)
    if(ios==0)then
        write(35) ens_xa_all(:,:,num)
    else
        write(*,*) "not created", fname
    end if
    close(35)
end do
write(*,*) "test,",ens_xa_all(1006,550:600,1)
write(*,*) "simmean",simmean(1006,550:600,1),global_wlv_all(1006,550:600,1),global_null(1006,550:600)

! save inflation parameter 
fname=trim(adjustl(expdir))//"inflation/parm_infl"//nxtyyyymmddhh//".bin"
open(36,file=fname,form="unformatted",access="direct",recl=4*lonpx*latpx,status="replace",iostat=ios)
if(ios==0)then
    write(36,rec=1) parm_infl
else
    write(*,*) "no parm_infl"
end if
close(36)

close(78)
close(89)
close(79)
close(84)
close(72)
close(73)
close(74)
close(82)
deallocate(global_xa,global_xa_all,global_wlv,global_wlv_all,globalx,global_sim,simmean,globaltrue,ens_xa,ens_xa_all,global_null)

end program data_assim
!*****************************************************************
subroutine lag_distance(i,j,x,y,nx,ny,nextX,nextY,nextdst,lag_dist)
implicit none 
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: nextdst
!--
real                                :: lag_dist
integer                             :: ix,iy,iix,iiy,tx,ty,ud
real                                :: length,rl
if (i==x .and. j==y) then
  ud=0
else
  ud=-1
end if
if (ud==-1) then
  tx=x
  ty=y
  ix=i
  iy=j
  length=0.0
  lag_dist=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy 
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix==-9 .or. iy==-9) then
      ud=+1
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=+1
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl!/2.0
  end do
end if
if (ud==+1) then
  tx=i
  ty=j
  ix=x
  iy=y
  length=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=-9999
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl
  end do
end if
if (ud==-9999) then
  lag_dist=-9999
elseif (ud==0) then
  lag_dist=0.0
else
  lag_dist=length
end if
return
end subroutine lag_distance
!**************************************************
function roundx(ix, nx)
implicit none
!-- for input -----------
integer                     ix, nx
!-- for output ----------
integer                     roundx
!------------------------
if (ix .ge. 1) then
  roundx = ix - int((ix -1)/nx)*nx
else
  roundx = nx - abs(mod(ix,nx))
end if 
return
end function roundx
!*****************************************************************
subroutine ixy2iixy(ix,iy, nx, ny, iix, iiy)
implicit none
!- for input -----------------
integer                   ix, iy, nx, ny
!- for output ----------------
integer                   iix, iiy,roundx
!-----------------------------
if (iy .lt. 1) then
  iiy = 2 - iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else if (iy .gt. ny) then
  iiy = 2*ny -iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else
  iiy = iy
  iix = roundx(ix, nx)
end if
return
end subroutine ixy2iixy
!*****************************************************************        
