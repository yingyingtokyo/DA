program make_restart
implicit none
integer                         :: i,j,ios,n,k
character(len=128)              :: fname,buf,camadir,expdir,mapname,camain_dir,camaout_dir,exp_name
!-map variables
real(kind=4)                    :: gsize,west, north, east, south,DFRCINC ! map boundries
integer                         :: latpx,lonpx,nflp    ! pixel size, calculated
integer,parameter               :: NLFP=10
real(kind=4),allocatable        :: rivsto(:,:),fldsto(:,:),resto(:,:,:),restoreC(:,:,:),rivstoC(:,:),fldstoC(:,:),dt_rivsto(:,:),dt_fldsto(:,:),fldfrac(:,:)
real(kind=4),allocatable        :: elevtn(:,:),rivlen(:,:),rivwth(:,:),rivsto_max(:,:),rivhgt(:,:)
real(kind=4),allocatable        :: fldhgt(:,:,:),fldstomax(:,:,:)
integer,allocatable             :: oceanmask(:,:),fldstage(:,:)
real(kind=4),allocatable        :: grid_area(:,:)
real(kind=4),allocatable        :: nextdst(:,:)
integer,allocatable             :: nextX(:,:),nextY(:,:)
integer,dimension(1500,1320,2)  :: nextXY
real(kind=4),parameter          :: g=9.80665,dt=86400.,man=0.03,man2=0.10,pdstmth=10000.
character(len=8)                :: yyyymmdd,onedaybef,onedayaft
real(kind=4)                    :: dhgtpre,DSTONOW,DSTOPRE,DWTHPRE,DWTHINC        !! for calculate flood storage 
character(len=3)                :: num_name, cal
character(len=10)               :: loop
character(len=1)                :: loopchar
integer                         :: ens_num,k
integer                         :: cor ! for parameter corruption
real(kind=4),allocatable        :: xa(:,:),rivdph_all(:,:,:),rivdph(:,:),flddph(:,:),rivdphC(:,:),rivdph_ori(:,:)
real(kind=4)                    :: hgt,pre,Across

call getarg(1,buf)
read(buf,*) yyyymmdd

call getarg(2,buf)
read(buf,*) onedaybef

call getarg(3,buf)
read(buf,*) onedayaft

call getarg(4,buf)
read(buf,*) loop

call getarg(5,buf)
read(buf,"(A)") camadir

call getarg(6,buf)
read(buf,"(A)") mapname

call getarg(7,buf)
read(buf,*) ens_num ! number of ensemble members

call getarg(8,buf)
read(buf,*) num_name

call getarg(9,buf)
read(buf,"(A)") expdir

call getarg(10,buf)
read(buf,"(A)") camaout_dir

call getarg(11,buf)
read(buf,"(A)") camain_dir

call getarg(12,buf)
read(buf,"(A)") exp_name

! call getarg(10,buf) ! not needed
! read(buf,*) cor ! for identifying the corruption 0: none, 1: rivhgt, 2: rivwth 
!                 ! 3:rivman 4: fldhgt, 5: rivhgt, rivwth, rivman, and fldhgt

!==
fname=trim(camadir)//"/map/"//trim(mapname)//"/params.txt"
open(11,file=fname,form='formatted')
read(11,*) lonpx
read(11,*) latpx
read(11,*) nflp
read(11,*) gsize
read(11,*) west
read(11,*) east
read(11,*) south
read(11,*) north
close(11)
!-------

if(trim(adjustl(loop))=="open") loopchar="C"
if(trim(adjustl(loop))=="assim") loopchar="A"
allocate(rivdph_ori(lonpx,latpx))
allocate(rivdph_all(lonpx,latpx,24))

! program for making restart file
fname=trim(adjustl(expdir))//"/logout/restartError_"//yyyymmdd//loopchar//num_name//".log"
open(82,file=fname,status='replace')
write(82,*) "make_restart.f90 Errors"

allocate(restoreC(lonpx,latpx,2))
allocate(rivstoC(lonpx,latpx))
allocate(fldstoC(lonpx,latpx))
if (loopchar=="C") then
    fname=trim(adjustl(expdir))//"/"//trim(adjustl(camain_dir))//"/restart/open/restart"//onedayaft//"C"//num_name//".bin"
else if (loopchar=="A") then
    fname=trim(adjustl(expdir))//"/"//trim(adjustl(camain_dir))//"/restart/assim/restart"//onedayaft//"A"//num_name//".bin"
end if
! read original restart file
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) restoreC
end if
rivstoC=restoreC(:,:,1)
fldstoC=restoreC(:,:,2)
close(34)

! read assimilated xa
allocate(xa(lonpx,latpx))
fname=trim(adjustl(expdir))//"/assim_"//exp_name//"/ens_xa/"//trim(adjustl(loop))//"/"//trim(adjustl(yyyymmdd))//trim(adjustl(num_name))//"A_xa.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) xa
else
    write(*,*) "no file at xa"
    write(82,*) "no file xa",fname
end if
close(34)

! read many parameters
! read CaMa-Flood parametes
allocate(rivlen(lonpx,latpx),rivwth(lonpx,latpx),rivhgt(lonpx,latpx),fldhgt(lonpx,latpx,nflp),fldstomax(lonpx,latpx,nflp))
allocate(elevtn(lonpx,latpx),nextX(lonpx,latpx),nextY(lonpx,latpx),nextdst(lonpx,latpx),grid_area(lonpx,latpx))

fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivlen.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivlen
    ! ocean is 1.39e4
else
    write(*,*) "no file rivlen"
    write(82,*) "no file rivlen at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivwth_gwdlr.bin"
! if (cor==2 .OR. cor==5) then
!     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivwth_corrupt.bin"
! end if
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth"
    write(82,*) "no file rivwth at:",fname 
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivhgt.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivhgt
    ! ocean is -9999
else
    write(*,*) "no file rivhgt"
    write(82,*) "no file rivhgt at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/fldhgt.bin"
! if (cor==4 .OR. cor==5) then
!     fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/fldhgt_corrupt.bin"
! end if
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx*10,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) fldhgt
    ! ocean is -9999
else
    write(*,*) "no file fldhgt"
    write(82,*) "no file fldhgt at:",fname
end if
close(34)

!    fname=trim(adjustl(camadir))//"map/"//trim(mapname)//"/elevtn.bin"
!    write(numch,'(i3.3)') num
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/elevtn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) elevtn
    ! ocean is -9999
else
    write(*,*) "no file elevtn"
    write(82,*) "no file elevtn at:",fname
end if
close(34)

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nextxy.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) nextXY
else
    write(*,*) "no file nextXY"
    write(82,*) "no file nextXY at:",fname
end if
nextX=nextXY(:,:,1)
nextY=nextXY(:,:,2)
close(34)

! read distance to next grid
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
    ! 海は
else
    write(*,*) "no file nextdst"
    write(82,*) "no file nextdst at:",fname
end if
close(34)

! read original water depth
fname="/work/a06/yingying/camada/HydroDA/src/"//trim(adjustl(camaout_dir))//"/"//yyyymmdd//"A"//num_name//"/rivdph"//yyyymmdd(1:4)//".bin"
write(*,*) 'test',fname
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0) then
    read(34) rivdph_all
else
    write(*,*) "no rivdph_ori"
end if
close(34)
rivdph_ori=rivdph_all(:,:,24)

! read grid area
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/ctmare.bin" ! after CaMa-Flood v3.9
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) grid_area
    ! 海は
else
    write(*,*) "no file grid_area"
    write(82,*) "no file grarea at:",fname
end if
close(34)

! =======================================================
! allocate
allocate(rivsto_max(lonpx,latpx),oceanmask(lonpx,latpx),fldstage(lonpx,latpx))
allocate(fldfrac(lonpx,latpx),rivdph(lonpx,latpx),rivdphC(lonpx,latpx),rivsto(lonpx,latpx),flddph(lonpx,latpx),fldsto(lonpx,latpx),resto(lonpx,latpx,2),dt_rivsto(lonpx,latpx),dt_fldsto(lonpx,latpx))

! make ocean mask
! 0:inland 1:mouth 2:ocean
do i=1,lonpx
    do j=1,latpx
        if(nextX(i,j)==-9999)then
            oceanmask(i,j)=2
        else if(nextX(i,j)==-9)then ! ocean
            oceanmask(i,j)=1
        else if(nextX(i,j)==-10)then ! inland termination
            oceanmask(i,j)=1
        else
            oceanmask(i,j)=0
        end if
    end do
end do

! =======================================================
! NEW judgement of flood stage from WSE
! 0: xa <= max(within river)   1: xa > max(out of river) 
rivdph = xa !+ rivhgt - elevtn  ! calculate river depth
fldstage=0
fldfrac=0
do i=1,lonpx
    do j=1,latpx
        if(oceanmask(i,j)<2)then
            if(rivdph(i,j)>rivhgt(i,j))then
                fldstage(i,j)=1
            end if
        end if
    end do
end do
write(*,*) "restart wlvori",num_name,rivdph_ori(1006,550:600)
! calc river storage max
rivsto_max = rivlen*rivwth*rivhgt

! [1] Assume all waters in river channel
rivsto = rivdph * rivlen * rivwth
fldsto = 0.
flddph = 0.
dhgtpre = 0. ! height of floodplain
fldstomax(:,:,:) = 0.
DFRCINC = 1./NLFP

! fldstage(i,j) == 0: no flood
! fldstage(i,j) 1~10: flood stage |  1  |  2  |  ...  |  10  |
!                          fldhgt ↑0m  ↑[1] ↑[2] ..↑[9]  ↑[10]
! fldstage(i,j) ==11: over fldmax

do i=1,lonpx
  do j=1,latpx
    ! calculate the river storage in different cases
    ! nearly no river water
    if(rivdph(i,j)<1e-8)then
       rivdph(i,j) = rivdph_ori(i,j)
       rivsto(i,j) = rivstoC(i,j)
    end if

    ! water out of river
    if(fldstage(i,j)==1)THEN
       flddph(i,j) == rivdph(i,j)-rivhgt(i,j)
       ! calculation start from only river storage 
       DSTOALL = rivsto(i,j) + fldsto(i,j) 
       DSTOPRE = rivsto_max(i,j)
       DWTHPRE = rivwth(i,j)
       ! DFRCINC: floodplain fraction increment (1/NLFP),NLFP: floodplain layer, when increase for a whole floodplain layer
       DWTHINC = grid_area(i,j)*rivlen**(-1)*DFRCINC 
       ! floodplain storage very small
       if (flddph(i,j)<1e-2)THEN
           flddph(i,j) = 0
           rivsto  = rivsto_max(i,j) 
           fldstp  = 0.
       else if
           ! calculate floodplain storage
           k = 1  !floodplain layer
           DO WHILE (fldght(i,j,k)<rivdph(i,j) .AND. k<=NLFP)
               ! formular from cmf_ctrl_maps_mod.F90
               DSTONOW = rivlen(i,j) * (rivwth(i,j) + DWTHINC*(dble(k)-0.5)) * (fldhgt(i,j,k)-dhgtpre)
               DSTOPRE = DSTOPRE + DSTONOW
               DWTHPRE = DWTHPRE + DWTHINC 
               dhgtpre = fldhgt(i,j,k)
               fldstage(i,j) = fldstage(i,j) + 1
               k = k+1
       hgt = hgt
       if (k<=NLFP)THEN
           ! [1] water within floodplain
           DSTONOW = rivlen(i,j) * (rivwth(i,j) + DWTHINC*(dble(k-1)-0.5)+DWTHINC*(1.)*((flddph(i,j)-dhgtpre)/(fldhgt(i,j,k)-dhgtpre+1e-20))) * (rivdph(i,j)-dhgtpre)
       else if
           ! [2] water out of floodplain
           DSTONOW = rivlen(i,j) * DWTHPRE * (flddph(i,j)-dhgtpre) 
       end if
       ! current whole water storage
       DSTOPRE = DSTOPRE + DSTONOW
       fldsto(i,j)  = DSTOPRE - rivsto(i,j)
    end if
  end do
end do
write(*,*) "restart wlv",num_name,rivdph(1006,550:600)
write(*,*) "restart ori",num_name,rivstoC(1006,550:600)
write(*,*) "restart sto",num_name,rivsto(1006,550:600)
write(*,*) "restart max",rivsto_max(1006,550:600)

! calculate flood inundation percentage
do i=1,lonpx
  do j=1,latpx
    if(fldstage(i,j)==1)then
      do n=1,10
        if(fldhgt(i,j,n)>flddph(i,j)) exit
        fldstage(i,j) = fldstage(i,j) + 1
      end do
    end if
    ! fldstage(i,j) == 0: no flood
    ! fldstage(i,j) 1~10: flood stage |  1  |  2  |  ...  |  10  |
    !                          fldhgt ↑0m  ↑[1] ↑[2] ..↑[9]  ↑[10]
    ! fldstage(i,j) ==11: over fldmax
  end do
end do

! calculate fldsto
fldsto = 0
do i=1,lonpx
    do j=1,latpx
        ! calculate flood storage until fldstage(i,j)-1
        if(fldstage(i,j)>1)then
          dhgtpre = 0
          do n=1,fldstage(i,j)-1
            fldsto(i,j) = fldsto(i,j) + grid_area(i,j)*0.1*(real(n)-0.5)*(fldhgt(i,j,n)-dhgtpre)
            dhgtpre = fldhgt(i,j,n)
            fldfrac(i,j) = fldfrac(i,j) + 0.1
          end do
        end if
        ! calculate flood storage at fldstage(i,j)
        if(fldstage(i,j)>0 .and. fldstage(i,j)<11)then
          if(fldstage(i,j)==1)then
            dhgtpre = 0
          end if
          if(fldstage(i,j)>1)then
            dhgtpre = fldhgt(i,j,fldstage(i,j)-1)
          end if
          n = fldstage(i,j)
          k = fldstage(i,j)-1
          hgt = fldhgt(i,j,n)
          pre = dhgtpre
          !Across = grid_area(i,j)*(k*(hgt-flddph(i,j))+(K+1)*(flddph(i,j)-pre))/(hgt-pre+1e-20)*0.1
          Across = 0.1*grid_area(i,j)*(k + ((flddph(i,j)-pre)/(hgt-pre+1e-20))*1.0)   
          fldfrac(i,j) = fldfrac(i,j) + 0.1/(hgt-pre+1e-20)*(0*(hgt-flddph(i,j))+1*(flddph(i,j)-pre))

          fldsto(i,j) = fldsto(i,j) + ((Across+0.1*grid_area(i,j)*k)*0.5)*(flddph(i,j)-pre)

        end if
        if (fldstage(i,j)==11)then
          pre = fldhgt(i,j,10)
          fldsto(i,j) = fldsto(i,j) + grid_area(i,j)*(flddph(i,j)-pre)
          fldfrac(i,j) = 1.0
        end if
    end do
end do

! ===================================================
! arrangement for ocean mask
do i=1,lonpx
    do j=1,latpx
        if(oceanmask(i,j)==2)then
            rivsto(i,j)=1e20
            fldsto(i,j)=1e20
        end if
    end do
end do

! =================================================
! save and store output & restart file
! make restart file
! rivsto,fldsto
resto(:,:,1) = rivsto
resto(:,:,2) = fldsto 
fname=trim(adjustl(expdir))//"/"//trim(adjustl(camain_dir))//"/restart/"//trim(adjustl(loop))//"/restart"//onedayaft//loopchar//num_name//".bin"
open(35,file=fname,form="unformatted",access="stream",status="replace",iostat=ios)
if(ios==0)then
    write(35) resto
end if
close(35)
write(82,*) "done restart file at:",fname
print*, "---> restart ", fname
close(82)
! deallocate
deallocate(restoreC,rivstoC,fldstoC,xa,rivlen,rivwth,rivhgt,fldhgt,fldstomax)
deallocate(elevtn,nextX,nextY,nextdst,grid_area)
deallocate(rivdph,rivdph_ori,rivdphC,rivsto,flddph,fldsto,resto,dt_rivsto,dt_fldsto)
end program make_restart
