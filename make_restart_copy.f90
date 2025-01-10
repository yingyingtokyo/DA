program make_restart
implicit none
integer                         :: i,j,ios,n
character(len=128)              :: fname,buf,camadir,expdir,mapname,camain_dir,camaout_dir,exp_name
!-map variables
real                            :: gsize,west, north, east, south ! map boundries
integer                         :: latpx,lonpx,nflp    ! pixel size, calculated
real,allocatable                :: rivsto(:,:),fldsto(:,:),resto(:,:,:),restoreC(:,:,:),rivstoC(:,:),fldstoC(:,:),dt_rivsto(:,:),dt_fldsto(:,:)

real,allocatable                :: elevtn(:,:)

real,allocatable                :: rivlen(:,:),rivwth(:,:),rivsto_max(:,:),rivhgt(:,:)
real,allocatable                :: fldhgt(:,:,:)
integer,allocatable             :: oceanmask(:,:),fldstage(:,:)
real,allocatable                :: grid_area(:,:)
integer,allocatable             :: nextX(:,:),nextY(:,:)
integer,dimension(1500,1320,2)  :: nextXY
real,parameter                  :: g=9.80665,dt=86400.,man=0.03,man2=0.10,pdstmth=10000.
character(len=8)                :: yyyymmdd,onedaybef,onedayaft
real                            :: dhgtpre        !! private
character(len=3)                :: num_name, cal
character(len=10)               :: loop

real,allocatable                :: fldfrac(:,:)

character(len=1)                :: loopchar
integer                         :: k
integer                         :: cor ! for parameter corruption

real,allocatable                :: rivdph_all(:,:,:),rivdph(:,:),flddph(:,:),rivdph_ori(:,:)
real                            :: hgt,pre,Across

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
read(buf,*) num_name

call getarg(8,buf)
read(buf,"(A)") expdir

call getarg(9,buf)
read(buf,"(A)") camaout_dir

call getarg(10,buf)
read(buf,"(A)") camain_dir

call getarg(11,buf)
read(buf,"(A)") exp_name

! call getarg(10,buf) ! not needed
! read(buf,*) cor ! for identifying the corruption 0: none, 1: rivhgt, 2: rivwth 
!                 ! 3:rivman 4: fldhgt, 5: rivhgt, rivwth, rivman, and fldhgt

!==
fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/params.txt"
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
fname=trim(adjustl(expdir))//"/logout/restartError_"//yyyymmdd//trim(adjustl(loopchar))//num_name//".log"
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

allocate(rivwth(lonpx,latpx),rivhgt(lonpx,latpx),fldhgt(lonpx,latpx,nflp),elevtn(lonpx,latpx),nextX(lonpx,latpx),nextY(lonpx,latpx),grid_area(lonpx,latpx),rivlen(lonpx,latpx))
if (trim(adjustl(exp_name))=='wid') then
     fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/rivwth_gwdlr"//num_name//".bin"
else
     fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/rivwth_gwdlr.bin"
end if
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth"
    write(82,*) "no file rivwth at:",fname 
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/rivhgt.bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0)then
    read(34) rivhgt
    ! ocean is -9999
else
    write(*,*) "no file rivhgt"
    write(82,*) "no file rivhgt at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/fldhgt.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx*10,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) fldhgt
    ! ocean is -9999
else
    write(*,*) "no file fldhgt"
    write(82,*) "no file fldhgt at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/elevtn.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) elevtn
    ! ocean is -9999
else
    write(*,*) "no file elevtn"
    write(82,*) "no file elevtn at:",fname
end if
close(34)

fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is -9999
else
    write(*,*) "no file rivlen"
    write(82,*) "no file rivlen at:",fname
end if
close(34)

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/nextxy.bin"
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

! read original water depth
fname="/work/a06/yingying/camada/HydroDA/src/"//trim(adjustl(camaout_dir))//"/"//trim(adjustl(yyyymmdd))//"A"//num_name//"/rivdph"//yyyymmdd(1:4)//".bin"
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0) then
    read(34) rivdph_all
else
    write(*,*) "no rivdph_ori"
end if
close(34)
rivdph_ori=rivdph_all(:,:,24)
deallocate(rivdph_all)

! read grid area
fname=trim(adjustl(camadir))//"/map/"//trim(adjustl(mapname))//"/ctmare.bin" ! after CaMa-Flood v3.9
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
allocate(fldfrac(lonpx,latpx),rivdph(lonpx,latpx),rivsto(lonpx,latpx),flddph(lonpx,latpx),fldsto(lonpx,latpx),resto(lonpx,latpx,2),dt_rivsto(lonpx,latpx),dt_fldsto(lonpx,latpx))

! rivdph after DA
fname="/work/a06/yingying/camada/HydroDA/src/assim_"//trim(adjustl(exp_name))//"/ens_xa/assim/"//trim(adjustl(yyyymmdd))//num_name//"A_xa.bin"
write(*,*) 'test',fname
open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
if(ios==0) then
    read(34) rivdph   !xa
else
    write(*,*) "no rivdph"
end if
close(34)

! calc river storage max
rivsto_max = rivlen*rivwth*rivhgt

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
! 0: xa <= max  1: xa > max
! WSE = elevtn - rivhgt + rivdph
! hence, max = elevtn
rivsto=0.
fldstage=0
fldfrac=0.
do i=1,lonpx
    do j=1,latpx
        if(oceanmask(i,j)<2)then
            write(*,*) "111"
            write(*,*) rivdph(i,j),rivhgt(i,j)
            if(rivdph(i,j)>rivhgt(i,j))then
                fldstage(i,j)=1
            write(*,*) "222"
            end if
        end if
    end do
end do

! calculate rivdph

! calculate flddph
flddph = 0
do i=1,lonpx
  do j=1,latpx
    write(*,*) rivdph(i,j)
    rivsto(i,j) = rivdph(i,j) * rivlen(i,j) * rivwth(i,j)
    write(*,*) "111"
    if(rivdph(i,j)<1e-8)then
       rivdph(i,j) = 0. 
       rivsto(i,j) = 0. 
    end if
    if(fldstage(i,j)>0)then
      flddph(i,j) = rivdph(i,j) - rivhgt(i,j)
      if (flddph(i,j)<1e-2)    flddph(i,j) =0
    end if
  end do
end do
write(*,*) "restart wlv",num_name,rivdph(994,606)
write(*,*) "restart ori",num_name,rivstoC(994,606)
write(*,*) "restart sto",num_name,rivsto(994,606)
write(*,*) "restart max",rivsto_max(994,606)

! calculate flood inundation percentage
do i=1,lonpx
  do j=1,latpx
    if(fldstage(i,j)==1)then
      do n=1,nflp
        if(fldhgt(i,j,n)>flddph(i,j)) exit
        fldstage(i,j) = fldstage(i,j) + 1
      end do
    end if
  end do
end do

! fldstage(i,j) == 0: no flood
! fldstage(i,j) 1~10: flood stage |  1  |  2  |  ...  |  10  |
!                          fldhgt ↑0m  ↑[1] ↑[2] ..↑[9]  ↑[10]
! fldstage(i,j) ==11: over fldmax
! calculate fldsto
fldsto = 0.
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
        if (fldstage(i,j)>nflp+1)then
          pre = fldhgt(i,j,nflp)
          fldsto(i,j) = fldsto(i,j) + grid_area(i,j)*(flddph(i,j)-pre)
          fldfrac(i,j) = 1.0
        end if
    end do
end do
write(*,*) "test fldsto",flddph(994,606),rivdph(994,606),rivhgt(994,606),fldsto(994,606)
write(*,*) "fld level",fldhgt(994,606,:)
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
deallocate(restoreC,rivstoC,fldstoC,oceanmask,rivlen,rivwth,rivhgt,fldhgt)
deallocate(elevtn,nextX,nextY,grid_area)
deallocate(rivdph,rivdph_ori,rivsto,fldstage,fldfrac,flddph,fldsto,resto,dt_rivsto,dt_fldsto)
end program make_restart
