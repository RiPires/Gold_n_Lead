C********************************************************************
C                            ULHISTOS                               *
C                                                                   *
C          Luis Peralta and Ana Farinha @ 2008 (version 1)          *
C                                                                   *
C          Luis Peralta   Version 6.9 - Copyright (c) 2021          * 
C                                                                   *
C            Universidade de Lisboa, Faculdade de Ciencias          *
C Laboratorio de Instrumentacao e Fisica Experimental de Particulas *
C                                                                   *
C           1-d, 2-d and 3-d histograms for Ulysses                 *
C********************************************************************


C    The ULHISTOS package provides histogramming capability to the 
c    Ulysses geometry package. Although a component of Ulysses, the
c    ULHISTOS package can be separatly used. The packages comes with
c    a manual, where all routines are explained.   


C    Copyright (C) 2008  Luis Peralta, Ana Farinha 
C    Copyright (C) 2012  Luis Peralta
C    Copyright (C) 2016  Luis Peralta
C    Copyright (C) 2017  Luis Peralta
C    Copyright (C) 2020  Luis Peralta
C
C    This program is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program.  If not, see <http://www.gnu.org/licenses/>.
C

C    This version is NOT fully compatible with previous versions


C===============================================================
c---------------- MODULES --------------------------------------
C===============================================================

         MODULE ulhglobal

!         Common variables for ULHISTOS
      
       character(len=5), parameter :: version_ulhistos='6.9.2'

       integer(kind=4) :: nh1d    ! Max number of 1-d histos
       integer(kind=4) :: nh2d    ! Max number of 2-d histos
       integer(kind=4) :: nh3d    ! Max number of 3-d histos

       integer(kind=4) :: nmax1d ! Max number of channels per 1-d histos
       integer(kind=4) :: nmax2d ! Max number of channels is nmax2d**2 for 2-d histos
       integer(kind=4) :: nmax3d ! Max number of channels is nmax3d**3 for 3-d histos

       integer(kind=4) :: nmax2dx,nmax2dy          ! Max number of channels is nmax2d**2 for 2-d histos
       integer(kind=4) :: nmax3dx,nmax3dy,nmax3dz  ! Max number of channels is nmax3d**3 for 3-d histos

!      ulhistos 1D histograms
       real(kind=8), dimension(:,:), allocatable :: histo1d,histo1dv   ! up to (nh1d,nmax1d)
       real(kind=8), dimension(:,:), allocatable :: histo1dc            ! up to (nh1d,nmax1d)
       real(kind=8), dimension(:), allocatable ::xmin1d,xmax1d,binx1d  ! up to (nh1d)
       real(kind=8), dimension(:), allocatable :: underflow1d          ! up to (nh1d)
       real(kind=8), dimension(:), allocatable :: overflow1d           ! up to (nh1d)
       integer(kind=4), dimension(:), allocatable :: nx1d,istatus1d    ! up to (nh1d)
       integer(kind=4), dimension(:), allocatable :: nentries1d        ! up to (nh1d)

!      ulhistos 2D histograms 
       real(kind=8), dimension(:,:,:), allocatable :: histo2d ! up to (nh2d,nmax2d,nmax2d)
       real(kind=8), dimension(:,:,:), allocatable :: histo2dv ! up to (nh2d,nmax2d,nmax2d)
       real(kind=8), dimension(:,:,:), allocatable :: histo2dc ! up to (nh2d,nmax2d,nmax2d)
       real(kind=8), dimension(:), allocatable :: xmin2d,xmax2d,binx2d ! up to (nh2d)
       real(kind=8), dimension(:), allocatable :: ymin2d,ymax2d,biny2d ! up to (nh2d)
       real(kind=8), dimension(:), allocatable :: underflow2d          ! up to (nh2d)
       real(kind=8), dimension(:), allocatable :: overflow2d           ! up to (nh2d)
       integer(kind=4), dimension(:), allocatable :: nx2d,ny2d         ! up to (nh2d)
       integer(kind=4), dimension(:), allocatable :: istatus2d         ! up to (nh2d)
       integer(kind=4), dimension(:), allocatable :: nentries2d        ! up to (nh2d)

!      ulhistos 3D histograms 
       real(kind=8), dimension(:,:,:,:), allocatable :: histo3d  ! up to (nh3d,nmax3d,nmax3d,nmax3d)
       real(kind=8), dimension(:,:,:,:), allocatable :: histo3dv ! up to (nh3d,nmax3d,nmax3d,nmax3d)
       real(kind=8), dimension(:,:,:,:), allocatable :: histo3dc ! up to (nh3d,nmax3d,nmax3d,nmax3d)
       real(kind=8), dimension(:), allocatable :: xmin3d,xmax3d,binx3d ! up to (nh3d)
       real(kind=8), dimension(:), allocatable :: ymin3d,ymax3d,biny3d ! up to (nh3d)
       real(kind=8), dimension(:), allocatable :: zmin3d,zmax3d,binz3d ! up to (nh3d)
       real(kind=8), dimension(:), allocatable :: underflow3d          ! up to (nh3d)
       real(kind=8), dimension(:), allocatable :: overflow3d           ! up to (nh3d)
       integer(kind=4), dimension(:), allocatable :: nx3d,ny3d,nz3d    ! up to (nh3d)
       integer(kind=4), dimension(:), allocatable :: istatus3d         ! up to (nh3d)
       integer(kind=4), dimension(:), allocatable :: nentries3d        ! up to (nh3d)

C      number of events (equal to number of entries in frequency histograms)
       integer(kind=4), dimension(:), allocatable :: nev1d     ! up to (nh1d)
       integer(kind=4), dimension(:), allocatable :: nev2d      ! up to (nh2d)
       integer(kind=4), dimension(:), allocatable :: nev3d       ! up to (nh3d)

!      partial counters
       real(kind=8), dimension(:,:), allocatable :: histo1dp    ! up to (nh1d,nmax1d)
       real(kind=8), dimension(:,:,:), allocatable :: histo2dp  ! up to (nh2d,nmax2d,nmax2d)
       real(kind=8), dimension(:,:,:,:), allocatable :: histo3dp ! up to (nh3d,nmax3d,nmax3d,nmax3d)


!      iddim : dimension of histo id
!      idtyp : type of histo id
       integer(kind=4), dimension(:), allocatable :: ids1d !(nh1d)
       integer(kind=4), dimension(:), allocatable :: ids2d !(nh2d)
       integer(kind=4), dimension(:), allocatable :: ids3d !(nh3d)

       integer(kind=4), parameter :: nmaxhistos=99999
       integer(kind=4), parameter :: nuserhistos=9999
       integer(kind=4), dimension(nmaxhistos) :: ids1d_1,ids2d_1,ids3d_1 ! (nmaxhistos)
       integer(kind=4), dimension(nmaxhistos) :: iddim,idtyp,istatus     ! (nmaxhistos)
       integer(kind=4) :: nids,nids1d,nids2d,nids3d

!      input/output units
       integer(kind=4) :: iue,ius,iuo,iui
!      Common events
       integer(kind=4) :: ntotevents

!      titles
       character(len=64), dimension(:), allocatable :: title1d !(nh1d)
       character(len=64), dimension(:), allocatable :: title2d !(nh2d)
       character(len=64), dimension(:), allocatable :: title3d !(nh3d)

!      2d hit-or-miss random generator
       real*8, dimension(:), allocatable :: histog2cmax ! histogram max. content value
       real*8, dimension(:), allocatable :: histog2xmin,histog2xmax
       real*8, dimension(:), allocatable :: histog2ymin,histog2ymax
       real*8, dimension(:), allocatable :: histog2binx,histog2biny

      END MODULE ulhglobal


      MODULE ULHRANDOM
!     random generator seeds
      integer(kind=4) :: iseed1,iseed2
      END MODULE ULHRANDOM


C================================================
C       Booking routines
C================================================


C************************************************
C       Initialize histogram database
C************************************************
       subroutine ulhinit(nh1,nmx1,nh2,nmx2,nmy2,nh3,nmx3,nmy3,nmz3)

       use ulhglobal
       use ulhrandom
       implicit none
       save
       integer(kind=4) :: nh1,nmx1,nh2,nmx2,nmy2,nh3,nmx3,nmy3,nmz3



!      default dimension values
       nmax1d=1000
       nmax2dx=100
       nmax2dy=100
       nmax3dx=50
       nmax3dy=50
       nmax3dz=50 
       nh1d=100
       nh2d=10 
       nh3d=0 

       write(*,*)' *** ULHISTOS Version: ',version_ulhistos,' ***'

       ISEED1=123456789 ! random generator seeds
       ISEED2=987654321

!      set no. of 1d, 2d and 3d histograms. If nh1,nh2 or nh3 =0 then ulhlimits value assumed
       if(nh1 >= 0) nh1d=nh1
       if(nh2 >= 0) nh2d=nh2
       if(nh3 >= 0) nh3d=nh3

       if(nmx1 > 0) nmax1d=nmx1

       if(nmx2 > 0) nmax2dx=nmx2
       if(nmy2 > 0) nmax2dy=nmy2

       if(nmx3 > 0) nmax3dx=nmx3
       if(nmy3 > 0) nmax3dy=nmy3
       if(nmz3 > 0) nmax3dz=nmz3

!      1D
       if(nh1d > 0) then
       allocate(histo1d(nh1d,nmax1d))  ! allocate 1d histograms memory space
       allocate(histo1dv(nh1d,nmax1d))
       allocate(histo1dc(nh1d,nmax1d))
       allocate(histo1dp(nh1d,nmax1d))
       allocate(nx1d(nh1d))
       allocate(xmin1d(nh1d))
       allocate(xmax1d(nh1d))
       allocate(binx1d(nh1d))
       allocate(underflow1d(nh1d))
       allocate(overflow1d(nh1d))
       allocate(nentries1d(nh1d))
       allocate(ids1d(nh1d))
       allocate(nev1d(nh1d))
       allocate(title1d(nh1d))

           histo1d=0.0d0
           histo1dv=0.0d0
           histo1dc=0.0d0
           histo1dp=0.0d0
         
         nx1d=0
         xmin1d=0.0d0
         xmax1d=0.0d0
         binx1d=0.0d0
         underflow1d=0.0d0
         overflow1d=0.0d0
         nentries1d=0
         ids1d=0
         nev1d=0
         title1d=" "

       nids1d=0
       endif

!      2D
       if(nh2d > 0) then
       allocate(histo2d(nh2d,nmax2dx,nmax2dy))
       allocate(histo2dv(nh2d,nmax2dx,nmax2dy))
       allocate(histo2dc(nh2d,nmax2dx,nmax2dy))
       allocate(histo2dp(nh2d,nmax2dx,nmax2dy))
       allocate(nx2d(nh2d))
       allocate(xmin2d(nh2d))
       allocate(xmax2d(nh2d))
       allocate(binx2d(nh2d))
       allocate(ny2d(nh2d))
       allocate(ymin2d(nh2d))
       allocate(ymax2d(nh2d))
       allocate(biny2d(nh2d))
       allocate(underflow2d(nh2d))
       allocate(overflow2d(nh2d))
       allocate(nentries2d(nh2d))
       allocate(ids2d(nh2d))
       allocate(nev2d(nh2d))
       allocate(title2d(nh2d))
c      2d hit-or-miss generator
       allocate(histog2cmax(nh2d))
       allocate(histog2xmin(nh2d))
       allocate(histog2xmax(nh2d))
       allocate(histog2ymin(nh2d))
       allocate(histog2ymax(nh2d))
       allocate(histog2binx(nh2d))
       allocate(histog2biny(nh2d))
         
             histo2d=0.0d0
             histo2dv=0.0d0
             histo2dc=0.0d0
             histo2dp=0.0d0

         nx2d=0
         xmin2d=0.0d0
         xmax2d=0.0d0
         binx2d=0.0d0
         ny2d=0
         ymin2d=0.0d0
         ymax2d=0.0d0
         biny2d=0.0d0
         underflow2d=0.0d0
         overflow2d=0.0d0
         nentries2d=0
         ids2d=0
         nev2d=0
         title2d=" "

       nids2d=0
       endif

!      3D
       if(nh3d > 0) then

        allocate(histo3d(nh3d,nmax3dx,nmax3dy,nmax3dz))
        allocate(histo3dv(nh3d,nmax3dx,nmax3dy,nmax3dz))
        allocate(histo3dc(nh3d,nmax3dx,nmax3dy,nmax3dz))
        allocate(histo3dp(nh3d,nmax3dx,nmax3dy,nmax3dz))
        allocate(nx3d(nh3d))
        allocate(xmin3d(nh3d))
        allocate(xmax3d(nh3d))
        allocate(binx3d(nh3d))
        allocate(ny3d(nh3d))
        allocate(ymin3d(nh3d))
        allocate(ymax3d(nh3d))
        allocate(biny3d(nh3d))
        allocate(nz3d(nh3d))
        allocate(zmin3d(nh3d))
        allocate(zmax3d(nh3d))
        allocate(binz3d(nh3d))
        allocate(underflow3d(nh3d))
        allocate(overflow3d(nh3d))
        allocate(nentries3d(nh3d))
        allocate(ids3d(nh3d))
        allocate(nev3d(nh3d))
        allocate(title3d(nh3d))

             histo3d=0.0d0
             histo3dv=0.0d0
             histo3dc=0.0d0
             histo3dp=0.0d0

         nx3d=0
         xmin3d=0.0d0
         xmax3d=0.0d0
         binx3d=0.0d0
         ny3d=0
         ymin3d=0.0d0
         ymax3d=0.0d0
         biny3d=0.0d0
         nz3d=0
         zmin3d=0.0d0
         zmax3d=0.0d0
         binz3d=0.0d0

         underflow3d=0.0d0
         overflow3d=0.0d0
         nentries3d=0
         ids3d=0
         nev3d=0
         title3d=" "

       endif

       nids3d=0


        ids1d_1=0
        ids2d_1=0
        ids3d_1=0
        iddim=0
        idtyp=1  ! histograms type default to 1: frequency
        istatus=0
      

       nids=0
       ntotevents=0 ! total number of events

C       set output units

       iue=6 ! output errors     
       ius=6 ! output summaries
       iuo=99 ! standard output for histograms
       iui=81 ! standard input for histograms

       end


C************************************************
C       Reset full database
C************************************************
!      this routine is kept for backward compatibility
       subroutine ulhstart

       use ulhglobal
       use ulhrandom
       implicit none
       save

!      initialize arrays with default
       call ulhinit(-1,0,-1,0,0,-1,0,0,0)

       end


C************************************************
C       Book 1-d  histograms                       
C************************************************
       subroutine ulhb1(id,tit,nx,xmin,xmax)

       use ulhglobal 
       implicit none     
       save
       integer(kind=4) :: id,nx,ip,idime,i
       character(len=64) :: tit
       real(kind=8) ::xmin,xmax

       if( -99999 <= id .and. id < -9999) then 
          id=abs(id)   ! 10000 - 99999 system reserved id
       else if(id <= 0 .or. id > nuserhistos) then
        write(iue,*)'ULHB1: id outside range 1-9999, ID=',id
        return
       endif

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.ne.0)then
         write(iue,*)'ULHB1: histogram already exist, ID=',id
         goto 99
        endif

       if(nx.le.0)then
         write(iue,*)'ULHB1: number channels must be positive, ID=',id
         goto 99
       endif

       if(nx.gt.nmax1d)then
         write(iue,*)'ULHB1: too many channels, ID=',id
         goto 99
       endif       

       call ulhputid(id,ip,1) ! put ID into database
       if(ip.eq.0) goto 99    ! an error occurred -> exit

c      reset all channels
       do i=1,nx
         histo1d(ip,i)=0.
         histo1dv(ip,i)=0.
         histo1dc(ip,i)=0.
         histo1dp(ip,i)=0.
       enddo

       nx1d(ip)=nx                 ! number of channels
       xmin1d(ip)=xmin             ! lower bound
       xmax1d(ip)=xmax             ! upper bound
       binx1d(ip)=(xmax-xmin)/nx   ! bin width
       underflow1d(ip)=0.
       overflow1d(ip)=0.
       nentries1d(ip)=0
       istatus(id)=0
       nev1d(ip)=0

       title1d(ip)=tit

99       end

C************************************************
C        Book 2-d  histograms                       
C************************************************
       subroutine ulhb2(id,tit,nx,xmin,xmax,ny,ymin,ymax)     
       use ulhglobal
       implicit none 
       save
       integer*4 :: id,ip,idime,nx,ny,i,j
       character(len=64) :: tit
       real(kind=8) ::xmin,xmax,ymin,ymax

       if( -99999 <= id .and. id < -9999) then 
          id=abs(id)   ! 10000 - 99999 system reserved id
       else if(id <= 0 .or. id > nuserhistos) then
        write(iue,*)'ULHB2: id outside range 1-9999, ID=',id
        return
       endif

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.ne.0)then
         write(iue,*)'ULHB2: histogram already exist, ID=',id
         goto 99
        endif

       if(nx.le.0 .or. ny.le.0)then
         write(iue,*)'ULHB2: number channels must be positive, ID=',id
         goto 99
       endif       

       if(nx.gt.nmax2dx .or. ny.gt.nmax2dy)then
         write(iue,*)'ULHB2: too many channels, ID=',id
         goto 99
       endif       

       call ulhputid(id,ip,2) ! put ID into database
       if(ip.eq.0) goto 99    ! an error occurred -> exit

c      reset all channels
       do i=1,nx
         do j=1,ny
            histo2d(ip,i,j)=0.0d0
            histo2dv(ip,i,j)=0.0d0
            histo2dc(ip,i,j)=0.0d0
            histo2dp(ip,i,j)=0.0d0
         enddo
       enddo

       nx2d(ip)=nx                 ! number of x channels
       ny2d(ip)=ny                 ! number of y channels
       xmin2d(ip)=xmin             ! lower x bound
       xmax2d(ip)=xmax             ! upper x bound
       ymin2d(ip)=ymin             ! lower y bound
       ymax2d(ip)=ymax             ! upper y bound
       binx2d(ip)=(xmax-xmin)/nx   ! x bin width       
       biny2d(ip)=(ymax-ymin)/ny   ! y bin width
       underflow2d(ip)=0.
       overflow2d(ip)=0.
       nentries2d(ip)=0
       nev2d(ip)=0
       istatus(id)=0
       title2d(ip)=tit

99       end

C************************************************
C        Book 3-d  histograms                       
C************************************************
       subroutine ulhb3(id,tit,nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax)
      
       use ulhglobal
       implicit none
       save
       integer(kind=4) :: id,nx,ny,nz,idime,ip,i,j,k
       character(len=64) :: tit
       real(kind=8) ::xmin,xmax,ymin,ymax,zmin,zmax

       if( -99999 <= id .and. id < -9999) then 
          id=abs(id)   ! 10000 - 99999 system reserved id
       else if(id <= 0 .or. id > nuserhistos) then
        write(iue,*)'ULHB3: id outside range 1-9999, ID=',id
        return
       endif

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.ne.0)then
         write(iue,*)'ULHB3: histogram already exist, ID=',id
         goto 99
        endif

       if(nx.le.0 .or. ny.le.0 .or. nz .le.0)then
         write(iue,*)'ULHB3: number channels must be positive, ID=',id
         goto 99
       endif       

       if(nx.gt.nmax3dx .or. ny.gt.nmax3dy .or. nz.gt.nmax3dz)then
         write(iue,*)'ULHB3: too many channels, ID=',id
         goto 99
       endif       

       call ulhputid(id,ip,3) ! put ID into database
       if(ip.eq.0) goto 99    ! an error occurred -> exit

c      reset all channels
       do i=1,nx
        do j=1,ny
         do k=1,nz
            histo3d(ip,i,j,k)=0.0d0
            histo3dv(ip,i,j,k)=0.0d0
            histo3dc(ip,i,j,k)=0.0d0
            histo3dp(ip,i,j,k)=0.0d0
         enddo
        enddo
       enddo

       nx3d(ip)=nx                 ! number of x channels
       ny3d(ip)=ny                 ! number of y channels
       nz3d(ip)=nz                 ! number of z channels
       xmin3d(ip)=xmin             ! lower x bound
       xmax3d(ip)=xmax             ! upper x bound
       ymin3d(ip)=ymin             ! lower y bound
       ymax3d(ip)=ymax             ! upper y bound
       zmin3d(ip)=zmin             ! lower z bound
       zmax3d(ip)=zmax             ! upper z bound
       binx3d(ip)=(xmax-xmin)/nx   ! x bin width       
       biny3d(ip)=(ymax-ymin)/ny   ! y bin width
       binz3d(ip)=(zmax-zmin)/nz   ! z bin width
       underflow3d(ip)=0.
       overflow3d(ip)=0.
       nentries3d(ip)=0
       nev3d(ip)=0
       istatus(id)=0
       title3d(ip)=tit

99       end

C************************************************
C           Copy histograms
C************************************************
       subroutine ulhcopy(id1,id2)

       use ulhglobal
       implicit none
       save
       integer*4 ::id1,id2,ip1,ip2,idime1,idime2
       integer*4 ::itype,i,j,k
       
       call ulhgetip(id1,ip1,idime1) ! get ID1 from database
        if(idime1.eq.0)then
         write(iue,*)'ULHCOPY: histogram ID1 does not exist, ID=',id1
         return
        endif

        call ulhgetip(id2,ip2,idime2) ! check ID2 in database
        select case(idime2)
         case(0)              ! id2 doesn't exist
           call ulhputid(id2,ip2,idime1) ! put ID2 into database
           if(ip2.eq.0) return     ! an error occurred -> exit
         case default
          write(iue,*)
     &    'Warning ULHCOPY: ID2 exists and will be replaced, ID=',id2
          if(idime2 .ne. idime1) then
           write(iue,*)
     &     'Error ULHCOPY: ID1, ID2 have diferent dimensions',id1,id2
           return
          endif
         end select

       itype=idtyp(id1)         ! histo type =1,2,3,4
       idtyp(id2)=idtyp(id1)         ! histo type

       if(idime1.eq.1)then

c      copy histo

       do i=1,nx1d(ip1)
         histo1d(ip2,i)=histo1d(ip1,i)
         histo1dv(ip2,i)=histo1dv(ip1,i)
         histo1dc(ip2,i)=histo1dc(ip1,i)
         histo1dp(ip2,i)=histo1dp(ip1,i)
       enddo

         nx1d(ip2)=nx1d(ip1)               
         xmin1d(ip2)=xmin1d(ip1)             
         xmax1d(ip2)=xmax1d(ip1)             
         binx1d(ip2)=binx1d(ip1)  
         underflow1d(ip2)=underflow1d(ip1)
         overflow1d(ip2)=overflow1d(ip1)
         nentries1d(ip2)=nentries1d(ip1)
         istatus(id2)=istatus(id1)
         nev1d(ip2)=nev1d(ip1)
         title1d(ip2)=title1d(ip1)
         ids1d(ip2)=id2
 
       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
        do j=1,ny2d(ip1)
         histo2d(ip2,i,j)=histo2d(ip1,i,j)
         histo2dv(ip2,i,j)=histo2dv(ip1,i,j)
         histo2dc(ip2,i,j)=histo2dc(ip1,i,j)
         histo2dp(ip2,i,j)=histo2dp(ip1,i,j)
        enddo
       enddo

         nx2d(ip2)=nx2d(ip1)
         ny2d(ip2)=ny2d(ip1)               
         xmin2d(ip2)=xmin2d(ip1)             
         xmax2d(ip2)=xmax2d(ip1) 
         ymin2d(ip2)=ymin2d(ip1)             
         ymax2d(ip2)=ymax2d(ip1)             
         binx2d(ip2)=binx2d(ip1)  
         biny2d(ip2)=biny2d(ip1)  
         underflow2d(ip2)=underflow2d(ip1)
         overflow2d(ip2)=overflow2d(ip1)
         nentries2d(ip2)=nentries2d(ip1)
         istatus(id2)=istatus(id1)
         nev2d(ip2)=nev2d(ip1)

         title2d(ip2)=title2d(ip1)

         ids2d(ip2)=id2

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
         histo3d(ip2,i,j,k)=histo3d(ip1,i,j,k)
         histo3dv(ip2,i,j,k)=histo3dv(ip1,i,j,k)
         histo3dc(ip2,i,j,k)=histo3dc(ip1,i,j,k)
         histo3dp(ip2,i,j,k)=histo3dp(ip1,i,j,k)
         enddo
        enddo
       enddo

         nx3d(ip2)=nx3d(ip1)
         ny3d(ip2)=ny3d(ip1)
         nz3d(ip2)=nz3d(ip1)    
         xmin3d(ip2)=xmin3d(ip1)             
         xmax3d(ip2)=xmax3d(ip1) 
         ymin3d(ip2)=ymin3d(ip1)             
         ymax3d(ip2)=ymax3d(ip1)
         zmin3d(ip2)=zmin3d(ip1)             
         zmax3d(ip2)=zmax3d(ip1)                 
         binx3d(ip2)=binx3d(ip1)  
         biny3d(ip2)=biny3d(ip1)  
         binz3d(ip2)=binz3d(ip1)  
         underflow3d(ip2)=underflow3d(ip1)
         overflow3d(ip2)=overflow3d(ip1)
         nentries3d(ip2)=nentries3d(ip1)
         istatus(id2)=istatus(id1)
         nev3d(ip2)=nev3d(ip1)
         title3d(ip2)=title3d(ip1)
         ids3d(ip2)=id2

       endif

       end

C************************************************
C           Copy histogram properties
C************************************************
       subroutine ulhcpp(id1,id2)

       use ulhglobal
       save
       integer(kind=4) ::id1,id2,ip1,ip2,idime1,idime2
       integer :: itype
       
       call ulhgetip(id1,ip1,idime1) ! get ID1 from database
        if(idime1.eq.0)then
         write(iue,*)'ULHCPP: histogram ID1 does not exist, ID=',id1
         return
        endif

        call ulhgetip(id2,ip2,idime2) ! check ID2 in database
        select case(idime2)
         case(0)              ! id2 doesn't exist
           call ulhputid(id2,ip2,idime1) ! put ID2 into database
           if(ip2.eq.0) return     ! an error occurred -> exit
         case default
          if(idime2 .ne. idime1) then
           write(iue,*)
     &     'Error ULHCPP: ID1, ID2 have diferent dimensions',id1,id2
           return
          endif
         end select

       itype=idtyp(id1)         ! histo type =1,2,3,4
       idtyp(id2)=idtyp(id1)    ! histo type


       if(idime1.eq.1)then

         nx1d(ip2)=nx1d(ip1)               
         xmin1d(ip2)=xmin1d(ip1)             
         xmax1d(ip2)=xmax1d(ip1)             
         binx1d(ip2)=binx1d(ip1)  
         istatus(id2)=istatus(id1)
         title1d(ip2)=title1d(ip1)
         ids1d(ip2)=id2
         

       elseif(idime1.eq.2)then
         nx2d(ip2)=nx2d(ip1)
         ny2d(ip2)=ny2d(ip1)               
         xmin2d(ip2)=xmin2d(ip1)             
         xmax2d(ip2)=xmax2d(ip1) 
         ymin2d(ip2)=ymin2d(ip1)             
         ymax2d(ip2)=ymax2d(ip1)             
         binx2d(ip2)=binx2d(ip1)  
         biny2d(ip2)=biny2d(ip1)  
         istatus(id2)=istatus(id1)
         title2d(ip2)=title2d(ip1)
         ids2d(ip2)=id2

       elseif(idime1.eq.3)then
         nx3d(ip2)=nx3d(ip1)
         ny3d(ip2)=ny3d(ip1)
         nz3d(ip2)=nz3d(ip1)    
         xmin3d(ip2)=xmin3d(ip1)             
         xmax3d(ip2)=xmax3d(ip1) 
         ymin3d(ip2)=ymin3d(ip1)             
         ymax3d(ip2)=ymax3d(ip1)
         zmin3d(ip2)=zmin3d(ip1)             
         zmax3d(ip2)=zmax3d(ip1)                 
         binx3d(ip2)=binx3d(ip1)  
         biny3d(ip2)=biny3d(ip1)  
         binz3d(ip2)=binz3d(ip1)  
         istatus(id2)=istatus(id1)
         title3d(ip2)=title3d(ip1)
         ids3d(ip2)=id2

       endif

       end

C***********************************************************
C       Fill id2 with no. of events from id1 in each channel
C***********************************************************
       subroutine ulhcpev(id1,id2)

       use ulhglobal
       save
       integer(kind=4) ::id1,id2,ip1,ip2,idime1,idime2
       integer :: itype
       
       call ulhgetip(id1,ip1,idime1) ! get ID1 from database
        if(idime1.eq.0)then
         write(iue,*)'ULHCPEV: histogram ID1 does not exist, ID=',id1
         return
        endif

        call ulhgetip(id2,ip2,idime2) ! check ID2 in database
        select case(idime2)
         case(0)              ! id2 doesn't exist
           call ulhputid(id2,ip2,idime1) ! put ID2 into database
           if(ip2.eq.0) return     ! an error occurred -> exit
         case default
          write(iue,*)
     &    'Warning ULHCPEV: ID2 exists and will be replaced, ID=',id2
          if(idime2 .ne. idime1) then
           write(iue,*)
     &     'Error ULHCPEV: ID1, ID2 have diferent dimensions',id1,id2
           return
          endif
         end select

       itype=idtyp(id1)         ! histo type =1,2,3,4
       idtyp(id2)=idtyp(id1)         ! histo type

       if(idime1.eq.1)then

c      copy histo

       do i=1,nx1d(ip1)
         histo1d(ip2,i)=histo1dc(ip1,i)
         histo1dv(ip2,i)=sqrt(histo1dc(ip1,i))
       enddo

         nx1d(ip2)=nx1d(ip1)               
         xmin1d(ip2)=xmin1d(ip1)             
         xmax1d(ip2)=xmax1d(ip1)             
         binx1d(ip2)=binx1d(ip1)  
         underflow1d(ip2)=underflow1d(ip1)
         overflow1d(ip2)=overflow1d(ip1)
         nentries1d(ip2)=nentries1d(ip1)
         istatus(id2)=istatus(id1)
         nev1d(ip2)=nev1d(ip1)
         title1d(ip2)=title1d(ip1)
         ids1d(ip2)=id2
 
       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
        do j=1,ny2d(ip1)
         histo2d(ip2,i,j)=histo2dc(ip1,i,j)
         histo2dv(ip2,i,j)=sqrt(histo2dc(ip1,i,j))
        enddo
       enddo

         nx2d(ip2)=nx2d(ip1)
         ny2d(ip2)=ny2d(ip1)               
         xmin2d(ip2)=xmin2d(ip1)             
         xmax2d(ip2)=xmax2d(ip1) 
         ymin2d(ip2)=ymin2d(ip1)             
         ymax2d(ip2)=ymax2d(ip1)             
         binx2d(ip2)=binx2d(ip1)  
         biny2d(ip2)=biny2d(ip1)  
         underflow2d(ip2)=underflow2d(ip1)
         overflow2d(ip2)=overflow2d(ip1)
         nentries2d(ip2)=nentries2d(ip1)
         istatus(id2)=istatus(id1)
         nev2d(ip2)=nev2d(ip1)

         title2d(ip2)=title2d(ip1)

         ids2d(ip2)=id2

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
         histo3d(ip2,i,j,k)=histo3dc(ip1,i,j,k)
         histo3dv(ip2,i,j,k)=sqrt(histo3dc(ip1,i,j,k))
         enddo
        enddo
       enddo

         nx3d(ip2)=nx3d(ip1)
         ny3d(ip2)=ny3d(ip1)
         nz3d(ip2)=nz3d(ip1)    
         xmin3d(ip2)=xmin3d(ip1)             
         xmax3d(ip2)=xmax3d(ip1) 
         ymin3d(ip2)=ymin3d(ip1)             
         ymax3d(ip2)=ymax3d(ip1)
         zmin3d(ip2)=zmin3d(ip1)             
         zmax3d(ip2)=zmax3d(ip1)                 
         binx3d(ip2)=binx3d(ip1)  
         biny3d(ip2)=biny3d(ip1)  
         binz3d(ip2)=binz3d(ip1)  
         underflow3d(ip2)=underflow3d(ip1)
         overflow3d(ip2)=overflow3d(ip1)
         nentries3d(ip2)=nentries3d(ip1)
         istatus(id2)=istatus(id1)
         nev3d(ip2)=nev3d(ip1)
         title3d(ip2)=title3d(ip1)
         ids3d(ip2)=id2

       endif

       istatus(id2)=3

       end

C************************************************
C       Book and fill an histogram with external
C       function fun(x)
C************************************************
       subroutine ulhbf1(id,tit,nx,xmin,xmax,fun)
       
       use ulhglobal 
       save      
       integer (kind=4) :: id,nx
       real (kind=8) :: x,xmin,xmax
       character(len=64) :: tit
       real*8 :: fun
       external fun

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.ne.0)then
         write(iue,*)'ULHBF1: histogram already exist, ID=',id
         goto 99
        endif

       if(nx.le.0)then
         write(iue,*)'ULHBF1: number channels must be positive, ID=',id
         goto 99
       endif

       if(nx.gt.nmax1d)then
         write(iue,*)'ULHBF1: too many channels, ID=',id
         goto 99
       endif       

       call ulhputid(id,ip,1) ! put ID into database
       if(ip.eq.0) goto 99    ! an error occurred -> exit

       nx1d(ip)=nx                 ! number of channels
       xmin1d(ip)=xmin             ! lower bound
       xmax1d(ip)=xmax             ! upper bound
       binx1d(ip)=(xmax-xmin)/nx   ! bin width
       underflow1d(ip)=0.
       overflow1d(ip)=0.
       nentries1d(ip)=0
       istatus(id)=4
       nev1d(ip)=0

       title1d(ip)=tit

c       set all channels
       do i=1,nx
         x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! compute channel center 
         histo1d(ip,i)=fun(x) ! compute function at channel center
         histo1dv(ip,i)=0.    ! error array is left empty for future use
         histo1dc(ip,i)=1.    ! set to 1 entry
         histo1dp(ip,i)=0.
       enddo

99       end


C************************************************
C       Reset contents of an histogram
C************************************************
       subroutine ulhreset(id)

       use ulhglobal
       save
       integer(kind=4) :: id
       
       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHRESET: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.eq.1)then

       do i=1,nx1d(ip)
         histo1d(ip,i)  =0.d0
         histo1dv(ip,i) =0.d0  
         histo1dc(ip,i) =0.d0 
         histo1dp(ip,i) =0.d0  
       enddo

         underflow1d(ip)=0.d0
         overflow1d(ip) =0.d0
         nentries1d(ip) =0
         istatus(id)  =0
         nev1d(ip)=0

       elseif(idime.eq.2)then

       do i=1,nx2d(ip)
        do j=1,ny2d(ip)
         histo2d(ip,i,j)=0.d0  
         histo2dv(ip,i,j)=0.d0 
         histo2dc(ip,i,j)=0.d0 
         histo2dp(ip,i,j)=0.d0  
        enddo
       enddo

         underflow2d(ip)=0.d0
         overflow2d(ip) =0.d0
         nentries2d(ip) =0
         istatus(id)  =0
         nev2d(ip)=0

       elseif(idime.eq.3)then

       do i=1,nx3d(ip)
        do j=1,ny3d(ip)
         do k=1,nz3d(ip)
         histo3d(ip,i,j,k)=0.d0  
         histo3dv(ip,i,j,k)=0.d0 
         histo3dc(ip,i,j,k)=0.d0 
         histo3dp(ip,i,j,k)=0.d0  
         enddo
        enddo
       enddo

         underflow3d(ip)=0.d0
         overflow3d(ip) =0.d0
         nentries3d(ip) =0
         istatus(id)  =0
         nev3d(ip)=0

       endif


99       end

C************************************************
C       Reset and redefine an histogram
C************************************************
       subroutine ulhreuse1(id,tit,nx,xmin,xmax)

       use ulhglobal
       implicit none
       save
       integer(kind=4) :: id,nx,idime,ip,i
       character(len=64) :: tit
       real(kind=8) ::xmin,xmax
       
       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHREUSE1: histogram does not exist, ID=',id
         stop
        endif

       if(idime /= 1)then
         write(iue,*)'ULHREUSE1: not 1d, id:',id
         stop
       endif

       if(nx.le.0)then
         write(iue,*)
     &   'ULHREUSE1: number channels must be positive, ID=',id
         goto 99
       endif

       if(nx.gt.nmax1d)then
         write(iue,*)'ULHREUSE1: too many channels, ID=',id
         goto 99
       endif       

       !call ulhputid(id,ip,1) ! put ID into database
       !if(ip.eq.0) goto 99    ! an error occurred -> exit

c      reset all channels
       do i=1,nx
         histo1d(ip,i)=0.
         histo1dv(ip,i)=0.
         histo1dc(ip,i)=0.
         histo1dp(ip,i)=0.
       enddo

       nx1d(ip)=nx                 ! number of channels
       xmin1d(ip)=xmin             ! lower bound
       xmax1d(ip)=xmax             ! upper bound
       binx1d(ip)=(xmax-xmin)/nx   ! bin width
       underflow1d(ip)=0.
       overflow1d(ip)=0.
       nentries1d(ip)=0
       istatus(id)=0
       nev1d(ip)=0

       title1d(ip)=tit

99       end

C************************************************
C       Reset and redefine an histogram
C************************************************
       subroutine ulhreuse2(id,tit,nx,xmin,xmax,ny,ymin,ymax)

       use ulhglobal
       implicit none
       save
       integer(kind=4) :: id,nx,ny,idime,ip,i,j
       character(len=64) :: tit
       real(kind=8) ::xmin,xmax,ymin,ymax
       
       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHREUSE2: histogram does not exist, ID=',id
         stop
        endif

       if(idime /= 2)then
         write(iue,*)'ULHREUSE2: not 2d, id:',id
         stop
       endif

       if(nx.le.0 .or. ny.le.0)then
         write(iue,*)
     &   'ULHREUSE2: number channels must be positive, ID=',id
         goto 99
       endif       

       if(nx.gt.nmax2dx .or. ny.gt.nmax2dy)then
         write(iue,*)'ULHREUSE2: too many channels, ID=',id
         goto 99
       endif       

c      reset all channels
       do i=1,nx
         do j=1,ny
            histo2d(ip,i,j)=0.0d0
            histo2dv(ip,i,j)=0.0d0
            histo2dc(ip,i,j)=0.0d0
            histo2dp(ip,i,j)=0.0d0
         enddo
       enddo

       nx2d(ip)=nx                 ! number of x channels
       ny2d(ip)=ny                 ! number of y channels
       xmin2d(ip)=xmin             ! lower x bound
       xmax2d(ip)=xmax             ! upper x bound
       ymin2d(ip)=ymin             ! lower y bound
       ymax2d(ip)=ymax             ! upper y bound
       binx2d(ip)=(xmax-xmin)/nx   ! x bin width       
       biny2d(ip)=(ymax-ymin)/ny   ! y bin width
       underflow2d(ip)=0.
       overflow2d(ip)=0.
       nentries2d(ip)=0
       nev2d(ip)=0
       istatus(id)=0
       title2d(ip)=tit

99       end

C************************************************
C       Reset and redefine an histogram
C************************************************
       subroutine ulhreuse3
     & (id,tit,nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax)


       use ulhglobal
       implicit none
       save
       integer(kind=4) :: id,nx,ny,nz,idime,ip,i,j,k
       character(len=64) :: tit
       real(kind=8) ::xmin,xmax,ymin,ymax,zmin,zmax
       
       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHREUSE3: histogram does not exist, ID=',id
         stop
        endif

       if(idime /= 3)then
         write(iue,*)'ULHREUSE3: not 2d, id:',id
         stop
       endif

       if(nx.le.0 .or. ny.le.0 .or. nz .le.0)then
         write(iue,*)
     &   'ULHREUSE3: number channels must be positive, ID=',id
         goto 99
       endif       

       if(nx.gt.nmax3dx .or. ny.gt.nmax3dy .or. nz.gt.nmax3dz)then
         write(iue,*)'ULHREUSE3: too many channels, ID=',id
         goto 99
       endif       

c      reset all channels
       do i=1,nx
        do j=1,ny
         do k=1,nz
            histo3d(ip,i,j,k)=0.0d0
            histo3dv(ip,i,j,k)=0.0d0
            histo3dc(ip,i,j,k)=0.0d0
            histo3dp(ip,i,j,k)=0.0d0
         enddo
        enddo
       enddo

       nx3d(ip)=nx                 ! number of x channels
       ny3d(ip)=ny                 ! number of y channels
       nz3d(ip)=nz                 ! number of z channels
       xmin3d(ip)=xmin             ! lower x bound
       xmax3d(ip)=xmax             ! upper x bound
       ymin3d(ip)=ymin             ! lower y bound
       ymax3d(ip)=ymax             ! upper y bound
       zmin3d(ip)=zmin             ! lower z bound
       zmax3d(ip)=zmax             ! upper z bound
       binx3d(ip)=(xmax-xmin)/nx   ! x bin width       
       biny3d(ip)=(ymax-ymin)/ny   ! y bin width
       binz3d(ip)=(zmax-zmin)/nz   ! z bin width
       underflow3d(ip)=0.
       overflow3d(ip)=0.
       nentries3d(ip)=0
       nev3d(ip)=0
       istatus(id)=0
       title3d(ip)=tit

99       end

C************************************************
C       Set histogram title
C************************************************
       subroutine ulhsetitle(id,tit)

       use ulhglobal
       save
       character(len=64) :: tit
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHSETITLE: histogram does not exist, ID=',id
         goto 99
        endif

       select case (idime)
        case(1)
        title1d(ip)=tit
        case(2)
        title2d(ip)=tit
        case(3) 
        title3d(ip)=tit
       end select

99       end

C************************************************
C       Set histogram type
C************************************************
c      itype=1  Frequency histogram (default). Error = sqrt(content)
c               Filling with ulhadd
c      itype=2  Event average histogram        content=sum/n, error = sigma/sqrt(n)
c               Filling with ulhave
c      itype=3  Bin average histogram          content=sum/n_i, error = sigma/sqrt(n_i)
c               Filling with ulhave
c      itype=4  Quantity histogram. Q*Q is also filled and kept in final array 
c               Filling with ulhadd
c      itype=5  Similar to itype=4 but after ulhend each bin is divided the number of entries in the bin
c               Filling with ulhadd
c      itype=6  Quantity and error are given. Errors are squared summed. After closing Error=sqrt(sum error^2)
c               Filling with ulhaddq  

       subroutine ulhsetype(id,itype)

       use ulhglobal
       save
       integer(kind=4) :: id,itype,ip,idime

          call ulhgetip(id,ip,idime) ! check ID in database
           if(idime.eq.0)then
            write(iue,*)'ULHSETYPE: histogram does not exist, ID=',id
            return
           endif

       idtyp(id)=itype

       select case (itype)

       case(1,2,3,4,5,6)
!       nothing to do
        return 
         case default
            write(iue,*)'ULHSETYPE: unknown itype - id,itype',id,itype
            return
        end select

        end

C************************************************
C       Set total number of events
C************************************************
       subroutine ulhsetnev(nev)

       use ulhglobal
       save
       integer(kind=4) :: nev

        if(nev .le. 0)then
         write(iue,*)'ULHSETNEV: invalid number of events'
         return
        endif

        ntotevents=nev

        end

C================================================
C       Filling routines
C================================================


C************************************************
C       1-d histogram filling 
C       Add a quantity q to histogram contents
c       No event weights are allowed
C************************************************
       subroutine ulhadd1(id,x,q)

       use ulhglobal
       save
       real(kind=8) :: x,q
       integer(kind=4) :: id
       integer :: i

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHADD1: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHADD1: Operation no long allowed'  
         write(iue,*)
     &   'ULHADD1: histo filling only before calling ULHEND'
         goto 99
       endif

       i=int((x-xmin1d(ip))/binx1d(ip))+1  
       nentries1d(ip)=nentries1d(ip)+1     ! total no. entries  
       nev1d(ip)=nev1d(ip)+1  ! no. events

       itype=idtyp(id)
       if(i.ge.1.and.i.le.nx1d(ip))then
         histo1d(ip,i)=histo1d(ip,i)+q
         histo1dc(ip,i)=histo1dc(ip,i)+1.d0 ! entries in channel

         select case(itype)
          case(1)
           histo1dv(ip,i)=histo1dv(ip,i)+q ! error=sqrt(q)
          case(4,5)
           histo1dv(ip,i)=histo1dv(ip,i)+q*q ! allow variance computation
          end select

       elseif(i.le.0)then
         underflow1d(ip)=underflow1d(ip)+1.
       else
          overflow1d(ip)=overflow1d(ip)+1.
       endif

99       end

C************************************************
C       2-d histogram filling 
C       Add a quantity q to histogram contents
c       No event weights are allowed
C************************************************
       subroutine ulhadd2(id,x,y,q)

       use ulhglobal
       save
       real(kind=8) :: x,y,q
       integer(kind=4) :: id
       integer :: i,j

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHADD2: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHADD2: Operation no long allowed'  
         write(iue,*)
     &  'ULHADD2: histo filling only before calling ULHEND'
         goto 99
       endif

       i=int((x-xmin2d(ip))/binx2d(ip))+1
       j=int((y-ymin2d(ip))/biny2d(ip))+1

       nentries2d(ip)=nentries2d(ip)+1 ! total no. entries
       nev2d(ip)=nev2d(ip)+1 ! no. events

       itype=idtyp(id)
       if( (i.ge.1.and.i.le.nx2d(ip)) .and.
     &      (j.ge.1.and.j.le.ny2d(ip)) )then
         histo2d(ip,i,j)=histo2d(ip,i,j)+q
         histo2dc(ip,i,j)=histo2dc(ip,i,j)+1.d0 ! entries in channel

         select case(itype)
          case(1)
           histo2dv(ip,i,j)=histo2dv(ip,i,j)+q ! error=sqrt(q)
          case(4,5)
           histo2dv(ip,i,j)=histo2dv(ip,i,j)+q*q ! allow variance computation
          end select

       elseif(i.le.0 .or. j.le.0)then
         underflow2d(ip)=underflow2d(ip)+1.
       else
          overflow2d(ip)=overflow2d(ip)+1.
       endif

99       end

C************************************************
C       3-d histogram filling 
C       Add a quantity q to histogram contents
c       No event weights are allowed
C************************************************
       subroutine ulhadd3(id,x,y,z,q)

       use ulhglobal
       save
       real(kind=8) :: x,y,z,q
       integer(kind=4) :: id
       integer :: i,j,k

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHADD3: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHADD3: Operation no long allowed'  
         write(iue,*)
     &  'ULHADD3: histo filling only before calling ULHEND'
         goto 99
       endif

       i=int((x-xmin3d(ip))/binx3d(ip))+1
       j=int((y-ymin3d(ip))/biny3d(ip))+1
       k=int((z-zmin3d(ip))/binz3d(ip))+1

       nentries3d(ip)=nentries3d(ip)+1 ! total no. entries
       nev3d(ip)=nev3d(ip)+1   ! no. valid entries

       itype=idtyp(id)
       if( (i.ge.1.and.i.le.nx3d(ip)) .and.
     &     (j.ge.1.and.j.le.ny3d(ip)) .and. 
     &     (k.ge.1.and.k.le.nz3d(ip))       )then

         histo3d(ip,i,j,k)=histo3d(ip,i,j,k)+q
         histo3dc(ip,i,j,k)=histo3dc(ip,i,j,k)+1.d0 ! entries in channel

         select case(itype)
          case(1)
           histo3dv(ip,i,j,k)=histo3dv(ip,i,j,k)+q ! error=sqrt(q)
          case(4,5)
           histo3dv(ip,i,j,k)=histo3dv(ip,i,j,k)+q*q ! ! allow variance computation
          end select

       elseif(i.le.0 .or. j.le.0 .or. k.le.0 )then
         underflow3d(ip)=underflow3d(ip)+1.
       else
          overflow3d(ip)=overflow3d(ip)+1.
       endif

99       end

C************************************************
C       1-d histogram filling for average type 
C       histogram. Accumulate a quantity q to 
C       histogram contents on a event by event
C       basis. No event weights are allowed.
C************************************************
       subroutine ulhave1(id,x,q) 

       use ulhglobal
       save
       real(kind=8) :: x,q
       integer(kind=4) :: id
       integer :: i

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHAVE1: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHAVE1: Operation no long allowed'  
         write(iue,*)
     &   'ULHAVE1: histo filling only before calling ULHEND'
         goto 99
       endif
       
c       idtyp(id)= 2 !   event average
                     !   no check is made on previous assignment 
       if(idtyp(id).ne.2 .and. idtyp(id).ne.3) then  
        write(iue,*)"ULHAVE1: histo type .ne. 2 .or. 3 - id, itype",
     &               id,idtyp(id)
        return
       endif

       i=int((x-xmin1d(ip))/binx1d(ip))+1

C      Accumulate q in partial counter
       nentries1d(ip)=nentries1d(ip)+1  ! total nb of entries in histogram . This number is greater than no. of events
       if( (i.ge.1.and.i.le.nx1d(ip)) )then
           histo1dp(ip,i)=histo1dp(ip,i)+q
         elseif(i.le.0)then
           underflow1d(ip)=underflow1d(ip)+1.
         else
            overflow1d(ip)=overflow1d(ip)+1.
         endif

99       end

C************************************************
C       2-d histogram filling for average type 
C       histogram. Accumulate a quantity q to 
C       histogram contents on a event by event
C       basis. No event weights are allowed
C************************************************
       subroutine ulhave2(id,x,y,q)  

       use ulhglobal
       save
       real(kind=8) :: x,y,q
       integer(kind=4) :: id
       integer :: i,j

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHAVE2: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHAVE2: Operation no long allowed'  
         write(iue,*)
     &   'ULHAVE2: histo filling only before calling ULHEND'
         goto 99
       endif


       if(idtyp(id).ne.2 .and. idtyp(id).ne.3) then  
        write(iue,*)"ULHAVE2: histo type .ne. 2 .or. 3 - id, itype",
     &               id,idtyp(id)
        return
       endif

       i=int((x-xmin2d(ip))/binx2d(ip))+1
       j=int((y-ymin2d(ip))/biny2d(ip))+1

C         Accumulate q in partial counter

       nentries2d(ip)=nentries2d(ip)+1 ! nb of entries in histogram
       if( (i.ge.1.and.i.le.nx2d(ip)) .and.
     &      (j.ge.1.and.j.le.ny2d(ip)) )then

           histo2dp(ip,i,j)=histo2dp(ip,i,j)+q
         elseif(i.le.0 .or. j.le.0)then
           underflow2d(ip)=underflow2d(ip)+1.
         else
            overflow2d(ip)=overflow2d(ip)+1.
         endif

99       end


C************************************************
C       3-d histogram filling for average type 
C       histogram. Accumulate a quantity q to 
C       histogram contents on a event by event
C       basis. No event weights are allowed
C************************************************
       subroutine ulhave3(id,x,y,z,q)

       use ulhglobal
       save
       real(kind=8) :: x,y,z,q
       integer(kind=4) :: id
       integer :: i,j,k

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHAVE3: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHAVE3: Operation no long allowed'  
         write(iue,*)
     &   'ULHAVE3: histo filling only before calling ULHEND'
         goto 99
       endif


       if(idtyp(id).ne.2 .and. idtyp(id).ne.3) then  
        write(iue,*)"ULHAVE3: histo type .ne. 2 .or. 3 - id, itype",
     &               id,idtyp(id)
        return
       endif
       i=int((x-xmin3d(ip))/binx3d(ip))+1
       j=int((y-ymin3d(ip))/biny3d(ip))+1
       k=int((z-zmin3d(ip))/binz3d(ip))+1
C         Accumulate q in partial counter

       nentries3d(ip)=nentries3d(ip)+1 ! nb of entries in histogram
       if( (i.ge.1.and.i.le.nx3d(ip)).and.
     &     (j.ge.1.and.j.le.ny3d(ip)).and. 
     &     (k.ge.1.and.k.le.nz3d(ip))        )then

           histo3dp(ip,i,j,k)=histo3dp(ip,i,j,k)+q
         elseif(i.le.0 .or. j.le.0 .or. k.le.0)then
           underflow3d(ip)=underflow3d(ip)+1.d0
         else
            overflow3d(ip)=overflow3d(ip)+1.d0
         endif

99       end

C**************************************************
C       1-d histogram filling 
C       Add a quantity q and its error to histogram
C**************************************************
       subroutine ulhaddq1(id,x,q,er)

       use ulhglobal
       save
       real(kind=8) :: x,q,er
       integer(kind=4) :: id
       integer :: i

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHADDQ1: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHADDQ1: Operation no long allowed'  
         write(iue,*)
     &   'ULHADDQ1: histo filling only before calling ULHEND'
         goto 99
       endif

       i=int((x-xmin1d(ip))/binx1d(ip))+1  
       nentries1d(ip)=nentries1d(ip)+1     ! total no. entries  
       nev1d(ip)=nev1d(ip)+1  ! no. events

       itype=idtyp(id)
       if(itype /= 6) then 
         write(iue,*)'ULHADDQ1: histo itype /=6 ID:',id
         stop
       endif

       if(i.ge.1.and.i.le.nx1d(ip))then
         histo1d(ip,i)=histo1d(ip,i)+q
         histo1dc(ip,i)=histo1dc(ip,i)+1.d0 ! entries in channel

         select case(itype)
          case(6)
           histo1dv(ip,i)=histo1dv(ip,i)+er*er ! variance
          end select

       elseif(i.le.0)then
         underflow1d(ip)=underflow1d(ip)+1.
       else
          overflow1d(ip)=overflow1d(ip)+1.
       endif

99       end

C************************************************
C       2-d histogram filling 
C       Add a quantity q to histogram contents
c       No event weights are allowed
C************************************************
       subroutine ulhaddq2(id,x,y,q,er)

       use ulhglobal
       save
       real(kind=8) :: x,y,q,er
       integer(kind=4) :: id
       integer :: i,j

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHADDQ2: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHADDQ2: Operation no long allowed'  
         write(iue,*)
     &  'ULHADDQ2: histo filling only before calling ULHEND'
         goto 99
       endif

       i=int((x-xmin2d(ip))/binx2d(ip))+1
       j=int((y-ymin2d(ip))/biny2d(ip))+1

       nentries2d(ip)=nentries2d(ip)+1 ! total no. entries
       nev2d(ip)=nev2d(ip)+1 ! no. events

       itype=idtyp(id)
       itype=idtyp(id)
       if(itype /= 6) then 
         write(iue,*)'ULHADDQ2: histo itype /=6 ID:',id
         stop
       endif
       if( (i.ge.1.and.i.le.nx2d(ip)) .and.
     &      (j.ge.1.and.j.le.ny2d(ip)) )then
         histo2d(ip,i,j)=histo2d(ip,i,j)+q
         histo2dc(ip,i,j)=histo2dc(ip,i,j)+1.d0 ! entries in channel

         select case(itype)
          case(6)
           histo2dv(ip,i,j)=histo2dv(ip,i,j)+er*er ! variance
          end select

       elseif(i.le.0 .or. j.le.0)then
         underflow2d(ip)=underflow2d(ip)+1.
       else
          overflow2d(ip)=overflow2d(ip)+1.
       endif

99       end

C************************************************
C       3-d histogram filling 
C       Add a quantity q to histogram contents
c       No event weights are allowed
C************************************************
       subroutine ulhaddq3(id,x,y,z,q,er)

       use ulhglobal
       save
       real(kind=8) :: x,y,z,q,er
       integer(kind=4) :: id
       integer :: i,j,k

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHADDQ3: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHADDQ3: Operation no long allowed'  
         write(iue,*)
     &  'ULHADDQ3: histo filling only before calling ULHEND'
         goto 99
       endif

       i=int((x-xmin3d(ip))/binx3d(ip))+1
       j=int((y-ymin3d(ip))/biny3d(ip))+1
       k=int((z-zmin3d(ip))/binz3d(ip))+1

       nentries3d(ip)=nentries3d(ip)+1 ! total no. entries
       nev3d(ip)=nev3d(ip)+1   ! no. valid entries

       itype=idtyp(id)
       itype=idtyp(id)
       if(itype /= 6) then 
         write(iue,*)'ULHADDQ3: histo itype /=6 ID:',id
         stop
       endif
       if( (i.ge.1.and.i.le.nx3d(ip)) .and.
     &     (j.ge.1.and.j.le.ny3d(ip)) .and. 
     &     (k.ge.1.and.k.le.nz3d(ip))       )then

         histo3d(ip,i,j,k)=histo3d(ip,i,j,k)+q
         histo3dc(ip,i,j,k)=histo3dc(ip,i,j,k)+1.d0 ! entries in channel

         select case(itype)
          case(6)
           histo3dv(ip,i,j,k)=histo3dv(ip,i,j,k)+er*er ! ! allow variance computation
          end select

       elseif(i.le.0 .or. j.le.0 .or. k.le.0 )then
         underflow3d(ip)=underflow3d(ip)+1.
       else
          overflow3d(ip)=overflow3d(ip)+1.
       endif

99       end



C************************************************
C       Event endding average type histogram.
C************************************************
       subroutine ulhevend(id)

       use ulhglobal
       save
       real(kind=8) :: q
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHEVEND: histogram does not exist, ID=',id
         goto 99
        endif

       if(istatus(id).ne.0) then
         write(iue,*)'ULHEVEND: Operation no long allowed'  
         write(iue,*)
     &   'ULHEVEND: histo filling only before calling ULHEND'
         goto 99
       endif


       if(idtyp(id).ne.2 .and. idtyp(id).ne.3) then  
        write(iue,*)"ULHEVEND: histo type .ne.2 .or.3: id, itype",
     &               id,idtyp(id)
        return
       endif

       select case (idime)

       case(1)
        do k=1,nx1d(ip)
         if(histo1dp(ip,k).ne.0) then
         q=histo1dp(ip,k)
         histo1d(ip,k)=histo1d(ip,k)+q ! transfer partial counter
         histo1dv(ip,k)=histo1dv(ip,k)+q*q 
         histo1dc(ip,k)=histo1dc(ip,k)+1.0d0 ! update event bin counter
         histo1dp(ip,k)=0.0D0 ! reset partial counters
         endif
        enddo
        nev1d(ip)=nev1d(ip)+1  ! total nb of events 

       case(2)
        do k=1,nx2d(ip)
         do l=1,ny2d(ip)
         if(histo2dp(ip,k,l).ne.0) then
         q=histo2dp(ip,k,l)
         histo2d(ip,k,l)= histo2d(ip,k,l)+q ! transfer partial counter
         histo2dv(ip,k,l)= histo2dv(ip,k,l)+q*q
         histo2dc(ip,k,l)= histo2dc(ip,k,l)+1.0d0
         histo2dp(ip,k,l)=0.0D0 ! reset partial counters
         endif
         enddo
        enddo
        nev2d(ip)=nev2d(ip)+1  ! total nb of events 

       case(3)
        do k1=1,nx3d(ip)
         do k2=1,ny3d(ip)
          do k3=1,nz3d(ip)
           if(histo3dp(ip,k1,k2,k3).ne.0) then
            q=histo3dp(ip,k1,k2,k3)
            histo3d(ip,k1,k2,k3)= histo3d(ip,k1,k2,k3)+q ! transfer partial counter
            histo3dv(ip,k1,k2,k3)= histo3dv(ip,k1,k2,k3)+q*q
            histo3dc(ip,k1,k2,k3)= histo3dc(ip,k1,k2,k3)+1.0d0     
            histo3dp(ip,k1,k2,k3)=0.0D0 ! reset partial counters
           endif
          enddo
         enddo
        enddo
        nev3d(ip)=nev3d(ip)+1  ! total nb of events that called ulhave3

       case default
        write(iue,*)"ULHAVE_EVEND: unknown histo dimension id:",id
       end select

99       end



C       Input/output routines
C================================================

C************************************************
C       Output an histogram to a file
C************************************************
       subroutine ulhout(id)
       
       use ulhglobal
       save
       real(kind=8) :: x,y,z
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHOUT: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.eq.1)then    ! 1-d histograms ----------------
       nx=nx1d(ip)
       write(iuo,'(''#nch:'',i6,'' id='',i6,6x,a64)')
     & nx,id,title1d(ip)
       write(iuo,21)'x      ','content  ','error    '
       do i=1,nx1d(ip)
         x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2
C         outputs channel content and standard error
         write(iuo,11)x,histo1d(ip,i),histo1dv(ip,i)

       enddo

        if(istatus(id).eq.0)then
         write(iue,*)'ULHOUT: uncertainties not correct, ID=',id
        endif

       elseif(idime.eq.2)then     ! 2-d histograms ----------------

       nx=nx2d(ip)
       ny=ny2d(ip)
       write(iuo,'(''#nch:'',i6,1x,i6,'' id='',i6,6x,a64)')
     & nx,ny,id,title2d(ip)
      ! write(iuo,25)'  id=',id,title2d(ip)
       write(iuo,22)'x      ','y      ','content  ','error    '
       do i=1,nx2d(ip)
         do j=1,ny2d(ip)
           x=xmin2d(ip)+i*binx2d(ip)-binx2d(ip)/2
           y=ymin2d(ip)+j*biny2d(ip)-biny2d(ip)/2
C           outputs channel content and standard error
           write(iuo,12)x,y,histo2d(ip,i,j),histo2dv(ip,i,j)
         enddo
        write(iuo,*) ! write a blank line at the end for gnuplot
       enddo

        if(istatus(id).eq.0)then
         write(iue,*)'ULHOUT: uncertainties not correct, ID=',id
        endif

       elseif(idime.eq.3)then     ! 3-d histograms ----------------

       nx=nx3d(ip)
       ny=ny3d(ip)
       nz=nz3d(ip)
       write(iuo,'(''#nch:'',i6,1x,i6,1x,i6,'' id='',i6,6x,a64)')
     & nx,ny,nz,id,title3d(ip)

       write(iuo,23)
     & 'x      ','y      ','z      ','content  ','error    '
       do i=1,nx3d(ip)
         do j=1,ny3d(ip)
          do k=1,nz3d(ip)
           x=xmin3d(ip)+i*binx3d(ip)-binx3d(ip)/2
           y=ymin3d(ip)+j*biny3d(ip)-biny3d(ip)/2
           z=zmin3d(ip)+k*binz3d(ip)-binz3d(ip)/2
C           outputs channel content and standard error
           write(iuo,13)x,y,z,histo3d(ip,i,j,k),histo3dv(ip,i,j,k)
          enddo
          write(iuo,*) ! write a blank line at the end for gnuplot
        enddo
        write(iuo,*) ! write a blank line at the end for gnuplot
       enddo

        if(istatus(id).eq.0)then
         write(iue,*)'ULHOUT: uncertainties not correct, ID=',id
        endif

       endif


11       format(1x,3(1x,e12.6))
12       format(1x,4(1x,e12.6))
13       format(1x,5(1x,e12.6))
21       format('#',3(1x,a12))
22       format('#',4(1x,a12))
23       format('#',5(1x,a12))


99       end

c*****************************************************************
c Histogram output routine with zero content channel suppression.
c This format is cannot be used if histogram is to read by ulhrin.
c*****************************************************************
       subroutine ulhoutnozero(id)
   
       use ulhglobal
       save
       real(kind=8) :: x,y,z
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHOUT: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.eq.1)then    ! 1-d histograms ----------------
       nx=nx1d(ip)
       write(iuo,'(''#nch:'',i6,'' id='',i6,6x,a64)')
     & nx,id,title1d(ip)
       write(iuo,21)'x      ','content  ','error    '
       do i=1,nx1d(ip)
         x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2
C         outputs channel content and standard error
         if(histo1d(ip,i) > 0.)then
           write(iuo,11)x,histo1d(ip,i),histo1dv(ip,i)
         endif

       enddo

        if(istatus(id).eq.0)then
         write(iue,*)'ULHOUT: uncertainties not correct, ID=',id
        endif

       elseif(idime.eq.2)then     ! 2-d histograms ----------------

       nx=nx2d(ip)
       ny=ny2d(ip)
       write(iuo,'(''#nch:'',i6,1x,i6,'' id='',i6,6x,a64)')
     & nx,ny,id,title2d(ip)
      ! write(iuo,25)'  id=',id,title2d(ip)
       write(iuo,22)'x      ','y      ','content  ','error    '
       do i=1,nx2d(ip)
         do j=1,ny2d(ip)
           x=xmin2d(ip)+i*binx2d(ip)-binx2d(ip)/2
           y=ymin2d(ip)+j*biny2d(ip)-biny2d(ip)/2
C           outputs channel content and standard error
         if(histo2d(ip,i,j) > 0.)then
           write(iuo,12)x,y,histo2d(ip,i,j),histo2dv(ip,i,j)
         endif
         enddo
        !write(iuo,*) ! write a blank line at the end for gnuplot
       enddo

        if(istatus(id).eq.0)then
         write(iue,*)'ULHOUT: uncertainties not correct, ID=',id
        endif

       elseif(idime.eq.3)then     ! 3-d histograms ----------------

       nx=nx3d(ip)
       ny=ny3d(ip)
       nz=nz3d(ip)
       write(iuo,'(''#nch:'',i6,1x,i6,1x,i6,'' id='',i6,6x,a64)')
     & nx,ny,nz,id,title3d(ip)

       write(iuo,23)
     & 'x      ','y      ','z      ','content  ','error    '
       do i=1,nx3d(ip)
         do j=1,ny3d(ip)
          do k=1,nz3d(ip)
           x=xmin3d(ip)+i*binx3d(ip)-binx3d(ip)/2
           y=ymin3d(ip)+j*biny3d(ip)-biny3d(ip)/2
           z=zmin3d(ip)+k*binz3d(ip)-binz3d(ip)/2
C           outputs channel content and standard error
         if(histo3d(ip,i,j,k) > 0.)then
           write(iuo,13)x,y,z,histo3d(ip,i,j,k),histo3dv(ip,i,j,k)
         endif
          enddo
          !write(iuo,*) ! write a blank line at the end for gnuplot
        enddo
        !write(iuo,*) ! write a blank line at the end for gnuplot
       enddo

        if(istatus(id).eq.0)then
         write(iue,*)'ULHOUT: uncertainties not correct, ID=',id
        endif

       endif


11       format(1x,3(1x,e12.6))
12       format(1x,4(1x,e12.6))
13       format(1x,5(1x,e12.6))
21       format('#',3(1x,a12))
22       format('#',4(1x,a12))
23       format('#',5(1x,a12))


99       end



c******************************************************
C      Special routine to output 2D VEUSZ histograms
c******************************************************
c      to import the data to VEUSZ
C      import-> 2D -> give a name (anyone) to the datset in the 
c      dataset field

       subroutine ulhoutvz2d(id)

       use ulhglobal
       save

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHOUT: histogram does not exist, ID=',id
         return
        endif

       if(idime.eq.1)then
         write(iue,*)'ULHOUTVZ2D: not a 2D histogram, ID=',id
         return
        endif
       
         write(iuo,*)'xrange',xmin2d(ip),xmax2d(ip)
         write(iuo,*)'yrange',ymin2d(ip),ymax2d(ip)
c        write x left (xmin) to right (xmax), y bottom (ymin) to top (ymax)
         do j=ny2d(ip),1,-1  
            write(iuo,*) ( histo2d(ip,i,j), i=1,nx2d(ip) )
         enddo

         end

C************************************************
C       List an histogram summary
C************************************************
       subroutine ulhlist(idin)

       use ulhglobal
       save
       integer(kind=4) :: idin

C       if idin=0 output all histograms summary
       if(idin.eq.0)then

        do i=1,nids1d ! 1-d histos
         ip=i
         id=ids1d(i)
         idime=iddim(id)
         call ulhoutls(id,ip,idime)
        enddo

        do i=1,nids2d ! 2-d histos
         ip=i
         id=ids2d(i)
         idime=iddim(id)
         call ulhoutls(id,ip,idime)
        enddo

        do i=1,nids3d ! 3-d histos
         ip=i
         id=ids3d(i)
         idime=iddim(id)
         call ulhoutls(id,ip,idime)
        enddo

       else

       call ulhgetip(idin,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHLIST: histogram does not exist, ID=',id
         goto 99
        endif

       call ulhoutls(idin,ip,idime)
       endif

99       end

C************************************************
C       Output histogram summary list
C************************************************
       subroutine ulhoutls(id,ip,idime)

       use ulhglobal
       save
       real(kind=8) :: sumt
       integer(kind=4) :: id,idime,nevent

       if(idime.eq.1)then           ! 1D histograms --------------------------
       
c        if(idtyp(id).eq.1)then      ! find number of events that filled histo
            nevent= nentries1d(ip)-int(underflow1d(ip))
     &                            -int(overflow1d(ip))
c         else
c           nevent=nev1d(ip)
c         endif


C       compute histo integral
       sumt=0.0D0
       do i=1,nx1d(ip)
         sumt=sumt+histo1d(ip,i)     ! sum of contents
       enddo

10       format(1x,a9,2x,a64)
11       format(1x,a22,i12)
12       format(1x,a22,g12.6)

        write(ius,*)'****************************************'
        write(ius,10)'*  title:',title1d(ip)
        write(ius,11)'*  histogram ID:      ', id
        write(ius,11)'*  histogram type:    ', idtyp(id)
        write(ius,11)'*  status:            ', istatus(id)
        write(ius,11)'*  number of channels:', nx1d(ip)
        write(ius,12)'*  xmin:              ', xmin1d(ip)             
        write(ius,12)'*  xmax:              ', xmax1d(ip)             
        write(ius,12)'*  bin size:          ', binx1d(ip)  
        write(ius,12)'*  nb of entries:     ', nentries1d(ip)*1.d0
        write(ius,12)'*  nb valid entries:  ', nevent*1.d0
        write(ius,12)'*  sum of contents:   ', sumt
        write(ius,12)'*  underflow:         ', underflow1d(ip)
        write(ius,12)'*  overflow:          ', overflow1d(ip)
        write(ius,*)'****************************************'

       elseif(idime.eq.2)then       ! 2D histograms -------------------------

c        if(idtyp(id).eq.1)then      ! find number of events that filled histo
            nevent= nentries2d(ip)-int(underflow2d(ip))
     &                            -int(overflow2d(ip))
c         else
c           nevent=nev2d(ip)
c         endif

C       compute histo integral
       sumt=0.0D0
       do i=1,nx2d(ip)
        do j=1,ny2d(ip)
         sumt=sumt+histo2d(ip,i,j)
        enddo
       enddo

        write(ius,*)'****************************************'
        write(ius,10)'*  title:',title2d(ip)
        write(ius,11)'*  histogram ID:      ', id
        write(ius,11)'*  histogram type:    ', idtyp(id)
        write(ius,11)'*  status:            ', istatus(id)
        write(ius,11)'*  x nb of channels:  ', nx2d(ip)
        write(ius,11)'*  y nb of channels:  ', ny2d(ip)
        write(ius,12)'*  xmin:              ', xmin2d(ip)             
        write(ius,12)'*  xmax:              ', xmax2d(ip)             
        write(ius,12)'*  x bin size:        ', binx2d(ip)  
        write(ius,12)'*  ymin:              ', ymin2d(ip)             
        write(ius,12)'*  ymax:              ', ymax2d(ip)             
        write(ius,12)'*  y bin size:        ', biny2d(ip)  
        write(ius,12)'*  nb of entries:     ', nentries2d(ip)*1.d0
        write(ius,12)'*  nb valid entries:  ', nevent*1.d0
        write(ius,12)'*  sum of contents:   ', sumt
        write(ius,12)'*  underflow:         ', underflow2d(ip)
        write(ius,12)'*  overflow:          ', overflow2d(ip)
        write(ius,*)'****************************************'

       elseif(idime.eq.3)then       ! 3D histograms -------------------------

c        if(idtyp(id).eq.1)then      ! find number of events that filled histo
            nevent= nentries3d(ip)-int(underflow3d(ip))
     &                            -int(overflow3d(ip))
c         else
c           nevent=nev3d(ip)
c         endif

C       compute histo integral
       sumt=0.0D0
       do i=1,nx3d(ip)
        do j=1,ny3d(ip)
         do k=1,nz3d(ip)
          sumt=sumt+histo3d(ip,i,j,k)
         enddo
        enddo
       enddo

        write(ius,*)'****************************************'
        write(ius,10)'*  title:',title3d(ip)
        write(ius,11)'*  histogram ID:      ', id
        write(ius,11)'*  histogram type:    ', idtyp(id)
        write(ius,11)'*  status:            ', istatus(id)
        write(ius,11)'*  x nb of channels:  ', nx3d(ip)
        write(ius,11)'*  y nb of channels:  ', ny3d(ip)
        write(ius,11)'*  z nb of channels:  ', nz3d(ip)
        write(ius,12)'*  xmin:              ', xmin3d(ip)             
        write(ius,12)'*  xmax:              ', xmax3d(ip)             
        write(ius,12)'*  x bin size:        ', binx3d(ip)  
        write(ius,12)'*  ymin:              ', ymin3d(ip)             
        write(ius,12)'*  ymax:              ', ymax3d(ip)             
        write(ius,12)'*  y bin size:        ', biny3d(ip)  
        write(ius,12)'*  zmin:              ', zmin3d(ip)             
        write(ius,12)'*  zmax:              ', zmax3d(ip)             
        write(ius,12)'*  z bin size:        ', binz3d(ip)  
        write(ius,12)'*  nb of entries:     ', nentries3d(ip)*1.d0
        write(ius,12)'*  nb valid entries:  ', nevent*1.d0
        write(ius,12)'*  sum of contents:   ', sumt
        write(ius,12)'*  underflow:         ', underflow3d(ip)
        write(ius,12)'*  overflow:          ', overflow3d(ip)
        write(ius,*)'****************************************'

       endif

       end

C************************************************
C       Open file for output/input
C************************************************
       subroutine ulhopen(iunit,filename,utype)

       use ulhglobal
       save
       character(len=1) :: utype
       character(len=64) :: filename
       integer(kind=4) :: iunit       

       if(utype.eq.'e' .or. utype.eq.'E')then
         if(iue.ne. 6) close(iue) ! close error unit 1st
         iue=iunit
       elseif(utype.eq.'s' .or. utype.eq.'S')then
         if(ius.ne. 6) close(ius) ! close summary unit 1st
         ius=iunit
       elseif(utype.eq.'o' .or. utype.eq.'O')then
        close(iuo) ! close output unit 1st
        iuo=iunit
       elseif(utype.eq.'i' .or. utype.eq.'I')then
        close(iui) ! close output unit 1st
        iui=iunit
       else
        write(iue,*)'ULHOPEN: unknown unit type, ',utype
        goto 99
       endif

       open(unit=iunit,file=filename,status='unknown')

99       end

C************************************************
C       Close standard units
C************************************************
       subroutine ulhclose

       use ulhglobal
       save

       close(iue)
       close(ius)
       close(iuo)
       close(iui)

       end    

C************************************************
C       Input 1-d histogram from file
C************************************************
       subroutine ulhin1(id)

       use ulhglobal
       save
       real(kind=8), dimension(nmax1d) :: getx,getc,gete
       real(kind=8) :: binx
       integer(kind=4) :: id,n,nx
       character(len=64) :: title,comment
       character(len=18) :: idtxt


       call ulhgetip(id,ip,idime) ! check ID in database

        if(idime.eq.2)then
         write(iue,*)'ULHIN1: histogram is not 1-d ID=',id,' Exit!'
         goto 99
        elseif(idime.eq.1)then
         write(iue,*)
     &   'ULHIN1: histogram exists. Data may be lost ID=',id
         call ulhreset(id) ! reset histogram
        else    
         call ulhputid(id,ip,1) ! put ID into database 
         if(ip.eq.0) then
           write(iue,*)
     &      'ULHIN1: error in ulhputid id=',id
           goto 99    ! an error occurred -> exit
         endif

        endif


       read(iui,11,err=999,end=999)idtxt,title ! read id and title line
11       format(a18,a64)

       read(iui,12,err=999,end=999)comment ! comment line
12       format(a64)

       n=1
10       read(iui,*,err=998,end=998)getx(n),getc(n),gete(n) ! x,content,error
       n=n+1
       goto 10
998       continue    ! end of file reading 

       nx=n-1
c       set channels
       do i=1,nx
         histo1d(ip,i)=getc(i)
         histo1dv(ip,i)=gete(i)
         histo1dc(ip,i)=1.        ! no. entries set to 1.
       enddo

       binx=abs(getx(2)-getx(1))
       
       nx1d(ip)=nx                 ! number of channels
       xmin1d(ip)=getx(1)-binx/2.  ! lower bound
       xmax1d(ip)=getx(nx)+binx/2. ! upper bound
       binx1d(ip)=binx             ! bin width
       underflow1d(ip)=0.
       overflow1d(ip)=0.
       nentries1d(ip)=0
       istatus(id)=3             ! errors given by user
       nev1d(ip)=0

       title1d(ip)=title

       goto 99 ! exit

999       write(iue,*)'ULHIN1: Error reading histogram id=',id,
     &              ' from unit iui=',iui
       goto 99

99       end

C************************************************
C       Input 2-d histogram from file
C************************************************
       subroutine ulhin2(id)

       use ulhglobal
       save
       real(kind=8), dimension(nmax2dx,nmax2dy) :: getx,gety
       real(kind=8) :: x1,xn,y1,yn
       real(kind=8) :: binx,biny
       integer(kind=4) :: id,nx,ny
       character(len=64) :: title,comment
       character(len=5) :: idtxt

       call ulhgetip(id,ip,idime) ! check ID in database

        if(idime.eq.1)then
         write(iue,*)'ULHIN2: histogram is not 2-d ID=',id,' Exit!'
         goto 99
        elseif(idime.eq.2)then
         write(iue,*)
     &   'ULHIN2: histogram exists. Data may be lost ID=',id
         call ulhreset(id) ! reset histogram
        else    
         call ulhputid(id,ip,2) ! put ID into database
         if(ip.eq.0) then
           write(iue,*)
     &      'ULHIN2: error in ulhputid id=',id
           goto 99    ! an error occurred -> exit
         endif
        endif

!       read(iui,11,err=999,end=999)idtxt,title ! read id and title line
!11       format(a18,a64)

       read(iui,*,err=999,end=999)idtxt,nx,ny,title ! read id and title line
       read(iui,12,err=999,end=999)comment ! comment line
12       format(a64)

       do i=1,nx
        do j=1,ny
          read(iui,*,err=999,end=999)
     &     getx(i,j),gety(i,j),histo2d(ip,i,j),histo2dv(ip,i,j) ! x,y,content,error
           
           histo2dc(ip,i,j)=1.d0   ! no. entries set to 1.
        enddo
          read(iui,12,err=999,end=999)comment ! blank line
       enddo

       x1=getx(1,1)
       xn=getx(nx,ny)
       y1=gety(1,1)
       yn=gety(nx,ny)

       binx=abs(getx(2,1)-getx(1,1) )
       biny=abs(gety(1,2)-gety(1,1) )
       
       nx2d(ip)=nx           ! number of x channels
       ny2d(ip)=ny           ! number of y channels
       xmin2d(ip)=x1-binx/2. ! lower x bound
       xmax2d(ip)=xn+binx/2. ! upper x bound
       ymin2d(ip)=y1-biny/2. ! lower y bound
       ymax2d(ip)=yn+biny/2. ! upper y bound
       binx2d(ip)=binx       ! x bin width       
       biny2d(ip)=biny       ! y bin width
       underflow2d(ip)=0.
       overflow2d(ip)=0.
       nentries2d(ip)=0
       nev2d(ip)=0
       istatus(id)=3   
       title2d(ip)=title

       goto 99 ! exit

999       write(iue,*)'ULHIN2: Error reading histogram id=',id,
     &              ' from unit iui=',iui
       goto 99

99       end

C************************************************
C       Input 3-d histogram from file
C************************************************
       subroutine ulhin3(id)

       use ulhglobal
       save
       real(kind=8),dimension(nmax3dx,nmax3dy,nmax3dz)::getx,gety,getz
       real(kind=8) :: x1,xn,y1,yn,z1,zn
       real(kind=8) :: binx,biny,binz
       integer(kind=4) :: id,nx,ny,nz
       character(len=64) :: title,comment
       character(len=5) :: idtxt

       call ulhgetip(id,ip,idime) ! check ID in database

        if(idime.eq.1 .or. idime .eq. 2)then
         write(iue,*)'ULHIN3: histogram is not 3-d ID=',id,' Exit!'
         goto 99
        elseif(idime.eq.3)then
         write(iue,*)
     &   'ULHIN3: histogram exists. Data may be lost ID=',id
         call ulhreset(id) ! reset histogram
        else    
         call ulhputid(id,ip,3) ! put ID into database
         if(ip.eq.0) then
           write(iue,*)
     &      'ULHIN3: error in ulhputid id=',id
           goto 99    ! an error occurred -> exit
         endif
        endif

!       read(iui,11,err=999,end=999)idtxt,title ! read id and title line
!11       format(a18,a64)

       read(iui,*,err=999,end=999)idtxt,nx,ny,nz,title ! read id and title line
       read(iui,12,err=999,end=999)comment ! comment line
12       format(a64)

       do i=1,nx
        do j=1,ny
         do k=1,nz
          read(iui,*,err=999,end=999)
     &     getx(i,j,k),gety(i,j,k),getz(i,j,k),
     &     histo3d(ip,i,j,k),histo3dv(ip,i,j,k) ! x,y,content,error

           histo3dc(ip,i,j,k)=1.   ! set no. entries to 1.
         enddo
          read(iui,12,err=999,end=999)comment ! blank line
        enddo 
         read(iui,12,err=999,end=999)comment ! blank line
       enddo

       x1=getx(1,1,1)
       xn=getx(nx,ny,nz)
       y1=gety(1,1,1)
       yn=gety(nx,ny,nz)
       z1=getz(1,1,1)
       zn=getz(nx,ny,nz)

       binx=abs(getx(2,1,1)-getx(1,1,1) )
       biny=abs(gety(1,2,1)-gety(1,1,1) )
       binz=abs(getz(1,1,2)-getz(1,1,1) )
       
       nx3d(ip)=nx           ! number of x channels
       ny3d(ip)=ny           ! number of y channels
       nz3d(ip)=nz           ! number of z channels
       xmin3d(ip)=x1-binx/2. ! lower x bound
       xmax3d(ip)=xn+binx/2. ! upper x bound
       ymin3d(ip)=y1-biny/2. ! lower y bound
       ymax3d(ip)=yn+biny/2. ! upper y bound
       zmin3d(ip)=z1-binz/2. ! lower z bound
       zmax3d(ip)=zn+binz/2. ! upper z bound
       binx3d(ip)=binx       ! x bin width       
       biny3d(ip)=biny       ! y bin width
       binz3d(ip)=binz       ! z bin width
       underflow3d(ip)=0.
       overflow3d(ip)=0.
       nentries3d(ip)=0
       nev3d(ip)=0
       istatus(id)=3   
       title3d(ip)=title

       goto 99 ! exit

999       write(iue,*)'ULHIN3: Error reading histogram id=',id,
     &              ' from unit iui=',iui
       goto 99

99       end


C************************************************
C       Input and define a spectrum from file
C       This histograms don't need to be booked
C       File must first be open with ulhopen!
C************************************************
       subroutine ulhinsp(id)

       use ulhglobal
       save
       real(kind=8), dimension(nmax1d):: getx,getc
       character(len=64) :: title
       real(kind=8) :: binx
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database

        if(idime.eq.2)then
         write(iue,*)'ULHINSP: histogram is not 1-d ID=',id,' Exit!'
         goto 99
        elseif(idime.eq.1)then
         write(iue,*)
     &   'ULHINSP: histogram exists. Data may be lost ID=',id
         call ulhreset(id) ! reset histogram
        else    
         call ulhputid(id,ip,1) ! put ID into database
         if(ip.eq.0) goto 99    ! an error occurred -> exit
        endif

       read(iui,12,err=999,end=999)title ! 1st line is title
12       format(a64)

       n=1
10       read(iui,*,err=998,end=998)getx(n),getc(n) ! x,content
       n=n+1
       goto 10
998       continue    ! end of file reading 

       nx=n-1
c       set channels
       do i=1,nx
         histo1d(ip,i)=getc(i)
          histo1dv(ip,i)=0.D0
       enddo

       binx=abs(getx(2)-getx(1))
       
       nx1d(ip)=nx                 ! number of channels
       xmin1d(ip)=getx(1)-binx/2.  ! lower bound
       xmax1d(ip)=getx(nx)+binx/2. ! upper bound
       binx1d(ip)=binx             ! bin width
       underflow1d(ip)=0.
       overflow1d(ip)=0.
       nentries1d(ip)=0
       istatus(id)=3             ! errors given by user
       nev1d(ip)=0

       title1d(ip)=title

       goto 99 ! exit

999       write(iue,*)'ULHINSP: Error reading histogram id=',id,
     &              ' from unit iui=',iui
       goto 99

99       end

C       Operations on histograms
C================================================

C************************************************
C       Sum of two histograms multiplied 
C       by constants 
C************************************************
       subroutine ulhsum(id1,id2,id3,a,b)

       use ulhglobal
       save
       real(kind=8) :: x,y,dx,dy,q,dq,a,b
       integer(kind=4) :: id1,id2,id3

       call ulhgetip(id1,ip1,idime1) ! check ID in database
        if(idime1.eq.0)then
         write(iue,*)'ULHSUM: histogram ID1 does not exist, ID=',id1
         return
        endif

       call ulhgetip(id2,ip2,idime2) ! check ID in database
        if(idime2.eq.0)then
         write(iue,*)'ULHSUM: histogram ID2 does not exist, ID=',id2
         return
        endif

       if(idime1.ne.idime2)then
         write(iue,*)
     &    'ULHSUM: histograms have diferent dimensions, IDs=',id1,id2
         return
        endif
       
        ist1=istatus(id1)
        ist2=istatus(id2)

        if(ist1.eq.0 .or. ist2.eq.0) then
         write(iue,*)
     &   'ULHSUM: Operation not allowed. Call ULHEND 1st, IDs=',id1,id2
         return
        endif

        call ulhgetip(id3,ip3,idime3) ! check ID in database
        select case(idime3)
         case(0)           ! id3 does not exist
          call ulhcpp(id1,id3) 
         case(1,2,3)
          if(idime3 .ne. idime1) then
           write(iue,*)
     &     'ERROR ULHSUM: ID1, ID3 with different dimensions',id1,id3
           return
          endif
          call ulhcpp(id1,id3) ! copy id1 properties to id3 (but not contents of id1)
         end select

       call ulhgetip(id3,ip3,idime3)

       if(idime1.eq.1)then        !-----------------------------------------

       do i=1,nx1d(ip1)
         x=histo1d(ip1,i)
         y=histo1d(ip2,i)
         dx=histo1dv(ip1,i)
         dy=histo1dv(ip2,i)
         q=x*a + y*b
         dq=sqrt(a*a*dx*dx + b*b*dy*dy)
         histo1d(ip3,i)=q
         histo1dv(ip3,i)=dq
         histo1dc(ip3,i)=histo1dc(ip1,i)+histo1dc(ip2,i)
       enddo
       
       nentries1d(ip3)=nentries1d(ip1)+nentries1d(ip2)
       overflow1d(ip3)=overflow1d(ip1)+overflow1d(ip2)
       underflow1d(ip3)=underflow1d(ip1)+underflow1d(ip2) 

       elseif(idime1.eq.2)then     !----------------------------------------

       do i=1,nx2d(ip1)
         do j=1,ny2d(ip1)

         x=histo2d(ip1,i,j)
         y=histo2d(ip2,i,j)
         dx=histo2dv(ip1,i,j)
         dy=histo2dv(ip2,i,j)
         q=x*a + y*b
         dq=sqrt(a*a*dx*dx + b*b*dy*dy)
         histo2d(ip3,i,j)=q
         histo2dv(ip3,i,j)=dq
         histo2dc(ip3,i,j)=histo2dc(ip1,i,j)+histo2dc(ip2,i,j)
         enddo
       enddo

       nentries2d(ip3)=nentries2d(ip1)+nentries2d(ip2)
       overflow2d(ip3)=overflow2d(ip1)+overflow2d(ip2)
       underflow2d(ip3)=underflow2d(ip1)+underflow2d(ip2)

       elseif(idime1.eq.3)then      !-----------------------------------------

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
          do k=1,nz3d(ip1)

         x=histo3d(ip1,i,j,k)
         y=histo3d(ip2,i,j,k)
         dx=histo3dv(ip1,i,j,k)
         dy=histo3dv(ip2,i,j,k)
         q=x*a + y*b
         dq=sqrt(a*a*dx*dx + b*b*dy*dy)
         histo3d(ip3,i,j,k)=q
         histo3dv(ip3,i,j,k)=dq
         histo3dc(ip3,i,j,k)=histo3dc(ip1,i,j,k)+histo3dc(ip2,i,j,k)
         enddo
        enddo
       enddo

       nentries3d(ip3)=nentries3d(ip1)+nentries3d(ip2)
       overflow3d(ip3)=overflow3d(ip1)+overflow3d(ip2)
       underflow3d(ip3)=underflow3d(ip1)+underflow3d(ip2)

       endif
       
       istatus(id3)=3
       end

C************************************************
C       Multiplication of two histograms 
C       multiplied by a constant 
C************************************************
       subroutine ulhmul(id1,id2,id3,a)

       use ulhglobal
       save
       real(kind=8) :: x,y,dx,dy,q,dq,a
       integer(kind=4) :: id1,id2,id3

       call ulhgetip(id1,ip1,idime1) ! check ID in database
        if(idime1.eq.0)then
         write(iue,*)'ULHMUL: histogram ID1 does not exist, ID=',id1
         goto 99
        endif

       call ulhgetip(id2,ip2,idime2) ! check ID in database
        if(idime2.eq.0)then
         write(iue,*)'ULHMUL: histogram ID2 does not exist, ID=',id2
         goto 99
        endif

       if(idime1.ne.idime2)then
         write(iue,*)
     &    'ULHMUL: histograms have diferent dimensions, IDs=',id1,id2
         goto 99
        endif

        ist1=istatus(id1)
        ist2=istatus(id2)

        if(ist1.eq.0 .or. ist2.eq.0) then
         write(iue,*)
     &   'ULHMUL: Operation not allowed. Call ULHEND 1st, IDs=',id1,id2
         goto 99
        endif

        call ulhgetip(id3,ip3,idime3) ! check ID in database
        select case(idime3)
         case(0)           ! id3 does not exist
          call ulhcpp(id1,id3) 
         case(1,2,3)
          if(idime3 .ne. idime1) then
           write(iue,*)
     &     'ERROR ULHSUM: ID1, ID3 with different dimensions',id1,id3
           return
          endif
          call ulhcpp(id1,id3) ! copy id1 properties to id3 (but not contents of id1)
         end select

       call ulhgetip(id3,ip3,idime3)

       if(idime1.eq.1)then

       do i=1,nx1d(ip1)
         x=histo1d(ip1,i)
         y=histo1d(ip2,i)
         dx=histo1dv(ip1,i)
         dy=histo1dv(ip2,i)
         q=x*y* a
         dq=a*sqrt(y*y*dx*dx + x*x*dy*dy)
         histo1d(ip3,i)=q
         histo1dv(ip3,i)=dq
         histo1dc(ip3,i)=histo1dc(ip1,i)+histo1dc(ip2,i)
       enddo

       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
        do j=1,ny2d(ip1)
         x=histo2d(ip1,i,j)
         y=histo2d(ip2,i,j)
         dx=histo2dv(ip1,i,j)
         dy=histo2dv(ip2,i,j)
         q=x*y* a
         dq=a*sqrt(y*y*dx*dx + x*x*dy*dy)
         histo2d(ip3,i,j)=q
         histo2dv(ip3,i,j)=dq
         histo2dc(ip3,i,j)=histo2dc(ip1,i,j)+histo2dc(ip2,i,j)
        enddo
       enddo

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
         x=histo3d(ip1,i,j,k)
         y=histo3d(ip2,i,j,k)
         dx=histo3dv(ip1,i,j,k)
         dy=histo3dv(ip2,i,j,k)
         q=x*y* a
         dq=a*sqrt(y*y*dx*dx + x*x*dy*dy)
         histo3d(ip3,i,j,k)=q
         histo3dv(ip3,i,j,k)=dq
         histo3dc(ip3,i,j,k)=histo3dc(ip1,i,j,k)+histo3dc(ip2,i,j,k)
        enddo
        enddo
       enddo

       endif

         istatus(id3)=3
99       end

C************************************************
C       Divison of two histograms multiplied 
C       by a constant 
C************************************************
       subroutine ulhdiv(id1,id2,id3,a)

       use ulhglobal
       save
       real(kind=8) :: x,y,dx,dy,q,dq,a
       integer(kind=4) :: id1,id2,id3

       call ulhgetip(id1,ip1,idime1) ! check ID in database
        if(idime1.eq.0)then
         write(iue,*)'ULHDIV: histogram ID1 does not exist, ID=',id1
         goto 99
        endif

       call ulhgetip(id2,ip2,idime2) ! check ID in database
        if(idime2.eq.0)then
         write(iue,*)'ULHDIV: histogram ID2 does not exist, ID=',id2
         goto 99
        endif

       if(idime1.ne.idime2)then
         write(iue,*)
     &    'ULHDIV: histograms have diferent dimensions, IDs=',id1,id2
         goto 99
        endif

        ist1=istatus(id1)
        ist2=istatus(id2)

        if(ist1.eq.0 .or. ist2.eq.0) then
         write(iue,*)
     &   'ULHDIV: Operation not allowed. Call ULHEND 1st, IDs=',id1,id2
         goto 99
        endif

        call ulhgetip(id3,ip3,idime3) ! check ID in database
        select case(idime3)
         case(0)           ! id3 does not exist
          call ulhcpp(id1,id3) 
         case(1,2,3)
          if(idime3 .ne. idime1) then
           write(iue,*)
     &     'ERROR ULHSUM: ID1, ID3 with different dimensions',id1,id3
           return
          endif
          call ulhcpp(id1,id3) ! copy id1 properties to id3 (but not contents of id1)
         end select

       call ulhgetip(id3,ip3,idime3)

       if(idime1.eq.1)then

       do i=1,nx1d(ip1)
         x=histo1d(ip1,i)
         y=histo1d(ip2,i)
         dx=histo1dv(ip1,i)
         dy=histo1dv(ip2,i)
         if(y.ne.0)then
          q=x/y * a
          dq=a* sqrt( (dx*dx)/(y*y) + (x*x*dy*dy)/(y*y*y*y) )          
         else
          q=0.
          dq=0.
         endif
         histo1d(ip3,i)=q
         histo1dv(ip3,i)=dq
         histo1dc(ip3,i)=histo1dc(ip1,i)+histo1dc(ip2,i)
       enddo

       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
         do j=1,ny2d(ip1)
          x=histo2d(ip1,i,j)
          y=histo2d(ip2,i,j)
          dx=histo2dv(ip1,i,j)
          dy=histo2dv(ip2,i,j)
          if(y.ne.0)then
           q=x/y * a
           dq=a* sqrt( (dx*dx)/(y*y) + (x*x*dy*dy)/(y*y*y*y) )          
          else
           q=0.
           dq=0.
          endif
          histo2d(ip3,i,j)=q
          histo2dv(ip3,i,j)=dq
          histo2dc(ip3,i,j)=histo2dc(ip1,i,j)+histo2dc(ip2,i,j)
         enddo
       enddo

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
         x=histo3d(ip1,i,j,k)
         y=histo3d(ip2,i,j,k)
         dx=histo3dv(ip1,i,j,k)
         dy=histo3dv(ip2,i,j,k)
          if(y.ne.0)then
           q=x/y * a
           dq=a* sqrt( (dx*dx)/(y*y) + (x*x*dy*dy)/(y*y*y*y) )          
          else
           q=0.
           dq=0.
          endif
          histo3d(ip3,i,j,k)=q
          histo3dv(ip3,i,j,k)=dq
          histo3dc(ip3,i,j,k)=histo3dc(ip1,i,j,k)+histo3dc(ip2,i,j,k)
         enddo
        enddo
       enddo

       endif

         istatus(id3)=3
99       end


C************************************************
C       Multiplication of an histogram
C       by a constant 
C************************************************
       subroutine ulhscale(id1,id3,a)

       use ulhglobal
       save
       real(kind=8) :: x,dx,q,dq,a
       integer(kind=4) :: id1,id3

       call ulhgetip(id1,ip1,idime1) ! check ID in database
        if(idime1.eq.0)then
         write(iue,*)'ULHSCALE: histogram ID1 does not exist, ID=',id1
         goto 99
        endif

        ist1=istatus(id1)

        if(ist1.eq.0) then
         write(iue,*)
     &   'ULHSCALE: Operation not allowed. Call ULHEND 1st, IDs=',
     &   id1,id2
         goto 99
        endif


        call ulhgetip(id3,ip3,idime3) ! check ID in database
        select case(idime3)
         case(0)           ! id3 does not exist
          call ulhcpp(id1,id3) 
         case(1,2,3)
          if(idime3 .ne. idime1) then
           write(iue,*)
     &     'ERROR ULHSUM: ID1, ID3 with different dimensions',id1,id3
           return
          endif
          call ulhcpp(id1,id3) ! copy id1 properties to id3 (but not contents of id1)
         end select

       call ulhgetip(id3,ip3,idime3)


       if(idime1.eq.1)then

       do i=1,nx1d(ip1)
         x=histo1d(ip1,i)
         dx=histo1dv(ip1,i)         
         q=x* a
         dq=a*dx
         histo1d(ip3,i)=q
         histo1dv(ip3,i)=dq
         histo1dc(ip3,i)=histo1dc(ip1,i)
       enddo

       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
        do j=1,ny2d(ip1)
         x=histo2d(ip1,i,j)
         dx=histo2dv(ip1,i,j)         
         q=x* a
         dq=a*dx
         histo2d(ip3,i,j)=q
         histo2dv(ip3,i,j)=dq
         histo2dc(ip3,i,j)=histo2dc(ip1,i,j)
        enddo
       enddo

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
           x=histo3d(ip1,i,j,k)
           dx=histo3dv(ip1,i,j,k)         
           q=x* a
           dq=a*dx
           histo3d(ip3,i,j,k)=q
           histo3dv(ip3,i,j,k)=dq
          histo3dc(ip3,i,j,k)=histo3dc(ip1,i,j,k)
         enddo
        enddo
       enddo

       endif

         istatus(id3)=3
99       end

C************************************************
C       Square root of an histogram
C************************************************
       subroutine ulhsqrt(id1,id3)

       use ulhglobal
       save
       real(kind=8) :: x,dx,q,dq
       integer(kind=4) :: id1,id3

       call ulhgetip(id1,ip1,idime1) ! check ID in database
        if(idime1.eq.0)then
         write(iue,*)'ULHSQRT: histogram ID1 does not exist, ID=',id1
         goto 99
        endif

        ist1=istatus(id1)

        if(ist1.eq.0) then
         write(iue,*)
     &   'ULHSQRT: Operation not allowed. Call ULHEND 1st, IDs=',
     &   id1,id2
         goto 99
        endif

        call ulhgetip(id3,ip3,idime3) ! check ID in database
        select case(idime3)
         case(0)           ! id3 does not exist
          call ulhcpp(id1,id3) 
         case(1,2,3)
          if(idime3 .ne. idime1) then
           write(iue,*)
     &     'ERROR ULHSUM: ID1, ID3 with different dimensions',id1,id3
           return
          endif
          call ulhcpp(id1,id3) ! copy id1 properties to id3 (but not contents of id1)
         end select

       call ulhgetip(id3,ip3,idime3)

       if(idime1.eq.1)then

       do i=1,nx1d(ip1)
         x=histo1d(ip1,i)
         dx=histo1dv(ip1,i)

         if(x.gt.0.)then         
          q=sqrt(x)
          dq=0.5*dx/q
         else
          q=0.
          dq=0.
         endif

         histo1d(ip3,i)=q
         histo1dv(ip3,i)=dq
         histo1dc(ip3,i)=histo1dc(ip1,i)
       enddo
       
       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
        do j=1,ny2d(ip1)

         x=histo2d(ip1,i,j)
         dx=histo2dv(ip1,i,j)

         if(x.gt.0.)then         
          q=sqrt(x)
          dq=0.5*dx/q
         else
          q=0.
          dq=0.
         endif

         histo2d(ip3,i,j)=q
         histo2dv(ip3,i,j)=dq
         histo2dc(ip3,i,j)=histo2dc(ip1,i,j)
        enddo
       enddo

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
           x=histo3d(ip1,i,j,k)
           dx=histo3dv(ip1,i,j,k)  
       
           if(x.gt.0.)then         
            q=sqrt(x)
            dq=0.5*dx/q
           else
            q=0.
            dq=0.
           endif

           histo3d(ip3,i,j,k)=q
           histo3dv(ip3,i,j,k)=dq
          histo3dc(ip3,i,j,k)=histo3dc(ip1,i,j,k)
         enddo
        enddo
       enddo

       endif

         istatus(id3)=3
99       end

C************************************************
C       Exponencial of an histogram 
C************************************************
       subroutine ulhexp(id1,id3,b)

       use ulhglobal
       save
       real(kind=8) :: x,dx,q,dq,b
       integer(kind=4) :: id1,id3

       call ulhgetip(id1,ip1,idime1) ! check ID in database
        if(idime1.eq.0)then
         write(iue,*)'ULHEXP: histogram ID1 does not exist, ID=',id1
         goto 99
        endif

        ist1=istatus(id1)

        if(ist1.eq.0) then
         write(iue,*)
     &   'ULHEXP: Operation not allowed. Call ULHEND 1st, IDs=',id1,id2
         goto 99
        endif

        call ulhgetip(id3,ip3,idime3) ! check ID in database
        select case(idime3)
         case(0)           ! id3 does not exist
          call ulhcpp(id1,id3) 
         case(1,2,3)
          if(idime3 .ne. idime1) then
           write(iue,*)
     &     'ERROR ULHSUM: ID1, ID3 with different dimensions',id1,id3
           return
          endif
          call ulhcpp(id1,id3) ! copy id1 properties to id3 (but not contents of id1)
         end select

       call ulhgetip(id3,ip3,idime3)

       if(idime1.eq.1)then

       do i=1,nx1d(ip1)
         x=histo1d(ip1,i)
         dx=histo1dv(ip1,i)
         q=exp(b*x)
         dq=b*q*dx

         histo1d(ip3,i)=q
         histo1dv(ip3,i)=dq
         histo1dc(ip3,i)=histo1dc(ip1,i)
       enddo

       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
        do j=1,ny2d(ip1)
         x=histo2d(ip1,i,j)
         dx=histo2dv(ip1,i,j)

         q=exp(b*x)
         dq=b*q*dx

         histo2d(ip3,i,j)=q
         histo2dv(ip3,i,j)=dq
         histo2dc(ip3,i,j)=histo2dc(ip1,i,j)
        enddo
       enddo

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
           x=histo3d(ip1,i,j,k)
           dx=histo3dv(ip1,i,j,k)    
     
            q=exp(b*x)
            dq=b*q*dx

           histo3d(ip3,i,j,k)=q
           histo3dv(ip3,i,j,k)=dq
           histo3dc(ip3,i,j,k)=histo3dc(ip1,i,j,k)
         enddo
        enddo
       enddo

       endif

         istatus(id3)=3
99       end

C************************************************
C       Natural logarithm of an histogram 
C************************************************
       subroutine ulhlog(id1,id3)

       use ulhglobal
       save
       real(kind=8) :: x,dx,q,dq
       integer(kind=4) :: id1,id3

       call ulhgetip(id1,ip1,idime1) ! check ID in database
        if(idime1.eq.0)then
         write(iue,*)'ULHLOG: histogram ID1 does not exist, ID=',id1
         goto 99
        endif

        ist1=istatus(id1)

        if(ist1.eq.0) then
         write(iue,*)
     &   'ULHLOG: Operation not allowed. Call ULHEND 1st, IDs=',id1,id2
         goto 99
        endif

        call ulhgetip(id3,ip3,idime3) ! check ID in database
        select case(idime3)
         case(0)           ! id3 does not exist
          call ulhcpp(id1,id3) 
         case(1,2,3)
          if(idime3 .ne. idime1) then
           write(iue,*)
     &     'ERROR ULHSUM: ID1, ID3 with different dimensions',id1,id3
           return
          endif
          call ulhcpp(id1,id3) ! copy id1 properties to id3 (but not contents of id1)
         end select

       call ulhgetip(id3,ip3,idime3)

       if(idime1.eq.1)then

       do i=1,nx1d(ip1)
         x=histo1d(ip1,i)
         dx=histo1dv(ip1,i)

         if(x.gt.0.)then         
          q=log(x)
          dq=dx/abs(x)
         else
          q=0.
          dq=0.
         endif

         histo1d(ip3,i)=q
         histo1dv(ip3,i)=dq
         histo1dc(ip3,i)=histo1dc(ip1,i)
       enddo

       elseif(idime1.eq.2)then

       do i=1,nx2d(ip1)
        do j=1,ny2d(ip1)
         x=histo2d(ip1,i,j)
         dx=histo2dv(ip1,i,j)

         if(x.gt.0.)then         
          q=log(x)
          dq=dx/abs(x)
         else
          q=0.
          dq=0.
         endif

         histo2d(ip3,i,j)=q
         histo2dv(ip3,i,j)=dq
         histo2dc(ip3,i,j)=histo2dc(ip1,i,j)
        enddo
       enddo

       elseif(idime1.eq.3)then

       do i=1,nx3d(ip1)
        do j=1,ny3d(ip1)
         do k=1,nz3d(ip1)
           x=histo3d(ip1,i,j,k)
           dx=histo3dv(ip1,i,j,k)    
     
           if(x.gt.0.)then         
            q=log(x)
            dq=dx/abs(x)
           else
            q=0.
            dq=0.
           endif

           histo3d(ip3,i,j,k)=q
           histo3dv(ip3,i,j,k)=dq
           histo3dc(ip3,i,j,k)=histo3dc(ip1,i,j,k)
         enddo
        enddo
       enddo

       endif

         istatus(id3)=3
99     end


C================================================
C       Close histograms for filling
C================================================



C************************************************
C       Get standard deviations for a 
C       frequency histogram. 
C************************************************
c      Errors are equal to the sqrt of channel content        itype=1
c      Errors are equal to the sqrt of sum err*err            itype=6
c      Error contains Q*Q                                     itype=4,5


       subroutine ulhaddEND(id)          

       use ulhglobal
       save
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHADDEND: histogram does not exist, ID=',id
         goto 99
        endif

        ist1=istatus(id)

       if(ist1.ne.0) then
         write(iue,*)'ULHADDEND: Histogram status not zero, ID=',id
         goto 99
       endif

         itype=idtyp(id)

         if(idime.eq.1)then

           select case(itype)
           case(1,6)
            do i=1,nx1d(ip)         
             histo1dv(ip,i)=sqrt(histo1dv(ip,i))
            enddo
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           case(4)         
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           case(5)
            do i=1,nx1d(ip)         
             if(histo1dc(ip,i)>0)then
              histo1d(ip,i)=histo1d(ip,i)/histo1dc(ip,i)
              histo1dv(ip,i)=histo1dv(ip,i)/histo1dc(ip,i)
             endif
            enddo
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           end select

         elseif(idime.eq.2)then

           select case(itype)
           case(1,6)
           do i=1,nx2d(ip)
            do j=1,ny2d(ip)         
             histo2dv(ip,i,j)=sqrt(histo2dv(ip,i,j))
            enddo
           enddo
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           case(4)
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           case(5)
           do i=1,nx2d(ip)
            do j=1,ny2d(ip)    
             if(histo2dc(ip,i,j)>0)then     
              histo2d(ip,i,j)=histo2d(ip,i,j)/histo2dc(ip,i,j)
              histo2dv(ip,i,j)=histo2dv(ip,i,j)/histo2dc(ip,i,j)
             endif
            enddo
           enddo
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           end select

         elseif(idime.eq.3)then

           select case(itype)
           case(1,6)
           do i=1,nx3d(ip)
            do j=1,ny3d(ip)
             do k=1,nz3d(ip)
               histo3dv(ip,i,j,k)=sqrt(histo3dv(ip,i,j,k))
             enddo
            enddo
           enddo
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           case(4)
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           case(5)
           do i=1,nx3d(ip)
            do j=1,ny3d(ip)
             do k=1,nz3d(ip)
              if(histo3dc(ip,i,j,k)>0)then
               histo3d(ip,i,j,k)=histo3d(ip,i,j,k)/histo3dc(ip,i,j,k)
               histo3dv(ip,i,j,k)=histo3dv(ip,i,j,k)/histo3dc(ip,i,j,k)
              endif
             enddo
            enddo
           enddo
           istatus(id)=1  ! the histogram status is changed to frequency type
           return

           end select

         endif

99       end

C************************************************
C       Get the average values and SD for  
C       Event average histograms
C************************************************
c
c
       subroutine ulhaveEND(id)                        

       use ulhglobal
       save
       real(kind=8) :: Q,Q2
       integer(kind=4) :: id,ist1,Nev
       integer :: i,j,k

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHAVEEND: histogram does not exist, ID=',id
         goto 99
        endif

        ist1=istatus(id)

       if(ist1.ne.0) then
         write(iue,*)'ULHAVEEND: Histogram status not zero, ID=',id
         goto 99
       endif

       itype=idtyp(id)

       if(idime.eq.1)then  ! 1D histograms ------------------------------------

          if(nev1d(ip) .eq. 0) then
             write(iue,*)'ULHAVEEND: No entries ID=',id
             return
          endif

 
          if(ntotevents .ne.0) then
            nev=ntotevents          ! no. events set
          else
            nev=nev1d(ip)           ! no. of events with valid entries
          endif


       select case(itype)
        case(2)
         do i=1,nx1d(ip)
          Q=histo1d(ip,i)/Nev   ! average value of physical quantity Q
          Q2=histo1dv(ip,i)/Nev ! average of Q*Q2
          histo1d(ip,i)=Q     ! contents are substituted by average
          if(Q2-Q*Q .ge.0)then
           histo1dv(ip,i)=sqrt((Q2-Q*Q)/Nev) ! standard deviation of Q
          else
           histo1dv(ip,i)=0.
          endif   
         enddo

        case(3)   
         do i=1,nx1d(ip)
          if(histo1dc(ip,i) .ne. 0)then
            Q=histo1d(ip,i)/histo1dc(ip,i)   ! bin average value of physical quantity Q
            Q2=histo1dv(ip,i)/histo1dc(ip,i) ! average of Q*Q
            histo1d(ip,i)=Q     ! contents are substituted by average
           if(Q2-Q*Q .ge.0)then
            histo1dv(ip,i)=sqrt((Q2-Q*Q)/histo1dc(ip,i)) ! standard deviation of Q
           else
            histo1dv(ip,i)=0.
           endif  
           else
            histo1d(ip,i)=0.
            histo1dv(ip,i)=0.
          endif 
         enddo

        end select


       elseif(idime.eq.2)then ! 2D histograms ---------------------------

          if(nev2d(ip) .eq. 0) then
             write(iue,*)'ULHAVEEND: No entries ID=',id
             return
          endif

          if(ntotevents .ne.0) then
            nev=ntotevents
          else
            nev=nev2d(ip)
          endif

       select case(itype)
        case(2)
         do i=1,nx2d(ip)
          do j=1,ny2d(ip)
           Q=histo2d(ip,i,j)/Nev   ! average value of physical quantity Q
           Q2=histo2dv(ip,i,j)/Nev ! average of Q*Q
           histo2d(ip,i,j)=Q     ! contents are substituted by average
           if(Q2-Q*Q .ge.0)then
            histo2dv(ip,i,j)=sqrt((Q2-Q*Q)/Nev) ! standard deviation of Q
           else
            histo2dv(ip,i,j)=0.
           endif
          enddo   
         enddo

        case(3)   
         do i=1,nx2d(ip)
          do j=1,ny2d(ip)
          if(histo2dc(ip,i,j) .ne. 0)then
            Q=histo2d(ip,i,j)/histo2dc(ip,i,j)   ! bin average value of physical quantity Q
            Q2=histo2dv(ip,i,j)/histo2dc(ip,i,j) ! average of Q*Q
            histo2d(ip,i,j)=Q     ! contents are substituted by average
           if(Q2-Q*Q .ge.0)then
            histo2dv(ip,i,j)=sqrt((Q2-Q*Q)/histo2dc(ip,i,j)) ! standard deviation of Q
           else
            histo2dv(ip,i,j)=0.
           endif  
           else
            histo2d(ip,i,j)=0.
            histo2dv(ip,i,j)=0.
          endif 
          enddo
         enddo 

        end select


       elseif(idime.eq.3)then ! 3D histograms ---------------------------

          if(nev3d(ip) .eq. 0) then
             write(iue,*)'ULHAVEEND: No entries ID=',id
             return
          endif

          if(ntotevents .ne.0) then
            nev=ntotevents
          else
            nev=nev3d(ip)
          endif

       select case(itype)
        case(2)
         do i=1,nx3d(ip)
          do j=1,ny3d(ip)
           do k=1,nz3d(ip)
           Q=histo3d(ip,i,j,k)/Nev   ! average value of physical quantity Q
           Q2=histo3dv(ip,i,j,k)/Nev ! average of Q*Q
           histo3d(ip,i,j,k)=Q     ! contents are substituted by average
           if(Q2-Q*Q .ge.0)then
            histo3dv(ip,i,j,k)=sqrt((Q2-Q*Q)/Nev) ! standard deviation of Q
           else
            histo3dv(ip,i,j,k)=0.
           endif
           enddo
          enddo   
         enddo

        case(3)   
         do i=1,nx3d(ip)
          do j=1,ny3d(ip)
           do k=1,nz3d(ip)
            if(histo3dc(ip,i,j,k) .ne. 0)then
               Q=histo3d(ip,i,j,k)/histo3dc(ip,i,j,k)   ! bin average value of physical quantity Q
               Q2=histo3dv(ip,i,j,k)/histo3dc(ip,i,j,k) ! average of Q*Q
               histo3d(ip,i,j,k)=Q     ! contents are substituted by average
               if(Q2-Q*Q .ge.0)then
                histo3dv(ip,i,j,k)=sqrt((Q2-Q*Q)/histo3dc(ip,i,j,k)) ! standard deviation of Q
               else
                histo3dv(ip,i,j,k)=0.
               endif  
            else
              histo3d(ip,i,j,k)=0.
              histo3dv(ip,i,j,k)=0.
            endif 
           enddo
          enddo 
         enddo

        end select

       endif

        istatus(id)=2   ! the histogram status is changed to average type

99     end


c------------------------------------------------------
        subroutine ulhend(idin)
c-------------------------------------------------------
c
c       End the filling of an histogram and comput errors
c       if ID=0 End all histograms


       use ulhglobal
       save
       integer (kind=4) :: idin

       if(idin.eq.0)then
          do i=1,nids1d ! 1-d histos
           ip=i
           id=ids1d(i)
           if(id <= nuserhistos) call ulhendID(id)
          enddo
          do i=1,nids2d ! 2-d histos
           ip=i
           id=ids2d(i)
           if(id <= nuserhistos) call ulhendID(id)
          enddo
          do i=1,nids3d ! 3-d histos
           ip=i
           id=ids3d(i)
           if(id <= nuserhistos) call ulhendID(id)
          enddo

        else ! only one ID

         call ulhgetip(idin,ip,idime) ! check ID in database
          if(idime.eq.0)then
            write(iue,*)'ULHEND: histogram does not exist, ID=',id
            return
          endif
           call ulhendID(idin)

        endif

        end

c------------------------------------------------------
        subroutine ulhendID(id)
c-------------------------------------------------------
c
c       End the filling of an histogram and comput errors
c
       use ulhglobal
       save
       integer(kind=4) :: id,itype

       itype=idtyp(id)
       select case (itype)
       case (1,4,5,6)        ! frequency type histograms
         call ulhaddEND(id)
       case(2,3)             ! average quantity
         call ulhaveEND(id)
       case default
          write(iue,*)'ULHEND: type unknown ID, type =',id,itype
       end select

       end

C================================================
C       Using information
C================================================


c************************************************
c       Retrives the information on the histogram
c       1-  nx
c       2-  xmin
c       3-  xmax
c       4-  binx
c       5-  ny
c       6-  ymin
c       7-  ymax
c       8-  biny
c       9-  nz
c       10- zmin
c       11- zmax
c       12- binz
c       13- nentries (number of entries) 
c       14- underflow
c       15- overflow
c       16- sum of contents
c       17- dimension
c       18- status
c       19- type
c       20- number of events ( <= no. entries)
c************************************************

       subroutine ulhinfo(id,hinfo)

       use ulhglobal
       save

       real(kind=8), dimension(20) :: hinfo
       real(kind=8) :: sumt
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHINFO: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.eq.1)then

C       compute histo integral
       sumt=0.0D0
       do i=1,nx1d(ip)
         sumt=sumt+histo1d(ip,i)     ! sum of contents
       enddo

        hinfo(1)=nx1d(ip)*1.0D0
        hinfo(2)=xmin1d(ip)             
        hinfo(3)=xmax1d(ip)             
        hinfo(4)=binx1d(ip)
        hinfo(5)=0.D0
        hinfo(6)=0.D0             
        hinfo(7)=0.D0             
        hinfo(8)=0.D0   
        hinfo(9)=0.D0
        hinfo(10)=0.D0             
        hinfo(11)=0.D0             
        hinfo(12)=0.D0   
        hinfo(13)=nentries1d(ip)*1.0D0
        hinfo(14)=underflow1d(ip)
        hinfo(15)=overflow1d(ip)
        hinfo(16)=sumt
        hinfo(17)=idime*1.0D0
        hinfo(18)=istatus(id)*1.0D0
        hinfo(19)=idtyp(id)*1.d0
        if(idtyp(id).eq.1)then
            hinfo(20)=hinfo(13)
         else
            hinfo(20)=nev1d(ip)*1.d0
         endif


       elseif(idime.eq.2)then

C       compute histo integral
       sumt=0.0D0
       do i=1,nx2d(ip)
        do j=1,ny2d(ip)
         sumt=sumt+histo2d(ip,i,j)
        enddo
       enddo

        hinfo(1)=nx2d(ip)*1.0D0
        hinfo(2)=xmin2d(ip)             
        hinfo(3)=xmax2d(ip)             
        hinfo(4)=binx2d(ip)
        hinfo(5)=ny2d(ip)*1.0D0
        hinfo(6)=ymin2d(ip)             
        hinfo(7)=ymax2d(ip)             
        hinfo(8)=biny2d(ip) 
        hinfo(9)=0.D0
        hinfo(10)=0.D0             
        hinfo(11)=0.D0             
        hinfo(12)=0.D0   
        hinfo(13)=nentries2d(ip)*1.0D0
        hinfo(14)=underflow2d(ip)
        hinfo(15)=overflow2d(ip)
        hinfo(16)=sumt
        hinfo(17)=idime*1.0D0
        hinfo(18)=istatus(id)*1.0D0
        hinfo(19)=idtyp(id)*1.d0   
        if(idtyp(id).eq.1)then
            hinfo(20)=hinfo(13)
         else
            hinfo(20)=nev2d(ip)*1.d0
         endif

       elseif(idime.eq.3)then

C       compute histo integral
       sumt=0.0D0
       do i=1,nx3d(ip)
        do j=1,ny3d(ip)
         do k=1,nz3d(ip)       
           sumt=sumt+histo3d(ip,i,j,k)
         enddo
        enddo
       enddo

        hinfo(1)=nx3d(ip)*1.0D0
        hinfo(2)=xmin3d(ip)             
        hinfo(3)=xmax3d(ip)             
        hinfo(4)=binx3d(ip)
        hinfo(5)=ny3d(ip)*1.0D0
        hinfo(6)=ymin3d(ip)             
        hinfo(7)=ymax3d(ip)             
        hinfo(8)=biny3d(ip) 
        hinfo(9)=nz3d(ip)*1.0D0
        hinfo(10)=zmin3d(ip)             
        hinfo(11)=zmax3d(ip)             
        hinfo(12)=binz3d(ip) 
        hinfo(13)=nentries3d(ip)*1.0D0
        hinfo(14)=underflow3d(ip)
        hinfo(15)=overflow3d(ip)
        hinfo(16)=sumt
        hinfo(17)=idime*1.0D0
        hinfo(18)=istatus(id)*1.0D0
        hinfo(19)=idtyp(id)*1.d0   
        if(idtyp(id).eq.1)then
            hinfo(20)=hinfo(13)
         else
            hinfo(20)=nev3d(ip)*1.d0
         endif

       endif

99     end

C************************************************
C       Computes the mean value of the x variable
C       using the contents as the (non-normalized)
C       probability distribution. The histogram
C        must be of the frequency type.
C************************************************

       subroutine ulhxmean(id,xmed,sigma,sx)

       use ulhglobal
       save
       real(kind=8) :: xmed,sigma,sx
       real(kind=8) :: x,xmed2, sumt, var
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHXMEAN: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHXMEAN: histogram must be 1-d, ID=',id
         goto 99
        endif

       sumt=0.0D0
       xmed=0.0D0
       xmed2=0.0D0
       do i=1,nx1d(ip)
        sumt=sumt+histo1d(ip,i)     ! sum of contents
        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 
        xmed=xmed+x*histo1d(ip,i)
        xmed2=xmed2+x*x*histo1d(ip,i)
       enddo
        if(sumt.gt.0.) then
              xmed=xmed/sumt     ! mean value
              xmed2=xmed2/sumt

              var=abs(xmed2-xmed*xmed) ! variance
              sigma=sqrt(var)     ! mean standard deviation
              n=nentries1d(ip)
              sx=sigma/sqrt(n*1.d0) ! mean value error 
                                    !error is underestimate in case of underflow or overflow
                                 
        else
              xmed=0.d0
              sigma=0.
              sx=0.
        endif


99       end

C************************************************
C       Computes the mean value,sigma and error
C       of the x variable between xmin and xmax, 
C       using the contents as the (non-normalized)
C       probability distribution. The histogram
C        must be of the frequency type.
C       Useful for peaks!
C************************************************

       subroutine ulhxm(id,xmin,xmax,xmed,sigma,sx)

       use ulhglobal
       save
       real(kind=8) :: xmin,xmax,xmed,sigma,sx
       real(kind=8) :: x,xmed2,sumt
       real(kind=8) :: var
       integer(kind=4) :: id


       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHXM: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHXM: histogram must be 1-d, ID=',id
         goto 99
        endif

       imin=int((xmin-xmin1d(ip))/binx1d(ip))+1
       if(imin.lt.1) imin=1
       imax=int((xmax-xmin1d(ip))/binx1d(ip))+1
       if(imax.gt.nx1d(ip)) imax=nx1d(ip)

       sumt=0.0D0
       xmed=0.0D0
       xmed2=0.0D0
       do i=imin,imax
        sumt=sumt+histo1d(ip,i)     ! sum of contents
        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 
        xmed=xmed+x*histo1d(ip,i)
        xmed2=xmed2+x*x*histo1d(ip,i)
       enddo
        if(sumt.gt.0.) then
              xmed=xmed/sumt     ! mean value
              xmed2=xmed2/sumt
              var=abs(xmed2-xmed*xmed) ! variance
               sigma=sqrt(var)     ! mean standard deviation
               n=nentries1d(ip)
               sx=sigma/sqrt(n*1.d0) ! mean value error
                                     !error is underestimate in case of underflow or overflow
        else
              xmed=0.d0
              sigma=0.
              sx=0.
        endif
       
99       end


C*************************************************************
C       Multi type Peak analysis. The mean (centroid), sigma, 
C       total and net area with background 
C       subtraction, are computed.
C*************************************************************

       subroutine ulhanapeak
     &(id,itype,xmin,xmax,m,em,b,eb,xmed,sigma,total,net,enet,ierro)

       use ulhglobal
       save
       real*8 ::xmin,xmax,m,em,b,eb,xmed,sigma,total,net,enet
       real*8 :: x,bgd,ebgd
       real*8 :: sx,sxx,xmed2
       real*8 :: counts,ex,enet2
       integer*4 :: id,itype,ierro
       integer*4 :: i
       real*8 :: exmed
c      itype  1: bgd=m*x+b, 2: bgd=b*exp(m*x),  3: bgd=b*x**m

       ierro=0
       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHANAPEAK: histogram does not exist, ID=',id
         ierro=1
         return
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHANAPEAK: histogram must be 1-d, ID=',id
         ierro=2
         return
        endif

       imin=int((xmin-xmin1d(ip))/binx1d(ip))+1
       if(imin.lt.1) imin=1
       imax=int((xmax-xmin1d(ip))/binx1d(ip))+1
       if(imax.gt.nx1d(ip)) imax=nx1d(ip)

       total=0.0D0
       net=0.0d0
       enet2=0.d0
       sx=0.0D0
       sxx=0.d0
       ex=binx1d(ip)/2.d0

       exmed=0.
       do i=imin,imax

        counts=histo1d(ip,i)
        total=total+counts     ! sum of contents
        x=xmin1d(ip)+i*ex*2.d0-ex ! comput channel center 

        select case(itype)
         case(1)          
          bgd=m*x+b
          ebgd=sqrt( (em*x)**2 + (m*ex)**2 + eb**2 )
         case(2)          
          bgd=b*exp(m*x)
          ebgd=sqrt( (exp(m*x)*eb)**2 + 
     &    (x*b*exp(m*x)*em)**2 + (m*b*exp(m*x)*ex)**2 )
         case(3)          
          bgd=b*(x**m)
          ebgd=sqrt( ((x**m)*eb)**2 +  (b*(x**m)*log(abs(x))*em)**2 + 
     &    (m*b*(x**(m-1.d0))*ex)**2 )
         case default
          ierro=3
          return
        end select

       if(counts-bgd >= 0)then
        sx=sx+x*(counts-bgd)
        sxx=sxx+x*x*(counts-bgd)
        net=net+counts-bgd
        enet2=enet2 + counts + ebgd

        exmed=exmed+(ex*(counts-bgd))**2+(x**2*counts)+ (x**2*ebgd**2)   
        !exmed=exmed+(x**2*counts)    
       endif
  
       enddo
       
        enet=sqrt(enet2)      ! NET uncertainty
        if(net.gt.0.) then
           xmed=sx/net     ! mean value 
           xmed2=sxx/net
           sigma=sqrt(abs(xmed2-xmed*xmed))   
           exmed=sqrt(exmed)/net
           !print *, 'erro media:', sigma/sqrt(net), exmed
        else
           xmed=0.d0
           sigma=0.d0
        endif
       
       end


C***********************************************************************
C       Get one Peak parameters based on smart peak limits determination 
C       and background subtraction
C***********************************************************************
c     input
c     id: histogram, itype: 1-linear, 2-exponential, 3-power
c     roi0, roi3: limits for analysis
c     m,em,b,eb: bgd fitted parameters and errors
c     output
c     xmin,xmax: peak limits where bgd=0
c     xmed,sigma,fwhm,total,net,enet: peak parameters

       subroutine ulhget1peak
     &(id,itype,roi0,roi3,mi,emi,bi,ebi,
     & xmin,xmax,xmed,sigma,fwhm,total,net,enet,ierro)

       use ulhglobal
       implicit none
       save
       real*8 :: mi,emi,bi,ebi,chi2,r
       real*8 :: m,em,b,eb
       real*8 :: roi0,roi3
       real*8 :: xmin,xmax,xmed,sigma,total,net,enet
       real*8 :: x,bgd,ebgd,xmed2
       real*8 :: sx,sxx
       real*8 :: counts,ex,enet2,cpeak,ecpeak,signal
       real*8 :: cmax,fwhm,xl,xr
       real*8, dimension(2) :: xfit,yfit,eyfit
       integer*4 :: ichmax,i1,i2
       integer*4 :: idime,ip,imed
       integer*4 :: id,itype,ierro,ierr
       integer*4 :: i

c      bgd parameters
       m=mi
       em=emi
       b=bi
       eb=ebi

       ierro=0
       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETPEAK: histogram does not exist, ID=',id
         ierro=1
         return
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHGETPEAK: histogram must be 1-d, ID=',id
         ierro=2
         return
        endif

c       make basic analysis
        call ulhanapeak
     &  (id,itype,roi0,roi3,
     &   m,em,b,eb,xmed,sigma,total,net,enet,ierr) 

       call ulhgeti(id,xmed,imed) ! get channel for centroid guess
       call ulhgetci(id,imed,cpeak,ecpeak)

       total=0.0D0
       net=0.0d0
       enet2=0.d0
       sx=0.0D0
       sxx=0.d0
       ex=binx1d(ip)/2.d0
       cmax=0. ! max counts

c      sweep peak right side
       i=imed 
       signal=cpeak 
       do while(signal > 0.d0)
        i=i+1 ! update 1 channel
        if(i > nx1d(ip) ) exit   ! exit loop if at histogram upper limit
        counts=histo1d(ip,i)
        x=xmin1d(ip)+i*ex*2.d0-ex ! comput channel center

        select case(itype)
         case(1)          
          bgd=m*x+b
          ebgd=sqrt( (em*x)**2 + (m*ex)**2 + eb**2 )
         case(2)          
          bgd=b*exp(m*x)
          ebgd=sqrt( (exp(m*x)*eb)**2 + 
     &    (x*b*exp(m*x)*em)**2 + (m*b*exp(m*x)*ex)**2 )
         case(3)          
          bgd=b*(x**m)
          ebgd=sqrt( ((x**m)*eb)**2 +  (b*(x**m)*log(abs(x))*em)**2 + 
     &    (m*b*(x**(m-1.d0))*ex)**2 )
         case default
          ierro=3
          return
        end select

       signal=counts-bgd
       if(signal > cmax ) then
         cmax=signal
         ichmax=i
       endif

       total=total+counts     ! sum of contents
       if(signal > 0)then
        xmax=x ! update peak upper limit
        sx=sx+x*(counts-bgd)
        sxx=sxx+x*x*(counts-bgd)
        net=net+counts-bgd
        enet2=enet2 + counts + ebgd
       endif
  
       end do ! end while


c      sweep peak left side
       i=imed+1 ! start loop at imed 
       signal=cpeak 
       do while(signal > 0.d0)
        i=i-1 ! update 1 channel
        if(i < 1 ) exit   ! exit loop if at histogram lower limit
        counts=histo1d(ip,i)
        x=xmin1d(ip)+i*ex*2.d0-ex ! comput channel center 

        select case(itype)
         case(1)          
          bgd=m*x+b
          ebgd=sqrt( (em*x)**2 + (m*ex)**2 + eb**2 )
         case(2)          
          bgd=b*exp(m*x)
          ebgd=sqrt( (exp(m*x)*eb)**2 + 
     &    (x*b*exp(m*x)*em)**2 + (m*b*exp(m*x)*ex)**2 )
         case(3)          
          bgd=b*(x**m)
          ebgd=sqrt( ((x**m)*eb)**2 +  (b*(x**m)*log(abs(x))*em)**2 + 
     &    (m*b*(x**(m-1.d0))*ex)**2 )
         case default
          ierro=3
          return
        end select

       signal=counts-bgd
       if(signal > cmax ) then
         cmax=signal
         ichmax=i
       endif

       total=total+counts     ! sum of contents
       if(signal > 0)then
        xmin=x ! update peak upper limit
        sx=sx+x*(counts-bgd)
        sxx=sxx+x*x*(counts-bgd)
        net=net+counts-bgd
        enet2=enet2 + counts + ebgd
       endif
  
       end do ! end while
 
        enet=sqrt(enet2)      ! NET uncertainty
        if(net.gt.0.) then
           xmed=sx/net     ! mean value 
           xmed2=sxx/net
           sigma=sqrt(abs(xmed2-xmed*xmed))   
        else
           xmed=0.d0
           sigma=0.d0
        endif

c----- find FWHM -----------------------------

c      sweep peak right side
       i=ichmax 
       signal=cmax 
       do while(signal > cmax/2.d0)
        i=i+1 ! update 1 channel
        if(i > nx1d(ip) ) exit   ! exit loop if at histogram upper limit
        counts=histo1d(ip,i)
        x=xmin1d(ip)+i*ex*2.d0-ex ! comput channel center

        select case(itype)
         case(1)          
          bgd=m*x+b
         case(2)          
          bgd=b*exp(m*x)
         case(3)          
          bgd=b*(x**m)
         case default
          ierro=3
          return
        end select

       signal=counts-bgd  
       end do ! end while
       i1=i-1
       i2=i
       xfit(2)=histo1d(ip,i2)
       xfit(1)=histo1d(ip,i1)
       yfit(2)=xmin1d(ip)+i2*ex*2.d0-ex
       yfit(1)=xmin1d(ip)+i1*ex*2.d0-ex
       eyfit(2)=0.
       eyfit(1)=0.

       call ulhfitvec
     &  (xfit,yfit,eyfit,2,3,1,m,em,b,eb,r,chi2,ierr)

       xr=b*(cmax/2.)**m  ! right side x

c      sweep peak left side
       i=ichmax+1 
       signal=cmax 
       do while(signal > cmax/2.d0)
        i=i-1 ! update 1 channel
        if(i < 1 ) exit   ! exit loop if at histogram upper limit
        counts=histo1d(ip,i)
        x=xmin1d(ip)+i*ex*2.d0-ex ! comput channel center

        select case(itype)
         case(1)          
          bgd=m*x+b
         case(2)          
          bgd=b*exp(m*x)
         case(3)          
          bgd=b*(x**m)
         case default
          ierro=3
          return
        end select

       signal=counts-bgd  
       end do ! end while
       i1=i-1
       i2=i
       xfit(2)=histo1d(ip,i2)
       xfit(1)=histo1d(ip,i1)
       yfit(2)=xmin1d(ip)+i2*ex*2.d0-ex
       yfit(1)=xmin1d(ip)+i1*ex*2.d0-ex
       eyfit(2)=0.
       eyfit(1)=0.
       
       call ulhfitvec
     &  (xfit,yfit,eyfit,2,3,1,m,em,b,eb,r,chi2,ierr)

       xl=b*(cmax/2.)**m  ! left side x

       fwhm=xr-xl

       end



C************************************************
C       Peak analysis. The mean (centroid), sigma, 
C       total and net area with background 
C       subtraction, are computed.
C       A linear background is assumed.
C************************************************

       subroutine ulhpeak(id,xmin,xmax,xmed,sigma,total,anet)

       use ulhglobal
       save
       real(kind=8) :: xmin,xmax,xmed,sigma,total,anet
       real(kind=8) :: x,x1,x2,y1,y2,am,b,bgd
       real(kind=8) :: sx,sxx,var
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPEAK: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHPEAK: histogram must be 1-d, ID=',id
         goto 99
        endif

       imin=int((xmin-xmin1d(ip))/binx1d(ip))+1
       if(imin.lt.1) imin=1
       imax=int((xmax-xmin1d(ip))/binx1d(ip))+1
       if(imax.gt.nx1d(ip)) imax=nx1d(ip)

C       comput backgroud line
       x1=xmin
       x2=xmax
       y1=histo1d(ip,imin)
       y2=histo1d(ip,imax)
       am=(y1-y2)/(xmin-xmax) ! line slope
       b=y1-am*xmin         ! intersection

C       compute xmed
       total=0.0D0
       anet=0.0d0
       sx=0.0D0
       do i=imin,imax
        total=total+histo1d(ip,i)     ! sum of contents
        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 
        bgd=am*x+b
       if(histo1d(ip,i)-bgd.ge.0)then
        sx=sx+x*(histo1d(ip,i)-bgd)
        anet=anet+histo1d(ip,i)-bgd
       endif

       enddo
        if(anet.gt.0.) then
              xmed=sx/anet     ! mean value
        else
              xmed=0.d0
        endif
       
c       comput the corrected width
       sxx=0.
       do i=imin,imax
        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 
        bgd=am*x+b
       if((histo1d(ip,i)-bgd).ge.0)then
        sxx=sxx+(x-xmed)**2*(histo1d(ip,i)-bgd) 
        endif
       enddo

       if(anet.gt.0.) then
          var=abs(sxx)/anet
          sigma=sqrt(var) 
       else
          sigma=0.
       endif

99       end

C*************************************************************
C       Peak analysis. The mean (centroid), sigma, total and 
C       net area with background subtraction, are computed.
C       A non-linear background is assumed.
C*************************************************************
C
C      OBSOLETA: usar ulhanapeak
C
       subroutine ulhpeaklog(id,xmin,xmax,xmed,sigma,total,net)

       use ulhglobal
       save
       real(kind=8) :: xmin,xmax,xmed,sigma,total,net
       real(kind=8) :: x,x1,x2,y1,y2,bgd
       real(kind=8) :: sx,sxx,var
       integer(kind=4) :: id
       real*8 :: lnx1,lnx2,lny1,lny2
       real*8 :: m,b
       integer :: icase

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPEAK: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHPEAK: histogram must be 1-d, ID=',id
         goto 99
        endif

       imin=int((xmin-xmin1d(ip))/binx1d(ip))+1
       if(imin.lt.1) imin=1
       imax=int((xmax-xmin1d(ip))/binx1d(ip))+1
       if(imax.gt.nx1d(ip)) imax=nx1d(ip)

C       comput backgroud line
       x1=xmin
       x2=xmax
       y1=histo1d(ip,imin)
       y2=histo1d(ip,imax)
       icase=1
       if(x1>0.)then
         lnx1=log(x1)
       else
         icase=2
       endif

       if(x2>0.)then
         lnx2=log(x2)
       else
         icase=2
       endif


       if(y1>0.)then
         lny1=log(y1)
       else
         icase=2
       endif
 
       if(y2>0.)then
         lny2=log(y2)
       else
         icase=2
       endif

       select case(icase)
        case(1) ! log interpolation
           m=(lny2-lny1)/(lnx2-lnx1)
           b=lny1-m*lnx1 
        case(2) ! use linear interpolation instead
        m=(y1-y2)/(x1-x2) ! line slope
        b=y1-m*x1         ! intersection
       end select

C       compute xmed
       total=0.0D0
       net=0.0d0
       sx=0.0D0
       do i=imin,imax
        total=total+histo1d(ip,i)     ! sum of contents
        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 

        select case(icase)
        case(1)
         bgd=exp( m*log(x) + b)
        case(2)
         bgd=m*x+b
        end select

       if(histo1d(ip,i)-bgd.ge.0)then
        sx=sx+x*(histo1d(ip,i)-bgd)
        net=net+histo1d(ip,i)-bgd
       endif

       enddo
        if(net.gt.0.) then
              xmed=sx/net     ! mean value
        else
              xmed=0.d0
        endif
       
c       comput the corrected width
       sxx=0.
       do i=imin,imax
        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 

        select case(icase)
        case(1)
         bgd=exp( m*log(x) + b)
        case(2)
         bgd=m*x+b
        end select

       if((histo1d(ip,i)-bgd).ge.0)then
        sxx=sxx+(x-xmed)**2*(histo1d(ip,i)-bgd) 
        endif
       enddo

       if(net.gt.0.) then
          var=abs(sxx)/net
          sigma=sqrt(var) 
       else
          sigma=0.
       endif

99       end


C*************************************************************
C       Fill signal and bgd spectrum in given ROI 
C*************************************************************

       subroutine ulhpeakfilllog(id,idsignal,idbgd,xmin,xmax)

       use ulhglobal
       save
       real(kind=8) :: xmin,xmax
       real(kind=8) :: x,x1,x2,y1,y2,signal,bgd
       integer(kind=4) :: id,idsignal,idbgd
       real*8 :: lnx1,lnx2,lny1,lny2
       real*8 :: m,b
       integer :: icase

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPEAK: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHPEAK: histogram must be 1-d, ID=',id
         goto 99
        endif

       call ulhgetip(idsignal,ip2,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPEAK: histogram does not exist, ID=',idsignal
         goto 99
        endif

       call ulhgetip(idbgd,ip3,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPEAK: histogram does not exist, ID=',idbgd
         goto 99
        endif


       imin=int((xmin-xmin1d(ip))/binx1d(ip))+1
       if(imin.lt.1) imin=1
       imax=int((xmax-xmin1d(ip))/binx1d(ip))+1
       if(imax.gt.nx1d(ip)) imax=nx1d(ip)

C       comput backgroud line
       x1=xmin
       x2=xmax
       y1=histo1d(ip,imin)
       y2=histo1d(ip,imax)
       icase=1
       if(x1>0.)then
         lnx1=log(x1)
       else
         icase=2
       endif

       if(x2>0.)then
         lnx2=log(x2)
       else
         icase=2
       endif


       if(y1>0.)then
         lny1=log(y1)
       else
         icase=2
       endif
 
       if(y2>0.)then
         lny2=log(y2)
       else
         icase=2
       endif

       select case(icase)
        case(1) ! log interpolation
           m=(lny2-lny1)/(lnx2-lnx1)
           b=lny1-m*lnx1 
        case(2) ! use linear interpolation instead
        m=(y1-y2)/(x1-x2) ! line slope
        b=y1-m*x1         ! intersection
       end select


       do i=imin,imax

        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 

        select case(icase)
        case(1)
         bgd=exp( m*log(x) + b)
        case(2)
         bgd=m*x+b
        end select

       signal=histo1d(ip,i)-bgd
       if(signal.ge.0)then
        call ulhadd1(idsignal,x,signal)
        call ulhadd1(idbgd,x,bgd)
       endif

       enddo

99       end

C******************************************************************
C       Computes the fwhm of a distribution.
C       The left-right hwhm are computed relative to local maximum
C       near the mean value.
C******************************************************************

       subroutine ulhfwhm(id,fwhm)

       use ulhglobal
       save
       real(kind=8) :: fwhm
       real(kind=8) :: xmed,sigma,sx
       real(kind=8) :: xup,xdow
       real(kind=8) :: cmax
       real(kind=8) :: hwhmR, hwhmL
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHFWHM: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHFWHM: histogram must be 1-d, ID=',id
         goto 99
        endif

       nx=nx1d(ip) ! no. of channels

       call ulhxmean(id,xmed,sigma,sx) ! find mean value and sigma

       call ulhgeti(id,xmed,imed)      ! channel of mean value
       xdow=xmed-sigma
       call ulhgeti(id,xdow,idow)
       xup= xmed+sigma
       call ulhgeti(id,xup,iup)

c      find max. value near xmed
       cmax=histo1d(ip,imed)
       imax=imed
       do i=idow,iup
        if(histo1d(ip,i) .gt. cmax) then
         cmax=histo1d(ip,i)
         imax=i
        endif
       enddo

        hwhmR=0.d0
        hwhmL=0.d0

        do j=imax,nx ! find right half-width
           if(histo1d(ip,j) .lt. cmax/2.d0) then
             call ulhgetx(id,j,hwhmR) 
             goto 10
           endif
        enddo
10      continue

        do j=imax,1,-1 ! find left half-width
           if(histo1d(ip,j) .lt. cmax/2.d0) then 
             call ulhgetx(id,j,hwhmL)
             goto 20
           endif
        enddo
20      continue

        fwhm=abs(hwhmR-hwhmL)
99      end


c================================================
c      Getting data out of histograms
c================================================


C************************************************
C       Get no. of events and entries
C************************************************
       subroutine ulhgetnev(id,nev,nent)

       use ulhglobal
       save
       integer(kind=4) :: id,nent,nev

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETNEV: histogram does not exist, ID=',id
         goto 99
        endif

        select case(idime)
         case(1)
         nent=nentries1d(ip)
         nev=nev1d(ip)
         case(2)
         nent=nentries2d(ip)
         nev=nev2d(ip)
         case(3)
         nent=nentries3d(ip)
         nev=nev3d(ip)
        end select

99       end


C************************************************
C       Get content and error for a given x
C************************************************
       subroutine ulhgetcx(id,x,c,e)

       use ulhglobal
       save
       real(kind=8) ::x,c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETCX: histogram ID does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHGETCX: histogram must be 1-d, ID=',id
         goto 99
        endif

       i=int((x-xmin1d(ip))/binx1d(ip))+1

        if(i.lt.1 .or. i .gt. nx1d(ip))then
         write(iue,*)'ULHGETCX: channel outside range, ID=',id
         goto 99
        endif

       c=histo1d(ip,i)
       e=histo1dv(ip,i)

99       end

C************************************************
C       Get content and error for a given (x,y)
C************************************************
       subroutine ulhgetcxy(id,x,y,c,e)

       use ulhglobal
       save
       real(kind=8) ::x,y,c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETCXY: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.2)then
         write(iue,*)'ULHGETCXY: histogram must be 2-d, ID=',id
         goto 99
        endif

       i=int((x-xmin2d(ip))/binx2d(ip))+1
       j=int((y-ymin2d(ip))/biny2d(ip))+1

        if(i.lt.1 .or. i .gt. nx2d(ip) .or.
     &      j.lt.1 .or. j .gt. ny2d(ip) )then
         write(iue,*)'ULHGETCXY: channel outside range, ID=',id
         goto 99
        endif

       c=histo2d(ip,i,j)
       e=histo2dv(ip,i,j)


99       end

C************************************************
C       Get content and error for a given (x,y,z)
C************************************************
       subroutine ulhgetcxyz(id,x,y,z,c,e)

       use ulhglobal
       save
       real(kind=8) ::x,y,z,c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETCXYZ: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.3)then
         write(iue,*)'ULHGETCXYZ: histogram must be 3-d, ID=',id
         goto 99
        endif

       i=int((x-xmin3d(ip))/binx3d(ip))+1
       j=int((y-ymin3d(ip))/biny3d(ip))+1
       k=int((z-zmin3d(ip))/binz3d(ip))+1

        if(i.lt.1 .or. i .gt. nx3d(ip) .or.
     &     j.lt.1 .or. j .gt. ny3d(ip) .or. 
     &     k.lt.1 .or. k .gt. nz3d(ip) )then
         write(iue,*)'ULHGETCXYZ: channel outside range, ID=',id
         goto 99
        endif

       c=histo3d(ip,i,j,k)
       e=histo3dv(ip,i,j,k)

99       end

C************************************************
C       Get content and error for given i channel  
C************************************************
       subroutine ulhgetci(id,i,c,e)

       use ulhglobal
       save
       real(kind=8) ::c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETCI: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHGETCI: histogram must be 1-d, ID=',id
         goto 99
        endif

        if(i.lt.1 .or. i .gt. nx1d(ip))then
         write(iue,*)'ULHGETCI: channel outside range, ID=',id
         goto 99
        endif

       c=histo1d(ip,i)
       e=histo1dv(ip,i)

99       end

C*******************************************************
C       Get content and error for a given (i,j) channel 
C*******************************************************
       subroutine ulhgetcij(id,i,j,c,e)

       use ulhglobal
       save
       real(kind=8) ::c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETCIJ: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.2)then
         write(iue,*)'ULHGETCIJ: histogram must be 2-d, ID=',id
         goto 99
        endif

        if(i.lt.1 .or. i .gt. nx2d(ip) .or.
     &      j.lt.1 .or. j .gt. ny2d(ip) )then
         write(iue,*)'ULHGETCIJ: channel outside range, ID=',id
         goto 99
        endif

       c=histo2d(ip,i,j)
       e=histo2dv(ip,i,j)

99       end

C*********************************************************
C       Get content and error for a given (i,j,k) channel 
C*********************************************************
       subroutine ulhgetcijk(id,i,j,k,c,e)

       use ulhglobal
       save
       real(kind=8) ::c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETCIJK: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.3)then
         write(iue,*)'ULHGETCIJK: histogram must be 3-d, ID=',id
         goto 99
        endif

        if(i.lt.1 .or. i .gt. nx3d(ip) .or.
     &     j.lt.1 .or. j .gt. ny3d(ip) .or. 
     &     k.lt.1 .or. k .gt. nz3d(ip) )then
         write(iue,*)'ULHGETCIJK: channel outside range, ID=',id
         goto 99
        endif

       c=histo3d(ip,i,j,k)
       e=histo3dv(ip,i,j,k)

99     end


C************************************************
C       Get channel i for a given x 
C************************************************
       subroutine ulhgeti(id,x,i)

       use ulhglobal
       save
       real(kind=8) :: x
       integer(kind=4) :: id
       
       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETI: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHGETI: histogram must be 1-d, ID=',id
         goto 99
        endif

       nx=nx1d(ip)
       i=int((x-xmin1d(ip))/binx1d(ip))+1
       if(i.lt.1) then
        i=1
        write(iue,*)'ULHGETI: x below histogram limits, ID=',id

       elseif(i.gt. nx) then
         i=nx
         write(iue,*)'ULHGETI: x above histogram limits, ID=',id
       endif

99       end

C************************************************
C       Get (i,j) for a given (x,y)
C************************************************
       subroutine ulhgetij(id,x,y,i,j)

       use ulhglobal
       save
       real(kind=8) :: x,y
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETIJ: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.2)then
         write(iue,*)'ULHGETIJ: histogram must be 2-d, ID=',id
         goto 99
        endif

       nx=nx2d(ip)
       ny=ny2d(ip)
       i=int((x-xmin2d(ip))/binx2d(ip))+1
       j=int((y-ymin2d(ip))/biny2d(ip))+1

       if(i.lt.1) then
        i=1
        write(iue,*)'ULHGETIJ: x below histogram limits, ID=',id

       elseif(i.gt. nx) then
         i=nx
         write(iue,*)'ULHGETIJ: x above histogram limits, ID=',id
       endif

       if(j.lt.1) then
        j=1
        write(iue,*)'ULHGETIJ: y below histogram limits, ID=',id

       elseif(j.gt. ny) then
         i=ny
         write(iue,*)'ULHGETIJ: y above histogram limits, ID=',id
       endif


99       end


C************************************************
C       Get (i,j,k) for a given (x,y,z)
C************************************************
       subroutine ulhgetijk(id,x,y,z,i,j,k)

       use ulhglobal
       save
       real(kind=8) :: x,y,z
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETIJK: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.3)then
         write(iue,*)'ULHGETIJK: histogram must be 3-d, ID=',id
         goto 99
        endif

       nx=nx3d(ip)
       ny=ny3d(ip)
       nz=nz3d(ip)
       i=int((x-xmin3d(ip))/binx3d(ip))+1
       j=int((y-ymin3d(ip))/biny3d(ip))+1
       k=int((z-zmin3d(ip))/binz3d(ip))+1

       if(i.lt.1) then
        i=1
        write(iue,*)'ULHGETIJK: x below histogram limits, ID=',id

       elseif(i.gt. nx) then
         i=nx
         write(iue,*)'ULHGETIJK: x above histogram limits, ID=',id
       endif

       if(j.lt.1) then
        j=1
        write(iue,*)'ULHGETIJK: y below histogram limits, ID=',id

       elseif(j.gt. ny) then
         j=ny
         write(iue,*)'ULHGETIJK: y above histogram limits, ID=',id
       endif

       if(k.lt.1) then
        k=1
        write(iue,*)'ULHGETIJK: z below histogram limits, ID=',id

       elseif(k.gt. nz) then
         k=nz
         write(iue,*)'ULHGETIJK: y above histogram limits, ID=',id
       endif

99       end


C************************************************
C       Get x for a given channel i 
C************************************************
       subroutine ulhgetx(id,i,x)

       use ulhglobal
       save
       real(kind=8) :: x
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETX: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHGETX: histogram must be 1-d, ID=',id
         goto 99
        endif

        x=xmin1d(ip)+i*binx1d(ip)-binx1d(ip)/2.d0

99       end

C************************************************
C       Get (x,y) for a given channel (i,j) 
C************************************************
       subroutine ulhgetxy(id,i,j,x,y)

       use ulhglobal
       save
       real(kind=8) :: x,y
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETXY: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.2)then
         write(iue,*)'ULHGETXY: histogram must be 2-d, ID=',id
         goto 99
        endif

        x=xmin2d(ip)+i*binx2d(ip)-binx2d(ip)/2.d0
        y=ymin2d(ip)+j*biny2d(ip)-biny2d(ip)/2.d0

99       end

C************************************************
C       Get (x,y,z) for a given channel (i,j,k) 
C************************************************
       subroutine ulhgetxyz(id,i,j,k,x,y,z)

       use ulhglobal
       save
       real(kind=8) :: x,y,z
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETXYZ: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.3)then
         write(iue,*)'ULHGETXYZ: histogram must be 3-d, ID=',id
         goto 99
        endif

        x=xmin3d(ip)+i*binx3d(ip)-binx3d(ip)/2.d0
        y=ymin3d(ip)+j*biny3d(ip)-biny3d(ip)/2.d0
        z=zmin3d(ip)+k*biny3d(ip)-biny3d(ip)/2.d0

99       end


C**************************************************
C      Get the channel, x, content value and error 
c      for the histogram max. value
C**************************************************
       subroutine ulhgetmax1(id,imax,xmax,cmax,emax)

       use ulhglobal
       save
       real(kind=8) :: xmax,cmax,emax
       real(kind=8) :: c
       integer(kind=4) :: id,imax,i

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETMAX1: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHGETMAX1: histogram must be 1-d, ID=',id
         goto 99
        endif

        imax=1
        cmax=histo1d(ip,1)
        
        do i=2,nx1d(ip)
           c=histo1d(ip,i)
           if(c.ge.cmax)then
            imax=i
            cmax=c
           endif
        enddo

        emax=histo1dv(ip,imax)
        call ulhgetx(id,imax,xmax)

99       end

C***********************************************************
C      Get the channel (i,j), (x,y), content value and error 
c      for the histogram max. value
C***********************************************************
       subroutine ulhgetmax2(id,imax,jmax,xmax,ymax,cmax,emax)

       use ulhglobal
       save
       real(kind=8) :: xmax,ymax,cmax,emax
       real(kind=8) :: c
       integer(kind=4) :: id,imax,jmax,i,j

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETMAX2: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.2)then
         write(iue,*)'ULHGETMAX2: histogram must be 2-d, ID=',id
         goto 99
        endif

        imax=1
        jmax=1
        cmax=histo2d(ip,1,1)
        
        do i=1,nx2d(ip)
         do j=1,ny2d(ip)
           c=histo2d(ip,i,j)
           if(c.ge.cmax)then
            imax=i
            jmax=j
            cmax=c
           endif
         enddo
        enddo

        emax=histo2dv(ip,imax,jmax)
        call ulhgetxy(id,imax,jmax,xmax,ymax)

99       end

C***************************************************************
C      Get the channel (i,j,k), (x,y,z), content value and error 
c      for the histogram max. value
C***************************************************************
       subroutine ulhgetmax3(id,imax,jmax,kmax,
     &                          xmax,ymax,zmax,cmax,emax)


       use ulhglobal
       save
       real(kind=8) :: xmax,ymax,zmax,cmax,emax
       real(kind=8) :: c
       integer(kind=4) :: id,imax,jmax,kmax,i,j,k

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETMAX3: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.3)then
         write(iue,*)'ULHGETMAX3: histogram must be 3-d, ID=',id
         goto 99
        endif

        imax=1
        jmax=1
        kmax=1
        cmax=histo3d(ip,1,1,1)
        
        do i=1,nx3d(ip)
         do j=1,ny3d(ip)
          do k=1,nz3d(ip)
           c=histo3d(ip,i,j,k)
           if(c.ge.cmax)then
            imax=i
            jmax=j
            kmax=k
            cmax=c
           endif
          enddo
         enddo
        enddo

        emax=histo3dv(ip,imax,jmax,kmax)
        call ulhgetxyz(id,imax,jmax,kmax,xmax,ymax,zmax)

99       end

C**************************************************
C      Get the channel, x, content value and error 
C      for the histogram min. value
C**************************************************
       subroutine ulhgetmin1(id,imin,xmin,cmin,emin)

       use ulhglobal
       save
       real(kind=8) :: xmin,cmin,emin
       real(kind=8) :: c
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETMIN1: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHGETMIN1: histogram must be 1-d, ID=',id
         goto 99
        endif

        imin=1
        cmin=histo1d(ip,1)
        
        do i=2,nx1d(ip)
           c=histo1d(ip,i)
           if(c.ge.cmin)then
            imin=i
            cmin=c
           endif
        enddo

        emin=histo1dv(ip,imin)
        call ulhgetx(id,imin,xmin)

99       end


C***********************************************************
C      Get the channel (i,j), (x,y), content value and error 
C      for the histogram min. value
C***********************************************************
       subroutine ulhgetmin2(id,imin,jmin,xmin,ymin,cmin,emin)

       use ulhglobal
       save
       real(kind=8) :: xmin,ymin,cmin,emin
       real(kind=8) :: c
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETMIN2: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.2)then
         write(iue,*)'ULHGETMIN2: histogram must be 2-d, ID=',id
         goto 99
        endif

        imin=1
        jmin=1
        cmin=histo2d(ip,1,1)
        
        do i=1,nx2d(ip)
         do j=1,ny2d(ip)
           c=histo2d(ip,i,j)
           if(c.ge.cmin)then
            imin=i
            jmin=j
            cmin=c
           endif
         enddo
        enddo

        emin=histo2dv(ip,imin,jmin)
        call ulhgetxy(id,imin,jmin,xmin,ymin)

99       end


C****************************************************************
C      Get the channel (i,j,k), (x,y,z), content value and error 
C      for the histogram min. value
C****************************************************************
       subroutine ulhgetmin3(id,imin,jmin,kmin,
     &                          xmin,ymin,zmin,cmin,emin)

       use ulhglobal
       save
       real(kind=8) :: xmin,ymin,zmin,cmin,emin
       real(kind=8) :: c
       integer(kind=4) :: id,imin,jmin,kmin,i,j,k

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHGETMIN3: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.3)then
         write(iue,*)'ULHGETMIN3: histogram must be 3-d, ID=',id
         goto 99
        endif

        imin=1
        jmin=1
        kmin=1
        cmin=histo3d(ip,1,1,1)
        
        do i=1,nx3d(ip)
         do j=1,ny3d(ip)
          do k=1,nz3d(ip)
           c=histo3d(ip,i,j,k)
           if(c.ge.cmin)then
            imin=i
            jmin=j
            kmin=k
            cmin=c
           endif
          enddo
         enddo
        enddo

        emin=histo3dv(ip,imin,jmin,kmin)
        call ulhgetxyz(id,imin,jmin,kmin,xmin,ymin,zmin)

99       end


C*****************************************************
C       Get for histo ID the pointer IP from database
C*****************************************************
       subroutine ulhgetip(id,ip,idime)

       use ulhglobal    
       save
       integer(kind=4) :: id,ip,idime

       ip=0
       idime=0
       if(id.lt.1 .or. id.gt.nmaxhistos) then 
         write(iue,*)'ULHGETIP: ID outside range, ID=',id       
         goto 99 
       endif

       idime=iddim(id)  !  get histo dimension
       if(idime.eq.1)then
        ip=ids1d_1(id)   !  pointer : if ip=0 error on determination
       elseif(idime.eq.2)then
        ip=ids2d_1(id)
       elseif(idime.eq.3)then
        ip=ids3d_1(id)
       endif

99       end

c================================================
c      Put data into histograms
c================================================


C************************************************
C       Put content and error for given i channel  
C************************************************
       subroutine ulhputci(id,i,c,e)

       use ulhglobal
       save
       real(kind=8) ::c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPUTCI: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHPUTCI: histogram must be 1-d, ID=',id
         goto 99
        endif

        if(i.lt.1 .or. i .gt. nx1d(ip))then
         write(iue,*)'ULHPUTCI: channel outside range, ID=',id
         goto 99
        endif

       histo1d(ip,i)=c
       histo1dv(ip,i)=e
       istatus(id)=3

99       end

C*******************************************************
C       Put content and error for a given (i,j) channel 
C*******************************************************
       subroutine ulhputcij(id,i,j,c,e)

       use ulhglobal
       save
       real(kind=8) ::c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPUTCIJ: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.2)then
         write(iue,*)'ULHPUTCIJ: histogram must be 2-d, ID=',id
         goto 99
        endif

        if(i.lt.1 .or. i .gt. nx2d(ip) .or.
     &     j.lt.1 .or. j .gt. ny2d(ip) )then
         write(iue,*)'ULHPUTCIJ: channel outside range, ID=',id
         goto 99
        endif

       histo2d(ip,i,j)=c
       histo2dv(ip,i,j)=e
       istatus(id)=3

99       end

C*********************************************************
C       Put content and error for a given (i,j,k) channel 
C*********************************************************
       subroutine ulhputcijk(id,i,j,k,c,e)

       use ulhglobal
       save
       real(kind=8) :: c,e
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHPUTCIJK: histogram does not exist, ID=',id
         goto 99
        endif

       if(idime.ne.3)then
         write(iue,*)'ULHPUTCIJK: histogram must be 3-d, ID=',id
         goto 99
        endif

        if(i.lt.1 .or. i .gt. nx3d(ip) .or.
     &     j.lt.1 .or. j .gt. ny3d(ip) .or.
     &     k.lt.1 .or. k .gt. nz3d(ip) )then

         write(iue,*)'ULHPUTCIJK: channel outside range, ID=',id
         goto 99
        endif

       histo3d(ip,i,j,k)=c
       histo3dv(ip,i,j,k)=e
       istatus(id)=3

99       end



C************************************************
C       Put an ID into the database
C************************************************
       subroutine ulhputid(id,ip,idime)

       use ulhglobal
       save
       integer(kind=4) :: id,ip,idime

       ip=0

       if(nids+1.gt.nh1d+nh2d+nh3d) then
        write(iue,*)'ULHPUTID: Too many histograms'
        goto 99
       endif

c       check if id is in the range
       if(id.lt.1 .or. id .gt. nmaxhistos)then
        write(iue,*)'ULHPUTID: histo ID=',id,' outside range!'
        goto 99
       endif

c       check if id exists
c       check in 1-d
       ip=ids1d_1(id)
       if(ip.ne.0)then
        write(iue,*)'ULHPUTID: ID=',id,' exists. Choose another ID.'
        goto 99
       endif
c       check in 2-d
       ip=ids2d_1(id)
       if(ip.ne.0)then
        write(iue,*)'ULHPUTID: ID=',id,' exists. Choose another ID.'
        goto 99
       endif
c       check in 3-d
       ip=ids3d_1(id)
       if(ip.ne.0)then
        write(iue,*)'ULHPUTID: ID=',id,' exists. Choose another ID.'
        goto 99
       endif

       if(idime.eq.1)then

         if((nids1d+1).gt.nh1d) then
          write(iue,*)'ULHPUTID: Too many 1-d histograms'
          goto 99
         endif

         nids=nids+1
         nids1d=nids1d+1    ! 1-d histo
         iddim(id)=1
         ids1d(nids1d)=id   ! ID array
         ids1d_1(id)=nids1d ! inverse ID array
         ip=nids1d

       elseif(idime.eq.2)then        

         if(nids2d.gt.nh2d) then
          write(iue,*)'ULHPUTID: Too many 2-d histograms'
          goto 99
         endif

         nids=nids+1
         nids2d=nids2d+1    ! 2-d histo
         iddim(id)=2
         ids2d(nids2d)=id   ! ID array
         ids2d_1(id)=nids2d ! inverse ID array
         ip=nids2d

       elseif(idime.eq.3)then  

         if(nids3d.gt.nh3d) then
          write(iue,*)'ULHPUTID: Too many 3-d histograms'
          goto 99
         endif

         nids=nids+1
         nids3d=nids3d+1    ! 3-d histo
         iddim(id)=3
         ids3d(nids3d)=id   ! ID array
         ids3d_1(id)=nids3d ! inverse ID array
         ip=nids3d

       endif

99       end
       


C=======================================================
C      Making slices 
C=======================================================


C*************************************************
C      Get a 1d slice from a 2d histogram
C*************************************************

       subroutine ulhslice2(id,id2,islice,cut)
       use ulhglobal
       save
       integer(kind=4) :: id,id2,islice
       real(kind=8) :: cut
       integer(kind=4) :: ip,ip2,idime,idime2,nch
       real(kind=8) :: x,y,xmin,xmax,ymin,ymax
       integer :: i,j
       character(len=64) :: tit

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHSLICE2: histogram does not exist, ID=',id
         return
        endif
        if(idime.ne.2)then
         write(iue,*)'ULHSLICE2: not a 2d histogram, ID=',id
         return
        endif


       call ulhgetip(id2,ip2,idime2) ! check ID2 in database
        if(idime2.ne.0)then
         write(iue,*)'ULHSLICE2: histogram exists, ID=',id2
         return
        endif

       tit=title2d(ip)       
       idtyp(id2)=idtyp(id)    ! histo type =1,2,3,4

       select case (islice)
       case (1)            ! x slice @ y= cut
         nch=nx2d(ip)
         xmin=xmin2d(ip)
         xmax=xmax2d(ip)
         call ulhb1(id2,tit,nch,xmin,xmax)
         call ulhgetip(id2,ip2,idime2)
         x=xmin
         y=cut
         call ulhgetij(id,x,y,i,j) ! get j at cut
         do i=1,nch
            histo1d(ip2,i)=histo2d(ip,i,j)
            histo1dv(ip2,i)=histo2dv(ip,i,j)
            histo1dc(ip2,i)=histo2dc(ip,i,j)
            histo1dp(ip2,i)=histo2dp(ip,i,j)
            nev1d(ip2)=nev1d(ip2)+int(histo2dc(ip,i,j)) ! total no. events in slice
         enddo
       case(2)             ! y slice @ x= cut
         nch=ny2d(ip)
         ymin=ymin2d(ip)
         ymax=ymax2d(ip)
         call ulhb1(id2,tit,nch,ymin,ymax)
         call ulhgetip(id2,ip2,idime2)
         x=cut
         y=ymin
         call ulhgetij(id,x,y,i,j) ! get i at cut
         do j=1,nch
            histo1d(ip2,j)=histo2d(ip,i,j)
            histo1dv(ip2,j)=histo2dv(ip,i,j)
            histo1dc(ip2,j)=histo2dc(ip,i,j)
            histo1dp(ip2,j)=histo2dp(ip,i,j)
            nev1d(ip2)=nev1d(ip2)+int(histo2dc(ip,i,j)) ! total no. events in slice
         enddo
       case default
          write(iue,*)'ULHSLICE2: unknown slice- ID, slice',id2,islice
          return
       end select

       istatus(id2)=istatus(id)

       end

C*************************************************
C      Get a 2d slice from a 3d histogram
C*************************************************

       subroutine ulhslice3(id,id2,islice,cut)
       use ulhglobal
       save
       integer(kind=4) :: id,id2,islice
       real(kind=8) :: cut
       integer(kind=4) :: ip,ip2,idime,idime2,nchx,nchy,nchz
       real(kind=8) :: x,y,z,xmin,xmax,ymin,ymax,zmin,zmax
       integer :: i,j,k
       character(len=64) :: tit

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHSLICE3: histogram does not exist, ID=',id
         return
        endif
        if(idime.ne.3)then
         write(iue,*)'ULHSLICE3: not a 3d histogram, ID=',id
         return
        endif

       call ulhgetip(id2,ip2,idime2) ! check ID2 in database
        if(idime2.ne.0)then
         write(iue,*)'ULHSLICE3: histogram exists, ID=',id2
         return
        endif

       tit=title2d(ip)
       idtyp(id2)=idtyp(id)    ! histo type =1,2,3,4

       select case (islice)
       case (1)            ! x-y slice @ z=cut
         nchx=nx3d(ip)
         xmin=xmin3d(ip)
         xmax=xmax3d(ip)
         nchy=ny3d(ip)
         ymin=ymin3d(ip)
         ymax=ymax3d(ip)
         call ulhb2(id2,tit,nchx,xmin,xmax,nchy,ymin,ymax)
         call ulhgetip(id2,ip2,idime2)
         x=xmin
         y=ymin
         z=cut
         call ulhgetijk(id,x,y,z,i,j,k) ! get k at cut
         do i=1,nchx
          do j=1,nchy
            histo2d(ip2,i,j)=histo3d(ip,i,j,k)
            histo2dv(ip2,i,j)=histo3dv(ip,i,j,k)
            histo2dc(ip2,i,j)=histo3dc(ip,i,j,k)
            histo2dp(ip2,i,j)=histo3dp(ip,i,j,k)
            nev2d(ip2)=nev2d(ip2)+int(histo3dc(ip,i,j,k)) ! total no. events in slice
          enddo
         enddo

       case(2)             ! x-z slice @ y=cut 
         nchx=nx3d(ip)
         xmin=xmin3d(ip)
         xmax=xmax3d(ip)
         nchz=nz3d(ip)
         zmin=zmin3d(ip)
         zmax=zmax3d(ip)
         call ulhb2(id2,tit,nchx,xmin,xmax,nchz,zmin,zmax)
         call ulhgetip(id2,ip2,idime2)
         x=xmin
         z=zmin
         y=cut
         call ulhgetijk(id,x,y,z,i,j,k) ! get j at cut
         do i=1,nchx
          do k=1,nchz
            histo2d(ip2,i,k)=histo3d(ip,i,j,k)
            histo2dv(ip2,i,k)=histo3dv(ip,i,j,k)
            histo2dc(ip2,i,k)=histo3dc(ip,i,j,k)
            histo2dp(ip2,i,k)=histo3dp(ip,i,j,k)
            nev2d(ip2)=nev2d(ip2)+int(histo3dc(ip,i,j,k)) ! total no. events in slice
          enddo
         enddo

       case(3)             ! y-z slice @ x=cut 
         nchy=ny3d(ip)
         ymin=ymin3d(ip)
         ymax=ymax3d(ip)
         nchz=nz3d(ip)
         zmin=zmin3d(ip)
         zmax=zmax3d(ip)
         call ulhb2(id2,tit,nchy,ymin,ymax,nchz,zmin,zmax)
         call ulhgetip(id2,ip2,idime2)
         x=cut
         y=ymin
         z=zmin
         call ulhgetijk(id,x,y,z,i,j,k) ! get i at cut
         do j=1,nchy
          do k=1,nchz
            histo2d(ip2,j,k)=histo3d(ip,i,j,k)
            histo2dv(ip2,j,k)=histo3dv(ip,i,j,k)
            histo2dc(ip2,j,k)=histo3dc(ip,i,j,k)
            histo2dp(ip2,j,k)=histo3dp(ip,i,j,k)
            nev2d(ip2)=nev2d(ip2)+int(histo3dc(ip,i,j,k)) ! total no. events in slice
          enddo
         enddo

       case default
          write(iue,*)'ULHSLICE3: unknown slice- ID, slice',id2,islice
          return
       end select

       istatus(id2)=istatus(id)

       end


C==================================================
C       Random number generators
C==================================================

c     compatibility routine
      subroutine ulhrndi1(id)
       use ulhglobal
       save
       call ulhrndi(id)
      end

C************************************************
C       Initialization routine to the inverse-
C       transform generator. Histogram id is 
C       replaced by its integral
C************************************************

       subroutine ulhrndi(id)

       use ulhglobal
       save
       real(kind=8) :: sumt,sumx
       integer(kind=4) :: id

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHRNDI: histogram does not exist, ID=',id
         goto 99
        elseif(idime.ge.2)then
         write(iue,*)'ULHRNDI:2-d or 3-d histograms not allowed,ID=',id
         goto 99
        endif

       if(istatus(id).eq.5)then
         write(iue,*)'ULHRNDI: histogram already transformed, ID=',id
         goto 99
       endif

       istatus(id)=5    ! set status to generator 

       sumt=0.D0              ! get histogram integral
       do i=1,nx1d(ip)
         if(histo1d(ip,i).lt.0.D0) then
          write(iue,*)
     &     'ULHRNDI: histogram contents must be positive, ID=',id
          goto 99
         endif
         sumt=sumt+histo1d(ip,i)
       enddo


       sumx=0.D0
       do i=1,nx1d(ip)
         sumx=sumx+histo1d(ip,i)/sumt
         histo1d(ip,i)=sumx
       enddo

99       end

C************************************************
C       Generate numbers according to distribution
C       given in histogram id using the inverse-
C       transform method. 
C************************************************
       subroutine ulhrnd1(id,x)

       use ulhglobal
       save
       real(kind=8) :: x,xmin,xmax,f1
       real(kind=8) :: ymin,ymax,ymed
       real(kind=8) :: ulhrand
       integer(kind=4) :: id
       integer :: imin,imax,imed
       integer :: j

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHRND1: histogram does not exist, ID=',id
         goto 99
        endif

        if(istatus(id).ne.5)then
         write(iue,*)'ULHRND1: call ULHRNDI first, ID=',id
         goto 99
        endif

       f1=ulhrand(1.d0)  ! random seek value

       imin=0              ! starting i seeking values
       imax=nx1d(ip)
       imed=(imin+imax)/2
       
c 
1       continue   ! seek channel dividing intervals in two halves

       if(imin.ne.0)then
         ymin=histo1d(ip,imin)
       else                    ! i=0 is a special case
         ymin=0.D0
       endif
       ymax=histo1d(ip,imax)
       ymed=histo1d(ip,imed)

C       comparison chain
       if(f1.ge.ymin .and. f1.le.ymed)then  ! lower half
          if(imed.eq.imin+1) then
           j=imed      ! found channel value
           goto 9      ! end comparison
         endif
        imax=imed
        imed=(imin+imax)/2

       else   ! upper half

        if(imed.eq.imax-1) then
          j=imax        ! found channel value
         goto 9       ! end comparison
        endif
        imin=imed
        imed=(imin+imax)/2

       endif

       goto 1    ! continue comparison
              
9       continue
       xmin=xmin1d(ip)+(j-1)*binx1d(ip) ! channel lower bound
       xmax=xmin+binx1d(ip)             ! channel upper bound
       x=xmax-(xmax-xmin)*ulhrand(2.d0)    ! choose x uniformely 
                                         ! distributed in the channel

C       x=xmin1d(ip)+j*binx1d(ip)-binx1d(ip)/2.d0 ! comput channel center 

99       end

C***************************************************************
C       Initialization routine to the hit-or-miss 2d generator. 
C       Routine find the histogram maximum
C***************************************************************

       subroutine ulhrndi2(id)

       use ulhglobal
       save
       real(kind=8) :: xmax,ymax,cmax,emax
       integer(kind=4) :: id,ip,idime,imax,jmax
       real(kind=8), dimension(20) :: hinfo

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHRNDI2: histogram does not exist, ID=',id
         goto 99
        elseif(idime.ne.2)then
         write(iue,*)'ULHRNDI2: histogram not 2D ID=',id
         goto 99
        endif

        istatus(id)=5    ! set status to generator 

c       get histo max value
        call ulhgetmax2(id,imax,jmax,xmax,ymax,cmax,emax)
        histog2cmax(ip)=cmax

c       get histo bounderies 
        call ulhinfo(id,hinfo)

        histog2binx(ip)=hinfo(1)
        histog2xmin(ip)=hinfo(2)
        histog2xmax(ip)=hinfo(3)

        histog2biny(ip)=hinfo(5)
        histog2ymin(ip)=hinfo(6)
        histog2ymax(ip)=hinfo(7)

99       end

C**********************************************************
C       Generate numbers according to 2d distribution
C       given in histogram id using the hit-or-miss method. 
C**********************************************************
       subroutine ulhrnd2(id,x,y)

       use ulhglobal
       implicit none
       save
       real*8 :: x,xmin,xmax,cmax
       real*8 :: y,ymin,ymax,c,e
       real*8 :: binx,biny,z
       real*8 :: ulhrand
       integer*4 :: id,ip,idime

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHRND2: histogram does not exist, ID=',id
         goto 99
        endif

        if(istatus(id).ne.5)then
         write(iue,*)'ULHRND2: call ULHRNDI2 first, ID=',id
         goto 99
        endif

        xmin=histog2xmin(ip)
        xmax=histog2xmax(ip)
        binx=histog2binx(ip)
        ymin=histog2ymin(ip)
        ymax=histog2ymax(ip)
        biny=histog2biny(ip)
        cmax=histog2cmax(ip)

10      continue

        x=xmin +(xmax-xmin)*ulhrand(1.d0)
        y=ymin +(ymax-ymin)*ulhrand(2.d0)
        z=cmax*ulhrand(3.d0)
        
        call ulhgetcxy(id,x,y,c,e)
        if(z > c) goto 10

99      end


      FUNCTION ULHRAND(DUMMY)
C
C  This is an adapted version of subroutine RANECU written by F. James
C  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
C  give a single random number at each call.
C
C  The 'seeds' ISEED1 and ISEED2 can be initialised in the main program
C  and transferred through the module ULHRANDOM.
C

C  Adapted by L. Peralta @2011 for ULHISTOS
!                        @2016

      USE ULHRANDOM
      real(kind=8) :: ulhrand,dummy,dum
      real(kind=8), PARAMETER :: USCALE=1.0D0/2.147483563D9
      integer(kind=4) :: I1,I2,IZ

C
      dum=dummy
      I1=ISEED1/53668
      ISEED1=40014*(ISEED1-I1*53668)-I1*12211
      IF(ISEED1.LT.0) ISEED1=ISEED1+2147483563
C
      I2=ISEED2/52774
      ISEED2=40692*(ISEED2-I2*52774)-I2*3791
      IF(ISEED2.LT.0) ISEED2=ISEED2+2147483399
C
      IZ=ISEED1-ISEED2
      IF(IZ.LT.1) IZ=IZ+2147483562
      ULHRAND=IZ*USCALE
C
      RETURN
      END

c------------------------------------------------------------------------
c fitting routines
c------------------------------------------------------------------------

c-----------------------------------------------------------
c       fit histogram 1D
c-----------------------------------------------------------
        subroutine ulhfith1
     &  (id,xmin,xmax,ifitype,ierrtype,npoints,m,em,b,eb,r,chi2,ierr)
        
        use ulhglobal
        implicit none
        save
        integer*4 :: id,ifitype,ierrtype,ierr
        integer*4 :: ip,idime,imin,imax,n,npoints,i
        real*8 :: xmin,xmax,m,em,b,eb,chi2,r,ndf
        real*8 :: temp,xn,yn,en
        integer*4, parameter :: nmaxpoints=1000000
        real*8, dimension(nmaxpoints) :: x,y,ey
        integer*4 :: nfit

c       input: id,xmin,xmax,ifitype,ierrtype
c       output: npoints,m,em,b,eb,r,chi2,ierr
c       ifitype 1: y=m*x+b  2: y=b*exp(m*x)   3: y=b*x**m
c
c       ierr: 0: no error, 1: histo doen't exist, 2:histo is not 1D, 
c             3: no.points< 2, 4: no. points> nmaxpoints,  
c             5: unknown type uncertainty, 6: unknown fit type
c            
c       ierrtype - uncertainty: 1: no errors, 2: given values, 3: sqrt(counts)

       ierr=0
       m=0.d0
       em=0.d0
       b=0.d0
       eb=0.d0
       chi2=0.d0

       call ulhgetip(id,ip,idime) ! check ID in database
        if(idime.eq.0)then
         write(iue,*)'ULHFITH1: histogram does not exist, ID=',id
         ierr=1
         return
        endif

       if(idime.ne.1)then
         write(iue,*)'ULHFITH1: histogram must be 1-d, ID=',id
         ierr=2
         return
        endif

        if(xmin > xmax) then ! inverted values
          temp=xmax
          xmax=xmin
          xmin=temp
        endif

        call ulhgeti(id,xmin,imin)
        call ulhgeti(id,xmax,imax)

        npoints=imax-imin+1
        if(npoints < 2 ) then
         ierr=3
         return
        elseif(npoints > nmaxpoints) then
         ierr=4
         return
        endif

c       read data from histogram ------------------------
        n=0
        select case(ifitype)

c       y=b+m*x
        case(1) ! straight line *************

        do i=imin,imax
         n=n+1
         call ulhgetx(id,i,xn)
         x(n)=xn
         call ulhgetci(id,i,yn,en)
         y(n)=yn

         select case(ierrtype)
          case(1)
          case(2)
           ey(n)=en
          case(3)
           ey(n)=sqrt(y(n))
          case default
            ierr=5
            return
         end select
        enddo


         select case(ierrtype)
           case(1)
            CALL ulhrecta1(X,Y,Npoints,m,em,b,eb,R,CHI2)
           case(2,3)
                CALL ulhrectaE(X,Y,EY,Npoints,m,em,b,eb,R,CHI2)
           case default
            ierr=5
            return
         end select

        NDF=Npoints-2
        chi2=chi2/ndf
c        end straight line fit

c       y=b*exp(m*x)
        case(2) ! exponential **************
        do i=imin,imax
         n=n+1
         call ulhgetx(id,i,xn)
         x(n)=xn
         call ulhgetci(id,i,yn,en)
         y(n)=yn
         if(yn <= 0.) then !discard negative or zero
            n=n-1
            goto 10
         endif

         select case(ierrtype)
          case(1)
          case(2)
           ey(n)=en
          case(3)
           ey(n)=sqrt(y(n))
          case default
            ierr=4
            return
         end select

         ey(n)=ey(n)/y(n) ! transform ey and y
         y(n)=log(y(n))

10      continue
        enddo

        nfit=n
        if(nfit < 2 ) then
         ierr=3
         return
        endif

         select case(ierrtype)
           case(1)
            CALL ulhrecta1(X,Y,nfit,m,Em,b,Eb,R,CHI2)
           case(2,3)
            CALL ulhrectaE(X,Y,EY,nfit,m,Em,b,Eb,R,CHI2)
           case default
            ierr=4
            return
         end select

        eb=exp(b)*eb
        b=exp(b)
        NDF=nfit-2
        chi2=chi2/ndf
c       end exponencial fit

c       y=b*x**m
        case(3) ! power **************
        do i=imin,imax
         n=n+1
         call ulhgetx(id,i,xn)
         x(n)=xn
         call ulhgetci(id,i,yn,en)
         y(n)=yn
         if(yn <= 0.) then !discard negative or zero
            n=n-1
            goto 20
         endif

         select case(ierrtype)
          case(1)
          case(2)
           ey(n)=en
          case(3)
           ey(n)=sqrt(y(n))
          case default
            ierr=4
            return
         end select

         ey(n)=ey(n)/y(n) ! transform ey and y
         y(n)=log(y(n))
         x(n)=log(x(n))

20      continue
        enddo

        nfit=n
        if(nfit < 2 ) then
         ierr=3
         return
        endif

         select case(ierrtype)
           case(1)
            CALL ulhrecta1(X,Y,nfit,m,Em,b,Eb,R,CHI2)
           case(2,3)
            CALL ulhrectaE(X,Y,EY,nfit,m,Em,b,Eb,R,CHI2)
           case default
            ierr=4
            return
         end select

        eb=exp(b)*eb
        b=exp(b)
        NDF=nfit-2
        chi2=chi2/ndf
c       end power fit

        case default !*******************
         ierr=6
         return
        end select
        
        end

c-----------------------------------------------------------
c       fit vector (x,y)
c-----------------------------------------------------------
        subroutine ulhfitvec
     &  (x,y,ey,npoints,ifitype,ierrtype,m,em,b,eb,r,chi2,ierr)
        
        use ulhglobal
        implicit none
        save
        integer*4 :: ifitype,ierrtype,ierr
        integer*4 :: n,npoints,i
        real*8 :: m,em,b,eb,chi2,r,ndf
        real*8, dimension(npoints) :: x,y,ey,xfit,yfit,eyfit
        integer*4 :: nfit

c       ifitype 1: y=m*x+b  2: y=b*exp(m*x)   3: y=b*x**m
c       ierrtype 0: no errors  1:ey   2: sqrt(y)
c
c       ierr: 0: no error,
c             3: no.points< 2, 
c             5: unknown error type, 6: unknown fit type
c            
c       ierrtype - uncertainty: 1: no errors, 2: given values, 3: sqrt(counts)

       ierr=0
       m=0.d0
       em=0.d0
       b=0.d0
       eb=0.d0
       chi2=0.d0

        if(npoints < 2 ) then
         ierr=3
         return
        endif

        do i=1,npoints ! reset fit variables
         xfit(i)=0.
         yfit(i)=0.
         eyfit(i)=0.
        enddo

        n=0
        select case(ifitype)

c       y=b+m*x
        case(1) ! straight line *************

        do i=1,npoints
         xfit(i)=x(i)
         yfit(i)=y(i)

         select case(ierrtype)
          case(1)
          case(2)
           eyfit(i)=ey(i)
          case(3)
           eyfit(i)=sqrt(y(i))
          case default
            ierr=5
            return
         end select

        enddo

         select case(ierrtype)
           case(1)
            CALL ulhrecta1(xfit,yfit,Npoints,m,em,b,eb,R,CHI2)
           case(2,3)
                CALL ulhrectaE(xfit,yfit,eyfit,Npoints,m,em,b,eb,R,CHI2)
           case default
            ierr=5
            return
         end select

        NDF=Npoints-2
        if(ndf >0) then
         chi2=chi2/ndf
        else
         chi2=0.
        endif
c        end straight line fit

c       y=b*exp(m*x)
        case(2) ! exponential **************
        n=0
        do i=1,npoints
        n=n+1
         if(y(i) <= 0.) then !discard negative or zero
            n=n-1
            goto 10
         endif

         select case(ierrtype)
          case(1)
          case(2)
           eyfit(n)=ey(i)
          case(3)
           eyfit(n)=sqrt(y(i))
          case default
            ierr=5
            return
         end select

         eyfit(n)=eyfit(n)/y(i) ! transform ey and y
         yfit(n)=log(y(i))
         xfit(n)=x(i)

10      continue
        enddo

        nfit=n
        if(nfit < 2 ) then
         ierr=3
         return
        endif

         select case(ierrtype)
           case(1)
            CALL ulhrecta1(Xfit,Yfit,nfit,m,Em,b,Eb,R,CHI2)
           case(2,3)
            CALL ulhrectaE(Xfit,Yfit,EYfit,nfit,m,Em,b,Eb,R,CHI2)
           case default
            ierr=5
            return
         end select

        eb=exp(b)*eb
        b=exp(b)
        NDF=nfit-2
        if(ndf >0) then
         chi2=chi2/ndf
        else
         chi2=0.
        endif
c       end exponencial fit

c       y=b*x**m
        case(3) ! power **************

        n=0
        do i=1,npoints
         n=n+1
         if(y(i) <= 0.) then !discard negative or zero
            n=n-1
            goto 20
         endif

         select case(ierrtype)
          case(1)
          case(2)
           eyfit(n)=ey(i)
          case(3)
           eyfit(n)=sqrt(y(i))
          case default
            ierr=5
            return
         end select

         eyfit(n)=eyfit(n)/y(i) ! transform ey and y
         yfit(n)=log(y(i))
         xfit(n)=log(x(i))

20      continue
        enddo
    
        nfit=n   
        if(nfit < 2 ) then
         ierr=3
         return
        endif

         select case(ierrtype)
           case(1)
            CALL ulhrecta1(Xfit,Yfit,nfit,m,Em,b,Eb,R,CHI2)
           case(2,3)
            CALL ulhrectaE(Xfit,Yfit,EYfit,nfit,m,Em,b,Eb,R,CHI2)
           case default
            ierr=5
            return
         end select

        eb=exp(b)*eb
        b=exp(b)
        NDF=Nfit-2
        if(ndf >0) then
         chi2=chi2/ndf
        else
         chi2=0.
        endif
c       end power fit

        case default !*******************
         ierr=6
         return
        end select
        
        end

c----------------------------------------------------------
        SUBROUTINE ulhrecta1(X,Y,N,RM,SM,B,SB,R,CHI2)
c----------------------------------------------------------
c        recta de regressao linear, sem erros experimentais

        implicit none
        REAL*8, DIMENSION(N) :: X,Y,W
        real*8 :: RM,SM,B,SB,R,CHI2
        integer*4 :: N,i
        real*8 :: sx,sy,sx2,sy2,sxy,delta,sd,alfa
        integer*4 :: k

        IF ( N .LT. 2 ) STOP

        SX=0.
        SY=0.
        SX2=0.
        SXY=0.
        SY2=0.

        DO 200 K=1,N
        SX=SX+X(K)
        SY=SY+Y(K)
        SX2=SX2+X(K)*X(K)
        SXY=SXY+X(K)*Y(K)
        SY2=SY2+Y(K)*Y(K)
200        CONTINUE
C
C        Calculo do declive M e ordenada na origem B  
C
        RM=(N*SXY-SX*SY)/(N*SX2-SX*SX)
        B=(SY*SX2-SX*SXY)/(N*SX2-SX*SX)
C
C        Calculo dos erros SM e SB
C
        IF ( N .EQ. 2 ) THEN
        SM=0.
        SB=0.
        GOTO 1000
        END IF
        SD=0.
        DO 300 K=1,N
        W(K)=RM*X(K)+B-Y(K)
        SD=SD+W(K)*W(K)
300        CONTINUE
        DELTA=N*SX2-SX*SX
        ALFA=SD/(N-2.)
        SM=SQRT(ALFA*N/DELTA)   ! valor corrigido em 2017
        SB=SQRT(ALFA*SX2/DELTA) ! valor corrigido em 2017

1000        CONTINUE
        R=(N*SXY-SX*SY)/SQRT((N*SX2-SX*SX)*(N*SY2-SY*SY))

c        calculo do chi2
        CHI2=0.
        IF(N .LE. 2) GOTO 9999
        DO I=1,N
        CHI2=CHI2+(Y(I)-RM*X(I)-B)**2./ABS(RM*X(I)+B)
        ENDDO

9999        RETURN
        END



c-------------------------------------------------------------------
        SUBROUTINE ulhrectaE (X,Y,SIGY,N,XM,SIGM,B,SIGB,R,CHI2)
c-------------------------------------------------------------------
c        Ajuste de uma recta com erros experimentais
C 
        implicit none
        REAL*8, DIMENSION(N) :: X,Y,SIGY
        real*8 :: XM,SIGM,B,SIGB,R,CHI2
        integer*4 :: N
        real*8 :: sx,sy,sxx,syy,sxy,sinv,s2xy
        real*8 :: xmed,ymed,delta,rsx,rsy
        integer*4 :: i,ierr

C
C
        XM=0.
        SIGM=0.
        B=0.
        SIGB=0.
        CHI2=0.
        R=0.


        SX=0.
        SXX=0.
        SY=0.
        SXY=0.
        SYY=0.
        SINV=0.

        DO I=1,N
         if( abs(SIGY(I)) > 0.) then
           SX=SX+X(I)/(SIGY(I)*SIGY(I))
           SXX=SXX+X(I)*X(I)/(SIGY(I)*SIGY(I))
           SY=SY+Y(I)/(SIGY(I)*SIGY(I))
           SXY=SXY+X(I)*Y(I)/(SIGY(I)*SIGY(I))
           SYY=SYY+Y(I)*Y(I)/(SIGY(I)*SIGY(I))
           SINV=SINV+1./(SIGY(I)*SIGY(I))
          endif
         enddo    

C
C        CALCULO DE DELTA
C 
        DELTA=SINV*SXX-SX*SX
c        todos os dados sao zero
        IF(DELTA .le. 0.)then
        ierr=1
        GOTO 999
        endif
C
C        PARAMETROS DA RECTA     Y=XM*X+B
C        CALCULO DE XM
C 
        XM=(SINV*SXY-SX*SY)/DELTA
C
C        CALCULO DE   B
C 
        B=(SXX*SY-SX*SXY)/DELTA
C
C        SIGMA B
C 
        SIGB=SQRT(SXX/DELTA)
C
C        SIGMA M
C 
        SIGM=SQRT(SINV/DELTA)
C
C        COEFICIENTE DE CORRELACAO LINEAR
C 
        XMED=SX/SINV
        YMED=SY/SINV
        S2XY=1./(N-1.)*(SXY/(SINV/N)-N*XMED*YMED)
        RSX=SQRT(1./(N-1.)*(SXX/(SINV/N)-N*XMED*XMED))
        RSY=SQRT(1./(N-1.)*(SYY/(SINV/N)-N*YMED*YMED))
        R=S2XY/(RSX*RSY)

*        calculo do chi2
        CHI2=0.
        IF(N .LE. 2) GOTO 999
        DO I=1,N
        IF(abs(SIGY(I)) > 0.) THEN
        CHI2=CHI2+(Y(I)-XM*X(I)-B)**2./SIGY(I)**2.
        ENDIF
        ENDDO

999        RETURN
        END
