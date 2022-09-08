C**********************************************************************
C       ULYSSES                                                       *
C                                                                     *
C       A geometry package for radiation transport                    *
C                                                                     *
C       Luis Peralta & Ana Farinha                                    *
C       2021    Version 5.7                                           *
C                                                                     *
C       Copyright (C) 2007-2021                                       *
C       Laboratorio de Instrumentacao e  Particulas                   *
C       Universidade de Lisboa                                        *
C                                                                     *
C**********************************************************************
C    Copyright (C) 2008  Ana Farinha, Luis Peralta
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

c    Modified to include optical components

c************************************************************************************
c      MODULES for Ulysses 5.x     
c************************************************************************************

       MODULE ullib_mod

       real*8, parameter :: pi=3.141592653589793
       character*5, parameter :: version_ullib='5.8' ! version

       integer*4, parameter :: nvmax=10000    ! max no. volumes
       integer*4, parameter :: nchildmax=100  ! max no. child volumes
       integer*4, parameter :: nparadim=6     ! max no. of volume parameters        
       real*8,    parameter :: delta_across=1.d-8  !  increment across border

c  init_arrays_flag ! initialization flag
c  nvolu            ! no. of defined volumes
       integer*4 :: init_arrays_flag
       integer*4 :: nvolu

c  mychild(nvmax,0:nchildmax) : id of child volumes. mychild(ivol,0)=no. childs
c  mytype(nvmax) : ivol type
c  myrot(nvmax)  : ivol kind of rotation 0-no rot, 1-xx, 2-yy, 3-zz
c  ivolmat(nvmax): ivol material
       integer*4, dimension(nvmax,0:nchildmax) :: mychild
       integer*4, dimension(nvmax) :: mytype,myrot,ivolmat

c  volcenter(nvmax,3): volume reference point
c  volpar(nvmax,nparadim): volume parameters
c  volangle(nvmax) : volume rotation angle
c  volcos(nvmax)   : cos rot angle
c  volsin(nvmax)   : sin rot angle

       real*8, dimension(nvmax,3) :: volcenter
       real*8, dimension(nvmax,nparadim) :: volpar
       real*8, dimension(nvmax) :: volangle,volcos,volsin

c  idall(serial no.) : id of volume 
c  ivolume(id)       : vol. serial no.
c  isearchflag(vol. serial no.): flag is set if vol. has been searched 

       integer*4, dimension(nvmax) :: idall,ivolume,isearchflag

c  involume : flag=true if particle inside volume
c  ivolcheck: volume to be checked for a hit
       logical   :: involume
       integer*4 :: ivolcheck

       integer*4, parameter :: nspc_ger=100     ! max number of generators
       integer*4, parameter :: nspc_ch=10000    ! max number of channels in a spectrum

c  spcX(nspc_ger,nspc_ch)= energy values
c  spcProb(nspc_ger,nspc_ch)= probability values
c  nspc(nspc_ger)= number of energy values
c  nspc_existing(id)= 1 if generator id exists
c  common for generator ulspc

       real*8, dimension(nspc_ger,nspc_ch) :: spcX,spcProb
       integer*4, dimension(nspc_ger) :: nspc,nspc_existing

c virtual volumes
c nvivolu= no. virtual volumes
c  idvirtual(serial no.) : id of volume 
c  ivivolume(id)       : vol. serial no.
c  myvitype(nvmax) : ivol type

       integer*4 :: nvivolu
       integer*4, dimension(nvmax) :: idvirtual,ivivolume,myvitype


       real*8, dimension(3):: xhitpoint  ! exact impact point
       integer*4, dimension(nvmax) :: ivolhitface ! volume id face hit point face. Used by ULOPTICS
       integer*4, dimension(10) :: ihitemp ! temp variable containing hit faces
       integer*4 :: ihitface ! temp variable containnig hit face

       SAVE
       END MODULE ullib_mod

!     random generator seeds
      MODULE ULRSEED
      integer(kind=4) :: iseed1,iseed2
      END MODULE ULRSEED

c***************************************************************************************

        subroutine ulversion_check
        
        use ullib_mod
        implicit double precision(a-h,o-z), integer*4(i-n)
        save

        write(*,*)' *** ULYSSES version: ',version_ullib,' ***'
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Initialization routines 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c************************************************************
c     Luis Peralta                                          *
c     2011                                                  *
c************************************************************
c Init common variables

       subroutine ulinit
       use ullib_mod
       use ulrseed
       implicit double precision(a-h,o-z), integer*4(i-n)
       save

       ISEED1=123456789 ! random generator seeds
       ISEED2=987654321

       dummy=0.
       call ulversion_check

       init_arrays_flag=1
       nvolu=0

c       geometry
              mytype=0
              myrot=0
              ivolume=0

              isearchflag=0

              volangle=0.
              volcos=1.
              volsin=0.

              mychild=0
              volcenter=0.
              volpar=0.

c       auxiliar
       involume=.false.
       ivolcheck=1

       idall=0
       ivolmat=0
       
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c    Input and setting routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c************************************************************
c     Ana Farinha                                           *
c     Luis Peralta                                          *
c     2006                                                  *
c************************************************************
c Positionning routine. Set volumes in space

      subroutine ulposi(id,xcenter)
      use ullib_mod
      implicit double precision(a-h,o-z), integer*4(i-n)    
      real*8, dimension(3) :: xcenter
      save

c     check if arrays are reset
      if(init_arrays_flag.ne.1)Stop 'ULPOSI: Call ULINIT first' 

      itype=mytype(id)
      if(itype.eq.0) then
       write(*,*)'ULPOSI: Unknown volume id:',id
       stop
      endif

      volcenter(id,1)=xcenter(1)
      volcenter(id,2)=xcenter(2)
      volcenter(id,3)=xcenter(3)

      end

c************************************************************
c     Luis Peralta                                          *
c     Dec 2011                                              *
c************************************************************
c     set rotation angle and rotation axis (irot) for volume ivol
c     angles are input in DEG
      subroutine ulrotate(ivol,irot,angle)
      use ullib_mod
      implicit double precision(a-h,o-z), integer*4(i-n)
      save
     
c     check if arrays are reset
      if(init_arrays_flag.ne.1)Stop 'ULROTATE: Call ULINIT first' 

      itype=mytype(ivol)
      if(itype.eq.0) then
       write(*,*)'ULROTATE: Unknown volume id:',ivol
       stop
      endif

      if(irot.eq.0 .or. irot.eq.1 .or. irot.eq.2 .or. irot.eq.3)then

        volangle(ivol)= angle

        theta=angle*pi/180.
        volcos(ivol) = cos(theta)
        volsin(ivol) = sin(theta)

        myrot(ivol)=irot     ! re-write rotation type

      else

       write(*,*)'ULROTATE: unknown rotation type:',irot
       stop
      endif

      end

c***************************************************
c       Luis Peralta                               *
c                                                  *
c       April. 2012                                *
c***************************************************
c     Input a new volume

       subroutine ulvolume(ivol,parin,itype,imother,mat)
       use ullib_mod
       implicit double precision(a-h,o-z), integer*4(i-n)
       real*8, dimension(*) :: parin
       data init/0/
       save
       real*8, dimension(nparadim) :: par

c     check if arrays are reset
      if(init_arrays_flag.ne.1)Stop 'ULVOLUME: Call ULINIT first' 

c     check if first volume is the universe
      if(init.eq.0)then
       if(ivol.ne.1) stop 'ULVOLUME: Define universe first'
       if(mat.ne.0) stop 'ULVOLUME: Universe matter must be zero'
       init=1
      endif

c     check if universe is a box
      if(ivol .eq. 1) then
       if(itype .ne. 100) stop 'ULVOLUME: Universe must be a box'
      endif

c     check if mother is no zero
      if(ivol.ne.1)then
       if(imother.le.0) then
        write(*,*)'ULVOLUME: mother zero or negative for volume:',ivol
        stop
       endif
      endif
      
c check if id already used
c nvolu: no. of existing volumes  

      inuse=0
      do i=1, nvolu
         if(ivol.eq.idall(i))then            
            write(*,*)'ULVOLUME: ID already in use:',ivol
            inuse=1
         endif
         if((ivol.le.0).or.(ivol.gt.nvmax))then
            write(*,*)'ULVOLUME: Invalid ID:',ivol  !modificar
            stop
         endif
      enddo


       call ulgetnpar(itype,npar) 
       do i=1,npar
        par(i)=parin(i)   ! move parameters to local array
       enddo

        if(inuse.eq.0)then
         nvolu=nvolu+1 !
         idall(nvolu)=ivol ! array with all volume ids 
         ivolume(ivol)=nvolu ! position of ivol in idall array 
          
         mytype(ivol)=itype

       call ulvcheck(ivol,par,itype) ! check parameter consistency

      do i=1, npar
         volpar(ivol,i)=par(i) 
      enddo

        ivolmat(ivol)=mat  ! material

c   set volume as a child of imother
       if(imother.gt.0)then
          nchild=mychild(imother,0) ! no. of childs
          nchild=nchild+1           ! add child
          mychild(imother,0)=nchild 
          mychild(imother,nchild)=ivol
       endif
      else    ! volume already defined. Add another mother 
       if(imother.gt.0)then
         nchild=mychild(imother,0) ! 
         nchild=nchild+1           ! 
         mychild(imother,0)=nchild 
         mychild(imother,nchild)=ivol
       endif

      endif
      end

c***************************************************
c       Luis Peralta                               *
c       Jan. 2011                                  *
c       May 2020                                   *
c***************************************************
c      for each volume type performe tests and checks on 
c      input volume parameters
       subroutine ulvcheck(ivol,par,itype)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)
       save
       real*8, dimension(*) :: par

        select case(itype)

        case(100) ! box
c        par(1)=Lx, par(2)=Ly, par(3)=Lz
         do i=1,3
            par(i)=abs(par(i)) ! set to positive values
         enddo
c-----------------------------------------------------------------------------
        case(200,201,202,203) ! cylinder
c       par(1)=r, par(2)=h   
       do i=1,2
         par(i)=abs(par(i)) ! set to positive values
       enddo
c------------------------------------------------------------------------------
        case(210,211,212,213)! tube
c       par(1)=rmin,  par(2)=rmax, par(3)=height
       do i=1,3
         par(i)=abs(par(i)) ! set to positive values
       enddo

         if(par(1).gt.par(2))then
           write(*,*)'ERROR: Rmin > Rmax  Tube id=',ivol
           stop
         endif
c------------------------------------------------------------------------------
        case(220,221,222,223) ! elliptic cylinder
c      par(1)=a, par(2)=b, par(3)=h
       do i=1,3
         par(i)=abs(par(i)) ! set to positive values
       enddo
       if(par(1).eq. 0. .or.par(2).eq. 0.) then
        write(*,*)'ERROR: elliptic cylinder with zero radius id=',ivol
       endif
c------------------------------------------------------------------------------
        case(300)! sphere
c        par(1)=r
         par(1)=abs(par(1)) ! set to positive values

c--------------------------------------------------------------------------------
        case(310,311,312,313,-310,-311,-312,-313)! cut sphere
c       par(1)=r, par(2)=hmin, par(3)=hmax    (zz axis)
       do i=1,3
         par(i)=abs(par(i)) ! set to positive values
       enddo

        if(par(2) .gt. par(3) )then
           write(*,*)'ERROR: hmin > hmax Cut Sphere id=',ivol
           write(*,*)par(2),par(3)
           stop
         endif
         if(abs(par(3)) .gt. par(1) )then
           write(*,*)'ERROR: Abs(hmax) >R Cut Sphere id=',ivol
           stop
         endif
         if(abs(par(2)) .gt. par(1) )then
           write(*,*)'ERROR: Abs(hmin) >R Cut Sphere id=',ivol
           stop
         endif
c--------------------------------------------------------------------------------
        case(400,401,402,403,-400,-401,-402,-403)! cone
c       par(1)=Rmax, par(2)=height of truncated cone = h2 
c       par(3)=h (total)
       do i=1,3
         par(i)=abs(par(i)) ! set to positive values
       enddo

         if(par(2).gt.par(3))then
           write(*,*)'ERROR: h2 > h  Cone id=',ivol
           stop
         endif
C----------------------------------------------------------------------------------
        case(500,501,502,503,-500,-501,-502,-503)! pyramid
c       par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
       do i=1,4
         par(i)=abs(par(i)) ! set to positive values
       enddo
         if(par(3).gt.par(4))then
           write(*,*) 'ERROR: hmin > hmax id=',ivol
           stop
         endif
c-----------------------------------------------------------------------------------
        case(600,601,602,603,-600,-601,-602,-603)! wedge
c       par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
       do i=1,4
         par(i)=abs(par(i)) ! set to positive values
       enddo
         if(par(3).ge.par(4))then
           write(*,*) 'ERROR: hmin >= hmax id=',ivol
           stop
         endif
c-----------------------------------------------------------------------------------
        case(700,701,702,703,-700,-701,-702,-703)! paraboloid
c       par(1)=a, par(2)=b, par(3)=zmax
       do i=1,3
         par(i)=abs(par(i)) ! set to positive values
       enddo
c-----------------------------------------------------------------------------------
        case(800)! ellipsoid
c       par(1)=a, par(2)=b, par(3)=c
       do i=1,3
         par(i)=abs(par(i)) ! set to positive values
       enddo

c-----------------------------------------------------------------------------------
        case default
          write(*,*)'Unknown volume type:',itype
          stop

        end select
        end

c***************************************************
c       Luis Peralta                               *
c       Mar. 2012                                  *
c       May 2020                                   *
c***************************************************
c      get the number of parameters of volume itype
       subroutine ulgetnpar(itype,npar) 
       implicit double precision(a-h,o-z), integer*4(i-n)

        select case(itype)
        case(100) ! box
c       par(1)=Lx, par(2)=Ly, par(3)=Lz
         npar=3
c--------------------------------------------------------------
        case(200,201,202,203) ! cylinder
c       par(1)=r, par(2)=h   
        npar=2
c--------------------------------------------------------------
        case(210,211,212,213)! tube
c       par(1)=rmin,  par(2)=rmax, par(3)=height
        npar=3
c--------------------------------------------------------------
        case(220,221,222,223) ! elliptic cylinder
c      par(1)=a, par(2)=b, par(3)=h
       npar=3
c--------------------------------------------------------------
        case(300)! sphere
c       par(1)=r
        npar=1
c--------------------------------------------------------------
        case(310,311,312,313,-310,-311,-312,-313)! cut sphere
c       par(1)=r, par(2)=hmin, par(3)=hmax    (zz axis)
        npar=3
c--------------------------------------------------------------
        case(400,401,402,403,-400,-401,-402,-403)! cone
c       par(1)=Rmax, par(2)=height of truncated cone = h2 
c       par(3)=h (total)
        npar=3
c--------------------------------------------------------------
        case(500,501,502,503,-500,-501,-502,-503)! pyramid
c       par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
        npar=4
c--------------------------------------------------------------
        case(600,601,602,603,-600,-601,-602,-603)! wedge
c       par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
        npar=4
c--------------------------------------------------------------
        case(700,701,702,703,-700,-701,-702,-703)! paraboloid
c       par(1)=a, par(2)=b, par(3)=zmax
        npar=3
c--------------------------------------------------------------
        case(800)! ellipsoid
c       par(1)=a, par(2)=b, par(3)=c
        npar=3

        case default
          write(*,*)'Unknown volume type:',itype
          stop

        end select
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Translaction and rotations
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c------------------------------------------
c        Luis Peralta
c        Dec. 2011  
c------------------------------------------
         subroutine ulrotxz(xpar)
c        rotate x <-> z  coordinates
         implicit double precision(a-h,o-z), integer*4(i-n)
         real*8, dimension(3) :: xpar

         x=xpar(1)
         y=xpar(2)
         z=xpar(3)
         xpar(1)=-z
         xpar(2)=y
         xpar(3)=x

         end
c------------------------------------------
c        Luis Peralta
c        Dec. 2011  
c------------------------------------------
         subroutine ulrotxz_1(xpar)
c        rotate z <-> x  coordinates
         implicit double precision(a-h,o-z), integer*4(i-n)
         real*8, dimension(3) :: xpar

         x=xpar(1)
         y=xpar(2)
         z=xpar(3)
         xpar(1)=z
         xpar(2)=y
         xpar(3)=-x

         end
c------------------------------------------
c        Luis Peralta
c        Dec. 2011  
c------------------------------------------
         subroutine ulrotyz(xpar)
c        rotate y <-> z  coordinates
         implicit double precision(a-h,o-z), integer*4(i-n)
         real*8, dimension(3) :: xpar

         x=xpar(1)
         y=xpar(2)
         z=xpar(3)
         xpar(1)=x
         xpar(2)=-z
         xpar(3)=y

         end
c------------------------------------------
c        Luis Peralta
c        Dec. 2011  
c------------------------------------------
         subroutine ulrotyz_1(xpar)
c        rotate z <-> y  coordinates
         implicit double precision(a-h,o-z), integer*4(i-n)
         real*8, dimension(3) :: xpar

         x=xpar(1)
         y=xpar(2)
         z=xpar(3)
         xpar(1)=x
         xpar(2)=z
         xpar(3)=-y

         end


c------------------------------------------
c        Luis Peralta
c        May 2020  
c------------------------------------------
         subroutine ulflipx(xpar)
c        x <-> -x  
         implicit double precision(a-h,o-z), integer*4(i-n)
         real*8, dimension(3) :: xpar

         x=-xpar(1)
         xpar(1)=x

         end


c------------------------------------------
c        Luis Peralta
c        May 2020  
c------------------------------------------
         subroutine ulflipy(xpar)
c        y <-> -y  
         implicit double precision(a-h,o-z), integer*4(i-n)
         real*8, dimension(3) :: xpar

         y=-xpar(2)
         xpar(2)=y

         end

c------------------------------------------
c        Luis Peralta
c        May 2020  
c------------------------------------------
         subroutine ulflipz(xpar)
c        z <-> -z  
         implicit double precision(a-h,o-z), integer*4(i-n)
         real*8, dimension(3) :: xpar

         z=-xpar(3)
         xpar(3)=z

         end

c************************************************************
c     Ana Farinha                                           *
c     Luis Peralta                                          *
c     2006, May 2020                                                  *
c************************************************************
c       This rotine makes inverse of the rotation of a volume

      subroutine ulrot_1(ivol,xpar,xrot)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xpar, xrot

      irot=myrot(ivol)
      cost=volcos(ivol)
      sint=volsin(ivol)

      x=xpar(1)
      y=xpar(2)
      z=xpar(3)
      select case(irot)

       case(1)
       xrot(1)=x
       xrot(2)=cost*y - sint*z
       xrot(3)=sint*y + cost*z

       case(2)
       xrot(1)=cost*x - sint*z
       xrot(2)=y
       xrot(3)=sint*x + cost*z

       case(3)
       xrot(1)=cost*x - sint*y
       xrot(2)=sint*x + cost*y
       xrot(3)=z

       case default
        write(*,*) 'Unknown rotation type'
        stop
      end select
       
       end
c************************************************************
c     Luis Peralta                                          *
c     Dec. 2011, May 2020                                             *
c************************************************************
c       This rotine makes the rotation of a volume

      subroutine ulrot(ivol,xpar,xrot)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xpar, xrot

      irot=myrot(ivol)
      cost=volcos(ivol)
      sint=volsin(ivol)

      x=xpar(1)
      y=xpar(2)
      z=xpar(3)

      select case(irot)

       case(1)
       xrot(1)=x
       xrot(2)=cost*y + sint*z
       xrot(3)=-sint*y + cost*z

       case(2)
       xrot(1)=cost*x + sint*z
       xrot(2)=y
       xrot(3)=-sint*x + cost*z

       case(3)
       xrot(1)=cost*x + sint*y
       xrot(2)=-sint*x + cost*y
       xrot(3)=z

       case default
        write(*,*) 'Unknown rotation type'
        stop
      end select
  
       end

c************************************************************
c     Luis Peralta                                          *
c     Nov. 2013, May 2020                                   *
c************************************************************
c     Rotate a point xpar about axis irot by theta
c     Rotation axis is at (0,0,0)

      subroutine ulrotp(irot,theta,xpar,xrot)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xpar, xrot

      !pi=acos(-1.d0)
      cost=cos(theta*pi/180.d0)
      sint=sin(theta*pi/180.d0)

      x=xpar(1)
      y=xpar(2)
      z=xpar(3)

      select case(irot)

       case(1)
       xrot(1)=x
       xrot(2)=cost*y - sint*z
       xrot(3)=sint*y + cost*z

       case(2)
       xrot(1)=cost*x - sint*z
       xrot(2)=y
       xrot(3)=sint*x + cost*z

       case(3)
       xrot(1)=cost*x - sint*y
       xrot(2)=sint*x + cost*y
       xrot(3)=z

       case default
        write(*,*) 'Unknown rotation type'
        stop
      end select
       
       end

c************************************************************
c     Ana Farinha                                           *
c     Luis Peralta                                          *
c     2006                                                  *
c************************************************************
c Inverse translaction

      subroutine ultrans_1(xin,xout,id)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xin, xout

      xout(1)=xin(1)+volcenter(id,1)
      xout(2)=xin(2)+volcenter(id,2)
      xout(3)=xin(3)+volcenter(id,3)    

      end

c************************************************************
c     Ana Farinha                                           *
c     Luis Peralta                                          *
c     2006                                                  *
c************************************************************
c Translaction

      subroutine ultrans(xin,xout,id)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xin, xout

      xout(1)=xin(1)-volcenter(id,1)
      xout(2)=xin(2)-volcenter(id,2)
      xout(3)=xin(3)-volcenter(id,3)

      end

c************************************************************
c     Luis Peralta                                          *
c     Dec. 2013                                                  *
c************************************************************
c Translaction of a point x to a new position x'

      subroutine ultransp(xin,xtr,xout)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xin, xout,xtr

      xout(1)=xin(1)+xtr(1)
      xout(2)=xin(2)+xtr(2)
      xout(3)=xin(3)+xtr(3)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Find the particle 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c***************************************************************
c       Luis Peralta                                           *
c       Nov. 2011                                              *
c***************************************************************
c       finds volume where particle is                         
c       searching all volumes                                  
      subroutine ulfind(x,y,z,ivnew,mat)      
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xpar
     
       do i=1,nvolu
        isearchflag(i)=0   ! reset search flags
       enddo

       xpar(1)=x
       xpar(2)=y
       xpar(3)=z

c       look for particle in all volumes
          do i=nvolu,1,-1
           ivol=idall(i)
             call ulvolsearch(xpar,ivol,i_flag)
             isearchflag(i)=1

             if(i_flag.eq.1)then
                ivnew=ivol
                goto 999
             endif
          enddo
        
         ivnew=1 ! Cannot find particle. Go to universe
         mat=0
         goto 9999 ! exit
 999        continue     
         
c        check if particle is inside one of the child volumes
         i_flag=0    ! reset flag
         ivaux=ivnew
 99      nchild=mychild(ivaux,0) ! no. of childs
         do i=1,nchild
            ivchild=mychild(ivaux,i) ! child volume id
            j=ivolume(ivchild)       ! volume serial number

              if(isearchflag(j).eq. 0)then ! volume has not been searched
                call ulvolsearch(xpar,ivchild,i_flag)
               isearchflag(j)=1    ! search flag set on
              endif

            if(i_flag.eq.1)then
               ivaux=ivchild
               goto 99  ! look into grand-son volume ....
            endif
        
         enddo

         ivnew=ivaux
         mat=ivolmat(ivnew)

9999    end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c    Tracking
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c************************************************************
c       Luis Peralta                                        *
c       Nov. 2011, Sep. 2021                                           *
c************************************************************
c     Tracking routine
      subroutine ultrack(x,y,z,u,v,w,ivold,step,ivnew,mat)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: ximp,xpar,upar

       xpar(1)=x   ! particle starting postition
       xpar(2)=y
       xpar(3)=z

       upar(1)=u  ! particle direction
       upar(2)=v
       upar(3)=w  

       ivolhitface=0 ! reset hit face vector for each volume

         call uldhit(xpar,upar,ivold,ximp,dmin) ! find impact point
         xhitpoint=ximp ! copy hit point array and store it in ullib_mod

      if(dmin .ge. step)then  ! particle stays in the same volume
        x=x+u*step
        y=y+v*step
        z=z+w*step

        xpar(1)=x   ! particle starting postition
        xpar(2)=y
        xpar(3)=z 
        call ulvolsearch(xpar,ivold,iflag) 
        if( iflag == 1 ) then
         ivnew=ivold ! particle stays in volume as assumed
         mat=ivolmat(ivnew)
        else      
         call ulfind (x,y,z,ivnew,mat) ! something went wrong. Find particle
        endif 

              
      else

c     particle gets out of volume

        x=ximp(1)+delta_across*u ! put particle across the border
        y=ximp(2)+delta_across*v
        z=ximp(3)+delta_across*w
        step=dmin+delta_across    ! update step
        call ulfind (x,y,z,ivnew,mat)
        !mat=ivolmat(ivnew)  ! set new material

      endif

       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Find particle in a given volume rotines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c*****************************************************************
c       Ana Farinha                                              
c       Luis Peralta                                             
c       Oct. 2006                                                
c*****************************************************************
c      Check if point(x,y,z) is inside box
c      itype=100
c      par(1)=Lx
c      par(2)=Ly
c      par(3)=Lz
      subroutine ulbox(xpar,par,i_flag)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xpar
      real*8, dimension(*) :: par
     
       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       xLx=par(1)/2.d0
       yLy=par(2)/2.d0
       zLz=par(3)/2.d0

       if ( (x .ge. -xLx) .and. (x .le. xLx) )then
          if ( (y .ge. -yLy) .and. (y .le. yLy) )then
             if ( (z .ge. -zLz) .and. (z .le. zLz))then

               i_flag=1
               
             endif
          endif
       endif

       end

c*****************************************************************
c       Ana Farinha                                              
c       Luis Peralta                                             
c       Oct. 2006 , Aug. 2010                                    
c      Check if point(x,y,z) is inside cylinder in zz
c      itype=200
c      par(1)=radius
c      par(2)=height
c*****************************************************************
      subroutine ulcyl(xpar,par,i_flag)
      use ullib_mod 
      implicit double precision(a-h,o-z), integer*4(i-n)  
      save
      real*8, dimension(3) :: xpar
      real*8, dimension(*) :: par
       
       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       r2=par(1)*par(1)
       h=par(2)/2.d0

       P=x*x + y*y

       if (P .le. r2)then
          if((z .ge. -h).and.(z .le. h))then

             i_flag=1 
             
          endif
       endif
       end
    
c***************************************************
c       Ana Farinha                                
c       Luis Peralta                               
c       Oct. 2006, Aug. 2010
c***************************************************
c      Check if point(x,y,z) is inside tube along zz
c      itype=210
c      par(1)=rmin
c      par(2)=rmax
c      par(3)=height
       subroutine ultub(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par

       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       rmin2=par(1)*par(1)
       rmax2=par(2)*par(2)
       h=par(3)/2.d0

       P=x*x+y*y

       if ((P .ge. rmin2).and. (P .le. rmax2))then
          if((z .ge. -h).and.(z .le. h))then
             
             i_flag=1 
             
          endif
       endif
       end

c*****************************************************************
c       Luis Peralta                                             
c                                                                
c       Mar. 2012                                                
c*****************************************************************
c      Check if point(x,y,z) is inside an elliptic cylinder in zz
c      itype=220
c      par(1)=a
c      par(2)=b
c      par(3)=h
       subroutine ulecyl(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par
       
       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       a=par(1)
       b=par(2)
       h=par(3)

c      check if the point is inside the volume
       P=(x*x)/(a*a) + (y*y)/(b*b) -1.d0

       if (P .le. 0)then
          if((z .ge. -h/2.).and.(z .le. +h/2.))then

             i_flag=1 
             
          endif
       endif
       end
       

c***************************************************
c       Ana Farinha                                 
c       Luis Peralta                                
c       Tiago Ribeiro                               
c       Jun. 2010, Aug. 2010                        
c***************************************************
c      Check if point(x,y,z) is inside a sphere 
c      itype=300
c      par(1)=radius
       subroutine ulsph(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par

       i_flag=0
       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       r=par(1)

       P=x*x+y*y+z*z

        if (P .le. r*r) then
          i_flag=1
        endif
       
       end

c***************************************************
c       Luis Peralta                                
c       Jul. 2011                           
c***************************************************
c      Check if point(x,y,z) is inside a cut sphere 
c      itype=310
c      par(1)=radius
c      par(2)=hmin    (zz axis)
c      par(3)=hmax    (zz axis)
       subroutine ulcsph(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par

       i_flag=0
       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       r=par(1)
       hmin=par(2)
       hmax=par(3)

       P=x*x+y*y+z*z

        if (P .le. r*r) then
         if(z .le. hmax) then 
          if(z .ge. hmin) then 
          i_flag=1
          endif
         endif
        endif
       
       end
c************************************************************
c     Luis  Peralta                                     
c     Jan. 2011                                            
c************************************************************
c     Check if point (x,y,z) is inside a cone or
c     inside a truncated cone with longitudinal axis in z   
c     Vertex is at (0,0,0) 
c
c       itype=400
c       par(1)=max radius               = Rmax
c       par(2)=height of truncated cone = h2 
c       par(3)=total height             = h
c
c      The reference point is at the vertex
       subroutine ulcone (xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par

       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)

       Rmax=par(1)
       h2=  par(2)
       h=   par(3)

       h1=h-h2
       tg=Rmax/h                 ! tan(theta)

       P2=x*x+y*y

       if((z .ge. h1).and.(z .le. h)) then
          if (P2 .le. tg*tg*z*z) then
            
             i_flag=1
             
          endif
       endif
       end     
c*****************************************************************
c       Luis Peralta   
c       Dec. 2011 
c*****************************************************************
c       Check if point(x,y,z) is inside pyramid
c
c       itype=500,501,502,503
c       par(1)=Lx
c       par(2)=Ly
c       par(3)=hmin
c       par(4)=hmax
c
c      The reference point is the base center 
       subroutine ulpyramid(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par
       real*8               :: Lx,Ly

       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)

       Lx=par(1)
       Ly=par(2)
       hmin=par(3)
       hmax=par(4)

       if ((z .ge. 0.).and.(z .le. hmin))then

         xmax= (hmax-z)* Lx /(2.*hmax)
         ymax= (hmax-z)* Ly /(2.*hmax)

          if((x .ge. -xmax).and.(x .le. xmax))then
             if((y .ge. -ymax).and.(y .le. ymax))then
               
               i_flag=1

             endif
          endif
       endif
       end

c*****************************************************************
c       Luis Peralta        
c       Dec. 2011      
c*****************************************************************
c       Check if point(x,y,z) is inside a wedge             
c                                              
c       itype=600
c       par(1)=Lx
c       par(2)=Ly
c       par(3)=hmin
c       par(4)=hmax
c
c      The reference point is the base left side 
       subroutine ulwedge(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par
       real*8               :: Lx,Ly


       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)

       Lx=par(1)
       Ly=par(2)
       hmin=par(3)
       hmax=par(4)

         Cx= Lx*hmax/(hmax-hmin)
         zmax=(Cx-x)*hmax/Cx

       if( x .ge. 0. .and. x .le. Lx)then
         if(y .ge. -Ly/2.d0 .and. y .le. +Ly/2.d0)then
           if(z .ge. 0. .and. z .le. zmax)then

            i_flag=1

           endif
         endif
       endif

       end
        
c************************************************************
c       Luis Peralta                                        
c       Dec. 2011     
c**********************************************************
c       Check if point(x,y,z) is inside paraboloid          
c
c       itype=700
c       par(1)=a
c       par(2)=b 
c       par(3)=zmax
c
c      The reference point is the vertex
       subroutine ulparaboloid(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       a=par(1)
       b=par(2)
       h=par(3)

       i_flag=0               
 
       if(z .ge.0 .and. z.le. h)then
        if ( (x/a)*(x/a) + (y/b)*(y/b) -z  .le. 0.d0) then
             
             i_flag=1

        endif      
       endif
       end

c**********************************************************
c       Ana Farinha
c       Luis Peralta  
c       Oct. 2006, Dec 2011   
c**********************************************************
c                                                          
c       Check if point(x,y,z) is inside ellipsoid          
c          (x/a)**2 + (y/b)**2 + (z/c)**2 =1               
c                                                          
c       itype=800
c       par(1)=a
c       par(2)=b 
c       par(3)=c
c       
c      The reference point is the figure center 
       subroutine ulellipsoid(xpar,par,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar
       real*8, dimension(*) :: par

       i_flag=0

       x=xpar(1)
       y=xpar(2)
       z=xpar(3)
       a=par(1)
       b=par(2)
       c=par(3)

       if ( (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c) -1.0 .le. 0.)then

              i_flag=1

       endif
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c    Computation of impact points and distances
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c**********************************************************
c       Luis Peralta  
c       Mar. 2012    
c**********************************************************
c       Determines the distance^2 from particle position xpar 
c       to impact point ximp
c       If impact point not along particle direction distance is 
c       set to infinity
c
        subroutine uldist2(xpar,upar,ximp,dist2)
        implicit double precision(a-h,o-z), integer*4(i-n)
        real*8, dimension(3) :: xpar,upar,ximp,xaux
        
          !xaux(1)=ximp(1)-xpar(1) 
          !xaux(2)=ximp(2)-xpar(2)
          !xaux(3)=ximp(3)-xpar(3)

          xaux=ximp-xpar

          if( ulvdot(upar,xaux) .ge. 0.d0 ) then ! particle is going in that direction
           dist2=ulvdot(xaux,xaux) ! squared distance to impact point

          else   ! particle is going the other direction
           dist2=1.0d+36  ! infinity
          endif
        
        end

c**********************************************************
c       Luis Peralta  
c       Mar. 2012
c       July 2021                                        
c**********************************************************
c       determines the minimum distance in array d2
c
c       nhit:  nb of impact points
c       d2  : array of distances^2
        subroutine ulfindmin(nhit,d2,imin)
        implicit double precision(a-h,o-z), integer*4(i-n)
        real*8, dimension(*) :: d2

        dmin2=d2(1)  ! start with first point
        imin=1

        do i=2,nhit    ! loop on impact points
        
          if(d2(i) .lt. dmin2) then   ! find minimum distance
           dmin2=d2(i)
           imin=i
          endif

        enddo
        end

c************************************************************
c     Luis Peralta                                          *
c     Mar 2012                                              *
c************************************************************
c  check if particle current path intersects a volume
c  iflag=0 no hit;  iflag=1 hit
      
       subroutine ulhit(x,y,z,u,v,w,ivol,iflag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp

        xpar(1)=x
        xpar(2)=y
        xpar(3)=z

        upar(1)=u
        upar(2)=v
        upar(3)=w

        iflag=0
        call uldvol(xpar,upar,ivol,ximp,dmin)
        if(dmin .lt. 1.0d+5) iflag=1         ! if distance lt 1 km 

        end

c************************************************************
c     Ana Farinha                                           *
c     Luis Peralta                                          *
c     2006, 2012, 2021                                      *
c************************************************************
c  Compute the hit point and minimum distance to boundary
c  following the particle current direction
      
      subroutine uldhit(xpar,upar,id,ximp,distmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp, ximp1
      
c     impact point in volume along line of flight
      involume=.true.     ! set flag for particle inside volume
      ivolcheck=id     ! volume under hit check. Debugging variable
      call uldvol(xpar,upar,id,ximp,distmin) 

      vnormal0=vnormaltmp ! normal vector

      nchild=mychild(id,0)

      do i=1,nchild
         ivchild=mychild(id,i)

c     impact point in child volume along line of flight
         involume=.false.     ! set flag for particle ouside volume
         ivolcheck=ivchild
         call uldvol(xpar,upar,ivchild,ximp1,distmin1) 

         sub=distmin-distmin1
         if(sub.gt.0.)then
           distmin=distmin1
            ximp=ximp1  ! hit point
            vnormal0=vnormaltmp ! normal vector

         endif
         
      enddo

      end

c************************************************************
c     Ana Farinha  2006                                     *
c     Luis Peralta  Dec. 2011                               *
c************************************************************
c      hit points on a box
       subroutine uldbox(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8, dimension(6) :: xlambda
       real*8, dimension(3) :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100) :: d2hit
       real*8, dimension(*) :: par
    
       nhit=0
       ihitemp=0

c      plane 1 and 2
       if(upar(1).ne.0)then
          xlambda(1)=(-par(1)/2.D0-xpar(1))/upar(1)
          xlambda(2)=(+par(1)/2.D0-xpar(1))/upar(1)

          nhit=nhit+1
          ihitemp(nhit)=1          
          hit(1)=-par(1)/2.D0
          hit(2)=xpar(2)+xlambda(1)*upar(2)
          hit(3)=xpar(3)+xlambda(1)*upar(3) 
    
            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

          nhit=nhit+1   
          ihitemp(nhit)=2       
          hit(1)=+par(1)/2.D0
          hit(2)=xpar(2)+xlambda(2)*upar(2)
          hit(3)=xpar(3)+xlambda(2)*upar(3)

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
         
       endif 
          
c      plane 3 and 4
       if(upar(2).ne.0)then
          xlambda(3)=(-par(2)/2.D0-xpar(2))/upar(2)
          xlambda(4)=(+par(2)/2.D0-xpar(2))/upar(2)

          nhit=nhit+1
          ihitemp(nhit)=3        
          hit(1)=xpar(1)+xlambda(3)*upar(1)
          hit(2)=-par(2)/2.D0
          hit(3)=xpar(3)+xlambda(3)*upar(3)

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

          nhit=nhit+1
          ihitemp(nhit)=4          
          hit(1)=xpar(1)+xlambda(4)*upar(1)
          hit(2)=+par(2)/2.D0
          hit(3)=xpar(3)+xlambda(4)*upar(3)

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif
          
c      plane 5 and 6
       if(upar(3).ne.0)then
          xlambda(5)=(-par(3)/2.D0-xpar(3))/upar(3)
          xlambda(6)=(+par(3)/2.D0-xpar(3))/upar(3)

          nhit=nhit+1
          ihitemp(nhit)=5          
          hit(1)=xpar(1)+xlambda(5)*upar(1)
          hit(2)=xpar(2)+xlambda(5)*upar(2)
          hit(3)=-par(3)/2.D0 

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
             
          nhit=nhit+1
          ihitemp(nhit)=6          
          hit(1)=xpar(1)+xlambda(6)*upar(1)
          hit(2)=xpar(2)+xlambda(6)*upar(2)
          hit(3)=+par(3)/2.D0

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c------------------------------------------------------------------
c       find closest hit point in particle direction  
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=ihitemp(imin)  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)



       call ulbox(hit,par,iflag)        ! check if hit point belongs to volume
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c----------------------------------------------------------------------------------

       end

c************************************************************
c     Ana Farinha   2006                                    *
c     Luis Peralta  Dec. 2011                               *
c************************************************************
c      hit points on a cylinder
       subroutine uldcyl(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par


       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       r=par(1)
       h=par(2)

       nhit=0
c      bottom side---------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        ihitemp(nhit)=1
        lambda=(-h/2.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = -h/2.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
       endif
       
c      top side ----------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        ihitemp(nhit)=2
        lambda=(+h/2.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = +h/2.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
       endif

c      round sides ---------------------------------------------
       a=ux*ux+uy*uy
       b=(2.*x0*ux) + (2.*y0*uy)
       c=x0*x0+y0*y0-r*r

       call ul2solver(a,b,c,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         ihitemp(nhit)=3
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

         nhit=nhit+1 
         ihitemp(nhit)=4
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c----------------------------------------------------------------
c       find closest hit point in particle direction 
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=ihitemp(imin)  ! copy hit face to common module
c -----------------------------------------------------

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)
       call ulcyl(hit,par,iflag)        ! check if hit point belongs to volume
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c----------------------------------------------------------------------

       end
       
c************************************************************
c     Ana Farinha 2006                                      *
c     Luis Peralta  Dec. 2011                               *
c************************************************************
c      hit points on a tube
       subroutine uldtub(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par


       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       rmin=par(1)
       rmax=par(2)
       h   =par(3)

       nhit=0
c      bottom side-------------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(-h/2.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = -h/2.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif
       
c      top side ----------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(+h/2.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = +h/2.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      round side inside -------------------------------------
       a=ux*ux+uy*uy
       b=(2.*x0*ux) + (2.*y0*uy)
       c=x0*x0+y0*y0-rmin*rmin

       call ul2solver(a,b,c,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)


         nhit=nhit+1 
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c      round side outside ------------------------------------
       a=ux*ux+uy*uy
       b=(2.*x0*ux) + (2.*y0*uy)
       c=x0*x0+y0*y0-rmax*rmax

       call ul2solver(a,b,c,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

         nhit=nhit+1 
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c-------------------------------------------------------------
c       find closest hit point in particle direction 
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ultub(hit,par,iflag)        ! check if hit point belongs to volume 
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c-------------------------------------------------------------------

       end


c************************************************************
c     Ana Farinha   2006                                    *
c     Luis Peralta  Dec. 2011                               *
c     2006                                                  *
c************************************************************
c      hit points on a elliptic cylinder
       subroutine uldecyl(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par

       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       aa=par(1)
       bb=par(2)
       h=par(3)

       nhit=0
c      bottom side-------------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(-h/2.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = -h/2.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
       endif
       
c      top side ----------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(+h/2.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = +h/2.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
       endif

c      round side --------------------------------------------
       a=(ux*ux)/(aa*aa) + (uy*uy)/(bb*bb)
       b=(2.*x0*ux)/(aa*aa) + (2.*y0*uy)/(bb*bb)
       c=(x0*x0)/(aa*aa)+(y0*y0)/(bb*bb)-1.d0

       call ul2solver(a,b,c,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

         nhit=nhit+1 
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c-----------------------------------------------------------------
c       find closest hit point in particle direction 
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)
       call ulecyl(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c--------------------------------------------------------------------

       end
       
c************************************************************
c     Ana Farinha 2006                                      *
c     Luis Peralta Dec. 2011                                *
c************************************************************
c      hit points on a sphere
       subroutine uldsph(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par


       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       r=par(1)

       nhit=0

       AA=1.d0  ! ux*ux+uy*uy+uz*uz
       BB=(2.*x0*ux) + (2.*y0*uy) + (2.*z0*uz)
       CC=x0*x0+y0*y0+z0*z0-r*r

       call ul2solver(AA,BB,CC,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

         nhit=nhit+1 
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c--------------------------------------------------------------
c       find closest hit point in particle direction  
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ulsph(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c------------------------------------------------------------------

       end

c*************************************************************
c       Luis Peralta                                          
c       Jul. 2011                                       
c*************************************************************
c      hit points on a truncated sphere  on both zz sides
       subroutine uldcsph(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par

       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

        R=par(1)
        hmin=par(2)
        hmax=par(3)

        nhit=0

c      intersection with sphere --------------------------------------

        a=1.d0 !ux*ux+uy*uy+uz*uz
        b=2.d0*(xpar(1)*upar(1)+xpar(2)*upar(2)+xpar(3)*upar(3))
        c=xpar(1)*xpar(1)+xpar(2)*xpar(2)+xpar(3)*xpar(3)-R*R

        call ul2solver(a,b,c,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

         nhit=nhit+1 
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c      intersection with planes ---------------------------------

        if(upar(3).ne.0.)then

          nhit=nhit+1
          xlambda=(hmax-xpar(3) ) / upar(3)
          hit(1)=xpar(1)+xlambda*upar(1)
          hit(2)=xpar(2)+xlambda*upar(2)
          hit(3)=xpar(3)+xlambda*upar(3)

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

          nhit=nhit+1
          xlambda=(hmin-xpar(3) ) / upar(3)
          hit(1)=xpar(1)+xlambda*upar(1)
          hit(2)=xpar(2)+xlambda*upar(2)
          hit(3)=xpar(3)+xlambda*upar(3)

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c-----------------------------------------------------------------
c       find closest hit point in particle direction  
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ulcsph(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c-----------------------------------------------------------------
       
       end

c************************************************************
c     Luis  Peralta                                     
c     Jan. 2011                                            
c************************************************************
c     Compute distance to and impact point for a truncated cone
c     along z. Vertex is at (0,0,0) 
c
c       itype=400
c       par(1)=max radius = Rmax
c       par(2)=h2 = height of truncated cone
c       par(3)=h = total height
c
c      the reference point is the vertex
c      hit points on a truncated cone
       subroutine uldcone(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par

       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       Rmax=par(1)
       h2=  par(2)
       h=   par(3)

       h1=h-h2
       tg=Rmax/h   ! tan(theta)
       hz=abs(z0)  ! height at point

       nhit=0   ! nb of impacts in the walls

c      cone smaller base  z=h1 -----------------------------
       if(uz .ne. 0.) then      ! no hit point if uz=0
        nhit=nhit+1
        alambda=(h1-z0)/uz
        hit(1)=x0+alambda*ux
        hit(2)=y0+alambda*uy
        hit(3)=h1

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      cone larger base  z=h=h1+h2 --------------------------
       if(uz .ne. 0.) then      ! no hit point if uz=0
        nhit=nhit+1
        alambda=(h-z0)/uz
        hit(1)=x0+alambda*ux
        hit(2)=y0+alambda*uy
        hit(3)=h

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      cone side --------------------------------------------
       a=ux*ux + uy*uy - uz*uz*tg*tg
       b=2.d0*(ux*x0 + uy*y0 - uz*z0*tg*tg)
       c=x0*x0 + y0*y0 - z0*z0*tg*tg

       delta=b*b-4.*a*c

       if( delta .gt. 0.) then

        nhit=nhit+1
        alambda= (-b -sqrt(delta) ) / (2.*a) ! -solution
        hit(1)=x0+alambda*ux
        hit(2)=y0+alambda*uy
        hit(3)=z0+alambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        nhit=nhit+1
        alambda= (-b +sqrt(delta) ) / (2.*a) ! +solution
        hit(1)=x0+alambda*ux
        hit(2)=y0+alambda*uy
        hit(3)=z0+alambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       elseif (delta .eq. 0.) then

        nhit=nhit+1
        alambda= (-b ) / (2.*a) ! only one solution
        hit(1)=x0+alambda*ux
        hit(2)=y0+alambda*uy
        hit(3)=z0+alambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c       find closest hit point in particle direction  ------------------------------
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ulcone(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c-------------------------------------------------------------------

       end
       
c************************************************************
c     Luis Peralta                                          *
c     Dec. 2011                                             *
c************************************************************
c      hit points on a pyramid
       subroutine uldpyramid(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par

       real*8 :: nx,ny,nz ! normal vector to side
       real*8 :: Lx,Ly,hmin,hmax
       real*8 :: ndotu

       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       Lx=par(1)
       Ly=par(2)
       hmin=par(3)
       hmax=par(4)

       costx= hmax/sqrt(hmax*hmax + Lx*Lx/4.)
       sintx= Lx/(2.*sqrt(hmax*hmax + Lx*Lx/4.))

       costy= hmax/sqrt(hmax*hmax + Ly*Ly/4.)
       sinty= Ly/(2.*sqrt(hmax*hmax + Ly*Ly/4.))
      
       nhit=0
c      top side-------------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(hmin-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = hmin

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif
       
c      bottom side ----------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(0.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = 0.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      left side  on xx -----------------------------------------
          nx=-costx
          ny=0.
          nz=sintx
          ndotu=nx*ux+ny*uy+nz*uz
        if(ndotu .ne. 0.) then
         nhit=nhit+1
         lambda= -(nx*x0+ny*y0+nz*(z0-hmax))/ndotu
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz  

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
        
        endif

c      right side  on xx --------------------------------------------
          nx=costx
          ny=0.
          nz=sintx
          ndotu=nx*ux+ny*uy+nz*uz
        if(ndotu .ne. 0.) then
         nhit=nhit+1
         lambda= -(nx*x0+ny*y0+nz*(z0-hmax))/ndotu
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz   

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
       
        endif

c      left side  on yy --------------------------------------------
          nx=0.
          ny=-costy
          nz=sinty
          ndotu=nx*ux+ny*uy+nz*uz
        if(ndotu .ne. 0.) then
         nhit=nhit+1
         lambda= -(nx*x0+ny*y0+nz*(z0-hmax))/ndotu
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz 

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
         
        endif

c      right side  on yy ------------------------------------------
          nx=0.
          ny=costy
          nz=sinty
          ndotu=nx*ux+ny*uy+nz*uz
        if(ndotu .ne. 0.) then
         nhit=nhit+1
         lambda= -(nx*x0+ny*y0+nz*(z0-hmax))/ndotu
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
      
        endif

c-------------------------------------------------------------------------
c       find closest hit point in particle direction  
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ulpyramid(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c--------------------------------------------------------------
     
       end

c************************************************************
c     Luis Peralta                                          *
c     Dec. 2011                                             *
c************************************************************
c      hit points on a wedge
       subroutine uldwedge(par,xpar,upar,ximp,dmin)  
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par

       real*8 :: nx,ny,nz ! normal vector to side
       real*8 :: Lx,Ly,hmin,hmax
       real*8 :: ndotu

       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       Lx=par(1)
       Ly=par(2)
       hmin=par(3)
       hmax=par(4)

       Cx= Lx*hmax/(hmax-hmin)
       cost= hmax/sqrt(hmax*hmax + Cx*Cx)
       sint=   Cx/sqrt(hmax*hmax + Cx*Cx)
      
       nhit=0
       
c      bottom side ----------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(0.-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = 0.

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      top slanted side -----------------------------------------
          nx=cost
          ny=0.
          nz=sint
          ndotu=nx*ux+ny*uy+nz*uz
        if(ndotu .ne. 0.) then
         nhit=nhit+1
         lambda= -(nx*x0+ny*y0+nz*(z0-hmax))/ndotu
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)
     
        endif

c      left side  on xx -----------------------------------
       if(ux.ne.0.)then 
        nhit=nhit+1 
        lambda=(0.-x0)/ux
        hit(1) = 0.
        hit(2) = y0 + lambda*uy
        hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      right side  on xx --------------------------------
       if(ux.ne.0.)then 
        nhit=nhit+1 
        lambda=(Lx-x0)/ux
        hit(1) = Lx
        hit(2) = y0 + lambda*uy
        hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      left side  on yy -------------------------------
       if(uy.ne.0.)then 
        nhit=nhit+1 
        lambda=(-Ly/2.-y0)/uy
        hit(1) = x0 + lambda*ux
        hit(2) = -Ly/2.
        hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      right side  on yy ----------------------------
       if(uy.ne.0.)then 
        nhit=nhit+1 
        lambda=(+Ly/2.-y0)/uy
        hit(1) = x0 + lambda*ux
        hit(2) = +Ly/2.
        hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c       --------------------------------------------------------------------
c       find closest hit point in particle direction  
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ulwedge(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c--------------------------------------------------------------
       
       end
       
c************************************************************
c     Luis Peralta                                          *
c     Dec 2011                                              *
c************************************************************
c      hit points on a paraboloid
       subroutine uldparaboloid(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par

       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       a=par(1)
       b=par(2)
       zmax=par(3)

       nhit=0
       
c      top side ----------------------------------------------
       if(uz.ne.0.)then 
        nhit=nhit+1 
        lambda=(zmax-z0)/uz
        hit(1) = x0 + lambda*ux
        hit(2) = y0 + lambda*uy
        hit(3) = zmax

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

       endif

c      paraboloid sides ----------------------------------------

       AA=(ux/a)*(ux/a)+(uy/b)*(uy/b)
       BB=(2.*x0*ux)/(a*a) + (2.*y0*uy)/(b*b)-uz
       CC=(x0/a)*(x0/a)+(y0/b)*(y0/b)-z0

       call ul2solver(AA,BB,CC,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

         nhit=nhit+1 
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c       -------------------------------------------------------------
c       find closest hit point in particle direction  
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ulparaboloid(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c---------------------------------------------------------------

       end

c************************************************************
c     Luis Peralta                                          *
c     Dec 2011                                              *
c************************************************************
c      hit points on a ellipsoid
       subroutine uldellipsoid(par,xpar,upar,ximp,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar, upar, ximp
       real*8               :: lambda
       real*8, dimension(3)     :: hit
       real*8, dimension(100,3) :: xhit
       real*8, dimension(100)   :: d2hit
       real*8, dimension(*)     :: par

       x0=xpar(1)
       y0=xpar(2)
       z0=xpar(3)

       ux=upar(1)
       uy=upar(2)
       uz=upar(3)

       a=par(1)
       b=par(2)
       c=par(3)

       nhit=0

       AA=(ux/a)*(ux/a)+(uy/b)*(uy/b)+(uz/c)*(uz/c)
       BB=(2.*x0*ux)/(a*a) + (2.*y0*uy)/(b*b) + (2.*z0*uz)/(c*c)
       CC=(x0/a)*(x0/a)+(y0/b)*(y0/b)+(z0/c)*(z0/c)-1.

       call ul2solver(AA,BB,CC,xlamb1,xlamb2,ierr)

        if(ierr.eq.0)then
         nhit=nhit+1 
         lambda=xlamb1
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

         nhit=nhit+1 
         lambda=xlamb2
         hit(1) = x0 + lambda*ux
         hit(2) = y0 + lambda*uy
         hit(3) = z0 + lambda*uz

            call uldist2(xpar,upar,hit,dist2)
            d2hit(nhit)=dist2
            xhit(nhit,1)=hit(1)
            xhit(nhit,2)=hit(2)
            xhit(nhit,3)=hit(3)

        endif

c----------------------------------------------------------------------
c       find closest hit point in particle direction 
c       belonging to a volume
        if(nhit .eq. 0) then    
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
        endif

       itry=0
10     continue
       call ulfindmin(nhit,d2hit,imin)  ! find closest hit point

       ihitface=imin  ! copy hit face to common module

       if(involume)then
        delta=-delta_across
       else
        delta=delta_across
       endif
       hit(1)=xhit(imin,1)+delta*upar(1) ! put point on the right boarder side
       hit(2)=xhit(imin,2)+delta*upar(2)
       hit(3)=xhit(imin,3)+delta*upar(3)

       call ulellipsoid(hit,par,iflag)        ! check if hit point belongs to volume    
       if(iflag.ne.0) then
        dmin=sqrt(d2hit(imin))
        ximp(1)=xhit(imin,1)
        ximp(2)=xhit(imin,2)
        ximp(3)=xhit(imin,3)
        return
       else     ! hit point is outside volume
        d2hit(imin)=1.0d+38
        itry=itry+1
       endif

       if(itry.le.nhit) then
        goto 10
       else
         dmin=1.D+38
         ximp(1)=1.d+36
         ximp(2)=1.d+36
         ximp(3)=1.d+36
         return
       endif
c-----------------------------------------------------------------------

       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c    Peak-a-volume routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c************************************************************
c     Ana Farinha & Luis Peralta    2006                    *
c     Luis Peralta revised Dec. 2011                        *
c                          May 2020                         *
c************************************************************
c     choose routine to compute hit point along particle direction 
c     and minimum distance      
      subroutine uldvol(xin,uin,id,ximpout,dmin)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xin,uin,ximp,ximpout,xpar,upar,xpar0
       real*8, dimension(3) :: xrot,urot
       real*8, dimension(nparadim) :: par
  
      itype=mytype(id)
      irot=myrot(id)
      angle=volangle(id)

      !xpar0(1)=xin(1)     ! keep input variables unchanged
      !xpar0(2)=xin(2)
      !xpar0(3)=xin(3)
      !upar(1)=uin(1)
      !upar(2)=uin(2)
      !upar(3)=uin(3)

      xpar0 = xin
      upar = uin

      call ulgetpar(id,par)
      call ultrans(xpar0,xpar,id)

      if(irot.ne.0)then
       call ulrot(id,xpar,xrot)
       !xpar(1)=xrot(1)
       !xpar(2)=xrot(2)
       !xpar(3)=xrot(3)

       xpar = xrot

       call ulrot(id,upar,urot)
       !upar(1)=urot(1)
       !upar(2)=urot(2)
       !upar(3)=urot(3)
       upar=urot

      endif
   
      select case(itype)

        case(100) ! box
         call uldbox(par,xpar,upar,ximp,dmin)


        case(200,203) ! cylinder
         call uldcyl(par,xpar,upar,ximp,dmin)

        case(201)
         call ulrotxz(xpar)
         call ulrotxz(upar)
         call uldcyl(par,xpar,upar,ximp,dmin)
         call ulrotxz_1(ximp)

        case(202)
         call ulrotyz(xpar)
         call ulrotyz(upar)
         call uldcyl(par,xpar,upar,ximp,dmin)
         call ulrotyz_1(ximp)

        case(210,213)! tube
         call uldtub(par,xpar,upar,ximp,dmin)

        case(211)
         call ulrotxz(xpar)
         call ulrotxz(upar)
         call uldtub(par,xpar,upar,ximp,dmin)
         call ulrotxz_1(ximp)

        case(212)
         call ulrotyz(xpar)
         call ulrotyz(upar)
         call uldtub(par,xpar,upar,ximp,dmin)
         call ulrotyz_1(ximp)

        case(220,223) ! elliptic cylinder
         call uldecyl(par,xpar,upar,ximp,dmin)
         
        case(221)
         call ulrotxz(xpar)
         call ulrotxz(upar)
         call uldecyl(par,xpar,upar,ximp,dmin)
         call ulrotxz_1(ximp)

        case(222)
         call ulrotyz(xpar)
         call ulrotyz(upar)
         call uldecyl(par,xpar,upar,ximp,dmin)
         call ulrotyz_1(ximp)

        case(300) ! sphere 
         call uldsph(par,xpar,upar,ximp,dmin)

        case(310,313)! cut sphere
         call uldcsph(par,xpar,upar,ximp,dmin)
        case(-310,-313)
          call ulflipz(xpar)
          call ulflipz(upar)
          call uldcsph(par,xpar,upar,ximp,dmin)
          call ulflipz(ximp)

        case(311)
         call ulrotxz(xpar)
         call ulrotxz(upar)
         call uldcsph(par,xpar,upar,ximp,dmin)
         call ulrotxz_1(ximp)
        case(-311)
          call ulflipx(xpar)
          call ulflipx(upar)      
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldcsph(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
          call ulflipx(ximp)

        case(312)
         call ulrotyz(xpar)
         call ulrotyz(upar)
         call uldcsph(par,xpar,upar,ximp,dmin)
         call ulrotyz_1(ximp)
        case(-312)
          call ulflipy(xpar)
          call ulflipy(upar) 
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldcsph(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
          call ulflipy(ximp)

        case(400,403)! cone
          call uldcone(par,xpar,upar,ximp,dmin)
        case(-400,-403)
          call ulflipz(xpar)
          call ulflipz(upar)
          call uldcone(par,xpar,upar,ximp,dmin)
          call ulflipz(ximp)

        case(401)
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldcone(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
        case(-401)    
          call ulflipx(xpar)
          call ulflipx(upar)      
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldcone(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
          call ulflipx(ximp)

        case(402)
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldcone(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
        case(-402)
          call ulflipy(xpar)
          call ulflipy(upar) 
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldcone(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
          call ulflipy(ximp)

        case(500,503)! pyramid
         call uldpyramid(par,xpar,upar,ximp,dmin)
        case(-500,-503)
          call ulflipz(xpar)
          call ulflipz(upar)
          call uldpyramid(par,xpar,upar,ximp,dmin)
          call ulflipz(ximp)


        case(501)
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldpyramid(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
        case(-501)
          call ulflipx(xpar)
          call ulflipx(upar)      
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldpyramid(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
          call ulflipx(ximp)

        case(502)
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldpyramid(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
        case(-502)
          call ulflipy(xpar)
          call ulflipy(upar) 
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldpyramid(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
          call ulflipy(ximp)

        case(600,603)! wedge
         call uldwedge(par,xpar,upar,ximp,dmin)
        case(-600,-603)
          call ulflipz(xpar)
          call ulflipz(upar)
          call uldwedge(par,xpar,upar,ximp,dmin)
          call ulflipz(ximp)

        case(601)  
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldwedge(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
        case(-601)  
          call ulflipx(xpar)
          call ulflipx(upar)      
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldwedge(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
          call ulflipx(ximp)

        case(602)       
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldwedge(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
        case(-602)
          call ulflipy(xpar)
          call ulflipy(upar) 
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldwedge(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
          call ulflipy(ximp)

        case(700,703)! paraboloid
         call uldparaboloid(par,xpar,upar,ximp,dmin)
        case(-700,-703)
          call ulflipz(xpar)
          call ulflipz(upar)
         call uldparaboloid(par,xpar,upar,ximp,dmin)
          call ulflipz(ximp)

        case(701)   
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldparaboloid(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
        case(-701) 
          call ulflipx(xpar)
          call ulflipx(upar)      
          call ulrotxz(xpar)
          call ulrotxz(upar)
          call uldparaboloid(par,xpar,upar,ximp,dmin)
          call ulrotxz_1(ximp)
          call ulflipx(ximp)

        case(702)
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldparaboloid(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
        case(-702)
          call ulflipy(xpar)
          call ulflipy(upar) 
          call ulrotyz(xpar)
          call ulrotyz(upar)
          call uldparaboloid(par,xpar,upar,ximp,dmin)
          call ulrotyz_1(ximp)
          call ulflipy(ximp)

        case(800) ! ellipsoid 
         call uldellipsoid(par,xpar,upar,ximp,dmin)

      case default
        write(*,*)'ULDVOL: volume unknown itype=',itype
        STOP  
      end select

      if(irot.ne.0)then
       call ulrot_1(id,ximp,xrot)
       !ximp(1)=xrot(1)
       !ximp(2)=xrot(2)
       !ximp(3)=xrot(3)

       ximp=xrot

      endif

      call ultrans_1(ximp,ximpout,id)

      ivolhitface(id)=ihitface   ! keep hit face for volume id

      end

c***************************************************
c       Ana Farinha                                *
c       Luis Peralta                               *
c       Oct. 2006  , Mar 2012                      *
c       May 2020                                   *
c***************************************************
c     Look for particle in volume id   
       subroutine ulvolsearch(xin,id,i_flag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xin,xpar,xpar0,xrot
       real*8, dimension(nparadim) :: par


      itype=mytype(id)
      irot=myrot(id)
      angle=volangle(id)

      !xpar0(1)=xin(1)   ! keep input variable unchanged
      !xpar0(2)=xin(2)
      !xpar0(3)=xin(3)

      xpar0 = xin

      call ulgetpar(id,par)

      call ultrans(xpar0,xpar,id)

      if(irot.ne.0)then
       call ulrot(id,xpar,xrot)
       !xpar(1)=xrot(1)
       !xpar(2)=xrot(2)
       !xpar(3)=xrot(3)

       xpar=xrot
      endif

      select case(itype)
        case(100) ! box
         call ulbox(xpar,par,i_flag)

        case(200,203) ! cylinder
         call ulcyl(xpar,par,i_flag)

        case(201)       
         call ulrotxz(xpar)
         call ulcyl(xpar,par,i_flag)

        case(202)
         call ulrotyz(xpar)
         call ulcyl(xpar,par,i_flag)

        case(210,213)! tube
         call ultub(xpar,par,i_flag)

        case(211)
         call ulrotxz(xpar)
         call ultub(xpar,par,i_flag)

        case(212)
         call ulrotyz(xpar)
         call ultub(xpar,par,i_flag)

        case(220,223) ! elliptic cylinder
         call ulecyl(xpar,par,i_flag)

        case(221)      
         call ulrotxz(xpar)
         call ulecyl(xpar,par,i_flag)

        case(222)
         call ulrotyz(xpar)
         call ulecyl(xpar,par,i_flag)

        case(300)! sphere 
         call ulsph(xpar,par,i_flag)

        case(310,313)! cut sphere
         call ulcsph(xpar,par,i_flag)
        case(-310,-313)
         call ulflipz(xpar)
         call ulcsph(xpar,par,i_flag)

        case(311)
         call ulrotxz(xpar)
         call ulcsph(xpar,par,i_flag)
        case(-311)
          call ulflipx(xpar)
          call ulrotxz(xpar)
          call ulcsph(xpar,par,i_flag)

        case(312)
         call ulrotyz(xpar)
         call ulcsph(xpar,par,i_flag)
        case(-312)
         call ulflipy(xpar)
         call ulrotyz(xpar)
         call ulcsph(xpar,par,i_flag)

        case(400,403)! cone
         call ulcone(xpar,par,i_flag)
        case(-400,-403)! cone
         call ulflipz(xpar)
         call ulcone(xpar,par,i_flag)

        case(401)
          call ulrotxz(xpar)
          call ulcone(xpar,par,i_flag)
        case(-401)
          call ulflipx(xpar)
          call ulrotxz(xpar)
          call ulcone(xpar,par,i_flag)

        case(402)
         call ulrotyz(xpar)
         call ulcone(xpar,par,i_flag)
        case(-402)
         call ulflipy(xpar)
         call ulrotyz(xpar)
         call ulcone(xpar,par,i_flag)

        case(500,503)! pyramid
         call ulpyramid(xpar,par,i_flag)
        case(-500,-503)
         call ulflipz(xpar)
         call ulpyramid(xpar,par,i_flag)

        case(501)
          call ulrotxz(xpar)
          call ulpyramid(xpar,par,i_flag)
        case(-501)
          call ulflipx(xpar)
          call ulrotxz(xpar)
          call ulpyramid(xpar,par,i_flag)

        case(502)
         call ulrotyz(xpar)
         call ulpyramid(xpar,par,i_flag)
        case(-502)
         call ulflipy(xpar)
         call ulrotyz(xpar)
         call ulpyramid(xpar,par,i_flag)

        case(600,603)! wedge
         call ulwedge(xpar,par,i_flag)
        case(-600,-603)
         call ulflipz(xpar)
         call ulwedge(xpar,par,i_flag)

        case(601)
          call ulrotxz(xpar)
          call ulwedge(xpar,par,i_flag)
        case(-601)
          call ulflipx(xpar)
          call ulrotxz(xpar)
          call ulwedge(xpar,par,i_flag)

        case(602)
         call ulrotyz(xpar)
         call ulwedge(xpar,par,i_flag)
        case(-602)
         call ulflipy(xpar)
         call ulrotyz(xpar)
         call ulwedge(xpar,par,i_flag)

        case(700,703)! paraboloid
         call ulparaboloid(xpar,par,i_flag)
        case(-700,-703)
         call ulflipz(xpar)
         call ulparaboloid(xpar,par,i_flag)

        case(701)
          call ulrotxz(xpar)
          call ulparaboloid(xpar,par,i_flag)
        case(-701)
          call ulflipx(xpar)
          call ulrotxz(xpar)
          call ulparaboloid(xpar,par,i_flag)

        case(702)
         call ulrotyz(xpar)
         call ulparaboloid(xpar,par,i_flag)
        case(-702)
         call ulflipy(xpar)
         call ulrotyz(xpar)
         call ulparaboloid(xpar,par,i_flag)

        case(800) ! ellipsoid 
         call ulellipsoid(xpar,par,i_flag)

        case default
         write(*,*)' ERROR: Unknown volume type -> STOP'
         stop
      end select

      end


c***************************************************************
c       Luis Peralta                                           *
c       Dec. 2013                                              *
c***************************************************************
c
c     Inflate a volume by volscale and look for a hit on 
c     new volume.
c     iflag=0 no hit
c     iflag=1 hit
c
      subroutine ultarget(id,volscale,x,y,z,u,v,w,iflag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xpar,upar,ximp
       real*8, dimension(nparadim) :: par0


c     scale volume id parameters by volscale factor
      do i=1, nparadim
         par0(i)=volpar(id,i)
         volpar(id,i)=par0(i)*volscale
      enddo  

         xpar(1)=x
         xpar(2)=y
         xpar(3)=z

         upar(1)=u
         upar(2)=v
         upar(3)=w

         involume=.false.                 ! check if point is inside volume id
         call ulvolsearch(xpar,id,iflag)
         if(iflag.ne.0) involume=.true.

       iflag=0
       call uldvol(xpar,upar,id,ximp,dmin)
      if(dmin .lt. 1.0d+5) iflag=1         ! if distance lt 1 km 

c     put back original values
      do i=1, nparadim
         volpar(id,i)=par0(i)
      enddo  
      
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c    Getting info
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This routines are meant to be used by the user
c
c***************************************************************
c       Luis Peralta                                           *
c       Oct. 2009                                              *
c***************************************************************
c      get the material type of a volume
       subroutine ulgetmat(id,mat)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save

      itype=mytype(id)
      if(itype.eq.0) then
       write(*,*)'ULGETMAT: Unknown volume id:',id
      endif

       mat = ivolmat(id)

      end

c***************************************************************
c       Luis Peralta                                           *
c       Oct. 2009                                              *
c***************************************************************
c      get the total number of volumes
       subroutine ulgetnvolu(ntot_vol)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save

       ntot_vol = nvolu

      end

c***************************************************************
c       Luis Peralta                                           *
c       2010                                                   *
c***************************************************************
c      get volume parameters par(i)
       subroutine ulgetpar(id,par)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(nparadim) :: par
      
      itype=mytype(id)
      if(itype.eq.0) then
       write(*,*)'ULGETPAR: Unknown volume id:',id
      endif

      do i=1, nparadim
         par(i)=volpar(id,i)
      enddo  
      
      end

c***************************************************************
c       Luis Peralta                                           *
c       2012                                                   *
c***************************************************************
c      change one volume parameter
      
       subroutine ulput1par(ivol,ipar,value)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save

      itype=mytype(ivol)
      if(itype.eq.0) then
       write(*,*)'ULPUT1PAR: Unknown volume id:',ivol
       stop
      endif
         if(ipar .lt.1 .or. ipar .gt. nparadim) 
     &   stop 'ULPUT1PAR: parameter unknown'
         volpar(ivol,ipar) = value
          
      end
c***************************************************************
c       Luis Peralta                                           *
c       Jan. 2014                                              *
c***************************************************************
c     change volume parameters
      
       subroutine ulputpar(ivol,parin,npar)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(*) :: parin

      itype=mytype(ivol)
      if(itype.eq.0) then
       write(*,*)'ULPUTPAR: Unknown volume id:',ivol
       stop
      endif
       do i=1,npar
         volpar(ivol,ipar) = parin(i)
       enddo   
      end

c************************************************************
c     Luis Peralta                                          *
c     Dec. 2013                                             *
c************************************************************
c     Get body center position coordinates

       subroutine ulgetcenter(id,xcenter)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xcenter

      itype=mytype(id)
      if(itype.eq.0) then
       write(*,*)'ULGETCENTER: Unknown volume id:',id
      endif

      xcenter(1)=volcenter(id,1)
      xcenter(2)=volcenter(id,2)
      xcenter(3)=volcenter(id,3)

      end

c***************************************************************
c       Luis Peralta                                           *
c       2010                                                   *
c***************************************************************
c
c      output geometry summary
       subroutine ulgsummary
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save

1     format(' Number of defined volumes:',I5)
      write(*,1)nvolu
      write(*,*)' '
      write(*,*)'-----------------------------------'
      write(*,*)'Volume parameters'
      write(*,*)'-----------------------------------'

      do i=1,nvolu
       ivol= idall(i)
       itype=mytype(ivol)
       call ulgetnpar(itype,npar)   
       write(*,*)' '
2      format(' Volume id:',I5,'    itype:',I5)
       write(*,2)ivol,itype
3             format(' Dimensions:      ',6(1x,g12.6))
              write(*,3)(volpar(ivol,j), j=1,npar)

4             format(' Center position: ',3(1x,g12.6))
              write(*,4)(volcenter(ivol,j), j=1,3) 

5      format(' Rot Angle:        ',1x,f8.2)
       write(*,5) volangle(ivol)

      enddo
      write(*,*)' '
      write(*,*)'-----------------------------------'
      write(*,*)'Children volumes'
      write(*,*)'-----------------------------------'
      write(*,*)' '

      do i=1,nvolu
       ivol= idall(i)
6        format(' Volume id:', I5)
         write(*,6)ivol
7        format(' Number of children:',I5)
         write(*,7)mychild(ivol,0)
         if((mychild(ivol,0)).ne.0)then 
8           format(' Child id: ',5I5)   
            write(*,8)(mychild(ivol,j), j=1,mychild(ivol,0)) 
         endif
         write(*,*)' '      
      enddo
      end

C*********************************************************
C        Luis Peralta
c        Dec. 2013
c        May 2020                                   
C*********************************************************
C       Computes a body volume
c       id           body id
c       v            body volume
  
       subroutine ulbodyvol(id,v)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(nparadim) :: par

      itype=mytype(id)
      if(itype.eq.0) then
       write(*,*)'ULBODYVOL: Unknown volume id:',id
      endif

        call ulgetpar(id,par)

        itype=mytype(id)        ! get volume type
        !pi=acos(-1.d0)


        select case(itype)

        case(100) ! box
c        par(1)=Lx, par(2)=Ly, par(3)=Lz         
         v= par(1)*par(2)*par(3)
c-----------------------------------------------------------------------------
        case(200,201,202,203) ! cylinder
c       par(1)=r, par(2)=h   
         v=pi * par(1)**2 * par(2)
c------------------------------------------------------------------------------
        case(210,211,212,213)! tube
c       par(1)=rmin,  par(2)=rmax, par(3)=height
        v1=pi * par(1)**2 * par(3)
        v2=pi * par(2)**2 * par(3)
        v=v2-v1
c------------------------------------------------------------------------------
        case(220,221,222,223) ! elliptic cylinder
c      par(1)=a, par(2)=b, par(3)=h
       v=pi*par(1)*par(2)*par(3)
c------------------------------------------------------------------------------
        case(300)! sphere
c        par(1)=r
         v=4.d0/3.d0*pi*par(1)**3
c--------------------------------------------------------------------------------
        case(310,311,312,313,-310,-311,-312,-313)! cut sphere
c       par(1)=r, par(2)=hmin, par(3)=hmax    (zz axis)

       r=par(1)
       hmin=par(2)
       hmax=par(3)

       h1=r-hmin
       a1=sqrt(r**2-hmin**2)
       v1=(pi*h1/6.d0)*(3.d0*a1**2+h1**2) ! dome volume = pi*h/6 * (3a**2 +h**2)
       h2=r-hmax
       a2=sqrt(r**2-hmax**2)
       v2=(pi*h2/6.d0)*(3.d0*a2**2+h2**2)
       v=v1-v2
c--------------------------------------------------------------------------------
        case(400,401,402,403,-400,-401,-402,-403)! cone
c       par(1)=Rmax, par(2)=height of truncated cone = h2 
c       par(3)=h (total)
       r=par(1)
       h2=par(2)
       h=par(3)
       v2=1.d0/3.d0 * pi*r**2 * h ! full cone
       h1=h-h2
       r1=r/h * h1
       v1=1.d0/3.d0 * pi*r1**2 * h1 ! cone upper part
       v=v2-v1
C----------------------------------------------------------------------------------
        case(500,501,502,503,-500,-501,-502,-503)! pyramid
c       par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
       xL=par(1)
       yL=par(2)
       hmin=par(3)
       hmax=par(4)

       v2=1.d0/3.d0 * hmax * xL * yL ! full pyramid
       h1=hmax-hmin
       xL1=xL/hmax * h1
       yL1=yL/hmax * h1
       v1=1.d0/3.d0 * h1 * xL1 * yL1 ! upper part
       v=v2-v1
c-----------------------------------------------------------------------------------
        case(600,601,602,603,-600,-601,-602,-603)! wedge
c       par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
       xL=par(1)
       yL=par(2)
       hmin=par(3)
       hmax=par(4)
       Cx=xL*hmax/(hmax-hmin)

       v2= Cx*hmax*(yL/3.+yL/6.) ! full wedge
       v1= (Cx-xL)*hmin*(yL/3.+yL/6.) ! small wedge
       v=v2-v1
c-----------------------------------------------------------------------------------
        case(700,701,702,703,-700,-701,-702,-703)! paraboloid
c       par(1)=a, par(2)=b, par(3)=zmax
       a=par(1)
       b=par(2)
       h=par(3)

       v=pi/2. * a*b*h
c-----------------------------------------------------------------------------------
        case(800)! ellipsoid
c       par(1)=a, par(2)=b, par(3)=c
       a=par(1)
       b=par(2)
       c=par(3)

       v=4./3. *pi * a*b*c
       
        case default
          write(*,*)'ULBODYVOL: Unknown volume type:',itype
          v=0.

        end select
        end

C*********************************************************
C        Luis Peralta
c        Dec. 2013
C*********************************************************
C       Checks if body id1 is inside body id2 using randomly selected points
c       id1,id2      body ids
c       iflag = 0    outside
c       iflag = 1    inside
c       iflag = 2    partially inside/outside
  
       subroutine ulinchk(id1,id2,iflag)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(3) :: xcenter, xin

      itype=mytype(id1)
      if(itype.eq.0) then
       write(*,*)'ULINCHK: Unknown volume id:',id1
      endif
      itype=mytype(id2)
      if(itype.eq.0) then
       write(*,*)'ULINCHK: Unknown volume id:',id2
      endif

        call ulbubbler(id1,rc)  ! get bubble radius

        ntry=1000000
        nfail=0

        call ulgetcenter(id1,xcenter)
        xc=xcenter(1)
        yc=xcenter(2)
        zc=xcenter(3)
        do i=1, ntry
          call ulrndvol(id1,xc,yc,zc,rc,x,y,z)
          xin(1)=x
          xin(2)=y
          xin(3)=z
          call ulvolsearch(xin,id2,iflag)
          if(iflag.eq.0) nfail=nfail+1
        enddo

        if(nfail.eq. ntry)then
          iflag=0    ! outside
        elseif(nfail.gt. 0)then
          iflag=2    ! partially inside/outside
        else
          iflag=1    ! inside
        endif

        end

C*********************************************************
C        Luis Peralta
c        Dec. 2013
c        May 2020                                   
C*********************************************************
c      get the best bubble radius to be used in ulrndvol
  
       subroutine ulbubbler(id,rc)
       use ullib_mod 
       implicit double precision(a-h,o-z), integer*4(i-n)  
       save
       real*8, dimension(nparadim) :: par

      itype=mytype(id)
      if(itype.eq.0) then
       write(*,*)'ULINCHK: Unknown volume id:',id1
      endif

        call ulgetpar(id,par)
        !pi=acos(-1.d0)

       select case(itype)

        case(100)  ! box   par(1)=Lx, par(2)=Ly, par(3)=Lz         
         rc= sqrt(par(1)**2 + par(2)**2 + par(3)**2)/2. ! diagonal= sqrt(Lx**2+Ly**2+Lz**2)
c-----------------------------------------------------------------------------
        case(200,201,202,203)  ! cylinder   par(1)=r, par(2)=h   
         r=par(1)
         h=par(2)
         rc=sqrt(h**2 + (2*r)**2)/2.
c------------------------------------------------------------------------------
        case(210,211,212,213)  ! tube   par(1)=rmin,  par(2)=rmax, par(3)=height
         r=par(2)
         h=par(3)
         rc=sqrt(h**2 + (2*r)**2)/2.
c------------------------------------------------------------------------------
        case(220,221,222,223)  ! elliptic cylinder   par(1)=a, par(2)=b, par(3)=h
         r=max(par(1),par(2))
         h=par(3)
         rc=sqrt(h**2 + (2*r)**2)/2.
c------------------------------------------------------------------------------
        case(300) ! sphere   par(1)=r
         rc=par(1)
c--------------------------------------------------------------------------------
        case(310,311,312,313,-310,-311,-312,-313)  ! cut sphere  par(1)=r, par(2)=hmin, par(3)=hmax (zz axis)
         rc=par(1)
c--------------------------------------------------------------------------------
        case(400,401,402,403,-400,-401,-402,-403)  ! cone par(1)=Rmax, par(2)= h2  par(3)=h (total)
         Rmax=par(1)
         h=par(3)
         rc=sqrt(Rmax**2+h**2)
C----------------------------------------------------------------------------------
        case(500,501,502,503,-500,-501,-502,-503) ! pyramid par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
         xL=par(1)
         yL=par(2)
         d=sqrt(xL**2+yL**2)/2.
         hmax=par(4)
         rc=max(d,hmax)
c-----------------------------------------------------------------------------------
        case(600,601,602,603,-600,-601,-602,-603) ! wedge par(1)=Lx, par(2)=Ly, par(3)=hmin, par(4)=hmax
         xL=par(1)
         yL=par(2)
         d1=sqrt(xL**2+(yL/2.)**2)
         d2=sqrt(yL**2+(xL/2.)**2)
         hmax=par(4)
         rc=max(d1,d2,hmax)
c-----------------------------------------------------------------------------------
        case(700,701,702,703,-700,-701,-702,-703)! paraboloid  par(1)=a, par(2)=b, par(3)=zmax
         a=par(1)
         b=par(2)
         h=par(3)

         d1=sqrt(h**2+h*a**2)
         d2=sqrt(h**2+h*b**2)
         rc=max(d1,d2)
c-----------------------------------------------------------------------------------
        case(800) ! ellipsoid  par(1)=a, par(2)=b, par(3)=c
         a=par(1)
         b=par(2)
         c=par(3)
         rc=max(a,b,c)

        case default
          rc=0.
          write(*,*)'ULBUBBLER: Unknown volume type:',id1

        end select

        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c    Functions and utility routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c***************************************************************
c       Luis Peralta                                           *
c       2012                                                   *
c***************************************************************
c      2nd degree equation solver
       subroutine ul2solver(a,b,c,x1,x2,ierr)
       implicit double precision(a-h,o-z), integer*4(i-n)

       ierr=0
       if(a.eq. 0.)then    ! not a 2nd deg. equation
         ierr = -1
         return
       endif

       delta= b*b - 4.d0*a*c

       if(delta .ge.0.) then    ! 1 or 2 solutions
        ierr=0
        sqrt_delta=sqrt(delta) 
        x1= (-b - sqrt_Delta) / (2.d0 * a)
        x2= (-b + sqrt_Delta) / (2.d0 * a)

       else
         ierr=-1    ! no real solutions

       endif

       end

c***************************************************************
c       Luis Peralta                                           *
c       2012                                                   *
c***************************************************************
c      compute the distance squared between points x1 and x2
       function uld2(x1,x2)
       implicit double precision(a-h,o-z), integer*4(i-n)
       real*8, dimension(3) :: x1,x2
       
       uld2=(x1(1)-x2(1))*(x1(1)-x2(1)) +
     &      (x1(2)-x2(2))*(x1(2)-x2(2)) +
     &      (x1(3)-x2(3))*(x1(3)-x2(3)) 

       end

c***************************************************
c       Ana Farinha                                *
c       Luis Peralta                               *
c       Oct. 2006                                  *
c***************************************************
c       Dot produt of two vectors
       double precision function ulvdot(x1,x2)
       implicit none
       real*8, dimension(3) :: x1,x2
       !real*8 :: ulvdot
       
       ulvdot=x1(1)*x2(1)+x1(2)*x2(2)+x1(3)*x2(3)
       
       end

c***************************************************                               *
c       Luis Peralta                               *
c       2021                                       *
c***************************************************
c      Modulus of a vector
       double precision function ulvmod(x)
       implicit none
       real*8, dimension(3) :: x
       !real*8 :: ulvmod
       
       ulvmod=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
       ulvmod=sqrt(ulvmod)
       
       end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c    Random numbers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C*********************************************************
C        Luis Peralta
c        Nov. 2013
C*********************************************************
C       generate a point inside an arbitrary volume
C        
c       generate points inside a bubble (sphere) an accept if inside volume id
c       WARNING:
c       The routine uses hit-or-miss logic, so it can be quite slow!!!
c

c       id           volume id
c       xc,yc,zc,rc  bubble center and radius
c       x,y,z        generated point coordinates
  
       subroutine ulrndvol(id,xc,yc,zc,rc,x,y,z)
       implicit double precision (a-h,o-z), integer*4 (i-n)

       external rand
       real*8, dimension(3) :: xin
       SAVE

       !pi=acos(-1.d0)
       one_third=1.d0/3.d0
       icount=0        ! loop counter 
c      generate inside a sphere os radius r

1       r3=rc*rc*rc*ulrand(1.d0)
        r=r3**one_third
        costhe=-1.d0 + 2.d0*ulrand(2.d0)
        sinthe=sqrt(abs(1.d0-costhe*costhe))
        phi=2.d0*pi*ulrand(3.d0)

        xin(1)=r*sinthe*cos(phi) + xc
        xin(2)=r*sinthe*sin(phi) + yc
        xin(3)=r*costhe          + zc

        call ulvolsearch(xin,id,iflag)
        if(iflag .eq. 0) then
         icount=icount+1
         if(mod(icount,1000000) .eq. 0) 
     &   write(*,*)"WARNING - Possible infinite loop in ULRNDVOL"
         goto 1
        endif
      
        x=xin(1)
        y=xin(2)
        z=xin(3)

        end

C*********************************************************
C        Luis Peralta
c        Dec. 2013
C*********************************************************
C       Computes the efficiency of ulrndvol random generator
C        
c       id           volume id
c       rc           bubble radius
c       eff          efficiency
  
       subroutine ulrndvoleff(id,rc,eff)
       use ullib_mod
       implicit double precision (a-h,o-z), integer*4 (i-n)
       save

      itype=mytype(id)
      if(itype.eq.0) then
       write(*,*)'ULRNDVOLEFF: Unknown volume id:',id
       stop
      endif

        !pi=acos(-1.d0)
        vsph=4.d0/3.d0 * pi * rc**3.d0 ! sphere volume

        call ulbodyvol(id,v)

        if(rc .gt. 0.)then
          eff=v/vsph
        else
          eff=0.
        endif
        end
C*********************************************************
C        Luis Peralta
c        Mar. 2008
C*********************************************************
C       Choose a case from a small list according to given
C       probabilities. 
C        
C       P = probability array
C       n = number of values
C       iout  = chosen case      
       subroutine ulrndls(P,n,iout)
       implicit double precision (a-h,o-z), integer*4 (i-n)
       real*8, dimension(n) ::  P,Pnew      
       external rand
       SAVE

        total=0.d0            ! sum of probabilities
        do i=1,n
         total=total+p(i)
        enddo

        partial=0.d0
        do i=1,n
         partial=partial+p(i)/total  ! integral probability
         Pnew(i)=partial
        enddo

       z=ulrand(1.d0)       ! generate a z~U(0,1)

       iout=1
       if(z.le.Pnew(1))goto 99
       do i=2,n
        if(z.gt.Pnew(i-1) .and. z.le.Pnew(i)) then  ! choose according to prob.
         iout=i
         goto 99
        endif
       enddo
       iout=n

99       end

C*********************************************************************
c       Luis Peralta
c       March 2008 / April 2019
C**********************************************************************
C       Generate an energy value from a given spectrum
C       with normalized intensities = spcProb(id,i)
C       Several diferent generators (with different ids are allowed.
C
C       Output values are discrete, chosen between spcX(i) values
C       
C       Routine ulspcini must have been called in the first
C       place, to get the integral probabilities.
C
C       The routine generates numbers using the inverse-
C       transform method. 

       subroutine ulspc(id,Eout)
       use ullib_mod
       implicit double precision(a-h,o-z), integer*4(i-n) 
       external rand
       save

        
       f1=ulrand(1.d0)  ! random seek value

       imin=0        ! starting i seeking values
       imax=nspc(id)
       imed=(imin+imax)/2
       
1       continue   ! seek channel dividing intervals in two halves

       if(imin.ne.0)then
         ymin=spcProb(id,imin)

       else       ! i=0 is a special case
         ymin=0.D0
       endif
       ymax=spcProb(id,imax)
       ymed=spcProb(id,imed)

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


       Eout=spcX(id,j)  ! random value output


       end

C*********************************************************************
c       Luis Peralta
c       April 2019
C**********************************************************************
c
c      ulspcu is similar to ulspc but
c      output values are randomize between specX array values
c      Intervals between spcX(i) values should be constant.
c      Non-contant interval can produce non-continuous output.
 
       subroutine ulspcu(id,xout)
       use ullib_mod
       implicit double precision(a-h,o-z), integer*4(i-n) 
       external rand
       save

       f1=ulrand(1.d0)  ! random seek value

       imin=0        ! starting i seeking values
       imax=nspc(id)
       n=imax
       imed=(imin+imax)/2
       iout=-1

       do while( iout < 0)       
! seek channel dividing intervals in two halves

       if(imin.ne.0)then
         ymin=spcProb(id,imin)
       else       ! i=0 is a special case
         ymin=0.D0
       endif

       ymax=spcProb(id,imax)
       ymed=spcProb(id,imed)

C       comparison chain

       if(f1.ge.ymin .and. f1.le.ymed)then  ! lower half
          if(imed.eq.imin+1) then
           j=imed      ! found channel value
           iout=j
           exit ! end comparison      
         endif

        imax=imed
        imed=(imin+imax)/2

       else   ! upper half

        if(imed.eq.imax-1) then
          j=imax        ! found channel value
          iout=j
          exit       ! end comparison
        endif

        imin=imed
        imed=(imin+imax)/2

       endif

       enddo ! continue comparison
      
       if(iout .ne. nspc(id)) then ! not last point
         x1=spcX(id,iout)
         x2=spcX(id,iout+1)
       else
         x1=spcX(id,iout)
         x2=spcX(id,iout-1)
       endif
         dx=abs(x2-x1)
         xout=(x1-dx/2.)+(dx)*ulrand(2.d0) ! randomize output value 


       end

c**********************************************************************
C       Luis Peralta
c       March 2008 
C**********************************************************************
C       Initialization of energy spectrum
C       id  = generator number
C       Ein = energy array        (input)
C       Pin = intensity array     (input)
C       n   = number of energy/intensity values (input) 
C
C       The probability array replaced by the integral probability
C       according to:
C       spcProb(n)=[pin(1)+pin(2)+ ....+pin(n) ] / sum pin(i) 


       subroutine ulspcini(id,Ein,Pin,n)
       use ullib_mod
       implicit double precision(a-h,o-z), integer*4(i-n)
       real*8, dimension(n) :: Pin,Ein
       save

       if(n.gt. nspc_ch)then
         stop 'ULSPCINI: array outside limits'
       endif

       if(id.gt. nspc_ger)then
         stop 'ULSPCINI: generator ID too BIG'
       endif

       if(nspc_existing(id).eq. 1)then
         stop 'ULSPCINI: existing generator ID'
       endif

       nspc_existing(id)=1

       nspc(id)=n       ! put value in common

       total=0.d0       ! comput sum of probabilities
       do i=1,n
        total=total+Pin(i)
       enddo

       partial=0.d0
       do i=1,n
        partial=partial+Pin(i)/total  
        spcProb(id,i)=partial    ! integral probability
        spcX(id,i)=Ein(i)
       enddo

       end

c---------------------------------------------------------
       subroutine ulspcres(id)
c---------------------------------------------------------
c      reset generator id in database
       use ullib_mod
       implicit none
       save
       integer*4 :: id

       nspc_existing(id)=0
       nspc(id)=0       
       end

c-----------------------------------------------------
c     set random generator seeds
c     May 2020
c-----------------------------------------------------
      subroutine ulsetseed(i1,i2)
      use ulrseed 
      implicit none
      integer*4 :: i1,i2
       iseed1=i1
       iseed2=i2
      end

c----------------------------------------------------
c     uniform random generator
c----------------------------------------------------
      FUNCTION ULRAND(DUMMY)
C
C  This is an adapted version of subroutine RANECU written by F. James
C  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
C  give a single random number at each call.
C
C  The 'seeds' ISEED1 and ISEED2 must be initialised in the main program
C  and transferred through the named common block /ULHRSEED/.
C

C  Adapted by L. Peralta @2011 for ullib
c                        @2019

      USE ULRSEED
      real(kind=8) :: ULRAND,dummy, dum
      real(kind=8), PARAMETER :: USCALE=1.0D0/2.147483563D9
      integer(kind=4) :: I1,I2,IZ

      dum=dummy
C
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
      ULRAND=IZ*USCALE
C
      RETURN
      END
