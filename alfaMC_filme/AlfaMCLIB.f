C***********************************************************************
C                            AlfaMC                                    *
C               A package for alpha particle transport                 *
C                    Luis Peralta  @ Oct 2010                          *
C                    Version 2.3   @ Dec. 2012                         *
C                    Version 3.0   @ Sep. 2013                         *
C                    Version 3.1   @ Mar. 2014                         *
C                    Version 3.2   @ Apr. 2014                         *
C                    Version 3.3   @ Jun. 2017                         *
C                    Version 4.0   @ Sep. 2017                         *
C                    Version 5.0   @ Feb. 2019                         *
c                    Version 6.0   @ Sep. 2021                         *
C                                                                      *
C            Universidade de Lisboa, Faculdade de Ciencias             *
C Laboratorio de Instrumentacao e Fisica Experimental de Particulas    *
C*************************************************************************
c
C   This program is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.
C
C   This program is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program.  If not, see http://www.gnu.org/licenses/
c

c       AlfaMC global variables module
        MODULE alfamcglobal

c       AlfaMC version 6.3
c       Luis Peralta - Copyright (c) 2020           

c
c   
            
        REAL*8, parameter :: pi=3.14159265358979d0
        REAL*8, parameter :: a0=  5.29177210903d4! Bohr radius (fm)   
        REAL*8, parameter :: alfa= 1.d0/137.035999084d0  ! fine structure constant 
        REAL*8, parameter :: hbar_c=197.32697        !hbar*c (MeV.fm)
        REAL*8, parameter :: u_to_MeV=931.454  ! MeV per u

        INTEGER*4, PARAMETER :: maxmat=50    ! max. nb of materials
        INTEGER*4, PARAMETER :: maxener=1000 ! max number of energies in database

c       number of computed energy values
c       proton (0.1 keV to 1000 MeV, 0.1 keV bins: 1E7 points)
c       alpha (0.1 keV to 100 MeV, 0.1 keV bins: 1E6 points)
c       muon (0.1 MeV to 10000 MeV, 0.1 MeV bins: 1E5 points)
        INTEGER*4,dimension(6) ::  NenergyPoints
        DATA NenergyPoints/0,0,0,1E7,1E6,1E5/
        INTEGER*4 :: npmaxE ! no. of energy points set in Alinit
        REAL*8, DIMENSION(6) :: Emax 
        DATA Emax/0.,0.,0.,100.d0,1000.d0,10000.d0/
        REAL*8, DIMENSION(6) :: Emin
        DATA Emin/0.,0.,0.,1.d-4,1.d-4,0.10d0/
        

        REAL*8, PARAMETER :: sigmastep=0.2d0 ! gaussian fluctuation on the step


c       particle
        integer :: idparticle=5 ! particle id: 4- proton, 5- alpha, 6- muon
        REAL*8  :: pmass        ! particle mass
        REAL*8  :: Zcharge      ! particle charge

c       particles mass
        REAL*8, PARAMETER :: muonmass=105.6583745d0   ! muon mass in MeV/c2
        REAL*8, PARAMETER :: protonmass=938.2720813d0 ! proton mass in MeV/c2
        REAL*8, PARAMETER :: alphamass=3727.379d0     ! alpha mass in MeV/c2

        REAL*8, PARAMETER :: emass=0.510998910d0 ! electron mass

c       particles charge
        REAL*8, PARAMETER :: Z_muon=1.d0    ! muon charge
        REAL*8, PARAMETER :: Z_proton=1.d0  ! proton charge
        REAL*8, PARAMETER :: Z_alpha=2.d0  ! alpha particle charge

c       particle step variables
        REAL*8 :: destep =0.01 ! default
        REAL*8, DIMENSION(:), allocatable :: stepmin,stepmax,Emin_part

c       header
        CHARACTER*80, DIMENSION(:), allocatable :: file_header

c       dEdx

         REAL*8, DIMENSION(:,:), allocatable :: dedxe,dedxn,dedxt
         REAL*8, DIMENSION(:), allocatable   :: Energy_table
         REAL*8, DIMENSION(:), allocatable ::Zmat,Amat
         REAL*8, DIMENSION(:), allocatable ::rho,X0mat,EnerImat

         REAL*8 :: dEdx ! dedx computed from tables in each step

c       kinematics
        REAL*8 :: beta_part,gamma_part,etotal_part,P_part,Tmax
        REAL*8 :: beta2,gamma2
        REAL*8 :: xi, MeanEloss, kappa
        REAL*8 :: z_eff  ! particle effective charge

c       Detector resolution
        REAL*8 :: adet,bdet

c      Multiple Scattering
c       iscat=-1: no scaterring, 0:Fermi model, 1: Single collision model
        INTEGER*4 :: iscat=0 ! default

c       Energy Straggling
c       istrag=-1: no straggling, istrag=0: only Gauss, istrag=1: Gauss/Vavilov/Landau
        INTEGER*4 :: istrag=0 ! default

c       Ionization/Excitation 
        integer*4 :: iresol=0 !iresol=0: no detector resolution 1: only Gauss, 2: Ionization/Excitation added

c       common input material
c       inputmat
        INTEGER*4 :: input_imat

c       Vavilov generator
        integer*4,parameter  :: nlambdavavi=535  ! max number of lambdaL bins
        integer*4,parameter  :: nkavavi=28
        integer*4,parameter  :: nbeta2vavi=11
        real*8, dimension(nkavavi,nbeta2vavi,nlambdavavi) :: datavavilov ! input from file distributions
        real*8, dimension(nkavavi,nbeta2vavi,nlambdavavi) :: xvavilov ! lambda L 
        real*8, dimension(nkavavi,nbeta2vavi,nlambdavavi) :: pvavilov ! Vavilov integral probabilities 
        integer*4, dimension(nkavavi,nbeta2vavi) :: npointvavilov ! number of lambda in Vavilov distribution
c       Landau generator
        integer*4, parameter :: nplandau=535
        real*8, dimension(nplandau) :: xlandau, plandau

        END MODULE alfamcglobal

        MODULE RAND_mod
        INTEGER*4 :: iseed1=123456789 ! Random number generator initialization
        INTEGER*4 :: iseed2=987654321
        END MODULE RAND_mod

c*-------------------------------------------------------
        SUBROUTINE Alsetparticle(id)
c*-------------------------------------------------------
        USE alfamcglobal
        IMPLICIT NONE  
        integer :: id
         idparticle=id
        end

c*-------------------------------------------------------
        SUBROUTINE Alinit(nmat,destepin,stmin,stmax,Emin_mat)
c*-------------------------------------------------------
c*      initialize variables

        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: nmat
        REAL*8 :: destepin
        REAL*8, DIMENSION(*) :: stmin,stmax,Emin_mat
        SAVE
        INTEGER*4 :: i
        
       write(*,*)' *** AlfaMC version 6.3 *** '

c       make some checks
        if(nmat.gt. maxmat) 
     &                 stop 'ERROR: Too many materials! change MAXMAT'
        if(nmat.lt. 1) stop 'ERROR: No materials to process'

        allocate(file_header(nmat))

        allocate(stepmin(nmat))
        allocate(stepmax(nmat))
        allocate(Emin_part(nmat))

        if(destepin > 0) destep=destepin
        do i=1,nmat
         stepmin(i)=stmin(i)
         stepmax(i)=stmax(i)
         Emin_part(i)=Emin_mat(i)
        enddo

c       common alphastep
        if(destep .gt. 1.d0) then
           destep=1.d0 ! set max fraction to 1
           write(*,*)'WARNING: DESTEP set to 1'
        endif

        if(destep .lt. 1.d-5) then
           destep=1.d-5 ! set min fraction
           write(*,*)'WARNING: DESTEP set to 1.D-5'
        endif

        do i=1,nmat
         if(stepmin(i).lt. 1.0D-8)then
            stepmin(i)=1.0D-8
            write(*,*)
     &     'WARNING: stepmin set to 1.D-8 cm for material ID=',i
         endif

         if(stepmax(i).gt. 1.0D+2)then
            stepmax(i)=1.0D+2
            write(*,*)
     &     'WARNING: stepmax set to 100 cm for material ID=',i
         endif


         if(Emin_part(i).lt. Emin(idparticle) ) then
          Emin_part(i)=Emin(idparticle)
          write(*,*)'WARNING: Emin set to ',Emin(idparticle),
     &    ' MeV for material ID=',i
         endif
        enddo

c       detector common
        adet=0.d0
        bdet=0.d0

c       scattering
c        iscat=0      ! init in alfamcglobal MODULE
c       straggling common
c        istrag=0     ! init in alfamcglobal MODULE

c       reset no. input materials
        input_imat=0

c       particle variables
        select case (idparticle)
          case (4) !proton
            pmass=protonmass
            Zcharge=Z_proton
            npmaxE=NenergyPoints(4)
            write(*,*)'           Set for Protons '
            write(*,*)' '
          case (5) !alpha
            pmass=alphamass
            Zcharge=Z_alpha
            npmaxE=NenergyPoints(5)
            write(*,*)'           Set for Alphas '
            write(*,*)' '
          case (6) !muon
            pmass=muonmass
            Zcharge=Z_muon
            npmaxE=NenergyPoints(6)
            write(*,*)'           Set for Muons '
            write(*,*)' '
          case default
            pmass=alphamass
            Zcharge=Z_alpha
            npmaxE=NenergyPoints(5)
            write(*,*)'           Set for Alphas '
            write(*,*)' '
        end select

c       dE/dx tables
        allocate(dedxe(nmat,npmaxE))
        allocate(dedxn(nmat,npmaxE))
        allocate(dedxt(nmat,npmaxE))
        allocate(Energy_table(npmaxE))
c       material properties
        allocate(Zmat(nmat))
        allocate(Amat(nmat))
        allocate(rho(nmat))
        allocate(X0mat(nmat))
        allocate(EnerImat(nmat))

c       straggling
         call alVavilovi
         call alLandaui

        end

c*------------------------------------------------
        subroutine Alimat(iuin,material)
c*------------------------------------------------
c        init material database
c
c        iuin:  material file unit
c        imat:  sequencial number of material to be read

        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: iuin
        CHARACTER*64 :: material
        REAL*8, DIMENSION(maxener) :: ek, stope,stopn,stopt
        SAVE
        INTEGER*4 :: i,j

        input_imat=input_imat+1    ! no. of input materials

        open(iuin,file=material,status='old') !Materials data file (input)
        i=input_imat        ! loop on the materials
         
          do j=1,maxener   ! reset
            ek(j)=0.d0
            stope(j)=0.d0
            stopn(j)=0.d0
            stopt(j)=0.d0
          enddo

          read(iuin,'(A80)',err=98,end=98)file_header(i)    ! material header
                      ! density, eff. atomic nb., eff. mass nb., Mean Exc. Energ
                
          read(iuin,*)rho(i),Zmat(i),Amat(i),X0mat(i),EnerImat(i) 

          do j=1,maxener
           read(iuin,*,end=99,err=99)
     &     ek(j),stope(j),stopn(j),stopt(j) ! energy, electronic,nuclear, total de/dx
          enddo

99        continue
           call AlmkdEdx(i,ek,stope,stopn,stopt)    ! make dEdx tables

          close(iuin)
          return

98        continue
          Stop "Empty data file"
          end

c--------------------------------------------------------------
        subroutine AlmkdEdx(imat,ek,stope,stopn,stopt) 
c--------------------------------------------------------------
c
c       make dEdx tables
c
        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: imat
        REAL*8, DIMENSION(maxener) :: ek, stope,stopn,stopt
        INTEGER*4 :: i,j
        REAL*8 :: energy, am, Ebin
        SAVE

       do i=1,npmaxE        ! reset dedx values
        dedxe(imat,i)=0.d0
        dedxn(imat,i)=0.d0
        dedxt(imat,i)=0.d0
       enddo


c       direct assignment from index to energy in 0.1 keV (proton/alpha) or 1 keV (muon)
        j=1
        Ebin=Emin(idparticle)
                        
        do i=1,npmaxE         

          !energy=(i+0.001d0)/1000.d0      ! energy in MeV . Max energy=maxmev/1000. MeV
           energy=Emin(idparticle) + (i-1)*Ebin +  Ebin/1000.d0
           Energy_table(i)=energy

          if(energy >= ek(1)) then

2          continue
           if( energy.ge. ek(j) .and. energy .le. ek(j+1) )then  ! energy is between two database values 

             if(stope(j+1) /=0. .and. stope(j) /=0.)
     &       am=log(stope(j+1)/stope(j))/log(ek(j+1)/ek(j))
             dedxe(imat,i)=stope(j)*(energy/ek(j))**am            ! logarithmic interpolation
           
             if(stopn(j+1) /=0. .and. stopn(j) /=0.)
     &       am=log(stopn(j+1)/stopn(j))/log(ek(j+1)/ek(j))
             dedxn(imat,i)=stopn(j)*(energy/ek(j))**am            ! logarithmic interpolation

             if(stopt(j+1) /=0. .and. stopt(j) /=0.)
     &       am=log(stopt(j+1)/stopt(j))/log(ek(j+1)/ek(j))
             dedxt(imat,i)=stopt(j)*(energy/ek(j))**am            ! logarithmic interpolation

           else
            j=j+1       ! update interval limits
            if(j == maxener) goto 99 ! energy greater than limit
            goto 2
           endif

          endif

        enddo
99      continue

        end

c-------------------------------------------------
        subroutine Alwrdedx(iout,imat)
c-------------------------------------------------
c       output interpolated dedx table
c
        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: iout,imat
        INTEGER*4 :: i
        SAVE

        write(iout,'(a80)')file_header(imat)
        write(iout,*)rho(imat),Zmat(imat),Amat(imat),X0mat(imat),
     &               EnerImat(imat)
        do i=1,npmaxE
         write(iout,*)
     &   energy_table(i),dedxe(imat,i),dedxn(imat,i),dedxt(imat,i)
        enddo

        end


c*------------------------------------------------
        subroutine Alkine(ek)
c*------------------------------------------------
c       compute particle kinematic variables
c       2019 version

        USE alfamcglobal
        IMPLICIT NONE  
        REAL*8 :: ek
        SAVE

          etotal_part=ek+pmass ! Particle total energy
          P_part=sqrt(etotal_part*etotal_part-pmass*pmass) !momentum
          beta_part=P_part/etotal_part
          beta2=beta_part*beta_part
          gamma_part=etotal_part/pmass
          gamma2=gamma_part*gamma_part

c         maximum energy transfer in one collision with an electron
          Tmax= 2.d0*emass*beta2 * gamma2 / 
     &   (1.d0 + 2*emass/pmass*gamma_part + (emass/pmass)**2)

        select case (idparticle)
         case(4)
         Z_eff=Zcharge
         case (5)
         Z_eff=Zcharge*(1.d0-exp(-125.d0*beta_part/2.d0**0.66666666d0))
         case (6)
         Z_eff=Zcharge
         case default
         Z_eff=Zcharge*(1.d0-exp(-125.d0*beta_part/2.d0**0.66666666d0))
        end select
        end

c-------------------------------------------------
        subroutine Alstep(mat,e,step)
c-------------------------------------------------
c       compute step
c       2021 version

        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: mat
        REAL*8 :: e,step,dxmin,dxmax
        REAL*8 :: dE,dx,sigma,dx0,de_dx
        SAVE
        REAL*8 :: alrndnor

        call alkine(e)  ! compute Lorentz beta and gamma and Tmax

        call aldedx(e,mat,de_dx) ! compute dedx
        dEdx=de_dx
        dE=destep*e  ! energy decrease according to destep= fraction of loss energy

        if(dedx .gt. 0.0d0) then
           dx=1.d0/dedx * dE ! step in cm
        else
           dx=stepmax(mat)
        endif

!       sigmastep is the flutuaction % allowed for step
        dxmax=dx+sigmastep*dx/2.
        dxmin=max(0.,dx-sigmastep*dx/2.) ! zero is the minimum value
        dx=dxmin +(dxmax-dxmin)*alrndnor(1.d0) ! uniform distribution

          ! sigma=sigmastep*dx
          ! dx0=alrndnor(1.d0)
          ! dx=dx0*sigma+dx ! apply a gaussian flutuation to step

         !step=abs(dx)

         if(dx.lt.stepmin(mat)) dx=stepmin(mat)
         if(dx.gt.stepmax(mat)) dx=stepmax(mat)
         step=dx
         MeanEloss=dedx*step ! Mean Energy loss in one step

c        energy straggling parameter
c        2*pi*Na*re^2mec^2 = 0.153537467 MeV cm2/g
         xi= 0.153537467 * rho(mat) * Zmat(mat)/Amat(mat) * z_eff**2 / 
     &          beta2 * step                 

         kappa=xi/Tmax

        end

c----------------------------------------------------
        subroutine aldedx(E,mat,dE_dx)
c----------------------------------------------------
c       compute dE/dx from data base
c       E: particle energy in MeV
c
        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: mat
        REAL*8 :: E,dE_dx
        INTEGER*4 :: i
        !REAL*8 :: Emax
        SAVE


        !i=int(E*1000.d0)
        i=int( (E-Emin(idparticle))/Emin(idparticle) ) +1
        if(i.ge.1 .and. i .le. npmaxE) then
          de_dx=dedxt(mat,i) * rho(mat)
        elseif(i .lt. 1) then
          de_dx=dedxt(mat,1) * rho(mat)
        else
          !Emax=1.d0*npmaxE/1000.d0
          write(*,'("Energy above limit! Emax:",F8.2," MeV")')
     &    Emax(idparticle)
          de_dx=dedxt(mat,npmaxE) * rho(mat)
          !stop
        endif

        end

c---------------------------------------------------
        subroutine alInter(step,mat,e,dE,ua,va,wa)
c---------------------------------------------------
c       OBSOLETE 
c       use alStragg

        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: mat
        REAL*8 :: step,e,dE,ua,va,wa,u0,v0,w0
        SAVE
          u0=ua
          v0=va
          w0=wa
          call alStragg(step,mat,e,dE)

        end
c---------------------------------------------------
        subroutine alStragg(step,mat,E,dE)
c---------------------------------------------------
c       Energy Loss and Straggling
c       2021 version

        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: mat
        REAL*8 :: step,E0,E,dE
        SAVE

         if(mat == 0) then
           dE=0.d0
           return
         endif

c     ultrack may change step
c     re-compute quantities after ultrack
         xi= 0.153537467 * rho(mat) * Zmat(mat)/Amat(mat) * z_eff**2 / 
     &          beta2 * step                 

         kappa=xi/Tmax
         MeanEloss=dedx*step
         dE=MeanEloss         ! mean energy loss
         E0=E

         select case(istrag)

          case(-1) ! no straggling
             E=E0-dE
             if(E < 0.d0)then
              dE=E0
              E=0.d0
             endif

          case(0) ! only Gauss
         
             call algauss(mat,step,dE)
             E=E0-dE
             if(E < 0.d0)then
              dE=E0
              E=0.d0
             endif

          case(1)! Gauss/Vavilov/Landau

            if(kappa > 10.d0) then       ! Gauss  kappa > 10
             call algauss(mat,step,dE)
             E=E0-dE
             if(E < 0.d0)then
              dE=E0
              E=0.d0
             endif

            elseif(kappa >= 0.01d0) then ! Vavilov 0.01 < kappa <= 10
             call alvavilov(dE)
             E=E0-dE
             if(E < 0.d0)then 
              dE=E0
              E=0.d0
             endif

            else

             call alLandau(dE) ! Landau  kappa <= 0.01
             E=E0-dE
             if(E < 0.d0)then
              dE=E0
              E=0.d0
             endif

            endif

        end select

        end

c--------------------------------------------------------
        subroutine alMScat(step,mat,ea,ua,va,wa) 
c-------------------------------------------------------- 
c       2019 version

        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: mat
        REAL*8 :: step,ea,ua,va,wa
        SAVE
        REAL*8 :: tstep

         select case(iscat)

         case(-1)  ! no multiple scattering
         
         case(0)   ! Fermi model
          call alFermiMScat(step,mat,ua,va,wa,tstep) !multiple scattering and step length correction
          step=tstep
         case(1)   ! single scattering model
          call alSingle_Col(Ea,step,mat,ua,va,wa,tstep)
          step=tstep 
         case default 
          call alFermiMScat(step,mat,ua,va,wa,tstep) !multiple scattering and step length correction
          step=tstep
         end select

         end
c----------------------------------------------------------
        subroutine alSingle_Col(E,step,mat,ua,va,wa,tstep)
c----------------------------------------------------------
        USE alfamcglobal

        IMPLICIT NONE  
        INTEGER*4 :: i,mat,Ncol
        REAL*8 :: E,step,ua,va,wa,tstep
        REAL*8 :: R,a,BigA,costhe,sinthe,Q,M
        REAL*8 :: u0,v0,w0,costheta0,theta0,phi
        REAL*8 :: cx,cy,cz
        REAL*8 :: alrand
        real*8 :: E0
        SAVE

         E0=E  ! not used
         
c        initial direction
         u0=ua
         v0=va
         w0=wa


c        compute mean energy loss per collision
         M=Amat(mat)*u_to_MeV
         Q=0.5*Tmax
         Ncol=int(dEdx*step/Q)
         if(Ncol == 0) Ncol=1 ! adicionado na versao 5.9

         R=0.885*Zmat(mat)**(-0.3333333333d0)*a0 ! R in fm
         BigA=(hbar_c/(2.d0*p_part*R))**2 * 
     &   (1.13+3.76*(alfa*Zmat(mat)/beta_part)**2)

         a=2.d0*BigA+1.d0

        do i=1,ncol
         call alwentzel(a,costhe)

         sinthe=sqrt(1.d0-costhe*costhe)
         phi=2.d0*pi*alrand(1.d0)
         cx=sinthe*cos(phi)
         cy=sinthe*sin(phi)
         cz=costhe

         call alRot(cx,cy,cz,ua,va,wa) ! rotate to Lab. ref. frame
        enddo

         costheta0=u0*ua+v0*va+w0*wa ! intern product = cos(ang)
         if(costheta0>1.d0) costheta0=1.d0
         theta0=acos(costheta0)
c        step length correction
         tstep=step + (theta0**2 / step) / 4.d0  * step*step

        end
c-------------------------------------------------------------
        subroutine alwentzel(a,x)
c-------------------------------------------------------------
c       generate a number according to Wentzel distribution
        real*8 :: y,x,a
        REAL*8 :: alrand
        y=alrand(1.d0)
        x=(2.d0*y*a-a+1.d0)/(2.d0*y+a-1.d0)

        end
c--------------------------------------------------------
        subroutine alFermiMScat(step,mat,ua,va,wa,tstep) 
c-------------------------------------------------------- 
c       Multiple scattering - Fermi model
c
        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: mat
        REAL*8 :: step,ua,va,wa,tstep
        REAL*8 :: theta0,thetax,thetay,cx,cy,cz  
        SAVE
        REAL*8 :: alrndnor
  
         call alTheta0(step,mat,theta0) ! comput the rms theta 

c        step length correction
         tstep=step + (theta0**2 / step) / 4.d0  * step*step

         thetax=alrndnor(1.d0)*theta0  ! angle in xz plane relative to z
         thetay=alrndnor(2.d0)*theta0  ! angle in yz plane relative to z

        cx=sin(thetax)! new incidence direction in particle ref. frame
        cy=sin(thetay)
        cz=sqrt(abs(1.d0-cx*cx-cy*cy))
c        c=sqrt(cx*cx+cy*cy+cz*cz)          ! renormalize cos
c        cx=cx/c
c        cy=cy/c
c        cz=cz/c

        call alRot(cx,cy,cz,ua,va,wa) ! rotate to Lab. ref. frame

        end

c--------------------------------------------------------
        subroutine altstep(mat,step,tstep) 
c-------------------------------------------------------- 
c       Compute the multiple scattering corrected step
c
        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: mat
        REAL*8 :: step,tstep
        REAL*8 :: theta0
        SAVE

         call alTheta0(step,mat,theta0) ! compute the rms theta 
c        step length correction
c        t approx =  s+ K/4 * s**2
c        the approximation theta0**2=K*s  is assumed
         tstep=step + (theta0**2 / step) / 4.d0  * step*step
         end

c--------------------------------------------------
        subroutine alRot(cx,cy,cz,u,v,w)
c--------------------------------------------------
c       rotation from beam ref.frame to Lab. ref. frame
c       use particular x,y beam axis choice. See manual!
        USE alfamcglobal
        IMPLICIT NONE  
        REAL*8 :: cx,cy,cz,u,v,w
        REAL*8 :: umod,sqruv,cxlab,cylab,czlab
        SAVE

c       renormalize u,v,w
        umod=sqrt(u*u+v*v+w*w)
        u=u/umod
        v=v/umod
        w=w/umod
        sqruv=sqrt(u*u+v*v)

        if(sqruv .ne. 0.d0) then
         cxlab= (v/sqruv) *cx + (u*w/sqruv)*cy + u*cz
         cylab= (-u/sqruv)*cx + (v*w/sqruv)*cy + v*cz
         czlab=       0.d0*cx -       sqruv*cy + w*cz

         u=cxlab   ! particle new direction
         v=cylab
         w=czlab

        else
         u=cx
         v=cy
         if(w .gt. 0.d0) then 
           w=cz
         else
           w=-cz
         endif
        endif

c       renormalize final result
        umod=sqrt(u*u+v*v+w*w)
        u=u/umod
        v=v/umod
        w=w/umod

        end

c*------------------------------------------------
        subroutine alGauss(mat,dX,dE)
c*------------------------------------------------
c       gaussian spreed of non-relativistic particles
c       2021 version

        USE alfamcglobal
        IMPLICIT NONE  
        REAL*8 :: dX,dE
        INTEGER*4 :: mat
        REAL*8 :: sigma0,sigma02,erand,dE0  
        SAVE
        REAL*8 :: alrndnor

c         sigma02= 0.1569 * rho * Z/A * zpart**2 * dX    ! Bohr

          sigma02= 0.1534d0 * rho(mat) * Zmat(mat)/Amat(mat)* 
     &              Z_eff**2 * dX / 
     &              beta2   * Tmax * (1.d0-beta2 /2.d0) ! see GEANT3 manual

          sigma0=sqrt(sigma02)
          erand=alrndnor(1.d0)
          dE0=dE
          dE=erand*sigma0+dE
          if(dE < 0.d0) dE=dE0  ! avoid negative Eloss. No straggling then.

          end
c*------------------------------------------------
        subroutine alVavilovi
c*------------------------------------------------
c       initialize Vavilov database
c       2021
        use alfamcglobal
        implicit none
        integer*4 :: n,ika,ib2,il,ib0
        integer*4 :: jka,jb2,jl
        real*8 :: ka,b2,l,p,ka0
        integer :: endof=0
        open(unit=1,file='vavilov-table.dat',status='old')

          ika=0  ! (1-28)
          il=1   ! (1-535)
          ib0=1 
          ka0=0

          xvavilov=0.  ! reset database
          datavavilov=0.  ! reset database
          
          do while(endof == 0)

          read(1,*,iostat=endof)ka,b2,l,p

          if( endof == 0) then

           if(ka > ka0) then ! see if kappa changes
             ka0=ka 
             ika=ika+1
             il=1
           endif

           ib2= int(b2*10)+1   ! (1-11)
           if(ib2 /= ib0) then ! see if beta2 changes
             ib0=ib2
             il=1
           endif

           xvavilov(ika,ib2,il)=l     ! lambda L
           datavavilov(ika,ib2,il)=p  ! probability
           il=il+1
           endif

           enddo

        npointvavilov=0  ! reset array

!       find number of non-zero probability values
        do ika=1,nkavavi
         do ib2=1,nbeta2vavi
          do il=1,nlambdavavi
           if(datavavilov(ika,ib2,il) > 0.) then
             npointvavilov(ika,ib2)=npointvavilov(ika,ib2)+1
           else
            exit
           endif
          enddo
         enddo
        enddo

        do ika=1,nkavavi
         do ib2=1,nbeta2vavi

          n=npointvavilov(ika,ib2)
          call algeni(datavavilov(ika,ib2,:),pvavilov(ika,ib2,:),n)

         enddo
        enddo
 

        close(1)
        end

c*------------------------------------------------
        subroutine alVavilov(dE)
c*------------------------------------------------
c       generate lambdaL according to Vavilov Distribution
c       2021
        use alfamcglobal
        implicit none
        real*8 :: dE
        integer*4 :: ika,ib2,n
        real*8 :: lambda
        real*8, parameter :: C=0.577215d0
        save


        if(kappa >=0.01d0 .and. kappa <0.1d0)then
          ika = int(kappa*100)
        elseif(kappa >=0.1d0 .and. kappa <1.d0)then
          ika=int(kappa*10)+9
        elseif(kappa >=1.d0 .and. kappa <=10.d0)then
          ika=int(kappa)+18
        else
          print *,'ERROR calling alVavilov'
          dE=0.
          return
        endif
    

        ib2= int(beta2*10+0.4999d0)+1 ! if beta2>0.95 use table for beta2=1
        n=npointvavilov(ika,ib2)
        call algen(Xvavilov(ika,ib2,:),Pvavilov(ika,ib2,:),lambda,n)

        dE=(lambda + (1.d0-C) + beta2 + log(kappa))*xi + MeanEloss


        if(dE < 0.d0)  dE=0. ! avoid negative Eloss. No straggling then.

        end

c*------------------------------------------------
        subroutine allandaui
c*------------------------------------------------
c        init Landau random generator
c        2021 version
         use alfamcglobal
         implicit none
         integer :: n
         real*8, dimension(nplandau) :: datalandau
         real*8 :: x0,step

       data datalandau/7.15d-06,2.17d-05,5.89d-05,1.45d-04,3.25d-04,
     & 0.000673,0.001295,0.002329,0.003943,0.006319,0.009637,0.014054,
     & 0.019682,0.026572,0.034703,0.043985,0.054259,0.065311,0.076892,
     & 0.088729,0.100551,0.112097,0.123135,0.133463,0.142922,0.151392,
     & 0.158793,0.165084,0.170256,0.174325,0.177334,0.17934,0.180414,
     & 0.180635,0.180087,0.178873,0.177032,0.174674,0.171875,0.168706,
     & 0.165234,0.161517,0.157611,0.153562,0.149415,0.145207,0.140969,
     & 0.136729,0.132512,0.128337,0.124221,0.120177,0.116217,0.112348,
     & 0.108579,0.104913,0.101355,0.097906,0.094568,0.091341,0.088224,
     & 0.085217,0.082318,0.079525,0.076836,0.074248,0.071758,0.069363,
     & 0.06706,0.064847,0.06272,0.060676,0.058712,0.056824,0.055011,
     & 0.053269,0.051594,0.049986,0.048439,0.046953,0.045525,0.044151,
     & 0.042831,0.041561,0.040339,0.039163,0.038032,0.036944,0.035895,
     & 0.034886,0.033914,0.032978,0.032075,0.031206,0.030367,0.029559,
     & 0.028779,0.028026,0.0273,0.0266,0.025923,0.025269,0.024638,
     & 0.024028,0.023438,0.022868,0.022317,0.021784,0.021269,0.02077,
     & 0.020286,0.019819,0.019366,0.018927,0.018502,0.01809,0.017691,
     & 0.017304,0.016928,0.016564,0.016211,0.015868,0.015535,0.015212,
     & 0.014898,0.014594,0.014298,0.01401,0.01373,0.013458,0.013194,
     & 0.012937,0.012687,0.012444,0.012207,0.011976,0.011752,0.011534,
     & 0.011321,0.011114,0.010912,0.010715,0.010524,0.010337,0.010155,
     & 0.009977,0.009804,0.009635,0.00947,0.009309,0.009152,0.008999,
     & 0.00885,0.008704,0.008561,0.008422,0.008286,0.008153,0.008023,
     & 0.007896,0.007772,0.007651,0.007532,0.007416,0.007303,0.007192,
     & 0.007083,0.006977,0.006873,0.006771,0.006672,0.006574,0.006479,
     & 0.006385,0.006294,0.006204,0.006116,0.00603,0.005946,0.005864,
     & 0.005783,0.005703,0.005626,0.005549,0.005475,0.005401,0.005329,
     & 0.005259,0.00519,0.005122,0.005055,0.00499,0.004926,0.004863,
     & 0.004801,0.004741,0.004681,0.004623,0.004565,0.004509,0.004454,
     & 0.004399,0.004346,0.004293,0.004242,0.004191,0.004141,0.004092,
     & 0.004044,0.003997,0.003951,0.003905,0.00386,0.003816,0.003772,
     & 0.00373,0.003688,0.003646,0.003606,0.003566,0.003526,0.003488,
     & 0.00345,0.003412,0.003375,0.003339,0.003303,0.003268,0.003233,
     & 0.003199,0.003165,0.003132,0.0031,0.003068,0.003036,0.003005,
     & 0.002974,0.002944,0.002914,0.002885,0.002856,0.002828,0.0028,
     & 0.002772,0.002745,0.002718,0.002691,0.002665,0.00264,0.002614,
     & 0.002589,0.002565,0.00254,0.002516,0.002493,0.002469,0.002446,
     & 0.002424,0.002401,0.002379,0.002357,0.002336,0.002315,0.002294,
     & 0.002273,0.002253,0.002232,0.002213,0.002193,0.002174,0.002155,
     & 0.002136,0.002117,0.002099,0.002081,0.002063,0.002045,0.002028,
     & 0.00201,0.001993,0.001977,0.00196,0.001944,0.001927,0.001911,
     & 0.001896,0.00188,0.001865,0.001849,0.001834,0.001819,0.001805,
     & 0.00179,0.001776,0.001762,0.001748,0.001734,0.00172,0.001707,
     & 0.001693,0.00168,0.001667,0.001654,0.001641,0.001629,0.001616,
     & 0.001604,0.001592,0.001579,0.001568,0.001556,0.001544,0.001533,
     & 0.001521,0.00151,0.001499,0.001488,0.001477,0.001466,0.001455,
     & 0.001445,0.001434,0.001424,0.001414,0.001403,0.001393,0.001383,
     & 0.001374,0.001364,0.001354,0.001345,0.001335,0.001326,0.001317,
     & 0.001308,0.001299,0.00129,0.001281,0.001272,0.001263,0.001255,
     & 0.001246,0.001238,0.00123,0.001221,0.001213,0.001205,0.001197,
     & 0.001189,0.001181,0.001174,0.001166,0.001158,0.001151,0.001143,
     & 0.001136,0.001128,0.001121,0.001114,0.001107,0.0011,0.001093,
     & 0.001086,0.001079,0.001072,0.001065,0.001059,0.001052,0.001046,
     & 0.001039,0.001033,0.001026,0.00102,0.001014,0.001007,0.001001,
     & 0.000995,0.000989,0.000983,0.000977,0.000971,0.000966,0.00096,
     & 0.000954,0.000948,0.000943,0.000937,0.000932,0.000926,0.000921,
     & 0.000915,0.00091,0.000905,0.0009,0.000894,0.000889,0.000884,
     & 0.000879,0.000874,0.000869,0.000864,0.000859,0.000854,0.00085,
     & 0.000845,0.00084,0.000835,0.000831,0.000826,0.000821,0.000817,
     & 0.000812,0.000808,0.000804,0.000799,0.000795,0.00079,0.000786,
     & 0.000782,0.000778,0.000773,0.000769,0.000765,0.000761,0.000757,
     & 0.000753,0.000749,0.000745,0.000741,0.000737,0.000733,0.000729,
     & 0.000726,0.000722,0.000718,0.000714,0.000711,0.000707,0.000703,
     & 0.0007,0.000696,0.000693,0.000689,0.000686,0.000682,0.000679,
     & 0.000675,0.000672,0.000669,0.000665,0.000662,0.000659,0.000655,
     & 0.000652,0.000649,0.000646,0.000642,0.000639,0.000636,0.000633,
     & 0.00063,0.000627,0.000624,0.000621,0.000618,0.000615,0.000612,
     & 0.000609,0.000606,0.000603,0.0006,0.000597,0.000595,0.000592,
     & 0.000589,0.000586,0.000583,0.000581,0.000578,0.000575,0.000573,
     & 0.00057,0.000567,0.000565,0.000562,0.000559,0.000557,0.000554,
     & 0.000552,0.000549,0.000547,0.000544,0.000542,0.000539,0.000537,
     & 0.000535,0.000532,0.00053,0.000527,0.000525,0.000523,0.00052,
     & 0.000518,0.000516,0.000513,0.000511,0.000509,0.000507,0.000504,
     & 0.000502,0.0005,0.000498,0.000496,0.000494,0.000491,0.000489,
     & 0.000487,0.000485,0.000483,0.000481,0.000479,0.000477,0.000475,
     & 0.000473,0.000471,0.000469,0.000467,0.000465,0.000463,0.000461,
     & 0.000459,0.000457,0.000455,0.000453,0.000452/
!         book Landau generator 
          step=0.1d0
          x0=-3.5d0
          do n=1,nplandau
           xlandau(n)= x0+(n-1)*step
          enddo
          !call ulspcini(1,xlandau,plandau,nplandau) ! init Landau generator
           call algeni(datalandau,plandau,nplandau)  ! init Landau generator
        end

c*------------------------------------------------
        subroutine alLandau(dE)
c*------------------------------------------------
c       Energy loss Landau model
c       2019 version

        USE alfamcglobal
        IMPLICIT NONE  
        REAL*8 :: dE  
        REAL*8 :: lambda, C
        data C/0.577215d0/
        SAVE


       ! call ulspcu(1,lambda)
        call algen(Xlandau,Plandau,lambda,nplandau)

        dE=(lambda + (1.d0-C) + beta2 + log(kappa))*xi + MeanEloss

        if(dE < 0.d0)  dE=0. ! avoid negative Eloss. No straggling then.

        end

c---------------------------------------------
         subroutine alTheta0(step,mat,theta0)
c---------------------------------------------
c       compute rms theta0 angle
        USE alfamcglobal
        IMPLICIT NONE  
        REAL*8 :: step,theta0
        INTEGER*4 :: mat
        REAL*8 :: X0
        SAVE

        X0=X0mat(mat)/rho(mat) ! X0mat in g cm-2   X0 in cm

c       Lynch and Dahl ! 4/5 is a correction for better comparison with SRIM
        theta0=4.d0/5.d0*13.6d0 /(beta_part*P_part)*2.d0*sqrt(step/X0)*  
     &  (1.d0+0.038d0*log((step*Z_eff**2)/(X0*beta2)))        
                       

        end

c*------------------------------------------------
        subroutine alSetCol(is)
c*------------------------------------------------
c       set multiscattering model
c       iscat= -1: no multiple-scattering 0: Fermi model 1: single collision model 
        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: is
        SAVE

        select case(is)
        case(-1)
         iscat=-1 ! no scattering
        case(0)
         iscat=0 ! Fermi model 
        case(1)
         iscat=1 !single collision model
        case  default
         iscat=0! Fermi
        end select
        end


c*------------------------------------------------
        subroutine alSetStrag(is)
c*------------------------------------------------
c       set type of straggling functions
c       istrag=0: only Gauss, istrag=1: Gauss/Vavilov/Landau
        USE alfamcglobal
        IMPLICIT NONE  
        INTEGER*4 :: is
        SAVE

        select case(is)
        case(-1)
         istrag=-1! no straggling 
        case(0)
         istrag=0! Gauss 
        case(1)
         istrag=1
        case  default
         istrag=0! Gauss
        end select
        end

c*------------------------------------------------
        subroutine alsetdet(a,b)
c*------------------------------------------------
c       set detector resolution parameters
c*      detector resolution computed as r=a/sqrt(E)+b

        USE alfamcglobal
        IMPLICIT NONE  
        REAL*8 :: a,b
        SAVE
          if( a > 0.)then
            adet=a
            bdet=b
          else
           adet=0.
           bdet=0.
          endif

        end

c*------------------------------------------------
        subroutine alsetRes(ir,a,b)
c*------------------------------------------------
        use alfamcglobal
        implicit none
        integer*4 :: ir   
        real*8 :: a,b
        save

c       set ionization/excitation function for energy resolution 
c       iresol=0: no detector resolution 1: only Gauss, 2: Ionization/Excitation added

          iresol=ir
          if( a > 0.)then
            adet=a
            bdet=b
          else
           adet=0.
           bdet=0.
          endif
        end
c*------------------------------------------------
        subroutine alsigmadet(E)
c*------------------------------------------------
c*      detector resolution
c       OBSOLETE : use alDetRes(E)
        USE alfamcglobal
        IMPLICIT NONE  
        REAL*8 :: E, Ein,sigma,erand
        real*8 :: alrndnor
        external alrandnor
        SAVE

         Ein=E
         if(adet > 0.) then
c        detector resolution computed as r=a/sqrt(E)+b
          sigma = adet*sqrt(E)+bdet*E    
          erand=alrndnor(1.d0)
          E=erand*sigma+E  
          if(E < 0.) E=0.

         endif

        end

c---------------------------------------------------
        subroutine alfunIE(Energy,Eout)
c---------------------------------------------------
c Initialize generator
c Ionization/Excitation function for silicon
c M.Jurado Vargas, Alfonso Fernandez Timon, J.F. Ziegler
c doi.org/10.1016/j.nima.2020.164134

       implicit none
       real*8 :: a,b,c,a1,b1,g1,a2,b2,g2,eta
       real*8 :: sigma, tau1,tau2
       real*8 :: Emin,Emax,dE
       real*8, parameter :: dE_def=0.1 ! Energy steps in keV for generator (default)
       real*8 :: Energy,E0,f1,f2,f,u,Eout

       integer   :: i
       integer*4 :: n  ! generator number of channels
       integer*4, parameter :: nmax=1000
       real*8, dimension(nmax) :: Ein,Pin,Pout

       Ein=0.d0
       Pin=0.d0
       E0= Energy*1000.d0  ! routine uses keV, while Energy in MeV
c sigma
       a=7.12d-1
       b=1.32d-3
       c=6.13d-1

c tau1
       a1=2.36d-1
       b1=1.56d-1
       g1=2.74d-4

c tau2
       a2=4.83E-2
       b2=5.60E-2
       g2=4.22E-4

c eta mean value
       eta=0.909 

         sigma=a+b*E0**c
         tau1=a1+b1*exp(-g1*E0)
         tau2=a2+b2*exp(-g2*E0)

         Emin=E0-20.*sigma
         if(Emin < 0.) Emin=0.
         Emax=E0+5.*sigma


         n=int((Emax-Emin)/dE_def)
         if( n <= nmax)then
            dE=dE_def
         else
            dE=(Emax-Emin)/(1.d0*nmax)
         endif

        do i=1,n

         u=Emin + (i-1)*dE

         f1=(     eta*tau1*exp((u-E0)*tau1 + sigma**2*tau1**2/2 ))*
     &           erfc(  1/sqrt(2.) *( (u-E0)/sigma +sigma*tau1 ))
         f2=( (1-eta)*tau2*exp((u-E0)*tau2 + sigma**2*tau2**2/2 ))*
     &           erfc(  1/sqrt(2.) *( (u-E0)/sigma +sigma*tau2 ))

         !f= sigma*sqrt(2.*pi)/2. *(f1+f2)
          f=f1+f2

          Ein(i)=u/1000. ! back to MeV
          Pin(i)=f
        enddo

        call algeni(Pin,Pout,n)
        call algen(Ein,Pout,Eout,n)

        end

c---------------------------------------------------
        subroutine algausstail(Energy,Eout)
c---------------------------------------------------

       implicit none

       real*8 :: sigma, tau1,tau2,eta
       real*8 :: Emin,Emax,dE
       real*8 :: Energy,E0,f1,f2,f,u,Eout

       integer   :: i
       integer*4 :: n  ! generator number of channels
       integer*4, parameter :: nmax=1000
       real*8, dimension(nmax) :: Ein,Pin,Pout

         E0=Energy
         sigma=12.E-3
         tau1=20.E-3
         tau2=30.E-3
         eta=0.9d0

         n=400
        ! Emin=7.4d0
        ! Emax=7.8d0
         Emin=E0-20.*sigma
         if(Emin < 0.) Emin=0.
         Emax=E0+10.*sigma
         dE=(Emax-Emin)/(1.d0*n)
         

        do i=1,n

         u=Emin + (i-1)*dE + dE/2.

         f1=(     eta/tau1*exp((u-E0)/tau1 + sigma**2/tau1**2/2 ))*
     &           erfc(  1/sqrt(2.) *( (u-E0)/sigma +sigma/tau1 ))
         f2=( (1-eta)/tau2*exp((u-E0)/tau2 + sigma**2/tau2**2/2 ))*
     &           erfc(  1/sqrt(2.) *( (u-E0)/sigma +sigma/tau2 ))

          f=f1+f2

          Ein(i)=u
          Pin(i)=f
        enddo

        call algeni(Pin,Pout,n)
        call algen(Ein,Pout,Eout,n)

        end

c*------------------------------------------------
        subroutine alDetRes(E,Eout)
c*------------------------------------------------
c*      detector resolution
        use alfamcglobal
        implicit none 
        real*8 :: E, Eout,sigma, erand
        save
        real*8 :: alrndnor
        external alrndnor
c       iresol=0: no detector resolution 1: only Gauss, 2: Ionization/Excitation added

        Eout=E
        if(iresol == 0) return
        if( adet <= 0.) return

        if( iresol == 4) then  ! gauss + tail
         call algausstail(E,Eout)
        endif

        if( iresol == 2 .or. iresol == 3) then  ! add ionization/excitation
         call alfunIE(E,Eout)
        endif

        if(iresol == 1 .or. iresol == 2) then
c        detector resolution computed as r=a/sqrt(E)+b
          sigma = adet*sqrt(Eout)+bdet*Eout    
          erand=alrndnor(1.d0)
          Eout=erand*sigma+Eout  
        endif

        if(E <=0.) then
           Eout=0.
           return
        endif

        end

c*------------------------------------------------
        double precision function alrndnor(dummy)
c*------------------------------------------------
c*      Normal distributed random numbers
        IMPLICIT NONE
        REAL*8 :: alrand
        external alrand
        REAL*8 :: dum,dummy,u1,u2
        REAL*8 :: pi=3.14159265358979d0
        SAVE

         dum=dummy
         u1=alrand(1.d0)
         u2=alrand(2.d0)
         alrndnor = sin(2.d0*pi*u1)*sqrt(-2.d0*log(u2))
        end

c-----------------------------------------------------
c     set random generator seeds
c     May 2020
c-----------------------------------------------------
      subroutine alsetseed(i1,i2)
      use RAND_mod
      implicit none
      integer*4 :: i1,i2
       iseed1=i1
       iseed2=i2
      end


c**********************************************************************
C       Luis Peralta
c       Oct. 2021 
C**********************************************************************
C       Generator initialization
CC 
C       Adapted from ULLIB

       subroutine algeni(Pin,Pout,n)
       use alfamcglobal
       implicit none 
       real*8, dimension(n) :: Pin,Pout
       real*8 :: total,partial
       integer*4 :: i,n

       total=0.d0       ! comput sum of probabilities
       Pout=0.d0        ! reset
       do i=1,n
        total=total+Pin(i)
       enddo

       partial=0.d0
       do i=1,n
        partial=partial+Pin(i)/total  
        Pout(i)=partial    ! integral probability
       enddo

       end


C*********************************************************************
c       Luis Peralta
c       Oct. 2021
C**********************************************************************
C      Random number according to given distribution
C
C      Adapted from ULLIB
       subroutine algen(Xin,P,Xout,n)
       implicit none
       real*8, dimension(n) :: Xin,P ! arrays containing values and probabilities
       real*8 :: Xout ! output random value
       integer*4 :: n
       real*8 :: f1,ymin,ymax,ymed,x1,x2,dx
       integer*4 :: imin,imax,imed,iout,j
       real*8 :: alrand
       save

       f1=alrand(1.d0)  ! random seek value

       imin=0        ! starting i seeking values
       imax=n
       imed=(imin+imax)/2
       iout=-1

       do while( iout < 0)       
! seek channel dividing intervals in two halves

       if(imin.ne.0)then
         ymin=P(imin)
       else       ! i=0 is a special case
         ymin=0.D0
       endif

       ymax=P(imax)
       ymed=P(imed)

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
      
       if(iout < n) then ! not last point
         x1=Xin(iout)
         x2=Xin(iout+1)
       else
         x1=Xin(iout)
         x2=Xin(iout-1)
       endif
         dx=abs(x2-x1)
         Xout=(x1-dx/2.)+(dx)*alrand(2.d0) ! randomize output value 

       end

c---------------------------------------------------
        double precision FUNCTION alRAND(DUMMY)
c---------------------------------------------------

C
C  This is an adapted version of subroutine RANECU written by F. James
C  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
C  give a single random number at each call.
C
C  The 'seeds' ISEED1 and ISEED2 are initialised in MODULE RAND_mod. 
C  It is advisable to declare alRAND as an external function in all sub-
C  programs that call it.
C
C     modified 2019

      use RAND_mod

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (USCALE=1.0D0/2.147483563D9)
      real*8 :: dummy,dum
      
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
      alRAND=IZ*USCALE
C
      RETURN
      END

