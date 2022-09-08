! Computation of film thickness
!       units: cm, MeV, degree
        module geometry_mod !--------------------------------------
        ! materials
          integer, parameter :: nmaterial=4 ! no. of materials

        ! cylindrical vacuum chamber
          real*8, parameter :: rcha=1.     ! radius
          real*8, parameter :: hcha=10.     ! height
          real*8, parameter :: zwcha=0.d0  ! chamber bottom
          real*8, parameter :: xcha=0.d0    ! x,y,z of chamber center
          real*8, parameter :: ycha=0.d0
          real*8, parameter :: zcha=zwcha+hcha/2.d0
          integer           :: idcha        ! chamber volume id

        ! cylindrical target/film
          real*8, parameter :: rtar=1.d0    ! radius
          real*8, parameter :: htar=50.d-7  ! thickness
          real*8, parameter :: zwtar=8.0     ! z of target entrance
          real*8, parameter :: xtar=0.      ! x,y,z of target center
          real*8, parameter :: ytar=0. 
          real*8, parameter :: ztar=zwtar+htar/2.
          integer           :: idtar        ! target volume id
          real*8, parameter :: sigmatar=0.  ! target thickness fluctuation in cm


        ! cylindrical gold detector window
          real*8, parameter :: rwin=0.57/2.d0 ! radius
          real*8, parameter :: hwin=2.07d-6   ! thickness ~ 40 ug/cm2 rho=19.32 g/cm3
          real*8, parameter :: zwwin=8.5     ! z of window entrance
          real*8, parameter :: xwin=0.d0      ! x,y,z of target center
          real*8, parameter :: ywin=0.d0
          real*8, parameter :: zwin=zwwin+hwin/2.
          integer           :: idwin        ! target volume id

        ! cylindrical Si detector
          real*8, parameter :: rdet=0.57/2.d0  ! radius
          real*8, parameter :: hdet=300.d-4  ! detector thickness/height
          real*8, parameter :: zwdet=zwwin+hwin       ! z of detector window
          real*8, parameter :: xdet=0.       ! x,y,z of detector center
          real*8, parameter :: ydet=0.       ! detector must be completly inside chamber
          real*8, parameter :: zdet=zwdet+hdet/2.
          integer           :: iddet         ! detector volume id
        save
        end module geometry_mod 

        module source_mod !---------------------------------------
        ! cylindrical source
          real*8, parameter :: rsrc=0.25d0     ! radius
          real*8, parameter :: hsrc=1.0d-5        ! height
          real*8, parameter :: xsrc=0.
          real*8, parameter :: ysrc=0.
          real*8, parameter :: zsrc=7.5
          real*8, parameter :: thetabeam=0.   ! theta angle for pencil beam / conic beam
          real*8, parameter :: phibeam=0.     ! phi angle for pencil beam
          integer,parameter :: idparticle=5   ! 5: alpha 
          integer,parameter :: isrcspace= 2   ! 1:point, 2:flat, 3:cylindrical 
          integer,parameter :: isrcdir=3      ! 1:pencil, 2:isotropic, 3:window, 4:conic beam  
          integer,parameter :: isrce=3        ! 1:mono, 2:some peaks, 3:spectrum 
          real*8,parameter  :: Epeak=5.48     ! E in MeV if monoenergetic
          integer*4 :: npeak                  ! no. of peaks/bins (isrce=2,3)
          real*8,dimension(:), allocatable :: Esrc    ! energy array (peak/spectrum) in keV
          real*8,dimension(:), allocatable :: Psrc    ! peak/spectrum intensity array
          character*64      :: src_file='U232_Calibrated_v1.in' ! peak/spectrum file name
          integer           :: idsrc          ! mc generator spectrum id
        save
        end module source_mod 

        module alfaMC_mod
         integer*4 :: iscol=0 ! scattering -1: no scaterring, 0:Fermi, 1: Single collision
         integer*4 :: istra=0 ! straggling -1: no stragglig, 0: Gauss, 1: Gauss/Vavilov/Landau

         integer*4 :: ntot= 1E7 ! Total number of events
         integer*4 :: ninfo=5E5 ! Information loop

         ! detector resolution R=sigma/E = FWHM/(2.355*E)
         real*8,parameter  :: adet=0.0005 !parameters for detector resolution
         real*8,parameter  :: bdet=0.0000   !R=a/sqrt(E)+b  - if adet=0 no resolution

         ! chamber step
         real*8,parameter  :: stepmincha= 100.d-4 !chamber min step in cm
         real*8,parameter  :: stepmaxcha= 1000.d-4 !chamber max step in cm

         ! target/film step
         real*8,parameter  :: stepmintar= 0.01d-4 !target min step in cm
         real*8,parameter  :: stepmaxtar= 9.d-4 !target max step in cm

         ! Gold window
         real*8,parameter  :: stepminwin= 1.0d-8 !chamber max step in cm
         real*8,parameter  :: stepmaxwin= 100.d-8 !target max step in cm

         ! detector step
         real*8,parameter  :: stepmindet= 1.0d-4 !detector min step in cm
         real*8,parameter  :: stepmaxdet= 50.d-4 !detector max step in cm

        save
        end module alfaMC_mod

        module myhistos
         ! histogramming

         ! Energy in detector
         integer*4, parameter :: nbindet=1000
         real*8,parameter     :: Ehdet=9.d0 ! Emax histogram in MeV

         ! Energy in target
         integer*4, parameter :: nbintar=1000
         real*8,parameter     :: Ehtar=10.d0 ! Emax histogram in MeV

        save
        end module myhistos


