C
C      	AlfaMC 
c       Main for cylindrical chamber+target+window+detector 
C
	PROGRAM MAIN

        use geometry_mod
        use source_mod
        use alfaMC_mod
        use myhistos

	IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N) 	
	character*64 tit,name_file,material,comment
        real*8 stepmin(20),stepmax(20),Emin_Alpha(20)
        real*8 tfoil
        integer :: endof
	!COMMON/RSEED/ISEED1,ISEED2 ! random number generator
C-----------------------------------------------------------------
	!iseed1=123456789 ! Random number generator initialization
	!iseed2=987654321
c-------------------------------------------------------
	call ulinit      ! start Ulysses

C       Histogram booking
C-------------------------------------------------------
c	 Start histo database. 
c        ulhinit: nh1,nmx1,nh2,nmx2,nmy2,nh3,nmx3,nmy3,nmz3
         call ulhinit(10,10000,10,100,100,0,0,0,0)  

         idhdet=1
         idhtar=2
         idhztar=3 

        tit='E in detector (MeV)'                
        call ulhb1(idhdet,tit,nbindet,0.D0,Ehdet) 

        tit='E in target (MeV)'                
        call ulhb1(idhtar,tit,nbintar,0.D0,Ehtar) 

c ------ input energy spectrum

         select case(isrce)
          case(2,3)
           open(unit=1,file=src_file,status='old')
           npeak=0
           endof=0
           do while(endof==0) ! find npeak
             read(1,*,IOSTAT=endof)dum1,dum2
             if(endof == 0) npeak=npeak+1
           enddo
           close(1)
           allocate(Esrc(npeak))
           allocate(Psrc(npeak))

           open(unit=1,file=src_file,status='old')
             do i=1,npeak   ! read energy, weight
              read(1,*)Esrc(i),Psrc(i)
             enddo
           close(1)

           idsrc=100
           call ulspcini(idsrc,Esrc,Psrc,npeak) ! init spectrum generator

          case default
          end select


C     Input AlfaMC parameters
C***************************************************************
c	ntot =1E4  ! number of events (alpha particles)
c        ninfo=1E3  ! info frequency

c       initialize AlfaMC

        nmat=nmaterial ! no. of materials

c       cuts
	destep=0.01 ! fraction energy loss in each step
        ! Inside AlfaMC energy in MeV 
	do i=1,nmat
          Emin_Alpha(i)= 1.0d-2  ! Minimum Energy in MeV
	enddo

          ! chamber
          stepmin(1)= stepmincha
          stepmax(1)= stepmaxcha
          ! target
          stepmin(2)= stepmintar
          stepmax(2)= stepmaxtar
          ! Gold window
          stepmin(3)= stepminwin
          stepmax(3)= stepmaxwin
          ! Si detector
          stepmin(4)= stepmindet
          stepmax(4)= stepmaxdet

        call Alinit(nmat,destep,stepmin,stepmax,Emin_Alpha)   ! initialize AlfaMC global variables
c -------detector resolution 
         if( adet .ne. 0) then
           call alsetdet(adet,bdet) 
         endif
        call alSetStrag(istra)! set straggling -1: no stragglig, 0: Gauss, 1: Gauss/Landau
        call alSetCol(iscol)  ! set scattering -1: no scaterring, 0:Fermi, 1: Single collision

c       input materials
        open(unit=2,file='material.in',status='old')

         iuin=1 ! input unit

         do i=1,nmat
          read(2,*)material, comment
          call Alimat(iuin,material)   ! input material
         enddo
        close(2)

c 	Define geometry for ULYSSES
C**********************************************************************
	write(*,*)'Reading Geometry ----------------------------------'	
	call ulgeo
!	call ulgsummary ! output geometry summary


C 	Start shower simulation
C**********************************************************************
	write(*,*)'start sending particles----------------------------'	

	!N=0
        do n=1,NTOT
10      continue   

	!N=N+1  ! nb of generated events
	if( mod(n,ninfo).eq.0) write(*,*)'Event: ',N ! user info

C       reset counters   
        Edet=0.   ! Energy in medium
        Etar=0.
       
C	Generates the primary particles
C***************************************************************
c       changing foil thickness
        if(sigmatar .gt.0) then
          rnd=alrndnor(1.d0) ! normal distribution N(0,1)
          tfoil=rnd*sigmatar+htar ! variable foil thickness 
          call ulput1par(idtar,2,tfoil) ! change foil thickness
        endif

19	call ulsource(ea,xa,ya,za,ua,va,wa,kpar)
        call ulfind(xa,ya,za,ibody,mat)
	if(ea .lt. Emin_alpha(mat)) goto 19  ! particle below cuts. 

20     continue
        call alstep(mat,ea,step) !Determines segment length
        mat0=mat                ! keep memory of current material

        call alMScat(step,mat,ea,ua,va,wa) ! multiple scattering

C	Track the particle through the geometry
        call ultrack(xa,ya,za,ua,va,wa,ibody,step,newbody,mat) ! a new step is found if a boundary is crossed

        call alInter(step,mat0,ea,dE,ua,va,wa) ! Energy loss and multiple scattering

C        Scoring
C*********************************************************
         if(ibody == iddet .or. ibody == idwin) Edet=Edet+dE ! Total energy in MeV
         if(ibody == idtar) Etar=Etar+dE ! Total energy in MeV
C**********************************************************

        ibody=newbody ! set the new body
	IF(MAT.EQ.0 .OR. IBODY.EQ.1) THEN !The particle left the material system
	  GOTO 40 !Exit
	ENDIF

	IF(Ea.LT.Emin_alpha(mat)) THEN !The particle has been absorbed
	     GOTO 40 !Exit
	ENDIF
	GOTO 20   ! another step

C>>>>>>>>>>>>>>>> The simulation of the track ends here.
40      CONTINUE

c       close average histograms     

c       filling histograms at end of event ---------------------
        if(Edet.ne.0) then 
          if(adet .ne. 0.) call alsigmadet(Edet)   ! resolucao do detector
          call ulhadd1(idhdet,Edet,1.d0) ! Energy in detector
        endif

        if(Etar.ne.0) call ulhadd1(idhtar,Etar,1.d0) ! Energy in target
c-----------------------------------------------------------

        enddo

c*******************************************************************
	call ulhend(0)  ! close histograms for filling
                               
c	Output histograms
C*************************************************************	
	name_file='summary.txt'
	call ulhopen(90,name_file,'s') !output file for histo 
	call ulhlist(0)   ! all histograms summary

	name_file='Edet.csv'
	call ulhopen(98,name_file,'o') !output file for histo 
	call ulhout(idhdet)    ! output histo 

	name_file='Etar.csv'
	call ulhopen(98,name_file,'o') !output file for histo 
	call ulhout(idhtar)    ! output histo  


	call ulhclose ! close all output units	

	END
C...+....1....+....2....+....3....+....4....+....5....+....6....+....7..
