
c**************************************************************
        subroutine ulsource(e,x,y,z,u,v,w,kpar) 
        use geometry_mod
        use source_mod
       
        IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N) 
        real*8, parameter :: pi=3.14159265358979323846
        external alrand

C       set particle type - 1-electrons;2-photons;3-positrons 
        kpar=idparticle


c       select space distribution ----------------------------
        select case(isrcspace) 
         case(1)
c        point source
          x=xsrc
          y=ysrc
          z=zsrc

         case(2)
c        flat circular source
          r2=rsrc**2 * alrand(1.d0)
          phi=2.*pi*alrand(2.d0)
          x=sqrt(r2)*cos(phi)+xsrc
          y=sqrt(r2)*sin(phi)+ysrc
          z=zsrc

         case(3)
c        cylindrical source
          r2=rsrc**2 *alrand(1.d0)
          phi=2.*pi*alrand(2.d0)
          x=sqrt(r2)*cos(phi)+xsrc
          y=sqrt(r2)*sin(phi)+ysrc
          z=hsrc*alrand(3.d0)+zsrc

         case default
c        point source
          x=xsrc
          y=ysrc
          z=zsrc
        end select


c       select direction distribution ----------------------------
        select case(isrcdir)

         case(1)
c        pencil beam
          u=cos(phibeam/180.d0*pi)*sin(thetabeam/180.d0*pi)
          v=sin(phibeam/180.d0*pi)*sin(thetabeam/180.d0*pi)
          w=cos(thetabeam/180.d0*pi)

         case(2)
c        isotropic 4pi source
          costhe=-1.+2.*alrand(4.d0) 
          sinthe=sqrt(1.d0-costhe*costhe)
          phi= 2.d0*pi*alrand(5.d0)
          u=sinthe*cos(phi)  
          v=sinthe*sin(phi)
          w=costhe

         case(3)
c        Detector dimensions and source-detector distance
          dsource=abs(zwdet-z)    ! source-detector distance (cm)
          tgtheta=rdet/dsource   ! aperture angle
          cosmin=1.d0/sqrt(1.d0+tgtheta*tgtheta)
          cosmax=1.d0  !(theta=0.)

c         isotropic source limited to [cosmin, cosmax]
          costhe=cosmin+(cosmax-cosmin)*alrand(4.d0) ! generate cos(theta)
          sinthe=sqrt(1.d0-costhe*costhe)
          phi= 2.d0*pi*alrand(5.d0)
          u=sinthe*cos(phi)  
          v=sinthe*sin(phi)
          w=costhe

         case(4)
c        Divergent beam 
          cosmax=1.d0  !(theta=0.)
          cosmin=cos(thetabeam/180.d0*pi) ! thetabeam: semi-aperture angle
c         isotropic source limited to [cosmin, cosmax]
          costhe=cosmin+(cosmax-cosmin)*alrand(4.d0) ! generate cos(theta)
          sinthe=sqrt(1.d0-costhe*costhe)
          phi= 2.d0*pi*alrand(5.d0)
          u=sinthe*cos(phi)  
          v=sinthe*sin(phi)
          w=costhe

        case default
c        pencil beam
          u=0.
          v=0.
          w=1.d0

       end select


c       select energy distribution ----------------------------
c       energy in MeV 

        select case(isrce)

         case(1)
         ! monoenergectic 
          e=Epeak

         case(2)
         ! some peaks
           call ulrndls(psrc,npeak,iout)
           e=Esrc(iout)

         case(3)
         ! spectrum
           call ulspc(idsrc,eout)
           e=eout

         case default
           e=Epeak
         end select

        end

