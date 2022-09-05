C     Chamber + target + window+Detector

      subroutine ulgeo

        use geometry_mod
        use source_mod

      implicit double precision(a-h,o-z), integer*4(i-n)
      dimension par(6), xcenter(3)

        idcha=2
        idtar=3
        idwin=4
        iddet=5

c Universe - 1st volume !Mandatory!                       
        mat=0   
        id=1
        par(1)=rcha*2.d0
        par(2)=rcha*2.d0 
        par(3)=hcha*2.
        itype=100     ! box
        imother=0     ! no mother volume to universe !
        call ulvolume(id,par,itype,imother,mat) ! define volume

        mat=1   ! chamber 
        id=idcha
        par(1)=rcha
        par(2)=hcha
        itype=200     ! cylinder
        imother=1   
        call ulvolume(id,par,itype,imother,mat)
        xcenter(1)=xcha
        xcenter(2)=ycha
        xcenter(3)=zcha 
        CALL ulposi(id,xcenter)

        mat=2   ! target
        id=idtar
        par(1)=rtar
        par(2)=htar
        itype=200     ! cylinder
        imother=idcha ! target inside chamber   
        call ulvolume(id,par,itype,imother,mat)
        xcenter(1)=xtar
        xcenter(2)=ytar
        xcenter(3)=ztar
        CALL ulposi(id,xcenter)

        mat=3   ! detector window
        id=idwin
        par(1)=rwin
        par(2)=hwin
        itype=200     ! cylinder
        imother=idcha ! window inside chamber   
        call ulvolume(id,par,itype,imother,mat)
        xcenter(1)=xwin
        xcenter(2)=ywin
        xcenter(3)=zwin
        call ulposi(id,xcenter)

        mat=4   ! detector
        id=iddet
        par(1)=rdet
        par(2)=hdet
        itype=200     ! cylinder
        imother=idcha ! detector inside chamber   
        call ulvolume(id,par,itype,imother,mat)
        xcenter(1)=xdet
        xcenter(2)=ydet
        xcenter(3)=zdet
        call ulposi(id,xcenter)

        end
