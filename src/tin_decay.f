      subroutine tin_decay
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      INCLUDE 'pmcomms.f'
      EXTERNAL rand
c---------------------------------------------------------------------c
c     Generating radiation mode; 
c     (1) Electron capture will either emit a k or l shell beta or a
c         391.698keV gamma. Half life of 115.09 days.  Assume the source 
c         used was a 1uCurie source of this gamma. Ref: NNDC
c     
c     (2) Internal conversion x-rays in the tin create auger electrons 
c         with energies ranging from ~3keV to 25keV. Ref: Table of 
c         Isotopes by R. B. Firestone (c1996). 
c---------------------------------------------------------------------c
c     Generating X-rays and Auger Electrons from K-Shell vacancies    c
c---------------------------------------------------------------------c
      bk  = 29.2001d3   ! Binding Energies for Atomic levels
      bL1 = 4.4647d3
      bL2 = 4.1561d3
      bL3 = 3.9288d3
      bM3 = 0.7144d3
      bM3p= 0.7656d3
c---------------------------------------------------------------------c
      goto 90

      if(mod(npar,2).eq.0)then
        goto 90
      elseif(mod(npar,2)-1.eq.0)then
        goto 80
      endif

80    continue
c---------------------------------------------------------------------c
      xka1=25.271d3  ! X-ray energies
      xka2=25.044d3
      xkb1=28.486d3
      xkb2=29.111d3
      xkb3=28.444d3
      xLa1=3.444d3
      xLa2=3.435d3
      xLb1=3.663d3
      xLb2=3.905d3
      xLb3=3.750d3
c----------------------------------------------------------------------c
      xe = 0.0d0
      E  = 0.0d0
      aprob=99.76d0*rand(1.0d0)
c
      if(aprob.le.1.05)then
         Xe=bL1
      elseif(aprob.gt.1.05.and.aprob.le.2.29)then
         Xe=bL1
      elseif(aprob.gt.2.29.and.aprob.le.3.90)then
         Xe=bL1
      elseif(aprob.gt.3.90.and.aprob.le.4.067)then
         Xe=bL2
      elseif(aprob.gt.4.067.and.aprob.le.7.667)then
         Xe=bL2
      elseif(aprob.gt.7.667.and.aprob.le.9.487)then
         Xe=bL3
      elseif(aprob.gt.9.487.and.aprob.le.55.187)then  
         E=xka1
      elseif(aprob.gt.55.187.and.aprob.le.79.887)then
         E=xka2
      elseif(aprob.gt.79.887.and.aprob.le.87.877)then
         E=xkb1
      elseif(aprob.gt.87.877.and.aprob.le.90.067)then
         E=xkb2
      elseif(aprob.gt.90.067.and.aprob.le.94.217)then
         E=xkb3
      elseif(aprob.gt.94.217.and.aprob.le.97.117)then
         E=xLa1
      elseif(aprob.gt.97.117.and.aprob.le.97.437)then
         E=xLa2
      elseif(aprob.gt.97.437.and.aprob.le.99.187)then
         E=xLb1
      elseif(aprob.gt.99.187.and.aprob.le.99.647)then
         E=xLb2
      elseif(aprob.gt.99.647.and.aprob.le.99.76)then
         E=xLb3
      endif
c
      if(aprob.gt.9.487)then
          kpar=2
      else
          kpar=1
          E=bK-Xe-bM3-0.75*(bM3p-bM3)
      endif
      goto 100
c---------------------------------------------------------------------c
c     Generating Auger Electrons from the 391.698keV gamma from electron
c     capture.
c---------------------------------------------------------------------c
90    continue
      econv= rand(1.d0)
      gprob= rand(1.d0)

      if(gprob.le.0.968545)then
            egamma=391.698d3
         else
            E   = 255.134d3
            kpar= 2
            goto 100
      endif

      if(econv.le.0.65239)then
         E=Egamma
         kpar=2
      else
         eprob=rand(1.0d0)
         if(eprob.le.0.80577)then
            Ebeta=363758.!egamma-bk
         else if(eprob.gt.0.80577.and.eprob.le.0.96245)then
            Ebeta=387461.!egamma-bL1
         else if(eprob.gt.0.96245.and.eprob.le.0.99426)then
            Ebeta=390872.!egamma-bL1
         else if(eprob.gt.0.99426)then
            Ebeta=391576.!
         endif
         E=Ebeta
         kpar=1
      endif
      
c---------------------------------------------------------------------c
c     Generating the angular distribution which will assumed to be    c
c     isotropic.                                                      c
c---------------------------------------------------------------------c
100   continue
      
1000  continue
      y = (0.15 - 0.3*rand(1.d0))  ! for background runs, tin source is
      x = (0.15 - 0.3*rand(1.d0))  ! approximately at (-5.5,0,155).
      !x = (4.68 - 0.3*rand(1.d0)) ! edge of decay trap test
      if(sqrt(x**2 + y**2) .gt. 0.15) goto 1000  
      !if(sqrt((x-4.53)**2 + y**2) .gt. 0.15) goto 1000
      z=1.0!1.55d2
      !z=1.87309d2!1.55d2 
      costh = 1.-2.*rand(1.d0)   ! sample from cos(theta)
      theta = dabs(dacos(costh)) !    pi*rand(1.0d0)
      phi   =  2.*pi*rand(1.0d0)
      w=dcos(theta)
      u=dsin(theta)*dcos(phi)
      v=dsin(theta)*dsin(phi)
      !w=1.
      !u=0.
      !v=0.
      win=w
      uin=u
      vin=v
      ein=e
c----------------------------------------------------------------------c
c       Start event time
c----------------------------------------------------------------------c
      time=3.37837*rand(1.d0)
c----------------------------------------------------------------------c

      return
      end
