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
      bk  = 27.9399d3   ! Binding Energies for Atomic levels
      bL1 = 4.2375d3
      bL2 = 3.9380d3
      bL3 = 3.7301d3
      bM3 = 0.6643d3
      bM3p= 0.7144d3
      xka1=24.210d3  ! X-ray energies
      xka2=24.002d3
      xkb1=27.276d3
      xkb2=27.863d3
      xkb3=27.238d3
      xLa1=3.287d3
      xLa2=3.279d3
      xLb1=3.487d3
      xLb2=3.714d3
      xLb3=3.573d3
      xLb4=3.535d3
c---------------------------------------------------------------------c
      goto 90  
c---------------------------------------------------------------------c
c     Generating Auger Electrons from the 391.698keV gamma from electron
c     capture.
c---------------------------------------------------------------------c
90    continue
      econv= 1.34761d0*rand(1.d0)
      gprob= rand(1.d0)

      if(gprob.le.0.968545)then
            Egamma=391.698d3
         else
            E   = 255.134d3
            kpar= 2
            goto 100
      endif

      if(econv.le.0.65239)then
         E=Egamma
         kpar=2
      else
         eprob=2.0d0*rand(1.0d0)
	 atot=99.185d0
         if(eprob.le.(1.11/atot))then
            Xe=bL1
         elseif(eprob.gt.(1.11/atot).and.eprob.le.(2.41/atot))then
            Xe=bL1
         elseif(eprob.gt.(2.41/atot).and.eprob.le.(4.14/atot))then
            Xe=bL1
         elseif(eprob.gt.(4.14/atot).and.eprob.le.(4.319/atot))then
            Xe=bL2
         elseif(eprob.gt.(4.319/atot).and.eprob.le.(8.119/atot))then
            Xe=bL2
         elseif(eprob.gt.(8.119/atot).and.eprob.le.(10.089/atot))then
            Xe=bL3
         elseif(eprob.gt.(10.089/atot).and.eprob.le.(55.389/atot))then
            Ex=xka1
         elseif(eprob.gt.(55.389/atot).and.eprob.le.(79.889/atot))then
            Ex=xka2
         elseif(eprob.gt.(79.889/atot).and.eprob.le.(87.739/atot))then
            Ex=xkb1
         elseif(eprob.gt.(87.739/atot).and.eprob.le.(89.829/atot))then
            Ex=xkb2
         elseif(eprob.gt.(89.829/atot).and.eprob.le.(93.899/atot))then
            Ex=xkb3
         elseif(eprob.gt.(93.899/atot).and.eprob.le.(96.699/atot))then
            Ex=xLa1
         elseif(eprob.gt.(96.699/atot).and.eprob.le.(97.009/atot))then
            Ex=xLa2
         elseif(eprob.gt.(97.009/atot).and.eprob.le.(98.639/atot))then
            Ex=xLb1
         elseif(eprob.gt.(98.639/atot).and.eprob.le.(99.079/atot))then
            Ex=xLb2
         elseif(eprob.gt.(99.079/atot).and.eprob.le.(99.144/atot))then
            Ex=xLb3
         elseif(eprob.gt.(99.144/atot).and.eprob.le.(99.185/atot))then
            Ex=xLb4
	 elseif(eprob.gt.1.0.and.eprob.le.1.80577)then
            Ebeta=363758.!egamma-bk
         elseif(eprob.gt.1.80577.and.eprob.le.1.96245)then
            Ebeta=387461.!egamma-bL1
         elseif(eprob.gt.1.96245.and.eprob.le.1.99426)then
            Ebeta=390872.!egamma-bL1
         elseif(eprob.gt.1.99426.and.eprob.le.2.0)then
            Ebeta=391576.!
         endif

         if(eprob.le.(10.089/atot))then
            kpar=1
            E=bK-Xe-bM3-0.75*(bM3p-bM3)
         elseif(eprob.gt.(10.089/atot).and.eprob.lt.1.0)then   
	    kpar=2
            E=Ex
      	 elseif(eprob.gt.1.0.and.eprob.lt.2.0)then
            kpar=1
            E=Ebeta
      	 endif
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
