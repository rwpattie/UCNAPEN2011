      subroutine bi_decay
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      INCLUDE 'pmcomms.f'
      EXTERNAL rand
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from Bi207
c    taken from :
c
c  http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=207BI&unc=nds
c------------------------------------------------------------------------
      GPROB = 1.029*RAND(1.0d0)

      ! Select the correct gamma branch
      ! there are 3 level 7/2- (7.03%), 13/2+ (84.0%), and 5/2- (8.9%)
      ! This routine selects the energy level of the initial decay
      ! and the routines randomally determine the cascade structure
      ! and whether electrons are ejected by the emitted gamma.
      if(GPROB.LE.0.84)THEN
         CALL BI1063K
         CALL BI569KS
      ELSE IF(GPROB.GT.0.84.AND.GPROB.LE.0.929)THEN
         CALL BI569K
      ELSE IF(GPROB.GT.0.929.AND.GPROB.LE.1.0)THEN
         GPROB = 7.0d0*RAND(1.0D0)
         IF(GPROB.LE.6.87)THEN
              CALL BIBRANCH1770
         ELSE 
              CALL BIBRANCH1442
         ENDIF
      ELSE IF(GPROB.GT.1.0)THEN
         CALL BIAUGER
      endif
  
      return
      end
C-----------------------------------------------------------------------------C
      SUBROUTINE BIBRANCH1442
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      EXTERNAL rand
      
      ILB(1) = 1
      ILB(2) = 0
      ILB(3) = 0
      ILB(4) = 0
      ILB(5) = 0
      
      phi   = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)
      u     = dsin(theta)*dcos(phi)
      v     = dsin(theta)*dsin(phi)
      w     = dcos(theta)
      kpar  = 1
      x     = 0.4*(0.5-1.0*rand(1.d0)) ! for background runs, tin source is
      y     = 0.4*(0.5-1.0*rand(1.d0)) !   approximately at (-5.5,0,155).
      z     = 1.
      
      XE = 1.310757D-1*RAND(1.0d0)
      IF(XE.LE.1.31D-1)THEN
         E    = 1442.2D3
         KPAR = 2
      ELSE IF(XE.GT.1.31D-1.AND.XE.LE.0.1310613)THEN
         E    = 1426.34D3
      ELSE IF(XE.GT.0.1310613)THEN
         E    = 1438.35D3
      ENDIF
      
      PTYPE = 5
      
      BPROB = 0.12869 * RAND(1.D0)
      IF(BPROB.LE.0.128)THEN
         CALL BI328KS
         CALL BI569KS
      ELSE 
         CALL BI897KS
      ENDIF
      
      RETURN
      END
C-----------------------------------------------------------------------------C
      SUBROUTINE BI897KS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      DIMENSION ILBH(5)
      EXTERNAL rand,energy1
      
      ILBH(1) = -1
      ILBH(2) = 0
      ILBH(3) = 0
      ILBH(4) = 0
      ILBH(5) = 0
      
      XS=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      YS=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      ZS = 1.!155.d0
      ! Create Second 
      theta  = 1.0- 2.0*rand(1.d0)
      psi    = 2*pi*rand(1.d0)
      us     = dsin(dacos(theta))*dcos(psi)
      vs     = dsin(dacos(theta))*dsin(psi)
      ws     = theta

      xe = 0.130962*RAND(1.0D0)
      KPAR2 = 1
      if(XE.LE.0.128)THEN
        E2    = 897.77D3
        KPAR2 = 2
      ELSE IF(XE.GT.0.128.AND.XE.LE.0.13046)then
        E2    = 809.77D3
      else if(xe.gt.0.13046.AND.XE.LE.0.130867)then
        E2    = 881.91D3
      else if(xe.gt.0.130867)then
        E2    = 893.93D3
      endif
      
      PTYPE = PTYPE + 30
      CALL STORES(E2,XS,YS,ZS,US,VS,WS,1,KPAR2,ILBH,0)
      
      RETURN       
      END
C-----------------------------------------------------------------------------C      
      SUBROUTINE BI328KS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      DIMENSION ILBH(5)
      EXTERNAL rand,energy1
      
      ILBH(1) = -1
      ILBH(2) = 0
      ILBH(3) = 0
      ILBH(4) = 0
      ILBH(5) = 0
      
      XS=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      YS=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      ZS = 1.
      ! Create Secondary
      theta  = 1.0- 2.0*rand(1.d0)
      psi    = 2*pi*rand(1.d0)
      us     = dsin(dacos(theta))*dcos(psi)
      vs     = dsin(dacos(theta))*dsin(psi)
      ws     = theta

      xe = 9.175E-4*RAND(1.0D0)
      KPAR2 = 1
      if(XE.LE.6.9D-4)THEN
        E2    = 328.1D3
        KPAR2 = 2
      ELSE IF(XE.GT.6.9D-4.AND.XE.LE.8.78D-4)then
        E2    = 240.10D3
      else if(e.gt.8.78D-4.AND.XE.LE.9.10D-4)then
        E2    = 312.24D3 
      else if(xe.gt.9.1D-4)then
        E2    = 324.25D3
      endif
      PTYPE = PTYPE + 20
      CALL STORES(E2,XS,YS,ZS,US,VS,WS,1,KPAR2,ILBH,0)
      RETURN 
      END
C-----------------------------------------------------------------------------C
      SUBROUTINE BIBRANCH1770
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      EXTERNAL rand
      
      ILB(1) = 1
      ILB(2) = 0
      ILB(3) = 0
      ILB(4) = 0
      ILB(5) = 0
      
      phi = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)
      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)
      kpar = 1
      x=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z = 1.
      
      XE = 6.8972D0*RAND(1.0D0)
      IF(XE.LE.6.87)THEN
         E    = 1770.228D3
         KPAR = 2
      ELSE IF(XE.GT.6.87.AND.XE.LE.6.8738)THEN
         E    = 1682.224D3
      ELSE IF(XE.GT.6.8738)THEN
         E    = 1754.367D3
      ENDIF
      
      PTYPE = 5
      
      CALL BI569KS
      
      RETURN
      END
C-----------------------------------------------------------------------------C
      SUBROUTINE BIAUGER
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      EXTERNAL rand
      
      ILB(1) = 1 ! THIS DEFINES THE EVENT AS BEING FROM THE 
      ILB(2) = 0  ! SOURCE PREVENTING LOCAL ENERGY SUBTRACTION
      ILB(3) = 0
      ILB(4) = 0
      ILB(5) = 0     
   
      kpar = 1
      x=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z = 1.!155.d0
        
      theta  = 1.0 - 2.0*rand(1.d0)
      PSI    = 2*PI*rand(1.d0)
      U      = dsin(dacos(theta))*dcos(psi)
      V      = dsin(dacos(theta))*dsin(psi)
      W      = theta
      E      = 56.7D3
      PTYPE  = 10
      RETURN
      END
C-----------------------------------------------------------------------------C      
      SUBROUTINE BI569KS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      DIMENSION ILBH(5)
      EXTERNAL rand
      
      ILBH(1) = -1
      ILBH(2) = 0
      ILBH(3) = 0
      ILBH(4) = 0
      ILBH(5) = 1
      
      XS=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      YS=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      ZS = 1.
      ! Create Second 
      theta  = 1.0- 2.0*rand(1.d0)
      psi    = 2*pi*rand(1.d0)
      us     = dsin(dacos(theta))*dcos(psi)
      vs     = dsin(dacos(theta))*dsin(psi)
      ws     = theta

      xe = 99.84*rand(1.0d0)
      KPAR2 = 1
      if(xe.le.1.537)then
        E2 = 4.816935e5
      Else if(xe.gt.1.537.and.xe.le.1.979)then
        E2 = 5.53872e5
      else if(xe.gt.1.979.and.xe.le.2.09)then
        E2 = 5.658473e5
      else if(xe.gt.2.09)then
        E2 = 5.69698e5
        KPAR2 = 2
      endif
      
      IF(KPAR2.EQ.2)THEN      
         PTYPE = PTYPE + 10
      ELSE IF(KPAR2.EQ.1)THEN
         PTYPE = PTYPE + 100
      ENDIF
      CALL STORES(E2,XS,YS,ZS,US,VS,WS,1.0D0,KPAR2,ILBH,0)
      
      RETURN 
      END
c-----------------------------------------------------------------------------c
      SUBROUTINE BI569K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      EXTERNAL rand
      
      ILB(1) = 1
      ILB(2) = 0
      ILB(3) = 0
      ILB(4) = 0
      ILB(5) = 0
      
      phi = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)
      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)
      kpar = 1
      x=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z = 1.
      
      xe = 99.84*rand(1.d0)
      if(xe.le.1.537)then
        E = 4.816935e5
      Else if(xe.gt.1.537.and.xe.le.1.979)then
        E = 5.53872e5
      else if(xe.gt.1.979.and.xe.le.2.09)then
        E = 5.658473e5
      else if(xe.gt.2.09)then
        E = 5.69698e5
        KPAR = 2
      endif
     
      PTYPE = 1
      
      RETURN 
      END
c-----------------------------------------------------------------------------c
      SUBROUTINE BI1063K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      EXTERNAL rand
      
      ILB(1) = 1
      ILB(2) = 0
      ILB(3) = 0
      ILB(4) = 0
      ILB(5) = 0
      
      phi   = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)
      u     = dsin(theta)*dcos(phi)
      v     = dsin(theta)*dsin(phi)
      w     = dcos(theta)
      kpar  = 1
      x     = 0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y     = 0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z     = 1.
      
      xe = 83.86*RAND(1.D0)
	    
      if(xe.le.7.08)then
          E = 975.651D3
      else if(xe.gt.7.08.and.xe.le.8.92)then
          E = 1047.795D3
      else if(xe.gt.8.92.and.xe.le.9.36)then
          E = 1.059805D6
      else if(xe.gt.9.36)then
          E = 1.063656e6
          kpar = 2
      endif
      
      IF(kpar.eq.1)then
        PTYPE = 3
      ELSEIF(KPAR.EQ.2)THEN
        PTYPE = 5
      ENDIF
      
      RETURN 
      END
C-----------------------------------------------------------------------------C
      subroutine ce_decay
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy1
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from Ce139
c    taken from :
c
c  http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=139CE&unc=nds
c------------------------------------------------------------------------

      eweight = 118.4231
      gweight = 169.43

      x = (eweight+gweight)*rand(1.d0)
        if(x.le.90.2)then
          kpar =1
          E = 3.8d3
        else if(x.gt.90.2.and.x.le.98.5)then
          kpar =1
          E = 27.4d3
        else if(x.gt.98.5.and.x.le.117.65)then
          kpar =1
          E = 126.9329d3
        else if(x.gt.117.65.and.x.le.117.948)then
          kpar =1
          E = 159.5912d3
        else if(x.gt.117.948.and.x.le.118.4231)then
          kpar =1
          E = 164.4962d3
        else if(x.gt.eweight.and.x.le.11.9+eweight)then
          kpar =2
          E = 4.65d3
        else if(x.gt.eweight+11.9.and.x.le.34.4+eweight)then
          kpar =2
          E = 33.034d3
        else if(x.gt.eweight+34.4.and.x.le.75.4+eweight)then
          kpar =2
          E = 33.442d3
        else if(x.gt.eweight+75.4.and.x.le.79.35+eweight)then
          kpar =2
          E = 37.72d3
        else if(x.gt.eweight+79.35.and.x.le.86.97+eweight)then
          kpar =2
          E = 37.801d3
        else if(x.gt.eweight+86.97.and.x.le.89.43+eweight)then
          kpar =2
          E = 38.726d3
        else if(x.gt.eweight+89.43.and.x.le.169.43+eweight)then
          kpar =2
          E = 165.8575d3
        endif

      phi = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)

      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)
1900  continue
      y = (0.15 - 0.3*rand(1.d0))
      x = (0.15 - 0.3*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 0.15) goto 1900
      z = 1.0

      return
      end
c-----------------------------------------------------------------------------C
      subroutine sr_decay
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy1
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from Sr85
c    taken from :
c
c  http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=85Sr&unc=nds
c------------------------------------------------------------------------

      x = 1.d0*rand(1.d0)
      if(x.le.0.89866)then
        E = 498.8070e3
      else if(x.gt.0.89866.and.x.le.1.0)then
        E = 511.9416e3
      endif

      phi = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)

      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)
      kpar = 1
      x=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z = 1.

      return
      end
c-----------------------------------------------------------------------------C
      subroutine xe_135_decay
c----------------------------------------------------------------------
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
c      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy3,energy5
      include 'pmcomms.f'
      DIMENSION ILBH(5)

      betaprob = rand(1.d0)
      if(betaprob.gt.0.25.and.betaprob.le.0.50)goto 1100
      if(betaprob.gt.0.50.and.betaprob.le.0.75)goto 1300
      if(betaprob.gt.0.75)goto 1500
      
      betaprob = rand(1.d0)
      nbetatype = 0
      qend = 0.d0

      ILBH(1) = -1
      ILBH(2) = 1
      ILBH(3) = 1
      ILBH(4) = 0
      ILBH(5) = 1

 
c     select the beta decay end point
c      betaprob = 0.50     
 
      if(betaprob.le.0.96)then
         qend = 915.0d3
         rej = 2.45d1
         rn = 4.7932d0
         nbetatype = 0
      else if(betaprob.gt.0.96.and.betaprob.le.0.9911)then
         qend = 557.0d3
         rej = 8.20d0
         rn = 4.7932d0
         nbetatype = 1
      else if(betaprob.gt.0.9911.and.betaprob.le.0.997)then
         qend = 757.0d3
         rej = 1.58d1
         rn = 4.7932d0
         nbetatype = 2
      else if(betaprob.gt.0.997.and.betaprob.le.0.99823)then
         qend = 103.0d3
         rej = 0.29d0
         rn = 4.7932d0
         nbetatype = 3
      else if(betaprob.gt.0.99823)then
         qend = 184.0e3
         rej = 0.90d0
         rn = 4.7932d0
         nbetatype = 4
      endif

      !rad   = 6.5*rand(1.d0)
      !psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
1000  continue 
      x     = 6.5*(1.0 - 2.0*rand(1.d0))
      y     = 6.5*(1.0 - 2.0*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 6.5) goto 1000
      !x     = rad*dcos(psi)
      !y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      
      Kpar  = 1
      e     = energy3(qend,rej,rn)
   
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      us     = dsin(dacos(theta))*dcos(psi)
      vs     = dsin(dacos(theta))*dsin(psi)
      ws     = theta
C     
C      AUGER ENERGIES
C 
      EE1 = 2.49794E5
      EE2 = 2.13809E5
      EE3 = 2.4408E5
      EE4 = 2.48577E5
      EE5 = 2.49563E5
      EE6 = 2.49778E5
C
C   AUGERS 2
C     
      EE7 = 6.08185E5
      EE8 = 5.72200E5
      EE9 = 6.02471E5
      E10 = 6.06968E5

      if(nbetatype.eq.0)then
        augerprob = 0.966388*rand(1.0d0)
        if(augerprob.le.0.90)then
          CALL STORES(EE1,X,Y,Z,US,VS,WS,1.d0,2,ILBH,0)  
        else if(augerprob.gt.0.90.and.augerprob.le.0.9561)then
          CALL STORES(EE2,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0)
        else if(augerprob.gt.0.9561.and.augerprob.le.0.9643)then
          CALL STORES(EE3,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0)
        else if(augerprob.gt.0.9643.and.augerprob.le.0.96599)then
          CALL STORES(EE4,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0)
        else if(augerprob.gt.0.96599.and.augerprob.le.0.9634)then
          CALL STORES(EE5,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0) 
        else if(augerprob.gt.0.9634)then
          CALL STORES(EE6,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0)
        endif
      else if(nbetatype.eq.1)then
        augerprob = 2.921005*rand(1.d0)
        if(augerprob.le.2.90)then
         CALL STORES(EE7,X,Y,Z,US,VS,WS,1.d0,2,ILBH,0)
         else if(augerprob.gt.2.90.and.augerprob.le.2.9182)then
         CALL STORES(EE8,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0)
        else if(augerprob.gt.2.9182.and.augerprob.le.2.92053)then
         CALL STORES(EE9,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0)
        else if(augerprob.gt.2.92053)then
         CALL STORES(E10,X,Y,Z,US,VS,WS,1.d0,1,ILBH,0)
        endif
      endif 
      ptype = 1
      return

1100  continue
      ptype = 2
      betaprob = 97.137*rand(1.d0)
 
      if(betaprob.lt.61.6)then
        E = 129.369d3
      else if(betaprob.ge.61.6.and.betaprob.lt.90.4)then
        E = 158.477d3
      else if(betaprob.ge.90.4.and.betaprob.lt.96.99)then
        E = 162.788d3
      else if(betaprob.ge.96.99.and.betaprob.lt.97.137)then
        E = 163.914d3
      endif
   
      !rad   = 6.5*rand(1.d0)
      !psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
1200  continue
      x     = 6.5*(1.0 - 2.0*rand(1.d0))
      y     = 6.5*(1.0 - 2.0*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 6.5) goto 1200
      !x     = rad*dcos(psi)
      !y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      Kpar  = 1
    
      return
  
1300  continue
      ptype = 3
      betaprob = rand(1.d0)
 
      if(betaprob.le.0.985)then
        qend = 346.4e3
        rej = 3.2d0
        rn = 4.7831d0
        e = energy5(qend,rej,rn)
      else if(betaprob.gt.0.985)then
        qend = 266.8e3
        rej = 1.9d0
        rn = 4.7831d0
        e = energy5(qend,rej,rn)
      !else if(betaprob.gt.0.999)then ! A lot of spectrum near rejection energy
        !qend = 43.6e3
        !rej = 0.058d0
        !rn = 4.7831d0
        !e = energy5(qend,rej,rn)
      endif 
      
      !rad   = 6.5*rand(1.d0)
      !psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
1400  continue
      x     = 6.5*(1.0 - 2.0*rand(1.d0))
      y     = 6.5*(1.0 - 2.0*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 6.5) goto 1400
      !x     = rad*dcos(psi)
      !y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      Kpar  = 1

      return
 
1500  continue
      ptype = 4
      betaprob = rand(1.d0)
 
      betaprob = 90.955*rand(1.d0)

      if(betaprob.lt.63.5)then
        E = 198.660d3
      else if(betaprob.ge.63.5.and.betaprob.lt.84.18)then
        E = 227.768e3
      else if(betaprob.ge.84.18.and.betaprob.lt.89.73)then
        E = 232.079e3
      else if(betaprob.ge.89.73)then
        E = 233.013e3
      endif

      !rad   = 6.5*rand(1.d0)
      !psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
1600  continue
      x     = 6.5*(1.0 - 2.0*rand(1.d0))
      y     = 6.5*(1.0 - 2.0*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 6.5) goto 1600
      !x     = rad*dcos(psi)
      !y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      Kpar  = 1

           
      return
      end      

c----------------------------------------------------------------------------C
      double precision function energy3(qend,rej,rn)
      implicit double precision(A-H,J-M,O-Z), integer*4(i,n)
      parameter(emass=510.998928d3)
      parameter(pi   =3.141592654d0)
      parameter(alphainv=1.37036d2)
      parameter(zed  =5.40d1)
      parameter(lambda=3.86159268d2)
      parameter(B    =1.0170d0)
      parameter(ga=1.73169d0)
      common/rseed/iseed1,iseed2
      external rand
c----------------------------------------------------------------------------c
c  Fermi factor from D.H. Wilkinson, Nucl. Phys. A 377, 474 (1982).
c  Nucl. radii from I. Angeli, K. P. Marinova, Atomic Data and Nuclear 
c  Tables 99, 69 (2013).  Linear extrapolation for Xe135 & appx. Xe133
c  F0L0 and L0 tables from Behrens and Janecke, Numerical Tables for 
c                          Beta Decay and Electron Capture (1969).
c----------------------------------------------------------------------------c

      EO = qend/emass + 1
50    E=(EO-1.0)*RAND(1.D0)
      if(E.lt.0.03) goto 50
      Y=rej*RAND(1.D0)
      R=rn/lambda
      C=zed/alphainv
      A=B-1.0052d0 
      S=A/EO
      P=DSQRT((E+1)**2-1)
      CALL ZGMFN(DSQRT(1.0d0-C**2),(C*(E+1))/P,1,Q,U)
      FERMI=2*(1.0d0+DSQRT(1.0d0-C**2))*(1/((2*P*R)**(2*
     1     (1-DSQRT(1.0d0-C**2)))))*DEXP((PI*C*(E+1))/P)*
     1     (DSQRT(Q**2+U**2)/ga)**2*(B-S*E) ! D.H. Wilkinson (also Behrens/Janecke)
      if (Q.lt.0.0) goto 50
      W=FERMI*P*(EO-(E+1))**2*(E+1)
      if (W.lt.Y) goto 50
      !print*,'FERMI = ', FERMI
      E = E*emass
      energy3 = E
c 
      return
      end

c----------------------------------------------------------------------------C
      double precision function energy5(qend,rej,rn)
      implicit double precision(A-H,J-M,O-Z), integer*4(i,n)
      parameter(emass=510.998928d3)
      parameter(pi   =3.141592654d0)
      parameter(alphainv=1.37036d2)
      parameter(zed  =5.40d1)
      parameter(lambda=3.86159268d2)
      parameter(B    =1.0170d0)
      parameter(ga=1.73169d0)
      common/rseed/iseed1,iseed2
      external rand
c----------------------------------------------------------------------------c
c  Fermi factor from D.H. Wilkinson, Nucl. Phys. A 377, 474 (1982).
c  Nucl. radii from I. Angeli, K. P. Marinova, Atomic Data and Nuclear 
c  Tables 99, 69 (2013).  Linear extrapolation for Xe135 & appx. Xe133
c  F0L0 and L0 tables from Behrens and Janecke, Numerical Tables for 
c                          Beta Decay and Electron Capture (1969).
c----------------------------------------------------------------------------c

      EO = qend/emass + 1
60    E=(EO-1.0)*RAND(1.D0)
      if(E.lt.0.03) goto 60
      Y=rej*RAND(1.D0)
      R=rn/lambda
      C=zed/alphainv
      A=B-1.0052d0
      S=A/EO
      P=DSQRT((E+1)**2-1)
      CALL ZGMFN(DSQRT(1.0d0-C**2),(C*(E+1))/P,1,Q,U)
      FERMI=2*(1.0d0+DSQRT(1.0d0-C**2))*(1/((2*P*R)**(2*
     1     (1-DSQRT(1.0d0-C**2)))))*DEXP((PI*C*(E+1))/P)*
     1     (DSQRT(Q**2+U**2)/ga)**2*(B-S*E) ! D.H. Wilkinson (also Behrens/Janecke)
      if (Q.lt.0.0) goto 60
      W=FERMI*P*(EO-(E+1))**2*(E+1)
      if (W.lt.Y) goto 60
      !print*,'FERMI2 = ', FERMI
      E = E*emass
      energy5 = E
c 
      return
      end
