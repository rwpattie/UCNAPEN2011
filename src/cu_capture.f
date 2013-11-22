c-----------------------------------------------------------------
      subroutine cu_capture
c-----------------------------------------------------------------
c
      implicit DOUBLE PRECISION(A-H,J-M,O-Z), integer*4(i,n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
      COMMON/CUCAP/CLINES(4,424),CLINES2(4,325)
c------------------------------------------------------------------
c     open and read into an array the gamma lins from neutron capture 
c     on cu65.
c     the line energies and intensities are taken from CapGam
c        http://www.nndc.bnl.gov/capgam/byn/page070.html 

      CUISOTOPE = RAND(1.D0)
      PICK = 1.D2*RAND(1.D0)

      IF(CUISOTOPE.LT.67.)THEN
100     CONTINUE
        NROW  = NINT(325*RAND(1.D0))
        IF(CLINES2(3,NROW) .LE. PICK) GOTO 100
        E = CLINES2(1,NROW)*1.D3 ! TABLE ENERGIES ARE IN keV
      ELSE
110     CONTINUE
        NROW  = NINT(424*RAND(1.D0))
        IF(CLINES(3,NROW) .LE. PICK) GOTO 110
        E = CLINES(1,NROW)*1.D3 ! TABLE ENERGIES ARE IN ke
      ENDIF

      KPAR  = 2 ! Set Particle Type as Photon
      theta = pi*rand(1.0D0) ! isotropic distribution
      phi   = 2.0*pi*rand(1.0D0)

      w   = dcos(theta)
      u   = dsin(theta)*dcos(phi)
      v   = dsin(theta)*dsin(phi)

      Z   = 1.5d2 - 3.0d2*rand(1.d0)
      RRO = 6.0001d0 ! SHOULD ORIGINATE THE PHOTON ON THE TUBE SURFACE
      PSI = 2.0*pi*rand(1.0d0)
      X   = RRO*DCOS(PSI)
      Y   = RRO*DSIN(PSI)

      RETURN
      END       
