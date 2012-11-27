      subroutine al_decay(npar)
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy1
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------c
      IF(MOD(NPAR,2).EQ.0)THEN
        KPAR = 2
        E = 1778969.0
        BETA = 1.
      ELSEIF(MOD(NPAR,2)-1.EQ.0)THEN
        KPAR = 1
100   CONTINUE
        E = energy1()
        E = E*ME
        IF(E.LT.1000)GOTO 100
      ENDIF
     
      RAD = 5.0+RAND(2.0D0)
      PHI = 2*PI*RAND(1.D0)
      X = RAD*DCOS(PHI)
      Y = RAD*DSIN(PHI)
      Z = 150.D0
      PHI2 = 2*PI*RAND(1.D0)
      THETA = PI*RAND(1.D0)
      U = DSIN(THETA)*DCOS(PHI2)
      V = DSIN(THETA)*DSIN(PHI2)
      W = DCOS(THETA)

      RETURN
      END
C----------------------------------------------------------------------C
      double precision function energy1()
      implicit double precision(a-h,J-M,o-z), integer*4(i,n)
      parameter(E0=5.603299497e0,pi=3.141592654d0)
      common/rseed/iseed1,iseed2
      external rand
c
50    E=(E0-1.0)*RAND(1.D0)
      y=100*RAND(1.D0)
      G=(-E/(DSQRT(E**2+2*E)*137))
      FERMI=(2.0*PI*G)/(DEXP(2*PI*G)-1)
      f=Fermi*DSQRT(E**2+2*E)*(E0-(E+1))**2*(E+1)
      if (f.lt.y) goto 50
      energy1 = E
c 
      return
      end

       

