C-------------------------------------------------------------------------C
      DOUBLE PRECISION FUNCTION POSIT(XX)
C-------------------------------------------------------------------------C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      DIMENSION XX(170)
      PARAMETER( MAXRAD = 8.683)
C------------------------------------------------------------------------C
      CENTERSPACE = (2.*MAXRAD)/(170.0)

      WS = 0
      ST = 0

      DO I=2,169
          IF(XX(I).GT.0)THEN
             WS = WS + ((I-1)*CENTERSPACE - MAXRAD)*XX(I)
             ST = ST + XX(I)
          ENDIF
      ENDDO

      IF(ST.GT.0)THEN
         POSIT = WS/ST
      ELSE
         POSIT = 0.
      ENDIF

      RETURN
      END
c----------------------------------------------------------------------c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Detrep attempts to simulate the imperfections in the energy 
c     resolutions of our plastic scintillator.  Given a gaussian
c     response function whose width is sigma=sqrt(2.5*E), where E is the
c     true energy, the recorded energy is determined.  Currently not
c     working with the best efficiency.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function detrep(E0)
      implicit double precision(a-h,o-z), integer*4(i-n)
      PARAMETER  (PI=3.1415926535897932D0)
      external rand
c----------------------------------------------------------------------c
      E0=E0/1000.0
      SIG = DSQRT(2.5*E0)

500   CONTINUE

      E=2*E0*rand(1.d0)

      Y=RAND(1.0D0)

      F=DEXP(-(E-E0)**2/(2*SIG*SIG))

      IF(F.LT.Y) then
        GOTO 500
      endif

      DETREP=E*1000.0
c----------------------------------------------------------------------c      
      RETURN
      END 
************************************************************************