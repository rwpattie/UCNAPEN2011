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
c----------------------------------------------------------------------c
      DOUBLE PRECISION FUNCTION BIRKS_LAW(E,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      PARAMETER (Sk = 0.42050, Bk = 0.01907)
      PARAMETER (AA = 116.7,BB= -0.7287,RHO = 1.032);
c      PARAMETER (Sk = 0.6685, Bk = 0.01907)

c     THIS FUNCTION IS AN APPROXIMATION OF BIRKS LAW USING PARAMETERS
C     MEASURED IN JUNHUA'S THESIS FOR THE SCINTILLATOR USED FOR THE
C     SUZUNO DETECTOR 
C       const G4double kb = 0.01907*cm/MeV;                     // Birk's law quenching constant
C       const G4double a = 116.7*MeV*cm*cm/g;           // dEdx fit parameter a*e^(b*E)
C       const G4double b = -0.7287;                                     // dEdx fit parameter a*e^(b*E)
C       const G4double rho = 1.032*g/cm3;                       // scintillator density
C       const G4double dEdx = a*rho*pow(E/keV,b);       // estimated dE/dx
C       return 1.0/(1+kb*dEdx);
c     The factor of e-6 takes the energy loss in eV to MeV
c     S is 1000 times greater since it takes npe/meV -> npe/keV 
      DE = D*1d-6
      dEdx = AA*RHO*(E*1d-3)**BB 
      BIRKS_LAW = DE / (1. + BK*dEdX)*1d6
      IF(BIRKS_LAW.EQ.0)Print*,DE,birks_law,DEDX,E,(E*1D-3),
     1 (E*1D-3)**BB
      
C      BIRKS_LAW = (Sk*(DE/XE)) / (1. + Bk*(DE/XE))
C      BIRKS_LAW = BIRKS_LAW*XE*1.0d6

      RETURN
      END
C-----------------------------------------------------------------------C
      DOUBLE PRECISION FUNCTION DELTAT(DS,E)
c-----------------------------------------------------------------------c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      PARAMETER(C=2.99792548D10,EMASS=510998.0D0)
      
        BETA = DSQRT(E*(E + 2*EMASS))/(E + EMASS)
        DELTAT = DS / (BETA*C)
        
      RETURN
      END
C-----------------------------------------------------------------------C
      SUBROUTINE RECORD_ENERGYLOSS(DTYPE,DE,EPRE)
C-----------------------------------------------------------------------C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'ucnapenmain.h'
      PARAMETER(NB = 5000)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/CNT1/TDEBO(NB),TDEBO2(NB),DEBO(NB)
      EXTERNAL BIRKS_LAW
C-----------------------------------------------------------------------C
       dtype = 1
c      IF(DTYPE.NE.2.AND.DTYPE.LE.8)THEN
c       PRINT*,X,Y,Z,IBODY,MAT,DTYPE,DS
         IF(MAT.EQ.1)CALL WIRECHECK(DE)
C        IF THE DECAY TRAP FOILS ENERGY COUNTERS WITH ENERGY LOST IN
C        BOTH THE FOIL AND THE COATING
         IF((MAT.EQ.8.OR.MAT.EQ.7.OR.MAT.EQ.3)
     1       .AND.Z.GT.140.AND.Z.LT.160)THEN
               EFLE = EFLE+DE
         ELSEIF((mat.eq.8.or.mat.eq.7.or.mat.eq.3)
     1          .and.z.lt.-140.and.z.gt.-160)then
               EFLW = EFLW+DE
         ENDIF
C        IF THE ELECTRON IS IN THE DEAD LAYER ADD ITS ENERGY TO THE EDEAD* 
C        COUNTERS
         IF(IBODY.EQ.416.)THEN
                                ! TURN THE ENERGY LOSS IN A PATH LENGTH INTO
                                ! NPE USING BIRKS LAW,
                                ! CURRENTLY TESTING WHETHER INTEGER OR CONTINUOUS 
                                ! MODEL
               PHTE  = PHTE  + BIRKS_LAW(EPRE,DE)
         ELSE IF(IBODY.EQ.398.)THEN
               PHTW  = PHTW  + BIRKS_LAW(EPRE,DE)
         END IF
c      ENDIF
c
c        check to see if the detectors are triggerred.
c
c      if(dtype .ne. 2) then
          CALL TRIGGERCHECK(PHTW,PHTE,DEBO(440),DEBO(441),TIME,TRGE,
     1 TRGW)
c      else if(dtype .eq. 2) then
c           if(DEBO(4).gt.0.and.trge(1).eq.0.0) trge(1)=time*1.0d9
c           if(DEBO(19).gt.0.and.trgw(1).eq.0.0)trgw(1)=time*1.0d9
c      endif
c
      RETURN 
      END
C-----------------------------------------------------------------------------C
c----------------------------------------------------------------------c
      SUBROUTINE WIRECHECK(DE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(NENTRIES=96)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/MWPC/ECX(170),ECY(170),WCX(170),WCY(170)
      COMMON/EDEP/EPE,EPW,EPER,EPED,EPWD,EPWR,EGEA,EGWA,
     1 EGEDF,EGEDB,EGWDF,EGWDB,EFLE,EFLW,ETUBE,EAL,EME1,EME2,
     1 EMW1,EMW2,PHW,PHEN,PHE,PHWN,EHOLDER,EEAN,EEC1,EEC2,EWAN,
     1 EWC1,EWC2,ECOL
      REAL DECS(NENTRIES)
      REAL COSTHETA(16)
      COMMON/HBOOKU/DECS,COSTHETA
      PARAMETER ( MAXRAD = 8.386)! , CENTERSPACE = 0.1016)
      PARAMETER ( ZCATH1 = 221.5 , ZCATH2 = 223.5 ) 
      PARAMETER ( ZCATHC = 222.5)
c     
c     routine for tallying the position of wirechamber 
c     energy loss.
c
c----------------------------------------------------------------------
      CENTERSPACE = (2.*MAXRAD)/(170.)
      IF(Z.GT.0)THEN
           IF(MAT.EQ.1.AND.Z.LT.ZCATHC)THEN
                   INTX = INT((X+MAXRAD)/CENTERSPACE)+1
                   IF(INTX.LT.0)THEN 
                      INTX = 1
                   ELSE IF(INTX.GT.169)THEN 
                      INTX = 170
                   ENDIF
                   ECX(INTX) = ECX(INTX)+DE
           ELSE IF(MAT.EQ.1.AND.Z.GT.ZCATHC)THEN
                   INTY = INT((Y+MAXRAD)/CENTERSPACE)+1
                   IF(INTY.LT.0)THEN
                     INTY = 1
                   ELSE IF(INTY.GT.169)THEN 
                     INTY = 170
                   ENDIF
                   ECY(INTY) = ECY(INTY)+DE
           ENDIF

           IF((IBODY.EQ.423.OR.IBODY.EQ.424).AND.Z.LT.ZCATH1)THEN
              DECS(16) = DECS(16)+REAL(DE)
              EGEDF = EGEDF + DE
           ELSE IF(IBODY.EQ.423.AND.Z.GT.ZCATH2)THEN
              DECS(46) = DECS(46)+REAL(DE)
              EGEDB = EGEDB + DE
           ELSE IF(IBODY.EQ.423.AND.DABS(Z-ZCATHC).LT.1)THEN
              DECS(15) = DECS(15)+REAL(DE)
              EGEA  = EGEA + DE
           ENDIF

      ELSEIF(Z.LT.0)THEN

           IF(MAT.EQ.1.AND.DABS(Z).GT.ZCATHC)THEN
                   INTX = INT((X+MAXRAD)/CENTERSPACE)+1
                   IF(INTX.LT.0)THEN
                     INTX  = 1
                   ELSE IF(INTX.GT.169)THEN
                     INTX = 170
                   ENDIF
                   WCX(INTX) = WCX(INTX)+DE
           ELSE IF(MAT.EQ.1.AND.DABS(Z).LT.ZCATHC)THEN
                   INTY = INT((Y+MAXRAD)/CENTERSPACE)+1
                   IF(INTY.LT.0)THEN
                     INTY = 1
                   ELSE IF(INTY.GT.169)THEN
                     INTY = 170
                   ENDIF
                   WCY(INTY) = WCY(INTY)+DE
           ENDIF

           IF((IBODY.EQ.406.OR.IBODY.EQ.405).AND.DABS(Z).LT.ZCATH1)THEN
              DECS(20)  =DECS(20)+REAL(DE)
              EGWA = EGWA + DE
           ELSE IF(IBODY.EQ.405.AND.DABS(Z).GT.ZCATH2)THEN 
              DECS(47)  =DECS(47)+REAL(DE)
              EGWDB = EGWDB + DE
           ELSE IF(IBODY.EQ.405.AND.DABS(DABS(Z)-ZCATHC).LT.1)THEN
              DECS(19)=DECS(19)+REAL(DE)
              EGWDF = EGWDF + DE 
           ENDIF
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------------C------------------------------------------------------------------------C
      SUBROUTINE TRIGGERCHECK(PEW,PEE,PEA,PEH,TIME,TRGE,TRGW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      DIMENSION TRGE(10),TRGW(10)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/MULTI/EMX,EMY,WMX,WMY,EPOSX,EPOSY,WPOSX,WPOSY,
     1 EXSCI,EYSCI,WXSCI,WYSCI,TAPD    
C------------------------------------------------------------------------C
      DO I = 1,4
C     CHECK IF THE DETECT HAS BEEN TRIGGERED AT THIS THRESHOLD
C     IF NOT TRIGGERED YET, CHECK THE ENERGY IN MWPC AND SCINTILLATOR
C     ONCE TRIGGERED RECORD THE TIME TO WRITE TO A FILE AT THE END OF
C     THE TRACK.
       IF(TRGE(I).EQ.0.0)THEN
         IF(PEE.GE.2.)THEN
           TRGE(I)=TIME
	   EXSCI=X
           EYSCI=Y
         ENDIF
       ENDIF
C
       IF(TRGW(I).EQ.0.0)THEN
         IF(PEW.GE.2.)THEN
           TRGW(I)=TIME
           WXSCI=X
	   WYSCI=Y
         ENDIF
       ENDIF
C
      ENDDO
      IF(TAPD.EQ.0.0)THEN
	IF(PEA.GE.2.)THEN
    	  TAPD=TIME
	ENDIF
      ENDIF
C
      RETURN
      END
C-------------------------------------------------------------------------C
c======================================================================c
      SUBROUTINE SETCOSTHETAS(IMAT,WO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(NENTRIES=96)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      REAL DECS(NENTRIES)
      REAL COSTHETA(16)
      COMMON/HBOOKU/DECS,COSTHETA
C----------------------------------------------------------------------C
      !ensure a boundary was crossed
!      if(mat.eq.imat)return

      IF(W*Z.GT.0)THEN
        IF(ABS(Z).GT.148.00.AND.ABS(Z).LT.150.0.AND.
     1        (MAT.EQ.3.OR.MAT.EQ.8))THEN
           IF(COSTHETA(1).EQ.0)COSTHETA(1) = ABS(WO)
        ELSE IF(ABS(Z).GT.150.00.AND.ABS(Z).LT.152.0.AND.
     1    (MAT.EQ.0.OR.MAT.EQ.7))THEN
          IF(COSTHETA(2).EQ.0)COSTHETA(2) = ABS(WO)
        ELSE IF(ABS(Z).GT.216.00.AND.ABS(Z).LT.220.0.AND.
     1      (MAT.EQ.3))THEN
          IF(COSTHETA(3).EQ.0)COSTHETA(3) = ABS(WO)
        ELSE IF(ABS(Z).GT.219.00.AND.ABS(Z).LT.222.0.AND.
     1       MAT.EQ.1)THEN
          IF(COSTHETA(4).EQ.0)COSTHETA(4) = ABS(WO)
        ELSE IF(ABS(Z).GT.222.00.and.MAT.EQ.3)THEN
          IF(COSTHETA(5).EQ.0)COSTHETA(5) = ABS(WO)
        ELSE IF(ABS(Z).GT.222.00.and.MAT.EQ.6)THEN
          IF(COSTHETA(6).EQ.0)COSTHETA(6) = ABS(WO)
        ELSE IF(ABS(Z).GT.224.00.AND.MAT.EQ.2)THEN
          IF(COSTHETA(7).EQ.0)COSTHETA(7) = ABS(WO)
        ELSE IF(ABS(Z).GT.224.00.AND.
     1    (MAT.EQ.0.OR.MAT.EQ.7))THEN
          IF(COSTHETA(8).EQ.0)COSTHETA(8) = ABS(WO)
        ENDIF

      ELSE IF(W*Z.LT.0)THEN
         IF(ABS(Z).GT.148.00.AND.ABS(Z).LT.150.01.AND.
     1        (IMAT.EQ.3.OR.IMAT.EQ.8))THEN
           IF(COSTHETA(9).EQ.0)COSTHETA(9) = ABS(WO)
        ELSE IF(ABS(Z).GT.150.00.AND.ABS(Z).LT.152.0.AND.
     1    (IMAT.EQ.0.OR.IMAT.EQ.7))THEN
          IF(COSTHETA(10).EQ.0)COSTHETA(10) = ABS(WO)
        ELSE IF(ABS(Z).GT.216.00.AND.ABS(Z).LT.220.0.AND.
     1      (IMAT.EQ.3))THEN
          IF(COSTHETA(11).EQ.0)COSTHETA(11) = ABS(WO)
        ELSE IF(ABS(Z).GT.219.00.AND.ABS(Z).LT.222.0.AND.
     1       IMAT.EQ.1)THEN
          IF(COSTHETA(12).EQ.0)COSTHETA(12) = ABS(WO)
        ELSE IF(ABS(Z).GT.222.00.AND.IMAT.EQ.3)THEN
          IF(COSTHETA(13).EQ.0)COSTHETA(13) = ABS(WO)
        ELSE IF(ABS(Z).GT.222.00.and.MAT.EQ.6)THEN
          IF(COSTHETA(14).EQ.0)COSTHETA(14) = ABS(WO)
        ELSE IF(ABS(Z).GT.224.00.AND.IMAT.EQ.2)THEN
          IF(COSTHETA(15).EQ.0)COSTHETA(15) = ABS(WO)
        ELSE IF(ABS(Z).GT.224.00.AND.
     1    (MAT.EQ.0.OR.MAT.EQ.7))THEN
          IF(COSTHETA(16).EQ.0)COSTHETA(16) = ABS(WO)
        ENDIF
      ENDIF


      RETURN
      END
C======================================================================C
      SUBROUTINE TRACKSTEPS
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f'
      INCLUDE 'ucnapenmain.h'
      
      IF(IBODY.EQ.432)THEN
           IDCE = IDCE + 1
      ELSE IF(IBODY.EQ.433)THEN
           IBEE = IBEE + 1
      ELSE IF(IBODY.EQ.434)THEN
           IBEW = IBEW + 1
      ELSE IF(IBODY.EQ.435)THEN
           IDCW = IDCW + 1
      ELSE IF(IBODY.EQ.429)THEN
           IDDE = IDDE + 1
      ELSE IF(IBODY.EQ.415)THEN
           IMYBE = IMYBE + 1
      ELSE IF(IBODY.EQ.414)THEN
           IMYFE = IMYFE + 1
      ELSE IF(IBODY.EQ.411)THEN
           IDDW = IDDW + 1
      ELSE IF(IBODY.EQ.397)THEN
           IMYBW = IMYBW + 1
      ELSE IF(IBODY.EQ.396)THEN
           IMYFW = IMYFW + 1
      ENDIF
      
      RETURN
      END
C=-===================================================================C
      SUBROUTINE WRITESTEPS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f'
      INCLUDE 'ucnapenmain.h'
      
      WRITE(47,'(1F10.3,1X,5I8,1X,1F10.3,1X,5I8)') 
     1 DEBO(416),IBEE,IDCE,IMYFE,IMYBE,IDDE,
     1 DEBO(398),IBEW,IDCW,IMYFW,IMYBW,IDDW
      
      RETURN
      END
