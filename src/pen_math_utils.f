C  *********************************************************************
C                       SUBROUTINE GCONE
C  *********************************************************************
      SUBROUTINE GCONE(UF,VF,WF)
C
C  This subroutine samples a random direction uniformly within a cone
C  with central axis in the direction (THETA,PHI) and aperture ALPHA.
C  Parameters are initialised by calling subroutine GCONE0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      EXTERNAL RAND
C  ****  Define a direction relative to the z-axis.
      WT=CAPER+(1.0D0-CAPER)*RAND(1.0D0)
      DF=TWOPI*RAND(2.0D0)
      SUV=SQRT(1.0D0-WT*WT)
      UT=SUV*COS(DF)
      VT=SUV*SIN(DF)
C  **** Rotate to the beam axis direction.
      UF=CPCT*UT-SPHI*VT+CPST*WT
      VF=SPCT*UT+CPHI*VT+SPST*WT
      WF=-STHE*UT+CTHE*WT
C  ****  Ensure normalisation.
      DXY=UF*UF+VF*VF
      DXYZ=DXY+WF*WF
      IF(ABS(DXYZ-1.0D0).GT.1.0D-14) THEN
        FNORM=1.0D0/SQRT(DXYZ)
        UF=FNORM*UF
        VF=FNORM*VF
        WF=FNORM*WF
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE GCONE0
C  *********************************************************************
      SUBROUTINE GCONE0(THETA,PHI,ALPHA)
C
C  This subroutine defines the parameters for sampling random directions
C  uniformly within a cone with axis in the direction (THETA,PHI) and
C  aperture ALPHA (in rad).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C  ****  Parameters for sampling directions within a cone.
      COMMON/CGCONE/CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER
C
      CPCT=COS(PHI)*COS(THETA)
      CPST=COS(PHI)*SIN(THETA)
      SPCT=SIN(PHI)*COS(THETA)
      SPST=SIN(PHI)*SIN(THETA)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      STHE=SIN(THETA)
      CTHE=COS(THETA)
      CAPER=COS(ALPHA)
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE RDPSF
C  *********************************************************************
      SUBROUTINE RDPSF(IUNIT,NSHI,ISEC,KODEPS)
C
C  This subroutine reads the phase-space file (psf). Blank lines and
C  lines starting with the pound sign (#) are considered as comments and
C  are ignored.
C
C  The reading of a new psf is initiated by calling subroutine RDPSF
C  with KODEPS=0. If the psf contains some particles, the oputput value
C  of KODEPS is set equal to 1. In subsequent calls to subroutine RDPSF,
C  the output value KODEPS=1 indicates that a valid record has been read
C  and that the state variables of a particle have been loaded in common
C  TRACK. The output value of ISEC is 1 if the psf contains more
C  particles from the current shower, and 0 otherwise. That is, when
C  ISEC=0 the next particle in the psf belongs to a new shower. The
C  output value KODEPS=-1 indicates that all particle records in the psf
C  have been read (or that the psf was empty). Note that when KODEPS=-1,
C  the contents of common TRACK remains unchanged.
C
C  To reduce disk-access time, particle records are read in bunches and
C  stored in a buffer. The buffer size (number of particles in a bunch)
C  is defined by the parameter NRBUFF.
C
C  To verify that the input psfs are read correctly, uncomment the lines
C  with 'CALL WRPSF' (which write the read particle records in an output
C  file, UNIT=13). You should also add a sentence 'CALL WRPSF(13,0,1)'
C  at the end of the main program, which forces the writting of particle
C  records that may remain in the buffer and closes the output file.
C  This file should contain all the records read from the psf files.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      PARAMETER (NRBUFF=1000)  ! Buffer size.
      CHARACTER PSFR(NRBUFF)*136
      SAVE PSFR
C  ****  Main-PENELOPE common.
      PARAMETER (MAXMAT=10)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C
      SAVE E0,X0,Y0,Z0,U0,V0,W0,WGHT0,KPAR0,ILB10,ILB20,ILB30,ILB40,
     1  NSHI0,ILAST,IRD,IEOF
C
C  ****  Initialisation (first call after opening the psf).
C
      IF(KODEPS.EQ.0) THEN
        NSHI=0
        ISEC=0
        IEOF=0
        DO I=1,NRBUFF
         PSFR(I)='#'  ! Identifies empty records.
        ENDDO
        READ(IUNIT,'(A)',END=10,ERR=10) (PSFR(I),I=1,NRBUFF)
        GO TO 11
   10   CONTINUE
        IEOF=1
   11   CONTINUE
C
        IRD=0
   20   CONTINUE
        IF(IRD.GE.NRBUFF) THEN
          IF(IEOF.EQ.1) THEN
            KODEPS=-1  ! End of file.
            RETURN
          ENDIF
          DO I=1,NRBUFF
           PSFR(I)='#'  ! Identifies empty records.
          ENDDO
          READ(IUNIT,'(A)',END=21,ERR=21) (PSFR(I),I=1,NRBUFF)
          GO TO 22
   21     CONTINUE
          IEOF=1
   22     CONTINUE
          IRD=0
        ENDIF
        IRD=IRD+1
        READ(PSFR(IRD),*,ERR=20,END=20) KPAR0,E0,X0,Y0,Z0,U0,V0,W0,
     1    WGHT0,ILB10,ILB20,ILB30,ILB40,NSHI0
        ILAST=0
        KODEPS=1
        RETURN
      ENDIF
C
C  ****  The last particle in the psf has already been loaded.
C
      IF(ILAST.EQ.1) THEN
        NSHI=0
        ISEC=0
        KODEPS=-1  ! End of file.
C       CALL WRPSF(13,NSHI,-1)  ! Uncomment to check psf reading.
        RETURN
      ENDIF
C
C  ****  A read particle is loaded.
C
      E=E0
      X=X0
      Y=Y0
      Z=Z0
      U=U0
      V=V0
      W=W0
      WGHT=WGHT0
      KPAR=KPAR0
      ILB(1)=ILB10
      ILB(2)=ILB20
      ILB(3)=ILB30
      ILB(4)=ILB40
      NSHI=NSHI0
C     CALL WRPSF(13,NSHI,0)  ! Uncomment to check psf reading.
      KODEPS=1
C
C  ****  Read a new particle and keep its state variables in memory.
C        This is needed to identify possible secondary particles that
C        belong to the same shower.
C
   30 CONTINUE
      IF(IRD.GE.NRBUFF) THEN
        IF(IEOF.EQ.1) THEN
          ILAST=1
          ISEC=0
          RETURN
        ENDIF
        DO I=1,NRBUFF
         PSFR(I)='#'  ! Identifies empty records.
        ENDDO
        READ(IUNIT,'(A)',END=31,ERR=31) (PSFR(I),I=1,NRBUFF)
        GO TO 32
   31   CONTINUE
        IEOF=1
   32   CONTINUE
        IRD=0
      ENDIF
      IRD=IRD+1
      READ(PSFR(IRD),*,ERR=30,END=30) KPAR0,E0,X0,Y0,Z0,U0,V0,W0,
     1    WGHT0,ILB10,ILB20,ILB30,ILB40,NSHI0
C
      IF(NSHI0.EQ.0) THEN
        ISEC=1
      ELSE
        ISEC=0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE WRPSF
C  *********************************************************************
      SUBROUTINE WRPSF(IPSFO,NSHI,ICLOSE)
C
C  This subroutine writes particle records in the output phase-space
C  file (psf). To reduce disk-access time, particle records are
C  collected in a buffer and moved to the psf when the buffer is full.
C  The buffer size (number of records) is set by the parameter NRBUFF.
C
C  Particles that remain in the buffer at the end of a simulation run
C  must be transferred to the psf by calling subroutine WRPSF with
C  ICLOSE=1.
C
C  Input parameters:
C     IPSFO ..... output unit.
C     NSHI ...... incremental shower number.
C     ICLOSE .... closing flag, acts only when its value is not 0.
C                 When ICLOSE .GT. 0, particles remaining in the buffer
C                 are transferred to the psf, and the output unit IPSFO
C                 is closed.
C                 If ICLOSE .LT. 0, particles in the buffer are moved to
C                 the psf, but the output unit remains open.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*10 LC10A,LC10B
C
      PARAMETER (NRBUFF=1000)  ! Buffer size.
      CHARACTER PSFR(NRBUFF)*136
      SAVE PSFR
C  ****  Main-PENELOPE common.
      PARAMETER (MAXMAT=10)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C
      DATA NREC/0/
      SAVE NREC
C
      IF(ICLOSE.NE.0) THEN
        IF(NREC.GT.0) THEN
          DO I=1,NREC
            NCHAR=INDEX(PSFR(I),'*')-1
            WRITE(IPSFO,'(A)') PSFR(I)(1:NCHAR)
          ENDDO
        ENDIF
        IF(ICLOSE.GT.0) CLOSE(IPSFO)
        NREC=0
        RETURN
      ENDIF
C
      NREC=NREC+1
      CALL N2CH10(ILB(4),LC10A,NDIGA)
      CALL N2CH10(NSHI,LC10B,NDIGB)
      WRITE(PSFR(NREC),'(I2,1P,8E13.5,I3,I2,I2,1X,A,1X,A,A)')
     1  KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),
     2  LC10A(1:NDIGA),LC10B(1:NDIGB),'*'
C
      IF(NREC.EQ.NRBUFF) THEN
        WRITE(IPSFO,'(A)') (PSFR(I)(1:INDEX(PSFR(I),'*')-1),I=1,NREC)
        NREC=0
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE N2CH10
C  *********************************************************************
      SUBROUTINE N2CH10(N,L,NDIG)
C
C  This subroutine writes an integer number N in a 10-character string
C  L. The number is written at the left, followed by unused blanks. NDIG
C  is the number of decimal digits of N.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*10 L,LT
C
      WRITE(LT,'(I10)') N
      DO I=1,10
        IF(LT(I:I).NE.' ') THEN
          IT=I-1
          GO TO 1
        ENDIF
      ENDDO
      IT=9
    1 CONTINUE
      L=LT(IT+1:10)
      NDIG=10-IT
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SOURCE
C  *********************************************************************
      SUBROUTINE SOURCE
C
C  Definition of the primary radiation source provided by the user. When
C  KPARP=0, the initial states of primary particles are set by this
C  subroutine, which may use or supersede other definitions in the input
C  file.
C
C  The first call to subroutine SOURCE (when the value of IDEF is equal
C  to 0) should set _all_ the source parameters that are not defined
C  through the input file of PENMAIN. In particular EPMAX (highest
C  energy of the primary particles) should be defined at the first call.
C  This parameter is needed for the initialization of PENELOPE.
C
C  Each subsequent call to SOURCE sets the initial state of a shower,
C  which may consist of several primary particles (e.g., emitted in a
C  cascade of nuclear transitions). The state variables E, X,Y,Z, U,V,W,
C  WGHT, KPAR and ILB(5) of each primary particle are loaded in common
C  TRACK. In the case of polarized photons, the polarization state
C  (Stokes parameters) and the IPOL flag are loaded in common STOKES.
C  Accompanying primary particles are sent to the secondary stack (by
C  calling subroutine STORES of the PENELOPE package); these particles
C  must be marked with ILB(1)=-1 to be treated as proper primary
C  particles by subroutine SHOWER. Since subroutine STORES gets
C  information through the common blocks TRACK and STOKES, the
C  accompanying particles should be defined first. The last particle
C  defined, just before returning control to the calling program, is
C  loaded directly in commons TRACK and STOKES.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f' ! Includes the following two common blocks:
C     COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
C     COMMON/STOKES/SP1,SP2,SP3,IPOL  ! Do not remove the comments!
      DIMENSION ILBS(5)
      EXTERNAL RAND
C  ****  Save IDEF and all parameters needed to set the initial states.
      SAVE IDEF,E1,E2
      DATA IDEF/0/
C
C  ****  Define the source characteristics.
C
      IF(IDEF.EQ.0) THEN  ! Example: 60-Co source, two gamma rays.
        E1=1.17E6
        E2=1.33E6
C        EPMAX=E1+E2+1.0D0  ! Note: for positrons add 5.12D5.
        IDEF=1
        RETURN
      ENDIF
C
C  ************  Set the initial states of a new shower.
C
C  ----  Initial position ...
      IF(LEXSRC) THEN
        IF(LEXBD) THEN
          NTRIAL=0
   10     CONTINUE
          X=SX0+(RAND(1.0D0)-0.5D0)*SSX
          Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
          Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
          CALL LOCATE
          NTRIAL=NTRIAL+1
          IF(NTRIAL.GT.200) THEN
            WRITE(26,'(3X,''WARNING: the sampling of initial '',
     1        ''positions may be very inefficient.'')')
            WRITE(6,'(3X,''WARNING: the sampling of initial '',
     1        ''positions may be very inefficient.'')')
          ENDIF
          IF(IXSBOD(IBODY).EQ.0) GO TO 10
        ELSE
          X=SX0+(RAND(1.0D0)-0.5D0)*SSX
          Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
          Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
        ENDIF
      ELSE
        X=SX0
        Y=SY0
        Z=SZ0
      ENDIF
C
C  ************  First gamma ray (stored in the secondary stack).
C
      KPARS=2
      ES=E1
C  ----  Initial direction ...
      IF(LSCONE) THEN
        CALL GCONE(US,VS,WS)  ! Conical beam.
      ELSE
        WS=CTHL+RAND(4.0D0)*DCTH  ! Rectangular beam.
        UV=SQRT(1.0D0-WS*WS)
        PHI=PHIL+RAND(5.0D0)*DPHI
        US=UV*COS(PHI)
        VS=UV*SIN(PHI)
      ENDIF
      WGHTS=1.0D0
C  ****  The flag ILBS(1)=-1 indicates that the particle is a primary
C  one, i.e., that its energy does not have to be subtracted from the
C  site.
      ILBS(1)=-1
      ILBS(2)=0
      ILBS(3)=0
      ILBS(4)=0
      ILBS(5)=0
      IPOL=0
      CALL STORES(ES,X,Y,Z,US,VS,WS,WGHTS,KPARS,ILBS,IPOL)
C
C  ************  Second gamma ray (ready to be tracked).
C
      KPAR=2
      E=E2
C  ----  Initial direction ...
      IF(LSCONE) THEN
        CALL GCONE(U,V,W)  ! Conical beam.
      ELSE
        W=CTHL+RAND(4.0D0)*DCTH  ! Rectangular beam.
        UV=SQRT(1.0D0-W*W)
        PHI=PHIL+RAND(5.0D0)*DPHI
        U=UV*COS(PHI)
        V=UV*SIN(PHI)
      ENDIF
      WGHT=1.0D0
      ILB(1)=1  ! Identifies primary particles.
      ILB(2)=0
      ILB(3)=0
      ILB(4)=0
      ILB(5)=0
      IPOL=0
      RETURN
      END
C **********************************************************************
C                      SUBROUTINE DOSEBOX
C **********************************************************************
      SUBROUTINE DOSEBOX(DEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f'
      
      IF(LDOSEM) THEN  ! The particle is inside the dose box.
        IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
     1     (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
     1     (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
          I1=INT(1.0D0+(X-DXL(1))*RBDOSE(1))
          I2=INT(1.0D0+(Y-DXL(2))*RBDOSE(2))
          I3=INT(1.0D0+(Z-DXL(3))*RBDOSE(3))
          DOSP=DEP*RHOI(MAT)
          IF(N.NE.LDOSE(I1,I2,I3)) THEN
            DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
            DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
            DOSEP(I1,I2,I3)=DOSP
            LDOSE(I1,I2,I3)=N
          ELSE
            DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DOSP
          ENDIF
C
          IF(N.NE.LDDOSE(I3)) THEN
            DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
            DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
            DDOSEP(I3)=DOSP
            LDDOSE(I3)=N
          ELSE
            DDOSEP(I3)=DDOSEP(I3)+DOSP
          ENDIF
        ENDIF
      ENDIF
      
      RETURN
      END
C *********************************************************************
C              SUBROUTINE IMPACT_DETECTOR2
C *********************************************************************
      SUBROUTINE IMPACT_DETECTOR2(IBODYL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f'
C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Impact detectors.
      IF(IBODYL.GT.NBODY)RETURN
      IDET=KDET(IBODY)
      IF(IDET.NE.0) THEN
        IF(KDET(IBODYL).NE.IDET.AND.KKDI(IDET,KPAR).EQ.1) THEN
C
          IF(IPSF(IDET).EQ.1) THEN
            NSHJ=INT(SHN-RLAST)
            CALL WRPSF(IPSFO,NSHJ,0)
            RWRITE=RWRITE+1.0D0
            RLAST=SHN
          ENDIF
C
          DEDI(IDET)=DEDI(IDET)+E*WGHT
          IF(LDILOG(IDET)) THEN
            IE=INT(1.0D0+(LOG(E)-EDILL(IDET))*RBDIEL(IDET))
          ELSE
            IE=INT(1.0D0+(E-EDIL(IDET))*RBDIE(IDET))
          ENDIF
          IF(IE.GT.0.AND.IE.LE.NDICH(IDET)) THEN
            IF(N.NE.LDIT(IDET,IE)) THEN
              DIT(IDET,IE)=DIT(IDET,IE)+DITP(IDET,IE)
              DIT2(IDET,IE)=DIT2(IDET,IE)+DITP(IDET,IE)**2
              DITP(IDET,IE)=WGHT
              LDIT(IDET,IE)=N
            ELSE
              DITP(IDET,IE)=DITP(IDET,IE)+WGHT
            ENDIF
            IF(N.NE.LDIP(IDET,IE,KPAR)) THEN
              DIP(IDET,IE,KPAR)=DIP(IDET,IE,KPAR)+DIPP(IDET,IE,KPAR)
              DIP2(IDET,IE,KPAR)=
     1            DIP2(IDET,IE,KPAR)+DIPP(IDET,IE,KPAR)**2
              DIPP(IDET,IE,KPAR)=WGHT
              LDIP(IDET,IE,KPAR)=N
            ELSE
              DIPP(IDET,IE,KPAR)=DIPP(IDET,IE,KPAR)+WGHT
            ENDIF
          ENDIF
!           IF(IDCUT(IDET).EQ.0) THEN
!             DEBO(IBODY)=DEBO(IBODY)+E*WGHT
!             IEXIT=3
!             GO TO 104
!           ENDIF
        ENDIF
      ENDIF
      
      RETURN
      END
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
C **********************************************************************
C        SUBROUTINE ENERGYFLUENCE
C **********************************************************************
      SUBROUTINE ENERGYFLUENCE(DSEF,IBODYL)
c---------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f'
C  ----  Energy distribution of fluence.
      IF(IBODYL.GT.NBODY)RETURN
      IDET=KDET(IBODYL)
      IF(IDET.NE.0) THEN
        IF(IDCUT(IDET).EQ.2) THEN
          IF(KKDI(IDET,KPAR).EQ.1) THEN
            IF(LDILOG(IDET)) THEN
              IE=INT(1.0D0+(LOG(E)-EDILL(IDET))*RBDIEL(IDET))
            ELSE
              IE=INT(1.0D0+(E-EDIL(IDET))*RBDIE(IDET))
            ENDIF
            IF(IE.GT.0.AND.IE.LE.NDICH(IDET)) THEN
              IF(N.NE.LFST(IDET,IE)) THEN
                FST(IDET,IE)=FST(IDET,IE)+FSTP(IDET,IE)
                FST2(IDET,IE)=FST2(IDET,IE)+FSTP(IDET,IE)**2
                FSTP(IDET,IE)=WGHT*DSEF
                LFST(IDET,IE)=N
              ELSE
                FSTP(IDET,IE)=FSTP(IDET,IE)+WGHT*DSEF
              ENDIF
              IF(N.NE.LFSP(IDET,IE,KPAR)) THEN
                FSP(IDET,IE,KPAR)=FSP(IDET,IE,KPAR)+FSPP(IDET,IE,KPAR)
                FSP2(IDET,IE,KPAR)=
     1              FSP2(IDET,IE,KPAR)+FSPP(IDET,IE,KPAR)**2
                FSPP(IDET,IE,KPAR)=WGHT*DSEF
                LFSP(IDET,IE,KPAR)=N
              ELSE
                FSPP(IDET,IE,KPAR)=FSPP(IDET,IE,KPAR)+WGHT*DSEF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      
      
      RETURN
      END
C **********************************************************************
      SUBROUTINE COLLECTEMERGE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      
C  ****  Energy distribution of emerging particles.
        K=1.0D0+(E-EMIN)*RBSE
        IF(K.GT.0.AND.K.LE.NBE) THEN
          IF(N.NE.LPDE(KPAR,IEXIT,K)) THEN
            PDE(KPAR,IEXIT,K)=PDE(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)
            PDE2(KPAR,IEXIT,K)=
     1        PDE2(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)**2
            PDEP(KPAR,IEXIT,K)=WGHT
            LPDE(KPAR,IEXIT,K)=N
          ELSE
            PDEP(KPAR,IEXIT,K)=PDEP(KPAR,IEXIT,K)+WGHT
          ENDIF
        ENDIF
C  ****  Angular distribution of emerging particles.
        THETA=ACOS(W)
        KTH=1.0D0+THETA*RA2DE*RBSTH
        IF(ABS(U).GT.1.0D-16) THEN  ! Azimuthal bin number corrected.
           PHI=ATAN2(V,U)
        ELSE IF(ABS(V).GT.1.0D-16) THEN
           PHI=ATAN2(V,U)
        ELSE
           PHI=0.0D0
        ENDIF
        IF(PHI.LT.0.0D0) PHI=TWOPI+PHI
        KPH=1.0D0+PHI*RA2DE*RBSPH
        IF(N.NE.LPDA(KPAR,KTH,KPH)) THEN
          PDA(KPAR,KTH,KPH)=PDA(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)
          PDA2(KPAR,KTH,KPH)=PDA2(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)**2
          PDAP(KPAR,KTH,KPH)=WGHT
          LPDA(KPAR,KTH,KPH)=N
        ELSE
          PDAP(KPAR,KTH,KPH)=PDAP(KPAR,KTH,KPH)+WGHT
        ENDIF
        
      RETURN
      END
C **********************************************************************
C        SUBROUTINE TALLYSPECTRA
C **********************************************************************
      SUBROUTINE TALLYSPECTRA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4(I-N)
      INCLUDE 'pmcomms.f'
      
      IF(NDEDEF.GT.0) THEN
        DO KD=1,NDEDEF
          DEDE(KD)=0.0D0
        ENDDO
        DO KB=1,NBODY
          IDET=KBDE(KB)
          IF(IDET.NE.0) THEN
            DEDE(IDET)=DEDE(IDET)+DEBO(KB)
          ENDIF
        ENDDO
C
        DO KD=1,NDEDEF
          TDED(KD)=TDED(KD)+DEDE(KD)
          TDED2(KD)=TDED2(KD)+DEDE(KD)**2
          IF(DEDE(KD).GT.1.0D-5) THEN
            IF(LDELOG(KD)) THEN
              KE=1.0D0+(LOG(DEDE(KD))-EDELL(KD))*RBDEEL(KD)
            ELSE
              KE=1.0D0+(DEDE(KD)-EDEL(KD))*RBDEE(KD)
            ENDIF
            IF(KE.GT.0.AND.KE.LE.NDECH(KD)) THEN
              DET(KD,KE)=DET(KD,KE)+1.0D0
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      DO KB=1,NBODY
        TDEBO(KB)=TDEBO(KB)+DEBO(KB)
        TDEBO2(KB)=TDEBO2(KB)+DEBO(KB)**2
      ENDDO
C  --  Average energies 'collected' by impact detectors.
      IF(NDIDEF.GT.0) THEN
        DO KD=1,NDIDEF
          TDID(KD)=TDID(KD)+DEDI(KD)
          TDID2(KD)=TDID2(KD)+DEDI(KD)**2
        ENDDO
      ENDIF
C  --  Final state counters.
      DO I=1,3
        PRIM(I)=PRIM(I)+DPRIM(I)
        PRIM2(I)=PRIM2(I)+DPRIM(I)**2
        DO K=1,3
          SEC(K,I)=SEC(K,I)+DSEC(K,I)
          SEC2(K,I)=SEC2(K,I)+DSEC(K,I)**2
        ENDDO
      ENDDO
      
      RETURN
      END
      
