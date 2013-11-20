C  *********************************************************************
C                       SUBROUTINE PMWRT
C  *********************************************************************
      SUBROUTINE PMWRT(ICLOSE)
C
C  Computes averages and writes results in output files.
C
C  ICLOSE is a closing flag, which is used only when the phase-space
C  file of an impact detector is being generated.
C  -- When ICLOSE .GT. 0, particles remaining in the buffer are
C     transferred to the psf, and the psf output unit is closed.
C  -- If ICLOSE .LT. 0, particles are moved to the psf, but the psf
C     output unit remains open.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f'
      DIMENSION WSEC(3,3),WSEC2(3,3)
C
C  ************  Transfer partial counters to global counters.
C
      DO KPAR=1,3
      DO IEXIT=1,2
      DO IE=1,NBE
        PDE(KPAR,IEXIT,IE)=PDE(KPAR,IEXIT,IE)+PDEP(KPAR,IEXIT,IE)
        PDE2(KPAR,IEXIT,IE)=PDE2(KPAR,IEXIT,IE)+PDEP(KPAR,IEXIT,IE)**2
        PDEP(KPAR,IEXIT,IE)=0.0D0
        LPDE(KPAR,IEXIT,IE)=0
      ENDDO
      ENDDO
      ENDDO
C
      DO KPAR=1,3
      DO KTH=1,NBTH
      DO KPH=1,NBPH
        PDA(KPAR,KTH,KPH)=PDA(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)
        PDA2(KPAR,KTH,KPH)=PDA2(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)**2
        PDAP(KPAR,KTH,KPH)=0.0D0
        LPDA(KPAR,KTH,KPH)=0
      ENDDO
      ENDDO
      ENDDO
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
        DO J=1,NDICH(ID)
          DIT(ID,J)=DIT(ID,J)+DITP(ID,J)
          DIT2(ID,J)=DIT2(ID,J)+DITP(ID,J)**2
          DITP(ID,J)=0.0D0
          LDIT(ID,J)=0
        ENDDO
        DO J=1,NDICH(ID)
        DO K=1,3
          DIP(ID,J,K)=DIP(ID,J,K)+DIPP(ID,J,K)
          DIP2(ID,J,K)=DIP2(ID,J,K)+DIPP(ID,J,K)**2
          DIPP(ID,J,K)=0.0D0
          LDIP(ID,J,K)=0
        ENDDO
        ENDDO
        ENDDO
      ENDIF
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
        DO J=1,NDICH(ID)
          FST(ID,J)=FST(ID,J)+FSTP(ID,J)
          FST2(ID,J)=FST2(ID,J)+FSTP(ID,J)**2
          FSTP(ID,J)=0.0D0
          LFST(ID,J)=0
        ENDDO
        DO J=1,NDICH(ID)
        DO K=1,3
          FSP(ID,J,K)=FSP(ID,J,K)+FSPP(ID,J,K)
          FSP2(ID,J,K)=FSP2(ID,J,K)+FSPP(ID,J,K)**2
          FSPP(ID,J,K)=0.0D0
          LFSP(ID,J,K)=0
        ENDDO
        ENDDO
        ENDDO
      ENDIF
C
      IF(LDOSEM) THEN
        DO I3=1,NDB(3)
        DO I2=1,NDB(2)
        DO I1=1,NDB(1)
          DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
          DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
          DOSEP(I1,I2,I3)=0.0D0
          LDOSE(I1,I2,I3)=0
        ENDDO
        ENDDO
          DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
          DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
          DDOSEP(I3)=0.0D0
          LDDOSE(I3)=0
        ENDDO
      ENDIF
C
      IF(NPSFO.GT.0) THEN
        CALL WRPSF(IPSFO,0,ICLOSE)
      ENDIF
C
C  ************  If 'DUMPTO' is active, write counters in a dump file.
C
      IF(LDUMP) THEN
        OPEN(9,FILE=PFILED)
        WRITE(9,*) SHN,TSIM
        WRITE(9,'(A65)') TITLE
        WRITE(9,*) ISEED1,ISEED2
        WRITE(9,*) NPSN,RLREAD
        WRITE(9,*) (SHIST(I),I=1,NSEB)
        WRITE(9,*) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        WRITE(9,*) ((SEC(K,I),I=1,3),K=1,3),((SEC2(K,I),I=1,3),K=1,3)
        WRITE(9,*) (TDEBO(I),I=1,NBODY), (TDEBO2(I),I=1,NBODY)
        WRITE(9,*) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        WRITE(9,*) (((PDA(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3),
     1             (((PDA2(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3)
        WRITE(9,*) (TDID(I),I=1,NIDM), (TDID2(I),I=1,NIDM)
        WRITE(9,*) (TDED(I),I=1,NEDM), (TDED2(I),I=1,NEDM)
        IF(NDIDEF.GT.0) THEN
          WRITE(9,*) RLAST,RWRITE
          DO ID=1,NDIDEF
            WRITE(9,*) (DIT(ID,J),J=1,NDICH(ID))
            WRITE(9,*) (DIT2(ID,J),J=1,NDICH(ID))
            WRITE(9,*) ((DIP(ID,J,K),J=1,NDICH(ID)),K=1,3)
            WRITE(9,*) ((DIP2(ID,J,K),J=1,NDICH(ID)),K=1,3)
            IF(IDCUT(ID).EQ.2) THEN
              WRITE(9,*) (FST(ID,J),J=1,NDICH(ID))
              WRITE(9,*) (FST2(ID,J),J=1,NDICH(ID))
              WRITE(9,*) ((FSP(ID,J,K),J=1,NDICH(ID)),K=1,3)
              WRITE(9,*) ((FSP2(ID,J,K),J=1,NDICH(ID)),K=1,3)
            ENDIF
          ENDDO
        ENDIF
        IF(NDEDEF.GT.0) THEN
          DO ID=1,NDEDEF
            WRITE(9,*) (DET(ID,J),J=1,NDECH(ID))
          ENDDO
        ENDIF
        IF(LDOSEM) THEN
          WRITE(9,*)
     1      (((DOSE(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1)),
     1      (((DOSE2(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
          WRITE(9,*) (DDOSE(I3),I3=1,NDB(3)),(DDOSE2(I3),I3=1,NDB(3))
        ENDIF
        WRITE(9,'(/3X,''*** END ***'')')
        CLOSE(9)
      ENDIF
C
C  ------------------------  Print simulation results.
C
C     IEXIT: 1=upbound, 2=downbound, 3=absorbed.
C

      OPEN(27,FILE='data/penmain-res.dat')
      WRITE(27,3000)
 3000 FORMAT(//3X,35('*')/3X,'**   Program PENMAIN. Results.   **',
     1  /3X,35('*'))
C
      WRITE(27,1001) DATE23
 1001 FORMAT(/3X,'Date and time: ',A23)
      WRITE(27,'(/3X,A65)') TITLE
C
      WRITE(27,3001) TSIM
 3001 FORMAT(/3X,'Simulation time ......................... ',
     1  1P,E13.6,' sec')
      TAVS=SHN/TSIM
      WRITE(27,3002) TAVS
 3002 FORMAT(3X,'Simulation speed ........................ ',
     1  1P,E13.6,' showers/sec')
      WRITE(27,3003) SHN
 3003 FORMAT(//3X,'Simulated primary showers ............... ',
     1  1P,E13.6)
C
      IF(KPARP.EQ.1) WRITE(27,1110)
 1110 FORMAT(/3X,'Primary particles: electrons')
      IF(KPARP.EQ.2) WRITE(27,1111)
 1111 FORMAT(/3X,'Primary particles: photons')
      IF(KPARP.EQ.3) WRITE(27,1112)
 1112 FORMAT(/3X,'Primary particles: positrons')
      IF(KPARP.EQ.0) WRITE(26,1113)
 1113 FORMAT(/3X,'Primary particles: set by the user subroutine SOURCE')
C
      WRITE(27,3004) PRIM(1)
 3004 FORMAT(/3X,'Upbound primary particles ............... ',
     1  1P,E13.6)
      WRITE(27,3005) PRIM(2)
 3005 FORMAT(3X,'Downbound primary particles ............. ',
     1  1P,E13.6)
      WRITE(27,3006) PRIM(3)
 3006 FORMAT(3X,'Absorbed primary particles .............. ',
     1  1P,E13.6)
C
      FNT=1.0D0/SHN
      FT=(PRIM(1)+SEC(KPARP,1))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(1)-PRIM(1)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,1)-SEC(KPARP,1)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(27,3007) FT,ERR
 3007 FORMAT(/3X,'Upbound fraction ................... ',
     1  1P,E13.6,' +-',E8.1)
      FB=(PRIM(2)+SEC(KPARP,2))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(2)-PRIM(2)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,2)-SEC(KPARP,2)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(27,3008) FB,ERR
 3008 FORMAT(3X,'Downbound fraction ................. ',
     1  1P,E13.6,' +-',E8.1)
      FA=PRIM(3)*FNT
      ERR=3.0D0*FNT*SQRT(ABS(PRIM2(3)-PRIM(3)**2*FNT))
      WRITE(27,3009) FA,ERR
 3009 FORMAT(3X,'Absorption fraction ................ ',
     1  1P,E13.6,' +-',E8.1)
C
      DO K=1,3
        DO I=1,3
          WSEC2(K,I)=3.0D0*FNT*SQRT(ABS(SEC2(K,I)-SEC(K,I)**2*FNT))
          WSEC(K,I)=SEC(K,I)*FNT
        ENDDO
      ENDDO
      WRITE(27,3010)
     1  WSEC(1,1),WSEC(2,1),WSEC(3,1),WSEC2(1,1),WSEC2(2,1),WSEC2(3,1),
     1  WSEC(1,2),WSEC(2,2),WSEC(3,2),WSEC2(1,2),WSEC2(2,2),WSEC2(3,2),
     1  WSEC(1,3),WSEC(2,3),WSEC(3,3),WSEC2(1,3),WSEC2(2,3),WSEC2(3,3)
 3010 FORMAT(/3X,'Secondary-particle generation probabilities:',
     1  /19X,46('-'),
     1  /19X,'|  electrons   |   photons    |  positrons   |',1P,
     1  /3X,62('-')/3X,'|   upbound     |',3(E13.6,1X,'|'),
     1  /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1  /3X,62('-')/3X,'|   downbound   |',3(E13.6,1X,'|'),
     1  /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1  /3X,62('-')/3X,'|   absorbed    |',3(E13.6,1X,'|'),
     1  /3X,'|               |',3('  +-',E8.1,2X,'|'),
     1  /3X,62('-'))
C
C  ****  Average energies deposited in bodies..
C
      DF=1.0D0/SHN
      WRITE(27,3011)
 3011 FORMAT(/3X,'Average deposited energies (bodies):')
      DO KB=1,NBODY
        IF(MATER(KB).NE.0) THEN
          QER=3.0D0*DF*SQRT(ABS(TDEBO2(KB)-TDEBO(KB)**2*DF))
          QAV=TDEBO(KB)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(27,3012) KB,QAV,QER,EFFIC
        ENDIF
      ENDDO
 3012 FORMAT(6X,'Body ',I4, ' ...... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
C
C  ****  Average energies 'collected' by impact detectors.
C
      IF(NDIDEF.GT.0) THEN
        WRITE(27,3013)
 3013   FORMAT(/3X,'Average incoming energies (impact detectors):')
        DO KD=1,NDIDEF
          QER=3.0D0*DF*SQRT(ABS(TDID2(KD)-TDID(KD)**2*DF))
          QAV=TDID(KD)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(27,3014) KD,QAV,QER,EFFIC
        ENDDO
      ENDIF
 3014 FORMAT(6X,'Detector #',I2,' ... ',1P,E13.6,' +-',E8.1,' eV',4X,
     1  '(effic. =',E9.2,')')
C
C  ****  Average deposited energies in energy-deposition detectors.
C
      IF(NDEDEF.GT.0) THEN
        WRITE(27,3015)
 3015   FORMAT(/3X,'Average deposited energies (energy detectors):')
        DO KD=1,NDEDEF
          QER=3.0D0*DF*SQRT(ABS(TDED2(KD)-TDED(KD)**2*DF))
          QAV=TDED(KD)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(27,3014) KD,QAV,QER,EFFIC
        ENDDO
      ENDIF
C
C  ****  Detector efficiencies.
C
      IF(NDIDEF.GT.0) THEN
        WRITE(27,3016)
 3016   FORMAT(/3X,'Detection efficiencies (impact detectors):')
        DO KD=1,NDIDEF
          SUM=0.0D0
          ESUM=0.0D0
          DO J=1,NDICH(KD)
            YERR=3.0D0*SQRT(ABS(DIT2(KD,J)-DIT(KD,J)**2*DF))
            YAV=MAX(DIT(KD,J)*DF,1.0D-35)
            YERR=MAX(YERR*DF,1.0D-35)
            SUM=SUM+YAV
            ESUM=ESUM+YERR
          ENDDO
          WRITE(27,3017) KD,SUM,ESUM
        ENDDO
      ENDIF
 3017 FORMAT(6X,'Detector #',I2,' ... ',1P,E13.6,' +-',E8.1)
C
      IF(NDEDEF.GT.0) THEN
        WRITE(27,3018)
 3018   FORMAT(/3X,'Detection efficiencies (energy detectors):')
        DO KD=1,NDEDEF
          SUM=0.0D0
          ESUM=0.0D0
          DO J=1,NDECH(KD)
            YAV=DET(KD,J)*DF
            YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV))*DF)
            SUM=SUM+YAV
            ESUM=ESUM+YERR
          ENDDO
          WRITE(27,3017) KD,SUM,ESUM
        ENDDO
      ENDIF
C
C  ************  Energy spectrum of the source.
C
      IF(LSPEC) THEN
        OPEN(9,FILE='data/psource.dat')
        WRITE(9,9000)
 9000 FORMAT(
     1  1X,'#  Results from PENMAIN. ',
     1 /1X,'#  Source energy spectrum.',
     1 /1X,'#  1st column: E (eV). 2nd column: spectrum (1/eV).',
     1 /1X,'#  3rd and 4th columns: simul. pdf limits (3SD, 1/eV).',/)
        PTOT=0.0D0
        WRITE(9,'(1X,1P,4E14.6)') ES(1),PTOT,PTOT,PTOT
        DO KE=1,NSEB
          PTOT=PTOT+PTS(KE)
        ENDDO
        IF(PTOT.GT.1.0D-35) THEN
          DO KE=1,NSEB
            YAV=SHIST(KE)*DF
            YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV)*DF))
            EINTL=ES(KE+1)-ES(KE)
            IF(EINTL.GT.1.0D-15) THEN
              FACT=1.0D0/EINTL
            ELSE
              FACT=1.0D15
            ENDIF
            WRITE(9,'(1X,1P,4E14.6)') ES(KE),PTS(KE)*FACT/PTOT,
     1        (YAV-YERR)*FACT,(YAV+YERR)*FACT
            WRITE(9,'(1X,1P,4E14.6)') ES(KE+1),PTS(KE)*FACT/PTOT,
     1        (YAV-YERR)*FACT,(YAV+YERR)*FACT
          ENDDO
          PTOT=0.0D0
          WRITE(9,'(1X,1P,4E14.6)') ES(NSEB+1),PTOT,PTOT,PTOT
        ELSE
          DO KE=1,NSEB
            EINTL=DBLE(KE-1)*DSHE
            YAV=SHIST(KE)*DF
            YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV)*DF))
            WRITE(9,'(1X,1P,4E14.6)') EINTL,YAV*RDSHE,
     1        (YAV-YERR)*RDSHE,(YAV+YERR)*RDSHE
            WRITE(9,'(1X,1P,4E14.6)') EINTL+DSHE,YAV*RDSHE,
     1        (YAV-YERR)*RDSHE,(YAV+YERR)*RDSHE
          ENDDO
        ENDIF
        CLOSE(9)
      ENDIF
C
C  ************  Energy distributions of emerging particles.
C
C  ****  Upbound electrons.
      OPEN(9,FILE='data/energy-el-up.dat')
      WRITE(9,9110)
 9110 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of upbound electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE
        YERR=3.0D0*SQRT(ABS(PDE2(1,1,K)-PDE(1,1,K)**2*DF))
        YAV=PDE(1,1,K)*DF*RBSE
        YERR=YERR*DF*RBSE
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Downbound electrons.
      OPEN(9,FILE='data/energy-el-down.dat')
      WRITE(9,9120)
 9120 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of downbound electrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE
        YERR=3.0D0*SQRT(ABS(PDE2(1,2,K)-PDE(1,2,K)**2*DF))
        YAV=PDE(1,2,K)*DF*RBSE
        YERR=YERR*DF*RBSE
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Upbound photons.
      OPEN(9,FILE='data/energy-ph-up.dat')
      WRITE(9,9130)
 9130 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of upbound photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE
        YERR=3.0D0*SQRT(ABS(PDE2(2,1,K)-PDE(2,1,K)**2*DF))
        YAV=PDE(2,1,K)*DF*RBSE
        YERR=YERR*DF*RBSE
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Downbound photons.
      OPEN(9,FILE='data/energy-ph-down.dat')
      WRITE(9,9140)
 9140 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of downbound photons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE
        YERR=3.0D0*SQRT(ABS(PDE2(2,2,K)-PDE(2,2,K)**2*DF))
        YAV=PDE(2,2,K)*DF*RBSE
        YERR=YERR*DF*RBSE
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Upbound positrons.
      OPEN(9,FILE='data/energy-po-up.dat')
      WRITE(9,9150)
 9150 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of upbound positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE
        YERR=3.0D0*SQRT(ABS(PDE2(3,1,K)-PDE(3,1,K)**2*DF))
        YAV=PDE(3,1,K)*DF*RBSE
        YERR=YERR*DF*RBSE
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C  ****  Downbound positrons.
      OPEN(9,FILE='data/energy-po-down.dat')
      WRITE(9,9160)
 9160 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Energy distribution of downbound positrons.',
     1 /1X,'#  1st column: E (eV).',
     1 /1X,'#  2nd and 3rd columns: probability density and STU',
     1         ' (1/(eV*particle)).',/)
      DO K=1,NBE
        XX=EMIN+(K-0.5D0)*BSE
        YERR=3.0D0*SQRT(ABS(PDE2(3,2,K)-PDE(3,2,K)**2*DF))
        YAV=PDE(3,2,K)*DF*RBSE
        YERR=YERR*DF*RBSE
        WRITE(9,'(1X,1P,3E14.6)')
     1     XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
      ENDDO
      CLOSE(9)
C
C  ************  Angular distributions of emerging particles.
C
      IF(NBPH.EQ.1) THEN
        OPEN(9,FILE='data/polar-angle-el.dat')
      ELSE
        OPEN(9,FILE='data/angle-el.dat')
      ENDIF
      WRITE(9,9210)
 9210 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Angular distribution of emerging electrons.',
     1 /1X,'#  1st and 2nd columns: THETA and PHI (deg).',
     1 /1X,'#  3rd and 4th columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*(BSPH*DE2RA)
        DO L=1,NBPH
          YY=(L-0.5D0)*BSPH
          YERR=3.0D0*SQRT(ABS(PDA2(1,K,L)-PDA(1,K,L)**2*DF))
          YAV=PDA(1,K,L)*DF/DSANG
          YERR=YERR*DF/DSANG
          WRITE(9,'(1X,1P,6E14.6)')
     1       XX,YY,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        IF(NBPH.GT.1) WRITE(9,*) '   '
      ENDDO
      CLOSE(9)
C
      IF(NBPH.EQ.1) THEN
        OPEN(9,FILE='data/polar-angle-ph.dat')
      ELSE
        OPEN(9,FILE='data/angle-ph.dat')
      ENDIF
      WRITE(9,9220)
 9220 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Angular distribution of emerging photons.',
     1 /1X,'#  1st and 2nd columns: THETA and PHI (deg).',
     1 /1X,'#  3rd and 4th columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*(BSPH*DE2RA)
        DO L=1,NBPH
          YY=(L-0.5D0)*BSPH
          YERR=3.0D0*SQRT(ABS(PDA2(2,K,L)-PDA(2,K,L)**2*DF))
          YAV=PDA(2,K,L)*DF/DSANG
          YERR=YERR*DF/DSANG
          WRITE(9,'(1X,1P,6E14.6)')
     1       XX,YY,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        IF(NBPH.GT.1) WRITE(9,*) '   '
      ENDDO
      CLOSE(9)
C
      IF(NBPH.EQ.1) THEN
        OPEN(9,FILE='data/polar-angle-po.dat')
      ELSE
        OPEN(9,FILE='data/angle-po.dat')
      ENDIF
      WRITE(9,9230)
 9230 FORMAT(
     1  1X,'#  Results from PENMAIN.',
     1 /1X,'#  Angular distribution of emerging positrons.',
     1 /1X,'#  1st and 2nd columns: THETA and PHI (deg).',
     1 /1X,'#  3rd and 4th columns: probability density and STU',
     1         ' (1/sr)',/)
      DO K=1,NBTH
        XX=(K-0.5D0)*BSTH
        XXR=(K-1.0D0)*BSTH*DE2RA
        DSANG=(COS(XXR)-COS(XXR+BSTH*DE2RA))*(BSPH*DE2RA)
        DO L=1,NBPH
          YY=(L-0.5D0)*BSPH
          YERR=3.0D0*SQRT(ABS(PDA2(3,K,L)-PDA(3,K,L)**2*DF))
          YAV=PDA(3,K,L)*DF/DSANG
          YERR=YERR*DF/DSANG
          WRITE(9,'(1X,1P,6E14.6)')
     1       XX,YY,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        IF(NBPH.GT.1) WRITE(9,*) '   '
      ENDDO
      CLOSE(9)
C
C  ************  Spectra from impact detectors.
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
          OPEN(9,FILE=SPCDIO(ID))
          WRITE(9,9310) ID
 9310 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from impact detector #',I3,
     1 /1X,'#  Energy spectra of incident particles.',
     1 /1X,'#  1st column: particle energy (eV).',
     1 /1X,'#  2nd column: probability density (1/(eV*particle)).',
     1 /1X,'#  3rd column: statistical uncertainty, STU (3 sigma).',
     1 /1X,'#  4th-5th columns: electron spectrum and STU (3 sigma).',
     1 /1X,'#  6th-7th columns: photon spectrum and STU (3 sigma).',
     1 /1X,'#  8th-9th columns: positron spectrum and STU (3 sigma).',/)
          DO J=1,NDICH(ID)
            IF(LDILOG(ID)) THEN
              XLOW=EXP(EDILL(ID)+(J-1)*BDIEL(ID))
              XUPP=EXP(EDILL(ID)+J*BDIEL(ID))
              XX=0.5D0*(XUPP+XLOW)
              BINS=XUPP-XLOW
            ELSE
              XX=EDIL(ID)+(J-0.5D0)*BDIE(ID)
              BINS=BDIE(ID)
            ENDIF
            YERR=3.0D0*SQRT(ABS(DIT2(ID,J)-DIT(ID,J)**2*DF))
            YAV=MAX(DIT(ID,J)*DF/BINS,1.0D-35)
            YERR=MAX(YERR*DF/BINS,1.0D-35)
C
            YERR1=3.0D0*SQRT(ABS(DIP2(ID,J,1)-DIP(ID,J,1)**2*DF))
            YAV1=MAX(DIP(ID,J,1)*DF/BINS,1.0D-35)
            YERR1=MAX(YERR1*DF/BINS,1.0D-35)
C
            YERR2=3.0D0*SQRT(ABS(DIP2(ID,J,2)-DIP(ID,J,2)**2*DF))
            YAV2=MAX(DIP(ID,J,2)*DF/BINS,1.0D-35)
            YERR2=MAX(YERR2*DF/BINS,1.0D-35)
C
            YERR3=3.0D0*SQRT(ABS(DIP2(ID,J,3)-DIP(ID,J,3)**2*DF))
            YAV3=MAX(DIP(ID,J,3)*DF/BINS,1.0D-35)
            YERR3=MAX(YERR3*DF/BINS,1.0D-35)
            WRITE(9,'(1X,1P,E14.6,4(E14.6,E10.2))')
     1        XX,YAV,YERR,YAV1,YERR1,YAV2,YERR2,YAV3,YERR3
          ENDDO
          CLOSE(9)
        ENDDO
      ENDIF
C
C  ************  Distributions of fluence with respect to energy
C                (integrated over the volume of the impact detectors).
C
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
          IF(IDCUT(ID).EQ.2) THEN
          OPEN(9,FILE=SPCFSO(ID))
          WRITE(9,9330) ID
 9330 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from impact detector #',I3,
     1 /1X,'#  Fluences integrated over the volume of the detector ',
     1   '(in cm/(eV*particle)).',
     1 /1X,'#  1st column: particle energy (eV).',
     1 /1X,'#  2nd column: energy distribution of total fluence.',
     1 /1X,'#  3rd column: statistical uncertainty, STU (3 sigma).',
     1 /1X,'#  4th-5th columns: electron fluence and STU (3 sigma).',
     1 /1X,'#  6th-7th columns: photon fluence and STU (3 sigma).',
     1 /1X,'#  8th-9th columns: positron fluence and STU (3 sigma).',/)
          DO J=1,NDICH(ID)
            IF(LDILOG(ID)) THEN
              XLOW=EXP(EDILL(ID)+(J-1)*BDIEL(ID))
              XUPP=EXP(EDILL(ID)+J*BDIEL(ID))
              XX=0.5D0*(XUPP+XLOW)
              BINS=XUPP-XLOW
            ELSE
              XX=EDIL(ID)+(J-0.5D0)*BDIE(ID)
              BINS=BDIE(ID)
            ENDIF
            YERR=3.0D0*SQRT(ABS(FST2(ID,J)-FST(ID,J)**2*DF))
            YAV=MAX(FST(ID,J)*DF/BINS,1.0D-35)
            YERR=MAX(YERR*DF/BINS,1.0D-35)
C
            YERR1=3.0D0*SQRT(ABS(FSP2(ID,J,1)-FSP(ID,J,1)**2*DF))
            YAV1=MAX(FSP(ID,J,1)*DF/BINS,1.0D-35)
            YERR1=MAX(YERR1*DF/BINS,1.0D-35)
C
            YERR2=3.0D0*SQRT(ABS(FSP2(ID,J,2)-FSP(ID,J,2)**2*DF))
            YAV2=MAX(FSP(ID,J,2)*DF/BINS,1.0D-35)
            YERR2=MAX(YERR2*DF/BINS,1.0D-35)
C
            YERR3=3.0D0*SQRT(ABS(FSP2(ID,J,3)-FSP(ID,J,3)**2*DF))
            YAV3=MAX(FSP(ID,J,3)*DF/BINS,1.0D-35)
            YERR3=MAX(YERR3*DF/BINS,1.0D-35)
            WRITE(9,'(1X,1P,E14.6,4(E14.6,E10.2))')
     1        XX,YAV,YERR,YAV1,YERR1,YAV2,YERR2,YAV3,YERR3
          ENDDO
          CLOSE(9)
          ENDIF
        ENDDO
      ENDIF
C
C  ************  Spectra from energy-deposition detectors.
C
      IF(NDEDEF.GT.0) THEN
        DO ID=1,NDEDEF
          OPEN(9,FILE=SPCDEO(ID))
          WRITE(9,9320) ID
 9320 FORMAT(
     1  1X,'#  Results from PENMAIN. Output from energy-deposition',
     1  ' detector #',I3,
     1 /1X,'#  Deposited energy spectrum.',
     1 /1X,'#  WARNING: May be strongly biased if interaction ',
     1  'forcing is used!',
     1 /1X,'#  1st column: deposited energy (eV).',
     1 /1X,'#  2nd column: probability density (1/(eV*particle)).',
     1 /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
          DO J=1,NDECH(ID)
            IF(LDELOG(ID)) THEN
              XLOW=EXP(EDELL(ID)+(J-1)*BDEEL(ID))
              XUPP=EXP(EDELL(ID)+J*BDEEL(ID))
              XX=0.5D0*(XUPP+XLOW)
              BINS=XUPP-XLOW
            ELSE
              XX=EDEL(ID)+(J-0.5D0)*BDEE(ID)
              BINS=BDEE(ID)
            ENDIF
            YAV=DET(ID,J)*DF
            YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV))*DF)
            YAV=YAV/BINS
            YERR=YERR/BINS
            WRITE(9,'(1X,1P,3E14.6)')
     1        XX,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
          ENDDO
          CLOSE(9)
        ENDDO
      ENDIF
C
C  ************  Dose distributions.
C
      IF(LDOSEM) THEN
C
C  ****  Depth-dose distribution.
C
        OPEN(9,FILE='data/depth-dose.dat')
        WRITE(9,9410)
 9410   FORMAT(
     1     1X,'#  Results from PENMAIN. Depth-dose distribution.',
     1    /1X,'#  (integrated over X and Y within the volume of the ',
     1      'dose-map box).',
     1    /1X,'#  1st column: z coordinate (cm).',
     1    /1X,'#  2nd column: depth-dose (eV/(g/cm**2)).',
     1    /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
        DO I3=1,NDB(3)
          ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          YAV=DDOSE(I3)
          YAV2=DDOSE2(I3)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF/BDOSE(3)
          YERR=YERR*DF/BDOSE(3)
          WRITE(9,'(1X,1P,3E14.6)')
     1      ZZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
C
        VOXEL=BDOSE(1)*BDOSE(2)*BDOSE(3)
        DMAX=0.0D0
        I1M=1
        I2M=1
        I3M=1
        DO I1=1,NDB(1)
          DO I2=1,NDB(2)
            DO I3=1,NDB(3)
              IF(DOSE(I1,I2,I3).GT.DMAX) THEN
                I1M=I1
                I2M=I2
                I3M=I3
                DMAX=DOSE(I1,I2,I3)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
        QAV=DOSE(I1M,I2M,I3M)
        QAV2=DOSE2(I1M,I2M,I3M)
        QER=3.0D0*SQRT(ABS(QAV2-QAV**2*DF))
        QAV=QAV*DF/VOXEL
        QER=QER*DF/VOXEL
        IF(QER.GT.1.0D-10*ABS(QAV)) THEN
          EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
        ELSE
          EFFIC=0.0D0
        ENDIF
        WRITE(27,3020) QAV,QER,EFFIC
 3020   FORMAT(/6X,'Maximum dose ... ',1P,E13.6,' +-',E8.1,' eV/g',2X,
     1    '(effic. =',E9.2,')')
C
        OPEN(9,FILE='data/3d-dose.dat')
        WRITE(9,9420)
 9420   FORMAT(1X,'#  Results from PENMAIN. 3D dose distribution.')
        WRITE(9,9421) DXL(1),DXU(1)
 9421   FORMAT(1X,'#  Dose-map box:  XL = ',1P,E13.6,
     1    ' cm,  XU = ',E13.6,' cm')
        WRITE(9,9422) DXL(2),DXU(2)
 9422   FORMAT(1X,'#',17X,'YL = ',1P,E13.6,' cm,  YU = ',E13.6,' cm')
        WRITE(9,9423) DXL(3),DXU(3)
 9423   FORMAT(1X,'#',17X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,' cm')
        WRITE(9,9424) NDB(1),NDB(2),NDB(3)
 9424   FORMAT(1X,'#  Numbers of bins:     NBX =',I4,', NBY =',I4,
     1        ', NBZ =',I4,/1X,'#')
        WRITE(9,9425)
 9425   FORMAT(1X,'#  columns 1 to 3: coordinates X,Y,Z of the bin',
     1    ' centres',/1X,'#  4th column: dose (eV/g).',
     1    /1X,'#  5th column: statistical uncertainty (3 sigma).',
     1    /1X,'#  columns 6 to 8: bin indices IX,IY,IZ.',/)
        DO I3=1,NDB(3)
          ZZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          DO I1=1,NDB(1)
            XX=DXL(1)+(I1-0.5D0)*BDOSE(1)
            DO I2=1,NDB(2)
              YY=DXL(2)+(I2-0.5D0)*BDOSE(2)
              YAV=DOSE(I1,I2,I3)
              YAV2=DOSE2(I1,I2,I3)
              YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
              YAV=MAX(YAV*DF/VOXEL,1.0D-35)
              YERR=MAX(YERR*DF/VOXEL,1.0D-35)
              WRITE(9,9426) XX,YY,ZZ,YAV,YERR,I1,I2,I3
            ENDDO
            WRITE(9,*) '   '
          ENDDO
          WRITE(9,*) '   '
        ENDDO
 9426   FORMAT(1P,3E11.3,1P,E13.5,E9.1,3I4)
        CLOSE(9)
C
C  ****  Dose distributions at the central axes.
C
        I1C=(NDB(1)/2)+1
        I2C=(NDB(2)/2)+1
        I3C=(NDB(3)/2)+1
C
        OPEN(9,FILE='data/x-dose.dat')
          WRITE(9,9440)
 9440   FORMAT(
     1     1X,'#  Results from PENMAIN.',
     1    /1X,'#  Dose distribution along the central X axis',
     1    /1X,'#  1st column: x (cm).',
     1    /1X,'#  2nd column: dose (eV/g).',
     1    /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
        DO I1=1,NDB(1)
          XYZ=DXL(1)+(I1-0.5D0)*BDOSE(1)
          YAV=DOSE(I1,I2C,I3C)
          YAV2=DOSE2(I1,I2C,I3C)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF/VOXEL
          YERR=YERR*DF/VOXEL
          WRITE(9,'(1X,1P,3E14.6)')
     1      XYZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
C
        OPEN(9,FILE='data/y-dose.dat')
          WRITE(9,9450)
 9450   FORMAT(
     1     1X,'#  Results from PENMAIN.',
     1    /1X,'#  Dose distribution along the central Y axis',
     1    /1X,'#  1st column: y (cm).',
     1    /1X,'#  2nd column: dose (eV/g).',
     1    /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
        DO I2=1,NDB(2)
          XYZ=DXL(2)+(I2-0.5D0)*BDOSE(2)
          YAV=DOSE(I1C,I2,I3C)
          YAV2=DOSE2(I1C,I2,I3C)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF/VOXEL
          YERR=YERR*DF/VOXEL
          WRITE(9,'(1X,1P,3E14.6)')
     1      XYZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
C
        OPEN(9,FILE='data/z-dose.dat')
          WRITE(9,9460)
 9460   FORMAT(
     1     1X,'#  Results from PENMAIN.',
     1    /1X,'#  Dose distribution along the central Z axis',
     1    /1X,'#  1st column: z (cm).',
     1    /1X,'#  2nd column: dose (eV/g).',
     1    /1X,'#  3rd column: statistical uncertainty (3 sigma).',/)
        DO I3=1,NDB(3)
          XYZ=DXL(3)+(I3-0.5D0)*BDOSE(3)
          YAV=DOSE(I1C,I2C,I3)
          YAV2=DOSE2(I1C,I2C,I3)
          YERR=3.0D0*SQRT(ABS(YAV2-YAV**2*DF))
          YAV=YAV*DF/VOXEL
          YERR=YERR*DF/VOXEL
          WRITE(9,'(1X,1P,3E14.6)')
     1      XYZ,MAX(YAV,1.0D-35),MAX(YERR,1.0D-35)
        ENDDO
        CLOSE(9)
      ENDIF
C
      WRITE(27,3030) ISEED1,ISEED2
 3030 FORMAT(/3X,'Last random seeds = ',I10,' , ',I10)
      WRITE(27,'(/3X,72(''-''))')
      CLOSE(27)
      IF(ICLOSE.LT.0) RETURN
C
C  ************  Write final average results in output file.
C
      WRITE(26,'(/3X,72(''-''))')
      WRITE(26,3000)
      WRITE(26,1001) DATE23
      WRITE(26,'(/3X,A65)') TITLE
C
      WRITE(26,3001) TSIM
      WRITE(26,3002) TAVS
      WRITE(26,3003) SHN
C
      IF(KPARP.EQ.1) WRITE(26,1110)
      IF(KPARP.EQ.2) WRITE(26,1111)
      IF(KPARP.EQ.3) WRITE(26,1112)
      IF(KPARP.EQ.0) WRITE(26,1113)
C
      WRITE(26,3004) PRIM(1)
      WRITE(26,3005) PRIM(2)
      WRITE(26,3006) PRIM(3)
C
      FNT=1.0D0/SHN
      FT=(PRIM(1)+SEC(KPARP,1))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(1)-PRIM(1)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,1)-SEC(KPARP,1)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(26,3007) FT,ERR
      FB=(PRIM(2)+SEC(KPARP,2))*FNT
      ERR1=3.0D0*FNT*SQRT(ABS(PRIM2(2)-PRIM(2)**2*FNT))
      ERR2=3.0D0*FNT*SQRT(ABS(SEC2(KPARP,2)-SEC(KPARP,2)**2*FNT))
      ERR=ERR1+ERR2
      WRITE(26,3008) FB,ERR
      FA=PRIM(3)*FNT
      ERR=3.0D0*FNT*SQRT(ABS(PRIM2(3)-PRIM(3)**2*FNT))
      WRITE(26,3009) FA,ERR
C
      WRITE(26,3010)
     1 WSEC(1,1),WSEC(2,1),WSEC(3,1),WSEC2(1,1),WSEC2(2,1),WSEC2(3,1),
     1 WSEC(1,2),WSEC(2,2),WSEC(3,2),WSEC2(1,2),WSEC2(2,2),WSEC2(3,2),
     1 WSEC(1,3),WSEC(2,3),WSEC(3,3),WSEC2(1,3),WSEC2(2,3),WSEC2(3,3)
C
C  ****  Average energies deposited in bodies..
C
      WRITE(26,3011)
      DO KB=1,NBODY
        IF(MATER(KB).NE.0) THEN
          QER=3.0D0*DF*SQRT(ABS(TDEBO2(KB)-TDEBO(KB)**2*DF))
          QAV=TDEBO(KB)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(26,3012) KB,QAV,QER,EFFIC
        ENDIF
      ENDDO
C
C  ****  Average energies 'collected' by impact detectors.
C
      IF(NDIDEF.GT.0) THEN
        WRITE(26,3013)
        DO KD=1,NDIDEF
          QER=3.0D0*DF*SQRT(ABS(TDID2(KD)-TDID(KD)**2*DF))
          QAV=TDID(KD)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(26,3014) KD,QAV,QER,EFFIC
        ENDDO
      ENDIF
C
C  ****  Average deposited energies in energy-deposition detectors.
C
      IF(NDEDEF.GT.0) THEN
        WRITE(26,3015)
        DO KD=1,NDEDEF
          QER=3.0D0*DF*SQRT(ABS(TDED2(KD)-TDED(KD)**2*DF))
          QAV=TDED(KD)*DF
          IF(QER.GT.1.0D-10*ABS(QAV)) THEN
            EFFIC=QAV**2/((QER/3.0D0)**2*TSIM)
          ELSE
            EFFIC=0.0D0
          ENDIF
          WRITE(26,3014) KD,QAV,QER,EFFIC
        ENDDO
      ENDIF
C
C  ****  Detector efficiencies.
C
      IF(NDIDEF.GT.0) THEN
        WRITE(26,3016)
        DO KD=1,NDIDEF
          SUM=0.0D0
          ESUM=0.0D0
          DO J=1,NDICH(KD)
            YERR=3.0D0*SQRT(ABS(DIT2(KD,J)-DIT(KD,J)**2*DF))
            YAV=MAX(DIT(KD,J)*DF,1.0D-35)
            YERR=MAX(YERR*DF,1.0D-35)
            SUM=SUM+YAV
            ESUM=ESUM+YERR
          ENDDO
          WRITE(26,3017) KD,SUM,ESUM
        ENDDO
      ENDIF
C
      IF(NDEDEF.GT.0) THEN
        WRITE(26,3018)
        DO KD=1,NDEDEF
          SUM=0.0D0
          ESUM=0.0D0
          DO J=1,NDECH(KD)
            YAV=DET(KD,J)*DF
            YERR=3.0D0*SQRT(ABS(YAV*(1.0D0-YAV))*DF)
            SUM=SUM+YAV
            ESUM=ESUM+YERR
          ENDDO
          WRITE(26,3017) KD,SUM,ESUM
        ENDDO
      ENDIF
C
      WRITE(26,3030) ISEED1,ISEED2
      WRITE(26,'(/3X,''***  END  ***'')')
      WRITE(26,'(/3X,72(''-''))')
      CLOSE(26)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE PMRDR
C  *********************************************************************
      SUBROUTINE PMRDR
C
C  Reads the input file and initializes PENELOPE and PENGEOM.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LIT
      CHARACTER*20 PMFILE,PFILE,PFILER
      CHARACTER*120 BUFFER,BUF2
C
      CHARACTER*6 KWORD,
     1  KWTITL,KWKPAR,KWSENE,KWSPEC,  KWSPOL,KWSPOS,KWSBOX,KWSBOD,
     1  KWSCON,KWSREC,KWPSFN,KWPSPL,  KWRRSP,KWEMAX,KWMATF,KWSIMP,
     1  KWGEOM,KWSMAX,KWEABS,KWIFOR,  KWNBE,KWNBAN,KWIDET,KWIPSF,
     1  KWISPC,KWIFLN,KWIBOD,KWIPAR,  KWEDET,KWESPC,KWEBOD,KGRDXX,
     1  KGRDYY,KGRDZZ,KGRDBN,KWRESU,  KWDUMP,KWDMPP,KWRSEE,KWNSIM,
     1  KWTIME,KWCOMM,KZMAXD,KRADMX,  KHBKFN,KSIMTP
      PARAMETER (
     1  KWTITL='TITLE ',KWKPAR='SKPAR ',KWSENE='SENERG',KWSPEC='SPECTR',
     1  KWSPOL='SGPOL ',KWSPOS='SPOSIT',KWSBOX='SBOX  ',KWSBOD='SBODY ',
     1  KWSCON='SCONE ',KWSREC='SPYRAM',KWPSFN='IPSFN ',KWPSPL='IPSPLI',
     1  KWRRSP='WGTWIN',KWEMAX='EPMAX ',KWMATF='MFNAME',KWSIMP='MSIMPA',
     1  KWGEOM='GEOMFN',KWSMAX='DSMAX ',KWEABS='EABSB ',KWIFOR='IFORCE',
     1  KWNBE ='NBE   ',KWNBAN='NBANGL',KWIDET='IMPDET',KWIPSF='IDPSF ',
     1  KWISPC='IDSPC ',KWIFLN='IDFLNC',KWIBOD='IDBODY',KWIPAR='IDKPAR',
     1  KWEDET='ENDETC',KWESPC='EDSPC ',KWEBOD='EDBODY',KGRDXX='GRIDX ',
     1  KGRDYY='GRIDY ',KGRDZZ='GRIDZ ',KGRDBN='GRIDBN',KWRESU='RESUME',
     1  KWDUMP='DUMPTO',KWDMPP='DUMPP ',KWRSEE='RSEED ',KWNSIM='NSIMSH',
     1  KWTIME='TIME  ',KWCOMM='      ',KZMAXD='ZMAX  ',KRADMX='RADMAX',
     1  KHBKFN='HBOOKF',KSIMTP='SIMTYP')
C
      INCLUDE 'pmcomms.f'
      DIMENSION PARINP(20)
      DIMENSION PMFILE(MAXMAT)
C
      OPEN(26,FILE='data/penmain.dat')  ! Global output/message file.
C
      DO I=1,3
        PRIM(I)=0.0D0
        PRIM2(I)=0.0D0
        DO J=1,3
          SEC(I,J)=0.0D0
          SEC2(I,J)=0.0D0
        ENDDO
      ENDDO
C
      DO I=1,NSEM
        SHIST(I)=0.0D0
      ENDDO
C
      DO I=1,NB
        TDEBO(I)=0.0D0
        TDEBO2(I)=0.0D0
        EABSB(1,I)=50.0D0
        EABSB(2,I)=50.0D0
        EABSB(3,I)=50.0D0
      ENDDO
C
      DO I=1,3
        DO J=1,2
          DO K=1,NBEM
            PDE(I,J,K)=0.0D0
            PDE2(I,J,K)=0.0D0
            PDEP(I,J,K)=0.0D0
            LPDE(I,J,K)=0
          ENDDO
        ENDDO
      ENDDO
C
      DO I=1,3
        DO J=1,NBTHM
          DO K=1,NBPHM
            PDA(I,J,K)=0.0D0
            PDA2(I,J,K)=0.0D0
            PDAP(I,J,K)=0.0D0
            LPDA(I,J,K)=0
          ENDDO
        ENDDO
      ENDDO
C
      DO I=1,NIDM
        LDILOG(I)=.FALSE.
        DO J=1,NIDCM
          DIT(I,J)=0.0D0
          DIT2(I,J)=0.0D0
          DITP(I,J)=0.0D0
          LDIT(I,J)=0
          DO K=1,3
            DIP(I,J,K)=0.0D0
            DIP2(I,J,K)=0.0D0
            DIPP(I,J,K)=0.0D0
            LDIP(I,J,K)=0
          ENDDO
        ENDDO
      ENDDO
C
      DO I=1,NIDM
        DO J=1,NIDCM
          FST(I,J)=0.0D0
          FST2(I,J)=0.0D0
          FSTP(I,J)=0.0D0
          LFST(I,J)=0
          DO K=1,3
            FSP(I,J,K)=0.0D0
            FSP2(I,J,K)=0.0D0
            FSPP(I,J,K)=0.0D0
            LFSP(I,J,K)=0
          ENDDO
        ENDDO
      ENDDO
C
      DO I=1,NEDM
        LDELOG(I)=.FALSE.
        DO J=1,NEDCM
          DET(I,J)=0.0D0
        ENDDO
      ENDDO
C
      DO K=1,NDZM
        DO J=1,NDYM
          DO I=1,NDXM
            DOSE(I,J,K)=0.0D0
            DOSE2(I,J,K)=0.0D0
            DOSEP(I,J,K)=0.0D0
            LDOSE(I,J,K)=0
          ENDDO
        ENDDO
        DDOSE(K)=0.0D0
        DDOSE2(K)=0.0D0
        DDOSEP(K)=0.0D0
        LDDOSE(K)=0
      ENDDO
C
      CTHL=0.0D0
      DCTH=0.0D0
      PHIL=0.0D0
      DPHI=0.0D0
      DO KB=1,NB
        IXSBOD(KB)=0
      ENDDO
C
C  ****  Time counter initiation.
C
      CALL TIME0
C
C  ------------------------  Read input data file.
C
      WRITE(26,1000)
 1000 FORMAT(//3X,61('*'),/3X,'**   Program PENMAIN. ',
     1 ' Input data and run-time messages.   **',/3X,61('*'))
C
      CALL PDATET(DATE23)
      WRITE(26,1001) DATE23
 1001 FORMAT(/3X,'Date and time: ',A23)
C
C  ****  Title.
C
      READ(5,'(A6,1X,A65)') KWORD,TITLE
      IF(KWORD.EQ.KWTITL) THEN
        WRITE(26,'(/3X,A65)') TITLE
      ELSE
        WRITE(26,*) 'The input file must begin with the TITLE line.'
        STOP 'The input file must begin with the TITLE line.'
      ENDIF
C
      KPARP=1
      NPSF=0
      NPSN=0
      NSEB=1
      NSPLIT=1
      RLREAD=0.0D0
      RWGMIN=1.0D35
      KBSMAX=0
      ISEC=0
      LDUMP=.FALSE.
      LPSF=.FALSE.
      LSPEC=.FALSE.
      LGPOL=.FALSE.
      LEXSRC=.FALSE.
      LEXBD=.FALSE.
      LSCONE=.TRUE.
      JOBEND=0
C
C  ************  Source description.
C
   11 CONTINUE
      READ(5,'(A6,1X,A120)') KWORD,BUFFER
      IF(KWORD.EQ.KWCOMM) GO TO 11
      IF(KWORD.EQ.KWPSFN) GO TO 21
C
      WRITE(26,1100)
 1100 FORMAT(/3X,72('-'),/3X,'>>>>>>  Source description.')
C
      IF(KWORD.EQ.KWKPAR) THEN
        READ(BUFFER,*) KPARP
   12   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 12
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ENDIF
      IF(KPARP.LT.0.OR.KPARP.GT.3) THEN
        WRITE(26,*) 'KPARP =',KPARP
        WRITE(26,*) 'Incorrect particle type.'
        STOP 'Incorrect particle type.'
      ENDIF
      IF(KPARP.EQ.1) WRITE(26,1110)
 1110 FORMAT(/3X,'Primary particles: electrons')
      IF(KPARP.EQ.2) WRITE(26,1111)
 1111 FORMAT(/3X,'Primary particles: photons')
      IF(KPARP.EQ.3) WRITE(26,1112)
 1112 FORMAT(/3X,'Primary particles: positrons')
      IF(KPARP.EQ.0) WRITE(26,1113)
 1113 FORMAT(/3X,'Primary particles: set by the user subroutine SOURCE')
 
      IF(KWORD.EQ.KSIMTP) THEN
        READ(BUFFER,*) NSIMTYPE
 1215   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 1215
        WRITE(26,1114)
 1114   FORMAT(/3X,'Source generation by option : ',I2)
      ENDIF
C
C  ****  Monoenergetic source.
C
      IF(KWORD.EQ.KWSENE) THEN
        READ(BUFFER,*) E0
        WRITE(26,1120) E0
 1120   FORMAT(3X,'Initial energy = ',1P,E13.6,' eV')
   13   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 13
        IF(KWORD.EQ.KWPSFN) GO TO 21
C
C  ****  Continuous energy spectrum.
C
      ELSE IF(KWORD.EQ.KWSPEC) THEN
        LSPEC=.TRUE.
        NSEB=0
   14   CONTINUE
        NSEB=NSEB+1
        IF(NSEB.GT.NSEM) THEN
          WRITE(26,*) 'Source energy spectrum.'
          WRITE(26,*) 'The number of energy bins is too large.'
          STOP 'The number of energy bins is too large.'
        ENDIF
        READ(BUFFER,*) ES(NSEB),PTS(NSEB)
        PTS(NSEB)=MAX(PTS(NSEB),0.0D0)
   15   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 15
        IF(KWORD.EQ.KWSPEC) GO TO 14
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        E0=1.0D9
        IF(KPARP.NE.0) WRITE(26,1120) E0
      ENDIF
      IF(LSPEC) THEN
        IF(NSEB.GT.1) THEN
          CALL SORT2(ES,PTS,NSEB)
          WRITE(26,1121)
 1121     FORMAT(/3X,'Spectrum:',7X,'I',4X,'E_low(eV)',4x,'E_high(eV)',
     1      5X,'P_sum(E)',/16X,45('-'))
          DO I=1,NSEB-1
            WRITE(26,'(16X,I4,1P,3E14.6)') I,ES(I),ES(I+1),PTS(I)
          ENDDO
          WRITE(26,*) '  '
          E0=ES(NSEB)
          NSEB=NSEB-1
          CALL IRND0(PTS,FS,IAS,NSEB)
        ELSE
          WRITE(26,*) 'The source energy spectrum is not defined.'
          STOP 'The source energy spectrum is not defined.'
        ENDIF
      ENDIF
      IF(E0.LT.50.0D0) THEN
        WRITE(26,*) 'The initial energy E0 is too small.'
        STOP 'The initial energy E0 is too small.'
      ENDIF
      EPMAX=E0
C  ----  Positrons eventually give annihilation gamma rays. The maximum
C        energy of annihilation photons is .lt. 1.21*(E0+me*c**2).
      IF(KPARP.EQ.3) EPMAX=1.21D0*(E0+5.12D5)
C
C  ****  Photon polarisation effects (only for primary photons).
C
      IF(KWORD.EQ.KWSPOL) THEN
        READ(BUFFER,*) SP10,SP20,SP30
        LGPOL=.TRUE.
   20   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 20
        IF(KWORD.EQ.KWPSFN) GO TO 21
        WRITE(26,1131) SP10,SP20,SP30
 1131   FORMAT(/3X,'Polarised primary photons. Stokes Parameters:',
     1    /6X,'P1 = ',E13.6,' (linear polarisation at 45 deg azimuth)',
     1    /6X,'P2 = ',E13.6,' (circular polarisation)',
     1    /6X,'P3 = ',E13.6,' (linear polarisation at zero azimuth)')
      ENDIF
C
C  ****  Position of the point source.
C
      IF(KWORD.EQ.KWSPOS) THEN
        READ(BUFFER,*) SX0,SY0,SZ0
   16   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 16
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        SX0=0.0D0
        SY0=0.0D0
        SZ0=0.0D0
      ENDIF
      WRITE(26,1132) SX0,SY0,SZ0
 1132 FORMAT(/3X,'Coordinates of centre:     SX0 = ',1P,E13.6,
     1  ' cm',/30X,'SY0 = ',E13.6,' cm',/30X,'SZ0 = ',E13.6,' cm')
C
C  ****  Extended source.
C
      IF(KWORD.EQ.KWSBOX) THEN
        LEXSRC=.TRUE.
        READ(BUFFER,*) SSX,SSY,SSZ
        SSX=ABS(SSX)
        SSY=ABS(SSY)
        SSZ=ABS(SSZ)
   17   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 17
        IF(KWORD.EQ.KWPSFN) GO TO 21
      ELSE
        LEXSRC=.FALSE.
        SSX=0.0D0
        SSY=0.0D0
        SSZ=0.0D0
      ENDIF
      IF(LEXSRC) WRITE(26,1135) SSX,SSY,SSZ
 1135 FORMAT(3X,'Source size:',15X,'SSX = ',1P,E13.6,' cm',
     1  /30X,'SSY = ',E13.6,' cm',/30X,'SSZ = ',E13.6,' cm')
C
C  ****  Active bodies of an extended source.
C
  717 CONTINUE
      IF(KWORD.EQ.KWSBOD) THEN
        READ(BUFFER,*) KB
        IF(KB.LE.0.OR.KB.GT.NB) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect body label.'
          STOP 'Incorrect body label.'
        ENDIF
        WRITE(26,7620) KB
 7620   FORMAT(3X,'Active body = ',I4)
        IXSBOD(KB)=1
        LEXBD=.TRUE.
        KBSMAX=MAX(KBSMAX,KB)
  766   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 766
        IF(KWORD.EQ.KWSBOD) GO TO 717
      ENDIF
C
C  ****  Angular distribution of primary particles.
C
      ISOURC=0
  777 CONTINUE
      IF(KWORD.EQ.KWSCON) THEN
        LSCONE=.TRUE.
        READ(BUFFER,*) STHETA,SPHI,SALPHA
        IF(STHETA.LT.-1.0D-9.OR.STHETA-180.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''THETA must be between 0 and 180 deg.'')')
          STOP 'THETA must be between 0 and 180 deg'
        ENDIF
        IF(SPHI.LT.-1.0D-9.OR.SPHI-360.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''PHI must be between 0 and 360 deg.'')')
          STOP 'PHI must be between 0 and 360 deg'
        ENDIF
        IF(SALPHA.LT.-1.0D-9.OR.SALPHA-180.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''ALPHA must be between 0 and 180 deg.'')')
          STOP 'ALPHA must be between 0 and 180 deg'
        ENDIF
        CALL GCONE0(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA)
        ISOURC=1
   18   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 18
        IF(KWORD.EQ.KWPSFN) GO TO 721
        GO TO 777
      ELSE IF(KWORD.EQ.KWSREC) THEN
        LSCONE=.FALSE.
        READ(BUFFER,*) THETLD,THETUD,PHILD,PHIUD
        IF(MIN(THETLD,THETUD).LT.-1.0D-9.OR.
     1     MAX(THETLD,THETUD)-180.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''THETA must be between 0 and 180 deg.'')')
          STOP 'THETA must be between 0 and 180 deg'
        ENDIF
        IF(MIN(PHILD,PHIUD).LT.-1.0D-9.OR.
     1     MAX(PHILD,PHIUD)-360.0D0.GT.1.0D-9) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,'(/3X,''PHI must be between 0 and 360 deg.'')')
          STOP 'PHI must be between 0 and 360 deg'
        ENDIF
        CTHL=COS(THETLD*DE2RA)
        DCTH=COS(THETUD*DE2RA)-CTHL
        PHIL=PHILD*DE2RA
        DPHI=PHIUD*DE2RA-PHIL
        ISOURC=1
   19   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 19
        IF(KWORD.EQ.KWPSFN) GO TO 721
        GO TO 777
      ELSE IF(ISOURC.EQ.0) THEN
        LSCONE=.TRUE.
        STHETA=0.0D0
        SPHI=0.0D0
        SALPHA=0.0D0
        CALL GCONE0(STHETA*DE2RA,SPHI*DE2RA,SALPHA*DE2RA)
      ENDIF
C
  721 CONTINUE
      IF(LSCONE) THEN
        WRITE(26,1133) STHETA,SPHI
 1133   FORMAT(/3X,'*** Conical beam:'
     1    /3X,'Beam axis direction:     THETA = ',1P,E13.6,' deg',
     1    /30X,'PHI = ',E13.6,' deg')
        WRITE(26,1134) SALPHA
 1134   FORMAT(3X,'Beam aperture:',11X,'ALPHA = ',1P,E13.6,' deg')
      ELSE
        WRITE(26,1733) THETLD,THETUD,PHILD,PHIUD
 1733   FORMAT(/3X,'*** Rectangular beam:'
     1    /3X,'Angular window: THETA = (',1P,E13.6,',',E13.6,') deg',
     1    /21X,'PHI = (',E13.6,',',E13.6,') deg')
      ENDIF
C
C  ************  Particle state variables read from a phase-space file.
C
   21 CONTINUE
      IF(KWORD.EQ.KWPSFN) THEN
        IF(KPARP.EQ.0) THEN
          WRITE(26,'(/3X,''With KPARP=0 (subroutine SOURCE activat'',
     1      ''ed),'')')
          WRITE(26,'(3X,''we cannot read particles from a phase-'',
     1      ''space file.'')')
          STOP 'Inconsistent definition of the primary source.'
        ENDIF
        NPSF=NPSF+1
        IF(NPSF.EQ.1) THEN
          WRITE(26,1200)
 1200     FORMAT(/3X,72('-'),/3X,'>>>>>>  Input phase-space files.'/)
        ENDIF
        IF(NPSF.GT.NPSFM) THEN
          WRITE(26,'(/3X,''Too many phase-space files.'')')
          STOP 'Too many phase-space files'
        ENDIF
        READ(BUFFER,'(A20)') PSFI(NPSF)
        WRITE(26,1201) NPSF,PSFI(NPSF)
 1201   FORMAT(3X,'Phase-space file #',I4,': ',A20)
        LPSF=.TRUE.
        LSPEC=.FALSE.
   22   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 22
        IF(KWORD.EQ.KWPSFN) GO TO 21
C
        IF(KWORD.EQ.KWPSPL) THEN
          READ(BUFFER,*) NSPLIT
          IF(NSPLIT.GT.1.AND.NSPLIT.LE.1000) THEN
            WRITE(26,1210) NSPLIT
 1210       FORMAT(/3X,'Particle splitting number = ',I3)
          ELSE
            NSPLIT=MIN(1000,ABS(NSPLIT))
            WRITE(26,1211) NSPLIT
 1211       FORMAT(/3X,'Particle splitting number = ',I3,
     1        ' (modified)')
          ENDIF
   23     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 23
        ENDIF
C
        IF(KWORD.EQ.KWRRSP) THEN
          READ(BUFFER,*) WGMIN,WGMAX
          WGMIN=ABS(WGMIN)
          WGMAX=MIN(ABS(WGMAX),1.0D10)
          IF(WGMIN.GT.WGMAX) THEN
            WRITE(26,*) 'WGMIN =',WGMIN
            WRITE(26,*) 'WGMAX =',WGMAX
            WRITE(26,'(/3X,''Inconsistent window end points.'')')
            STOP 'Inconsistent window end points.'
          ENDIF
          WRITE(26,1291) WGMIN,WGMAX
 1291     FORMAT(/3X,'Initial weight window = (',1P,E13.6,',',
     1      E13.6,')')
   24     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 24
        ENDIF
        RWGMIN=1.0D0/WGMIN
      ENDIF
C
      IF(KPARP.EQ.0) THEN
        LSPEC=.TRUE.  ! Generate the energy spectrum from the source.
!        CALL GENEVENT(1,1,PTYPE,DECAYPAR,NCNT,NINCREASE)
        CALL SOURCE
        IF(LSPEC) THEN
          NSEB=200
          DSHE=1.000000001D0*EPMAX/DBLE(NSEB)
          RDSHE=1.0D0/DSHE
        ENDIF
      ENDIF
C
C  ****  Maximum particle energy.
C
      IF(KWORD.EQ.KWEMAX) THEN
        READ(BUFFER,*) EPMAXR
        IF(KPARP.NE.0) EPMAX=EPMAXR
        WRITE(26,1220) EPMAX
 1220   FORMAT(/3X,'Maximum particle energy = ',1P,E13.6,' eV')
   25   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 25
      ELSE
        IF(LPSF) THEN
          EPMAX=1.0D9
          WRITE(26,1220) EPMAX
          WRITE(26,'(3X,''WARNING: You should have specified the '',
     1      ''maximum energy EPMAX.'')')
        ENDIF
        WRITE(26,1220) EPMAX
      ENDIF
C
C  ************  Material data and simulation parameters.
C
      WRITE(26,1300)
 1300 FORMAT(/3X,72('-'),/
     1  3X,'>>>>>>  Material data and simulation parameters.')
C
C  ****  Simulation parameters.
C
      DO M=1,MAXMAT
        EABS(1,M)=0.010D0*EPMAX
        EABS(2,M)=0.001D0*EPMAX
        EABS(3,M)=0.010D0*EPMAX
        C1(M)=0.10D0
        C2(M)=0.10D0
        WCC(M)=EABS(1,M)
        WCR(M)=EABS(2,M)
      ENDDO
      DO IB=1,NB
        DSMAX(IB)=1.0D0
      ENDDO
C
      NMAT=0
   31 CONTINUE
      IF(KWORD.EQ.KWMATF) THEN
        NMAT=NMAT+1
        READ(BUFFER,'(A20)') PMFILE(NMAT)
   32   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 32
        IF(KWORD.EQ.KWMATF) GO TO 31
      ENDIF
C
      IF(KWORD.EQ.KWSIMP) THEN
        READ(BUFFER,*) EABS(1,NMAT),EABS(2,NMAT),EABS(3,NMAT),
     1    C1(NMAT),C2(NMAT),WCC(NMAT),WCR(NMAT)
   33   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 33
        IF(KWORD.EQ.KWMATF) GO TO 31
      ENDIF
C
      IF(NMAT.EQ.0) THEN
        WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
        WRITE(26,*) 'You have to specify a material file (line MFNAME).'
        STOP 'You have to specify a material file (line MFNAME).'
      ENDIF
      IF(NMAT.GT.MAXMAT) THEN
        WRITE(26,*) 'Wrong number of materials.'
        WRITE(26,'(''NMAT ='',I4,'' is larger than MAXMAT ='',I4)')
     1    NMAT,MAXMAT
        STOP 'Wrong number of materials.'
      ENDIF
C
      DO M=1,NMAT
        IF(M.EQ.1) LIT='st'
        IF(M.EQ.2) LIT='nd'
        IF(M.EQ.3) LIT='rd'
        IF(M.GT.3) LIT='th'
        WRITE(26,1320) M,LIT
 1320   FORMAT(/3X,'**** ',I2,A2,' material')
        WRITE(26,1325) PMFILE(M)
 1325   FORMAT(3X,'Material data file: ',A)
        IF(EABS(1,M).LT.5.0D1) EABS(1,M)=5.0D1
        IF(EABS(2,M).LT.5.0D1) EABS(2,M)=5.0D1
        IF(EABS(3,M).LT.5.0D1) EABS(3,M)=5.0D1
        WRITE(26,1321) EABS(1,M)
 1321   FORMAT(3X,'Electron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1322) EABS(2,M)
 1322   FORMAT(3X,'  Photon absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1323) EABS(3,M)
 1323   FORMAT(3X,'Positron absorption energy = ',1P,E13.6,' eV')
        WRITE(26,1324) C1(M),C2(M),WCC(M),WCR(M)
 1324   FORMAT(3X,'Electron-positron simulation parameters:',
     1    /4X,'C1 =',1P,E13.6,',      C2 =',E13.6,/3X,'Wcc =',E13.6,
     1    ' eV,  Wcr =',E13.6,' eV')
      ENDDO
C
C  ****  Initialisation of PENELOPE.
C
      WRITE(6,*) '  Initialising PENELOPE ...'
      IWR=16
      OPEN(IWR,FILE='data/material.dat')
        INFO=1
        CALL PEINIT(EPMAX,NMAT,IWR,INFO,PMFILE)
      CLOSE(IWR)
C  ----  Inverse densities are used to score the local dose.
      DO M=1,NMAT
        RHOI(M)=1.0D0/RHO(M)
      ENDDO
C
C  ************  Geometry definition.
C
C  Define here the geometry parameters that are to be altered, if any.
C     PARINP(1)=
C     PARINP(2)=  ...
      NPINP=0
C
      WRITE(6,*) '  Initialising PENGEOM ...'
      IF(KWORD.EQ.KWGEOM) THEN
        READ(BUFFER,'(A20)') PFILE
        WRITE(26,1340) PFILE
 1340   FORMAT(/3X,72('-'),/3X,'>>>>>>  Geometry definition.',
     1    /3X,'PENGEOM''s geometry file: ',A20)
        OPEN(15,FILE=PFILE,IOSTAT=KODE)
        IF(KODE.NE.0) THEN
          WRITE(26,'(''File '',A20,'' could not be opened.'')') PFILE
          STOP
        ENDIF
        OPEN(16,FILE='data/geometry.rep')
        CALL GEOMIN(PARINP,NPINP,NMATG,NBODY,15,16)
        CLOSE(15)
        CLOSE(16)
        IF(NMATG.LT.1) THEN
          WRITE(26,*) 'NMATG must be greater than 0.'
          STOP 'NMATG must be greater than 0.'
        ENDIF
C
        IF(NBODY.GT.NB) THEN
          WRITE(26,'(/6X,''Too many bodies.'')')
          STOP 'Too many bodies.'
        ENDIF
C
        IF(NMATG.GT.NMAT) THEN
          WRITE(26,'(/6X,''Too many different materials.'')')
          STOP 'Too many different materials.'
        ENDIF
C
        IF(KBSMAX.GT.NBODY) THEN
          WRITE(26,'(/6X,''KBSMAX = '',I4)') KBSMAX
          WRITE(26,'(6X,'' NBODY = '',I4)') NBODY
          WRITE(26,'(6X,''Some source bodies are undefined. STOP.'')')
          STOP 'Some source bodies are undefined.'
        ENDIF
C
   34   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 34
      ELSE
        WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
        WRITE(26,*) 'You have to specify a geometry file.'
        STOP 'You have to specify a geometry file.'
      ENDIF
C
C  ****  Maximum step lengths of electrons and positrons.
C
      IF(KWORD.EQ.KWSMAX) THEN
        READ(BUFFER,*) IB
        IF(IB.LT.1.OR.IB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect body number.'
          STOP 'Incorrect body number.'
        ENDIF
        READ(BUFFER,*) IB,DSMAX(IB)
        IF(DSMAX(IB).LT.1.0D-7) DSMAX(IB)=1.0D20
   35   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 35
        IF(KWORD.EQ.KWSMAX) THEN
          READ(BUFFER,*) IB
          IF(IB.LT.1.OR.IB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body number.'
            STOP 'Incorrect body number.'
          ENDIF
          READ(BUFFER,*) IB,DSMAX(IB)
          IF(DSMAX(IB).LT.1.0D-7) DSMAX(IB)=1.0D20
          GO TO 35
        ENDIF
      ENDIF
C
C  ****  Local absorption energies (useful to reduce simulation work
C        in regions of lesser interest).
C
      DO IB=1,NBODY
        M=MATER(IB)
        IF(M.GT.0) THEN
          EABSB(1,IB)=EABS(1,M)
          EABSB(2,IB)=EABS(2,M)
          EABSB(3,IB)=EABS(3,M)
        ENDIF
      ENDDO
C
      IF(KWORD.EQ.KWEABS) THEN
   36   CONTINUE
        READ(BUFFER,*) IB
        IF(IB.LT.1.OR.IB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect body number.'
          STOP 'Incorrect body number.'
        ENDIF
        READ(BUFFER,*) JB,EAB1,EAB2,EAB3
        IF(MATER(IB).GT.0) THEN
          EABSB(1,IB)=MAX(EABSB(1,IB),EAB1)
          EABSB(2,IB)=MAX(EABSB(2,IB),EAB2)
          EABSB(3,IB)=MAX(EABSB(3,IB),EAB3)
        ENDIF
   37   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 37
        IF(KWORD.EQ.KWEABS) GO TO 36
      ENDIF
C
      WRITE(26,1350)
 1350 FORMAT(/9X,'Maximum allowed step lengths of',
     1  ' electrons and positrons',/9X,'and local absorption ',
     2  'energies (non-void bodies).',
     3  //3X,'Body',3X,'DSMAX(IB)',4X,'EABSB(1,IB)',3X,
     3  'EABSB(2,IB)',3X,'EABSB(3,IB)',/4X,'IB',7X,'(cm)',10X,'(eV)',
     4  10X,'(eV)',10X,'(eV)')
      DO IB=1,NBODY
        IF(MATER(IB).GT.0) WRITE(26,1351) IB,DSMAX(IB),EABSB(1,IB),
     1    EABSB(2,IB),EABSB(3,IB)
 1351   FORMAT(3X,I4,1P,4E14.6)
      ENDDO
C
C  ************  Variance reduction (only interaction forcing).
C
      DO KB=1,NB
        DO ICOL=1,8
          DO KPAR=1,3
            FORCE(KB,KPAR,ICOL)=1.0D0
          ENDDO
        ENDDO
        DO KPAR=1,3
          LFORCE(KB,KPAR)=.FALSE.
          WLOW(KB,KPAR)=0.0D0
          WHIG(KB,KPAR)=1.0D6
        ENDDO
      ENDDO
C
      IF(KWORD.EQ.KWIFOR) THEN
        WRITE(26,1400)
 1400   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  Interaction forcing: FORCE(IBODY,KPAR,ICOL)')
   41   CONTINUE
        READ(BUFFER,*) KB,KPAR,ICOL,FORCER,WWLOW,WWHIG
C  ****  Negative FORCER values are re-interpreted, as described in the
C        heading comments above.
        IF(FORCER.LT.-1.0D-6) THEN
          MM=MATER(KB)
          EVENTS=MAX(ABS(FORCER),1.0D0)
          PLT=PRANGE(E0,KPAR,MM)
          RMFP=PHMFP(E0,KPAR,MM,ICOL)
          IF(RMFP.LT.1.0D25) THEN
            FORCER=EVENTS*RMFP/PLT
          ELSE
            FORCER=1.0D0
          ENDIF
        ENDIF
        IF(WWLOW.LT.1.0D-6) WWLOW=1.0D-6
        IF(WWHIG.GT.1.0D6) WWHIG=1.0D6
        IF(KB.LT.1.OR.KB.GT.NBODY) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect KB value.'
          STOP 'Incorrect KB value.'
        ENDIF
        IF(KPAR.LT.1.OR.KPAR.GT.3) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of KPAR.'
          STOP 'Incorrect value of KPAR.'
        ENDIF
        IF(ICOL.LT.1.OR.ICOL.GT.8) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect value of ICOL.'
          STOP 'Incorrect value of ICOL.'
        ENDIF
        WLOW(KB,KPAR)=MAX(WLOW(KB,KPAR),WWLOW)
        WHIG(KB,KPAR)=MIN(WHIG(KB,KPAR),WWHIG)
        IF(WLOW(KB,KPAR).GT.WHIG(KB,KPAR)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect weight window limits.'
          STOP 'Incorrect weight window limits.'
        ENDIF
        IF(FORCER.GT.1.0001D0) THEN
          LFORCE(KB,KPAR)=.TRUE.
          FORCE(KB,KPAR,ICOL)=FORCER
          WRITE(26,1410) KB,KPAR,ICOL,FORCER,WLOW(KB,KPAR),WHIG(KB,KPAR)
 1410     FORMAT(3X,'FORCE(',I4,',',I1,',',I1,') =',1P,E13.6,
     1      ',  weight window = (',E9.2,',',E9.2,')')
        ENDIF
   42   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 42
        IF(KWORD.EQ.KWIFOR) GO TO 41
      ENDIF
C
C  ************  Energy and angular distributions of emerging
C                particles.
C
      WRITE(26,1500)
 1500 FORMAT(/3X,72('-'),/
     1  3X,'>>>>>>  Energy and angular distributions of emerging',
     1  ' particles.')
C
      IF(KWORD.EQ.KWNBE) THEN
        READ(BUFFER,*) EMIN,EMAX,NBE
   51   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 51
      ELSE
        EMIN=0.0D0
        EMAX=EPMAX
        NBE=100
      ENDIF
      EMIN=MAX(EMIN,0.0D0)
      EMAX=MAX(EMIN+50.0D0,EMAX)
      WRITE(26,1510) NBE,EMIN,EMAX
 1510 FORMAT(3X,'E:       NBE = ',I3,
     1  ',  EMIN =',1P,E13.6,' eV,  EMAX =',E13.6,' eV')
      IF(NBE.LT.1) THEN
        WRITE(26,*) 'Wrong number of energy bins.'
        STOP 'Wrong number of energy bins.'
      ELSE IF (NBE.GT.NBEM) THEN
        WRITE(26,*) 'NBE is too large.'
        WRITE(26,*) 'Set the parameter NBEM equal to ',NBE
        STOP 'NBE is too large.'
      ENDIF
C
      IF(KWORD.EQ.KWNBAN) THEN
        READ(BUFFER,*) NBTH,NBPH
   52   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 52
      ELSE
        NBTH=90
        NBPH=1
      ENDIF
      WRITE(26,1520) NBTH
 1520 FORMAT(3X,'Theta:  NBTH = ',I3)
      WRITE(26,1530) NBPH
 1530 FORMAT(3X,'Phi:    NBPH = ',I3)
      IF(NBTH.LT.1) THEN
        WRITE(26,*) 'Wrong number of THETA bins.'
        STOP 'Wrong number of THETA bins.'
      ELSE IF (NBTH.GT.NBTHM) THEN
        WRITE(26,*) 'NBTH is too large.'
        WRITE(26,*) 'Set the parameter NBTHM equal to ',NBTH
        STOP 'NBTH is too large.'
      ENDIF
      IF(NBPH.LT.1) THEN
        WRITE(26,*) 'Wrong number of PHI bins.'
        STOP 'Wrong number of PHI bins.'
      ELSE IF (NBPH.GT.NBPHM) THEN
        WRITE(26,*) 'NBPH is too large.'
        WRITE(26,*) 'Set the parameter NBPHM equal to ',NBPH
        STOP 'NBPH is too large.'
      ENDIF
C
C  ****  Bin sizes.
C  ----  The factor FSAFE=1.000000001 serves to ensure that the upper
C  limit of the tallied interval is within the last bin (otherwise, the
C  array dimensions could be exceeded).
C
      FSAFE=1.000000001D0
      BSE=FSAFE*(EMAX-EMIN)/DBLE(NBE)
      RBSE=1.0D0/BSE
      BSTH=FSAFE*180.0D0/DBLE(NBTH)
      RBSTH=1.0D0/BSTH
      BSPH=FSAFE*360.0D0/DBLE(NBPH)
      RBSPH=1.0D0/BSPH
C
C  ************  Impact detectors.
C
      DO KD=1,NIDM
        KKDI(KD,1)=0
        KKDI(KD,2)=0
        KKDI(KD,3)=0
        TDID(KD)=0.0D0
        TDID2(KD)=0.0D0
      ENDDO
      NDBOD=0
      NPSFO=0
C
      NDIDEF=0
   61 CONTINUE
      IF(KWORD.EQ.KWIDET) THEN
        NDIDEF=NDIDEF+1
        NDBOD=0
        IF(NDIDEF.GT.NIDM) THEN
          WRITE(26,'(3X,''NDIDEF = '',I4)') NDIDEF
          WRITE(26,*) 'Too many detectors.'
          STOP 'Too many detectors.'
        ENDIF
        WRITE(26,1600) NDIDEF
 1600   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  Impact detector #', I2)
        READ(BUFFER,*) EDIL(NDIDEF),EDIU(NDIDEF),NDICH(NDIDEF),
     1    IPSF(NDIDEF),IDCUT(NDIDEF)
        IF(EDIL(NDIDEF).LT.0.0D0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'EDIL must be positive.'
          STOP 'EDIL must be positive.'
        ENDIF
        IF(EDIU(NDIDEF).LT.EDIL(NDIDEF)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect energy limits.'
          STOP 'Incorrect energy limits.'
        ENDIF
        IF(NDICH(NDIDEF).LT.0) THEN
          LDILOG(NDIDEF)=.TRUE.  ! Logarithmic bins.
          NDICH(NDIDEF)=ABS(NDICH(NDIDEF))
        ENDIF
        IF(NDICH(NDIDEF).LT.10.OR.NDICH(NDIDEF).GT.NIDCM) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect number of energy bins.'
          STOP 'Incorrect number of energy bins.'
        ENDIF
C
        IF(LDILOG(NDIDEF)) THEN
          WRITE(26,1614) EDIL(NDIDEF),EDIU(NDIDEF),NDICH(NDIDEF)
 1614     FORMAT(3X,'Energy window = (',1P,E12.5,',',E12.5,') eV',
     1      /3X,'Number of energy bins = ',I4,'  (logarithmic mesh)')
          EDIUL(NDIDEF)=LOG(EDIU(NDIDEF))
          EDILL(NDIDEF)=LOG(EDIL(NDIDEF))
          BDIEL(NDIDEF)=FSAFE*(EDIUL(NDIDEF)-EDILL(NDIDEF))
     1      /DBLE(NDICH(NDIDEF))
          RBDIEL(NDIDEF)=1.0D0/BDIEL(NDIDEF)
        ELSE
          WRITE(26,1610) EDIL(NDIDEF),EDIU(NDIDEF),NDICH(NDIDEF)
 1610     FORMAT(3X,'Energy window = (',1P,E12.5,',',E12.5,') eV',
     1      /3X,'Number of energy bins = ',I4,'  (uniform mesh)')
          BDIE(NDIDEF)=FSAFE*(EDIU(NDIDEF)-EDIL(NDIDEF))
     1      /DBLE(NDICH(NDIDEF))
          RBDIE(NDIDEF)=1.0D0/BDIE(NDIDEF)
        ENDIF
C
        IF(IPSF(NDIDEF).LT.0.OR.ABS(IPSF(NDIDEF)).GT.1) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Wrong IPSF value.'
          STOP 'Wrong IPSF value.'
        ENDIF
C
        IF(IDCUT(NDIDEF).LT.0.OR.ABS(IDCUT(NDIDEF)).GT.2) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Wrong IDCUT value.'
          STOP 'Wrong IDCUT value.'
        ENDIF
C
   62   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 62
C
        IF(KWORD.EQ.KWISPC) THEN
          IF(NDIDEF.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,'(A20)') SPCDIO(NDIDEF)
          WRITE(26,1612) SPCDIO(NDIDEF)
 1612     FORMAT(3X,'Output energy spectrum: ',A20)
   64     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 64
        ELSE
          WRITE(BUF2,'(I5)') 1000+NDIDEF
          SPCDIO(NDIDEF)='data/spc-impdet-'//BUF2(4:5)//'.dat'
          WRITE(26,1612) SPCDIO(NDIDEF)
        ENDIF
C
        IF(KWORD.EQ.KWIPSF) THEN
          IF(NDIDEF.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,'(A20)') PSFDIO(NDIDEF)
          IF(IPSF(NDIDEF).GT.0) THEN
            WRITE(26,1611) PSFDIO(NDIDEF)
 1611       FORMAT(3X,'Output phase-space file: ',A20)
            NPSFO=NPSFO+1
            IF(NPSFO.GT.1) THEN
              WRITE(26,*) 'You cannot generate more than one PSF in ',
     1          'a single run.'
              STOP 'Only one PSF can be generated in a each run.'
            ENDIF
          ELSE
            WRITE(26,1613)
 1613       FORMAT(3X,'No phase-space file is generated.')
          ENDIF
   63     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 63
        ELSE
          IF(IPSF(NDIDEF).GT.0) THEN
            WRITE(BUF2,'(I5)') 1000+NDIDEF
            PSFDIO(NDIDEF)='data/psf-impdet-'//BUF2(4:5)//'.dat'
            WRITE(26,1611) PSFDIO(NDIDEF)
            NPSFO=NPSFO+1
            IF(NPSFO.GT.1) THEN
              WRITE(26,*) 'You cannot generate more than one PSF in ',
     1          'a single run.'
              STOP 'Only one PSF can be generated in each run.'
            ENDIF
          ENDIF
        ENDIF
C
        IF(ABS(IDCUT(NDIDEF)).EQ.0) THEN
          WRITE(26,1601)
 1601     FORMAT(3X,'Detected particles are absorbed')
        ELSE
          WRITE(26,1602)
 1602     FORMAT(3X,'Particles are transported through this detector')
        ENDIF
C
        IF(KWORD.EQ.KWIFLN) THEN
          IF(NDIDEF.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,'(A20)') SPCFSO(NDIDEF)
   83     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 83
          IF(IDCUT(NDIDEF).EQ.2) THEN
            WRITE(26,1662) SPCFSO(NDIDEF)
 1662       FORMAT(3X,'Output fluence distribution: ',A20)
          ENDIF
        ELSE
          IF(IDCUT(NDIDEF).EQ.2) THEN
            WRITE(BUF2,'(I5)') 1000+NDIDEF
            SPCFSO(NDIDEF)='data/fln-impdet-'//BUF2(4:5)//'.dat'
            WRITE(26,1662) SPCFSO(NDIDEF)
          ENDIF
        ENDIF
C
   65   CONTINUE
        IF(KWORD.EQ.KWIBOD) THEN
          IF(NDIDEF.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          READ(BUFFER,*) KB
          IF(KB.LE.0.OR.KB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body label.'
            STOP 'Incorrect body label.'
          ENDIF
          IF(KDET(KB).NE.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A body cannot be part of two detectors.'
            STOP 'A body cannot be part of two detectors.'
          ENDIF
          IF(MATER(KB).EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A void body cannot be part of a detector.'
            STOP 'A void body cannot be part of a detectors.'
          ENDIF
          WRITE(26,1620) KB
 1620     FORMAT(3X,'Active body = ',I4)
          KDET(KB)=NDIDEF
          NDBOD=NDBOD+1
   66     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 66
          IF(KWORD.EQ.KWIBOD) GO TO 65
          IF(KWORD.EQ.KWIDET) THEN
            ITST=MAX(KKDI(NDIDEF,1),KKDI(NDIDEF,2),KKDI(NDIDEF,3))
            IF(ITST.EQ.0) THEN
              KKDI(NDIDEF,1)=1
              KKDI(NDIDEF,2)=1
              KKDI(NDIDEF,3)=1
              WRITE(26,1630)
 1630         FORMAT(3X,'Detected particles = electrons, photons and ',
     1          'positrons')
            ENDIF
            IF(NDBOD.EQ.0) THEN
              WRITE(26,*) 'This detector has no active bodies.'
              STOP 'This detector has no active bodies.'
            ENDIF
            GO TO 61
          ENDIF
        ENDIF
C
   67   CONTINUE
        IF(KWORD.EQ.KWIPAR) THEN
          IF(NDIDEF.EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'No impact detector has been defined yet.'
            STOP 'No impact detector has been defined yet.'
          ENDIF
          IF(NDBOD.EQ.0) THEN
            WRITE(26,*) 'This detector has no active bodies.'
            STOP 'This detector has no active bodies.'
          ENDIF
C
          READ(BUFFER,*) KPARD
          IF(KPARD.EQ.1) THEN
            KKDI(NDIDEF,1)=1
            WRITE(26,1631)
 1631       FORMAT(3X,'Detected particles = electrons')
          ELSE IF(KPARD.EQ.2) THEN
            KKDI(NDIDEF,2)=1
            WRITE(26,1632)
 1632       FORMAT(3X,'Detected particles = photons')
          ELSE IF(KPARD.EQ.3) THEN
            KKDI(NDIDEF,3)=1
            WRITE(26,1633)
 1633       FORMAT(3X,'Detected particles = positrons')
          ENDIF
   68     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 68
          IF(KWORD.EQ.KWIPAR) GO TO 67
          IF(KWORD.EQ.KWIBOD) GO TO 65
          IF(KWORD.EQ.KWIDET) THEN
            ITST=MAX(KKDI(NDIDEF,1),KKDI(NDIDEF,2),KKDI(NDIDEF,3))
            IF(ITST.EQ.0) THEN
              KKDI(NDIDEF,1)=1
              KKDI(NDIDEF,2)=1
              KKDI(NDIDEF,3)=1
              WRITE(26,1630)
            ENDIF
            GO TO 61
          ENDIF
        ENDIF
      ENDIF
C
      IF(NDIDEF.GT.0) THEN
        IF(NDBOD.EQ.0) THEN
          WRITE(26,*) 'This detector has no active bodies.'
          STOP 'This detector has no active bodies.'
        ENDIF
        ITST=MAX(KKDI(NDIDEF,1),KKDI(NDIDEF,2),KKDI(NDIDEF,3))
        IF(ITST.EQ.0) THEN
          KKDI(NDIDEF,1)=1
          KKDI(NDIDEF,2)=1
          KKDI(NDIDEF,3)=1
          WRITE(26,1630)
        ENDIF
      ENDIF
C
C  ************  Energy-deposition detectors.
C
      DO KB=1,NBODY
        KBDE(KB)=0
      ENDDO
      DO KD=1,NEDM
        TDED(KD)=0.0D0
        TDED2(KD)=0.0D0
      ENDDO
      NDBOD=0
C
      NDEDEF=0
   43 CONTINUE
      IF(KWORD.EQ.KWEDET) THEN
        IF(NDEDEF.GT.0) THEN
          IF(NDBOD.EQ.0) THEN
            WRITE(26,*) 'This detector has no active bodies.'
            STOP 'This detector has no active bodies.'
          ENDIF
        ENDIF
        NDEDEF=NDEDEF+1
        NDBOD=0
        IF(NDEDEF.GT.NEDM) THEN
          WRITE(26,'(3X,''NDEDEF = '',I4)') NDEDEF
          WRITE(26,*) 'Too many energy-deposition detectors.'
          STOP 'Too many energy-deposition detectors.'
        ENDIF
        WRITE(26,1650) NDEDEF
 1650   FORMAT(/3X,72('-'),/
     1    3X,'>>>>>>  Energy-deposition detector #', I2)
        READ(BUFFER,*) EDEL(NDEDEF),EDEU(NDEDEF),NDECH(NDEDEF)
        IF(EDEL(NDEDEF).LT.0.0D0) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'EDEL must be positive.'
          STOP 'EDEL must be positive.'
        ENDIF
        IF(EDEU(NDEDEF).LT.EDEL(NDEDEF)) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect energy limits.'
          STOP 'Incorrect energy limits.'
        ENDIF
        IF(NDECH(NDEDEF).LT.0) THEN
          LDELOG(NDEDEF)=.TRUE.  ! Logarithmic bins.
          NDECH(NDEDEF)=ABS(NDECH(NDEDEF))
        ENDIF
        IF(NDECH(NDEDEF).LT.10.OR.NDECH(NDEDEF).GT.NEDCM) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Incorrect number of bins.'
          STOP 'Incorrect number of bins.'
        ENDIF
C
        IF(LDELOG(NDEDEF)) THEN
          WRITE(26,1614) EDEL(NDEDEF),EDEU(NDEDEF),NDECH(NDEDEF)
          EDEUL(NDEDEF)=LOG(EDEU(NDEDEF))
          EDELL(NDEDEF)=LOG(EDEL(NDEDEF))
          BDEEL(NDEDEF)=FSAFE*(EDEUL(NDEDEF)-EDELL(NDEDEF))
     1      /DBLE(NDECH(NDEDEF))
          RBDEEL(NDEDEF)=1.0D0/BDEEL(NDEDEF)
        ELSE
          WRITE(26,1610) EDEL(NDEDEF),EDEU(NDEDEF),NDECH(NDEDEF)
          BDEE(NDEDEF)=FSAFE*(EDEU(NDEDEF)-EDEL(NDEDEF))
     1      /DBLE(NDECH(NDEDEF))
          RBDEE(NDEDEF)=1.0D0/BDEE(NDEDEF)
        ENDIF
C
   44   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 44
        IF(KWORD.EQ.KWESPC) THEN
          READ(BUFFER,'(A20)') SPCDEO(NDEDEF)
          WRITE(26,1651) SPCDEO(NDEDEF)
 1651     FORMAT(3X,'Output spectrum: ',A20)
   45     CONTINUE
          READ(5,'(A6,1X,A120)') KWORD,BUFFER
          IF(KWORD.EQ.KWCOMM) GO TO 45
        ELSE
          WRITE(BUF2,'(I5)') 1000+NDEDEF
          SPCDEO(NDEDEF)='data/spc-enddet-'//BUF2(4:5)//'.dat'
          WRITE(26,1651) SPCDEO(NDEDEF)
        ENDIF
C
   46   CONTINUE
        IF(KWORD.EQ.KWEBOD) THEN
          READ(BUFFER,*) KB
          IF(KB.LT.0.OR.KB.GT.NBODY) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'Incorrect body label.'
            STOP 'Incorrect body label.'
          ENDIF
          IF(KBDE(KB).NE.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A body cannot be part of two detectors.'
            STOP 'A body cannot be part of two detectors.'
          ENDIF
          IF(MATER(KB).EQ.0) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'A void body cannot be part of a detector.'
            STOP 'A void body cannot be part of a detectors.'
          ENDIF
          WRITE(26,1652) KB
 1652     FORMAT(3X,'Active body = ',I4)
          IF(LFORCE(KB,1).OR.LFORCE(KB,2).OR.LFORCE(KB,3)) THEN
            WRITE(26,'(3X,''#  WARNING: Spectrum may be strongly'',
     1        '' biased'',/15X,''when interaction forcing is used!'')')
          ENDIF
          KBDE(KB)=NDEDEF
          NDBOD=NDBOD+1
        ENDIF
   47   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 47
        IF(KWORD.EQ.KWEDET) GO TO 43
        IF(KWORD.EQ.KWEBOD) GO TO 46
      ENDIF
C
      IF(NDEDEF.GT.0) THEN
        IF(NDBOD.EQ.0) THEN
          WRITE(26,*) 'This detector has no active bodies.'
          STOP 'This detector has no active bodies.'
        ENDIF
      ENDIF
C
C  ************  Dose distribution.
C
      LDOSEM=.FALSE.
      IF(KWORD.EQ.KGRDXX) THEN
        LDOSEM=.TRUE.
        WRITE(26,1700)
 1700   FORMAT(/3X,72('-'),/3X,'>>>>>>  3D Dose distribution.')
        READ(BUFFER,*) DXL(1),DXU(1)
        IF(DXL(1).GT.DXU(1)) THEN
          SAVE=DXL(1)
          DXL(1)=DXU(1)
          DXU(1)=SAVE
        ENDIF
        IF(DXU(1).LT.DXL(1)+1.0D-6) THEN
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'XU must be greater than XL+1.0E-6.'
          STOP 'XU must be greater than XL+1.0E-6.'
        ENDIF
        WRITE(26,1710) DXL(1),DXU(1)
 1710   FORMAT(3X,'Dose-map box:  XL = ',1P,E13.6,' cm,  XU = ',
     1    E13.6,' cm')
   71   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 71
        IF(KWORD.EQ.KGRDYY) THEN
          READ(BUFFER,*) DXL(2),DXU(2)
          IF(DXL(2).GT.DXU(2)) THEN
            SAVE=DXL(2)
            DXL(2)=DXU(2)
            DXU(2)=SAVE
          ENDIF
          IF(DXU(2).LT.DXL(2)+1.0D-6) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'YU must be greater than YL+1.0E-6.'
            STOP 'YU must be greater than YL+1.0E-6.'
          ENDIF
          WRITE(26,1711) DXL(2),DXU(2)
 1711     FORMAT(18X,'YL = ',1P,E13.6,' cm,  YU = ',E13.6,' cm')
        ELSE
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Unrecognized keyword.'
          STOP 'Unrecognized keyword.'
        ENDIF
   72   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 72
        IF(KWORD.EQ.KGRDZZ) THEN
          READ(BUFFER,*) DXL(3),DXU(3)
          IF(DXL(3).GT.DXU(3)) THEN
            SAVE=DXL(3)
            DXL(3)=DXU(3)
            DXU(3)=SAVE
          ENDIF
          IF(DXU(3).LT.DXL(3)+1.0D-6) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,*) 'ZU must be greater than ZL+1.0E-6.'
            STOP 'ZU must be greater than ZL+1.0E-6.'
          ENDIF
          WRITE(26,1712) DXL(3),DXU(3)
 1712     FORMAT(18X,'ZL = ',1P,E13.6,' cm,  ZU = ',E13.6,' cm')
        ELSE
          WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
          WRITE(26,*) 'Unrecognized keyword.'
          STOP 'Unrecognized keyword.'
        ENDIF
   73   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 73
        IF(KWORD.EQ.KGRDBN) THEN
          READ(BUFFER,*) NDBX,NDBY,NDBZ
          IF(NDBX.LT.0.OR.NDBX.GT.NDXM) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,'(''NDBX must be .GT.0. and .LT.'',I4)') NDXM
            WRITE(26,*) 'Increase the value of the parameter NDXM.'
            STOP 'NDBX must be .GT.0. and .LE.NDXM'
          ENDIF
          IF(NDBY.LT.0.OR.NDBY.GT.NDYM) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,'(''NDBY must be .GT.0. and .LT.'',I4)') NDYM
            WRITE(26,*) 'Increase the value of the parameter NDYM.'
            STOP 'NDBY must be .GT.0. and .LE.NDYM'
          ENDIF
          IF(NDBZ.LT.0.OR.NDBZ.GT.NDZM) THEN
            WRITE(26,'(A6,1X,A120)') KWORD,BUFFER
            WRITE(26,'(''NDBZ must be .GT.0. and .LT.'',I4)') NDZM
            WRITE(26,*) 'Increase the value of the parameter NDZM.'
            STOP 'NDBZ must be .GT.0. and .LE.NDZM'
          ENDIF
          NDB(1)=NDBX
          NDB(2)=NDBY
          NDB(3)=NDBZ
        ELSE
          NDB(1)=25
          NDB(2)=25
          NDB(3)=25
        ENDIF
        WRITE(26,1713) NDB(1),NDB(2),NDB(3)
 1713   FORMAT(3X,'Numbers of bins:   NDBX =',I4,', NDBY =',I4,
     1    ', NDBZ =',I4)
        DO I=1,3
          BDOSE(I)=FSAFE*(DXU(I)-DXL(I))/DBLE(NDB(I))
          RBDOSE(I)=1.0D0/BDOSE(I)
        ENDDO
   74   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 74
      ENDIF
C
C  ************  Job characteristics.
C
      WRITE(26,1800)
 1800 FORMAT(/3X,72('-'),/
     1  3X,'>>>>>>  Job characteristics.')
C
      IRESUM=0
      IF(KWORD.EQ.KWRESU) THEN
        READ(BUFFER,'(A20)') PFILER
        WRITE(26,1810) PFILER
 1810   FORMAT(3X,'Resume simulation from previous dump file: ',A20)
        IRESUM=1
   75   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 75
      ENDIF
C
      DUMPP=1.0D15
      IF(KWORD.EQ.KWDUMP) THEN
        READ(BUFFER,'(A20)') PFILED
        WRITE(26,1820) PFILED
 1820   FORMAT(3X,'Write final counter values on the dump file: ',A20)
        LDUMP=.TRUE.
   76   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 76
      ENDIF
C
      IF(KWORD.EQ.KWDMPP) THEN
        READ(BUFFER,*) DUMPP
        IF(LDUMP) THEN
          IF(DUMPP.LT.15.0D0) DUMPP=15.0D0
          IF(DUMPP.GT.86400.0D0) DUMPP=86400.0D0
          WRITE(26,1830) DUMPP
 1830     FORMAT(3X,'Dumping period: DUMPP =',1P,E13.6)
        ENDIF
   77   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 77
      ENDIF
C
      IF(KWORD.EQ.KWRSEE) THEN
        READ(BUFFER,*) ISEED1,ISEED2
   82   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 82
      ELSE
        ISEED1=1
        ISEED2=1
      ENDIF
      WRITE(26,1850) ISEED1,ISEED2
 1850 FORMAT(3X,'Random-number generator seeds = ',I10,', ',I10)
C
      IF(KWORD.EQ.KWNSIM) THEN
        READ(BUFFER,*) DSHN
        IF(DSHN.LT.1.0D0) DSHN=2.0D9
   78   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 78
      ELSE
        DSHN=2.0D9
      ENDIF
      WRITE(26,1840) DSHN
 1840 FORMAT(/3X,'Number of showers to be simulated =',1P,E13.6)
C
      IF(KWORD.EQ.KWTIME) THEN
        READ(BUFFER,*) TIMEA
   79   CONTINUE
        READ(5,'(A6,1X,A120)') KWORD,BUFFER
        IF(KWORD.EQ.KWCOMM) GO TO 79
      ELSE
        TIMEA=2.0D9
      ENDIF
      IF(TIMEA.LT.1.0D0) TIMEA=100.0D0
      WRITE(26,1860) TIMEA
 1860 FORMAT(3X,'Computation time available = ',1P,E13.6,' sec')
C
      CALL TIMER(TSECIN)
      TSECA=TIMEA+TSECIN
      TSECAD=TSECIN
      WRITE(26,1870)
 1870 FORMAT(/3X,72('-'))
 
      IF(KWORD.EQ.KHBKFN)THEN
         READ(BUFFER,'(A20)') HBKFN
         WRITE(6,*)'OUPUT FILE IS: ',HBKFN
   85    CONTINUE
         READ(5,'(A6,1X,A120)') KWORD,BUFFER
         IF(KWORD.EQ.KWCOMM) GO TO 85
      ELSE
         HBKFN = 'test.hbook'
      ENDIF
      WRITE(26,1975) HBKFN
 1975 FORMAT(3X,'HBOOK file name = ',A20)
 
      IF(KWORD.EQ.KZMAXD)THEN
         READ(BUFFER,*) ZEND
   80    CONTINUE
         READ(5,'(A6,1X,A120)') KWORD,BUFFER
         IF(KWORD.EQ.KWCOMM) GO TO 80
      ELSE
         ZEND = 1E30
      ENDIF
      WRITE(26,1960) ZEND
 1960 FORMAT(3X,'Maximum Z distance = ',1P,E13.6,' cm')
 
      IF(KWORD.EQ.KRADMX)THEN
         READ(BUFFER,*) RMAX
   81    CONTINUE
         READ(5,'(A6,1X,A120)') KWORD,BUFFER
         IF(KWORD.EQ.KWCOMM) GO TO 81
      ELSE
         RMAX = 1E30
      ENDIF
      WRITE(26,1970) RMAX
 1970 FORMAT(3X,'Maximum radial distance = ',1P,E13.6,' cm')
C     
C  ************  If 'RESUME' is active, read previously generated
C                counters...
C
      SHNA=0.0D0
      CPUTA=0.0D0
      N=0
C
      RLAST=0.0D0
      RWRITE=0.0D0
C
      IF(IRESUM.EQ.1) THEN
        WRITE(6,*) '  Reading the DUMP file ...'
        OPEN(9,FILE=PFILER)
        READ (9,*,ERR=91,END=91) SHNA,CPUTA
        READ (9,'(A65)',ERR=90) TITLE2
        IF(TITLE2.NE.TITLE) THEN
          WRITE(26,*)
     1      'The dump file is corrupted (the TITLE does not match).'
          STOP 'The dump file is corrupted (the TITLE does not match).'
        ENDIF
        READ (9,*,ERR=90) ISEED1,ISEED2
        READ (9,*,ERR=90) NPSN,RLREAD
        READ (9,*,ERR=90) (SHIST(I),I=1,NSEB)
        READ (9,*,ERR=90) (PRIM(I),I=1,3),(PRIM2(I),I=1,3)
        READ (9,*,ERR=90) ((SEC(K,I),I=1,3),K=1,3),
     1             ((SEC2(K,I),I=1,3),K=1,3)
        READ (9,*,ERR=90) (TDEBO(I),I=1,NBODY), (TDEBO2(I),I=1,NBODY)
        READ (9,*,ERR=90) (((PDE(I,J,K),K=1,NBE),J=1,2),I=1,3),
     1             (((PDE2(I,J,K),K=1,NBE),J=1,2),I=1,3)
        READ (9,*,ERR=90) (((PDA(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3),
     1             (((PDA2(I,J,K),K=1,NBPH),J=1,NBTH),I=1,3)
        READ (9,*,ERR=90) (TDID(I),I=1,NIDM), (TDID2(I),I=1,NIDM)
        READ (9,*,ERR=90) (TDED(I),I=1,NEDM), (TDED2(I),I=1,NEDM)
        IF(NDIDEF.GT.0) THEN
          READ (9,*,ERR=90) RLAST,RWRITE
          DO ID=1,NDIDEF
            READ (9,*,ERR=90) (DIT(ID,J),J=1,NDICH(ID))
            READ (9,*,ERR=90) (DIT2(ID,J),J=1,NDICH(ID))
            READ (9,*,ERR=90) ((DIP(ID,J,K),J=1,NDICH(ID)),K=1,3)
            READ (9,*,ERR=90) ((DIP2(ID,J,K),J=1,NDICH(ID)),K=1,3)
            IF(IDCUT(ID).EQ.2) THEN
              READ (9,*,ERR=90) (FST(ID,J),J=1,NDICH(ID))
              READ (9,*,ERR=90) (FST2(ID,J),J=1,NDICH(ID))
              READ (9,*,ERR=90) ((FSP(ID,J,K),J=1,NDICH(ID)),K=1,3)
              READ (9,*,ERR=90) ((FSP2(ID,J,K),J=1,NDICH(ID)),K=1,3)
            ENDIF
          ENDDO
        ENDIF
        IF(NDEDEF.GT.0) THEN
          DO ID=1,NDEDEF
            READ (9,*,ERR=90) (DET(ID,J),J=1,NDECH(ID))
          ENDDO
        ENDIF
        IF(LDOSEM) THEN
          READ (9,*,ERR=90)
     1      (((DOSE(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1)),
     1      (((DOSE2(I1,I2,I3),I3=1,NDB(3)),I2=1,NDB(2)),I1=1,NDB(1))
          READ (9,*,ERR=90)
     1      (DDOSE(I3),I3=1,NDB(3)),(DDOSE2(I3),I3=1,NDB(3))
        ENDIF
        CLOSE(9)
        WRITE(26,1880) PFILER
 1880   FORMAT(3X,'Simulation has been resumed from dump file: ',A20)
        GO TO 92
   90   CONTINUE
        WRITE(26,*) 'The dump file is empty or corrupted.'
        STOP 'The dump file is empty or corrupted.'
   91   CONTINUE
        WRITE(26,1890)
 1890   FORMAT(3X,'WARNING: Could not resume from dump file...')
        CLOSE(9)
        IRESUM=0
      ENDIF
   92 CONTINUE
C
      IPSFO=21
      IF(NDIDEF.GT.0) THEN
        DO ID=1,NDIDEF
          IF(IPSF(ID).GT.0) THEN
            OPEN(IPSFO,FILE=PSFDIO(ID),IOSTAT=KODE)
            IF(KODE.NE.0) THEN
              WRITE(26,'(''File '',A20,'' could not be opened.'')')
     1          PSFDIO(ID)
              STOP 'File could not be opened.'
            ENDIF
            RWR=0.0D0
            IF(RWRITE.GT.0) THEN
   93         CONTINUE
              KODEPS=0
              CALL RDPSF(IPSFO,NSHI,ISEC,KODEPS)
              IF(KODEPS.EQ.-1) THEN
                GO TO 94
              ELSE
                RWR=RWR+1.0D0
                IF(RWR.LT.RWRITE-0.5D0) GO TO 93
                GO TO 94
              ENDIF
            ENDIF
   94       CONTINUE
            IF(RWR.LT.0.5D0) THEN
              WRITE(IPSFO,1901) ID
 1901         FORMAT(1X,'#  Results from PENMAIN. Phase-space fi',
     1          'le of detector no.',I3,/1X,'#')
              WRITE(IPSFO,1902)
 1902         FORMAT(1X,'#/KPAR',2X,'E',12X,'X',12X,'Y',12X,'Z',12X,
     1          'U',12X,'V',12X,'W',11X,'WGHT',5X,'ILB(1:4)',7X,'NSHI',
     1          /1X,'#',125('-'))
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
      IPSFI=20
      IF(LPSF) THEN
        IF(IRESUM.EQ.1) THEN
          IF(NPSN.GT.NPSF) THEN
            WRITE(6,*) '  **** The simulation was already completed.'
            WRITE(26,*) '  **** The simulation was already completed.'
            SHN=SHNA
            TSIM=CPUTA
            CPUT0=CPUTIM()
            JOBEND=1
            RETURN
          ENDIF
          IF(NPSN.GT.1) THEN
            DO I=1,NPSN-1
              WRITE(6,1903) PSFI(I)
              WRITE(26,1903) PSFI(I)
 1903         FORMAT(/3X,'+ The phase-space file ',A20,
     1          ' was read in previous runs.')
            ENDDO
          ENDIF
   95     CONTINUE
          OPEN(IPSFI,FILE=PSFI(NPSN))
          KODEPS=0
          CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
 1904     FORMAT(/3X,'+ The phase-space file ',A20,' is opened.')
          IF(KODEPS.EQ.-1) THEN
            WRITE(26,*) ' STOP. The file is empty or corrupted.'
            STOP 'The file is empty or corrupted.'
          ENDIF
          IF(RLREAD.GT.0.5D0) THEN
            RI=0.0D0
   96       CONTINUE
            RI=RI+1.0D0
            CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
            IF(KODEPS.EQ.-1) THEN
              NPSN=NPSN+1
              IF(NPSN.GT.NPSF) THEN
                WRITE(6,*)
     1            '  **** The simulation was already completed.'
                WRITE(26,*)
     1            '  **** The simulation was already completed.'
                SHN=SHNA
                TSIM=CPUTA
                CPUT0=CPUTIM()
                JOBEND=2
                RETURN
              ELSE
                WRITE(6,1905) PSFI(NPSN)
                WRITE(26,1905) PSFI(NPSN)
 1905           FORMAT(/3X,'+ The phase-space file ',A20,
     1            ' was completed in the last run.')
                CLOSE(IPSFI)
                RLREAD=0.0D0
                GO TO 95
              ENDIF
            ENDIF
            IF(RI.LT.RLREAD-0.5D0) GO TO 96
          ELSE
            RLREAD=0.0D0
          ENDIF
        ELSE
          NPSN=1
          RLREAD=0.0D0
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
          OPEN(IPSFI,FILE=PSFI(NPSN),IOSTAT=KODE)
          IF(KODE.NE.0) THEN
            WRITE(26,'(''File '',A20,'' could not be opened.'')')
     1        PSFI(NPSN)
            STOP 'File could not be opened.'
          ENDIF
          KODEPS=0
          CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
          IF(KODEPS.EQ.-1) THEN
            WRITE(26,*) ' STOP. The file is empty or corrupted.'
            STOP 'The file is empty or corrupted.'
          ENDIF
        ENDIF
      ENDIF
C
C  ************  Initialise constants.
C
      SHN=SHNA  ! Shower counter, particles from the dump file.
      N=MOD(SHN,2.0D9)+0.5D0
      TSIM=CPUTA
      CPUT0=CPUTIM()
      IF(SHN.GT.DSHN-0.5D0) THEN
        WRITE(6,*) '  **** The simulation was already completed.'
        WRITE(26,*) '  **** The simulation was already completed.'
        JOBEND=3
      ELSE
        WRITE(6,*) '  The simulation is started ...'
      ENDIF
C
      RETURN
      END
