c----------------------------------------------------------------------c
      SUBROUTINE BUILDHBOOK(NFILE)
c----------------------------------------------------------------------c
      INCLUDE 'pmcomms.f'
      INCLUDE 'ucnapenmain.h'
c----------------------------------------------------------------------c
      CHARACTER*26 HBOOKFILE
      CHARACTER*3  FILENUMBER,FORMSTR
      CALL HLIMIT(NWPAWC) ! A FEW PAW FORUMS SUGGEST THIS FOR MAKING
      IQUEST(10) = 65000  ! BIGGER NTUPLES. KINDA WORKS BUT STILL 
      
      IF(NFILE.GT.0)THEN
          CALL HROUT(0,ICYCLE,'T')
          CALL HREND('EVENT')
      ENDIF
      
      WRITE(FILENUMBER,'(I0.3)')NFILE
      HBOOKFILE = TRIM(HBKFN)//'_'//TRIM(FILENUMBER)//'.hbook'
      
      CALL HROPEN(1,'EVENT',HBOOKFILE,'NQE',1024,IERR)
      CALL HBOOKN(34,'DECAYS',NENTRIES,'//EVENT',1024,DTAGS) 
      CALL HBOOK1(100,'INTIAL BETA DECAY SPECTRUM',100,0.,1000.,0.)
      CALL HBOOK1(110,'THROWN ASYMMETRY',100,0.,1000.,0.)
      CALL HBOOK1(120,'MEASURED ASYMMETRY (RAW)',100,0.,1000.,0.)
      CALL HBOOK1(130,'RAW MWPC SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(140,'RAW MWCP SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(150,'RAW PMT SPECTRUM - EAST',100,0.,1500.,0.)
      CALL HBOOK1(160,'RAW PMT SPECTRUM - WEST',100,0.,1500.,0.)
      CALL HBOOK1(170,'RAW MYLAR SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(180,'RAW MYLAR SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(190,'SPECTRUM HITS',100,0.,1000.,0.)
      CALL HBOOK1(195,'PHOTON ELECTRONS',36,0.,360.,0.)
      CALL HBOOK1(200,'BACKSCATTER ENERGY SPECTRUM',100,0.,1000.,0.)
      CALL HBOOK2(210,'BACKSCATTER ENERGY VS ANGLE SPECTRUM',100,0.,
     1                1000.,100,0.,1.,0.)
      CALL HBOOK1(220,'ANGLUAR DISTRIBUTION OF BACKSCATTERS',100,0.,
     1                1.,0.)
      CALL HBOOK1(225,'COMPUTER TIME DISTRIBUTION',100,0.,3.,0.)

      RETURN
      END
C----------------------------------------------------------------------C     
      SUBROUTINE FILLBETATREE(EFOILE,EFOILW,DTYPE)
CC---------------------------------------------------------------------C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'ucnapenmain.h'
      PARAMETER(NB=5000)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/CNT1/TDEBO(NB),TDEBO2(NB),DEBO(NB)
      COMMON/CSOUR0/CTHL,DCTH,PHIL,DPHI,KPARP,LSCONE,LGPOL,LPSF,JOBEND
     1 ,NSIMTYPE,PTYPE
      EXTERNAL POSIT,DETREP
C----------------------------------------------------------------------C
      IF(DTYPE.EQ.1)THEN
        dtype = 1
      ENDIF
         
          ! SCINTILLATOR      
           DECS(14) = REAL(DEBO(398))
           DECS(13) = REAL(DEBO(416))
          ! DEAD LAYERS
           DECS(50) = REAL(DEBO(429))
           DECS(51) = REAL(DEBO(411))
          ! DECS(50) = REAL(DEBO(418))
           DECS(18) = REAL(DEBO(396)) ! East Detector Foils
           DECS(49) = REAL(DEBO(397)) 
           DECS(17) = REAL(DEBO(414)) ! West Detector Foils 
           DECS(48) = REAL(DEBO(415))
           DECS(21) = REAL(DEBO(434)+DEBO(432))
           DECS(22) = REAL(DEBO(433)+DEBO(435))
           DECS(23) = REAL(DEBO(431))  ! Decay Tube
           DECS(24) = REAL(DEBO(430))
           DECS(25) = REAL(POSIT(ECX))
           DECS(26) = REAL(POSIT(ECY))
           DECS(27) = REAL(POSIT(WCX))
           DECS(28) = REAL(POSIT(WCY))
           DECS(29) = REAL(TRGWEST(1))
           DECS(30) = REAL(TRGEAST(1))
           DECS(31) = REAL(DETREP(DEBO(416)))
           DECS(32) = REAL(DETREP(DEBO(398)))
           DECS(77) = REAL(DEBO(434))

        DO I = 1,14
            IF(ECX(I).GT.0)DECS(33)=REAL(DECS(33)+1.)
            IF(ECY(I).GT.0)DECS(34)=REAL(DECS(34)+1.)
            IF(WCX(I).GT.0)DECS(35)=REAL(DECS(35)+1.)
            IF(WCY(I).GT.0)DECS(36)=REAL(DECS(36)+1.)
         ENDDO
      
         DECS(37) = REAL(TRACKTIME)
         DECS(38) = REAL(abor)
         DECS(39) = REAL(ae)
         DECS(40) = REAL(AX)
         DECS(41) = REAL(AY)
         DECS(42) = REAL(AZ)
         DECS(43) = REAL(AU)
         DECS(44) = REAL(AV)
         DECS(45) = REAL(AW)
      
         PCONVERT = 420.5/1.0E6

         DECS(52) = REAL(DETREP(PHTE))
         DECS(53) = REAL(DETREP(PHTW))
         DECS(54) = REAL(PHTE)
         DECS(55) = REAL(PHTW)
         DECS(56) = REAL(EPr)
         DECS(57) = REAL(XPr)
         DECS(58) = REAL(YPr)
         DECS(59) = REAL(ZPr)
         DECS(60) = REAL(UPr)
         DECS(61) = REAL(VPr)
         DECS(62) = REAL(WPr)
         DECS(63) = REAL(PTOF)
         DECS(64) = REAL(PTYPE)
         DECS(65) = REAL(DEBO(437)+DEBO(438))
         DECS(70) = REAL(DEBO(440))

         DO I = 1,64
             DECS(71) = DECS(71) + REAL(DEBO(I))
             DECS(72) = DECS(72) + REAL(DEBO(I+65))
             DECS(73) = DECS(73) + REAL(DEBO(I+132))
             DECS(74) = DECS(74) + REAL(DEBO(I+197))
             DECS(75) = DECS(75) + REAL(DEBO(I+262))
             DECS(76) = DECS(76) + REAL(DEBO(I+328))
        ENDDO 
        DO I = 1,12
             DECS(77+I) = REAL(COSTHETA(I))
        ENDDO 
   
        CALL HFN(34,DECS) ! FILL EVENT RECORD
        IF(DEBO(398).GT.0)CALL HFILL(150,REAL(DEBO(398))/1E3,0,1.)
        IF(DEBO(416).GT.0)CALL HFILL(160,REAL(DEBO(416))/1E3,0,1.)
        CALL HFILL(225,REAL(TRACKTIME),0.,1.)
        
      RETURN
      END
C======================================================================C
C======================================================================c
      SUBROUTINE INITIALIZE_EVENT(NPAR)
c----------------------------------------------------------------------c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(NB=5000)
      INCLUDE 'ucnapenmain.h'
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/CNT1/TDEBO(NB),TDEBO2(NB),DEBO(NB)
c----------------------------------------------------------------------c
      DECS(1) = REAL(NPAR) ! Filling the initial event ntuple
      DECS(2) = REAL(KPAR) ! Particle number,type,position,energy,and
      DECS(3) = REAL(E)    ! momentum direction are stored in this
      DECS(4) = REAL(X)    ! ntuple.
      DECS(5) = REAL(Y)
      DECS(6) = REAL(Z)
      DECS(7) = REAL(U)
      DECS(8) = REAL(V)
      DECS(9) = REAL(W)
  
      DECS(10) = 0.       ! BACKSCATTERED ENERGY
      DECS(11) = 0.       ! BACKSCATTERED POLAR COSINE
      DECS(12) = 0.       ! LOGICAL SWITCH 1 IF MIRRORED
      
      DO I = 13,NENTRIES  ! CLEAR DATA ARRAY
         DECS(I)=0.
      ENDDO
      DO I = 1,170        ! CLEAR CATHODE CHANNELS
         ECX(I) = 0.
         ECY(I) = 0.
         WCX(I) = 0.
         WCY(I) = 0.
      ENDDO
c    
      TIME  = 0.0D0
      PHTE  = 0.
      PHTW  = 0.
      PHTEN = 0.
      PHTWN = 0.
      PTOF  = 0.
      wi  = w ! INITIAL PITCH ANGLE
      wim = w
c
      DO I=1,12
           IF(I.LE.10)THEN
             TRGEAST(I)  = 0.00D0  ! Set trigger times to 0
             TRGWEST(I)  = 0.00D0
           ENDIF
           COSTHETA(I) = 0.0D0
      ENDDO
c
      DO I=1,NB
          DEBO(I)=0.0 ! Set Energy depostion to 0
      ENDDO

      ILB(1) = 1   ! ILB(1) = 1 TELLS PENELOPE THIS IS A PRIMARY
      ILB(5) = 0   ! ILB(1) > 1 ARE SECONDARIES, TERTIARIES...
      
      IBEE = 0
      IBEW = 0 
      IDCE = 0 
      IDCW = 0 
      IMYFE= 0
      IMYBE= 0
      IMYFW= 0
      IMYBW= 0
      IDDE = 0
      IDDW = 0      

      RETURN
      END
C---------------------------------------------------------------------------C
