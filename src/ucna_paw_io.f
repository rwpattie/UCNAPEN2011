c----------------------------------------------------------------------c
      SUBROUTINE BUILDHBOOK()
c----------------------------------------------------------------------c
      INCLUDE 'ucnapenmain.h'
c----------------------------------------------------------------------c
      CALL HLIMIT(NWPAWC) ! A FEW PAW FORUMS SUGGEST THIS FOR MAKING
      IQUEST(10) = 65000  ! BIGGER NTUPLES. KINDA WORKS BUT STILL 
                          ! LIMITED.
      
      CALL HROPEN(1,'EVENT','test.hbook','NQE',1024,IERR)
! c      if(NTYPE.NE.10)THEN
      CALL HBOOKN(34,'DECAYS',nentries,'//EVENT',1024,DTAGS) 
c      ELSE 
c        CALL HBOOKN(34,'DECAYS',ngentries,'//EVENT',1024,GETAGS)
c      ENDIF
C      ELSE IF(DTYPE.EQ.2)THEN
c        CALL HBOOKN(34,'DECAYS',nCDentries,'//EVENT',1024,CDTAGS)
c      ENDIF
      CALL HBOOK1(100,'INTIAL BETA DECAY SPECTRUM',300,0.,3000.,0.)
      CALL HBOOK1(110,'THROWN ASYMMETRY',100,0.,1000.,0.)
      CALL HBOOK1(120,'MEASURED ASYMMETRY (RAW)',100,0.,1000.,0.)
      CALL HBOOK1(130,'RAW MWPC SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(140,'RAW MWCP SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(150,'RAW PMT SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(160,'RAW PMT SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(170,'RAW MYLAR SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(180,'RAW MYLAR SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(190,'SPECTRUM HITS',100,0.,1000.,0.)
      CALL HBOOK1(195,'PHOTON ELECTRONS',36,0.,360.,0.)
      CALL HBOOK1(200,'BACKSCATTER ENERGY SPECTRUM',100,0.,1000.,0.)
      CALL HBOOK2(210,'BACKSCATTER ENERGY VS ANGLE SPECTRUM',100,0.,
     1                1000.,100,0.,1.,0.)
      CALL HBOOK1(220,'ANGLUAR DISTRIBUTION OF BACKSCATTERS',100,0.,
     1                1.,0.)

      RETURN
      END
C----------------------------------------------------------------------C     
      SUBROUTINE FILLBETATREE(DECS,EDEP,EFOILE,EFOILW,EAL,PTYPE
     1                       ,NENTRIES,DTYPE,COSTHETA)
CC---------------------------------------------------------------------C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/TRACK/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/MWPC/ECX(170),ECY(170),WCX(170),WCY(170)
      COMMON/NPEBL/PHTE,PHTW,PHTEN,PHTWN,PTOF,wi,wim,TIME
      COMMON/TRGT/TRGEAST(10),TRGWEST(10)
      COMMON/ABORT/TRACKTIME,ABOR,AX,AY,AZ,AU,AV,AW
      COMMON/PROTON/EP,XP,YP,ZP,UP,VP,WP
      REAL DECS,COSTHETA
      DIMENSION EDEP(500)
      DIMENSION DECS(nentries)
      DIMENSION COSTHETA(12)
      EXTERNAL POSIT,DETREP
C----------------------------------------------------------------------C

      ! SCINTILLATOR
           DECS(14) = REAL(EDEP(398))
           DECS(13) = REAL(EDEP(416))
        ! DEAD LAYERS
           DECS(50) = REAL(EDEP(429))
           DECS(51) = REAL(EDEP(411))
!           DECS(50) = REAL(EDEP(418))
           DECS(18) = REAL(EDEP(396)) ! East Detector Foils
           DECS(49) = REAL(EDEP(397)) 
           DECS(17) = REAL(EDEP(414)) ! West Detector Foils 
           DECS(48) = REAL(EDEP(415))
           DECS(21) = REAL(efoile)
           DECS(22) = REAL(efoilw)
           DECS(23) = REAL(EDEP(431))  ! Decay Tube
           DECS(24) = REAL(EDEP(430))
           DECS(25) = REAL(POSIT(ECX))
           DECS(26) = REAL(POSIT(ECY))
           DECS(27) = REAL(POSIT(WCX))
           DECS(28) = REAL(POSIT(WCY))
           DECS(29) = REAL(TRGWEST(1))
           DECS(30) = REAL(TRGEAST(1))
           DECS(31) = REAL(DETREP(EDEP(416)))
           DECS(32) = REAL(DETREP(EDEP(398)))
           DECS(77) = REAL(EDEP(434))
         
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
         DECS(56) = REAL(EP)
         DECS(57) = REAL(XP)
         DECS(58) = REAL(YP)
         DECS(59) = REAL(ZP)
         DECS(60) = REAL(UP)
         DECS(61) = REAL(VP)
         DECS(62) = REAL(WP)
         DECS(63) = REAL(PTOF)
         DECS(64) = REAL(PTYPE)
         DECS(65) = REAL(EDEP(435))
         DECS(70) = REAL(EDEP(440))

         DO I = 1,64
             DECS(71) = DECS(71) + REAL(EDEP(I))
             DECS(72) = DECS(72) + REAL(EDEP(I+65))
             DECS(73) = DECS(73) + REAL(EDEP(I+132))
             DECS(74) = DECS(74) + REAL(EDEP(I+197))
             DECS(75) = DECS(75) + REAL(EDEP(I+262))
             DECS(76) = DECS(76) + REAL(EDEP(I+328))
        ENDDO 
        DO I = 1,12
             DECS(77+I) = COSTHETA(I)
        ENDDO 
              
      RETURN
      END
C======================================================================C