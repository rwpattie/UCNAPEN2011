      PARAMETER(EMASS=510998.0,QE=-1.0,NENTRIES=92) 
      PARAMETER(WIREDIV = 170)
      PARAMETER(NCDENTRIES=18)
      PARAMETER(NGENTRIES=12)
c      PARAMETER(DLAYER=6D-4)               ! DEAD LAYER THICKNESS
      PARAMETER(NWPAWC=500000)              ! PARAMETER TO MAKE NTUPLES BIGGER
      PARAMETER(DTI = 1E-11)                ! INITIAL TIME STEP GUESS
      CHARACTER*8 DTAGS(nentries)           ! NTUPLE VARIABLE ARRAY
      CHARACTER*8 CDTAGS(NCDENTRIES)
      CHARACTER*8 GETAGS(NGENTRIES)
      REAL DECAYPAR
      DIMENSION DECAYPAR(12)
      DIMENSION paw(nwpawc)
      
      DATA DTAGS/'NPAR','KPAR','E','X','Y','Z','U','V','W','Eb','Wb',
     1 'Mt','Epe','Epw','Egea','Egedf','Eme1','Emw1',
     1 'Egwa','Egwdf','Efle','Eflw','Etube','Eal',
     1 'eposx','eposy','wposx','wposy','tw','te','Eper','Epwr','EMX',
     1 'EMY','WMX','WMY','TIME','abort','ae','ax','ay','az','au',
     1 'av','aw','Egedb','Egwdb','Eme2','Emw2','Eped','Epwd','Phe',
     1 'Phw','nPhe','nPhw','EP','XP','YP','ZP','UP','VP','WP','PTOF',
     1 'Ptype','Eholder','Exsci','Eysci','Wxsci','Wysci','BckV','Eean',
     1 'Eec1','Eec2','Ewan','Ewc1','Ewc2','Ecol','W1','W2','W3','W4',
     1 'W5','W6','W7','W8','W9','W10','W11','W12','Eapd','Epadh','tapd'/

      DATA CDTAGS/'NPAR','KPAR','E','X','Y','Z','U','V','W',
     1 'ESI','EDEAD','EBACKING','ERING','EFOIL','AX','AY','AZ','AMAT'/

      DATA GETAGS/'E','U','V','W','EP','UP','VP','WP','NE','NU','NV',
     1 'NW'/
     
      COMMON/PAWC/PAW,H(NWPAWC)
      COMMON/QUEST/IQUEST(100)
      REAL DECS(NENTRIES)
      REAL COSTHETA(12)
      COMMON/CUCAP/CLINES(4,424),CLINES2(4,325)
      COMMON/PROTON/EPr,XPr,YPr,ZPr,UPr,VPr,WPr,PTOF
      COMMON/NPEBL/PHTE,PHTW,PHTEN,PHTWN,wi,wim,TIME
      COMMON/TRGT/TRGE(10),TRGW(10)
      COMMON/INIT/NPAR0,KPAR0,Ei,X0,Y0,Z0,U0,V0,W0,PTYPE0
      COMMON/EDEP/EPE,EPW,EPER,EPED,EPWD,EPWR,EGEA,EGWA,
     1 EGEDF,EGEDB,EGWDF,EGWDB,EFLE,EFLW,ETUBE,EAL,EME1,EME2,
     1 EMW1,EMW2,PHW,PHEN,PHE,PHWN,EHOLDER,EEAN,EEC1,EEC2,EWAN,
     1 EWC1,EWC2,ECOL,EAPD,EAPDH
      COMMON/COS/W1(12)
      COMMON/MULTI/EMX,EMY,WMX,WMY,EPOSX,EPOSY,WPOSX,WPOSY,
     1 EXSCI,EYSCI,WXSCI,WYSCI,TAPD
      COMMON/HBOOKU/DECS,COSTHETA
      COMMON/INTACT/IBEE,IBEW,IDCE,IDCW,IMYFE,IMYFW,IMYBE,IMYBW,IDDE,
     1 IDDW
      COMMON/MWPC/ECX(170),ECY(170),WCX(170),WCY(170)
      COMMON/ABORT/TRACKTIME,ABOR,AE,AX,AY,AZ,AU,AV,AW
