      SUBROUTINE PROTONS(DECAYPAR)
c---------------------------------------------------------------------c
c     This program generates an event distribution from neutron beta
c     decay using functions from F. Gluck, "Measurable Distributions of 
c     Unpolarized Neutron Beta Decay", Phys. Rev. D., V 47 N 7 2840 (1993).
c---------------------------------------------------------------------c
      IMPLICIT DOUBLE PRECISION(A-H,M-Z), INTEGER*4(I-L)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      external energy2,aprob2
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/PROTON/EPr,XPr,YPr,ZPr,UPr,VPr,WPr
      COMMON/RSEED/ISEED1,ISEED2
      DIMENSION DECAYPAR(12)
      REAL DECAYPAR

      PARAMETER( MN = 939.5656D6, MP = 938.2723D6)
      PARAMETER( ELAMBDA = -1.267,alpha =1./137.)
      
      IGLUCHECK = 2
c      PRINT*,'SET NEUTRON POLARIZATION AND EXPERIMENT: (POL,AEXP,BEXP)'
c      READ*,POL,AEXP,BEXP
      pol = 1.0
      aexp = 1
      bexp = 1
      REJECT = 100
c     DEFINE MASSES IN TERMS OF THE ELECTRON REST MASS
      MN1 = MN/ME 
      MP1 = MP/ME
      ME1 = 1.
      DELTA = MN1 - MP1
C     ELECTRON CRITICAL ENERGY
      XCRIT = 0.5*(DELTA + ME1*ME1/DELTA)
C     MAXIMUM ELECTRON ENERGY
 
      EEMAX = DELTA - (DELTA*DELTA-1.)/(2.*MN1)

C     CALCULATE STANDARD MODEL VALUES FOR a,b,A,B
C
      A = -2.*elambda*(1.+elambda)/(1.+3.*elambda**2)
      B = -2.*elambda*(1.-elambda)/(1.+3.*elambda**2)
      alit = (1.-elambda**2)/(1.+3.*elambda**2)
      el1 = Elambda
      BLIT = 0
c
C     SELECT ELECTRON ENERGY
C
 10   CONTINUE
      EE = ENERGY2()+1.
      PE = SQRT(EE*EE-1.)
C
C     SELECT DETERMINE NEUTRON SPIN DIRECTION ANGLE
C
      check15 = 1.0*rand(1.d0)
C
      if(abs(pol).gt.check15.and..not.pol.eq.1.0)then
         ntheta = 0.
         NEPHI = 0.
      elseif(abs(pol).lt.check15.and..not.pol.eq.1.0)then
         ntheta = pi
         NEPHI = 0.
      endif
c
      if(pol.eq.0.)THEN 
           ntheta=pi*rand(1.d0)
           NEPHI =2*PI*RAND(1.d0)
      ENDIF
c
      if(pol.eq.1.0)then
        nephi = 0.0
        ntheta= 0.0
      endif
c
      NEU = COS(NEPHI)*SIN(NTHETA)
      NEV = SIN(NEPHI)*SIN(NTHETA)
      NEW = COS(NTHETA)
C
C     OFF SET ANGLES BY THE NEUTRON POLARIZATION
C
 9    continue
      THETAE = PI*RAND(1.d0)
      PHIE = 2.*pi*RAND(1.d0)
c
      if(pol.eq.1.0)then
        ff = aprob2(A,el1,thetae,Ee)
        gg = 1.08*rand(1.d0)
        if(gg.gt.ff) goto 10
      endif
C
C     DEFINE ELECTRON DIRECTION COSINES
C
      EU = COS(PHIE)*SIN(THETAE)
      EV = SIN(PHIE)*SIN(THETAE)
      EW = COS(THETAE)

11    continue

      THETAP = DACOS(2*RAND(1.d0)-1)
      PHIP = 2.*PI*RAND(1.d0)
C
C     DEFINE PROTON DIRECTION COSINES
C
      PU = COS(PHIP)*SIN(THETAP)
      PV = SIN(PHIP)*SIN(THETAP)
      PW = COS(THETAP)

      THETAEP = EU*PU+PV*EV+EW*PW
C
C     CHECK TO SEE IF THETAEP IS INSIDE THE ALLOW KINEMATIC REGION
C 
      IF(EE.GE.XCRIT)THEN
            THETAMAX = -SQRT(1.-((MP1/MN1)*(EEMAX-EE)/PE)**2)
            IF(THETAEP.GT.THETAMAX)GOTO 11
      ENDIF
     
C
C     EVAULATE DELTA FUNCTION TO FIND PROTON ENERGY
C
      XM = MN1 - EE
      X  = XM*XM-PE*PE+MP1*MP1
      Y  = PE*THETAEP
C 
C     SOLUTIONS TO THE QUADRATIC EQUATION FOR THE PROTON MOMENTA
C 
      aa = 4.*(xm*xm-y*y)
      bb = 4.*y*x
      hbig = 0.5*(mp1+mn1+1./(mp1+mn1))
      hlit = 4.*(mn1*mn1-mp1*mp1)*(ee-xcrit)*(hbig-ee)
      s = sqrt(bb*bb - 4.*aa*hlit)
c
      ppp = (-bb + s)/(2.*aa)
      ppm = (-bb - s)/(2.*aa)
c 
      RRP = SQRT(PPP*PPP+MP1*MP1)
      RRM = SQRT(PPM*PPM+MP1*MP1)
C
C     CHECK KINEMATIC LIMITS OF PROTON ENERGY
C
c
      EPMAX = SQRT((PE+(MN1*(DELTA-EE)/(MN1-EE-PE)))**2 + MP1**2)
c
      IF(EE.LT.XCRIT)THEN
           EPMIN=SQRT(((MN1*(DELTA-EE)/(MN1-EE-PE))-PE)**2 + MP1**2)
      ELSE
           EPMIN=SQRT((PE-(MN1*(DELTA-EE)/(MN1-EE-PE)))**2 + MP1**2)
      ENDIF
c
      if(ee.gt.xcrit)then
           if(rrp.lt.epmax.and.rrp.gt.epmin)jblah2=1
           if(rrm.lt.epmax.and.rrm.gt.epmin)jblah3=1
      endif

      IF(IGLUCHECK.EQ.2)THEN
C
C     DEFINE THE JACOBIAN FOR THE DELTA FUNCTION 
C
       FDP = ABS(PPP/RRP + (PPP+Y)/(MN1-EE-RRP))
       FDM = ABS(PPM/RRM + (PPM+Y)/(MN1-EE-RRM))
C
C     DEFINE THE MATRIX ELEMENT
C
       e1p = sqrt(ppp*ppp + pe*pe + 2*pe*ppp*thetaep)
       e1m = sqrt(ppm*ppm + pe*pe + 2*pe*ppm*thetaep)
c
      CHI1 =( 1. - ALIT*(PE*PE+Y*PPP)/(EE*E1P) + BLIT/EE 
     1  + POL*(A*PE/EE*EW-B*(PE*EW+PPP*PW)/E1P))
     1  * EE*PE*PPP*PPP*SIN(THETAp)*SIN(THETAe)
c
      CHI2 =( 1. - ALIT*(PE*PE+Y*PPM)/(EE*E1M) + BLIT/EE 
     1  + POL*(A*PE/EE*EW-B*(PE*EW+PPM*PW)/E1M))
     1  *EE*PE*PPM*PPM*SIN(THETAp)*SIN(THETAe)
C
C     DEFINE THE DIFFERENTIALS
C
      DIFFP = 1.!SQRT(1.-THETAEP**2)!EE*PE*PPP*PPP*SIN(THETAE)*SIN(THETAP)
      DIFFM = 1.!SQRT(1.-THETAEP**2)!EE*PE*PPM*PPM*SIN(THETAE)*SIN(THETAP)

      ELSEIF(IGLUCKCHECK.EQ.1)THEN
c
c     now for the gluck way
C 
          da = -8.*y*pe
          db = 4.*pe*x
c
         fdp=1./dabs((ppp/rrp)*((bb-s)/(2*aa*aa)*da+(bb-s)/(2.*aa*s)*db-
     1      hlit/(aa*s)*da))
c
         fdm=1./dabs((ppm/rrm)*((bb+s)/(2*aa*aa)*da-(bb+s)/(2.*aa*s)*db+
     1      hlit/(aa*s)*da))
cC
         el1 = mn1 - ee - rrp
         el2 = mn1 - ee - rrm
c
         fkap = 1.855
c
         elmax = delta - (delta*delta+1.)/(2.*mn1)
         emax  = mp1 + (delta*delta-1.)/(2.*mn1)
c
         dv1 = ee*(eemax-ee)+el1*(elmax-el1)-mp1*(emax-rrp)
         dv2 = ee*(eemax-ee)+el2*(elmax-el2)-mp1*(emax-rrm)
c
         da1 = ee*(eemax-ee)+el1*(elmax-el1)+mp1*(emax-rrp)
         da2 = ee*(eemax-ee)+el2*(elmax-el2)+mp1*(emax-rrm)
c
         di1 = 2.*(ee*(eemax-ee)-el1*(elmax-el1))
         di2 = 2.*(ee*(eemax-ee)-el2*(elmax-el2))
c
         chi1 = dv1+elambda**2*da1+elambda*(1.+2.*fkap)*di1
         chi2 = dv1+elambda**2*da1+elambda*(1.+2.*fkap)*di2
c
         DIFFP = 1.
         DIFFM = 1.
      ENDIF
C
C     SET INTEGRATION LIMITS
C
      IF(EE.GE.XCRIT)THEN
         H = 1.*RAND(1.d0)
         F1 = (CHI1*DIFFP/FDP)/(CHI1*DIFFP/FDP + CHI2*DIFFM/FDM)
         F2 = (CHI2*DIFFM/FDM)/(CHI1*DIFFP/FDP + CHI2*DIFFM/FDM)
           IF(H.GT.F2)THEN
                EP = RRp
                PP = PPp
           ELSE
                EP = RRm
                PP = PPm
           ENDIF
      ELSEIF(EE.LT.XCRIT)THEN
         EP = RRP
         PP = PPP
         CHI2 = 0.
      ENDIF
      IF(EP.GT.EPMAX.OR.EP.LT.EPMIN)GOTO 11
C
C     DO THE MONTE CARLO INTEGRATION 
C 
      beta = pe/ee !electron's beta
c
      betar = abs(beta-(1.0-beta**2)*pp/ep*thetaep) ! relative beta between proton and electron
c
      gam = 0.5772 ! a function not defined See Gluck
c
      FERMI = 1. + alpha*pi/betar + alpha**2*(11./4.- gam -
     1        log(2.*betar*ee*(0.01)/4.)+pi*pi/(3*betar**2)) 

      DGAMMA1 =FERMI*(CHI1*DIFFP/FDP+CHI2*DIFFM/FDM)
c
      if(ee.gt.xcrit)then 
         reject=700.
      elseif(ee.le.xcrit)then
         reject = 10.
      endif
      G = REJECT*RAND(1.d0)
c
      IF(G.GT.DGAMMA1)GOTO 11
C
C     DETERMINE THE NEUTRINO VARIALBES
C
      EN = MN1 - EE - EP
      NU = -(PE*EU+PP*PU)/EN
      NV = -(PE*EV+PP*PV)/EN
      NW = -(PE*EW+PP*PW)/EN
      THETANE = EU*NU+EV*NV+EW*NW
C     Set the electron energy and dirction
      E = (EE -1.) * ME
c
      RADIUS = 6.0
c      PPSI   = 2*PI*RAND(1.D0)
c
320   CONTINUE
      X = RADIUS*(1.0 - 2.0*RAND(1.0D0))
      Y = RADIUS*(1.0 - 2.0*RAND(1.0D0))
      IF(SQRT(X*X + Y*Y) .GT. RADIUS) GOTO 320
      Z = 1.5d2 - 3.0e2*RAND(1.D0)
c 
c     quick hack for unpolarized decay

      U = EU
      V = EV
      W = EW

c     set the proton energy and direction
      EP = EP*ME - MP
      XP = X
      YP = Y
      ZP = Z
      UP = PU
      VP = PV
      WP = PW

      DECAYPAR(1) = REAL(E)
      DECAYPAR(2) = REAL(U)
      DECAYPAR(3) = REAL(V)
      DECAYPAR(4) = REAL(W)
      DECAYPAR(5) = REAL(EP)
      DECAYPAR(6) = REAL(UP)
      DECAYPAR(7) = REAL(VP)
      DECAYPAR(8) = REAL(WP)
      DECAYPAR(9) = REAL(EN)
      DECAYPAR(10) = REAL(NU)
      DECAYPAR(11) = REAL(NV)
      DECAYPAR(12) = REAL(NW)
      
      RETURN
      END
c----------------------------------------------------------------c
      double precision function energy2()
      implicit double precision(a-h,J-M,o-z), integer*4(i,n)
      external rand
      parameter(E0=2.5309874e0, b=0.0)
c
50    E=(E0-1.0)*RAND(1.d0)
      y=1.80*RAND(1.d0)
      GAM=(-E/(DSQRT(E*E+2*E)*137.))
      FERMI=2.*3.14*GAM/(EXP(2.*3.14*GAM)-1.)
      f=FERMI*DSQRT(E**2+2*E)*(E0-(E+1))**2*(E+1)*(1+b*1/(e+1))
      if (f.lt.y) goto 50
      energy2 = E
c
      return
      end
c---------------------------------------------------------------------c
      double precision function aprob2(A,lambda,theta,E)
      implicit double precision(a-h,J-M,o-z), integer*4(i,n)
      PARAMETER (PI = 3.141592654d0)
      external rand
c----------------------------------------------------------------------c
      beta = dsqrt(E*E - 1.d0)/E
      lambda=-1.0*lambda
      W     = 2.529476
      alpha = 1./137.036
      MM    = 1837.0
      MU    = 3.71
c----------------------------------------------------------------------c
      AM = (LAMBDA+MU)/(LAMBDA*MM*(1.0-LAMBDA)*(1.0+3.*LAMBDA**2))
      A1 = LAMBDA**2+(2./3.)*LAMBDA-(1./3.)
      A2 = -LAMBDA**3-3.*LAMBDA**2-(5./3.)*LAMBDA+(1./3.)
      A3 = (1.0-LAMBDA)*2.*LAMBDA**2
c----------------------------------------------------------------------c
      Abb   = A*(1.0+AM*(A1*W+A2*E+A3/E))
      aprob2 = dsin(theta)*(1.0+Abb*beta*dcos(theta))
c----------------------------------------------------------------------c
      return
      end
