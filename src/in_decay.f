      subroutine in_decay 
c----------------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand,energy6
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from In114m
c    taken from two sources:
c
c    1) C. Wrede et al., NIM B 269, 1113 (2011)
c    2) http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=114IN&unc=nds
c----------------------------------------------------------------------------c

      gweight = 21.96
      eweight = 175.682

      x = (gweight+eweight)*rand(1.d0)
   
        if(x.le.40.1)then
          kpar =1 
          E = 162.33e3
        else if(x.gt.40.1.and.x.le.72.0)then
          kpar = 1
          E = 186.03e3
        else if(x.gt.72.0.and.x.le.78.71)then
          kpar = 1
          E = 189.44e3
        else if(x.gt.78.71.and.x.le.80.061)then
          kpar = 1 
          E =190.15e3
        else if(x.gt.80.061.and.x.le.175.682)then
          kpar = 1
          qend =198.86e4
          E = energy6(qend)
        else if(x.gt.eweight.and.x.le.eweight+15.56)then
          kpar = 2
          E = 190.27e3
        else if(x.gt.15.56+eweight.and.x.le.18.76+eweight)then
          kpar =2 
          E = 558.43e3
        else if(x.gt.18.76+eweight.and.x.le.21.96+eweight)then
          kpar =2 
          E = 725.24e3
        endif

      phi   = 2.*pi*rand(1.d0)
      costh = 1. - 2.*rand(1.d0)
      theta = dacos(costh)
 
      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)      

1000  continue
      y = (0.15 - 0.3*rand(1.d0))
      x = (0.15 - 0.3*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 0.15) goto 1000      
      z = 1.0

      return
      end

c-----------------------------------------------------------------------------C
      double precision function energy6(qend)
      implicit double precision(A-H,J-M,O-Z), integer*4(i,n)
      parameter(emass=510.998928d3)
      parameter(pi   =3.141592654d0)
      parameter(alphainv=1.37036d2)
      parameter(zed  =4.90d1)
      parameter(lambda=3.86159268d2)
      parameter(B    =1.0148d0)
      parameter(rn   =4.6056d0)
      parameter(ga=1.7765d0)
      common/rseed/iseed1,iseed2
      external rand
c----------------------------------------------------------------------------c
c  Fermi factor from D.H. Wilkinson, Nucl. Phys. A 377, 474 (1982).
c  Nucl. radii from I. Angeli, K. P. Marinova, Atomic Data and Nuclear 
c  Tables 99, 69 (2013).
c  F0L0 and L0 tables from Behrens and Janecke, Numerical Tables for 
c                          Beta Decay and Electron Capture (1969).
c----------------------------------------------------------------------------c

      EO = qend/emass + 1
50    E=(EO-1.0)*RAND(1.D0)
      if(E.lt.0.03) goto 50
      Y=170.0*RAND(1.D0)
      R=rn/lambda
      C=zed/alphainv
      A=B-0.9911d0 ! Endpt. of In114 at p = 4.9
      S=A/EO
      P=DSQRT((E+1)**2-1)
      CALL ZGMFN(DSQRT(1.0d0-C**2),(C*(E+1))/P,1,Q,U)
      FERMI=2*(1.0d0+DSQRT(1.0d0-C**2))*(1/((2*P*R)**(2*
     1     (1-DSQRT(1.0d0-C**2)))))*DEXP((PI*C*(E+1))/P)*
     1     (DSQRT(Q**2+U**2)/ga)**2*(B-S*E) ! D.H. Wilkinson (also Behrens/Janecke)
      if (Q.lt.0.0) goto 50
      W=FERMI*P*(EO-(E+1))**2*(E+1)
      if (W.lt.Y) goto 50
      !print*,'FO = ', FO
      E = E*emass
      energy6 = E
c 
      return
      end
