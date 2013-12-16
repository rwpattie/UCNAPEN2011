c---------------------------------------------------------------------c
      subroutine cd_113m_decay
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      EXTERNAL rand,energy4
      include 'pmcomms.f'
c----------------------------------------------------------------------c
c    Endpoint and mean energy of Cd113m can be found here:
c
c    http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=113CD&unc=nds
c------------------------------------------------------------------------

      qend  = 585.7e3
      Kpar  = 1
      e     = energy4(qend)
      phi   = 2.*pi*rand(1.d0)
      costh = 1. - 2.*rand(1.d0)
      theta = dacos(costh)

      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)

1900  continue
      y = (0.15 - 0.3*rand(1.d0))
      x = (0.15 - 0.3*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 0.15) goto 1900
      z = 1.0

      return
      end

c-----------------------------------------------------------------------------C
      double precision function energy4(qend)
      implicit double precision(a-h,J-M,o-z), integer*4(i,n)
      parameter(emass=510.998928d3)
      parameter(pi   =3.141592654d0)
      parameter(alphainv=1.37036d2)
      parameter(zed  =4.80d1)
      parameter(lambda=3.86159268d2)
      parameter(b    =1.0144d0)
      parameter(rn   =4.6012d0)
      parameter(ga=1.78506d0)
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
      y=6.55*RAND(1.D0)
      R=rn/lambda
      C=zed/alphainv
      a=b-1.0086d0 ! Endpt. of Cd113m at p = 0.9
      s=a/EO
      p=DSQRT((E+1)**2-1)
      CALL ZGMFN(DSQRT(1.0d0-C**2),(C*(E+1))/p,1,Q,M)
      LO=b-s*E ! Linear extrapolation of Behrens and Janecke L0
      FO=2*(1.0d0+DSQRT(1.0d0-C**2))*(1/((2*p*R)**(2*
     1     (1-DSQRT(1.0d0-C**2)))))*DEXP((PI*C*(E+1))/p)*
     1     (DSQRT(Q**2+M**2)/ga)**2 ! D.H. Wilkinson (also Behrens/Janecke)
      FERMI=FO*LO 
      if (Q.lt.0.0.and.M.lt.0.0) goto 50
      j=Fermi*p*(EO-(E+1))**2*(E+1)
      if (j.lt.y) goto 50
      !print*,'FO = ', FO
      E = E*emass
      energy4 = E
c 
      return
      end
