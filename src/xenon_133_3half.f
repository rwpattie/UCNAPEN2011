c-----------------------------------------------------------------------------C
      subroutine xenon_133_3half
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand,energy5
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      btxeprob = rand(1.d0)

      if(btxeprob.lt.0.986)then
        qend = 346.4e3
        rej = 3.2d0
        rn = 4.7831d0
        e = energy5(qend,rej,rn)
      else if(btxeprob.ge.0.986)then
        qend = 266.8e3
        rej = 1.9d0
        rn = 4.7831d0
        e = energy5(qend,rej,rn)
      endif

      z     = 220.0*(1.0 - 2.0*rand(1.d0))
1200  continue
      x     = 6.23189*(1.0 - 2.0*rand(1.d0))
      y     = 6.23189*(1.0 - 2.0*rand(1.d0))
      if(sqrt(x**2 + y**2) .ge. 6.23189) goto 1200 ! IR for 2011/2012
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      Kpar  = 1

      return
      end

c----------------------------------------------------------------------------C
      double precision function energy5(qend,rej,rn)
      implicit double precision(A-H,J-M,O-Z), integer*4(i,n)
      parameter(emass=510.998928d3)
      parameter(pi   =3.141592654d0)
      parameter(alphainv=1.37036d2)
      parameter(zed  =5.40d1)
      parameter(lambda=3.86159268d2)
      parameter(B    =1.0170d0)
      parameter(ga=1.73169d0)
      common/rseed/iseed1,iseed2
      external rand
c----------------------------------------------------------------------------c
c  Fermi factor from D.H. Wilkinson, Nucl. Phys. A 377, 474 (1982).
c  Nucl. radii from I. Angeli, K. P. Marinova, Atomic Data and Nuclear 
c  Tables 99, 69 (2013).  Linear extrapolation for Xe135 & appx. Xe133
c  F0L0 and L0 tables from Behrens and Janecke, Numerical Tables for 
c                          Beta Decay and Electron Capture (1969).
c----------------------------------------------------------------------------c

      EO = qend/emass + 1
50    E=(EO-1.0)*RAND(1.D0)
      if(E.lt.0.03) goto 50
      Y=rej*RAND(1.D0)
      R=rn/lambda
      C=zed/alphainv
      A=B-1.0052d0
      S=A/EO
      P=DSQRT((E+1)**2-1)
      CALL ZGMFN(DSQRT(1.0d0-C**2),(C*(E+1))/P,1,Q,U)
      FERMI=2*(1.0d0+DSQRT(1.0d0-C**2))*(1/((2*P*R)**(2*
     1     (1-DSQRT(1.0d0-C**2)))))*DEXP((PI*C*(E+1))/P)*
     1     (DSQRT(Q**2+U**2)/ga)**2*(B-S*E) ! D.H. Wilkinson (also Behrens/Janecke)
      if (Q.lt.0.0) goto 50
      W=FERMI*P*(EO-(E+1))**2*(E+1)
      if (W.lt.Y) goto 50
      !print*,'FERMI = ', FERMI
      E = E*emass
      energy5 = E
c 
      return
      end
