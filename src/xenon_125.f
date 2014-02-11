c-----------------------------------------------------------------------------C
      subroutine xenon_125
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      gweight = 96.36
      eweight = 10.324

      xprob = (gweight+eweight)*rand(1.d0)

      if(xprob.lt.6.78)then
        Kpar  = 2
        E = 54.968e3
      else if(xprob.ge.6.78.and.xprob.lt.60.58)then
        Kpar  = 2
        E = 188.418e3
      else if(xprob.ge.60.58.and.xprob.lt.90.58)then
        Kpar  = 2
        E = 243.378e3
      else if(xprob.ge.90.58.and.xprob.lt.95.25)then
        Kpar  = 2
        E = 453.796e3
      else if(xprob.ge.95.25.and.xprob.lt.96.36)then
        Kpar  = 2
        E = 846.511e3
      else if(xprob.ge.gweight.and.xprob.lt.3.25+gweight)then
        Kpar  = 1
        E = 49.780e3
      else if(xprob.ge.3.25+gweight.and.xprob.lt.9.45+gweight)then
        Kpar  = 1
        E = 155.249e3
      else if(xprob.ge.9.45+gweight)then
        Kpar  = 1
        E = 183.230e3
      endif

      z     = 220.0*(1.0 - 2.0*rand(1.d0))
1200  continue
      x     = 6.23189*(1.0 - 2.0*rand(1.d0))
      y     = 6.23189*(1.0 - 2.0*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 6.23189) goto 1200 ! IR for 2011/2012
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta

      return
      end
