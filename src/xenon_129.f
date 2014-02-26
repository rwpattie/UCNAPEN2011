c-----------------------------------------------------------------------------C
      subroutine xenon_129
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      eweight = 201.745
      gweight = 12.09
      zprob = (eweight+gweight)*rand(1.d0)

      if(zprob.lt.78.85)then
        Kpar  = 1
        E = 5.017e3
      else if(zprob.ge.78.85.and.zprob.lt.94.85)then
        Kpar  = 1
        E = 24.600e3
      else if(zprob.ge.94.85.and.zprob.lt.105.58)then
        Kpar  = 1
        E = 34.125e3
      else if(zprob.ge.105.58.and.zprob.lt.107.745)then
        Kpar  = 1
        E = 38.436e3
      else if(zprob.ge.107.745.and.zprob.lt.171.745)then
        Kpar  = 1
        E = 162.00e3
      else if(zprob.ge.171.745.and.zprob.lt.196.245)then
        Kpar  = 1
        E = 191.11e3
      else if(zprob.ge.196.245.and.zprob.lt.201.745)then
        Kpar  = 1
        E = 195.42e3
      else if(zprob.ge.eweight.and.zprob.lt.7.50+eweight)then
        Kpar  = 2
        E = 39.578e3
      else if(zprob.ge.7.50+eweight)then
        Kpar  = 2
        E = 196.56e3
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

      return
      end
