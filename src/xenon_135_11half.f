c-----------------------------------------------------------------------------C
      subroutine xenon_135_11half
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      gweight = 80.4
      eweight = 18.879

      yprob = (gweight+eweight)*rand(1.d0)

      if(yprob.le.80.4)then
        Kpar  = 2
        E = 526.561e3
      else if(yprob.gt.gweight.and.yprob.le.15.34+gweight)then
        Kpar  = 1
        E = 492.000e3
      else if(yprob.gt.15.34+gweight.and.yprob.le.18.26+gweight)then
        Kpar  = 1
        E = 521.108e3
      else if(yprob.gt.18.26+gweight)then
        Kpar  = 1
        E = 525.419e3
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
