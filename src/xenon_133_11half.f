c-----------------------------------------------------------------------------C
      subroutine xenon_133_11half
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      convxprob = 89.56*rand(1.d0)

      if(convxprob.lt.62.90)then
        E = 198.66e3
      else if(convxprob.ge.62.90.and.convxprob.lt.83.90)then
        E = 227.768e3
      else if(convxprob.ge.83.90.and.convxprob.lt.88.60)then
        E = 232.079e3
      else if(convxprob.ge.88.60)then
        E = 233.013e3
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
      Kpar  = 1

      return
      end
