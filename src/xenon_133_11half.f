c-----------------------------------------------------------------------------C
      subroutine xenon_133_11half
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      convxprob = 106.787*rand(1.d0)

      if(convxprob.lt.7.00)then
        Kpar  = 1
        E = 24.60e3
      else if(convxprob.ge.7.00.and.convxprob.lt.69.90)then
        Kpar  = 1
        E = 198.66e3
      else if(convxprob.ge.69.90.and.convxprob.lt.90.90)then
        Kpar  = 1
        E = 227.768e3
      else if(convxprob.ge.90.90.and.convxprob.lt.95.60)then
        Kpar  = 1
        E = 232.079e3
      else if(convxprob.ge.95.60.and.convxprob.lt.96.56)then
        Kpar  = 1
        E = 233.013e3
      else if(convxprob.ge.96.56.and.convxprob.lt.96.667)then
        Kpar  = 1
        E = 233.205e3
      else if(convxprob.ge.96.667)then
        Kpar  = 2
        E = 233.221e3
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
