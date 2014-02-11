c-----------------------------------------------------------------------------C
      subroutine xenon_131
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      cnvprob = 97.137*rand(1.d0)

      if(cnvprob.lt.61.6)then
        E = 129.369e3
      else if(cnvprob.ge.61.6.and.cnvprob.lt.90.4)then
        E = 158.477e3
      else if(cnvprob.ge.90.4.and.cnvprob.lt.96.99)then
        E = 162.788e3
      else if(cnvprob.ge.96.99)then
        E = 163.914e3
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
