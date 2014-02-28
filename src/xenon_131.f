c-----------------------------------------------------------------------------C
      subroutine xenon_131
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      cnvprob = 107.327*rand(1.d0)

      if (cnvprob.lt.6.9)then
        Kpar  = 1
        E = 24.600e3
      else if(cnvprob.ge.6.9.and.cnvprob.lt.68.5)then
        Kpar  = 1
        E = 129.369e3
      else if(cnvprob.ge.68.5.and.cnvprob.lt.82.4104)then
        Kpar  = 1
        E = 158.480e3
      else if(cnvprob.ge.82.4104.and.cnvprob.lt.85.6072)then
        Kpar  = 1
        E = 158.830e3
      else if(cnvprob.ge.85.6072.and.cnvprob.lt.97.3)then
        Kpar  = 1
        E = 159.150e3
      else if(cnvprob.ge.97.3.and.cnvprob.lt.103.89)then
        Kpar  = 1
        E = 162.788e3
      else if(cnvprob.ge.103.89.and.cnvprob.lt.105.23)then
        Kpar  = 1
        E = 162.720e3
      else if(cnvprob.ge.105.23.and.cnvprob.lt.105.377)then
        Kpar  = 1
        E = 163.914e3
      else if(cnvprob.ge.105.377)then
        Kpar  = 2
        E = 163.930e3
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
