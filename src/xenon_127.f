c-----------------------------------------------------------------------------C
      subroutine xenon_127
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      gweight = 117.2643
      eweight = 17.298

      crudprob = (gweight+eweight)*rand(1.d0)

      if(crudprob.lt.1.24)then
        Kpar  = 2
        E = 57.61e3
      else if(crudprob.ge.1.24.and.crudprob.lt.5.55)then
        Kpar  = 2
        E = 145.252e3
      else if(crudprob.ge.5.55.and.crudprob.lt.31.25)then
        Kpar  = 2
        E = 172.132e3
      else if(crudprob.ge.31.25.and.crudprob.lt.99.95)then
        Kpar  = 2
        E = 202.860e3
      else if(crudprob.ge.99.95.and.crudprob.lt.117.25)then
        Kpar  = 2
        E = 374.991e3
      else if(crudprob.ge.117.25.and.crudprob.lt.117.2643)then
        Kpar  = 2
        E = 618.400e3
      else if(crudprob.ge.gweight.and.crudprob.lt.3.91+gweight)then
        Kpar  = 1
        E = 24.441e3
      else if(crudprob.ge.3.91+gweight.and.crudprob.lt.5.45+gweight)then
        Kpar  = 1
        E = 112.083e3
      else if(crudprob.ge.5.45+gweight.and.crudprob.lt.9.10+gweight)then
        Kpar  = 1
        E = 138.963e3
      else if(crudprob.ge.9.10+gweight.and.crudprob.lt.15.73
     1       +gweight)then
        Kpar  = 1
        E = 169.691e3
      else if(crudprob.ge.15.73+gweight.and.crudprob.lt.16.71
     1       +gweight)then
        Kpar  = 1
        E = 197.672e3
      else if(crudprob.ge.16.71+gweight.and.crudprob.lt.17.009
     1       +gweight)then
        Kpar  = 1
        E = 201.788e3
      else if(crudprob.ge.17.009+gweight)then
        Kpar  = 1
        E = 341.822e3
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
