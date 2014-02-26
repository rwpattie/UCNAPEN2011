c-----------------------------------------------------------------------------C
      subroutine xenon_125
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      gweight = 99.6859
      eweight = 52.9057

      xprob = (gweight+eweight)*rand(1.d0)

      if(xprob.lt.6.78)then
        Kpar  = 2
        E = 54.968e3
      else if(xprob.ge.6.78.and.xprob.lt.6.898)then
        Kpar  = 2
        E = 74.875e3
      else if(xprob.ge.6.898.and.xprob.lt.7.377)then
        Kpar  = 2
        E = 113.551e3
      else if(xprob.ge.7.377.and.xprob.lt.61.177)then
        Kpar  = 2
        E = 188.418e3
      else if(xprob.ge.61.177.and.xprob.lt.61.252)then
        Kpar  = 2
        E = 210.418e3
      else if(xprob.ge.61.252.and.xprob.lt.91.252)then
        Kpar  = 2
        E = 243.378e3
      else if(xprob.ge.91.252.and.xprob.lt.91.423)then
        Kpar  = 2
        E = 372.081e3
      else if(xprob.ge.91.423.and.xprob.lt.96.093)then
        Kpar  = 2
        E = 453.796e3
      else if(xprob.ge.96.093.and.xprob.lt.96.1237)then
        Kpar  = 2
        E = 553.690e3
      else if(xprob.ge.96.1237.and.xprob.lt.96.2417)then
        Kpar  = 2
        E = 635.382e3
      else if(xprob.ge.96.2417.and.xprob.lt.96.3547)then
        Kpar  = 2
        E = 636.110e3
      else if(xprob.ge.96.3547.and.xprob.lt.96.4096)then
        Kpar  = 2
        E = 727.096e3
      else if(xprob.ge.96.4096.and.xprob.lt.97.5196)then
        Kpar  = 2
        E = 846.511e3
      else if(xprob.ge.97.5196.and.xprob.lt.98.0976)then
        Kpar  = 2
        E = 901.510e3
      else if(xprob.ge.98.0976.and.xprob.lt.98.2486)then
        Kpar  = 2
        E = 937.494e3
      else if(xprob.ge.98.2486.and.xprob.lt.98.3503)then
        Kpar  = 2
        E = 992.430e3
      else if(xprob.ge.98.3503.and.xprob.lt.98.5113)then
        Kpar  = 2
        E = 1007.431e3   
      else if(xprob.ge.98.5113.and.xprob.lt.98.5726)then
        Kpar  = 2
        E = 1075.540e3
      else if(xprob.ge.98.5726.and.xprob.lt.98.6377)then
        Kpar  = 2
        E = 1089.860e3
      else if(xprob.ge.98.6377.and.xprob.lt.98.9367)then
        Kpar  = 2
        E = 1138.23e3
      else if(xprob.ge.98.9367.and.xprob.lt.99.6197)then
        Kpar  = 2
        E = 1180.838e3
      else if(xprob.ge.99.6197.and.xprob.lt.99.6859)then
        Kpar  = 2
        E = 1193.230e3
      else if(xprob.ge.gweight.and.xprob.lt.24.4+gweight)then
        Kpar  = 1
        E = 21.799e3
      else if(xprob.ge.24.4+gweight.and.xprob.lt.38.4+gweight)then
        Kpar  = 1
        E = 23.600e3
      else if(xprob.ge.38.4+gweight.and.xprob.lt.38.74+gweight)then
        Kpar  = 1
        E = 41.706e3
      else if(xprob.ge.38.74+gweight.and.xprob.lt.41.99+gweight)then
        Kpar  = 1
        E = 49.780e3
      else if(xprob.ge.41.99+gweight.and.xprob.lt.42.644+gweight)then
        Kpar  = 1
        E = 53.896e3
      else if(xprob.ge.42.644+gweight.and.xprob.lt.42.776+gweight)then
        Kpar  = 1
        E = 54.782e3
      else if(xprob.ge.42.776+gweight.and.xprob.lt.42.966+gweight)then
        Kpar  = 1
        E = 69.687e3
      else if(xprob.ge.42.966+gweight.and.xprob.lt.43.007+gweight)then
        Kpar  = 1
        E = 73.803e3
      else if(xprob.ge.43.007+gweight.and.xprob.lt.43.224+gweight)then
        Kpar  = 1
        E = 80.382e3
      else if(xprob.ge.43.224+gweight.and.xprob.lt.49.424+gweight)then
        Kpar  = 1
        E = 155.249e3
      else if(xprob.ge.49.424+gweight.and.xprob.lt.50.298+gweight)then
        Kpar  = 1
        E = 183.230e3
      else if(xprob.ge.50.298+gweight.and.xprob.lt.50.475+gweight)then
        Kpar  = 1
        E = 187.346e3
      else if(xprob.ge.50.475+gweight.and.xprob.lt.52.425+gweight)then
        Kpar  = 1
        E = 210.209e3
      else if(xprob.ge.52.425+gweight.and.xprob.lt.52.780+gweight)then
        Kpar  = 1
        E = 238.190e3
      else if(xprob.ge.52.780+gweight.and.xprob.lt.52.8531+gweight)then
        Kpar  = 1
        E = 242.306e3
      else if(xprob.ge.52.8531+gweight)then
        Kpar  = 1
        E = 420.627e3
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
