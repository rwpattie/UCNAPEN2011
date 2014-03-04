c-----------------------------------------------------------------------------C
      subroutine xenon_125
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      gweight = 99.84
      eweight = 52.9851

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
      else if(xprob.ge.91.252.and.xprob.lt.91.2655)then
        Kpar  = 2
        E = 258.360e3
      else if(xprob.ge.91.2655.and.xprob.lt.91.2854)then
        Kpar  = 2
        E = 340.220e3
      else if(xprob.ge.91.2854.and.xprob.lt.91.4564)then
        Kpar  = 2
        E = 372.081e3
      else if(xprob.ge.91.4564.and.xprob.lt.96.1264)then
        Kpar  = 2
        E = 453.796e3
      else if(xprob.ge.96.1264.and.xprob.lt.96.1571)then
        Kpar  = 2
        E = 553.690e3
      else if(xprob.ge.96.1571.and.xprob.lt.96.2751)then
        Kpar  = 2
        E = 635.382e3
      else if(xprob.ge.96.2751.and.xprob.lt.96.3881)then
        Kpar  = 2
        E = 636.110e3
      else if(xprob.ge.96.3881.and.xprob.lt.96.4016)then
        Kpar  = 2
        E = 717.900e3
      else if(xprob.ge.96.4016.and.xprob.lt.96.4565)then
        Kpar  = 2
        E = 727.096e3
      else if(xprob.ge.96.4565.and.xprob.lt.96.464)then
        Kpar  = 2
        E = 764.170e3
      else if(xprob.ge.96.464.and.xprob.lt.96.471)then
        Kpar  = 2
        E = 809.180e3
      else if(xprob.ge.96.471.and.xprob.lt.96.4952)then
        Kpar  = 2
        E = 819.020e3
      else if(xprob.ge.96.4952.and.xprob.lt.97.6052)then
        Kpar  = 2
        E = 846.511e3
      else if(xprob.ge.97.6052.and.xprob.lt.97.6212)then
        Kpar  = 2
        E = 894.420e3
      else if(xprob.ge.97.6212.and.xprob.lt.98.1992)then
        Kpar  = 2
        E = 901.510e3
      else if(xprob.ge.98.1992.and.xprob.lt.98.3502)then
        Kpar  = 2
        E = 937.494e3
      else if(xprob.ge.98.3502.and.xprob.lt.98.4519)then
        Kpar  = 2
        E = 992.430e3
      else if(xprob.ge.98.4519.and.xprob.lt.98.6129)then
        Kpar  = 2
        E = 1007.431e3   
      else if(xprob.ge.98.6129.and.xprob.lt.98.6366)then
        Kpar  = 2
        E = 1020.550e3
      else if(xprob.ge.98.6366.and.xprob.lt.98.6544)then
        Kpar  = 2
        E = 1070.850e3
      else if(xprob.ge.98.6544.and.xprob.lt.98.7157)then
        Kpar  = 2
        E = 1075.540e3
      else if(xprob.ge.98.7157.and.xprob.lt.98.7808)then
        Kpar  = 2
        E = 1089.860e3
      else if(xprob.ge.98.7808.and.xprob.lt.99.0798)then
        Kpar  = 2
        E = 1138.230e3
      else if(xprob.ge.99.0798.and.xprob.lt.99.7628)then
        Kpar  = 2
        E = 1180.838e3
      else if(xprob.ge.99.7628.and.xprob.lt.99.829)then
        Kpar  = 2
        E = 1193.230e3
      else if(xprob.ge.99.829.and.xprob.lt.99.833)then
        Kpar  = 2
        E = 1199.670e3
      else if(xprob.ge.99.833.and.xprob.lt.99.84)then
        Kpar  = 2
        E = 1442.700e3
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
      else if(xprob.ge.43.224+gweight.and.xprob.lt.43.2534+gweight)then
        Kpar  = 1
        E = 108.363e3
      else if(xprob.ge.43.2534+gweight.and.xprob.lt.49.4534+gweight)then
        Kpar  = 1
        E = 155.249e3
      else if(xprob.ge.49.4534+gweight.and.xprob.lt.50.3274+gweight)then
        Kpar  = 1
        E = 183.230e3
      else if(xprob.ge.50.3274+gweight.and.xprob.lt.50.5044+gweight)then
        Kpar  = 1
        E = 187.346e3
      else if(xprob.ge.50.5044+gweight.and.xprob.lt.50.54+gweight)then
        Kpar  = 1
        E = 188.232e3
      else if(xprob.ge.50.54+gweight.and.xprob.lt.52.49+gweight)then
        Kpar  = 1
        E = 210.209e3
      else if(xprob.ge.52.49+gweight.and.xprob.lt.52.845+gweight)then
        Kpar  = 1
        E = 238.190e3
      else if(xprob.ge.52.845+gweight.and.xprob.lt.52.9181+gweight)then
        Kpar  = 1
        E = 242.306e3
      else if(xprob.ge.52.9181+gweight.and.xprob.lt.52.9325+gweight)then
        Kpar  = 1
        E = 243.192e3
      else if(xprob.ge.52.9325+gweight)then
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
