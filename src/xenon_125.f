c-----------------------------------------------------------------------------C
      subroutine xenon_125
c-----------------------------------------------------------------------------C
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2

      gweight = 99.84789
      eweight = 53.0504096

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
      else if(xprob.ge.98.7808.and.xprob.lt.98.7834)then
        Kpar  = 2
        E = 1108.710e3
      else if(xprob.ge.98.7834.and.xprob.lt.99.0824)then
        Kpar  = 2
        E = 1138.230e3
      else if(xprob.ge.99.0824.and.xprob.lt.99.7654)then
        Kpar  = 2
        E = 1180.838e3
      else if(xprob.ge.99.7654.and.xprob.lt.99.8316)then
        Kpar  = 2
        E = 1193.230e3
      else if(xprob.ge.99.8316.and.xprob.lt.99.8356)then
        Kpar  = 2
        E = 1199.670e3
      else if(xprob.ge.99.8356.and.xprob.lt.99.8372)then
        Kpar  = 2
        E = 1254.350e3
      else if(xprob.ge.99.8372.and.xprob.lt.99.8383)then
        Kpar  = 2
        E = 1318.910e3
      else if(xprob.ge.99.8383.and.xprob.lt.99.83981)then
        Kpar  = 2
        E = 1381.000e3
      else if(xprob.ge.99.83981.and.xprob.lt.99.84681)then
        Kpar  = 2
        E = 1442.700e3
      else if(xprob.ge.99.84681.and.xprob.lt.99.84789)then
        Kpar  = 2
        E = 1562.400e3
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
      else if(xprob.ge.42.776+gweight.and.xprob.lt.42.7914+gweight)then
        Kpar  = 1
        E = 54.961e3
      else if(xprob.ge.42.7914+gweight.and.xprob.lt.42.9814+gweight)then
        Kpar  = 1
        E = 69.687e3
      else if(xprob.ge.42.9814+gweight.and.xprob.lt.43.0224+gweight)then
        Kpar  = 1
        E = 73.803e3
      else if(xprob.ge.43.0224+gweight.and.xprob.lt.43.0303+gweight)then
        Kpar  = 1
        E = 74.689e3
      else if(xprob.ge.43.0303+gweight.and.xprob.lt.43.03101
     1       +gweight)then
        Kpar  = 1
        E = 74.868e3
      else if(xprob.ge.43.03101+gweight.and.xprob.lt.43.24801
     1       +gweight)then
        Kpar  = 1
        E = 80.382e3
      else if(xprob.ge.43.24801+gweight.and.xprob.lt.43.27741
     1       +gweight)then
        Kpar  = 1
        E = 108.363e3
      else if(xprob.ge.43.27741+gweight.and.xprob.lt.43.28335
     1       +gweight)then
        Kpar  = 1
        E = 112.479e3
      else if(xprob.ge.43.28335+gweight.and.xprob.lt.43.28455
     1       +gweight)then
        Kpar  = 1
        E = 113.365e3
      else if(xprob.ge.43.28455+gweight.and.xprob.lt.43.284689
     1       +gweight)then
        Kpar  = 1
        E = 113.544e3
      else if(xprob.ge.43.284689+gweight.and.xprob.lt.49.484689
     1       +gweight)then
        Kpar  = 1
        E = 155.249e3
      else if(xprob.ge.49.484689+gweight.and.xprob.lt.49.491689
     1       +gweight)then
        Kpar  = 1
        E = 177.249e3
      else if(xprob.ge.49.491689+gweight.and.xprob.lt.50.365689
     1       +gweight)then
        Kpar  = 1
        E = 183.230e3
      else if(xprob.ge.50.365689+gweight.and.xprob.lt.50.542689
     1       +gweight)then
        Kpar  = 1
        E = 187.346e3
      else if(xprob.ge.50.542689+gweight.and.xprob.lt.50.578289
     1       +gweight)then
        Kpar  = 1
        E = 188.232e3
      else if(xprob.ge.50.578289+gweight.and.xprob.lt.50.582339
     1       +gweight)then
        Kpar  = 1
        E = 188.411e3
      else if(xprob.ge.50.582339+gweight.and.xprob.lt.50.583539
     1       +gweight)then
        Kpar  = 1
        E = 205.230e3
      else if(xprob.ge.50.583539+gweight.and.xprob.lt.50.583779
     1       +gweight)then
        Kpar  = 1
        E = 209.346e3
      else if(xprob.ge.50.583779+gweight.and.xprob.lt.52.533779
     1       +gweight)then
        Kpar  = 1
        E = 210.209e3
      else if(xprob.ge.52.533779+gweight.and.xprob.lt.52.533826
     1       +gweight)then
        Kpar  = 1
        E = 210.232e3
      else if(xprob.ge.52.533826+gweight.and.xprob.lt.52.888826
     1       +gweight)then
        Kpar  = 1
        E = 238.190e3
      else if(xprob.ge.52.888826+gweight.and.xprob.lt.52.961926
     1       +gweight)then
        Kpar  = 1
        E = 242.306e3
      else if(xprob.ge.52.961926+gweight.and.xprob.lt.52.976326
     1       +gweight)then
        Kpar  = 1
        E = 243.192e3
      else if(xprob.ge.52.976326+gweight.and.xprob.lt.52.977836
     1       +gweight)then
        Kpar  = 1
        E = 243.371e3
      else if(xprob.ge.52.977836+gweight.and.xprob.lt.52.978286
     1       +gweight)then
        Kpar  = 1
        E = 307.050e3
      else if(xprob.ge.52.978286+gweight.and.xprob.lt.52.978357
     1       +gweight)then
        Kpar  = 1
        E = 335.030e3
      else if(xprob.ge.52.978357+gweight.and.xprob.lt.52.981397
     1       +gweight)then
        Kpar  = 1
        E = 338.912e3
      else if(xprob.ge.52.981397+gweight.and.xprob.lt.52.9814115
     1       +gweight)then
        Kpar  = 1
        E = 339.150e3
      else if(xprob.ge.52.9814115+gweight.and.xprob.lt.52.9818315
     1       +gweight)then
        Kpar  = 1
        E = 366.893e3
      else if(xprob.ge.52.9818315+gweight.and.xprob.lt.52.9819165
     1       +gweight)then
        Kpar  = 1
        E = 371.009e3
      else if(xprob.ge.52.9819165+gweight.and.xprob.lt.52.9819337
     1       +gweight)then
        Kpar  = 1
        E = 371.895e3
      else if(xprob.ge.52.9819337+gweight.and.xprob.lt.53.0345337
     1       +gweight)then
        Kpar  = 1
        E = 420.627e3
      else if(xprob.ge.53.0345337+gweight.and.xprob.lt.53.0411337
     1       +gweight)then
        Kpar  = 1
        E = 448.608e3
      else if(xprob.ge.53.0411337+gweight.and.xprob.lt.53.0424637
     1       +gweight)then
        Kpar  = 1
        E = 452.724e3
      else if(xprob.ge.53.0424637+gweight.and.xprob.lt.53.0427327
     1       +gweight)then
        Kpar  = 1
        E = 453.610e3
      else if(xprob.ge.53.0427327+gweight.and.xprob.lt.53.0427644
     1       +gweight)then
        Kpar  = 1 
        E = 453.789e3
      else if(xprob.ge.53.0427644+gweight.and.xprob.lt.53.0432844
     1       +gweight)then
        Kpar  = 1
        E = 602.213e3
      else if(xprob.ge.53.0432844+gweight.and.xprob.lt.53.0437844
     1       +gweight)then
        Kpar  = 1
        E = 602.941e3
      else if(xprob.ge.53.0437844+gweight.and.xprob.lt.53.0438514
     1       +gweight)then
        Kpar  = 1
        E = 630.194e3
      else if(xprob.ge.53.0438514+gweight.and.xprob.lt.53.0439154
     1       +gweight)then
        Kpar  = 1
        E = 630.922e3
      else if(xprob.ge.53.0439154+gweight.and.xprob.lt.53.0439289
     1       +gweight)then
        Kpar  = 1
        E = 634.310e3
      else if(xprob.ge.53.0439289+gweight.and.xprob.lt.53.0439418
     1       +gweight)then
        Kpar  = 1
        E = 635.038e3
      else if(xprob.ge.53.0439418+gweight.and.xprob.lt.53.0441218
     1       +gweight)then
        Kpar  = 1
        E = 693.927e3
      else if(xprob.ge.53.0441218+gweight.and.xprob.lt.53.0441438
     1       +gweight)then
        Kpar  = 1
        E = 721.908e3
      else if(xprob.ge.53.0441438+gweight.and.xprob.lt.53.0469338
     1       +gweight)then
        Kpar  = 1
        E = 813.342e3
      else if(xprob.ge.53.0469338+gweight.and.xprob.lt.53.0472778
     1       +gweight)then
        Kpar  = 1
        E = 841.323e3
      else if(xprob.ge.53.0472778+gweight.and.xprob.lt.53.0473468
     1       +gweight)then
        Kpar  = 1
        E = 845.439e3
      else if(xprob.ge.53.0473468+gweight.and.xprob.lt.53.0473608
     1       +gweight)then
        Kpar  = 1
        E = 846.325e3
      else if(xprob.ge.53.0473608+gweight.and.xprob.lt.53.0484608
     1       +gweight)then
        Kpar  = 1
        E = 868.340e3
      else if(xprob.ge.53.0484608+gweight.and.xprob.lt.53.0485998
     1       +gweight)then
        Kpar  = 1
        E = 892.320e3
      else if(xprob.ge.53.0485998+gweight.and.xprob.lt.53.0486278
     1       +gweight)then
        Kpar  = 1
        E = 900.440e3
      else if(xprob.ge.53.0486278+gweight.and.xprob.lt.53.0488978
     1       +gweight)then
        Kpar  = 1
        E = 904.325e3
      else if(xprob.ge.53.0486278+gweight.and.xprob.lt.53.0489308
     1       +gweight)then
        Kpar  = 1
        E = 932.306e3
      else if(xprob.ge.53.0489308+gweight.and.xprob.lt.53.0491708
     1       +gweight)then
        Kpar  = 1
        E = 974.260e3
      else if(xprob.ge.53.0491708+gweight.and.xprob.lt.53.0492008
     1       +gweight)then
        Kpar  = 1
        E = 1002.240e3
      else if(xprob.ge.53.0492008+gweight.and.xprob.lt.53.0495408
     1       +gweight)then
        Kpar  = 1
        E = 1105.060e3
      else if(xprob.ge.53.0495408+gweight.and.xprob.lt.53.0495828
     1       +gweight)then
        Kpar  = 1
        E = 1133.040e3
      else if(xprob.ge.53.0495828+gweight.and.xprob.lt.53.0503028
     1       +gweight)then
        Kpar  = 1
        E = 1147.670e3
      else if(xprob.ge.53.0503028+gweight.and.xprob.lt.53.0503918
     1       +gweight)then
        Kpar  = 1
        E = 1175.650e3
      else if(xprob.ge.53.0503918+gweight)then
        Kpar  = 1
        E = 1179.770e3
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
