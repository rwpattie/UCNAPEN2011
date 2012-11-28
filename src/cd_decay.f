      subroutine cd_decay 
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from Cd109
c    taken from :
c
c    http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=109Cd&unc=nds
c------------------------------------------------------------------------

      gweight = 116.11
      eweight = 284.1

      x = (gweight+eweight)*rand(1.d0)
   
        if(x.le.167.1)then
          kpar =1 
          E = 2.61e3
        else if(x.gt.167.1.and.x.le.187.9)then
          kpar = 1
          E = 18.5e3
        else if(x.gt.187.9.and.x.le.229.6)then
          kpar = 1
          E = 62.5196e3
        else if(x.gt.229.6.and.x.le.273.6)then
          kpar = 1 
          E =84.2278e3
        else if(x.gt.273.6.and.x.le.282.5)then
          kpar = 1
          E = 87.3161e3
        else if(x.gt.282.5.and.x.le.284.1)then
          kpar = 1 
          E = 87.9384e3        
        else  if(x.le.10.3+eweight)then
          kpar = 2
          E = 2.98e3
        else if(x.gt.10.3+eweight.and.x.le.40.1+eweight)then
          kpar =2 
          E = 21.99e3
        else if(x.gt.40.1+eweight.and.x.le.96.1+eweight)then
          kpar =2 
          E = 22.163e3
        else if(x.gt.96.1+eweight.and.x.le.100.9+eweight)then
          kpar =2 
          E = 24.912e3
        else if(x.gt.100.9+eweight.and.x.le.110.1+eweight)then
          kpar =2 
          E = 24.943e3
        else if(x.gt.110.1+eweight.and.x.le.112.41+eweight)then
          kpar =2
          E = 25.455e3
        else if(x.gt.112.41+eweight.and.x.le.116.11+eweight)then
          kpar =2
          E = 88.0336e3
        endif

      phi   = 2.*pi*rand(1.d0)
      costh = 1. - 2.*rand(1.d0)
      theta = dacos(costh)
 
      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)      

      y = (0.1 - 0.2*rand(1.d0))
      x = (0.1 - 0.2*rand(1.d0))
      z = -9.9

      return
      end
