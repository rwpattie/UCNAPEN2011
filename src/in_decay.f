      subroutine in_decay 
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from In114
c    taken from :
c
c    http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=114IN&unc=nds
c------------------------------------------------------------------------

      gweight = 53.96
      eweight = 149.69

      x = (gweight+eweight)*rand(1.d0)
   
        if(x.le.65.0)then
          kpar =1 
          E = 2.84e3
        else if(x.gt.65.0.and.x.le.70.98)then
          kpar = 1
          E = 20.1e3
        else if(x.gt.70.98.and.x.le.111.08)then
          kpar = 1
          E = 162.33e3
        else if(x.gt.111.08.and.x.le.142.98)then
          kpar = 1 
          E =186.03e3
        else if(x.gt.142.98.and.x.le.149.69)then
          kpar = 1
          E = 189.44e3
        else  if(x.le.4.78+eweight)then
          kpar = 2
          E = 3.29e3
        else if(x.gt.4.78+eweight.and.x.le.14.58+eweight)then
          kpar =2 
          E = 24.002e3
        else if(x.gt.14.58+eweight.and.x.le.32.78+eweight)then
          kpar =2 
          E = 24.21e3
        else if(x.gt.32.78+eweight.and.x.le.34.41+eweight)then
          kpar =2 
          E = 27.238e3
        else if(x.gt.34.41+eweight.and.x.le.37.56+eweight)then
          kpar =2 
          E = 27.276e3
        else if(x.gt.37.56+eweight.and.x.le.38.4+eweight)then
          kpar =2
          E = 27.863e3
        else if(x.gt.38.4+eweight.and.x.le.53.96+eweight)then
          kpar =2
          E = 190.27e3
        endif

      phi   = 2.*pi*rand(1.d0)
      costh = 1. - 2.*rand(1.d0)
      theta = dacos(costh)
 
      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)      

1000  continue
      y = (0.15 - 0.3*rand(1.d0))
      x = (0.15 - 0.3*rand(1.d0))
      if(sqrt(x**2 + y**2) .gt. 0.15) goto 1000      
      z = 1.0

      return
      end
