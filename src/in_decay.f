      subroutine in_decay 
c----------------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from In114m
c    taken from two sources:
c
c    1) C. Wrede et al., NIM B 269, 1113 (2011)
c    2) http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=114IN&unc=nds
c----------------------------------------------------------------------------c

      gweight = 21.96
      eweight = 80.061

      x = (gweight+eweight)*rand(1.d0)
   
        if(x.le.40.1)then
          kpar =1 
          E = 162.33e3
        else if(x.gt.40.1.and.x.le.72.0)then
          kpar = 1
          E = 186.03e3
        else if(x.gt.72.0.and.x.le.78.71)then
          kpar = 1
          E = 189.44e3
        else if(x.gt.78.71.and.x.le.80.061)then
          kpar = 1 
          E =190.15e3
        else if(x.gt.eweight.and.x.le.eweight+15.56)then
          kpar = 2
          E = 190.27e3
        else if(x.gt.15.56+eweight.and.x.le.18.76+eweight)then
          kpar =2 
          E = 558.43e3
        else if(x.gt.18.76+eweight.and.x.le.21.96+eweight)then
          kpar =2 
          E = 725.24e3
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
