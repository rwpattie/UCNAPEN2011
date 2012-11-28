      subroutine sr_decay 
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy1
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from Sr85
c    taken from :
c
c  http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=85Sr&unc=nds
c------------------------------------------------------------------------

      x = 1.d0*rand(1.d0)
      if(x.le.0.89866)then
        E = 498.8070e3
      else if(x.gt.0.89866.and.x.le.1.0)then
        E = 511.9416e3
      endif
   
      phi = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)
 
      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)
      kpar = 1
      y = -.5
      x = .1
      z = 1.

      return
      end
