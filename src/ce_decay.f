      subroutine ce_decay
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy1
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from Bi207
c    taken from :
c
c  http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=207BI&unc=nds
c------------------------------------------------------------------------

      x = 20.0503*rand(1.d0)
      if(x.le.17.146)then
        E = 126.932d3
      else if(x.gt.17.146.and.x.le.19.444)then
        E = 159.591d3
      else if(x.gt.19.444.and.x.le.19.9191)then
        E = 164.496d3
      else if(x.gt.19.9191.and.x.le.20.0503)then
        E = 165.587d3
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
      z = 1.!155.d0

      return
      end
            
