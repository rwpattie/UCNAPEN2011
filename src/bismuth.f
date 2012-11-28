      subroutine bi_decay 
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy1
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
      DIMENSION ILBH(5)
c----------------------------------------------------------------------c
c    This give the approximate level scheme for IC beta's from Bi207
c    taken from :
c
c  http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=207BI&unc=nds
c------------------------------------------------------------------------
      gprob = 1.0*rand(1.d0)
      
      if(grob.le.0.745)then
         ngamma = 2
      else 
         ngamma = 1
      endif
   
      phi = 2.*pi*rand(1.d0)
      costh = 1.-2.*rand(1.d0)
      theta = dacos(costh)

      ILBH(1) = 2
      ILBH(2) = 1
      ILBH(3) = 1
      ILBH(4) = 0
      ILBH(5) = 1
 
      u = dsin(theta)*dcos(phi)
      v = dsin(theta)*dsin(phi)
      w = dcos(theta)
      kpar = 1
      x=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z = 1.!155.d0

      if(ngamma.eq.1)then
           xe = 100.0*rand(1.d0)          
          if(xe.le.1.537)then
            E = 481693.5
           else if(xe.gt.1.537.and.xe.le.1.979)then
            E = 5.53872e5
           else if(xe.gt.1.979.and.xe.le.2.190)then
            E = 5.658473e5
           else if(xe.gt.2.105)then 
            E    = 5.69698e5
            kpar = 2 
          endif
      elseif(ngamma.eq.2)then
          xe = 100.0*rand(1.d0)+2.105
          ! Create first particle
          if(xe.gt.2.105.and.xe.le.3.945)then
             E = 1.047795e6
          else if(xe.gt.3.945.and.xe.le.10.975)then
             E = 9.75651e5
          else if(xe.gt.10.975.and.xe.le.11.515)then
             E = 1.0598e6
          else if(xe.gt.11.515)then
             E = 1.063656e6
             kpar = 2
          endif
          ! Create Second 
            theta  = 1.0- 2.0*rand(1.d0)
            psi    = 2*pi*rand(1.d0)
            us     = dsin(dacos(theta))*dcos(psi)
            vs     = dsin(dacos(theta))*dsin(psi)
            ws     = theta

           xe = 100.*rand(1.d0)
           KPAR2 = 1
           if(xe.le.1.52)then
            E2 = 4.816935e5
           else if(xe.gt.1.52.and.xe.le.1.958)then
            E2 = 5.53872e5
           else if(xe.gt.1.958.and.xe.le.2.105)then
            E2 = 5.658473e5
           else if(xe.gt.2.105)then
            E2 = 5.69698e5
            KPAR2 = 2
           endif
           CALL STORES(E2,X,Y,Z,US,VS,WS,1,KPAR2,ILBH)
      endif
        
      gprob = 1.0*rand(1.d0)
      if(grob.gt.0.029)then
         theta   = 1.0- 2.0*rand(1.d0)
          psi    = 2*pi*rand(1.d0)
          us     = dsin(dacos(theta))*dcos(psi)
          vs     = dsin(dacos(theta))*dsin(psi)
          ws     = theta
          E3     = 56.7e3
          CALL STORES(E3,X,Y,Z,US,VS,WS,1,1,ILBH)
      endif
   
      return
      end
c-----------------------------------------------------------------------------c
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
      x = 0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y = 0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z = 1

      return
      end
c-----------------------------------------------------------------
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
      x=0.4*(0.5-1*rand(1.d0)) ! for background runs, tin source is
      y=0.4*(0.5-1*rand(1.d0)) !   approximately at (-5.5,0,155).
      z = 1.

      return
      end
c----------------------------------------------------------------------
      subroutine xe_135_decay(ptype)
c----------------------------------------------------------------------
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0,ME=510.9870d3)
      EXTERNAL rand,energy3
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
      DIMENSION ILBH(5)

      betaprob = rand(1.d0)
      if(betaprob.gt.0.25.and.betaprob.le.0.50)goto 1100
      if(betaprob.gt.0.50.and.betaprob.le.0.75)goto 1200
      if(betaprob.gt.0.75)goto 1300
      
      betaprob = rand(1.d0)
      nbetatype = 0
      qend = 0.d0

      ILBH(1) = 2
      ILBH(2) = 1
      ILBH(3) = 1
      ILBH(4) = 0
      ILBH(5) = 1

 
c     select the beta decay end point
c      betaprob = 0.50     
 
      if(betaprob.le.0.96)then
         qend = 915.0d3
         nbetatype = 0
      else if(betaprob.gt.0.96.and.betaprob.le.0.9911)then
         qend = 557.0d3
         nbetatype = 1
      else if(betaprob.gt.0.9911.and.betaprob.le.0.997)then
         qend = 757.0d3
         nbetatype = 2
      else if(betaprob.gt.0.997.and.betaprob.le.0.99823)then
         qend = 103.0d3
         nbetatype = 3
      else if(betaprob.gt.0.99823)then
         qend = 184.0e3
         nbetatype = 4
      endif

      rad   = 6.5*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
      x     = rad*dcos(psi)
      y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      
      Kpar  = 1
      e     = energy3(qend)
   
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      us     = dsin(dacos(theta))*dcos(psi)
      vs     = dsin(dacos(theta))*dsin(psi)
      ws     = theta
     
      if(nbetatype.eq.0)then
        augerprob = 0.966388*rand(1.0d0)
        if(augerprob.le.0.90)then
          CALL STORES(249.794e3,X,Y,Z,US,VS,WS,1,2,ILBH)  
        else if(augerprob.gt.0.90.and.augerprob.le.0.9561)then
          CALL STORES(213.809e3,X,Y,Z,US,VS,WS,1,1,ILBH)
        else if(augerprob.gt.0.9561.and.augerprob.le.0.9643)then
          CALL STORES(244.08e3,X,Y,Z,US,VS,WS,1,1,ILBH)
        else if(augerprob.gt.0.9643.and.augerprob.le.0.96599)then
          CALL STORES(248.577e3,X,Y,Z,US,VS,WS,1,1,ILBH)
        else if(augerprob.gt.0.96599.and.augerprob.le.0.9634)then
          CALL STORES(249.563e3,X,Y,Z,US,VS,WS,1,1,ILBH) 
        else if(augerprob.gt.0.9634)then
          CALL STORES(249.778e3,X,Y,Z,US,VS,WS,1,1,ILBH)
        endif
      else if(nbetatype.eq.1)then
        augerprob = 2.921005*rand(1.d0)
        if(augerprob.le.2.90)then
         CALL STORES(608.185e3,X,Y,Z,US,VS,WS,1,2,ILBH)
        else if(augerprob.gt.2.90.and.augerprob.le.2.9182)then
         CALL STORES(572.200e3,X,Y,Z,US,VS,WS,1,1,ILBH)
        else if(augerprob.gt.2.9182.and.augerprob.le.2.92053)then
         CALL STORES(602.471e3,X,Y,Z,US,VS,WS,1,1,ILBH)
        else if(augerprob.gt.2.92053)then
         CALL STORES(606.968e3,X,Y,Z,US,VS,WS,1,1,ILBH)
        endif
      endif 
      ptype = 1
      return

1100  continue
      ptype = 2
      betaprob = 97.137*rand(1.d0)
 
      if(betaprob.lt.61.6)then
        E = 129.369d3
      else if(betaprob.ge.61.6.and.betaprob.lt.90.4)then
        E = 158.477d3
      else if(betaprob.ge.90.4.and.betaprob.lt.96.99)then
        E = 162.788d3
      else if(betaprob.ge.96.99.and.betaprob.lt.97.137)then
        E = 163.914d3
      endif
   
      rad   = 6.5*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
      x     = rad*dcos(psi)
      y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      Kpar  = 1
    
      return
  
1200  continue
      ptype = 3
      betaprob = rand(1.d0)
 
      if(betaprob.le.0.99)then
        qend = 346.4e3
        e = energy3(qend)
      else if(betaprob.gt.0.99.and.betaprob.le.0.9981)then
        qend = 266.8e3
        e = energy3(qend)
      else if(betaprob.gt.0.9981)then
        qend = 43.5e3
        e = energy3(qend)
      endif 
      
      rad   = 6.5*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
      x     = rad*dcos(psi)
      y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      Kpar  = 1

      return
 
1300  continue
      ptype = 4
      betaprob = rand(1.d0)
 
      betaprob = 90.955*rand(1.d0)

      if(betaprob.lt.63.5)then
        E = 198.660d3
      else if(betaprob.ge.63.5.and.betaprob.lt.84.18)then
        E = 227.768e3
      else if(betaprob.ge.84.18.and.betaprob.lt.89.73)then
        E = 232.079e3
      else if(betaprob.ge.89.73)then
        E = 233.013e3
      endif

      rad   = 6.5*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      z     = 220.0*(1.0 - 2.0*rand(1.d0))
      x     = rad*dcos(psi)
      y     = rad*dsin(psi)
      theta = 1.0- 2.0*rand(1.d0)
      psi   = 2*pi*rand(1.d0)
      u     = dsin(dacos(theta))*dcos(psi)
      v     = dsin(dacos(theta))*dsin(psi)
      w     = theta
      Kpar  = 1

           
      return
      end      
c--------------------------------------------------------------------------

      double precision function energy3(qend)
      implicit double precision(a-h,J-M,o-z), integer*4(i,n)
      parameter(emass=510.991e3)
      parameter(pi   =3.141592654d0)
      common/rseed/iseed1,iseed2
      external rand
c
      E0 = qend/emass + 1
50    E=(E0-1.0)*RAND(1.D0)
      y=2.80*RAND(1.D0)
      G=(-E/(DSQRT(E**2+2*E)*137))
      FERMI=(2.0*PI*G)/(DEXP(2*PI*G)-1)
      f=Fermi*DSQRT(E**2+2*E)*(E0-(E+1))**2*(E+1)
      if (f.lt.y) goto 50

      E = E*emass
      energy3 = E
c 
      return
      end
      


