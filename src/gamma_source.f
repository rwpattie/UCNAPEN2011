      subroutine gamma_source(EG,XG,YG,ZG)
c---------------------------------------------------------------------c
      implicit DOUBLE PRECISION(A-H,O-Z), integer*4(i-n)
      parameter(pi=3.141592654d0)
      EXTERNAL rand,energy1
      COMMON/track/E,X,Y,Z,U,V,W,WGHT,KPAR,IBODY,MAT,ILB(5)
      COMMON/RSEED/ISEED1,ISEED2
c-------------------------------------------------------------------c
      KPAR = 2
      E=EG
      X=XG 
      Y=YG
      Z=ZG
   
      thetamax = ATAN(7.5/220.0)      
      THETA = thetamax*rand(1.d0)
      PHI =  2.*pi*rand(1.d0)

      V = sin(theta)*sin(phi)
      W = cos(theta)
      U = sin(theta)*cos(phi)
    
      RETURN
      END
       
      
      

