      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     Plotter for analytic solution
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN, IP, IPM
      DIMENSION YA(10), EK(4,10), Y(10)
      DOUBLE PRECISION LP,LAM
      
      
      
      PRINT*,"Hello world"
      
      ! opening the file for reading
      open (2, file = 'parameters.dat', status = 'old')

      read(2,*) LP, LAM, PCC, XP, PU
   
      close(2)
   
      write(*,*) LP, LAM, PCC, XP, PU
       
      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      
      ! opening the file for writing
      open (1, file = 'analytic.dat', status = 'unknown')
      
      XA=XP
      H=-PU/XP
      IPM=1000000000
      
      PRINT*,PU,XP,H
      !STOP
      
      IP=1
      DO IP=1,IPM,1
        XA=XA+H*(IP-1)
        PRESS=  -1.D0/(8.D0*PI*GS*XA*XA)*(1.D0-(XA/XP)**2*
     &          EXP((1.D0+LAM)*(XA*XA-XP*XP)/(LP*LP)))
        EDEN =   1.D0/(8.D0*PI*GS*XA*XA)*(1.D0+(XA/XP)**2*
     &          EXP((1.D0+LAM)*(XA*XA-XP*XP)/(LP*LP)))
        MASST=  (XA**3*(-1.D0 + LAM**2 - 8.D0*GS*LP**2*PI*PRESS))/
     &          (2.D0*GS*(LP**2 + XA**2*(-1.D0 + LAM**2)))/MSS
        WRITE(1,*)(XA/1.D3),PRESS/1.D4,MASST,EDEN/1.D4
        WRITE(*,*)IP,(XA/1.D3),PRESS/1.D4
        IF (XA.LE.PU) STOP
      ENDDO
      
      
      END