      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     TEST1 linear EOS+corr initial boundary 28 Mei 2020
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN    
      DIMENSION YA(10), EK(4,10), Y(10)
      DOUBLE PRECISION LP,LAM


       OPEN (unit=2,STATUS='unknown',FILE='parameters.dat')
       OPEN (unit=3,STATUS='unknown',FILE='radmass-vTOV.dat')
       OPEN (unit=1,STATUS='unknown',FILE='profil-vTOV.dat')

C     IM = NUMBER OF EQUATIONS Y(I)=Pressure, Y(2)=NS Mass and Y(3)=E density

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      IM=3
      IN=IM-1
      
      !DO 10 IL=1,300,1
      FIXEDIL=1
      IL=FIXEDIL
      PC=-1.D-9*IL
       
       XP=10.0D3!*1.616255D-35
       !PC=-1.D9!*2.89598D81
       LP=3.D0!*1.616255D-35 
       XPA=XP
       LAM=10.0D0
       
       EDC=FED(PC)
       
       PCC=PC 
       MCC=1.D-24 !4.D0*PI*XC*XC*XC*EDC/3.D0 
       
      Y(1)=PCC
       XA=XP
       MSUR=(XA**3*(-1.D0 + LAM**2 - 8.D0*GS*LP**2*PI*Y(1)))/
     &     (2.D0*GS*(LP**2 + XA**2*(-1.D0 + LAM**2)))
      Y(2)=MSUR
c-----------------------------------------------------------------------
c     Y(1)=PC 
c      Y(2)=1.0D-8
      
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      
      P0=Y(1)
 
      Y(3)=FED(P0)
      
      
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D0
      NS=40000
      !XL=30.0D3

      H=-PU/NS
      HH=H/(2.0D0)

      !IF (IL.EQ.FIXEDIL) THEN
      ! WRITE(1,*)(XP/1.D3),Y(1),Y(2),Y(3)!,(XP-XPA)/1.D3*1.D1
      !ENDIF
      WRITE(*,*)(XP/1.D3),Y(1),Y(2),Y(3)
       
C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH
         
         AWAL=Y(1)

C    COMPUTE K1, L1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H,LP,LAM)

C    COMPUTE K2, L2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
            
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      


         YA(3)=FED(P0)

         XA=XM

         CALL FUNCT(EK,J,YA,XA,H,LP,LAM)

C    COMPUTE K3, L3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  
            YA(3)=FED(P0)          

         XA=XM


         CALL FUNCT(EK,J,YA,XA,H,LP,LAM)

C    COMPUTE K4, L4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)
         
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

 

            YA(3)=FED(P0)
     

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H,LP,LAM)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)
         
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  


          Y(3)=FED(P0)
          Y(2)=(XA**3*(-1.D0 + LAM**2 - 8.D0*GS*LP**2*PI*Y(1)))/
     &         (2.D0*GS*(LP**2 + XA**2*(-1.D0 + LAM**2)))/MSS
          
       END DO
       
       XA=XP



       !WRITE(1,*)(XP/1.D3),Y(1),Y(2),Y(3)
       WRITE(*,*)(XP/1.D3),(Y(1)/1.D4),Y(2),(Y(3)/1.D4)

       IF (IL.EQ.FIXEDIL) THEN
        WRITE(1,*)(XP/1.D3),(Y(1)/1.D4),Y(2),(Y(3)/1.D4)!,(XP-XPA)/1.D3*1.D1
       ENDIF
     

       PS=Y(1)
       PMIN=1.0D-9
c       PMIN=2.0D-5
      
      ! IF (MCC.LT.0.0D0) THEN
            ! WRITE(*,*) "MCC negative"
            ! WRITE(*,*) Y(2)
            ! GOTO 10
      ! ELSE IF (ABS(PS).GT.10*ABS(PCC)) THEN
            ! WRITE(*,*) "ABS(P) =",ABS(PS),"> 10 ABS(PCC) =",10*ABS(PCC)
            ! GOTO 10
      ! ELSE IF (LI.GE.1000000) THEN
            ! WRITE(*,*) "Max iteration = 1000000 reached"
            ! WRITE(*,*)"p=",Y(1),"m=",Y(2),"rho=",Y(3)
            ! GOTO 10
      ! END IF

      !IF (ABS(PS) .GT. PMIN  ) GOTO 28
      !IF (XA/XPA .GT.0) GOTO 28
      IF (XA .GT. PU) GOTO 28

    

 50     PRINT*,""

        WRITE(3,*)(Y(1)/1.D4),(Y(3)/1.D4),(XPA/1.D3),MSUR,
     &  GS*Y(2)*MSS/XP

        WRITE(*,*)(Y(1)/1.D4),(Y(3)/1.D4),(XPA/1.D3),MSUR,
     &  GS*Y(2)*MSS/XP
     

       
       IF (IL.EQ.FIXEDIL) THEN
        WRITE(2,*)LP, LAM, PCC, XPA, PU
       ENDIF
 
   

  10   CONTINUE
      
      
      END

      SUBROUTINE FUNCT(EK,J,YA,XA,H,LP,LAM)
C     *********************************************************
C     DEFINES TOV EQUATIONS 
C
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(4,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
c--------------------------------------------------------------
      EDEN=YA(3)
      PRESS=YA(1)
      MASST=YA(2)

      !  XX=(2.D0*GS*MASST*MSS/XA)
      !&   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
      !&   /(1.D0-2.D0*GS*MASST*MSS/XA)
c-------------------------------------------------------------      
c       TPA=0.0D0
c       TMA=1.0D0
c      OPEN (unit=7,STATUS='unknown',FILE='cekTMAXL10mm.dat')
c      write(7,*)TPA,FF,GG,TMA
c-------------------------------------------------------------      
  
      !EK(J,1)=-EDEN*(1.D0+PRESS/EDEN)*XX/(2.D0*XA)*H
      
      !EK(J,2)=4.D0*PI*XA*XA*EDEN*H/MSS 
      
      
      
      PERSP=((1.D0 + LAM)**2*(LP**4 + LP**2*XA**2*(-2.D0 + LAM)
     &      - XA**4*(-1.D0 + LAM**2)))/
     &      (4.D0*GS*LP**2*PI*XA*(LP**2 - XA**2*(1.D0 + LAM))
     &      *(LP**2 + XA**2*(-1.D0 + LAM**2))) + 
     &      (2.D0*XA*(1.D0 + LAM)*(-(XA**4*(-1.D0 + LAM)
     &      *(1.D0 + LAM)**2) + LP**4*(1.D0 + 2*LAM) + 
     &      LP**2*XA**2*(-2.D0 - 2.D0*LAM + LAM**2 
     &      + LAM**3))*PRESS)/
     &      (LP**2*(LP**2 - XA**2*(1.D0 + LAM))
     &      *(LP**2 + XA**2*(-1.D0 + LAM**2)))
      EK(J,1)=H*PERSP
      EK(J,2)=0.D0
      
      ! PRINT*,XA,EK(J,1),EK(J,2),"tes"
     
      
      

      RETURN
      END

 
C-------------------------------------------------
C     MIT BAG B=145 MeV^4
c-------------------------------------------------
  
      
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      HC  = 197.327D0
      B=(145.D0)**4  ! PC until 800
      !B=(185.D0)**4  ! PC until 1600
      BMF=B/(HC*HC*HC)
      FED= -1.D0*P0+240.D0  ! 4.D0*BMF
   
       RETURN
       END
c-----------------------------------------------------------------------      
