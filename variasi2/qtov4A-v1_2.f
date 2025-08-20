      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     TEST1 linear EOS+corr initial boundary 28 Mei 2020
C                  
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IL, IN    
      DIMENSION YA(10), EK(4,10), Y(10)


c      OPEN (unit=2,STATUS='unknown',FILE='CatatanB145SC.dat')
       OPEN (unit=3,STATUS='unknown',FILE='radmass-vfortran.dat')
       OPEN (unit=1,STATUS='unknown',FILE='profil-vfortran.dat')

C     IM = NUMBER OF EQUATIONS Y(I)=Pressure, Y(2)=NS Mass and Y(3)=E density

      HC  = 197.327D0
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
      IM=4   ! 3 degrees of freedom plus 1
      IN=IM-1  ! 3 degrees of freedom
c---------------------------------------------------------------------
      XL=1.0D-3
      ALPHA=4.D0
C     XL=1.0D-2
C     ALPHA=5.0D-3
      XC=XL*ALPHA      
c      XC=1.0D-3
      
c---------------------------------------------------------------------      
       !DO 10 IL=2,1600,2
       FIXEDIL=300
       IL=FIXEDIL
       PC=1.D0*IL
 
       YA(10)=XL
c       PC=800.D0
       EDC=FED(PC)
       FIXP=1.D0-2.D0*XL*XL*PI*GS*(3.D0*PC+EDC)/3.D0
c     FIXM=1.D0+XL*XL/(XC*XC)-XL*XL*XL/(XC*XC*XC)*DATANH(XC/XL)
c     we used series expansion above because in fortran77 no function arctanh.
c     In fortran 90/95 the function is exist, but because the code is written
c     with 77 version, if we compile with gfortran, the code is not always
c     stable                   !!!       
       !FIXM=2.D0/3.D0-(XC/XL)**2/5.D0-(XC/XL)**4/7.D0-(XC/XL)**6/9.D0
       FIXM=1.D0
       IF (ABS(ALPHA)<1)  FIXM=1.D0 + (XC/XL)**(-2)  
     &     - LOG((1 + (XC/XL))/(1 - (XC/XL)))/(2.D0*(XC/XL)**3)

c       FIXP=1.0D0
c       FIXM=1.0D0
       
       PCC=PC-2.D0*PI*GS*XC*XC*(PC+EDC)*(3.D0*PC+EDC)*FIXP/3.D0
       MCC=4.D0*PI*XC*XC*XC*EDC*FIXM/3.D0
       NP=0
       XNUC=0.0D0
 40    XNUC=XNUC  ! metric at center

      !WRITE(*,*) "XL <",SQRT(3/(2*PI*GS*(EDC+3.D0*PCC)))/1.D3
      IF ((EDC+3.D0*PCC).GE.0.D0) THEN
            IF (XL.GE.SQRT(3/(2*PI*GS*(EDC+3.D0*PCC)))) THEN
                  WRITE(*,*) "XL is too large"
                  GOTO 10
            END IF
      ELSE
            WRITE(*,*) "SEC is violated"
            GOTO 10
      END IF
       
      Y(1)=PCC

      Y(2)=MCC
      
      Y(3)=XNUC
c-----------------------------------------------------------------------
c     Y(1)=PC 
c      Y(2)=1.0D-8
      
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      
      P0=Y(1)
 
      Y(4)=FED(P0)
    

     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D0
      NS=32
      !XL=30.0D3

      H=PU/NS
c     XP should be larger than XC in order to avoid the unphysical behavior
c      near center!!    
      XP=1.0D0
      HH=H/(2.0D0)

      IF (XP.LT.XC) THEN
            WRITE(*,*) "XP=",XP," is NOT larger than XC=",XC
            STOP
      END IF

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH

C    COMPUTE K1, L1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K2, L2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
            
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      


         YA(4)=FED(P0)

         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K3, L3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  
            YA(4)=FED(P0)          

         XA=XM


         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K4, L4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)
         
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

 

            YA(4)=FED(P0)
     

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)
         
          
c--------------------------------------------------------------------------
c Presure vs energy density relation
c---------------------------------------------------------------------      

  


          Y(4)=FED(P0)
          
       END DO
       
       XA=XP
       PRESS=Y(1)
       MASST=Y(2)
       EDEN=Y(4)
       XX=(2.D0*GS*MASST*MSS/XA)
     &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
     &   /(1.D0-2.D0*GS*MASST*MSS/XA)
       TPA=SQRT(1.D0+(XL/XA)**2*XX)
       
       
    

      
      IF (NP .EQ. 1) THEN
      IF (IL.EQ.FIXEDIL) THEN 
C       WRITE(1,*)IL,LI,(XP/1.D3),Y(1),Y(2),Y(4),TPA,EXP(Y(3)),
C      & (1.D0-2.D0*GS*MASST*MSS/XA)
       WRITE(1,*)(XP/1.D3),Y(1),Y(2)
     &      ,(1.D0-2.D0*GS*MASST*MSS/XA),EXP(Y(3))
       WRITE(*,*)(XP/1.D3),Y(1),Y(2)
      ENDIF
      ENDIF
     
      


       PS=Y(1)
       PMIN=1.0D-8
c       PMIN=2.0D-5
      
      IF (MCC.LT.0.0D0) THEN
            WRITE(*,*) "MCC negative"
            WRITE(*,*) Y(2)
            GOTO 10
      ELSE IF (ABS(PS).GT.10*ABS(PCC)) THEN
            WRITE(*,*) "ABS(P) =",ABS(PS),"> 10 ABS(PCC) =",10*ABS(PCC)
            GOTO 10
      ELSE IF (LI.GE.1000000) THEN
            WRITE(*,*) "Max iteration = 1000000 reached"
            WRITE(*,*)"p=",Y(1),"m=",Y(2),"rho=",Y(4)
            GOTO 10
      END IF

      IF (PS .GT. PMIN  ) GOTO 28
      
       
       IF(NP.EQ.0) THEN
            TPB=1.D0-2.D0*GS*Y(2)*MSS/XP
            !WRITE(*,*)LOG(TPB)-Y(3)
            IF(ABS(LOG(TPB)-Y(3)).GT.1.D-15) THEN
                  XNUC=XNUC+LOG(TPB)-Y(3)
                  GOTO 40
            ELSE
                  NP=NP+1
                  GOTO 40
            ENDIF
       ENDIF
    

       XLUP=SQRT(3/(2*PI*GS*(EDC+3.D0*PCC))) 

c      WRITE(2,*)IL,(XP/1.D3), (Y(I),I=1,IM)
        WRITE(3,*)PC,(EDC/1.D3),(XP/1.D3),Y(2),
     &  2.0D0*GS*Y(2)*MSS/XP,(XLUP/1.D3)
c       write(*,*)X
        WRITE(*,*)PC,(EDC/1.D3),(XP/1.D3),Y(2),
     &  2.0D0*GS*Y(2)*MSS/XP,(XLUP/1.D3)
 
   

  10   CONTINUE
      
      
      END

      SUBROUTINE FUNCT(EK,J,YA,XA,H)
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
c     XL=1.0D-4
      XL=YA(10)
      EDEN=YA(4)
      PRESS=YA(1)
      MASST=YA(2)

      XX=(2.D0*GS*MASST*MSS/XA)
     &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
     &   /(1.D0-2.D0*GS*MASST*MSS/XA)

      TPA=SQRT(1.D0+(XL/XA)**2*XX)

C        AA=3.D0*XL*XL*PRESS/(EDEN*XA*XA*TPA)
C      &  - 3.D0*XL*XL*MASST*MSS
C      &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
C      &   /(4.D0*PI*EDEN*XA*XA*XA*XA*XA*TPA)
C      &  - GS*XL*XL*(MASST*MSS)**2
C      &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
C      &   /(2.D0*PI*EDEN*XA*XA*XA*XA*XA*XA*TPA)
C      &   /(1.D0-2.D0*GS*MASST*MSS/XA)
C      &  - MASST*MSS*(1.D0-TPA)
C      &   /(4.D0*PI*EDEN*XA*XA*XA)
C      &  - (1.D0-2.D0*GS*MASST*MSS/XA)*(1.D0-TPA)
C      &   /(4.D0*PI*GS*EDEN*XA*XA)
C      &  + (1.D0-2.D0*GS*MASST*MSS/XA)*(1.D0-TPA)**2
C      &   /(8.D0*PI*GS*EDEN*XL*XL)
C      &  + XL*XL/(EDEN*XA*TPA)
C      &   *(-(PRESS+EDEN)*XX/XA/(1.D0+TPA))

C        BB=4.D0*PI*XL*XL*PRESS*XA/(MASST*MSS*TPA)
C      &  - XL*XL/(XA*XA*TPA)
C      &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
C      &  - 2.D0*GS*XL*XL*MASST*MSS
C      &   *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MASST*MSS))
C      &   /(XA*XA*XA*TPA)
C      &   /(1.D0-2.D0*GS*MASST*MSS/XA)
C      &  - (1.D0-TPA)

      AA1=(3.D0*XL**2*PRESS)/
     &     (XA**2*TPA*EDEN)
	  AA2=(-3.D0*XL**2*MSS*MASST
     &     *(1.D0 + (4.D0*PI*XA**3*PRESS)/(MSS*MASST)))/
     &     (4.D0*PI*XA**5*TPA*EDEN)
      AA3=-(GS*XL**2*MSS**2*MASST**2
     &     *(1.D0 + (4.D0*PI*XA**3*PRESS)/(MSS*MASST)))/
     &     (2.D0*PI*XA**6*(1.D0-2.D0*GS*MASST*MSS/XA)*
     &      TPA*EDEN)
	  AA4=-(MSS*MASST*(1.D0 - TPA))
     &     /(4.D0*PI*XA**3*EDEN)
	  AA5=-((1.D0-2.D0*GS*MASST*MSS/XA)*(1.D0 - TPA))
     &     /(4.D0*GS*PI*XA**2*EDEN)
	  AA6=((1.D0-2.D0*GS*MASST*MSS/XA)*(1.D0 - TPA)**2)
     &     /(8.D0*GS*XL**2*PI*EDEN)
	  AA7= -(((1.D0 - TPA)*(-PRESS - EDEN))/
     &    (TPA*EDEN))
	 
	 
	  BB1=(4.D0*XL**2*PI*XA*PRESS)/
     &     (MSS*MASST*TPA)
	  BB2=-((XL**2
     &     *(1.D0 + (4.D0*PI*XA**3*PRESS)/(MSS*MASST)))
     &     /(XA**2*TPA))
	  BB3=(-2.D0*GS*XL**2*MSS*MASST
     &     *(1.D0 + (4.D0*PI*XA**3*PRESS)/(MSS*MASST)))
     &     /(XA**3*(1.D0 - (2.D0*GS*MSS*MASST)/(XA))*TPA)
	  BB4=-(1.D0-TPA)
	 
	  AA=AA1+AA2+AA3+AA4+AA5+AA6+AA7
	  BB=BB1+BB2+BB3+BB4
	   
	   
	  TMA=(1.D0+AA)/(1.D0+BB)
c-------------------------------------------------------------      
c       TPA=0.0D0
c       TMA=1.0D0
c      OPEN (unit=7,STATUS='unknown',FILE='cekTMAXL10mm.dat')
c      write(7,*)TPA,FF,GG,TMA
c-------------------------------------------------------------      
  
      EK(J,1)=-EDEN*(1.D0+PRESS/EDEN)*XX/((1.D0+TPA)*XA)*H
      EK(J,2)=4.D0*PI*XA*XA*EDEN*TMA*H/MSS
      EK(J,3)=2.D0*XX/((1.D0+TPA)*XA)*H
	  
c	  EK(J,1)=H*(-((XA*(1.D0 - TPA)*(-PRESS - EDEN))/XL**2))
c	  EK(J,2)=H*4.D0*PI*XA*XA*EDEN*TMA/MSS
c	  EK(J,3)=H*(-2.D0*((XA*(1.D0 - TPA))/XL**2))

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
	  W = 1.D0
      FED= P0/W+4.D0*BMF
   
       RETURN
       END
c-----------------------------------------------------------------------      
 