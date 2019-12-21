! OJO, COMO QUE ESTE MAN TAMBIEN REBOTA? SHAPE AND VALUE DEPENDS ON DELTAT, MEANING IT COULD BE DUE TO NUMERICAL NOISE.

MODULE geometrical_constants
  INTEGER, PARAMETER :: iterations = 1024
  INTEGER, PARAMETER :: frames     = 1
  INTEGER, PARAMETER :: D          = 3
  INTEGER, PARAMETER :: N_x        = 128
  INTEGER, PARAMETER :: N_y        = 128
  INTEGER, PARAMETER :: N_z        = 128
  INTEGER, PARAMETER :: NN         = 2097152

  DOUBLE PRECISION, PARAMETER:: BOX_SIZE_X = 100E-6
  DOUBLE PRECISION, PARAMETER:: BOX_SIZE_Y = 100E-6
  DOUBLE PRECISION, PARAMETER:: BOX_SIZE_Z = 200E-6
  DOUBLE PRECISION, PARAMETER:: length_scale = 10.0E-6
END MODULE geometrical_constants

MODULE physical_constants
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER:: e       = 1.602176462E-19
  DOUBLE PRECISION, PARAMETER:: hbar    = 1.054E-34
  DOUBLE PRECISION, PARAMETER:: mu_B    = 9.2740154E-24
  DOUBLE PRECISION, PARAMETER:: k_B     = 1.38E-23
  DOUBLE PRECISION, PARAMETER:: pi      = 3.1415926
  DOUBLE PRECISION, PARAMETER:: mu_cero = 12.566370614E-7
  DOUBLE PRECISION, PARAMETER:: epsilon_cero = 8.854187817E-12
  DOUBLE PRECISION, PARAMETER:: deltat  = 0.5E-3
END MODULE physical_constants

MODULE atomic_properties
  USE physical_constants
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: mass = 87.0*1.6604E-27
  DOUBLE PRECISION, PARAMETER :: a    = 5.2E-9
  DOUBLE PRECISION, PARAMETER :: mu   = 9.2740154E-24
  DOUBLE PRECISION, PARAMETER :: N_0  = 2.5E5
  DOUBLE PRECISION            :: U_0  = 4.0*pi*hbar*hbar*a*N_0/mass
END MODULE atomic_properties

MODULE potential_constants
  USE physical_constants
  DOUBLE PRECISION, PARAMETER :: w_x  = 2.0*pi*200.0
  DOUBLE PRECISION, PARAMETER :: w_y  = 2.0*pi*200.0
  DOUBLE PRECISION, PARAMETER :: w_z  = 2.0*pi*200.0  
  DOUBLE PRECISION, PARAMETER :: width  = 1.0
  DOUBLE PRECISION, PARAMETER :: height = 1.0
  DOUBLE PRECISION            :: mu_TF
  INTEGER :: channels
END MODULE potential_constants


program tridiagonal      
  
  USE geometrical_constants
  USE physical_constants
  USE atomic_properties
  USE potential_constants
  
  IMPLICIT NONE
  
  !------------------- I. a.  PHYSICAL PARAMETERS -------------------------------- 
  
  INTEGER          :: N(3)
  DOUBLE PRECISION :: DELTA(3),BOX_SIZE(3),w(3)
  
  !------------------- I. b. SCALAR ARGUMENTS ----------------------------------
  
  INTEGER          :: info,nrhs,ldb,r,p,s,i,j,k,l,m,t, N_med(3)
  CHARACTER*1      :: trans
  COMPLEX*16       :: alpha,beta(3)
  DOUBLE PRECISION :: gamma,x_i,x_j,y_k,z_l,norma
  
  !------------------- I. c. ARRAY ARGUMENTS -----------------------------------
  
  COMPLEX*16 :: diagM(NN),diagU(NN-1),diagL(NN-1),diagU2(NN-2)
  COMPLEX*16 :: b(NN,1),x_aux(NN,1),V_ext(NN)
  INTEGER, DIMENSION(NN) :: ipiv
  
  N(1) = N_x
  N(2) = N_y
  N(3) = N_z
  
  N_med(1) = INT(0.5*N(1))
  N_med(2) = INT(0.5*N(2))
  N_med(3) = INT(0.5*N(3))
  
  info=0
  nrhs=1
  ldb = NN
  trans = "N" 
  
  w(1) = w_x
  w(2) = w_y
  w(3) = w_z
  
  BOX_SIZE(1) = BOX_SIZE_X 
  BOX_SIZE(2) = BOX_SIZE_Y
  BOX_SIZE(3) = BOX_SIZE_Z
  
  DELTA(1) =  BOX_SIZE(1)/(N(1)-1)
  DELTA(2) =  BOX_SIZE(2)/(N(2)-1)
  DELTA(3) =  BOX_SIZE(3)/(N(3)-1)
  
  beta(1)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(1)**2)),0.)
  beta(2)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(2)**2)),0.)
  beta(3)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(3)**2)),0.)
  alpha    = DCMPLX(-deltat/(2*D*hbar),0.)

  DO channels=0,0

  CALL EXTERNAL_POTENTIAL(V_ext,delta,w,BOX_SIZE,N)
  CALL INITIAL_STATE(x_aux,delta,w,BOX_SIZE,N)	
  CALL NORMALIZATION(x_aux,delta)
  CALL WRITE_FUN(x_aux,delta,BOX_SIZE)

!  CALL i_TIME_EVOLUTION(N,delta,V_ext,x_aux,BOX_SIZE)
  
  DO r = 1,iterations
  DO s = 1,3
         
     CALL RHS(N,alpha,beta,V_ext,x_aux,b,s)		
     DO p=1,2
        
        CALL MATRIX(N,alpha,beta,V_ext,x_aux,diagL,diagU,diagM,s)                   	  	  	
        CALL ZGTTRF(NN,diagL,diagM,diagU,diagU2,ipiv,info)  
        x_aux = b
        CALL ZGTTRS(trans,NN,nrhs,diagL,diagM,diagU,diagU2,ipiv,x_aux,ldb,info)		   
        CALL NORMALIZATION(x_aux,delta)
        
     ENDDO  
     
     CALL ROTATION_1(N,V_ext,s)		   	  
     CALL ROTATION_2(N,x_aux,s)
     
  END DO
  
!   IF(MOD(r,int(iterations/frames)).EQ.0.) THEN  
!      CALL WRITE_FUN(x_aux,delta,BOX_SIZE)
!   END IF
  
  END DO

  CALL WRITE_FUN(x_aux,delta,BOX_SIZE)

  END DO
  write(6,*)'----- Programa para la evolucion temporal de un BEC-2D -----'
  write(6,*)'-------- Nottingham, March 16th 2007, PhD studies ----------'
  write(6,*)'------------------ ppxgs@nottingham.ac.uk ------------------'
  
201 format(4E15.6)
202 format(5E15.6)
  
END PROGRAM tridiagonal


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


SUBROUTINE INITIAL_STATE(x,delta,w,BOX_SIZE,N)
  
  USE geometrical_constants
  USE physical_constants
  USE atomic_properties
  USE potential_constants
 
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: N(3)
  DOUBLE PRECISION, INTENT (IN)  :: delta(3),w(3),BOX_SIZE(3)
  COMPLEX*16, INTENT (OUT) :: x(NN,1)  

  INTEGER          :: i,j,k,l, N_med(3)
  DOUBLE PRECISION :: x_j,y_k,z_l,w_bar,a_bar,partereal
  DOUBLE PRECISION :: sigma2(3),length_scale_(3)

  !open(4,file="init_test_2.dat",action="read") 
  open(2,file="ShellSmallWeak_128x128x128.dat",action="read")
  
  N_med(3) = INT(0.5*N(3))
 
  length_scale_(1) = SQRT(60.0*hbar/(mass*w(1)))
  length_scale_(2) = SQRT(1.0*hbar/(mass*w(2)))
  length_scale_(3) = SQRT(1.0*hbar/(mass*w(3)))
  
  sigma2(1) =   length_scale_(1)*length_scale_(1)
  sigma2(2) =   length_scale_(2)*length_scale_(2)
  sigma2(3) =   length_scale_(3)*length_scale_(3)
  
  w_bar = (w(1)*w(2)*w(3))**(0.333333)
  a_bar = sqrt(hbar/(mass*w_bar))
  mu_TF = 0.5*(15*N_0*a/a_bar)**(0.6)*hbar*w_bar 

  read(2,*) mu_TF
  x = 0.0

  DO i=1,NN
     j = MOD(i-1,N(1)) + 1  
     l = (i-1)/(N(1)*N(2)) + 1
     k = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1
     
     x_j =  -0.5*BOX_SIZE(1) + (j-1)*DELTA(1)
     y_k =  -0.5*BOX_SIZE(2) + (k-1)*DELTA(2)  
     z_l =  -0.5*BOX_SIZE(3) + (l-1)*DELTA(3)
     
  !   x(i,1) = DCMPLX(((pi*pi*pi*sigma2(1)*sigma2(2)*sigma2(3))**(-0.25))*exp(-0.5*((x_j)**2)/sigma2(1)),0.0) 
  !   x(i,1) = x(i,1)*DCMPLX(exp(- 0.5*((y_k)**2)/sigma2(2) - 0.5*((z_l)**2)/sigma2(3)),0.0)

     read(2,*) partereal
     IF(partereal.LT.mu_TF) THEN
       x(i,1) = sqrt( (1.0/U_0)*(mu_TF - partereal)*((hbar*hbar/(mass*length_scale*length_scale))) )
     END IF	

  ENDDO
  CLOSE(2)
  !write(*,*) BOX_SIZE(1),BOX_SIZE(2),BOX_SIZE(3),DELTA(1),DELTA(2),DELTA(3),length_scale(1),length_scale(2),length_scale(3),sigma2(1),sigma2(2),sigma2(3) 
201 format(1E15.6)
  
END SUBROUTINE INITIAL_STATE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


SUBROUTINE EXTERNAL_POTENTIAL(V_ext,delta,w,BOX_SIZE,N)
  USE physical_constants
  USE geometrical_constants
  USE potential_constants
  USE atomic_properties
	
  IMPLICIT NONE
  INTEGER,INTENT (IN) :: N(3)
  DOUBLE PRECISION, INTENT (IN)  :: delta(3),w(3),BOX_SIZE(3)
  COMPLEX*16, INTENT (OUT) :: V_ext(NN)

  INTEGER          :: i,j,k,l
  DOUBLE PRECISION :: x_j,y_k,z_l,energy

  open(2,file="ShellSmallWeak_128x128x128.dat",action="read")
  

 ! read(2,201) mu_TF
  read(2,*) mu_TF
 ! mu_TF = mu_TF * ((hbar*hbar/(mass*lenght_scale)))

  !------------- FREE EXPANSTION ---------------
  !DO i=1,N
  !	V_ext(i) = 0.0
  !END DO
  !------------- --------------- ---------------
  
  !------------- HARMONIC POTENTIAL ---------------
  DO i=1,NN
     j = MOD(i-1,N(1)) + 1  
     l = (i-1)/(N(1)*N(2)) + 1
     k = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1
     
      x_j =  -0.5*BOX_SIZE(1) + (j-1)*DELTA(1)
      y_k =  -0.5*BOX_SIZE(2) + (k-1)*DELTA(2)  
      z_l =  -0.5*BOX_SIZE(3) + (l-1)*DELTA(3)
      
!      V_ext(i) = DCMPLX(0.5*mass*((w(1)*x_j)**2  + (w(2)*y_k)**2 + (w(3)*z_l)**2),0.0)   

!      read(2,201) energy
!	write(*,*)i,j,k,l
     read(2,*) energy
      V_ext(i) = DCMPLX(energy,0.0)*((hbar*hbar/(mass*length_scale*length_scale)))
  END DO
  CLOSE(2)
  !201 format(1E15.6)  
END SUBROUTINE EXTERNAL_POTENTIAL


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

 
 SUBROUTINE RHS(N,alpha,beta,V_ext,x_aux,b,s)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   INTEGER,INTENT (IN)     :: s,N(3)
   COMPLEX*16, INTENT (IN) :: alpha,beta(3)
   COMPLEX*16, INTENT (IN) :: V_ext(NN),x_aux(NN,1)
   COMPLEX*16, INTENT (OUT) :: b(NN,1)
   
   INTEGER :: i,j,k,l,i_prima_r,i_prima_l,s_,i_, i_prima
   
   DO i=1,NN	

     
!      IF(s.EQ.1) THEN
!         j = MOD(i-1,N(1)) + 1  
!!$         l = (i-1)/(N(1)*N(2)) + 1
!!$         k = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1
!!$         s_ = s
!!$         i_ = j     
!!$         i_prima_r =  j + 1 + (k+1)*N(1)   + (l-1)*N(1)*N(2)
!!$         i_prima_l =  j - 1 + (k-1)*N(1) + (l-1)*N(1)*N(2)
!!$      END IF
!!$      
!!$      IF(s.EQ.2) THEN
!!$         k = MOD(i-1,N(2)) + 1  
!!$         j = (i-1)/(N(2)*N(3)) + 1
!!$         l = (i-k-(j-1)*(N(2)*N(3)))/N(2) + 1
!!$         s_ = s
!!$         i_ = k   
!!$         i_prima_r =   k + 1 + (l-1)*N(2) + (j-1)*N(2)*N(3)
!!$         i_prima_l =   k - 1 + (l-1)*N(2) + (j-1)*N(2)*N(3)
!!$      END IF
!!$            
!!$      IF(s.EQ.3) THEN
!!$         l = MOD(i-1,N(3)) + 1  
!!$         k = (i-1)/(N(3)*N(1)) + 1
!!$         j = (i-l-(k-1)*(N(1)*N(3)))/N(3) + 1
!!$         s_ = s
!!$         i_ = l     
!!$         i_prima_r =  l + 1 + (j-1)*N(3) + (k-1)*N(3)*N(1)
!!$         i_prima_l =  l - 1 + (j-1)*N(3) + (k-1)*N(3)*N(1)
!!$      END IF
				

!      i_prima =  N(s_)

      s_ = s
      i_prima = 1
      b(i,1) = x_aux(i,1)*(1. - 2*beta(s_) + alpha*(V_ext(i) + U_0*x_aux(i,1)*CONJG(x_aux(i,1))))     			
      IF(MOD(i,N(s_)).EQ.1) b(i,1) = b(i,1) + beta(s_)*x_aux(i+i_prima,1)
      IF(MOD(i,N(s_)).EQ.0) b(i,1) = b(i,1) + beta(s_)*x_aux(i-i_prima,1)
      IF((MOD(i,N(s_)).NE.0).AND.(MOD(i,N(s_)).NE.1)) b(i,1) = b(i,1) + beta(s_)*x_aux(i+i_prima,1)+ beta(s_)*x_aux(i-i_prima,1)	
 
!      IF(i_.EQ.1)     b(i,1) = b(i,1) + beta(s_)*x_aux(i_prima_r,1)
!      IF(i_.EQ.N(s_)) b(i,1) = b(i,1) + beta(s_)*x_aux(i_prima_l,1)
!      IF((i_.NE.0).AND.(i_.NE.N(s_))) b(i,1) = b(i,1) + beta(s_)*x_aux(i_prima_r,1)+ beta(s_)*x_aux(i_prima_l,1)	

   END DO
   
 END SUBROUTINE RHS


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


 SUBROUTINE MATRIX(N,alpha,beta,V_ext,x_aux,diagL,diagU,diagM,s)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   INTEGER,INTENT (IN)     :: s,N(3)
   COMPLEX*16, INTENT (IN) :: alpha,beta(3)
   COMPLEX*16, INTENT (IN) :: V_ext(NN),x_aux(NN,1)
   COMPLEX*16, INTENT (OUT) :: diagL(NN-1),diagU(NN-1),diagM(NN)
   
   INTEGER :: i,j,k,l,i_prima
   
   DO i=1,NN-1
      IF (MOD(i,N(s)).EQ.0) THEN
         diagL(i) = 0. 
         diagU(i) = 0.
      ELSE
         diagL(i) = -beta(s)
         diagU(i) = -beta(s)
      ENDIF
   ENDDO
   
   DO i=1,NN			 
      i_prima = i		
      diagM(i)= 1. + 2.0*beta(s) - alpha*(V_ext(i_prima) + U_0*x_aux(i_prima,1)*CONJG(x_aux(i_prima,1)))
   ENDDO
   
 END SUBROUTINE MATRIX
 

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

 SUBROUTINE NORMALIZATION(x,delta)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT (IN) :: delta(3)
   COMPLEX*16, INTENT (INOUT) :: x(NN,1)
	
   DOUBLE PRECISION :: norma
   INTEGER :: i
   
   norma = 0.
   DO i=1,NN
      norma = norma + x(i,1)*CONJG(x(i,1))*delta(1)*delta(2)*delta(3)
   ENDDO

   x = x/sqrt(norma)   
   
!   write(*,*) norma
   
 END SUBROUTINE NORMALIZATION
 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


SUBROUTINE ROTATION_1(N,V,s)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: N(3),s
   COMPLEX*16, INTENT (INOUT) :: V(NN)
   
   INTEGER :: i,j,k,l,i_prima
   COMPLEX*16 :: x(NN,1)
   
   do i = 1,NN
      
      IF(s.EQ.1) THEN
         j = MOD(i-1,N(1)) + 1  
         l = (i-1)/(N(1)*N(2)) + 1
         k = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1
         i_prima = (j-1)*N(3)*N(2) + (l-1)*N(2) + k
      END IF
      
      IF(s.EQ.2) THEN
         k = MOD(i-1,N(2)) + 1  
         j = (i-1)/(N(2)*N(3)) + 1
         l = (i-k-(j-1)*(N(2)*N(3)))/N(2) + 1
         i_prima = (k-1)*N(3)*N(1) + (j-1)*N(3) + l
      END IF
            
      IF(s.EQ.3) THEN
         l = MOD(i-1,N(3)) + 1  
         k = (i-1)/(N(3)*N(1)) + 1
         j = (i-l-(k-1)*(N(1)*N(3)))/N(3) + 1
         i_prima = (l-1)*N(1)*N(2) + (k-1)*N(1) + j
      END IF
 
      x(i_prima,1) = V(i) 
      
   enddo
   
   do i = 1,NN
      V(i) = x(i,1)
   enddo
   
 END SUBROUTINE ROTATION_1


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

 
 SUBROUTINE ROTATION_2(N,x_aux,s)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: N(3),s
   COMPLEX*16, INTENT (INOUT) :: x_aux(NN,1)
   COMPLEX*16 :: x(NN,1)
   
   INTEGER :: i,j,k,l,i_prima
   
   DO i = 1,NN
   
   IF(s.EQ.1) THEN
         j = MOD(i-1,N(1)) + 1  
         l = (i-1)/(N(1)*N(2)) + 1
         k = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1
         i_prima = (j-1)*N(3)*N(2) + (l-1)*N(2) + k
      END IF
      
      IF(s.EQ.2) THEN
         k = MOD(i-1,N(2)) + 1  
         j = (i-1)/(N(2)*N(3)) + 1
         l = (i-k-(j-1)*(N(2)*N(3)))/N(2) + 1
         i_prima = (k-1)*N(3)*N(1) + (j-1)*N(3) + l
      END IF
            
      IF(s.EQ.3) THEN
         l = MOD(i-1,N(3)) + 1  
         k = (i-1)/(N(3)*N(1)) + 1
         j = (i-l-(k-1)*(N(1)*N(3)))/N(3) + 1
         i_prima = (l-1)*N(1)*N(2) + (k-1)*N(1) + j
      END IF
 

      
      !write(*,*) s,i,j,k,l,i_prima
  
      x(i_prima,1) = x_aux(i,1) 
      
   ENDDO
   
   DO i = 1,NN
      x_aux(i,1) = x(i,1)
  END DO
  
  
END SUBROUTINE ROTATION_2


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


SUBROUTINE WRITE_FUN(x_aux,delta,BOX_SIZE)
  USE physical_constants
  USE geometrical_constants
  USE atomic_properties
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT (IN) :: delta(3),BOX_SIZE(3)
  COMPLEX*16, INTENT (IN) :: x_aux(NN,1)
  INTEGER :: i,j,k,l
  DOUBLE PRECISION :: density_1D(N_x),density_2D(N_x,N_y),density_2D_sec(N_x,N_y)
  DOUBLE PRECISION :: x_j,y_k
  
  open(2,file="gs1D_Shell_SW.dat",access='append')
  open(3,file="gs3D_Shell_SW.dat",access='append')
  
  DO i = 1,NN						  
     write(3,201) REAL(x_aux(i,1))
  END DO
  
  do j=1,N_x
	density_1D(j) = 0.0
  end do 
  

  DO l=1,N_z
  DO k=1,N_y
  DO j=1,N_x
     i = j + (k-1)*N_x + (l-1)*N_x*N_y
     density_1D(j) = density_1D(j) + real(x_aux(i,1)*conjg(x_aux(i,1)))
  END DO
  END DO
  END DO

  DO j = 1,N_x	
     x_j =  -0.5*BOX_SIZE(1) + (j-1)*DELTA(1)	
     density_1D(j) = density_1D(j)*DELTA(2)*DELTA(3)
     write(2,202) 1E6*x_j,N_0*density_1D(j)/1E6
  END DO

  write(2,*)  
  write(2,*)

  DO l=1,N_z
  DO k=1,N_y
  DO j=1,N_x
     i = j + (k-1)*N_x + (l-1)*N_x*N_y
     density_2D(j,k) = density_2D(j,k) + real(x_aux(i,1)*conjg(x_aux(i,1)))
     IF(l.eq.N_z/2) density_2D_sec(j,k) = real(x_aux(i,1)*conjg(x_aux(i,1)))
  END DO
  END DO
  END DO

  DO k=1,N_y
  y_k =  -0.5*BOX_SIZE(2) + (k-1)*DELTA(2)	
  DO j = 1,N_x	
     x_j =  -0.5*BOX_SIZE(1) + (j-1)*DELTA(1)	
     density_2D(j,k) = density_2D(j,k)*DELTA(3)
     write(2,204) 1E6*x_j,1E6*y_k,N_0*density_2D(j,k)/1E12,N_0*density_2D_sec(j,k)/1E12
  END DO
  WRITE(2,*)
  END DO

  write(2,*)  
  write(2,*)

  close(2)
  close(3)


201 format(1E15.6)
202 format(2E15.6)
203 format(3E15.6)
204 format(4E15.6)
  
END SUBROUTINE WRITE_FUN


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


SUBROUTINE WRITE_FUN_SEC(N,N_med,ax,delta,BOX_SIZE)
  
  USE physical_constants
  USE geometrical_constants
  USE atomic_properties

  INTEGER, INTENT (IN) :: N_med(3),N(3)
  DOUBLE PRECISION, INTENT (IN) :: delta(3),BOX_SIZE(3)
  COMPLEX*16, INTENT (IN) ::ax(NN,1)
  
  INTEGER :: i,j,k,l
  DOUBLE PRECISION :: x_j,y_k,z_l
  
  open(2,file="gs1.dat", position = "append")
  
  DO i = 1,NN
 	
     j = MOD(i-1,N(1)) + 1  
     l = (i-1)/(N(1)*N(2)) + 1
     k = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1       
     
     x_j =  -0.5*BOX_SIZE(1) + (j-1)*delta(1)
     y_k =  -0.5*BOX_SIZE(2) + (k-1)*delta(2)  
     z_l =  -0.5*BOX_SIZE(3) + (l-1)*delta(3)
     
     IF(l.EQ.N_med(3)) THEN						
     !   write(2,201) 1E6*x_j,1E6*y_k,1E6*z_l,REAL(ax(i,1))
        write(*,*) i,REAL(ax(i,1))
     ENDIF
     IF(j.EQ.N(1).AND.l.EQ.N_med(3)) write(2,*)     
  END DO
  write(2,*) 
  write(2,*)
  
201 format(4E15.6)
  
END SUBROUTINE WRITE_FUN_SEC


!!$SUBROUTINE i_TIME_EVOLUTION(N,delta,V_ext,x_aux,BOX_SIZE)
!!$
!!$  USE geometrical_constants
!!$  USE physical_constants
!!$  USE atomic_properties
!!$  USE potential_constants
!!$  
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER, INTENT(IN)          :: N(3)
!!$  DOUBLE PRECISION,INTENT (IN) :: DELTA(3),BOX_SIZE(3)
!!$  COMPLEX*16, INTENT (IN)      :: x_aux(NN,1),V_ext(NN)
!!$
!!$  
!!$  INTEGER          :: info,nrhs,ldb,r,p,s,i,j,k,l,m,t, N_med(3)
!!$  CHARACTER*1      :: trans
!!$  DOUBLE PRECISION :: gamma,x_i,x_j,y_k,z_l,norma,w(3)
!!$  COMPLEX*16       :: alpha,beta(3)
!!$  
!!$  
!!$  COMPLEX*16 :: diagM(NN),diagU(NN-1),diagL(NN-1),diagU2(NN-2)
!!$  COMPLEX*16 :: b(NN,1)
!!$  INTEGER, DIMENSION(NN) :: ipiv
!!$
!!$  N_med(1) = INT(0.5*N(1)+1)
!!$  N_med(2) = INT(0.5*N(2)+1)
!!$  N_med(3) = INT(0.5*N(3)+1)
!!$  
!!$  info=0
!!$  nrhs=1
!!$  ldb = NN
!!$  trans = "N" 
!!$  
!!$  w(1) = w_x
!!$  w(2) = w_y
!!$  w(3) = w_z
!!$
!!$  beta(1)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(1)**2)),0.)
!!$  beta(2)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(2)**2)),0.)
!!$  beta(3)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(3)**2)),0.)
!!$  alpha    = DCMPLX(-deltat/(2*D*hbar),0.)
!!$  
!!$  DO r = 1,iterations
!!$  DO s = 1,3
!!$         
!!$     CALL RHS(N,alpha,beta,V_ext,x_aux,b,s)		
!!$     !CALL WRITE_FUN_SEC(N,N_med,b,delta,BOX_SIZE)
!!$     DO p=1,2
!!$        
!!$        CALL MATRIX(N,alpha,beta,V_ext,x_aux,diagL,diagU,diagM,s)                   	  	  	
!!$        !CALL TRIDIAGONAL_SOLUTION(diagL,diagM,diagU,b,x_aux)
!!$        CALL ZGTTRF(NN,diagL,diagM,diagU,diagU2,ipiv,info)  
!!$        x_aux = b
!!$        CALL ZGTTRS(trans,NN,nrhs,diagL,diagM,diagU,diagU2,ipiv,x_aux,ldb,info)		   
!!$        CALL NORMALIZATION(x_aux,delta)
!!$        
!!$     ENDDO   ! ENDDO FOR HALF AND WHOLE TIME STEP IN ONE DIMENSION
!!$     
!!$     CALL ROTATION_1(N,V_ext,s)		   	  
!!$     CALL ROTATION_2(N,x_aux,s)
!!$     
!!$  END DO ! ENDDO FOR THE three DIMENSIONS
!!$  
!!$  IF(MOD(r,int(iterations/frames)).EQ.0.) THEN  
!!$     !CALL WRITE_FUN_SEC(N,N_med,x_aux,delta,BOX_SIZE)  
!!$     CALL WRITE_FUN(x_aux,delta,BOX_SIZE)
!!$  END IF
!!$  
!!$END DO !ENDDO FOR WHOLE TIME EVOLUTION
!!$
!!$!CALL WRITE_FUN(x_aux,delta,BOX_SIZE)
!!$!CALL WRITE_FUN_SEC(N,N_med,x_aux,delta,BOX_SIZE)
!!$!CALL WRITE_nD(x_aux,2)
!!$
!!$write(6,*)'----- Programa para la evolucion temporal de un BEC-2D -----'
!!$write(6,*)'-------- Nottingham, March 16th 2007, PhD studies ----------'
!!$write(6,*)'------------------ ppxgs@nottingham.ac.uk ------------------'
!!$
!!$END SUBROUTINE i_TIME_EVOLUTION

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


! SUBROUTINE TRIDIAGONAL_SOLUTION(a,b,c,d,x)
!   !CALL TRIDIAGONAL_SOLUTION(diagL,diagM,diagU,b,x_aux)
!   USE physical_constants
!   USE geometrical_constants
!   USE atomic_properties
! 
!   COMPLEX*16, INTENT (IN)    :: a(NN-1),b(NN)
!   COMPLEX*16, INTENT (INOUT) :: c(NN-1),d(NN),x(NN,1)
! 
!   INTEGER i
!   COMPLEX*16 id
! 
!   c(1) = c(1)/b(1)
!   d(1) = d(1)/b(1)
! 
!   DO i=2,NN,1
!      id   = 1.0/(b(i) - c(i-1)*a(i-1))
!      c(i) = c(i)*id
!      d(i) = (d(i) - a(i-1)*d(i-1))*id
!   END DO
!   
!   x(NN,1) = d(NN)
! 
!   DO i = NN-1,-1,1
!      x(i) = d(i) - c(i)*x(i+1)
!   END DO

! END SUBROUTINE TRIDIAGONAL_SOLUTION
