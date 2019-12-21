

MODULE geometrical_constants
  INTEGER, PARAMETER :: iterations = 10000
  INTEGER, PARAMETER :: frames     = 1
  INTEGER, PARAMETER :: D          = 3
  INTEGER, PARAMETER :: N_x        = 400
  INTEGER, PARAMETER :: N_y        = 40
  INTEGER, PARAMETER :: N_z        = 40
  INTEGER, PARAMETER :: NN         = 640000!432000

  DOUBLE PRECISION, PARAMETER:: BOX_SIZE_X = 200.0E-6
  DOUBLE PRECISION, PARAMETER:: BOX_SIZE_Y = 20.0E-6
  DOUBLE PRECISION, PARAMETER:: BOX_SIZE_Z = 20.0E-6
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
  DOUBLE PRECISION, PARAMETER:: deltat  = 0.5E-7
END MODULE physical_constants

MODULE atomic_properties
  USE physical_constants
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: mass = 87.0*1.6604E-27
  DOUBLE PRECISION, PARAMETER :: a    = 5.4E-9
  DOUBLE PRECISION, PARAMETER :: mu   = 9.2740154E-24
  DOUBLE PRECISION, PARAMETER :: N_0  = 60.0E3
  DOUBLE PRECISION            :: U_0  = 4.0*pi*hbar*hbar*a*N_0/mass
END MODULE atomic_properties

MODULE potential_constants
  USE physical_constants
  DOUBLE PRECISION, PARAMETER :: w_x  = 2.0*pi*15.0
  DOUBLE PRECISION, PARAMETER :: w_y  = 2.0*pi*200.0
  DOUBLE PRECISION, PARAMETER :: w_z  = 2.0*pi*200.0    
  DOUBLE PRECISION, PARAMETER :: width  = 1.0
  DOUBLE PRECISION, PARAMETER :: height = 1.0
  INTEGER :: channels
END MODULE potential_constants
!--------------------------------------------------------------------------------------------------------------------------
!s_m: large vector of how many elements each procesor has to send to each other processor (send matric) 0:0, 0:1.... 1:0,1:1,1:2 ....
!s_c: vector of how many elements each procesor have to send (sent counter): i:0,i:1,i:2....
!r_m: large vector of how many elements each procesor has to recieve to (recieve matrix): 0:0, 0:1.... 1:0,1:1,1:2 ....
!r_c: vector of how many elements each procesor recieve from (recive count): i:0,i:1,i:2....
!s_i: index used to create the message.
!old_proc: The origin of the element i
!new_proc: Where to sent the element i
!i_prima:  new index of the i-th element under the new arragment (global)
!i_message:  new index of the i-th element under the new arragment (local).
!i_abs: The index corresponding to the recieved message(?)!
!i_abs_: i_abs repartido entre los procesadores
!i_abs_sent: i_abs_ reorganizado para enviar
!i_final: Es el indice correspondiente al mensaje que se recibe, para reorganizar lo que se recibe.
!--------------------------------------------------------------------------------------------------------------------------
program CN_PARALELO

  USE geometrical_constants
  USE physical_constants
  USE atomic_properties
  USE potential_constants
  IMPLICIT NONE
  INCLUDE 'mpif.h'

!--------- VARIABLES AND CONSTANTS FOR MPI -------------

  integer MPI_Status(MPI_STATUS_SIZE)
  integer myproc, numprocs, ierr

  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: x
  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: x_,x_sent,x_recv

  COMPLEX*16, DIMENSION(:),ALLOCATABLE :: v_ext
  COMPLEX*16, DIMENSION(:),ALLOCATABLE :: v_ext_,v_ext_sent,v_ext_recv

  INTEGER, DIMENSION(:),ALLOCATABLE     :: s_m,r_m
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: s_c,r_c
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: s_d,r_d
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: s_i,i_final
  INTEGER, DIMENSION(:),ALLOCATABLE     :: i_abs_sent
  INTEGER, DIMENSION(:), ALLOCATABLE    :: i_prima,i_message,i_abs,i_abs_,index_
  INTEGER old_proc,new_proc
  INTEGER i,j,k,l,s,i_,s_
  INTEGER n_p_p

  CHARACTER*40 file_name
  CHARACTER*10 sub

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

!-----------VARIABLES AND CONSTANTS FOR TIME EVOLUTION-------------------

  !------------------- I. a.  PHYSICAL PARAMETERS -------------------------------- 
  
  INTEGER          :: N(3)
  DOUBLE PRECISION :: DELTA(3),BOX_SIZE(3),w(3)
  
  !------------------- I. b. SCALAR ARGUMENTS ----------------------------------
  
  INTEGER          :: info,nrhs,ldb,r,p,m,t, N_med(3)
  CHARACTER*1      :: trans
  COMPLEX*16       :: alpha,beta(3)
  COMPLEX*16	   :: gamma,x_i,x_j,y_k,z_l
  DOUBLE PRECISION :: norma_l,norma
  
  !------------------- I. c. ARRAY ARGUMENTS -----------------------------------

  COMPLEX*16, DIMENSION(:),  ALLOCATABLE :: diagM,diagU,diagL,diagU2
  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: b
  INTEGER,DIMENSION(:),ALLOCATABLE       :: ipiv
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!------------------------------------------------------------------------

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myproc, ierr)

  n_p_p = NN/numprocs


  ALLOCATE(s_c(3,numprocs),r_c(3,numprocs),s_d(3,numprocs),r_d(3,numprocs),s_i(3,n_p_p),i_final(3,n_p_p))
  
  ALLOCATE(i_abs_(n_p_p),i_abs_sent(n_p_p))

  IF(myproc.eq.0) THEN
     ALLOCATE(s_m(numprocs*numprocs),r_m(numprocs*numprocs))
     ALLOCATE(i_prima(NN))
     ALLOCATE(i_message(NN))
     ALLOCATE(i_abs(NN))     
     ALLOCATE(index_(NN))
  END IF
  
  write(sub,'(i2)') myproc
  file_name ='test_'//adjustl(trim(sub))
  file_name = adjustl(file_name)

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

  beta(1)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(1)**2)),0.0)
  beta(2)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(2)**2)),0.0)
  beta(3)  = DCMPLX(0.25*deltat*hbar/(mass*(DELTA(3)**2)),0.0)
  alpha    = DCMPLX(-deltat/(2*D*hbar),0.0)

! !   IF(myproc.eq.0) THEN
! ! ! 	open(111,file="sorted_index.dat",action="read")  
! !   	open(110,file="sorted_index.dat",action="write")
! !   END IF	
  
  DO s=1,3
     IF(myproc.eq.0) then
        !CALL init_parallel_com(numprocs,NN,s,N)
        CALL init_parallel_com(s_m,r_m,i_prima,i_abs,i_message,numprocs,NN,s,N)
! !        DO i=1,NN
! !           write(110,105) s,i,i_prima(i),i_abs(i),i_message(i)
! !        END DO 
! !        write(110,*) 
! !        write(110,*)
! !        DO i=1,numprocs*numprocs
! !           write(110,104) s,i,s_m(i),r_m(i)
! !        END DO
!!       DO i_=1,NN
!!          read(111,105) s_,i,i_prima(i_),i_abs(i_),i_message(i_)
!!       END DO
!!       DO i_=1,numprocs*numprocs
!!          read(111,104) s_,i,s_m(i_),r_m(i_)
!!       END DO
            !write(*,*) (s_m(i), i=1,numprocs*numprocs)
           !write(*,*) (r_m(i), i=1,numprocs*numprocs) 
     END IF
     !write(*,*)  n_p_p,numprocs

      CALL mpi_scatter(i_message, n_p_p,MPI_INTEGER,s_i(s,:),n_p_p,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL mpi_scatter(s_m,numprocs,MPI_INTEGER,s_c(s,:),numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL mpi_scatter(r_m,numprocs,MPI_INTEGER,r_c(s,:),numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL mpi_scatter(i_abs,n_p_p,MPI_INTEGER,i_abs_,n_p_p,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL init_parallel_local(i_abs_,s_i,s_c,r_c,i_abs_sent,s_d,r_d,s,numprocs,n_p_p)
      CALL MPI_AlltoallV(i_abs_sent,s_c(s,:),s_d(s,:),MPI_INTEGER,i_final(s,:),r_c(s,:),r_d(s,:), &
           & MPI_INTEGER,MPI_COMM_WORLD,ierr)
  END DO
  
  IF(myproc.eq.0) THEN
!     close(111)
!     close(110)
     DEALLOCATE(s_m,r_m,i_message,i_prima,i_abs)
  END IF

  DEALLOCATE(i_abs_,i_abs_sent)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  

  DO channels=0,0
     
     IF(myproc.eq.0) THEN
        ALLOCATE(x(NN,1))
        ALLOCATE(v_ext(NN))
        CALL INITIAL_STATE(x,delta,w,BOX_SIZE,N)	
        CALL EXTERNAL_POTENTIAL(v_ext,delta,w,BOX_SIZE,N)
        CALL NORMALIZATION(norma_l,x,delta,NN)
        x = x/sqrt(norma_l)	
     END IF
     
     ALLOCATE(x_(n_p_p,1))
     ALLOCATE(v_ext_(n_p_p))
     
     CALL mpi_scatter(x(:,1),n_p_p,MPI_DOUBLE_COMPLEX,x_(:,1),n_p_p,MPI_DOUBLE_COMPLEX,0, &
          & MPI_COMM_WORLD,ierr)
     CALL mpi_scatter(v_ext,n_p_p,MPI_DOUBLE_COMPLEX,v_ext_,n_p_p,MPI_DOUBLE_COMPLEX,0, &
          & MPI_COMM_WORLD,ierr)
     IF(myproc.eq.0) DEALLOCATE(x,v_ext)
     CALL NORMALIZATION(norma_l,x_,delta,n_p_p)
     CALL MPI_REDUCE(norma_l,norma,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Bcast(norma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
     x_ = x_/sqrt(norma)	
     
     ALLOCATE(x_sent(n_p_p,1))  
     ALLOCATE(x_recv(n_p_p,1))  
     ALLOCATE(v_ext_sent(n_p_p))  
     ALLOCATE(v_ext_recv(n_p_p))  
     
     ALLOCATE(diagM(n_p_p),diagU(n_p_p-1),diagL(n_p_p-1),diagU2(n_p_p-2))
     ALLOCATE(b(n_p_p,1))
     ALLOCATE(ipiv(n_p_p))
     
     !CALL WRITE_FUN(x_,delta,BOX_SIZE,N,n_p_p,file_name)
     !CALL WRITE_FUN(v_ext_,delta,BOX_SIZE,N,n_p_p,file_name)
     
     DO r = 1,1000,1!iterations
        DO s = 1,3	        
           CALL RHS(N,alpha,beta,v_ext_,x_,b,s,n_p_p)		
           DO p=1,2
              CALL MATRIX(N,alpha,beta,v_ext_,x_,diagL,diagU,diagM,s,n_p_p)                   	  	  	
              CALL ZGTTRF(n_p_p,diagL,diagM,diagU,diagU2,ipiv,info)  
              x_ = b
              CALL ZGTTRS(trans,n_p_p,nrhs,diagL,diagM,diagU,diagU2,ipiv,x_,ldb,info)		   
              CALL NORMALIZATION(norma_l,x_,delta,n_p_p)
              CALL MPI_REDUCE(norma_l,norma,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
              CALL MPI_Bcast(norma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
              x_ = x_/sqrt(norma)	
           ENDDO   ! ENDDO FOR HALF AND WHOLE TIME STEP IN ONE DIMENSION        
           DO i=1,n_p_p
              x_sent(i,1) = x_(s_i(s,i),1)        ! Creating the message to sent
           END DO
           CALL MPI_AlltoallV(x_sent(:,1),s_c(s,:),s_d(s,:),MPI_DOUBLE_COMPLEX,x_recv(:,1),r_c(s,:),r_d(s,:), &
                & MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)        
           DO i = 1,n_p_p
              x_(i_final(s,i),1) = x_recv(i,1)   ! Reorganizing the recieved buffer
           END DO
           DO i=1,n_p_p
              x_sent(i,1) = v_ext_(s_i(s,i))        ! Creating the message to sent
           END DO
           CALL MPI_AlltoallV(x_sent(:,1),s_c(s,:),s_d(s,:),MPI_DOUBLE_COMPLEX,x_recv(:,1),r_c(s,:),r_d(s,:), &
                & MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)     
           DO i = 1,n_p_p
              v_ext_(i_final(s,i)) = x_recv(i,1)   ! Reorganizing the recieved buffer
           END DO
        END DO
     END DO
     CALL WRITE_FUN(x_,delta,BOX_SIZE,N,n_p_p,file_name)
     DEALLOCATE(x_sent,x_recv,v_ext_sent,v_ext_recv,diagM,diagU,diagL,diagU2,b,ipiv)
     
  END DO
  
  CALL mpi_finalize(ierr)
  
  
101 format(1i8)
103 format(3i8)
104 format(4i8)
105 format(5i8)
106 format(6i8)
109 format(9i8)
  
203 format(3E15.6)
201 format(1E15.6)

3021 format(2I5,1E15.6)
3022 format(2I5,2E15.6)
3041 format(4i4,1E15.6)
3074 format(7i4,4E15.6)
  
30121 format(1i4,2E15.6,1i4)
  
END program CN_PARALELO

! SUBROUTINE SORTING(i_prima,s_m,n_p_p,index,NN,N,numprocs,s)
!   
!   INTEGER, INTENT (IN) :: n_p_p,NN,numprocs,s,N(3)
!   INTEGER, DIMENSION(numprocs*numprocs) :: s_m
!   INTEGER, DIMENSION(NN),INTENT(IN)  :: i_prima
!   INTEGER, DIMENSION(NN),INTENT(OUT) :: index
! 
!   INTEGER, DIMENSION(n_p_p) :: c
!   INTEGER i,j,l
! 
!   INTEGER i_,sort_(numprocs),sort(numprocs)
!   INTEGER, DIMENSION(NN) :: indice
!     
! 
!   DO l=0,numprocs-1
! 
!      DO j=1,n_p_p
!         c(j) = i_prima(j+l*n_p_p)
! !        write(*,*) j,j+l*n_p_p,c(j)
!      END DO
! 
!      DO j=2,n_p_p
!         y = c(j)
!         i = j
!         DO WHILE((i.GT.1) .AND. (c(i-1).GT.y)) 
!            c(i) = c(i-1)
!            i = i-1           
!         END DO
!         c(i) = y        
!      END DO 
!  
!      DO j=1,n_p_p
!         i = 1
!         DO WHILE ((i.LE.n_p_p) .AND.(c(j).NE.i_prima(i+l*n_p_p)))
!            i = i + 1           
!         END DO
!         index(j+l*n_p_p) = i  
!      END DO
!  
!   END DO
! 
! 
! END SUBROUTINE SORTING
! 

SUBROUTINE init_parallel_com(s_m,r_m,i_prima,i_abs,i_message,numprocs,NN,s,N)
	
    IMPLICIT NONE
    INTEGER, DIMENSION(numprocs*numprocs), INTENT (OUT)  :: s_m,r_m
    INTEGER, DIMENSION(NN), INTENT (OUT)  :: i_prima,i_abs,i_message
    INTEGER, DIMENSION(3),INTENT (IN)     :: N
    INTEGER, INTENT (IN)                  :: s,NN,numprocs



    INTEGER old_proc,new_proc
    INTEGER i,j,k,l
    INTEGER n_p_p,myproc
    INTEGER new_p,n_x,n_y,n_z
  

    n_p_p = INT(NN/numprocs)
    
    !DO i=1,numprocs
    !   DO j=1,numprocs
    !      s_m((i-1)*numprocs+j) = 0
    !   END DO
    !END DO
    s_m = 0
    
     DO i=1,NN

       IF(s.EQ.1) THEN
           j  = MOD(i-1,N(1)) + 1
           l  = (i-1)/(N(1)*N(2)) + 1
           k  = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1
           i_prima(i) =  k + (l-1)*N(2) + (j-1)*N(2)*N(3)
        END IF

        IF(s.EQ.2) THEN
           k  = MOD(i-1,N(2)) + 1
           j  = (i-1)/(N(2)*N(3)) + 1
           l  = (i-k-(j-1)*(N(2)*N(3)))/N(2) + 1
           i_prima(i) =  l + (j-1)*N(3) + (k-1)*N(3)*N(1)
        END IF

       IF(s.EQ.3) THEN
           l  = MOD(i-1,N(3)) + 1
           k  = (i-1)/(N(3)*N(1)) + 1
           j  = (i-l-(k-1)*(N(3)*N(1)))/N(3) + 1
           i_prima(i) =  j + (k-1)*N(1) + (l-1)*N(1)*N(2)
       END IF

       
        old_proc = int((i-1)/n_p_p)
        new_proc = int((i_prima(i)-1)/n_p_p)
        s_m((old_proc+1-1)*numprocs + new_proc + 1) = s_m((old_proc+1-1)*numprocs + new_proc + 1) + 1
        i_abs(i) = i_prima(i) - INT((i_prima(i)-1)/n_p_p)*n_p_p

     END DO
 
     DO i=1,numprocs
        DO j=1,numprocs
           r_m((j-1)*numprocs +i) = s_m((i-1)*numprocs +j)
        END DO
     END DO
     
!     !CALL sorting(i_prima,s_m,n_p_p,i_message,NN,N,numprocs,s)	
     n_x = N(1)/numprocs
     n_y = N(2)/numprocs
     n_z = N(3)/numprocs
     
     IF(s.EQ.1) THEN
        DO myproc=0,numprocs-1
           i = 1 + myproc*n_p_p
           DO new_p=0,numprocs-1
             DO j=1+new_p*n_x,(new_p+1)*n_x
                DO l=1+myproc*n_z,(myproc+1)*n_z
                   DO k=1,N(2)                      	             
                      i_message(i) = j + (k-1)*N(1)+(l-1)*N(1)*N(2) - myproc*n_p_p
                      i = i+1
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF
    
    IF(s.EQ.2) THEN
       DO myproc=0,numprocs-1
          i = 1 + myproc*n_p_p
          DO new_p=0,numprocs-1
             DO k=1+new_p*n_y,(new_p+1)*n_y
                DO j=1+myproc*n_x,(myproc+1)*n_x
                   DO l=1,N(3)                      
                      i_message(i) = k + (l-1)*N(2)+(j-1)*N(2)*N(3) - myproc*n_p_p
                      i = i+1
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF  

    IF(s.EQ.3) THEN
       DO myproc=0,numprocs-1
          i = 1 + myproc*n_p_p
          DO new_p=0,numprocs-1
             DO l=1+new_p*n_z,(new_p+1)*n_z
                DO k=1+myproc*n_y,(myproc+1)*n_y
                   DO j=1,N(1)                      
                      i_message(i) = l + (j-1)*N(3)+(k-1)*N(3)*N(1) - myproc*n_p_p
                      i = i+1
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF


END SUBROUTINE init_parallel_com
 
SUBROUTINE init_parallel_local(i_abs_,s_i,s_c,r_c,i_abs_sent,s_d,r_d,s,numprocs,n_p_p)

    INTEGER, DIMENSION(n_p_p), INTENT (IN)   :: i_abs_
    INTEGER, DIMENSION(3,n_p_p), INTENT (IN) :: s_i
    INTEGER, DIMENSION(3,numprocs), INTENT (IN)  :: s_c,r_c

    INTEGER, DIMENSION(n_p_p), INTENT (OUT)  :: i_abs_sent
    INTEGER, DIMENSION(3,numprocs), INTENT (OUT)  :: s_d,r_d
  
    INTEGER, INTENT (IN) :: s,numprocs,n_p_p

    INTEGER i

    DO i=1,n_p_p
       i_abs_sent(i) = i_abs_(s_i(s,i))
    END DO
 
    r_d(s,1) = 0
    s_d(s,1) = 0
 
    DO i=2,numprocs
        s_d(s,i) = s_d(s,i-1) + s_c(s,i-1)
        r_d(s,i) = r_d(s,i-1) + r_c(s,i-1)
     END DO

END SUBROUTINE

SUBROUTINE INITIAL_STATE(x,delta,w,BOX_SIZE,N)

  USE geometrical_constants
  USE physical_constants
  USE atomic_properties

  IMPLICIT NONE
  INTEGER,INTENT (IN) :: N(3)
  DOUBLE PRECISION, INTENT (IN)  :: delta(3),w(3),BOX_SIZE(3)
  COMPLEX*16, INTENT (OUT) :: x(NN,1)

  INTEGER          :: i,j,k,l,i_
  DOUBLE PRECISION :: x_j,y_k,z_l,w_bar,a_bar,mu_TF,partereal
  DOUBLE PRECISION :: sigma2(3),length_scale(3)

!  open(2,file="init_gs_QPC_2.dat",action="read")

  length_scale(1) = SQRT(1.0*hbar/(mass*w(1)))
  length_scale(2) = SQRT(1.0*hbar/(mass*w(2)))
  length_scale(3) = SQRT(1.0*hbar/(mass*w(3)))

  sigma2(1) =   length_scale(1)*length_scale(1)
  sigma2(2) =   length_scale(2)*length_scale(2)
  sigma2(3) =   length_scale(3)*length_scale(3)

  DO i=1,NN

     j = MOD(i-1,N(1)) + 1  
     l = (i-1)/(N(1)*N(2)) + 1
     k = (i-j-(l-1)*(N(1)*N(2)))/N(1) + 1
     
     x_j =  -0.5*BOX_SIZE(1) + (j-1)*DELTA(1)
     y_k =  -0.5*BOX_SIZE(2) + (k-1)*DELTA(2)  
     z_l =  -0.5*BOX_SIZE(3) + (l-1)*DELTA(3)
     
     x(i,1) = 1.0!DCMPLX(((pi*pi*pi*sigma2(1)*sigma2(2)*sigma2(3))**(-0.25))*exp(-0.5*((x_j)**2)/sigma2(1)),0.0) 
     !x(i,1) = x(i,1)*DCMPLX(exp(- 0.5*((y_k)**2)/sigma2(2) - 0.5*((z_l)**2)/sigma2(3)),0.0)
!!$!
!!$!    read(2,201) partereal
!!$!    x(i,1)  = DCMPLX(partereal,0.0)
  ENDDO 
!  CLOSE(2)
	
!201 format(1E15.6)

END SUBROUTINE

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
      
      V_ext(i) = DCMPLX(0.5*mass*((w(1)*x_j)**2  + (w(2)*y_k)**2 + (w(3)*z_l)**2),0.0) 
  END DO 
  !------------- --------------- ---------------
  
 
  !------------- EXTERNAL POTENTIAL ---------------
!!$!  open(4,file="energy_QPC_2.dat",action="read")
!!$!  DO i=1,NN
!!$!     read(4,201) energy
!!$!     V_ext(i) = DCMPLX(energy,0.0)
!!$!  END DO
!!$!  CLOSE(4)
  !------------- --------------- ---------------
  
201 format(1E15.6)
END SUBROUTINE EXTERNAL_POTENTIAL

SUBROUTINE RHS(N,alpha,beta,V_ext,x_aux,b,s,n_p_p)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   INTEGER,INTENT (IN)     :: s,N(3),n_p_p
   COMPLEX*16, INTENT (IN) :: alpha,beta(3)
   COMPLEX*16, INTENT (IN) :: V_ext(n_p_p),x_aux(n_p_p,1)
   COMPLEX*16, INTENT (OUT) :: b(n_p_p,1)
   
   INTEGER :: i,j,k,l,i_prima_r,i_prima_l,s_,i_, i_prima
   
   DO i=1,n_p_p

      s_ = s
      i_prima = 1
      b(i,1) = x_aux(i,1)*(1. - 2*beta(s_) + alpha*(V_ext(i) + U_0*x_aux(i,1)*CONJG(x_aux(i,1))))     			
      IF(MOD(i,N(s_)).EQ.1) b(i,1) = b(i,1) + beta(s_)*x_aux(i+i_prima,1)
      IF(MOD(i,N(s_)).EQ.0) b(i,1) = b(i,1) + beta(s_)*x_aux(i-i_prima,1)
      IF((MOD(i,N(s_)).NE.0).AND.(MOD(i,N(s_)).NE.1)) b(i,1) = b(i,1) + beta(s_)*x_aux(i+i_prima,1)+ beta(s_)*x_aux(i-i_prima,1)	
 
   END DO

 END SUBROUTINE RHS
! 
SUBROUTINE MATRIX(N,alpha,beta,V_ext,x_aux,diagL,diagU,diagM,s,n_p_p)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   INTEGER,INTENT (IN)     :: s,N(3),n_p_p
   COMPLEX*16, INTENT (IN) :: alpha,beta(3)
   COMPLEX*16, INTENT (IN) :: V_ext(n_p_p),x_aux(n_p_p,1)
   COMPLEX*16, INTENT (OUT) :: diagL(n_p_p-1),diagU(n_p_p-1),diagM(n_p_p)
   
   INTEGER :: i,j,k,l,i_prima
   
   DO i=1,n_p_p-1
      IF (MOD(i,N(s)).EQ.0) THEN
         diagL(i) = 0. 
         diagU(i) = 0.
      ELSE
         diagL(i) = -beta(s)
         diagU(i) = -beta(s)
      ENDIF
   ENDDO
   
   DO i=1,n_p_p			 
      i_prima = i		
      diagM(i)= 1. + 2.0*beta(s) - alpha*(V_ext(i_prima) + U_0*x_aux(i_prima,1)*CONJG(x_aux(i_prima,1)))
   ENDDO
   
 END SUBROUTINE MATRIX

SUBROUTINE NORMALIZATION(norma_l,x,delta,n_p_p)
   
   USE physical_constants
   USE geometrical_constants
   USE atomic_properties
   
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: n_p_p   
   DOUBLE PRECISION, INTENT (IN) :: delta(3)
   COMPLEX*16, INTENT (IN) :: x(n_p_p,1)
	
   DOUBLE PRECISION, INTENT (OUT) :: norma_l	
 
   INTEGER :: i
   
   norma_l = 0.
   DO i=1,n_p_p
      norma_l = norma_l + REAL(x(i,1)*CONJG(x(i,1))*delta(1)*delta(2)*delta(3))
   ENDDO
      
 END SUBROUTINE NORMALIZATION

SUBROUTINE WRITE_FUN(x_aux,delta,BOX_SIZE,N,num_datos,file_name)
  
  DOUBLE PRECISION, INTENT (IN) :: delta(3),BOX_SIZE(3)
  COMPLEX*16, INTENT (IN) :: x_aux(num_datos,1)
  INTEGER, INTENT (IN)    :: N(3),num_datos
  CHARACTER*40,INTENT (IN) :: file_name
  INTEGER :: i
  
  write(*,*) file_name,num_datos
  open(2,file=file_name, position = "append")
  DO i = 1,num_datos						  
      write(2,201) REAL(x_aux(i,1))
  END DO
  close(2)

201 format(1E15.6)
  
END SUBROUTINE WRITE_FUN







