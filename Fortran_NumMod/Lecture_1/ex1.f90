!========================================================
!   gfortran ex1.f90 -o ex1 -lblas -llapack && ./ex1 
!   ./ex1
!========================================================
PROGRAM ylesanded_nr_1

    IMPLICIT NONE

    CHARACTER(len=8)   :: title_A, title_b, title_C, title_D
    CHARACTER(len=8)   :: title_2AT, title_invA
    !INTEGER, parameter :: dp = selected_real_kind(15, 307)
    INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)
    INTEGER            :: sizev, size1, size2
    REAL(KIND=SP)      :: x, y, z, w 
    REAL(KIND=SP)      :: u, v, ans
    REAL(KIND=SP), DIMENSION(:,:), ALLOCATABLE :: A, B, C, D, twoAtrans, inv_A
    REAL(KIND=SP), DIMENSION(:),   ALLOCATABLE :: vec_b
    REAL(KIND=SP), EXTERNAL :: func_y, func_z, F
    
    WRITE(*, "(A)") "Exercise 1."
    x = 13.0
    y =  7.0
    z = SQRT(x)
    w = 6.0*x*y - x/y
    WRITE(*,"(A, F8.4)") "z = ", z
    WRITE(*,"(A, F8.4)") "w = ", w

    WRITE(*, "(A)") "Exercise 2."
    x =  2.0;    y = -1.0;    z =  5.0
    !x = -6.0;    y =  5.0;    z = -2.0
    u = (x*y + z/(x-4.0*y*y) + z*z) / (x*x*x + y*y - 3.0*z*x)
    v = (x*y*z + x*x - (2*y)**2)**(5.0/(x+y))
    WRITE(*,"(A, ES11.4)") "u = ", u
    WRITE(*,"(A, ES11.4)") "v = ", v

    WRITE(*, "(A)") "Exercise 3."
    !x = -1.0
    !x = 4.0
    x = 3.0
    
    y = func_y(x)
    WRITE(*,"(A, F3.1, A, ES10.4)") "y(", x, ") = ", y

    WRITE(*, "(A)") "Exercise 4."
    !u = 2.0
    u = 8.0
    z = func_z(u)
    WRITE(*,"(A, F3.1, A, ES11.4)") "z(", u, ") = ", z

    WRITE(*, "(A)") "Exercise 5."
    x =  50.0
    y = -30.0
    ans = F(x, y)
    WRITE(*,"(A, F4.1, A, F5.1, A, ES10.4)") "z(", x, ",", y, ") = ", ans

    WRITE(*, "(A)") "Exercise 6."
    size1 = 3
    size2 = 3
    sizev = 3
    ALLOCATE(A(size1, size2))
    ALLOCATE(vec_b(sizev))
    ALLOCATE(twoAtrans(size2, size1))
    ALLOCATE(inv_A(size1, size2))
    ALLOCATE(B(size1, size2))
    ALLOCATE(C(size1, size2))
    ALLOCATE(D(size1, size2))

    A = TRANSPOSE(RESHAPE(        &
                [3.0, 12.0, 52.0, &
                 4.0, 6.0, -11.0, &
                 -2.0, 7.0, 2.0], &
                 [size1, size2]))
    vec_b = [13.0, -2.0, 5.0]

    title_A = "Matrix A"
    CALL print_matrix(A, title_A)
    
    title_b = "Vector b"
    CALL print_vector(vec_b, title_b)   ! technically the same as title_B

    twoAtrans = 2.0*A
    title_2AT = "Mat 2A^T"
    CALL print_matrix(twoAtrans, title_2AT)

    ! BLAS/LAPACK boosted inverse function utilizing LU-decomp
    inv_A = inv(A)
    title_invA = "Mat A^-1"
    CALL print_matrix(inv_A, title_invA)

    B = twoAtrans + inv_A
    title_B = "Matrix B"    ! notice how we don't need to define it separately
    CALL print_matrix(B, title_B)

    ! BLAS/LAPACK boosted matmul()
    C = matfunc(A, B)
    title_C = "Matrix C"
    CALL print_matrix(C, title_C)

    ! Element-wise multiplication
    D = A * B
    title_D = "Matrix D"
    CALL print_matrix(D, title_D)
    
    ! Free allocated memory
    DEALLOCATE(A)
    DEALLOCATE(vec_b)
    DEALLOCATE(twoAtrans)
    DEALLOCATE(inv_A)
    DEALLOCATE(B)
    DEALLOCATE(C)
    DEALLOCATE(D)

    ! Start of functions which require EXPLICIT definition   
    CONTAINS

!=============================================================================
    
    SUBROUTINE print_matrix(matrix, title)
        IMPLICIT NONE
        !INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)
        CHARACTER(len=*),INTENT(IN), OPTIONAL :: title   
        REAL(KIND=SP), DIMENSION(:,:)         :: matrix
        !INTENT(IN) - information for the programmer and the 
        ! compiler that this is an INPUT; OPTIONAL - not ne-
        ! cessary for the user to input it to WORK correctly.
        CHARACTER(len=15) :: mformat = '(100(F14.6,1x))'
        INTEGER           :: i, j
    
        WRITE(*, *) title
        DO i = 1, MIN(5, SIZE(matrix, 1))
            WRITE(*, mformat) (matrix(i, j), j = 1, MIN(5, SIZE(matrix, 2)))
        END DO
    END SUBROUTINE print_matrix

!-----------------------------------------------------------------------------

    SUBROUTINE print_vector(vector, title)
        IMPLICIT NONE
        REAL(KIND=SP), DIMENSION(:),INTENT(IN) :: vector  
        CHARACTER(len=*),INTENT(IN), OPTIONAL  :: title  ! Optional title
        INTEGER :: i
        CHARACTER(len=15) :: format = '(F14.6, 1X)'  ! Print format for each element
    
        ! If title is provided, print it
        IF (PRESENT(title)) THEN
            WRITE(*, *) title
        END IF
    
        ! Loop through the vector and print the elements
        DO i = 1, MIN(5, SIZE(vector))
            WRITE(*, format) vector(i)
        END DO
    END SUBROUTINE print_vector

!-----------------------------------------------------------------------------

    FUNCTION inv(MAT) RESULT(MAT_inv)
        IMPLICIT NONE
        
        REAL(KIND=SP), DIMENSION(:, :), INTENT(IN)  :: MAT  ! Input matrix A
        REAL(KIND=SP), DIMENSION(:, :), ALLOCATABLE :: MAT_inv  ! Output matrix (inverse of A)
        REAL(KIND=SP), ALLOCATABLE, DIMENSION(:) :: WORK
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV  ! Array for pivoting during LU decomposition
        INTEGER :: INFO, LDA, LWORK                 ! Error status for LAPACK routines
        INTEGER :: ROWS, COLS
        
        ROWS = SIZE(MAT, 1)
        COLS = SIZE(MAT, 2)
        
        ! Ensure the matrix is square
        IF (ROWS /= COLS) THEN
            PRINT *, 'Matrix is not square. Cannot compute inverse.'
            STOP
        END IF

        LDA = ROWS
        LWORK = ROWS*COLS
        
        ! Allocate the inverse matrix and IPIV array
        ALLOCATE(MAT_inv(LDA, COLS))
        ALLOCATE(WORK(LWORK))
        ALLOCATE(IPIV(ROWS))
        
        ! Copy the input matrix MAT into MAT_inv
        MAT_inv = MAT
        
        ! Perform LU decomposition of MAT (MAT = L*U)
        CALL SGETRF(ROWS, COLS, MAT_inv, LDA, IPIV, INFO)
        IF (INFO .NE. 0) THEN
            PRINT *, 'Error in LU decomposition'
            DEALLOCATE(MAT_inv, IPIV)
            STOP
        END IF
        
        ! Compute the inverse of MAT using the LU decomposition
        CALL SGETRI(ROWS, MAT_inv, COLS, IPIV, WORK, LWORK, INFO)
        IF (INFO .NE. 0) THEN
            PRINT *, 'Error in computing the inverse of MAT'
            DEALLOCATE(MAT_inv, IPIV)
            STOP
        END IF
        
        ! DeALLOCATE the IPIV and WORK array
        DEALLOCATE(IPIV)
        DEALLOCATE(WORK)
        
    END FUNCTION inv

!-----------------------------------------------------------------------------

    FUNCTION matfunc(matA, matB) result(matC)
        REAL(KIND=SP), dimension(:,:),INTENT(IN)   :: matA, matB    !INTENT(IN) - information for the programmer and the compiler that this is an INPUT
        REAL(KIND=SP), allocatable, dimension(:,:) :: matC
        REAL(KIND=SP) :: ALPHA, BETA
        INTEGER       :: M, N, K

        ! Get dimensions of matrices
        M = SIZE(matA, 1)
        K = SIZE(matA, 2)
        N = SIZE(matB, 2)

        ! Set default values
        ALPHA = 1.0_sp
        BETA  = 0.0_sp

        ! Check dimensions
        IF (SIZE(matB, 1) /= K) then
            WRITE(*,*) 'Error: Matrix dimensions do not match for multiplication.'
            matC = 0.0_sp
            RETURN
        END IF

        ! Perform matrix multiplication
        ALLOCATE(matC(M, N))
        CALL SGEMM('N', 'N', M, N, K, ALPHA, matA, M, matB, K, BETA, matC, M)

        ! You can use the resulting matrix C as needed
    END FUNCTION matfunc
    
END PROGRAM ylesanded_nr_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FORTRAN 77 style of writing a function
FUNCTION func_y(x)
    IMPLICIT NONE
    INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)
    REAL(KIND=SP)      :: func_y
    REAL(KIND=SP)      :: x
    func_y = (x**2 - 3.0)*(2.0 + x)**4 - 5.0*exp(x) + 2.0*cos(x + 1.0)
    RETURN  ! Optional, but not mandatory any more
END FUNCTION func_y

! FORTRAN 90 style of writing a function
FUNCTION func_z(u) RESULT(z)
    IMPLICIT NONE
    INTEGER, parameter :: sp = selected_real_kind(6, 37)
    REAL(KIND=SP), INTENT(IN) :: u ! input
    REAL(KIND=SP)             :: z ! output
    z = u**5 - 3.0_sp*u**4 + u**3 + 1.0_sp
END FUNCTION func_z

FUNCTION F(x,y) RESULT(ans)
    IMPLICIT NONE
    INTEGER, parameter :: sp = selected_real_kind(6, 37)
    REAL(KIND=SP), INTENT(IN) :: x, y
    REAL(KIND=SP)             :: ans ! output
    ans = sin(x-y)/(x*x) + (cos(2.0_sp*x + y)/(x-y)**4)**(1.0_sp/3.0_sp)
END FUNCTION F

