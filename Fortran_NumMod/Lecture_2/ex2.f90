!========================================================
!   gfortran ex2.f90 -o ex2 
!   ./ex2
!========================================================
program main
    implicit none
    real           :: x, y
    real, external :: f
    
    x = -1.0
    !x = 4.0
    y = f(x)
    write(*,"(a, f4.1, a, f7.4)") "y(", x, ") = ", y

end program main

function f(x)
    implicit none
    real :: f, x
    f = 3.0*x**4 - x**3 + 5.0*x**2 + 7.0
    return
end function f