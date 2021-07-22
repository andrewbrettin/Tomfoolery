! Program for generating random numbers
program generate_random_numbers
    IMPLICIT NONE

    ! Declarations
    integer, parameter :: NIGLOBAL = 10, NJGLOBAL = 8
    real :: x1(NIGLOBAL, NJGLOBAL) ! Generate grids of uniform random vars
    real :: x2(NIGLOBAL, NJGLOBAL)
    real :: z1(NIGLOBAL, NJGLOBAL)
    real :: z2(NIGLOBAL, NJGLOBAL)
    real :: s(NIGLOBAL, NJGLOBAL)
    

    real, parameter :: PI = 4 * atan(1.0)

    call random_number(x1)
    call random_number(x2)
    
    ! Implement the Box-Muller Transformation to get standard normal vars
    z1 = SQRT(-2 * LOG(x1)) * COS(2 * PI * x2)
    z2 = SQRT(-2 * LOG(x1)) * SIN(2 * PI *x2)

    s = z1 + z2

    print *, s(1:10, 1)

end program generate_random_numbers