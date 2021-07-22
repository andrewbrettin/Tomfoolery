
program fortran_variables
    implicit none
    ! Variable declarations
    real, parameter :: PI = 3.1
    real :: r_num= 0.0, r_num2 = 0.0
    double precision :: dbl_num = 1.111111111d+0
    integer :: int1 = 0
    logical :: can_vote = .TRUE.
    character (len=10) :: month
    complex :: compnum = (2.0, 4.9)
    print *, "Biggest real: ", huge(r_num)



end program fortran_variables