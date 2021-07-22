! This program is based on a tutorial from: https://www.youtube.com/watch?v=__2UgFNYgf8.


program fortran_tut
    implicit none ! Explicit variable declaration
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !  Variable declarations ! !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real                                :: rand
    integer                             :: n = 1
    real, parameter                     :: PI = 3.14159
    real, dimension(1:50)               :: arr_1
    integer, dimension(5,5)             :: arr_2
    integer, dimension(:), allocatable  :: a5
    integer                             :: numvals = 0
    integer, dimension(1:5)             :: a6 = (/1,2,3,4,5/)
    integer, dimension(1:3,1:3)         :: a7
    integer :: i, j
    
    ! Structures
    type Customer
       character (len=40) :: name
       integer :: age
       real :: balance
    end type Customer

    type(Customer) :: cust1
    cust1%name    = "Sally Wally"
    cust1%age     = 24
    cust1%balance = 133.25

    



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ! ! ! ! Statements ! ! ! ! ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, mod(5,4)
    print *, (5**4)

    call random_number(rand)
    print *, rand

    ! Functions
    print *, "Square root: ", sqrt(2.0)
    print *, "PI", 4*atan(1.0)
    print *, "Exp: ", exp(1.0)
    print *, "Trig: ", cos(PI), cosh(1.0)
    print *, "Logarithm: ", log(2.71)

    ! Loops: loop from n = 1 to 10 (inclusive), skip by 2
    print *, "Loops"
    do n = 1, 10
        if (mod(n,2) == 0) then
            print *, n
        end if
    end do

    ! Arrays
    do i = 1, 5
        do j = 1, 5
            call random_number(rand)
            arr_2(i,j) = rand
        end do
    end do

    print *, "functions",  get_sum2(3, 4)

    ! Functions
    contains
        integer function get_sum(n1, n2)
            implicit NONE
            integer :: n1, n2, sum
            sum = n1 + n2
        end function get_sum

        function get_sum2(n1, n2) result(sum)
            implicit NONE
            integer, intent(in) :: n1, n2
            integer :: sum
            sum = n1 + n2
        end function get_sum2

end program fortran_tut