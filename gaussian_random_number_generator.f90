subroutine gaussian_random_number_generator( mean, standard_deviation, &
                                             the_random_nbr_we_re_looking_for )

implicit none

! two random numbers with uniform distribution
double precision :: random_number1 , random_number2

! mean of the gaussian distribution we want a random number from
double precision, intent(in) :: mean

! standard_deviation of the gaussian distribution
double precision, intent(in) :: standard_deviation

double precision, intent ( out ) :: the_random_nbr_we_re_looking_for

double precision, parameter :: pi = 3.14159265358979323846D+00


!-----------------------------------------------------------------------------------------------

! first step, get two random numbers with uniform distribution.
! This is done easily with the fortran 90 intrinsic function random_number

call random_number( random_number1 )
call random_number( random_number2 )

! convert to normal distribution, which is the purpose of this code :)

the_random_nbr_we_re_looking_for = &
       mean + &
       standard_deviation * &
       sqrt ( -2.0d0 * log ( random_number1 ) ) * &
       cos ( 2.0d0 * pi * random_number2 )

end subroutine gaussian_random_number_generator
