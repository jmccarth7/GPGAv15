subroutine median_calc( median )


! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module


implicit none

real(kind=8), allocatable, dimension(:)  :: temp_SSE

real(kind=8) :: median

integer(kind=4) :: k1, k2


!----------------------------------------------------------------------------------


! calculate the median of the individual_SSE

! first, sort the individual_SSE in ascending order
! then the median value is the value at (or near) index =  n_GA_individuals/2

allocate( temp_SSE( n_GA_individuals ) )

temp_SSE = individual_SSE

call sort( n_GA_individuals, temp_SSE )

if( mod( n_GA_individuals,2 ) == 0 ) then

    k1 = n_GA_individuals / 2
    k2 = k1 + 1
    median = 0.5d0 * ( temp_SSE(k1) + temp_SSE(k2) )

else

    k1 = ( n_GA_individuals + 1 ) / 2
    median = temp_SSE(k1)

endif

write(GA_print_unit,'(/A,1x,E24.16)') 'medcalc: median = ', median

deallocate( temp_SSE )

return


end subroutine median_calc

!234567890123456789012345678901234567890123456789012345678901234567890
