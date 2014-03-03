subroutine sse0_calc()


use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module
use GP_data_module
use GA_variables_module

implicit none


real(kind=8) :: fvec(n_time_steps)

integer(kind=4) :: i_CODE_equation
integer(kind=4) :: i_time_step

!-------------------------------------------------------------------------------


! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.




write(GP_print_unit,*) ' '

do  i_CODE_equation=1,n_CODE_equations
    write(GP_print_unit,'(A,1x,I6,2(1x,E15.7))') &
          'ssec: i_eqn, data_variance, data_variance_inv ', &
                 i_CODE_equation, data_variance(i_CODE_equation), &
                                  data_variance_inv(i_CODE_equation)
enddo !  i_CODE_equation

write(GP_print_unit,'(/A/)') 'ssec: using data_variance inv   '




SSE0 = 0.0D+0
do  i_time_step=1,n_time_steps

    fvec(i_time_step)=0.

    !write(GP_print_unit,'(/A,1x,I6, 1x,I10)') &
    !'ssec: myid, i_time_step ', myid, i_time_step

    do  i_CODE_equation=1,n_CODE_equations


        !write(GP_print_unit,'(A,2(1x,I6), 3(1x,E15.7))') &
        !      'ssec: myid, i_eqn, RK_soln, data_array, var ', &
        !            myid, i_CODE_equation,                    &
        !            Runge_Kutta_Solution(i_time_step,i_CODE_equation), &
        !            Data_Array(i_time_step,i_CODE_equation), &
        !            data_variance(i_CODE_equation)

        !write(GP_print_unit,'(A,2(1x,I6), 1x,E15.7)') &
        !      'ssec: myid, i_eqn, data_variance ', &
        !            myid, i_CODE_equation, data_variance(i_CODE_equation)

        fvec(i_time_step) = fvec(i_time_step)  +                   &
             Data_Array(i_time_step,i_CODE_equation)**2  *         &
                                Data_Variance_inv(i_CODE_equation)

        !fvec(i_time_step) = fvec(i_time_step)  +                 &
        !     Data_Array(i_time_step,i_CODE_equation)**2  /       &
        !                        Data_Variance(i_CODE_equation)

    enddo ! i_CODE_equation


    SSE0 = SSE0 + fvec(i_time_step)

    !write(GP_print_unit,'(A,1x,I6, 1x,I6, 1x, E15.7)')&
    !      'ssec: myid, i_time_step, fvec ', &
    !             myid, i_time_step, fvec(i_time_step)

enddo ! i_time_step


write(GP_print_unit,'(A,1x,I6,2x,E15.7/)') 'ssec: myid, SSE0 = ',myid, SSE0

return
end subroutine sse0_calc
