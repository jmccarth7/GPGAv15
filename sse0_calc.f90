subroutine sse0_calc()


use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module
use GP_data_module
use GA_variables_module

implicit none


real (kind=8) :: fvec(n_time_steps)

real (kind=8) :: ssum, ssum2, totobs, cff

!-------------------------------------------------------------------------------


! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.

do i_CODE_equation=1,n_CODE_equations

  ssum  = 0.0D+0
  ssum2 = 0.0D+0

  do  i_time_step=0,n_time_steps
      ssum  = ssum  +  Data_Array(i_time_step,i_CODE_equation)
      ssum2 = ssum2 +  Data_Array(i_time_step,i_CODE_equation)**2

      !write(6,*) 'ssec: i_time_step, Data_Array(i_time_step,1:n_CODE_equations) ', &
      !                  i_time_step, Data_Array(i_time_step,1:n_CODE_equations) 
  enddo ! i_time_step

  totobs=dble(n_time_steps+1)

  cff=( (totobs*ssum2)-(ssum**2) ) / ( totobs*(totobs-1.0D+0) )

  if( cff .gt. 0.0D+0) then  ! set variance to observed variance for normalize by the s.d.
      Data_Variance(i_CODE_equation)=cff
  else ! set variance to 1.0 for normalization to be 'unaltered'
      Data_Variance(i_CODE_equation)=1.0D+0
  endif !   cff .gt. 0.0D+0

enddo !  i_CODE_equation



write(6,*) ' '

do  i_CODE_equation=1,n_CODE_equations
    write(6,'(A,1x,I6, 1x,E15.7)') &
          'ssec: i_eqn, data_variance ', &
                 i_CODE_equation, data_variance(i_CODE_equation)
enddo !  i_CODE_equation

write(6,*) ' '




SSE0 = 0.0D+0
do i_time_step=1,n_time_steps

  fvec(i_time_step)=0.

  !write(6,'(/A,1x,I6, 1x,I10)')'ssec: myid, i_time_step ', myid, i_time_step

  do i_CODE_equation=1,n_CODE_equations

    if( abs(  data_variance(i_CODE_equation) ) > 1.0D-30  )then

        !write(6,'(A,2(1x,I6), 3(1x,E15.7))') &
        !      'ssec: myid, i_eqn, RK_soln, data_array, var ', &
        !            myid, i_CODE_equation,                    &
        !            Runge_Kutta_Solution(i_time_step,i_CODE_equation), &
        !            Data_Array(i_time_step,i_CODE_equation), &
        !            data_variance(i_CODE_equation)

        !write(6,'(A,2(1x,I6), 1x,E15.7)') &
        !      'ssec: myid, i_eqn, data_variance ', &
        !            myid, i_CODE_equation, data_variance(i_CODE_equation)

        fvec(i_time_step) = fvec(i_time_step)  +                                  &
             Data_Array(i_time_step,i_CODE_equation)**2  /     &
                                Data_Variance(i_CODE_equation)

    else

        !write(6,'(/A,1x,I6, 1x,I10)')'ssec: bad variance myid, i_time_step ', &
        !                                                 myid, i_time_step
        !write(6,'(A,2(1x,I6), 1x,E15.7)') &
        !      'ssec: myid, i_CODE_equation, data_variance(i_CODE_equation)              ', &
        !             myid, i_CODE_equation, data_variance(i_CODE_equation)

    endif ! abs(  data_variance(i_CODE_equation) ) > 1.0D-30



  enddo ! i_CODE_equation


  SSE0 = SSE0 + fvec(i_time_step)

  !write(6,'(A,1x,I6, 1x,I6, 1x, E15.7)')&
  !      'ssec: myid, i_time_step, fvec ', &
  !             myid, i_time_step, fvec(i_time_step)

enddo ! i_time_step

!write(6,'(A,1x,I6,2x,E24.16)') 'ssec: myid, SSE0 = ',myid, SSE0

return
end subroutine sse0_calc