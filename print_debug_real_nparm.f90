subroutine print_debug_real_nparm( label, input_array  )



! print REAL arrays of the form:

!  input_array(1:n_Maximum_Number_Parameters, 1:n_GP_Individuals)



use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module

implicit none


character(*) :: label


integer(kind=4) :: i_GP_Individual


integer(kind=4) :: i_parm


real(kind=8), dimension(1:n_Maximum_Number_Parameters, 1:n_GP_Individuals) :: &
                         input_array
!--------------------------------------------------------------------------------

write(6,'(/A)') 'pd2: entry print_debug1'


!!! debug
write(6,'(/A,1x,A)') 'pd2: print ', label
write(6,'(A)') &
   'pd2: i_parm, input_array(i_parm, i_GP_individual )'
do  i_GP_individual = 1, n_GP_individuals
    do  i_parm = 1, n_Maximum_Number_Parameters
        if( abs( input_array(i_parm, i_GP_individual ) ) > 0.0d0 )then
            write(6,'(I6,1x,I6, 10x, E15.7)') &
                  i_GP_Individual, i_parm, &
                    input_array(i_parm, i_GP_individual )
        endif ! abs( input_array(i_parm, i_GP_individual ) ) > 0.0d0
    enddo
enddo ! i_GP_individual

!-------------------------------------------------------------------------------------------------

write(6,'(A//)') 'pd2: at return   '

return

end subroutine print_debug_real_nparm