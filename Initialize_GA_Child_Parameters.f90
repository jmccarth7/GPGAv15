subroutine Initialize_GA_Child_Parameters(Child_Parameters)


use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none


real(kind=8) :: Child_Parameters(n_maximum_number_parameters,n_GA_Individuals)
real(kind=4) :: cff
real(kind=8) :: dff

integer(kind=4) :: i_parameter
!integer(kind=4) :: jj

!----------------------------------------------------------------------------

Run_GA_lmdif=.true.

!if( L_ga_print )then                                                                                
!    write(GA_print_unit,'(/A,1x, I6/)')  'Init: n_parameters ', n_Parameters
!    write(GA_print_unit,'(/A/)') &
!          'Init:  i_GA_individual  child parameters(:,i_GA_individual)  '
!endif ! L_ga_print  


do  i_GA_Individual=1,n_GA_Individuals

    do  i_Parameter=1,n_Parameters


        call random_real(cff) ! random real number generator
        dff = cff
        Child_Parameters(i_Parameter,i_GA_Individual) = dff


        !if( L_ga_print )then 
        !    write(GA_print_unit,'(A,2(1x, I6),1x,E24.16 )') &
        !         'Init: ', i_GA_individual, i_Parameter, &
        !                   child_parameters(i_parameter, i_GA_individual)
        !endif ! L_ga_print 

    enddo ! i_parameter

    !if( L_ga_print )then  
    !    write(GA_print_unit,'(A,1x, I6,1x,12(1x,E15.7))') &
    !         'Init: ', i_GA_individual, &
    !          ( child_parameters(jj, i_GA_individual), jj = 1,n_parameters )
    !         !'Init: ', i_GA_individual, child_parameters(1:n_parameters, i_GA_individual)
    !endif ! L_ga_print  

enddo ! i_ga_individual


!!--------------------------------------------------------------------------------
!!do  i_GA_Individual=1,n_GA_Individuals                           ! debug_only
!!    read(5,*) Child_Parameters(1:n_parameters,i_GA_Individual)   ! debug_only
!!enddo ! i_ga_individual                                          ! debug_only
!!--------------------------------------------------------------------------------


return


end subroutine Initialize_GA_Child_Parameters
