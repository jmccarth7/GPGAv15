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

integer(kind=4) :: i
integer(kind=4) :: i_parameter
integer(kind=4) :: jj

!----------------------------------------------------------------------------

Run_GA_lmdif=.true.

write(GA_print_unit,'(/A,1x, I6/)')  'Init: n_parameters ', n_Parameters

!write(GA_print_unit,'(/A/)') &
!      'Init:  i_GA_individual  child parameters(:,i_GA_individual)  '


do  i_GA_Individual=1,n_GA_Individuals

    do  i_Parameter=1,n_Parameters


        call random_real(cff) ! random real number generator
        dff = cff
        Child_Parameters(i_Parameter,i_GA_Individual) = dff


         !----------------------------------------------------------------------------

!!         if( n_linked_parms > 0 )then
!!
!!             if( any( i_parameter == linked_parms(2,:) )    ) cycle
!!
!!             if( any( i_parameter == linked_parms(1,:) )    )then
!!
!!                 do  i = 1, n_linked_parms
!!
!!                     if( i_parameter == linked_parms(1,i) )then
!!
!!                         child_parameters(linked_parms(2,i),i_GA_Individual) = &
!!                         child_parameters(linked_parms(1,i),i_GA_Individual)
!!
!!                         exit
!!
!!                     endif ! i_parameter == linked_parms(i,1)
!!
!!                 enddo ! i
!!
!!             endif !  any( i_parameter == linked_parms(1,:) )
!!
!!         endif !  n_linked_parms > 0

         !----------------------------------------------------------------------------


        !write(GA_print_unit,'(A,2(1x, I6),1x,E24.16 )') &
        !     'Init: ', i_GA_individual, i_Parameter, &
        !               child_parameters(i_parameter, i_GA_individual)

    enddo ! i_parameter

    write(GA_print_unit,'(A,1x, I6,1x,12(1x,E15.7))') &
         'Init: ', i_GA_individual, &
          ( child_parameters(jj, i_GA_individual), jj = 1,n_parameters )
         !'Init: ', i_GA_individual, child_parameters(1:n_parameters, i_GA_individual)

enddo ! i_ga_individual


!!--------------------------------------------------------------------------------
!!do  i_GA_Individual=1,n_GA_Individuals                           ! debug only
!!    read(5,*) Child_Parameters(1:n_parameters,i_GA_Individual)   ! debug only
!!enddo ! i_ga_individual                                          ! debug only
!!--------------------------------------------------------------------------------


return


end subroutine Initialize_GA_Child_Parameters
