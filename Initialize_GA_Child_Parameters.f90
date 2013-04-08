subroutine Initialize_GA_Child_Parameters(Child_Parameters)


use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none


real(kind=8) :: Child_Parameters(n_GA_Individuals,n_maximum_number_parameters)
real(kind=4) :: cff
real(kind=8) :: dff

integer(kind=4) :: i

!----------------------------------------------------------------------------

Run_GA_lmdif=.true.

write(6,'(/A,1x, I6/)')  'Init: n_parameters ', n_Parameters

!write(6,'(/A/)') &
!      'Init:  i_GA_individual  child parameters(i_GA_individual,:)  '


do  i_GA_Individual=1,n_GA_Individuals

    do  i_Parameter=1,n_Parameters


        call random_real(cff) ! random real number generator
        dff = cff
        Child_Parameters(i_GA_Individual,i_Parameter) = dff


!        !  parameter 11 == parameter 7
!        !  parameter 12 == parameter 8
!
!        if( i_parameter == 11 )then
!            Child_Parameters(i_GA_Individual,i_Parameter) =  &
!                            Child_Parameters(i_GA_Individual,7)
!        endif
!        if( i_parameter == 12 )then
!            Child_Parameters(i_GA_Individual,i_Parameter) =  &
!                            Child_Parameters(i_GA_Individual,8)
!        endif


         !----------------------------------------------------------------------------

         if( n_linked_parms > 0 )then

             if( any( i_parameter == linked_parms(2,:) )    ) cycle

             if( any( i_parameter == linked_parms(1,:) )    )then

                 do  i = 1, n_linked_parms

                     if( i_parameter == linked_parms(1,i) )then

                         child_parameters(i_GA_Individual, linked_parms(2,i) ) = &
                         child_parameters(i_GA_Individual, linked_parms(1,i) )

                         exit

                     endif ! i_parameter == linked_parms(i,1)

                 enddo ! i

             endif !  any( i_parameter == linked_parms(1,:) )

         endif !  n_linked_parms > 0

         !----------------------------------------------------------------------------



        !write(6,'(A,2(1x, I6),1x,E24.16 )') &
        !     'Init: ', i_GA_individual, i_Parameter, &
        !               child_parameters(i_GA_individual,i_parameter)

    enddo ! i_parameter

    !write(6,'(A,1x, I6,1x,12(1x,E15.7))') &
    !     'Init: ', i_GA_individual, child_parameters(i_GA_individual,1:n_parameters)

enddo ! i_ga_individual


!!--------------------------------------------------------------------------------
!!do  i_GA_Individual=1,n_GA_Individuals                           ! debug only
!!    read(5,*) Child_Parameters(i_GA_Individual,1:n_parameters)   ! debug only
!!enddo ! i_ga_individual                                          ! debug only
!!--------------------------------------------------------------------------------


return


end subroutine Initialize_GA_Child_Parameters
!234567890123456789012345678901234567890123456789012345678901234567890
