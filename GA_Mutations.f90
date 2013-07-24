subroutine GA_Mutations(Child_Parameters, individual_quality )


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=8) :: child_parameters(n_maximum_number_parameters,n_GA_Individuals)
real(kind=4) :: cff
real(kind=8) :: dff

integer (kind=4) :: i_GA_Mutation
integer (kind=4) :: i_GA_Individual_Mutation, i_Parameter_Mutation

integer(kind=4) :: individual_quality(n_GA_individuals)

integer (kind=4) :: n_mutated

!---------------------------------------------------------------------

if( n_GA_Mutations < 1 ) return


!write(GA_print_unit,'(//A,1x,I6/)') 'gam: n_GA_Mutations ', n_GA_Mutations


n_mutated  = 0

do i_GA_Mutation=1,n_GA_Mutations


  !---------------------------------------------------------------------


  ! randomly pick an individual to mutate [presently a child]

  ! if the index i_GA_Mutation is in the array individual_elites,
  ! do not replace this individual - it is an elite individual

  ! check_for_elite generates random numbers for the individual number
  ! until it finds one not in the list of elite individuals


  call check_for_elite( i_GA_Individual_mutation )


  !--------------------------------------------------------------------

  !write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6)') &
  !      'gam: i_GA_Mutation, dff, i_GA_Individual_mutation ', &
  !            i_GA_Mutation, dff, i_GA_Individual_mutation
  !write(GA_print_unit,'(/A/I6,12(1x,E15.7))') &
  !      'gam: before i_GA_Individual_mutation,  &
  !  &child_parameters(1:n_parameters, i_GA_Individual_mutation ) ', &
  !                    i_GA_Individual_mutation,  &
  !   child_parameters(1:n_parameters, i_GA_Individual_mutation)

  !--------------------------------------------------------------------

  !  randomly pick which parameter will be replaced

  call random_number(cff)   ! uniform random number generator
  dff = cff

  i_Parameter_Mutation=1+int( dff*dble(n_parameters-1) )
  i_Parameter_Mutation = min( i_Parameter_Mutation , n_parameters )

  !write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6)') &
  !      'gam: i_GA_Mutation, dff, i_Parameter_Mutation     ', &
  !            i_GA_Mutation, dff, i_Parameter_Mutation

  !--------------------------------------------------------------------

  !  randomly pick a new real number for this parameter

  call random_real(cff)
  dff = cff

  child_parameters(i_Parameter_Mutation, i_GA_Individual_Mutation) = dff

  !----------------------------------------------------------------------------

  !write(GA_print_unit,'(A/I6,12(1x,E15.7))') &
  !      'gam: after ', &
  !      i_GA_Individual_mutation,  &
  !      child_parameters(1:n_parameters, i_GA_Individual_mutation )

  !write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6/)') &
  !      'gam: i_GA_Individual_Mutation, child_parameters(i_Parm_Mut, i_GA_Ind_Mut) ', &
  !            i_GA_Individual_Mutation, &
  !       child_parameters(i_Parameter_Mutation, i_GA_Individual_Mutation)

  !--------------------------------------------------------------------

  ! set the flag to do the RK integration on this parameter

  Run_GA_lmdif(i_GA_Individual_Mutation)=.true.


  ! I don't think this is needed,
  ! since the individual_quality will be set to 1 later

  individual_quality(i_GA_Individual_Mutation) = 1   


  n_mutated  = n_mutated  + 1

enddo

write(GA_print_unit,'(A,1x,I6,1x,I10/)') &
      'gam: i_GA_generation, n_mutated ',  &
            i_GA_generation, n_mutated

return
end subroutine GA_Mutations
