subroutine GA_Mutations(Child_Parameters, individual_quality )


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=8) :: child_parameters(n_GA_Individuals,n_maximum_number_parameters)
real(kind=4) :: cff
real(kind=8) :: dff

integer (kind=4) :: i
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

  call check_for_elite( i_GA_Individual_mutation )

  !--------------------------------------------------------------------

  !write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6)') &
  !      'gam: i_GA_Mutation, dff, i_GA_Individual_mutation ', &
  !            i_GA_Mutation, dff, i_GA_Individual_mutation
  !write(GA_print_unit,'(/A/I6,12(1x,E15.7))') &
  !      'gam: before i_GA_Individual_mutation,  &
  !  &child_parameters(i_GA_Individual_mutation, 1:n_parameters) ', &
  !                    i_GA_Individual_mutation,  &
  !   child_parameters(i_GA_Individual_mutation, 1:n_parameters)


  !--------------------------------------------------------------------

  !  randomly pick which parameter will be replaced 

  call random_number(cff)   ! uniform random number generator
  dff = cff

  i_Parameter_Mutation=1+int( dff*dble(n_parameters-1) ) 

  !write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6)') &
  !      'gam: i_GA_Mutation, dff, i_Parameter_Mutation     ', &
  !            i_GA_Mutation, dff, i_Parameter_Mutation

  !--------------------------------------------------------------------

  !  randomly pick a new real number for this parameter 

  call random_real(cff)   
  dff = cff

  child_parameters(i_GA_Individual_Mutation,i_Parameter_Mutation) = dff 

!  if( i_parameter_mutation == 7  )then
!      child_parameters(i_GA_Individual_Mutation,11) = &
!      child_parameters(i_GA_Individual_Mutation,i_Parameter_Mutation) 
!  endif !  i_parameter_mutation == 7
!
!  if( i_parameter_mutation == 8  )then
!      child_parameters(i_GA_Individual_Mutation,12) = &
!      child_parameters(i_GA_Individual_Mutation,i_Parameter_Mutation) 
!  endif !  i_parameter_mutation == 8
!
!  if( i_parameter_mutation == 11 )then
!      child_parameters(i_GA_Individual_Mutation,7)  = &
!      child_parameters(i_GA_Individual_Mutation,i_Parameter_Mutation) 
!  endif !  i_parameter_mutation == 11 
!
!  if( i_parameter_mutation == 12 )then
!      child_parameters(i_GA_Individual_Mutation,8)  = &
!      child_parameters(i_GA_Individual_Mutation,i_Parameter_Mutation) 
!  endif !  i_parameter_mutation == 12 

  !----------------------------------------------------------------------------

  !write(GA_print_unit,'(A/I6,12(1x,E15.7))') &
  !      'gam: after ', &
  !      i_GA_Individual_mutation,  child_parameters(i_GA_Individual_mutation, 1:n_parameters)

  !write(GA_print_unit,'(A,1x,I6,1x,E15.7,1x,I6/)') &
  !      'gam: i_GA_Individual_Mutation, child_parameters(i_GA_Ind_Mut,i_Parm_Mut) ', &
  !            i_GA_Individual_Mutation, child_parameters(i_GA_Individual_Mutation,i_Parameter_Mutation)

  !--------------------------------------------------------------------

  ! set the flag to do the RK integration on this parameter

  Run_GA_lmdif(i_GA_Individual_Mutation)=.true.  


  ! I don't think this is needed, since the individual_quality will be set to 1 later

  individual_quality(i_GA_Individual_Mutation) = 1   ! reset quality since this has been mutated


  n_mutated  = n_mutated  + 1


  !----------------------------------------------------------------------------

  if( n_linked_parms > 0 )then

      if( any( i_parameter_mutation == linked_parms(2,:) )    ) cycle

      if( any( i_parameter_mutation == linked_parms(1,:) )    )then

          do  i = 1, n_linked_parms 
    
              if( i_parameter_mutation == linked_parms(1,i) )then
    
                  child_parameters(i_GA_Individual_Mutation,linked_parms(2,i) ) = &
                  child_parameters(i_GA_Individual_Mutation,i_Parameter_Mutation) 
    
                  exit
    
              endif ! i_parameter_mutation == linked_parms(i,1)
    
          enddo ! i 

      endif !  any( i_parameter_mutation == linked_parms(1,:) )  

  endif !  n_linked_parms > 0 

  !----------------------------------------------------------------------------

enddo

write(GA_print_unit,'(A,1x,I6,1x,I10/)') &
      'gam: i_GA_generation, n_mutated ',  &
            i_GA_generation, n_mutated

return
end subroutine GA_Mutations
!234567890123456789012345678901234567890123456789012345678901234567890
