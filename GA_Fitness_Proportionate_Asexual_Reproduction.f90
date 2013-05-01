subroutine GA_Fitness_Proportionate_Reproduction( &
                            Parent_Parameters,Child_Parameters, &
                            individual_quality )
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=8)    :: parent_parameters(n_GA_Individuals,n_maximum_number_parameters)
real(kind=8)    :: child_parameters(n_GA_Individuals,n_maximum_number_parameters)
integer(kind=4) :: individual_quality(n_GA_individuals)

real(kind=4) :: cff
real(kind=8) :: dff
real(kind=8) :: mean_fit_before
real(kind=8) :: mean_fit_after
integer (kind=4) :: icff

integer (kind=4) :: icount
integer (kind=4) :: n_replaced

integer(kind=4) :: i_GP_individual                                                                           
integer(kind=4) :: i_GP_Best_Parent                                                                          
integer(kind=4) :: i_GP_Generation                                                                           
integer(kind=4) :: j_GP_Individual                                                                           
integer(kind=4) :: i_Tree                                                                                    
integer(kind=4) :: i_Node            

!--------------------------------------------------------------------------


!write(GA_print_unit,'(/A/)') 'gafp: at entry'


! for each individual,i,  choose a random number in  [0.0, 1.0]

! the range of the integrated_ranked_fitness is also [0.0, 1.0]

! cycle through all individuals until one, j,  is found such that:
!     the integrated_ranked_fitness(j) > random number

! then replace child parameters of i with child parameters of j




mean_fit_before = 0.0d0
icount = 0
do  i_GA_individual = 1, n_GA_individuals

    !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
    !      'gafp: before  i_GA_individual , Indiv_Fit ', &
    !                     i_GA_individual , Individual_Ranked_Fitness(i_GA_Individual)

    if( individual_quality(i_GA_individual) > 0  .and.  &
        Individual_Ranked_Fitness(i_GA_Individual) > 1.0d0 )then

        mean_fit_before =  mean_fit_before + Individual_Ranked_Fitness(i_GA_Individual)
        icount = icount + 1

    endif ! individual_quality...

enddo ! i_GA_individual


if( icount > 0 )then
    mean_fit_before =  mean_fit_before / real( icount, kind=8 )
else
    mean_fit_before = 0.0d0
endif ! icount > 0

!-------------------------------------------------------------------------------

n_replaced = 0
icff = 0

i_loop:&
do i_GA_Individual=1,n_GA_Individuals

  Run_GA_lmdif(i_GA_Individual)=.false.

  call Random_Number(cff) ! uniform random number generator

  dff = cff

  !write(GA_print_unit,'(/A,1x,I6,1x,E15.7)') 'gafp: i_GA_Individual, dff', &
  !                                                  i_GA_Individual, dff

  !--------------------------------------------------------------------------

  ! if the index i_GA_individual is in the array individual_elites,
  ! do not replace this individual - it is an elite individual


  if( any( individual_elites == i_GA_individual ) )then
      !write(GA_print_unit,'(A,1x,I6)') &
      !'gafp: skip -- index is elite i_GA_individual = ', &
      !                              i_GA_individual
      cycle
  endif   ! any( individual_elites == i_GA_individual )

  !--------------------------------------------------------------------------


  ! find an individual to replace the current i_GA_individual


  j_loop:&
  do  j_GA_Individual=1,n_GA_Individuals ! normalize to the maximum values
                                         ! so that the range is [0. , 1.]

      !write(GA_print_unit,'(A,1x,I6,2(1x,E15.7))') &
      !      'gafp: j_GA_Individual, Integrated_Ranked_Fitness(j_GA_Individual), dff ', &
      !             j_GA_Individual, Integrated_Ranked_Fitness(j_GA_Individual), dff

      !----------------------------------------------------------------------------------

      ! don't replace with this individual since it is bad

      if( individual_quality( j_GA_Individual ) < 0 )then
          cycle
      endif

      !----------------------------------------------------------------------------------

      ! set icff to -1 so that, if no j_GA_individual satisfies the test on cff,
      ! then this i_GA_individual will be skipped, and run_ga_lmdif will be false for it.

      icff = -1

      !----------------------------------------------------------------------------------


      if( dff .le. Integrated_Ranked_Fitness(j_GA_Individual) ) then

          icff=j_GA_Individual
          exit j_loop

      endif  ! dff .le. ...


  enddo j_loop ! j_GA_Individual



  if( icff > 0 )then
      j_GA_Individual=icff ! index to move over both 1) the parent parameters and
                           !                         2) the individual fitness levels
  else
      cycle i_loop   ! skip replacing this individual
  endif

  !write(GA_print_unit,'(/A,2(1x,I6))') &
  !'gafp:aft loop  j_GA_Individual, icff', j_GA_Individual, icff

  !-----------------------------------------------------------------------------------------

  ! do the replacements here

  Individual_Ranked_Fitness(i_GA_Individual) = Individual_Ranked_Fitness(j_GA_Individual)
  individual_quality( i_GA_individual )      = individual_quality( j_GA_individual )

  n_replaced = n_replaced + 1


  do i_Parameter=1,n_Parameters

    Child_Parameters(i_GA_Individual,i_Parameter)=Parent_Parameters(j_GA_Individual,i_Parameter)

  enddo ! i_Parameter


  !-----------------------------------------------------------------------------------------


  !write(GA_print_unit,'(A,2(1x,I6))') 'gafp: i_GA_Individual replaced by j_GA_individual', &
  !                               i_GA_Individual, j_GA_Individual

  !write(GA_print_unit,'(A,1x,I6, 1x, E15.7/)') &
  !      'gafp:aft loop i_GA_Individual, Individual_Ranked_Fitness', &
  !                     i_GA_Individual, Individual_Ranked_Fitness(i_GA_Individual)

  !write(GA_print_unit,'(I6,12(1x,E15.7))') &
  !      i_GA_Individual,  Child_Parameters(i_GA_Individual,1:n_Parameters)

enddo i_loop  ! i_GA_Individual



!write(GA_print_unit,'(/A)') &
!'gafp: i_GA_Individual,  Child_Parameters(i_GA_Individual,1:n_Parameters) '

!do  i_GA_Individual=1,n_GA_Individuals
!    write(GA_print_unit,'(I6,12(1x,E15.7))') &
!         i_GA_Individual,  Child_Parameters(i_GA_Individual,1:n_Parameters)
!enddo  ! i_GA_Individual



!write(GA_print_unit,'(A)') ' '


mean_fit_after = 0.0d0
icount = 0
do  i_GA_individual = 1, n_GA_individuals

    !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
    !      'gafp: after i_GA_individual, Indiv_Fit ', &
    !                   i_GA_individual, Individual_Ranked_Fitness(i_GA_Individual)

    if( individual_quality(i_GA_individual) > 0  .and.  &
        Individual_Ranked_Fitness(i_GA_Individual) >  1.0d0 )then

        mean_fit_after =  mean_fit_after + Individual_Ranked_Fitness(i_GA_Individual)
        icount = icount + 1

    endif ! individual_quality...

enddo ! i_GA_individual

if( icount > 0 )then
    mean_fit_after =  mean_fit_after / real( icount , kind=8 )
else
    mean_fit_after =  0.0D0
endif ! icount > 0


write(GA_print_unit,'(/A,1x,I6,1x, E20.10)') &
      'gafp: i_ga_generation, mean_fit_before      ', &
             i_ga_generation, mean_fit_before
write(GA_print_unit,'(A,1x,I6,1x, E20.10)')  &
      'gafp: i_ga_generation, mean_fit_after       ', &
             i_ga_generation, mean_fit_after
write(GA_print_unit,'(A,1x,I6,2(1x,I10)/)')  &
      'gafp: i_ga_generation, n_replaced, icount   ', &
             i_ga_generation, n_replaced, icount

return


end subroutine GA_Fitness_Proportionate_Reproduction
!234567890123456789012345678901234567890123456789012345678901234567890
