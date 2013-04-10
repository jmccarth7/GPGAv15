subroutine GA_Tournament_Style_Sexual_Reproduction(&
              Parent_Parameters,Child_Parameters, individual_quality  )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!  randomly choose two 'parents' using the Tournament-Style Selection and
!  cross the parameter strings to create two new 'children' parameter strings
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=8) :: parent_parameters(n_GA_Individuals,n_maximum_number_parameters)
real(kind=8) ::  child_parameters(n_GA_Individuals,n_maximum_number_parameters)

integer(kind=4) :: individual_quality(n_GA_individuals)

integer (kind=4) :: i_GA_Crossover,i_GA_Crossover_Point
integer (kind=4) :: k_GA_Individual_Male(2),k_GA_Individual_Female(2)

real(kind=8) :: child_one_parameters(n_parameters)
real(kind=8) :: child_two_parameters(n_parameters)

real(kind=8) :: temp_male_parameters(n_parameters)
real(kind=8) :: temp_female_parameters(n_parameters)

real(kind=4) :: cff
real(kind=8) :: dff

integer (kind=4) :: n_replaced
integer (kind=4) :: i            

!---------------------------------------------------------------------------

!write(GA_print_unit,'(/A,1x,I6)')'gato: n_GA_Crossovers', n_GA_Crossovers

n_replaced = 0

do i_GA_Crossover=1,n_GA_Crossovers


  !--------------------------------------------------------------------

  ! pick male parent 1 for sexual crossing of parameter strings


  ! check_for_elite generates random numbers and computes an individual
  ! number until it finds one which is not in the elite set of individuals
  ! this non-elite number is returned to the calling program

  call check_for_elite( k_GA_Individual_Male(1) )


  !write(GA_print_unit,'(/A,1x,I6)')      'gato: k_GA_Individual_Male(1)  ', k_GA_Individual_Male(1)
  !write(GA_print_unit,'(A/(15(1x,I6)))') 'gato: individual_elites ', individual_elites


  !--------------------------------------------------------------------

  ! pick male parent 2 for sexual crossing of parameter strings


  call check_for_elite( k_GA_Individual_Male(2) )

  !write(GA_print_unit,'(//A,3(1x,I6))')&
  !      'gato: aft random i_GA_Crossover, k_GA_Individual_Male(1:2)  ', &
  !                        i_GA_Crossover, k_GA_Individual_Male(1:2)

  !--------------------------------------------------------------------

  ! you picked the same individual for both male parents, so choose another

  if( k_GA_Individual_Male(2) .eq. k_GA_Individual_Male(1)) then

      if( k_GA_Individual_Male(1) .ne. N_GA_Individuals) then
          k_GA_Individual_Male(2)= min( k_GA_Individual_Male(1)+1, N_GA_Individuals )
      else
          k_GA_Individual_Male(2)= max( k_GA_Individual_Male(1)-1, 1 )
      endif !   k_GA_Individual_Male(1) .ne. N_GA_Individuals

  endif ! k_GA_Individual_Male(2) .eq. k_GA_Individual_Male(1)


  !--------------------------------------------------------------------


  ! select the individual of the two with the best fitness
  ! best fitness means Individual_Ranked_Fitness is largest


  !write(GA_print_unit,'(A,1x,E24.16)')&
  !      'gato: Individual_Ranked_Fitness(k_GA_Individual_Male(1)) ', &
  !             Individual_Ranked_Fitness(k_GA_Individual_Male(1))
  !write(GA_print_unit,'(A,1x,E24.16)')&
  !      'gato: Individual_Ranked_Fitness(k_GA_Individual_Male(2)) ', &
  !             Individual_Ranked_Fitness(k_GA_Individual_Male(2))

  if( Individual_Ranked_Fitness(k_GA_Individual_Male(1)) .lt. &
      Individual_Ranked_Fitness(k_GA_Individual_Male(2))        ) then

      k_GA_Individual_Male(1)=k_GA_Individual_Male(2)

  endif !   Individual_Ranked_Fitness(k_GA_Individual_Male(1)) .lt. ...

  !write(GA_print_unit,'(A,3(1x,I6))')&
  !      'gato: selected male i_GA_Crossover, k_GA_Individual_Male(1) ', &
  !                           i_GA_Crossover, k_GA_Individual_Male(1)

  !---------------------------------------------------------------------------------

  ! pick female parent 1 for sexual crossing of parent parameter strings

  call check_for_elite( k_GA_Individual_Female(1) )

  !---------------------------------------------------------------------------------

  ! pick female parent 2 for sexual crossing of parent parameter strings

  call check_for_elite( k_GA_Individual_Female(2) )

  !write(GA_print_unit,'(/A,3(1x,I6))')&
  !      'gato: aft random i_GA_Crossover, k_GA_Individual_Female(1:2)', &
  !                        i_GA_Crossover, k_GA_Individual_Female(1:2)

  !---------------------------------------------------------------------------------

  ! you picked the same individual for both female parents, so choose another


  if( k_GA_Individual_Female(2) .eq. k_GA_Individual_Female(1)) then

      if( k_GA_Individual_Female(1) .ne. N_GA_Individuals) then
          k_GA_Individual_Female(2) =  min( k_GA_Individual_Female(1)+1, N_GA_Individuals )
      else
          k_GA_Individual_Female(2) =  max( k_GA_Individual_Female(1)-1, 1 )
      endif !   k_GA_Individual_Female(1) .ne. N_GA_Individuals)

  endif !   k_GA_Individual_Female(2) .eq. k_GA_Individual_Female(1)


  !---------------------------------------------------------------------------------


  ! select the individual of the two with the best fitness
  ! best fitness means Individual_Ranked_Fitness is largest

  !write(GA_print_unit,'(A,1x,E24.16)')&
  !      'gato: Individual_Ranked_Fitness(k_GA_Individual_Female(1)) ', &
  !             Individual_Ranked_Fitness(k_GA_Individual_Female(1))
  !write(GA_print_unit,'(A,1x,E24.16)')&
  !      'gato: Individual_Ranked_Fitness(k_GA_Individual_Female(2)) ', &
  !             Individual_Ranked_Fitness(k_GA_Individual_Female(2))


  if( Individual_Ranked_Fitness(k_GA_Individual_Female(1)) .lt. &
      Individual_Ranked_Fitness(k_GA_Individual_Female(2))         ) then

      k_GA_Individual_Female(1)=k_GA_Individual_Female(2)

  endif ! Individual_Ranked_Fitness(k_GA_Individual_Female(1)) .lt. ...

  !---------------------------------------------------------------------------------

  !  save parameters for selected male and female parents before crossover
  !  (just for comparison )

  temp_male_parameters(1:n_parameters)   = Child_Parameters(k_GA_Individual_Male(1),  1:n_parameters)
  temp_female_parameters(1:n_parameters) = Child_Parameters(k_GA_Individual_Female(1),1:n_parameters)

  !---------------------------------------------------------------------------------

  ! choose the location along the parameter string for the crossover to occur

  call Random_Number(cff) ! uniform random number generator
  dff = cff

  i_GA_Crossover_Point=1+int( dff*dble(n_Parameters-2) )  ! pick a location from 1 to n_parameters-1

  !write(GA_print_unit,'(/A,2(1x,I6))')&
  !      'gato: i_GA_Crossover, i_GA_Crossover_Point ', &
  !             i_GA_Crossover, i_GA_Crossover_Point

  !--------------------------------------------------------------------------------

  ! do the crossover at the selected location

  do  i_Parameter=1,n_Parameters

      if( i_parameter .le. i_GA_Crossover_Point) then

          Child_One_Parameters(i_Parameter)=Parent_Parameters(k_GA_Individual_Male(1),i_parameter)
          Child_Two_Parameters(i_Parameter)=Parent_Parameters(k_GA_Individual_Female(1),i_parameter)

      else

          Child_One_Parameters(i_Parameter)=Parent_Parameters(k_GA_Individual_Female(1),i_parameter)
          Child_Two_Parameters(i_Parameter)=Parent_Parameters(k_GA_Individual_Male(1),i_parameter)

      endif !   i_parameter .le. i_GA_Crossover_Point

  enddo ! i_parameter


  !--------------------------------------------------------------------------------

  ! replace the crossover point parameter value
  ! with a new random number in each child

  call random_real(cff)
  dff = cff
  Child_One_Parameters(i_GA_Crossover_Point) = dff


  !----------------------------------------------------------------------------

  if( n_linked_parms > 0 )then

      if( any( i_GA_Crossover_Point == linked_parms(:,:) )    )then

              do  i = 1, n_linked_parms

                  if( i_GA_Crossover_Point == linked_parms(1,i) )then

                      child_one_parameters(linked_parms(2,i) ) = &
                      child_one_parameters(linked_parms(1,i) )

                  endif ! i_GA_Crossover_Point == linked_parms(i,1)

                  if( i_GA_Crossover_Point == linked_parms(2,i) )then

                      child_one_parameters(linked_parms(1,i) ) = &
                      child_one_parameters(linked_parms(2,i) )

                  endif ! i_GA_Crossover_Point == linked_parms(i,1)


              enddo ! i

      endif !  any( i_GA_Crossover_Point == linked_parms(1,:) )

  endif !  n_linked_parms > 0

  !----------------------------------------------------------------------------


  call random_real(cff)
  dff = cff
  Child_Two_Parameters(i_GA_Crossover_Point) = dff


  !----------------------------------------------------------------------------

  if( n_linked_parms > 0 )then

      if( any( i_GA_Crossover_Point == linked_parms(:,:) )    )then

              do  i = 1, n_linked_parms

                  if( i_GA_Crossover_Point == linked_parms(1,i) )then

                      child_two_parameters(linked_parms(2,i) ) = &
                      child_two_parameters(linked_parms(1,i) )

                  endif ! i_GA_Crossover_Point == linked_parms(i,1)

                  if( i_GA_Crossover_Point == linked_parms(2,i) )then

                      child_two_parameters(linked_parms(1,i) ) = &
                      child_two_parameters(linked_parms(2,i) )

                  endif ! i_GA_Crossover_Point == linked_parms(i,1)


              enddo ! i

      endif !  any( i_GA_Crossover_Point == linked_parms(1,:) )

  endif !  n_linked_parms > 0



  !--------------------------------------------------------------------------------

  ! replace the mating pool with the newly crossed parameter strings

  do  i_parameter=1,n_parameters

      Child_Parameters(k_GA_Individual_Male(1),  i_parameter) = Child_One_Parameters(i_Parameter)
      Child_Parameters(k_GA_Individual_Female(1),i_parameter) = Child_Two_Parameters(i_Parameter)

  enddo ! i_parameter

  !------------------------------------------------------------------------------------------

  !  print parameters for selected male and female parents before and after  crossover

  !write(GA_print_unit,'(A,3(1x,I6))')&
  !      'gato: selected i_GA_Crossover, k_GA_Individual_Male(1), k_GA_Individual_Female(1) ', &
  !                      i_GA_Crossover, k_GA_Individual_Male(1), k_GA_Individual_Female(1)

  !write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
  ! 'gato: before k_GA_Individual_Male(1), &
  !  &Child_Parameters(k_GA_Individual_Male(1), 1:n_parameters ) ', &
  !  k_GA_Individual_Male(1),   temp_male_parameters(1:n_parameters)
  !write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
  ! 'gato: after ', &
  !  k_GA_Individual_Male(1), Child_Parameters(k_GA_Individual_Male(1), 1:n_parameters)

  !write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
  ! 'gato: before k_GA_Individual_Female(1), &
  !  &Child_Parameters(k_GA_Individual_Female(1), 1:n_parameters)', &
  !  k_GA_Individual_Female(1), temp_female_parameters(1:n_parameters)
  !write(GA_print_unit,'(A/I6,12(1x,E15.7))')&
  ! 'gato: after  ', &
  !  k_GA_Individual_Female(1), &
  !  &Child_Parameters(k_GA_Individual_Female(1), 1:n_parameters )


  Run_GA_lmdif( k_GA_Individual_Male(1) )   = .true.
  Run_GA_lmdif( k_GA_Individual_Female(1) ) = .true.

  individual_quality( k_GA_Individual_Male(1) )   = 1
  individual_quality( k_GA_Individual_Female(1) ) = 1

  n_replaced = n_replaced + 2

enddo

write(GA_print_unit,'(A,1x,I6,1x,I10)')&
      'gato: i_ga_generation, n_replaced in tournament ', &
             i_ga_generation, n_replaced


return
end subroutine GA_Tournament_Style_Sexual_Reproduction
!234567890123456789012345678901234567890123456789012345678901234567890
