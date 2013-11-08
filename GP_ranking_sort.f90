subroutine GP_ranking_sort


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module

implicit none

real(kind=8) :: cff

integer(kind=4),dimension(n_GP_Individuals)  :: Ranked_Fitness_Index

integer(kind=4) :: i_GP_Individual
integer(kind=4) :: j_GP_Individual

integer(kind=4) :: icff
integer(kind=4) :: i_tree
integer(kind=4) :: i_node

integer(kind=4) :: i_parm

real(kind=8), dimension( 1:n_Nodes,1:n_Trees, 1:n_GP_individuals ) :: &
                         GP_population_node_parameters_temp

real(kind=8), dimension( 1:n_CODE_equations, 1:n_GP_individuals ) :: &
                         GP_Population_Initial_Conditions_temp

!!!!!real(kind=8), dimension(1:n_Maximum_Number_Parameters, 1:n_GP_Individuals) :: &
!!!!!                         GP_Adult_Population_Parameter_Solution_temp


!--------------------------------------------------------------------------------


write(6,'(/A)') 'gprs: entry GP_ranking_sort '



! Set up a simple 'index' array

do  i_GP_Individual=1,n_GP_Individuals
    Ranked_Fitness_Index(i_GP_Individual)=i_GP_Individual
enddo

!write(6,'(/A)') 'gprs: before sort '
!write(6,'(A)')                    &
!      'gprs:i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
!                            &GP_Child_Individual_SSE(i_GP_Individual)'
!do  i_GP_Individual=1,n_GP_Individuals
!    write(6,'(I6,1x, I6, 1x, E15.7)') &
!          i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
!                           GP_Child_Individual_SSE(i_GP_Individual)
!enddo


!-------------------------------------------------------------------------------------------------

! Now, rank the Individual SSE so that
! the Individual with the lowest (highest) SSE is First (Last)

do  i_GP_Individual=1,n_GP_Individuals

    do  j_GP_Individual=1,n_GP_Individuals-1

        if( GP_Child_Individual_SSE(j_GP_Individual+1) .lt. &
              GP_Child_Individual_SSE(j_GP_Individual)) then

            !     Swap the two ranked fitness and index array values around
            cff=GP_Child_Individual_SSE(j_GP_Individual)
            GP_Child_Individual_SSE(j_GP_Individual)=GP_Child_Individual_SSE(j_GP_Individual+1)
            GP_Child_Individual_SSE(j_GP_Individual+1)=cff

            icff=Ranked_Fitness_Index(j_GP_Individual)
            Ranked_Fitness_Index(j_GP_Individual)=Ranked_Fitness_Index(j_GP_Individual+1)
            Ranked_Fitness_Index(j_GP_Individual+1)=icff

        endif !GP_Child_Individual_SSE(j_GP_Individual+1) .lt. ...

    enddo ! j_GP_Individual

enddo  ! i_GP_Individual

!-------------------------------------------------------------------------------------------------

write(6,'(/A)') 'gprs: after  sort '
write(6,'(A)')  'gprs:i_GP_Individual, Ranked_Fitness_Index, &
                            &GP_Child_Individual_SSE'
do  i_GP_Individual=1,n_GP_Individuals
    write(6,'(5x,I10,1x, I10, 16x, E15.7)') &
          i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
                           GP_Child_Individual_SSE(i_GP_Individual)
enddo

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------


GP_Adult_Population_SSE=GP_Child_Individual_SSE
GP_Adult_Individual_SSE=GP_Child_Individual_SSE

!GP_Child_Individual_SSE=GP_Child_Population_SSE


!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
! Re-rank ALL of the Individuals to keep the code simple and not replicate copies of children
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

!! debug
!write(6,'(/A)') 'gprs: before applying  sort to GP_Child_Population_Node_Type '
!call print_debug_integer_node_tree( &
!      'from GP_ranking_sort before sort GP_Child_Population_Node_Type', &
!      GP_Child_Population_Node_Type ) 

!-------------------------------------------------------------------------------------------------

GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,  1:n_GP_Individuals) = &
GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,  Ranked_Fitness_Index(1:n_GP_Individuals) )

! Copy this back across to the Child Population values
! to allow the Elite codes to propagate along in the next generations

GP_Child_Population_Node_Type=GP_Adult_Population_Node_Type

!-------------------------------------------------------------------------------------------------

!write(6,'(/A)') 'gprs: after applying  sort to GP_Child_Population_Node_Type '
!call print_debug_integer_node_tree( &
!      'from GP_ranking_sort after  sort GP_Child_Population_Node_Type', &
!      GP_Child_Population_Node_Type ) 
!
!write(6,'(/A)') 'gprs: after applying  sort to GP_Adult_Population_Node_Type '
!call print_debug_integer_node_tree( &
!      'from GP_ranking_sort after  sort GP_Adult_Population_Node_Type', &
!      GP_Adult_Population_Node_Type ) 

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

!!! debug
!!!!!!write(6,'(/A)') 'gprs: before applying  sort to GP_Adult_Population_Parameter_Solution '
!!!!!!call print_debug_real_nparm( 'from GP_ranking_sort before sort GP_Adult_Population_Parameter_Solution', &
!!!!!!                    GP_Adult_Population_Parameter_Solution ) 

!!!!!!-------------------------------------------------------------------------------------------------
!!!!!
!!!!!GP_Adult_Population_Parameter_Solution_temp(1:n_Maximum_Number_Parameters,  1:n_GP_Individuals) = &
!!!!!     GP_Adult_Population_Parameter_Solution(1:n_Maximum_Number_Parameters,  &
!!!!!                                                      Ranked_Fitness_Index(1:n_GP_Individuals) )
!!!!!
!!!!!GP_Adult_Population_Parameter_Solution = GP_Adult_Population_Parameter_Solution_temp
!!!!!
!!!!!!-------------------------------------------------------------------------------------------------

! debug
!!!!!!write(6,'(/A)') 'gprs: after  applying  sort to GP_Adult_Population_Parameter_Solution '
!!!!!
!!!!!!call print_debug_real_nparm( &
!!!!!!      'from GP_ranking_sort AFTER sort GP_Adult_Population_Parameter_Solution', &
!!!!!!      GP_Adult_Population_Parameter_Solution ) 

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

! GP_Population_Initial_Conditions(1:n_CODE_Equations, 1:n_GP_Individuals )


! debug
!write(6,'(/A)') 'gprs: before applying  sort to GP_Population_Initial_Conditions       '
!call print_debug_real_nparm( &
!       'from GP_ranking_sort before sort GP_Population_Initial_Conditions ', &
!       GP_Population_Initial_Conditions ) 

!-------------------------------------------------------------------------------------------------

! sort the GP_population_initial_conditions

do  i_GP_individual = 1, n_GP_individuals

    GP_Population_Initial_Conditions_temp(1:n_CODE_Equations, i_GP_individual ) = &
         GP_Population_Initial_Conditions(1:n_CODE_Equations, &
                                                  Ranked_Fitness_Index(i_GP_individual) )

enddo ! i_GP_individual

GP_Population_Initial_Conditions = GP_Population_Initial_Conditions_temp

!-------------------------------------------------------------------------------------------------


!write(6,'(/A)') 'gprs: after applying  sort to GP_Population_Initial_Conditions       '

!call print_debug_real_nparm( &
!      'from GP_ranking_sort after sort GP_Population_Initial_Conditions ', &
!      GP_Population_Initial_Conditions ) 


!-------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

! debug
!write(6,'(/A)') 'gprs: before applying  sort to GP_population_node_parameters '

!call print_debug_real_node_tree( &
!      'from GP_ranking_sort before  sort GP_population_node_parameters', &
!      GP_population_node_parameters) 

!-------------------------------------------------------------------------------------------------

! sort the GP_population_node_parameters

do  i_GP_individual = 1, n_GP_individuals

    GP_population_node_parameters_temp(1:n_Nodes,1:n_Trees, i_GP_individual ) = &
         GP_population_node_parameters(1:n_Nodes,1:n_Trees, Ranked_Fitness_Index(i_GP_individual) )

enddo ! i_GP_individual

GP_population_node_parameters = GP_population_node_parameters_temp

!-------------------------------------------------------------------------------------------------

! debug
!write(6,'(/A)') 'gprs: after applying  sort to GP_population_node_parameters '

!call print_debug_real_node_tree( &
!  'from GP_ranking_sort after sort GP_population_node_parameters', &
!  GP_population_node_parameters) 


!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!!
!!! original calculation
!!
!!! Calculate the Adult Population's Total SSE
!!
!!cff=0.0d0
!!do  i_GP_Individual=1,n_GP_Individuals
!!    cff=cff+GP_Child_Individual_SSE(i_GP_Individual)
!!enddo
!!
!!write(6,'(/A, 1x, E15.7)') &
!!      'gprs: after: sum GP_Child_Individual_SSE ', cff
!!
!!! Calculate a simple 'normalized' ranking of the SSE as an estimate of fitness
!!
!!! [Fitness = (Total-SSE)/Total
!!!  ==> higher individual SSE == lower value/ranking; Ranging from 0-1]
!!
!!
!!GP_Population_Ranked_Fitness = 0.0D0
!!do  i_GP_Individual=1,n_GP_Individuals
!!
!!    if( cff > 0.0D0 )then
!!        GP_Population_Ranked_Fitness(i_GP_Individual) = &
!!                ( cff - GP_Child_Individual_SSE(i_GP_Individual) ) / cff
!!    else
!!        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0d0
!!    endif ! cff > 0.0d0
!!
!!enddo  ! i_GP_Individual
!!
!!
!------------------------------------------------------------------------------------

! Calculate the Adult Population's Total SSE

cff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals
    if( GP_Child_Individual_SSE(i_GP_Individual) < 1.0e12 )then
        cff=cff+GP_Child_Individual_SSE(i_GP_Individual)
    endif
enddo

write(6,'(/A, 1x, E15.7)') &
      'gprs: after: sum GP_Child_Individual_SSE ', cff

! Calculate a simple 'normalized' ranking of the SSE as an estimate of fitness

! [Fitness = (Total-SSE)/Total
!  ==> higher individual SSE == lower value/ranking; Ranging from 0-1]


GP_Population_Ranked_Fitness = 0.0D0
do  i_GP_Individual=1,n_GP_Individuals

    if( cff > 0.0D0 )then !.and. GP_Child_Individual_SSE(i_GP_Individual) < 1.0e12 )then
        GP_Population_Ranked_Fitness(i_GP_Individual) = &
                ( cff - GP_Child_Individual_SSE(i_GP_Individual) ) / cff
    else
        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0d0
    endif ! cff > 0.0d0

enddo  ! i_GP_Individual


!------------------------------------------------------------------------------------


! Calculate the Integrated Ranked Fitness values for creating the next generation

GP_Integrated_Population_Ranked_Fitness = 0.0D0

cff=0.0d0

do  i_GP_Individual=1,n_GP_Individuals
    cff = cff + GP_Population_Ranked_Fitness(i_GP_individual)
    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = cff
enddo ! i_GP_Individual

! Normalize to the integrated ranking values so that the ranking integration ranges from [0. to 1.]

write(6,'(/A, 1x, E15.7)') &
      'gprs: GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) ', &
             GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

if( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) > 0.0d0 )then

    do  i_GP_Individual=1,n_GP_Individuals

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) =  &
             GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) /  &
                       GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

    enddo  ! i_GP_Individual

endif ! GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) > 0.0d0

!------------------------------------------------------------------------------------------

write(6,'(/A)') 'gprs: after  sort '

write(6,'(A)')                    &
    'gprs:i_GP_Individual GP_Integ_Pop_Ranked_Fitness GP_Pop_Ranked_Fitness GP_Child_Indiv_SSE'

do  i_GP_Individual=1,n_GP_Individuals
    write(6,'(I10, 20x,3(3x, E15.7))') &
          i_GP_Individual, GP_Integrated_Population_Ranked_Fitness(i_GP_Individual), &
                           GP_Population_Ranked_Fitness(i_GP_Individual), &
                           GP_Child_Individual_SSE(i_GP_Individual)
enddo   ! i_GP_Individual



write(6,'(/A)') 'gprs: at return   '

return

end subroutine GP_ranking_sort
