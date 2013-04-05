subroutine GP_fitness_reset( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module


implicit none

!real (kind=4) :: GP_Node_Parameters_Answer(n_Nodes,n_Trees)
!integer (kind=4) :: GP_Node_Type_Answer(n_Nodes,n_Trees)
!integer (kind=4) :: GP_Node_Type_for_Plotting(9,n_Nodes,n_Trees)

!integer(kind=4) :: i_GP_Best_Parent
!integer(kind=4) :: i_GP_Generation
!integer(kind=4) :: i_GP_Individual
!integer(kind=4) :: j_GP_Individual
!integer(kind=4) :: i_Tree
!integer(kind=4) :: i_Node

logical bad_node


integer :: i
integer :: nop
integer :: message_len

logical :: Lplot


real(kind=8), allocatable, dimension(:) :: answer
real(kind=8), allocatable, dimension(:) :: output_array

real (kind=8) :: ssum, ssum2, totobs, dff

real (kind=8) :: xcount

!----------------------------------------------------------------------------------------

allocate( output_array( n_maximum_number_parameters ) )


!write(*,*) GP_Child_Individual_SSE
!if( myid == 0 )then
!    write(6,'(A,1x,I6)') 'gpfr: i_GP_generation ',  i_GP_generation
!
!    write(6,'(/A/(5(1x,E15.7)))') 'gpfr: GP_Child_Individual_SSE = ',&
!                                      GP_Child_Individual_SSE
!
!endif ! myid == 0
!---------------------------------------------------------------------------------------------------

! fitness reset region

! calculate the total population's SSE

dff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals
    dff=dff+GP_Child_Individual_SSE(i_GP_Individual)
enddo

write(6,'(A,1x,I6)') 'gpfr: i_GP_generation ',  i_GP_generation
write(6,'(/A,1x,E15.7)') 'gpfr: sum of all GP_Child_Individual_SSE = ', dff


! calculate a normalized ranking of the errors (higher individual SSE == lower value/ranking)
do  i_GP_Individual=1,n_GP_Individuals
    GP_Individual_Ranked_Fitness(i_GP_Individual) = &
         ( dff - GP_Child_Individual_SSE(i_GP_Individual) ) / dff
enddo

write(6,'(/A)') &
          'gpfr: i_GP_individual, GP_Individual_Ranked_Fitness '
do  i_GP_Individual=1,n_GP_Individuals
    write(6,'(I6,1x,E15.7)') &
          i_GP_individual, GP_Individual_Ranked_Fitness(i_GP_Individual)
enddo


! calculate the sum of the rankings

dff=0.0
do  i_GP_Individual=1,n_GP_Individuals
    dff=dff+GP_Individual_Ranked_Fitness(i_GP_Individual)
    GP_Integrated_Ranked_Fitness(i_GP_Individual)=dff
enddo

!do  i_GP_Individual=1,n_GP_Individuals
!    write(6,'(A,1x,I6,1x,E15.7)') &
!          'gpfr: i_GP_individual, GP_Integrated_Ranked_Fitness ', &
!                 i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual)
!enddo

! normalize to the integrated ranking values so that the ranking integration ranges from [0. to 1.]

write(6,'(/A)') &
      'gpfr: norm  i_GP_individual, GP_Integrated_Ranked_Fitness '

do  i_GP_Individual=1,n_GP_Individuals

    GP_Integrated_Ranked_Fitness(i_GP_Individual) = &
    GP_Integrated_Ranked_Fitness(i_GP_Individual)/GP_Integrated_Ranked_Fitness(n_GP_Individuals)

    write(6,'(I6,1x,E15.7)') &
          i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual)
enddo


!write(*,*) i_GP_Generation,'MAIN',GP_Individual_Ranked_Fitness(1)
write(*,'(//A,1x,I6,1x,E15.7/)') &
      'gpfr: MAIN, i_GP_Generation, GP_Individual_Ranked_Fitness(1) ', &
                   i_GP_Generation, GP_Individual_Ranked_Fitness(1)



!---------------------------------------------------------------------------------------------------

! find GP_best_parent


i_GP_Best_Parent=1
dff=GP_Individual_Ranked_Fitness(1)
do  i_GP_Individual=2,n_GP_individuals
    if( GP_Individual_Ranked_Fitness(i_GP_individual) .gt. dff) then
        dff=GP_Individual_Ranked_Fitness(i_GP_individual)
        i_GP_Best_Parent=i_GP_Individual
    endif
enddo

write(*,'(/A,2(1x,I6)/)') 'gpfr: i_GP_Generation, i_GP_Best_Parent ', &
                                i_GP_Generation, i_GP_Best_Parent


!---------------------------------------------------------------------------

! fill output array of parameters for best individual and write on GP_output_unit

do  i_CODE_equation=1,n_CODE_equations

    write(6,'(I6,1x, E20.10, 4x, E20.10)') &
          i_CODE_equation, &
          Runge_Kutta_Initial_Conditions(i_CODE_equation),  &
          GP_individual_initial_conditions(i_CODE_equation)

    output_array( i_CODE_equation ) = &
            GP_individual_initial_conditions(i_CODE_equation)

    ! ?? write(GA_output_unit,'(E24.16)') &
    ! ??       GP_individual_initial_conditions(i_CODE_equation)

enddo ! i_CODE_equation


nop = n_CODE_equations
do  i_tree=1,n_trees
    do  i_node=1,n_nodes
        if( abs( GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree) ) >  &
                                                                     1.0d-20   )then

            write(6,'(2(1x,I6), 1x, E20.10, 4x, E20.10)') &
                  i_node, i_tree, &
                  GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree)

             nop = nop + 1
             output_array(nop) = GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree)
        endif

    enddo ! i_node
enddo  ! i_tree


write( GP_output_unit, '(I6,1x,I6, 12(1x,E15.7))') &
       i_GP_Generation, i_GP_best_parent, &
       GP_individual_ranked_fitness(i_GP_Best_Parent), &
       output_array(1:nop) 

!---------------------------------------------------------------------------



GP_Adult_Individual_SSE  =  GP_Child_Individual_SSE

write(*,'(/A)') &
      'gpfr: i_GP_Individual, GP_Adult_Individual_SSE, GP_Individual_Ranked_Fitness '

do  i_GP_Individual=1,n_GP_individuals

    write(*,'(I6,2(1x,E15.7))') &
          i_GP_Individual, GP_Adult_Individual_SSE, GP_Individual_Ranked_Fitness

enddo ! i_GP_individual


!off if( i_GP_Generation .eq. 3) Stop

do i_GP_Individual=1,n_GP_individuals

   GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
   GP_Adult_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

   !off  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees)=GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)

enddo ! i_GP_individual



write(unit_gp_out) GP_Node_Type_for_Plotting


return


end subroutine GP_fitness_reset
