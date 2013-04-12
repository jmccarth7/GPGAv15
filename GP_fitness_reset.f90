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


!write(GP_print_unit,*) GP_Child_Individual_SSE
!if( myid == 0 )then
!    write(GP_print_unit,'(A,1x,I6)') 'gpfr: i_GP_generation ',  i_GP_generation
!
!    write(GP_print_unit,'(/A/(5(1x,E15.7)))') 'gpfr: GP_Child_Individual_SSE = ',&
!                                      GP_Child_Individual_SSE
!
!endif ! myid == 0
!---------------------------------------------------------------------------------------------------

! fitness reset region

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(A,1x,I6/)')     'gpfr: i_GP_generation ',  i_GP_generation
    
    do  i_GP_Individual=1,n_GP_Individuals
        write(GP_print_unit,'(A,2(1x,I6))') &
              'gpfr: i_GP_Individual, GP_Individual_N_GP_param(i_GP_Individual) ', &
                     i_GP_Individual, GP_Individual_N_GP_param(i_GP_Individual) 
    enddo
    
    write(GP_print_unit,'(/A)')' '

endif ! i_GP_generation ...

!---------------------------------------------------------------------------------------------------


! calculate the total population's SSE

dff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) <= 0 ) cycle

    dff=dff+GP_Child_Individual_SSE(i_GP_Individual)

enddo ! i_gp_individual

!---------------------------------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    do  i_GP_Individual=1,n_GP_Individuals
            write(GP_print_unit,'(A,1x,I6,1x,E15.7)') &
                  'gpfr: i_GP_Individual, GP_Child_Individual_SSE(i_GP_Individual) ', &
                         i_GP_Individual, GP_Child_Individual_SSE(i_GP_Individual) 
    enddo ! i_gp_individual

endif ! i_GP_generation ...

write(GP_print_unit,'(/A,1x,E15.7)') 'gpfr: sum of all GP_Child_Individual_SSE = ', dff

!---------------------------------------------------------------------------------------------------

! calculate a normalized ranking of the errors (higher individual SSE == lower value/ranking)

GP_Individual_Ranked_Fitness = 0.0d0

do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) <= 0 ) cycle

    if( abs( dff ) > 1.0D-20 )then

        GP_Individual_Ranked_Fitness(i_GP_Individual) = &
             ( dff - GP_Child_Individual_SSE(i_GP_Individual) ) / dff
    else
        GP_Individual_Ranked_Fitness(i_GP_Individual) = 0.0D0
    endif ! abs( dff ) > 1.0D-20
enddo

!---------------------------------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(/A)') &
          'gpfr: i_GP_individual, GP_Individual_Ranked_Fitness '
    
    do  i_GP_Individual=1,n_GP_Individuals
        write(GP_print_unit,'(I6,1x,E15.7)') &
              i_GP_individual, GP_Individual_Ranked_Fitness(i_GP_Individual)
    enddo

endif ! i_GP_generation == 1 .or. ...

!---------------------------------------------------------------------------------------------------

! calculate the sum of the rankings

GP_Integrated_Ranked_Fitness = 0.0d0

dff=0.0
do  i_GP_Individual=1,n_GP_Individuals

    dff = dff + GP_Individual_Ranked_Fitness(i_GP_Individual)

    GP_Integrated_Ranked_Fitness(i_GP_Individual)=dff
enddo

!do  i_GP_Individual=1,n_GP_Individuals
!    write(GP_print_unit,'(A,1x,I6,1x,E15.7)') &
!    'gpfr: i_GP_individual, GP_Integrated_Ranked_Fitness ', &
!           i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual)
!enddo

!---------------------------------------------------------------------------------------------------

! normalize to the integrated ranking values so that 
! the ranking integration ranges from [0. to 1.]


do  i_GP_Individual=1,n_GP_Individuals

    if( abs( GP_Integrated_Ranked_Fitness(n_GP_Individuals) ) > 1.0D-20 )then
        GP_Integrated_Ranked_Fitness(i_GP_Individual) = &
        GP_Integrated_Ranked_Fitness(i_GP_Individual) / &
                        GP_Integrated_Ranked_Fitness(n_GP_Individuals)
    else
        GP_Integrated_Ranked_Fitness(i_GP_Individual) = 0.0D0
    endif ! abs( GP_Integrated_Ranked_Fitness(n_GP_Individuals) ) > 1.0D-20

enddo

!---------------------------------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(/A)') &
          'gpfr: norm  i_GP_individual, GP_Integrated_Ranked_Fitness '
    
    do  i_GP_Individual=1,n_GP_Individuals
    
        write(GP_print_unit,'(I6,1x,E15.7)') &
              i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual)
    enddo

endif ! i_GP_generation == 1 .or. ...


!---------------------------------------------------------------------------------------------------

!write(GP_print_unit,*) i_GP_Generation,'MAIN',GP_Individual_Ranked_Fitness(1)
write(GP_print_unit,'(//A,1x,I6,1x,E15.7/)') &
      'gpfr: MAIN, i_GP_Generation, GP_Individual_Ranked_Fitness(1) ', &
                   i_GP_Generation, GP_Individual_Ranked_Fitness(1)

!---------------------------------------------------------------------------------------------------

! find GP_best_parent


i_GP_Best_Parent=1

dff=GP_Individual_Ranked_Fitness(1)

do  i_GP_Individual=2,n_GP_individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) <= 0 ) cycle

    if( GP_Individual_Ranked_Fitness(i_GP_individual) .gt. dff) then
        dff=GP_Individual_Ranked_Fitness(i_GP_individual)
        i_GP_Best_Parent=i_GP_Individual
    endif

enddo ! i_GP_Individual

write(GP_print_unit,'(/A,2(1x,I6),1x,E15.7/)') &
      'gpfr: i_GP_Gen, i_GP_Best_Parent, GP_Indiv_Ranked_Fit(i_GP_Best_Parent) ', &
             i_GP_Generation, i_GP_Best_Parent, &
             GP_Individual_Ranked_Fitness(i_GP_Best_Parent)


!---------------------------------------------------------------------------

! fill output array of parameters for best individual and write on GP_output_unit

write(GP_print_unit,'(/A)') &
          'gpfr: i_CODE_eq  R_K_Init_Cond(i_CODE_eq)  &
          &GP_indiv_init_cond(i_CODE_eq)'
do  i_CODE_equation=1,n_CODE_equations

    write(GP_print_unit,'(6x,I6,7x, E20.10, 10x, E20.10)') &
          i_CODE_equation, &
          Runge_Kutta_Initial_Conditions(i_CODE_equation),  &
          GP_individual_initial_conditions(i_CODE_equation)

    output_array( i_CODE_equation ) = &
            GP_individual_initial_conditions(i_CODE_equation)

    ! ?? write(GA_output_unit,'(E24.16)') &
    ! ??       GP_individual_initial_conditions(i_CODE_equation)

enddo ! i_CODE_equation


write(GP_print_unit,'(/A)') &
          'gpfr: i_node  itree  nop  GP_pop_node_params(i_GP_Best_Parent,i_node,i_tree)'
nop = n_CODE_equations
do  i_tree=1,n_trees
    do  i_node=1,n_nodes
        if( abs( GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree) ) >  &
                                                                     1.0d-20   )then
            nop = nop + 1
            write(GP_print_unit,'(2x,3(1x,I6), 1x, E20.10, 4x, E20.10)') &
                  i_node, i_tree, nop, &
                  GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree)

            output_array(nop) = GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree)

        endif !   abs( GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree) ) >...

    enddo ! i_node
enddo  ! i_tree


write( GP_output_unit, '(I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
       i_GP_Generation, i_GP_best_parent, &
       GP_individual_ranked_fitness(i_GP_Best_Parent), &
       nop, output_array(1:nop) 

write(GP_print_unit, '(//A,1x,I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
       'gpfr: i_GP_gen, i_GP_best_parent, GP_indiv_ranked_fit, output_array) ', &
       i_GP_Generation, i_GP_best_parent, &
       GP_individual_ranked_fitness(i_GP_Best_Parent), &
       nop, output_array(1:nop) 

!---------------------------------------------------------------------------



GP_Adult_Individual_SSE  =  GP_Child_Individual_SSE



!---------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(/A)') &
          'gpfr: i_GP_Indiv, GP_Adult_Indiv_SSE, GP_Indiv_Ranked_Fitness '
    
    do  i_GP_Individual=1,n_GP_individuals
    
        write(GP_print_unit,'(6x,I6,2(7x,E15.7))') &
              i_GP_Individual, &
              GP_Adult_Individual_SSE(i_GP_Individual), &
              GP_Individual_Ranked_Fitness(i_GP_Individual)
    
    enddo ! i_GP_individual

endif ! i_GP_generation == 1 .or. ...

!off if( i_GP_Generation .eq. 3) Stop

!-----------------------------------------------------------------------------------------

do i_GP_Individual=1,n_GP_individuals

   GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
        GP_Adult_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

   !off  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
   !off  GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)

enddo ! i_GP_individual



write(unit_gp_out) GP_Node_Type_for_Plotting


!-----------------------------------------------------------------------------------------

return


end subroutine GP_fitness_reset
