subroutine GP_calc_fitness( output_array )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations

! inputs:

! GP_Child_Individual_SSE


! outputs:

! GP_Adult_Individual_SSE
! GP_Population_Ranked_Fitness
! GP_Integrated_Population_Ranked_Fitness
! i_GP_Best_Parent
! output_array


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!use mpi
!use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module


implicit none



integer :: nop



real(kind=8), dimension(n_maximum_number_parameters) :: output_array

real (kind=8) ::  dff


!----------------------------------------------------------------------------------------

! this routine is only called by processor 0


!allocate( output_array( n_maximum_number_parameters ) )

output_array = 0.0d0

!!!!write(GP_print_unit,*) GP_Child_Individual_SSE

write(GP_print_unit,'(A,1x,I6)') 'gpcf: i_GP_generation ',  i_GP_generation

write(GP_print_unit,'(/A/(5(1x,E15.7)))') 'gpcf: GP_Child_Individual_SSE = ',&
                                                 GP_Child_Individual_SSE


!---------------------------------------------------------------------------------------------------

! fitness reset region (??)

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(A,1x,I6/)') &
          'gpcf: i_GP_generation ',  i_GP_generation
    
    write(GP_print_unit,'(A)') &
          'gpcf: i_GP_Indiv, GP_Indiv_N_GP_param(i_GP_Indiv) '

    do  i_GP_Individual=1,n_GP_Individuals
        write(GP_print_unit,'(5x,2(1x,I6))') &
                     i_GP_Individual, GP_Individual_N_GP_param(i_GP_Individual) 
    enddo
    
    write(GP_print_unit,'(/A)')' '

endif ! i_GP_generation ...

!-------------------------------------------------------------------------------------------------


! calculate the total population's SSE

dff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) <= 0 ) cycle

    dff=dff+GP_Child_Individual_SSE(i_GP_Individual)

enddo ! i_gp_individual

!-------------------------------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(A)') &
              'gpcf: i_GP_Indiv, GP_Child_Indiv_SSE(i_GP_Indiv) '
    do  i_GP_Individual=1,n_GP_Individuals
        write(GP_print_unit,'(6x,I6,6x,E15.7)') &
                     i_GP_Individual, GP_Child_Individual_SSE(i_GP_Individual) 
    enddo ! i_gp_individual

endif ! i_GP_generation ...

write(GP_print_unit,'(/A,1x,E15.7)') 'gpcf: sum of all GP_Child_Indiv_SSE = ', dff

!-------------------------------------------------------------------------------------------------

! calculate a normalized ranking of the errors (higher individual SSE == lower value/ranking)

GP_Population_Ranked_Fitness = 0.0d0

do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) <= 0 ) cycle

    if( abs( dff ) > 1.0D-20 )then

        GP_Population_Ranked_Fitness(i_GP_Individual) = &
             ( dff - GP_Child_Individual_SSE(i_GP_Individual) ) / dff
    else
        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0
    endif ! abs( dff ) > 1.0D-20

enddo ! i_GP_Individual

!-------------------------------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(/A)') &
          'gpcf: i_GP_indiv, GP_Indiv_Ranked_Fitness '
    
    do  i_GP_Individual=1,n_GP_Individuals
        write(GP_print_unit,'(5x,I6,5x,E15.7)') &
              i_GP_individual, GP_Population_Ranked_Fitness(i_GP_Individual)
    enddo

endif ! i_GP_generation == 1 .or. ...

!-------------------------------------------------------------------------------------------------

! calculate the sum of the rankings

GP_Integrated_Population_Ranked_Fitness = 0.0d0

dff=0.0
do  i_GP_Individual=1,n_GP_Individuals

    dff = dff + GP_Population_Ranked_Fitness(i_GP_Individual)

    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)=dff
enddo

write(GP_print_unit,'(/A)') &
    'gpcf: i_GP_indiv, GP_Integrated_Population_Ranked_Fitness '
do  i_GP_Individual=1,n_GP_Individuals
    write(GP_print_unit,'(5x,I6,5x,E15.7)') &
           i_GP_individual, GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)
enddo

!-------------------------------------------------------------------------------------------------

! normalize to the integrated ranking values so that 
! the ranking integration ranges from [0. to 1.]


do  i_GP_Individual=1,n_GP_Individuals

    if( abs( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) ) > 1.0D-20 )then

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = &
        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) / &
                        GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

    else

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0

    endif ! abs( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) ) > 1.0D-20

enddo ! i_GP_Individual

!-------------------------------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(/A)') &
          'gpcf: i_GP_indiv, GP_Integrated_Population_Ranked_Fitness (norm) '
    
    do  i_GP_Individual=1,n_GP_Individuals
    
        write(GP_print_unit,'(5x,I6,5x,E15.7)') &
              i_GP_individual, GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)
    enddo

endif ! i_GP_generation == 1 .or. ...


!-------------------------------------------------------------------------------------------------

!write(GP_print_unit,*) i_GP_Generation,'MAIN',GP_Population_Ranked_Fitness(1)
write(GP_print_unit,'(/A,17x,A,15x,I6,8x,E15.7)') &
      'gpcf: i_GP_Gen,', 'GP_Pop_Ranked_Fit(1) ', &
             i_GP_Generation, GP_Population_Ranked_Fitness(1)

!-------------------------------------------------------------------------------------------------

! find GP_best_parent


i_GP_Best_Parent=1

dff=GP_Population_Ranked_Fitness(1)

do  i_GP_Individual=2,n_GP_individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) <= 0 ) cycle

    if( GP_Population_Ranked_Fitness(i_GP_individual) .gt. dff) then
        dff=GP_Population_Ranked_Fitness(i_GP_individual)
        i_GP_Best_Parent=i_GP_Individual
    endif

enddo ! i_GP_Individual

write(GP_print_unit,'(A,2(1x,I6),1x,E15.7/)') &
      'gpcf: i_GP_Gen,i_GP_Best_Parent,&
            &GP_Pop_Ranked_Fit(i_GP_Best_Parent)', &
             i_GP_Generation, i_GP_Best_Parent, &
             GP_Population_Ranked_Fitness(i_GP_Best_Parent)


!---------------------------------------------------------------------------

! fill output array of parameters for best individual and write on GP_output_unit


write(GP_print_unit,'(/A)') &
          'gpcf: i_CODE_eq  R_K_Init_Cond(i_CODE_eq)  &
          &GP_Pop_init_cond(i_GP_Best_Parent,i_CODE_eq)'
do  i_CODE_equation=1,n_CODE_equations

    write(GP_print_unit,'(6x,I6,7x, E20.10, 10x, E20.10)') &
          i_CODE_equation, &
          Runge_Kutta_Initial_Conditions(i_CODE_equation),  &
          GP_Population_Initial_Conditions(i_GP_Best_Parent, i_CODE_equation )
          !!!GP_individual_initial_conditions(i_CODE_equation)

    output_array( i_CODE_equation ) = &
            GP_Population_Initial_Conditions(i_GP_Best_Parent, i_CODE_equation )
            !GP_individual_initial_conditions(i_CODE_equation)

    ! ?? write(GA_output_unit,'(E24.16)') &
    ! ??       GP_individual_initial_conditions(i_CODE_equation)

enddo ! i_CODE_equation


nop = n_CODE_equations

write(GP_print_unit,'(/A,2(1x,I6)/)') &
      'gpcf: before tree loop n_code_equations, nop ', &
                              n_code_equations, nop
write(GP_print_unit,'(/A)') &
     'gpcf: i_node  itree  nop  GP_pop_node_params(i_GP_Best_Parent,i_node,i_tree)'

nop = nop + 1

tree_loop:&
do  i_tree=1,n_trees
    node_loop:&
    do  i_node=1,n_nodes

        !if( abs( GP_population_node_parameters( &
        !              i_GP_Best_Parent,i_node,i_tree) ) >  1.0d-20   )then

        if( GP_Adult_Population_Node_Type(i_GP_Best_Parent,i_Node,i_Tree) == 0 )then

            write(GP_print_unit,'(2x,3(1x,I6), 1x, E20.10, 4x, E20.10)') &
                  i_node, i_tree, nop, &
                  GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree)

        !endif ! abs( GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree) ) >...
            output_array(nop) = &
                   GP_population_node_parameters(i_GP_Best_Parent,i_node,i_tree)

            nop = nop + 1

        endif ! GP_Adult_Pop_Node_Type(i_GP_Best_Parent,i_Node,i_Tree) == 0

        !write(GP_print_unit,'(3(1x,I6))') i_tree, i_node, nop

        if( nop > n_maximum_number_parameters ) then
            write(GP_print_unit,'(A)') &
                  'gpcf: nop >  n_maximum_number_parameters  '
            write(GP_print_unit,'(A,3(1x,I6))') &
                  'gpcf: nop, n_maximum_number_parameters  ', &
                         nop, n_maximum_number_parameters 
            nop = min( nop, n_maximum_number_parameters )
            exit tree_loop
        endif  ! nop > ...



    enddo node_loop ! i_node
enddo tree_loop ! i_tree


write( GP_output_unit, '(I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
       i_GP_Generation, i_GP_best_parent, &
       GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
       nop, output_array(1:nop) 

!write(GP_print_unit, '(//A,1x,I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
!       'gpcf: i_GP_gen,i_GP_best_parent,GP_indiv_ranked_fit, nop, output_array', &
!       i_GP_Generation, i_GP_best_parent, &
!       GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
!       nop, output_array(1:nop) 

!---------------------------------------------------------------------------



GP_Adult_Individual_SSE  =  GP_Child_Individual_SSE



!---------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then 

    write(GP_print_unit,'(/A)') &
          'gpcf: i_GP_Indiv, GP_Adult_Indiv_SSE, GP_Pop_Ranked_Fitness '
    
    do  i_GP_Individual=1,n_GP_individuals
    
        write(GP_print_unit,'(6x,I6,2(7x,E15.7))') &
              i_GP_Individual, &
              GP_Adult_Individual_SSE(i_GP_Individual), &
              GP_Population_Ranked_Fitness(i_GP_Individual)
    
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


end subroutine GP_calc_fitness
