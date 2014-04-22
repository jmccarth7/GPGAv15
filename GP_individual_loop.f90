subroutine GP_individual_loop( new_group, new_comm, i_GP_generation )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module

use fasham_variables_module
use Tree_Node_Factory_module
use class_Tree_Node


implicit none




integer :: i
!integer :: i_diversity
integer :: message_len

integer(kind=4),intent(in) :: new_group
integer(kind=4),intent(in) :: new_comm 
integer(kind=4),intent(in) :: i_GP_Generation

integer(kind=4) :: i_GP_individual
integer(kind=4) :: GP_minSSE_Individual
integer(kind=4) :: GP_minSSE_generation
integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: jj
!!integer(kind=4) :: nn

integer(kind=4) :: i_CODE_equation
integer(kind=4) :: max_n_gp_params

integer(kind=4) :: n_GP_vars
integer(kind=4) :: nop

integer(kind=4) :: i_GP_best_parent
integer(kind=4) :: ierror
integer(kind=4) :: ierror_t
integer(kind=4) :: ierror_m
integer(kind=4) :: ierror_tb
integer(kind=4) :: i_start_generation

!integer(kind=4) :: new_rank 
integer(kind=4) :: n_procs  

integer(kind=4) :: i_part
integer(kind=4) :: i_gp_1
integer(kind=4) :: i_gp_2
integer(kind=4) :: ind1  
integer(kind=4) :: ind2  


!real(kind=8) :: t1
!real(kind=8) :: t2


!---------------------------------------------------------------------------------------

divider = ( numprocs - 1 ) / n_partitions

call mpi_comm_rank( new_comm, new_rank, ierr )  
call mpi_comm_size( new_comm, n_procs,  ierr )  

if( myid == 0 )then
!if( new_rank == 0 )then
    write(GP_print_unit,'(/A,4(1x,I6))')&
     'gil: before loop myid, new_rank, n_code_equations,  n_GP_params ',&
                       myid, new_rank, n_code_equations,  n_GP_parameters
    write(GP_print_unit,'(A,5(1x,I6))')&
     'gil: before loop myid, new_rank, numprocs, divider, n_partitions', &
                       myid, new_rank, numprocs, divider, n_partitions
    write(GP_print_unit,'(A,3(1x,I6))')&
     'gil: before loop myid, new_rank, n_procs                        ', &
                       myid, new_rank, n_procs              
endif !  myid == 0

!!!!!!!!!!!!!!!!!!!!!!!!!if( myid == 0 )return

! do the loop over the GP individuals in n_partitions chunks
part_loop:&
do  i_part = 1, n_partitions

    ! i_gp_1 and i_gp_2 are limits on the processor number for this partition

    i_gp_1 = divider * (i_part - 1) + 1
    i_gp_2 = divider *  i_part

    if( i_gp_1 <= myid  .and.   &
        myid   <= i_gp_2         )then

        write(GP_print_unit,'(A,5(1x,I6))')&
         'gil:in loop myid, new_rank, i_part, i_gp_1, i_gp_2 ',&
                      myid, new_rank, i_part, i_gp_1, i_gp_2

        ! ind1 and ind2 are limits on the i_GP_individuals processed in this partition

        ind1 =  (n_GP_individuals / n_partitions)  * (i_part-1) + 1
        ind2 =  (n_GP_individuals / n_partitions)  *  i_part

        ! get any remaining individuals in the last partition

        if( i_part == n_partitions )then
            ind2 = n_GP_individuals
        endif ! i_part == n_partitions 

        write(GP_print_unit,'(A,5(1x,I6))')&
         'gil:in loop myid, new_rank, i_part, ind1, ind2     ',&
                      myid, new_rank, i_part, ind1, ind2  


        gp_ind_loop:&
        do  i_GP_individual= ind1, ind2    ! 1,n_GP_individuals
        
        
            !------------------------------------------------------------------------------
            ! calculate how many parameters total to fit for the specific individual CODE
            ! and save this number in GP_Individual_N_GP_param(i_GP_individual)
        
            n_GP_Parameters = n_code_equations
        
            if( new_rank == 0 )then
                write(GP_print_unit,'(A,2(1x,I6))')&
                 'gil: before loop n_code_equations,  n_GP_params',&
                                   n_code_equations,  n_GP_parameters
            endif !  new_rank == 0
        
            do  i_Tree=1,n_Trees
        
                do  i_Node=1,n_Nodes
        
                    if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0) then
        
                        n_GP_Parameters = n_GP_Parameters+1
        
                        !if( new_rank == 0 )then
                        !    write(GP_print_unit,'(A,5(1x,I6))')&
                        !    'gil: i_GP_indiv, i_tree, i_node, GP_Adult_Pop_Node_Type, n_GP_params',&
                        !        i_GP_individual, i_tree, i_node, &
                        !        GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual), &
                        !        n_GP_parameters
                        !endif !  new_rank == 0
        
                    endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)
        
                    !if( new_rank == 0 )then
                    !    if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) > -9999)then
                    !        write(GP_print_unit,'(A,5(1x,I6))')&
                    !        'gil: i_GP_indiv, i_tree, i_node, GP_Adult_Pop_Node_Type, n_GP_params',&
                    !            i_GP_individual, i_tree, i_node, &
                    !            GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual), &
                    !            n_GP_parameters
                    !    endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) > -9999
                    !endif !  new_rank == 0
        
                enddo ! i_node
        
            enddo ! i_tree
        
        
            if( new_rank == 0 )then
                write(GP_print_unit,'(A,2(1x,I6))')&
                 'gil: after loop n_code_equations,  n_GP_params',&
                                  n_code_equations,  n_GP_parameters
                !flush(GP_print_unit)
            endif !  new_rank == 0
        
            !------------------------------------------------------------------------
        
            ! if n_GP_parameters > n_maximum_number_parameters, this is an error
            ! for now, just print a warning and set n_GP_parameters = n_max...
        
            if( n_GP_parameters > n_maximum_number_parameters )then
        
                if( new_rank == 0 )then
                    write(GP_print_unit,'(/A,3(1x,I5),2(1x,I10))')&
                      'gil:WARNING  new_rank, i_GP_generation, i_GP_Individual, &
                          &n_GP_parameters, n_maximum_number_parameters', &
                          new_rank, i_GP_generation, i_GP_Individual, &
                           n_GP_parameters, n_maximum_number_parameters
                    !flush(GP_print_unit)
                endif !  new_rank == 0
        
                !call MPI_FINALIZE(ierr)
                !stop 'too many parms'
        
        
                GP_Child_Individual_SSE(i_GP_Individual) = 1.0D13
        
                GP_Adult_Individual_SSE(i_GP_Individual) = 1.0D13
                GP_Adult_Population_SSE(i_GP_Individual) = 1.0D13
        
                cycle gp_ind_loop
        
            endif  ! n_GP_parameters > n_maximum_number_parameters
        
            !------------------------------------------------------------------------
        
            GP_Individual_N_GP_param(i_GP_individual) = n_GP_parameters
        
            !------------------------------------------------------------------------------
        
        
            ! run GPCODE_... to evaluate this individual  if Run_GP_Calculate_Fitness is true
        
        
        
            if( Run_GP_Calculate_Fitness(i_GP_Individual) ) then
        
                !-----------------------------------------------------------------------------------
        
                ! these get set randomly in the GA-lmdif search algorithm ( in GPCODE* )
                !GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) = 0.0d0               ! 20131209
        
                !-----------------------------------------------------------------------------------
                if( new_rank == 0 )then
                    write(GP_print_unit,'(/A)')&
                    'gil:----------------------------------------------------------------------'
                    write(GP_print_unit,'(/A,2(1x,I6),3x,L1)') &
                          'gil: i_GP_Gen, i_GP_indiv, Run_GP_Calculate_Fitness', &
                                i_GP_Generation, i_GP_individual, &
                                       Run_GP_Calculate_Fitness(i_GP_Individual)
                    write(GP_print_unit,'(A)')&
                    'gil:----------------------------------------------------------------------'
                    !!flush(GP_print_unit)
                endif !  new_rank == 0
                !-----------------------------------------------------------------------------------
        
                if( new_rank == 0 )then
                    write(GP_print_unit,'(/A,4(1x,I6))') &
                      'gil: i_GP_individual, n_trees, n_nodes, n_GP_parameters ', &
                            i_GP_individual, n_trees, n_nodes, n_GP_parameters
                    !flush(GP_print_unit)
                endif !  new_rank == 0
        
                !-----------------------------------------------------------------------------------
        
                do  i_Tree=1,n_Trees
                    do  i_Node=1,n_Nodes
        
                        GP_Individual_Node_Type(i_Node,i_Tree) = &
                           GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)
        
                        GP_Individual_Node_Parameters(i_Node,i_Tree) = 0.0d0
                    enddo ! i_node
                enddo ! i_tree
        
                !-----------------------------------------------------------------------------------
        
                ! calculate how many variables are in the tree
        
                n_GP_vars = 0
                do  i_Tree=1,n_Trees
                    do  i_Node=1,n_Nodes
        
                        if( GP_Individual_Node_Type(i_Node,i_Tree) < 0  .and. &
                            GP_Individual_Node_Type(i_Node,i_Tree) > -9999  ) then
                            n_GP_vars = n_GP_vars + 1
                        endif ! GP_Individual_Node_Type(i_Node,i_Tree) > 0 ....
        
                        !if( new_rank == 0 .and. GP_Individual_Node_Type(i_Node,i_Tree) > -9999 )then
                        !    write(GP_print_unit,'(A,5(1x,I6))')&
                        !    'gil: i_GP_indiv, i_tree, i_node, GP_Indiv_Node_Type, n_GP_vars ', &
                        !          i_GP_individual, i_tree, i_node, &
                        !          GP_Individual_Node_Type(i_Node,i_Tree), n_GP_vars
                        !endif !  new_rank == 0 .and. GP_Individual_Node_Type > -9999
        
                    enddo ! i_node
                enddo ! i_tree
        
                !------------------------------------------------------------------------------
        
                if( new_rank == 0 )then
                    write(GP_print_unit,'(A,1x,I6,5x,A,2(1x,I6))')&
                          'gil: for i_GP_Individual', i_GP_Individual, &
                          'number of parameters, variables =', &
                                n_GP_parameters, n_GP_vars
                    !flush(GP_print_unit)
                endif !  new_rank == 0
        
        
                !-------------------------------------------------------------------
        
                ! cycle the i_GP_individual loop if there are no GP parameters
                ! or if n_GP_parameters <=  n_code_equations
        
                if( n_GP_parameters == 0 .or. &
                    n_GP_parameters > n_maximum_number_parameters .or.  &
                    n_GP_parameters <=  n_code_equations                 ) then   ! new jjm 20130814
        
                    if( new_rank == 0 )then
                        write(GP_print_unit,'(A,1x,I6)')&
                              'gil: skipping this i_GP_Individual --&
                              &  the number of parameters is ', n_GP_parameters
                        !flush(GP_print_unit)
                    endif !  new_rank == 0
        
        
                    ! set SSE values for this rejected individual so that
                    ! its fitness will be very small
        
                    GP_Child_Individual_SSE(i_GP_Individual) = 1.0D13
        
                    GP_Adult_Individual_SSE(i_GP_Individual) = 1.0D13
                    GP_Adult_Population_SSE(i_GP_Individual) = 1.0D13
        
        
                    if( new_rank == 0 )then
                        write(GP_print_unit,'(/A/9x,A/A,2(1x,I5), 1x, E15.7)')&
                              'gil: rejected for n_GP_parameters <=  n_code_equations',&
                              'or for n_GP_parameters >   n_maximum_number_parameters',&
                               'i_GP_gen,i_GP_indiv,GP_Child_Pop_SSE(i_GP_Indiv)  ', &
                                   i_GP_generation, i_GP_individual, &
                                   GP_Child_Individual_SSE(i_GP_Individual)
                    endif !  new_rank == 0
        
                    cycle gp_ind_loop
        
                endif ! n_GP_parameters == 0
        
                !-------------------------------------------------------------------
        
        
                if( new_rank == 0 )then
        
                    write(GP_print_unit,'(/A,2(1x,I6))') &
                     'gil: call GPCODE_GA_lmdif_Param_Opt         i_GP_Gen, i_GP_indiv', &
                                               i_GP_Generation, i_GP_individual
                    !!flush(GP_print_unit)
        
                    if( L_ga_print )then
                        write(GA_print_unit,'(//A/A)') &
                         '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>&
                         &>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
                         'gil: call GPCODE_GA_lmdif_Parameter_Optimization routine'
                        write(GA_print_unit,'(A,2(1x,I6)/)') &
                         'gil: i_GP_Generation, i_GP_individual',&
                               i_GP_Generation, i_GP_individual
                    endif ! L_ga_print
        
                endif ! new_rank == 0
        
        
                !----------------------------------------------------------------------------
        
        
                ! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND
                ! LINK THE SSE OUTPUT TO THE ARRAY AT THE END
                ! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEED TO BE PULLED OUT
        
                !      individual_fitness
                !   GP_Individual_Initial_Conditions
                !   GP_Individual_Node_Parameters
        
                ! these arrays are broadcast in GPCODE_GA...
        
        
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
                call GPCODE_GA_lmdif_Parameter_Optimization( &
                                  i_GP_Generation,i_GP_individual, &
                                  new_group, new_comm, ind2-ind1+1  )
        
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
                !if( myid == 0 )then
                    write(GP_print_unit,'(A,1x,I6)') &
                     'gil: aft call GPCODE_GA_lmdif_Parameter_Optimization routine'
                    !flush(GP_print_unit)
                !endif ! myid == 0


        
                GP_population_fitness(i_GP_individual) = individual_fitness
        


                !if( myid == 0 )then

                    write(GP_print_unit,'(/A,4(1x,I6), 1x, E15.7)')&
                      'gil: myid, new_rank, i_GP_gen, i_GP_indiv, GP_pop_fit(i_GP_indiv) ', &
                            myid, new_rank, i_GP_generation, i_GP_individual, &
                            GP_population_fitness(i_GP_individual)
        
                !endif ! myid == 0 

                !if( myid == 0 )then
                    write(GP_print_unit,'(/A,2(1x,I6), 1x, E15.7)')&
                          'gil:1 i_GP_gen, i_GP_indiv, GP_pop_fit(i_GP_indiv) ', &
                               i_GP_generation, i_GP_individual, &
                               GP_population_fitness(i_GP_individual)
                !endif !  myid == 0

    
                if( myid == 0 )then
    
                    ! print side-by-side comparisons of
                    ! starting values and values from optimization
        
                    !write(GP_print_unit,'(/A/)') &
                    !  'gil:           truth value           GP_individual_initial_conditions'
        
                    do  i_CODE_equation=1,n_CODE_equations
        
                        !write(GP_print_unit,'(I6,1x, E20.10, 4x, E20.10)') &
                        !      i_CODE_equation, &
                        !      Runge_Kutta_Initial_Conditions(i_CODE_equation),  &
                        !      GP_individual_initial_conditions(i_CODE_equation)
        
                        if( L_GA_output_parameters )then
                            write(GA_output_unit,'(E15.7)') &
                                  GP_individual_initial_conditions(i_CODE_equation)
                        endif ! L_GA_output_parameters )then
        
                    enddo ! i_CODE_equation
        
                endif !  myid == 0
        
                !--------------------------------------------------------------------------------
        
                ! compute GP_Child_Individual_SSE(i_GP_Individual)
        
                ! use the best sse from the GPCODE subroutine
        
                GP_Child_Individual_SSE(i_GP_Individual) = Individual_SSE_best_parent
        
                GP_Adult_Individual_SSE(i_GP_Individual) = Individual_SSE_best_parent
                GP_Adult_Population_SSE(i_GP_Individual) = Individual_SSE_best_parent
        
                if( myid == 0 )then
                    !if( i_GP_generation == 1                                  .or. &
                    !    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
                    !    i_GP_generation == n_GP_generations                          )then
        
                        write(GP_print_unit,'(A,4(1x,I5), 1x, E12.5,2x,G10.3,2x,F6.3)') &
                          'gil:3 i_GP_gen, i_indiv, nparms, nvars, SSE, SSE/SSE0, GP_pop_fit', &
                                 i_GP_generation, i_GP_individual, &
                                 n_GP_parameters, n_GP_vars, &
                                 GP_Child_Individual_SSE(i_GP_Individual), &
                                 GP_Child_Individual_SSE(i_GP_Individual)/SSE0, &
                                 GP_population_fitness(i_GP_individual)
                        !!flush(GP_print_unit)
        
                    !endif ! i_GP_generation == 1 .or. ...
                endif !  myid == 0
        
                !-----------------------------------------------------------------------------------
        
        
                ! set the GA_lmdif-optimized initial condition array
        
                do  i_code_equation = 1, n_code_equations  ! 20131209
        
                    GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual) = &
                        GP_Individual_Initial_Conditions(i_CODE_Equation)
        
                enddo ! i_code_equation
        
                !-----------------------------------------------------------------------------------
        
        
                ! set the GA_lmdif-optimized CODE parameter set array
        
                do  i_Tree=1,n_Trees    ! 20131209
                    do  i_Node=1,n_Nodes
        
                        GP_Population_Node_Parameters(i_node, i_tree, i_GP_Individual) = &
                               GP_Individual_Node_Parameters(i_node, i_tree)
        
                        GP_Adult_Population_Node_Type(i_node, i_tree, i_GP_Individual) = &
                                         GP_Individual_Node_Type(i_node, i_tree)
        
                        GP_Child_Population_Node_Type(i_node, i_tree, i_GP_Individual) = &
                                         GP_Individual_Node_Type(i_node, i_tree)
        
                    enddo ! i_node
                enddo ! i_tree
        
                !------------------------------------------------------------------------------
        
        
        
                !if( myid == 0 )then
        
                    !!----------------------------------------------------------------------------
                    !! print the node parameters (if there are any)
                    !write(GP_print_unit,'(/A/)') &
                    !       'gil:  tree  node  Runge_Kutta_Node_Params&
                    !                  &   GP_population_node_params'
                    !do  i_tree=1,n_trees
                    !    do  i_node=1,n_nodes
                    !        if( GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0 ) then
                    !            write(GP_print_unit,'(2(1x,I6), 1x, E20.10, 4x, E20.10)') &
                    !             i_tree, i_node, &
                    !             GP_individual_Node_Parameters(i_node,i_tree), &
                    !             GP_population_node_parameters(i_node,i_tree,i_GP_individual)
                    !        endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0
                    !    enddo ! i_node
                    !enddo  ! i_tree
                    !!---------------------------------------------------------------------------------
        
        
                    !write(GP_print_unit,'(/A)')  &
                    !      'gil: after loading GP_Pop arrays with GP_indiv array values '
        
        
                !endif !  myid == 0
        
            endif !   Run_GP_Calculate_Fitness(i_GP_Individual)
        
        
            if( myid == 0 )then
        
                if( L_GP_all_summary )then
        
                    ! this prints a summary of the initial conditions,
                    ! parameters,  and node types for this individual,
                    ! after being optimized in GPCODE*opt
                    ! and writes the tree to the summary file
        
                    !if( i_GP_generation == 1                                  .or. &
                    !    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
                    !    i_GP_generation == n_GP_generations                          )then
                    !
                        write(GP_print_unit,'(/A)') &
                          'gil:------------------------------------------&
                           &-----------------------------'
                        write(GP_print_unit,'(A,2(1x,I6))') &
                        'gil: call summary_GP_indiv i_GP_generation, i_GP_individual ', &
                                                    i_GP_generation, i_GP_individual
                     
                    !endif ! i_GP_generation == 1 .or. ...
        
                    call summary_GP_indiv( i_GP_generation, i_GP_individual, 0 )
        
                endif ! L_GP_all_summary
        
            endif !  myid == 0
        
        
        enddo  gp_ind_loop    !   i_GP_individual

    endif ! i_gp_1 <= myid  .and. ...

call MPI_FINALIZE(ierr)  !debug only                                                                                            
stop                     !debug only                                                                                            

enddo  part_loop

                                                                                                                                
call MPI_BARRIER( MPI_COMM_WORLD, ierr )  
!call MPI_BARRIER( new_comm , ierr )  
write(6,'(/A,4(1x,I6), 1x, E15.7/)')&
              'gil:  AFT BARRIER myid, new_rank ', &
                                 myid, new_rank

if( myid == 0 )then
    do  i_GP_individual = 1, n_GP_individuals
        write(GP_print_unit,'(A,4(1x,I6), 1x, E15.7)')&
              'gil:1 myid, new_rank, i_GP_gen, i_GP_indiv, GP_pop_fit(i_GP_indiv) ', &
                     myid, new_rank, i_GP_generation, i_GP_individual, &
                                GP_population_fitness(i_GP_individual)
    enddo ! i_GP_individual
endif !  myid == 0


return



end subroutine GP_individual_loop
