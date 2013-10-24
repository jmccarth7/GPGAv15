subroutine GP_calc_fitness( i_GP_generation, output_array, &
                            i_GP_best_parent, nop )

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

integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_GP_Best_Parent
integer(kind=4),intent(in) :: i_GP_Generation
integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node


integer(kind=4) :: i_CODE_equation

real(kind=8), dimension(n_maximum_number_parameters), intent(out) :: output_array

real (kind=8) ::  dff

real (kind=8) ::  mean_fit
real (kind=8) ::  rms_fit
real (kind=8) ::  std_dev_fit

!-------------------------------------------------------------------------------

! this routine is only called by processor 0

! fitness reset region (??)

output_array = 0.0d0


write(GP_print_unit,'(/A,1x,I6/)') &
      'gpcf: i_GP_generation ',  i_GP_generation
write(GP_print_unit,'(/A/)') &
  'gpcf: i_GP_individual, GP_Child_Individual_SSE(i_GP_individual) '

do  i_GP_individual = 1, n_GP_individuals

    write(GP_print_unit,'(I6,1x,E15.7)')  &
          i_GP_individual, GP_Child_Individual_SSE(i_GP_individual)

enddo ! i_GP_individual

!-------------------------------------------------------------------------------


if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then

    write(GP_print_unit,'(/A,1x,I6/)') &
          'gpcf: i_GP_generation ',  i_GP_generation

    write(GP_print_unit,'(A)') &
          'gpcf: i_GP_Indiv, GP_Indiv_N_GP_param(i_GP_Indiv) '

    !do  i_GP_Individual=1,n_GP_Individuals
    !    write(GP_print_unit,'(5x,2(1x,I6))') &
    !          i_GP_Individual, GP_Individual_N_GP_param(i_GP_Individual)
    !enddo

    write(GP_print_unit, '(5(1x,I8," :", I3))') &
           ( i_GP_Individual, GP_Individual_N_GP_param(i_GP_Individual), &
             i_GP_Individual=1,n_GP_Individuals )


    write(GP_print_unit,'(/A)')' '

endif ! i_GP_generation ...

!-------------------------------------------------------------------------------


! calculate the total population's SSE

dff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) cycle

    dff=dff+GP_Child_Individual_SSE(i_GP_Individual)

enddo ! i_gp_individual

!write(GP_print_unit,'(/A,1x,E15.7)') 'gpcf: sum of all GP_Child_Indiv_SSE = ', dff

!-------------------------------------------------------------------------------

!if( i_GP_generation == 1                                 .or. &
!    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
!    i_GP_generation == n_GP_generations                          ) then
!
!    write(GP_print_unit,'(A)') &
!              'gpcf: i_GP_Indiv, GP_Child_Indiv_SSE(i_GP_Indiv) '
!    do  i_GP_Individual=1,n_GP_Individuals
!        write(GP_print_unit,'(6x,I6,6x,E15.7)') &
!               i_GP_Individual, GP_Child_Individual_SSE(i_GP_Individual)
!    enddo ! i_gp_individual
!
!endif ! i_GP_generation ...


!-------------------------------------------------------------------------------

! calculate a normalized ranking of the errors
! (higher individual SSE == lower value/ranking)

! calculate fitness as 1/individual sse

GP_Population_Ranked_Fitness = 0.0d0

do  i_GP_Individual=1,n_GP_Individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) cycle


    if( abs( GP_Child_Individual_SSE(i_GP_Individual) ) > 1.0D-30 )then
        GP_Population_Ranked_Fitness(i_GP_Individual) = &
             sse0  /  GP_Child_Individual_SSE(i_GP_Individual)
             !1.0d0 /  GP_Child_Individual_SSE(i_GP_Individual)
    else
        GP_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0
    endif ! abs( dff ) > 1.0D-30

enddo ! i_GP_Individual


!-------------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then

    write(GP_print_unit,'(/A)') &
          'gpcf: i_GP_indiv  GP_Child_Indiv_SSE  GP_Pop_Ranked_Fitness   '

    do  i_GP_Individual=1,n_GP_Individuals
        write(GP_print_unit,'(5x,I6,2x, 2(5x,E15.7))') &
              i_GP_individual, GP_Child_Individual_SSE(i_GP_Individual), &
             GP_Population_Ranked_Fitness(i_GP_Individual)
    enddo

endif ! i_GP_generation == 1 .or. ...

!-------------------------------------------------------------------------------

! calculate the sum of the rankings

GP_Integrated_Population_Ranked_Fitness = 0.0d0

dff=0.0
do  i_GP_Individual=1,n_GP_Individuals

    dff = dff + GP_Population_Ranked_Fitness(i_GP_Individual)

    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)=dff
enddo

!write(GP_print_unit,'(/A)') &
!    'gpcf: i_GP_indiv  GP_Child_Indiv_SSE  GP_Integ_Pop_Ranked_Fitness '
!do  i_GP_Individual=1,n_GP_Individuals
!    write(GP_print_unit,'(5x,I6,2(5x,E15.7))') &
!           i_GP_individual, &
!           GP_Child_Individual_SSE(i_GP_Individual), &
!           GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)
!enddo

!-------------------------------------------------------------------------------

! normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]


do  i_GP_Individual=1,n_GP_Individuals

    if( abs( GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals) ) &
                                                            > 1.0D-30 )then

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = &
        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) / &
                      GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)

    else

        GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) = 0.0D0

    endif ! abs(GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals))...

enddo ! i_GP_Individual

!-------------------------------------------------------------------------------

!if( i_GP_generation == 1                                 .or. &
!    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
!    i_GP_generation == n_GP_generations                          ) then
!
!    write(GP_print_unit,'(/A)') &
!      'gpcf: i_GP_indiv  GP_Child_Indiv_SSE  GP_Integ_Pop_Ranked_Fitness (norm)'
!
!    do  i_GP_Individual=1,n_GP_Individuals
!
!        write(GP_print_unit,'(5x,I6,2x,2(5x,E15.7))') &
!              i_GP_individual, &
!              GP_Child_Individual_SSE(i_GP_Individual), &
!              GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)
!    enddo
!
!endif ! i_GP_generation == 1 .or. ...


!-------------------------------------------------------------------------------

!write(GP_print_unit,*) i_GP_Generation,'MAIN',GP_Population_Ranked_Fitness(1)
!write(GP_print_unit,'(/A,17x,A,21x,I6,8x,E15.7)') &
!      'gpcf: i_GP_Gen,', 'GP_Pop_Ranked_Fit(1) ', &
!             i_GP_Generation, GP_Population_Ranked_Fitness(1)

!-------------------------------------------------------------------------------

! find GP_best_parent


i_GP_Best_Parent=1

dff=GP_Population_Ranked_Fitness(1)

do  i_GP_Individual=2,n_GP_individuals

    if(  GP_Individual_N_GP_param( i_GP_Individual ) < min_N_param ) cycle

    if( GP_Population_Ranked_Fitness(i_GP_individual) .gt. dff) then
        dff=GP_Population_Ranked_Fitness(i_GP_individual)
        i_GP_Best_Parent=i_GP_Individual
    endif

enddo ! i_GP_Individual

write(GP_print_unit,'(/A,2(1x,I6),2(1x,E15.7))') &
      'gpcf: i_GP_Gen,i_GP_Best_Parent,&
            &GP_Pop_Ranked_Fit(), GP_Child_Indiv_SSE()', &
             i_GP_Generation, i_GP_Best_Parent, &
             GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
             GP_Child_Individual_SSE(i_GP_Best_Parent)


!---------------------------------------------------------------------------

! fill output array of parameters for best individual
! and write on GP_print_unit


write(GP_print_unit,'(/A)') &
          'gpcf: i_CODE_eq  R_K_Init_Cond(i_CODE_eq)  &
          &GP_Pop_init_cond(i_CODE_eq,i_GP_Best_Parent)'
do  i_CODE_equation=1,n_CODE_equations

    write(GP_print_unit,'(6x,I6,7x, E20.10, 10x, E20.10)') &
          i_CODE_equation, &
          Runge_Kutta_Initial_Conditions(i_CODE_equation),  &
          GP_Population_Initial_Conditions( i_CODE_equation, i_GP_Best_Parent )

    output_array( i_CODE_equation ) = &
            GP_Population_Initial_Conditions( i_CODE_equation,i_GP_Best_Parent )


enddo ! i_CODE_equation


nop = n_CODE_equations

write(GP_print_unit,'(/A)') &
      'gpcf: count number of parameters, nop, in tree'

!write(GP_print_unit,'(A,2(1x,I6))') &
!      'gpcf: before tree loop n_code_equations, nop ', &
!                              n_code_equations, nop

write(GP_print_unit,'(/A)') &
     'gpcf: i_tree  i_node  nop  &
     &GP_pop_node_params(i_node,i_tree,i_GP_Best_Parent)'


tree_loop:&
do  i_tree=1,n_trees

    node_loop:&
    do  i_node=1,n_nodes

        if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0 )then

            nop = nop + 1

            write(GP_print_unit,'(2x,3(1x,I6), 1x, E20.10, 4x, E20.10)') &
                  i_tree, i_node, nop, &
                  GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent)

            output_array(nop) = &
                   GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent)


        endif ! GP_Adult_Pop_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0

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

write(GP_print_unit,'(/A,2(1x,I6))') &
      'gpcf: after tree loop n_code_equations, nop ', &
                             n_code_equations, nop

!---------------------------------------------------------------------------

!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        write(GP_print_unit,'(A,3(1x,I6))') &
!        'gpsbrl:1 i_tree, i_node, GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent)', &
!                  i_tree, i_node, GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent)
!    enddo ! i_node
!enddo ! i_tree



!---------------------------------------------------------------------------

if( L_GP_output_parameters )then

    write( GP_output_unit, '(I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
           i_GP_Generation, i_GP_best_parent, &
           GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
           nop, output_array(1:nop)

    !write(GP_print_unit, '(//A,1x,I6,1x,I6,1x,E15.7,1x,I6, 12(1x,E15.7))') &
    ! 'gpcf: i_GP_gen,i_GP_best_parent,GP_indiv_ranked_fit, nop, output_array', &
    !        i_GP_Generation, i_GP_best_parent, &
    !       GP_Population_Ranked_Fitness(i_GP_Best_Parent), &
    !       nop, output_array(1:nop)

endif ! L_GP_output_parameters

!---------------------------------------------------------------------------

if( i_GP_generation == 1                                 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
    i_GP_generation == n_GP_generations                          ) then


    write(GP_print_unit,'(/A)') &
         'gpcf: i_tree  i_node   GP_Adult_Pop_Node_Type(,,i_GP_Best_Parent)'

    tree_loop2:&
    do  i_tree=1,n_trees
        node_loop2:&
        do  i_node=1,n_nodes

            if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) &
                                                                > -9999  )then

                write(GP_print_unit,'(2x,3(1x,I6))') &
                  i_tree, i_node, &
                  GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent)

            endif !GP_Adult_Pop_Node_Type(i_Node,i_Tree,i_GP_Best_Parent)> -9999

        enddo node_loop2 ! i_node
    enddo tree_loop2 ! i_tree

    write(GP_print_unit,'(/A,1x,I6)') &
         'gpcf: print the tree for the best individual =', i_GP_Best_Parent

    call print_trees( i_GP_generation, i_GP_Best_Parent, i_GP_Best_Parent, &
                      GP_Adult_Population_Node_Type, 'best parent' )

endif ! i_GP_generation == 1 .or. ...


!---------------------------------------------------------------------------

GP_Adult_Individual_SSE  =  GP_Child_Individual_SSE

GP_Adult_Population_SSE  =  GP_Child_Individual_SSE
GP_Child_Population_SSE  =  GP_Child_Individual_SSE
!---------------------------------------------------------------------------
!
!if( i_GP_generation == 1                                 .or. &
!    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
!    i_GP_generation == n_GP_generations                          ) then
!
    write(GP_print_unit,'(/A)') &
          'gpcf: i_GP_Indiv   GP_Adult_Indiv_SSE    GP_Pop_Ranked_Fitness'

    do  i_GP_Individual=1,n_GP_individuals

        write(GP_print_unit,'(6x,I6,2(7x,E15.7))') &
              i_GP_Individual, &
              GP_Adult_Individual_SSE(i_GP_Individual), &
              GP_Population_Ranked_Fitness(i_GP_Individual)

    enddo ! i_GP_individual
!
!endif ! i_GP_generation == 1 .or. ...
!
!!off if( i_GP_Generation .eq. 3) Stop

!-------------------------------------------------------------------------------

call calc_stats( n_GP_individuals, GP_Population_Ranked_Fitness,  &
                 mean_fit, rms_fit, std_dev_fit )

write(GP_print_unit,'(/A,1x,I6,3(1x,E15.7)/)') &
   'gpcf: GP_Gen, GP_Pop_Rank_Fit mean, rms, std_dev', &
          i_GP_Generation, mean_fit, rms_fit, std_dev_fit

!-------------------------------------------------------------------------------

if( L_GP_log )then

    ! write information to a GP log file giving:
    ! generation, individual, SSE, individual_fitness

    do  i_GP_Individual=1,n_GP_individuals

        write(GP_log_unit) &
              i_GP_generation, &
              i_GP_Individual, &
              GP_Adult_Individual_SSE(i_GP_Individual), &
              GP_Population_Ranked_Fitness(i_GP_Individual)

    enddo ! i_GP_individual

endif ! L_GP_log

!-------------------------------------------------------------------------------

if( L_unit50_output )then

    ! calculate array for writing on unit50.txt ( unit_gp_out )

    do i_GP_Individual=1,n_GP_individuals

       GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
            GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual)

       !off  GP_Node_Type_for_Plotting(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
       !off  GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)

    enddo ! i_GP_individual



    write(unit_gp_out) GP_Node_Type_for_Plotting

endif ! L_unit50_output



!------------------------------------------------------------------------------



! don't call GP_ranking_sort on last generation since after the generation loop
! GP_select_best...      is called and uses the arrays and the i_GP_best parent

! GP_ranking re-orders all these arrays, so that the best parent is no longer at
! the index it was in GP_calc_fitness

if( i_GP_generation < n_GP_generations )then

    ! re-sort based on rankings

    ! uses:
    !  GP_Child_Population_SSE
    !  GP_population_node_parameters
    !  GP_Adult_Population_Parameter_Solution
    !  GP_Child_Population_Parameter_Solution
    !  GP_Child_Population_Node_Type

    ! sets:
    !  GP_Child_Population_SSE
    !  GP_Adult_Population_Parameter_Solution
    !  GP_Child_Population_Parameter_Solution
    !  GP_Adult_Population_Node_Type
    !  GP_Child_Population_Node_Type
    !  GP_Population_Initial_Conditions
    !  GP_Adult_Population_SSE
    !  GP_Adult_Individual_SSE
    !  GP_Child_Individual_SSE
    !  GP_population_node_parameters
    !  GP_Population_Ranked_Fitness
    !  GP_Integrated_Population_Ranked_Fitness


    write(GP_print_unit,'(/A,1x,I6/)') &
       'gpcf: call GP_ranking GP_Gen ',  i_GP_Generation

    call GP_ranking_sort()


endif ! i_GP_generation < n_GP_generations



!-------------------------------------------------------------------------------


return


end subroutine GP_calc_fitness
