program GPCODE_GA_lmdif_parameter_optimization_test

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
integer :: i_diversity 
integer :: nop
integer :: message_len

logical :: Lplot


real(kind=8), allocatable, dimension(:) :: answer
real(kind=8), allocatable, dimension(:) :: output_array

real (kind=8) :: ssum, ssum2, totobs, dff

real (kind=8) :: xcount

!----------------------------------------------------------------------------------------



!---------------------------------------
! current setup
! lmdif runs only on best individual of each generation
! no replacement for bad points in GP*n
! no retry in setup_run_fcn or setup_run_lmdif
! in random_real, boundary is 0.1
! --  0.1 for random range [0,1], 0.9 for random range [0,50]
! max interations in lmdif is 100
!---------------------------------------


call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)




Lprint_lmdif = .TRUE.

!------------------------------------------------------------------

if( myid == 0 )then

    ! read the control input from file  "GPCODE_cntl"
    ! and broadcast the values read to the other cpu's

    call read_cntl_stuff( )

    open( unit_gp_out, file = 'unit50.txt', &
          form = 'unformatted', access='sequential', &
          status = 'unknown' )

endif !   myid == 0

!------------------------------------------------------------------



! broadcast the values read in by cpu 0 to others

call MPI_BCAST( GA_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_save_elites_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_GA_Generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_GA_Individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_time_steps, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( dt, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( model, 30,     &
                MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( GP_Tree_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( GP_Elitist_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( GP_Asexual_Reproduction_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( GP_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( GP_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_gp_individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_gp_generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
!------------------------------------------------------------------

! wait until everybody has the values

call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------


! set the scalar values for the model


call init_values( 0 )



!------------------------------------------------------------------

if( myid == 0 )then

    write(GP_print_unit,'(A,1x,I10)')    '0: n_GA_individuals           ', &
                                             n_GA_individuals
    write(GP_print_unit,'(A,1x,I10)')    '0: n_time_steps               ', &
                                             n_time_steps
    write(GP_print_unit,'(A,1x,I10)')    '0: n_GA_Generations           ', &
                                             n_GA_Generations
    write(GP_print_unit,'(A,1x, E15.7)') '0: GA_Crossover_Probability   ', &
                                             GA_Crossover_Probability
    write(GP_print_unit,'(A,1x, E15.7)') '0: GA_Mutation_Probability    ', &
                                             GA_Mutation_Probability
    write(GP_print_unit,'(A,1x, E15.7)') '0: GA_save_elites_Probability ', &
                                             GA_save_elites_Probability

    write(GP_print_unit,'(/A,1x, E15.7)')'0: GP_Tree_Probability        ', &
                                             GP_Tree_Probability

    write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Elitist_Probability     ', &
                                             GP_Elitist_Probability

    write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Crossover_Probability   ', &
                                             GP_Crossover_Probability

    write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Asexual_Reproduction_Probability ', &
                                             GP_Asexual_Reproduction_Probability

    write(GP_print_unit,'(A,1x, E15.7)') '0: GP_Mutation_Probability    ', &
                                             GP_Mutation_Probability
    write(GP_print_unit,'(/A,1x,I10)')   '0: n_gp_individuals           ', &
                                             n_gp_individuals
    write(GP_print_unit,'(A,1x,I10)')    '0: n_gp_generations           ', &
                                             n_gp_generations

endif ! myid == 0

!------------------------------------------------------------------

! allocate variable dimension arrays

call allocate_arrays1( )

allocate( answer( n_maximum_number_parameters ) )
allocate( output_array( n_maximum_number_parameters ) )


!------------------------------------------------------------------


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the twin experiment 'nature'
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]



GP_Population_Node_Type=-9999

! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]

GP_Individual_Node_Parameters=0.                        ! Matrix Operation
GP_Individual_Node_Type=-9999                           ! Matrix Operation
GP_Population_Node_Parameters=0.                        ! Matrix Operation
GP_Adult_Population_Node_Type=-9999                     ! Matrix Operation
GP_Child_Population_Node_Type=-9999                     ! Matrix Operation

!------------------------------------------------------------------

! fill the model arrays

! sets Runge_Kutta_Initial_Conditions,
!      GP_Individual_Node_Type and
!      GP_Individual_Node_Parameters

call init_values( 1 )


!------------------------------------------------------------------
! set the desired 'twin experiment' population node type
! and parameter using the info from the set up file

! run the Runge-Kutta model only once with proc 0

call set_answer_arrays( )


!------------------------------------------------------------------------

! then broadcast the R-K result to all processors

message_len = ( n_time_steps + 1 ) * n_CODE_equations

if( myid == 0 )then
    write(GP_print_unit,'(/A/)') &
          '0: time_step   Runge_Kutta_Solution(time_step,1:n_CODE_equations)'
    do  i = 0, n_time_steps
        write(GP_print_unit,'(I6,8(1x,E15.7))') &
              i, Runge_Kutta_Solution(i,1:n_CODE_equations)
    enddo ! i
    !write(GP_print_unit,'(A,1x,I10)') &
    !      '0: Runge_Kutta message length message_len = ', message_len
endif ! myid == 0


call MPI_BCAST( Runge_Kutta_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



Data_Array=Runge_Kutta_Solution          ! Matrix Operation


!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

call comp_data_variance( )


!--------------------------------------------------------------------------------

! put the desired model parameters in the answer array


answer = 0.0d0 ! set all to zero

n_parameters = 0

do  i_CODE_equation=1,n_CODE_equations
    n_parameters=n_parameters+1
    answer(n_parameters)=Runge_Kutta_Initial_Conditions(i_CODE_equation)
enddo ! i_CODE_equation



! calculate how many parameters total to fit for the specific individual CODE

do  i_tree=1,n_trees
    do  i_node=1,n_nodes

        if( GP_individual_node_type(i_node,i_tree) .eq. 0) then  ! there is a set parameter
            n_parameters=n_parameters+1
            answer(n_parameters)=GP_Individual_Node_Parameters(i_node,i_tree)
        endif ! GP_individual_node_type(i_node,i_tree) .eq. 0

    enddo ! i_node
enddo ! i_tree


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if( myid == 0 )then

    write(GP_print_unit,'(/A,1x,I6)') 'Total Parameters for this run = ',n_parameters
    write(GP_print_unit,'(A,1x,I6)')  '0: n_GA_individuals ', n_GA_individuals

    !write(GP_print_unit,'(A,1x,I6)') 'number of initial conditions   = ',n_initial_conditions
    write(GP_print_unit,'(A,1x,I6)')   '0: n_time_steps     ', n_time_steps
    write(GP_print_unit,'(A,1x,I6)')   '0: n_GA_Generations ', n_GA_Generations
    write(GP_print_unit,'(A,1x,I6)')   '0: n_levels         ', n_levels
    write(GP_print_unit,'(A,1x,I6)')   '0: n_functions      ', n_functions
    write(GP_print_unit,'(A,1x,I6)')   '0: n_CODE_equations ', n_CODE_equations
    write(GP_print_unit,'(A,1x,I6)')   '0: n_parameters     ', n_parameters
    write(GP_print_unit,'(A,1x,I6)')   '0: n_trees          ', n_trees
    write(GP_print_unit,'(A,1x,I6)')   '0: n_nodes          ', n_nodes
    write(GP_print_unit,'(A,1x,E15.7)')'0: dt (days)        ', dt
    write(GP_print_unit,'(A,1x,E15.7)')'0: dt (min)         ', dt * 1440.d0
    !write(GP_print_unit,*) individual_fitness
    !write(GP_print_unit,*) data_array
    write(GP_print_unit,'(A,1x,E15.7)') '0: GA_Crossover_Probability    ', &
                                            GA_Crossover_Probability
    write(GP_print_unit,'(A,1x,E15.7)') '0: GA_Mutation_Probability     ', &
                                            GA_Mutation_Probability

    write(GP_print_unit,'(A,1x,E15.7)') '0: GA_save_elites_Probability  ', &
                                            GA_save_elites_Probability

    !write(GP_print_unit,'(A)') '0: code calls lmdif for all individuals and generations '
    write(GP_print_unit,'(A)') '0: code calls lmdif only for best individual on last generation'

    !write(GP_print_unit,'(A)') '0: regular random_real initialization of child parameters '
    write(GP_print_unit,'(A)') '0: 2-range random_real initialization of child parameters '

    write(GP_print_unit,'(A//)') '0: replacement of bad individuals used'
    !write(GP_print_unit,'(A//)') '0: NO replacement of bad individuals '


    ! calculate the generation interval for printing the list of children

    GA_child_print_interval = n_GA_generations /  number_GA_child_prints

    GP_child_print_interval = n_GP_generations /  number_GP_child_prints

    write(GP_print_unit,'(A,1x,I6)')   '0: GP_child_print_interval ', &
                                           GP_child_print_interval
    write(GP_print_unit,'(A,1x,I6)')   '0: GA_child_print_interval ', &
                                           GA_child_print_interval

    !-----------------------------------------------------------------------------

    ! this call calculates the sse value for RK solution = 0 for all time steps
    ! note:  sse0 is only used by cpu 0 which does all fitness calculations

    call sse0_calc( )

    !-----------------------------------------------------------------------------

    write(GP_print_unit,'(/A)') ' '
    do  i = 1, n_parameters
        write(GP_print_unit,'(A,1x,I6,2x,E24.16)') &
              '0: i, answer(i) ', i, answer(i)
    enddo ! i
    write(GP_print_unit,'(/A)') ' '
    !-----------------------------------------------------------------------------


    open( GA_output_unit, file='output_parameters', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

    open( GP_output_unit, file='GP_output_parameters', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

    open( GP_summary_output_unit, file='GP_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )


    write(GP_summary_output_unit, '(I10)') n_code_equations
     
endif ! myid == 0


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the Genetic Programming optimization routine for the Binary Tree Evolution
! with the embedded GA_lmdif parameter optimization scheme
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! Number of Carry-Over Elitists

n_GP_Elitists = nint(GP_Elitist_Probability*n_GP_individuals)

! Number of GP Fitness Proportionate Reproduction

n_GP_Asexual_Reproductions = nint(GP_Asexual_Reproduction_Probability*n_GP_individuals)

! Number of GP Sexual Crossovers

n_GP_Crossovers = nint(GP_Crossover_Probability*n_GP_individuals)

! Number of GP Mutations

n_GP_Mutations = n_GP_Individuals - &
                 ( n_GP_Elitists + n_GP_Crossovers + n_GP_Asexual_Reproductions )

!n_GP_Mutations = &
!min( nint( GP_Mutation_Probability * n_GP_individuals ) , n_GP_Mutations )
if( myid == 0 )then
    write(GP_print_unit,'(A,1x,i6,1x,I6)') &
      '0: n_GP_Mut, min(nint(GP_Mut_Prob*n_GP_indiv),n_GP_Mut) ', &
          n_GP_Mutations, &
          min(nint(GP_Mutation_Probability*n_GP_individuals),n_GP_Mutations)
endif ! myid == 0


if( myid == 0 )then

    write(GP_print_unit,'(/A,1x,I6)')   '0: n_gp_individuals           ', &
                                            n_gp_individuals
    write(GP_print_unit,'(A,1x,I6/)')   '0: n_gp_generations           ', &
                                            n_gp_generations

    write(GP_print_unit,'(A,1x,F10.6)') '0: GP_Elitist_Probability     ', &
                                            GP_Elitist_Probability
    write(GP_print_unit,'(A,1x,F10.6)') '0: GP_Crossover_Probability   ', &
                                            GP_Crossover_Probability
    write(GP_print_unit,'(A,1x,F10.6)') '0: GP_Mutation_Probability    ', &
                                            GP_Mutation_Probability
    write(GP_print_unit,'(A,1x,F10.6)') '0: GP_Asexual_Reproduction_Probability ', &
                                            GP_Asexual_Reproduction_Probability

    write(GP_print_unit,'(/A,1x,I6)')   '0: n_GP_Elitists              ', &
                                            n_GP_Elitists
    write(GP_print_unit,'(A,1x,I6)')    '0: n_GP_Crossovers            ', &
                                            n_GP_Crossovers
    write(GP_print_unit,'(A,1x,I6)')    '0: n_GP_Mutations             ', &
                                            n_GP_Mutations
    write(GP_print_unit,'(A,1x,I6/)')   '0: n_GP_Asexual_Reproductions ',  &
                                            n_GP_Asexual_Reproductions

endif ! myid == 0


!  make sure numbers add up to total number of individuals

if( n_GP_Elitists              + &
    n_GP_Asexual_Reproductions + &
    n_GP_Crossovers            + &
    n_GP_Mutations                 .gt. n_GP_Individuals) then

    write(GP_print_unit,'(/A/)') &
          '0:Sum of n_GP_Elitists + n_Asexual_Reproduction + &
          &n_GP_Crossovers + n_GP_Mutations is too high'

    call MPI_FINALIZE(ierr)
    stop '0:sum too big'

elseif( n_GP_Elitists              + &
        n_GP_Asexual_Reproductions + &
        n_GP_Crossovers            + &
        n_GP_Mutations                 .lt. n_GP_Individuals) then

    write(GP_print_unit,'(/A/)') &
          '0: Sum of n_GP_Elitists + n_Asexual_Reproduction + &
          &n_GP_Crossovers + n_GP_Mutations is too low'

    call MPI_FINALIZE(ierr)
    stop '0:sum too small'

endif !   n_GP_Elitists + ...



!------------------------------------------------------------------------------



do  i_GP_Generation=1,n_GP_Generations


    if( myid == 0 )then
        write(GP_print_unit,'(/A,1x,I6,1x,A,1x,I6/)') &
              '0: GP Generation # ',i_GP_Generation,&
                ' is underway.   n_Nodes * n_Trees = ', n_Nodes*n_Trees
    endif ! myid == 0


    GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type  ! Matrix Operation


    ! Run_GP_Calculate_Fitness determines if the new GP child
    ! has to be sent to GA_lmdif for parameter optimization

    Run_GP_Calculate_Fitness=.false.


    ! randomly create the initial tree arrays for each individual and
    ! send them all to GA_lmdif for parameter optimization

    if( i_GP_Generation .eq. 1) then

        ! determines if the new GP child
        !has to be sent to GA_lmdif for parameter optimization

        Run_GP_Calculate_Fitness=.true.

        if( myid == 0 )then

            write(GP_print_unit,'(A,1x,I6/)') &
              '0: call GP_Tree_Build   Generation = ',i_GP_Generation

            ! initialize the GP_Adult_Population_Node_Type array with random trees

            call GP_Tree_Build

            !debug only do  i_GP_individual = 1, n_GP_Individuals
            !debug only     GP_Adult_Population_Node_Type(i_GP_individual,:,:) = &
            !debug only                       GP_Node_Type_Answer(:,:) ! debug only
            !debug only enddo


        endif ! myid == 0

        ! send values for GP_Adult_Population_Node_Type to all processors

        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6)') &
        !      '0: broadcast  GP_Adult_Population_Node_Type Generation = ',i_GP_Generation
        !endif ! myid == 0

        message_len = n_GP_Individuals * n_Nodes * n_Trees

        call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        ! wait until everybody has the values

        !call MPI_BARRIER( MPI_COMM_WORLD, ierr )

        if( myid == 0 )then
            call GP_calc_diversity_index( n_GP_individuals,  &
                                          GP_Adult_Population_Node_Type, i_diversity )
        endif ! myid == 0 

        !-----------------------------------------------------------------------------


    else !  i_GP_Generation > 1

        if( myid == 0 )then

            ! create the next 'generation' of tree structures using either:
            !    i)  Copy the top n_GP_Elitists individuals into the next generation
            !   ii)  GP Fitness-Proportionate Asexual Reproduction;
            !  iii)  GP Tournament-Style Sexual Reproduction, and;
            !   iv)  GP Mutation


            !write(GP_print_unit,'(/A/)')  &
            !      '0: ############################################################################'
            !write(GP_print_unit,'(A/)')  '0: GP_Adult trees before call selection routines '

            !call print_trees( 1, n_GP_individuals, GP_Adult_Population_Node_Type )

            !write(GP_print_unit,'(/A/)') &
            !      '0: ############################################################################'


            !   i) move the top n_GP_Elitists into the next generation

            if( n_GP_Elitists .gt. 0 ) then
                write(GP_print_unit,'(A)')'0:  call GP_Elitists '
                call GP_Elitists
            endif !  n_GP_Elitists .gt. 0


            !write(GP_print_unit,'(/A/)')  &
            !      '0: ############################################################################'
            !write(GP_print_unit,'(A/)')  '0: GP_Child trees after  call to GP_Elitists'
            !call print_trees( 1, n_GP_individuals, GP_Child_Population_Node_Type )
            !write(GP_print_unit,'(/A/)') &
            !      '0: ############################################################################'


            !   ii) Carry out "GP Fitness-Proportionate Reproduction"

            if( n_GP_Asexual_Reproductions .gt. 0 )then
                write(GP_print_unit,'(A)') &
                      '0:  call GP_Fitness_Proportionate_Asexual_Reproduction '
                call GP_Fitness_Proportionate_Asexual_Reproduction
            endif !  n_GP_Asexual_Reproductions .gt. 0

            !write(GP_print_unit,'(/A/)')  &
            !      '0: ############################################################################'
            !write(GP_print_unit,'(A/)')  '0: GP_Child trees after call to &
            !                  &GP_Fitness_Proportionate_Asexual_Reproduction'
            !call print_trees( 1, n_GP_individuals, GP_Child_Population_Node_Type )
            !write(GP_print_unit,'(/A/)') &
            !      '0: ############################################################################'



            !  iii) Carry out "GP Tree Crossover" Operations
            !       Using Tournament-Style Sexual Reproduction Selection
            !       and randomly use it to replace the new children

            if( n_GP_Crossovers .gt. 0 )then
                write(GP_print_unit,'(A)') &
                      '0:  call GP_Tournament_Style_Sexual_Reproduction '
                call GP_Tournament_Style_Sexual_Reproduction
            endif !  n_GP_Crossovers .gt. 0



            !write(GP_print_unit,'(/A/)')  &
            !      '0: ############################################################################'
            !write(GP_print_unit,'(A/)')  '0: GP_Child trees after call to &
            !                  &GP_Tournament_Style_Sexual_Reproduction'      
            !call print_trees( 1, n_GP_individuals, GP_Child_Population_Node_Type )
            !write(GP_print_unit,'(/A/)') &
            !      '0: ############################################################################'


            !   iv) Carry out "GP Parameter Mutation" Operations

            if( n_GP_Mutations .gt. 0 )then
                write(GP_print_unit,'(A)')'0:  call GP_Mutations '
                call GP_Mutations
            endif !  n_GP_Mutations .gt. 0



            !write(GP_print_unit,'(/A/)')  &
            !      '0: ############################################################################'
            !write(GP_print_unit,'(A/)')  '0: GP_Child trees after call to GP_Mutations'
            !call print_trees( 1, n_GP_individuals, GP_Child_Population_Node_Type )
            !write(GP_print_unit,'(/A/)') &
            !      '0: ############################################################################'

            call GP_calc_diversity_index( n_GP_individuals, &
                                          GP_Child_Population_Node_Type, i_diversity )

            
        endif ! myid == 0

        !------------------------------------------------------------------------------------
        ! broadcast


        ! GP_Child_Population_Node_Type


        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6)') &
        !      '0: broadcast  GP_Child_Pop_Node_Type  Generation = ',&
        !                                        i_GP_Generation
        !endif ! myid == 0

        message_len = n_GP_Individuals * n_Nodes * n_Trees
        call MPI_BCAST( GP_Child_Population_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        ! GP_Adult_Population_Node_Type

        message_len = n_GP_Individuals * n_Nodes * n_Trees
        call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


        ! Parent_Tree_Swap_Node_Type

        message_len = 2 * n_Nodes
        call MPI_BCAST( Parent_Tree_Swap_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        ! GP_Child_Individual_SSE

        call MPI_BCAST( GP_Child_Individual_SSE, n_GP_individuals,    &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

        ! GP_Individual_Ranked_Fitness

        call MPI_BCAST( GP_Individual_Ranked_Fitness, n_GP_individuals,    &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

        ! GP_Integrated_Ranked_Fitness

        call MPI_BCAST( GP_Integrated_Ranked_Fitness, n_GP_individuals,    &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

        ! Run_GP_Calculate_Fitness array

        call MPI_BCAST( Run_GP_Calculate_Fitness , n_GP_Individuals,    &
                        MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

        !------------------------------------------------------------------------------------

        ! wait until everybody has the values

        call MPI_BARRIER( MPI_COMM_WORLD, ierr )

        !------------------------------------------------------------------------------------

        !   move the children into adulthood

        if( myid == 0 )then
            write(GP_print_unit,'(A,1x,I6)') &
              '0: move GP_child_pop_node_type to GP_adult_pop_node_type    Generation = ',&
                                                                      i_GP_Generation
        endif ! myid == 0

        GP_Adult_Population_Node_Type=GP_Child_Population_Node_Type      ! Matrix Operation

    endif ! i_GP_Generation .eq. 1


    !-----------------------------------------------------------------------------------------

    ! print trees before tree clean

    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A/)')  &
    !          '0: ############################################################################'
    !    write(GP_print_unit,'(A/)')  '0: trees before call to GP_Clean_Tree_Nodes'
    !    call print_trees( 1, n_GP_individuals, GP_Adult_Population_Node_Type )
    !    write(GP_print_unit,'(/A/)') &
    !          '0: ############################################################################'
    !endif !  myid == 0


    !-----------------------------------------------------------------------------------------

    ! sweep through all the GP_Adult_Population_Node_Type
    ! to replace function nodes that have both terminals set as parameters
    ! and set the node to a parameter itself

    if( myid == 0 )then
        write(GP_print_unit,'(/A,1x,I6/)') &
              '0: call GP_Clean_Tree_Nodes   Generation = ', &
                                        i_GP_Generation
    endif ! myid == 0


    call GP_Clean_Tree_Nodes


    !-----------------------------------------------------------------------------------------

    ! print trees after call to GP_Clean_Tree_Nodes

    if( myid == 0 )then
        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. & 
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)')  &
            '0: ############################################################################'
            write(GP_print_unit,'(A/)')  '0: trees after call to GP_Clean_Tree_Nodes'
    
            call print_trees( 1, n_GP_individuals, GP_Adult_Population_Node_Type )
    
            write(GP_print_unit,'(/A/)')  &
            '0: ############################################################################'
            endif ! i_GP_generation == 1 .or. ...

    endif !  myid == 0

    !-----------------------------------------------------------------------------------------


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! GA_lmdif subroutine segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    gp_ind_loop:&
    do  i_GP_individual=1,n_GP_individuals

        if( myid == 0 )then
            write(GP_print_unit,'(/A)')&
                  '0:-----------------------------------------------------------------'
            write(GP_print_unit,'(A,1x,I6,4x,L1   )') &
                  '0: i_GP_individual, Run_GP_Calculate_Fitness ', &
                      i_GP_individual, Run_GP_Calculate_Fitness(i_GP_Individual)
        endif !  myid == 0

        GP_Individual_Node_Parameters = 0.0 ! these get set randomly in the GA-lmdif search algorithm

!>>>>>>>>>> jjm 20130417
!        if( myid == 0 )then
!            write(GP_print_unit,'(/A,3(1x,I6)/)') &
!                  '0: i_GP_individual, n_nodes, n_trees ', &
!                      i_GP_individual, n_nodes, n_trees
!            do  i_Tree=1,n_Trees
!                do  i_Node=1,n_Nodes
!                    if( GP_Individual_Node_Type(i_Node,i_Tree) /= -9999 )then
!                        write(GP_print_unit,'(A,3(1x,I6))') &
!                          '0: i_node, i_tree, GP_Individual_Node_Type(i_Node,i_Tree)', &
!                              i_node, i_tree, GP_Individual_Node_Type(i_Node,i_Tree)
!                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) /= -9999
!                enddo ! i_node
!            enddo ! i_tree
!            write(GP_print_unit,'(/A)') ' '
!            do  i_Tree=1,n_Trees
!                do  i_Node=1,n_Nodes
!                    if( GP_Individual_Node_Type(i_Node,i_Tree) /= -9999 )then
!                        write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
!                          '0: i_node, i_tree, GP_Individual_Node_Parameters(i_Node,i_Tree)', &
!                              i_node, i_tree, GP_Individual_Node_Parameters(i_Node,i_Tree)
!                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) /= -9999
!                enddo ! i_node
!            enddo ! i_tree
!            write(GP_print_unit,'(/A)') ' '
!        endif !  myid == 0
!<<<<<<<<<<< jjm 20130417

        if( Run_GP_Calculate_Fitness(i_GP_Individual) ) then

            GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) = &
              GP_Adult_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

            ! calculate how many parameters total to fit for the specific individual CODE

            n_GP_Parameters=0
            do  i_Tree=1,n_Trees

                do  i_Node=1,n_Nodes

                    if( GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0) then  ! a set parameter
                        n_GP_Parameters=n_GP_Parameters+1
                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

                    !if( myid == 0 )then
                    !    write(GP_print_unit,'(A,4(1x,I6))')&
                    !   '0: i_GP_individual, i_node, i_tree, n_GP_parameters ', &
                    !       i_GP_individual, i_node, i_tree, n_GP_parameters
                    !endif !  myid == 0

                enddo ! i_node

            enddo ! i_tree

            if( myid == 0 )then
                write(GP_print_unit,'(/A,3(1x,I6)/)') &
                      '0: i_GP_individual, n_nodes, n_trees ', &
                          i_GP_individual, n_nodes, n_trees
                write(GP_print_unit,'(A,1x,I6,3x,A,1x,I6/)')&
                      '0: for i_GP_Individual', i_GP_Individual, &
                      'the number of parameters is:  n_GP_parameters =', n_GP_parameters
            endif !  myid == 0

            GP_Individual_N_GP_param(i_GP_individual) = n_GP_parameters


            ! cycle the i_GP_individual loop if there are no GP parameters


            if( n_GP_parameters == 0 ) then
                if( myid == 0 )then
                    write(GP_print_unit,'(A,1x,I6)')&
                          '0: skipping this i_GP_Individual --&
                          &  the number of parameters is ', n_GP_parameters
                endif !  myid == 0
                cycle gp_ind_loop
            endif ! n_GP_parameters == 0

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND
            ! LINK THE SSE OUTPUT TO THE ARRAY AT THE END
            ! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEED TO BE PULLED OUT
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !--------------------------------------------------------------------------------

            if( myid == 0 )then
                write(GP_print_unit,'(/A,1x,I6)') &
                 '0: call GPCODE_GA_lmdif_Parameter_Optimization routine', myid
                write(GP_print_unit,'(A,2(1x,I6)/)') &
                 '0: i_GP_Generation, i_GP_individual ', & 
                     i_GP_Generation, i_GP_individual 
                write(GA_print_unit,'(/A,1x,I6)') &
                 '0: call GPCODE_GA_lmdif_Parameter_Optimization routine', myid
                write(GA_print_unit,'(A,2(1x,I6)/)') &
                 '0: i_GP_Generation, i_GP_individual ', & 
                     i_GP_Generation, i_GP_individual 
            endif ! myid == 0

            !-------------------------------------------------

            call GPCODE_GA_lmdif_Parameter_Optimization()

            ! wait until everybody has the values

            call MPI_BARRIER( MPI_COMM_WORLD, ierr )

            !-------------------------------------------------

            if( myid == 0 )then
                write(GP_print_unit,'(//A,1x,I6//)') &
                 '0: aft call GPCODE_GA_lmdif_Parameter_Optimization routine'
            endif ! myid == 0


            GP_population_fitness(i_GP_individual) = individual_fitness


            !if( myid == 0 )then
            !    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7/)') &
            !          '0:1 i_GP_gen, i_GP_Indiv, GP_pop_fit(i_GP_indiv)', &
            !               i_GP_generation, i_GP_Individual, &
            !               GP_population_fitness(i_GP_individual)
            !endif ! myid == 0

            !--------------------------------------------------------------------------------

            if( myid == 0 )then


                write(GP_print_unit,'(/A,2(1x,I6), 1x, E15.7)')&
                      '0:1 i_GP_gen, i_GP_indiv, GP_pop_fit(i_GP_indiv) ', &
                           i_GP_generation, i_GP_individual, &
                           GP_population_fitness(i_GP_individual)

                !----------------------------------------------------------------------------

                ! print side-by-side comparisons of starting values
                ! and values from optimization

                write(GP_print_unit,'(/A/)') &
                      '0:           truth value           output value '

                do  i_CODE_equation=1,n_CODE_equations

                    write(GP_print_unit,'(I6,1x, E20.10, 4x, E20.10)') &
                          i_CODE_equation, &
                          Runge_Kutta_Initial_Conditions(i_CODE_equation),  &
                          GP_individual_initial_conditions(i_CODE_equation)

                          output_array( i_CODE_equation ) = &
                                   GP_individual_initial_conditions(i_CODE_equation)

                    !write(GA_output_unit,'(E24.16)') &
                    !      GP_individual_initial_conditions(i_CODE_equation)

                enddo ! i_CODE_equation

            endif !  myid == 0

            !--------------------------------------------------------------------------------

            ! compute GP_Child_Individual_SSE(i_GP_Individual)

            call comp_GP_child_indiv_sse()

            !-------------------------------------------------------------------------------------


            ! set the GA_lmdif-optimized initial condition array

            GP_Population_Initial_Conditions(i_GP_Individual,1:n_CODE_Equations) = &
                GP_Individual_Initial_Conditions(1:n_CODE_Equations) ! Matrix Operation

            if( myid == 0 )then
                write(GP_print_unit,'(/A,2(1x,I6), 2(1x, E15.7))')&
                      '0:3 i_GP_gen, i_GP_indiv, &
                      &GP_Pop_Init_Cond(i_GP_Indiv,1:n_CODE_Eq) ', &
                          i_GP_generation, i_GP_individual, &
                       GP_Population_Initial_Conditions(i_GP_Individual,1:n_CODE_Equations)
            endif !  myid == 0


            !-------------------------------------------------------------------------------------

            ! set the GA_lmdif-optimized CODE parameter set array

            GP_Population_Node_Parameters(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
                GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) ! Matrix Operation


            if( myid == 0 )then

                !----------------------------------------------------------------------------

                ! print the node parameters (if there are any)

                !write(GP_print_unit,'(/A)')  &
                !   '0: GP_population_node_parameters(i_GP_individual,:,:) '
                !write(GP_print_unit,'(5(1x,E15.7))')  &
                !       GP_population_node_parameters(i_GP_individual,:,:)

                if( any( abs( GP_population_node_parameters(i_GP_individual,:,:) ) &
                                                                            > 1.0d-20 ) )then

                    write(GP_print_unit,'(/A/)') &
                       '0:  node  tree  Runge_Kutta_Node_Parameters&
                       &   GP_population_node_parameters'

                    do  i_tree=1,n_trees
                        do  i_node=1,n_nodes
                            if( abs( GP_population_node_parameters( &
                                      i_GP_individual,i_node,i_tree) ) > 1.0d-20   )then

                                write(GP_print_unit,'(2(1x,I6), 1x, E20.10, 4x, E20.10)') &
                                      i_node, i_tree, &
                                      Runge_Kutta_Node_Parameters(i_node,i_tree), &
                                      GP_population_node_parameters(i_GP_individual,i_node,i_tree)
                            endif
                        enddo ! i_node
                    enddo  ! i_tree

                endif ! any( abs( GP_population_node_parameters(i_GP_individual,:,:) )> 1.0d-20 )

                !------------------------------------------------------------------------------------

                write(GP_print_unit,'(/A)') &
                      '0:---------------------------------------------------------------'
                write(GP_print_unit,'(/A/)')  &
                      '0: after loading GP_Pop arrays with GP_indiv array values '

                write(GP_print_unit,'(A,1x,I6,4x,L1)') &
                      '0: i_GP_individual, Run_GP_Calculate_Fitness(i_GP_Individual) ', &
                          i_GP_individual, Run_GP_Calculate_Fitness(i_GP_Individual)


                ! this prints a summary of the initial conditions, parameters and node types
                ! for this individual after being optimized in GPCODE*opt

                call summary_GP_indiv( i_GP_individual )

            endif !  myid == 0

        endif !   Run_GP_Calculate_Fitness(i_GP_Individual)

    enddo  gp_ind_loop    !   i_GP_individual


    !-------------------------------------------------------------------------------------

    ! do fitness calculations for this GP generation

    if( myid == 0 )then

        write(GP_print_unit,'(/A)')&
              '0:-----------------------------------------------------------------'
        write(GP_print_unit,'(/A/)') '0: call GP_fitness_reset '

        call GP_fitness_reset

        write(GP_print_unit,'(/A/)') '0: aft call GP_fitness_reset '
        write(GP_print_unit,'(A/)')&
              '0:-----------------------------------------------------------------'

    endif ! myid == 0

    !-------------------------------------------------------------------------------------



    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
    !          '0: broadcast GP_Adult_Individual_SSE = ',&
    !                        GP_Adult_Individual_SSE
    !endif ! myid == 0

    call MPI_BCAST( GP_Adult_Individual_SSE, n_GP_individuals,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
    !          '0: broadcast GP_Individual_Ranked_Fitness ',&
    !                        GP_Individual_Ranked_Fitness
    !endif ! myid == 0

    call MPI_BCAST( GP_Individual_Ranked_Fitness, n_GP_individuals,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A/(5(1x,E15.7)))') &
    !          '0: broadcast GP_Integrated_Ranked_Fitness ',&
    !                        GP_Integrated_Ranked_Fitness
    !endif ! myid == 0

    call MPI_BCAST( GP_Integrated_Ranked_Fitness, n_GP_individuals,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )




enddo !  i_GP_Generation




if( myid == 0 )then
    write(GP_print_unit,'(/A/)') '0: after i_GP_generation loop  '
endif ! myid == 0

!------------------------------------------------------------------

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! plot results

!if( myid == 0 )then
!    Lplot = .true.
!    if( Lplot ) call plot_results(Runge_Kutta_Solution,n_time_steps, n_CODE_equations )
!endif ! myid == 0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!------------------------------------------------------------------

! deallocate variable dimension arrays

call deallocate_arrays1( )

deallocate( answer )
deallocate( output_array )


!------------------------------------------------------------------
close( GA_output_unit )
close( GP_output_unit )
close( GP_summary_output_unit )
close( unit_gp_out )
!------------------------------------------------------------------

call MPI_FINALIZE(ierr)


!write(GP_print_unit,*)'0: aft mpi_finalize   ierr = ', ierr



end program GPCODE_GA_lmdif_parameter_optimization_test
