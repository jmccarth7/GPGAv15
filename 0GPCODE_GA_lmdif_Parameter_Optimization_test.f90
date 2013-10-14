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




integer :: i
!integer :: i_diversity
integer :: message_len

integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_GP_Generation
integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: jj
integer(kind=4) :: nn

integer(kind=4) :: i_CODE_equation

integer(kind=4) :: n_GP_vars
integer(kind=4) :: nop

integer(kind=4) :: i_GP_best_parent

real(kind=8), allocatable, dimension(:) :: answer
real(kind=8), allocatable, dimension(:) :: output_array

!real(kind=8) :: t1
!real(kind=8) :: t2

!character(200) :: tree_descrip

character(10),parameter :: program_version   = '201308.001'
character(10),parameter :: modification_date = '20131004'
character(50),parameter :: branch  =  'old_elite_parallel_lmdif'


!---------------------------------------------------------------------------------------

GP_para_flag = .FALSE.


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

    write(6,'(/A/)') '0: run lmdif in parallel on each GP generation'
    write(6,'(/A/)') '0: changed RK sub to make it faster'
    write(6,'(/A/)') '0: using the old_elite_scheme in GP_Fit* GP_Tou*, GP_Mut*'
    write(6,'(/A/)') '0: fast mod 1: remove GP diversity and tree printout     '
    write(6,'(/A/)') '0: fast mod 1: remove GP_calc_fit, GP_ranking, summary printout  '
    write(6,'(A)')'0:  run with clean tree call for only processor 0'
    write(6,'(A)')'0:  run with no barrier before call GPCODE  '
    write(6,'(A)')'0:  and with no barrier after  call GPCODE  '
    write(6,'(A)')'0: removed barrier in GPCODE aft bcast of L_stop'

    !------------------------------------------------------
    write(GP_print_unit, '(/3(A,1x,A,1x)//)') &
    '0: GPGACODE program version', trim(program_version), &
    '  branch:', trim( branch ) , &
    '  Last modified on:', trim( modification_date )
    !------------------------------------------------------


    ! read the control input from file  "GPCODE_cntl"

    call read_cntl_stuff( )

    !------------------------------------------------------

    if( L_unit50_output )then
        open( unit_gp_out, file = 'unit50.txt', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_unit50_output


    if( L_GP_log )then
        open( GP_log_unit, file = 'GP_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_GP_log


    if( L_GA_log )then
        open( GA_log_unit, file = 'GA_log', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_GA_log


endif !   myid == 0

!------------------------------------------------------------------

sum_lmdif = 0.0d0

! broadcast the values read in by cpu 0 to others

call bcast1()

!------------------------------------------------------------------

CALL RANDOM_SEED(size = n_seed)

if( .not. allocated( seed ) )then
    ALLOCATE(seed(n_seed))
endif ! .not. allocated( seed )


if( user_input_random_seed > 0 )then

    clock = user_input_random_seed

    if( myid == 0 )then
        write(6,'(A,1x,I12)') &
              '0: user input random seed       clock = ', clock
    endif !   myid == 0

    seed = user_input_random_seed + &
              37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)
else

    CALL SYSTEM_CLOCK(COUNT=clock)

    if( myid == 0 )then
        write(6,'(A,1x,I12)')&
              '0: random seed input clock = ', clock
    endif !   myid == 0

    seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)

endif ! user_input_random_seed > 0

CALL RANDOM_SEED(PUT = seed)

!------------------------------------------------------------------

! wait until everybody has the values

call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!------------------------------------------------------------------

! load table of 2**ilevel - 1  for RK process

call load_pow2_level(  )

!------------------------------------------------------------------

! set the scalar values for the model
                                                                                                                      
! sets:
  n_levels   
! n_functions    
! n_CODE_equations   
! n_trees   
! n_nodes
  

call init_values( 0 )

!------------------------------------------------------------------

if( myid == 0 )then

    call print_values1()

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


GP_Population_Node_Type=-9999

! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]

GP_Individual_Node_Parameters=0.0D0              ! Matrix Operation
GP_Individual_Node_Type=-9999                    ! Matrix Operation
GP_Population_Node_Parameters=0.0D0              ! Matrix Operation

GP_Adult_Population_Parameter_Solution=0.0D0     ! Matrix Operation
GP_Child_Population_Parameter_Solution=0.0D0     ! Matrix Operation

GP_Adult_Population_Node_Type=-9999              ! Matrix Operation
GP_Child_Population_Node_Type=-9999              ! Matrix Operation

!------------------------------------------------------------------

! fill the model arrays

! sets: 
!      Runge_Kutta_Initial_Conditions
!      GP_Individual_Node_Type 
!      GP_Individual_Node_Parameters
!      tree_evaluation
!      Node_Probability

call init_values( 1 )

!------------------------------------------------------------------

! fill a string used to number nodes in print_trees

! this needs to have a variable size since the number of nodes
! may be different in different runs

! create_tree_node_string makes it long enough for n_nodes

call create_tree_node_string()

!------------------------------------------------------------------

! set the desired 'twin experiment' population node type
! and parameter using the info from the set up file

! run the Runge-Kutta model only once with proc 0

! sets:

! GP_Node_Type_Answer
! GP_Node_Parameters_Answer
! GP_Population_Node_Type
! Runge_Kutta_Solution
! Runge_Kutta_Node_Parameters
! Runge_Kutta_Node_Type

! GP_Node_Type_for_Plotting (if L_unit50_output true)
  

call set_answer_arrays( )


!------------------------------------------------------------------------

! then broadcast the R-K result: Runge_Kutta_Solution 


if( myid == 0 )then
    write(GP_print_unit,'(/A/)') &
          '0: time_step   Runge_Kutta_Solution(time_step,1:n_CODE_equations)'
    do  i = 0, n_time_steps
        write(GP_print_unit,'(I6,8(1x,E15.7))') &
              i, (Runge_Kutta_Solution(i,jj), jj = 1,n_CODE_equations )
    enddo ! i
endif ! myid == 0

call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

message_len = ( n_time_steps + 1 ) * n_CODE_equations
call MPI_BCAST( Runge_Kutta_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?


Data_Array=Runge_Kutta_Solution          ! Matrix Operation


!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

! sets:
! Data_Variance
! Data_Variance_inv

call comp_data_variance( )


!--------------------------------------------------------------------------------

! put the desired model parameters in the array:  answer


answer = 0.0d0 ! set all to zero

n_parameters = 0

do  i_CODE_equation=1,n_CODE_equations
    n_parameters=n_parameters+1
    answer(n_parameters)=Runge_Kutta_Initial_Conditions(i_CODE_equation)
enddo ! i_CODE_equation


!--------------------------------------------------------------------------------


! calculate how many parameters total to fit for the specific individual CODE

do  i_tree=1,n_trees
    do  i_node=1,n_nodes

        if( GP_individual_node_type(i_node,i_tree) .eq. 0) then  
            n_parameters=n_parameters+1
            answer(n_parameters)=GP_Individual_Node_Parameters(i_node,i_tree)
        endif ! GP_individual_node_type(i_node,i_tree) .eq. 0

    enddo ! i_node
enddo ! i_tree


!--------------------------------------------------------------------------------

if( myid == 0 )then

    call print_values2( answer )

    !-----------------------------------------------------------------------------

    ! this call calculates sse0,  the sse value for RK solution = 0 for all time steps
    ! note:  sse0 is only used by cpu 0 which does all fitness calculations

    call sse0_calc( )

    !---------------------------------------------------------------------------


    ! open output files

    if( L_GA_output_parameters )then
        open( GA_output_unit, file='GA_output_parameters', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )
    endif ! L_GA_output_parameters

    if( L_GP_output_parameters )then
        open( GP_output_unit, file='GP_output_parameters', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )
    endif ! L_GP_output_parameters


    open( GP_summary_output_unit, file='GP_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

endif ! myid == 0


!---------------------------------------------------------------------------

! calculate n_GP_Asexual_Reproductions, n_GP_Crossovers,  etc.
! from the number of GP individuals and the probabilities such as:
! GP_Asexual_Reproduction_Probability, GP_Crossover_Probability, etc.

call set_modified_indiv( )


!---------------------------------------------------------------------------

! begin the GP generation loop


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the Genetic Programming optimization routine for the Binary Tree Evolution
! with the embedded GA_lmdif parameter optimization scheme
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


generation_loop:&
do  i_GP_Generation=1,n_GP_Generations


    if( myid == 0 )then
        write(GP_print_unit,'(/A/A,1x,I6,1x,A,1x,I6/A/)') &
          '===============================================================================', &
          '0: GP Generation # ',i_GP_Generation,&
          ' is underway.   n_Nodes * n_Trees = ', n_Nodes*n_Trees, &
          '==============================================================================='
    endif ! myid == 0


    GP_Child_Population_Node_Type = GP_Adult_Population_Node_Type  


    ! Run_GP_Calculate_Fitness determines if the new GP child
    ! has to be put through the R-K process for parameter optimization

    ! Run_GP_Calculate_Fitness is true for all individuals in generation 1

    ! Run_GP_Calculate_Fitness will be FALSE for generations > 1 
    ! if the individual did not change on the last generation (so it needs no recalculation)

    !Run_GP_Calculate_Fitness= .TRUE.  ! debug only .false.

    Run_GP_Calculate_Fitness= .false.


    ! randomly create the initial tree arrays for each individual and
    ! send them all to GA_lmdif for parameter optimization on generation 1

    if( i_GP_Generation .eq. 1) then

        ! determines if the new GP child
        !has to be sent to GA_lmdif for parameter optimization

        Run_GP_Calculate_Fitness=.true.

        !---------------------------------------------------------------------------------

        if( myid == 0 )then

            write(GP_print_unit,'(/A,1x,I6/)') &
              '0: call GP_Tree_Build   Generation = ',i_GP_Generation

            ! initializes the GP_Child_Population_Node_Type array with random trees

            ! sets GP_Adult_Population_Node_Type to GP_Child_Population_Node_Type


            !t1 = MPI_Wtime()

            call GP_Tree_Build

            !t2 = MPI_Wtime()


            !write(GP_print_unit,'(A,1x,E15.7)') &
            !  '0: time spent in GP_Tree_Build  = ', t2 - t1

            !debug only! set all GP tree models to the "truth" model            
            !debug only do  i_GP_individual = 1, n_GP_Individuals
            !debug only     GP_Adult_Population_Node_Type(:,:,i_GP_individual) = &
            !debug only                       GP_Node_Type_Answer(:,:) ! debug only
            !debug only enddo

        endif ! myid == 0


        !---------------------------------------------------------------------------------

        ! broadcast GP_Adult_Population_Node_Type 

        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6)') &
        !      '0: broadcast  GP_Adult_Population_Node_Type Generation = ',i_GP_Generation
        !endif ! myid == 0


        call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

        message_len = n_GP_Individuals * n_Nodes * n_Trees
        call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

        call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?



        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,E15.7)') &
        !      '0: time spent in bcast GP_Adult_Pop_Node_Type 2 = ', t2 - t1
        !endif ! myid == 0

        !---------------------------------------------------------------------------------

        ! compute a "diversity index" which characterizes each individual with a
        ! number derived from the number of nodes, etc.

        !if( myid == 0 )then
        !    call GP_calc_diversity_index( n_GP_individuals,  &
        !                                  GP_Adult_Population_Node_Type, &
        !                                  i_diversity, i_gp_generation )
        !endif ! myid == 0

        !-----------------------------------------------------------------------------


    else !  i_GP_Generation > 1

        ! create the next 'generation' of tree structures using either:

        !    i)  GP Fitness-Proportionate Asexual Reproduction;
        !   ii)  GP Tournament-Style Sexual Reproduction, and;
        !  iii)  GP Mutation

        ! then broadcast the arrays to all processors

        if( myid == 0 )then

            !tree_descrip =  ' GP_Adult trees before call selection routines '
            !call print_trees( i_GP_generation, 1, n_GP_individuals, GP_Adult_Population_Node_Type, &
            !                  trim( tree_descrip )  )

            !----------------------------------------------------------------------------------

            !   i) Carry out "GP Fitness-Proportionate Reproduction"

            !      randomly replaces values of individuals in child arrays 
            !      with values from the adult arrays of fitter individuals 

            !   uses: 
            !   GP_Integrated_Population_Ranked_Fitness
            !   GP_Adult_Population_Node_Type
            !   GP_Adult_Population_Parameter_Solution
            !   GP_Child_Population_SSE

            !   sets:
            !   GP_Child_Population_Node_Type
            !   GP_Child_Population_Parameter_Solution
            !   GP_Child_Population_SSE

            if( n_GP_Asexual_Reproductions .gt. 0 )then

                write(GP_print_unit,'(/A,1x,I6)') &
                      '0: call GP_Fitness_Proportionate_Asexual_Reproduction  &
                      &n_GP_Asexual_Reproductions =', n_GP_Asexual_Reproductions


                !t1 = MPI_Wtime()

                call GP_Fitness_Proportionate_Asexual_Reproduction

                !t2 = MPI_Wtime()

                !write(GP_print_unit,'(A,1x,E15.7)') &
                !  '0: time spent in GP_Fitness_Proportionate_Asexual_Reproduction  = ', t2 - t1

                !write(GP_print_unit,'(/A)') &
                !      '0:aft  call GP_Fitness_Proportionate_Asexual_Reproduction '

                !tree_descrip=' GP_Child trees after call to &
                !              &GP_Fitness_Proportionate_Asexual_Reproduction'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Child_Population_Node_Type, &
                !                    trim( tree_descrip )  )

            endif !  n_GP_Asexual_Reproductions .gt. 0

            !----------------------------------------------------------------------------------

            !  ii) Carry out "GP Tree Crossover" Operations
            !      Using Tournament-Style Sexual Reproduction Selection
            !      and randomly use it to replace the new children


            ! uses:
            !    GP_Adult_Population_Node_Type
            !    GP_Adult_Population_SSE
 
            ! sets:
            !    GP_Child_Population_Node_Type
            !    Run_GP_Calculate_Fitness ( to true for modified individuals )
   

            if( n_GP_Crossovers .gt. 0 )then

                write(GP_print_unit,'(A,1x,I6)') &
                      '0: call GP_Tournament_Style_Sexual_Reproduction n_GP_Crossovers =', &
                                                                       n_GP_Crossovers

                !t1 = MPI_Wtime()

                call GP_Tournament_Style_Sexual_Reproduction

                !t2 = MPI_Wtime()


                !write(GP_print_unit,'(A,1x,E15.7)') &
                !  '0: time spent in GP_Tournament_Style_Sexual_Reproduction = ', t2 - t1

                !write(GP_print_unit,'(/A)') &
                !      '0: aft  call GP_Tournament_Style_Sexual_Reproduction '

                !tree_descrip = ' GP_Child trees after call to &
                !                  &GP_Tournament_Style_Sexual_Reproduction'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Child_Population_Node_Type, &
                !                    trim( tree_descrip )  )

            endif !  n_GP_Crossovers .gt. 0


            !----------------------------------------------------------------------------------

            !   iii) Carry out "GP Parameter Mutation" Operations

            ! uses:

            !  GP_Adult_Population_Node_Type

            ! sets:
                                                                                                              
            !  GP_Child_Population_Node_Type
                                                                                                              
            !  Run_GP_Calculate_Fitness  ( to true for modified individuals )
      

            if( n_GP_Mutations .gt. 0 )then

                write(GP_print_unit,'(A,1x,I6)')'0: call GP_Mutations n_GP_Mutations =', &
                                                                      n_GP_Mutations

                !t1 = MPI_Wtime()

                call GP_Mutations

                !t2 = MPI_Wtime()

                !write(GP_print_unit,'(A,1x,E15.7)') &
                !  '0: time spent in GP_Mutations = ', t2 - t1

                !write(GP_print_unit,'(/A)')'0:aft call GP_Mutations '

                !tree_descrip =  ' GP_Child trees after call to GP_Mutations'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Child_Population_Node_Type, &
                !                    trim( tree_descrip )  )

            endif !  n_GP_Mutations .gt. 0


            !---------------------------------------------------------------------------

            !   Move over any newly created children into the adult arrays                 

            GP_Adult_Population_Node_Type = GP_Child_Population_Node_Type
            GP_Adult_Population_SSE       = GP_Child_Population_SSE

            write(GP_print_unit,'(/A)')&
                  '0:aft  move Child_Node_Type and SSE to Adult'

            !---------------------------------------------------------------------------

            ! calculate the diversity index for each individual for generations > 1

            !write(GP_print_unit,'(/A)')&
            !      '0: call GP_calc_diversity_index '

            !t1 = MPI_Wtime()
            !call GP_calc_diversity_index( n_GP_individuals, &
            !                              GP_Child_Population_Node_Type, &
            !                              i_diversity, i_gp_generation )
            !t2 = MPI_Wtime()

            !write(GP_print_unit,'(A,1x,E15.7)') &
            !  '0: time spent in GP_calc_diversity_index = ', t2 - t1

            !write(GP_print_unit,'(/A)')&
            !      '0: aft call GP_calc_diversity_index '


            !-----------------------------------------------------------------------------------------
        
            ! print trees before tree clean
        
            !if( i_GP_generation == 1                                  .or. &
            !    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            !    i_GP_generation == n_GP_generations                          )then
            !    tree_descrip =  ' trees before call to GP_Clean_Tree_Nodes'
            !    call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !                      GP_Adult_Population_Node_Type, &
            !                      trim( tree_descrip )  )
            !endif ! i_GP_generation == 1
        
        
            !-----------------------------------------------------------------------------------------
        
            ! GP_Clean_Tree_Nodes sweeps through the GP_Adult_Population_Node_Type array
            ! to replace function nodes that have both terminals set as parameters
            ! and to set the replaced node to a parameter itself
        
            ! sets:
            ! GP_Adult_Population_Node_Type
        
            write(GP_print_unit,'(/A,1x,I6)') &
                  '0: call GP_Clean_Tree_Nodes   Generation =', i_GP_Generation
            write(GP_print_unit,'(A)') &
                  '0: GP_Clean_Tree_Nodes  modifies GP_Adult_Population_Node_Type'
    
    
            call GP_Clean_Tree_Nodes
    
    
            !write(GP_print_unit,'(/A,1x,I6/)') &
            !         '0: aft call GP_Clean_Tree_Nodes   Generation =', i_GP_Generation
        
            !-----------------------------------------------------------------------------------------
        
            ! print trees after call to GP_Clean_Tree_Nodes
        
            !!if( i_GP_generation == 1                                  .or. &
            !!    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            !!    i_GP_generation == n_GP_generations                          )then
        
            !    !t1 = MPI_Wtime()
        
            !    tree_descrip =  ' trees after call to GP_Clean_Tree_Nodes'
            !    call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !                      GP_Adult_Population_Node_Type, &
            !                      trim( tree_descrip )  )
        
            !    !t2 = MPI_Wtime()
        
            !    !write(GP_print_unit,'(A,1x,E15.7)') &
            !    !  '0: time spent in print_trees = ', t2 - t1
        
            !!endif ! i_GP_generation == 1 .or. ...
        
            !-----------------------------------------------------------------------------------------
            !>>>>>>>>>> jjm 20130417
            !if( i_GP_generation == 1                                 .or. &
            !    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. &
            !    i_GP_generation == n_GP_generations                         )then
            !
            !    write(GP_print_unit,'(/A, 1x, I6/)') &
            !          '0: after call to GP_Clean_Tree_Nodes i_GP_generation =',i_GP_generation
        
            ! print node type information for each GP individual
        
            !    call print_gp_node_type_parm( )
            !endif !  i_GP_generation == 1 ...
            !<<<<<<<<<<< jjm 20130417
        

        endif ! myid == 0


        !------------------------------------------------------------------------------------


        ! broadcast:

        ! GP_Child_Population_Node_Type
        ! GP_Adult_Population_Node_Type
        ! Parent_Tree_Swap_Node_Type (removed from broadcast)
        ! GP_Child_Individual_SSE
        ! GP_Adult_Individual_SSE
        ! GP_Integrated_Population_Ranked_Fitness
        ! GP_Population_Ranked_Fitness
        ! Run_GP_Calculate_Fitness 

        !if( myid == 0 )then
        !    write(GP_print_unit,'(/A)')&
        !          '0: call bcast2 '
        !endif ! myid == 0

        call bcast2()

        !if( myid == 0 )then
        !    write(GP_print_unit,'(/A)')&
        !          '0: aft call bcast2 '
        !endif ! myid == 0

    endif ! i_GP_Generation .eq. 1


    !------------------------------------------------------------------------------------

    ! removed because GP_Tree_Clean now done before call to bcast2 
    ! which broadcasts GP_Adult_Population_Node_Type


    !! broadcast GP_Adult_Population_Node_Type
    !call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?
    !message_len = n_GP_Individuals * n_Nodes * n_Trees
    !call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
    !                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
    !call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?
    !!if( myid == 0 )then
    !!    write(GP_print_unit,'(A,1x,I6)') &
    !!      '0: aft broadcast  GP_Adult_Pop_Node_Type  Generation = ',&
    !!                                            i_GP_Generation
    !!endif ! myid == 0


    !-----------------------------------------------------------------------------------------


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! GA_lmdif subroutine segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    if( myid == 0 )then
        write(GP_print_unit,'(/A/A)')&
              '0:-----------------------------------------------------------------', &
              '0:  before starting loop on GP individuals '
        write(GP_print_unit,'(A/A,4x,L1)')&
              '0: Are there any individuals to calculate fitness for? ', &
              '0: any( Run_GP_Calculate_Fitness ) = ', any( Run_GP_Calculate_Fitness )
    endif !  myid == 0


    !-----------------------------------------------------------------------------------

    ! exit the generation loop since 
    ! there are no more individuals to evaluate fitness for

    if( .not.  any( Run_GP_Calculate_Fitness ) ) exit generation_loop

    !-----------------------------------------------------------------------------------


    gp_ind_loop:&
    do  i_GP_individual=1,n_GP_individuals

        if( myid == 0 )then
            write(GP_print_unit,'(//A)')&
            '0:----------------------------------------------------------------------'
            write(GP_print_unit,'(A,2(1x,I6),3x,L1)') &
                  '0: i_GP_Gen, i_GP_indiv, Run_GP_Calculate_Fitness', &
                      i_GP_Generation, i_GP_individual, Run_GP_Calculate_Fitness(i_GP_Individual)
            write(GP_print_unit,'(A/)')&
            '0:----------------------------------------------------------------------'

        endif !  myid == 0


        ! run GPCODE_... to evaluate this individual  if Run_GP_Calculate_Fitness is true 

        if( Run_GP_Calculate_Fitness(i_GP_Individual) ) then

            GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) = &
               GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual)


            ! these get set randomly in the GA-lmdif search algorithm ( in GPCODE* )

            GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) = 0.0d0


            !------------------------------------------------------------------------------

            ! calculate how many parameters total to fit for the specific individual CODE
            ! and save this number in GP_Individual_N_GP_param(i_GP_individual)

            n_GP_Parameters = n_code_equations

            do  i_Tree=1,n_Trees

                do  i_Node=1,n_Nodes

                    if( GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0) then  ! a set parameter
                        n_GP_Parameters=n_GP_Parameters+1
                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

                    !if( myid == 0 )then
                    !    if( GP_Individual_Node_Type(i_Node,i_Tree) > -9999 )then
                    !        write(GP_print_unit,'(A,5(1x,I6))')&
                    !        '0: i_GP_indiv, i_tree, i_node, &
                    !          &GP_Indiv_Node_Type, n_GP_params',&
                    !           i_GP_individual, i_tree, i_node, &
                    !           GP_Individual_Node_Type(i_Node,i_Tree), n_GP_parameters
                    !    endif ! GP_Individual_Node_Type(i_Node,i_Tree) > -9999
                    !endif !  myid == 0

                enddo ! i_node

            enddo ! i_tree

            if( myid == 0 )then
                write(GP_print_unit,'(/A,3(1x,I6))') &
                  '0: i_GP_individual, n_trees, n_nodes ', &
                      i_GP_individual, n_trees, n_nodes
                write(GP_print_unit,'(A,1x,I6,3x,A,1x,I6/)')&
                  '0: for i_GP_Individual', i_GP_Individual, &
                  'the number of parameters is:  n_GP_parameters =', n_GP_parameters
            endif !  myid == 0


            !------------------------------------------------------------------------
            ! if n_GP_parameters > n_maximum_number_parameters, this is an error
            ! for now, just print a warning and set n_GP_parameters = n_max...

            if( n_GP_parameters > n_maximum_number_parameters )then

                write(GP_print_unit,'(/A,5(1x,I10)/)')&
                  '0:WARNING  myid, i_GP_generation, i_GP_Individual, &
                      &n_GP_parameters, n_maximum_number_parameters', &
                      myid, i_GP_generation, i_GP_Individual, &
                       n_GP_parameters, n_maximum_number_parameters

                 !!!n_GP_parameters = n_maximum_number_parameters
            endif
            !------------------------------------------------------------------------


            GP_Individual_N_GP_param(i_GP_individual) = n_GP_parameters


            !------------------------------------------------------------------------------


            ! calculate how many variables are in the tree


            n_GP_vars = 0

            do  i_Tree=1,n_Trees

                do  i_Node=1,n_Nodes

                    if( GP_Individual_Node_Type(i_Node,i_Tree) < 0  .and. &
                        GP_Individual_Node_Type(i_Node,i_Tree) > -9999  ) then
                        n_GP_vars = n_GP_vars + 1
                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) > 0 ....

                    !if( myid == 0 )then
                    !    if( GP_Individual_Node_Type(i_Node,i_Tree) > -9999 )then
                    !        write(GP_print_unit,'(A,5(1x,I6))')&
                    !        '0: i_GP_indiv, i_tree, i_node, &
                    !          &GP_Indiv_Node_Type, n_GP_vars ', &
                    !           i_GP_individual, i_tree, i_node, &
                    !           GP_Individual_Node_Type(i_Node,i_Tree), n_GP_vars
                    !    endif ! GP_Individual_Node_Type(i_Node,i_Tree) > -9999
                    !endif !  myid == 0

                enddo ! i_node

            enddo ! i_tree

            !------------------------------------------------------------------------------

            if( myid == 0 )then
                !write(GP_print_unit,'(/A,3(1x,I6))') &
                !      '0: i_GP_individual, n_trees, n_nodes ', &
                !          i_GP_individual, n_trees, n_nodes
                write(GP_print_unit,'(A,1x,I6,2(/3x,A,1x,I6)/)')&
                      '0: for i_GP_Individual', i_GP_Individual, &
                      'the number of parameters is: ', n_GP_parameters, &
                      'the number of variables  is: ', n_GP_vars
            endif !  myid == 0


            !-------------------------------------------------------------------

            ! cycle the i_GP_individual loop if there are no GP parameters
            ! or if n_GP_parameters <=  n_code_equations

            !orig if( n_GP_parameters == 0 .or. &
            !orig     n_GP_parameters < n_code_equations ) then

            if( n_GP_parameters == 0 .or. &
                n_GP_parameters > n_maximum_number_parameters .or.  &
                n_GP_parameters <=  n_code_equations ) then   ! new jjm 20130814

                if( myid == 0 )then
                    write(GP_print_unit,'(A,1x,I6)')&
                          '0: skipping this i_GP_Individual --&
                          &  the number of parameters is ', n_GP_parameters
                endif !  myid == 0

                !------------------------------------------------------------------------

                ! set SSE values for this rejected individual so that
                ! its fitness will be very small

                GP_Child_Individual_SSE(i_GP_Individual) = 1.0D13
                GP_Child_Population_SSE(i_GP_Individual) = 1.0D13

                GP_Adult_Individual_SSE(i_GP_Individual) = 1.0D13
                GP_Adult_Population_SSE(i_GP_Individual) = 1.0D13

                GP_Individual_N_GP_param(i_GP_individual) = 0
                n_GP_parameters                           = 0


                if( myid == 0 )then
                    write(GP_print_unit,'(/A/A/A,2(1x,I5), 1x, E15.7)')&
                          '0: rejected for n_GP_parameters <=  n_code_equations',&
                          'or for n_GP_parameters > n_maximum_number_parameters',&
                           'i_GP_gen,i_GP_indiv,GP_Child_Pop_SSE(i_GP_Indiv)  ', &
                               i_GP_generation, i_GP_individual, &
                               GP_Child_Population_SSE(i_GP_Individual)
                endif !  myid == 0

                !------------------------------------------------------------------------

                cycle gp_ind_loop

                !!call GP_tree_reset( i_GP_individual )

            endif ! n_GP_parameters == 0

            !-------------------------------------------------------------------


            if( myid == 0 )then
                write(GP_print_unit,'(//A)') &
                 '0: call GPCODE_GA_lmdif_Parameter_Optimization routine'
                write(GP_print_unit,'(A,2(1x,I6))') &
                 '0: i_GP_Generation, i_GP_individual', &
                     i_GP_Generation, i_GP_individual

                if( L_ga_print )then
                    write(GA_print_unit,'(//A/A)') &
                     '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>&
                     &>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
                     '0: call GPCODE_GA_lmdif_Parameter_Optimization routine'
                    write(GA_print_unit,'(A,2(1x,I6)/)') &
                     '0: i_GP_Generation, i_GP_individual',&
                         i_GP_Generation, i_GP_individual
                endif ! L_ga_print
            endif ! myid == 0


            !----------------------------------------------------------------------------

            !call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

            ! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND
            ! LINK THE SSE OUTPUT TO THE ARRAY AT THE END
            ! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEED TO BE PULLED OUT

            ! GPCODE_GA_lmdif_Parameter_Optimization sets:

            !      individual_fitness
            !   GP_Individual_Initial_Conditions
            !   GP_Individual_Node_Parameters

            ! these arrays are broadcast in GPCODE_GA...


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call GPCODE_GA_lmdif_Parameter_Optimization( i_GP_Generation,i_GP_individual )

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !if( myid == 0 )then
            !    write(GP_print_unit,'(A,1x,I6)') &
            !     '0: aft1 call GPCODE_GA_lmdif_Parameter_Optimization routine'
            !endif ! myid == 0

            ! wait until everybody has the values

            !call MPI_BARRIER( MPI_COMM_WORLD, ierr )

            !----------------------------------------------------------------------------

            if( myid == 0 )then
                write(GP_print_unit,'(A,1x,I6/)') &
                 '0: aft2 call GPCODE_GA_lmdif_Parameter_Optimization routine'
            endif ! myid == 0

            GP_population_fitness(i_GP_individual) = individual_fitness

            !--------------------------------------------------------------------------------

            if( myid == 0 )then

                write(GP_print_unit,'(/A,2(1x,I6), 1x, E15.7)')&
                      '0:1 i_GP_gen, i_GP_indiv, GP_pop_fit(i_GP_indiv) ', &
                           i_GP_generation, i_GP_individual, &
                           GP_population_fitness(i_GP_individual)

                ! print side-by-side comparisons of
                ! starting values and values from optimization

                write(GP_print_unit,'(/A/)') &
                  '0:           truth value           GP_individual_initial_conditions'

                do  i_CODE_equation=1,n_CODE_equations

                    write(GP_print_unit,'(I6,1x, E20.10, 4x, E20.10)') &
                          i_CODE_equation, &
                          Runge_Kutta_Initial_Conditions(i_CODE_equation),  &
                          GP_individual_initial_conditions(i_CODE_equation)

                    if( L_GA_output_parameters )then
                        write(GA_output_unit,'(E24.16)') &
                              GP_individual_initial_conditions(i_CODE_equation)
                    endif ! L_GA_output_parameters )then



                    !output_array( i_CODE_equation ) = &
                    !           GP_individual_initial_conditions(i_CODE_equation)


                enddo ! i_CODE_equation

            endif !  myid == 0

            !--------------------------------------------------------------------------------

            ! compute GP_Child_Individual_SSE(i_GP_Individual)
            ! and     GP_Child_Population_SSE(i_GP_Individual)

            ! use the best sse from the GPCODE subroutine

            GP_Child_Individual_SSE(i_GP_Individual) = Individual_SSE_best_parent
            GP_Child_Population_SSE(i_GP_Individual) = Individual_SSE_best_parent

            GP_Adult_Individual_SSE(i_GP_Individual) = GP_Child_Individual_SSE(i_GP_Individual)
            GP_Adult_Population_SSE(i_GP_Individual) = GP_Child_Population_SSE(i_GP_Individual)


            if( myid == 0 )then
                write(GP_print_unit,'(/A,2(1x,I5), 1x, E15.7)')&
                      '0:3 i_GP_gen,i_GP_indiv,GP_Child_Pop_SSE(i_GP_Indiv)  ', &
                           i_GP_generation, i_GP_individual, &
                           GP_Child_Population_SSE(i_GP_Individual)
            endif !  myid == 0

            !-----------------------------------------------------------------------------------


            ! set the GA_lmdif-optimized initial condition array

            GP_Population_Initial_Conditions(1:n_CODE_Equations, i_GP_Individual) = &
                GP_Individual_Initial_Conditions(1:n_CODE_Equations) ! Matrix Operation



            !if( myid == 0 )then
            !    write(GP_print_unit,'(A,2(1x,I5), 5(1x, E15.7))')&
            !          '0:3 i_GP_gen,i_GP_indiv,&
            !          &GP_Pop_Init_Cond(i_GP_Indiv,1:n_CODE_Eq)', &
            !               i_GP_generation, i_GP_individual, &
            !           GP_Population_Initial_Conditions(1:n_CODE_Equations, i_GP_Individual)
            !endif !  myid == 0

            !------------------------------------------------------------------------------

            ! load the array GP_Adult_Population_Parameter_Solution 


            !  first, load  the initial conditions

            GP_Adult_Population_Parameter_Solution(1:n_CODE_equations, i_GP_Individual) = &
                GP_Individual_Initial_Conditions(1:n_CODE_Equations) 


            ! then the rest of the parameters 

            nn = n_code_equations

            do  i_Tree=1,n_Trees
                do  i_Node=1,n_Nodes

                    if( GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0) then 

                        nn = nn + 1

                        GP_Adult_Population_Parameter_Solution( nn, i_GP_Individual ) = &
                                GP_Individual_Node_Parameters(i_Node,i_Tree) 

                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

                enddo ! i_node
            enddo ! i_tree

            !-------------------------------------------------------------------------------------

            ! set the GA_lmdif-optimized CODE parameter set array

            GP_Population_Node_Parameters(1:n_Nodes,1:n_Trees, i_GP_Individual) = &
                GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) ! Matrix Operation

            !------------------------------------------------------------------------------

            GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
                             GP_Individual_Node_Type(1:n_Nodes,1:n_Trees)

            GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
                             GP_Individual_Node_Type(1:n_Nodes,1:n_Trees)

            !------------------------------------------------------------------------------

            if( myid == 0 )then

                !----------------------------------------------------------------------------

                ! print the node parameters (if there are any)

                write(GP_print_unit,'(/A/)') &
                       '0:  tree  node  Runge_Kutta_Node_Params&
                                  &   GP_population_node_params'
                do  i_tree=1,n_trees
                    do  i_node=1,n_nodes
                        if( GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0 ) then
                            write(GP_print_unit,'(2(1x,I6), 1x, E20.10, 4x, E20.10)') &
                             i_tree, i_node, &
                             Runge_Kutta_Node_Parameters(i_node,i_tree), &
                             GP_population_node_parameters(i_node,i_tree,i_GP_individual)
                        endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

                    enddo ! i_node
                enddo  ! i_tree

                !---------------------------------------------------------------------------------

                write(GP_print_unit,'(/A)') &
                      '0:---------------------------------------------------------------'

                !write(GP_print_unit,'(/A)')  &
                !      '0: after loading GP_Pop arrays with GP_indiv array values '

                !write(GP_print_unit,'(A,2(1x,I6),4x,L1)') &
                !  '0: i_GP_gen,i_GP_indiv, Run_GP_Calculate_Fitness(i_GP_Indiv) ', &
                !      i_GP_generation,i_GP_individual, &
                !      Run_GP_Calculate_Fitness(i_GP_Individual)

            endif !  myid == 0

        endif !   Run_GP_Calculate_Fitness(i_GP_Individual)


        if( myid == 0 )then

            ! this prints a summary of the initial conditions,
            ! parameters,  and node types for this individual,
            ! after being optimized in GPCODE*opt
            ! and writes the tree to the summary file

            call summary_GP_indiv( i_GP_generation, i_GP_individual )

        endif !  myid == 0


    enddo  gp_ind_loop    !   i_GP_individual


    !-------------------------------------------------------------------------------------

    if( myid == 0 )then
        write(GP_print_unit,'(/A,1x,I6)') &
          '0: before GP_para_lmdif_process   i_GP_generation =', i_GP_Generation

        write(GP_print_unit, '(/A )') &
              '0:i_GP_Individual,  GP_Individual_N_GP_param '
        do  i_GP_individual = 1, n_GP_individuals
            write(GP_print_unit, '(I6,10x,I6 )') &
                  i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual)
        enddo


        write(GP_print_unit, '(/A )') &
              '0:i_GP_Individual,  GP_Adult_Population_SSE'
        do  i_GP_individual = 1, n_GP_individuals
            write(GP_print_unit, '(I6,10x,E15.7)') &
                  i_GP_Individual,  GP_Adult_Population_SSE(i_GP_Individual)
        enddo
    endif ! myid == 0

    !-------------------------------------------------------------------------------------

    !  call routine to run lmdif in parallel on all the GP individuals

    !  GP_para_lmdif_process returns arrays to be used in GP_calc_fitness:

    ! uses:
    !  GP_Population_Initial_Conditions
    !  GP_Adult_population_Node_Type
    !  GP_population_node_parameters

    ! sets:
    !  GP_Population_Initial_Conditions
    !  GP_population_node_parameters
    !  child_parameters
    !  GP_child_individual_SSE
    !  individual_quality
    !  GP_n_parms
    !  GP_Adult_population_parameter_solution 

    if( myid == 0 )then
        write(GP_print_unit,'(/A)')&
          '0:-----------------------------------------------------------------'
        write(GP_print_unit,'(A,1x,I6)') &
          '0: call GP_para_lmdif_process i_GP_generation =', &
                                         i_GP_Generation
        write(GP_print_unit,'(A/)')&
          '0:-----------------------------------------------------------------'
    endif ! myid == 0


    GP_para_flag = .TRUE. 

    call  GP_para_lmdif_process( i_GP_generation )

    GP_para_flag = .FALSE.

    !-------------------------------------------------------------------------------------

                                                                                                                        
    GP_Child_Population_Parameter_Solution(1:n_Maximum_Number_Parameters, 1:n_GP_Individuals) = &
    GP_Adult_Population_Parameter_Solution(1:n_Maximum_Number_Parameters, 1:n_GP_Individuals)


    !-------------------------------------------------------------------------------------

    ! do fitness calculations for this GP generation

    if( myid == 0 )then

        write(GP_print_unit,'(/A)')&
          '0:#################################################################'
        write(GP_print_unit,'(A,1x,I6)') &
          '0: call GP_calc_fitness i_GP_generation =', &
                                   i_GP_Generation
        write(GP_print_unit,'(A/)')&
          '0:#################################################################'

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  debug
        ! print the node parameters (if there are any)

        write(GP_print_unit,'(/A/)') &
           '0:  i_GP_indiv tree  node  GP_population_node_params'
        do  i_GP_individual = 1, n_GP_individuals
            do  i_tree=1,n_trees
                do  i_node=1,n_nodes
                    if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) .eq. 0 ) then
                        if( GP_population_node_parameters(i_node,i_tree,i_GP_individual) > 0.0d0 )then
                            write(GP_print_unit,'(3(1x,I6),  4x, E20.10)') &
                             i_GP_individual, i_tree, i_node, &
                             GP_population_node_parameters(i_node,i_tree,i_GP_individual)
                        endif ! GP_population_node_parameters(i_node,i_tree,i_GP_individual) > 0.0d0
                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

                enddo ! i_node
            enddo  ! i_tree
        enddo ! i_GP_individual
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  debug

        !----------------------------------------------------------------------------------


        ! do fitness calculations for this GP generation

        !t1 = MPI_Wtime()

        call GP_calc_fitness( i_GP_generation, output_array, &
                              i_GP_best_parent, nop )

        !t2 = MPI_Wtime()
        !----------------------------------------------------------------------------------


        !write(GP_print_unit,'(A,1x,E15.7)') &
        !        '0: time spent in GP_calc_fitness = ', t2 - t1

        write(GP_print_unit,'(/A)')&
          '0:#################################################################'
        write(GP_print_unit,'(A,1x,I6)') &
          '0: aft call GP_calc_fitness i_GP_generation =', &
                                       i_GP_Generation
        write(GP_print_unit,'(A/)')&
          '0:#################################################################'


        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  debug

        ! print the node parameters (if there are any)

        write(GP_print_unit,'(/A/)') &
           '0:  i_GP_indiv tree  node  GP_population_node_params'
        do  i_GP_individual = 1, n_GP_individuals
            do  i_tree=1,n_trees
                do  i_node=1,n_nodes
                    if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) .eq. 0 ) then
                        if( GP_population_node_parameters(i_node,i_tree,i_GP_individual) > 0.0d0 )then
                            write(GP_print_unit,'(3(1x,I6),  4x, E20.10)') &
                             i_GP_individual, i_tree, i_node, &
                             GP_population_node_parameters(i_node,i_tree,i_GP_individual)
                        endif ! GP_population_node_parameters(i_node,i_tree,i_GP_individual) > 0.0d0
                    endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

                enddo ! i_node
            enddo  ! i_tree
        enddo ! i_GP_individual
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  debug


    endif ! myid == 0

    !---------------------------------------------------------------------------


    !if( myid == 0 )then
    !    call print_gp_node_type_parm( )
    !endif ! myid == 0

    !---------------------------------------------------------------------------

    ! broadcast:

    !  GP_Adult_Individual_SSE
    !  GP_population_node_parameters
    !  GP_Population_Ranked_Fitness
    !  GP_Integrated_Population_Ranked_Fitness

    call bcast3( )




enddo generation_loop !  i_GP_Generation


call MPI_BARRIER( MPI_COMM_WORLD, ierr )


if( myid == 0 )then
    write(GP_print_unit,'(/A/)') '0: after i_GP_generation loop  '


    !---------------------------------------------------------------------------


    !tree_loop:&
    !do  i_tree=1,n_trees
    !    node_loop:&
    !    do  i_node=1,n_nodes
    !        if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0 )then
    !            write(GP_print_unit,'(2x,2(1x,I6), 1x, E20.10, 4x, E20.10)') &
    !                  i_tree, i_node,  &
    !                  GP_population_node_parameters(i_node,i_tree,i_GP_Best_Parent)
    !        endif ! GP_Adult_Pop_Node_Type(i_Node,i_Tree,i_GP_Best_Parent) == 0
    !        write(GP_print_unit,'(3(1x,I6))') i_tree, i_node, nop
    !    enddo node_loop ! i_node
    !enddo tree_loop ! i_tree

    write(GP_print_unit,'(A,1x,I6,1x,E15.7/)') &
          '0: i_GP_best_parent, GP_child_individual_sse( i_GP_best_parent ) ', &
              i_GP_best_parent, GP_child_individual_sse( i_GP_best_parent )

    write(GP_print_unit,'(A/)') &
          '0: call GP_select_best_RK_lmdif_result to run lmdif for best parent'


    ! uses:
    !  output_array
    !  GP_Adult_Population_Node_Type

    ! sets:

    !  individual_quality
    !  Parent_Parameters
    !  GP_individual_node_type
    !  parent_parameters
    !  GP_child_individual_SSE
    !  GP_individual_ranked_fitness
    !  child_parameters
    !  GP_Individual_Initial_Conditions
    !  GP_Individual_Node_Parameters


    call GP_select_best_RK_lmdif_result( i_GP_best_parent, output_array , nop )


    !---------------------------------------------------------------------------


    write(GP_print_unit,'(/A,1x,E15.7/)') '0: sum of time spent in lmdif = ', sum_lmdif

endif ! myid == 0


call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! plot results

!if( myid == 0 )then
!    Lplot = .true.
!    if( Lplot ) call plot_results(Runge_Kutta_Solution,n_time_steps, n_CODE_equations )
!endif ! myid == 0

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


if( L_unit50_output )then
    close( unit_gp_out )
endif ! L_unit50_output

if( L_GP_log )then
    close( GP_log_unit )
endif ! L_GP_log


if( L_GA_log )then
    close( GA_log_unit )
endif ! L_GA_log


if( L_GA_output_parameters )then
    close( GA_output_unit )
endif ! L_GA_output_parameters

if( L_GP_output_parameters )then
    close( GP_output_unit )
endif ! L_GP_output_parameters


close( GP_summary_output_unit )

!------------------------------------------------------------------

! deallocate variable dimension arrays

call deallocate_arrays1( )

deallocate( answer )
deallocate( output_array )

!------------------------------------------------------------------

call MPI_FINALIZE(ierr)


!write(GP_print_unit,*)'0: aft mpi_finalize   ierr = ', ierr


stop

end program GPCODE_GA_lmdif_parameter_optimization_test
