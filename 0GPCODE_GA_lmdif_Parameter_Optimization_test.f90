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

integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_GP_Generation
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

real(kind=8), allocatable, dimension(:) :: answer
real(kind=8), allocatable, dimension(:) :: output_array

!real(kind=8) :: t1
!real(kind=8) :: t2

character(200) :: tree_descrip

character(10),parameter :: program_version   = '201401.002'
character(10),parameter :: modification_date = '20140303'
character(50),parameter :: branch  =  'master'
!!type(Tree_Node_Pointer) :: t1, t2, t3, t4
!!type(Tree_Node_Pointer), dimension(:), allocatable :: T1_Nodes, T2_Nodes
!!!real(kind=8) :: r
!!integer(kind=4) :: index, swap_node_1, swap_node_2


!---------------------------------------------------------------------------------------

GP_para_flag = .FALSE.
ierror = 0

GP_minSSE_Individual = 0
GP_minSSE_generation = 0

!--------------------------------------------------------------
! current setup
! lmdif runs only on best individual of each generation
! no replacement for bad points in GP*n
! no retry in setup_run_fcn or setup_run_lmdif
! in random_real, boundary is 0.1
! --  0.1 for random range [0,1], 0.9 for random range [0,50]
! max interations in lmdif is 100
!--------------------------------------------------------------


call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


if( myid == 0 )then
    write(6,'(/A,1x,I2/)') '0: numprocs = ', numprocs
endif ! myid == 0


!------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(/A)') '0: call Global_Setup '
endif ! myid == 0

call Global_Setup()

if( myid == 0 )then
    write(6,'(A)') '0: aft call Global_Setup '
endif ! myid == 0

!------------------------------------------------------------------

GP_para_flag = .FALSE.  ! .True.
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
    write(6,'(A)')'0: removed several barriers in 0*f90 and GPCODE*f90'

    !------------------------------------------------------
    write(GP_print_unit, '(/3(A,1x,A,1x)//)') &
    '0: GPGACODE program version', trim(program_version), &
    '  branch:', trim( branch ) , &
    '  Last modified on:', trim( modification_date )
    !------------------------------------------------------


    ! read the control input from file  "GPCODE_cntl"

    call read_cntl_stuff( ierror  )

    !------------------------------------------------------

    ! open output units

    if( L_unit50_output )then
        open( unit_gp_out, file = 'unit50.txt', &
              form = 'unformatted', access='sequential', &
              status = 'unknown' )
    endif ! L_unit50_output


    if( L_GPSSE_log )then

        open( GPSSE_log_unit, file = 'GPSSE_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

        open( GPSSE_best_log_unit, file = 'GPSSE_best_log', &
              form = 'formatted', access='sequential', &
              status = 'unknown' )

    endif ! L_GP_log



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

!----------------------------------------------------------

! ierror > 0 if read_cntl_stuff has encountered a problem
! stop all processes in this case

if( ierror > 0 ) then

    call MPI_FINALIZE(ierr)
    stop 'ierror > 0'

endif  ! ierror > 0

!----------------------------------------------------------

if( n_input_vars > 0 )then

    ! read in the data number of points, number of vars
    
    if( myid == 0 )then
    
        write(6, '(/A)') '0: call read_input_data_size '
        call read_input_data_size( )
    
    endif !   myid == 0
    
    
    
    !---------------------------------------------------------------------
    
    ! broadcast number of data points, number of input variables
    
    call MPI_BCAST( n_input_data_points, 1,    &
                    MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
    
    call MPI_BCAST( n_input_vars, 1,    &
                    MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
    
    
    call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?
    
    !---------------------------------------------------------------------
    
    !write(6, '(/A,2(1x,I6))') '0: myid, n_input_data_points', &
    !                              myid, n_input_data_points
    !write(6, '(/A,2(1x,I6))') '0: myid, n_input_vars', myid, n_input_vars
    
    
    ! allocate input data names
    
    !write(6, '(/A)') '0: allocate input_data_names'
    allocate( input_data_names( 0:n_input_vars ) )
    
    
    ! allocate input data array
    
    !write(6, '(/A)') '0: allocate input_data_array'
    allocate( input_data_array( 0:n_input_vars, n_input_data_points) )
    
    
    call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?
    
    !---------------------------------------------------------------------
    
    ! read in the data into the arrays  input_data_names, input_data_array
    
    
    if( myid == 0 )then
    
        write(6, '(/A)') '0: call read_input_data'
    
        call read_input_data( )
    
        write(6,'(/A,2(1x,I6))') '0: n_input_data_points', &
                                     n_input_data_points
        write(6,'(/A,2(1x,I6))') '0: n_input_vars      =', n_input_vars
        write(6,'(/A,2(1x,I6))') '0: n_functions_input =', n_functions_input
    
    endif !   myid == 0


endif  !  n_input_vars > 0 

!---------------------------------------------------------------------

!sum_lmdif = 0.0d0

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
        write(6,'(/A,1x,I12)') &
              '0: user input random seed       clock = ', clock
    endif !   myid == 0

    seed = user_input_random_seed + &
              37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)
else

    CALL SYSTEM_CLOCK(COUNT=clock)

    if( myid == 0 )then
        write(6,'(/A,1x,I12)')&
              '0: random seed input clock = ', clock
    endif !   myid == 0

    seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)

endif ! user_input_random_seed > 0

CALL RANDOM_SEED(PUT = seed)


if( myid == 0 )then

    write(6,'(/A,1x,I12)') '0: n_seed ', n_seed
    write(6,'(/A)') '0: seed array '
    do  i = 1, n_seed

        write(6,'(I12,1x,I12)')  i, seed(i)

    enddo ! i
    write(6,'(A)') ' '

endif ! myid == 0

!------------------------------------------------------------------

! load table of 2**ilevel - 1  for RK process

call load_pow2_level(  )

!------------------------------------------------------------------

! set the scalar values for the model

! sets:
! n_levels
! n_functions
! n_CODE_equations
! n_trees
! n_nodes


call init_values( 0 )

n_Variables = n_CODE_equations

!---------------------------------------------------------------------

! n_inputs is used in deser*2 to point to input values in rk_data_array

n_inputs = n_input_vars

!------------------------------------------------------------------


if( myid == 0 )then

    write(6,'(/A,1x,I6)')     '0: n_code_equations ', n_code_equations
    write(6,'(A,1x,I6)')      '0: n_variables      ', n_variables
    write(6, '(/A,2(1x,I6))') '0: n_input_vars     ', n_input_vars
    write(6, '(A,2(1x,I6))')  '0: n_code_equations ', n_code_equations
    write(6, '(A,2(1x,I6)/)') '0: n_inputs         ', n_inputs

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



! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]

GP_Individual_Node_Parameters=0.0D0              ! Matrix Operation
GP_Individual_Node_Type=-9999                    ! Matrix Operation
GP_Population_Node_Parameters=0.0D0              ! Matrix Operation


GP_Adult_Population_Node_Type=-9999              ! Matrix Operation
GP_Child_Population_Node_Type=-9999              ! Matrix Operation

GP_minSSE_Individual_SSE = 1.0d99

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
! and parameter using the info from the setup file

! in set_answer_arrays, run the Runge-Kutta model only once with proc 0

! sets:

! GP_Node_Type_Answer
! GP_Node_Parameters_Answer
! Runge_Kutta_Solution
! Runge_Kutta_Node_Parameters
! Runge_Kutta_Node_Type

! GP_Node_Type_for_Plotting (if L_unit50_output true)


if( myid == 0 )then

    call set_answer_arrays( )

endif ! myid == 0

!------------------------------------------------------------------------

! then broadcast the R-K result: Runge_Kutta_Solution


if( myid == 0 )then    ! 20131209
    write(GP_print_unit,'(/A/)') &
          '0: time_step   Numerical_Code_Solution(time_step,1:n_CODE_equations)'
    do  i = 1, n_time_steps
        write(GP_print_unit,'(I6,2x,10(1x,E14.7))') &
              i, (Numerical_Code_Solution(i,jj), jj = 1,n_CODE_equations )
    enddo ! i
endif ! myid == 0


!message_len = ( n_time_steps + 1 ) * n_CODE_equations
!call MPI_BCAST( Runge_Kutta_Solution, message_len,    &
!                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


message_len = ( n_time_steps + 1 ) * n_CODE_equations
call MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



!Data_Array=Runge_Kutta_Solution          ! Matrix Operation
Data_Array=Numerical_CODE_Solution        ! Matrix Operation


!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

! compute the data variance with cpu 0 only, then broadcast results

! sets:
! Data_Variance
! Data_Variance_inv

if( myid == 0 )then    ! 20131209
    call comp_data_variance( )
endif ! myid == 0



message_len =  n_CODE_equations
call MPI_BCAST( Data_Variance_inv, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!--------------------------------------------------------------------------------


! put the desired model parameters in the array:  answer

answer = 0.0d0 ! set all to zero

n_parameters = 0

if( myid == 0 )then    ! 20131209
    write(GP_print_unit,'(/A,1x,I6/)') '0: n_code_equations ', n_code_equations
endif ! myid == 0

do  i_CODE_equation=1,n_CODE_equations
    n_parameters=n_parameters+1
    answer(n_parameters)=Numerical_CODE_Initial_Conditions(i_CODE_equation)
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

! calculate the generation interval for printing the list of children

GA_child_print_interval = n_GA_generations /  number_GA_child_prints

if( GA_child_print_interval == 0) then
    GA_child_print_interval = max( 1, n_GA_generations / 2 )
endif


GP_child_print_interval = n_GP_generations /  number_GP_child_prints

if( GP_child_print_interval == 0) then
    GP_child_print_interval = max( 1, n_GP_generations / 2 )
endif


call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

message_len = 1
call MPI_BCAST( GA_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


message_len = 1
call MPI_BCAST( GP_child_print_interval, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


!--------------------------------------------------------------------------------


if( myid == 0 )then

    call print_values2( answer )

    !-----------------------------------------------------------------------------

    ! this call calculates sse0,  the sse value obtained when the
    ! RK solution = 0 for all time steps

    ! note:  sse0 is only used by cpu 0 which does all fitness calculations

    call sse0_calc( )

    !---------------------------------------------------------------------------


    ! open more output files

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


    if( L_GP_all_summary  )then

        open( GP_summary_output_unit, file='GP_ALL_summary_file', &
              form = 'formatted', access = 'sequential', &
              status = 'unknown' )

    endif ! L_GP_all_summary  


    open( GP_best_summary_output_unit, file='GP_summary_file', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

    open( GP_minSSE_summary_output_unit, file='GP_minSSE_summary_file', &
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
    ! if the individual did not change on the last generation
    ! (so it needs no recalculation)


    Run_GP_Calculate_Fitness= .false.


    ! randomly create the initial tree arrays for each individual and
    ! send them all to GA_lmdif for parameter optimization on generation 1

    if( i_GP_Generation .eq. 1) then

        ! determines if the new GP child
        !has to be sent to GA_lmdif for parameter optimization

        Run_GP_Calculate_Fitness=.true.

        !---------------------------------------------------------------------------------

        if( myid == 0 )then

            write(GP_print_unit,'(/A,1x,I6)') &
              '0: call GP_Tree_Build        Generation =',i_GP_Generation

            ! initialize the GP_Adult_Population_Node_Type array with random trees

            !t1 = MPI_Wtime()


            call GP_Tree_Build


            !t2 = MPI_Wtime()


            !write(GP_print_unit,'(A,1x,E15.7)') &
            !  '0: time spent in GP_Tree_Build  = ', t2 - t1

            !! debug only >>>>>>>>>>>>>>>>
            !! set all GP tree models to the "truth" model
            !do  i_GP_individual = 1, n_GP_Individuals
            !    GP_Adult_Population_Node_Type(:,:,i_GP_individual) = &
            !    GP_Node_Type_Answer(:,:) ! debug only
            !enddo
            !! debug only <<<<<<<<<<<<<<<<<
            

            !! debug only >>>>>>>>>>>>>>>>
            ! fasham model


            i_gp_individual = 1
            GP_Adult_Population_Node_Type(:,:,: ) = -9999
            GP_Individual_Node_Type(:,:) =  -9999

        
!!            GP_Individual_Node_Type(9,1) =  -2
!!            !GP_Individual_Node_Parameters(8,1) = am  
!!            GP_Individual_Node_Type(8,1) =  0
!!            GP_Individual_Node_Type(5,1) =  -1
!!            GP_Individual_Node_Type(4,1) =  1
!!            !GP_Individual_Node_Parameters(3,1) = aN0
!!            GP_Individual_Node_Type(3,1) =  0
!!            GP_Individual_Node_Type(2,1) =  4
!!            GP_Individual_Node_Type(1,1) =  3
!!
!!
!!            GP_Individual_Node_Type(9,8) =  -2
!!            !GP_Individual_Node_Parameters(8,8) = am  
!!            GP_Individual_Node_Type(8,8) =  0
!!            GP_Individual_Node_Type(5,8) =  -1
!!            GP_Individual_Node_Type(4,8) =  1
!!            GP_Individual_Node_Type(3,8) =  -1
!!            GP_Individual_Node_Type(2,8) =  4
!!            GP_Individual_Node_Type(1,8) =  3
!!
!!
!!        
!!            GP_Individual_Node_Type(19,13) = -1 
!!            !GP_Individual_Node_Parameters(18,13) = aK1 
!!            GP_Individual_Node_Type(18,13) =  0
!!            GP_Individual_Node_Type(17,13) = -2 
!!            !GP_Individual_Node_Parameters(16,13) = phi 
!!            GP_Individual_Node_Type(16,13) =  0
!!            GP_Individual_Node_Type(9,13) =  6
!!            GP_Individual_Node_Type(8,13) =  9
!!            GP_Individual_Node_Type(5,13) = -6 
!!            GP_Individual_Node_Type(4,13) =  3
!!            GP_Individual_Node_Type(3,13) = -4 
!!            GP_Individual_Node_Type(2,13) =  3
!!            GP_Individual_Node_Type(1,13) =  3
!!
!!        
!!            GP_Individual_Node_Type(9,15) = -2
!!            !GP_Individual_Node_Parameters(8,15) = am  
!!            GP_Individual_Node_Type(8,15) =  0
!!            GP_Individual_Node_Type(5,15) = -1
!!            GP_Individual_Node_Type(4,15) =  1
!!            GP_Individual_Node_Type(3,15) = -2
!!            GP_Individual_Node_Type(2,15) =  4
!!            GP_Individual_Node_Type(1,15) =  3
!!
!!        
!!            GP_Individual_Node_Type(25,19) =  -3
!!            !GP_Individual_Node_Parameters(24,19) = eta 
!!            GP_Individual_Node_Type(24,19) =  0
!!            GP_Individual_Node_Type(17,19) =  -3
!!            !GP_Individual_Node_Parameters(16,19) = eta 
!!            GP_Individual_Node_Type(16,19) =  0
!!            GP_Individual_Node_Type(15,19) =  -3
!!            !GP_Individual_Node_Parameters(14,19) = aK4 
!!            GP_Individual_Node_Type(14,19) =  0
!!            GP_Individual_Node_Type(13,19) =  -2
!!            GP_Individual_Node_Type(12,19) =  3
!!            GP_Individual_Node_Type(11,19) =  -5
!!            !GP_Individual_Node_Parameters(10,19) = Vb  
!!            GP_Individual_Node_Type(10,19) =  0
!!            GP_Individual_Node_Type(9,19) =  -2
!!            GP_Individual_Node_Type(8,19) =  3
!!            GP_Individual_Node_Type(7,19) =  1
!!            GP_Individual_Node_Type(6,19) =  10
!!            GP_Individual_Node_Type(5,19) =  3
!!            GP_Individual_Node_Type(4,19) =  10
!!            GP_Individual_Node_Type(3,19) =  1
!!            GP_Individual_Node_Type(2,19) =  3
!!            GP_Individual_Node_Type(1,19) =  4
!!        
!!
!!        
!!            GP_Individual_Node_Type(9,20) =  -2
!!            !GP_Individual_Node_Parameters(8,20) = aK2 
!!            GP_Individual_Node_Type(8,20) =  0
!!            GP_Individual_Node_Type(5,20) =  -4
!!            GP_Individual_Node_Type(4,20) =  6
!!            GP_Individual_Node_Type(3,20) =  -6
!!            GP_Individual_Node_Type(2,20) =  3
!!            GP_Individual_Node_Type(1,20) =  3
!!
!!
!!        
!!            GP_Individual_Node_Type(9,22) =  -2
!!            !GP_Individual_Node_Parameters(8,22) = am  
!!            GP_Individual_Node_Type(8,22) =  0
!!            GP_Individual_Node_Type(5,22) =  -1
!!            GP_Individual_Node_Type(4,22) =  1
!!            GP_Individual_Node_Type(3,22) =  -3
!!            GP_Individual_Node_Type(2,22) =  4
!!            GP_Individual_Node_Type(1,22) =  3
!!        
!!        
!!            GP_Individual_Node_Type(49,26) =  -3
!!            !GP_Individual_Node_Parameters(48,26) = eta 
!!            GP_Individual_Node_Type(48,26) =  0
!!            GP_Individual_Node_Type(25,26) =  -2
!!            GP_Individual_Node_Type(24,26) =  3
!!            !GP_Individual_Node_Parameters(13,26) = aK4 
!!            GP_Individual_Node_Type(13,26) =  0
!!            GP_Individual_Node_Type(12,26) =  10
!!            GP_Individual_Node_Type(9,26) =  -5
!!            !GP_Individual_Node_Parameters(8,26) = Vb  
!!            GP_Individual_Node_Type(8,26) =  0
!!            GP_Individual_Node_Type(7,26) =  -3
!!            GP_Individual_Node_Type(6,26) =  1
!!            GP_Individual_Node_Type(5,26) =  -3
!!            GP_Individual_Node_Type(4,26) =  3
!!            GP_Individual_Node_Type(3,26) =  1
!!            GP_Individual_Node_Type(2,26) =  3
!!            GP_Individual_Node_Type(1,26) =  4
!!
!!
!!
!!            GP_Individual_Node_Type(17,29) =  -2
!!            !GP_Individual_Node_Parameters(16,29) = am 
!!            GP_Individual_Node_Type(16,29) =  0
!!            !GP_Individual_Node_Parameters(9,29) = V  
!!            GP_Individual_Node_Type(9,29) =  0
!!            GP_Individual_Node_Type(8,29) =  1
!!            GP_Individual_Node_Type(5,29) =  -1
!!            GP_Individual_Node_Type(4,29) =  1
!!            GP_Individual_Node_Type(3,29) =  -4
!!            GP_Individual_Node_Type(2,29) =  4
!!            GP_Individual_Node_Type(1,29) =  3
!!
!!        
!!            GP_Individual_Node_Type(3,32) =  -4
!!            !GP_Individual_Node_Parameters(9,32) = amu4
!!            GP_Individual_Node_Type(9,32) =  0
!!            GP_Individual_Node_Type(1,32) =  3
!!
!!        
!!            GP_Individual_Node_Type( 96,35) = 0                                                              
!!            GP_Individual_Node_Type( 98,35) = 0                                                              
!!            GP_Individual_Node_Type( 50,35) = 0                                                              
!!            GP_Individual_Node_Type( 13,35) = 0                                                              
!!            GP_Individual_Node_Type( 57,35) = 0                                                              
!!            GP_Individual_Node_Type( 29,35) = 0                                                              
!!            GP_Individual_Node_Type( 121,35) = 0                                                             
!!            GP_Individual_Node_Type( 61,35) = 0                                                              
!!            GP_Individual_Node_Type( 125,35) = 0                                                             
!!            GP_Individual_Node_Type( 63,35) = 0       
!!
!!!!            GP_Individual_Node_Type(17,35) =  -7
!!!!            !GP_Individual_Node_Parameters(16,35) = g         
!!!!            GP_Individual_Node_Type(16,35) =  0
!!!!            !GP_Individual_Node_Parameters(11,35) = 2.0d0 
!!!!            GP_Individual_Node_Type(11,35) =  0
!!!!            GP_Individual_Node_Type(10,35) =  -4
!!!!            !GP_Individual_Node_Parameters(9,35) = p3    
!!!!            GP_Individual_Node_Type(9,35) =  0
!!!!            GP_Individual_Node_Type(8,35) =  3
!!!!            GP_Individual_Node_Type(5,35) =  8
!!!!            GP_Individual_Node_Type(4,35) =  3
!!!!            GP_Individual_Node_Type(3,35) =  21
!!!!            GP_Individual_Node_Type(2,35) =  3
!!!!            GP_Individual_Node_Type(1,35) =  4
!!
!!            ! n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
!!        
!!            GP_Individual_Node_Type(17,35) =  -7
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n16 => GetParameterNode(g)
!!        
!!            GP_Individual_Node_Parameters(16,35) = g
!!            GP_Individual_Node_Type(16,35) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n11 => GetParameterNode(2.D+0)
!!        
!!            GP_Individual_Node_Parameters(11,35) = 2.0d0
!!            GP_Individual_Node_Type(11,35) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n10 => GetVariableNode(btmp(abs(SPECIES_DETRITUS)),SPECIES_DETRITUS)
!!        
!!            GP_Individual_Node_Type(10,35) =  -4
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n9 => GetParameterNode(p3)
!!        
!!            GP_Individual_Node_Parameters(9,35) = p3
!!            GP_Individual_Node_Type(9,35) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n8 => GetMathNode(Multiply, n16, n17)
!!        
!!            GP_Individual_Node_Type(8,35) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n5 => GetMathNode(Power, n10, n11)
!!        
!!            GP_Individual_Node_Type(5,35) =  8
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n4 => GetMathNode(Multiply, n8, n9)
!!        
!!            GP_Individual_Node_Type(4,35) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n3 => f_G_Lower()    ! make function 21
!!        
!!            GP_Individual_Node_Type(3,35) =  21
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n2 => GetMathNode(Multiply, n4, n5)
!!        
!!            GP_Individual_Node_Type(2,35) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n1 => GetMathNode(ProtectedDivide, n2, n3)
!!        
!!            GP_Individual_Node_Type(1,35) =  4
!!        
!!            !---------------------------------------------------------------------------
!!        
!!        
!!            GP_Individual_Node_Type(1,35)= 4             ! "[1] /"];
!!            GP_Individual_Node_Type(2,35)= 3             ! "[2] *"];
!!            GP_Individual_Node_Type(4,35)= 3             ! "[4] *"];
!!            GP_Individual_Node_Type(8,35)= 3             ! "[8] *"];
!!            GP_Individual_Node_Type(16,35)= 0
!!            GP_Individual_Node_Parameters(16,35)= g             ! "[16] (P)   1.00000000"];
!!            GP_Individual_Node_Type(17,35)= -7            ! ZOO  "[17] (V)       0.00"];
!!            GP_Individual_Node_Type(9,35)= 0 
!!            GP_Individual_Node_Parameters(9,35)= p3             ! "[9] (P)   1.00000000"];
!!            GP_Individual_Node_Type(5,35)= 8             !"[5] pow"];
!!            GP_Individual_Node_Type(10,35)= -4            ! DET  "[10] (V)       0.00"];
!!            GP_Individual_Node_Type(11,35)= 0
!!            GP_Individual_Node_Parameters(11,35) = 2.0d0              ! "[11] (P)   2.00000000"];
!!            GP_Individual_Node_Type(3,35)= 1             ! "[3] +"];
!!            GP_Individual_Node_Type(6,35)= 3             ! "[6] *"];
!!            GP_Individual_Node_Type(12,35)= 1             ! "[12] +"];
!!            GP_Individual_Node_Type(24,35)= 1             ! "[24] +"];
!!            GP_Individual_Node_Type(48,35)= 3             ! "[48] *"];
!!            GP_Individual_Node_Type(96,35)= 0
!!            GP_Individual_Node_Parameters(96,35)= p1            !   "[96] (P)   1.00000000"];
!!            GP_Individual_Node_Type(97,35)= -6             ! PHY  "[97] (V)       0.00"];
!!            GP_Individual_Node_Type(49,35)= 3             ! "[49] *"];
!!            GP_Individual_Node_Type(98,35)= 0
!!            GP_Individual_Node_Parameters(98,35)= p2             ! "[98] (P)   1.00000000"];
!!            GP_Individual_Node_Type(99,35)= -5             !  BACT  "[99] (V)       0.00"];
!!            GP_Individual_Node_Type(25,35)= 3             !  "[25] *"];
!!            GP_Individual_Node_Type(50,35)= 0
!!            GP_Individual_Node_Parameters(50,35)= p3            !  "[50] (P)   1.00000000"];
!!            GP_Individual_Node_Type(51,35)=  -4              ! DET "[51] (V)       0.00"];
!!            GP_Individual_Node_Type(13,35)= 0
!!            GP_Individual_Node_Parameters(13,35)= ak3            !"[13] (P)   1.00000000"];
!!            GP_Individual_Node_Type(7,35)= 1              ! "[7] +"];
!!            GP_Individual_Node_Type(14,35)= 3             ! "[14] *"];
!!            GP_Individual_Node_Type(28,35)= 8             ! "[28] pow"];
!!            GP_Individual_Node_Type(56,35)= -6             !  PHY "[56] (V)       0.00"];
!!            GP_Individual_Node_Type(57,35)= 0
!!            GP_Individual_Node_Parameters(57,35)=  2.0d0            ! "[57] (P)   2.00000000"];
!!            GP_Individual_Node_Type(29,35)= 0
!!            GP_Individual_Node_Parameters(29,35)= p1            ! "[29] (P)   1.00000000"];
!!            GP_Individual_Node_Type(15,35)= 1             ! "[15] +"];
!!            GP_Individual_Node_Type(30,35)= 3             ! "[30] *"];
!!            GP_Individual_Node_Type(60,35)= 8             ! "[60] pow"];
!!            GP_Individual_Node_Type(120,35)= -5            ! BACT  "[120] (V)       0.00"];
!!            GP_Individual_Node_Type(121,35)= 0
!!            GP_Individual_Node_Parameters(121,35)= 2.0d0           ! "[121] (P)   2.00000000"];
!!            GP_Individual_Node_Type(61,35)= 0
!!            GP_Individual_Node_Parameters(61,35)= p2            !  "[61] (P)   1.00000000"];
!!            GP_Individual_Node_Type(31,35)= 3              ! "[31] *"];
!!            GP_Individual_Node_Type(62,35)= 8              ! "[62] pow"];
!!            GP_Individual_Node_Type(124,35)= -4            !  DET "[124] (V)       0.00"];
!!            GP_Individual_Node_Type(125,35)= 0
!!            GP_Individual_Node_Parameters(125,35)= 2.0d0            ! "[125] (P)   2.00000000"];
!!            GP_Individual_Node_Type(63,35)= 0
!!            GP_Individual_Node_Parameters(63,35)=  p3             !  "[63] (P)   1.00000000"];
!!
!!
!!
!!
!!            GP_Individual_Node_Type(9,36) =  -2
!!            !GP_Individual_Node_Parameters(8,36) = am     
!!            GP_Individual_Node_Type(8,36) =  0
!!            GP_Individual_Node_Type(5,36) =  -1
!!            GP_Individual_Node_Type(4,36) =  1
!!            GP_Individual_Node_Type(3,36) =  -5
!!            GP_Individual_Node_Type(2,36) =  4
!!            GP_Individual_Node_Type(1,36) =  3
!!        
!!        
!!            GP_Individual_Node_Type(3,38) = -5
!!            !GP_Individual_Node_Parameters(2,38) = amu3  
!!            GP_Individual_Node_Type(2,38) =  0
!!            GP_Individual_Node_Type(1,38) =  3
!!        
!!
!!!!            GP_Individual_Node_Type( 96,42) = 0                                                              
!!!!            GP_Individual_Node_Type( 98,42) = 0                                                              
!!!!            GP_Individual_Node_Type( 50,42) = 0                                                              
!!!!            GP_Individual_Node_Type( 13,42) = 0                                                              
!!!!            GP_Individual_Node_Type( 57,42) = 0                                                              
!!!!            GP_Individual_Node_Type( 29,42) = 0                                                              
!!!!            GP_Individual_Node_Type( 121,42) = 0                                                             
!!!!            GP_Individual_Node_Type( 61,42) = 0                                                              
!!!!            GP_Individual_Node_Type( 125,42) = 0                                                             
!!!!            GP_Individual_Node_Type( 63,42) = 0       
!!!!
!!!!            GP_Individual_Node_Type(17,42) = -7
!!!!            !GP_Individual_Node_Parameters(16,42) = g
!!!!            GP_Individual_Node_Type(16,42) =  0
!!!!            !GP_Individual_Node_Parameters(11,42) = 2.0d0 
!!!!            GP_Individual_Node_Type(11,42) =  0
!!!!            GP_Individual_Node_Type(10,42) = -5  
!!!!            !GP_Individual_Node_Parameters(9,42) = p2    
!!!!            GP_Individual_Node_Type(9,42) =  0
!!!!            GP_Individual_Node_Type(8,42) = 3  
!!!!            GP_Individual_Node_Type(5,42) = 8  
!!!!            GP_Individual_Node_Type(4,42) = 3  
!!!!            GP_Individual_Node_Type(3,42) =  21
!!!!            GP_Individual_Node_Type(2,42) =  3
!!!!            GP_Individual_Node_Type(1,42) =  4
!!
!!            ! n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
!!        
!!            GP_Individual_Node_Type(17,42) = -7
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n16 => GetParameterNode(g)
!!        
!!            GP_Individual_Node_Parameters(16,42) = g
!!            GP_Individual_Node_Type(16,42) =  0
!!            !---------------------------------------------------------------------------
!!        
!!            ! n11 => GetParameterNode(2.D+0)
!!        
!!            GP_Individual_Node_Parameters(11,42) = 2.0d0
!!            GP_Individual_Node_Type(11,42) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n10 => GetVariableNode(btmp(abs(SPECIES_BACTERIA)),SPECIES_BACTERIA)
!!        
!!            GP_Individual_Node_Type(10,42) = -5
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n9 => GetParameterNode(p2)
!!        
!!            GP_Individual_Node_Parameters(9,42) = p2
!!            GP_Individual_Node_Type(9,42) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n8 => GetMathNode(Multiply, n16, n17)
!!        
!!            GP_Individual_Node_Type(8,42) = 3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n5 => GetMathNode(Power, n10, n11)
!!        
!!            GP_Individual_Node_Type(5,42) = 8
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n4 => GetMathNode(Multiply, n8, n9)
!!        
!!            GP_Individual_Node_Type(4,42) = 3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n3 => f_G_Lower()
!!        
!!            GP_Individual_Node_Type(3,42) =  21
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n2 => GetMathNode(Multiply, n4, n5)
!!        
!!            GP_Individual_Node_Type(2,42) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n1 => GetMathNode(ProtectedDivide, n2, n3)
!!        
!!            GP_Individual_Node_Type(1,42) =  4
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            GP_Individual_Node_Type(1,42)= 4            !  "[1] /"];
!!            GP_Individual_Node_Type(2,42)= 3            !  "[2] *"];
!!            GP_Individual_Node_Type(4,42)= 3            !  "[4] *"];
!!            GP_Individual_Node_Type(8,42)= 3            !  "[8] *"];
!!            GP_Individual_Node_Type(16,42)= 0
!!            GP_Individual_Node_Parameters(16,42)= g            ! "[16] (P)   1.00000000"];
!!            GP_Individual_Node_Type(17,42)= -7            ! ZOO  "[17] (V)       0.00"];
!!            GP_Individual_Node_Type(9,42)= 0
!!            GP_Individual_Node_Parameters(9,42)= p2            ! "[9] (P)   1.00000000"];
!!            GP_Individual_Node_Type(5,42)= 8            !  "[5] pow"];
!!            GP_Individual_Node_Type(10,42)= -5             ! BACT "[10] (V)       0.00"];
!!            GP_Individual_Node_Type(11,42)= 0
!!            GP_Individual_Node_Parameters(11,42)= 2.0d0            ! "[11] (P)   2.00000000"];
!!            GP_Individual_Node_Type(3,42)= 1            !  "[3] +"];
!!            GP_Individual_Node_Type(6,42)= 3            !  "[6] *"];
!!            GP_Individual_Node_Type(12,42)= 1            !  "[12] +"];
!!            GP_Individual_Node_Type(24,42)= 1            !  "[24] +"];
!!            GP_Individual_Node_Type(48,42)= 3            !  "[48] *"];
!!            GP_Individual_Node_Type(96,42)= 0
!!            GP_Individual_Node_Parameters(96,42)= p1           !  "[96] (P)   1.00000000"];
!!            GP_Individual_Node_Type(97,42)= -6            ! PHY  "[97] (V)       0.00"];
!!            GP_Individual_Node_Type(49,42)= 3            !  "[49] *"];
!!            GP_Individual_Node_Type(98,42)= 0
!!            GP_Individual_Node_Parameters(98,42)= p2            !  "[98] (P)   1.00000000"];
!!            GP_Individual_Node_Type(99,42)=  -5            !  BACT "[99] (V)       0.00"];
!!            GP_Individual_Node_Type(25,42)= 3            !  "[25] *"];
!!            GP_Individual_Node_Type(50,42)= 0
!!            GP_Individual_Node_Parameters(50,42)= p3            !  "[50] (P)   1.00000000"];
!!            GP_Individual_Node_Type(51,42)=  -4            !  DET "[51] (V)       0.00"];
!!            GP_Individual_Node_Type(13,42)= 0
!!            GP_Individual_Node_Parameters(13,42)= ak3            !  "[13] (P)   1.00000000"];
!!            GP_Individual_Node_Type(7,42)= 1            !  "[7] +"];
!!            GP_Individual_Node_Type(14,42)= 3            !  "[14] *"];
!!            GP_Individual_Node_Type(28,42)= 8            !  "[28] pow"];
!!            GP_Individual_Node_Type(56,42)= -5           ! BACT "[56] (V)       0.00"];
!!            GP_Individual_Node_Type(57,42)= 0
!!            GP_Individual_Node_Parameters(57,42)= 2.0d0            ! "[57] (P)   2.00000000"];
!!            GP_Individual_Node_Type(29,42)= 0
!!            GP_Individual_Node_Parameters(29,42)= p2            ! "[29] (P)   1.00000000"];
!!            GP_Individual_Node_Type(15,42)= 1            !  "[15] +"];
!!            GP_Individual_Node_Type(30,42)= 3            !  "[30] *"];
!!            GP_Individual_Node_Type(60,42)= 8            !  "[60] pow"];
!!            GP_Individual_Node_Type(120,42)= -5            ! BACT "[120] (V)       0.00"];
!!            GP_Individual_Node_Type(121,42)= 0
!!            GP_Individual_Node_Parameters(121,42)= 2.0d0            ! "[121] (P)   2.00000000"];
!!            GP_Individual_Node_Type(61,42)= 0
!!            GP_Individual_Node_Parameters(61,42)= 2.0d0            ! "[61] (P)   1.00000000"];
!!            GP_Individual_Node_Type(31,42)= 3            !  "[31] *"];
!!            GP_Individual_Node_Type(62,42)= 8            !  "[62] pow"];
!!            GP_Individual_Node_Type(124,42)= -4            ! DET "[124] (V)       0.00"];
!!            GP_Individual_Node_Type(125,42)= 0
!!            GP_Individual_Node_Parameters(125,42)= 2.0d0           ! "[125] (P)   2.00000000"];
!!            GP_Individual_Node_Type(63,42)= 0
!!            GP_Individual_Node_Parameters(63,42)= p3            ! "[63] (P)   1.00000000"];
!!
!!
!!        
!!        
!!            GP_Individual_Node_Type(9,43) =  -2
!!            !GP_Individual_Node_Parameters(8,43) = am    
!!            GP_Individual_Node_Type(8,43) =  0
!!            GP_Individual_Node_Type(5,43) =  -1
!!            GP_Individual_Node_Type(4,43) =  1
!!            GP_Individual_Node_Type(3,43) =  -6
!!            GP_Individual_Node_Type(2,43) =  4
!!            GP_Individual_Node_Type(1,43) =  3
!!
!!        
!!            GP_Individual_Node_Type(67,46) =  -1 
!!            !GP_Individual_Node_Parameters(66,46) = aK1    
!!            GP_Individual_Node_Type(66,46) =  0
!!            GP_Individual_Node_Type(65,46) =  -2 
!!            !GP_Individual_Node_Parameters(64,46) = phi    
!!            GP_Individual_Node_Type(64,46) =  0
!!            GP_Individual_Node_Type(35,46) =  -2 
!!            !GP_Individual_Node_Parameters(34,46) = ak2    
!!            GP_Individual_Node_Type(34,46) =  0
!!            GP_Individual_Node_Type(33,46) =  6
!!            GP_Individual_Node_Type(32,46) =  9
!!            GP_Individual_Node_Type(17,46) =  6
!!            GP_Individual_Node_Type(16,46) =  3
!!            !GP_Individual_Node_Parameters(9,46) = gamma1
!!            GP_Individual_Node_Type(9,46) =  0
!!            GP_Individual_Node_Type(8,46) =  1
!!            GP_Individual_Node_Type(3,46) =  -6
!!            GP_Individual_Node_Type(4,46) =  3
!!            GP_Individual_Node_Type(3,46) =  -4
!!            GP_Individual_Node_Type(2,46) =  3
!!            GP_Individual_Node_Type(1,46) =  3
!!        
!!        
!!            GP_Individual_Node_Type(3,47) = -6
!!            !GP_Individual_Node_Parameters(2,47) = amu1        
!!            GP_Individual_Node_Type(2,47) =  0
!!            GP_Individual_Node_Type(1,47) =  3
!!        
!!
!!!!            GP_Individual_Node_Type( 96,49) = 0                                                              
!!!!            GP_Individual_Node_Type( 98,49) = 0                                                              
!!!!            GP_Individual_Node_Type( 50,49) = 0                                                              
!!!!            GP_Individual_Node_Type( 13,49) = 0                                                              
!!!!            GP_Individual_Node_Type( 57,49) = 0                                                              
!!!!            GP_Individual_Node_Type( 29,49) = 0                                                              
!!!!            GP_Individual_Node_Type( 121,49) = 0       
!!!!            GP_Individual_Node_Type( 61,49) = 0                                                              
!!!!            GP_Individual_Node_Type( 125,49) = 0                                                             
!!!!            GP_Individual_Node_Type( 63,49) = 0         
!!!!
!!!!            GP_Individual_Node_Type(17,49) = -7 
!!!!            !GP_Individual_Node_Parameters(16,49) = g           
!!!!            GP_Individual_Node_Type(16,49) =  0
!!!!            !GP_Individual_Node_Parameters(11,49) = 2.0d0       
!!!!            GP_Individual_Node_Type(11,49) =  0
!!!!            GP_Individual_Node_Type(10,49) = -6 
!!!!            !GP_Individual_Node_Parameters(9,49) = p1          
!!!!            GP_Individual_Node_Type(9,49) =  0
!!!!            GP_Individual_Node_Type(8,49) = 3  
!!!!            GP_Individual_Node_Type(3,49) = 8  
!!!!            GP_Individual_Node_Type(4,49) = 3  
!!!!            GP_Individual_Node_Type(3,49) =  21
!!!!            GP_Individual_Node_Type(2,49) =  3
!!!!            GP_Individual_Node_Type(1,49) =  4
!!
!!            ! n17 => GetVariableNode(btmp(abs(SPECIES_ZOOPLANKTON)),SPECIES_ZOOPLANKTON)
!!        
!!            GP_Individual_Node_Type(17,49) = -7
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n16 => GetParameterNode(g)
!!        
!!            GP_Individual_Node_Parameters(16,49) = g
!!            GP_Individual_Node_Type(16,49) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n11 => GetParameterNode(2.D+0)
!!        
!!            GP_Individual_Node_Parameters(11,49) = 2.0d0
!!            GP_Individual_Node_Type(11,49) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n10 => GetVariableNode(btmp(abs(SPECIES_PHYTOPLANKTON)),SPECIES_PHYTOPLANKTON)
!!        
!!            GP_Individual_Node_Type(10,49) = -6
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n9 => GetParameterNode(p1)
!!        
!!            GP_Individual_Node_Parameters(9,49) = p1
!!            GP_Individual_Node_Type(9,49) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n8 => GetMathNode(Multiply, n16, n17)
!!        
!!            GP_Individual_Node_Type(8,49) = 3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n5 => GetMathNode(Power, n10, n11)
!!        
!!            GP_Individual_Node_Type(3,49) = 8
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n4 => GetMathNode(Multiply, n8, n9)
!!        
!!            GP_Individual_Node_Type(4,49) = 3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n3 => f_G_Lower()
!!        
!!            GP_Individual_Node_Type(3,49) =  21
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n2 => GetMathNode(Multiply, n4, n5)
!!        
!!            GP_Individual_Node_Type(2,49) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n1 => GetMathNode(ProtectedDivide, n2, n3)
!!        
!!            GP_Individual_Node_Type(1,49) =  4
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            GP_Individual_Node_Type(1,49)= 1            !  "[1] /"];
!!            GP_Individual_Node_Type(2,49)= 3            !  "[2] *"];
!!            GP_Individual_Node_Type(4,49)= 3            !  "[4] *"];
!!            GP_Individual_Node_Type(8,49)= 3            !  "[8] *"];
!!            GP_Individual_Node_Type(16,49)= 0
!!            GP_Individual_Node_Parameters(16,49)= g            ! "[16] (P)   1.00000000"];
!!            GP_Individual_Node_Type(17,49)= -7            ! ZOO "[17] (V)       0.00"];
!!            GP_Individual_Node_Type(9,49)= 0
!!            GP_Individual_Node_Parameters(9,49)= p1     ! "[9] (P)   1.00000000"];
!!            GP_Individual_Node_Type(5,49)= 8            !  "[5] pow"];
!!            GP_Individual_Node_Type(10,49)= -6            ! PHY "[10] (V)       0.00"];
!!            GP_Individual_Node_Type(11,49)= 0
!!            GP_Individual_Node_Parameters(11,49)= 2.0d0            ! "[11] (P)   2.00000000"];
!!            GP_Individual_Node_Type(3,49)= 1            !  "[3] +"];
!!            GP_Individual_Node_Type(6,49)= 3            !  "[6] *"];
!!            GP_Individual_Node_Type(12,49)= 1            !  "[12] +"];
!!            GP_Individual_Node_Type(24,49)= 1            !  "[24] +"];
!!            GP_Individual_Node_Type(48,49)= 3            !  "[48] *"];
!!            GP_Individual_Node_Type(96,49)= 0
!!            GP_Individual_Node_Parameters(96,49)= p1            ! "[96] (P)   1.00000000"];
!!            GP_Individual_Node_Type(97,49)= -6            ! PHY  "[97] (V)       0.00"];
!!            GP_Individual_Node_Type(49,49)= 3            !  "[49] *"];
!!            GP_Individual_Node_Type(98,49)= 0
!!            GP_Individual_Node_Parameters(98,49)= p2            ! "[98] (P)   1.00000000"];
!!            GP_Individual_Node_Type(99,49)= -5            ! BACT "[99] (V)       0.00"];
!!            GP_Individual_Node_Type(25,49)= 3            !  "[25] *"];
!!            GP_Individual_Node_Type(50,49)= 0
!!            GP_Individual_Node_Parameters(50,49)= p3            ! "[50] (P)   1.00000000"];
!!            GP_Individual_Node_Type(51,49)= -4           ! DET "[51] (V)       0.00"];
!!            GP_Individual_Node_Type(13,49)= 0
!!            GP_Individual_Node_Parameters(13,49)= ak3            !  "[13] (P)   1.00000000"];
!!            GP_Individual_Node_Type(7,49)= 1            !  "[7] +"];
!!            GP_Individual_Node_Type(14,49)= 3            !  "[14] *"];
!!            GP_Individual_Node_Type(28,49)= 8            !  "[28] pow"];
!!            GP_Individual_Node_Type(56,49)= -5            ! BACT "[56] (V)       0.00"];
!!            GP_Individual_Node_Type(57,49)= 0
!!            GP_Individual_Node_Parameters(57,49)= 2.0d0            ! "[57] (P)   2.00000000"];
!!            GP_Individual_Node_Type(29,49)= 0
!!            GP_Individual_Node_Parameters(29,49)= p2            ! "[29] (P)   1.00000000"];
!!            GP_Individual_Node_Type(15,49)= 1            !  "[15] +"];
!!            GP_Individual_Node_Type(30,49)= 3            !  "[30] *"];
!!            GP_Individual_Node_Type(60,49)= 8            !  "[60] pow"];
!!            GP_Individual_Node_Type(120,49)=  -5            ! BACT "[120] (V)       0.00"];
!!            GP_Individual_Node_Type(121,49)= 0
!!            GP_Individual_Node_Parameters(121,49)=2.0d0            ! "[121] (P)   2.00000000"];
!!            GP_Individual_Node_Type(61,49)= 0
!!            GP_Individual_Node_Parameters(61,49)= p2           ! "[61] (P)   1.00000000"];
!!            GP_Individual_Node_Type(31,49)= 3            !  "[31] *"];
!!            GP_Individual_Node_Type(62,49)= 8            !  "[62] pow"];
!!            GP_Individual_Node_Type(124,49)= -4            ! DET "[124] (V)       0.00"];
!!            GP_Individual_Node_Type(125,49)= 0
!!            GP_Individual_Node_Parameters(125,49)= 2.0d0            ! "[125] (P)   2.00000000"];
!!            GP_Individual_Node_Type(63,49)= 0
!!            GP_Individual_Node_Parameters(63,49)= p3            ! "[63] (P)   1.00000000"];
!!
!!
!!
!!        
!!
!!            !GP_Individual_Node_Parameters(13,50) = omega       
!!            GP_Individual_Node_Type(13,50) =  0
!!            !GP_Individual_Node_Parameters(12,50) = amu5        
!!            GP_Individual_Node_Type(12,50) =  0
!!            GP_Individual_Node_Type(9,50) =  -7
!!            GP_Individual_Node_Type(8,50) =  -3
!!            GP_Individual_Node_Type(7,50) =  -7
!!            GP_Individual_Node_Type(6,50) =  3
!!            GP_Individual_Node_Type(5,50) =  -1
!!            GP_Individual_Node_Type(4,50) =  3
!!            GP_Individual_Node_Type(3,50) =  3
!!            GP_Individual_Node_Type(2,50) =  4
!!            GP_Individual_Node_Type(1,50) =  1
!!        
!!
!!            !GP_Individual_Node_Parameters(17,52) = omega
!!            GP_Individual_Node_Type(17,52) =  0
!!            !GP_Individual_Node_Parameters(16,52) = 1.0d0
!!            GP_Individual_Node_Type(16,52) =  0
!!            !GP_Individual_Node_Parameters(13,52) = amu2        
!!            GP_Individual_Node_Type(13,52) =  0
!!            !GP_Individual_Node_Parameters(12,52) = epsilon     
!!            GP_Individual_Node_Type(12,52) =  0
!!            !GP_Individual_Node_Parameters(9,52) = amu5        
!!            GP_Individual_Node_Type(9,52) =  0
!!            GP_Individual_Node_Type(8,52) =  3
!!            GP_Individual_Node_Type(7,52) = -7
!!            GP_Individual_Node_Type(6,52) =  3
!!            GP_Individual_Node_Type(5,52) = -7
!!            GP_Individual_Node_Type(4,52) =  3
!!            GP_Individual_Node_Type(3,52) =  3
!!            GP_Individual_Node_Type(2,52) =  3
!!            GP_Individual_Node_Type(1,52) =  1
!!        
!!
!!        
!!
!!        
!!            !GP_Individual_Node_Parameters(9,53) = epsilon     
!!            GP_Individual_Node_Type(9,53) =  0
!!            !GP_Individual_Node_Parameters(8,53) = 1.0d0       
!!            GP_Individual_Node_Type(8,53) =  0
!!            !GP_Individual_Node_Parameters(5,53) = amu2        
!!            GP_Individual_Node_Type(5,53) =  0
!!            GP_Individual_Node_Type(4,53) = 2
!!            GP_Individual_Node_Type(3,53) = -7
!!            GP_Individual_Node_Type(2,53) =  3
!!            GP_Individual_Node_Type(1,53) =  3
!!
!!
!!
!!!!            GP_Individual_Node_Type( 128,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 65,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 67,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 544,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 546,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 274,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 69,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 281,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 141,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 569,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 285,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 573,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 287,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 160,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 81,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 83,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 672,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 674,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 338,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 85,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 345,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 173,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 697,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 349,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 701,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 351,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 96,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 49,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 51,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 416,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 418,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 210,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 53,54) = 0                                                              
!!!!            GP_Individual_Node_Type( 217,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 109,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 441,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 221,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 445,54) = 0                                                             
!!!!            GP_Individual_Node_Type( 223,54) = 0             
!!!!            !GP_Individual_Node_Parameters(23,54) = beta2
!!!!            GP_Individual_Node_Type(23,54) =  0
!!!!            !GP_Individual_Node_Parameters(22,54) = 1.0d0 
!!!!            GP_Individual_Node_Type(22,54) =  0
!!!!            !GP_Individual_Node_Parameters(19,54) = beta1 
!!!!            GP_Individual_Node_Type(19,54) =  0
!!!!            !GP_Individual_Node_Parameters(18,54) = 1.0d0 
!!!!            GP_Individual_Node_Type(18,54) =  0
!!!!            !GP_Individual_Node_Parameters(15,54) = beta3
!!!!            GP_Individual_Node_Type(15,54) =  0
!!!!            !GP_Individual_Node_Parameters(14,54) = 1.0d0       
!!!!            GP_Individual_Node_Type(14,54) =  0
!!!!            GP_Individual_Node_Type(11,54) =  2
!!!!            GP_Individual_Node_Type(10,54) =  23
!!!!            GP_Individual_Node_Type(9,54) =  2
!!!!            GP_Individual_Node_Type(8,54) =  22
!!!!            GP_Individual_Node_Type(7,54) =  2
!!!!            GP_Individual_Node_Type(6,54) =  24
!!!!            GP_Individual_Node_Type(5,54) =  3
!!!!            GP_Individual_Node_Type(4,54) =  3
!!!!            GP_Individual_Node_Type(3,54) =  3
!!!!            GP_Individual_Node_Type(2,54) =  1
!!!!            GP_Individual_Node_Type(1,54) =  1
!!
!!            ! n23 => GetParameterNode(beta2)
!!        
!!            GP_Individual_Node_Parameters(23,54) = beta2
!!            GP_Individual_Node_Type(23,54) =  0
!!            !---------------------------------------------------------------------------
!!        
!!            ! n22 => GetParameterNode(1.D+0)
!!        
!!            GP_Individual_Node_Parameters(22,54) = 1.0d0
!!            GP_Individual_Node_Type(22,54) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n19 => GetParameterNode(beta1)
!!        
!!            GP_Individual_Node_Parameters(19,54) = beta1
!!            GP_Individual_Node_Type(19,54) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n18 => GetParameterNode(1.D+0)
!!        
!!            GP_Individual_Node_Parameters(18,54) = 1.0d0
!!            GP_Individual_Node_Type(18,54) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n15 => GetParameterNode(beta3)
!!        
!!            GP_Individual_Node_Parameters(15,54) = beta3
!!            GP_Individual_Node_Type(15,54) =  0
!!            !---------------------------------------------------------------------------
!!        
!!            ! n14 => GetParameterNode(1.D+0)
!!        
!!            GP_Individual_Node_Parameters(14,54) = 1.0d0
!!            GP_Individual_Node_Type(14,54) =  0
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n11 => GetMathNode(Subtract, n22, n23)
!!        
!!            GP_Individual_Node_Type(11,54) =  2
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n10 => f_G2()
!!        
!!            GP_Individual_Node_Type(10,54) =  23
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n9 => GetMathNode(Subtract, n18, n19)
!!        
!!            GP_Individual_Node_Type(9,54) =  2
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n8 => f_G1()
!!        
!!            GP_Individual_Node_Type(8,54) =  22
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n7 => GetMathNode(Subtract, n14, n15)
!!        
!!            GP_Individual_Node_Type(7,54) =  2
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n6 => f_G3()
!!        
!!            GP_Individual_Node_Type(6,54) =  24
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n5 => GetMathNode(Multiply, n10, n11)
!!        
!!            GP_Individual_Node_Type(5,54) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n4 => GetMathNode(Multiply, n8, n9)
!!        
!!            GP_Individual_Node_Type(4,54) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n3 => GetMathNode(Multiply, n6, n7)
!!        
!!            GP_Individual_Node_Type(3,54) =  3
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n2 => GetMathNode(Add, n4, n5)
!!        
!!            GP_Individual_Node_Type(2,54) =  1
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            ! n1 => GetMathNode(Add, n2, n3)
!!        
!!            GP_Individual_Node_Type(1,54) =  1
!!        
!!            !---------------------------------------------------------------------------
!!        
!!            GP_Individual_Node_Type(1,54)= 1            !  "[1] +"];
!!            GP_Individual_Node_Type(2,54)= 1            !  "[2] +"];
!!            GP_Individual_Node_Type(4,54)= 3            !  "[4] *"];
!!            GP_Individual_Node_Type(8,54)= 4            !  "[8] /"];
!!            GP_Individual_Node_Type(16,54)= 3           !   "[16] *"];
!!            GP_Individual_Node_Type(32,54)= 3           !   "[32] *"];
!!            GP_Individual_Node_Type(64,54)= 3           !   "[64] *"];
!!            GP_Individual_Node_Type(128,54)= 0
!!            GP_Individual_Node_Parameters(128,54)= g            !  "[128] (P)   1.00000000"];
!!            GP_Individual_Node_Type(129,54)= -7            ! ZOO  "[129] (V)       0.00"];
!!            GP_Individual_Node_Type(65,54)= 0
!!            GP_Individual_Node_Parameters(65,54)= p1            !  "[65] (P)   1.00000000"];
!!            GP_Individual_Node_Type(33,54)= 8            !  "[33] pow"];
!!            GP_Individual_Node_Type(66,54)= -6            ! PHY "[66] (V)       0.00"];
!!            GP_Individual_Node_Type(67,54)= 0
!!            GP_Individual_Node_Parameters(67,54)=  2.0d0             ! "[67] (P)   2.00000000"];
!!            GP_Individual_Node_Type(17,54)= 1            !  "[17] +"];                               ****1
!!            GP_Individual_Node_Type(34,54)= 3           !   "[34] *"];                               ****1
!!            GP_Individual_Node_Type(68,54)= 1            !  "[68] +"];                               ****1
!!            GP_Individual_Node_Type(136,54)= 1            !  "[136] +"];                             ****1
!!            GP_Individual_Node_Type(272,54)= 3           !   "[272] *"];                             ****1
!!            GP_Individual_Node_Type(544,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(544,54)= p1           ! "[544] (P)   1.00000000"];         ****1
!!            GP_Individual_Node_Type(545,54)= -6           ! PHY "[545] (V)       0.00"];             ****1
!!            GP_Individual_Node_Type(273,54)= 3           !   "[273] *"];                             ****1
!!            GP_Individual_Node_Type(546,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(546,54)= p2           ! "[546] (P)   1.00000000"];         ****1
!!            GP_Individual_Node_Type(547,54)= -5            ! BACT "[547] (V)       0.00"];           ****1
!!            GP_Individual_Node_Type(137,54)= 3           !   "[137] *"];                             ****1
!!            GP_Individual_Node_Type(274,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(274,54)= p3            !  "[274] (P)   1.00000000"];       ****1
!!            GP_Individual_Node_Type(275,54)= -4            ! DET "[275] (V)       0.00"];            ****1
!!            GP_Individual_Node_Type(69,54)= 0                                         !               ****1
!!            GP_Individual_Node_Parameters(69,54)= ak3            !  "[69] (P)   1.00000000"];        ****1
!!            GP_Individual_Node_Type(35,54)= 1            !  "[35] +"];                               ****1
!!            GP_Individual_Node_Type(70,54)= 3           !   "[70] *"];                               ****1
!!            GP_Individual_Node_Type(140,54)= 8            !  "[140] pow"];                           ****1
!!            GP_Individual_Node_Type(280,54)= -6            ! PHY "[280] (V)       0.00"];            ****1
!!            GP_Individual_Node_Type(281,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(281,54)= 2.0d0           !  "[281] (P)   2.00000000"];     ****1
!!            GP_Individual_Node_Type(141,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(141,54)= p1            !  "[141] (P)   1.00000000"];       ****1
!!            GP_Individual_Node_Type(71,54)= 1            !  "[71] +"];                               ****1
!!            GP_Individual_Node_Type(142,54)= 3           !   "[142] *"];                             ****1
!!            GP_Individual_Node_Type(284,54)= 8            !  "[284] pow"];                           ****1
!!            GP_Individual_Node_Type(568,54)= -5            ! BACT "[568] (V)       0.00"];           ****1
!!            GP_Individual_Node_Type(569,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(569,54)= 2.0d0           !  "[569] (P)   2.00000000"];     ****1
!!            GP_Individual_Node_Type(285,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(285,54)= p2            !  "[285] (P)   1.00000000"];       ****1
!!            GP_Individual_Node_Type(143,54)= 3           !   "[143] *"];                             ****1
!!            GP_Individual_Node_Type(286,54)= 8            !  "[286] pow"];                           ****1
!!            GP_Individual_Node_Type(572,54)= -4            ! DET  "[572] (V)       0.00"];           ****1
!!            GP_Individual_Node_Type(573,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(573,54)= 2.0d0            !  "[573] (P)   2.00000000"];    ****1
!!            GP_Individual_Node_Type(287,54)= 0                                        !               ****1
!!            GP_Individual_Node_Parameters(287,54)= p3            !  "[287] (P)   1.00000000"];       ****1
!!            GP_Individual_Node_Type(9,54)= 2            !  "[9] -"];
!!            GP_Individual_Node_Type(18,54)= 0
!!            GP_Individual_Node_Parameters(18,54)= 1.0d0           ! "[18] (P)   1.00000000"];
!!            GP_Individual_Node_Type(19,54)= 0
!!            GP_Individual_Node_Parameters(19,54)= beta1            ! "[19] (P)   0.75000000"];
!!            GP_Individual_Node_Type(5,54)= 3           !   "[5] *"];
!!            GP_Individual_Node_Type(10,54)= 4            !  "[10] /"];
!!            GP_Individual_Node_Type(20,54)= 3           !   "[20] *"];
!!            GP_Individual_Node_Type(40,54)= 3           !   "[40] *"];
!!            GP_Individual_Node_Type(80,54)= 3           !   "[80] *"];
!!            GP_Individual_Node_Type(160,54)= 0
!!            GP_Individual_Node_Parameters(160,54)= g            !  "[160] (P)   1.00000000"];
!!            GP_Individual_Node_Type(161,54)= -7            ! ZOO "[161] (V)       0.00"];
!!            GP_Individual_Node_Type(81,54)= 0
!!            GP_Individual_Node_Parameters(81,54)= p2            ! "[81] (P)   1.00000000"];
!!            GP_Individual_Node_Type(41,54)= 8            !  "[41] pow"];
!!            GP_Individual_Node_Type(82,54)= -5            ! BACT "[82] (V)       0.00"];
!!            GP_Individual_Node_Type(83,54)= 0
!!            GP_Individual_Node_Parameters(83,54)=  2.0d0            !  "[83] (P)   2.00000000"];
!!            GP_Individual_Node_Type(21,54)= 1            !  "[21] +"];                               ****2
!!            GP_Individual_Node_Type(42,54)= 3           !   "[42] *"];                               ****2
!!            GP_Individual_Node_Type(84,54)= 1            !  "[84] +"];                               ****2
!!            GP_Individual_Node_Type(168,54)= 1            !  "[168] +"];                             ****2
!!            GP_Individual_Node_Type(336,54)= 3           !   "[336] *"];                             ****2
!!            GP_Individual_Node_Type(672,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(672,54)= p1           !  "[672] (P)   1.00000000"];        ****2
!!            GP_Individual_Node_Type(673,54)= -6            ! PHY "[673] (V)       0.00"];            ****2
!!            GP_Individual_Node_Type(337,54)= 3           !   "[337] *"];                             ****2
!!            GP_Individual_Node_Type(674,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(674,54)= p2            ! "[674] (P)   1.00000000"];        ****2
!!            GP_Individual_Node_Type(675,54)= -5           ! BACT "[675] (V)       0.00"];            ****2
!!            GP_Individual_Node_Type(169,54)= 3           !   "[169] *"];                             ****2
!!            GP_Individual_Node_Type(338,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(338,54)= p3           !  "[338] (P)   1.00000000"];        ****2
!!            GP_Individual_Node_Type(339,54)= -4            ! DET "[339] (V)       0.00"];            ****2
!!            GP_Individual_Node_Type(85,54)= 0                                           !             ****2
!!            GP_Individual_Node_Parameters(85,54)= ak3           ! "[85] (P)   1.00000000"];          ****2
!!            GP_Individual_Node_Type(43,54)= 1            !  "[43] +"];                               ****2
!!            GP_Individual_Node_Type(86,54)= 3           !   "[86] *"];                               ****2
!!            GP_Individual_Node_Type(172,54)= 8            !  "[172] pow"];                           ****2
!!            GP_Individual_Node_Type(344,54)= -6           ! PHY "[344] (V)       0.00"];             ****2
!!            GP_Individual_Node_Type(345,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(345,54)= 2.0d0           ! "[345] (P)   2.00000000"];      ****2
!!            GP_Individual_Node_Type(173,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(173,54)= p1           !  "[173] (P)   1.00000000"];        ****2
!!            GP_Individual_Node_Type(87,54)= 1            !  "[87] +"];                               ****2
!!            GP_Individual_Node_Type(174,54)= 3           !   "[174] *"];                             ****2
!!            GP_Individual_Node_Type(348,54)= 8            !  "[348] pow"];                           ****2
!!            GP_Individual_Node_Type(696,54)= -5            ! BACT "[696] (V)       0.00"];           ****2
!!            GP_Individual_Node_Type(697,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(697,54)= 2.0d0           !  "[697] (P)   2.00000000"];     ****2
!!            GP_Individual_Node_Type(349,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(349,54)= p2            !  "[349] (P)   1.00000000"];       ****2
!!            GP_Individual_Node_Type(175,54)= 3           !   "[175] *"];                             ****2
!!            GP_Individual_Node_Type(350,54)= 8            !  "[350] pow"];                           ****2
!!            GP_Individual_Node_Type(700,54)= -4            ! DET "[700] (V)       0.00"];            ****2
!!            GP_Individual_Node_Type(701,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(701,54)= 2.0d0            !  "[701] (P)   2.00000000"];    ****2
!!            GP_Individual_Node_Type(351,54)= 0                                          !             ****2
!!            GP_Individual_Node_Parameters(351,54)= p3            !  "[351] (P)   1.00000000"];       ****2
!!            GP_Individual_Node_Type(11,54)= 2            !  "[11] -"];
!!            GP_Individual_Node_Type(22,54)= 0
!!            GP_Individual_Node_Parameters(22,54)=  1.0d0            ! "[22] (P)   1.00000000"];
!!            GP_Individual_Node_Type(23,54)= 0
!!            GP_Individual_Node_Parameters(23,54)= beta2             ! "[23] (P)   0.75000000"];
!!            GP_Individual_Node_Type(3,54)= 3           !   "[3] *"];
!!            GP_Individual_Node_Type(6,54)= 4            !  "[6] /"];
!!            GP_Individual_Node_Type(12,54)= 3           !   "[12] *"];
!!            GP_Individual_Node_Type(24,54)= 3           !   "[24] *"];
!!            GP_Individual_Node_Type(48,54)= 3           !   "[48] *"];
!!            GP_Individual_Node_Type(96,54)= 0
!!            GP_Individual_Node_Parameters(96,54)= g            !  "[96] (P)   1.00000000"];
!!            GP_Individual_Node_Type(97,54)= -7            ! ZOO  "[97] (V)       0.00"];
!!            GP_Individual_Node_Type(49,54)= 0
!!            GP_Individual_Node_Parameters(49,54)= p3            !  "[49] (P)   1.00000000"];
!!            GP_Individual_Node_Type(25,54)= 8            !  "[25] pow"];
!!            GP_Individual_Node_Type(50,54)= -4            ! DET "[50] (V)       0.00"];
!!            GP_Individual_Node_Type(51,54)= 0
!!            GP_Individual_Node_Parameters(51,54)= 2.0d0           ! "[51] (P)   2.00000000"];
!!            GP_Individual_Node_Type(13,54)= 1            !  "[13] +"];                               ****3
!!            GP_Individual_Node_Type(26,54)= 3           !   "[26] *"];                               ****3
!!            GP_Individual_Node_Type(52,54)= 1            !  "[52] +"];                               ****3
!!            GP_Individual_Node_Type(104,54)= 1            !  "[104] +"];                             ****3
!!            GP_Individual_Node_Type(208,54)= 3           !   "[208] *"];                             ****3
!!            GP_Individual_Node_Type(416,54)= 0                                        !               ****3
!!            GP_Individual_Node_Parameters(416,54)= p1            !  "[416] (P)   1.00000000"];       ****3
!!            GP_Individual_Node_Type(417,54)=  -6            !  PHY  "[417] (V)       0.00"];         ****3
!!            GP_Individual_Node_Type(209,54)= 3           !   "[209] *"];                             ****3
!!            GP_Individual_Node_Type(418,54)= 0                                        !               ****3
!!            GP_Individual_Node_Parameters(418,54)= p2            !  "[418] (P)   1.00000000"];       ****3
!!            GP_Individual_Node_Type(419,54)= -5            ! BACT "[419] (V)       0.00"];           ****3
!!            GP_Individual_Node_Type(105,54)= 3           !   "[105] *"];                             ****3
!!            GP_Individual_Node_Type(210,54)= 0                                        !               ****3
!!            GP_Individual_Node_Parameters(210,54)= p3            !  "[210] (P)   1.00000000"];       ****3
!!            GP_Individual_Node_Type(211,54)= -4            ! DET  "[211] (V)       0.00"];           ****3
!!            GP_Individual_Node_Type(53,54)= 0                                         !               ****3
!!            GP_Individual_Node_Parameters(53,54)= ak3            ! "[53] (P)   1.00000000"];         ****3
!!            GP_Individual_Node_Type(210,54)= 0                                        !               ****3
!!            GP_Individual_Node_Parameters(210,54)= p3           !  "[210] (P)   1.00000000"];        ****3
!!            GP_Individual_Node_Type(27,54)= 1            !  "[27] +"];                               ****3
!!            GP_Individual_Node_Type(54,54)= 3           !   "[54] *"];                               ****3
!!            GP_Individual_Node_Type(108,54)= 8            !  "[108] pow"];                           ****3
!!            GP_Individual_Node_Type(216,54)= -6            ! PHY "[216] (V)       0.00"];            ****3
!!            GP_Individual_Node_Type(217,54)= 0                                        !               ****3
!!            GP_Individual_Node_Parameters(217,54)= 2.0d0           !  "[217] (P)   2.00000000"];     ****3
!!            GP_Individual_Node_Type(109,54)= 0                                        !               ****3
!!            GP_Individual_Node_Parameters(109,54)=  p1            ! "[109] (P)   1.00000000"];       ****3
!!            GP_Individual_Node_Type(55,54)= 1            !  "[55] +"];                               ****3
!!            GP_Individual_Node_Type(110,54)= 3           !   "[110] *"];                             ****3
!!            GP_Individual_Node_Type(220,54)= 8            !  "[220] pow"];                           ****3
!!            GP_Individual_Node_Type(440,54)= -5            ! BACT  "[440] (V)       0.00"];          ****3
!!            GP_Individual_Node_Type(441,54)= 0                      !                                 ****3
!!            GP_Individual_Node_Parameters(441,54)= 2.0d0           ! "[441] (P)   2.00000000"];      ****3
!!            GP_Individual_Node_Type(221,54)= 0                      !                                 ****3
!!            GP_Individual_Node_Parameters(221,54)= p2           !  "[221] (P)   1.00000000"];        ****3
!!            GP_Individual_Node_Type(111,54)= 3           !   "[111] *"];                             ****3
!!            GP_Individual_Node_Type(222,54)= 8            !  "[222] pow"];                           ****3
!!            GP_Individual_Node_Type(444,54)= -4            ! DET "[444] (V)       0.00"];            ****3
!!            GP_Individual_Node_Type(445,54)= 0                      !                                 ****3
!!            GP_Individual_Node_Parameters(445,54)= 2.0d0            !  "[445] (P)   2.00000000"];    ****3
!!            GP_Individual_Node_Type(223,54)= 0                      !                                 ****3
!!            GP_Individual_Node_Parameters(223,54)= p3            !  "[223] (P)   1.00000000"];       ****3
!!            GP_Individual_Node_Type(7,54)= 2            !  "[7] -"];
!!            GP_Individual_Node_Type(14,54)= 0
!!            GP_Individual_Node_Parameters(14,54)=  1.0d0            !  "[14] (P)   1.00000000"];
!!            GP_Individual_Node_Type(15,54)= 0
!!            GP_Individual_Node_Parameters(15,54)=  beta3            ! "[15] (P)   0.75000000"];


    

            do  i_GP_individual = 1, n_GP_Individuals
                GP_Adult_Population_Node_Type(:,:,i_GP_individual) = &
                      GP_Individual_Node_Type(:,:)
                GP_Population_Node_parameters(:,:,i_GP_individual) = &
                      GP_Individual_Node_parameters(:,:) 
            enddo

            !! debug only <<<<<<<<<<<<<<<<<

        endif ! myid == 0


        !---------------------------------------------------------------------------------

        ! broadcast GP_Adult_Population_Node_Type

        if( myid == 0 )then
            write(GP_print_unit,'(A,1x,I6)') &
              '0: broadcast  GP_Adult_Population_Node_Type Generation = ',i_GP_Generation
        endif ! myid == 0



        message_len = n_GP_Individuals * n_Nodes * n_Trees
        call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                        MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6)') &
        !      '0: aft broadcast  GP_Adult_Population_Node_Type  ierr = ',ierr
        !endif ! myid == 0



        !if( myid == 0 )then
        !    call print_debug_integer_node_tree( GP_print_unit, &
        !         'aft bcast 528  GP_Adult_Population_Node_Type ', &
        !         GP_Adult_Population_Node_Type )
        !endif ! myid == 0


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

            !----------------------------------------------------------------------------------

            !tree_descrip =  ' GP_Adult trees before call selection routines '
            !call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !                  GP_Adult_Population_Node_Type, &
            !                  trim( tree_descrip )  )

            if( i_GP_generation == 1                                  .or. &
                mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
                i_GP_generation == n_GP_generations                          )then

                write(GP_print_unit,'(//A)') '0:3 before modifications'                             
                write(GP_print_unit,'(A)')&
                      '0:3 i_GP_gen i_GP_indiv    GP_Child_Indiv_SSE&
                      &   GP_Child_Indiv_SSE/SSE0'

                do  i_GP_individual = 1, n_GP_individuals
                    write(GP_print_unit,'(2(1x,I10), 2(1x, E20.10))') &
                               i_GP_generation, i_GP_individual, &
                               GP_Child_Individual_SSE(i_GP_Individual), &
                               GP_Child_Individual_SSE(i_GP_Individual)/SSE0
                enddo ! i_GP_individual

            endif ! i_GP_generation == 1 .or. ...

            !----------------------------------------------------------------------------------

            !   i) Carry out "GP Fitness-Proportionate Reproduction"

            !      randomly replaces values of individuals in child arrays
            !      with values from the adult arrays of fitter individuals

            !   uses:
            !   GP_Integrated_Population_Ranked_Fitness
            !   GP_Adult_Population_Node_Type
            !   GP_Adult_Population_SSE
            !   GP_Population_Node_Parameters
            !   GP_Population_Initial_Conditions

            !   sets:
            !   GP_Child_Population_Node_Type
            !   GP_Child_Individual_SSE
            !   GP_Population_Node_Parameters
            !   GP_Population_Initial_Conditions


            if( n_GP_Asexual_Reproductions .gt. 0 )then

                write(GP_print_unit,'(/A,1x,I6)') &
                      '0: call GP_Fit_Prop_Asexual_Repro &
                      &n_GP_Asexual_Reproductions =', n_GP_Asexual_Reproductions


                call GP_Fitness_Proportionate_Asexual_Reproduction


                ! don't replace adult with child until
                ! after all 3 modification routines have been called

                !write(GP_print_unit,'(/A)') &
                !      '0:aft  call GP_Fitness_Proportionate_Asexual_Reproduction '
                !write(GP_print_unit,'(A/)') &
                !      '0: Adult NT, SSE NOT REPLACED with Child NT, SSE'

                !tree_descrip=' GP_Child trees after call to &
                !              &GP_Fitness_Proportionate_Asexual_Reproduction'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Child_Population_Node_Type, &
                !                    trim( tree_descrip )  )

                !tree_descrip=' GP_Adult trees after call to &
                !              &GP_Fitness_Proportionate_Asexual_Reproduction'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Adult_Population_Node_Type, &
                !                    trim( tree_descrip )  )

                !call print_debug_real_node_tree( GP_print_unit, &
                !         'aft GP_Fit  print GP_population_node_parameters ', &
                !         GP_population_node_parameters )


            endif !  n_GP_Asexual_Reproductions .gt. 0

            if( i_GP_generation == 1                                  .or. &
                mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
                i_GP_generation == n_GP_generations                          )then

                write(GP_print_unit,'(//A)') '0:3 after Asexual_Reproduction' 
                write(GP_print_unit,'(A)')&
                      '0:3 i_GP_gen i_GP_indiv    GP_Child_Indiv_SSE&
                      &   GP_Child_Indiv_SSE/SSE0'

                do  i_GP_individual = 1, n_GP_individuals
                    write(GP_print_unit,'(2(1x,I10), 2(1x, E20.10))') &
                               i_GP_generation, i_GP_individual, &
                               GP_Child_Individual_SSE(i_GP_Individual), &
                               GP_Child_Individual_SSE(i_GP_Individual)/SSE0
                enddo ! i_GP_individual

            endif ! i_GP_generation == 1 .or. ...

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

                write(GP_print_unit,'(/A,1x,I6)') &
                      '0: call GP_Tour_Style_Sexual_Repro n_GP_Crossovers =', &
                                                          n_GP_Crossovers

                call GP_Tournament_Style_Sexual_Reproduction


                !write(GP_print_unit,'(/A)') &
                !      '0: aft  call GP_Tournament_Style_Sexual_Reproduction '

                !tree_descrip = ' GP_Child trees after call to &
                !                  &GP_Tournament_Style_Sexual_Reproduction'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Child_Population_Node_Type, &
                !                    trim( tree_descrip )  )

                !call print_debug_real_node_tree( GP_print_unit, &
                !         'aft GP_Tour print GP_population_node_parameters ', &
                !         GP_population_node_parameters )

            endif !  n_GP_Crossovers .gt. 0

            if( i_GP_generation == 1                                  .or. &
                mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
                i_GP_generation == n_GP_generations                          )then

                write(GP_print_unit,'(//A)') '0:3 after Crossover ' 
                write(GP_print_unit,'(A)')&
                      '0:3 i_GP_gen i_GP_indiv    GP_Child_Indiv_SSE&
                      &   GP_Child_Indiv_SSE/SSE0'

                do  i_GP_individual = 1, n_GP_individuals
                    write(GP_print_unit,'(2(1x,I10), 2(1x, E20.10))') &
                               i_GP_generation, i_GP_individual, &
                               GP_Child_Individual_SSE(i_GP_Individual), &
                               GP_Child_Individual_SSE(i_GP_Individual)/SSE0
                enddo ! i_GP_individual

            endif ! i_GP_generation == 1 .or. ...


            !----------------------------------------------------------------------------------

            !   iii) Carry out "GP Parameter Mutation" Operations

            ! uses:
            !  GP_Adult_Population_Node_Type

            ! sets:
            !  GP_Child_Population_Node_Type

            !  Run_GP_Calculate_Fitness  ( to true for modified individuals )


            if( n_GP_Mutations .gt. 0 )then

                write(GP_print_unit,'(/A,13x,I6, 1x, E15.7)')&
                      '0: call GP_Mutations n_GP_Mutations, prob_no_elite', &
                                            n_GP_Mutations, prob_no_elite

                !tree_descrip =  ' GP_Adult trees BEFORE call to GP_Mutations'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Adult_Population_Node_Type, &
                !                    trim( tree_descrip )  )

                call GP_Mutations

                !tree_descrip =  ' GP_Child trees after call to GP_Mutations'
                !call print_trees( i_GP_generation, 1, n_GP_individuals, &
                !                    GP_Child_Population_Node_Type, &
                !                    trim( tree_descrip )  )

            endif !  n_GP_Mutations .gt. 0



            !tree_descrip =  ' GP_Adult trees after call to GP_Mutations'
            !call print_trees( i_GP_generation, 1, n_GP_individuals, &
            !                    GP_Adult_Population_Node_Type, &
            !                    trim( tree_descrip )  )


            !---------------------------------------------------------------------------

            !   Move over any newly created children into the adult arrays

            GP_Adult_Population_Node_Type = GP_Child_Population_Node_Type
            GP_Adult_Population_SSE       = GP_Child_Individual_SSE


            !write(GP_print_unit,'(/A)')&
            !      '0:aft  move Child_Node_Type and SSE to Adult'

            !write(GP_print_unit,'(/A/(10(3x,L1)))')&
            !      '0: Run_GP_Calculate_Fitness ', Run_GP_Calculate_Fitness

            !call print_debug_integer_node_tree( GP_print_unit, &
            !         'aft mutation print GP_Adult_Population_Node_Type ', &
            !         GP_Adult_Population_Node_Type )

            ! call print_debug_real_node_tree( GP_print_unit, &
            !          'aft mutation print GP_population_node_parameters ', &
            !          GP_population_node_parameters )

            if( i_GP_generation == 1                                  .or. &
                mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
                i_GP_generation == n_GP_generations                          )then


                write(GP_print_unit,'(//A)') '0:3 after Mutations' 

                write(GP_print_unit,'(A)')&
                      '0:3 i_GP_gen i_GP_indiv    GP_Child_Indiv_SSE&
                      &   GP_Child_Indiv_SSE/SSE0'

                do  i_GP_individual = 1, n_GP_individuals
                    write(GP_print_unit,'(2(1x,I10), 2(1x, E20.10))') &
                               i_GP_generation, i_GP_individual, &
                               GP_Child_Individual_SSE(i_GP_Individual), &
                               GP_Child_Individual_SSE(i_GP_Individual)/SSE0
                enddo ! i_GP_individual

                write(GP_print_unit,'(/A/(10(3x,L1)))')&
                      '0: Run_GP_Calculate_Fitness ', Run_GP_Calculate_Fitness

            endif ! i_GP_generation == 1 .or. ...
            !---------------------------------------------------------------------------

            ! calculate the diversity index for each individual for generations > 1

            !write(GP_print_unit,'(/A)')&
            !      '0: call GP_calc_diversity_index '

            !call GP_calc_diversity_index( n_GP_individuals, &
            !                              GP_Child_Population_Node_Type, &
            !                              i_diversity, i_gp_generation )

            !write(GP_print_unit,'(/A)')&
            !      '0: aft call GP_calc_diversity_index '

            !---------------------------------------------------------------------------

        endif ! myid == 0


        !------------------------------------------------------------------------------------


        ! broadcast:
        ! GP_Child_Population_Node_Type
        ! GP_Adult_Population_Node_Type
        ! GP_Child_Individual_SSE
        ! GP_Integrated_Population_Ranked_Fitness
        ! GP_Population_Ranked_Fitness
        ! Run_GP_Calculate_Fitness

        !if( myid == 0 )then
        !    write(GP_print_unit,'(/A)')&
        !          '0: call bcast2 '
        !endif ! myid == 0

        call bcast2()


    endif ! i_GP_Generation .eq. 1


    !-----------------------------------------------------------------------------------------

    ! print trees before tree clean

    if( myid == 0 )then
        !if( i_GP_generation == 1                                  .or. &
        !    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
        !    i_GP_generation == n_GP_generations                          )then
            tree_descrip =  ' trees before call to GP_Clean_Tree_Nodes'
            call print_trees( i_GP_generation, 1, n_GP_individuals, &
                 GP_Adult_Population_Node_Type, trim( tree_descrip )  )
        !endif ! i_GP_generation == 1
    endif !  myid == 0


    !-----------------------------------------------------------------------------------------

    ! GP_Clean_Tree_Nodes sweeps through the GP_Adult_Population_Node_Type array
    ! to replace function nodes that have both terminals set as parameters
    ! and to set the replaced node to a parameter itself

    if( myid == 0 )then

        !write(GP_print_unit,'(/A,1x,I6/)') &
        !      '0: call GP_Clean_Tree_Nodes  Generation =', i_GP_Generation

        !tree_descrip =  ' trees BEFORE call to GP_Clean_Tree_Nodes'
        !call print_trees( i_GP_generation, 1, n_GP_individuals, &
        !         GP_Adult_Population_Node_Type, trim( tree_descrip )  )


        call GP_Clean_Tree_Nodes


        !tree_descrip =  ' trees after call to GP_Clean_Tree_Nodes'
        !call print_trees( i_GP_generation, 1, n_GP_individuals, &
        !         GP_Adult_Population_Node_Type, trim( tree_descrip )  )

    endif ! myid == 0


    ! broadcast GP_Adult_Population_Node_Type changed by GP_Clean_Tree_Nodes


    message_len = n_GP_Individuals * n_Nodes * n_Trees
    call MPI_BCAST( GP_Adult_Population_Node_Type, message_len,    &
                    MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )


    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A,2(1x,I6))') &
    !      '0: aft broadcast  GP_Adult_Pop_Node_Type  Generation, ierr = ', &
    !                                            i_GP_Generation, ierr
    !endif ! myid == 0

    GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type


    !-----------------------------------------------------------------------------------------

    ! print trees after call to GP_Clean_Tree_Nodes

    !if( myid == 0 )then
    !    if( i_GP_generation == 1                                  .or. &
    !        mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
    !        i_GP_generation == n_GP_generations                          )then
    !
    !        tree_descrip =  ' trees after call to GP_Clean_Tree_Nodes'
    !        call print_trees( i_GP_generation, 1, n_GP_individuals, &
    !                   GP_Adult_Population_Node_Type, trim( tree_descrip )  )
    !
    !        write(GP_print_unit,'(/A, 1x, I6/)') &
    !          '0: after call to GP_Clean_Tree_Nodes i_GP_generation =',i_GP_generation
    !        print node type information for each GP individual
    !           call print_gp_node_type_parm( )
    !
    !    endif ! i_GP_generation == 1 .or. ...
    !
    !endif !  myid == 0

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


        !------------------------------------------------------------------------------
        ! calculate how many parameters total to fit for the specific individual CODE
        ! and save this number in GP_Individual_N_GP_param(i_GP_individual)

        n_GP_Parameters = n_code_equations

        if( myid == 0 )then
            write(GP_print_unit,'(/A,2(1x,I6)/)')&
             '0: before loop n_code_equations,  n_GP_params',&
                             n_code_equations,  n_GP_parameters
        endif !  myid == 0

        do  i_Tree=1,n_Trees

            do  i_Node=1,n_Nodes

                if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0) then

                    n_GP_Parameters = n_GP_Parameters+1

                    !if( myid == 0 )then
                    !    write(GP_print_unit,'(A,5(1x,I6))')&
                    !    '0: i_GP_indiv, i_tree, i_node, GP_Adult_Pop_Node_Type, n_GP_params',&
                    !        i_GP_individual, i_tree, i_node, &
                    !        GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual), &
                    !        n_GP_parameters
                    !endif !  myid == 0

                endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

                if( myid == 0 )then
                    if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) > -9999 )then
                        write(GP_print_unit,'(A,5(1x,I6))')&
                        '0: i_GP_indiv, i_tree, i_node, GP_Adult_Pop_Node_Type, n_GP_params',&
                            i_GP_individual, i_tree, i_node, &
                            GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual), &
                            n_GP_parameters
                    endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) > -9999
                endif !  myid == 0

            enddo ! i_node

        enddo ! i_tree


        if( myid == 0 )then
            write(GP_print_unit,'(/A,2(1x,I6)/)')&
             '0: after loop n_code_equations,  n_GP_params',&
                            n_code_equations,  n_GP_parameters
        endif !  myid == 0

        !------------------------------------------------------------------------

        ! if n_GP_parameters > n_maximum_number_parameters, this is an error
        ! for now, just print a warning and set n_GP_parameters = n_max...

        if( n_GP_parameters > n_maximum_number_parameters )then

            if( myid == 0 )then
                write(GP_print_unit,'(/A,5(1x,I10)/)')&
                  '0:WARNING  myid, i_GP_generation, i_GP_Individual, &
                      &n_GP_parameters, n_maximum_number_parameters', &
                      myid, i_GP_generation, i_GP_Individual, &
                       n_GP_parameters, n_maximum_number_parameters
            endif !  myid == 0

            call MPI_FINALIZE(ierr)
            stop 'too many parms'

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
            if( myid == 0 )then
                write(GP_print_unit,'(/A)')&
                '0:----------------------------------------------------------------------'
                write(GP_print_unit,'(A,2(1x,I6),3x,L1)') &
                      '0: i_GP_Gen, i_GP_indiv, Run_GP_Calculate_Fitness', &
                          i_GP_Generation, i_GP_individual, &
                                   Run_GP_Calculate_Fitness(i_GP_Individual)
                write(GP_print_unit,'(A)')&
                '0:----------------------------------------------------------------------'
            endif !  myid == 0
            !-----------------------------------------------------------------------------------

            if( myid == 0 )then
                write(GP_print_unit,'(/A,4(1x,I6))') &
                  '0: i_GP_individual, n_trees, n_nodes, n_GP_parameters ', &
                      i_GP_individual, n_trees, n_nodes, n_GP_parameters
            endif !  myid == 0

            !-----------------------------------------------------------------------------------

            do  i_Tree=1,n_Trees
                do  i_Node=1,n_Nodes

                    GP_Individual_Node_Type(i_Node,i_Tree) = &
                       GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

                            !!debug only  !GP_Individual_Node_Parameters(i_Node,i_Tree) = 0.0d0

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
                        write(GP_print_unit,'(A,1x,I6,5x,A,2(1x,I6))')&
                              '0: for i_GP_Individual', i_GP_Individual, &
                              'number of parameters, variables =', &
                                    n_GP_parameters, n_GP_vars
                    endif !  myid == 0


                    !-------------------------------------------------------------------

                    ! cycle the i_GP_individual loop if there are no GP parameters
                    ! or if n_GP_parameters <=  n_code_equations

                    if( n_GP_parameters == 0 .or. &
                        n_GP_parameters > n_maximum_number_parameters .or.  &
                        n_GP_parameters <=  n_code_equations                 ) then   ! new jjm 20130814

                        !if( myid == 0 )then
                        !    write(GP_print_unit,'(A,1x,I6)')&
                        !          '0: skipping this i_GP_Individual --&
                        !          &  the number of parameters is ', n_GP_parameters
                        !endif !  myid == 0


                        ! set SSE values for this rejected individual so that
                        ! its fitness will be very small

                        GP_Child_Individual_SSE(i_GP_Individual) = 1.0D13

                        GP_Adult_Individual_SSE(i_GP_Individual) = 1.0D13
                        GP_Adult_Population_SSE(i_GP_Individual) = 1.0D13


                        if( myid == 0 )then

                            write(GP_print_unit,'(/A/9x,A/A,2(1x,I5), 1x, E15.7)')&
                                  '0: rejected for n_GP_parameters <=  n_code_equations',&
                                  'or for n_GP_parameters >   n_maximum_number_parameters',&
                                   'i_GP_gen,i_GP_indiv,GP_Child_Pop_SSE(i_GP_Indiv)  ', &
                                       i_GP_generation, i_GP_individual, &
                                       GP_Child_Individual_SSE(i_GP_Individual)

                        endif !  myid == 0

                        cycle gp_ind_loop

                    endif ! n_GP_parameters == 0

                    !-------------------------------------------------------------------


                    if( myid == 0 )then

                        !write(GP_print_unit,'(/A,2(1x,I6))') &
                        ! '0: call GPCODE_GA_lmdif_Param_Opt         i_GP_Gen, i_GP_indiv', &
                        !                           i_GP_Generation, i_GP_individual

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
                    !     '0: aft call GPCODE_GA_lmdif_Parameter_Optimization routine'
                    !    flush(GP_print_unit)
                    !endif ! myid == 0


                    !--------------------------------------------------------------------------------

                    GP_population_fitness(i_GP_individual) = individual_fitness


                    if( myid == 0 )then

                        !write(GP_print_unit,'(/A,2(1x,I6), 1x, E15.7)')&
                        !      '0:1 i_GP_gen, i_GP_indiv, GP_pop_fit(i_GP_indiv) ', &
                        !           i_GP_generation, i_GP_individual, &
                        !           GP_population_fitness(i_GP_individual)

                        ! print side-by-side comparisons of
                        ! starting values and values from optimization

                        !write(GP_print_unit,'(/A/)') &
                        !  '0:           truth value           GP_individual_initial_conditions'

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

                                write(GP_print_unit,'(A,2(1x,I5), 2(1x, E15.7))')&
                                '0:3 i_GP_gen,  i_GP_indiv,  GP_Child_Indiv_SSE,  GP_pop_fit', &
                                     i_GP_generation, i_GP_individual, &
                                     GP_Child_Individual_SSE(i_GP_Individual), &
                                     GP_population_fitness(i_GP_individual)

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



            if( myid == 0 )then

                !----------------------------------------------------------------------------

                ! print the node parameters (if there are any)

                !write(GP_print_unit,'(/A/)') &
                !       '0:  tree  node  Runge_Kutta_Node_Params&
                !                  &   GP_population_node_params'
                !do  i_tree=1,n_trees
                !    do  i_node=1,n_nodes
                !        if( GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0 ) then
                !            write(GP_print_unit,'(2(1x,I6), 1x, E20.10, 4x, E20.10)') &
                !             i_tree, i_node, &
                !             Runge_Kutta_Node_Parameters(i_node,i_tree), &
                !             GP_population_node_parameters(i_node,i_tree,i_GP_individual)
                !        endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

                !    enddo ! i_node
                !enddo  ! i_tree

                !---------------------------------------------------------------------------------


                !write(GP_print_unit,'(/A)')  &
                !      '0: after loading GP_Pop arrays with GP_indiv array values '


            endif !  myid == 0

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
                !    write(GP_print_unit,'(/A)') &
                !      '0:------------------------------------------&
                !       &-----------------------------'
                !    write(GP_print_unit,'(A,2(1x,I6))') &
                !    '0: call summary_GP_indiv i_GP_generation, i_GP_individual ', &
                !                              i_GP_generation, i_GP_individual
                !
                !endif ! i_GP_generation == 1 .or. ...
    
                call summary_GP_indiv( i_GP_generation, i_GP_individual, 0 )

            endif ! L_GP_all_summary 

        endif !  myid == 0


    enddo  gp_ind_loop    !   i_GP_individual


    !-------------------------------------------------------------------------------------

    if( myid == 0 )then

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)') &
            '==================================&
            &==============================================='
            write(GP_print_unit,'(A,1x,I6)') &
            '0: aft indiv loop and before GP_para_lmdif_process   i_GP_generation =',&
                                                                  i_GP_Generation
            write(GP_print_unit,'(A/)') &
            '=================================&
            &================================================'

            write(GP_print_unit, '(A )') &
                 '0:i_GP_Indiv  GP_Indiv_N_param    GP_Adult_Pop_SSE   &
                  &  GP_Child_Indiv_SSE   GP_Child_Indiv_SSE/SSE0'

            do  i_GP_individual = 1, n_GP_individuals
                write(GP_print_unit, '(5x,I6,6x,I6,6x,3(1x, E20.10) )') &
                      i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual), &
                                        GP_Adult_Population_SSE(i_GP_Individual), &
                                        GP_Child_Individual_SSE(i_GP_Individual), &
                                        GP_Child_Individual_SSE(i_GP_Individual)/SSE0
            enddo

        endif ! i_GP_generation == 1 .or. ...

    endif ! myid == 0

    !-------------------------------------------------------------------------------------


    !  call GP_para_lmdif_process routine to run lmdif
    !  in parallel on all the GP individuals

    !  GP_para_lmdif_process returns arrays to be used in GP_calc_fitness

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

    GP_para_flag = .TRUE.

    max_n_gp_params = maxval( GP_Individual_N_GP_param )

    if( myid == 0 )then
        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            max_n_gp_params = maxval( GP_Individual_N_GP_param )

            write(GP_print_unit,'(/A)')&
            '0:-----------------------------------------------------------------'
            write(GP_print_unit,'(A,2(1x,I6))') &
            '0: call GP_para_lmdif_process i_GP_generation, max_n_gp_params', &
                                           i_GP_Generation, max_n_gp_params
            write(GP_print_unit,'(A/)')&
            '0:-----------------------------------------------------------------'

            flush( GP_print_unit )

        endif ! i_GP_generation == 1 .or. ...

    endif ! myid == 0


    !---------------------------------------------------------------

    call GP_para_lmdif_process( i_GP_generation, max_n_gp_params  )

    !---------------------------------------------------------------

    !GP_para_flag = .FALSE.

    !-------------------------------------------------------------------------------------

    if( myid == 0 )then

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then


            write(GP_print_unit,'(/A)') &
            '====================================&
            &============================================='
            write(GP_print_unit,'(A,1x,I6)') &
            '0: aft indiv loop and AFTER  GP_para_lmdif_process   i_GP_generation =',&
                                                                  i_GP_Generation
            write(GP_print_unit,'(A/)') &
            '====================================&
            &============================================='

            write(GP_print_unit, '(/A )') &
                 '0:i_GP_Indiv  GP_Indiv_N_param    GP_Adult_Pop_SSE    &
                  & GP_Child_Indiv_SSE   GP_Child_Indiv_SSE/SSE0'

            do  i_GP_individual = 1, n_GP_individuals
                write(GP_print_unit, '(5x,I6,6x,I6,6x,3(1x, E20.10) )') &
                 i_GP_Individual,  GP_Individual_N_GP_param(i_GP_individual), &
                                   GP_Adult_Population_SSE(i_GP_Individual), &
                                   GP_Child_Individual_SSE(i_GP_Individual), &
                                   GP_Child_Individual_SSE(i_GP_Individual)/SSE0
            enddo

        endif ! i_GP_generation == 1 .or. ...

    endif ! myid == 0

    !-------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------


    ! do fitness calculations for this GP generation

    if( myid == 0 )then

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)')&
            '0:#################################################################'
            write(GP_print_unit,'(A,1x,I6)') &
            '0: call GP_calc_fitness i_GP_generation =', &
                                     i_GP_Generation
            write(GP_print_unit,'(A/)')&
            '0:#################################################################'

        endif ! i_GP_generation == 1 .or. ...

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  debug
        ! print the node parameters (if there are any)
        !write(GP_print_unit,'(/A/)') &
        !   '0:  i_GP_indiv tree  node  GP_population_node_params'
        !call print_debug_real_node_tree( GP_print_unit, &
        !              'GP_population_node_parameters before GP_calc_fitness', &
        !               GP_population_node_parameters )
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  debug

        !----------------------------------------------------------------------------------

        ! do fitness calculations for this GP generation

        call GP_calc_fitness( i_GP_generation, output_array, &
                              i_GP_best_parent, nop )

        !----------------------------------------------------------------------------------

        !write(GP_print_unit,'(/A,1x,I6/)') &
        !  '0:1  n_GP_individuals =', n_GP_individuals

        if( i_GP_generation == 1                                  .or. &
            mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
            i_GP_generation == n_GP_generations                          )then

            write(GP_print_unit,'(/A)')&
              '0:#################################################################'
            write(GP_print_unit,'(A,3(1x,I6))') &
              '0: aft call GP_calc_fitness n_GP_indiv, i_GP_gen, i_GP_best_parent =', &
                                           n_GP_individuals, i_GP_Generation, i_GP_best_parent
            write(GP_print_unit,'(A)')&
              '0:#################################################################'

        endif ! i_GP_generation == 1 .or. ...

        !----------------------------------------------------------------------------------

        ! whenever the SSE for the best parent is less than GP_minSSE_Individual_SSE, 
        ! load the GP_minSSE* arrays with the corresponding arrays for the best parent

        ! thus, at the end of the run, the GP_minSSE* arrays will have the values of the 
        ! best individual over all generations

        ! this is needed only if the prob_no_elite parameter > 0, so that it is possible that
        ! the best individual on the final generation is not the best found in the run 


        if( GP_Child_Individual_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE  ) then

            GP_minSSE_Individual_SSE = GP_Child_Individual_SSE(i_GP_best_parent)

            GP_minSSE_Individual_Initial_Conditions(1:n_CODE_equations)  = &
                      GP_Population_Initial_Conditions(1:n_CODE_equations, i_GP_best_parent)

            do  i_tree=1,n_trees
                do  i_node=1,n_nodes

                    if( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_best_parent) == 0 ) then

                        GP_minSSE_Individual_Node_Parameters(i_node,i_tree) = &
                               GP_population_node_parameters(i_node,i_tree,i_GP_best_parent) 

                    endif ! GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_best_parent) == 0

                enddo ! i_node
            enddo  ! i_tree

            GP_minSSE_Individual_Node_Type(1:n_nodes,1:n_trees)  = &
                    GP_Adult_population_Node_Type(1:n_Nodes, 1:n_Trees, i_GP_best_parent)

            GP_minSSE_Individual_N_GP_param =  GP_Individual_N_GP_param(i_GP_best_parent)

            GP_minSSE_Individual =  i_GP_best_parent
            GP_minSSE_generation =  i_GP_generation                                      

        endif  !  GP_Child_Individual_SSE(i_GP_best_parent) <  GP_minSSE_Individual_SSE



    endif ! myid == 0

    !---------------------------------------------------------------------------

    ! broadcast results of GP_calc_fitness:

    !  GP_Adult_Individual_SSE
    !  GP_population_node_parameters
    !  GP_Population_Ranked_Fitness
    !  GP_Integrated_Population_Ranked_Fitness

    call MPI_BARRIER( MPI_COMM_WORLD, ierr )


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
    !---------------------------------------------------------------------------

    write(GP_print_unit,'(A,1x,I5,2(1x,E15.7)/)') &
    '0: i_GP_best_parent, GP_child_indiv_sse(), SSE/SSE0', &
        i_GP_best_parent, GP_child_individual_sse(i_GP_best_parent), &
                          GP_child_individual_sse(i_GP_best_parent)/SSE0


    ! GP_select_best_RK_lmdif_result runs lmdif on the best parent

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


    write(GP_print_unit,'(A/)') &
       '0: DO NOT call GP_select_best_RK_lmdif_result to run lmdif for best parent'

    !!!call GP_select_best_RK_lmdif_result( i_GP_best_parent, output_array , nop )


    !---------------------------------------------------------------------------


endif ! myid == 0

flush(GP_print_unit)


!----------------------------------------------------------------------------------------

! plot results

!if( myid == 0 )then
!    Lplot = .true.
!    if( Lplot ) call plot_results(Runge_Kutta_Solution,n_time_steps, n_CODE_equations )
!endif ! myid == 0


if( myid == 0 )then

    !------------------------------------------------------------------------------------

    max_n_gp_params = maxval( GP_Individual_N_GP_param )

    write(GP_print_unit,'(/A,3(1x,I5))') &
    '0: call print_time_series  i_GP_best_parent, max_n_gp_params, nop ', &
        i_GP_best_parent, max_n_gp_params, nop

    call print_time_series( i_GP_best_parent, output_array, nop )


    !------------------------------------------------------------------------------------

    ! this prints a summary of the initial conditions,
    ! parameters,  and node types for the individual with the minimum SSE
    ! and writes the tree to the summary file


        write(GP_print_unit,'(/A)') &
          '0:------------------------------------------&
           &-----------------------------'
        write(GP_print_unit,'(A,2(1x,I6))') &
        '0: call summary_GP_minSSE_indiv GP_minSSE_generation, GP_minSSE_Individual ', &
                                         GP_minSSE_generation, GP_minSSE_Individual

    call summary_GP_minSSE_indiv( GP_minSSE_generation, GP_minSSE_Individual )


    call print_time_series_minSSE( )

    !------------------------------------------------------------------------------------

endif ! myid == 0


!---------------------------------------------------------------------------------------

!  close output units

if( myid == 0 )then


    if( L_unit50_output )then
        close( unit_gp_out )
    endif ! L_unit50_output


    if( L_GP_log )then
        close( GP_log_unit )
    endif ! L_GP_log


    if( L_GPSSE_log )then
        close( GPSSE_log_unit )
        close( GPSSE_best_log_unit )
    endif ! L_GPSSE_log


    if( L_GA_log )then
        close( GA_log_unit )
    endif ! L_GA_log


    if( L_GA_output_parameters )then
        close( GA_output_unit )
    endif ! L_GA_output_parameters


    if( L_GP_output_parameters )then
        close( GP_output_unit )
    endif ! L_GP_output_parameters


    if( L_GP_all_summary )then
        close( GP_summary_output_unit )
    endif ! L_GP_all_summary

    close( GP_best_summary_output_unit )
    close( GP_minSSE_summary_output_unit )


endif ! myid == 0


!------------------------------------------------------------------

! deallocate variable dimension arrays

call deallocate_arrays1( )

deallocate( answer )
deallocate( output_array )


!------------------------------------------------------------------
if( myid == 0 )then

    write(GP_print_unit,'(//A)')  '0: NORMAL TERMINATION'
    write(GP_print_unit, '(3(A,1x,A,1x)//)') &
        '0: GPGACODE program version', trim(program_version), &
        '  branch:', trim( branch ) , &
        '  Last modified on:', trim( modification_date )

endif ! myid == 0

call MPI_FINALIZE(ierr)

!write(GP_print_unit,*)'0: aft mpi_finalize   ierr = ', ierr

stop

end program GPCODE_GA_lmdif_parameter_optimization_test
