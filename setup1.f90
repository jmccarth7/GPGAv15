subroutine setup1( ) 

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
integer(kind=4) :: ierror_t
integer(kind=4) :: ierror_m
integer(kind=4) :: ierror_tb
integer(kind=4) :: i_start_generation

!real(kind=8), allocatable, dimension(:) :: answer
!real(kind=8), allocatable, dimension(:) :: output_array

!real(kind=8) :: t1
!real(kind=8) :: t2

character(200) :: tree_descrip

!character(15),parameter :: program_version   = '201402.001_v12'
!character(10),parameter :: modification_date = '20140416'
!character(50),parameter :: branch  =  'master'


!---------------------------------------------------------------------------------------

!ierror_t  = 0
!ierror_m  = 0
!ierror_tb = 0
!
!GP_para_flag = .FALSE.
!ierror = 0
!
!GP_minSSE_Individual = 0
!GP_minSSE_generation = 0
!
!!--------------------------------------------------------------
!! current setup
!! lmdif runs only on best individual of each generation
!! no replacement for bad points in GP*n
!! no retry in setup_run_fcn or setup_run_lmdif
!! in random_real, boundary is 0.1
!! --  0.1 for random range [0,1], 0.9 for random range [0,50]
!! max interations in lmdif is 100
!!--------------------------------------------------------------
!
!
!call MPI_INIT(ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
!call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
!
!
!if( myid == 0 )then
!    write(6,'(/A,1x,I2/)') 'set1: numprocs = ', numprocs
!endif ! myid == 0
!
!
!!------------------------------------------------------------------
!
!if( myid == 0 )then
!    write(6,'(/A)') 'set1: call Global_Setup '
!endif ! myid == 0
!
!call Global_Setup()
!
!if( myid == 0 )then
!    write(6,'(A)') 'set1: aft call Global_Setup '
!endif ! myid == 0
!
!!------------------------------------------------------------------
!
!GP_para_flag = .FALSE.  ! .True.
!Lprint_lmdif = .TRUE.
!
!
!!------------------------------------------------------------------
!
!CALL RANDOM_SEED(size = n_seed)
!
!!------------------------------------------------------------------
!
!
!if( myid == 0 )then
!
!    write(6,'(/A/)') 'set1: run lmdif in parallel on each GP generation'
!    write(6,'(/A/)') 'set1: changed RK sub to make it faster'
!    write(6,'(/A/)') 'set1: using the old_elite_scheme in GP_Fit* GP_Tou*, GP_Mut*'
!    write(6,'(/A/)') 'set1: fast mod 1: remove GP diversity and tree printout     '
!    write(6,'(/A/)') 'set1: fast mod 1: remove GP_calc_fit, GP_ranking, summary printout  '
!    write(6,'(A)')'set1:  run with clean tree call for only processor 0'
!    write(6,'(A)')'set1:  run with no barrier before call GPCODE  '
!    write(6,'(A)')'set1:  and with no barrier after  call GPCODE  '
!    write(6,'(A)')'set1: removed barrier in GPCODE aft bcast of L_stop'
!    write(6,'(A)')'set1: removed several barriers in 0*f90 and GPCODE*f90'
!
!    !------------------------------------------------------
!    write(GP_print_unit, '(/3(A,1x,A,1x)//)') &
!    'set1: GPGACODE program version', trim(program_version), &
!    '  branch:', trim( branch ) , &
!    '  Last modified on:', trim( modification_date )
!    !------------------------------------------------------
!
!
!    ! read the control input from file  "GPCODE_cntl"
!
!    call read_cntl_stuff( ierror  )
!
!    !------------------------------------------------------
!
!    ! open output units
!
!    if( L_unit50_output )then
!        open( unit_gp_out, file = 'unit50.txt', &
!              form = 'unformatted', access='sequential', &
!              status = 'unknown' )
!    endif ! L_unit50_output
!
!
!    if( L_GPSSE_log )then
!
!        open( GPSSE_log_unit, file = 'GPSSE_log', &
!              form = 'formatted', access='sequential', &
!              status = 'unknown' )
!
!        open( GPSSE_best_log_unit, file = 'GPSSE_best_log', &
!              form = 'formatted', access='sequential', &
!              status = 'unknown' )
!
!    endif ! L_GP_log
!
!
!
!    if( L_GP_log )then
!        open( GP_log_unit, file = 'GP_log', &
!              form = 'unformatted', access='sequential', &
!              status = 'unknown' )
!    endif ! L_GP_log
!
!
!    if( L_GA_log )then
!        open( GA_log_unit, file = 'GA_log', &
!              form = 'unformatted', access='sequential', &
!              status = 'unknown' )
!    endif ! L_GA_log
!
!
!endif !   myid == 0
!
!!----------------------------------------------------------
!
!! ierror > 0 if read_cntl_stuff has encountered a problem
!! stop all processes in this case
!
!message_len =  1                                                                                         
!call MPI_BCAST( ierror, message_len,    &                                                             
!                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )                                                  
!
!if( ierror > 0 ) then
!
!    call MPI_FINALIZE(ierr)
!    stop 'ierror > 0'
!
!endif  ! ierror > 0
!
!!----------------------------------------------------------
!
!if( n_input_vars > 0 )then
!
!    ! read in the data number of points, number of vars
!
!    if( myid == 0 )then
!
!        write(6, '(/A)') 'set1: call read_input_data_size '
!        call read_input_data_size( )
!
!    endif !   myid == 0
!
!
!
!    !---------------------------------------------------------------------
!
!    ! broadcast number of data points, number of input variables
!
!    call MPI_BCAST( n_input_data_points, 1,    &
!                    MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
!
!    call MPI_BCAST( n_input_vars, 1,    &
!                    MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
!
!
!    call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?
!
!    !---------------------------------------------------------------------
!
!    !write(6, '(/A,2(1x,I6))') 'set1: myid, n_input_data_points', &
!    !                              myid, n_input_data_points
!    !write(6, '(/A,2(1x,I6))') 'set1: myid, n_input_vars', myid, n_input_vars
!
!
!    ! allocate input data names
!
!    !write(6, '(/A)') 'set1: allocate input_data_names'
!    allocate( input_data_names( set1:n_input_vars ) )
!
!
!    ! allocate input data array
!
!    !write(6, '(/A)') 'set1: allocate input_data_array'
!    allocate( input_data_array( set1:n_input_vars, n_input_data_points) )
!
!
!    call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?
!
!    !---------------------------------------------------------------------
!
!    ! read in the data into the arrays  input_data_names, input_data_array
!
!
!    if( myid == 0 )then
!
!        write(6, '(/A)') 'set1: call read_input_data'
!
!        call read_input_data( )
!
!        write(6,'(/A,2(1x,I6))') 'set1: n_input_data_points', &
!                                     n_input_data_points
!        write(6,'(/A,2(1x,I6))') 'set1: n_input_vars      =', n_input_vars
!        write(6,'(/A,2(1x,I6))') 'set1: n_functions_input =', n_functions_input
!
!    endif !   myid == 0
!

!endif  !  n_input_vars > 0

!---------------------------------------------------------------------

!sum_lmdif = 0.0d0

! broadcast the values read in by cpu 0 to others
!if( myid == 0 )then
!    write(6, '(/A/)') 'set1: call bcast1 '                      
!    flush(6)
!endif !   myid == 0


call bcast1()

!if( myid == 0 )then
!    write(6, '(/A/)') 'set1: after call bcast1 '                      
!    flush(6)
!endif !   myid == 0
!
!!------------------------------------------------------------------
!
!! call to random_seed moved to before read_cntl_stuff
!
!!! CALL RANDOM_SEED(size = n_seed)
!
!if( .not. allocated( seed ) )then
!    ALLOCATE(seed(n_seed))
!endif ! .not. allocated( seed )
!
!if( .not. allocated( current_seed ) )then
!    ALLOCATE(current_seed(n_seed))
!endif ! .not. allocated( current_seed )
!
!!if( myid == 0 )then
!!    write(6, '(/A,1x,I6,5x,L1/)') 'set1: myid, L_restart ', myid, L_restart             
!!    write(6,'(A,2(1x,I12))') 'set1: myid, n_seed ', myid, n_seed
!!endif !   myid == 0
!
!
!if( L_restart )then
!
!    seed(1:n_seed) = temp_seed(1:n_seed)
!
!    if( myid == 0 )then
!        write(6,'(A)') 'set1: temp_seed array '
!        do  i = 1, n_seed
!            write(6,'(I12,1x,I12)')  i, temp_seed(i)
!        enddo ! i
!        write(6,'(A)') ' '
!    endif !   myid == 0
!
!else
!
!    if( user_input_random_seed > 0 )then
!    
!        clock = user_input_random_seed
!    
!        if( myid == 0 )then
!            write(6,'(/A,1x,I12)') &
!                  'set1: user input random seed       clock = ', clock
!        endif !   myid == 0
!    
!        seed = user_input_random_seed + &
!                  37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)
!    else
!    
!        CALL SYSTEM_CLOCK(COUNT=clock)
!    
!        if( myid == 0 )then
!            write(6,'(/A,1x,I12)')&
!                  'set1: random seed input clock = ', clock
!        endif !   myid == 0
!    
!        seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)
!    
!    endif ! user_input_random_seed > 0
!
!endif ! L_restart
!
!
!CALL RANDOM_SEED(PUT = seed)
!
!
!if( myid == 0 )then
!
!    write(6,'(A,1x,I12)') 'set1: n_seed ', n_seed
!    write(6,'(A)') 'set1: seed array '
!    do  i = 1, n_seed
!
!        write(6,'(I12,1x,I12)')  i, seed(i)
!
!    enddo ! i
!    write(6,'(A)') ' '
!
!endif ! myid == 0


!write(6,'(A,1x,I3,1x,I12)') 'set1: myid,   n_seed ', myid, n_seed
!write(6,'(A,1x,I3)') 'set1: seed array ', myid
!do  i = 1, n_seed
!
!    write(6,'(I3,1x, I12,1x,I12)')  i, seed(i)
!
!enddo ! i
!write(6,'(A)') ' '


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

!write(6,'(A,1x,I3,2(1x,I12))') 'set1: myid,   n_seed, n_code_equations ', myid, n_seed, n_code_equations

!---------------------------------------------------------------------

! n_inputs is used in deser*2 to point to input values in rk_data_array

n_inputs = n_input_vars

!------------------------------------------------------------------


if( myid == 0 )then

    write(6,'(/A,1x,I6)')     'set1: n_code_equations ', n_code_equations
    write(6,'(A,1x,I6)')      'set1: n_variables      ', n_variables
    write(6, '(/A,2(1x,I6))') 'set1: n_input_vars     ', n_input_vars
    write(6, '(A,2(1x,I6))')  'set1: n_code_equations ', n_code_equations
    write(6, '(A,2(1x,I6)/)') 'set1: n_inputs         ', n_inputs

    call print_values1()

endif ! myid == 0

!------------------------------------------------------------------

! allocate variable dimension arrays



if( myid == 0 )then
    write(6, '(A,2(1x,I6)/)') 'set1: call allocate_arrays1 '            
    flush(6)
endif ! myid == 0 


call allocate_arrays1( )

!if( myid == 0 )then
!    write(6, '(A,2(1x,I6)/)') 'set1: after allocate_arrays1 '            
!    flush(6)
!endif ! myid == 0 


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
    write(6, '(A,2(1x,I6)/)') 'set1: before set_answer_arrays'            
    !flush(6)
endif ! myid == 0 


if( myid == 0 )then

    call set_answer_arrays( )

endif ! myid == 0


!if( myid == 0 )then
!    write(6, '(A,2(1x,I6)/)') 'set1: after set_answer_arrays'            
!    !flush(6)
!endif ! myid == 0 
!------------------------------------------------------------------------

! then broadcast the R-K result: Runge_Kutta_Solution


if( myid == 0 )then    ! 20131209
    write(GP_print_unit,'(/A/)') &
          'set1: time_step   Numerical_Code_Solution(time_step,1:n_CODE_equations)'
    do  i = 1, n_time_steps
        write(GP_print_unit,'(I6,2x,10(1x,E14.7))') &
              i, (Numerical_Code_Solution(i,jj), jj = 1,n_CODE_equations )
    enddo ! i
endif ! myid == 0




message_len = ( n_time_steps + 1 ) * n_CODE_equations
call MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )



Data_Array=Numerical_CODE_Solution        ! Matrix Operation

!--------------------------------------------------------------------------------

! zero so that later solutions don't have answer results in array

Numerical_CODE_Solution(1:n_time_steps, 1:n_code_equations) = 0.0d0 

message_len = ( n_time_steps + 1 ) * n_CODE_equations
call MPI_BCAST( Numerical_CODE_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

! compute the data variance with cpu 0 only, then broadcast results

! sets:
! Data_Variance
! Data_Variance_inv
if( myid == 0 )then
    write(6, '(A,2(1x,I6)/)') 'set1: before comp_data_variance'            
    !flush(6)
endif ! myid == 0 

if( myid == 0 )then    ! 20131209
    call comp_data_variance( )
endif ! myid == 0

!if( myid == 0 )then
!    write(6, '(A,2(1x,I6)/)') 'set1: after comp_data_variance'            
!    !flush(6)
!endif ! myid == 0 


message_len =  n_CODE_equations
call MPI_BCAST( Data_Variance_inv, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


!--------------------------------------------------------------------------------


! put the desired model parameters in the array:  answer

answer = 0.0d0 ! set all to zero

n_parameters = 0

if( myid == 0 )then    ! 20131209
    write(GP_print_unit,'(/A,1x,I6/)') 'set1: n_code_equations ', n_code_equations
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

    !if( myid == 0 )then
    !    write(6, '(A,2(1x,I6)/)') 'set1: after sse0_calc'
    !    !flush(6)
    !endif ! myid == 0 


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


!if( myid == 0 )then
!    write(6, '(A,2(1x,I6)/)') 'set1: after set_modified_indiv'            
!    !flush(6)
!endif ! myid == 0 



return 

endsubroutine setup1
