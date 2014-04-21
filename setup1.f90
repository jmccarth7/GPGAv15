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


!real(kind=8) :: t1
!real(kind=8) :: t2

character(200) :: tree_descrip


!---------------------------------------------------------------------------------------

!------------------------------------------------
!write(6,'(/A/)') &
!       'set1: myid,  call load_pow2_level      '

! load table of 2**ilevel - 1  for RK process

call load_pow2_level(  )

!------------------------------------------------

! set the scalar values for the model

! sets:
! n_levels
! n_functions
! n_CODE_equations
! n_trees
! n_nodes


call init_values( 0 )


n_Variables = n_CODE_equations

write(6,'(A,1x,I3,2(1x,I12))') &
       'set1: myid,   n_seed, n_code_equations ', &
              myid,   n_seed, n_code_equations

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

if( myid == 0 )then
    write(6, '(A,2(1x,I6)/)') 'set1: after allocate_arrays1 '            
    flush(6)
endif ! myid == 0 



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


if( myid == 0 )then
    write(6, '(A,2(1x,I6)/)') 'set1: after set_answer_arrays'            
    flush(6)
endif ! myid == 0 
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

! broadcast SSE0

message_len = 1 
call MPI_BCAST( SSE0, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


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
