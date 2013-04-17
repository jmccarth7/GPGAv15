program GPCODE_GA_lmdif_parameter_optimization_test

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module


implicit none


!logical bad_node


!integer :: i
!integer :: nop
!integer :: message_len

!logical :: Lplot


real(kind=8), allocatable, dimension(:) :: answer
real(kind=8), allocatable, dimension(:) :: output_array

!real (kind=8) :: ssum, ssum2, totobs, dff

!real (kind=8) :: xcount

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




!------------------------------------------------------------------


    ! read the control input from file  "GPCODE_cntl"
    ! and broadcast the values read to the other cpu's

    call read_cntl_stuff( )

    open( unit_gp_out, form = 'unformatted', access='sequential', &
          status = 'unknown' )


!------------------------------------------------------------------


! set the scalar values for the model


call init_values( 0 )



!------------------------------------------------------------------


    write(6,'(A,1x,I10)')    '0: n_GA_individuals           ', &
                                 n_GA_individuals
    write(6,'(A,1x,I10)')    '0: n_time_steps               ', &
                                 n_time_steps
    write(6,'(A,1x,I10)')    '0: n_GA_Generations           ', &
                                 n_GA_Generations

    write(6,'(/A,1x,I10)')   '0: n_gp_individuals           ', &
                                 n_gp_individuals
    write(6,'(A,1x,I10)')    '0: n_gp_generations           ', &
                                 n_gp_generations

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

!off include 'Lotka_Volterra_Example_Set_Up.f901'  ! replaced with routine init_values
!off include 'Franks_etal_NPZ_Mar_Bio_Example_Set_Up.f901'

! set the desired 'twin experiment' population node type
! and parameter using the info from the set up file

GP_Node_Type_Answer       = GP_Individual_Node_Type       ! Matrix Operation
GP_Node_Parameters_Answer = GP_Individual_Node_Parameters ! Matrix Operation


do i_GP_Individual=1,n_GP_individuals
  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
                        GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)
enddo

    write(unit_gp_out) GP_Node_Type_for_Plotting

! set the initial population node type using the info obtained from the set up file
! set the Initial Conditions, Model Parameters and Node Type for the 'twin experiment case'

GP_Population_Node_Type(1,1:n_nodes,1:n_trees)=GP_Individual_Node_Type(1:n_nodes,1:n_trees)


! initialize the biological data fields

Runge_Kutta_Solution(0,1:n_CODE_equations)=Runge_Kutta_Initial_Conditions ! Array Assignment

Runge_Kutta_Node_Parameters = GP_Individual_Node_Parameters  ! Matrix Operation
Runge_Kutta_Node_Type=GP_Individual_Node_Type                ! Matrix Operation


!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

call Runge_Kutta_Box_Model


!------------------------------------------------------------------------

! then broadcast the R-K result to all processors




!Data_Array=Runge_Kutta_Solution          ! Matrix Operation


!--------------------------------------------------------------------------------

! compute the data_variance  -- to be used in computing SSE

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!    write(6,'(A)') ' '
!
!do i_CODE_equation=1,n_CODE_equations
!
!  ssum  = 0.0D+0
!  ssum2 = 0.0D+0
!
!  do  i_time_step=0,n_time_steps
!      ssum  = ssum  +  Data_Array(i_time_step,i_CODE_equation)
!      ssum2 = ssum2 +  Data_Array(i_time_step,i_CODE_equation)**2
!  enddo !   i_time_step
!
!  totobs=dble(n_time_steps+1)
!
!  dff=( (totobs*ssum2)-(ssum**2) ) / ( totobs*(totobs-1.0D+0) )
!
!  if( dff .gt. 0.0D+0) then  ! set variance to observed variance for normalize by the s.d.
!      Data_Variance(i_CODE_equation)=dff
!  else ! set variance to 1.0 for normalization to be 'unaltered'
!      Data_Variance(i_CODE_equation)=1.0D+0
!  endif !   dff .gt. 0.0D+0
!
!  if(  abs( Data_Variance(i_CODE_equation) ) < 1.0D-30 )then
!       write(6,'(/A,1x,I6,2x,E15.7)') &
!       '0: i_CODE_equation, Data_Variance(i_CODE_equation) ', &
!           i_CODE_equation, Data_Variance(i_CODE_equation)
!       write(6,'(A/)') '0: bad data variance -- stopping program '
!       stop 'bad data var'
!
!  endif ! abs( Data_Variance(i_CODE_equation) ) < 1.0D-30
!
!      write(6,'(A,1x,I6,2x,E15.7)') &
!           '0: i_CODE_equation, Data_Variance(i_CODE_equation) ', &
!               i_CODE_equation, Data_Variance(i_CODE_equation)
!
!enddo !  i_CODE_equation
!
!    write(6,'(A)') ' '
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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


    write(*,'(/A,1x,I6)') 'Total Parameters for this run = ',n_parameters
    write(6,'(A,1x,I6)')  '0: n_GA_individuals ', n_GA_individuals

    !write(*,'(A,1x,I6)') 'number of initial conditions   = ',n_initial_conditions
    write(*,'(A,1x,I6)')   '0: n_time_steps     ', n_time_steps
    write(*,'(A,1x,I6)')   '0: n_GA_Generations ', n_GA_Generations
    write(*,'(A,1x,I6)')   '0: n_levels         ', n_levels
    write(*,'(A,1x,I6)')   '0: n_functions      ', n_functions
    write(*,'(A,1x,I6)')   '0: n_CODE_equations ', n_CODE_equations
    write(*,'(A,1x,I6)')   '0: n_parameters     ', n_parameters
    write(*,'(A,1x,I6)')   '0: n_trees          ', n_trees
    write(*,'(A,1x,I6)')   '0: n_nodes          ', n_nodes
    write(*,'(A,1x,E15.7)')'0: dt (days)        ', dt
    write(*,'(A,1x,E15.7)')'0: dt (min)         ', dt * 1440.d0
    !write(*,*) individual_fitness
    !write(*,*) data_array

    !write(6,'(/A,1x,I10)')  '0: n_gp_individuals            ', &
    !                            n_gp_individuals
    !write(6,'(A,1x,I10)')   '0: n_gp_generations            ', &
    !                            n_gp_generations

    !write(*,'(A)') '0: code calls lmdif for all individuals and generations '
    write(*,'(A)') '0: code calls lmdif only for best individual on last generation'

    !write(*,'(A)') '0: regular random_real initialization of child parameters '
    write(*,'(A)') '0: 2-range random_real initialization of child parameters '

    write(*,'(A//)') '0: replacement of bad individuals used'
    !write(*,'(A//)') '0: NO replacement of bad individuals '


    ! calculate the generation interval for printing the list of children

    child_print_interval = n_time_steps / number_child_prints


    !-----------------------------------------------------------------------------

    ! this call calculates the sse value for RK solution = 0 for all time steps
    ! note:  sse0 is only used by cpu 0 which does all fitness calculations

    !call sse0_calc( )

    !-----------------------------------------------------------------------------

    !write(6,'(/A)') ' '
    !do  i = 1, n_parameters
    !    write(6,'(A,1x,I6,2x,E24.16)') '0: i, answer(i) ', i, answer(i)
    !enddo ! i
    !write(6,'(/A)') ' '
    !-----------------------------------------------------------------------------


    open( GA_output_unit, file='output_parameters', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )

    open( GP_output_unit, file='GP_output_parameters', &
          form = 'formatted', access = 'sequential', &
          status = 'unknown' )



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the Genetic Programming optimization routine for the Binary Tree Evolution
! with the embedded GA_lmdif parameter optimization scheme
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    write(6,'(/A,1x,I6)')   '0: n_gp_individuals           ', &
                                n_gp_individuals
    write(6,'(A,1x,I6/)')   '0: n_gp_generations           ', &
                                n_gp_generations





!------------------------------------------------------------------------------



do  i_GP_Generation=1,n_GP_Generations


        write(*,'(/A,1x,I6,1x,A,1x,I6/)') &
              '0: GP Generation # ',i_GP_Generation,&
                ' is underway.   n_Nodes * n_Trees = ', n_Nodes*n_Trees


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


        write(6,'(A,1x,I6/)') &
          '0: call GP_Tree_Build   Generation = ',i_GP_Generation

        ! initialize the GP_Adult_Population_Node_Type array with random trees

        call GP_Tree_Build


    endif ! i_GP_Generation .eq. 1


    !-----------------------------------------------------------------------------------------

    call Initialize_GA_Child_Parameters( GP_Individual_Node_Parameters )

    !-----------------------------------------------------------------------------------------

    ! print trees before tree clean


    write(6,'(/A/)')  '0: ############################################################################'
    write(6,'(A/)')  '0: trees before call to GP_Clean_Tree_Nodes'


    do  i_GP_individual = 1, n_GP_individuals

        write(6,'(/A,1x,I6)')  '0: i_GP_indiv ', i_GP_individual
        write(6,'(A)') '0: i_tree nodes '
        write(6,'(A)') '            1  2  3  4  5  6  7  8  9 10 11 12 13 14 15'

        do  i_Tree=1,n_Trees

            write(6,'(I6,4x,20(1x,I2))' ) &
                 i_tree, GP_Adult_Population_Node_Type(i_GP_individual, 1:n_nodes, i_tree)

        enddo ! i_tree
    enddo  ! i_GP_individual


    write(6,'(/A/)')  '0: ############################################################################'



    i_GP_individual = 1
    call show_tree( i_GP_individual, GP_Individual_Node_Parameters )



    !-----------------------------------------------------------------------------------------

    ! sweep through all the GP_Adult_Population_Node_Type
    ! to replace function nodes that have both terminals set as parameters
    ! and set the node to a parameter itself

    write(6,'(/A,1x,I6/)') &
          '0: call GP_Clean_Tree_Nodes   Generation = ', &
                                    i_GP_Generation


    call GP_Clean_Tree_Nodes


    !-----------------------------------------------------------------------------------------

    ! print trees after call to GP_Clean_Tree_Nodes


    write(6,'(/A)')  '0: ############################################################################'
    write(6,'(A/)')  '0: trees after call to GP_Clean_Tree_Nodes'



    do  i_GP_individual = 1, n_GP_individuals

        write(6,'(/A,1x,I6)')  '0: i_GP_indiv ', i_GP_individual
        write(6,'(A)') '0: i_tree nodes '
        write(6,'(A)') '            1  2  3  4  5  6  7  8  9 10 11 12 13 14 15'

        do  i_Tree=1,n_Trees

            write(6,'(I6,4x,20(1x,I2))' ) &
                 i_tree, GP_Adult_Population_Node_Type(i_GP_individual, 1:n_nodes, i_tree)

        enddo ! i_tree
    enddo  ! i_GP_individual


    write(6,'(/A/)')  '0: ############################################################################'



    !-----------------------------------------------------------------------------------------


enddo !  i_GP_Generation



! plot results

!    Lplot = .true.
!    if( Lplot ) call plot_results(Runge_Kutta_Solution,n_time_steps, n_CODE_equations )


!------------------------------------------------------------------

! deallocate variable dimension arrays

                                                                                                                          
! deallocate variable dimension arrays                                                                                    
                                                                                                                          
call deallocate_arrays1( )                                                                                                
                                                                                                                          
deallocate( answer )                                                                                                      
deallocate( output_array )                                                                                                
                                                                                                                          
                                                                                                                          
!------------------------------------------------------------------
!close( GA_output_unit )
!close( GP_output_unit )
!close( unit_gp_out )
!------------------------------------------------------------------






end program GPCODE_GA_lmdif_parameter_optimization_test
