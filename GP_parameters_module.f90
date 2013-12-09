module GP_parameters_module

! set the GP related parameters

!!!!!use GP_model_parameters_module

implicit none


!--------------------------------------------------------------------------------------

!!!!!! module GP_model_parameters_module



integer :: n_levels
integer :: n_functions
integer(kind=4) :: n_CODE_equations


! number of possible node functions
! [needed for tree generations/alterations]

integer(kind=4) :: n_Node_Functions   ! =7

! n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
integer(kind=4) :: n_trees

! n_nodes = int(2**n_levels)-1
integer(kind=4) :: n_nodes 

! n_maximum_number_parameters = n_CODE_equations +  n_nodes
integer(kind=4) :: n_maximum_number_parameters


integer(kind=4) :: n_Variables ! = 7

integer(kind=4), parameter :: n_CODE_Forcing = 4

integer(kind=4), parameter :: n_Tracked_Resources = 1 ! number of different resources being tracked


!real(kind=8), parameter :: Delta_Time_in_Days = 1.0D+0/(24.0D+0*2.0D+0) ! 1/2 hour
!integer(kind=4), parameter :: n_Time_Steps=5*365*24*2          ! 5 years in half-hour increments

real(kind=8) :: Delta_Time_in_Days 

character (len=*), parameter :: Model_Name = 'Fasham'


!--------------------------------------------------------------------------------------



! if this value is non-zero, then the random number generator uses
! this instead of the system clock value
! used for debugging the program so that multiple runs will have the same
! set of random numbers

integer(kind=4) :: user_input_random_seed

! this value is the minimum number of parameters for a valid model
integer(kind=4),parameter :: min_N_param = 2

integer(kind=4),parameter :: GP_print_unit = 6


integer(kind=4),parameter :: GP_summary_output_unit = 40

integer(kind=4),parameter :: unit_gp_out = 50
logical ::   L_unit50_output

integer(kind=4),parameter :: GP_log_unit = 80
logical ::   L_GP_log


integer(kind=4),parameter :: GP_output_unit = 30
logical ::          L_GP_output_parameters

logical ::   L_print_equations

integer(kind=4) :: n_GP_individuals

integer(kind=4) :: n_GP_Generations


integer(kind=4) :: n_GP_parameters

character(30) :: model



integer(kind=4) :: n_parameters

! this decides what the tree's shape needs to be like, i.e. bush or logpole pine
! Suggestion: Later this might be modulated with a phase, amplitude, etc. f-n]
! NOTE: the last value must be set to zero to set the last level as terminals only.

!-------------------------------------------------------------------
integer, parameter :: str_len = 500

character(160) :: tree_node_string
character(3) :: node_element_string
!-------------------------------------------------------------------

! n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
! n_nodes = int(2**n_levels)-1
! n_maximum_number_parameters = n_CODE_equations +  n_nodes



! GP Probability of a Tree being assigned
! Estimated from previous work by Joel Cohen

!real (kind=4), parameter :: GP_Tree_Probability=0.5
real(kind=8) :: GP_Tree_Probability !=0.5 ! Estimated from previous work by Joel Cohen


! set the parameters for the Lotka Volterra Example
!real(kind=8), dimension(n_Levels), parameter :: &
!    Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

! Node_Probability = !(/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]
real(kind=8), allocatable, dimension(:) :: Node_Probability

!---------------------------------------------------------------------------------------

! Note: The next 4 parameters must add up to 1.0

! Keeps the top n_GP_Elitists of the
! Best Fit Individuals from Generation to Generation

!real (kind=4), parameter :: GP_Elitist_Probability = 0.1
real(kind=8) :: GP_Elitist_Probability

!real(kind=4),parameter :: &
!  GP_Asexual_Reproduction_Probability =0.4 ! prob of asexual reproduction
real(kind=8) :: GP_Asexual_Reproduction_Probability

!real(kind=4),parameter :: &
!  GP_Crossover_Probability=0.4 ! prob of sexual crossing of binary string
real(kind=8) :: GP_Crossover_Probability

!real (kind=4), parameter :: &
!  GP_Mutation_Probability = 0.1 ! prob of mutation in binary string
real(kind=8) :: GP_Mutation_Probability


!---------------------------------------------------------------------------------------


!character(4), dimension( n_nodes, n_trees )       ::  node_type_string
!character(4), dimension( n_nodes, n_trees )       ::  node_parameters_string
!character(str_len), dimension( n_nodes, n_trees ) ::  tree_evaluation_string

character(4), allocatable, dimension( : , : )       ::  node_type_string
character(4), allocatable, dimension( : , : )       ::  node_parameters_string
character(str_len), allocatable, dimension( : , : ) ::  tree_evaluation_string

character(str_len), allocatable, dimension( : )     ::  tree_value_string

!-------------------------------------------------------------------

! with 2500 steps, the LV predator-prey cycle curve is closed
!integer(kind=4), parameter :: n_time_steps= 2500 ! 8 ! 10

integer(kind=4) :: n_time_steps


!-------------------------------------------------------------------


! this is the interval to use to determine if the child values
! should be printed

! print if   mod( i_ga_generation, child_print_interval ) == 0

integer(kind=4)  :: GA_child_print_interval
integer(kind=4)  :: GP_child_print_interval



! this is the total number of child printouts for the entire run

!integer(kind=4),parameter :: number_GA_child_prints = 10 ! 20
!integer(kind=4),parameter :: number_GP_child_prints = 10 ! 20

integer(kind=4) :: number_GA_child_prints ! = 10 ! 20
integer(kind=4) :: number_GP_child_prints ! = 10 ! 20

!-------------------------------------------------------------------

logical :: L_bad_result


!----------------------------------------------------------------------------------------



integer (kind=4) :: ier_file,idummy,iwkid,iwktype  ! NCAR Graphics

! orig !character (len=*), parameter :: output_dir = 'Output/'//Model_Name
character (len=*), parameter :: output_dir = '.' 

character (len=*), parameter :: input_dir = 'Input'


!real (kind=4), dimension(0:n_Time_Steps) :: x_Time_Steps = 0.D+0 !The temporal unit depends on the delta time in days
real (kind=4), dimension(:),allocatable :: x_Time_Steps 



! real (kind=8), parameter :: dt = 1.0D+0 /(24.0D+0*60.0D+0)   ! [d^-1; 1 minute time step]
!!real (kind=8), parameter :: dt = 10.0D+0/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

real(kind=8) :: dt

!------------------------------------------------------------------------------

! used in random_real


real(kind=8) ::  random_scale_large
real(kind=8) ::  random_scale_small
real(kind=8) ::  random_scale_fraction



!--------------------------------------------------------------------

! table to store 2**level - 1  for use in RK integration

integer, parameter :: max_level = 10
integer, dimension(0:max_level) :: pow2_table

!--------------------------------------------------------------------

! number of parameters for each GP individual

integer, allocatable, dimension(:) :: GP_n_parms
logical :: GP_para_flag


end module GP_parameters_module
