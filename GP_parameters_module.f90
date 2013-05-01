module GP_parameters_module

! set the GP related parameters

use GP_model_parameters_module

implicit none


integer(kind=4),parameter :: GA_output_unit = 20
integer(kind=4),parameter :: GP_output_unit = 30

integer(kind=4),parameter :: GP_summary_output_unit = 40

integer(kind=4),parameter :: GP_print_unit = 6

integer(kind=4),parameter :: unit_gp_out = 50

!integer(kind=4), parameter :: n_GP_individuals=1
integer(kind=4) :: n_GP_individuals

integer(kind=4) :: n_GP_Generations

!integer(kind=4) :: i_GP_individual
!integer(kind=4) :: i_GP_Best_Parent
!integer(kind=4) :: i_GP_Generation
!integer(kind=4) :: j_GP_Individual
!integer(kind=4) :: i_Tree
!integer(kind=4) :: i_Node

integer(kind=4) :: n_GP_parameters

character(30) :: model

integer(kind=4) :: n_linked_parms
integer(kind=4) :: n_linked_parms_dim

integer(kind=4),allocatable,dimension(:,:)  :: linked_parms

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

integer(kind=4) :: i_level
integer(kind=4) :: i_function
integer(kind=4) :: i_CODE_equation
integer(kind=4) :: j_CODE_equation
integer(kind=4) :: i_node_left
integer(kind=4) :: i_node_right


integer(kind=4) :: n_parameters

! this decides what the tree's shape needs to be like, i.e. bush or logpole pine
! Suggestion: Later this might be modulated with a phase, amplitude, etc. f-n]
! NOTE: the last value must be set to zero to set the last level as terminals only.

!-------------------------------------------------------------------
integer, parameter :: str_len = 80
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

!old >>>>>>>>
!  NOTE: in the next 2 parameters:
!        GP_Crossover_Probability+GP_Mutation_Probability must be <= 1.

!real(kind=8) :: GP_Crossover_Probability =0.5d0 ! prob of sexual crossing of binary string
!real(kind=8) :: GP_Mutation_Probability  =0.1d0 ! prob of mutation in binary string
!old <<<<<<<<


!---------------------------------------------------------------------------------------


!character(4), dimension( n_nodes, n_trees )       ::  node_type_string
!character(4), dimension( n_nodes, n_trees )       ::  node_parameters_string
!character(str_len), dimension( n_nodes, n_trees ) ::  tree_evaluation_string

character(4), allocatable, dimension( : , : )       ::  node_type_string
character(4), allocatable, dimension( : , : )       ::  node_parameters_string
character(str_len), allocatable, dimension( : , : ) ::  tree_evaluation_string

!-------------------------------------------------------------------

! with 2500 steps, the LV predator-prey cycle curve is closed
!integer(kind=4), parameter :: n_time_steps= 2500 ! 8 ! 10

integer(kind=4) :: n_time_steps


!-------------------------------------------------------------------


! this is the interval to use to determine if the child values
! should be printed

! print if   mod( i_ga_generation, child_print_interval ) == 0

integer(kind=4)           :: GA_child_print_interval
integer(kind=4)           :: GP_child_print_interval



! this is the total number of child printouts for the entire run

integer(kind=4),parameter :: number_GA_child_prints = 2 ! 20
integer(kind=4),parameter :: number_GP_child_prints = 2 ! 20

!-------------------------------------------------------------------


integer (kind=4) :: i_time_step

logical :: L_bad_result

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



integer :: i_parameter


!----------------------------------------------------------------------------------------



! real (kind=8), parameter :: dt = 1.0D+0 /(24.0D+0*60.0D+0)   ! [d^-1; 1 minute time step]
!!real (kind=8), parameter :: dt = 10.0D+0/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

real(kind=8) :: dt

!------------------------------------------------------------------------------

! used in random_real


real(kind=8) ::  random_scale_large
real(kind=8) ::  random_scale_small
real(kind=8) ::  random_scale_fraction







end module GP_parameters_module
