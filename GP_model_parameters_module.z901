module GP_model_parameters_module

! set the GP related parameters


implicit none



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

integer(kind=4), parameter :: n_CODE_Forcing=4

integer(kind=4), parameter :: n_Tracked_Resources = 1 ! number of different resources being tracked


!real(kind=8), parameter :: Delta_Time_in_Days = 1.0D+0/(24.0D+0*2.0D+0) ! 1/2 hour
!integer(kind=4), parameter :: n_Time_Steps=5*365*24*2          ! 5 years in half-hour increments

real(kind=8) :: Delta_Time_in_Days 

character (len=*), parameter :: Model_Name = 'Fasham'



end module GP_model_parameters_module
