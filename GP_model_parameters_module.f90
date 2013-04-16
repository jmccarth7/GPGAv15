module GP_model_parameters_module

! set the GP related parameters


implicit none



integer :: n_levels
integer :: n_functions
integer(kind=4) :: n_CODE_equations


! number of possible node functions
! [needed for tree generations/alterations

integer(kind=4) :: n_Node_Functions   ! =7

! n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)
integer(kind=4) :: n_trees

! n_nodes = int(2**n_levels)-1
integer(kind=4) :: n_nodes

! n_maximum_number_parameters = n_CODE_equations +  n_nodes
integer(kind=4) :: n_maximum_number_parameters



end module GP_model_parameters_module
