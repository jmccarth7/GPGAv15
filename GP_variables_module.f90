module GP_variables_module

use GP_Parameters_module

implicit none

real(kind=8) :: Individual_Fitness
real(kind=8) :: Individual_SSE_best_parent

integer(kind=4) :: n_GP_Elitists
integer(kind=4) :: n_GP_Asexual_Reproductions
integer(kind=4) :: n_GP_Crossovers
integer(kind=4) :: n_GP_Mutations


! GP_Node_Parameters_Answer(n_Nodes,n_Trees)
real(kind=8), allocatable, dimension(:,:)  :: GP_Node_Parameters_Answer

! GP_Node_Type_Answer(n_Nodes,n_Trees)
integer(kind=4), allocatable, dimension(:,:)  :: GP_Node_Type_Answer

! GP_Node_Type_for_Plotting(9,n_Nodes,n_Trees)
integer(kind=4), allocatable, dimension(:,:,:)  :: GP_Node_Type_for_Plotting

!real(kind=8) :: GP_Population_Node_Parameters(n_nodes,n_trees,n_GP_individuals)
real(kind=8),allocatable,dimension(:,:,:)     :: GP_Population_Node_Parameters

!real(kind=8) :: GP_Individual_Node_Parameters(n_nodes,n_trees)
real(kind=8),allocatable,dimension(:,:)       :: GP_Individual_Node_Parameters

!real(kind=8) :: GP_Population_Parameter_Solution(n_maximum_number_parameters,n_GP_individuals)
!real(kind=8),allocatable,dimension(:,:)       :: GP_Population_Parameter_Solution


! in GPCODE*n.f90
real(kind=8),allocatable, dimension(:)  :: buffer
real(kind=8),allocatable, dimension(:)  :: buffer_recv

! in GP_para*lmdif*f90
real(kind=8),allocatable, dimension(:)  :: buffer2
real(kind=8),allocatable, dimension(:)  :: buffer2_recv


! GP_Solution(0:n_Time_Steps,n_CODE_Equations)
!real(kind=8),allocatable,dimension(:,:) :: GP_Solution

! GP_diversity_index(n_GP_Individuals)
integer(kind=4),allocatable,dimension(:) :: GP_diversity_index

! GP_Adult_Population_Node_Type(n_Nodes,n_Trees,n_GP_Individuals)
integer(kind=4),allocatable,dimension(:,:,:) :: GP_Adult_Population_Node_Type

! GP_Child_Population_Node_Type(n_Nodes,n_Trees,n_GP_Individuals)
integer(kind=4),allocatable,dimension(:,:,:) :: GP_Child_Population_Node_Type

! Parent_Tree_Swap_Node_Type(n_Nodes,2)
integer(kind=4),allocatable,dimension(:,:) :: Parent_Tree_Swap_Node_Type

! GP_Individual_Node_Type(n_Nodes,n_Trees)
integer(kind=4),allocatable,dimension(:,:) :: GP_Individual_Node_Type

!real(kind=8) :: GP_Individual_Parameter_Solution(n_maximum_number_parameters)
!!real(kind=8),allocatable,dimension(:)         :: GP_Individual_Parameter_Solution

! GP_Adult_Population_Parameter_Solution(n_Maximum_Number_Parameters,n_GP_Individuals)
real (kind=8),allocatable,dimension(:,:)  :: GP_Adult_Population_Parameter_Solution

! GP_Child_Population_Parameter_Solution(n_Maximum_Number_Parameters,n_GP_Individuals)
real (kind=8),allocatable,dimension(:,:)  :: GP_Child_Population_Parameter_Solution

!real(kind=8) :: GP_Population_Initial_Conditions(n_CODE_equations,n_GP_individuals)
real(kind=8),allocatable,dimension(:,:)       :: GP_Population_Initial_Conditions

!real(kind=8) :: GP_Individual_Initial_Conditions(n_CODE_equations)
real(kind=8),allocatable,dimension(:)         :: GP_Individual_Initial_Conditions

!real(kind=8) :: GP_Population_Fitness(n_GP_individuals)
real(kind=8),allocatable,dimension(:)         :: GP_Population_Fitness

!integer(kind=4) :: GP_Population_Node_Type(n_Nodes,n_Trees,n_GP_Individuals)
integer(kind=4),allocatable,dimension(:,:,:) :: GP_Population_Node_Type


! GP_Individual_N_GP_param
integer(kind=4),allocatable,dimension(:) :: GP_Individual_N_GP_param
!------------------------------------------------------------------------------

! Runge_Kutta_Node_Type(n_Nodes,n_Trees)
integer(kind=4),allocatable,dimension(:,:) :: Runge_Kutta_Node_Type
integer(kind=4),allocatable,dimension(:,:) :: RK_Node_Type

!real(kind=8) :: Runge_Kutta_Node_Parameters(n_nodes,n_trees)
real(kind=8),allocatable,dimension(:,:)       :: Runge_Kutta_Node_Parameters
real(kind=8),allocatable,dimension(:,:)       :: RK_Node_Parameters

!real (kind=8) :: Runge_Kutta_Initial_Conditions(n_CODE_equations)
real (kind=8),allocatable,dimension(:)        :: Runge_Kutta_Initial_Conditions
real (kind=8),allocatable,dimension(:)        :: RK_Initial_Conditions

!real (kind=8) :: Runge_Kutta_Solution(0:n_time_steps,n_CODE_equations)
!
real(kind=8),allocatable, dimension(:,:) :: Runge_Kutta_Solution
real(kind=8),allocatable, dimension(:,:) :: RK_Solution

!------------------------------------------------------------------------------

real(kind=8) :: GP_Individual_Lowest_SSE

!---------------------------------------------------------------------------

! must be kept for re-evaluations of next generations >>>

! GP_Adult_Individual_SSE(n_GP_Individuals)
! GP_Child_Individual_SSE(n_GP_Individuals)
! GP_Individual_Ranked_Fitness(n_GP_Individuals)
! GP_Integrated_Ranked_Fitness(n_GP_Individuals)

real(kind=8),allocatable, dimension(:) :: GP_Adult_Individual_SSE
real(kind=8),allocatable, dimension(:) :: GP_Child_Individual_SSE

real(kind=8),allocatable, dimension(:) :: GP_Adult_Population_SSE  ! ???
real(kind=8),allocatable, dimension(:) :: GP_Child_Population_SSE  ! ???

real(kind=8),allocatable, dimension(:) :: GP_Individual_Ranked_Fitness
real(kind=8),allocatable, dimension(:) :: GP_Integrated_Ranked_Fitness

real(kind=8),allocatable, dimension(:) :: GP_Population_Ranked_Fitness               ! ???
real(kind=8),allocatable, dimension(:) :: GP_Integrated_Population_Ranked_Fitness    ! ???


! must be kept for re-evaluations of next generations <<<

!---------------------------------------------------------------------------

!logical, dimension(n_GP_Individuals) :: Run_GP_Calculate_Fitness
logical, allocatable, dimension(:) :: Run_GP_Calculate_Fitness

! random number routine variables

integer values(1:8), i_seed
integer, dimension(:), allocatable :: seed
real (kind=8) :: r
integer clock,n_seed

!----------------------------------------------------

end module GP_variables_module
