!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
module GP_Parameters

implicit none

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the integer parameters for the Lotka Volterra Example
integer (kind=4), parameter :: n_Levels=4
integer (kind=4), parameter :: n_Node_Functions=7              ! number of possible node functions [needed for tree generations/alterations
integer (kind=4), parameter :: n_CODE_Equations=2         ! number of coupled equations
integer (kind=4), parameter :: n_Trees=((n_CODE_Equations+1)**2)-(n_CODE_Equations+1)
integer (kind=4), parameter :: n_Nodes=int(2**n_Levels)-1
integer (kind=4), parameter :: n_Time_Steps=2000     ! 400 days
real (kind=4), dimension(n_Levels), parameter :: Node_Probability = (/0.8,0.6,0.4,0./)  ! NOTE: Last value MUST BE 0.0!!!]
!off character*1, dimension(n_Node_Functions), parameter :: Node_Function_Type = (/'+','-','*','/','e','g','p'/)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the integer parameters for the Franks et al. Mar Bio NPZ Example
!off integer (kind=4), parameter :: n_Levels=5
!off integer (kind=4), parameter :: n_Functions=7              ! number of possible node functions [needed for tree generations/alterations
!off integer (kind=4), parameter :: n_CODE_Equations=3         ! number of coupled equations
!off integer (kind=4), parameter :: n_Trees=((n_CODE_Equations+1)**2)-(n_CODE_Equations+1)
!off integer (kind=4), parameter :: n_Nodes=int(2**n_Levels)-1
!off integer (kind=4), parameter :: n_Time_Steps=40*24*60     ! 400 days
!off real (kind=4), dimension(n_Levels), parameter :: Node_Probability = (/0.95,0.4,0.8,0.6,0.0/)
!off character*1, dimension(n_Node_Functions), parameter :: Node_Function_Type = (/'+','-','*','/','e','g','p'/)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

integer (kind=4), parameter :: n_Maximum_Number_Parameters = n_CODE_Equations+n_Nodes
integer (kind=4) :: n_Parameters

! this decides what the tree's shape needs to be like, i.e. bush or logpole pine
! Suggestion: Later this might be modulated with a phase, amplitude, etc. f-n]
! NOTE: the last value must be set to zero to set the last level as terminals only.

! GP routine-specific variables that can be modified
integer (kind=4), parameter :: n_GP_Generations = 1000
integer (kind=4), parameter :: n_GP_Individuals = 10

! GP Probability of a Tree being assigned
real (kind=4), parameter :: GP_Tree_Probability=0.5 ! Estimated from previous work by Joel Cohen

! Note: The next 4 parameters must add up to 1.0
real (kind=4), parameter :: GP_Elitist_Probability = 0.1                ! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation
real (kind=4), parameter :: GP_Asexual_Reproduction_Probability = 0.4   ! probability of asexual reproduction
real (kind=4), parameter :: GP_Crossover_Probability = 0.4              ! probability of sexual crossing of binary string 
real (kind=4), parameter :: GP_Mutation_Probability = 0.1               ! probability of mutation in binary string

real (kind=4), parameter :: GP_Set_Terminal_to_Parameter_Probability = 0.6  ! probability of setting a terminal node to a parameter

!off real (kind=8), parameter :: dt = 1.0D+0/(24.0D+0*60.0D+0)   ! [d^-1; 1 minute time step]
real (kind=8), parameter :: dt = 10.0D+0/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

end module GP_Parameters


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


module GA_Parameters

implicit none

! GA routine-specific variables that can be modified
integer (kind=4), parameter :: n_GA_Generations = 200
integer (kind=4), parameter :: n_GA_Individuals = 9

! Note: The next 4 parameters must add up to 1.0
real (kind=4), parameter :: GA_Elitist_Probability=0.1   ! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation
real (kind=4), parameter :: GA_Asexual_Reproduction_Probability=0.4   ! probability of asexual reproduction
real (kind=4), parameter :: GA_Crossover_Probability = 0.4  ! probability of sexual crossing of binary string 
real (kind=4), parameter :: GA_Mutation_Probability = 0.1   ! probability of mutation in binary string

end module GA_Parameters


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


module GP_Variables

use GP_Parameters

implicit none

integer (kind=4) :: n_GP_Elitists
integer (kind=4) :: n_GP_Asexual_Reproductions
integer (kind=4) :: n_GP_Crossovers
integer (kind=4) :: n_GP_Mutations

real (kind=4) :: GP_Population_Node_Parameters(n_GP_Individuals,n_Nodes,n_Trees)
real (kind=4) :: GP_Individual_Node_Parameters(n_Nodes,n_Trees)
real (kind=8) :: Runga_Kutta_Node_Parameters(n_Nodes,n_Trees)

real (kind=4) :: GP_Population_Parameter_Solution(n_GP_Individuals,n_Maximum_Number_Parameters)
real (kind=4) :: GP_Individual_Parameter_Solution(n_Maximum_Number_Parameters)

real (kind=4) :: GP_Population_Initial_Conditions(n_GP_Individuals,n_CODE_Equations)
real (kind=4) :: GP_Individual_Initial_Conditions(n_CODE_Equations)

real (kind=4) :: GP_Solution(0:n_Time_Steps,n_CODE_Equations)

integer (kind=4) :: GP_Adult_Population_Node_Type(n_GP_Individuals,n_Nodes,n_Trees)
integer (kind=4) :: GP_Child_Population_Node_Type(n_GP_Individuals,n_Nodes,n_Trees)
integer (kind=4) :: Parent_Tree_Swap_Node_Type(n_Nodes,2)
integer (kind=4) :: GP_Individual_Node_Type(n_Nodes,n_Trees)
integer (kind=4) :: Runga_Kutta_Node_Type(n_Nodes,n_Trees)

real (kind=8) :: Runga_Kutta_Initial_Conditions(n_CODE_Equations)
real (kind=8) :: Runga_Kutta_Solution(0:n_Time_Steps,n_CODE_Equations)

real (kind=4) :: GP_Individual_Lowest_SSE
real (kind=4) :: GP_Adult_Individual_SSE(n_GP_Individuals)       ! must be kept for re-evaluations of next generations
real (kind=4) :: GP_Child_Individual_SSE(n_GP_Individuals)       ! must be kept for re-evaluations of next generations
real (kind=4) :: GP_Individual_Ranked_Fitness(n_GP_Individuals)  ! must be kept for re-evaluations of next generations
real (kind=4) :: GP_Integrated_Ranked_Fitness(n_GP_Individuals)  ! must be kept for re-evaluations of next generations

logical, dimension(n_GP_Individuals) :: Run_GP_Calculate_Fitness

! random number routine variables
integer values(1:8), i_seed
integer, dimension(:), allocatable :: seed
real (kind=8) :: r
integer clock,n_seed

end module GP_Variables


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


module GA_Variables

use GA_Parameters

implicit none

integer (kind=4) :: n_GA_Elitists
integer (kind=4) :: n_GA_Asexual_Repoductions
integer (kind=4) :: n_GA_Crossovers
integer (kind=4) :: n_GA_Mutations

real (kind=4) :: GA_Individual_Lowest_SSE
real (kind=4) :: GA_Adult_Individual_SSE(n_GA_Individuals)       ! must be kept for re-evaluations of next generations
real (kind=4) :: GA_Child_Individual_SSE(n_GA_Individuals)       ! must be kept for re-evaluations of next generations
real (kind=4) :: GA_Integrated_SSE(n_GA_Individuals)             ! needed to support sexual and "tournament-style" reproduction
real (kind=4) :: GA_Individual_Ranked_Fitness(n_GA_Individuals)  ! must be kept for re-evaluations of next generations
real (kind=4) :: GA_Integrated_Ranked_Fitness(n_GA_Individuals)  ! must be kept for re-evaluations of next generations

logical, dimension(n_GA_Individuals) :: Run_GA_lmdif

end module GA_Variables


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


module GP_Data

use GP_Parameters

implicit none

real (kind=8) :: Data_Array(0:n_Time_Steps,n_CODE_Equations)

end module GP_Data


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
