module Runge_Kutta_Variables_module

use GP_Parameters_module

implicit none

!real(kind=8) :: Node_Values(n_nodes,n_trees)
real(kind=8),allocatable, dimension( : , : ) :: Node_Values

!real(kind=8) :: Tree_Evaluation(n_nodes,n_trees)
real(kind=8),allocatable, dimension( : , : ) :: Tree_Evaluation


!real(kind=8) :: Tree_Value(n_trees)
real(kind=8),allocatable, dimension( : ) :: Tree_Value


!integer(kind=4) :: Node_Eval_Type(n_nodes,n_trees)
integer(kind=4),allocatable, dimension( : , : ) :: Node_Eval_Type
!
!real(kind=8) :: bioflo(0:n_CODE_equations,0:n_CODE_equations)
real(kind=8),allocatable, dimension( : , : ) :: bioflo

character(str_len),allocatable, dimension( : , : ) :: bioflo_string

!real(kind=8) :: b_tmp(n_CODE_equations)
real(kind=8),allocatable, dimension( : ) :: b_tmp


!--------------------------------------------------------------------

! Runge-Kutta specific work arrays

!real(kind=8) :: kval(4,n_CODE_equations)
real(kind=8),allocatable, dimension( : , : ) :: kval

!real(kind=8) :: btmp(n_CODE_equations)
real(kind=8),allocatable, dimension( : ) :: btmp

!real(kind=8) :: fbio(n_CODE_equations)
real(kind=8),allocatable, dimension( : ) :: fbio


!! table to store 2**level - 1  for use in RK integration
!integer, parameter :: max_level = 10                                                                          
!integer, dimension(0:max_level) :: pow2_table                                                                 
 

!--------------------------------------------------------------------


!real(kind=8) :: left_node_value,right_node_value
!real(kind=8) :: cff

!integer (kind=4) :: icff
!integer (kind=4) :: iter

real(kind=8),dimension(4) :: Runge_Kutta_Time_Step

data Runge_Kutta_Time_Step /0.0D+0,0.5D+0,0.5D+0,1.0D+0/  ! fraction of R-K time step


end module Runge_Kutta_Variables_module
