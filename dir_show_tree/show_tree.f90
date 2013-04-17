subroutine show_tree( i_GA_indiv,  child_parameters )

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!use mpi
!use mpi_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none


integer, intent(in)  ::  i_GA_indiv

integer ::  iflag


! lmdif arrays and variables

!real (kind=8) :: x_LMDIF(n_maximum_number_parameters)
!real (kind=8) :: fvec(n_time_steps)
real(kind=8), allocatable, dimension(:) :: x_LMDIF
real(kind=8), allocatable, dimension(:) :: fvec


real(kind=8), parameter :: tol = 1.0d-30
real(kind=8), parameter :: epsfcn = 1.0d-6    ! original
real(kind=8), parameter :: factor=1.0D+0
real(kind=8), parameter :: zero = 0.0d0

integer (kind=4) ::   info


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

!integer(kind=4) :: individual_quality(n_GA_individuals)


!real(kind=8) :: child_parameters(n_GA_individuals,n_maximum_number_parameters)
real(kind=8) :: child_parameters(n_GA_individuals,1:n_trees, 1:n_nodes )

external :: fcn

!--------------------------------------------------------------------------------------------

write(6,'(/A,1(1x,I6)/)') 'shotr: i_GA_indiv ', i_GA_indiv
write(6,'(/A,1x,I6/)') 'shotr: n_maximum_number_parameters ', n_maximum_number_parameters

allocate( x_LMDIF(n_maximum_number_parameters) ) 
allocate( fvec(n_time_steps) ) 

x_LMDIF(1:n_maximum_number_parameters) = 0.0D0

!do  i_parameter=1,n_parameters
i_parameter = 0
                                                                                                                          
do  i_CODE_equation=1,n_CODE_equations                                                                                    
    i_parameter =  i_parameter + 1                                                                                           
    x_LMDIF(i_parameter)=Runge_Kutta_Initial_Conditions(i_CODE_equation)                                                  
enddo ! i_CODE_equation                                                                                                   
  
do  i_tree =1, n_trees
    do  i_node =1, n_nodes

        if( GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) == -9999 ) cycle 

        i_parameter = i_parameter + 1

        X_LMDIF(i_parameter) = child_parameters(i_GA_indiv,i_tree, i_node )
    
        !write(6,'(A,4(1x,I6),1x,E20.10)') &
        !      'shotr:1 i_GA_indiv,i_tree, i_node,i_parameter, child_parameters ', &
        !               i_GA_indiv,i_tree, i_node,i_parameter, &
        !               child_parameters(i_GA_indiv,i_tree, i_node)
        write(6,'(A,2(1x,I6),1x,E20.10)') &
              'shotr:1 i_GA_indiv, i_parameter,  X_LMDIF', &
                       i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)

    enddo ! i_node =1, n_nodes
enddo !i_tree =1, n_trees


!----------------------------------------------------------------------------------------

! call fcn



!write(6,'(/A,3(1x,I10))') &
!      'shotr: call fcn, i_GA_indiv, n_time_steps, n_parameters     ', &
!                        i_GA_indiv, n_time_steps, n_parameters


iflag = 1

call fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag )

info = iflag


! write(6,'(A,4(1x,I10)/)') &
!           'shotr: aft call fcn i_GA_indiv, n_time_steps, n_parameters, info ', &
!                                i_GA_indiv, n_time_steps, n_parameters, info


!----------------------------------------------------------------------------------------



!if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!write(6,'(/A/ 1(1x, I6), 12( 1x,E12.5))') &
!      'shotr:3 i_GA_indiv, X_LMDIF', &
!               i_GA_indiv, X_LMDIF(1:n_parameters)

deallocate( x_LMDIF ) 
deallocate( fvec ) 



return

end subroutine show_tree
