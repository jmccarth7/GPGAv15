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

real (kind=8) :: x_LMDIF(n_maximum_number_parameters)
real (kind=8) :: fvec(n_time_steps)


real(kind=8), parameter :: tol = 1.0d-30
real(kind=8), parameter :: epsfcn = 1.0d-6    ! original
real(kind=8), parameter :: factor=1.0D+0
real(kind=8), parameter :: zero = 0.0d0

integer (kind=4) ::   info


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(n_GA_individuals)


!real(kind=8) :: child_parameters(n_GA_individuals,n_maximum_number_parameters)
real(kind=8) :: child_parameters(n_GA_individuals,1:n_trees, 1:n_nodes )

external :: fcn

!--------------------------------------------------------------------------------------------

!write(6,'(/A,1(1x,I6)/)') 'shotr: i_GA_indiv ', i_GA_indiv

x_LMDIF(1:n_maximum_number_parameters) = 0.0D0

!do  i_parameter=1,n_parameters
i_parameter = 0
do  i_tree =1, n_trees
    do  i_node =1, n_nodes

        i_parameter = i_parameter + 1

        X_LMDIF(i_parameter) = child_parameters(i_GA_indiv,i_tree, i_node )
    
        !if( myid == 1 )then
            write(6,'(A,4(1x,I6),1x,E20.10)') &
                  'shotr:1 i_GA_indiv,i_tree, i_node,i_parameter, child_parameters ', &
                           i_GA_indiv,i_tree, i_node,i_parameter, &
                           child_parameters(i_GA_indiv,i_tree, i_node)
            write(6,'(A,2(1x,I6),1x,E20.10)') &
                  'shotr:1 i_GA_indiv, i_parameter,  X_LMDIF', &
                           i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
        !endif ! myid == 1

    enddo ! i_node =1, n_nodes
enddo !i_tree =1, n_trees


!----------------------------------------------------------------------------------------

! call fcn



!if( myid == 1 )then
!    write(6,'(/A,3(1x,I10))') &
!          'shotr: call fcn, i_GA_indiv, n_time_steps, n_parameters     ', &
!                            i_GA_indiv, n_time_steps, n_parameters
!endif ! myid == 1


iflag = 1

call fcn( n_time_steps, n_parameters, x_LMDIF, fvec, iflag )

info = iflag


!if( myid == 1 )then
!    write(6,'(A,4(1x,I10)/)') &
!              'shotr: aft call fcn i_GA_indiv, n_time_steps, n_parameters, info ', &
!                                   i_GA_indiv, n_time_steps, n_parameters, info
!endif ! myid == 1


!----------------------------------------------------------------------------------------



!if (info .eq. 8) info = 4

!-----------------------------------------------------------------------------------

!write(6,'(/A/ 1(1x, I6), 12( 1x,E12.5))') &
!      'shotr:3 i_GA_indiv, X_LMDIF', &
!               i_GA_indiv, X_LMDIF(1:n_parameters)


!do  i_parameter=1,n_parameters
!    child_parameters(i_GA_indiv,i_parameter) = &
!                            dabs( x_LMDIF(i_parameter) )
!    !write(6,'(A,2(1x,I6),1x,E20.10)') &
!    !      'shotr:3 aft RK  i_GA_indiv, i_parameter,  X_LMDIF', &
!    !                       i_GA_indiv, i_parameter,  X_LMDIF(i_parameter)
!enddo ! i_parameter

!write(6,'(/A/ 1(1x, I6), 12( 1x,E24.16))') &
!      'shotr:4  i_GA_indiv, child_parameters(i_GA_indiv,:)', &
!                i_GA_indiv, child_parameters(i_GA_indiv,1:n_parameters)

!-----------------------------------------------------------------------------------

!  calculate the individual SSE values by summing fvec over all time steps
!  fvec(i) = ( fcn(i) - truth(i) )**2
!  so SSE is calculated by summing fvec, not fvec**2

!write(6,'(/A/)')'shotr: calculate the individual SSE values '


!individual_SSE(i_GA_indiv)=0.0D+0
!
!if( individual_quality( i_GA_indiv ) > 0 ) then
!
!    !write(20,'(A,1x,I6)') 'shotr: i_GA_indiv ', i_GA_indiv
!    !write(6,'(A,1x,I6)') 'shotr: i_GA_indiv ', i_GA_indiv
!
!    do i_time_step=1,n_time_steps
!
!       if( isnan(fvec(i_time_step)) )    fvec(i_time_step) = 0.0d0
!       if( abs(fvec(i_time_step)) >  1.0d20 ) fvec(i_time_step) =  1.0d20
!
!       !write(20,'(A,1x,I6,1x,E24.16)' ) &
!       !      'shotr: i_time_step, fvec(i_time_step) ', &
!       !              i_time_step, fvec(i_time_step)
!       !write(6,'(A,1x,I6,1x,E24.16)' ) &
!       !      'shotr: i_time_step, fvec(i_time_step) ', &
!       !              i_time_step, fvec(i_time_step)
!
!       individual_SSE(i_GA_indiv) = individual_SSE(i_GA_indiv) + &
!                                    fvec(i_time_step)
!
!    enddo ! i_time_step
!
!endif !  individual_quality( i_GA_indiv ) > 0
!
!!write(6,'(A,2(1x,I6), 1x, E24.16)') &
!!      'shotr:  i_GA_indiv, individual_quality, individual_SSE', &
!!               i_GA_indiv, &
!!              individual_quality( i_GA_indiv ), &
!!              individual_SSE(i_GA_indiv)



return

end subroutine show_tree
