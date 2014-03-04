subroutine print_time_series_minSSE(  )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module


use GP_Parameters_module
use GA_Parameters_module
use GP_variables_module
use GA_Variables_module
use GP_Data_module


use Tree_Node_Factory_module
use class_Tree_Node


implicit none


!integer(kind=4),intent(in) :: i_GP_best_parent
!integer(kind=4),intent(in) :: nop
!integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_tree
integer(kind=4) :: i_node
integer(kind=4) :: ii
integer(kind=4) :: i
!integer(kind=4) :: jj

!logical :: buildTrees


real(kind=8), dimension( n_input_data_points ) :: resid
real(kind=8), dimension( n_input_data_points ) :: temp_data_array

!real(kind=8), intent(in),dimension( nop ) :: output_array


real(kind=8) :: RKmean, RKrms, RKstddev
real(kind=8) :: data_mean, data_rms, data_stddev
real(kind=8) :: resid_mean, resid_rms, resid_stddev

real(kind=8) :: RK_min
real(kind=8) :: data_min
real(kind=8) :: resid_min
real(kind=8) :: RK_max
real(kind=8) :: data_max
real(kind=8) :: resid_max

real(kind=8) :: resid_SSE
real(kind=8) :: r_corr
real(kind=8) :: prob_r
real(kind=8) :: fisher_z

integer, parameter :: plotMS_unit = 187

!------------------------------------------------------------------------------

GP_individual_Initial_Conditions = GP_minSSE_Individual_Initial_Conditions 
GP_Individual_Node_Parameters    = GP_minSSE_Individual_Node_Parameters
GP_Individual_Node_Type          = GP_minSSE_Individual_Node_Type


Numerical_CODE_Solution(0,1:n_CODE_equations)         = GP_individual_Initial_Conditions
Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = GP_individual_Initial_Conditions






!--------------------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(/A)') 'ptsMS: call Initialize_Model  '
endif ! myid == 0

call Initialize_Model( .true., .true., 6 )

!if( myid == 0 )then
!    write(6,'(/A/)') 'ptsMS: aft call Initialize_Model  '
!endif ! myid == 0

!------------------------------------------------------------------------------

! Generate PDF representation of trees


if( myid == 0 )then

    write(6,'(/A)') 'ptsMS: call Generate_Dot_Graph'

    call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, output_dir )

    !write(6,'(/A/)') 'ptsMS: aft call Generate_Dot_Graph'

endif ! myid == 0


! Write trees to disk

!if( myid == 0 )then
!    write(6,'(/A/)') 'ptsMS: call Serialize_Trees   '
!    call Serialize_Trees( GP_Trees(:,:), &
!                          n_Trees, n_Tracked_Resources, output_dir )
!    write(6,'(/A/)') 'ptsMS: aft call Serialize_Trees   '
!endif ! myid == 0


!------------------------------------------------------------------------------

! set the initial population node type using the info obtained
! from the setup file


Numerical_CODE_Solution(0,1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions
Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = &
                                 GP_individual_Initial_Conditions



if( myid == 0 )then

    write(6,'(A)') ' '
    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E24.16)') &
              'ptsMS: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                      ii, Numerical_CODE_Initial_Conditions(ii)
    enddo ! ii

    !write(6,'(A)') ' '

    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E24.16)') &
              'ptsMS: ii, Numerical_CODE_Solution(0,ii)         ', &
                      ii, Numerical_CODE_Solution(0,ii)
    enddo ! ii


    !write(6,'(A)') ' '
    write(6,'(/A,2(1x,I6))') 'ptsMS: n_trees, n_nodes ', n_trees, n_nodes

    write(6,'(/A)') &
          'ptsMS: i_tree  i_node  &
          &GP_Individual_Node_Parameters( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) == 0     )then

                write(6,'(2(1x,I8),6x,E24.16)') &
                      i_tree, i_node, &
                      GP_Individual_Node_Parameters( i_node, i_tree )

            endif ! GP_Individual_Node_Type( i_node, i_tree ) == 0

        enddo ! i_node
    enddo ! i_tree

    write(6,'(//A)') &
          'ptsMS: i_tree  i_node  &
          &GP_Individual_Node_Type( i_node, i_tree ) '

    do  i_tree = 1, n_trees
        do  i_node = 1, n_nodes

            if( GP_Individual_Node_Type( i_node, i_tree ) /= -9999 )then
                write(6,'(3(1x,I8))') &
                        i_tree, i_node, &
                        GP_Individual_Node_Type( i_node, i_tree )
            endif ! GP_Individual_Node_Type( i_node, i_tree ) /= -9999

        enddo ! i_node
    enddo ! i_tree

    write(6,'(A)') ' '

endif ! myid == 0



!------------------------------------------------------------------------


! run the Runge-Kutta model only once with proc 0

if( myid == 0 )then

    ! RK_Box_Model now put the time series in Numerical_CODE_Solution

    !call Runge_Kutta_Box_Model( .true. )  ! print
    call Runge_Kutta_Box_Model( .false. )   ! don't print



    !write(GP_print_unit,'(/A/)') &
    ! 'ptsMS: data_point   Numerical_CODE_Solution(data_point,1:n_CODE_equations)'
    !do  i = 0, n_input_data_points
    !    write(GP_print_unit,'(I6,2x,10(1x,E14.7))') &
    !          i, (Numerical_CODE_Solution(i,jj), jj = 1,n_CODE_equations )
    !enddo ! i

    open( plotMS_unit, file = 'plotMS.txt', status = 'unknown', &
          form = 'formatted', access = 'sequential' )


    write(GP_print_unit,'(/A/)') &
      'ptsMS:point   RK_Soln          input_data       &
       &diff = input_data_array - RK_Soln'

    write(plotMS_unit,'(A)') &
     '#ptsMS:point   RK_Soln          input_data       &
       &diff = input_data_array - RK_Soln'

    resid_SSE = 0.0d0
    do  i = 1, n_input_data_points

        resid(i) = input_data_array(0,i) -  Numerical_CODE_Solution(i,1)

        resid_SSE = resid_SSE + &
                   ( Data_Array(i,1) - Numerical_CODE_Solution(i,1) )**2  * &
                                                     Data_Variance_inv(1)


        write(GP_print_unit,'(I6,2x,3(2x,E15.7))') &
              i, Numerical_CODE_Solution(i,1),  input_data_array(0,i), &
                 input_data_array(0,i) - Numerical_CODE_Solution(i,1)

        write(plotMS_unit,'(I6,2x,3(2x,E15.7))') &
              i, Numerical_CODE_Solution(i,1),  input_data_array(0,i), &
                 input_data_array(0,i) - Numerical_CODE_Solution(i,1)
    enddo ! i


    call calc_stats( n_input_data_points, Numerical_CODE_Solution(1,1), &
                     RKmean, RKrms, RKstddev )


    temp_data_array = 0.0d0
    do  i = 1, n_input_data_points
        temp_data_array(i) = input_data_array(0,i)
    enddo

    call calc_stats( n_input_data_points, temp_data_array, &
                     data_mean, data_rms, data_stddev )


    call calc_stats( n_input_data_points, resid,                 &
                     resid_mean, resid_rms, resid_stddev )


    !call pearsn( Numerical_CODE_Solution(1,1), temp_data_array, &
    !             n_input_data_points, r_corr, prob_r, fisher_z )

    call corr( Numerical_CODE_Solution(1,1), temp_data_array, &
               n_input_data_points, 0, r_corr )

    write(GP_print_unit, '(//A,1x, I6,1x,E24.16/)') &
          'ptsMS: n_input_data_points, resid_SSE', &
                  n_input_data_points, resid_SSE

    RK_min = minval( Numerical_CODE_Solution(:,1) )
    RK_max = maxval( Numerical_CODE_Solution(:,1) )


    data_min = minval( temp_data_array(:) )
    data_max = maxval( temp_data_array(:) )

    resid_min = minval( resid(:) )
    resid_max = maxval( resid(:) )

    write(GP_print_unit, '(A)') &
          'ptsMS:            mean            rms             &
          &stddev            min            max'
    write(GP_print_unit, '(A,1x,5(1x,E15.7))') &
          'ptsMS: RK_Soln', RKmean, RKrms, RKstddev, RK_min, RK_max
    write(GP_print_unit, '(A,1x,5(1x,E15.7))') &
          'ptsMS: data   ', data_mean, data_rms, data_stddev, data_min, data_max
    write(GP_print_unit, '(A,1x,5(1x,E15.7)/)') &
          'ptsMS: resid  ', resid_mean, resid_rms, resid_stddev, resid_min, resid_max
    write(GP_print_unit, '(A,1x,5(1x,E15.7))') &
          'ptsMS: correlation coef. ', r_corr

    !write(GP_print_unit, '(A,1x,5(1x,E15.7))') &
    !      'ptsMS: R probability     ', prob_r
    !write(GP_print_unit, '(A,1x,5(1x,E15.7)/)') &
    !      'ptsMS: Fisher''s Z        ', fisher_z

    write(plotMS_unit, '(A)') &
          '#ptsMS:            mean            rms             &
          &stddev            min            max'
    write(plotMS_unit, '(A,1x,5(1x,E15.7))') &
          '#ptsMS: RK_Soln', RKmean, RKrms, RKstddev, RK_min, RK_max
    write(plotMS_unit, '(A,1x,5(1x,E15.7))') &
          '#ptsMS: data   ', data_mean, data_rms, data_stddev, data_min, data_max
    write(plotMS_unit, '(A,1x,5(1x,E15.7))') &
          '#ptsMS: resid  ', resid_mean, resid_rms, resid_stddev, resid_min, resid_max
    write(plotMS_unit, '(A,1x,5(1x,E15.7))') &
          '#ptsMS: correlation coef. ', r_corr

    !write(plotMS_unit, '(A,1x,5(1x,E15.7))') &
    !      '#ptsMS: R probability     ', prob_r
    !write(plotMS_unit, '(A,1x,5(1x,E15.7))') &
    !      '#ptsMS: Fisher''s Z        ', fisher_z

    close( plotMS_unit )

endif ! myid == 0


!--------------------------------------------------------------------------------

!if( allocated(  GP_Trees ) )then
!    do  i = 1, n_trees
!        call GP_Trees(i,1)%n%delete()
!        deallocate( GP_Trees(i,1)%n )
!    enddo
!    deallocate( GP_Trees )
!endif ! allocated( GP_Trees )

do  i = 1, n_trees
    call GP_Trees(i,1)%n%delete()
    deallocate( GP_Trees(i,1)%n )
enddo

!--------------------------------------------------------------------------------

return

end subroutine print_time_series_minSSE
