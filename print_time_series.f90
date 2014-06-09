subroutine print_time_series( i_GP_best_parent,  nop )

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

character(1000) :: title_string

integer(kind=4),intent(in) :: i_GP_best_parent
integer(kind=4),intent(in) :: nop
integer(kind=4) :: i_tree
integer(kind=4) :: i_node
integer(kind=4) :: ii
integer(kind=4) :: i
integer(kind=4) :: j
integer(kind=4) :: jj



real(kind=8) :: x_time_step


!real(kind=8), dimension( n_input_data_points, n_code_equations ) :: resid
real(kind=8), dimension( n_time_steps, n_code_equations ) :: resid

!real(kind=8), dimension( n_input_data_points ) :: temp_data_array

real(kind=8), dimension( n_time_steps ) :: temp_data_array

!real(kind=8), intent(in),dimension( nop ) :: output_array


real(kind=8),dimension(n_code_equations)  :: RKmean
real(kind=8),dimension(n_code_equations)  :: RKrms
real(kind=8),dimension(n_code_equations)  :: RKstddev
real(kind=8),dimension(n_code_equations)  :: data_mean
real(kind=8),dimension(n_code_equations)  :: data_rms
real(kind=8),dimension(n_code_equations)  :: data_stddev
real(kind=8),dimension(n_code_equations)  :: resid_mean
real(kind=8),dimension(n_code_equations)  :: resid_rms
real(kind=8),dimension(n_code_equations)  :: resid_stddev
real(kind=8),dimension(n_code_equations)  :: RK_min
real(kind=8),dimension(n_code_equations)  :: data_min
real(kind=8),dimension(n_code_equations)  :: resid_min
real(kind=8),dimension(n_code_equations)  :: RK_max
real(kind=8),dimension(n_code_equations)  :: data_max
real(kind=8),dimension(n_code_equations)  :: resid_max
real(kind=8),dimension(n_code_equations)  :: r_corr
real(kind=8),dimension(n_code_equations)  :: prob_r
real(kind=8),dimension(n_code_equations)  :: fisher_z

real(kind=8) :: resid_SSE
real(kind=8) :: y_min
real(kind=8) :: y_max
integer, parameter :: plot_unit = 177

!------------------------------------------------------------------------------

GP_individual_Initial_Conditions = GP_Population_Initial_Conditions(:, i_GP_best_parent)
GP_Individual_Node_Parameters    = GP_population_node_parameters(:,:,i_GP_best_parent)
GP_Individual_Node_Type          = GP_Adult_Population_Node_Type(:,:,i_GP_best_parent)


Numerical_CODE_Solution(0,1:n_CODE_equations)         = GP_individual_Initial_Conditions
Numerical_CODE_Initial_Conditions(1:n_CODE_equations) = GP_individual_Initial_Conditions

Numerical_CODE_Solution(1:n_time_steps,1:n_CODE_equations) = 0.0d0


if( myid == 0 )then
    write(6,'(/A,2(1x,I6)/)') 'pts: nop',  nop
endif ! myid == 0


!--------------------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(/A)') 'pts: call Initialize_Model  '
endif ! myid == 0

call Initialize_Model( .true., .true., 6 )

!if( myid == 0 )then
!    write(6,'(/A/)') 'pts: aft call Initialize_Model  '
!endif ! myid == 0

!------------------------------------------------------------------------------

! Generate PDF representation of trees


if( myid == 0 )then

    write(6,'(/A)') 'pts: call Generate_Dot_Graph'

    call Generate_Dot_Graph( GP_Trees(:,1), n_Trees, './pts')

    write(6,'(/A/)') 'pts: aft call Generate_Dot_Graph'

endif ! myid == 0


!------------------------------------------------------------------------------

! Write trees to disk

!if( myid == 0 )then
!    write(6,'(/A/)') 'pts: call Serialize_Trees   '
!    call Serialize_Trees( GP_Trees(:,:), &
!                          n_Trees, n_Tracked_Resources, output_dir )
!    write(6,'(/A/)') 'pts: aft call Serialize_Trees   '
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
              'pts: ii, Numerical_CODE_Initial_Conditions(ii) ', &
                    ii, Numerical_CODE_Initial_Conditions(ii)
    enddo ! ii

    write(6,'(A)') ' '

    do  ii = 1, n_CODE_equations
        write(6,'(A,1x,I6,1x,E24.16)') &
              'pts: ii, Numerical_CODE_Solution(0,ii)         ', &
                    ii, Numerical_CODE_Solution(0,ii)
    enddo ! ii


    write(6,'(/A,2(1x,I6))') 'pts: n_trees, n_nodes ', n_trees, n_nodes

    write(6,'(/A)') &
          'pts: i_tree  i_node  &
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
          'pts: i_tree  i_node  &
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


    ! RK_Box_Model now puts the time series in Numerical_CODE_Solution


    call Runge_Kutta_Box_Model( .true. )  ! print
    !call Runge_Kutta_Box_Model( .false. )   ! don't print



    write(GP_print_unit,'(/A/)') &
     'pts: data_point   Numerical_CODE_Solution(data_point,1:n_CODE_equations)'
    do  i = 1, n_time_steps   !n_input_data_points
        write(GP_print_unit,'(I6,2x,10(1x,E14.7))') &
              i, (Numerical_CODE_Solution(i,jj), jj = 1,n_CODE_equations )
    enddo ! i



    open( plot_unit, file = 'plot.txt', status = 'unknown', &
          form = 'formatted', access = 'sequential' )


    title_string = '#pts: pt'
    title_string = trim( title_string ) // &
                       '   RK_Soln      input_data  resid'
    do  j = 2, n_code_equations
        title_string = trim( title_string ) // &
                       '         RK_Soln      input_data  resid'
    enddo



    write(GP_print_unit,'(/A/)')  trim( title_string )
    write(plot_unit,'(A)')        trim( title_string )


    !------------------------------------------------------------------------------------

    ! calculate the resid_SSE only for times between sse_min_time and sse_max_time

    resid_SSE = 0.0d0

    do  i = 1, n_time_steps   !  n_input_data_points

        x_time_step = real( i, kind=8 ) * dt

        if( x_time_step < sse_min_time ) then
            sse_wt = sse_low_wt
        else
            sse_wt = 1.0d0      
        endif ! x_time_step < sse_min_time

        if( x_time_step > sse_max_time ) exit

        do  j = 1, n_code_equations
            resid_SSE = resid_SSE + &
                       ( Data_Array(i,j) - Numerical_CODE_Solution(i,j) )**2  * &
                                                     Data_Variance_inv(j) * &
                                                     sse_wt
            write(GP_print_unit,'(A,1x,I6,1x,I3,1x,5(1x,E15.7))') &
                  'pts: i, j, Soln, Data, Data_Var_inv, sse_wt, SSE', &
                  i, j, Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                       Data_Variance_inv(j), sse_wt, resid_SSE

        enddo ! j

    enddo ! i

    !------------------------------------------------------------------------------------

    do  i = 1, n_time_steps   !  n_input_data_points

        do  j = 1, n_code_equations

            resid(i,j) = Data_Array(i,j) -  Numerical_CODE_Solution(i,j)

        enddo ! j

        write(GP_print_unit,'(I6,2x,50(1x,E12.5))') &
              i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                                  j = 1, n_code_equations )

        write(plot_unit, '(I6,2x,50(1x,E12.5))') &
              i, ( Numerical_CODE_Solution(i,j),  Data_Array(i,j), &
                   Data_Array(i,j) - Numerical_CODE_Solution(i,j), &
                                                  j = 1, n_code_equations )

    enddo ! i



    !--------------------------------------------------------------------------------

    do  j = 1, n_code_equations

        call calc_stats( n_time_steps,  Numerical_CODE_Solution(1,j), &
                         RKmean(j), RKrms(j), RKstddev(j) , &
                         1.0d0, 0.0d0, 1.0d9, 1.0d0 ) 


        call calc_stats( n_time_steps, Data_Array(1,j), &
                         data_mean(j), data_rms(j), data_stddev(j), &
                         1.0d0, 0.0d0, 1.0d9, 1.0d0 ) 


        call calc_stats( n_time_steps, resid(1,j) ,              &
                         resid_mean(j), resid_rms(j), resid_stddev(j), &
                         1.0d0, 0.0d0, 1.0d9, 1.0d0 ) 


        !call pearsn( Numerical_CODE_Solution(1,1), temp_data_array, &
        !             n_input_data_points, r_corr, prob_r, fisher_z )

        call corr( Numerical_CODE_Solution(1,j), Data_Array(1,j), &
                   n_time_steps, 0, r_corr(j) , &
                   1.0d0, 0.0d0, 1.0d9, 1.0d0 ) 


        !RK_min(j) = minval( Numerical_CODE_Solution(:,j) )
        !RK_max(j) = maxval( Numerical_CODE_Solution(:,j) )


        !data_min(j) = minval( Data_Array(:,j) )
        !data_max(j) = maxval( Data_Array(:,j) )

        !resid_min(j) = minval( resid(:,j) )
        !resid_max(j) = maxval( resid(:,j) )

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts: j, RK_min(j), RK_max(j)      ', &
        !            j, RK_min(j), RK_max(j)

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts: j, data_min(j), data_max(j)  ', &
        !            j, data_min(j), data_max(j)

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts: j, resid_min(j), resid_max(j)', &
        !            j, resid_min(j), resid_max(j)

    enddo ! j

    !--------------------------------------------------------------------------------

    do  j = 1, n_code_equations

        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = Numerical_CODE_Solution(i,j)
        enddo ! i

        rk_min(j) =  minval( temp_data_array )
        rk_max(j) =  maxval( temp_data_array )

        !write(GP_print_unit,'(A,1x,I3,2x,2(1x,E12.5))') &
        !      'pts:2 j, RK_min(j), RK_max(j)      ', &
        !             j, RK_min(j), RK_max(j)


        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = data_array(i,j)
        enddo ! i

        data_min(j) =  minval( temp_data_array )
        data_max(j) =  maxval( temp_data_array )

        !write(GP_print_unit,'(A,1x,I3,2x,2(1x,E12.5))') &
        !      'pts:2 j, data_min(j), data_max(j)  ', &
        !             j, data_min(j), data_max(j)



        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = resid(i,j)
        enddo ! i

        resid_min(j) =  minval( temp_data_array )
        resid_max(j) =  maxval( temp_data_array )

        !write(GP_print_unit,'(A,1x,I3,2x,2(1x,E12.5))') &
        !      'pts:2 j, resid_min(j), resid_max(j)', &
        !             j, resid_min(j), resid_max(j)


        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts: j, RK_min(j), RK_max(j)      ', &
        !            j, RK_min(j), RK_max(j)

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts: j, data_min(j), data_max(j)  ', &
        !            j, data_min(j), data_max(j)

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts: j, resid_min(j), resid_max(j)', &
        !            j, resid_min(j), resid_max(j)

    enddo ! j

    !--------------------------------------------------------------------------------

    do  j = 1, n_code_equations

        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = Numerical_CODE_Solution(i,j)
        enddo ! i

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts:2 j, RK_min(j), RK_max(j)      ', &
        !             j, minval( temp_data_array ), maxval( temp_data_array )


        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = data_array(i,j)
        enddo ! i

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts:2 j, data_min(j), data_max(j)      ', &
        !             j, minval( temp_data_array ), maxval( temp_data_array )



        temp_data_array = 0.0d0
        do  i = 1, n_time_steps
            temp_data_array(i) = resid(i,j)
        enddo ! i

        !write(GP_print_unit,'(A,1x,I6,2x,2(1x,E12.5))') &
        !      'pts:2 j, resid_min(j), resid_max(j)      ', &
        !             j, minval( temp_data_array ), maxval( temp_data_array )


    enddo ! j

    !--------------------------------------------------------------------------------

    ! calculate overall y_min and y_max for plotting

    y_min =  1.0d99
    y_max = -1.0d99

    do  j = 1, n_code_equations

        y_min = min( y_min, RK_min(j) )
        y_max = max( y_max, RK_max(j) )

        y_min = min( y_min, data_min(j) )
        y_max = max( y_max, data_max(j) )

    enddo ! j

    !--------------------------------------------------------------------------------

    !write(GP_print_unit, '(//A,1x, I6,1x,E24.16/)') &
    !     'pts: n_input_data_points, resid_SSE', &
    !           n_input_data_points, resid_SSE

    write(GP_print_unit, '(//A,1x, I6,1x,E24.16/)') &
         'pts: n_time_steps, resid_SSE', &
               n_time_steps, resid_SSE

    !--------------------------------------------------------------------------------

    ! print results

    do  j = 1, n_code_equations

        write(GP_print_unit, '(/A)') &
              'pts: i_code_eq           mean            rms             &
              &stddev            min            max'
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'pts: RK_Soln', &
              j, RKmean(j), RKrms(j), RKstddev(j), RK_min(j), RK_max(j)
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'pts: data   ', &
              j, data_mean(j), data_rms(j), data_stddev(j), data_min(j), data_max(j)
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7)/)') &
              'pts: resid  ', &
              j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min(j), resid_max(j)
        write(GP_print_unit, '(A,1x,I2, 5(1x,E15.7))') &
              'pts: corr coef. ', j, r_corr(j)

    enddo ! j

    write(GP_print_unit, '(/A,1x,E15.7)') 'pts: y_min', y_min
    write(GP_print_unit, '(A,1x,E15.7/)') 'pts: y_max', y_max

    !write(GP_print_unit, '(A,1x,5(1x,E15.7))') &
    !      'pts: R probability     ', prob_r
    !write(GP_print_unit, '(A,1x,5(1x,E15.7)/)') &
    !      'pts: Fisher''s Z        ', fisher_z

    !flush(GP_print_unit) 


    !--------------------------------------------------------------------------------

    !  write results to output file

    do  j = 1, n_code_equations

        write(plot_unit, '(A)') &
              '#pts:  i_code_eq          mean            rms             &
              &stddev            min            max'
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: RK_Soln', &
              j, RKmean(j), RKrms(j), RKstddev(j), RK_min(j), RK_max(j)
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: data   ', &
              j, data_mean(j), data_rms(j), data_stddev(j), data_min(j), data_max(j)
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: resid  ', &
              j, resid_mean(j), resid_rms(j), resid_stddev(j), resid_min(j), resid_max(j)
        write(plot_unit, '(A,1x,I2,5(1x,E15.7))') &
              '#pts: corr coef. ', j, r_corr(j)

    enddo ! j

    write(plot_unit, '(A,1x,E15.7)')  '#pts: y_min', y_min
    write(plot_unit, '(A,1x,E15.7)')  '#pts: y_max', y_max

    !write(plot_unit, '(A,1x,5(1x,E15.7))') &
    !      '#pts: R probability     ', prob_r
    !write(plot_unit, '(A,1x,5(1x,E15.7))') &
    !      '#pts: Fisher''s Z        ', fisher_z

    !flush(plot_unit)

    close( plot_unit )

endif ! myid == 0


!--------------------------------------------------------------------------------

do  i = 1, n_trees
    if( associated( GP_Trees(i,1)%n ) ) then
        call GP_Trees(i,1)%n%delete()
        deallocate( GP_Trees(i,1)%n )
    endif !  associated( GP_Trees(i,1)%n )
enddo ! i

!--------------------------------------------------------------------------------

return

end subroutine print_time_series
