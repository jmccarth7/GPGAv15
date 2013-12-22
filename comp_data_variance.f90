subroutine comp_data_variance()

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
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none

real(kind=8) :: ssum, ssum2, totobs, dff
real(kind=8) :: totobs_m1

integer(kind=4) :: i_CODE_equation
integer(kind=4) :: i_time_step

integer(kind=4) :: nn

real(kind=8) :: sum1, sum2, sum3, mean, variance
real(kind=8) :: variance2
real(kind=8) :: M2, delta

!----------------------------------------------------------------------------------------



! compute the data_variance  -- to be used in computing SSE

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! there is some art form wiggle room to this calculation
! for instance, one can normalize by the variance of the
! individual observation types to give each observation
! equal weight, and there are other options that can be considered.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if( myid == 0 )then
    write(GP_print_unit,'(A)') ' '
endif

!write(GP_print_unit,'(A,3(1x,I6))') &
!         'cdv: myid, n_CODE_equations, n_time_steps ', &
!               myid, n_CODE_equations, n_time_steps

do i_CODE_equation=1,n_CODE_equations

    !-------------------------------------------------------------------------------
    ! original
    ssum  = 0.0D+0
    ssum2 = 0.0D+0

    do  i_time_step=0,n_time_steps
        ssum  = ssum  +  Data_Array(i_time_step,i_CODE_equation)
        ssum2 = ssum2 +  Data_Array(i_time_step,i_CODE_equation)**2
    enddo !   i_time_step

    !write(GP_print_unit,'(A,2(1x,I6),2(1x,E15.7) )') &
    !      'cdv: myid, i_CODE_equation, ssum, ssum2', &
    !            myid, i_CODE_equation, ssum, ssum2

    totobs    = dble(n_time_steps+1)
    totobs_m1 = dble(n_time_steps)

    !write(GP_print_unit,'(A,1x,I6, 2(1x,E15.7) )') &
    !      'cdv: myid, totobs, totobs_m1', &
    !            myid, totobs, totobs_m1


    dff=( (totobs*ssum2)-(ssum**2) ) / totobs / totobs_m1


    !write(GP_print_unit,'(/A,2x,E24.16)') 'cdv: original dff ', dff

    !write(GP_print_unit,'(A,1x, I6, 2(1x,E15.7) )') &
    !     'cdv: myid, dff', myid, dff

    !!-------------------------------------------------------------------------------
    !nn = n_time_steps + 1
    !sum1 = 0.0D0
    !do  i_time_step=0,n_time_steps
    !    sum1 = sum1 + Data_Array(i_time_step,i_CODE_equation)
    !enddo
    !
    !mean = sum1/real(nn, kind = 8 )
    !
    !sum2 = 0.0D0
    !sum3 = 0.0D0
    !do  i_time_step=0,n_time_steps
    !    sum2 = sum2 + (Data_Array(i_time_step,i_CODE_equation) - mean)**2
    !    sum3 = sum3 + (Data_Array(i_time_step,i_CODE_equation) - mean)
    !enddo
    !
    !variance = ( sum2 - sum3**2/real(nn, kind=8) ) / real(nn-1, kind=8)
    !
    !write(GP_print_unit,'(/A,2x,E24.16)') 'cdv: 1 variance     ', variance
    !
    !!-------------------------------------------------------------------------------
    !nn = 0
    !mean = 0.0d0
    !M2 = 0.0d0
    !
    !
    !do  i_time_step=0,n_time_steps
    !    nn = nn + 1
    !    delta = Data_Array(i_time_step,i_CODE_equation) - mean
    !    mean = mean + delta / real(nn, kind=8)
    !    M2 = M2 + delta*(Data_Array(i_time_step,i_CODE_equation) - mean)
    !enddo
    !
    !variance2 = M2/real(nn - 1, kind=8)
    !
    !write(GP_print_unit,'(/A,2x,E24.16)') 'cdv: 2 variance2    ', variance2
    !
    !!-------------------------------------------------------------------------------

    if( dff .gt. 0.0D+0) then

        ! set variance to observed variance for normalize by the s.d.

        Data_Variance(i_CODE_equation)=dff

    else

        ! set variance to 1.0 for normalization to be 'unaltered'

        Data_Variance(i_CODE_equation)=1.0D+0

    endif !   dff .gt. 0.0D+0



    ! compute data_variance_inv

    if( abs( Data_Variance(i_CODE_equation) ) > 0.0D0 ) then
        Data_Variance_inv(i_CODE_equation) = 1.0D0 / Data_Variance(i_CODE_equation)
    endif ! abs( data_var...


    if( abs( Data_Variance(i_CODE_equation) ) < 1.0D-30 )then

        write(GP_print_unit,'(/A,1x,I6,2x,E15.7)') &
         'cdv: i_CODE_equation, Data_Variance(i_CODE_equation) ', &
               i_CODE_equation, Data_Variance(i_CODE_equation)

        write(GP_print_unit,'(A/)') 'cdv: bad data variance -- stopping program '
        call MPI_FINALIZE(ierr)
        stop 'bad data var'

    endif ! abs( Data_Variance(i_CODE_equation) ) < 1.0D-30

    if( abs( Data_Variance_inv(i_CODE_equation) ) <= 0.0D0  )then

        write(GP_print_unit,'(/A,1x,I6,2x,E15.7)') &
         'cdv: i_CODE_equation, Data_Variance_inv(i_CODE_equation) ', &
               i_CODE_equation, Data_Variance_inv(i_CODE_equation)

        write(GP_print_unit,'(A/)') 'cdv: bad data variance inv -- stopping program '
        call MPI_FINALIZE(ierr)
        stop 'bad data var_inv'

    endif ! abs( Data_Variance_inv(i_CODE_equation) ) <=0.0D0


    if( myid == 0 )then

        write(GP_print_unit,'(A,1x,I4,2(1x,E23.16))') &
             'cdv: i_CODE_eq, Data_Variance, Data_Variance_inv ', &
                   i_CODE_equation, Data_Variance(i_CODE_equation), &
                                  Data_Variance_inv(i_CODE_equation)

    endif ! myid == 0

enddo !  i_CODE_equation

if( myid == 0 )then
    write(GP_print_unit,'(A)') ' '
endif !  myid == 0


return

end subroutine comp_data_variance
