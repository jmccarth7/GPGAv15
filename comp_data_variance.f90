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
use Runge_Kutta_Variables_module


implicit none

real (kind=8) :: ssum, ssum2, totobs, dff
!integer(kind=4) :: i_GP_individual                                                                           
!integer(kind=4) :: i_GP_Best_Parent                                                                          
!integer(kind=4) :: i_GP_Generation                                                                           
!integer(kind=4) :: j_GP_Individual                                                                           
!integer(kind=4) :: i_Tree                                                                                    
!integer(kind=4) :: i_Node            

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

do i_CODE_equation=1,n_CODE_equations

  ssum  = 0.0D+0
  ssum2 = 0.0D+0

  do  i_time_step=0,n_time_steps
      ssum  = ssum  +  Data_Array(i_time_step,i_CODE_equation)
      ssum2 = ssum2 +  Data_Array(i_time_step,i_CODE_equation)**2
  enddo !   i_time_step

  totobs=dble(n_time_steps+1)

  dff=( (totobs*ssum2)-(ssum**2) ) / ( totobs*(totobs-1.0D+0) )

  if( dff .gt. 0.0D+0) then  ! set variance to observed variance for normalize by the s.d.
      Data_Variance(i_CODE_equation)=dff
  else ! set variance to 1.0 for normalization to be 'unaltered'
      Data_Variance(i_CODE_equation)=1.0D+0
  endif !   dff .gt. 0.0D+0

  if(  abs( Data_Variance(i_CODE_equation) ) < 1.0D-30 )then
       write(GP_print_unit,'(/A,1x,I6,2x,E15.7)') &
       '0: i_CODE_equation, Data_Variance(i_CODE_equation) ', &
           i_CODE_equation, Data_Variance(i_CODE_equation)
       write(GP_print_unit,'(A/)') '0: bad data variance -- stopping program '
       stop 'bad data var'

  endif ! abs( Data_Variance(i_CODE_equation) ) < 1.0D-30

  if( myid == 0 )then
      write(GP_print_unit,'(A,1x,I6,2x,E15.7)') &
           '0: i_CODE_equation, Data_Variance(i_CODE_equation) ', &
               i_CODE_equation, Data_Variance(i_CODE_equation)
  endif ! myid == 0

enddo !  i_CODE_equation

if( myid == 0 )then
    write(GP_print_unit,'(A)') ' '
endif !  myid == 0


return

end subroutine comp_data_variance
