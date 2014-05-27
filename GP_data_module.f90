module GP_data_module

use GP_parameters_module

implicit none

!real(kind=8) :: Data_Array(0:n_time_steps,n_CODE_equations)
real(kind=8),allocatable, dimension(:,:) :: Data_Array


real(kind=8),allocatable, dimension( : )  :: Data_Variance_inv

real(kind=8),allocatable, dimension( : )  :: ratio_Data_Variance_inv


end module GP_data_module
