subroutine read_input_data_size( )


use kinds_mod

use mpi
use mpi_module


use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module


IMPLICIT NONE


integer(kind=i4b) :: istat

!integer(kind=i4b), parameter :: data_unitnum  = 501
integer(kind=i4b), parameter :: line_length   = 250

CHARACTER(line_length) :: Aline


integer(kind=i4b) ::  ncount               
!integer(kind=i4b) ::  i                    
!integer(kind=i4b) ::  j                    




!----------------------------------------------------------------------


! START OF EXECUTABLE CODE

! open the control input file

open( unit = data_unitnum, file = 'GPGACODE_data', form = 'formatted',&
      status = 'old' )


rewind(data_unitnum)


!-------------------------------------------------------

! count number of data points

rewind(data_unitnum)

ncount = 0
do 
    read( data_unitnum, '(A)', iostat = istat ) Aline
    if( istat /= 0 ) exit

    ncount = ncount + 1
enddo

!-------------------------------------------------------


! subtract 1 because line 1 contains labels

n_input_data_points = ncount - 1

write(6, '(/A,1x,I7)') 'ris: n_input_data_points = ', n_input_data_points        
write(6, '(/A,1x,I7)') 'ris: n_input_vars        = ', n_input_vars               

return 

END subroutine read_input_data_size
