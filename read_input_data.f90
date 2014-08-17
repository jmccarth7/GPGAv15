subroutine read_input_data( )


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

!CHARACTER(line_length) :: Aline


integer(kind=i4b) ::  ncount               
integer(kind=i4b) ::  i                    
integer(kind=i4b) ::  j                    


real(kind=r8b), allocatable, dimension(:) ::  temp_array


!----------------------------------------------------------------------


allocate( temp_array( 0:n_input_vars ) ) 


!---------------------------------------------------------------------

! read data names and values

rewind(data_unitnum)


read( data_unitnum, *, iostat = istat ) input_data_names(0:n_input_vars) 

write(6, '(A)') 'rid:  input data names '

write(6, '(/A,1x,A)') 'rid: dependent variable = ', trim(input_data_names(0))

write(6, '(/A)') 'rid: independent variable i, input_data_names(i)  '

do  i = 1, n_input_vars
    write(6, '(I2,1x,A)') i, trim(input_data_names(i))
enddo

flush(6)

!---------------------------------------------------------------------

ncount = 0
do

    read( data_unitnum, *, iostat = istat ) temp_array(0:n_input_vars)
    if( istat /= 0 )exit

    ncount = ncount + 1

    input_data_array(0:n_input_vars, ncount ) = temp_array(0:n_input_vars) 

enddo

!---------------------------------------------------------------------

! echo input data

write(6,'(//A/)') 'rid:  input data '


!title_string = ' '
!do  i = 0, n_input_vars
!    title_string = title_string // input_data_names(i) 
!enddo

write(6,'(10(1x,A20))') ( trim( input_data_names(i) ), i = 0, n_input_vars )

do  j = 1, n_input_data_points 
    write(6,'(10(1x,E20.10))') &
            ( input_data_array(i,j), i = 0, n_input_vars ) 
enddo 

flush(6)


deallocate( temp_array ) 


return 

END subroutine read_input_data
