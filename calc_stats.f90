subroutine calc_stats( n_array, array, mean, rms, std_dev )

implicit none

integer, intent(in)   :: n_array
real(kind=8),intent(in), dimension( n_array ) :: array

real(kind=8) :: mean
real(kind=8) :: rms
real(kind=8) :: std_dev
real(kind=8) :: sum1
real(kind=8) :: sum2

integer :: i
integer :: icount

real(kind=8) :: xcount

!-------------------------------------------------------------------

mean = 0.0d0
rms  = 0.0d0
std_dev  = 0.0d0
sum1     = 0.0d0
sum2     = 0.0d0

icount = 0
do  i = 1, n_array

    sum1 = sum1 + array(i)
    sum2 = sum2 + array(i)**2
    icount = icount + 1

enddo

if( icount == 0 )then

    mean = 0.0d0
    rms  = 0.0d0
    std_dev = 0.0d0
    return

endif ! icount == 0

xcount = real( icount, kind = 8 )

mean = sum1 / xcount
rms  = sum2 / xcount

std_dev =  sqrt( abs( rms   - mean**2 ) )

rms  = sqrt( rms )

!write(6,'(A,3(1x,E24.16))') 'cs: xcount, sum1, sum2', &
!                                xcount, sum1, sum2
!write(6,'(A,3(1x,E24.16))') 'cs: mean, rms, std_dev', &
!                                mean, rms, std_dev
return

end subroutine calc_stats
