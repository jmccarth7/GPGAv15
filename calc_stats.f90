subroutine calc_stats( n_array, array, mean, rms, std_dev, &
                       dt, sse_min_time, sse_max_time, sse_low_wt  )

implicit none

integer, intent(in)   :: n_array
real(kind=8),intent(in), dimension( n_array ) :: array

real(kind=8) :: mean
real(kind=8) :: rms
real(kind=8) :: std_dev
real(kind=8) :: sum1
real(kind=8) :: sum2
real(kind=8) :: arr 
real(kind=8),intent(in) :: dt
real(kind=8),intent(in) :: sse_min_time
real(kind=8),intent(in) :: sse_max_time
real(kind=8),intent(in) :: sse_low_wt   

integer :: i
integer :: icount

real(kind=8) :: xcount

real(kind=8) :: xi       

!-------------------------------------------------------------------

mean = 0.0d0
rms  = 0.0d0
std_dev  = 0.0d0
sum1     = 0.0d0
sum2     = 0.0d0

!write(6,'(A,3(1x,E15.7))') &
!  'cs: dt, sse_min_time, sse_max_time', &
!       dt, sse_min_time, sse_max_time

icount = 0
do  i = 1, n_array

    xi = real( i, kind=8) * dt
    if( xi < sse_min_time )then
        arr = array(i) * sse_low_wt
    else
        arr = array(i) 
    endif 
    if( xi > sse_max_time ) exit

    
    !sum1 = sum1 + array(i)
    !sum2 = sum2 + array(i)**2

    sum1 = sum1 + arr
    sum2 = sum2 + arr**2
    icount = icount + 1

    !write(6,'(A,1x,F10.2,2(1x,I6),3(1x,E15.7))') &
    !  'cs: xi, i, icount, array(i), arr, sum1', &
    !       xi, i, icount, array(i), arr, sum1

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

!write(6,'(A,3(1x,E15.7))') 'cs: xcount, sum1, sum2', &
!                                xcount, sum1, sum2
!write(6,'(A,3(1x,E15.7))') 'cs: mean, rms, std_dev', &
!                                mean, rms, std_dev

return

end subroutine calc_stats
