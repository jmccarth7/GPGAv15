subroutine random_real(bff)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This code generates a random number ranging between the full range of
! numbers that can be represented in single precision and the numbers
! are generated "uniformly over a log10 scale"
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

implicit none

!real (kind=4) :: aff,bff,cff,dff
real (kind=4) ::     bff,cff
!real (kind=4) :: log_huge,log_tiny,log_range

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! uniform random choice over log scale of full range of +/- real*4 values
! aff=1.
! call random_number(bff) ! uniform random number generator
! bff=bff-0.5
! aff=sign(aff,bff)
! call random_number(bff) ! uniform random number generator
! log_huge=log10(huge(dff))
! log_tiny=log10(tiny(dff))
! log_range=log_huge-log_tiny
! cff=aff*(10.**((bff*log_range)+log_tiny))

!-----------------------------------------------------------------
!call random_number(bff) ! uniform random number generator

!bff=50.*bff
!-----------------------------------------------------------------

! code below attempts to get a better balance of random numbers

! choose a random number between 0 and 1

! if <  0.3,  get another random number in the range [0.,1.]
! if >= 0.3,  get another random number in the range [0.,50.]

!-----------------------------------------------------------------

call random_number(cff) ! uniform random number generator

if( cff <= 0.2 )then

    call random_number(bff) ! uniform random number generator

    bff = 1.0 * bff

else

    call random_number(bff) ! uniform random number generator

    bff = 10. * bff

endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end
!234567890123456789012345678901234567890123456789012345678901234567890
