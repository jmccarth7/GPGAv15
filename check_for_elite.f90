subroutine check_for_elite( index0  )

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

integer(kind=4) :: index0
integer(kind=4) :: ksafe

real(kind=4) :: cff
real(kind=8) :: dff


!---------------------------------------------------------------------------

!  generate indices for the tournament selection

!  keep generating numbers until one is found
!  which is not an index of an elite individual which must not be replaced

!  ksafe is used to prevent infinite loops

ksafe = 0

do

    ksafe = ksafe + 1

    if( ksafe > n_GA_individuals ) then
        write(GA_print_unit,'(A,2(1x,I6))') &
              'cfe: no good index found  ksafe, n_GA_individuals ', &
                                         ksafe, n_GA_individuals
        stop 'check_elite bad'
    endif

    call random_number(cff) ! uniform random number generator
    dff = cff

    index0  = 1 + int(  dff * real( n_GA_Individuals-1, kind=8 )  )


    !write(GA_print_unit,'(A,2(1x,I6))')    'cfe: ksafe, index0 ', ksafe, index0
    !write(GA_print_unit,'(A/(15(1x,I6)))') 'cfe: individual_elites ',  &
    !                                 individual_elites(1:n_GA_save_elites)

    if( any( individual_elites == index0 ) )then

        !write(GA_print_unit,'(A,1x,I6)') 'cfe: index is elite       ', index0
        cycle

    endif   ! any( individual_elites == index0 )


    if( .not. any( individual_elites == index0 ) ) exit

enddo

!write(GA_print_unit,'(A,1x,I6)') 'cfe: at return index0     ', index0


return
end subroutine check_for_elite
