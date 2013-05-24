subroutine GA_replace_bad_individuals( Child_Parameters, &
                                       individual_quality  )


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module

implicit none

real(kind=8)    :: child_parameters(n_GA_Individuals,n_maximum_number_parameters)
integer(kind=4) :: individual_quality(n_GA_individuals)

real(kind=4) :: cff
real(kind=8) :: dff

integer(kind=4) :: n_replaced
integer(kind=4) :: i
integer(kind=4) :: i_parameter


!----------------------------------------------------------------------------


! for each individual, i, which has a quality < 0,
! generate a new set of random numbers for the parameters, and
! set the run_ga_lmdif to .true. so that this individual will
! be processed by fcn later.


n_replaced  = 0

i_loop:&
do  i_GA_Individual = 1, n_GA_Individuals

    !write(GA_print_unit,'(A,1x,I6,1x,I6)') &
    !      'grbi: i_GA_Individual, individual_quality', &
    !             i_GA_Individual, individual_quality( i_GA_Individual )

    if( individual_quality( i_GA_Individual ) < 0 )then

        do  i_Parameter = 1, n_Parameters

            call random_real(cff) ! random real number generator
            dff = cff
            Child_Parameters(i_GA_Individual,i_Parameter) = dff


            !----------------------------------------------------------------------------

!!            if( n_linked_parms > 0 )then
!!
!!                if( any( i_parameter == linked_parms(2,:) )    ) cycle
!!
!!                if( any( i_parameter == linked_parms(1,:) )    )then
!!
!!                    do  i = 1, n_linked_parms
!!
!!                        if( i_parameter == linked_parms(1,i) )then
!!
!!                            child_parameters(i_GA_Individual, linked_parms(2,i) ) = &
!!                            child_parameters(i_GA_Individual, linked_parms(1,i) )
!!
!!                            exit
!!
!!                        endif ! i_parameter == linked_parms(i,1)
!!
!!                    enddo ! i
!!
!!                endif !  any( i_parameter == linked_parms(1,:) )
!!
!!            endif !  n_linked_parms > 0

            !----------------------------------------------------------------------------



            !write(GA_print_unit,'(A,2(1x,I6),1x,E24.16)') &
            !'grbi: i_GA_Individual, i_parameter, Child_Parameters', &
            !       i_GA_Individual, i_parameter, &
            !       Child_Parameters(i_GA_Individual,i_Parameter)

        enddo ! i_Parameter


        n_replaced = n_replaced + 1

        Run_GA_lmdif(i_GA_Individual) = .true.

        ! don't set quality here since this will cause the fitness means computed in
        ! the fitness proportionate subroutine to be wrong since the
        ! replaced parameters will not yet have a fitness value

        ! individual_quality will be set to 1 for all individuals before the next RK run

        !!!!individual_quality(i_GA_Individual) = 1

    endif ! individual_quality( i_GA_Individual ) < 0

enddo i_loop  ! i_GA_Individual

write(GA_print_unit,'(A,1x,I6,1x,i10/)') 'grbi: i_GA_generation, n_replaced ', &
                                                i_GA_generation, n_replaced

return


end subroutine GA_replace_bad_individuals
