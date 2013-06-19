subroutine calc_fitness( child_parameters, individual_quality, &
                         i_GA_Best_Parent, parent_parameters,  &
                         L_stop_run  )


! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module


implicit none


integer,parameter ::  itag = 1


real(kind=8) :: parent_parameters(n_GA_individuals,n_maximum_number_parameters)
real(kind=8) :: child_parameters(n_GA_individuals,n_maximum_number_parameters)


real(kind=8) :: dble_cff

integer (kind=4) ::    i_GA_Best_Parent
integer (kind=4) ::    n_counted
integer (kind=4) ::    index_min_sse
integer (kind=4) ::    icount

real(kind=8), parameter :: max_err = 1.0d8  !100.0d0

real(kind=8) :: edit_level
real(kind=8) :: mean_individual_fit
real(kind=8) :: mean_individual_SSE
real(kind=8) :: mean_fitness

real(kind=8) :: xn
real(kind=8) :: var_fitness
real(kind=8) :: sigma_fitness


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(n_GA_individuals)


!integer(kind=4) :: n_retry

external :: fcn

real(kind=8), external :: indiv_fitness

logical :: L_stop_run

integer(kind=4) :: jj
integer(kind=4) :: i_parameter

!----------------------------------------------------------------------------------

L_stop_run = .FALSE.

!write(GA_print_unit,'(/A)') 'calcfit:  finished the loop on i_GA_individual  '

do  i_parameter=1,n_parameters
    do  i_GA_individual=1,n_GA_individuals

        parent_parameters(i_GA_individual,i_parameter) = &
         child_parameters(i_GA_individual,i_parameter)

    enddo !  i_GA_individual
enddo ! i_parameter

!-----------------------------------------------------------------------------------

!write(GA_print_unit,'(/A)') 'calcfit: i_GA_individual, parent params'

!do  i_GA_individual=1,n_GA_individuals
!    write(GA_print_unit,'(I6,(12(1x,E15.7)))') i_GA_individual, &
!              parent_parameters( i_GA_individual, 1:n_parameters)
!enddo !  i_GA_individual

!-----------------------------------------------------------------------------------

!write(GA_print_unit,*) ' '
!do  i_GA_individual=1,n_GA_individuals
!    write(GA_print_unit,'(A,1x,I6,1x,E24.16, 1x, I6)') &
!        'calcfit: i_GA_individual, individual_SSE, individual_quality', &
!                  i_GA_individual, individual_SSE(i_GA_individual), &
!                                   individual_quality( i_GA_individual )
!enddo !  i_GA_individual

!-----------------------------------------------------------------------------------


! compute the edit level = n_time_steps *  max_err**2


! that is, the edit level is reached
! if the residuals on all time steps are equal
! to the maximum, user-specified, error


edit_level = real(n_time_steps,kind=8) * max_err**2

write(GA_print_unit,'(/A,1x,I6, 2(1x, E15.7))') &
      'calcfit: n_time_steps, max_err, edit_level ', &
                n_time_steps, max_err, edit_level

!-----------------------------------------------------------------------------------

! if the indiv_sse > sse0,
! then this is a bad result, so give individual_quality = -1

! only do this if individual_quality > 0, since
! if the RK process was bad, individual_quality  was set to -10

! edit the individual_SSE by removing values > edit_level
! also remove values with sse >= sse0 since these functions are not good either

do  i_GA_individual=1,n_GA_individuals  ! calculate the total populations SSE

    if( individual_quality(i_GA_individual) > 0 ) then

        if( individual_SSE(i_GA_individual) > edit_level  ) then
        !if( individual_SSE(i_GA_individual) >  edit_level .or.  &
        !    individual_SSE(i_GA_individual) >= sse0      ) then

            individual_quality( i_GA_individual ) = -1

        endif !   individual_SSE(i_GA_individual) >  edit_level

    endif !  individual_quality > 0

enddo ! i_GA_individual

!-----------------------------------------------------------------------------------

! calculate the integrated ranked fitness levels
! to support the "Fitness Proportionate Reproduction" events

write(GA_print_unit,'(/A/)')'calcfit: calculate the integrated ranked fitness levels'

!----------------------------------------------------------------------------------


! calculate a normalized ranking of the errors
!     (higher individual SSE == lower value/ranking)


! calculate the individual fitness

write(GA_print_unit,'(A,1x,E24.16)')  'calcfit: sse0 ', sse0


!write(GA_print_unit,'(A)')  'calcfit: i_GA_ind, ind_SSE, ind_ranked_fitness   ind_quality'

do  i_GA_individual=1,n_GA_individuals


    if( individual_quality( i_GA_individual ) > 0 ) then

        !write(GA_print_unit,'(A,1x,I6,2(1x,E24.16))') &
        !      'calcfit: i_GA_individual, SSE , sse0/sse', &
        !                i_GA_individual, individual_SSE(i_GA_individual), &
        !                            sse0/individual_SSE(i_GA_individual)

        ! indiv_fitness is a function

        individual_ranked_fitness(i_GA_individual) = indiv_fitness( i_GA_individual )


    else

        individual_ranked_fitness(i_GA_individual) = 0.0d0

    endif ! individual_quality( i_GA_individual ) > 0


    !write(GA_print_unit,'(I6,2(1x,E24.16),1x,I6)') &
    !      i_GA_individual, individual_SSE(i_GA_individual), &
    !            individual_ranked_fitness(i_GA_individual), &
    !                  individual_quality( i_GA_individual )

enddo ! i_GA_individual


!-------------------------------------------------------------------------

! calculate the total populations sum of SSE, mean and minimum of SSE
! calculate the sum, mean and minimum of fitness for all individuals

dble_cff=0.0D+0
n_counted = 0
min_sse = 1.0D20
index_min_sse = 0
sum_individual_SSE = 0.0D0

do i_GA_individual=1,n_GA_individuals

  !write(GA_print_unit,'(A,2(1x,I6))') &
  !      'calcfit: i_GA_individual, individual_quality( i_GA_individual )', &
  !                i_GA_individual, individual_quality( i_GA_individual )

  if( individual_quality( i_GA_individual ) > 0 ) then


      ! indiv_fitness is a function

      !dble_cff = dble_cff +  indiv_fitness( i_GA_individual )    ! new
      dble_cff = dble_cff +  individual_ranked_fitness(i_GA_individual)

      integrated_SSE(i_GA_individual)=dble_cff   ! integrated_SSE never used

      sum_individual_SSE = sum_individual_SSE +  individual_SSE(i_GA_individual)

      n_counted = n_counted + 1

      !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
      !      'calcfit: i_GA_individual, integrated_SSE(i_GA_individual)  ', &
      !                i_GA_individual, integrated_SSE(i_GA_individual)


      if( individual_SSE(i_GA_individual) < min_sse )then

          min_sse = individual_SSE(i_GA_individual)
          index_min_sse = i_GA_individual

      endif !    individual_SSE(i_GA_individual) < min_sse


  endif !   individual_quality( i_GA_individual ) > 0

enddo ! i_GA_individual


! dble_cff is now the sum of the "good" individual_SSE's

sum_individual_fit = dble_cff


mean_individual_fit = 0.0d0
if( n_counted > 0 )then
    mean_individual_fit = sum_individual_fit / real( n_counted, kind=8)
endif ! n_counted > 0


write(GA_print_unit,'(/A,1x,I6,1x,E20.10)')&
      'calcfit: generation,   sum_individual_fit  = ', i_GA_generation, sum_individual_fit
write(GA_print_unit,'(A,1x,I6,1x,E20.10)')&
      'calcfit: generation,   mean_individual_fit = ', i_GA_generation, mean_individual_fit
write(GA_print_unit,'(A,1x,I6,1x,E20.10, 1x, I6)')&
      'calcfit: generation,   min_sse, index      = ', i_GA_generation, min_sse, index_min_sse
write(GA_print_unit,'(A,2(1x,I6)/)')&
      'calcfit: generation, number good           = ', i_GA_generation, n_counted

mean_individual_SSE = 0.0D0

if( n_counted > 0 )then
    mean_individual_SSE = sum_individual_SSE / real( n_counted, kind=8)
endif ! n_counted > 0


write(GA_print_unit,'(/A,1x,I6,1x,E20.10)')&
      'calcfit: generation,   sum_individual_SSE  = ', i_GA_generation, sum_individual_SSE
write(GA_print_unit,'(A,1x,I6,1x,E20.10)')&
      'calcfit: generation,   mean_individual_SSE = ', i_GA_generation, mean_individual_SSE


!---------------------------------------------------------------------------------

!write(GA_print_unit,'(A)')  'calcfit: i_GA_ind, ind_SSE, ind_ranked_fitness   ind_quality'

mean_fitness = 0.0d0
sigma_fitness = 0.0d0
var_fitness = 0.0d0
icount = 0

do  i_GA_individual=1,n_GA_individuals

    if( individual_quality( i_GA_individual ) > 0 ) then
    !if( individual_quality( i_GA_individual ) > 0 .and.  &
    !    individual_ranked_fitness(i_GA_individual) > 1.0d0  ) then

        mean_fitness = mean_fitness +  individual_ranked_fitness(i_GA_individual)
        var_fitness  = var_fitness  +  individual_ranked_fitness(i_GA_individual)**2
        icount = icount + 1

    endif !  individual_quality( i_GA_individual ) > 0

enddo ! i_GA_individual

xn =  real( icount, kind=8 )

if( icount > 0  ) then
    mean_fitness = mean_fitness / xn
    sigma_fitness = sqrt( abs( var_fitness / xn  - mean_fitness**2 )   )
else

    mean_fitness  =  0.0d0
    sigma_fitness =  0.0d0

endif

!var_fitness   = sqrt( var_fitness  / xn )

write(GA_print_unit,'(/A,2(1x,I6))') 'calcfit: n_GA_individuals, icount       ',&
                                               n_GA_individuals, icount

write(GA_print_unit,'(/A,1x,I6,1x,E15.7)') 'calcfit: i_GA_generation, mean_fitness  ',&
                                                     i_GA_generation, mean_fitness
write(GA_print_unit,'(A,1x,I6,1x,E15.7/)') 'calcfit: i_GA_generation, sigma_fitness ',&
                                                     i_GA_generation, sigma_fitness

!---------------------------------------------------------------------------------


!write(GA_print_unit,'(/A)') 'calcfit: i_GA_individual, integrated_ranked_fitness'

dble_cff=0.0D+0
do  i_GA_individual=1,n_GA_individuals  ! calculate the sum of the rankings


    dble_cff = dble_cff + individual_ranked_fitness(i_GA_individual)

    integrated_ranked_fitness(i_GA_individual)=dble_cff

    !write(GA_print_unit,'(I6,1x,E24.16)') &
    !      i_GA_individual, integrated_ranked_fitness(i_GA_individual)


enddo ! i_GA_individual

!---------------------------------------------------------------------------------

! normalize to the integrated ranking values so that
! the ranking integration ranges from [0. to 1.]


!write(GA_print_unit,'(/A)') 'calcfit: i_GA_individual, normalized integrated_ranked_fitness'

do  i_GA_individual=1,n_GA_individuals

    if( abs( integrated_ranked_fitness(n_GA_individuals) ) > 1.0D-20 )then

        integrated_ranked_fitness(i_GA_individual) = &
        integrated_ranked_fitness(i_GA_individual) / &
                          integrated_ranked_fitness(n_GA_individuals)

    else
        integrated_ranked_fitness(i_GA_individual) = 0.0D0

    endif ! abs( integrated_ranked_fitness(n_GA_individuals) ) > 1.0D-20

!    write(GA_print_unit,'(I6,1x,E24.16)') &
!          i_GA_individual, integrated_ranked_fitness(i_GA_individual)

enddo ! i_GA_individual

!-------------------------------------------------------------------------------

if( i_GA_generation == 1                                 .or. &
    mod(i_GA_generation, GA_child_print_interval ) == 0  .or. &
    i_GA_generation == n_GA_generations       )then

    write(GA_print_unit,'(/A)')&
     'i_GA_ind   ind_SSE            ind_ranked_fitness    integ_rank_fitness  ind_quality'

    do  i_GA_individual=1,n_GA_individuals

        write(GA_print_unit,'(I6,3(1x,E20.12),1x,I6)') &
              i_GA_individual, individual_SSE(i_GA_individual), &
                    individual_ranked_fitness(i_GA_individual), &
                    integrated_ranked_fitness(i_GA_individual), &
                          individual_quality( i_GA_individual )

    enddo ! i_GA_individual

endif !  i_GA_generation == 1 ...

!-------------------------------------------------------------------------------


! select best parent
! best fit individual has largest individual_ranked_fitness


i_GA_Best_Parent=1

dble_cff=individual_ranked_fitness(1)

!write(GA_print_unit,'(/A, 2(1x,E24.16)/)') &
!      'calcfit: individual_ranked_fitness(1) , dble_cff ', &
!                individual_ranked_fitness(1) , dble_cff


do  i_GA_individual=2,n_GA_individuals

    if( individual_ranked_fitness(i_GA_individual) .gt. dble_cff) then

        dble_cff=individual_ranked_fitness(i_GA_individual)

        i_GA_Best_Parent=i_GA_Individual

    endif !   individual_ranked_fitness(i_GA_individual) .gt. dble_cff

    !write(GA_print_unit,'(A,1x,I6,2(1x,E24.16))') &
    !      'calcfit: i_GA_ind, individual_ranked_fitness(i_GA_ind), dble_cff ', &
    !         i_GA_individual, individual_ranked_fitness(i_GA_individual), dble_cff

enddo ! i_GA_individual

!------------------------------------------------------------------------------

write(GA_print_unit,'(//A,2(1x, I6))') 'calcfit: Generation, i_GA_Best_Parent ', &
                                            i_GA_Generation, i_GA_Best_Parent
write(GA_print_unit,'(A,1x,I6,1x,E24.16)') &
      'calcfit: i_GA_Best_Parent, individual_ranked_fitness(i_GA_Best_Parent) ', &
                i_GA_Best_Parent, individual_ranked_fitness(i_GA_Best_Parent)
write(GA_print_unit,'(A,1x,I6,12(1x,E15.7)/(12(1x,E15.7)))') &
      'calcfit: i_GA_Best_Parent, parent_parameters ', &
                i_GA_Best_Parent, &
                (parent_parameters( i_GA_Best_Parent, jj),jj=1,n_parameters)
                !i_GA_Best_Parent, parent_parameters( i_GA_Best_Parent, 1:n_parameters)


!-----------------------------------------------------------------------

! if fitness >= 100, then sse <= 0.01 * sse0
! if this is true, stop the run

!if( individual_ranked_fitness(i_GA_Best_Parent) >= 100.0d0 ) then
!if( individual_ranked_fitness(i_GA_Best_Parent) >= 200.0d0 ) then
!    L_stop_run = .TRUE.
!endif ! individual_ranked_fitness(i_GA_Best_Parent) >= 100.0d0

!-----------------------------------------------------------------------

return


end subroutine calc_fitness
