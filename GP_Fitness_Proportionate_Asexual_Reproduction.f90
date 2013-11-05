subroutine GP_Fitness_Proportionate_Asexual_Reproduction

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff

integer(kind=4) :: i
integer(kind=4) :: icff
integer(kind=4) :: i_GP_individual
integer(kind=4) :: j_GP_Individual
integer(kind=4) :: i_GP_Asexual_Reproduction

integer(kind=4) :: i_parm
integer(kind=4) :: iii

real(kind=8) :: sse_ind

!-----------------------------------------------------------------------------

!do  i = 1, n_GP_individuals
!    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
!         'gpfpar: i, GP_Child_Individual_SSE(i)',&
!                  i, GP_Child_Individual_SSE(i)
!enddo



i_GP_Individual = n_GP_Elitists

write(GP_print_unit,'(/A,1x,I6)' ) &
      'gpfpar: n_GP_Asexual_Reproductions ', n_GP_Asexual_Reproductions
write(GP_print_unit,'(A,1x,I6)' ) &
      'gpfpar: n_GP_Elitists         ', n_GP_Elitists
write(GP_print_unit,'(A,1x,I6)' ) &
      'gpfpar: start i_GP_individual ', n_GP_Elitists  + 1


do  i_GP_Asexual_Reproduction=1,n_GP_Asexual_Reproductions


    !--------------------------------------------------------------------------

    i_GP_Individual=i_GP_Individual+1

    write(6,'(A,2(1x,I6))') &
          'gpfpar: i_GP_Asexual_Reproduction, i_GP_individual ', &
                   i_GP_Asexual_Reproduction, i_GP_individual

    sse_ind = GP_Child_Individual_SSE(i_GP_Individual)

    !write(6,'(A,1x,I6,1x,E16.7)') &
    !      'gpfpar: i_GP_individual, GP_Child_Individual_SSE(i_GP_Individual)  ', &
    !               i_GP_individual, GP_Child_Individual_SSE(i_GP_Individual)
    !--------------------------------------------------------------------------


    call Random_Number(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    ! GP_Integrated_Population_Ranked_Fitness is normalized so that
    ! the range is from [0. to 1.]

    icff = -1

    do  j_GP_Individual=1,n_GP_Individuals

        !write(GP_print_unit,'(A,1x,I6,2(1x,E15.7))' ) &
        !  'gpfpar: j_GP_Indiv, cff, &
        !  &GP_Integ_Pop_Ranked_Fitness(j_GP_Indiv) ', &
        !           j_GP_Individual, cff, &
        !   GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)

        if( cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)) then

            icff=j_GP_Individual

            !write(GP_print_unit,'(A,2(1x,I6),1(1x,E15.7))' ) &
            !  'gpfpar: j_GP_Indiv, icff, &
            !  &GP_Integ_Pop_Ranked_Fitness(j_GP_Indiv) ', &
            !           j_GP_Individual, icff, &
            !   GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)

            exit

        endif !   cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)

    enddo ! j_GP_Individual

    ! index to move over both the parent parameters and the individual fitness levels

    j_GP_Individual=icff

    !write(GP_print_unit,'(/A,2(1x,I6)/)' ) 'gpfpar: j_GP_Individual, icff ', &
    !                                                j_GP_Individual, icff

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
       GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,j_GP_Individual)


    if( i_GP_individual /= j_GP_individual )then
        GP_Adult_Population_Parameter_Solution(1:n_Maximum_Number_Parameters,i_GP_Individual) = &
        GP_Adult_Population_Parameter_Solution(1:n_Maximum_Number_Parameters,j_GP_Individual)
    endif !  i_GP_individual /= j_GP_individual

    GP_Population_Node_Parameters(1:n_Nodes,1:n_Trees, i_GP_Individual) = &             ! 20131030
            GP_Population_Node_Parameters(1:n_Nodes,1:n_Trees, j_GP_Individual)         ! 20131030

    GP_Population_Initial_Conditions(1:n_CODE_Equations, i_GP_Individual) = &           ! 20131030
            GP_Population_Initial_Conditions(1:n_CODE_Equations, j_GP_Individual)       ! 20131030


    ! give the child the adult's SSE value
    GP_Child_Individual_SSE(i_GP_Individual) = GP_Adult_Population_SSE(j_GP_Individual)



    !do  iii = 1, n_maximum_number_parameters
    !    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)' ) &
    !      'gpfpar: j_GP_individual, iii,GP_Adult_Population_Parameter_Solution(iii,j_GP_Individual)',&
    !               j_GP_individual, iii,GP_Adult_Population_Parameter_Solution(iii,j_GP_Individual)
    !enddo

    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
       'gpfpar: j_GP_individual, GP_Adult_Population_SSE(j_GP_Individual)',&
                j_GP_individual, GP_Adult_Population_SSE(j_GP_Individual)
    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
          'gpfpar: PREVIOUS i_GP_individual, GP_Child_Individual_SSE(i_GP_Individual)',&
                            i_GP_individual, sse_ind
    write(GP_print_unit,'(A)' ) &
          'gpfpar: j_GP_individual replaces i_GP_individual'



enddo ! i_GP_Asexual_Reproduction




! debug
!write(6,'(/A)') 'gpfpar: GP_Adult_Population_Parameter_Solution  at end of GP_Fit...'
!call print_debug_real_nparm( &
!       'GP_Adult_Population_Parameter_Solution  at end of GP_Fit...' , &
!        GP_Adult_Population_Parameter_Solution )


return

end subroutine GP_Fitness_Proportionate_Asexual_Reproduction
