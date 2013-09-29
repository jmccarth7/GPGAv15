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


real(kind=8) :: sse_ind

!-----------------------------------------------------------------------------

!do  i = 1, n_GP_individuals
!    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
!         'gpfpar: i, GP_Child_Population_SSE(i)',&
!                  i, GP_Child_Population_SSE(i)
!enddo
  


i_GP_Individual = n_GP_Elitists

!if( myid == 0 )then
!    write(GP_print_unit,'(/A,1x,I6)' ) &
!          'gpfpar: n_GP_Asexual_Reproductions ', n_GP_Asexual_Reproductions
!    write(GP_print_unit,'(A,1x,I6)' ) &
!          'gpfpar: n_GP_Elitists         ', n_GP_Elitists
!    write(GP_print_unit,'(A,1x,I6)' ) &
!          'gpfpar: start i_GP_individual ', n_GP_Elitists  + 1
!endif ! myid == 0


do  i_GP_Asexual_Reproduction=1,n_GP_Asexual_Reproductions


    !--------------------------------------------------------------------------

    i_GP_Individual=i_GP_Individual+1

    !call Random_Number(cff) ! uniform random number generator

    !i_GP_individual = min( nint( cff * n_GP_Individuals ) , n_GP_Individuals ) 
    !i_GP_individual = max( 1, i_GP_individual )                                

    !if( myid == 0 )then
    !    write(6,'(A,2(1x,I6))') &
    !          'gpfpar: i_GP_Asexual_Reproduction, i_GP_individual ', &
    !                   i_GP_Asexual_Reproduction, i_GP_individual 
    !endif ! myid == 0

    sse_ind = GP_Child_Population_SSE(i_GP_Individual)

    !--------------------------------------------------------------------------


    call Random_Number(cff) ! uniform random number generator

    ! the range of cff is [0. to 1.]

    ! GP_Integrated_Population_Ranked_Fitness is normalized so that
    ! the range is from [0. to 1.]

    icff = -1

    do  j_GP_Individual=1,n_GP_Individuals

        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6,2(1x,E15.7))' ) &
        !      'gpfpar: j_GP_Indiv, cff, &
        !      &GP_Integ_Pop_Ranked_Fitness(j_GP_Indiv) ', &
        !               j_GP_Individual, cff, &
        !       GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)
        !endif ! myid == 0

        if( cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)) then

            icff=j_GP_Individual

            !if( myid == 0 )then
            !    write(GP_print_unit,'(A,2(1x,I6),1(1x,E15.7))' ) &
            !      'gpfpar: j_GP_Indiv, icff, &
            !      &GP_Integ_Pop_Ranked_Fitness(j_GP_Indiv) ', &
            !               j_GP_Individual, icff, &
            !       GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)
            !endif ! myid == 0

            exit

        endif !   cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)

    enddo ! j_GP_Individual

    ! index to move over both the parent parameters and the individual fitness levels

    j_GP_Individual=icff

    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A,2(1x,I6)/)' ) 'gpfpar: j_GP_Individual, icff ', &
    !                                                    j_GP_Individual, icff
    !endif ! myid == 0

    GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,i_GP_Individual) = &
       GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,j_GP_Individual)

    GP_Child_Population_Parameter_Solution(1:n_Maximum_Number_Parameters,i_GP_Individual) = &
       GP_Adult_Population_Parameter_Solution(1:n_Maximum_Number_Parameters,j_GP_Individual)

    ! give the child the adult's SSE value
    GP_Child_Population_SSE(i_GP_Individual) = GP_Adult_Population_SSE(j_GP_Individual)


    !if( myid == 0 )then
    !    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
    !          'gpfpar: j_GP_individual, GP_Adult_Population_SSE(j_GP_Individual)',&
    !                   j_GP_individual, GP_Adult_Population_SSE(j_GP_Individual)
    !    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
    !          'gpfpar: i_GP_individual, GP_Child_Population_SSE(i_GP_Individual)',&
    !                   i_GP_individual, sse_ind
    !    write(GP_print_unit,'(A)' ) &
    !          'gpfpar: j_GP_individual replaces i_GP_individual'
    !endif ! myid == 0



    !!!Run_GP_Calculate_Fitness(i_GP_Individual)=.false.


enddo ! i_GP_Asexual_Reproduction


return

end subroutine GP_Fitness_Proportionate_Asexual_Reproduction
