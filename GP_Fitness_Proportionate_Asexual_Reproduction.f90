subroutine GP_Fitness_Proportionate_Asexual_Reproduction
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) cff
real(kind=8) dff

integer(kind=4) icff
integer(kind=4) i_Individual
integer(kind=4) j_Individual
integer(kind=4) i_GP_Asexual_Reproduction

logical Carry_On

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



i_GP_Individual = n_GP_Elitists ! + n_GP_Asexual_Reproductions


write(GP_print_unit,'(A,1x,I6)' ) 'gpfpar: n_GP_Elitists         ', n_GP_Elitists
write(GP_print_unit,'(A,1x,I6)' ) 'gpfpar: start i_GP_individual ', &
                               n_GP_Elitists  + 1


do i_GP_Asexual_Reproduction=1,n_GP_Asexual_Reproductions

  i_GP_Individual=i_GP_Individual+1
  call Random_Number(cff) ! uniform random number generator
  Carry_On=.true.

  ! normalize to the maximum values so that the range is from [0. to 1.]

  do j_GP_Individual=1,n_GP_Individuals

    if( Carry_On ) then

      if( cff .le. GP_Integrated_Ranked_Fitness(j_GP_Individual)) then

        icff=j_GP_Individual
        Carry_On=.false.

      endif !   cff .le. GP_Integrated_Ranked_Fitness(j_GP_Individual)

    endif ! Carry_On

  enddo ! j_GP_Individual

  ! index to move over both the parent parameters and the individual fitness levels
  j_GP_Individual=icff


  GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
  GP_Adult_Population_Node_Type(j_GP_Individual,1:n_Nodes,1:n_Trees)

  GP_Child_Individual_SSE(i_GP_Individual)=GP_Adult_Individual_SSE(j_GP_Individual)

  Run_GP_Calculate_Fitness(i_GP_Individual)=.false.


enddo ! i_GP_Asexual_Reproduction

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Fitness_Proportionate_Asexual_Reproduction
!234567890123456789012345678901234567890123456789012345678901234567890
