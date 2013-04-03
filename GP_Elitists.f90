subroutine GP_Elitists
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters
use GA_Parameters
use GP_Variables
use GA_Variables

implicit none

real (kind=4) :: Fitness(n_GP_Individuals)
real (kind=4) :: tmp(n_GP_Individuals)
real (kind=4) :: cff

integer (kind=4) :: Ranked_Fitness(n_GP_Individuals)
integer (kind=4) :: icff
integer (kind=4) :: i_GP_Individual
integer (kind=4) :: j_GP_Individual

logical Carry_On

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Run_GP_Calculate_Fitness=.false.  ! just as a backcheck

Fitness=GP_Individual_Ranked_Fitness

do i_GP_Individual=1,n_GP_Individuals
  Ranked_Fitness(i_GP_Individual)=i_GP_Individual
enddo

! calculate the fitness rankings
do i_GP_Individual=1,n_GP_Individuals
  do j_GP_Individual=1,n_GP_Individuals-1
    if (Fitness(j_GP_Individual+1) .gt. Fitness(j_GP_Individual)) then
      cff=Fitness(j_GP_Individual+1)
      Fitness(j_GP_Individual+1)=Fitness(j_GP_Individual)
      Fitness(j_GP_Individual)=cff
      icff=Ranked_Fitness(j_GP_Individual+1)
      Ranked_Fitness(j_GP_Individual+1)=Ranked_Fitness(j_GP_Individual)
      Ranked_Fitness(j_GP_Individual)=icff
    endif
  enddo
enddo

! rerank ALL of the Individuals to keep the code simple and not replicate copies of children
do i_GP_Individual=1,n_GP_Individuals
  GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
  GP_Adult_Population_Node_Type(Ranked_Fitness(i_GP_Individual),1:n_Nodes,1:n_Trees)
  GP_Child_Individual_SSE(i_GP_Individual)=GP_Adult_Individual_SSE(Ranked_Fitness(i_GP_Individual))
  tmp(i_GP_Individual)=GP_Individual_Ranked_Fitness(Ranked_Fitness(i_GP_Individual))
enddo
do i_GP_Individual=1,n_GP_Individuals
  GP_Individual_Ranked_Fitness(i_GP_Individual)=tmp(i_GP_Individual)
enddo

cff=0.0
do i_GP_Individual=1,n_GP_Individuals  ! calculate the sum of the rankings
  cff=cff+GP_Individual_Ranked_Fitness(i_GP_individual)
  GP_Integrated_Ranked_Fitness(i_GP_Individual)=cff
enddo
do i_GP_Individual=1,n_GP_Individuals ! normalize to the integrated ranking values so that the ranking integration ranges from [0. to 1.]
  GP_Integrated_Ranked_Fitness(i_GP_Individual) = &
  GP_Integrated_Ranked_Fitness(i_GP_Individual)/GP_Integrated_Ranked_Fitness(n_GP_Individuals)
enddo

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Elitists
!234567890123456789012345678901234567890123456789012345678901234567890
