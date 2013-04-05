subroutine GP_Elitists
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=8) :: Fitness(n_GP_Individuals)
real(kind=8) :: tmp(n_GP_Individuals)
real(kind=8) :: cff

integer(kind=4) :: Ranked_Fitness(n_GP_Individuals)
integer(kind=4) :: icff
!integer (kind=4) :: i_GP_Individual
!integer (kind=4) :: j_GP_Individual

logical Carry_On

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Run_GP_Calculate_Fitness=.false.  ! just as a backcheck

Fitness = GP_Individual_Ranked_Fitness

do  i_GP_Individual=1,n_GP_Individuals
    Ranked_Fitness(i_GP_Individual)=i_GP_Individual
enddo

! calculate the fitness rankings

! replace with sort routine later

do  i_GP_Individual=1,n_GP_Individuals

    do  j_GP_Individual=1,n_GP_Individuals-1

        if( Fitness(j_GP_Individual+1) .gt. Fitness(j_GP_Individual)) then

            
            cff=Fitness(j_GP_Individual+1)
            Fitness(j_GP_Individual+1)=Fitness(j_GP_Individual)
            Fitness(j_GP_Individual)=cff
            icff=Ranked_Fitness(j_GP_Individual+1)
            Ranked_Fitness(j_GP_Individual+1)=Ranked_Fitness(j_GP_Individual)
            Ranked_Fitness(j_GP_Individual)=icff

        endif !   Fitness(j_GP_Individual+1) .gt....

    enddo ! j_GP_Individual

enddo ! i_GP_Individual


! rerank ALL of the Individuals to keep the code simple and not replicate copies of children

do  i_GP_Individual=1,n_GP_Individuals
    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
    GP_Adult_Population_Node_Type(Ranked_Fitness(i_GP_Individual),1:n_Nodes,1:n_Trees)
    GP_Child_Individual_SSE(i_GP_Individual)=GP_Adult_Individual_SSE(Ranked_Fitness(i_GP_Individual))
    tmp(i_GP_Individual)=GP_Individual_Ranked_Fitness(Ranked_Fitness(i_GP_Individual))
enddo !  i_GP_Individual

do  i_GP_Individual=1,n_GP_Individuals
    GP_Individual_Ranked_Fitness(i_GP_Individual)=tmp(i_GP_Individual)
enddo ! i_GP_Individual


cff=0.0
do  i_GP_Individual=1,n_GP_Individuals  ! calculate the sum of the rankings
    cff=cff+GP_Individual_Ranked_Fitness(i_GP_individual)
    GP_Integrated_Ranked_Fitness(i_GP_Individual)=cff
enddo ! i_GP_Individual


! normalize to the integrated ranking values so that the ranking integration ranges from [0. to 1.]

do  i_GP_Individual=1,n_GP_Individuals 
    GP_Integrated_Ranked_Fitness(i_GP_Individual) = &
    GP_Integrated_Ranked_Fitness(i_GP_Individual)/GP_Integrated_Ranked_Fitness(n_GP_Individuals)
enddo ! i_GP_Individual

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Elitists
!234567890123456789012345678901234567890123456789012345678901234567890
