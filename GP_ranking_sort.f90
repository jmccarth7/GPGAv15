subroutine GP_ranking_sort


use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module

implicit none

real(kind=8) :: cff

integer(kind=4),dimension(n_GP_Individuals)  :: Ranked_Fitness_Index

integer(kind=4) :: i_GP_Individual
integer(kind=4) :: j_GP_Individual

integer(kind=4) :: icff              


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!write(6,'(/A/)') 'gprs: entry GP_ranking_sort '

! Set up a simple 'index' array

do  i_GP_Individual=1,n_GP_Individuals
    Ranked_Fitness_Index(i_GP_Individual)=i_GP_Individual
enddo

!write(6,'(/A)') 'gprs: before sort '

!write(6,'(A)')                    &
!      'gprs:i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
!                            &GP_Child_Population_SSE(i_GP_Individual)'

!do  i_GP_Individual=1,n_GP_Individuals
!    write(6,'(I6,1x, I6, 1x, E15.7)') &
!          i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
!                           GP_Child_Population_SSE(i_GP_Individual)
!enddo



! Now, rank the Individual SSE so that 
! the Individual with the lowest (highest) SSE is First (Last)

do  i_GP_Individual=1,n_GP_Individuals

    do  j_GP_Individual=1,n_GP_Individuals-1

        if( GP_Child_Population_SSE(j_GP_Individual+1) .lt. &
              GP_Child_Population_SSE(j_GP_Individual)) then

            !     Swap the two ranked fitness and index array values around
            cff=GP_Child_Population_SSE(j_GP_Individual)
            GP_Child_Population_SSE(j_GP_Individual)=GP_Child_Population_SSE(j_GP_Individual+1)
            GP_Child_Population_SSE(j_GP_Individual+1)=cff

            icff=Ranked_Fitness_Index(j_GP_Individual)
            Ranked_Fitness_Index(j_GP_Individual)=Ranked_Fitness_Index(j_GP_Individual+1)
            Ranked_Fitness_Index(j_GP_Individual+1)=icff

        endif !GP_Child_Population_SSE(j_GP_Individual+1) .lt. ...

    enddo ! j_GP_Individual

enddo  ! i_GP_Individual


!write(6,'(/A)') 'gprs: after  sort '

!write(6,'(A)')                    &
!      'gprs:i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
!                            &GP_Child_Population_SSE(i_GP_Individual)'

!do  i_GP_Individual=1,n_GP_Individuals
!    write(6,'(I6,1x, I6, 1x, E15.7)') &
!          i_GP_Individual, Ranked_Fitness_Index(i_GP_Individual), &
!                           GP_Child_Population_SSE(i_GP_Individual)
!enddo


! Re-rank ALL of the Individuals to keep the code simple and not replicate copies of children

GP_Adult_Population_Node_Type(1:n_Nodes,1:n_Trees,  1:n_GP_Individuals) = &
GP_Child_Population_Node_Type(1:n_Nodes,1:n_Trees,  Ranked_Fitness_Index(1:n_GP_Individuals) )

GP_Adult_Population_Parameter_Solution(1:n_Maximum_Number_Parameters,  1:n_GP_Individuals) = &
GP_Child_Population_Parameter_Solution(1:n_Maximum_Number_Parameters,  &
                                                      Ranked_Fitness_Index(1:n_GP_Individuals) )

GP_Adult_Population_SSE=GP_Child_Population_SSE

! Copy this back across to the Child Population values
! to allow the Elite codes to propagate along in the next generations

GP_Child_Population_Node_Type=GP_Adult_Population_Node_Type
GP_Child_Population_Parameter_Solution=GP_Adult_Population_Parameter_Solution


!------------------------------------------------------------------------------------------

! Calculate the Adult Population's Total SSE

cff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals
    cff=cff+GP_Child_Population_SSE(i_GP_Individual)
enddo

! Calculate a simple 'normalized' ranking of the SSE as an estimate of fitness

! [Fitness = (Total-SSE)/Total ==> higher individual SSE == lower value/ranking; Ranging from 0-1]

do  i_GP_Individual=1,n_GP_Individuals
    GP_Population_Ranked_Fitness(i_GP_Individual) = (cff-GP_Child_Population_SSE(i_GP_Individual))/cff
enddo  ! i_GP_Individual

! Calculate the Integrated Ranked Fitness values for creating the next generation

cff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals
    cff=cff+GP_Population_Ranked_Fitness(i_GP_individual)
    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual)=cff
enddo ! i_GP_Individual

! Normalize to the integrated ranking values so that the ranking integration ranges from [0. to 1.]

do  i_GP_Individual=1,n_GP_Individuals
    GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) =  &
         GP_Integrated_Population_Ranked_Fitness(i_GP_Individual) /  &
                   GP_Integrated_Population_Ranked_Fitness(n_GP_Individuals)
enddo  ! i_GP_Individual

!------------------------------------------------------------------------------------------

!write(6,'(/A)') 'gprs: after  sort '

!write(6,'(A)')                    &
!      'gprs:i_GP_Individual, GP_Integ_Pop_Ranked_Fitness, GP_Pop_Ranked_Fitness, GP_Child_Pop_SSE'

!do  i_GP_Individual=1,n_GP_Individuals
!    write(6,'(I6, 3(1x, E15.7))') &
!          i_GP_Individual, GP_Integrated_Population_Ranked_Fitness(i_GP_Individual), &
!                           GP_Population_Ranked_Fitness(i_GP_Individual), &
!                           GP_Child_Population_SSE(i_GP_Individual)
!enddo   ! i_GP_Individual



!write(6,'(/A)') 'gprs: at return   '

return

end subroutine GP_ranking_sort
