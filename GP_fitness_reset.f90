subroutine GP_fitness_reset( )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module


implicit none

!real (kind=4) :: GP_Node_Parameters_Answer(n_Nodes,n_Trees)
!integer (kind=4) :: GP_Node_Type_Answer(n_Nodes,n_Trees)
!integer (kind=4) :: GP_Node_Type_for_Plotting(9,n_Nodes,n_Trees)

!integer(kind=4) :: i_GP_Best_Parent
!integer(kind=4) :: i_GP_Generation
!integer(kind=4) :: i_GP_Individual
!integer(kind=4) :: j_GP_Individual
!integer(kind=4) :: i_Tree
!integer(kind=4) :: i_Node

logical bad_node


integer :: i
integer :: nop
integer :: message_len

logical :: Lplot


real(kind=8), allocatable, dimension(:) :: answer
real(kind=8), allocatable, dimension(:) :: output_array

real (kind=8) :: ssum, ssum2, totobs, cff

real (kind=8) :: xcount
!----------------------------------------------------------------------------------------



!write(*,*) GP_Child_Individual_SSE
!if( myid == 0 )then
!    write(6,'(A,1x,I6)') '0: i_GP_generation ',  i_GP_generation
!
!    write(6,'(/A/(5(1x,E15.7)))') '0: GP_Child_Individual_SSE = ',&
!                                      GP_Child_Individual_SSE
!
!endif ! myid == 0 
!---------------------------------------------------------------------------------------------------

! fitness reset region

! calculate the total population's SSE

cff=0.0d0
do  i_GP_Individual=1,n_GP_Individuals
    cff=cff+GP_Child_Individual_SSE(i_GP_Individual)
enddo

write(6,'(A,1x,I6)') '0: i_GP_generation ',  i_GP_generation
write(6,'(/A,1x,E15.7)') '0: sum of all GP_Child_Individual_SSE = ', cff  


! calculate a normalized ranking of the errors (higher individual SSE == lower value/ranking)
do  i_GP_Individual=1,n_GP_Individuals
    GP_Individual_Ranked_Fitness(i_GP_Individual) = &
         ( cff - GP_Child_Individual_SSE(i_GP_Individual) ) / cff
enddo

do  i_GP_Individual=1,n_GP_Individuals
    write(6,*) 'gpfr: i_GP_individual, GP_Individual_Ranked_Fitness(i_GP_Individual) ', &
                      i_GP_individual, GP_Individual_Ranked_Fitness(i_GP_Individual) 
enddo


! calculate the sum of the rankings

cff=0.0
do  i_GP_Individual=1,n_GP_Individuals
    cff=cff+GP_Individual_Ranked_Fitness(i_GP_Individual)
    GP_Integrated_Ranked_Fitness(i_GP_Individual)=cff
enddo

do  i_GP_Individual=1,n_GP_Individuals
    write(6,*) 'gpfr: i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual) ', &
                      i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual) 
enddo

! normalize to the integrated ranking values so that the ranking integration ranges from [0. to 1.]

do  i_GP_Individual=1,n_GP_Individuals
    GP_Integrated_Ranked_Fitness(i_GP_Individual) = &
    GP_Integrated_Ranked_Fitness(i_GP_Individual)/GP_Integrated_Ranked_Fitness(n_GP_Individuals)
    write(6,*) 'gpfr: norm  i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual) ', &
                            i_GP_individual, GP_Integrated_Ranked_Fitness(i_GP_Individual) 
enddo


!write(*,*) i_GP_Generation,'MAIN',GP_Individual_Ranked_Fitness(1)
write(*,'(A,1x,I6,1x,E15.7)') &
      '0: i_GP_Generation, MAIN GP_Individual_Ranked_Fitness(1) ', &
          i_GP_Generation,      GP_Individual_Ranked_Fitness(1)



!---------------------------------------------------------------------------------------------------




i_GP_Best_Parent=1
cff=GP_Individual_Ranked_Fitness(1)
do  i_GP_Individual=2,n_GP_individuals
    if( GP_Individual_Ranked_Fitness(i_GP_individual) .gt. cff) then
        cff=GP_Individual_Ranked_Fitness(i_GP_individual)
        i_GP_Best_Parent=i_GP_Individual
    endif
enddo

write(*,'(A,2(1x,I6))') '0: i_GP_Generation, i_GP_Best_Parent ', &                        
                            i_GP_Generation, i_GP_Best_Parent


GP_Adult_Individual_SSE  =  GP_Child_Individual_SSE

write(*,'(/A/(5(1x,E15.7)))') '0: GP_Adult_Individual_SSE = ',&
                                  GP_Adult_Individual_SSE

write(*,'(/A/(5(1x,E15.7)))') '0: GP_Individual_Ranked_Fitness ',&
                                  GP_Individual_Ranked_Fitness



!off if( i_GP_Generation .eq. 3) Stop

do i_GP_Individual=1,n_GP_individuals

   GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
   GP_Adult_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

   !off  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees)=GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)

enddo ! i_GP_individual



write(unit_gp_out) GP_Node_Type_for_Plotting


return


end program GPCODE_GA_lmdif_parameter_optimization_test
