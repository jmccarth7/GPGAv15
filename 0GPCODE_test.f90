!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
program GPCODE_test
! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of 
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters
use GA_Parameters
use GP_Variables
use GA_Variables
use GP_Data

implicit none

real (kind=4) :: GP_Node_Parameters_Answer(n_Nodes,n_Trees)

real (kind=4) :: icount
real (kind=4) :: cff

integer (kind=4) :: GP_Node_Type_Answer(n_Nodes,n_Trees)
integer (kind=4) :: GP_Node_Type_for_Plotting(9,n_Nodes,n_Trees)

integer (kind=4) :: i_GP_Best_Parent
integer (kind=4) :: i_GP_Generation
integer (kind=4) :: i_GP_Individual
integer (kind=4) :: j_GP_Individual
integer (kind=4) :: i_Tree
integer (kind=4) :: i_Node

logical bad_node

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

CALL RANDOM_SEED(size = n_seed)
ALLOCATE(seed(n_seed))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)
CALL RANDOM_SEED(PUT = seed)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! set the twin experiment 'nature' 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]
GP_Individual_Node_Parameters=0.                                          ! Matrix Operation
GP_Individual_Node_Type=-9999                                             ! Matrix Operation

include 'Lotka_Volterra_Example_Set_Up.f901'
!off include 'Franks_etal_NPZ_Mar_Bio_Example_Set_Up.f901'

! set the desired 'twin experiment' population node type and parameter using the info from the set up file
GP_Node_Type_Answer=GP_Individual_Node_Type                               ! Matrix Operation
GP_Node_Parameters_Answer=GP_Individual_Node_Parameters                   ! Matrix Operation

do i_GP_Individual=1,9
  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees)=GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)
enddo
write(50) GP_Node_Type_for_Plotting

! set the Initial Conditions, Model Parameters and Node Type for the 'twin experiment case'
Runga_Kutta_Solution(0,1:n_CODE_Equations)=Runga_Kutta_Initial_Conditions(1:n_CODE_Equations) ! Matrix Operation
Runga_Kutta_Node_Parameters=dble(GP_Individual_Node_Parameters)           ! Matrix Operation
Runga_Kutta_Node_Type=GP_Individual_Node_Type                             ! Matrix Operation

!off for testing call Runga_Kutta_Box_Model

Data_Array=Runga_Kutta_Solution                                           ! Matrix Operation

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! blank/set the values [0. = zero-valued parameter; -9999 = blank node type]
GP_Individual_Node_Parameters=0.                                          ! Matrix Operation
GP_Individual_Node_Type=-9999                                             ! Matrix Operation
GP_Population_Node_Parameters=0.                                          ! Matrix Operation
GP_Adult_Population_Node_Type=-9999                                       ! Matrix Operation
GP_Child_Population_Node_Type=-9999                                       ! Matrix Operation

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! run the Genetic Programming optimization routine for the Binary Tree Evolution
! with the embedded GA_lmdif parameter optimization scheme
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

n_GP_Elitists=int(GP_Elitist_Probability*n_GP_individuals)                           ! Number of Carry-Over Elitists
n_GP_Asexual_Reproductions=int(GP_Asexual_Reproduction_Probability*n_GP_individuals) ! Number of GP Fitness Proportionate Reproduction
n_GP_Crossovers=int(GP_Crossover_Probability*n_GP_individuals)                              ! Number of GP Sexual Crossovers
n_GP_Mutations=n_GP_Individuals-(n_GP_Elitists+n_GP_Crossovers+n_GP_Asexual_Reproductions)  ! Number of GP Mutations

if (n_GP_Elitists+n_GP_Asexual_Reproductions+n_GP_Crossovers+n_GP_Mutations .gt. n_GP_Individuals) then
  write(*,*) 'Sum of n_GP_Elitists+n_Asexual_Reproduction+n_GP_Crossovers+n_GP_Mutations is too high'
  stop
elseif (n_GP_Elitists+n_GP_Asexual_Reproductions+n_GP_Crossovers+n_GP_Mutations .lt. n_GP_Individuals) then
  write(*,*) 'Sum of n_GP_Elitists+n_Asexual_Reproduction+n_GP_Crossovers+n_GP_Mutations is too low'
  stop
endif

do i_GP_Generation=1,n_GP_Generations
  write(*,*) 'GP Generation # ',i_GP_Generation,' is underway',n_Nodes*n_Trees

  GP_Child_Population_Node_Type=GP_Adult_Population_Node_Type                               ! Matrix Operation

  Run_GP_Calculate_Fitness=.false.  ! determines if the new GP child has to be sent to GA_lmdif for parameter optimization

  if (i_GP_Generation .eq. 1) then  ! randomly create the initial tree arrays for each individual and send them all to GA_lmdif for parameter optimization

    Run_GP_Calculate_Fitness=.true. ! determines if the new GP child has to be sent to GA_lmdif for parameter optimization
    call GP_Tree_Build              ! initialize the GP_Adult_Population_Node_Type array with random trees

  else ! create the next 'generation' of tree structures using either: 
!         i)  Copy the top n_GP_Elitists individuals into the next generation
!        ii)  GP Fitness-Proportionate Asexual Reproduction;
!       iii)  GP Tournament-Style Sexual Reproduction, and;
!        iv)  GP Mutation

!   i) move the top n_GP_Elitists into the next generation
    if (n_GP_Elitists .gt. 0) call GP_Elitists
    
!   ii) Carry out "GP Fitness-Proportionate Reproduction"
    if (n_GP_Asexual_Reproductions .gt. 0) call GP_Fitness_Proportionate_Asexual_Reproduction

!   iii) Carry out "GP Tree Crossover" Operations Using Tournament-Style Sexual Reproduction Selection and randomly use it to replace the new children
    if (n_GP_Crossovers .gt. 0) call GP_Tournament_Style_Sexual_Reproduction

!   iv) Carry out "GP Parameter Mutation" Operations
    if (n_GP_Mutations .gt. 0)  call GP_Mutations

!   move the children into adulthood
    GP_Adult_Population_Node_Type=GP_Child_Population_Node_Type      ! Matrix Operation

  endif

! sweep through all the GP_Adult_Population_Node_Type to replace function nodes that have both terminals set as parameters with a parameter itself
  call GP_Clean_Tree_Nodes

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! GA_lmdif subroutine segment
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  do i_GP_Individual=1,n_GP_Individuals

    if (Run_GP_Calculate_Fitness(i_GP_Individual)) then

      GP_Individual_Node_Type(1:n_Nodes,1:n_Trees)=GP_Adult_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

!     calculate how many parameters total to fit for the specific individual CODE
      n_Parameters=0
      do i_Tree=1,n_Trees
        do i_Node=1,n_Nodes
           if (GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0) then  ! there is a set parameter
             n_Parameters=n_Parameters+1
          endif
        enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND LINK THE SSE OUTPUT TO THE ARRAY AT THE END
! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEEDS TO BE PULLED OUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cff=0.
      icount=0.
      do i_Node=1,n_Nodes
        do i_Tree=1,n_Trees

            if (GP_Individual_Node_Type(i_Node,i_Tree) .ne. GP_Node_Type_Answer(i_Node,i_Tree)) then
              if (GP_Node_Type_Answer(i_Node,i_Tree) .eq. -9999 .or. GP_Individual_Node_Type(i_Node,i_Tree) .eq. -9999.) then
                cff=1.
                icount=icount+cff
              else
                cff=float(GP_Individual_Node_Type(i_Node,i_Tree)-GP_Node_Type_Answer(i_Node,i_Tree))**2
                icount=icount+cff
              endif
            endif

        enddo
      enddo
      GP_Child_Individual_SSE(i_GP_Individual)=icount

!off      GP_Child_Individual_SSE(i_GP_individual)=GA_Individual_Lowest_SSE  ! from the GA_lmdif routine output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     set the GA_lmdif-optimized initial condition array
      GP_Population_Initial_Conditions(i_GP_Individual,1:n_CODE_Equations)=GP_Individual_Initial_Conditions(1:n_CODE_Equations) ! Matrix Operation

!     set the GA_lmdif-optimized CODE parameter set array
      GP_Population_Node_Parameters(i_GP_Individual,1:n_Nodes,1:n_Trees)=GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) ! Martix Operation

    endif

  enddo
  write(*,*) GP_Child_Individual_SSE

! fitness reset region

  cff=0.0
  do i_GP_Individual=1,n_GP_Individuals  ! calculate the total population's SSE
    cff=cff+GP_Child_Individual_SSE(i_GP_Individual)
  enddo
  
  do i_GP_Individual=1,n_GP_Individuals  ! calculate a normalized ranking of the errors (higher individual SSE == lower value/ranking)
    GP_Individual_Ranked_Fitness(i_GP_Individual)=(cff-GP_Child_Individual_SSE(i_GP_Individual))/cff
  enddo
  cff=0.0
  do i_GP_Individual=1,n_GP_Individuals  ! calculate the sum of the rankings
    cff=cff+GP_Individual_Ranked_Fitness(i_GP_Individual)
    GP_Integrated_Ranked_Fitness(i_GP_Individual)=cff
  enddo
  do i_GP_Individual=1,n_GP_Individuals ! normalize to the integrated ranking values so that the ranking integration ranges from [0. to 1.]
    GP_Integrated_Ranked_Fitness(i_GP_Individual) = &
    GP_Integrated_Ranked_Fitness(i_GP_Individual)/GP_Integrated_Ranked_Fitness(n_GP_Individuals)
  enddo

  write(*,*) i_GP_Generation,'MAIN',GP_Individual_Ranked_Fitness(1)

  i_GP_Best_Parent=1
  cff=GP_Individual_Ranked_Fitness(1)
  do i_GP_Individual=2,n_GP_individuals
    if (GP_Individual_Ranked_Fitness(i_GP_individual) .gt. cff) then
      cff=GP_Individual_Ranked_Fitness(i_GP_individual)
      i_GP_Best_Parent=i_GP_Individual
    endif
  enddo

  GP_Adult_Individual_SSE=GP_Child_Individual_SSE

! write(*,*) 'SSE = ',GP_Adult_Individual_SSE
! write(*,*)
! write(*,*) 'Fitness out= ',GP_Individual_Ranked_Fitness
! stop

!off if (i_GP_Generation .eq. 3) Stop

do i_GP_Individual=1,9
  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
  GP_Adult_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

!off  GP_Node_Type_for_Plotting(i_GP_Individual,1:n_Nodes,1:n_Trees)=GP_Node_Type_Answer(1:n_Nodes,1:n_Trees)

enddo

write(50) GP_Node_Type_for_Plotting

enddo

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
stop
end program GPCODE_test
!234567890123456789012345678901234567890123456789012345678901234567890
