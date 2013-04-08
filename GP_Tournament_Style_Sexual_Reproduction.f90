subroutine GP_Tournament_Style_Sexual_Reproduction

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! randomly choose two 'parents' using the Tournament-Style Selection 
! and cross the parameter strings to create two new 'children' parameter strings
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real (kind=4) :: cff
real (kind=8) :: dff

integer(kind=4) :: i_GP_Crossover
integer(kind=4) :: i_GP_Crossover_Point
integer(kind=4),dimension(2) :: k_GP_Individual_Male
integer(kind=4),dimension(2) :: k_GP_Individual_Female

!integer(kind=4),dimension(n_GP_Parameters) :: Child_One_Tree
!integer(kind=4),dimension(n_GP_Parameters) :: Child_Two_Tree

integer(kind=4) :: i_Male_Tree
integer(kind=4) :: i_Female_Tree
!integer(kind=4) :: i_GP_Individual
integer(kind=4) :: i_Node_Count
!integer(kind=4) :: i_Tree
integer(kind=4) :: icff

logical CROSS

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_GP_Individual=n_GP_Elitists+n_GP_Asexual_Reproductions


do i_GP_Crossover=1,n_GP_Crossovers


  i_GP_Individual=i_GP_Individual+1


  !----------------------------------------------------------------------

  ! pick the male parent for sexual crossing of parameter strings

  call random_number(cff) ! uniform random number generator
  k_GP_Individual_Male(1)=1+int(cff*float(n_GP_Individuals-1))

  call random_number(cff) ! uniform random number generator
  k_GP_Individual_Male(2)=1+int(cff*float(n_GP_Individuals-1))

  if( k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)) then

    if( k_GP_Individual_Male(1) .ne. n_GP_Individuals) then
        k_GP_Individual_Male(2)=k_GP_Individual_Male(1)+1
    else
        k_GP_Individual_Male(2)=k_GP_Individual_Male(1)-1
    endif !   k_GP_Individual_Male(1) .ne. n_GP_Individuals

  endif !   k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)

! select the individual of the two with the best fitness

  if( GP_Individual_Ranked_Fitness(k_GP_Individual_Male(2)) .gt.          &
      GP_Individual_Ranked_Fitness(k_GP_Individual_Male(1))      ) then

      k_GP_Individual_Male(1)=k_GP_Individual_Male(2)

  endif !   GP_Individual_Ranked_Fitness(k_GP_Individual_Male(2)) .gt....


  !----------------------------------------------------------------------

! pick the female parent for sexual crossing of parent parameter strings

  call Random_Number(cff) ! uniform random number generator
  k_GP_Individual_Female(1)=1+int(cff*float(n_GP_Individuals-1))

  call Random_Number(cff) ! uniform random number generator
  k_GP_Individual_Female(2)=1+int(cff*float(n_GP_Individuals-1))

  if( k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)  ) then

      if( k_GP_Individual_Female(1) .ne. N_GP_Individuals) then
          k_GP_Individual_Female(2)=k_GP_Individual_Female(1)+1
      else
          k_GP_Individual_Female(2)=k_GP_Individual_Female(1)-1
      endif !   k_GP_Individual_Female(1) .ne. N_GP_Individuals

  endif !   k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)

  ! select the individual of the two with the best fitness

  if( GP_Individual_Ranked_Fitness(k_GP_Individual_Female(2)) .gt.        &
      GP_Individual_Ranked_Fitness(k_GP_Individual_Female(1))      ) then

      k_GP_Individual_Female(1)=k_GP_Individual_Female(2)

  endif !   GP_Individual_Ranked_Fitness(k_GP_Individual_Female(2)) ...


  !----------------------------------------------------------------------

  ! randomly choose which tree structures from the male and female GP_CODEs 
  ! will participate in the genetic crossovers
  ! Find out how many trees there are in each GP_CODE

  i_Node_Count=0
  do i_Tree=1,n_Trees
    if( GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree) .ne. -9999)then

        i_Node_Count=i_Node_Count+1

    endif !   GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree)...
  enddo ! i_tree

  !write(6,'(A,1x,I6)') 'gpto: i_Node_Count ', i_Node_Count


  if( i_Node_Count .gt. 0) then

    CROSS=.true.  ! there is at least one Tree structure to cross with

    call Random_Number(cff) ! uniform random number generator
    icff=1+int(cff*float(i_Node_Count-1))  ! pick a tree

    i_Node_Count=0

    do i_Tree=1,n_Trees

      if( GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree) .ne. -9999) then

          i_Node_Count=i_Node_Count+1
          if( i_Node_Count .eq. icff) i_Male_Tree = i_Tree

      endif !   GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree)...

    enddo  ! i_tree

  else

    CROSS=.false. ! there are no Trees on this GP_CODE

  endif !   i_Node_Count .gt. 0
    
  if( CROSS) then

    i_Node_Count=0
    do  i_Tree=1,n_Trees

        if( GP_Adult_Population_Node_Type(k_GP_Individual_FeMale(1),1,i_Tree) .ne. -9999)then
            i_Node_Count=i_Node_Count+1
        endif ! GP_Adult_Population_Node_Type(k_GP_Individual_FeMale(1),1,i_Tree) .ne. -9999

    enddo ! i_tree


    if( i_Node_Count .gt. 0) then

      CROSS=.true.

      call Random_Number(cff) ! uniform random number generator
      icff=1+int(cff*float(i_Node_Count-1))  ! pick tree

      i_Node_Count=0
      do  i_Tree=1,n_Trees

          if( GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1,i_Tree) .ne. -9999) then

              i_Node_Count=i_Node_Count+1
              if( i_Node_Count .eq. icff) i_Female_Tree=i_Tree

          endif !GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1,i_Tree) ...

      enddo ! i_tree

    else

      CROSS=.false.

    endif !   i_Node_Count .gt. 0

  endif ! CROSS

  if( CROSS) then

    call Random_Number(cff) ! uniform random number generator

    ! i_Male_Tree   =1+int(cff*float(n_Trees-1))  ! pick a location from 1 to n_Trees
    ! i_Female_Tree =1+int(cff*float(n_Trees-1))  ! pick a location from 1 to n_Trees

    ! move the selected trees from the selected 
    ! male and female individuals that are to be randomly swapped

    Parent_Tree_Swap_Node_Type(1:n_Nodes,1) = &
            GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,i_Male_Tree)
    Parent_Tree_Swap_Node_Type(1:n_Nodes,2) = &
            GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1:n_Nodes,i_Female_Tree)


    !   perform the random tree swap
    !off    call GP_Tree_Swap

    !   move one of the swapped trees into GP_Child_Population_Node_Type

    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,i_Male_Tree) = &
                       Parent_Tree_Swap_Node_Type(1:n_Nodes,1)

    Run_GP_Calculate_Fitness(i_GP_Individual)=.true.

  endif ! CROSS

enddo !  i_GP_Crossover

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Tournament_Style_Sexual_Reproduction
!234567890123456789012345678901234567890123456789012345678901234567890
