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

real(kind=4) :: cff
real(kind=8) :: dff

integer(kind=4) :: i_GP_Crossover
!integer(kind=4) :: i_GP_Crossover_Point
integer(kind=4),dimension(2) :: k_GP_Individual_Male
integer(kind=4),dimension(2) :: k_GP_Individual_Female


integer(kind=4) :: i_Male_Tree
integer(kind=4) :: i_Female_Tree
integer(kind=4) :: i_Node_Count
integer(kind=4) :: icff
integer(kind=4) :: i_Error

logical CROSS

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_GP_Individual = n_GP_Elitists + n_GP_Asexual_Reproductions


write(GP_print_unit,'(A,2(1x,I6))' ) &
      'gptssr: n_GP_Elitists, n_GP_Asexual_Reproductions ', &
               n_GP_Elitists, n_GP_Asexual_Reproductions
write(GP_print_unit,'(A,1x,I6)' ) &
      'gptssr: start i_GP_individual  =  ', &
               n_GP_Elitists + n_GP_Asexual_Reproductions +1

do  i_GP_Crossover = 1,n_GP_Crossovers


    i_GP_Individual = i_GP_Individual+1


    !----------------------------------------------------------------------

    ! pick the male parent for sexual crossing of parameter strings

    call random_number(cff) ! uniform random number generator
    k_GP_Individual_Male(1) = 1+int(cff*float(n_GP_Individuals-1))

    call random_number(cff) ! uniform random number generator
    k_GP_Individual_Male(2) = 1+int(cff*float(n_GP_Individuals-1))

    ! Check to make sure that the two males are not the same

    if( k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)) then

        if( k_GP_Individual_Male(1) .ne. n_GP_Individuals) then
            k_GP_Individual_Male(2) = k_GP_Individual_Male(1) + 1
        else
            k_GP_Individual_Male(2) = k_GP_Individual_Male(1) - 1
        endif !   k_GP_Individual_Male(1) .ne. n_GP_Individuals

    endif !   k_GP_Individual_Male(2) .eq. k_GP_Individual_Male(1)

    ! select the individual of the two with the best fitness
    ! select the individual with the least SSE level between the two chosen males


    if( GP_Adult_Population_SSE(k_GP_Individual_Male(2)) .lt.  &
        GP_Adult_Population_SSE(k_GP_Individual_Male(1))         ) then
  
        write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
              'gptssr: k_GP_Indiv_Male(1), GP_Adult_Pop_SSE(k_GP_Indiv_Male(1)) ', &
                       k_GP_Individual_Male(1), GP_Adult_Population_SSE(k_GP_Individual_Male(1))
        write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
              'gptssr: k_GP_Indiv_Male(2), GP_Adult_Pop_SSE(k_GP_Indiv_Male(2)) ', &
                       k_GP_Individual_Male(2), GP_Adult_Population_SSE(k_GP_Individual_Male(2))
  
        !write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !      'gptssr: k_GP_Indiv_Male(1), GP_Indiv_Ranked_Fitness(k_GP_Indiv_Male(1)) ', &
        !               k_GP_Individual_Male(1), GP_Individual_Ranked_Fitness(k_GP_Individual_Male(1))
        !write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !      'gptssr: k_GP_Indiv_Male(2), GP_Indiv_Ranked_Fitness(k_GP_Indiv_Male(2)) ', &
        !               k_GP_Individual_Male(2), GP_Individual_Ranked_Fitness(k_GP_Individual_Male(2))
  
    !!if( GP_Individual_Ranked_Fitness(k_GP_Individual_Male(2)) .gt.          &
    !!    GP_Individual_Ranked_Fitness(k_GP_Individual_Male(1))      ) then
  
        k_GP_Individual_Male(1) = k_GP_Individual_Male(2)
  
    endif !   GP_Adult_Population_SSE(k_GP_Individual_Male(2)) .lt....
  
  
    !----------------------------------------------------------------------
  
    ! pick the female parent for sexual crossing of parent parameter strings
  
    call Random_Number(cff) ! uniform random number generator
    k_GP_Individual_Female(1)  =  1+int(cff*float(n_GP_Individuals-1))
  
    call Random_Number(cff) ! uniform random number generator
    k_GP_Individual_Female(2)  =  1+int(cff*float(n_GP_Individuals-1))
  
    ! Check to make sure that the two females are not the same
  
    if( k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)  ) then
  
        if( k_GP_Individual_Female(1) .ne. N_GP_Individuals) then
            k_GP_Individual_Female(2)  =  k_GP_Individual_Female(1)+1
        else
            k_GP_Individual_Female(2)  =  k_GP_Individual_Female(1)-1
        endif !   k_GP_Individual_Female(1) .ne. N_GP_Individuals
  
    endif !   k_GP_Individual_Female(2) .eq. k_GP_Individual_Female(1)
  
  
    ! select the individual with the lowest SSE level between the two chosen females
    if( GP_Adult_Population_SSE(k_GP_Individual_Female(2)) .lt.  &
        GP_Adult_Population_SSE(k_GP_Individual_Female(1))          ) then
        write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
              'gptssr: k_GP_Indiv_Female(1), GP_Adult_Pop_SSE(k_GP_Indiv_Female(1)) ', &
                       k_GP_Individual_Female(1), GP_Adult_Population_SSE(k_GP_Individual_Female(1))
        write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
              'gptssr: k_GP_Indiv_Female(2), GP_Adult_Pop_SSE(k_GP_Indiv_Female(2)) ', &
                       k_GP_Individual_Female(2), GP_Adult_Population_SSE(k_GP_Individual_Female(2))
  
        !write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !      'gptssr: k_GP_Indiv_Female(1), GP_Indiv_Ranked_Fitness(k_GP_Indiv_Female(1)) ', &
        !               k_GP_Individual_Female(1), GP_Individual_Ranked_Fitness(k_GP_Individual_Female(1))
        !write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !      'gptssr: k_GP_Indiv_Female(2), GP_Indiv_Ranked_Fitness(k_GP_Indiv_Female(2)) ', &
        !               k_GP_Individual_Female(2), GP_Individual_Ranked_Fitness(k_GP_Individual_Female(2))
  
        ! select the individual of the two with the best fitness
  
    !!if( GP_Individual_Ranked_Fitness(k_GP_Individual_Female(2)) .gt.        &
    !!    GP_Individual_Ranked_Fitness(k_GP_Individual_Female(1))      ) then
  
        k_GP_Individual_Female(1)  =  k_GP_Individual_Female(2)
  
    endif !   GP_Adult_Population_SSE(k_GP_Individual_Female(2)) ...
  
  
    !----------------------------------------------------------------------
  
    ! Randomly choose the tree structure location from the best male 
    ! to participate in the genetic crossovers
  
    ! randomly choose which tree structures from the male and female GP_CODEs
    ! will participate in the genetic crossovers
    ! Find out how many trees there are in each GP_CODE
  
    i_Node_Count = 0
    do  i_Tree = 1,n_Trees
        if( GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree) .ne. -9999)then
  
            i_Node_Count = i_Node_Count+1
  
        endif !   GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree)...
    enddo ! i_tree
  
    write(GP_print_unit,'(A,1x,I6)') 'gpto: i_Node_Count ', i_Node_Count
  
  
    if( i_Node_Count .gt. 0) then
  
        CROSS = .true.  ! there is at least one Tree structure to cross with
  
        call Random_Number(cff) ! uniform random number generator
        icff = 1+int(cff*float(i_Node_Count-1))  ! pick a tree
  
        i_Node_Count = 0
  
        do  i_Tree = 1,n_Trees
  
            if( GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree) .ne. -9999) then
  
                i_Node_Count = i_Node_Count+1
                if( i_Node_Count .eq. icff) i_Male_Tree  =  i_Tree
  
            endif !   GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1,i_Tree)...
  
        enddo  ! i_tree
  
    else
  
        CROSS = .false. ! there are no Trees on this GP_CODE
  
    endif !   i_Node_Count .gt. 0
  



    ! Randomly choose the tree structure location from the best female 
    ! to participate in the genetic crossovers
  
    if( CROSS) then
  
        i_Node_Count = 0
        do  i_Tree = 1,n_Trees
  
            if( GP_Adult_Population_Node_Type(k_GP_Individual_FeMale(1),1,i_Tree) .ne. -9999)then
                i_Node_Count = i_Node_Count+1
            endif ! GP_Adult_Population_Node_Type(k_GP_Individual_FeMale(1),1,i_Tree) .ne. -9999
  
        enddo ! i_tree
  
  
        if( i_Node_Count .gt. 0) then
  
            CROSS = .true.
  
            call Random_Number(cff) ! uniform random number generator
            icff = 1+int(cff*float(i_Node_Count-1))  ! pick tree
  
            i_Node_Count = 0
            do  i_Tree = 1,n_Trees
  
                if( GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1,i_Tree) .ne. -9999) then
  
                    i_Node_Count = i_Node_Count+1
                    if( i_Node_Count .eq. icff) i_Female_Tree = i_Tree
  
                endif !GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1,i_Tree) ...
  
            enddo ! i_tree
  
        else
  
            CROSS = .false.
  
        endif !   i_Node_Count .gt. 0
  
    endif ! CROSS


    ! stick the entire chosen male node/tree set into the new child node/tree set
    
    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees) =  &
           GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,1:n_Trees)
    
    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
          'gptssr: k_GP_Indiv_Male(1), &
                  &GP_Adult_Pop_Node_Type(k_GP_Indiv_Male(1),1:n_Nodes,1:n_Trees)', &
                   k_GP_Individual_Male(1), &
                   GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,1:n_Trees)
    
    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
          'gptssr: k_GP_Indiv_Male(1), &
                  &GP_Child_Pop_Node_Type(k_GP_Indiv_Male(1),1:n_Nodes,1:n_Trees)', &
                   k_GP_Individual_Male(1), &
                   GP_Child_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,1:n_Trees)
    
    !???! Do the genetic crossovers but only keep the solution 
    !???! from one (the male) of the two (male and female) generated child tree


    if( CROSS) then

        call Random_Number(cff) ! uniform random number generator

        ! i_Male_Tree    = 1+int(cff*float(n_Trees-1))  ! pick a location from 1 to n_Trees
        ! i_Female_Tree  = 1+int(cff*float(n_Trees-1))  ! pick a location from 1 to n_Trees

        ! move the selected trees from the selected
        ! male and female individuals that are to be randomly swapped

        Parent_Tree_Swap_Node_Type(1:n_Nodes,1)  =  &
            GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,i_Male_Tree)
        Parent_Tree_Swap_Node_Type(1:n_Nodes,2)  =  &
            GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1:n_Nodes,i_Female_Tree)

        GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) =  &
            GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,1:n_Trees)

        do  i_Tree = 1,n_Trees
            do  i_Node = 1,n_Nodes

                write(GP_print_unit,'(A,3(1x,I6),1x,E15.7)' ) &
                      'gptssr: k_GP_Indiv_Male(1),i_node,i_tree, &
                              &GP_Adult_Pop_Node_Type(k_GP_Indiv_Male(1),i_Node,i_Tree)', &
                               k_GP_Individual_Male(1),i_node,i_tree, &
                               GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),i_Node,i_Tree)
    
                write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)' ) &
                      'gptssr:i_node,i_tree, &
                              &GP_Indiv_Node_Type(i_Node,i_Tree)', &
                              i_node,i_tree,  &
                               GP_Individual_Node_Type(i_Node,i_Tree)

            enddo !  i_Node
        enddo ! i_Tree

        call GP_Check_Terminals(i_Error)

        if( i_Error .eq. 1) then
            write(*,*) 'Pre-GP_Check_Error [Male] in GP_Tournament_Style_Sexual_Reproduction', &
                     i_GP_Individual,i_Error
          stop
        endif

        GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) =  &
            GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1:n_Nodes,1:n_Trees)

        call GP_Check_Terminals(i_Error)

        if( i_Error .eq. 1) then
            write(*,*) 'Pre-GP_Check_Error [Female] in GP_Tournement_Style_Sexual_Reproduction', &
                  i_GP_Individual,i_Error
            stop
        endif



        call GP_Tree_Swap    !   perform the random tree swap

        !   move one of the swapped trees into the new child GP_Child_Population_Node_Type

        GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,i_Male_Tree)  =  &
                      Parent_Tree_Swap_Node_Type(1:n_Nodes,1)

        GP_Individual_Node_Type(1:n_Nodes,1:n_Trees)  =  &
        GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

        call GP_Check_Terminals(i_Error)

        if( i_Error .eq. 1) then
            write(*,*) 'Post-GP_Check_Error in GP_Tournament_Style_Sexual_Reproduction', &
                    i_GP_Individual,i_Error
            stop
        endif


        Run_GP_Calculate_Fitness(i_GP_Individual) = .true.

    endif ! CROSS

enddo !  i_GP_Crossover

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Tournament_Style_Sexual_Reproduction
!234567890123456789012345678901234567890123456789012345678901234567890
