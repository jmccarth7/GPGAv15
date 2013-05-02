subroutine GP_Tournament_Style_Sexual_Reproduction

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! randomly choose two 'parents' using the Tournament-Style Selection
! and cross the parameter strings to create two new 'children' parameter strings

! modifies    GP_Child_Population_Node_Type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use mpi
use mpi_module
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff

integer(kind=4) :: i_GP_Crossover
integer(kind=4),dimension(2) :: k_GP_Individual_Male
integer(kind=4),dimension(2) :: k_GP_Individual_Female


integer(kind=4) :: i_Male_Tree
integer(kind=4) :: i_Female_Tree
integer(kind=4) :: i_Error

integer(kind=4) :: i_GP_individual                                                                           

!----------------------------------------------------------------------------------

i_GP_Individual = n_GP_Elitists + n_GP_Asexual_Reproductions

if( myid == 0 )then

    write(GP_print_unit,'(/A,1x,I6)' ) &
          'gptssr: n_GP_Crossovers ', n_GP_Crossovers

    write(GP_print_unit,'(A,2(1x,I6))' ) &
          'gptssr: n_GP_Elitists, n_GP_Asexual_Reproductions ', &
                   n_GP_Elitists, n_GP_Asexual_Reproductions
    write(GP_print_unit,'(A,1x,I6)' ) &
          'gptssr: start i_GP_individual  =  ', &
                   n_GP_Elitists + n_GP_Asexual_Reproductions +1

endif ! myid == 0


do  i_GP_Crossover = 1,n_GP_Crossovers


    i_GP_Individual = i_GP_Individual+1

    !if( myid == 0 )then
    !    write(GP_print_unit,'(/A,2(1x,I6)/)' ) &
    !          'gptssr: i_GP_Crossover, i_GP_Individual ', &
    !                   i_GP_Crossover, i_GP_Individual
    !endif ! myid == 0

    !----------------------------------------------------------------------

    ! pick the male parent for sexual crossing of parameter strings

    call random_number(cff) ! uniform random number generator
    k_GP_Individual_Male(1) = 1+int(cff*float(n_GP_Individuals))  ! -1))

    call random_number(cff) ! uniform random number generator
    k_GP_Individual_Male(2) = 1+int(cff*float(n_GP_Individuals)) !  -1))

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

        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !          'gptssr: k_GP_Indiv_Male(1), GP_Adult_Pop_SSE(k_GP_Indiv_Male(1)) ', &
        !                   k_GP_Individual_Male(1), &
        !                   GP_Adult_Population_SSE(k_GP_Individual_Male(1))
        !    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !          'gptssr: k_GP_Indiv_Male(2), GP_Adult_Pop_SSE(k_GP_Indiv_Male(2)) ', &
        !                   k_GP_Individual_Male(2), &
        !                   GP_Adult_Population_SSE(k_GP_Individual_Male(2))
        !endif ! myid == 0

        k_GP_Individual_Male(1) = k_GP_Individual_Male(2)

    endif !   GP_Adult_Population_SSE(k_GP_Individual_Male(2)) .lt....


    !----------------------------------------------------------------------

    ! pick the female parent for sexual crossing of parent parameter strings

    call Random_Number(cff) ! uniform random number generator
    k_GP_Individual_Female(1)  =  1+int(cff*float(n_GP_Individuals))  ! -1))

    call Random_Number(cff) ! uniform random number generator
    k_GP_Individual_Female(2)  =  1+int(cff*float(n_GP_Individuals)) !  -1))

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

        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !          'gptssr: k_GP_Indiv_Female(1), GP_Adult_Pop_SSE(k_GP_Indiv_Female(1)) ', &
        !                   k_GP_Individual_Female(1), &
        !                   GP_Adult_Population_SSE(k_GP_Individual_Female(1))
        !    write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !          'gptssr: k_GP_Indiv_Female(2), GP_Adult_Pop_SSE(k_GP_Indiv_Female(2)) ', &
        !                   k_GP_Individual_Female(2), &
        !                   GP_Adult_Population_SSE(k_GP_Individual_Female(2))
        !endif ! myid == 0

        k_GP_Individual_Female(1)  =  k_GP_Individual_Female(2)

    endif !   GP_Adult_Population_SSE(k_GP_Individual_Female(2)) ...


    !----------------------------------------------------------------------

    ! Randomly choose the tree structure location from the best male
    ! to participate in the genetic crossovers

    ! randomly choose which tree structures from the male and female GP_CODEs
    ! will participate in the genetic crossovers
    ! Find out how many trees there are in each GP_CODE


    call Random_Number(cff) ! uniform random number generator
    i_Male_Tree=1+int(cff*float(n_Trees))  ! pick a tree

    call Random_Number(cff) ! uniform random number generator
    i_Female_Tree=1+int(cff*float(n_Trees))  ! pick a tree


    ! stick the entire chosen male node/tree set into the new child node/tree set

    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees) =  &
           GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,1:n_Trees)

    !if( myid == 0 )then
    !    write(GP_print_unit,'(A,1x,I6/(10(1x,I6)))' ) &
    !          'gptssr: k_GP_Indiv_Male(1), &
    !          &GP_Adult_Pop_Node_Type(k_GP_Indiv_Male(1),1:n_Nodes,1:n_Trees)', &
    !          k_GP_Individual_Male(1), &
    !          GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,1:n_Trees)
    !    write(GP_print_unit,'(A,1x,I6/(10(1x,I6)))' ) &
    !          'gptssr: i_GP_Individual, &
    !          &GP_Child_Pop_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)', &
    !          i_GP_Individual, &
    !          GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)
    !endif ! myid == 0


    !???! Do the genetic crossovers but only keep the solution
    !???! from one (the male) of the two (male and female) generated child tree


    !!?? if( CROSS) then

    !!??  call Random_Number(cff) ! uniform random number generator

    !!??  ! i_Male_Tree    = 1+int(cff*float(n_Trees-1))  ! pick a location from 1 to n_Trees
    !!??  ! i_Female_Tree  = 1+int(cff*float(n_Trees-1))  ! pick a location from 1 to n_Trees

    ! move the selected trees from the selected
    ! male and female individuals that are to be randomly swapped

    Parent_Tree_Swap_Node_Type(1:n_Nodes,1)  =  &
        GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,i_Male_Tree)

    Parent_Tree_Swap_Node_Type(1:n_Nodes,2)  =  &
        GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1:n_Nodes,i_Female_Tree)

    !-----------------------------------------------------------------------------------------

    GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) =  &
        GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),1:n_Nodes,1:n_Trees)

    !-----------------------------------------------------------------------------------------
    !do  i_Tree = 1,n_Trees
    !    do  i_Node = 1,n_Nodes
    !        if( GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),i_Node,i_Tree) /= -9999 )then
    !            if( myid == 0 )then
    !                write(GP_print_unit,'(A,4(1x,I6))' ) &
    !                      'gptssr: k_GP_Indiv_Male(1),i_node,i_tree, &
    !                      &GP_Adult_Pop_Node_Type(k_GP_Indiv_Male(1),i_Node,i_Tree)', &
    !                      k_GP_Individual_Male(1),i_node,i_tree, &
    !                      GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),i_Node,i_Tree)
    !            endif ! myid == 0
    !        endif ! GP_Adult_Population_Node_Type(k_GP_Individual_Male(1),i_Node,i_Tree) /= -9999
    !    enddo !  i_Node
    !enddo ! i_Tree
    !-----------------------------------------------------------------------------------------

    call GP_Check_Terminals(i_Error)

    if( i_Error .eq. 1) then

        if( myid == 0 )then
            write(6,'(A)')&
                  'gptssr: Pre-GP_Check_Error [Male] in GP_Tournament_Style_Sexual_Reproduction'
            write(6,'(A,3(1x,I6)/)') &
                  'gptssr: i_GP_Individual, k_GP_Individual_Male(1), i_Error  ', &
                           i_GP_Individual, k_GP_Individual_Male(1), i_Error
        endif ! myid == 0

      stop 'GP_Tou check error'
    endif

    !-----------------------------------------------------------------------------------------

    GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) =  &
        GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),1:n_Nodes,1:n_Trees)

    !-----------------------------------------------------------------------------------------
    !do  i_Tree = 1,n_Trees
    !    do  i_Node = 1,n_Nodes
    !        if( GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),i_Node,i_Tree) /= -9999 )then
    !            if( myid == 0 )then
    !                write(GP_print_unit,'(A,4(1x,I6))' ) &
    !                      'gptssr: k_GP_Indiv_Female(1),i_node,i_tree, &
    !                      &GP_Adult_Pop_Node_Type(k_GP_Indiv_Female(1),i_Node,i_Tree)', &
    !                      k_GP_Individual_Female(1),i_node,i_tree, &
    !                      GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),i_Node,i_Tree)
    !            endif ! myid == 0
    !        endif ! GP_Adult_Population_Node_Type(k_GP_Individual_Female(1),i_Node,i_Tree) /= -9999
    !    enddo !  i_Node
    !enddo ! i_Tree
    !-----------------------------------------------------------------------------------------

    call GP_Check_Terminals(i_Error)

    if( i_Error .eq. 1) then
        if( myid == 0 )then
            write(6,'(A)')&
                  'gptssr: Pre-GP_Check_Error [Female] in GP_Tournament_Style_Sexual_Reproduction'
            write(6,'(A,3(1x,I6)/)') &
                  'gptssr: i_GP_Individual, k_GP_Individual_Female(1), i_Error  ', &
                           i_GP_Individual, k_GP_Individual_Female(1), i_Error
        endif ! myid == 0
        stop 'GP_Tou stop error 2'
    endif


    !-----------------------------------------------------------------------------------

    !if( myid == 0 )then
    !    write(6,'(/A/)') 'gptssr: call GP_Tree_Swap '
    !endif ! myid == 0

    call GP_Tree_Swap    !   perform the random tree swap

    !if( myid == 0 )then
    !    write(6,'(/A/)') 'gptssr: aft call GP_Tree_Swap '
    !endif ! myid == 0

    !-----------------------------------------------------------------------------------

    !   move one of the swapped trees into the new child GP_Child_Population_Node_Type

    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,i_Male_Tree)  =  &
                  Parent_Tree_Swap_Node_Type(1:n_Nodes,1)


    !-----------------------------------------------------------------------------------

    GP_Individual_Node_Type(1:n_Nodes,1:n_Trees)  =  &
            GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

    !-----------------------------------------------------------------------------------

    call GP_Check_Terminals(i_Error)

    if( i_Error .eq. 1) then
        if( myid == 0 )then
            write(6,'(A)')&
                  'gptssr: Post-GP_Check_Error in GP_Tournament_Style_Sexual_Reproduction'
            write(6,'(A,3(1x,I6)/)') 'gptssr: i_GP_Individual, i_Male_Tree, i_Error  ', &
                                              i_GP_Individual, i_Male_Tree, i_Error
        endif ! myid == 0
        stop 'GP_Tou stop error 3'
    endif

    !-----------------------------------------------------------------------------------

    Run_GP_Calculate_Fitness(i_GP_Individual) = .true.

    !!?? endif ! CROSS

enddo !  i_GP_Crossover


return


end subroutine GP_Tournament_Style_Sexual_Reproduction
