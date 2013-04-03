subroutine GP_Mutations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Mutations in this subroutine are targeted to the nodes only.  
! The terminals are optimized later on using GA_lmdif.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters
use GA_Parameters
use GP_Variables
use GA_Variables

implicit none

real(kind=8) :: cff

integer(kind=4) :: i_GP_Mutation
integer(kind=4) :: i_Tree_Mutation
integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node
integer(kind=4) :: i_Swap_Node
integer(kind=4) :: icnt_Nodes
integer(kind=4) :: icnt
integer(kind=4) :: i_GP_Individual_Mutation
integer(kind=4) :: i_Node_Mutation
integer(kind=4) :: Node_Function
integer(kind=4) :: i_GP_Individual

logical Node_Not_Found

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

i_GP_Individual=n_GP_Elitists+n_GP_Asexual_Reproductions+n_GP_Crossovers


do i_GP_Mutation=1,n_GP_Mutations

  i_GP_Individual=i_GP_Individual+1

! fill in the chosen Child with the Parent's tree information
  do i_Node=1,n_Nodes
    do i_Tree=1,n_Trees
      GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)=GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)
    enddo ! i_tree
  enddo ! i_node


  call Random_Number(cff) ! uniform random number generator
  i_GP_Individual_Mutation=1+int(cff*float(n_GP_Individuals-1))       ! randomly pick one of the n_GP_Individuals to mutate

  call random_number(cff)                                             ! uniform random number generator
  i_Tree_Mutation=1+int(cff*float(n_Trees-1))                         ! randomly pick one of the equation trees


! count the number of nodes on the tree selected for a mutation.  Only nodes are mutated.

  icnt_Nodes=0
  do i_Node=1,n_Nodes
    if( GP_Adult_Population_Node_Type(i_GP_Individual_Mutation,i_Node,i_Tree_Mutation) .gt. 0) then
      icnt_Nodes=icnt_Nodes+1
    endif ! GP_Adult_Population_Node_Type...
  enddo ! i_node

! look to see if there are actually any nodes to mutate


  if( icnt_Nodes .gt. 0) then

    !   randomly choose a node location to mutate
    call random_number(cff) ! uniform random number generator
    i_Swap_Node=1+int(float(icnt_Nodes-1)*cff)


    icnt=0
    Node_Not_Found=.true.
    do i_Node=1,n_Nodes
      if (Node_Not_Found) then
        if (GP_Adult_Population_Node_Type(i_GP_Individual_Mutation,i_Node,i_Tree_Mutation) .ge. 0) then
          icnt=icnt+1
          if (icnt .eq. i_Swap_Node) then
            i_Swap_Node=i_Node
            Node_Not_Found=.false.
          endif !   icnt .eq. i_Swap_Node
        endif !   GP_Adult_Population_Node_Type...
      endif !   Node_Not_Found

    enddo ! i_node


!   fill in the child node with the mutation

    call random_number(cff) ! uniform random number generator
    node_function=1+int(float(n_Node_Functions-1)*cff)

    GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)=Node_Function

    Run_GP_Calculate_Fitness(i_GP_Individual)=.true.

  endif !   icnt_Nodes .gt. 0

enddo !  i_GP_Mutation

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Mutations
!234567890123456789012345678901234567890123456789012345678901234567890
