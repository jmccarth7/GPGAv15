subroutine GP_Mutations

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! Mutations in this subroutine are targeted to the nodes only.
! The terminals are optimized later on using GA_lmdif.

! Modifies  GP_Child_Population_Node_Type

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff

integer(kind=4) :: i_GP_Mutation
integer(kind=4) :: i_GP_Individual_Mutation
integer(kind=4) :: i_Tree_Mutation
integer(kind=4) :: icnt_Nodes
integer(kind=4) :: icnt


integer(kind=4) :: Node_to_Mutate
integer(kind=4) :: Node_Function
integer(kind=4) :: i_Error

logical Node_Not_Found

!------------------------------------------------------------------------



i_GP_Individual_Mutation = 0

i_GP_Individual = n_GP_Elitists + n_GP_Asexual_Reproductions + n_GP_Crossovers

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,I6)' ) &
          'gpmut: n_GP_Mutations ', n_GP_Mutations                                     
    write(GP_print_unit,'(A,3(1x,I6))' ) &
          'gpmut: n_GP_Elitists, n_GP_Asexual_Reproductions, n_GP_Crossovers ', &
                  n_GP_Elitists, n_GP_Asexual_Reproductions, n_GP_Crossovers
    write(GP_print_unit,'(A,1x,I6)' ) &
          'gpmut: start i_GP_individual  =  ', &
             n_GP_Elitists + n_GP_Asexual_Reproductions + n_GP_Crossovers +1
endif ! myid == 0 

do  i_GP_Mutation = 1,n_GP_Mutations

    i_GP_Individual = i_GP_Individual+1

!    ! fill in the chosen Child with the Parent's tree information
!    do  i_Node = 1,n_Nodes
!        do  i_Tree = 1,n_Trees
!            GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)  =  &
!            GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)
!        enddo ! i_tree
!    enddo ! i_node


    ! randomly pick one of the n_GP_Individuals to mutate

    call Random_Number(cff) ! uniform random number generator


    ! randomly choose from the population pool
    ! randomly pick one of the n_GP_Individuals Adults to mutate

    i_GP_Individual_Mutation = 1+int(cff*float(n_GP_Individuals))  ! -1))

! choose sequentially from the best of the population  [SHOWN TO CONVERGE FASTER THAN RANDOMLY CHOSEN]
!off i_GP_Individual_Mutation=i_GP_Individual_Mutation+1

    ! Fill in the Child nodes with the chosen Parent's node/tree information

    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees) =  &
        GP_Adult_Population_Node_Type(i_GP_Individual_Mutation,1:n_Nodes,1:n_Trees)

    GP_Individual_Node_Type(1:n_Nodes,1:n_Trees)  =  &
      GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

    !----------------------------------------------------------------------------------

    call GP_Check_Terminals(i_Error)

    if( i_Error .eq. 1) then
        if( myid == 0 )then
            !write(*,*) 'Pre-GP_Check_Error in GP_Mutation',&
            !         i_GP_Individual,i_GP_Mutation,i_Error
            write(6,'(/A)') 'gpm: Pre-GP_Check_Error in GP_Mutation'
            write(6,'(A,2(1x,I6)/)') 'gpm: i_GP_Individual, i_GP_Mutation, i_Error  ', &
                                           i_GP_Individual, i_GP_Mutation, i_Error
        endif ! myid == 0 
        stop 'GP Mut check error'
    endif

    !----------------------------------------------------------------------------------


    ! Randomly choose which tree to mutate

    call random_number(cff) ! uniform random number generator

    i_Tree_Mutation = 1+int(cff*float(n_Trees)) ! -1))   ! randomly pick one of the equation trees


    ! count the number of function nodes and terminals
    ! on the tree selected for a mutation

    ! Only function nodes are mutated.

    icnt_Nodes = 0

    do  i_Node = 1, n_Nodes

        !!if( GP_Adult_Population_Node_Type( &
        !!    i_GP_Individual_Mutation,i_Node,i_Tree_Mutation) .gt. 0) then

        if( GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree_Mutation) .ne. -9999) then

            icnt_Nodes = icnt_Nodes+1

        endif ! GP_Child_Population_Node_Type...

    enddo ! i_node


    ! look to see if there are actually any nodes to mutate

    ! if there are nodes to mutate (i.e. icnt_Nodes > 0),
    ! randomly pick one and give it a randomly chosen (mutated) new node

    if( icnt_Nodes .gt. 0) then

        !   randomly choose a node to mutate

        call random_number(cff) ! uniform random number generator
        !!??Node_to_Mutate = 1+int(float(icnt_Nodes-1)*cff)
        Node_to_Mutate = 1+int(float(icnt_Nodes)*cff)


        icnt = 0
        Node_Not_Found = .true.
        do  i_Node = 1,n_Nodes

            !!!if( GP_Adult_Population_Node_Type( i_GP_Individual_Mutation,i_Node,i_Tree_Mutation) .ge. 0)then

            if( GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree_Mutation) .ne. -9999) then

                ! this is a node with a function value
                icnt = icnt+1
                if( icnt .eq. Node_to_Mutate) then
                    Node_to_Mutate = i_Node
                    exit
                endif !   icnt .eq. Node_to_Mutate
            endif !   GP_Adult_Population_Node_Type...
        enddo ! i_node


        !   fill in the child node with the randomly chosen node function mutation

        call random_number(cff) ! uniform random number generator

        !!node_function = 1+int(float(n_Node_Functions-1)*cff)

        if( GP_Child_Population_Node_Type(i_GP_Individual,Node_to_Mutate,i_Tree_Mutation) .le. 0 ) then                  

            ! [Ranges from: -n_CODE_Equations to 0]              

            Node_Function = -1 *  int(cff*float(n_CODE_Equations + 1 ))  

            Node_Function = max( Node_Function , -n_CODE_Equations )  

        else                                                                                                            

            ! [Ranges from: 1 to n_Node_Functions]               

            Node_Function =  1 +  int(cff*float( n_Node_Functions ))       

            Node_Function =  min( Node_Function, n_Node_Functions )       

        endif                                                                    

        GP_Child_Population_Node_Type( i_GP_Individual, Node_to_Mutate, i_Tree_Mutation) = Node_Function

        !orig Run_GP_Calculate_Fitness(i_GP_Individual) = .true.

    endif !   icnt_Nodes .gt. 0


    Run_GP_Calculate_Fitness(i_GP_Individual) = .true.




    GP_Individual_Node_Type(1:n_Nodes,1:n_Trees)  =  &
    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)


    !----------------------------------------------------------------------------------

    call GP_Check_Terminals(i_Error)

    if( i_Error .eq. 1) then
        if( myid == 0 )then
            write(6,'(A)') 'gpm: Post-GP_Check_Error in GP_Mutation'
            write(6,'(A,2(1x,I6)/)') 'gpm: i_GP_Individual, i_GP_Mutation, i_Error  ', &
                                           i_GP_Individual, i_GP_Mutation, i_Error
        endif ! myid == 0
        stop 'GP_Mut check error'
    endif

    !----------------------------------------------------------------------------------


enddo !  i_GP_Mutation


return


end subroutine GP_Mutations
