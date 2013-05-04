subroutine GP_Tree_Build
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module
use GP_Parameters_module
use GP_model_parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff
!real(kind=8) :: dff

integer(kind=4) :: n_Nodes_at_Level
integer(kind=4) :: i_Level_Node
integer(kind=4) :: Node_Function
integer(kind=4) :: Node_Variable
integer(kind=4) :: i_Error

!-----------------------------------------------------------------------------

GP_Child_Population_Node_Type=-9999 ! set all to null [-9999]



!write(GP_print_unit,'(/A,1x,I6/)')    'gtb:  n_GP_Individuals ', n_GP_Individuals
!write(GP_print_unit,'(/A,2(1x,I6)/)') 'gtb:  n_trees, n_levels', n_trees, n_levels


do  i_GP_Individual=1,n_GP_Individuals  ! for each GP individual

    do  i_Tree=1,n_Trees                ! for each GPCODE tree

        call random_number(cff) ! uniform random number generator

        if( cff .le. GP_Tree_Probability ) then  ! go ahead and put in an equation

            ! always set the first node to zero
            GP_Child_Population_Node_Type(i_GP_Individual,1,i_Tree)=0

            i_Node=0
            do  i_Level=1,n_Levels-1

                n_Nodes_at_Level=int(2**(i_Level-1))

                !write(GP_print_unit,'(/A,2(1x,I6)/)') 'gtb: i_level, n_Nodes_at_Level ', &
                !                                            i_level, n_Nodes_at_Level

                do  i_Level_Node=1,n_Nodes_at_Level

                    i_Node=i_Node+1

                    if( GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)  .eq. 0 ) then

                        ! randomly decide function or terminal

                        call random_number(cff) ! uniform random number generator


                        !write(GP_print_unit,'(A,1x,I6,2(1x,F10.4))') &
                        !      'gtb: i_Level, cff, Node_Probability(i_Level)', &
                        !            i_Level, cff, Node_Probability(i_Level)

                        if( cff .lt. Node_Probability(i_Level) ) then  ! set as a terminal

                            call random_number(cff) ! uniform random number generator

                            node_function=1+int(cff*float(n_Node_Functions))

                            !write(GP_print_unit,'(A,1x,I6,1x,F10.4)') 'gtb: Node_Function, cff ', &
                            !                                                Node_Function, cff
                            !write(GP_print_unit,'(A,1x,I6)') 'gtb: n_Node_Functions', &
                            !                                       n_Node_Functions

                            Node_Function = min( Node_Function, n_Node_Functions )

                            !write(GP_print_unit,'(A,1x,I6)') 'gtb: Node_Function', Node_Function


                            GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) =  Node_Function

                            !  set the node vs terminal selection capability
                            !  for the node inputs at the next level

                            !write(GP_print_unit,'(/A,3(1x,I6)/)') &
                            !      'gtb: i_node , i_Level, N_Levels ', &
                            !            i_node , i_Level, N_Levels

                            if( i_Level .lt. N_Levels-1 ) then

                                ! set the node lowel level inputs to open

                                GP_Child_Population_Node_Type(i_GP_Individual, 2*i_Node  ,i_Tree)=0
                                GP_Child_Population_Node_Type(i_GP_Individual, 2*i_Node+1,i_Tree)=0
                            else

                                ! complete setting the node lowest level nodes with terminals

                                GP_Child_Population_Node_Type(i_GP_Individual,2*i_Node   ,i_Tree)=-1
                                GP_Child_Population_Node_Type(i_GP_Individual,2*i_Node+1 ,i_Tree)=-1

                            endif !   i_Level .lt. N_Levels-1


                        else ! set it as a Parameter or Variable at a later point in the code

                            ! set a parameter or variable later

                            GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)=-1


                        endif !   cff .lt. Node_Probability(i_Level)

                    endif !  GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) .eq. 0

                enddo !  i_Level_Node

            enddo !  i_Level

        endif !   cff .le. GP_Tree_Probability

    enddo !  i_Tree

    !write(GP_print_unit,'(/A,3(1x,I6)/)') &
    !      'gtb: i_GP_individual, n_nodes, n_trees ',  i_GP_individual, n_nodes, n_trees
    !do  i_Tree=1,n_Trees
    !    do  i_Node=1,n_Nodes
    !        if( GP_Child_Population_Node_Type(i_GP_Individual, i_Node,i_Tree) /= -9999 )then
    !            write(GP_print_unit,'(A,3(1x,I6))') &
    !                  'gtb: i_node, i_tree, GP_Child_Population_Node_Type(i_GP_Individual, i_Node,i_Tree)', &
    !                        i_node, i_tree, GP_Child_Population_Node_Type(i_GP_Individual, i_Node,i_Tree)
    !        endif ! GP_Child_Population_Node_Type(i_Node,i_Tree) /= -9999
    !    enddo ! i_node
    !enddo ! i_tree
    !write(GP_print_unit,'(/A)') ' '

enddo !  i_GP_Individual

!--------------------------------------------------------------------------------------------------

! randomly fill the terminals of the GP_Child_Population_Node_Type array
! with parameter or variable 'types'


do  i_GP_Individual=1,n_GP_Individuals

    do  i_Tree=1,n_Trees

        i_Node=0
        do  i_Level=1,n_Levels

            n_Nodes_at_Level = int(2**(i_Level-1))

            do  i_Level_Node = 1,n_Nodes_at_Level

                i_Node=i_Node+1

                if( GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) .eq. -1) then

                    call random_number(cff)   ! uniform random number generator

                    if( cff .le. GP_Set_Terminal_to_Parameter_Probability ) then

                        ! Set the Terminal to a Parameter

                        call random_number(cff) ! uniform random number generator

                        ! One of the OBSERVATIONS, one for each equations N, P, Z, etc.

                        Node_Variable=1+int(cff*float(n_CODE_Equations))

                        Node_Variable = min( Node_Variable, n_CODE_Equations )

                        !write(GP_print_unit,'(A,2(1x,I6))') &
                        !      'gtb: Node_Variable, n_CODE_Equations', &
                        !            Node_Variable, n_CODE_Equations


                        GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) = -Node_Variable


                    else  ! set as a random parameter

                        ! The setting to zero allows the parameters to be set in GA_lmdif

                        GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) = 0

                    endif !   cff .le. GP_Set_Terminal_to_Parameter_Probability
                endif !   GP_Child_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) .eq. -1
            enddo !  i_Level_Node
        enddo !  i_Level
    enddo !  i_Tree



    !------------------------------------------------------------------------
    ! loading GP_Individual_Node_Type is done here only because the
    ! subroutine GP_Check_Terminals uses GP_Individual_Node_Type

    ! GP_Individual_Node_Type is over-written each time you go through the
    ! loop on i_GP_Individual

    ! in the main program, later, GP_Individual_Node_Type is loaded with
    ! GP_Adult_Population_Node_Type in the loop in i_GP_Individual
    !------------------------------------------------------------------------

    !write(GP_print_unit,'(/A,3(1x,I6)/)') &
    !      'gtb: i_GP_individual, n_nodes, n_trees ',  i_GP_individual, n_nodes, n_trees
    !do  i_Tree=1,n_Trees
    !    do  i_Node=1,n_Nodes
    !        if( GP_Child_Population_Node_Type(i_GP_Individual, i_Node,i_Tree) /= -9999 )then
    !            write(GP_print_unit,'(A,3(1x,I6))') &
    !                  'gtb: i_node, i_tree, GP_Child_Population_Node_Type(i_GP_Individual, i_Node,i_Tree)', &
    !                        i_node, i_tree, GP_Child_Population_Node_Type(i_GP_Individual, i_Node,i_Tree)
    !        endif ! GP_Child_Population_Node_Type(i_Node,i_Tree) /= -9999
    !    enddo ! i_node
    !enddo ! i_tree
    !write(GP_print_unit,'(/A)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> '


    GP_Individual_Node_Type(1:n_Nodes,1:n_Trees) = &
    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees)

!    !---------------------------------------------------------------------------------
!    GP_Individual_Node_Type = -9999   ! debug only
!    GP_Individual_Node_Type(1,1)=4    ! debug only
!    GP_Individual_Node_Type(2,1)=6    ! debug only
!    GP_Individual_Node_Type(3,1)=4    ! debug only
!    GP_Individual_Node_Type(4,1)=2    ! debug only
!    GP_Individual_Node_Type(5,1)=-2    ! debug only
!    GP_Individual_Node_Type(6,1)=3    ! debug only
!    GP_Individual_Node_Type(7,1)=1    ! debug only
!    GP_Individual_Node_Type(8,1)=-2    ! debug only
!    GP_Individual_Node_Type(9,1)=-1    ! debug only
!    GP_Individual_Node_Type(12,1)=-1    ! debug only
!    GP_Individual_Node_Type(13,1)=-2    ! debug only
!    GP_Individual_Node_Type(14,1)=0    ! debug only
!    GP_Individual_Node_Type(15,1)=-1    ! debug only
!    GP_Individual_Node_Type(1,2)=4    ! debug only
!    GP_Individual_Node_Type(2,2)=0    ! debug only
!    GP_Individual_Node_Type(3,2)=7    ! debug only
!    GP_Individual_Node_Type(6,2)=1    ! debug only
!    GP_Individual_Node_Type(7,2)=0    ! debug only
!    GP_Individual_Node_Type(12,2)=-1    ! debug only
!    GP_Individual_Node_Type(13,2)=-1    ! debug only
!    GP_Individual_Node_Type(1,5)=-2    ! debug only
!
!
!    !---------------------------------------------------------------------------------

    !write(GP_print_unit,'(/A,2(1x,I6)/)') &
    !      'gtb: n_nodes, n_trees ', n_nodes, n_trees
    !do  i_Tree=1,n_Trees
    !    do  i_Node=1,n_Nodes
    !        if( GP_Individual_Node_Type(i_Node,i_Tree) /= -9999 )then
    !            write(GP_print_unit,'(A,3(1x,I6))') &
    !                  'gtb: i_node, i_tree, GP_Individual_Node_Type(i_Node,i_Tree)', &
    !                        i_node, i_tree, GP_Individual_Node_Type(i_Node,i_Tree)
    !        endif ! GP_Individual_Node_Type(i_Node,i_Tree) /= -9999
    !    enddo ! i_node
    !enddo ! i_tree
    !write(GP_print_unit,'(/A)') ' '


    call GP_Check_Terminals(i_Error)


    if( i_Error .eq. 1 ) then
        if( myid == 0 )then
            write(6,'(/A)') 'gtb: GP_Check_Error in GP_Tree_Build'
            write(6,'(A,2(1x,I6)/)') 'gtb: i_GP_Individual, i_Error  ', &
                                           i_GP_Individual, i_Error
        endif ! myid == 0 
        stop  'GP Tree Build error'
    endif

    !call MPI_FINALIZE(ierr)  ! debug only
    !stop                     ! debug only


enddo !  i_GP_Individual



GP_Adult_Population_Node_Type=GP_Child_Population_Node_Type



return

end subroutine GP_Tree_Build

