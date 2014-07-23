subroutine GP_Tree_Build( i_Error )
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
use kinds_mod 
use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff

integer(kind=i4b) :: i_GP_individual
integer(kind=i4b) :: i_Error
integer(kind=i4b) :: i_Node
integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Level
integer(kind=i4b) :: n_Nodes_at_Level
integer(kind=i4b) :: i_Level_Node
integer(kind=i4b) :: Node_Function
integer(kind=i4b) :: Node_Variable
integer(kind=i4b) :: test_function_index
integer(kind=i4b) :: n_parms
integer(kind=i4b) :: n_parms_per_tree


!real(kind=i4b),parameter :: prob_forcing = 0.10  ! 0.20   ! 0.045
!real(kind=i4b) :: prob_forcing = 0.10  ! 0.20   ! 0.045
real(kind=i4b),parameter :: prob_choose_forcing_type = 0.25  ! DO NOT CHANGE

real(kind=i4b),dimension(2,4) :: frac_forcing_type

data frac_forcing_type(1,1) / 0.578947 /
data frac_forcing_type(2,1) / 1.00000  /

data frac_forcing_type(1,2) / 0.210526 /
data frac_forcing_type(2,2) / 0.578947 /

data frac_forcing_type(1,3) / 0.0      /
data frac_forcing_type(2,3) / 0.052632 /

data frac_forcing_type(1,4) / 0.052632 /
data frac_forcing_type(2,4) / 0.210526 /

integer(kind=i4b) :: iforce

!-----------------------------------------------------------------------------



! allow user to turn off forcing with input card "no_forcing"   ! OLD VERSION
!prob_forcing = 0.10
!if( L_no_forcing ) then
!    prob_forcing = 0.00
!endif !  L_no_forcing 


write(GP_print_unit,'(A,1x,E15.7 )') 'gtb: prob_forcing ', prob_forcing



GP_Child_Population_Node_Type=-9999 ! set all to null [-9999]



!write(GP_print_unit,'(/A,1x,I6)')    'gtb: n_GP_Individuals ', n_GP_Individuals
!write(GP_print_unit,'(A,2(1x,I6)/)') 'gtb: n_trees, n_levels', n_trees, n_levels

!write(GP_print_unit,'(A,1x,E15.7 )') &
!      'gtb: GP_Set_Terminal_to_Parameter_Probability', &
!            GP_Set_Terminal_to_Parameter_Probability

do  i_GP_Individual=1,n_GP_Individuals  ! for each GP individual

    !write(GP_print_unit,'(/A,1(1x,I6))') &
    !             'gtb: i_GP_individual ', i_GP_individual


    do  i_Tree=1,n_Trees                ! for each GPCODE tree


        !write(GP_print_unit,'(A,1(1x,I6))') 'gtb: i_Tree ', i_Tree


        call random_number(cff) ! uniform random number generator

        !write(GP_print_unit,'(/A,2(1x,E15.7))') 'gtb: cff, GP_Tree_Probability ', &
        !                                              cff, GP_Tree_Probability

        if( cff .le. GP_Tree_Probability ) then  ! go ahead - put in an equation

            ! always set the first node to zero
            GP_Child_Population_Node_Type(1,i_Tree,i_GP_Individual)=0

            i_Node=0
            level_loop:&
            do  i_Level=1,n_Levels-1                    !original

                n_Nodes_at_Level= pow2_table( i_level-1 ) + 1 ! int(2**(i_Level-1))


                !write(GP_print_unit,'(/A,2(1x,I6)/)') &
                ! 'gtb: i_level, n_Nodes_at_Level ', &
                !       i_level, n_Nodes_at_Level

                do  i_Level_Node=1,n_Nodes_at_Level

                    i_Node=i_Node+1

                    !write(GP_print_unit,'(/A,2(1x,I6)/)') &
                    !      'gtb: i_node , n_Nodes ', i_node , n_Nodes

                    if( i_node > n_nodes ) exit level_loop


                    if( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)  .eq. 0 ) then


                        ! randomly decide function or terminal

                        call random_number(cff) ! uniform random number generator


                        !write(GP_print_unit,'(A,1x,I6,2(1x,F10.4))') &
                        !      'gtb: i_Level, cff, Node_Probability(i_Level)', &
                        !            i_Level, cff, Node_Probability(i_Level)


                        if( cff .lt. Node_Probability(i_Level) ) then  ! set as a function


                            ! new random number to choose the function

                            call random_number(cff) ! uniform random number generator


                            if( L_node_functions )then

                                node_function=1+int(cff*float(n_Node_Functions))

                                !write(GP_print_unit,'(A,1x,I6,1x,F10.4)') &
                                !      'gtb:1 Node_Function, cff ', &
                                !             Node_Function, cff
                                !write(GP_print_unit,'(A,1x,I6)') 'gtb:1 n_Node_Functions', &
                                !                                        n_Node_Functions

                                Node_Function = min( Node_Function, n_Node_Functions )

                                !write(GP_print_unit,'(A,1x,I6)') &
                                !       'gtb:1 Node_Function', Node_Function

                            else


                                test_function_index = 1+int(cff*float(n_functions_input))
                                test_function_index = max( 1, test_function_index  )
                                test_function_index = min( n_functions_input, test_function_index  )

                                !write(GP_print_unit,'(A,1x,I6,1x,F10.4)') &
                                !      'gtb:2 n_functions_input, cff ', &
                                !             n_functions_input, cff
                                !write(GP_print_unit,'(A,1x,I6)')  &
                                  !  'gtb:2 test_function_index',  test_function_index


                                node_function = selected_functions( test_function_index )

                                !write(GP_print_unit,'(A,1x,I6)') &
                                !       'gtb:2 Node_Function', Node_Function

                            endif ! L_node_functions

                            !--------------------------------------------------------------------

                            GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) =  &
                                                                             Node_Function
                            !--------------------------------------------------------------------

                            !  set the node vs terminal selection capability
                            !  for the node inputs at the next level

                            !write(GP_print_unit,'(/A,3(1x,I6)/)') &
                            !  'gtb: i_node, i_Level, N_Levels ', &
                            !        i_node, i_Level, N_Levels

                            if( i_Level .lt. N_Levels-1 ) then

                                ! set the node lowel level inputs to open

                                GP_Child_Population_Node_Type( min(2*i_Node, n_nodes)  , &
                                                                       i_Tree,i_GP_Individual) = 0
                                GP_Child_Population_Node_Type( min(2*i_Node+1, n_nodes), &
                                                                       i_Tree,i_GP_Individual) = 0

                            else

                                ! complete setting the node lowest level nodes with terminals

                                GP_Child_Population_Node_Type( min(2*i_Node, n_nodes)  , &
                                                                       i_Tree,i_GP_Individual) = -1
                                GP_Child_Population_Node_Type( min(2*i_Node+1, n_nodes), &
                                                                       i_Tree,i_GP_Individual) = -1


                            endif !   i_Level .lt. N_Levels-1



                        else


                            ! set it as a Parameter or Variable at a later point in the code

                            !  cff >=   Node_Probability(i_Level)
                            !  so set a parameter or variable later

                            GP_Child_Population_Node_Type(i_Node,i_Tree, i_GP_Individual)=-1


                        endif !   cff .lt. Node_Probability(i_Level)

                    endif !  GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0

                enddo !  i_Level_Node

            enddo level_loop  !  i_Level

        endif !   cff .le. GP_Tree_Probability

    enddo !  i_Tree

    !write(GP_print_unit,'(/A,3(1x,I6)/)') &
    !      'gtb:1 i_GP_individual, n_trees, n_nodes ',  i_GP_individual, n_trees, n_nodes
    !write(GP_print_unit,'(A)') &
    !      'gtb: i_tree, i_node, GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual)'
    !do  i_Tree=1,n_Trees
    !    do  i_Node=1,n_Nodes
    !        if( GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual) /= -9999 )then
    !            write(GP_print_unit,'(3(1x,I8))') &
    !                i_tree, i_node, GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual)
    !        endif ! GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) /= -9999
    !    enddo ! i_node
    !enddo ! i_tree
    !write(GP_print_unit,'(/A)') ' '

enddo !  i_GP_Individual


!------------------------------------------------------------------------------------------------


! randomly fill the terminals of the GP_Child_Population_Node_Type array
! with parameter or variable 'types'


do  i_GP_Individual=1,n_GP_Individuals

    !write(GP_print_unit,'(/A,1(1x,I6))') &
    !             'gtb:2 i_GP_individual ', i_GP_individual

    n_parms = 0

    do  i_Tree=1,n_Trees

        n_parms_per_tree = 0

        i_Node=0

        level_loop2:&
        do  i_Level=1,n_Levels

            n_Nodes_at_Level = pow2_table( i_level - 1 ) + 1  ! int(2**(i_Level-1))


            do  i_Level_Node = 1,n_Nodes_at_Level

                i_Node=i_Node+1

                if( i_node > n_nodes ) exit level_loop2


                if( GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1) then

                    call random_number(cff)   ! uniform random number generator

                    !write(GP_print_unit,'(A,1x,I6,2(1x,F10.4))') &
                    !          'gtb:2 i_Level, cff, GP_Set_Terminal_to_Parameter_Probability', &
                    !                 i_Level, cff, GP_Set_Terminal_to_Parameter_Probability

                    if( cff .le. GP_Set_Terminal_to_Parameter_Probability ) then


                        ! Set the Terminal to a Variable

                        call random_number(cff) ! uniform random number generator

                        ! One of the OBSERVATIONS, one for each equations N, P, Z, etc.

                        Node_Variable=1+int(cff*float(n_CODE_Equations))

                        Node_Variable = min( Node_Variable, n_CODE_Equations )

                        !write(GP_print_unit,'(A,1x,E15.7, 2(1x,I6))') &
                        !      'gtb:2 cff, Node_Variable, n_CODE_Equations', &
                        !             cff, Node_Variable, n_CODE_Equations

                        GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = &
                                                                          -Node_Variable

                        !----------------------------------------------------------------------

                        if( model == 'fasham' )then

                            !  set some variables to the forcing functions -5001 -> -5004

                            node_variable = 0
                            call set_forcing_node( node_variable )

                            !    call random_number(cff)

                            !!write(GP_print_unit,'(A,2(1x,E15.7))') &
                            !!      'gtb:3 cff, prob_forcing', cff, prob_forcing

                            !if( cff < prob_forcing )then

                                !call random_number(cff)

                                !!write(GP_print_unit,'(A,2(1x,E15.7))') &
                                !!      'gtb:4 cff, prob_choose_forcing_type', &
                                !!             cff, prob_choose_forcing_type


                                !node_variable = 0
                                !!if( cff < prob_choose_forcing_type ) node_variable = -5001
                                !!
                                !!do  iforce=1,3
                                !!
                                !!if( prob_choose_forcing_type * float(iforce) < cff .and. &
                                !!prob_choose_forcing_type * float(iforce+1) >= cff  )then
                                !!node_variable = -1 * (5000 + iforce + 1 )
                                !!endif ! prob_choose_forcing_type * float(iforce) < cff...
                                !!
                                !!enddo ! iforce
                                !!
                                !!! in case cff is very close to 1.0000
                                !!if( node_variable == 0 ) node_variable = -5004

                                !!if( cff < frac_forcing_type(3) )then
                                !!    node_variable = -5003
                                !!elseif( cff < frac_forcing_type(4) ) then
                                !!    node_variable = -5004
                                !!elseif( cff < frac_forcing_type(2) ) then
                                !!    node_variable = -5002
                                !!else
                                !!    node_variable = -5001
                                !!endif ! cff < frac_forcing_type(3)

                                !if( cff >  frac_forcing_type(1, 3) .and.  &
                                !    cff <= frac_forcing_type(2, 3)         )then

                                !    node_variable = -5003

                                !elseif( cff >  frac_forcing_type(1, 4) .and.  &
                                !        cff <= frac_forcing_type(2, 4)         )then

                                !    node_variable = -5004

                                !    ! turn off the daily forcing only
                                !    if( L_no_forcing ) then
                                !        node_variable = 0
                                !    endif ! L_no_forcing 

                                !elseif( cff >  frac_forcing_type(1, 2) .and.  &
                                !        cff <= frac_forcing_type(2, 2)         )then

                                !    node_variable = -5002

                                !elseif( cff >  frac_forcing_type(1, 1) .and.  &
                                !        cff <= frac_forcing_type(2, 1)         )then

                                !    node_variable = -5001

                                !endif ! cff < frac_forcing_type(1,3) ...



                                !write(GP_print_unit,'(A,2(1x,I6))') &
                                !      'gtb:4 node_variable', node_variable

                            !endif !  cff < prob_forcing

                            GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = &
                                                                               Node_Variable
                            !----------------------------------------------------------------------

                            !write(GP_print_unit,'(A,4(1x,I6))') &
                            !    'gtb:5 i_GP_Individual, i_Tree, i_Node, &
                            !        &GP_Child_Population_Node_Type', &
                            !           i_GP_Individual, i_Tree, i_Node, &
                            !         GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

                            !----------------------------------------------------------------------

                        endif ! model == 'fasham'

                    else  !   cff > GP_Set_Terminal_to_Parameter_Probability


                        ! set as a random parameter


                        ! Setting GP_Child_Population_Node_Type to zero
                        ! allows the parameters to be set in GA_lmdif


                        GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = 0
                        n_parms = n_parms + 1
                        n_parms_per_tree = n_parms_per_tree + 1


                        ! if there are too many parameters, set subsequent parameter nodes to undefined

                        if( n_parms > n_maximum_number_parameters )then
                            GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) = -9999
                        endif !   n_parms > n_maximum_number_parameters

                    endif !   cff .le. GP_Set_Terminal_to_Parameter_Probability

                endif !   GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. -1

            enddo !  i_Level_Node
        enddo level_loop2 !  i_Level

        !write(GP_print_unit,'(A,5(1x,I6))') &
        !      'gtb:2 i_GP_individual, i_tree, n_parms_per_tree, n_parms, max_parms ',  &
        !             i_GP_individual, i_tree, n_parms_per_tree, n_parms, n_maximum_number_parameters

    enddo !  i_Tree



    !------------------------------------------------------------------------

    !write(GP_print_unit,'(/A,3(1x,I6))') &
    !      'gtb:2 i_GP_individual, n_trees, n_nodes ',  i_GP_individual, n_trees, n_nodes

    !write(GP_print_unit,'(A)') &
    !      'gtb: i_tree, i_node, &
    !              &GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual)'
    !do  i_Tree=1,n_Trees
    !    do  i_Node=1,n_Nodes
    !        if( GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual) /= -9999 )then
    !            write(GP_print_unit,'(3(1x,I8))') &
    !                i_tree, i_node, GP_Child_Population_Node_Type( i_Node,i_Tree,i_GP_Individual)
    !        endif ! GP_Child_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) /= -9999
    !    enddo ! i_node
    !enddo ! i_tree
    !write(GP_print_unit,'(/A)') ' '


    !---------------------------------------------------------------------------------

    call GP_Check_Terminals(&
         GP_Child_Population_Node_Type( 1, 1, i_GP_Individual) , i_Error )


    if( i_Error .eq. 1 ) then
        write(6,'(/A)') 'gtb: GP_Check_Error in GP_Tree_Build'
        write(6,'(A,2(1x,I6)/)') 'gtb: i_GP_Individual, i_Error  ', &
                                       i_GP_Individual, i_Error
        return
    endif !   i_Error .eq. 1

enddo !  i_GP_Individual



GP_Adult_Population_Node_Type=GP_Child_Population_Node_Type




return

end subroutine GP_Tree_Build
