subroutine comp_GP_child_indiv_sse( )

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


real (kind=8) :: dff
real (kind=8) :: xcount

!----------------------------------------------------------------------------------------


! compute GP_Child_Individual_SSE(i_GP_Individual)

if( myid == 0 )then
    write(GP_print_unit,'(//A)') &
    'cgcsse: i_node, i_tree, GP_Indiv_Node_Type, GP_Node_Type_Answer  dff           xcount'
endif ! myid == 0

dff = 0.0d0
xcount = 0.0d0

do  i_Node=1,n_Nodes
    do  i_Tree=1,n_Trees

        !if( myid == 0 )then
        !    write(GP_print_unit,'(A,4(1x,I6),1x,E15.7)') &
        !     'cgcsse: i_node, i_tree, GP_Indiv_Node_Type, GP_Node_Type_Answer ',&
        !         i_node, i_tree, GP_Individual_Node_Type(i_Node,i_Tree), &
        !         GP_Node_Type_Answer(i_Node,i_Tree)
        !endif ! myid == 0

        if( GP_Individual_Node_Type(i_Node,i_Tree) .ne. &
            GP_Node_Type_Answer(i_Node,i_Tree)             ) then

            if( GP_Node_Type_Answer(i_Node,i_Tree)     .eq. -9999  .or. &
                GP_Individual_Node_Type(i_Node,i_Tree) .eq. -9999         ) then

                dff = 1.0d2        ! NOTE: larger  jjm 20130415
                xcount=xcount+dff

                !if( myid == 0 )then
                !    write(GP_print_unit,'(A,4(1x,I6),2(1x,E15.7))') &
                !     'cgcsse: i_node, i_tree, &GP_Indiv_Node_Type, &
                !     &GP_Node_Type_Answer  dff, xcount ',&
                !     i_node, i_tree, GP_Individual_Node_Type(i_Node,i_Tree), &
                !     GP_Node_Type_Answer(i_Node,i_Tree), dff, xcount 
                !endif ! myid == 0

            else

                dff = dble( ( GP_Individual_Node_Type(i_Node,i_Tree) -             &
                                  GP_Node_Type_Answer(i_Node,i_Tree)    )**2   )
                xcount=xcount+dff

                !if( myid == 0 )then
                !    write(GP_print_unit,'(A,4(1x,I6),2(1x,E15.7))') &
                !     'cgcsse: i_node, i_tree, GP_Indiv_Node_Type, &
                !        &GP_Node_Type_Answer  dff, xcount ',&
                !         i_node, i_tree, &
                !         GP_Individual_Node_Type(i_Node,i_Tree), &
                !         GP_Node_Type_Answer(i_Node,i_Tree), dff, xcount 
                !endif ! myid == 0

                if( myid == 0 )then
                    write(GP_print_unit,'(5x,3(1x,I6),15X,I6,10x,2(1x,E15.7))') &
                          i_node, i_tree, &
                          GP_Individual_Node_Type(i_Node,i_Tree), &
                          GP_Node_Type_Answer(i_Node,i_Tree), dff, xcount 
                endif ! myid == 0


                !if( myid == 0 )then
                !    write(GP_print_unit,'(A,2(1x,I6))') &
                !          'cgcsse: i_GP_generation, i_GP_Individual ', &
                !              i_GP_generation, i_GP_Individual
                !    write(GP_print_unit,'(A,4(1x,I6),1x,E15.7)') &
                !      'cgcsse: i_node, i_tree, GP_Indiv_Node_Type, &
                !         &GP_Node_Type_Answer, dff ',&
                !          i_node, i_tree, &
                !          GP_Individual_Node_Type(i_Node,i_Tree), &
                !          GP_Node_Type_Answer(i_Node,i_Tree), dff
                !endif ! myid == 0

            endif !   GP_Node_Type_Answer(i_Node,i_Tree) .eq. -9999 ...

        endif !   GP_Individual_Node_Type(i_Node,i_Tree) .ne. ...

    enddo ! i_tree
enddo  ! i_node


GP_Child_Individual_SSE(i_GP_Individual) = xcount
GP_Child_Population_SSE(i_GP_Individual) = xcount   ! ???

!---------------------------------------------------------------------------------

if( myid == 0 )then
    write(GP_print_unit,'(/A)') &
    'cgcsse:------------------------------------------------------------------------------'

    write(GP_print_unit,'(/A,2(1x,I6),1x,E15.7)') &
          'cgcsse: i_GP_gen, i_GP_Indiv, GP_Child_Indiv_SSE(i_GP_Indiv)', &
                   i_GP_generation, i_GP_Individual, &
                   GP_Child_Individual_SSE(i_GP_Individual)
    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
          'cgcsse: i_GP_gen, i_GP_Indiv, GP_Child_Pop_SSE(i_GP_Indiv)  ', &
                   i_GP_generation, i_GP_Individual, &
                   GP_Child_Population_SSE(i_GP_Individual)
endif ! myid == 0

!off GP_Child_Individual_SSE(i_GP_individual)=GA_Individual_Lowest_SSE  ! from the GA_lmdif routine output




return

end subroutine comp_GP_child_indiv_sse
