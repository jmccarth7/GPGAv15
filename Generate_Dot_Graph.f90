! Generate_Dot_Graph
! Written by Erik Wisuri [14 June, 2013]
! These subroutines generate DOT output from a tree array,
! and use the dot executable to create PDF files.

subroutine Generate_Dot_Graph ( Graph_Name, Trees, Parameter_Values, n_Nodes, n_Trees )
implicit none

! Input
character (len=*) :: Graph_Name ! _<tree number> will be appended to this for each file
integer (kind=4),intent(in) :: n_Nodes, n_Trees ! The number of trees and number of nodes in each tree
integer (kind=4), dimension(n_Nodes, n_Trees), intent(in) :: Trees ! The array of trees and nodes
real (kind=4), intent(in) :: Parameter_Values(n_Nodes, n_Trees) ! The values for the tree parameters

! Local variables
integer (kind=4) :: i, j, gFile, Node_Type
character (len=80) :: Graph_File

!-----------------------------------------------------------------------------------------

gFile = 80

do  i = 1,n_Trees

    write(Graph_File, '(A,I0.0)') 'trees/'//trim(Graph_Name)//'_', i  ! orig

    open(gFile, FILE=trim(Graph_File)//'.dot')

    write(gFile,*) 'digraph '//trim(Graph_Name)//' {'
    write(gFile,*) 'splines=false;'

    do  j = 1,n_Nodes

        Node_Type = Trees(j,i)

        if( Node_Type .ne. -9999) then
            if( Node_Type .gt. 0) then
                call Dot_Graph_Function ( gFile, Node_Type, j )
            elseif( Node_Type .lt. 0) then
                write (gFile,'(I0.0,A,I0.1,A)') j, '[label="', Node_Type, '"];'
            else
                write (gFile,*) j, '[label="', Parameter_Values(j,i), '"];'
            endif ! Node_Type .gt. 0

            call Dot_Graph_Hierarchy( gFile, j )

        endif ! Node_Type .ne. -9999

    enddo ! j

    write(gFile,*) '}'
    close(gFile)

    call system( 'dot ' // trim(Graph_File) // &
                 &'.dot -T pdf -Nheight=0.5 -Nwidth=0.02 -o' // trim(Graph_File) // '.pdf' )

enddo ! i 

end subroutine Generate_Dot_Graph


subroutine Dot_Graph_Function ( File, Function_Index, Node_ID )
implicit none

! Input
integer (kind=4), intent(in) :: File, Function_Index, Node_ID

!------------------------------------------------------------------------------

select case (Function_Index)
    case (1)
        write(File,'(I0.0,A)') Node_ID, '[label="+"];'
    case (2)
        write(File,'(I0.0,A)') Node_ID, '[label="-"];'
    case (3)
        write(File,'(I0.0,A)') Node_ID, '[label="*"];'
    case (4)
        write(File,'(I0.0,A)') Node_ID, '[label="/"];'
    case (5)
        write(File,'(I0.0,A)') Node_ID, '[label="IGF"];'
    case (6)
        write(File,'(I0.0,A)') Node_ID, '[label="MMT"];'
    case (7)
        write(File,'(I0.0,A)') Node_ID, '[label="MPGF"];'
    case (8)
        write(File,'(I0.0,A)') Node_ID, '[label="pow"];'
    case (9)
        write(File,'(I0.0,A)') Node_ID, '[label="min"];'
end select

end subroutine Dot_Graph_Function


subroutine Dot_Graph_Hierarchy ( File, Node_Index )
implicit none

! Input
integer (kind=4), intent(in) :: File, Node_Index

! Locals
integer (kind=4) :: parent_Node

!--------------------------------------------------------------------------------

parent_Node = Node_Index / 2
if( parent_Node .gt. 0) then
    write (File,*) parent_Node, ' -> ', Node_Index, ';'
endif

end subroutine Dot_Graph_Hierarchy
