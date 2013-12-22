module Tree_Helper_module

use class_Tree_Node


contains

!===========================================================================

recursive function GetNodeCount(node) result(count)

implicit none
integer(kind=4) :: count

type(Tree_Node), pointer :: node

!--------------------------------------------

!write(6,*)'gnc:0 node%node_type, count ', node%node_type, count

if( node%node_type .eq. 1 ) then


    node%node_count = 1 + GetNodeCount(node%left) + &
                          GetNodeCount(node%right)

    !write(6,*)'gnc:1 GetNodeCount(node%left) ', GetNodeCount(node%left)
    !write(6,*)'gnc:1 GetNodeCount(node%right)', GetNodeCount(node%right)
    !write(6,*)'gnc:1 node%node_count, count ', node%node_count, count

    count = node%node_count

    !write(6,*)'gnc:2 node%node_count, count ', node%node_count, count
endif

count = node%node_count

!write(6,*)'gnc:3 node%node_type, count ', node%node_type, count

end function

!===========================================================================

function GetMaxHeight(Trees, Tree_count)  result(maxHeight)

implicit none
integer(kind=4) :: Tree_count, currentHeight, maxHeight, i
type(Tree_Node_Pointer), dimension(Tree_count), intent(in) :: Trees ! The array of trees

!--------------------------------------------

maxHeight = 0
currentHeight = 0

do  i = 1,Tree_count

    currentHeight = GetTreeHeight(Trees(i)%n)
    if( currentHeight .gt. maxHeight) then
        maxHeight = currentHeight
    endif

enddo

end function

!===========================================================================

recursive function GetTreeHeight(node) result(height)

implicit none
type(Tree_Node), pointer :: node
integer(kind=4) :: height

!--------------------------------------------

height = 0
! Sanity Check
if( associated(node)) then
    height = max(GetTreeHeight(node%left), GetTreeHeight(node%right)) + 1
endif

end function

!===========================================================================

end module  Tree_Helper_module
