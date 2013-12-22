subroutine Deserialize_Trees2( Trees, num_Tracked_resources, Tree_count )

use mpi
use mpi_module

use GP_variables_module

use class_Serialization_Visitor
use class_Tree_Node
use Tree_Helper_module

use Tree_Node_Factory_module

implicit none

! Input
integer (kind=4):: Tree_count, num_Tracked_resources

! Input/Output
type(Tree_Node_Pointer), &
     dimension(Tree_count, num_Tracked_resources), intent(inout) :: Trees

! Local variables

integer (kind=4):: node_count, left, right, node_type
integer (kind=4):: node_id, node_operation, variable_index
integer (kind=4):: i, j, k, l
integer (kind=4):: inode
integer (kind=4):: inodec
integer (kind=4):: counter
integer (kind=4):: n_inode

real(kind=8) :: parameter_value

type(Tree_Node_Pointer), dimension(:), allocatable :: Nodes
integer(kind=4), dimension(:), allocatable :: Node_IDs

character (len=80) :: file_name
logical :: file_exists
type(Tree_Node), pointer :: parent, root

logical,save :: first  = .TRUE.
integer (kind=4):: temp_myid

!------------------------------------------------------------------------------------

if( first ) then        ! debug_only
    temp_myid = 0
    first = .false.
else
    temp_myid = 1
endif

! MathNodeType      = 1
! VariableNodeType  = 2
! ParameterNodeType = 3


!if( myid == temp_myid )then
!
!     write(60,'(//A,2(1x,I6))') 'DsT2: at entry n_nodes, tree_count ', n_nodes, tree_count
!     write(60,'(A,2(1x,I6))')   'DsT2: num_Tracked_resources   ', num_Tracked_resources
!
!     write(6,'(//A,2(1x,I6))') 'DsT2: at entry n_nodes, tree_count ', n_nodes, tree_count
!     write(6,'(A,2(1x,I6))')   'DsT2: num_Tracked_resources   ', num_Tracked_resources
!
!    call print_debug_integer_node_tree( 6,  'GP_Individual_Node_Type in Deser_tree2 ', &
!                                          GP_Individual_Node_Type )
!    call print_debug_integer_node_tree( 60,  &
!                                       'GP_Individual_Node_Type in Deser_tree2 ', &
!                                        GP_Individual_Node_Type )
!
!    call print_debug_real_node_tree( 6, 'GP_Individual_Node_parameters in Deser_tree2 ', &
!                                      GP_Individual_Node_parameters  )
!    call print_debug_real_node_tree( 60, &
!                                    'GP_Individual_Node_parameters in Deser_tree2 ', &
!                                     GP_Individual_Node_parameters  )
!endif !  myid == temp_myid



do i = 1, Tree_count

    !if( myid == temp_myid )then
    !    write(60,'(//A,1x,I6,1x,A/)') 'DsT2: Tree  i = ', i, &
    !       '#########################################################################'
    !    write(6,'(//A,1x,I6,1x,A/)')  'DsT2: Tree  i = ', i, &
    !       '#########################################################################'
    !endif !  myid == temp_myid


    do j = 1,  num_Tracked_resources

        Trees(i, j)%n => null()


        counter = 0
        do  inodec = 1, n_nodes

            if( GP_Individual_Node_Type( inodec, i ) > -9999 )then
                counter =  counter + 1
            endif ! GP_Individual_Node_Type...

        enddo ! inodec

        node_count = counter

        !if( myid == temp_myid )then
        !    write(60,'(//A,3(1x,I6))') 'DsT2: tree i, node_count, counter ', &
        !                                           i, node_count, counter
        !    write(6,'(//A,3(1x,I6))') 'DsT2: tree i, node_count, counter ', &
        !                                          i, node_count, counter
        !endif !  myid == temp_myid


        if( node_count <= 0 )  exit  ! j loop



        ! Dimension arrays that will hold nodes and node ids

        !if( myid == temp_myid )then
        !    write(6,'(A,2(1x,I6))')  'DsT2: allocate Nodes, Node_IDS for tree = ', i
        !    write(60,'(A,2(1x,I6))') 'DsT2: allocate Nodes, Node_IDS for tree = ', i
        !endif !  myid == temp_myid

        allocate( Nodes(node_count), Node_IDs(node_count) )



        counter = 0
        do  inode = 1,  n_nodes

            node_type        = 0
            node_id          = inode
            node_operation   = 0
            variable_index   = 0
            parameter_value  = 0.0D0


            if( GP_Individual_Node_Type( inode, i ) > -9999 )then

                !if( myid == temp_myid )then
                !    write(60,'(//A,3(1x,I6))') &
                !      'DsT2: tree i, inode, GP_Individual_Node_Type( inode, i )', &
                !                  i, inode, GP_Individual_Node_Type( inode, i )
                !    write(6,'(//A,3(1x,I6))') &
                !      'DsT2: tree i, inode, GP_Individual_Node_Type( inode, i )', &
                !                  i, inode, GP_Individual_Node_Type( inode, i )
                !endif !  myid == temp_myid

                counter =  counter + 1

                !!node_type =  GP_Individual_Node_Type( inode, i )
                !node_id   =  counter

                !write(60,'(A,3(1x,I6))') 'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                !                                i, inode, GP_Individual_Node_Type( inode, i )
                !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                !                               i, inode, GP_Individual_Node_Type( inode, i )



                if( GP_Individual_Node_Type( inode, i ) == 0 )then

                    parameter_value =  GP_Individual_Node_parameters( inode, i )
                    node_type = ParameterNodeType

                    !write(60,'(A,2(1x,I6),1x,E15.7)') &
                    !      'DsT2: i, inode, GP_Individual_Node_parameters( inode, i )', &
                    !             i, inode, GP_Individual_Node_parameters( inode, i )
                    !write(60,'(A,2(1x,I6),1x,E15.7)') 'DsT2: i, inode,  parameter_value', &
                    !                                         i, inode,  parameter_value

                    !write(6,'(A,2(1x,I6),1x,E15.7)') &
                    !      'DsT2: i, inode, GP_Individual_Node_parameters( inode, i )', &
                    !             i, inode, GP_Individual_Node_parameters( inode, i )
                    !write(6,'(A,2(1x,I6),1x,E15.7)') 'DsT2: i, inode,  parameter_value', &
                    !                                        i, inode,  parameter_value

                endif


                if( GP_Individual_Node_Type( inode, i ) <  0 )then

                    variable_index =  GP_Individual_Node_type( inode, i )
                    node_type = VariableNodeType

                    !write(60,'(A,3(1x,I6))') 'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                    !                                i, inode, GP_Individual_Node_Type( inode, i )
                    !write(60,'(A,4(1x,I6))') 'DsT2: i, inode, node_type, variable_index ', &
                    !                                i, inode, node_type, variable_index

                    !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                    !                               i, inode, GP_Individual_Node_Type( inode, i )
                    !write(6,'(A,4(1x,I6))') 'DsT2: i, inode, node_type, variable_index ', &
                    !                               i, inode, node_type, variable_index

                endif

                if( GP_Individual_Node_Type( inode, i ) >  0 )then

                    node_operation =  GP_Individual_Node_type( inode, i )
                    node_type = MathNodeType

                    !write(60,'(A,3(1x,I6))') 'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                    !                                i, inode, GP_Individual_Node_Type( inode, i )
                    !write(60,'(A,4(1x,I6))') 'DsT2: i, inode, node_type, node_operation ', &
                    !                                i, inode, node_type, node_operation

                    !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, GP_Individual_Node_Type( inode, i )', &
                    !                               i, inode, GP_Individual_Node_Type( inode, i )
                    !write(6,'(A,4(1x,I6))') 'DsT2: i, inode, node_type, node_operation ', &
                    !                               i, inode, node_type, node_operation

                endif



                !write(60,'(/A,4(1x,I6))') 'DsT2: i, inode, node_id, node_type', &
                !                                 i, inode, node_id, node_type
                !write(60,'(A,2(1x,I6),1x,E15.7)') &
                !                         'DsT2: i, inode, parameter_value   ', &
                !                                i, inode, parameter_value
                !write(60,'(A,3(1x,I6))') 'DsT2: i, inode, variable_index     ', &
                !                                i, inode, variable_index
                !write(60,'(A,3(1x,I6))') 'DsT2: i, inode, node_operation     ', &
                !                                i, inode, node_operation

                !write(6,'(/A,4(1x,I6))') 'DsT2: i, inode, node_id, node_type', &
                !                                i, inode, node_id, node_type
                !write(6,'(A,2(1x,I6),1x,E15.7)') &
                !                         'DsT2: i, inode, parameter_value   ', &
                !                                i, inode, parameter_value
                !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, variable_index     ', &
                !                               i, inode, variable_index
                !write(6,'(A,3(1x,I6))') 'DsT2: i, inode, node_operation     ', &
                !                               i, inode, node_operation


                Nodes(counter)%n => GetNodeFromRead( node_type, node_operation, &
                                                    variable_index, parameter_value)
                Node_IDs(counter) = node_id


                !if( myid == temp_myid )then
                !    write(60,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%variable_index', &
                !                 i, inode, counter, Nodes(counter)%n%variable_index
                !    write(60,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%operation     ', &
                !                 i, inode, counter, Nodes(counter)%n%operation
                !    write(60,'(A,3(1x,I6),1x,E15.7)') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%param         ', &
                !                 i, inode, counter, Nodes(counter)%n%param
                !    write(60,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Node_IDs(counter)              ', &
                !                 i, inode, counter, Node_IDs(counter)
                !    write(6,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%variable_index', &
                !                 i, inode, counter, Nodes(counter)%n%variable_index
                !    write(6,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%operation     ', &
                !                 i, inode, counter, Nodes(counter)%n%operation
                !    write(6,'(A,3(1x,I6),1x,E15.7)') &
                !          'DsT2: i, inode, counter, Nodes(counter)%n%param         ', &
                !                 i, inode, counter, Nodes(counter)%n%param
                !    write(6,'(A,4(1x,I6))') &
                !          'DsT2: i, inode, counter, Node_IDs(counter)              ', &
                !                 i, inode, counter, Node_IDs(counter)
                !endif !  myid == temp_myid

            endif ! GP_Individual_Node_Type...

        enddo ! inode



        ! First node is always the root
        root => Nodes(1)%n

        !if( myid == temp_myid )then
        !    write(60,'(//A,2(1x,I8))') 'DsT2: tree i, root%node_type', i, root%node_type
        !    write(60,'(A,2(1x,I8))')   'DsT2: tree i, node_count    ', i, node_count
        !    write(6,'(//A,2(1x,I8))')  'DsT2: tree i, root%node_type', i, root%node_type
        !    write(6,'(A,2(1x,I8))')    'DsT2: tree i, node_count    ', i, node_count
        !endif !  myid == temp_myid



        ! At this point, we have a collection of nodes, but there are:

        !   1) No associations between Math nodes and their children
        !   2) No associations between Variable nodes, and
        !      the array positions where their value will be calculated

        ! Make the above two associations

        !write(60,'(//A/)')  'DsT2: k loop  '
        !write(6, '(//A/)')  'DsT2: k loop  '

        do  k = 1, node_count

            !write(60,'(A,3(1x,I6))') 'DsT2: i, k, Nodes(k)%n%Node_Type ', &
            !                                i, k, Nodes(k)%n%Node_Type
            !write(6,'(A,3(1x,I6))')  'DsT2: i, k, Nodes(k)%n%Node_Type ', &
            !                                i, k, Nodes(k)%n%Node_Type


            if( Nodes(k)%n%Node_Type .eq. MathNodeType ) then

                ! Grab the parent & calculate children indexes

                parent => Nodes(k)%n
                left = Node_IDs(k)*2
                right = Node_IDs(k)*2+1

                !write(60,'(/A,3(1x,I6))') 'DsT2: i, k, MathNodeType     ', &
                !                                 i, k, MathNodeType
                !write(60,'(A,3(1x,I6))')  'DsT2: i, k, parent node type ', &
                !                                 i, k, parent%node_type
                !write(60,'(A,4(1x,I6))')  'DsT2: i, k, left, right      ', &
                !                                 i, k, left, right

                !write(6,'(/A,3(1x,I6))') 'DsT2: i, k, MathNodeType     ', &
                !                                i, k, MathNodeType
                !write(6,'(A,3(1x,I6))')  'DsT2: i, k, parent node type ', &
                !                                i, k, parent%node_type
                !write(6,'(A,4(1x,I6))')  'DsT2: i, k, left, right      ', &
                !                                i, k, left, right


                ! Grab the children and associate

                do  l = 1,node_count

                    if( Node_IDs(l) .eq. left ) then

                        !write(60,'(/A,3(1x,I6))') 'DsT2:left i, l, Node_IDs(l)', &
                        !                                     i, l, Node_IDs(l)
                        !write(6,'(/A,3(1x,I6))') 'DsT2:left i, l, Node_IDs(l)', &
                        !                                    i, l, Node_IDs(l)

                        parent%left => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    endif
                enddo

                do  l = 1, node_count

                    if( Node_IDs(l) .eq. right ) then

                        !write(60,'(/A,3(1x,I6))') 'DsT2:right i, l, Node_IDs(l)  ', &
                        !                                      i, l, Node_IDs(l)

                        !write(6,'(/A,3(1x,I6))') 'DsT2:right i, l, Node_IDs(l)  ', &
                        !                                     i, l, Node_IDs(l)

                        parent%right => Nodes(l)%n
                        Nodes(l)%n%parent => parent
                        exit

                    endif

                enddo


            elseif( Nodes(k)%n%Node_Type .eq. VariableNodeType) then

                ! If the index is in the -5000 range,
                ! this is a forcing function variable.
                ! Associate it with the correct array

                !write(60,'(/A,3(1x,I6))') &
                !      'DsT2: i, k, VariableNodeType         ', i, k, VariableNodeType
                !write(60,'(A,3(1x,I6))')  &
                !      'DsT2: i, k, Nodes(k)%n%variable_index', &
                !             i, k, Nodes(k)%n%variable_index

                !write(6,'(/A,3(1x,I6))') &
                !      'DsT2: i, k, VariableNodeType         ', i, k, VariableNodeType
                !write(6,'(A,3(1x,I6))')  &
                !      'DsT2: i, k, Nodes(k)%n%variable_index', &
                !             i, k, Nodes(k)%n%variable_index


                if( Nodes(k)%n%variable_index < -5000) then

                    Nodes(k)%n%variable => &
                     Numerical_CODE_Forcing_Functions(  abs( 5000 + Nodes(k)%n%variable_index )  )

                    !write(60,'(A,1x,E15.7)') &
                    !   'DsT2: Numerical_CODE_Forcing_Functions(abs(5000 + Nodes(k)%n%variable_index)) ',&
                    !          Numerical_CODE_Forcing_Functions(abs(5000 + Nodes(k)%n%variable_index))
                    !write(6,'(A,1x,E15.7)') &
                    !   'DsT2: Numerical_CODE_Forcing_Functions(abs(5000 + Nodes(k)%n%variable_index)) ',&
                    !          Numerical_CODE_Forcing_Functions(abs(5000 + Nodes(k)%n%variable_index))

                else

                    Nodes(k)%n%variable =>   btmp( abs( Nodes(k)%n%variable_index ) )

                    !write(60,'(A,3(1x,I6))') &
                    !    'DsT2: i, k, Nodes(k)%n%variable_index', &
                    !           i, k, Nodes(k)%n%variable_index
                    !write(60,'(A,2(1x,I6),1x,E15.7)') &
                    !    'DsT2: i, k, btmp(  abs( Nodes(k)%n%variable_index )  ) ', &
                    !           i, k, btmp(  abs( Nodes(k)%n%variable_index )  )

                    !write(6,'(A,3(1x,I6))') &
                    !    'DsT2: i, k, Nodes(k)%n%variable_index', &
                    !           i, k, Nodes(k)%n%variable_index
                    !write(6,'(A,2(1x,I6),1x,E15.7)') &
                    !    'DsT2: i, k, btmp(  abs( Nodes(k)%n%variable_index )  ) ', &
                    !           i, k, btmp(  abs( Nodes(k)%n%variable_index )  )

                endif


            elseif( Nodes(k)%n%Node_Type .eq. ParameterNodeType) then

                !write(60,'(/A,3(1x,I6))') &
                !      'DsT2: i, k, ParameterNodeType    ', i, k, ParameterNodeType
                !write(60,'(A,3(1x,I6))') &
                !      'DsT2: i, k, Nodes(k)%n%Node_Type ', i, k, Nodes(k)%n%Node_Type

                !write(6,'(/A,3(1x,I6))') &
                !      'DsT2: i, k, ParameterNodeType    ', i, k, ParameterNodeType
                !write(6,'(A,3(1x,I6))') &
                !      'DsT2: i, k, Nodes(k)%n%Node_Type ', i, k, Nodes(k)%n%Node_Type

            endif

        enddo ! k


        ! Finally, compute the node count for each node and
        ! assign the root node to the position in the Tree matrix

        root%node_count = GetNodeCount(root)

        Trees(i, j)%n => root

        !if( myid == temp_myid )then
        !    write(60,'(A,2(1x,I6))') 'DsT2: i, root%node_count ', i, root%node_count
        !    write(6,'(A,2(1x,I6))')  'DsT2: i, root%node_count ', i, root%node_count
        !endif !  myid == temp_myid



        ! Clean up
        deallocate( Nodes, Node_IDs )

    enddo ! j

    !write(60,'(A,1x,I6)') 'DsT2: aft J loop i =  ', i
    !write(6,'(A,1x,I6)')  'DsT2: aft J loop i =  ', i

enddo ! i


end subroutine Deserialize_Trees2
