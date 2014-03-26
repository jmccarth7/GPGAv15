!
! File:   Tree_Node_Factory.f03
! Author: Dave
!
! Created on August 5, 2013, 11:18 AM
!

module Tree_Node_Factory_module

    use mpi
    use mpi_module

    use class_Tree_Node

    integer(kind=4), parameter :: Add = 1
    integer(kind=4), parameter :: Subtract = 2
    integer(kind=4), parameter :: Multiply = 3
    integer(kind=4), parameter :: ProtectedDivide = 4
    integer(kind=4), parameter :: IvlevGrazingFunction = 5
    integer(kind=4), parameter :: MichealisMenton = 6
    integer(kind=4), parameter :: MayzaudPouletGrazingFunction = 7
    integer(kind=4), parameter :: Power = 8
    integer(kind=4), parameter :: ExponentialDecay = 9
    integer(kind=4), parameter :: Minimize = 10
    integer(kind=4), parameter :: Maximize = 11
!orig    integer(kind=4), parameter :: Minimize = 9
!orig    integer(kind=4), parameter :: Maximize = 10
!orig    integer(kind=4), parameter :: ExponentialDecay = 11
    integer(kind=4), parameter :: IfThen = 12
    integer(kind=4), parameter :: IfGt = 13
    integer(kind=4), parameter :: IfGte = 14
    integer(kind=4), parameter :: IfLt = 15
    integer(kind=4), parameter :: IfLte = 16
    integer(kind=4), parameter :: ExponentialLeftPlus  = 17
    integer(kind=4), parameter :: ExponentialRightPlus = 18
    integer(kind=4), parameter :: ExponentialLeftMinus  = 19
    integer(kind=4), parameter :: ExponentialRightMinus = 20

    integer(kind=4), parameter :: MathNodeType = 1
    integer(kind=4), parameter :: VariableNodeType = 2
    integer(kind=4), parameter :: ParameterNodeType = 3

    contains

    !--------------------------------------------------------------------------------

    recursive function GetMathNode( OperationIndex, Left, Right ) result(node)

        implicit none
        type(Tree_Node), pointer :: node
        type(Tree_Node), pointer :: Left, Right
        integer(kind=4) :: OperationIndex

        !--------------------------------------------------------------------

        ! Sanity Check
        if( .not. associated(Left) .or. &
            .not. associated(Right)       ) then

            write(*,'(//A//)') &
                       'GetMathNode:  MathNode with NULL child(ren)'
            call MPI_FINALIZE(ierr)
            stop 1 ! Stop program

        endif

        nullify(node)

        allocate(node)

        write(*,'(A)')       'GMN: entry GetMathNode '
        write(*,'(A,1x,I6)') 'GMN: OperationIndex ', OperationIndex

        ! Constructor: Node Type, Node Count, Parent,
        !              LeftChild, RightChild,
        !              Operation, Variable_Index, VariableValue, ParameterValue

        node = Tree_Node( MathNodeType, 1 + Left%node_count + Right%node_count, NULL(), &
                          Left, Right,    &
                          OperationIndex, 0, Null(), 0.D+0, &
                          Tree_Math_Node_Val,               &
                          Tree_Math_Node_Delete,            &
                          Tree_Math_Node_Accept_Visitor,    &
                          Tree_Math_Node_Randomize,         &
                          Tree_Math_Node_Get_Pointers,      &
                          Tree_Node_Swap  )

        Left%parent => node
        Right%parent => node

        !write(*,'(A)') 'GMN: leave GetMathNode '

    end function


    !--------------------------------------------------------------------------------

    recursive function GetVariableNode(VariableValue, Variable_Index) result (node)

        implicit none
        type(Tree_Node), pointer :: node
        real(kind=8), target:: VariableValue
        integer(kind=4) :: Variable_Index

        !--------------------------------------------------------------------------

        write(*,'(A)')           'GVN: entry GetVariableNode'
        write(*,'(A,1x,E24.16)') 'GVN: VariableValue  ', VariableValue
        write(*,'(A,1x,I6   )')  'GVN: Variable_Index ', Variable_Index

        allocate(node)

        ! Constructor: Node Type, Node Count, Parent, LeftChild, RightChild,
        !              Operation, Variable_Index, VariableValue, ParameterValue

        node = Tree_Node( VariableNodeType, 1, Null(), Null(), Null(), 0, &
                          Variable_Index, VariableValue, 0.D+0, &
                          Tree_Variable_Node_Val,  &
                          Tree_Node_Delete,  &
                          Tree_Variable_Node_Accept_Visitor,  &
                          Tree_Variable_Node_Randomize,  &
                          Tree_Node_Get_Pointers,  &
                          Tree_Node_Swap)

        write(*,'(A)') 'GVN: leave GetVariableNode'

    end function


    !--------------------------------------------------------------------------------

    recursive function GetParameterNode(ParameterValue) result (node)

        implicit none
        type(Tree_Node), pointer :: node
        real(kind=8) :: ParameterValue

        !----------------------------------------------------------------

        allocate(node)

        write(*,'(A)')           'GPN:  enter GetParameterNode '
        write(*,'(A,1x,E24.16)') 'GPN:  ParameterValue ', ParameterValue

        ! Constructor: Node Type, Node Count, Parent,
        !              LeftChild, RightChild, Operation,
        !              Variable_Index, VariableValue, ParameterValue

        node = Tree_Node( ParameterNodeType, 1, Null(), &
                          Null(), Null(), 0, 0, Null(), ParameterValue, &
                          Tree_Parameter_Node_Val, &
                          Tree_Node_Delete, &
                          Tree_Parameter_Node_Accept_Visitor, &
                          Tree_Parameter_Node_Randomize, &
                          Tree_Node_Get_Pointers, &
                          Tree_Node_Swap)

        write(*,'(A)') 'GPN:  leave GetParameterNode '

    end function


    !--------------------------------------------------------------------------------

    function GetNodeFromRead( node_type, node_operation, &
                              variable_index, parameter_value)     result(node)
        implicit none
        type(Tree_Node), pointer :: node
        integer(kind=4) :: node_type, node_operation, variable_index
        real(kind=8) :: parameter_value

        allocate(node)

        select case( node_type )

            case( MathNodeType )

               ! Left and Right are NULL, Count is 1

               node = Tree_Node( MathNodeType, 1, NULL(), NULL(), NULL(), &
                                 node_operation, 0, Null(), 0.D+0, &
                                 Tree_Math_Node_Val, &
                                 Tree_Math_Node_Delete, &
                                 Tree_Math_Node_Accept_Visitor, &
                                 Tree_Math_Node_Randomize, &
                                 Tree_Math_Node_Get_Pointers, &
                                 Tree_Node_Swap  )

            case( VariableNodeType )

               ! Variable Value is NULL

               node = Tree_Node( VariableNodeType, 1, Null(), Null(), Null(), &
                                 0, variable_index, NULL(), 0.D+0, &
                                 Tree_Variable_Node_Val, &
                                 Tree_Node_Delete, &
                                 Tree_Variable_Node_Accept_Visitor, &
                                 Tree_Variable_Node_Randomize, &
                                 Tree_Node_Get_Pointers, &
                                 Tree_Node_Swap )

            case( ParameterNodeType )

               node = Tree_Node( ParameterNodeType, 1, Null(), Null(), Null(), &
                                 0, 0, Null(), parameter_value, &
                                 Tree_Parameter_Node_Val, &
                                 Tree_Node_Delete, &
                                 Tree_Parameter_Node_Accept_Visitor, &
                                 Tree_Parameter_Node_Randomize, &
                                 Tree_Node_Get_Pointers, &
                                 Tree_Node_Swap)

            case default

               deallocate(node)

               node => NULL()

               write(*,'(//A,1x,I10//)') &
                  'GetNodeFromRead:  Bad Node Type:  node_type = ', node_type
               call MPI_FINALIZE(ierr)
               stop 1 ! Stop program

        end select

    end function


    !--------------------------------------------------------------------------------

    recursive function GetRandomNode() result(node)

        implicit none
        type(Tree_Node), pointer :: node
        real(kind = 8) :: rrnd

        call random_number(rrnd)

!        rrnd = rrnd * 3

        if( rrnd .ge. 0.66D+0 ) then

            node => GetMathNode( 1, GetRandomNode(), GetRandomNode() )

        elseif( rrnd .ge. 0.33D+0 ) then

            node => GetParameterNode(0.D+0)
        else

            ! -1 is a dummy variable standing in for the variable index

            node => GetVariableNode( 0.D+0, -1 )

        endif

        call node%Randomize()

    end function

    !--------------------------------------------------------------------------------

end module Tree_Node_Factory_module
