! Generate_Dot_Graph
! Written by Erik Wisuri [14 June, 2013]
! These subroutines generate DOT output from a tree array,
! and use the dot executable to create PDF files.

module class_Dot_Graph_Visitor

    use class_Tree_Node

    type, public, extends(Tree_Node_Visitor) :: Dot_Graph_Visitor
        integer(kind=4) :: funit, node_id
        contains
        procedure :: Visit_Tree_Node => Dot_Visit_Tree_Node
        procedure :: Visit_Tree_Math_Node => Dot_Visit_Math_Node
        procedure :: Visit_Tree_Parameter_Node => Dot_Visit_Parameter_Node
        procedure :: Visit_Tree_Variable_Node => Dot_Visit_Variable_Node
    end type Dot_Graph_Visitor

contains


    !-------------------------------------------------------------------------------------


    subroutine Dot_Visit_Tree_Node(this, node)
        class(Dot_Graph_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node

        write(6,'(//A//)') &
             'Generate_Dot_Graph: Error: generic type Tree_Node encountered in tree traversal.'

        !!stop 1 ! Stop program
        return 

    end subroutine Dot_Visit_Tree_Node


    !-------------------------------------------------------------------------------------


    subroutine Dot_Visit_Math_Node(this, node)

        class(Dot_Graph_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node
        integer(kind=4) :: myid

        myid = this%node_id

        write(this%funit,'(I0.0,A,I0.0,A)',advance='no') &
                               this%node_id, '[label="[', this%node_id, '] '
        call Dot_Graph_Function(this%funit, node%operation)

        write(this%funit,'(A)') '"];'
        call Dot_Graph_Hierarchy(this%funit, myid)

        this%node_id = myid*2
        call node%left%accept(this)

        this%node_id = myid*2 + 1
        call node%right%accept(this)

    end subroutine Dot_Visit_Math_Node


    !-------------------------------------------------------------------------------------


    subroutine Dot_Visit_Parameter_Node(this, node)
        class(Dot_Graph_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node

        write(this%funit,'(I0.0,A,I0.0,A)',advance='no') &
                            this%node_id, '[label="[', this%node_id, '] '
        write(this%funit,'(A,E12.5)',advance='no') '(P) ', node%val()
        write(this%funit,'(A)') '"];'

        call Dot_Graph_Hierarchy(this%funit, this%node_id)

    end subroutine Dot_Visit_Parameter_Node


    !-------------------------------------------------------------------------------------


    subroutine Dot_Visit_Variable_Node(this, node)
        class(Dot_Graph_Visitor), intent(inout) :: this
        class(Tree_Node), intent(in) :: node

        write(this%funit,'(I0.0,A,I0.0,A)',advance='no') &
                         this%node_id, '[label="[', this%node_id, '] '
        write(this%funit,'(A,E12.5)',advance='no') '(V) ', node%val()
        write(this%funit,'(A)') '"];'

        call Dot_Graph_Hierarchy(this%funit, this%node_id)

    end subroutine Dot_Visit_Variable_Node


end module class_Dot_Graph_Visitor


!---------------------------------------------------------------------------------------------------


subroutine Generate_Dot_Graph( Trees, Tree_count, output_dir )
    use class_Tree_Node
    use class_Dot_Graph_Visitor
    implicit none

    ! Input
    character(len=*), intent(in) :: output_dir
    integer(kind=4), intent(in) :: Tree_count
    type(Tree_Node_Pointer), dimension(Tree_count), intent(in) :: Trees ! The array of trees
    type(Dot_Graph_Visitor) :: grapher

    ! Local variables
    integer(kind=4) :: i, gFile
    character(len=80) :: Graph_File

!------------------------------------------------------------------------------------------------

    gFile = 85;

    do  i = 1,Tree_count
        if( associated(Trees(i)%n) ) then

            write(Graph_File, '(A,I0.0)') output_dir//'/Trees/', i

            open(gFile, FILE=trim(Graph_File)//'.dot')

            write(gFile,*) 'digraph g {'
            write(gFile,*) 'splines=false;'
            grapher = Dot_Graph_Visitor(gFile, 1)

            call Trees(i)%n%accept(grapher)

            write(gFile,*) '}'

            close(gFile)

            call system('dot '//trim(Graph_File)// &
                        '.dot -T pdf -Nheight=0.5 -Nwidth=0.02 -o'//trim(Graph_File)//'.pdf')
        endif ! associated

    enddo ! i 

end subroutine Generate_Dot_Graph



subroutine Dot_Graph_Function( File, Function_Index)
    implicit none

    ! Input
    integer(kind=4), intent(in) :: File, Function_Index

    select case (Function_Index)
        case (1)
            write(File,'(A)',advance='no') '+'
        case (2)
            write(File,'(A)',advance='no') '-'
        case (3)
            write(File,'(A)',advance='no') '*'
        case (4)
            write(File,'(A)',advance='no') '/'
        case (5)
            write(File,'(A)',advance='no') 'IGF'
        case (6)
            write(File,'(A)',advance='no') 'MMT'
        case (7)
            write(File,'(A)',advance='no') 'MPGF'
        case (8)
            write(File,'(A)',advance='no') 'pow'
        case (9)
            write(File,'(A)',advance='no') 'exp'
            !orig write(File,'(A)',advance='no') 'min'
        case (10)
            write(File,'(A)',advance='no') 'min'
            !orig write(File,'(A)',advance='no') 'max'
        case (11)
            write(File,'(A)',advance='no') 'max'
            !orig write(File,'(A)',advance='no') 'exp'
        case (12)
            write(File,'(A)',advance='no') 'if'
        case (13)
            write(File,'(A)',advance='no') '>'
        case (14)
            write(File,'(A)',advance='no') '>='
        case (15)
            write(File,'(A)',advance='no') '<'
        case (16)
            write(File,'(A)',advance='no') '<='
        case (17)
            write(File,'(A)',advance='no') 'expLP'
        case (18)
            write(File,'(A)',advance='no') 'expRP'
        case (19)
            write(File,'(A)',advance='no') 'expLM'
        case (20)
            write(File,'(A)',advance='no') 'expRM'
    end select

end subroutine Dot_Graph_Function



subroutine Dot_Graph_Hierarchy( File, Node_Index )
    implicit none

    ! Input
    integer(kind=4), intent(in) :: File, Node_Index

    ! Locals
    integer(kind=4) :: parent_Node

    parent_Node = Node_Index / 2
    if( parent_Node .gt. 0) then
        write(File,*) parent_Node, ' -> ', Node_Index, ';'
    endif

end subroutine Dot_Graph_Hierarchy
