subroutine print_trees( i_gen, n_indiv_start, n_indiv_stop, &
                        tree_type, tree_descrip )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!use mpi
!use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none


integer,intent(in) :: n_indiv_start
integer,intent(in) :: n_indiv_stop
character(*),intent(in) :: tree_descrip

character(2),dimension( 1:n_nodes ) :: tree_type_string

integer(kind=4), parameter :: string_length  = 1000
integer(kind=4), parameter :: element_length = 5   ! 7

character(1)  ::  element_fmt2     
character(4)  ::  element_format   

character(string_length)   :: node_string
character(string_length)   :: value_string
character(element_length)  :: node_element_string2
character(element_length)  :: value_element_string2

character(7)  :: tree_type_fmt1 
character(2)  :: tree_type_fmt2
character(7)  :: tree_type_fmt3 

character(16)  :: tree_type_fmt
integer(kind=4), intent(in), &
        dimension( 1:n_nodes, 1:n_trees, n_GP_individuals) :: tree_type

integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_Tree
integer(kind=4) :: i_node
integer(kind=4) :: i_gen
integer(kind=4) :: i
integer(kind=4) :: jj

integer(kind=4) :: nodes_filled
integer(kind=4) :: isub1         
integer(kind=4) :: isub2         
integer(kind=4),parameter :: node_boundary = 32
integer(kind=4),parameter :: nodes_on_line = 15
!----------------------------------------------------------------------------------------

write( element_fmt2, '(I1)') element_length

element_format= '(I' // element_fmt2 // ')'

!-----------------------------------------------

! format for small number of nodes

if( n_nodes < node_boundary ) then

    tree_type_fmt1 = '(I6,4x,'
    write(tree_type_fmt2, '(I2)')  n_nodes
    tree_type_fmt3 = '(1x,A))'

    tree_type_fmt = tree_type_fmt1 // &
                    trim( tree_type_fmt2 ) // tree_type_fmt3

    !write(6,'(A,A)') 'pt: tree_type_fmt1 ', tree_type_fmt1
    !write(6,'(A,A)') 'pt: tree_type_fmt2 ', tree_type_fmt2
    !write(6,'(A,A)') 'pt: tree_type_fmt3 ', tree_type_fmt3
    !write(6,'(A,A)') 'pt: tree_type_fmt  ', tree_type_fmt 

endif ! n_nodes < node_boundary 

!-----------------------------------------------

write(6,'(/A,2(1x,I6)/)' ) 'pt: n_nodes ', n_nodes

! print trees

write(GP_print_unit,'(/A)')  &
 'pt: ############################################################################'
if( len( trim(tree_descrip) ) > 0 )write(GP_print_unit,'(A)')  tree_descrip

write(6,'(A,2(1x,I6))' ) &
     'pt: n_indiv_start, n_indiv_stop ', n_indiv_start, n_indiv_stop  ! debug

do  i_GP_individual = n_indiv_start, n_indiv_stop

    if( n_indiv_stop - n_indiv_start > 0 )then

        write(GP_print_unit,'(/A/A,2(1x,I6)/)') &
          '============================================================================', &
          'pt: i_generation, i_GP_indiv ', i_gen, i_GP_individual
    else

        write(GP_print_unit,'(/A,2(1x,I6)/)') &
          'pt: i_generation, i_GP_indiv ', i_gen, i_GP_individual

    endif !  n_indiv_stop - n_indiv_start > 1


    if( n_nodes < node_boundary ) then
        write(GP_print_unit,'(A)') 'pt: i_tree                    nodes '
        write(GP_print_unit,'(10x,A)') trim( tree_node_string )
    endif ! n_nodes < node_boundary 


    write(6,'(/A,2(1x,I6)/)') &
          'pt: begin tree loop n_trees = ', n_trees   ! debug

    do  i_Tree=1,n_Trees

        !write(6,'(I6,4x,20(1x,I2))' ) &
        !     i_tree, Tree_Type( 1:n_nodes, i_tree,i_GP_individual)  ! debug

        if( .not. any( Tree_Type > -9999 )  ) cycle 

        if( n_nodes < node_boundary ) then

            tree_type_string = '  '
            do  i_node = 1, n_nodes
    
                write(6,'(A,3(1x,I6))' ) &
                     'pt: i_tree, i_node, Tree_Type( i_node, i_tree,i_GP_individual ) ', &
                          i_tree, i_node, Tree_Type( i_node, i_tree,i_GP_individual )  ! debug
    
                if( Tree_Type( i_node, i_tree, i_GP_individual) == -9999 )then
                    tree_type_string(i_node)  = ' .'
                else
                    write(tree_type_string(i_node), '(I2)') &
                          tree_type( i_node, i_tree, i_GP_individual )
                endif ! tree_type(...-9999
    
                write(6,'(A,2(1x,I6), 1x, A)' ) &
                     'pt: i_tree, i_node, Tree_Type_string( i_node ) ', &
                          i_tree, i_node, Tree_Type_string( i_node )  ! debug
    
            enddo ! i_node
    
    
            write(6,'(I6,4x,20(1x,I5))' ) &
                 i_tree, ( Tree_Type(jj, i_tree,i_GP_individual), jj = 1, n_nodes )  ! debug
    

            !write(GP_print_unit,'(I6,4x,50(1x,A))' ) &
            write(GP_print_unit, tree_type_fmt  ) &
                 i_tree, Tree_Type_string(1:n_nodes)



        else  ! n_nodes > node_boundary


            write(GP_print_unit,'(/A,1x,I6)' )  'pt: i_tree ', i_tree

            node_string = ''
            value_string = ''

            nodes_filled = 0
            do  i_node = 1, n_nodes

                if( Tree_Type( i_node, i_tree, i_GP_individual) > -9999 ) then

                    !write(node_element_string2, '(I7)') i_node
                    !write(value_element_string2, '(I7)') &

                    write( node_element_string2, element_format ) i_node
                    write( value_element_string2, element_format ) &
                                  Tree_Type( i_node, i_tree, i_GP_individual)

                    node_string  = trim(node_string) // node_element_string2
                    value_string = trim(value_string) // value_element_string2
                    nodes_filled =  nodes_filled + 1

                endif ! Tree_Type( i_node, i_tree, i_GP_individual) > -9999

            enddo ! i_node 

            write(6,'(/A, 2(1x,I6)/)')  'pt; i_tree, nodes_filled ', i_tree, nodes_filled  ! debug

            if( nodes_filled > 0 ) then 
                do  i = 1, nodes_filled/nodes_on_line + 1
                    !write(GP_print_unit,'(A, 1x,I6 )')  'pt; i ', i

                    !isub1 = 1 + (i-1) * 7 * nodes_on_line
                    !isub2 =         i * 7 * nodes_on_line
                    isub1 = 1 + (i-1) * element_length * nodes_on_line
                    isub2 =         i * element_length * nodes_on_line

                    if( nodes_filled - (i-1)*nodes_on_line > 0 )then
                        write(GP_print_unit,'(A, A      )') &
                                  'node: ', trim( node_string( isub1: isub2    ) )
                        write(GP_print_unit,'(A, A     /)') &
                                  'value:', trim( value_string( isub1: isub2   ) )
                    endif ! nodes_filled - (i-1)*nodes_on_line > 0 
                enddo 
            endif ! nodes_filled > 0 


        endif ! n_nodes < node_boundary

        !!write(GP_print_unit,'(A)') ' '

    enddo ! i_tree


    !---------------------------------------------------------------

    ! print equations corresponding to the tree

    if( L_print_equations )then

        write(6,'(/A/)' ) &
                 'pt: call create_equations '  ! debug

        call create_equations( i_gen, i_GP_individual, tree_type )
    endif ! L_print_equations

    !---------------------------------------------------------------


enddo  ! i_GP_individual

write(GP_print_unit,'(A/)')  &
      'pt: ############################################################################'


return

end subroutine print_trees
