subroutine read_all_summary_file( i_GP_generation,  icall  )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv )
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv  )
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv  )


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use kinds_mod

use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none



integer(kind=i4b) :: i_code_eq
integer(kind=i4b) :: i
integer(kind=i4b) :: istat


integer(kind=i4b),intent(in)  :: i_GP_Generation
!!integer(kind=i4b),intent(in)  :: i_GP_individual
integer(kind=i4b),intent(in)  :: icall       

integer(kind=i4b)             :: i_GP_Gen
integer(kind=i4b)             :: i_GP_indiv

integer(kind=i4b) :: i_Tree
integer(kind=i4b) :: i_Node
!integer(kind=i4b) :: nparm_temp

logical :: Lprint

character(2) :: arrow1
character(2) :: arrow2
character(200) :: Aline 

!----------------------------------------------------------------------------------------

! after GPCODE*opt*  these arrays represent the arrays
! for the  best individual  found in the GA process

!  GP_Population_Initial_Conditions( 1:n_CODE_Equations, i_GP_indiv)
!  GP_Population_Node_Parameters(    1:n_Nodes,1:n_Trees,i_GP_indiv)
!  GP_Adult_Population_Node_Type(    1:n_Nodes,1:n_Trees,i_GP_indiv)

!---------------------------------------------------
! assume this subroutine is called only by cpu 0
!---------------------------------------------------


!--------------------------------------------------------------------------------

open( GP_restart_file_input_unit, file='GP_restart_file', &                                        
      form = 'formatted', access = 'sequential', &                                                  
      status = 'old' )                                 

! set Lprint so printing is done only under the conditions in the if-test

!!Lprint = .TRUE.   ! debug only
Lprint = .FALSE. 

if( i_GP_generation == 1                                  .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0  .or. &
    i_GP_generation == n_GP_generations                          )then
    Lprint = .TRUE.
endif ! i_GP_generation == 1 .or. ...

!--------------------------------------------------------------------------------

readloop:&
do
! write the summary file header for each individual
! which has n_GP_parameters >= n_code_equations


!if( Lprint )then
!
!    write(GP_print_unit, '(/A/7(1x,I10))') &
!      'rasf: i_GP_gen i_GP_indiv   n_code_eq  &
!             &n_trees    n_nodes  n_levels    n_parms', &
!             i_GP_generation, i_GP_indiv, &
!             n_code_equations, n_trees, n_nodes, n_levels, &
!             GP_Individual_N_GP_param(i_GP_indiv)
!
!endif ! Lprint



!  icall is 0 if called from the main program
!  icall is 1 if called from GP_calc_fitness for the best individual

!if( icall == 0 )then

    read(GP_restart_file_input_unit, *, iostat=istat) &
             i_GP_Gen, i_GP_indiv, &
             n_code_equations, n_trees, n_nodes, n_levels
    if( istat /= 0 ) exit readloop

!else
!
!    write(GP_best_summary_output_unit, '(2x,6(1x,I6))') &
!             i_GP_generation, i_GP_indiv, &
!             n_code_equations, n_trees, n_nodes, n_levels
!
!
!endif ! icall == 0 

!--------------------------------------------------------------------------------

! initial conditions

!if( icall == 0 )then

    if( Lprint )then
        write(GP_print_unit,'(/A)')&
          'rasf: i_GP_gen  i_GP_indiv  i_code_eq  &
                &GP_Pop_Init_Cond(i_code_eq, i_GP_Indiv) '
    endif ! Lprint
    
    do  !  i = 1, n_CODE_Equations
    
        read(GP_restart_file_input_unit, '(A)', iostat=istat) Aline

        if( istat /= 0 )then
            exit readloop
        endif ! istat /=0 
        write(6,'(A,1x,A)')'rasf: Aline: ', trim(Aline)

        if( Aline(1:2) == '> '   ) exit

        read(Aline, *)&
              i_GP_Gen, i_GP_indiv, i_code_eq, &
              GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv ) !, &
              !'gen_indiv_eq'
    
        if( Lprint )then
            write(GP_print_unit,'(3(1x,I10), 7x, E24.16)')&
            i_GP_Gen, i_GP_indiv, i_code_eq, &
            GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv )
        endif ! Lprint
    
    enddo  ! i


    !read(GP_restart_file_input_unit, * ) &
    !      arrow1, i_GP_Gen, i_GP_indiv


!else 
!
!
!
!    if( Lprint )then
!        write(GP_print_unit,'(/A)')&
!          'rasf: i_GP_gen  i_GP_indiv  i_code_eq  &
!                &GP_Pop_Init_Cond(i_code_eq, i_GP_Indiv) '
!    endif ! Lprint
!    
!    do  i_code_eq = 1, n_CODE_Equations
!    
!        if( Lprint )then
!            write(GP_print_unit,'(3(1x,I10), 7x, E24.16)')&
!            i_GP_generation, i_GP_indiv, i_code_eq, &
!            GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv )
!        endif ! Lprint
!    
!        write(GP_best_summary_output_unit, '(2x,2(1x,I6),1x,I3, 1x, E24.16,2x,A)')&
!              i_GP_generation, i_GP_indiv, i_code_eq, &
!              GP_Population_Initial_Conditions( i_code_eq, i_GP_indiv ), &
!              'gen_indiv_eq'
!    
!    enddo  ! i_code_eq
!
!
!    write(GP_best_summary_output_unit, '(A,2(1x,I6))') &
!          '> ', i_GP_generation, i_GP_indiv
!
!endif ! icall == 0 

!!--------------------------------------------------------------------------------


! print the node types if node /= -9999


!!if( Lprint )then
!    write(GP_print_unit,'(A)')  &
!          'rasf: i_GP_gen  i_GP_indiv   i_tree     i_node    GP_Indiv_Node_Type'
!!endif ! Lprint

!  write node types to summary file

!if( icall == 0 )then


    !do  i_Tree=1,n_Trees
    !    do  i_Node=1,n_Nodes
    
    do 
            !if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999 ) then
    
        read(GP_restart_file_input_unit, '(A)',iostat=istat ) Aline


        if( istat /= 0 )then
            exit readloop
        endif ! istat /=0 

        write(6,'(A,1x,A)')'rasf: Aline: ', trim(Aline)

        if( Aline(1:2) == '> '   ) exit

           
        read(Aline, * ) &
             i_GP_Gen, i_GP_indiv,i_tree, i_node, &
             GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv)
    

                !if( Lprint )then
                !    write(GP_print_unit,'(2(1x,I10),3(1x,I10))') &
                !          i_GP_generation, i_GP_indiv,i_tree, i_node, &
                !          GP_Individual_Node_Type(i_Node,i_Tree)
                !endif ! Lprint

            !endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999
    
    
        !enddo  ! i_node
    !enddo ! i_tree

    enddo 

    !read(GP_restart_file_input_unit, * ) &
    !      arrow1, i_GP_Gen, i_GP_indiv

!else 
!
!    do  i_Tree=1,n_Trees
!        do  i_Node=1,n_Nodes
!    
!            if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999 ) then
!    
!                write(GP_best_summary_output_unit, '(2x,2(1x,I6),3(1x,I6))') &
!                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
!                      GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv)
!    
!            endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999
!    
!        enddo  ! i_node
!    enddo ! i_tree
!
!    write(GP_best_summary_output_unit, '(A,2(1x,I6))') &
!          '> ', i_GP_generation, i_GP_indiv
!
!endif ! icall == 0 

!!------------------------------------------------------------------------------
!
!
!! print the node parameters (if there are any)
!

!if( Lprint )then
!    write(GP_print_unit,'(/A/)') &
!    'rasf: i_GP_gen i_GP_indiv     tree        node   &
!    &GP_population_node_parameters'
!
!    do  i_tree=1,n_trees
!        do  i_node=1,n_nodes
!
!            if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv ) == 0  ) then
!
!                write(GP_print_unit,'(2(1x,I10), 2(1x,I10),2x, E15.7)') &
!                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
!                      GP_population_node_parameters(i_node,i_tree,i_GP_indiv)
!
!            endif  ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0
!
!        enddo ! i_node
!    enddo  ! i_tree
!
!endif ! Lprint


!------------------------------------------------------------------------------

!if( Lprint )then
!    call print_trees( i_GP_Gen, i_GP_indiv, i_GP_indiv, &
!                      GP_Adult_Population_Node_Type, ' ' )
!endif ! Lprint


!------------------------------------------------------------------------------


! write all non-zero parameters to output file

!if( icall == 0 )then

    !do  i_tree=1,n_trees
    !    do  i_node=1,n_nodes

    do 
            !if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) .ne. -9999 ) then
    
        read(GP_restart_file_input_unit, '(A)',iostat=istat ) Aline


        if( istat /= 0 )then
            exit readloop
        endif ! istat /=0 

        write(6,'(A,1x,A)')'rasf: Aline: ', trim(Aline)

        if( Aline(1:2) == '> ' ) exit
        if( Aline(1:2) == '>>' ) exit readloop
           
    
        !if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0 ) then
    
        read(Aline ,*) &
              i_GP_Gen, i_GP_indiv,i_tree, i_node, &
              GP_population_node_parameters( i_node,i_tree, i_GP_indiv)
    
            !endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0
    
        !enddo ! i_node
    !enddo  ! i_tree

    enddo  

    !write(GP_restart_file_input_unit, '(A,2(1x,I6))') '>>', i_GP_generation, i_GP_indiv
    !read(GP_restart_file_input_unit, * ) arrow2, i_GP_Gen, i_GP_indiv


!else 
!
!    do  i_tree=1,n_trees
!        do  i_node=1,n_nodes
!    
!            if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0        ) then
!    
!                write(GP_best_summary_output_unit,'(2x,2(1x,I6),2(1x,I6), 1x,E24.16)') &
!                      i_GP_generation, i_GP_indiv,i_tree, i_node, &
!                      GP_population_node_parameters( i_node,i_tree, i_GP_indiv)
!    
!            endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_indiv) == 0
!    
!        enddo ! i_node
!    enddo  ! i_tree
!
!    write(GP_best_summary_output_unit, '(A,2(1x,I6))') '>>', i_GP_generation, i_GP_indiv
!
!endif ! icall == 0 

enddo readloop 


!! write for each indiv.  first write in 0*.f90
!!write(GP_restart_file_input_unit, '(4(1x,I10))') n_code_equations, n_trees, n_nodes, n_levels


!---------------------------------------------------------------------------------

close( GP_best_summary_output_unit )  

return

end subroutine read_all_summary_file
