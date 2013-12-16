subroutine Runge_Kutta_Box_Model( )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

use class_Tree_Node
use class_Serialization_Visitor
use Tree_Helper_module
use Tree_Node_Factory_module


use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module



implicit none


!--------------------------------------------------------------------------------------------


real(kind=8),dimension(4) :: Runge_Kutta_Time_Step

data Runge_Kutta_Time_Step /0.0D+0,0.5D+0,0.5D+0,1.0D+0/  ! fraction of R-K time step

!------------------------------------------------------------------------------------------


!Forcing functions are used in computations, so are included here for book keeping purposes



real(kind=8) :: cff

integer(kind=4) :: iter
integer(kind=4) :: i_Time_Step, i_Track, i_Tree
integer(kind=4) :: i_CODE_Equation, j_CODE_Equation, i_Variable
integer(kind=4) :: i_node

integer(kind=4) :: tree_node_count


!--------------------------------------------------------------------------------------

tree_node_count = 0

if( myid == 0 )then
    write(6,'(/A,1x,I6/)') 'rkbm: entry Runge_Kutta_Box_Model myid = ', myid
    write(6,'(A,1x,I6/)')  'rkbm: n_Variables ', n_Variables
endif ! myid == 0

if( L_ga_print )then ! .and. myid == 1 )then
    write(GA_print_unit,'(/A,1x,I6/)') 'rkbm: entry Runge_Kutta_Box_Model myid = ', myid
    write(GA_print_unit,'(A,1x,I6/)')  'rkbm: n_Variables ', n_Variables
    !!!write(GA_print_unit,'(A,1x,I6/)')  'rkbm: tree_value modified'               
endif ! L_ga_print .and. myid == 1


if( dt <= 0.0d0 )then
    if( myid == 0 .or. myid == 1 )then
        write(6,'(/A/)') 'rkbm: BAD VALUE for dt'
        write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt
    endif ! myid == 0 .or. myid == 1
    call MPI_FINALIZE(ierr)
    stop 'bad delta_time'
endif ! dt <= 0.0D0


!if( myid == 0 )then
!    write(6,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop Numerical_CODE_Solution(0,:)', &
!                         Numerical_CODE_Solution(0,:)
!    write(6,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)
!endif ! myid == 0

!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop Numerical_CODE_Solution(0,:)', &
!                         Numerical_CODE_Solution(0,:)
!    write(GA_print_unit,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)
!endif ! L_ga_print .and. myid == 1




! debug >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_type(i_node, i_tree) > -9999 )then
!            write(GA_print_unit,'(A,3(1x,I4))')&
!                  'RuKbm: i_tree, i_node, GP_indiv_node_type', &
!                          i_tree, i_node, GP_individual_node_type(i_node, i_tree)
!        endif !  GP_individual_node_type(i_node, i_tree) > -9999
!    enddo ! i_node
!enddo ! i_tree
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_parameters(i_node, i_tree) > 0.0d0 )then
!            write(GA_print_unit,'(A,2(1x,I4),1x,E15.7)')&
!                  'RuKbm: i_tree, i_node, GP_indiv_node_parms', &
!                          i_tree, i_node, GP_individual_node_parameters(i_node, i_tree)
!        endif !  GP_individual_node_parameters(i_node, i_tree) > 0.0d0
!    enddo ! i_node
!enddo ! i_tree
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_type(i_node, i_tree) > -9999 )then
!            write(6,'(A,3(1x,I4))')&
!                  'RuKbm: i_tree, i_node, GP_indiv_node_type', &
!                          i_tree, i_node, GP_individual_node_type(i_node, i_tree)
!        endif !  GP_individual_node_type(i_node, i_tree) > -9999
!    enddo ! i_node
!enddo ! i_tree
!
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_parameters(i_node, i_tree) > 0.0d0 )then
!            write(6,'(A,2(1x,I4),1x,E15.7)')&
!                  'RuKbm: i_tree, i_node, GP_indiv_node_parms', &
!                          i_tree, i_node, GP_individual_node_parameters(i_node, i_tree)
!        endif !  GP_individual_node_parameters(i_node, i_tree) > 0.0d0
!    enddo ! i_node
!enddo ! i_tree
!
! debug <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



! start the time stepping loop


do  i_Time_Step = 1, n_Time_Steps



    b_tmp(:) = Numerical_CODE_Solution(i_Time_Step-1,:)  ! Array Assignment

    btmp = b_tmp  


    !if( myid == 0 )then
    !    write(6,'(A,1x,I6,10(1x,E15.7) )') &
    !          'rkbm: i_time_step, Num_CODE_Sol(i_Time_Step-1,:) ', &
    !                 i_time_step, Numerical_CODE_Solution(i_Time_Step-1,:)
    !    write(6,'(A,1x,I6,10(1x,E15.7) )') &
    !          'rkbm: i_time_step, b_tmp(:)                      ', &
    !                 i_time_step, b_tmp(:)
    !    write(6,'(A,1x,I6,10(1x,E15.7) )') &
    !          'rkbm: i_time_step, btmp(:)                       ', &
    !                 i_time_step, btmp(:)
    !endif ! myid == 0


    ! carry out a Runge-Kutta time step
    do  iter=1,4

        !if( myid == 0 )then
        !    write(6,'(//A,2(1x,I6)/)') 'rkbm: i_time_step, iter ', &
        !                                      i_time_step, iter
        !    write(GA_print_unit,'(//A,2(1x,I6)/)') 'rkbm: i_time_step, iter ', &
        !                                                  i_time_step, iter
        !endif ! myid == 0

        ! Call forcing functions for the box model

        call DoForcing(btmp, Runge_Kutta_Time_Step(iter), i_Time_Step)

        fbio = 0.0D+0

        do  i_Track = 1,n_Tracked_Resources

            !!!! Call Model_Diagnostics()

            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Evaluate the trees
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            Tree_Value = 0.0D+0                            ! Matrix Assignment

            do  i_Tree=1,n_Trees

                !if( myid == 0 )then
                !    write(6,'(//A,1x,I6,5x,L1)') &
                !          'rkbm: i_tree, associated(GP_Trees(i_Tree,i_Track)%n)  ', &
                !                 i_tree, associated(GP_Trees(i_Tree,i_Track)%n)
                !    write(GA_print_unit,'(//A,1x,I6,5x,L1)') &
                !          'rkbm: i_tree, associated(GP_Trees(i_Tree,i_Track)%n)  ', &
                !                 i_tree, associated(GP_Trees(i_Tree,i_Track)%n)
                !endif ! myid == 0

                if( associated( GP_Trees(i_Tree,i_Track)%n) ) then


                    Tree_Value(i_Tree) = GP_Trees( i_Tree, i_Track )%n%val()

                    !write(6,'(/A,22x,I6,1x,E15.7 )') &
                    !          'rkbm: i_tree, Tree_Value(i_tree) ', &
                    !                 i_tree, Tree_Value(i_tree)
                    write(GA_print_unit,'(/A,22x,I6,1x,E15.7)') &
                              'rkbm: i_tree, Tree_Value(i_tree) ', &
                                     i_tree, Tree_Value(i_tree)

                    !----------------------------------------------------------------------------------------
                    tree_node_count = GetNodeCount( GP_Trees( i_Tree, i_Track )%n )
                    !debug only !if( tree_node_count <= 1 )   Tree_Value(i_Tree) = 0.0d0   ! jjm 20131213
                    !----------------------------------------------------------------------------------------

                    !write(6,'(A,22x,I6,1x,E15.7 )') &
                    !'rkbm: i_tree,Tree_Value(i_tree)',i_tree,Tree_Value(i_tree)
                    write(GA_print_unit,'(A,22x,I6,1x,E15.7)') &
                    'rkbm: i_tree,Tree_Value(i_tree)',i_tree,Tree_Value(i_tree)


                    !if( myid == 0 )then
                    !write(6,'(A,22x,I6,1x,I6//)') &
                    !'rkbm: i_tree,tree_node_count   ',i_tree,tree_node_count
                    write(GA_print_unit,'(A,22x,I6,1x,I6//)') &
                    'rkbm: i_tree,tree_node_count   ',i_tree,tree_node_count
                    !endif ! myid == 0

                    

                endif ! associated(GP_Trees...

            enddo ! i_Trees

            !stop ! debug only

            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Calculate the flow terms from the determined tree_value terms
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            !if( myid == 0 )then
            !    write(6,'(/A/)') 'rkbm: Calculate the flow terms from the determined tree_value terms'
            !endif ! myid == 0

            !write(GA_print_unit,'(/A)') ' '
            !write(6,'(/A)') ' '



            i_Tree=0
            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material
                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .ne. j_CODE_Equation) then

                        i_Tree=i_Tree+1

                        ! 'abs' forces flow of material in one direction

                        bioflo(i_CODE_Equation,j_CODE_Equation)=abs(Tree_Value(i_Tree))


                        !if( L_ga_print )then
                        !    write(GA_print_unit,'(A,3(1x,I6),2(2x, E15.7))') &
                        !    'RuKbm:1 i_tree, i_eq, j_eq, tree_value, bioflo(i_eq,j_eq) ', &
                        !             i_tree, i_code_equation, j_code_equation, &
                        !     tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation)
                        !endif ! L_ga_print

                        !write(6,'(A,3(1x,I6),2(2x, E15.7))') &
                        !    'RuKbm:1 i_tree, i_eq, j_eq, tree_value, bioflo(i_eq,j_eq) ', &
                        !             i_tree, i_code_equation, j_code_equation, &
                        !     tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation)


                    else

                        ! never flow to/from same component
                        bioflo(i_CODE_Equation,j_CODE_Equation)=0.0D+0

                    endif ! i_CODE_Equation .ne. j_CODE_Equation

                enddo ! j_CODE_Equation
            enddo ! i_CODE_Equation

            !write(GA_print_unit,'(/A)') ' '
            !write(6,'(/A)') ' '

            !write(6,'(A,3(1x,I6))') &
            !      'RuKbm: bioflo_map(1:n_eqs, 1) ', &
            !              bioflo_map(1:n_eqs, 1)
            !write(60,'(A,3(1x,I6))') &
            !      'RuKbm: bioflo_map(1:n_eqs, 1) ', &
            !              bioflo_map(1:n_eqs, 1)

            ! bring in the component flow sources and sinks

            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material

                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .gt. 0 ) then

                        !write(6,'(A,3(1x,I6))') &
                        !      'RuKbm: i_eq, bioflo_map(i_eq, 1) ', &
                        !              i_code_equation, bioflo_map(i_code_equation, 1) 
                        !write(60,'(A,3(1x,I6))') &
                        !      'RuKbm: i_eq, bioflo_map(i_eq, 1) ', &
                        !              i_code_equation, bioflo_map(i_code_equation, 1) 

                        if( bioflo_map(i_CODE_Equation,i_Track) .gt. 0 ) then

                            !write(6,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: bef i_eq, fbio(bioflo_map(i_eq, 1)) ', &
                            !                  i_code_equation, fbio(bioflo_map(i_code_equation, 1)) 
                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: bef i_eq, fbio(bioflo_map(i_eq, 1)) ', &
                            !                  i_code_equation, fbio(bioflo_map(i_code_equation, 1)) 

                            fbio(bioflo_map(i_CODE_Equation,i_Track)) = &
                                fbio(bioflo_map(i_CODE_Equation,i_Track)) -  &
                                         bioflo(i_CODE_Equation,j_CODE_Equation)

                            !write(6,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: aft i_eq, fbio(bioflo_map(i_eq,1)),  ', &
                            !                  i_code_equation, fbio(bioflo_map(i_code_equation, 1)) 
                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: aft i_eq, fbio(bioflo_map(i_eq, 1)) ', &
                            !                  i_code_equation, fbio(bioflo_map(i_code_equation, 1)) 
                            !write(6,'(A,2(1x,I6),1x,E20.10)') &
                            !      'RuKbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !              i_code_equation, j_code_equation, &
                            !             bioflo(i_CODE_Equation,j_CODE_Equation) 
                            !write(60,'(A,2(1x,I6),1x,E20.10)') &
                            !      'RuKbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !              i_code_equation, j_code_equation, &
                            !             bioflo(i_CODE_Equation,j_CODE_Equation) 

                        endif ! bioflo_map(i_CODE_Equation,i_Track) .gt. 0

                    endif ! i_CODE_Equation .gt. 0

                    if( j_CODE_Equation .gt. 0 ) then

                        !write(6,'(A,3(1x,I6))') &
                        !      'RuKbm: j_eq, bioflo_map(j_eq, 1) ', &
                        !              j_code_equation, bioflo_map(j_code_equation, 1) 
                        !write(60,'(A,3(1x,I6))') &
                        !      'RuKbm: j_eq, bioflo_map(j_eq, 1) ', &
                        !              j_code_equation, bioflo_map(j_code_equation, 1) 

                        if( bioflo_map(j_CODE_Equation,i_Track) .gt. 0 ) then

                            !write(6,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: bef j_eq, fbio(bioflo_map(j_eq, 1)) ', &
                            !                  j_code_equation, fbio(bioflo_map(j_code_equation, 1)) 
                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: bef j_eq, fbio(bioflo_map(j_eq, 1)) ', &
                            !                  j_code_equation, fbio(bioflo_map(j_code_equation, 1)) 

                            fbio(bioflo_map(j_CODE_Equation,i_Track)) = &
                                 fbio(bioflo_map(j_CODE_Equation,i_Track)) + &
                                          bioflo(i_CODE_Equation,j_CODE_Equation)

                            !write(6,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: aft j_eq, fbio(bioflo_map(j_eq,1)),  ', &
                            !                  j_code_equation, fbio(bioflo_map(j_code_equation, 1)) 
                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'RuKbm: aft j_eq, fbio(bioflo_map(j_eq, 1)) ', &
                            !                  j_code_equation, fbio(bioflo_map(j_code_equation, 1)) 
                            !write(6,'(A,2(1x,I6),1x,E20.10)') &
                            !      'RuKbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !              i_code_equation, j_code_equation, &
                            !             bioflo(i_CODE_Equation,j_CODE_Equation) 
                            !write(60,'(A,2(1x,I6),1x,E20.10)') &
                            !      'RuKbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !              i_code_equation, j_code_equation, &
                            !             bioflo(i_CODE_Equation,j_CODE_Equation) 

                        endif ! bioflo_map(j_CODE_Equation,i_Track) .gt. 0

                    endif ! j_CODE_Equation .gt. 0

                enddo ! j_CODE_Equation
            enddo ! i_CODE_Equation

        enddo ! End Tracked Resources loop


        !write(GA_print_unit,'(/A)') ' '
        !write(6,'(/A)') ' '

        !do  i_CODE_equation=1,n_CODE_equations   ! source of material
        !    write(6,'(A,1x,I6,1x,E15.7)') &
        !          'RuKbm: i_eq, fbio(i_eq) ', &
        !                  i_code_equation, fbio(i_CODE_equation)
        !    write(GA_print_unit,'(A,1x,I6,1x,E15.7)') &
        !          'RuKbm: i_eq, fbio(i_eq) ', &
        !                  i_code_equation, fbio(i_CODE_equation)
        !enddo ! i_CODE_equation


        !write(GA_print_unit,'(/A)') ' '
        !write(6,'(/A)') ' '


        ! Capture any export terms, or boundary conditions after bio flow is calculated
        ! If the mode does not contain any of these, SecondaryForcing() should do nothing

        !!!!call SecondaryForcing(fbio)


        do  i_Variable=1,n_Variables

            !kval(iter,i_Variable) = Delta_Time_in_Days * fbio(i_Variable)
            kval(iter,i_Variable) = dt * fbio(i_Variable)

            if( iter .eq. 1) then

                btmp(i_Variable) = b_tmp(i_Variable) + (kval(iter,i_Variable)/2.0D+0)

            elseif( iter .eq. 2) then

                btmp(i_Variable) = b_tmp(i_Variable) + (kval(iter,i_Variable)/2.0D+0)

            elseif( iter .eq. 3) then

                btmp(i_Variable) = b_tmp(i_Variable) + kval(iter,i_Variable)

            elseif( iter .eq. 4) then

                cff = (kval(1,i_Variable)/6.0D+0) + &
                      (kval(2,i_Variable)/3.0D+0) + &
                      (kval(3,i_Variable)/3.0D+0) + &
                      (kval(4,i_Variable)/6.0D+0)

                b_tmp(i_Variable) = b_tmp(i_Variable)+cff

            endif

        enddo ! End Kval loop

        !write(6,'(/A,1x,I6,3(1x,E15.7)/)') 'rkbm: iter, btmp(1:n_eqs)' , &
        !                                          iter, btmp(1:n_code_equations)

        !write(6,'(/A,1x,I6,3(1x,E15.7)/)') 'rkbm: iter, b_tmp(1:n_eqs)' , &
        !                                          iter, b_tmp(1:n_code_equations)

        !write(GA_print_unit,'(/A,1x,I6,3(1x,E15.7)/)') 'rkbm: iter, btmp(1:n_eqs)' , &
        !                                                      iter, btmp(1:n_code_equations)

        !write(GA_print_unit,'(/A,1x,I6,3(1x,E15.7)/)') 'rkbm: iter, b_tmp(1:n_eqs)' , &
        !                                                      iter, b_tmp(1:n_code_equations)
    enddo ! End iter loop

    !---------------------------------------------------------------------------

    ! if b_tmp is bad on any time step, then return with a bad result

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > 1.0d20 ) ) then

        L_bad_result = .TRUE.

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'RuKbm: bad result myid, i_time_step, b_tmp ', &
        !                             myid, i_time_step, b_tmp(1:n_CODE_equations)
        !endif ! L_GP_print

        !write(6,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'rkbm: bad result myid, i_time_step, b_tmp ', &
        !                            myid, i_time_step, b_tmp(1:n_CODE_equations)
        !if( L_GA_print )then
        !    write(GA_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'rkbm: bad result myid, i_time_step, b_tmp ', &
        !                            myid, i_time_step, b_tmp(1:n_CODE_equations)
        !endif ! L_GA_print

        return

    endif !   any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > 1.0d20

    !---------------------------------------------------------------------------

    Numerical_CODE_Solution(i_Time_Step,1:n_Variables)=max(b_tmp(1:n_Variables),0.0D+0)


    !write(6,'(//A,2(1x,I6),12(1x,E15.7))') &
    !      'rkbm: myid, i_time_step, b_tmp ', &
    !             myid, i_time_step, b_tmp(1:n_CODE_equations)
    !write(GA_print_unit,'(//A,2(1x,I6),12(1x,E15.7))') &
    !      'rkbm: myid, i_time_step, b_tmp ', &
    !             myid, i_time_step, b_tmp(1:n_CODE_equations)

    if( myid == 0 )then
        write(6,'(A,1x,I6,1x,6(1x,E15.7))') 'rkbm: i_time_step, solution ', &
                       i_Time_Step, Numerical_CODE_Solution(i_Time_Step,1:n_Variables)
    endif ! myid == 0

    if( L_ga_print )then ! .and. myid == 1 )then
        write(GA_print_unit,'(A,1x,I6,1x,6(1x,E15.7))') &
              'rkbm: i_time_step, solution ', &
                     i_Time_Step, Numerical_CODE_Solution(i_Time_Step,1:n_Variables)
    endif ! L_ga_print .and. myid == 1




enddo ! End Time step loop


!if( myid == 0 )then
!    write(6,'(A)') 'rkbm: leave Runge_Kutta_Box_Model '
!endif ! myid == 0

!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(A)') 'rkbm: leave Runge_Kutta_Box_Model '
!endif ! L_ga_print .and. myid == 1


return

end subroutine Runge_Kutta_Box_Model
