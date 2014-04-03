subroutine Runge_Kutta_Box_Model( L_print_RK ) 

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
!integer(kind=4),intent(in) :: i_GA_indiv

logical :: L_GP_print
logical,intent(in) :: L_print_RK 

!--------------------------------------------------------------------------------------
L_GP_print = .true.

tree_node_count = 0

!if( L_ga_print )then 
!    write(GA_print_unit,'(/A,1x,I6/)') 'rkbm: entry Runge_Kutta_Box_Model myid = ', myid
!    write(GA_print_unit,'(A,1x,I6/)')  'rkbm: n_Variables ', n_Variables
!!    !!!write(GA_print_unit,'(A,1x,I6/)')  'rkbm: tree_value modified'
!endif ! L_ga_print

!write(6,'(/A,1x,I6)')  'rkbm: n_Variables     ', n_Variables
!write(6,'(A,1x,I6/)')  'rkbm: n_code_equations', n_code_equations
!flush(6)
!flush(GA_print_unit)

!--------------------------------------------------------------------------------------
!! debug >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! debug only  - put in discover problem tree
!GP_individual_node_type(:, :) =  -9999
!
!i_tree = 1
!i_node = 1
!GP_individual_node_type(i_node, i_tree) =  6
!i_node = 2
!GP_individual_node_type(i_node, i_tree) =  0
!i_node = 3
!GP_individual_node_type(i_node, i_tree) =  4
!i_node = 6
!GP_individual_node_type(i_node, i_tree) = -1
!i_node = 7
!GP_individual_node_type(i_node, i_tree) =  7
!i_node = 14
!GP_individual_node_type(i_node, i_tree) = -2
!i_node = 15
!GP_individual_node_type(i_node, i_tree) = -1
!
!i_tree = 5
!i_node = 1
!GP_individual_node_type(i_node, i_tree) = -2
!
!
!GP_individual_node_parameters(:, :) = 0.0D0
!i_tree = 1
!i_node = 2
!GP_individual_node_parameters(i_node, i_tree) = 0.7191251516342163D+02
!
!Numerical_CODE_Solution( 0 , 1) = 0.6718252785503864D-02
!Numerical_CODE_Solution( 0 , 2) = 0.8888030052185059D+02
!
!!--------------------------------------------------------------------------------------
!


if( dt <= 0.0d0 )then
    if( myid == 0 .or. myid == 1 )then
        write(6,'(/A/)') 'rkbm: BAD VALUE for dt'
        write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt
    endif ! myid == 0 .or. myid == 1
    call MPI_FINALIZE(ierr)
    stop 'bad delta_time'
endif ! dt <= 0.0D0

!write(6,'(/A,1x,E20.10/)') 'rkbm: dt', dt

!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop Numerical_CODE_Solution(0,:)', &
!                         Numerical_CODE_Solution(0,:)
!    write(GA_print_unit,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)
!endif ! L_ga_print .and. myid == 1

!write(6,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop Numerical_CODE_Solution(0,:)', &
!                         Numerical_CODE_Solution(0,:)
!write(6,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)

!write(6,'(A,1x,E20.10/)') 'rkbm: dt', dt


!!! debug >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_type(i_node, i_tree) > -9999 )then
!            write(GA_print_unit,'(A,3(1x,I4))')&
!                  'rkbm: i_tree, i_node, GP_indiv_node_type', &
!                          i_tree, i_node, GP_individual_node_type(i_node, i_tree)
!        endif !  GP_individual_node_type(i_node, i_tree) > -9999
!    enddo ! i_node
!enddo ! i_tree
!
!do  i_tree = 1, n_trees
!    do  i_node = 1, n_nodes
!        if( GP_individual_node_parameters(i_node, i_tree) > 0.0d0 )then
!            write(GA_print_unit,'(A,2(1x,I4),1x,E15.7)')&
!                  'rkbm: i_tree, i_node, GP_indiv_node_parms', &
!                          i_tree, i_node, GP_individual_node_parameters(i_node, i_tree)
!        endif !  GP_individual_node_parameters(i_node, i_tree) > 0.0d0
!    enddo ! i_node
!enddo ! i_tree
!! debug <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!flush(6)
!flush(GA_print_unit)


! start the time stepping loop


do  i_Time_Step = 1, n_Time_Steps


    !------------------------------------------------------------------------------

    !!!RK_data_array(1:n_input_vars) = input_data_array(1:n_input_vars, i_data_point )

    b_tmp(:) = Numerical_CODE_Solution(i_Time_Step-1,:)  ! Array Assignment


    !write(GA_print_unit,'(/A,1x,I6,10(1x,E15.7)/)') &
    !      'rkbm: i_time_step, b_tmp(1:n_eqs)' , &
    !             i_time_step, b_tmp(1:n_code_equations)
    !write(6,'(/A,1x,I6,10(1x,E15.7)/)') &
    !      'rkbm: i_time_step, b_tmp(1:n_eqs)' , &
    !             i_time_step, b_tmp(1:n_code_equations)

    !flush(6)
    !flush(GA_print_unit)

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > 1.0d20 ) ) then

        !write(6,'(/A,1x,I6/)') &
        !     'rkbm: bad b_tmp  i_time_step', i_time_step 
        !flush(6)

        L_bad_result = .TRUE.

        return
    endif !  any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > 1.0d20 )


    btmp = b_tmp

    !write(6,'(/A,10(1x,E15.7)/)') &
    !      'rkbm: btmp(1:n_eqs)' , &
    !             btmp(1:n_code_equations)
    !flush(6)
    !flush(GA_print_unit)
    !------------------------------------------------------------------------------

    ! carry out a Runge-Kutta time step
    do  iter=1,4

        !if( myid == 0 )then
        !if( L_ga_print )then
        !    write(GA_print_unit,'(A,2(1x,I6))') 'rkbm: i_time_step, iter ', &
        !                                               i_time_step, iter
        !endif ! L_ga_print
        !!endif ! myid == 0
        !write(6,'(A,2(1x,I6))') 'rkbm: i_time_step, iter ', &
        !                               i_time_step, iter
        !flush(6)
        !flush(GA_print_unit)

        !--------------------------------------------------------------------------
        ! Call forcing functions for the box model

        !write(6,'(/A)') 'rkbm: call DoForcing'
        !flush(6)
        !flush(GA_print_unit)

        call DoForcing( btmp, Runge_Kutta_Time_Step(iter), i_Time_Step )

        !write(6,'(A/)') 'rkbm: aft call DoForcing'
        !flush(6)
        !flush(GA_print_unit)
        !--------------------------------------------------------------------------

        fbio = 0.0D+0

        do  i_Track = 1,n_Tracked_Resources

            !!!! Call Model_Diagnostics()

            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Evaluate the trees
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            Tree_Value = 0.0D+0                            ! Matrix Assignment

            do  i_Tree=1,n_Trees



                if( associated( GP_Trees(i_Tree, i_Track)%n) ) then

                    !write(6,'(//A,1x,I6,5x,L1)') &
                    !      'rkbm: i_tree, associated(GP_Trees(i_Tree, i_Track)%n)  ', &
                    !             i_tree, associated(GP_Trees(i_Tree, i_Track)%n)

                    !write(GA_print_unit,'(/A,2x,I6)') &
                    !          'rkbm: bef size( GP_Trees ) ', size( GP_Trees )

                    Tree_Value(i_Tree) = GP_Trees( i_Tree,  i_Track )%n%val()

                    !if( i_time_step < 251 ) then
                    !    write(GA_print_unit,'(A,22x,I6,1x,I6,1x,E15.7)') &
                    !              'rkbm: iter, i_tree, Tree_Value(i_tree)', &
                    !                     iter, i_tree, Tree_Value(i_tree)
                    !    !if( myid == 0 )then

                    !if( abs( Tree_Value(i_tree) ) > 0.0d0 )then
                    !    write(6,'(A,22x,I6,1x,I6,1x,E15.7)') &
                    !          'rkbm: iter, i_tree, Tree_Value(i_tree)', &
                    !                 iter, i_tree, Tree_Value(i_tree)
                    !endif ! abs( Tree_Value(i_tree) ) > 0.0d0 
                    !flush(6)
                    !flush(GA_print_unit)

                    !    !endif ! myid == 0
                    !endif ! i_time_step < 251

                    if( isnan( Tree_Value(i_Tree) )          .or.   &
                          abs( Tree_Value(i_Tree) )  > 1.0d20    ) then

                        L_bad_result = .TRUE.

                        !if( L_ga_print )then
                        !    write(GA_print_unit,'(A,1x,I6,1x,I6,1x,E24.16)') &
                        !      'rkbm: bad value i_time_step, i_tree, Tree_Value(i_tree)', &
                        !                       i_time_step, i_tree, Tree_Value(i_tree)
                        !endif ! L_ga_print 
                        !write(6,'(A,1x,I6,1x,I6,1x,E24.16)') &
                        !      'rkbm: bad value i_time_step, i_tree, Tree_Value(i_tree)', &
                        !                       i_time_step, i_tree, Tree_Value(i_tree)
                        !flush(6)
                        !flush(GA_print_unit)

                        return
                    endif !  isnan( Tree_Value(i_Tree) ) .or. abs(Tree_Value(i_Tree)) > 1.0d20 

                    !write(GA_print_unit,'(/A,2x,I6)') &
                    !          'rkbm: aft size( GP_Trees ) ', size( GP_Trees )

                    !----------------------------------------------------------------------------
                    !tree_node_count = GetNodeCount( GP_Trees( i_Tree, i_Track )%n )
                    !debug only !if( tree_node_count <= 1 ) Tree_Value(i_Tree) = 0.0d0 !jjm 20131213
                    !---------------------------------------------------------------------------

                    !write(GA_print_unit,'(A,22x,I6,1x,E15.7)') &
                    !'rkbm: i_tree,Tree_Value(i_tree) ',i_tree,Tree_Value(i_tree)


                    !if( myid == 0 )then
                    !write(GA_print_unit,'(A,22x,I6,1x,I6/)') &
                    !'rkbm: i_tree,tree_node_count    ',i_tree,tree_node_count
                    !endif ! myid == 0



                endif ! associated(GP_Trees...

            enddo ! i_Trees


            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Calculate the flow terms from the determined tree_value terms
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            !if( myid == 0 )then
            ! write(6,'(/A/)') 'rkbm: Calculate the flow terms from the determined tree_value terms'
            !endif ! myid == 0

            !write(GA_print_unit,'(/A)') ' '



            i_Tree=0
            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material
                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .ne. j_CODE_Equation) then

                        i_Tree=i_Tree+1

                        ! 'abs' forces flow of material in one direction

                        bioflo(i_CODE_Equation,j_CODE_Equation)=abs(Tree_Value(i_Tree))


                        !if( L_ga_print )then
                        !    write(GA_print_unit,'(A,3(1x,I6),2(2x, E15.7))') &
                        !    'rkbm:1 i_tree, i_eq, j_eq, tree_value, bioflo(i_eq,j_eq) ', &
                        !            i_tree, i_code_equation, j_code_equation, &
                        !     tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation)
                        !endif ! L_ga_print
                        !write(6,'(A,3(1x,I6),2(2x, E15.7))') &
                        !    'rkbm:1 i_tree, i_eq, j_eq, tree_value, bioflo(i_eq,j_eq) ', &
                        !            i_tree, i_code_equation, j_code_equation, &
                        !     tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation)


                    else

                        ! never flow to/from same component
                        bioflo(i_CODE_Equation,j_CODE_Equation)=0.0D+0

                    endif ! i_CODE_Equation .ne. j_CODE_Equation

                enddo ! j_CODE_Equation
            enddo ! i_CODE_Equation

            !write(GA_print_unit,'(/A)') ' '

            !if( i_time_step < 251 ) then
            !    write(6,'(A,3(1x,I6))') &
            !          'rkbm: bioflo_map(1:n_eqs, 1) ', &
            !                 bioflo_map(1:n_code_equations, 1)
            !endif ! i_time_step < 251 

            ! bring in the component flow sources and sinks

            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material

                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .gt. 0 ) then

                        !write(60,'(A,3(1x,I6))') &
                        !      'rkbm: i_eq, bioflo_map(i_eq, 1) ', &
                        !              i_code_equation, bioflo_map(i_code_equation, 1)

                        if( bioflo_map(i_CODE_Equation,i_Track) .gt. 0 ) then

                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'rkbm: bef i_eq, fbio(bioflo_map(i_eq, 1)) ', &
                            !                  i_code_equation, fbio(bioflo_map(i_code_equation, 1))

                            fbio(bioflo_map(i_CODE_Equation,i_Track)) = &
                                fbio(bioflo_map(i_CODE_Equation,i_Track)) -  &
                                         bioflo(i_CODE_Equation,j_CODE_Equation)

                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'rkbm: aft i_eq, fbio(bioflo_map(i_eq, 1)) ', &
                            !                  i_code_equation, fbio(bioflo_map(i_code_equation, 1))
                            !if( abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0 )then 
                            !    write(6,'(A,2(1x,I6),1x,E20.10)') &
                            !          'rkbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !                 i_code_equation, j_code_equation, &
                            !                 bioflo(i_CODE_Equation,j_CODE_Equation)
                            !endif ! abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0 

                        endif ! bioflo_map(i_CODE_Equation,i_Track) .gt. 0

                    endif ! i_CODE_Equation .gt. 0

                    if( j_CODE_Equation .gt. 0 ) then

                        !write(60,'(A,3(1x,I6))') &
                        !      'rkbm: j_eq, bioflo_map(j_eq, 1) ', &
                        !              j_code_equation, bioflo_map(j_code_equation, 1)

                        if( bioflo_map(j_CODE_Equation,i_Track) .gt. 0 ) then

                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'rkbm: bef j_eq, fbio(bioflo_map(j_eq, 1)) ', &
                            !                  j_code_equation, fbio(bioflo_map(j_code_equation, 1))

                            fbio(bioflo_map(j_CODE_Equation,i_Track)) = &
                                 fbio(bioflo_map(j_CODE_Equation,i_Track)) + &
                                          bioflo(i_CODE_Equation,j_CODE_Equation)

                            !write(60,'(A,1x,I6,1x,E20.10)') &
                            !      'rkbm: aft j_eq, fbio(bioflo_map(j_eq, 1)) ', &
                            !                  j_code_equation, fbio(bioflo_map(j_code_equation, 1))
                            !if( abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0 )then
                            !    write(6,'(A,2(1x,I6),1x,E20.10)') &
                            !          'rkbm: i_eq, j_eq, bioflo(i_eq,j_eq) ', &
                            !                 i_code_equation, j_code_equation, &
                            !                 bioflo(i_CODE_Equation,j_CODE_Equation)
                            !endif !  abs(bioflo(i_CODE_Equation,j_CODE_Equation)) > 0.0d0 

                        endif ! bioflo_map(j_CODE_Equation,i_Track) .gt. 0

                    endif ! j_CODE_Equation .gt. 0

                enddo ! j_CODE_Equation
            enddo ! i_CODE_Equation

        enddo ! End Tracked Resources loop


        !if( i_time_step < 251 ) then
        !    write(6,'(A)') ' '
        !    do  i_CODE_equation=1,n_CODE_equations   ! source of material
        !        write(6,'(A,1x,I6,1x,E15.7)') &
        !              'rkbm: i_eq, fbio(i_eq) ', &
        !                     i_code_equation, fbio(i_CODE_equation)
        !    enddo ! i_CODE_equation
        !    write(6,'(A)') ' '
        !endif ! i_time_step < 251

        !write(GA_print_unit,'(/A)') ' '
        !do  i_CODE_equation=1,n_CODE_equations   ! source of material
        !    write(GA_print_unit,'(A,1x,I6,1x,E15.7)') &
        !          'rkbm: i_eq, fbio(i_eq) ', &
        !                  i_code_equation, fbio(i_CODE_equation)
        !enddo ! i_CODE_equation
        !write(GA_print_unit,'(/A)') ' '


        ! Capture any export terms, or boundary conditions after bio flow is calculated
        ! If the mode does not contain any of these, SecondaryForcing() should do nothing

        !!!!call SecondaryForcing(fbio)


        do  i_Variable=1,n_Variables

            !kval(iter,i_Variable) = Delta_Time_in_Days * fbio(i_Variable)
            kval(iter,i_Variable) = dt * fbio(i_Variable)


            !if( i_time_step < 251 ) then
            !    write(6,'(A,1x,I1,1x,I1,1x,i1,2(1x,E15.7))') &
            !              'rkbm:myid, iter, i_eq, kval(iter,i_eq, fbio(i_eq)', &
            !                    myid, iter, i_variable, &
            !                      kval(iter,i_variable), &
            !                       fbio(i_variable)
            !    write(GA_print_unit,'(A,1x,I1,1x,I6,1x,i1,3(1x,E15.7))') &
            !              'rkbm:myid, iter, i_eq, kval(iter,i_eq, fbio(i_eq)', &
            !                    myid, iter, i_variable, &
            !                      kval(iter,i_variable), &
            !                       fbio(i_variable)
            !endif ! i_time_step < 251

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


                !if( i_time_step < 251 ) then
                !            write(6,'(A,1x,I1,1x,I6,1x,i1,1x,I1,5(1x,E15.7))') &
                !                    'rkbm:', myid, i_time_step, iter, i_variable, &
                !                             kval(1:4,i_variable),cff
                !            write(6,'(/A,4(1x,I6))') &
                !                    'rkbm: myid, i_time_step, iter, i_variable  ', &
                !                           myid, i_time_step, iter, i_variable
                !            write(6,'(A,2(1x,I6),5(1x,E15.7))') &
                !                    'rkbm: myid, i_variable, kval(1:4,i_variable),cff  ', &
                !                           myid, i_variable, kval(1:4,i_variable),cff
                !endif ! i_time_step < 251

                b_tmp(i_Variable) = b_tmp(i_Variable)+cff

            endif

        enddo ! End Kval loop

        !write(GA_print_unit,'(/A,1x,I6,3(1x,E15.7)/)') 'rkbm: iter, btmp(1:n_eqs)' , &
        !                                                      iter, btmp(1:n_code_equations)

        !write(GA_print_unit,'(/A,1x,I6,3(1x,E15.7)/)') 'rkbm: iter, b_tmp(1:n_eqs)' , &
        !                                                      iter, b_tmp(1:n_code_equations)
        !write(6,'(A,1x,I6,7(1x,E15.7))') 'rkbm: iter, btmp(1:n_eqs) ' , &
        !                                        iter, btmp(1:n_code_equations)

        !write(6,'(A,1x,I6,7(1x,E15.7))') 'rkbm: iter, b_tmp(1:n_eqs)' , &
        !                                        iter, b_tmp(1:n_code_equations)

        !if( i_time_step < 251 ) then
        !    write(6,'(A, 1x,i1,3(1x,E15.7) )') &
        !         'rkbm: iter, btmp(  1:n_CODE_equations )', &
        !                iter, btmp(  1:n_CODE_equations )
        !
        !    write(6,'(A, 1x,i1,3(1x,E15.7) )') &
        !         'rkbm: iter, b_tmp( 1:n_CODE_equations )', &
        !                iter, b_tmp( 1:n_CODE_equations )
        !
        !    write(GA_print_unit,'(A, 1x,i1,3(1x,E15.7) )') &
        !         'rkbm: iter, btmp( 1:n_CODE_equations )', &
        !                iter, btmp( 1:n_CODE_equations )
        !
        !    write(GA_print_unit,'(A, 1x,i1,3(1x,E15.7) )') &
        !         'rkbm: iter, b_tmp( 1:n_CODE_equations )', &
        !                iter, b_tmp( 1:n_CODE_equations )
        !endif ! i_time_step < 251

    enddo ! End iter loop

    !---------------------------------------------------------------------------

    ! if b_tmp is bad on any time step, then return with a bad result

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > 1.0d20 ) ) then

        L_bad_result = .TRUE.

        !if( L_GP_print )then
        !    write(6,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'rkbm: bad result myid, i_time_step, b_tmp ', &
        !                            myid, i_time_step, b_tmp(1:n_CODE_equations)
        !endif ! L_GP_print

        !if( L_GA_print )then
        !    write(GA_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'rkbm: bad result myid, i_time_step, b_tmp ', &
        !                            myid, i_time_step, b_tmp(1:n_CODE_equations)
        !endif ! L_GA_print

        !flush(6)
        !flush(GA_print_unit)
        return

    endif !   any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > 1.0d20

    !---------------------------------------------------------------------------

    Numerical_CODE_Solution(i_Time_Step,1:n_Variables)=max(b_tmp(1:n_Variables),0.0D+0)


    !write(GA_print_unit,'(//A,2(1x,I6),12(1x,E15.7))') &
    !      'rkbm: myid, i_time_step, b_tmp ', &
    !             myid, i_time_step, b_tmp(1:n_CODE_equations)
    !write(6,'(//A,2(1x,I6),12(1x,E15.7))') &
    !      'rkbm: myid, i_time_step, b_tmp ', &
    !             myid, i_time_step, b_tmp(1:n_CODE_equations)

    !if( myid == 0 )then
    !if( L_ga_print )then 
    !    write(GA_print_unit,'(/A,1x,I6,1x,10(1x,E15.7)/)') 'rkbm: i_time_step, solution ', &
    !                  i_Time_Step, Numerical_CODE_Solution(i_Time_Step,1:n_Variables)
    !endif ! L_ga_print
    !     write(6,'(A,1x,I6,1x,10(1x,E15.7))') 'rkbm: i_time_step, solution ', &
    !                     i_Time_Step, Numerical_CODE_Solution(i_Time_Step,1:n_Variables)
    !!endif ! myid == 0

    !flush(6)
    !flush(GA_print_unit)

    !if( L_ga_print )then ! .and. myid == 1 )then
    !    write(GA_print_unit,'(A,1x,I6,1x,6(1x,E15.7)/)') &
    !          'rkbm: i_time_step, solution ', &
    !                 i_Time_Step, Numerical_CODE_Solution(i_Time_Step,1:n_Variables)
    !endif ! L_ga_print .and. myid == 1

    !if( i_time_step < 251 ) then
    !!if( i_time_step == 250 .or. i_time_step == 1 ) then
    !    if( myid == 0 )then
    !        write(6,'(A,2(1x,I6),12(1x,E15.7))') &
    !        'rkbm:6 myid, i_time_step, RK_Soln ', &
    !                myid, i_time_step, Numerical_CODE_Solution(i_time_step,1:n_CODE_equations)
    !    endif ! myid == 0

    !    if( L_GA_print )then
    !        write(GA_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
    !        'rkbm:g myid, i_time_step, RK_Soln ', &
    !                myid, i_time_step, Numerical_CODE_Solution(i_time_step,1:n_CODE_equations)
    !    endif ! L_ga_print
    !!endif !  i_time_step == 250 .or. i_time_step == 1
    !endif ! i_time_step < 251




enddo ! End Time step loop


!if( myid == 0 )then
!    write(6,'(A)') 'rkbm: leave Runge_Kutta_Box_Model '
!endif ! myid == 0

!if( L_ga_print )then ! .and. myid == 1 )then
!    write(GA_print_unit,'(/A/)') 'rkbm: leave Runge_Kutta_Box_Model '
!endif ! L_ga_print .and. myid == 1

!flush(6)

return

end subroutine Runge_Kutta_Box_Model
