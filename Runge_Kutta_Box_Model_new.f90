subroutine Runge_Kutta_Box_Model_new

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use mpi
use mpi_module

use GA_Parameters_module
use GP_Parameters_module
use GP_Variables_module
use Runge_Kutta_Variables_module

implicit none

real(kind=8),parameter  :: one_third = 1.0d0 / 3.0d0
real(kind=8),parameter  :: one_half  = 0.5d0
real(kind=8) :: cff

integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node
integer(kind=4) :: i_level
integer(kind=4) :: i_function

integer(kind=4) :: i_node_left
integer(kind=4) :: i_node_right

integer(kind=4) :: i_CODE_equation
integer(kind=4) :: j_CODE_equation

integer(kind=4) :: i_time_step

integer(kind=4) :: icff
integer(kind=4) :: iter
real(kind=8) :: left_node_value
real(kind=8) :: right_node_value

logical,parameter :: L_GP_print = .TRUE.

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

left_node_value  = 0.0D0
right_node_value = 0.0D0


! start the time stepping loop

do  i_time_step=1,n_time_steps

    b_tmp = Runge_Kutta_Solution( i_time_step-1 , 1:n_CODE_equations )  ! Array Assignment

    ! carry out a Runge-Kutta time step

    do  iter=1,4

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !   Evaluate the trees
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        tree_value=0.0D+0

        !jjm note:  could you initialize tree_evaluationa and node_eval_type only once
        !  and use another array in place of tree_evaluation after the select case statement?

        ! jjm note: would it save a small amount of time to use node instead of i_node_left and
        ! use node+1 instead of i_node_right ?


        !Tree_Evaluation = Runge_Kutta_Node_Parameters       ! Matrix Assignment
        !Node_Eval_Type  = Runge_Kutta_Node_Type             ! Matrix Assignment
        !  array( 1:n_elements)

        !type treex
        !     integer :: i_node
        !     integer :: i_tree
        !     integer :: i_node_type
        !     read(kind=8) :: value
        !end type treex


        do  i_tree = 1, n_trees

            if( .not. any( i_tree == array%i_tree ) ) cycle

            do  i_level=n_levels-1,1,-1   ! move up the tree structure from level "n_level-1" to level "1"

                ! the function number at the right end of the upper level

                i_function= pow2_table(i_level-1) !  2**(i_level-1) - 1


                ! run through each function at the level

                do  i_node = pow2_table(i_level)+1 , pow2_table(i_level+1) , 2  ! 2**i_level,(2*(2**i_level))-1,2

                    if( .not. any(  i_node ==  array%i_node ) ) cycle


                    i_function=i_function+1       ! sets the 'function' node's index

                    i_node_left=i_node            ! sets the 'left terminal' node's index;
                                                  !    i_node_left=i_function*2 would also work

                    i_node_right=i_node+1         ! sets the 'right terminal' node's index;
                                                  !    i_node_right=(i_function*2)+1 would also work

                    call get_node_type( i_function, i_tree, i_node_type1, value1 )

                    if( i_node_type1         .gt. 0) then ! run the calculation
                    !if( Node_Eval_Type(i_function,i_tree) .gt. 0) then ! run the calculation

                        icff= i_node_type1    ! Node_Eval_Type(i_node_left,i_tree)

                        if( icff .eq. 0) then

                            left_node_value= value1 ! Tree_Evaluation(i_node_left,i_tree)

                        elseif( icff .lt. 0 .and. icff .ne. -9999) then

                            if( iter .eq. 1) then
                                left_node_value=b_tmp(abs(icff))
                            else
                                left_node_value=btmp(abs(icff))
                            endif ! iter .eq. 1

                        endif ! icff .eq. 0

                        call get_node_type( i_node_right, i_tree, i_node_type2, value2 )

                        icff= i_node_type2

                        if( icff .eq. 0) then

                            right_node_value=  value2

                        elseif( icff .lt. 0 .and. icff .ne. -9999) then

                            if( iter .eq. 1) then
                                right_node_value=b_tmp(abs(icff))
                            else
                                right_node_value=btmp(abs(icff))
                            endif ! iter .eq. 1

                        endif ! icff .eq. 0

                        ! evaluate the node set

                        ! Function types used
                        ! Type 1: ==> Addition  left + right
                        ! Type 2: ==> Subtraction  left - right
                        ! Type 3: ==> Multiply  left * right
                        ! Type 4: ==> Divide (protected) left / right
                        ! Type 5: ==> Ivlev Grazing Function ==> (1 - e^-abs(left*right))
                        ! Type 6: ==> Michaelis-Menton Term (modified for Forward-Backward)
                        !                                    (1 / (abs(LHS) + abs(RHS)))
                        ! Type 7: ==> Mayzaud-Poulet Grazing Function ==>
                        !                     abs(left*right)*(1 -e^-abs(left*right))


                        SELECT CASE( i_node_type1  )


                           CASE(1)  ! LHS + RHS

                             Tree_Evaluation(i_function,i_tree)=left_node_value + right_node_value


                           CASE(2)  ! LHS - RHS

                             Tree_Evaluation(i_function,i_tree)=left_node_value-right_node_value


                           CASE(3)  ! LHS * RHS

                             Tree_Evaluation(i_function,i_tree)=left_node_value*right_node_value


                           CASE(4)  ! protected: LHS/RHS


                             if( right_node_value .ne. 0.0D+0) then
                             !if( abs(right_node_value) > 1.0d-50 ) then

                                 Tree_Evaluation(i_function,i_tree)=left_node_value/right_node_value

                             else

                                 Tree_Evaluation(i_function,i_tree)=0.0d0

                             endif


                           CASE(5)  ! '1.0D+0 - dexp(-dabs(LHS*RHS))'


                             cff=dabs(left_node_value*right_node_value)

                             Tree_Evaluation(i_function,i_tree) = 1.0D+0 - dexp(-1.0D+0*cff)



                           CASE(6)  ! 'Michealis-Menton (abs(RHS) / (abs(LHS) + abs(RHS)))'

                             cff = abs(left_node_value) + abs(right_node_value)

                             if( cff .gt. 0.0D+0 ) then
                                 Tree_Evaluation(i_function,i_tree) = abs(right_node_value) / cff
                             else
                                 Tree_Evaluation(i_function,i_tree)=0.0D+0
                             endif



                           CASE(7)  ! 'abs(LHS*RHS)*(1.0D+0 - dexp(-dabs(LHS*RHS)))'

                             cff=dabs(left_node_value*right_node_value)

                             Tree_Evaluation(i_function,i_tree)=cff*(1.0D+0-dexp(-1.0D+0*cff))

                             !if( L_ga_print )then
                             !    write(GA_print_unit,'(8x, A, 1x,I2,3(1x,E15.7) )') &
                             !     'icff, left, right, tree_eval ', &
                             !      icff, left_node_value, right_node_value, tree_evaluation(i_function,i_tree)
                             !endif ! L_ga_print

                           CASE DEFAULT

                             if( L_ga_print )then
                                 write(GA_print_unit,*) 'wrong case number chosen in Runge Kutta evaluations'
                             endif ! L_ga_print
                             stop 'RK bad case number'

                        END SELECT

                        Node_Eval_Type(i_function,i_tree)=0  ! now there is a value/variable

                        !node_eval_type(i_node_left,i_tree)  = 0
                        !node_eval_type(i_node_right,i_tree) = 0
                        !node_eval_type(i_function,i_tree)   = 0

                    endif !  Node_Eval_Type(i_function,i_tree) .gt. 0

                enddo !  i_node

            enddo !  i_level

            Tree_Value(i_tree)=Tree_Evaluation(1,i_tree)

            !if( L_ga_print )then
            !    write(GA_print_unit,'(8x, A, 1x,E15.7)') &
            !          'RuKbm: tree_eval(1, i_tree )  = ', &
            !           tree_evaluation(1, i_tree )
            !endif ! L_ga_print

        enddo !  i_tree

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !   Calculate the flow terms from the determined tree_value terms
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        i_tree=0
        do  i_CODE_equation=0,n_CODE_equations   ! source of material

            do  j_CODE_equation=0,n_CODE_equations ! sink of material

                if( i_CODE_equation .ne. j_CODE_equation) then

                    i_tree=i_tree+1

                    ! 'dabs' forces flow of material in one direction

                    bioflo(i_CODE_equation,j_CODE_equation)=dabs(tree_value(i_tree))

                    !if( L_ga_print )then
                    !    write(GA_print_unit,'(A,3(1x,I6),2(2x, E15.7))') &
                    !    'RuKbm:1 i_tree, i_code_equation, j_code_equation, &
                    !    &tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation) ', &
                    !            i_tree, i_code_equation, j_code_equation, &
                    !     tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation)
                    !endif ! L_ga_print

                else

                    bioflo(i_CODE_equation,j_CODE_equation)=0.0D+0   ! never flow to/from same component

                endif !   i_CODE_equation .ne. j_CODE_equation

            enddo !  j_CODE_equation

        enddo !  i_CODE_equation

        do  i_CODE_equation=1,n_CODE_equations   ! source of material
            fbio(i_CODE_equation)=0.0D+0
        enddo ! i_CODE_equation

        !   bring in the component flow sources and sinks

        do  i_CODE_equation=0,n_CODE_equations   ! source of material

            do  j_CODE_equation=0,n_CODE_equations ! sink of material

                if( i_CODE_equation .gt. 0) then

                    fbio(i_CODE_equation) = &
                      fbio(i_CODE_equation) - bioflo(i_CODE_equation,j_CODE_equation)

                endif !   i_CODE_equation .gt. 0

                if( j_CODE_equation .gt. 0) then

                    fbio(j_CODE_equation) = &
                      fbio(j_CODE_equation) + bioflo(i_CODE_equation,j_CODE_equation)

                endif !   j_CODE_equation .gt. 0

            enddo !  j_CODE_equation
        enddo !  i_CODE_equation


        do  i_CODE_equation=1,n_CODE_equations

            kval(iter,i_CODE_equation) = dt * fbio(i_CODE_equation)

            !write(6,'(A,1x,I1,1x,I6,1x,i1,3(1x,E24.16))') &
            !          'RuKbm:1', myid, iter, i_CODE_equation, &
            !                   kval(iter,i_code_equation),dt, &
            !                    fbio(i_CODE_equation)

            if( iter .eq. 1) then

                btmp(i_CODE_equation) = b_tmp(i_CODE_equation) + &
                                        (kval(iter,i_CODE_equation)/2.0D+0)
                !btmp(i_CODE_equation) = b_tmp(i_CODE_equation) + &
                !                        kval(iter,i_CODE_equation) * one_half

            elseif( iter .eq. 2) then

                btmp(i_CODE_equation) = b_tmp(i_CODE_equation) + &
                                        (kval(iter,i_CODE_equation)/2.0D+0)
                !btmp(i_CODE_equation) = b_tmp(i_CODE_equation) + &
                !                        kval(iter,i_CODE_equation) * one_half

            elseif( iter .eq. 3) then

                btmp(i_CODE_equation) = b_tmp(i_CODE_equation) + &
                                        kval(iter,i_CODE_equation)

            elseif( iter .eq. 4) then

                cff =  kval(1,i_CODE_equation)/6.0D+0  + &
                       kval(2,i_CODE_equation)/3.0D+0  + &
                       kval(3,i_CODE_equation)/3.0D+0  + &
                       kval(4,i_CODE_equation)/6.0D+0

                !write(6,'(A,1x,I1,1x,I6,1x,i1,1x,I1,5(1x,E24.16))') &
                !        'RuKbm:', myid, i_time_step, iter, i_CODE_equation, &
                !                 kval(1:4,i_code_equation),cff

                !write(6,'(/A,4(1x,I6))') &
                !        'RuKbm: myid, i_time_step, iter, i_CODE_equation  ', &
                !                myid, i_time_step, iter, i_CODE_equation
                !write(6,'(A,2(1x,I6),5(1x,E15.7))') &
                !        'RuKbm: myid, i_CODE_equation, kval(1:4,i_code_equation),cff  ', &
                !                myid, i_CODE_equation, kval(1:4,i_code_equation),cff

                !cff =  ( kval(1,i_CODE_equation) * one_half  + &
                !         kval(2,i_CODE_equation)             + &
                !         kval(3,i_CODE_equation)             + &
                !         kval(4,i_CODE_equation) * one_half ) * one_third

                !cff = ( ( kval(1,i_CODE_equation) + kval(4,i_CODE_equation) ) * one_half + &
                !          kval(2,i_CODE_equation) + kval(3,i_CODE_equation)             ) * one_third

                !write(6,'(A,1x,I1,1x,I6,1x,i1,1x,I1,5(1x,E24.16))') &
                !        'RuKbm:', myid, i_time_step, iter, i_CODE_equation, &
                !                kval(1:4,i_code_equation),cff

                b_tmp(i_CODE_equation)=b_tmp(i_CODE_equation)  +  cff

            endif !   iter .eq. 1

        enddo !  i_CODE_equation

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    enddo ! iter



    !---------------------------------------------------------------------------

    ! if b_tmp is bad on any time step, then return with a bad result

    if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > 1.0d20 ) ) then

        L_bad_result = .TRUE.

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
        !          'RuKbm: bad result myid, i_time_step, b_tmp ', &
        !                             myid, i_time_step, b_tmp(1:n_CODE_equations)
        !endif ! L_ga_print

        return

    endif !   any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > 1.0d20

    !---------------------------------------------------------------------------

    !   Array Assignment with no negative values

    !Runge_Kutta_Solution(i_time_step,1:n_CODE_equations) = max(b_tmp,0.0d+0)

    ! NPZ model

    Runge_Kutta_Solution(i_time_step,1:n_CODE_equations) = max( b_tmp, 0.0d0 )

    !---------------------------------------------------------------------------

    !if( L_GP_print .and. GP_para_flag .and. myid == 1 .and. mod( i_time_step, 50 ) == 0 )then
    !    write(GP_print_unit,'(A,2(1x,I6),12(1x,E15.7))') &
    !          'RuKbm: myid, i_time_step, RK_Soln ', &
    !                  myid, i_time_step, Runge_Kutta_Solution(i_time_step,1:n_CODE_equations)
    !endif ! L_ga_print



enddo ! i_time_step


return

end subroutine Runge_Kutta_Box_Model_new
