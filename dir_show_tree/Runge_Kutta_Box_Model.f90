subroutine Runge_Kutta_Box_Model
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


!use mpi
!use mpi_module

use GP_Parameters_module
use GP_Variables_module
use Runge_Kutta_Variables_module

implicit none

real(kind=8) :: cff



character(str_len) ::  left_node_value_string
character(str_len) ::  right_node_value_string
character(str_len) ::  cff_string

character(4), allocatable,  dimension( : , : )      ::  node_eval_type_string
character(str_len),  allocatable, dimension( : )    ::  tree_value_string

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!write(6,*)'RKBM: allocate node_eval_type_string '
allocate( node_eval_type_string( n_nodes,n_trees) ) 

!write(6,*)'RKBM: allocate tree_value_string '
allocate( tree_value_string( n_trees) ) 


left_node_value_string    = ''
right_node_value_string    = ''
cff_string    = ''
!tree_evaluation_string = ''

node_eval_type_string    = ''
tree_value_string    = ''

! start the time stepping loop

do i_time_step=1,1  !n_time_steps

  b_tmp = Runge_Kutta_Solution(i_time_step-1,1:n_CODE_equations)  ! Array Assignment

  ! carry out a Runge-Kutta time step

  do iter=1, 1 !4

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !   Evaluate the trees
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    tree_value=0.0D+0

    Tree_Evaluation = Runge_Kutta_Node_Parameters       ! Matrix Assignment
    Node_Eval_Type  = Runge_Kutta_Node_Type             ! Matrix Assignment

    tree_evaluation_string = node_parameters_string 

    do i_tree=1,n_trees

      !write(6,'(/A,1x,I6)')'RKBM: i_tree = ', i_tree

      do i_level=n_levels-1,1,-1   ! move up the tree structure from level "n_level-1" to level "1"

        i_function=(2**(i_level-1))-1 ! the function number at the right end of the upper level

        !write(6,'(A,1x,I6)')'    RKBM: i_level    = ', i_level
        !write(6,'(A,1x,I6)')'    RKBM: i_function = ', i_function

        do i_node=2**i_level,(2*(2**i_level))-1,2     ! run through each function at the level

          !write(6,'(A,1x,I6)')'        RKBM: i_node     = ', i_node    

          i_function=i_function+1       ! sets the 'function' node's index
          i_node_left=i_node            ! sets the 'left terminal' node's index;
                                        !    i_node_left=i_function*2 would also work
          i_node_right=i_node+1         ! sets the 'right terminal' node's index;
                                        !    i_node_right=(i_function*2)+1 would also work

          !write(6,'(A,1x,I6)')'    RKBM: i_function = ', i_function

          !write(6,'(A,1x,I6)')'        RKBM: i_node_left  = ', i_node_left
          !write(6,'(A,1x,I6)')'        RKBM: i_node_right = ', i_node_right

          !write(6,'(A,3(1x,I6))') 'RKBM: i_function, i_tree, node_eval_type(i_function,i_tree) ', &
          !                               i_function, i_tree, node_eval_type(i_function,i_tree)

          if( Node_Eval_Type(i_function,i_tree) .gt. 0) then ! run the calculation

            icff=Node_Eval_Type(i_node_left,i_tree)

            !write(6,'(A,3(1x,I6))')'        RKBM: icff, i_node_left, i_tree = ', &
            !                                      icff, i_node_left, i_tree

            if( icff .eq. 0) then

                left_node_value=Tree_Evaluation(i_node_left,i_tree)
                left_node_value_string = &
                             trim( tree_evaluation_string(i_node_left,i_tree) )

            elseif( icff .lt. 0 .and. icff .ne. -9999) then

                if( iter .eq. 1) then
                    left_node_value=b_tmp(abs(icff))
                    !write(6,*) '        left_node_value ', left_node_value 
                    write(left_node_value_string, '(F7.3)') left_node_value
                    !write(6,*) '        left_node_value_string ', left_node_value_string
                else
                    left_node_value=btmp(abs(icff))
                    !write(6,*) '        left_node_value ', left_node_value 
                    write(left_node_value_string, '(F7.3)') left_node_value
                    !write(6,*) '        left_node_value_string ', left_node_value_string
                endif ! iter .eq. 1

            endif ! icff .eq. 0

            icff=node_eval_type(i_node_right,i_tree)

            if( icff .eq. 0) then

                right_node_value=Tree_Evaluation(i_node_right,i_tree)
                right_node_value_string = trim( tree_evaluation_string(i_node_right,i_tree) )

            elseif( icff .lt. 0 .and. icff .ne. -9999) then

                if( iter .eq. 1) then
                    right_node_value=b_tmp(abs(icff))
                    write(right_node_value_string, '(F7.3)') right_node_value
                else
                    right_node_value=btmp(abs(icff))
                    write(right_node_value_string, '(F7.3)') right_node_value
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


            !write(6,'(A,3(1x,I6))') 'RKBM: i_function, i_tree, node_eval_type(i_function,i_tree) ', &
            !                               i_function, i_tree, node_eval_type(i_function,i_tree)

            SELECT CASE( node_eval_type(i_function,i_tree) )


               CASE(1)  ! LHS + RHS

                 Tree_Evaluation(i_function,i_tree)=left_node_value + right_node_value

                 tree_evaluation_string(i_function,i_tree) =  &
                          trim( left_node_value_string ) // ' + ' // &
                          trim( right_node_value_string )

                 !write(6,'(A)')'        RKBM: call print4  case 1'
                 !call print4( i_time_step, icff, &
                 !             left_node_value,  left_node_value_string, &
                 !             right_node_value, right_node_value_string, &
                 !             tree_evaluation )

               CASE(2)  ! LHS - RHS

                 Tree_Evaluation(i_function,i_tree)=left_node_value-right_node_value

                 tree_evaluation_string(i_function,i_tree) =  &
                                 trim( left_node_value_string ) // ' - ' // &
                                 trim( right_node_value_string )

                 !write(6,'(A)')'        RKBM: call print4  case 2'
                 !call print4( i_time_step, icff, &
                 !             left_node_value,  left_node_value_string, &
                 !             right_node_value, right_node_value_string, &
                 !             tree_evaluation )

               CASE(3)  ! LHS * RHS

                 Tree_Evaluation(i_function,i_tree)=left_node_value*right_node_value

                 tree_evaluation_string(i_function,i_tree) =  &
                                 trim( left_node_value_string ) // ' * ' // &
                                 trim( right_node_value_string )

                 !write(6,'(A)')      '        RKBM: call print4  case 3'
                 !write(6,'(A,1x,A)') '        RKBM; left_node_value_string ', &
                 !                              trim(left_node_value_string)

                 !call print4( i_time_step, icff, &
                 !             left_node_value,  left_node_value_string, &
                 !             right_node_value, right_node_value_string, &
                 !             tree_evaluation )


               CASE(4)  ! protected: LHS/RHS


                 if( right_node_value .ne. 0.0D+0) then

                     Tree_Evaluation(i_function,i_tree)=left_node_value/right_node_value

                     tree_evaluation_string(i_function,i_tree) =  &
                                 trim( left_node_value_string) // ' / ' // &
                                 trim( right_node_value_string)

                 else

                     Tree_Evaluation(i_function,i_tree)=0.

                     tree_evaluation_string(i_function,i_tree) =  '0.0'

                 endif


               CASE(5)  ! '1.0D+0 - dexp(-dabs(LHS*RHS))'


                 cff=dabs(left_node_value*right_node_value)
                 Tree_Evaluation(i_function,i_tree)=1.0D+0-dexp(-1.0D+0*cff)

                 tree_evaluation_string(i_function,i_tree) = &
                       '( 1.0-exp(-1.0* abs( ' // trim( left_node_value_string )// &
                       ' * ' // trim( right_node_value_string) // ' ) )'

                 !write(6,'(A)')'        RKBM: call print4  case 5'

                 !call print4( i_time_step, icff, &
                 !             left_node_value,  left_node_value_string, &
                 !             right_node_value, right_node_value_string, &
                 !             tree_evaluation )


               CASE(6)  ! 'Michealis-Menton (abs(RHS) / (abs(LHS) + abs(RHS)))'


                 cff=abs(left_node_value)+abs(right_node_value)
                 if( cff .gt. 0.0D+0) then
                     Tree_Evaluation(i_function,i_tree)=abs(right_node_value)/cff

                     tree_evaluation_string(i_function,i_tree) = &
                       'abs( '   // trim( right_node_value_string)  // &
                        ' ) / ( ' // trim( left_node_value_string ) // &
                        ' + '     // trim( right_node_value_string ) // ' )'

                 else
                     Tree_Evaluation(i_function,i_tree)=0.0D+0

                     tree_evaluation_string(i_function,i_tree) = '0.0'

                 endif

                !write(6,'(A)')'        RKBM: call print4  case 6'

                !call print4( i_time_step, icff, &
                !                 left_node_value,  left_node_value_string, &
                !                 right_node_value, right_node_value_string, &
                !                 tree_evaluation )



               CASE(7)  ! 'abs(LHS*RHS)*(1.0D+0 - dexp(-dabs(LHS*RHS)))'


                 cff=dabs(left_node_value*right_node_value)

                 cff_string = 'abs( ' // trim( left_node_value_string ) // &
                     ' * ' // trim( right_node_value_string) // ' )'
                 cff_string = trim( cff_string )

                 Tree_Evaluation(i_function,i_tree)=cff*(1.0D+0-dexp(-1.0D+0*cff))

                 tree_evaluation_string(i_function,i_tree) = trim( cff_string ) // ' * ' //  &
                      '( 1.0-exp(-1.0* ' // trim( cff_string )  // ' ) )'

                 !write(6,'(8x, A, 1x,I2,3(1x,E15.7) )') 'icff, left, right, tree_eval ', &
                 !      icff, left_node_value, right_node_value, tree_evaluation(i_function,i_tree)

                 !write(6,'(8x, A, A )') 'left_node_value_string  ', &
                 !      trim( left_node_value_string )
                 !write(6,'(8x, A, A )') 'right_node_value_string ', &
                 !      trim( right_node_value_string )
                 !write(6,'(8x, A, A )') 'tree_evaluation_string  ', &
                 !      trim( tree_evaluation_string(i_function,i_tree) )

                 !write(6,'(A)')'        RKBM: call print4  case 7'

                 !call print4( i_time_step, icff, &
                 !             left_node_value,  left_node_value_string, &
                 !             right_node_value, right_node_value_string, &
                 !             tree_evaluation )

               CASE DEFAULT

                 write(*,*) 'wrong case number chosen in Runge Kutta evaluations'
                 stop 'RK bad case number'

            END SELECT

            Node_Eval_Type(i_function,i_tree)=0  ! now there is a value/variable

            node_eval_type(i_node_left,i_tree)  = 0
            node_eval_type(i_node_right,i_tree) = 0
            node_eval_type(i_function,i_tree)   = 0
            node_eval_type_string(i_node_left,i_tree)  = ' '
            node_eval_type_string(i_node_right,i_tree) = ' '
            node_eval_type_string(i_function,i_tree)   = ' '

          endif !  Node_Eval_Type(i_function,i_tree) .gt. 0

        enddo !  i_node

      enddo !  i_level

      Tree_Value(i_tree)=Tree_Evaluation(1,i_tree)

      tree_value_string(i_tree) = tree_evaluation_string(1,i_tree)

      !write(6,'(8x, A, 1x,E15.7)') &
      !      'RKBM: tree_eval(1, i_tree )  = ', &
      !       tree_evaluation(1, i_tree )
      write(6,'(A, 1x,I6,1x,A,5x,A)') &
            'i_tree, tree_eval(1, i_tree )  = ', &
             i_tree, ':',  trim( tree_evaluation_string( 1, i_tree ) )

    enddo !  i_tree

!    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!    !   Calculate the flow terms from the determined tree_value terms
!    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!    i_tree=0
!    do i_CODE_equation=0,n_CODE_equations   ! source of material   ! orig
!
!      do j_CODE_equation=0,n_CODE_equations ! sink of material   ! orig
!
!        if( i_CODE_equation .ne. j_CODE_equation) then
!
!          i_tree=i_tree+1
!
!          ! 'dabs' forces flow of material in one direction
!          bioflo(i_CODE_equation,j_CODE_equation)=dabs(tree_value(i_tree))
!
!          write(6,'(A,3(1x,I6),2(2x, E15.7))') &
!          'RKBM:1 i_tree, i_code_equation, j_code_equation, &
!          &tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation) ', &
!                  i_tree, i_code_equation, j_code_equation, &
!           tree_value(i_tree), bioflo(i_CODE_equation,j_CODE_equation)
!
!        else
!
!          bioflo(i_CODE_equation,j_CODE_equation)=0.0D+0   ! never flow to/from same component
!
!        endif !   i_CODE_equation .ne. j_CODE_equation
!
!      enddo !  j_CODE_equation
!
!    enddo !  i_CODE_equation
!
!    do i_CODE_equation=1,n_CODE_equations   ! source of material
!      fbio(i_CODE_equation)=0.0D+0
!    enddo ! i_CODE_equation
!
!    !   bring in the component flow sources and sinks
!
!    do i_CODE_equation=0,n_CODE_equations   ! source of material
!
!      do j_CODE_equation=0,n_CODE_equations ! sink of material
!
!        if( i_CODE_equation .gt. 0) then
!          fbio(i_CODE_equation)=fbio(i_CODE_equation)-bioflo(i_CODE_equation,j_CODE_equation)
!        endif !   i_CODE_equation .gt. 0
!
!        if( j_CODE_equation .gt. 0) then
!          fbio(j_CODE_equation)=fbio(j_CODE_equation)+bioflo(i_CODE_equation,j_CODE_equation)
!        endif !   j_CODE_equation .gt. 0
!
!      enddo !  j_CODE_equation
!    enddo !  i_CODE_equation
!
!
!    do i_CODE_equation=1,n_CODE_equations
!
!      kval(iter,i_CODE_equation)=dt*fbio(i_CODE_equation)
!
!      if( iter .eq. 1) then
!
!        btmp(i_CODE_equation)=b_tmp(i_CODE_equation)+(kval(iter,i_CODE_equation)/2.0D+0)
!
!      elseif( iter .eq. 2) then
!
!        btmp(i_CODE_equation)=b_tmp(i_CODE_equation)+(kval(iter,i_CODE_equation)/2.0D+0)
!
!      elseif( iter .eq. 3) then
!
!        btmp(i_CODE_equation)=b_tmp(i_CODE_equation)+kval(iter,i_CODE_equation)
!
!      elseif( iter .eq. 4) then
!
!        cff = (kval(1,i_CODE_equation)/6.0D+0) + (kval(2,i_CODE_equation)/3.0D+0) + &
!              (kval(3,i_CODE_equation)/3.0D+0) + (kval(4,i_CODE_equation)/6.0D+0)
!        b_tmp(i_CODE_equation)=b_tmp(i_CODE_equation)+cff
!
!      endif !   iter .eq. 1
!
!    enddo !  i_CODE_equation
!
!  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
    !exit
  enddo ! iter
!
!
!
!  ! if b_tmp is bad on any time step, then return with a bad result
!
!  if( any( isnan( b_tmp ) ) .or.  any( abs(b_tmp)  > 1.0d20 ) ) then
!
!      L_bad_result = .TRUE.
!      write(6,'(A,1(1x,I6),12(1x,E15.7))') &
!            'RKBM: bad result i_time_step, b_tmp ', &
!                              i_time_step, b_tmp(1:n_CODE_equations)
!      return
!
!  endif !   any( isnan( b_tmp ) ) .or.  any( abs(b_tmp) > 1.0d20
!
!  !---------------------------------------------------------------------------
!
!  !   Array Assignment with no negative values
!
!  !Runge_Kutta_Solution(i_time_step,1:n_CODE_equations) = max(b_tmp,0.0d+0)
!
!  ! NPZ model
!
!  Runge_Kutta_Solution(i_time_step,1:n_CODE_equations) = max( b_tmp, 0.0d0 )
!
!  !---------------------------------------------------------------------------
!
!  write(6,'(A,1(1x,I6),12(1x,E15.7))') &
!        'RKBM: i_time_step, RK_Soln ', &
!               i_time_step, Runge_Kutta_Solution(i_time_step,1:n_CODE_equations)
!
  !exit

enddo ! i_time_step

!write(6,*)'RKBM: deallocate node_eval_type_string '
deallocate( node_eval_type_string ) 

!write(6,*)'RKBM: deallocate tree_value_string '
deallocate( tree_value_string ) 

return
end subroutine Runge_Kutta_Box_Model
!234567890123456789012345678901234567890123456789012345678901234567890
