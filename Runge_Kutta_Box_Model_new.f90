subroutine Runge_Kutta_Box_Model( )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! carry out a prescribed Runge-Kutta numerical integration
! using the GP architecture to solve a coupled system of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!use GP_Parameters
!use GA_Parameters
!use GP_Variables

use GP_model_parameters_module                                                                              
use GP_Parameters_module                                                                                    
use GA_Parameters_module                                                                                    
use GP_variables_module                                                                                     
use GA_Variables_module                                                                                     
use GP_Data_module                                                                                          
use GP_variables_module

!use Model_Parameters

! Runge-Kutta specific work arrays

implicit none

!real(kind=8) :: Tree_Value(n_Trees)

!real(kind=8) :: bioflo(0:n_CODE_Equations,0:n_CODE_Equations) !Forcing function are not part of the bio-flow

!!real(kind=8) :: Runge_Kutta_Time_Step(4) = (/ 0.0D+0, 0.5D+0, 0.5D+0, 1.0D+0 /)

!Forcing functions are used in computations, so are included here for book keeping purposes
!real(kind=8) :: b_tmp(n_Variables),fbio(n_Variables)

!real(kind=8) :: kval(4,n_Variables)
real(kind=8) :: cff

integer(kind=4) :: iter
integer(kind=4) :: i_Time_Step, i_Track, i_Tree
integer(kind=4) :: i_CODE_Equation, j_CODE_Equation, i_Variable


!real(kind=8) :: left_node_value,right_node_value
!real(kind=8) :: cff

!integer(kind=4) :: icff
!integer(kind=4) :: iter



!--------------------------------------------------------------------------------------


write(6,'(/A/)') 'rkbm: entry Runge_Kutta_Box_Model '

!write(6,'(A,10(1x,E15.7)/ )') &
!      'rkbm: before loop  btmp(:)', btmp(:)                                 


! start the time stepping loop
do  i_Time_Step = 1, n_Time_Steps
    
    if( mod(i_Time_Step, 100) .eq. 0 ) then
        !write (6,'(A,I5,A,I5,A)') 'Time step'//char(9), i_Time_Step, &
        !            char(9)//'of'//char(9), n_Time_Steps, char(13)
        write (6,'(A,I5,A,I5)') 'Time step'//char(9), i_Time_Step, &
                    char(9)//'of'//char(9), n_Time_Steps
    endif
    
    b_tmp(:) = Numerical_CODE_Solution(i_Time_Step-1,:)  ! Array Assignment

    !write(6,'(A,1x,I6,10(1x,E15.7) )') &
    !      'rkbm: i_time_step, Num_CODE_Sol(i_Time_Step-1,:) ', &
    !               i_time_step, Numerical_CODE_Solution(i_Time_Step-1,:) 
    !write(6,'(A,1x,I6,10(1x,E15.7) )') &
    !      'rkbm: i_time_step, b_tmp(:)                      ', &
    !               i_time_step, b_tmp(:)                                 

    btmp = b_tmp

    !write(6,'(A,1x,I6,10(1x,E15.7) )') &
    !      'rkbm: i_time_step, btmp(:)                       ', &
    !               i_time_step, btmp(:)                                 

    ! carry out a Runge-Kutta time step
    do  iter=1,4
        
        write(6,'(/A,2(1x,I6)/)') 'rkbm: i_time_step, iter ', &
                                         i_time_step, iter 

        ! Call forcing functions for the box model
        !call DoForcing(btmp, Runge_Kutta_Time_Step(iter), i_Time_Step)
        
        fbio = 0.0D+0
        do  i_Track = 1,n_Tracked_Resources

            !!!! Call Model_Diagnostics()

            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Evaluate the trees 
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            Tree_Value=0.0D+0                                 ! Matrix Assignment

            do  i_Tree=1,n_Trees

                write(6,'(A,1x,I6,5x,L1)') &
                      'rkbm: i_tree, associated(GP_Trees(i_Tree,i_Track)%n)  ', &
                             i_tree, associated(GP_Trees(i_Tree,i_Track)%n)  

                if( associated(GP_Trees(i_Tree,i_Track)%n) ) then
                    Tree_Value(i_Tree) = GP_Trees(i_Tree,i_Track)%n%val()
                    write(6,'(A,22x,I6,1x,E15.7 )') &
                          'rkbm: i_tree, Tree_Value(i_tree) ', &
                                 i_tree, Tree_Value(i_tree)                      
                endif ! associated(GP_Trees...

            enddo ! i_Trees


            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !   Calculate the flow terms from the determined tree_value terms
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            !write(6,'(/A/)') 'rkbm: Calculate the flow terms from the determined tree_value terms'

            i_Tree=0
            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material
                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material

                    if( i_CODE_Equation .ne. j_CODE_Equation) then

                        i_Tree=i_Tree+1

                        ! 'abs' forces flow of material in one direction
                        bioflo(i_CODE_Equation,j_CODE_Equation)=abs(Tree_Value(i_Tree))  
                        
                        !write(*,*) 'Tree: ', i_Tree, ' - ', Tree_Value(i_Tree)
                        !write(,'(A,1x,I6,2x,E20.10)') 'rkbm: i_Tree, Tree_Value(i_Tree)', &
                        !                                       i_Tree, Tree_Value(i_Tree)

                    else

                        ! never flow to/from same component
                        bioflo(i_CODE_Equation,j_CODE_Equation)=0.0D+0   

                    endif ! i_CODE_Equation .ne. j_CODE_Equation

                enddo ! j_CODE_Equation
            enddo ! i_CODE_Equation

            
            ! bring in the component flow sources and sinks
            do  i_CODE_Equation=0,n_CODE_Equations   ! source of material

                do  j_CODE_Equation=0,n_CODE_Equations ! sink of material
  
                    if( i_CODE_Equation .gt. 0 ) then
    
                        if( bioflo_map(i_CODE_Equation,i_Track) .gt. 0 ) then
                            
                            fbio(bioflo_map(i_CODE_Equation,i_Track)) = &
                                fbio(bioflo_map(i_CODE_Equation,i_Track)) -  &
                                         bioflo(i_CODE_Equation,j_CODE_Equation)
                              
                        endif ! bioflo_map(i_CODE_Equation,i_Track) .gt. 0 
    
                    endif ! i_CODE_Equation .gt. 0 
    
                    if( j_CODE_Equation .gt. 0 ) then

                        if( bioflo_map(j_CODE_Equation,i_Track) .gt. 0 ) then
                            
                            fbio(bioflo_map(j_CODE_Equation,i_Track)) = &
                                 fbio(bioflo_map(j_CODE_Equation,i_Track)) + &
                                          bioflo(i_CODE_Equation,j_CODE_Equation)
                              
                        endif ! bioflo_map(j_CODE_Equation,i_Track) .gt. 0 

                    endif ! j_CODE_Equation .gt. 0 
  
                enddo ! j_CODE_Equation
            enddo ! i_CODE_Equation

        enddo ! End Tracked Resources loop

        
        ! Capture any export terms, or boundary conditions after bio flow is calculated
        ! If the mode does not contain any of these, SecondaryForcing() should do nothing

        !!!!call SecondaryForcing(fbio)

        
        do  i_Variable=1,n_Variables

            kval(iter,i_Variable)=Delta_Time_in_Days*fbio(i_Variable)

            if( iter .eq. 1) then

                btmp(i_Variable)=b_tmp(i_Variable)+(kval(iter,i_Variable)/2.0D+0)

            elseif( iter .eq. 2) then

                btmp(i_Variable)=b_tmp(i_Variable)+(kval(iter,i_Variable)/2.0D+0)

            elseif( iter .eq. 3) then

                btmp(i_Variable)=b_tmp(i_Variable)+kval(iter,i_Variable)

            elseif( iter .eq. 4) then

                cff=(kval(1,i_Variable)/6.0D+0)+(kval(2,i_Variable)/3.0D+0)+ & 
                    (kval(3,i_Variable)/3.0D+0)+(kval(4,i_Variable)/6.0D+0)
                b_tmp(i_Variable)=b_tmp(i_Variable)+cff

            endif

        enddo ! End Kval loop
        
    enddo ! End iter loop
     
    Numerical_CODE_Solution(i_Time_Step,1:n_Variables)=max(b_tmp(1:n_Variables),0.0D+0)  


    write(6,'(A,1x,I6,1x,6(1x,E15.7))') 'rkbm: i_time_step, solution ', &
                   i_Time_Step, Numerical_CODE_Solution(i_Time_Step,1:n_Variables)

                                                                                                                                   
    !write(rk_output_unit,'(I10,1x,8(1x,E18.8E3))') & 
    !    i_time_step, Numerical_CODE_Solution(i_time_step,1:n_CODE_equations)
                                                                                                                                   

    
enddo ! End Time step loop

write(6,'(A)') 'rkbm: leave Runge_Kutta_Box_Model '


return

end subroutine Runge_Kutta_Box_Model
