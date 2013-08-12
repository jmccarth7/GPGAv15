subroutine select_best_RK_lmdif_result( &
                i_GA_best_parent, parent_parameters, &
                child_parameters, &
                individual_quality, L_stop_run )

! written by: Dr. John R. Moisan [NASA/GSFC] 5 December, 2012
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum parameter set for a coupled set of equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module
use clock_module

use GP_parameters_module
use GA_parameters_module
use GP_variables_module
use GA_variables_module
use GP_data_module


implicit none

integer :: child_number

integer ::  start_individ
integer ::  stop_individ
integer ::  start_limit
integer ::  stop_limit
integer ::  isource
integer ::  iproc
integer ::  chunk
integer ::  message_len
integer ::  numsent          
integer ::  sender           
integer ::  nsafe
integer ::  i_dummy
integer ::  i_individual
integer ::  i_2_individual

integer,parameter ::  itag  = 1
integer,parameter ::  itag2 = 2
integer,parameter ::  itag3 = 3


real(kind=8),&
 dimension(n_maximum_number_parameters,n_GA_individuals) ::  parent_parameters
real(kind=8),&
 dimension(n_maximum_number_parameters,n_GA_individuals) ::  child_parameters


real(kind=8) :: individual_SSE_best_1
real(kind=8) :: individual_ranked_fitness_best_1
real(kind=8) :: Individual_Fitness_best_1

real(kind=8),dimension(n_maximum_number_parameters) :: parent_parameters_best_1



real(kind=8) :: buffer(n_maximum_number_parameters + 2)
real(kind=8) :: buffer_recv(n_maximum_number_parameters + 2)


integer (kind=4) ::      i
integer (kind=4) :: i_GA_Best_Parent
integer (kind=4) :: i_GA_Best_Parent_1

integer (kind=4) :: i_GA_generation_last

real(kind=8), parameter :: tol = 1.0d-30

real(kind=8),parameter :: zero = 0.0d0


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4),intent(in) :: individual_quality(n_GA_individuals)

real(kind=8), external :: indiv_fitness

logical :: L_stop_run

logical :: Lprint_GA

integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: jj
integer(kind=4) :: i_parameter

!----------------------------------------------------------------------


n_parameters = n_GP_parameters



if( myid == 0 )then
    if( L_ga_print )then 
    write(GA_print_unit,'(//A)') 'sbrl: at entry  '
    write(GA_print_unit,'(A,1x,E15.7)') 'sbrl: dt ', dt
    write(GA_print_unit,'(A,1x,I10)') &
          'sbrl: n_parameters    =   ', n_parameters
    write(GA_print_unit,'(A,1x,I10)') &
          'sbrl: n_GP_parameters =   ', n_GP_parameters
    endif ! L_ga_print  

endif ! myid == 0




! save best parent parameters from the RK process, 
! then run lmdif to try to improve the best parent


i_GA_best_parent_1                                         = i_GA_best_parent
individual_SSE_best_1                       = individual_SSE(i_GA_best_parent)
individual_ranked_fitness_best_1 = individual_ranked_fitness(i_GA_best_parent)

Individual_Fitness        = Individual_Ranked_Fitness(i_GA_Best_Parent)
Individual_Fitness_best_1 = Individual_Fitness

parent_parameters_best_1(1:n_parameters) =  &
                        Parent_Parameters(1:n_parameters, i_GA_Best_Parent)

if( L_ga_print )then                                                                                
write(GA_print_unit,'(/A,1x,I10, 2(1x,E24.16))') &
      'sbrl: i_GA_best_parent_1, individual_SSE_best_1, &
                     &individual_ranked_fitness_best_1', &
             i_GA_best_parent_1, individual_SSE_best_1, &
                      individual_ranked_fitness_best_1
write(GA_print_unit,'(/A,1x,E24.16/)') &
      'sbrl: individual_fitness_best_1 ', individual_fitness_best_1

!write(GA_print_unit,'(/A/1x,I6, 12(1x,E15.7))') &
!      'sbrl: i_GA_best_parent_1, parent_parameters_best_1(1:n_parameters) ', &
!             i_GA_best_parent_1, parent_parameters_best_1(1:n_parameters)

write(GA_print_unit,'(/A)') 'sbrl: i, parent_parameters_best_1(i) '
do  i = 1, n_parameters
    write(GA_print_unit,'(I6, 1x,E24.16 )') i, parent_parameters_best_1(i)
enddo

endif ! L_ga_print  

!-------------------------------------------------------------------------------

!  run lmdif on best individual from the RK process

if( L_ga_print )then                                                                                
if( L_stop_run )then
    write(GA_print_unit,'(A,2(1x,I6))') &
      'sbrl: i_GA_Generation_last, i_GA_Best_Parent  call setup_run_lmdif ', &
             i_GA_Generation_last, i_GA_Best_Parent
else
    write(GA_print_unit,'(A,2(1x,I6))') &
      'sbrl: n_GA_Generations, i_GA_Best_Parent  call setup_run_lmdif ', &
             n_GA_Generations, i_GA_Best_Parent
endif ! L_stop_run
endif ! L_ga_print  


call setup_run_lmdif( i_GA_Best_Parent, parent_parameters, individual_quality )


if( L_ga_print )then                                                                                
if( L_stop_run )then
    write(GA_print_unit,'(A,2(1x,I6))') &
      'sbrl: i_GA_Generation_last, i_GA_Best_Parent  aft call setup_run_lmdif ', &
                  i_GA_Generation_last, i_GA_Best_Parent
else
    write(GA_print_unit,'(A,2(1x,I6))') &
      'sbrl: n_GA_Generations, i_GA_Best_Parent  aft call setup_run_lmdif ', &
                  n_GA_Generations, i_GA_Best_Parent
endif ! L_stop_run
endif ! L_ga_print  


!-------------------------------------------------------------------------------


!  compute fitness for parameters of the best parent after lmdif has been run


if( L_ga_print )then                                                                                
write(GA_print_unit,'(/A)')  &
      'sbrl: i_GA_ind, ind_SSE, ind_ranked_fitness  ----  aft lmdif '


write(GA_print_unit,'(/A,1x,I6)') &
      'sbrl: lmdif i_GA_best_parent ', i_GA_best_parent
write(GA_print_unit,'(/A)') &
      'sbrl: lmdif i, parent_parameters(i, i_GA_best_parent) '

do  i = 1, n_parameters
    write(GA_print_unit,'(I6, 1x,E24.16 )') &
          i, parent_parameters(i, i_GA_best_parent)
enddo ! i



write(GA_print_unit,'(A,1x,I6,1x,E24.16)') &
      'sbrl: i_GA_best_parent, individual_SSE', &
             i_GA_best_parent, individual_SSE(i_GA_best_parent)

endif ! L_ga_print  


individual_ranked_fitness(i_GA_best_parent) = indiv_fitness( i_GA_best_parent )


Individual_Fitness = Individual_Ranked_Fitness(i_GA_Best_Parent)


if( L_ga_print )then                                                                                
write(GA_print_unit,'(/A,1x,I6, 2(1x,E24.16))') &
      'sbrl: lmdif i_GA_best_parent, &
      &individual_SSE, individual_ranked_fitness', &
                   i_GA_best_parent, &
                   individual_SSE(i_GA_best_parent), &
                   individual_ranked_fitness(i_GA_best_parent)

write(GA_print_unit,'(/A,1x,E24.16/)') &
      'sbrl: lmdif individual_fitness ', individual_fitness
endif ! L_ga_print  


!------------------------------------------------------------------------------

!  test if lmdif has improved the best parent parameters


!  compare the fitness of the parameter set from the RK integrations
!  with    the fitness of the parameter set returned by lmdif 


!  select the set of parameters with the best fitness


if( L_ga_print )then                                                                                
write(GA_print_unit,'(/A, 1x,E24.16)') &
      'sbrl: fcn   individual_ranked_fitness_best_1            ', &
                   individual_ranked_fitness_best_1
write(GA_print_unit,'(A, 1x,E24.16/)') &
      'sbrl: lmdif individual_ranked_fitness(i_GA_best_parent) ', &
                   individual_ranked_fitness(i_GA_best_parent)
endif ! L_ga_print  


if( individual_ranked_fitness(i_GA_best_parent) <= &
                        individual_ranked_fitness_best_1 )then


    ! the fitness of the RK process output is better


    if( L_ga_print )then
    write(GA_print_unit,'(/A/)')&
          'sbrl:  the fitness of the RK process output is better '
    endif ! L_ga_print  

    individual_fitness         = individual_ranked_fitness_best_1
    Individual_SSE_best_parent = individual_SSE_best_1

!new 20130713
    child_parameters(1:n_parameters,i_GA_Best_Parent) =  &
                        parent_parameters_best_1(1:n_parameters)
!new 20130713

    ! choose the parameters of the best parent from the RK fcn integration

    if( L_ga_print )then
    write(GA_print_unit,'(/A)')&
          'sbrl: set the GA-optimized initial condition array '

    write(GA_print_unit,'(/A/1x,I6, 6(1x,E24.16))') &
          'sbrl: i_GA_best_parent_1, parent_parameters_best_1(1:n_CODE_Equations) ', &
                 i_GA_best_parent_1, &
                 ( parent_parameters_best_1(jj), jj = 1, n_CODE_Equations )
    endif ! L_ga_print  

    GP_Individual_Initial_Conditions(1:n_CODE_Equations) = &
            parent_parameters_best_1(1:n_CODE_Equations)

    if( L_ga_print )then  
    write(GA_print_unit,'(/A/ 6(1x,E24.16))') &
          'sbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
               ( GP_Individual_Initial_Conditions(jj), jj = 1,n_CODE_Equations )
    endif ! L_ga_print  

    if( L_stop_run )then
        write( GA_output_unit, '(I6,1x,I6, 12(1x,E15.7))') &
          i_GA_Generation_last, i_GA_best_parent_1, &
          individual_ranked_fitness_best_1, &
          (parent_parameters_best_1(jj),jj = 1,n_parameters)
    else
        write( GA_output_unit, '(I6,1x,I6, 12(1x,E15.7))') &
          n_GA_Generations, i_GA_best_parent_1, &
          individual_ranked_fitness_best_1, &
          (parent_parameters_best_1(jj),jj = 1,n_parameters)
    endif ! L_stop_run


    !-------------------------------------------------------------------------------


    ! load the parameters from the RK process into GP_Individual_Node_Parameters


    if( L_ga_print )then  
    write(GA_print_unit,'(/a/)') &
          'sbrl: set the GA-optimized CODE parameter array'
    endif ! L_ga_print  

    i_parameter = n_CODE_equations ! start at this number because of the
                                   ! initial conditions (n_CODE_Equations of them)

    do i_tree=1,n_trees
      do i_node=1,n_nodes

        if( GP_individual_node_type(i_node,i_tree) .eq. 0 ) then  ! there is a set parameter

            !if( L_ga_print )then 
            !write(GA_print_unit,'(A,3(1x,I6))') &
            !      'sbrl:1 i_tree, i_node,  GP_individual_node_type(i_node,i_tree) ', &
            !              i_tree, i_node,  GP_individual_node_type(i_node,i_tree)
            !endif ! L_ga_print  

            i_parameter=i_parameter+1

            GP_Individual_Node_Parameters(i_node,i_tree) = &
                          parent_parameters_best_1( i_parameter )

            !if( L_ga_print )then            
            !write(GA_print_unit,'(A,1x,I6,1x,E20.10)') &
            !      'sbrl:1 i_Parameter, parent_parameters_best_1( i_parameter ) ', &
            !              i_Parameter, parent_parameters_best_1( i_parameter )

            !write(GA_print_unit,'(A,2(1x,I6),1x,E20.10)') &
            !      'sbrl:1 i_tree, i_node, GP_indiv_node_params', &
            !              i_tree, i_node, GP_individual_node_parameters(i_node,i_tree)
            !endif ! L_ga_print  

        endif !   GP_individual_node_type(i_node,i_tree) .eq. 0

     enddo ! i_node

    enddo ! i_tree


!--------------------------------------------------------------------------------------


else  ! lmdif is best


    ! the fitness of the lmdif output is better


    if( L_ga_print )then  
    write(GA_print_unit,'(/A/)')&
          'sbrl:  the fitness of the lmdif output is better '
    endif ! L_ga_print  


    individual_fitness         = individual_ranked_fitness(i_GA_best_parent)
    Individual_SSE_best_parent = individual_SSE(i_GA_best_parent)

!new 20130713
    child_parameters(1:n_parameters,i_GA_Best_Parent) =  &
                        Parent_Parameters(1:n_parameters, i_GA_Best_Parent)
!new 20130713

    ! choose the parameters from the lmdif output for the best parent  

    if( L_ga_print )then   
    write(GA_print_unit,'(/A,1x,I6, 12(1x,E15.7))') &
          'sbrl: i_GA_best_parent, Parent_Parameters ', &
                 i_GA_best_parent, &
                 (Parent_Parameters(jj, i_GA_Best_Parent),jj= 1,n_parameters)
    endif ! L_ga_print  

    GP_Individual_Initial_Conditions(1:n_CODE_Equations) = &
                    Parent_Parameters( 1:n_CODE_Equations, i_GA_Best_Parent )

    if( L_ga_print )then    
    write(GA_print_unit,'(/A/ 6(1x,E24.16))') &
          'sbrl: GP_Individual_Initial_Conditions(1:n_CODE_Equations) ', &
                (GP_Individual_Initial_Conditions(jj),jj=1,n_CODE_Equations)
    endif ! L_ga_print  


    if( L_stop_run )then
        write( GA_output_unit , '(I6,1x,I6, 12(1x,E15.7))') &
          i_GA_Generation_last, i_GA_best_parent, &
          individual_ranked_fitness(i_GA_best_parent), &
          (parent_parameters(jj, i_GA_best_parent), jj=1,n_parameters)
    else
        write( GA_output_unit , '(I6,1x,I6, 12(1x,E15.7))') &
          n_GA_Generations, i_GA_best_parent, &
          individual_ranked_fitness(i_GA_best_parent), &
          (parent_parameters(jj, i_GA_best_parent), jj=1,n_parameters)
    endif ! L_stop_run



    !--------------------------------------------------------------------------------------

    ! load the parameters output by lmdif into GP_Individual_Node_Parameters


    if( L_ga_print )then      
    write(GA_print_unit,'(/a/)')&
          'sbrl: set the GA-optimized CODE parameter array'
    endif ! L_ga_print  

    i_parameter = n_CODE_equations ! start at this number because of the
                                   ! initial conditions (n_CODE_Equations of them)

    do i_tree=1,n_trees
      do i_node=1,n_nodes


        if( GP_individual_node_type(i_node,i_tree) .eq. 0 ) then  ! there is a set parameter

            !if( L_ga_print )then   
            !write(GA_print_unit,'(A,4(1x,I6))') &
            !  'sbrl:2 myid, i_tree, i_node,  GP_individual_node_type(i_node,i_tree) ', &
            !          myid, i_tree, i_node,  GP_individual_node_type(i_node,i_tree)
            !endif ! L_ga_print  

            i_parameter=i_parameter+1

            GP_Individual_Node_Parameters(i_node,i_tree) = &
                        Parent_Parameters(i_Parameter, i_GA_Best_Parent)

            !if( L_ga_print )then              
            !write(GA_print_unit,'(A,2(1x,I6),1x,E20.10)') &
            !      'sbrl:2 i_GA_Best_Parent,i_Parameter, &
            !             &Parent_Parameters(i_Parameter,i_GA_Best_Parent) ', &
            !              i_GA_Best_Parent,i_Parameter, &
            !              Parent_Parameters(i_Parameter, i_GA_Best_Parent)

            !write(GA_print_unit,'(A,2(1x,I6),1x,E20.10)') &
            !      'sbrl:2 i_tree, i_node, GP_indiv_node_params', &
            !              i_tree, i_node, GP_individual_node_parameters(i_node,i_tree)
            !endif ! L_ga_print  

        endif ! GP_individual_node_type(i_node,i_tree) .eq. 0

      enddo ! i_node

    enddo ! i_tree



endif ! individual_ranked_fitness...



return


end subroutine select_best_RK_lmdif_result