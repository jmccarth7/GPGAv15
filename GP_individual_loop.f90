subroutine GP_individual_loop( new_group, new_comm, i_GP_generation )

! program written by: Dr. John R. Moisan [NASA/GSFC] 31 January, 2013

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! program to use a twin experiment to test the effectiveness of
! a finding the optimum equation and parameter sets for a system of
! coupled ordinary differential equations
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module

use GP_Parameters_module
use GP_variables_module
use GA_Parameters_module
use GA_Variables_module
use GP_Data_module

use fasham_variables_module
use Tree_Node_Factory_module
use class_Tree_Node


implicit none



integer(kind=4),intent(in) :: new_group
integer(kind=4),intent(in) :: new_comm
integer(kind=4),intent(in) :: i_GP_Generation


integer :: i

integer(kind=4) :: i_GP_individual
integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: ii


integer(kind=4) :: n_GP_vars

integer(kind=4) :: n_procs

integer(kind=4) :: i_part
integer(kind=4) :: i_gp_1
integer(kind=4) :: i_gp_2
integer(kind=4) :: ind1
integer(kind=4) :: ind2
integer(kind=4) :: n_indiv

integer(kind=4),parameter :: tag_ind_sse = 200000
integer(kind=4),parameter :: tag_ind_fit = 100000

integer(kind=4) :: tag_fit_r
integer(kind=4) :: tag_fit_s
integer(kind=4) :: tag_sse_r
integer(kind=4) :: tag_sse_s


real(kind=8) ::   individual_fit_rec
real(kind=8) ::   individual_SSE_rec

real(kind=8),allocatable,dimension(:)  ::   fit_buffer_send
real(kind=8),allocatable,dimension(:)  ::   sse_buffer_send
real(kind=8),allocatable,dimension(:)  ::   sse_buffer_recv

!---------------------------------------------------------------------------------------

GP_population_fitness   = 0.0d0
GP_Child_Individual_SSE = 1.0d13
GP_Adult_Individual_SSE = 1.0d13
GP_Adult_Population_SSE = 1.0d13



allocate( fit_buffer_send( n_GP_individuals ) )
allocate( sse_buffer_send( n_GP_individuals ) )
allocate( sse_buffer_recv( n_GP_individuals ) )



divider = ( numprocs - 1 ) / n_partitions

call mpi_comm_rank( new_comm, new_rank, ierr )
call mpi_comm_size( new_comm, n_procs,  ierr )

if( myid == 0 )then
!if( new_rank == 0 )then
    write(GP_print_unit,'(/A,4(1x,i3))')&
     'gil: before loop myid, new_rank, n_code_equations,  n_GP_params ',&
                       myid, new_rank, n_code_equations,  n_GP_parameters
    write(GP_print_unit,'(A,5(1x,i3))')&
     'gil: before loop myid, new_rank, numprocs, divider, n_partitions', &
                       myid, new_rank, numprocs, divider, n_partitions
    write(GP_print_unit,'(A,3(1x,i3))')&
     'gil: before loop myid, new_rank, n_procs                        ', &
                       myid, new_rank, n_procs
endif !  myid == 0


! do the loop over the GP individuals in n_partitions chunks

part_loop:&
do  i_part = 1,  n_partitions

    ! i_gp_1 and i_gp_2 are limits on the processor number for this partition

    i_gp_1 = divider * (i_part - 1) + 1
    i_gp_2 = divider *  i_part

    !---------------------------------------------------------------------------------
    ! ind1 and ind2 are limits on the i_GP_individuals processed in this partition

    ind1 =  (n_GP_individuals / n_partitions)  * (i_part-1) + 1
    ind2 =  (n_GP_individuals / n_partitions)  *  i_part

    ! get any remaining individuals in the last partition

    if( i_part == n_partitions )then
        ind2 = n_GP_individuals
    endif ! i_part == n_partitions

    !write(GP_print_unit,'(A,5(1x,i3))')&
    !    'gil:in loop myid, new_rank, i_part, ind1, ind2     ',&
    !                 myid, new_rank, i_part, ind1, ind2


    !---------------------------------------------------------------------------------

    if( myid == 0 ) then


        ! receive the fitness information 

        tag_fit_r =  tag_ind_fit 

        write(GP_print_unit,'(A,5(1x,I7)/)')&
              'gil:31t myid, new_rank, tag_fit_r, ind1, ind2', &
                       myid, new_rank, tag_fit_r, ind1, ind2     

        n_indiv = ind2 - ind1 + 1

        call MPI_RECV( GP_population_fitness(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_fit_r,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )



        ! receive the SSE information 

        tag_sse_r = tag_ind_sse   

        !write(GP_print_unit,'(A,4(1x,I7)/)')&
        !      'gil:33t myid, new_rank, tag_sse_r, i_GP_individual', &
        !               myid, new_rank, tag_sse_r, i_GP_individual

        n_indiv = ind2 - ind1 + 1
        call MPI_RECV( sse_buffer_recv(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_sse_r,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        !write(GP_print_unit,'(A,3(1x,i3),1x,E15.7)')&
        ! 'gil:in loop RECV myid, new_rank, i_GP_individual, individual_SSE_rec ', &
        !                   myid, new_rank, i_GP_individual, individual_SSE_rec


        GP_Child_Individual_SSE(ind1:ind2) =  sse_buffer_recv(ind1:ind2)
        GP_Adult_Individual_SSE(ind1:ind2) =  sse_buffer_recv(ind1:ind2)
        GP_Adult_Population_SSE(ind1:ind2) =  sse_buffer_recv(ind1:ind2)



    elseif( i_gp_1 <= myid  .and.   &
                      myid   <= i_gp_2         )then

        gp_ind_loop:&
        do  i_GP_individual= ind1, ind2    ! 1,n_GP_individuals

            !write(GP_print_unit,'(A,8(1x,I4)/)')&
            !      'gil: myid, new_rank, i_part, i_gp_1, i_gp_2, ind1, ind2, i_GP_individual', &
            !            myid, new_rank, i_part, i_gp_1, i_gp_2, ind1, ind2, i_GP_individual


            !------------------------------------------------------------------------------

            ! calculate how many parameters total to fit for the specific individual CODE
            ! and save this number in GP_Individual_N_GP_param(i_GP_individual)

            n_GP_Parameters = n_code_equations


            do  i_Tree=1,n_Trees
                do  i_Node=1,n_Nodes
                    if( GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual) .eq. 0) then
                        n_GP_Parameters = n_GP_Parameters+1
                    endif ! GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)
                enddo ! i_node
            enddo ! i_tree

            if( new_rank == 0 )then
                write(GP_print_unit,'(A,2(1x,i3))')&
                 'gil: after loop n_code_equations,  n_GP_params',&
                                  n_code_equations,  n_GP_parameters
                !flush(GP_print_unit)
            endif !  new_rank == 0

            !------------------------------------------------------------------------

            GP_Individual_N_GP_param(i_GP_individual) = n_GP_parameters

            !------------------------------------------------------------------------------


            ! run GPCODE_... to evaluate this individual  if Run_GP_Calculate_Fitness is true


            if( Run_GP_Calculate_Fitness(i_GP_Individual) ) then

                !-----------------------------------------------------------------------------------

                ! these get set randomly in the GA-lmdif search algorithm ( in GPCODE* )
                !GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) = 0.0d0               ! 20131209

                !-----------------------------------------------------------------------------------

                if( new_rank == 0 )then
                    write(GP_print_unit,'(/A)')&
                    'gil:----------------------------------------------------------------------'
                    write(GP_print_unit,'(/A,2(1x,i3),3x,L1)') &
                          'gil: i_GP_Gen, i_GP_indiv, Run_GP_Calculate_Fitness', &
                                i_GP_Generation, i_GP_individual, &
                                       Run_GP_Calculate_Fitness(i_GP_Individual)
                    write(GP_print_unit,'(A)')&
                    'gil:----------------------------------------------------------------------'
                    !!flush(GP_print_unit)
                endif !  new_rank == 0

                !-----------------------------------------------------------------------------------

                if( new_rank == 0 )then
                    write(GP_print_unit,'(/A,4(1x,i3))') &
                      'gil: i_GP_individual, n_trees, n_nodes, n_GP_parameters ', &
                            i_GP_individual, n_trees, n_nodes, n_GP_parameters
                    !flush(GP_print_unit)
                endif !  new_rank == 0

                !-----------------------------------------------------------------------------------

                do  i_Tree=1,n_Trees
                    do  i_Node=1,n_Nodes

                        GP_Individual_Node_Type(i_Node,i_Tree) = &
                           GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

                        GP_Individual_Node_Parameters(i_Node,i_Tree) = 0.0d0
                    enddo ! i_node
                enddo ! i_tree

                !-----------------------------------------------------------------------------------

                ! calculate how many variables are in the tree

                n_GP_vars = 0
                do  i_Tree=1,n_Trees
                    do  i_Node=1,n_Nodes

                        if( GP_Individual_Node_Type(i_Node,i_Tree) < 0  .and. &
                            GP_Individual_Node_Type(i_Node,i_Tree) > -9999  ) then
                            n_GP_vars = n_GP_vars + 1
                        endif ! GP_Individual_Node_Type(i_Node,i_Tree) > 0 ....

                    enddo ! i_node
                enddo ! i_tree

                !------------------------------------------------------------------------------

                if( new_rank == 0 )then
                    write(GP_print_unit,'(A,1x,i3,5x,A,2(1x,i3))')&
                          'gil: for i_GP_Individual', i_GP_Individual, &
                          'number of parameters, variables =', &
                                n_GP_parameters, n_GP_vars
                    !flush(GP_print_unit)
                endif !  new_rank == 0


                !-------------------------------------------------------------------

                ! cycle the i_GP_individual loop if there are no GP parameters
                ! or if n_GP_parameters <=  n_code_equations

                if( n_GP_parameters == 0 .or. &
                    n_GP_parameters > n_maximum_number_parameters .or.  &
                    n_GP_parameters <=  n_code_equations                 ) then 

                    if( new_rank == 0 )then
                        write(GP_print_unit,'(A,1x,i3)')&
                              'gil: skipping this i_GP_Individual --&
                              &  the number of parameters is ', n_GP_parameters
                        flush(GP_print_unit)
                    endif !  new_rank == 0


                    individual_fitness = 0.0d0

                    if( new_rank == 0 )then
                        write(GP_print_unit,'(A,7(1x,I3), 1x, E15.7)')&
                         'gil: myid, new_rank, i_part, i_gp_1, i_gp_2, &
                               &i_GP_gen, i_GP_indiv, indiv_fit', &
                               myid, new_rank, i_part, i_gp_1, i_gp_2, &
                               i_GP_generation, i_GP_individual, &
                               individual_fitness
                    endif !  new_rank == 0

                    cycle gp_ind_loop

                endif ! n_GP_parameters == 0

                !-------------------------------------------------------------------


                if( new_rank == 0 )then

                    write(GP_print_unit,'(/A,2(1x,i3))') &
                     'gil: call GPCODE_GA_lmdif_Param_Opt         i_GP_Gen, i_GP_indiv', &
                                               i_GP_Generation, i_GP_individual
                    !flush(GP_print_unit)

                    if( L_ga_print )then
                        write(GA_print_unit,'(//A/A)') &
                         '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>&
                         &>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
                         'gil: call GPCODE_GA_lmdif_Parameter_Optimization routine'
                        write(GA_print_unit,'(A,2(1x,I3))') &
                         'gil: i_GP_Generation, i_GP_individual',&
                               i_GP_Generation, i_GP_individual
                    endif ! L_ga_print

                endif ! new_rank == 0


                !----------------------------------------------------------------------------


                ! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND
                ! LINK THE SSE OUTPUT TO THE ARRAY AT THE END
                ! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEED TO BE PULLED OUT

                !      individual_fitness
                !   GP_Individual_Initial_Conditions
                !   GP_Individual_Node_Parameters

                ! these arrays are broadcast in GPCODE_GA...


                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                call GPCODE_GA_lmdif_Parameter_Optimization( &
                                  i_GP_Generation,i_GP_individual, &
                                  new_group, new_comm  ) ! , ind2-ind1+1  )

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                if( new_rank == 0 )then
                    write(GP_print_unit,'(A,1x,i3)') &
                     'gil: aft call GPCODE_GA_lmdif_Parameter_Optimization routine'
                    !flush(GP_print_unit)
                endif ! new_rank == 0


                fit_buffer_send(i_GP_individual) = individual_fitness

                if( new_rank == 0 )then
                    write(GP_print_unit,'(/A,3(1x,I7),1x,E15.7/)')&
                          'gil:9t myid, new_rank, i_GP_individual, fit_buf_send', &
                                  myid, new_rank, i_GP_individual, &
                                  fit_buffer_send(i_GP_individual)
                endif ! new_rank == 0


                !--------------------------------------------------------------------------------


                sse_buffer_send( i_GP_individual ) = Individual_SSE_best_parent


            endif !   Run_GP_Calculate_Fitness(i_GP_Individual)

        enddo  gp_ind_loop    !   i_GP_individual


        !--------------------------------------------------------------------------------
        !  AFTER LOOP ON GP INDIVIDUALS  --  still in partition loop
        !--------------------------------------------------------------------------------


        !--------------------------------------------------------------------------------

   
        if( new_rank == 0 )then

            !  send the fitness buffer for the GP individuals already completed

            tag_fit_s = tag_ind_fit

            write(GP_print_unit,'(A,5(1x,I7)/)')&
                  'gil:9t myid, new_rank, tag_fit_s, ind1, ind2', &
                          myid, new_rank, tag_fit_s, ind1, ind2     

            n_indiv = ind2 - ind1 + 1
            call MPI_SEND( fit_buffer_send(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                          0,  tag_fit_s, MPI_COMM_WORLD, ierr )


           !--------------------------------------------------------------------------------

            !  send the SSE buffer for the GP individuals already completed

            tag_sse_s = tag_ind_sse

            write(GP_print_unit,'(/A,4(1x,I7)/)')&
                  'gil:10t myid, new_rank, tag_sse_s, i_GP_individual', &
                           myid, new_rank, tag_sse_s, i_GP_individual

            n_indiv = ind2 - ind1 + 1
            call MPI_SEND( sse_buffer_send(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                           0, tag_sse_s, MPI_COMM_WORLD, ierr )

        endif ! new_rank == 0

        !-----------------------------------------------------------------------------------


    endif ! i_gp_1 <= myid  .and. ...


    !write(6,'(A,4(1x,i3))')&
    !      'gil:12 AFT BARRIER myid, new_rank, ind1, ind2 ', &
    !                          myid, new_rank, ind1, ind2
    
enddo  part_loop





call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!write(6,'(/A,4(1x,i3), 1x, E15.7/)')&
!              'gil:  AFT BARRIER myid, new_rank ', &
!                                 myid, new_rank

if( myid == 0 )then
    write(GP_print_unit,'(//A//)' )&
         '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    do  ii = 1, n_GP_individuals
        write(GP_print_unit,'(A,4(1x,i3), 2(1x, E15.7))')&
         'gil:5 myid, new_rank, i_GP_gen, ii, GP_pop_fit, child_indiv_SSE', &
                myid, new_rank, i_GP_generation, ii, &
                           GP_population_fitness(ii), &
                         GP_Child_Individual_SSE(ii)
   enddo ! i_GP_individual
endif !  myid == 0

deallocate( fit_buffer_send )
deallocate( sse_buffer_send )
deallocate( sse_buffer_recv )



return



end subroutine GP_individual_loop
