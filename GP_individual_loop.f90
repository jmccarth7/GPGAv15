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
integer(kind=4) :: i_code_equation

integer(kind=4) :: ii
integer(kind=4) :: ii2
integer(kind=4) :: ii3
integer(kind=4) :: jj
integer(kind=4) :: jj2

integer(kind=4) :: message_len

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
integer(kind=4),parameter :: tag_parm    = 500000
integer(kind=4),parameter :: tag_node_type    = 600000
integer(kind=4),parameter :: tag_init_cond    = 700000
integer(kind=4),parameter :: tag_node_parm    = 800000


integer(kind=4) :: tag_fit_r
integer(kind=4) :: tag_fit_s
integer(kind=4) :: tag_sse_r
integer(kind=4) :: tag_sse_s

real(kind=8) ::   individual_fit_rec
real(kind=8) ::   individual_SSE_rec


real(kind=8),allocatable,dimension(:)    ::   init_cond_buff

real(kind=8),allocatable,dimension(:,:)  ::   node_parm_buff2

integer,allocatable,dimension(:,:)       ::   node_type_buff2

real(kind=8),allocatable,dimension(:)  ::   fit_buffer_send
real(kind=8),allocatable,dimension(:)  ::   sse_buffer_send
integer,allocatable,dimension(:)       ::   buff_parm_send


!---------------------------------------------------------------------------------------


if( myid /= 0 )then

!    if( myid == 1 )then
!        write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!         'gil: before allocate myid, (n_GP_individuals / n_partitions)', &
!                               myid, (n_GP_individuals / n_partitions)
!    endif !  myid == 1

    if( .not. allocated( fit_buffer_send ) )then
        allocate( fit_buffer_send( n_GP_individuals / n_partitions ) )
        allocated_memory = allocated_memory + &
                           real( n_GP_Individuals / n_partitions * 8, kind=8 )
    endif ! .not. allocated( fit_buffer_send )
    
!    if( myid == 1 )then
!        write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!         'gil: aft allocate  fit_buf_send  myid n_GP_Individuals / n_partitions * 8 ', &
!                               myid,  n_GP_Individuals / n_partitions * 8
!    endif !  myid == 1


    if( .not. allocated( sse_buffer_send ) )then
        allocate( sse_buffer_send( n_GP_individuals / n_partitions ) )
        allocated_memory = allocated_memory + &
                           real( n_GP_Individuals / n_partitions * 8, kind=8 )
    endif ! .not. allocated( sse_buffer_send )
!    if( myid == 1 )then
!        write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!         'gil: aft allocate  sse_buf_send  myid n_GP_Individuals / n_partitions * 8 ', &
!                               myid, n_GP_Individuals / n_partitions * 8
!    endif !  myid == 1
    
    if( .not. allocated( buff_parm_send ) )then
        allocate( buff_parm_send(  n_GP_individuals / n_partitions ) )
        allocated_memory = allocated_memory + &
                           real( n_GP_Individuals / n_partitions * 4, kind=8 )
    endif ! .not. allocated( buff_parm_send ) 

!    if( myid == 1 )then
!        write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!         'gil: aft allocate  buff_parm_send  myid n_GP_Individuals / n_partitions * 4 ', &
!                               myid, n_GP_Individuals / n_partitions * 4
!    endif !  myid == 1

endif ! myid /= 0


if( .not. allocated( init_cond_buff ) )then 
    allocate( init_cond_buff(  n_code_equations ) ) 
    allocated_memory = allocated_memory + real( n_code_equations * 8, kind=8 )
endif !  .not. allocated( init_cond_buff ) 

!if( myid == 0 )then
!    write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!     'gil: aft allocate  init_cond_buff myid n_code_equations*8 ', &
!                                         myid,  n_code_equations*8
!endif !  myid == 0

if( .not. allocated( node_parm_buff2 ) )then 
    allocate( node_parm_buff2(  n_nodes, n_trees ) )
    allocated_memory = allocated_memory + real( n_nodes * n_trees  * 8, kind=8 )
endif ! .not. allocated( node_parm_buff2 ) 

!if( myid == 0 )then
!    write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!     'gil: aft allocate  node_parm_buff2 myid  n_nodes * n_trees * 8 ', &
!                                         myid,  n_nodes * n_trees * 8
!endif !  myid == 0


if( .not. allocated( node_type_buff2 ) )then 
    allocate( node_type_buff2(  n_nodes, n_trees ) )
    allocated_memory = allocated_memory + real( n_nodes * n_trees * 4, kind=8 )
endif ! .not. allocated( node_type_buff2 ) 


!if( myid == 0 )then
!    write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!     'gil: aft allocate  node_type_buff2 myid  n_nodes * n_trees * 4 ', &
!                                         myid,  n_nodes * n_trees * 4
!endif !  myid == 0



divider = ( numprocs - 1 ) / n_partitions

call mpi_comm_rank( new_comm, new_rank, ierr )
call mpi_comm_size( new_comm, n_procs,  ierr )


!if( myid == 0 )then
!    write(GP_print_unit,'(/A,4(1x,i4))')&
!     'gil: before loop myid, new_rank, n_code_equations,  n_GP_params ',&
!                       myid, new_rank, n_code_equations,  n_GP_parameters
!    write(GP_print_unit,'(A,5(1x,i4))')&
!     'gil: before loop myid, new_rank, numprocs, divider, n_partitions', &
!                       myid, new_rank, numprocs, divider, n_partitions
!    write(GP_print_unit,'(A,3(1x,i4))')&
!     'gil: before loop myid, new_rank, n_procs                        ', &
!                       myid, new_rank, n_procs
!    write(GP_print_unit,'(A,1x,i4, 1x,I6)')&
!     'gil: before loop myid, (n_GP_individuals / n_partitions)', &
!                       myid, (n_GP_individuals / n_partitions)
!endif !  myid == 0


!if( myid == 0 )then
!    write(GP_print_unit,'(A,1x,i4, 1x,E15.7)')&
!     'gil: before loop myid, allocated_memory ', myid, allocated_memory
!    flush(GP_print_unit)
!endif !  myid == 0


!if( myid == 1 )then
!    write(GP_print_unit,'(A,1x,i4, 1x,E15.7)')&
!     'gil: before loop myid, allocated_memory ', myid, allocated_memory
!    flush(GP_print_unit)
!endif !  myid == 1

!-------------------------------------------------------------------------------


! do the loop over the GP individuals in n_partitions chunks

part_loop:&
do  i_part = 1,  n_partitions

    ! i_gp_1 and i_gp_2 are limits on the processor number for this partition

    i_gp_1 = divider * (i_part - 1) + 1
    i_gp_2 = divider *  i_part

    !---------------------------------------------------------------------------------
    ! ind1 and ind2 are limits on the i_GP_individuals processed in this partition

    ind1 =  (n_GP_individuals / n_partitions) * (i_part-1)  +  1
    ind2 =  (n_GP_individuals / n_partitions) *  i_part

    ! get any remaining individuals in the last partition

    if( i_part == n_partitions )then
        ind2 = n_GP_individuals
    endif ! i_part == n_partitions

    !write(GP_print_unit,'(A,7(1x,I5))')&
    !    'gil:in loop myid, new_rank, i_part, ind1, ind2, i_gp_1, i_gp_2',&
    !                 myid, new_rank, i_part, ind1, ind2, i_gp_1, i_gp_2

    !---------------------------------------------------------------------------------

    ! initialize the send buffers to values of the GP arrays
    ! this is necessary because if Run_GP_Fitness is .FALSE., for an individual
    ! then no value will be put into the place of that individual
    ! This code will then retain the value for that individual from the previous
    ! generation.

    !init_cond_buff  = GP_Population_Initial_Conditions(:, ind1:ind2 )

    if( myid /= 0 )then
        fit_buffer_send(1:ind2-ind1+1)  = GP_Population_Ranked_Fitness(ind1:ind2)
        sse_buffer_send(1:ind2-ind1+1)  = GP_Child_Individual_SSE(ind1:ind2)
        buff_parm_send(1:ind2-ind1+1)   = GP_Individual_N_GP_param(ind1:ind2)
    endif ! myid /= 0

!    write(GP_print_unit,'(A,8(1x,I4))')&
!        'gil:in loop myid, new_rank, i_part, ind1, ind2, i_gp_1, i_gp_2, ind2-ind1+1',&
!                     myid, new_rank, i_part, ind1, ind2, i_gp_1, i_gp_2, ind2-ind1+1

    !---------------------------------------------------------------------------------

    if( myid == 0 ) then

        !-------------------------------------------------------------------

        ! receive the number of GP parameters

        !write(GP_print_unit,'(A,5(1x,I7)/)')&
        !      'gil:30r myid, new_rank, tag_parm, ind1, ind2', &
        !               myid, new_rank, tag_parm, ind1, ind2


        n_indiv = ind2 - ind1 + 1
        call MPI_RECV( GP_Individual_N_GP_param(ind1), n_indiv, MPI_INTEGER, &
                       MPI_ANY_SOURCE, tag_parm,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        !write(GP_print_unit,'(A,4(1x,I7)/(10(1x,i4)))')&
        !      'gil:30r myid, n_indiv, ind1, ind2, GP_Individual_N_GP_param(ind1:ind2)', &
        !               myid, n_indiv, ind1, ind2, GP_Individual_N_GP_param(ind1:ind2)

        !-------------------------------------------------------------------

        ! receive the fitness information

        tag_fit_r =  tag_ind_fit

        !write(GP_print_unit,'(A,5(1x,I7)/)')&
        !      'gil:31r myid, new_rank, tag_fit_r, ind1, ind2', &
        !               myid, new_rank, tag_fit_r, ind1, ind2

        n_indiv = ind2 - ind1 + 1

        call MPI_RECV( GP_Population_Ranked_Fitness(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_fit_r,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

        !write(GP_print_unit,'(A,4(1x,I7)/(5(1x,E15.7)))')&
        !      'gil:31r myid, n_indiv, ind1, ind2, GP_Population_Ranked_Fitness(ind1:ind2)', &
        !               myid, n_indiv, ind1, ind2, GP_Population_Ranked_Fitness(ind1:ind2)


        !-------------------------------------------------------------------


        ! receive the SSE information

        tag_sse_r = tag_ind_sse

        n_indiv = ind2 - ind1 + 1

!        write(GP_print_unit,'(/A,7(1x,I7))')&
!              'gil:32r myid, new_rank, i_part, tag_sse_r, ind1, ind2, n_indiv', &
!                       myid, new_rank, i_part, tag_sse_r, ind1, ind2, n_indiv

        call MPI_RECV( GP_Child_Individual_SSE(ind1), n_indiv, MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, tag_sse_r,                       &
                       MPI_COMM_WORLD, MPI_STAT, ierr )

!        write(GP_print_unit,'(A,3(1x,I4)/(5(1x,E15.7)))')&
!              'gil:32r myid, new_rank, i_part, GP_Child_Individual_SSE(ind1:ind2)', &
!                       myid, new_rank, i_part, GP_Child_Individual_SSE(ind1:ind2)

        !GP_Child_Individual_SSE(ind1:ind2) =  sse_buffer_recv(1:n_indiv)
        GP_Adult_Individual_SSE(ind1:ind2) =  GP_Child_Individual_SSE(ind1:ind2)
        GP_Adult_Population_SSE(ind1:ind2) =  GP_Child_Individual_SSE(ind1:ind2)


        !--------------------------------------------------------------------------------

        n_indiv = ind2 - ind1 + 1
        message_len = n_Nodes  * n_Trees

        !write(GP_print_unit,'(A,7(1x,I7)/)')&
        !      'gil:33r myid,new_rank,tag_node_type,ind1,ind2,n_indiv,message_len', &
        !               myid,new_rank,tag_node_type,ind1,ind2,n_indiv,message_len

        do  jj = ind1, ind2 

            call MPI_RECV( GP_Adult_Population_Node_Type(:,:,jj), message_len,    &
                           MPI_INTEGER,  MPI_ANY_SOURCE, tag_node_type+jj, &
                           MPI_COMM_WORLD, MPI_STAT, ierr )
        enddo

        !write(GP_print_unit,'(A)')&
        !      'gil:33r indiv i_tree i_node GP_Adult_Pop_Node_Type(i_node, i_tree, indiv)'

        !do  ii = ind1, ind2
        !    do  i_tree = 1, n_trees
        !        do  i_node = 1, n_nodes
        !            if( GP_Adult_Population_Node_Type(i_node,i_tree,ii) > -9999 )then
        !                write(GP_print_unit,'(4(1x,I6))' )&
        !                ii,i_tree, i_node, GP_Adult_Population_Node_Type(i_node,i_tree,ii)
        !            endif ! GP_Adult_Population_Node_Type(i_node,i_tree,ii) > -9999
        !        enddo !  i_node
        !    enddo !  i_tree
        !enddo ! ii

        GP_Child_Population_Node_Type =  GP_Adult_Population_Node_Type

        !--------------------------------------------------------------------------------


        n_indiv = ind2 - ind1 + 1
        message_len = n_code_equations

        !write(GP_print_unit,'(A,7(1x,I7)/)')&
        !  'gil:34r myid,new_rank,tag_init_cond,ind1,ind2,n_indiv, message_len', &
        !           myid,new_rank,tag_init_cond,ind1,ind2,n_indiv, message_len

        do  ii = ind1, ind2 

            call MPI_RECV( GP_Population_Initial_Conditions(1,ii), message_len,    &
                           MPI_double_precision,  MPI_ANY_SOURCE, tag_init_cond+ii, &
                           MPI_COMM_WORLD, MPI_STAT, ierr )

        enddo ! ii

        !write(GP_print_unit,'(A,4(1x,I7)/(5(1x,E15.7)))')&
        !  'gil:34r myid, n_indiv, ind1, ind2, GP_Population_Initial_Conditions(:,ind1:ind2)', &
        !           myid, n_indiv, ind1, ind2, GP_Population_Initial_Conditions(:,ind1:ind2)


        !--------------------------------------------------------------------------------


        n_indiv = ind2 - ind1 + 1
        message_len = n_Nodes * n_Trees

        !write(GP_print_unit,'(A,4(1x,I7)/)')&
        !      'gil:35r myid, new_rank, n_nodes, n_trees', &
        !               myid, new_rank, n_nodes, n_trees

        !write(GP_print_unit,'(A,7(1x,I7)/)')&
        !      'gil:35r myid, new_rank, tag_init_cond, ind1, ind2, n_indiv, message_len', &
        !               myid, new_rank, tag_init_cond, ind1, ind2, n_indiv, message_len

        do  jj2 = ind1, ind2 

            call MPI_RECV( GP_Population_Node_Parameters(1,1,jj2), message_len,   &
                        MPI_double_precision,  MPI_ANY_SOURCE, tag_node_parm+jj2, &
                        MPI_COMM_WORLD, MPI_STAT, ierr )
        enddo ! jj2 

        !write(GP_print_unit,'(A)')&
        !  'gil:35r indiv i_tree i_node GP_Pop_Node_Params(i_node, i_tree,,indiv)'

        !do  ii = ind1, ind2
        !    do  i_tree = 1, n_trees
        !        do  i_node = 1, n_nodes
        !            if( isnan( GP_Population_Node_Parameters(i_node,i_tree,ii) ) )then
        !                write(GP_print_unit,'(A,5(1x,I7)/)')&
        !                  'gil:35r ii, i_tree, i_node, GP_Population_Node_Parameters is NAN', &
        !                           ii, i_tree, i_node
        !            endif !  isnan( GP_Population_Node_Parameters(i_node,i_tree,ii)
        !            if( GP_Population_Node_Parameters(i_node,i_tree,ii) > 0.0d0 )then
        !                write(GP_print_unit,'(3(1x,I6),1x,E15.7)' )&
        !                ii, i_tree, i_node, GP_Population_Node_Parameters(i_node,i_tree,ii)
        !            endif ! GP_Population_Node_Parameters(i_node,i_tree,ii) > 0.0d0
        !        enddo !  i_node
        !    enddo !  i_tree
        !enddo ! ii

        !--------------------------------------------------------------------------------



    elseif( i_gp_1 <= myid  .and.   &
                      myid   <= i_gp_2         )then

        gp_ind_loop:&
        do  i_GP_individual= ind1, ind2    ! 1,n_GP_individuals

            !write(GP_print_unit,'(A,8(1x,I4))')&
            !   'gil:2 myid, new_rank, i_part, i_gp_1, i_gp_2, ind1, ind2, i_GP_individual ', &
            !          myid, new_rank, i_part, i_gp_1, i_gp_2, ind1, ind2, i_GP_individual


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


            !if( new_rank == 0 )then
            !    write(GP_print_unit,'(A,2(1x,I5))')&
            !     'gil: after loop n_code_equations,  n_GP_params',&
            !                      n_code_equations,  n_GP_parameters
            !    !flush(GP_print_unit)
            !endif !  new_rank == 0


            !------------------------------------------------------------------------

            GP_Individual_N_GP_param(i_GP_individual) = n_GP_parameters

            buff_parm_send(i_GP_individual-ind1+1) = n_GP_parameters



            !  send the number of parameters for the GP individual
            !call MPI_SEND( n_GP_parameters, 1, MPI_INTEGER,        &
            !              0,  i_GP_individual, MPI_COMM_WORLD, ierr )

            !------------------------------------------------------------------------------


            ! run GPCODE_... to evaluate this individual  if Run_GP_Calculate_Fitness is true


            if( Run_GP_Calculate_Fitness(i_GP_Individual) ) then

                !-----------------------------------------------------------------------------------

                ! these get set randomly in the GA-lmdif search algorithm ( in GPCODE* )
                !GP_Individual_Node_Parameters(1:n_Nodes,1:n_Trees) = 0.0d0               ! 20131209

                !-----------------------------------------------------------------------------------

                !if( new_rank == 0 )then
                !    write(GP_print_unit,'(/A)')&
                !    'gil:----------------------------------------------------------------------'
                !    write(GP_print_unit,'(A,2(1x,i5),3x,L1)') &
                !          'gil: i_GP_Gen, i_GP_indiv, Run_GP_Calculate_Fitness', &
                !                i_GP_Generation, i_GP_individual, &
                !                       Run_GP_Calculate_Fitness(i_GP_Individual)
                !    write(GP_print_unit,'(A)')&
                !    'gil:----------------------------------------------------------------------'
                !    !!flush(GP_print_unit)
                !endif !  new_rank == 0

                !-----------------------------------------------------------------------------------

                !if( new_rank == 0 )then
                !    write(GP_print_unit,'(/A,4(1x,i6))') &
                !      'gil: i_GP_individual, n_trees, n_nodes, n_GP_parameters ', &
                !            i_GP_individual, n_trees, n_nodes, n_GP_parameters
                !    !flush(GP_print_unit)
                !endif !  new_rank == 0

                !-----------------------------------------------------------------------------------

                do  i_Tree=1,n_Trees
                    do  i_Node=1,n_Nodes

                        GP_Individual_Node_Type(i_Node,i_Tree) = &
                           GP_Adult_Population_Node_Type(i_Node,i_Tree,i_GP_Individual)

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

                !if( new_rank == 0 )then
                !    write(GP_print_unit,'(A,1x,i5,5x,A,2(1x,i5))')&
                !          'gil: for i_GP_Individual', i_GP_Individual, &
                !          'number of parameters, variables =', &
                !                n_GP_parameters, n_GP_vars
                !    !flush(GP_print_unit)
                !endif !  new_rank == 0


                !-------------------------------------------------------------------

                ! cycle the i_GP_individual loop if there are no GP parameters
                ! or if n_GP_parameters <=  n_code_equations

                if( n_GP_parameters == 0 .or. &
                    n_GP_parameters > n_maximum_number_parameters .or.  &
                    n_GP_parameters <=  n_code_equations                 ) then

                    !if( new_rank == 0 )then
                    !    write(GP_print_unit,'(A,1x,i5)')&
                    !          'gil: skipping this i_GP_Individual --&
                    !          &  the number of parameters is ', n_GP_parameters
                    !    flush(GP_print_unit)
                    !endif !  new_rank == 0


                    individual_fitness = 0.0d0

                    !if( new_rank == 0 )then
                    !    write(GP_print_unit,'(A,7(1x,I5), 1x, E15.7)')&
                    !     'gil: myid, new_rank, i_part, i_gp_1, i_gp_2, &
                    !           &i_GP_gen, i_GP_indiv, indiv_fit', &
                    !           myid, new_rank, i_part, i_gp_1, i_gp_2, &
                    !           i_GP_generation, i_GP_individual, &
                    !           individual_fitness
                    !endif !  new_rank == 0

                    cycle gp_ind_loop

                endif ! n_GP_parameters == 0

                !-------------------------------------------------------------------


                !if( new_rank == 0 )then

                    !write(GP_print_unit,'(/A,2(1x,I5))') &
                    ! 'gil: call GPCODE_GA_lmdif_Param_Opt         i_GP_Gen, i_GP_indiv', &
                    !                           i_GP_Generation, i_GP_individual
                    !flush(GP_print_unit)

                    !if( L_ga_print )then
                    !    write(GA_print_unit,'(//A/A)') &
                    !     '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>&
                    !     &>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
                    !     'gil: call GPCODE_GA_lmdif_Parameter_Optimization routine'
                    !    write(GA_print_unit,'(A,2(1x,I5))') &
                    !     'gil: i_GP_Generation, i_GP_individual',&
                    !           i_GP_Generation, i_GP_individual
                    !endif ! L_ga_print

                !endif ! new_rank == 0


                !----------------------------------------------------------------------------


                ! THIS IS WHERE YOU NEED TO INSERT THE GA_LMDIF CALL AND
                ! LINK THE SSE OUTPUT TO THE ARRAY AT THE END
                ! ALSO, THE OPTIMAL PARAMETER SETS FROM THE BEST CHILD NEED TO BE PULLED OUT

                !      individual_fitness
                !   GP_Individual_Initial_Conditions
                !   GP_Individual_Node_Parameters

                ! these arrays are broadcast in GPCODE_GA...


                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !write(GP_print_unit,'(A,4(1x,I5))') &
                !      'gil: call GPCODE myid, new_rank, i_GP_individual',&
                !                        myid, new_rank, i_GP_individual

                call GPCODE_GA_lmdif_Parameter_Optimization( &
                                  i_GP_Generation,i_GP_individual, &
                                  new_group, new_comm  )

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                !if( new_rank == 0 )then
                !    write(GP_print_unit,'(A,1x,I5)') &
                !     'gil: aft call GPCODE_GA_lmdif_Parameter_Optimization routine'
                !    !flush(GP_print_unit)
                !endif ! new_rank == 0


                fit_buffer_send(i_GP_individual - ind1 + 1) = individual_fitness

                !if( new_rank == 0 )then
                !    write(GP_print_unit,'(/A,3(1x,I7),1x,E15.7/)')&
                !          'gil:9t myid, new_rank, i_GP_individual, fit_buf_send', &
                !                  myid, new_rank, i_GP_individual, &
                !                  fit_buffer_send(i_GP_individual - ind1 + 1)
                !endif ! new_rank == 0


                !--------------------------------------------------------------------------------


                sse_buffer_send( i_GP_individual-ind1+1 ) = Individual_SSE_best_parent

                !if( new_rank == 0 )then
                !    write(GP_print_unit,'(A,3(1x,I7),1x,E15.7)')&
                !     'gil:9t myid, new_rank, i_GP_indiv, sse_buf_send(i_GP_indiv-ind1+1)', &
                !             myid, new_rank, i_GP_individual, &
                !             sse_buffer_send(i_GP_individual-ind1+1)
                !endif ! new_rank == 0


                !--------------------------------------------------------------------------------


                ! set the GA_lmdif-optimized initial condition array

                do  i_code_equation = 1, n_code_equations  ! 20131209

                    GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual) = &
                        GP_Individual_Initial_Conditions(i_CODE_Equation)

                enddo ! i_code_equation

                !-----------------------------------------------------------------------------------


                ! set the GA_lmdif-optimized CODE parameter set array

                do  i_Tree=1,n_Trees    ! 20131209
                    do  i_Node=1,n_Nodes

                        GP_Population_Node_Parameters(i_node, i_tree, i_GP_Individual) = &
                               GP_Individual_Node_Parameters(i_node, i_tree)

                        GP_Adult_Population_Node_Type(i_node, i_tree, i_GP_Individual) = &
                                         GP_Individual_Node_Type(i_node, i_tree)

                        GP_Child_Population_Node_Type(i_node, i_tree, i_GP_Individual) = &
                                         GP_Individual_Node_Type(i_node, i_tree)

                    enddo ! i_node
                enddo ! i_tree

                !------------------------------------------------------------------------------



            endif !   Run_GP_Calculate_Fitness(i_GP_Individual)

            !write(GP_print_unit,'(A,4(1x,I6))') &
            !    '0: call summary_GP_indiv myid, new_rank, i_GP_generation, i_GP_individual ', &
            !                              myid, new_rank, i_GP_generation, i_GP_individual
            !call summary_GP_indiv( i_GP_generation, i_GP_individual, 0 )





        enddo  gp_ind_loop    !   i_GP_individual


        !--------------------------------------------------------------------------------
        !  AFTER LOOP ON GP INDIVIDUALS  --  still in partition loop
        !--------------------------------------------------------------------------------


        !--------------------------------------------------------------------------------


        if( new_rank == 0 )then



            !  send the number of parameters for the GP individual

            n_indiv = ind2 - ind1 + 1
            call MPI_SEND( buff_parm_send, n_indiv, MPI_INTEGER,        &
                           0,  tag_parm, MPI_COMM_WORLD, ierr )

            !--------------------------------------------------------------------------------

            !  send the fitness buffer for the GP individuals already completed

            tag_fit_s = tag_ind_fit

            !write(GP_print_unit,'(A,5(1x,I7)/)')&
            !      'gil:9t myid, new_rank, tag_fit_s, ind1, ind2', &
            !              myid, new_rank, tag_fit_s, ind1, ind2

            n_indiv = ind2 - ind1 + 1
            call MPI_SEND( fit_buffer_send, n_indiv, MPI_DOUBLE_PRECISION, &
                           0,  tag_fit_s, MPI_COMM_WORLD, ierr )


            !--------------------------------------------------------------------------------

            !  send the SSE buffer for the GP individuals already completed

            tag_sse_s = tag_ind_sse

            !write(GP_print_unit,'(/A,4(1x,I7)/)')&
            !      'gil:10t myid, new_rank, tag_sse_s, i_GP_individual', &
            !               myid, new_rank, tag_sse_s, i_GP_individual

            !write(GP_print_unit,'(/A,5(1x,I7)/(5(1x,E15.7)))')&
            !      'gil:10t myid, new_rank, i_part, tag_sse_s, n_indiv, sse_buffer_send ', &
            !               myid, new_rank, i_part, tag_sse_s, n_indiv, sse_buffer_send 

            n_indiv = ind2 - ind1 + 1
            call MPI_SEND( sse_buffer_send, n_indiv, MPI_DOUBLE_PRECISION, &
                           0, tag_sse_s, MPI_COMM_WORLD, ierr )

            !--------------------------------------------------------------------------------

            n_indiv = ind2 - ind1 + 1
            message_len = n_nodes * n_Trees

            do  ii = ind1 , ind2 

                node_type_buff2(:,:) = GP_Adult_Population_Node_Type(:,:,ii)

                call MPI_SEND( node_type_buff2, message_len,    &
                            MPI_INTEGER,  0, tag_node_type+ii, MPI_COMM_WORLD, ierr )

            enddo  ! ii
            !--------------------------------------------------------------------------------


            n_indiv = ind2 - ind1 + 1
            message_len = n_code_equations

            !write(GP_print_unit,'(/A,6(1x,I7)/)')&
            !      'gil: myid, new_rank, ind1, ind2, n_indiv, message_len', &
            !            myid, new_rank, ind1, ind2, n_indiv, message_len

            do  ii3 = ind1, ind2

                init_cond_buff(:) = GP_Population_Initial_Conditions(:, ii3)

                call MPI_SEND( init_cond_buff, message_len,    &
                               MPI_double_precision,  0, tag_init_cond+ii3, &
                               MPI_COMM_WORLD, ierr )

            enddo  ! ii3 
            !--------------------------------------------------------------------------------


            n_indiv = ind2 - ind1 + 1
            message_len = n_Nodes * n_Trees

            do  ii2 = ind1, ind2

                node_parm_buff2(:, :) = GP_Population_Node_Parameters(:,:,ii2)

                call MPI_SEND( node_parm_buff2, message_len,                &
                               MPI_double_precision,  0, tag_node_parm+ii2, &
                               MPI_COMM_WORLD, ierr )

            enddo ! ii2

            !--------------------------------------------------------------------------------



        endif ! new_rank == 0

        !-----------------------------------------------------------------------------------


    endif ! i_gp_1 <= myid  .and. ...


    call MPI_BARRIER( MPI_COMM_WORLD, ierr )

    !write(GP_print_unit,'(A,6(1x,I5))')&
    !      'gil:12 AFT BARRIER myid, new_rank, i_part, ind1, ind2 ', &
    !                          myid, new_rank, i_part, ind1, ind2

enddo  part_loop






call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!write(GP_print_unit,'(/A,4(1x,I5), 1x, E15.7/)')&
!      'gil:  AFT BARRIER myid, new_rank ', &
!                         myid, new_rank

!---------------------------------------------------------------------------------

message_len =  n_GP_individuals
call MPI_BCAST( GP_Individual_N_GP_param, message_len,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    do  ii = 1, n_GP_individuals
!        write(GP_print_unit,'(A,5(1x,i4))')&
!         'gil:5 myid, new_rank, i_GP_gen, ii, GP_Individual_N_GP_param(ii) ', &
!                myid, new_rank, i_GP_generation, ii, &
!                        GP_Individual_N_GP_param(ii)
!   enddo ! ii
!endif !  myid == 0


!---------------------------------------------------------------------------------

!if( myid == 0 )then
!
!    write(GP_print_unit,'(//A//)' )&
!         '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!    do  ii = 1, n_GP_individuals
!        write(GP_print_unit,'(A,4(1x,I5), 2(1x, E15.7))')&
!         'gil:5 myid, new_rank, i_GP_gen, ii, GP_pop_fit, child_indiv_SSE', &
!                myid, new_rank, i_GP_generation, ii, &
!                           GP_Population_Ranked_Fitness(ii), &
!                         GP_Child_Individual_SSE(ii)
!            !write(GP_print_unit,'(A,4(1x,I6))') &
!            !    'gil: call summary_GP_indiv myid, new_rank, i_GP_generation, ii', &
!            !                                myid, new_rank, i_GP_generation, ii
!            !call summary_GP_indiv( i_GP_generation, ii, 0 )
!
!   enddo ! i_GP_individual
!
!endif !  myid == 0

if( myid /= 0 )then

    if( allocated( fit_buffer_send ) )then 
        deallocate( fit_buffer_send )
        allocated_memory = allocated_memory - &
                           real( n_GP_Individuals / n_partitions * 8, kind=8 )
    endif ! allocated( fit_buffer_send ) 

    if( allocated( sse_buffer_send ) )then 
        deallocate( sse_buffer_send )
        allocated_memory = allocated_memory - &
                           real( n_GP_Individuals / n_partitions * 8, kind=8 )
    endif ! allocated( sse_buffer_send )
    
    if( allocated( buff_parm_send ) )then 
        deallocate( buff_parm_send  )
        allocated_memory = allocated_memory - &
                           real( n_GP_Individuals / n_partitions * 8, kind=8 )
    endif ! allocated( buff_parm_send ) 
    
endif ! myid /= 0 

if( allocated( init_cond_buff ) )then 
    deallocate( init_cond_buff  )
    allocated_memory = allocated_memory - &
       real( n_code_equations * (n_GP_Individuals / n_partitions) * 8, kind=8 )
endif ! allocated( init_cond_buff ) 

if( allocated( node_parm_buff2 ) )then 
    deallocate( node_parm_buff2  )
    allocated_memory = allocated_memory - real( n_nodes * n_trees * 8, kind=8 )
endif ! allocated( node_parm_buff2 ) 

if( allocated( node_type_buff2 ) )then 
    deallocate( node_type_buff2  )
    allocated_memory = allocated_memory -  real( n_nodes * n_trees * 4, kind=8 )
endif ! allocated( node_type_buff2 ) 

call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!if( myid == 0 )then
!    write(GP_print_unit,'(A,1x,i4, 1x,E15.7)')&
!     'gil: at GP_ind return myid, allocated_memory ', myid, allocated_memory
!    flush(GP_print_unit)
!endif !  myid == 0



return



end subroutine GP_individual_loop
