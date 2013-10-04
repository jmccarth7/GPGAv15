subroutine GP_para_lmdif_process( i_GP_Generation )

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

integer :: i_GP_Generation
integer :: i_GP_individual

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
 dimension(n_maximum_number_parameters,n_GP_individuals) ::  parent_parameters
real(kind=8),&
 dimension(n_maximum_number_parameters,n_GP_individuals) ::  child_parameters


real(kind=8) :: individual_SSE_best_1
real(kind=8) :: individual_ranked_fitness_best_1
real(kind=8) :: Individual_Fitness_best_1

real(kind=8),dimension(n_maximum_number_parameters) :: parent_parameters_best_1



real(kind=8) :: buffer(n_maximum_number_parameters + 2)
real(kind=8) :: buffer_recv(n_maximum_number_parameters + 2)


integer (kind=4) ::      i
integer (kind=4) :: i_GP_Best_Parent
integer (kind=4) :: i_GP_Best_Parent_1

integer (kind=4) :: i_GP_generation_last

real(kind=8), parameter :: tol = 1.0d-30

real(kind=8),parameter :: zero = 0.0d0


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(n_GP_individuals)

real(kind=8), external :: indiv_fitness

logical :: L_stop_run


integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: jj
integer(kind=4) :: i_parameter

real(kind=8) :: t1
real(kind=8) :: t2


logical :: L_GP_print

!----------------------------------------------------------------------

L_GP_print = .TRUE.

i_dummy = 0


buffer(1:n_maximum_number_parameters+2)      = 0.0D0
buffer_recv(1:n_maximum_number_parameters+2) = 0.0D0



n_parameters = n_GP_parameters



if( myid == 0 )then
    if( L_GP_print )then
        write(GP_print_unit,'(//A)') 'gplp: at entry  '
        write(GP_print_unit,'(A,1x,E15.7)') 'gplp: dt ', dt
        write(GP_print_unit,'(A,1x,I10)') &
              'gplp: n_parameters    =   ', n_parameters
        write(GP_print_unit,'(A,1x,I10)') &
              'gplp: n_GP_parameters =   ', n_GP_parameters
    endif ! L_GP_print

endif ! myid == 0


! jjm 20130417 >>>>>>>>>>>>>>>
!if( myid == 0) then
!
!    if( L_GP_print )then
!        !write(GP_print_unit,'(A)')' '
!        !write(GP_print_unit,'(A)') &
!        !              'gplp: i_node, i_tree, GP_Indiv_Node_Param'
!        !do  i_tree=1,n_trees
!        !    do  i_node=1,n_nodes
!        !        !if( abs( GP_Individual_Node_Parameters(i_node,i_tree) ) > 1.0e-20 )then
!        !            !write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
!        !            !  'gplp: i_node, i_tree, GP_Indiv_Node_Param', &
!        !            write(GP_print_unit,'(8x,2(1x,I6),1x,E15.7)') &
!        !                          i_node, i_tree, GP_Individual_Node_Parameters(i_node,i_tree)
!        !        !endif ! abs( GP_Indiv_Node_Param(i_node,i_tree) ) > 1.0e-20
!        !    enddo ! i_node
!        !enddo  ! i_tree
!
!        write(GP_print_unit,'(/A)') &
!                     'gplp: i_node, i_tree, GP_Indiv_Node_Type'
!
!
!        do  i_tree=1,n_trees
!            do  i_node=1,n_nodes
!                if( GP_Individual_Node_Type(i_node,i_tree) > -9999 )then
!                    !write(GP_print_unit,'(A,3(1x,I6))') &
!                    !'gplp: i_node, i_tree, GP_Indiv_Node_Type', &
!                    write(GP_print_unit,'(8x,3(1x,I6))') &
!                                 i_node, i_tree, GP_Individual_Node_Type(i_node,i_tree)
!                endif ! GP_Indiv_Node_Type(i_node,i_tree) > -9999
!            enddo ! i_node
!        enddo  ! i_tree
!
!        write(GP_print_unit,'(A)')' '
!    endif ! L_GP_print
!
!endif ! myid == 0
! jjm 20130417 <<<<<<<<<<<<<<<




if( n_parameters .le. 0) then
    write(GP_print_unit,'(A)')        &
          'gplp: ERROR: n_parameters </= 0'
    write(GP_print_unit,'(A,1x,I10)') &
          'gplp: n_parameters =   ', n_parameters
    stop 'n_par<=0'
endif

!-----------------------------------------------------------------------------


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  debug
! print the node parameters (if there are any)

if( myid == 0 )then
    write(GP_print_unit,'(/A/)') &
       'gplp:  i_GP_indiv node  tree  GP_population_node_params'
    do  i_GP_individual = 1, n_GP_individuals
        do  i_tree=1,n_trees
            do  i_node=1,n_nodes
                if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) .eq. 0 ) then
                    if( GP_population_node_parameters(i_node,i_tree,i_GP_individual) > 0.0d0 )then
                        write(GP_print_unit,'(3(1x,I6),  4x, E20.10)') &
                         i_GP_individual, i_node, i_tree, &
                         GP_population_node_parameters(i_node,i_tree,i_GP_individual)
                    endif ! GP_population_node_parameters(i_node,i_tree,i_GP_individual) > 0.0d0
                endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0
    
            enddo ! i_node
        enddo  ! i_tree
    enddo ! i_GP_individual
endif ! myid == 0

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  debug


! load the population node parameters into the child parameters


if( myid == 0 )then

    child_parameters( 1:n_maximum_number_parameters, 1:n_GP_individuals) = 0.0d0
    do  i_GP_individual = 1, n_GP_individuals
    
        n_parameters = 0
    
    
        do  i_tree=1,n_trees
            do  i_node=1,n_nodes
    
                if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) .eq. 0 ) then
    
                    n_parameters = n_parameters + 1 
                   
                    child_parameters( n_parameters, i_GP_individual) =  &
                         GP_population_node_parameters(i_node,i_tree,i_GP_individual)
    
                endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0
    
            enddo ! i_node
        enddo  ! i_tree
    
    
    enddo ! i_GP_individual

endif ! myid == 0

!-----------------------------------------------------------------------------

! set up MPI process


if( myid == 0 )then

    if( L_GP_print )then
        write(GP_print_unit,'(/A,2(1x,I6))')&
          'gplp: n_GP_individuals, numprocs ', &
                 n_GP_individuals, numprocs
        write(GP_print_unit,'(A,1x,I6/)')&
          'gplp: n_parameters ', n_parameters
    endif ! L_GP_print

endif ! myid == 0

!-----------------------------------------------------------------------------

L_stop_run  = .FALSE.
!L_stop_run  = .TRUE.




!------------------------------------------------------------------------

!  broadcast child parameters


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,2(1x,I6))') &
!    'gplp:  broadcast child parameters myid, i_GP_generation ', &
!                                       myid, i_GP_generation
!endif ! L_GP_print

child_number =  n_GP_individuals * n_maximum_number_parameters

!if( L_GP_print )then
!    write(GP_print_unit,'(A,4(1x,I6)/)') &
!     'gplp:  myid, n_GP_individuals, n_maximum_number_parameters, child_number = ', &
!             myid, n_GP_individuals, n_maximum_number_parameters, child_number
!endif ! L_GP_print

call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

call MPI_BCAST( Child_Parameters,  child_number,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,2(1x,I10)/)') &
!     'gplp: child  broadcast myid, ierr = ', myid, ierr
!    write(GP_print_unit,'(/A,2(1x,I10)/)') &
!     'gplp: myid, n_GP_individuals = ', myid, n_GP_individuals
!    do  i_GP_individual = 1, n_GP_individuals
!        write(GP_print_unit,'(I3,1x,I3,1x,12(1x,E15.7))') &
!              myid, i_GP_individual, &
!              ( child_parameters(jj,i_GP_individual), jj = 1,n_parameters )
!    enddo ! i_GP_individual
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast Run_GA_lmdif

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!      'gplp:  broadcast Run_GA_lmdif i_GP_generation ', i_GP_generation
!endif ! L_GP_print

!call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?
!call MPI_BCAST( Run_GA_lmdif,  n_GP_individuals,    &
!                    MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
!call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?
!
!if( L_GP_print )then
!    write(GP_print_unit,'(A,1x,I10/)') &
!     'gplp: Run_GA_lmdif  broadcast ierr = ', ierr
!    write(GP_print_unit,'(/A,1x,I6,/(10(2x,L1)))') &
!     'gplp: myid, Run_GA_lmdif  ', myid, Run_GA_lmdif
!endif ! L_GP_print


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   begin RK fcn integration segment
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


if( myid == 0 )then
    if( L_GP_print )then
        write(GP_print_unit,'(/A,1x,I6/)') &
              'gplp: begin parallel lmdif segment i_GP_generation ', &
                                                  i_GP_generation
        write(GP_print_unit,'(A,1x,I6/)') 'gplp: n_GP_individuals ', &
                                                 n_GP_individuals
    endif ! L_GP_print
endif !  myid == 0


individual_quality( 1: n_GP_individuals ) = 1

nsafe = 0

if( myid == 0  )then

    ! processor 0 sends job assignments to N  processors
    ! where N is the smaller of (number of processors -1)
    !                    and the number of individuals

    ! numsent is the number of messages sent up to now

    numsent = 0
    i_GP_individual = 0

    do  isource = 1, min( numprocs-1, n_GP_individuals )


        i_GP_individual = i_GP_individual + 1

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,1x,I6, 4x, L1)') &
        !     'gplp:1 494 i_GP_individual, Run_GA_lmdif(i_GP_individual)  ', &
        !                 i_GP_individual, Run_GA_lmdif(i_GP_individual)
        !endif ! L_GP_print


        call MPI_SEND( i_dummy,  1, MPI_INTEGER,    &
                       isource, isource,  MPI_COMM_WORLD, ierr )
        numsent = numsent + 1

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,4(1x,I6))') &
        !     'gplp:1 504 myid, isource, i_GP_individual, numsent ', &
        !                 myid, isource, i_GP_individual, numsent
        !endif ! L_GP_print

    enddo ! isource


    ! at this point i_GP_individual = numsent


    !if( L_GP_print )then
    !    write(GP_print_unit,'(A,4(1x,I6))') &
    !         'gplp: aft source loop 1 myid, i_GP_individual, numsent ', &
    !                                  myid, i_GP_individual, numsent
    !endif ! L_GP_print

    !-------------------------------------------------------------------------------------

    ! processor 0 loops over the number of individuals and waits for a message
    ! from the other processors

    do  isource = 1, n_GP_individuals

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 521 myid, isource, numsent ', &
        !                 myid, isource, numsent
        !endif ! L_GP_print


        buffer_recv = 0.0d0
        call MPI_RECV( buffer_recv, n_parameters+2, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                       MPI_COMM_WORLD, MPI_STAT,  ierr )

        sender       = MPI_STAT( MPI_SOURCE )
        i_individual = MPI_STAT( MPI_TAG )


        ! received a message from processor "sender" which processed
        ! individual "i_individual"


        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,5(1x,I6))') &
        !     'gplp:2 529 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )', &
        !                 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )
        !    write(GP_print_unit,'(A,5(1x,I6))') &
        !     'gplp:2 529 myid, isource, numsent, sender, i_individual ', &
        !                 myid, isource, numsent, sender, i_individual
        !    write(GP_print_unit,'(A,2(1x,I6),4x,L1)') &
        !      'gplp:2 myid, i_individual, Run_GA_lmdif(i_individual)', &
        !              myid, i_individual, Run_GA_lmdif(i_individual)
        !endif ! L_GP_print


        ! store the information received in the above message

        !if( Run_GA_lmdif(i_individual) ) then

            child_parameters(1:n_parameters,i_individual) =  &
                                 buffer_recv(1:n_parameters)

            GP_Child_Individual_SSE(i_individual)     =       buffer_recv( n_parameters+1)
            individual_quality(i_individual) = nint( buffer_recv( n_parameters+2) )

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,3(1x,I6))') &
            !     'gplp:2 554 myid, n_parameters, i_individual ', &
            !                 myid, n_parameters, i_individual
            !    write(GP_print_unit,'(A/(5(1x,E15.7)))') &
            !     'gplp:2 child_parameters(1:n_parameters,i_individual) ', &
            !             child_parameters(1:n_parameters,i_individual)
            !    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
            !     'gplp:2 myid, i_individual, GP_Child_Individual_SSE(i_individual)', &
            !             myid, i_individual, GP_Child_Individual_SSE(i_individual)
            !    write(GP_print_unit,'(A,3(1x,I6))') &
            !     'gplp:2 myid, i_individual, individual_quality(i_individual) ', &
            !             myid, i_individual,  individual_quality(i_individual)
            !endif ! L_GP_print

        !endif ! Run_GA_lmdif(i_individual)

        !-------------------------------------------------------------------------------------

        ! check to see if all individuals have been processed

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,1x,I6, 4x,L1)') &
        !     'gplp:2 542 myid, numsent < n_GP_individuals ', &
        !                 myid, numsent < n_GP_individuals
        !endif ! L_GP_print

        if( numsent <  n_GP_individuals )then

            ! numsent <  n_GP_individuals    means not all individuals have been processed

            ! send a message to the processor "sender" which just sent a message saying it
            ! completed an individual, and tell it to process
            ! the individual "i_GP_individual" as the  "numsent+1"  task

            i_GP_individual = i_GP_individual + 1

            call MPI_SEND( i_GP_individual, 1, MPI_INTEGER,    &
                           sender, numsent+1,  MPI_COMM_WORLD, ierr )



            !write(GP_print_unit,'(A,4(1x,I6))') &
            !     'gplp:2 554 send   myid, sender, numsent, i_GP_individual ', &
            !                        myid, sender, numsent, i_GP_individual

            ! just sent a new task, so increment the number sent

            numsent = numsent + 1

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,4(1x,I6))') &
            !     'gplp:2 556  myid, sender, numsent, i_GP_individual ', &
            !                  myid, sender, numsent, i_GP_individual
            !endif ! L_GP_print


        else

            ! DONE !

            ! number of tasks sent out is >= number of individuals, so
            ! all the work has been completed

            ! tell the "sender" processor that it is done and
            ! send it a message to stop

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,3(1x,I6))') &
            !      'gplp:2 send msg to stop  myid, numsent, i_GP_individual ', &
            !                                myid, numsent, i_GP_individual
            !endif ! L_GP_print

            call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                           sender, 0,  MPI_COMM_WORLD, ierr )


        endif ! numsent

    enddo ! isource


    !----------------------------------------------------------------------

    ! this section takes care of the case where there are fewer GA individuals
    ! than (number of procs) - 1

    ! without the code below,  the program hangs because the processors
    ! with numbers  (n_GP_individuals+1)  to (numprocs-1)
    ! are waiting for a signal to stop
    ! and that is never going to be sent from the loop above.

    ! so when the above loop is finished, send a stop signal to the unused
    ! processors so the program can continue

    if( n_GP_individuals < numprocs -1 )then

        do  i = n_GP_individuals+1, numprocs-1

            call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                           i , 0,  MPI_COMM_WORLD, ierr )

        enddo ! i

    endif ! n_GP_individuals < numprocs -1

    !----------------------------------------------------------------------



else  ! not myid == 0


    ! code for processors 1 - ( numprocs - 1 )


    ! these processors wait until a message is received from
    ! processor 0 telling it to process
    ! the individual named in the message tag = MPI_STAT( MPI_TAG )

    !if( L_GP_print )then
    !    write(GP_print_unit,'(A,1x,I6)') &
    !          'gplp:3  myid, before recv_loop    ', myid
    !endif ! L_GP_print

    recv_loop:&
    do

        !if( myid == 1 )then
        !    write(GP_print_unit,'(A,1x,I6)') &
        !          'gplp: myid  proc1 in recv_loop', myid
        !endif



        call MPI_RECV( i_dummy, 1, MPI_INTEGER,    &
                       0, MPI_ANY_TAG,  MPI_COMM_WORLD, MPI_STAT, ierr )

        !---------------------------------------------------------------

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6))') &
        !      'gplp:3  myid, MPI_STAT( MPI_TAG ) ', &
        !                    myid, MPI_STAT( MPI_TAG )
        !endif ! L_GP_print

        ! was a stop signal received ?

        ! if the tag is <= 0, this is a stop signal

        if( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop

        !---------------------------------------------------------------

        ! process the individual named in the message tag

        i_2_individual = MPI_STAT( MPI_TAG )

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,4(1x,I6))') &
        !      'gplp:3  myid, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual ', &
        !               myid, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual
        !endif ! L_GP_print

        buffer = 0.0D0

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6),4x,L1)') &
        !      'gplp:3 myid, i_2_individual, Run_GA_lmdif(i_2_individual)', &
        !                   myid, i_2_individual, Run_GA_lmdif(i_2_individual)
        !endif ! L_GP_print

        !if( Run_GA_lmdif(i_2_individual)) then

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,2(1x,I6))') &
            !      'gplp:3 call setup_run_fcn  myid, i_2_individual   ', &
            !                                  myid, i_2_individual
            !endif ! L_GP_print

            ! do the Runge-Kutta integration for individual i_2_individual

            !t1 = MPI_Wtime()


            !call setup_run_fcn( i_2_individual, child_parameters, individual_quality )
            call setup_run_para_lmdif( i_2_individual, child_parameters, &
                                       individual_quality, &
                                       n_GP_individuals, GP_child_individual_SSE,  &
                                       L_GP_print, GP_print_unit )

            !t2 = MPI_Wtime()

            !write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
            !  'gplp: myid, i_2_individual, time spent in setup_run_fcn = ', &
            !         myid, i_2_individual, t2 - t1

            !if( L_GP_print )then
            !    write(GP_print_unit,'(A,3(1x,I6))') &
            !      'gplp:3 AFTER call setup_run_fcn  myid, i_2_individual   ', &
            !                                        myid, i_2_individual
            !endif ! L_GP_print

            !-------------------------------------------------------------------------

            buffer(1:n_parameters) = child_parameters(1:n_parameters,i_2_individual)
            buffer(n_parameters+1) = GP_Child_Individual_SSE(i_2_individual)
            buffer(n_parameters+2) = real( individual_quality(i_2_individual), kind=8 )

        !endif !  Run_GA_lmdif(i_2_individual)

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !    'gplp:3 705 myid, n_parameters, i_2_individual', &
        !                myid, n_parameters, i_2_individual
        !    write(GP_print_unit,'(A/(5(1x,E15.7)))') &
        !    'gplp:3 child_parameters(1:n_parameters,i_2_individual)', &
        !            child_parameters(1:n_parameters,i_2_individual)
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !    'gplp:3 myid, i_2_individual, individual_quality(i_2_individual)', &
        !            myid, i_2_individual, individual_quality(i_2_individual)
        !    write(GP_print_unit,'(A,2(1x,I6), 1x, E15.7)') &
        !    'gplp:3 send results myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)', &
        !                         myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)
        !endif ! L_GP_print



        ! send the R-K integration results for individual i_2_individual to processor 0

        call MPI_SEND( buffer, n_parameters+2, &
                       MPI_DOUBLE_PRECISION, 0, i_2_individual, MPI_COMM_WORLD, ierr )


        !---------------------------------------------------------------

        ! code to ensure that an error does not allow this loop to run forever

        nsafe = nsafe + 1

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I10))') &
        !      'gplp: myid,   nsafe = ', myid, nsafe
        !endif ! L_GP_print

        if( nsafe > 10 * n_GP_individuals ) then
            write(GP_print_unit,'(A,1x,I10)') &
              'gplp: too many iterations  nsafe = ', nsafe
            call MPI_FINALIZE(ierr)
            stop 'bad nsafe'
        endif ! nsafe

        !---------------------------------------------------------------

     enddo  recv_loop


endif ! myid == 0

!if( L_GP_print )then
!    write(GP_print_unit,'(A,2(1x,I6))') &
!      'gplp: after recv_loop  myid = ', myid
!endif ! L_GP_print


!-------------------------------------------------------------------

! wait until all n_GP_individuals have been processed

call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(A,2(1x,I6))') &
!      'gplp: after barrier 2 i_GP_generation, myid = ', &
!                             i_GP_generation, myid
!endif ! L_GP_print

!-------------------------------------------------------------------


! load the new child parameters into the population node parameters 


if( myid == 0 )then

    do  i_GP_individual = 1, n_GP_individuals
    
        n_parameters = 0
    
    
        do  i_tree=1,n_trees
            do  i_node=1,n_nodes
    
                if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) .eq. 0 ) then
    
                    n_parameters = n_parameters + 1 
                   
                    GP_population_node_parameters(i_node,i_tree,i_GP_individual) = &
                                child_parameters( n_parameters, i_GP_individual) 
    
                endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0
    
            enddo ! i_node
        enddo  ! i_tree
    
    
    enddo ! i_GP_individual

endif ! myid == 0

!-----------------------------------------------------------------------------

!
!
!!  calculate the fitness for this generation
!
!
!if( myid == 0  )then
!
!    if( L_GP_print )then
!
!        write(GP_print_unit,'(A)')  'gplp: GP_Child_Individual_SSE  '
!        write(GP_print_unit,'(5(1x,E15.7))')  GP_Child_Individual_SSE(1:n_GP_individuals)
!
!        write(GP_print_unit,'(/A,1x,I6)') &
!              'gplp: call calc_fitness i_GP_generation ', &
!                                       i_GP_generation
!    endif ! L_GP_print
!
!
!    call calc_fitness( child_parameters, individual_quality, &
!                       i_GA_Best_Parent, Parent_Parameters, L_stop_run, &
!                       i_GP_Generation, i_GP_individual )
!
!
!    if( L_GP_print )then
!        write(GP_print_unit,'(/A,1x,I6/)') &
!              'gplp: aft call calc_fitness i_GP_generation ', &
!                                                i_GP_generation
!    endif ! L_GP_print
!
!    !---------------------------------------------------------------------
!
!    ! this subroutine, if called, will run lmdif on the best individual
!    ! for this generation, instead of waiting to run lmdif after the
!    ! generation loop
!
!    !if( mod( i_GP_generation, 10 ) == 0 )then
!    !    if( L_GP_print )then
!    !        write(GP_print_unit,'(/A)') &
!    !          'gplp: call select_best_RK_lmdif_result '
!    !    endif ! L_GP_print
!    !    call select_best_RK_lmdif_result( &
!    !                i_GP_Generation,i_GP_individual, &
!    !                i_GA_best_parent, parent_parameters, &
!    !                child_parameters, &
!    !                individual_quality, L_stop_run  )
!    !endif ! mod( i_GP_generation, 10 ) == 0
!
!    !if( L_GP_print )then
!    !    write(GP_print_unit,'(A/)') &
!    !      'gplp: aft call select_best_RK_lmdif_result '
!    !endif ! L_GP_print
!
!    !---------------------------------------------------------------------
!
!
!
!endif ! myid == 0

!if( i_GP_generation > 1 )then
!    if( L_GP_print )then
!        write(GP_print_unit,'(A,2(1x,I6))') &
!          'gplp: at stop  i_GP_generation, myid = ', &
!                               i_GP_generation, myid
!        call MPI_FINALIZE(ierr) ! debug only
!        stop ! debug only
!    endif ! L_GP_print
!endif ! i_GP_generation > 1

!-------------------------------------------------------------------

!call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! necessary ?

!call MPI_BCAST( L_stop_run,  1,    &
!                MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )

!call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! necessary ?

!if( L_stop_run )then
!
!    !if( L_GP_print )then
!    !    write(GP_print_unit,'(/A,1x,I6/)') &
!    !    'gplp: L_stop_fitness is true &
!    !     &so last generation is i_GP_generation ', &
!    !                            i_GP_generation
!    !endif ! L_GP_print
!
!    i_GP_generation_last = i_GP_generation
!
!    exit
!
!endif ! L_stop_run



!if( L_GP_print )then
!    write(GP_print_unit,'(//A,1x,I6/)') &
!      'gplp:  finished the loop on i_GP_generation  myid =  ', myid
!endif ! L_GP_print

!----------------------------------------------------------------------

! wait until all processors have finished the generation loop

call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

!if( L_GP_print )then
!    write(GP_print_unit,'(A,1x,I6)') &
!          'gplp: after barrier 3 myid = ', myid
!endif ! L_GP_print

!----------------------------------------------------------------------
!
!
!! finished all generations,
!
!! now call lmdif on the best individual of the last generation
!! and determine if lmdif has improved the fitness of this individual
!! then save the parameters of the fitter of the two results, the RK result
!! and the lmdif result
!
!! GP_Individual_Node_Parameters loaded from the best parent parameters (RK or lmdif)
!! after the RK process and lmdif
!
!
!if( myid == 0  )then
!
!    if( L_GP_print )then
!        write(GP_print_unit,'(/A)') &
!          'gplp: call select_best_RK_lmdif_result '
!    endif ! L_GP_print
!
!    !t1 = MPI_Wtime()
!
!    call select_best_RK_lmdif_result( &
!                i_GP_Generation,i_GP_individual, &
!                i_GA_best_parent, parent_parameters, &
!                child_parameters, &
!                individual_quality, L_stop_run  )
!
!    !t2 = MPI_Wtime()
!
!    !write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
!    !      'gplp: i_GP_Generation,i_GP_individual, &
!    !       &time spent in select_best_RK_lmdif_result = ', &
!    !          i_GP_Generation,i_GP_individual, t2 - t1
!
!
!    if( L_GP_print )then
!        write(GP_print_unit,'(A/)') &
!          'gplp: aft call select_best_RK_lmdif_result '
!    endif ! L_GP_print
!
!
!endif ! myid == 0
!
!
!call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?
!
!------------------------------------------------------------------------

! broadcast individual_fitness


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast individual_fitness myid = ', myid
!endif ! L_GP_print

message_len = 1
call MPI_BCAST( individual_fitness, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast individual_fitness  ierr = ', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast Individual_SSE_best_parent


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast Individual_SSE_best_parent myid = ', myid
!endif ! L_GP_print

message_len = 1
call MPI_BCAST( Individual_SSE_best_parent, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast Individual_SSE_best_parent  ierr = ', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_Individual_Node_Parameters


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast GP_Individual_Node_Parameters  myid = ', myid
!endif ! L_GP_print

message_len = n_trees * n_nodes

call MPI_BCAST( GP_Individual_Node_Parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast GP_Individual_Node_Parameters  ierr = ', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_Individual_Initial_Conditions


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast GP_Individual_Initial_Conditions myid = ', myid
!endif ! L_GP_print

message_len = n_CODE_equations

call MPI_BCAST( GP_Individual_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast GP_Individual_Initial_Conditions ierr = ', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------



call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

!------------------------------------------------------------------------


if( myid == 0  )then
    if( L_GP_print )then
        write(GP_print_unit,'(//A/)') 'gplp:  final parent parameters  '
        write(GP_print_unit,'(A)') 'i_GP_individual                  parent_parameters '
        do  i_GP_individual = 1, n_GP_individuals
            write(GP_print_unit,'(I6,12(1x,E15.7 ))') &
              i_GP_individual, parent_parameters(1:n_parameters,i_GP_individual)
        enddo !  i_GP_individual
    endif ! L_GP_print
endif ! myid == 0


return


end subroutine GP_para_lmdif_process
