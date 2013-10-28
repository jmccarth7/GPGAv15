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

integer ::  isource
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
 dimension(n_maximum_number_parameters,n_GP_individuals) ::  child_parameters


real(kind=8) :: buffer2(n_maximum_number_parameters + 2)
real(kind=8) :: buffer2_recv(n_maximum_number_parameters + 2)


integer(kind=4) :: i
integer(kind=4) :: ii
integer(kind=4) :: jj

integer(kind=4) :: n_parms
integer(kind=4) :: n_parms_dim
integer(kind=4) :: nn
integer(kind=4) :: i_CODE_equation


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4), dimension( n_GP_individuals ) :: individual_quality

real(kind=8), external :: indiv_fitness


integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

logical :: L_GP_print

real(kind=8) ::  temp_SSE

!----------------------------------------------------------------------

L_GP_print = .FALSE.
if( i_GP_generation == 1 .or. &
    mod( i_GP_generation, GP_child_print_interval ) == 0 .or. & 
    i_GP_generation == n_GP_generations )then

    L_GP_print = .TRUE. 

endif ! i_GP_generation...

i_dummy = 0


buffer2(1:n_maximum_number_parameters+2)      = 0.0D0 
buffer2_recv(1:n_maximum_number_parameters+2) = 0.0D0 




if( myid == 0 )then
    if( L_GP_print )then
        write(GP_print_unit,'(//A)') 'gplp: at entry  '
        write(GP_print_unit,'(A,1x,E15.7)') 'gplp: dt ', dt
        !write(GP_print_unit,'(A,1x,I10)') &
        !      'gplp: n_parms    =   ', n_parms
        write(GP_print_unit,'(A,1x,I10)') &
              'gplp: n_GP_parameters =   ', n_GP_parameters
    endif ! L_GP_print

endif ! myid == 0


! jjm 20130417 >>>>>>>>>>>>>>>
!if( myid == 0) then
!
!    if( L_GP_print )then
!        write(GP_print_unit,'(A)')' '
!        write(GP_print_unit,'(A)') &
!                      'gplp: i_GP_indiv, i_tree, i_node, GP_pop_Node_Param'
!        do  i_GP_individual = 1, n_GP_individuals
!            do  i_tree=1,n_trees
!                do  i_node=1,n_nodes
!                    if( abs(GP_population_Node_Parameters(i_node,i_tree,i_GP_individual)) > 0.0d0 )then
!                        write(GP_print_unit,'(8x,3(1x,I6),3x,E15.7)') &
!                        i_GP_individual, i_tree, i_node, &
!                          GP_population_Node_Parameters(i_node,i_tree, i_GP_individual)
!                    endif ! abs( GP_pop_Node_Param(i_node,i_tree) ) > 0.0d0
!                enddo ! i_node
!            enddo  ! i_tree
!        enddo  ! i_GP_individual
!
!
!
!        write(GP_print_unit,'(/A)') &
!                     'gplp: i_GP_indiv, i_tree, i_node, GP_pop_Node_Type'
!
!        do  i_GP_individual = 1, n_GP_individuals
!            do  i_tree=1,n_trees
!                do  i_node=1,n_nodes
!                    if( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) > -9999 )then
!                        write(GP_print_unit,'(8x,4(1x,I6))') &
!                                i_GP_individual, i_tree, i_node, &
!                                GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual )
!                    endif ! GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) > -9999
!                enddo ! i_node
!            enddo  ! i_tree
!        enddo  ! i_GP_individual
!
!        write(GP_print_unit,'(A)')' '
!    endif ! L_GP_print
!
!endif ! myid == 0
!if( myid == 0) then
!
!    if( L_GP_print )then
!        write(GP_print_unit,'(A)')' '
!        write(GP_print_unit,'(A)') &
!                      'gplp: i_GP_indiv, i_tree, i_node, GP_pop_Node_Param'
!        do  i_GP_individual = 1, n_GP_individuals
!            do  i_tree=1,n_trees
!                do  i_node=1,n_nodes
!                    if( abs(GP_population_Node_Parameters(i_node,i_tree,i_GP_individual)) > 0.0d0 )then
!                        write(GP_print_unit,'(8x,3(1x,I6),3x,E15.7)') &
!                        i_GP_individual, i_tree, i_node, &
!                          GP_population_Node_Parameters(i_node,i_tree, i_GP_individual)
!                    endif ! abs( GP_pop_Node_Param(i_node,i_tree) ) > 0.0d0
!                enddo ! i_node
!            enddo  ! i_tree
!        enddo  ! i_GP_individual
!
!
!
!        write(GP_print_unit,'(/A)') &
!            'gplp: i_GP_indiv  i_tree  i_node  Node_Type    GP_pop_Node_Param'
!
!        do  i_GP_individual = 1, n_GP_individuals
!            do  i_tree=1,n_trees
!                do  i_node=1,n_nodes
!                    if( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) == 0 )then
!
!                        if( abs(GP_population_Node_Parameters(i_node,i_tree,i_GP_individual)) > 0.0d0 )then
!                            write(GP_print_unit,'(8x,4(1x,I6),9x,E15.7)') &
!                            i_GP_individual, i_tree, i_node, &
!                              GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ), &
!                              GP_population_Node_Parameters(i_node,i_tree, i_GP_individual)
!                        else
!
!                            write(GP_print_unit,'(8x,4(1x,I6))') &
!                                    i_GP_individual, i_tree, i_node, &
!                                    GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual )
!
!                        endif ! abs( GP_pop_Node_Param(i_node,i_tree) ) > 0.0d0
!
!                    elseif( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) > -9999 )then
!                        write(GP_print_unit,'(8x,4(1x,I6))') &
!                                i_GP_individual, i_tree, i_node, &
!                                GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual )
!                    endif ! GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) == 0       
!
!                enddo ! i_node
!            enddo  ! i_tree
!        enddo  ! i_GP_individual
!
!        write(GP_print_unit,'(A)')' '
!    endif ! L_GP_print
!
!endif ! myid == 0
! jjm 20130417 <<<<<<<<<<<<<<<


!-----------------------------------------------------------------------------


! load the population node parameters into the child parameters


child_parameters( 1:n_maximum_number_parameters, 1:n_GP_individuals) = 0.0d0

nn = 0

if( myid == 0 )then

    !write(GP_print_unit,'(/A)' ) &
    ! 'gplp: myid, i_GP_individual, n_parms, child_parameters, GP_pop_init_cond/GP_pop_node_params'

    do  i_GP_individual = 1, n_GP_individuals

        nn = 0

        do  i_CODE_equation=1,n_CODE_equations

            nn = nn + 1
            child_parameters( nn, i_GP_individual) =  &
                GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual)

            !write(GP_print_unit,'(I10,1x,I10,1x,I10,2(6x,E15.7))') &
            !      myid, i_GP_individual, nn, &
            !      child_parameters(nn,i_GP_individual), &
            !      GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual)

        enddo  ! i_CODE_equation

        do  i_tree=1,n_trees
            do  i_node=1,n_nodes

                if( GP_Adult_population_Node_Type(i_Node,i_Tree, i_GP_individual ) == 0 ) then

                    nn = nn + 1
                    child_parameters( nn, i_GP_individual) =  &
                         GP_population_node_parameters(i_node,i_tree,i_GP_individual)

                    !write(GP_print_unit,'(I10,1x,I10,1x,I10,2(6x,E15.7))') &
                    !      myid, i_GP_individual, nn, &
                    !      child_parameters(nn,i_GP_individual), &
                    !      GP_population_node_parameters(i_node,i_tree,i_GP_individual)

                endif ! GP_Individual_Node_Type(i_Node,i_Tree) == 0

            enddo ! i_node
        enddo  ! i_tree

        GP_n_parms( i_GP_individual ) = nn

    enddo ! i_GP_individual

endif ! myid == 0


!if( myid == 0  )then
!    if( L_GP_print )then
!        write(GP_print_unit,'(//A/)') 'gplp:  GP_n_parms '
!        write(GP_print_unit,'(A)') &
!              'i_GP_individual        GP_n_parms '
!        do  i_GP_individual = 1, n_GP_individuals
!            write(GP_print_unit,'(I10,10x,I10)') &
!                 i_GP_individual,  GP_n_parms(i_GP_individual) 
!        enddo !  i_GP_individual
!    endif ! L_GP_print
!endif ! myid == 0

!-----------------------------------------------------------------------------

if( myid == 0  )then
    if( L_GP_print )then

        write(GP_print_unit,'(//A/)') 'gplp:  initial child parameters  '
        !write(GP_print_unit,'(A)') &
        !      'i_GP_individual                  child_parameters '

        do  i_GP_individual = 1, n_GP_individuals

            write(GP_print_unit,'(/A)') &
              'gplp:  indiv   equation   GP_Population_Initial_Conditions'

            do  ii=1,n_CODE_equations
                write(GP_print_unit,'(I10,1x,I10,3x, E20.10)') &
                 i_GP_individual,  ii, &
                 GP_Population_Initial_Conditions(ii, i_GP_individual)
            enddo ! ii

            write(GP_print_unit,'(/A)') &
              'gplp:  indiv   parameter  child_parameters'

            n_parms = GP_n_parms( i_GP_individual )

            do  nn= n_code_equations + 1, n_parms
                write(GP_print_unit,'(I10,1x,I10,3x,E20.10)') &
                 i_GP_individual, nn, &
                 child_parameters(nn,i_GP_individual)
            enddo ! nn

        enddo !  i_GP_individual
    endif ! L_GP_print
endif ! myid == 0



!-----------------------------------------------------------------------------

! set up MPI process



!L_stop_run  = .FALSE.
!L_stop_run  = .TRUE.




!------------------------------------------------------------------------

!  broadcast GP_n_parms


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,2(1x,I6))') &
!    'gplp:  broadcast GP_N_parms  myid, i_GP_generation ', &
!                                  myid, i_GP_generation
!endif ! L_GP_print


!if( L_GP_print )then
!    write(GP_print_unit,'(A,4(1x,I6)/)') &
!     'gplp:  myid, n_GP_individuals =', &
!             myid, n_GP_individuals
!endif ! L_GP_print

call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

call MPI_BCAST( GP_n_parms,  n_GP_individuals,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?

if( L_GP_print )then
    if( myid == 0 )then

        write(GP_print_unit,'(/A,2(1x,I10)/)') &
         'gplp: child  broadcast myid, ierr = ', myid, ierr

        write(GP_print_unit,'(/A)') &
         'gplp: myid, i_GP_individual  GP_N_parms(i_GP_individual)  '

        do  i_GP_individual = 1, n_GP_individuals
            write(GP_print_unit,'(I10,1x,I10,1x,I10)') &
              myid, i_GP_individual, GP_N_parms(i_GP_individual)
        enddo ! i_GP_individual
    endif ! myid == 0
endif ! L_GP_print



!------------------------------------------------------------------------

!  broadcast child parameters


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,2(1x,I6))') &
!    'gplp:  broadcast child parameters myid, i_GP_generation', &
!                                       myid, i_GP_generation
!endif ! L_GP_print

child_number =  n_GP_individuals * n_maximum_number_parameters

!if( L_GP_print )then
!    write(GP_print_unit,'(A,4(1x,I6)/)') &
!    'gplp: myid, n_GP_individuals, n_maximum_number_parameters, child_number =', &
!           myid, n_GP_individuals, n_maximum_number_parameters, child_number
!endif ! L_GP_print


call MPI_BCAST( Child_Parameters,  child_number,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

call MPI_BARRIER( MPI_COMM_WORLD, ierr )  ! necessary ?



!------------------------------------------------------------------------


if( myid == 0 )then
    if( L_GP_print )then
        write(GP_print_unit,'(/A,1x,I6)') &
         'gplp: begin parallel lmdif segment i_GP_generation', &
                                             i_GP_generation
        write(GP_print_unit,'(A,1x,I6/)') &
         'gplp: n_GP_individuals', &
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
    !     'gplp: aft source loop 1 myid, i_GP_individual, numsent', &
    !                              myid, i_GP_individual, numsent
    !endif ! L_GP_print

    !-------------------------------------------------------------------------------------

    ! processor 0 loops over the number of individuals and waits for a message
    ! from the other processors

    do  isource = 1, n_GP_individuals

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 521 myid, isource, numsent', &
        !                 myid, isource, numsent
        !endif ! L_GP_print


        buffer2_recv = 0.0d0  
        call MPI_RECV( buffer2_recv, n_maximum_number_parameters+2, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                       MPI_COMM_WORLD, MPI_STAT,  ierr )

        sender       = MPI_STAT( MPI_SOURCE )
        i_individual = MPI_STAT( MPI_TAG )


        n_parms = GP_n_parms( i_individual )


        ! received a message from processor "sender" which processed
        ! individual "i_individual"


        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,5(1x,I6))') &
        !     'gplp:2 529 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )', &
        !                 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )
        !    write(GP_print_unit,'(A,5(1x,I6))') &
        !     'gplp:2 529 myid, isource, numsent, sender, i_individual ', &
        !                 myid, isource, numsent, sender, i_individual
        !endif ! L_GP_print


        ! store the information received in the above message


        child_parameters(1:n_maximum_number_parameters,i_individual) =  &
                             buffer2_recv(1:n_maximum_number_parameters)

        GP_Child_Individual_SSE(i_individual) =  &
                             buffer2_recv( n_maximum_number_parameters+1)
        individual_quality(i_individual) = &
                       nint( buffer2_recv( n_maximum_number_parameters+2) )


        !if( L_GP_print .and. i_individual == 3 )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 554 myid, n_parms, i_individual', &
        !                 myid, n_parms, i_individual
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 554 myid, i_individual, GP_n_parms( i_individual )', &
        !                 myid, i_individual, GP_n_parms( i_individual )
        !    write(GP_print_unit,'(A/(5(1x,E15.7)))') &
        !     'gplp:2 child_parameters(1:n_parms,i_individual)', &
        !             child_parameters(1:n_parms,i_individual)
        !    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
        !     'gplp:2 myid, i_individual, GP_Child_Individual_SSE(i_individual)', &
        !             myid, i_individual, GP_Child_Individual_SSE(i_individual)
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:2 myid, i_individual, individual_quality(i_individual)', &
        !             myid, i_individual, individual_quality(i_individual)
        !endif ! L_GP_print


        !-------------------------------------------------------------------------------------

        ! check to see if all individuals have been processed

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,1x,I6, 4x,L1)') &
        !     'gplp:2 542 myid, numsent < n_GP_individuals', &
        !                 myid, numsent < n_GP_individuals
        !endif ! L_GP_print

        if( numsent <  n_GP_individuals )then

            ! numsent <  n_GP_individuals
            ! means not all individuals have been processed

            ! send a message to the processor "sender"
            ! which just sent a message saying it has
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

            !     DONE !

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

    ! this section takes care of the case where there are fewer GP individuals
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



    recv_loop:&
    do


        call MPI_RECV( i_dummy, 1, MPI_INTEGER,    &
                       0, MPI_ANY_TAG,  MPI_COMM_WORLD, MPI_STAT, ierr )

        !---------------------------------------------------------------

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6))') &
        !      'gplp:3  myid, MPI_STAT( MPI_TAG ) ', &
        !               myid, MPI_STAT( MPI_TAG )
        !endif ! L_GP_print


        ! was a stop signal received ?

        ! if the tag is <= 0, this is a stop signal

        if( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop


        !---------------------------------------------------------------


        ! process the individual named in the message tag


        i_2_individual = MPI_STAT( MPI_TAG )

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,4(1x,I6))') &
        !      'gplp:3 myid, i_dummy, MPI_STAT(MPI_TAG), i_2_individual', &
        !              myid, i_dummy, MPI_STAT(MPI_TAG), i_2_individual
        !endif ! L_GP_print

        buffer2 = 0.0D0


        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,2(1x,I6))') &
        !      'gplp:3 call setup_run_para_lmdif  myid, i_2_individual', &
        !                                         myid, i_2_individual
        !endif ! L_GP_print


        !-------------------------------------------------------------------------


        ! do the Runge-Kutta integration for individual i_2_individual


        n_parms = GP_n_parms( i_2_individual )

        n_parms_dim = max( 1, n_parms )

        !if( L_GP_print .and. myid == 1  )then
        !    write(GP_print_unit,'(A,2(1x,I6),1x,I10)') &
        !      'gplp:3 myid, i_2_individual, n_parms', &
        !              myid, i_2_individual, n_parms
        !endif ! L_GP_print

        !if( L_GP_print .and. i_2_individual == 3 )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:6 554 myid, n_parms, i_2_individual', &
        !                 myid, n_parms, i_2_individual
        !    write(GP_print_unit,'(A/(5(1x,E15.7)))') &
        !     'gplp:6 child_parameters(1:n_parms,i_2_individual)', &
        !             child_parameters(1:n_parms,i_2_individual)
        !    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
        !     'gplp:6 myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)', &
        !             myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:6 myid, i_2_individual, individual_quality(i_2_individual)', &
        !             myid, i_2_individual, individual_quality(i_2_individual)
        !endif ! L_GP_print



        temp_SSE = GP_child_individual_SSE(i_2_individual)

        call setup_run_para_lmdif( i_2_individual, &
                                   child_parameters(1,i_2_individual), &
                                   individual_quality(i_2_individual), &
                                   n_GP_individuals, &
                                   temp_SSE,  &
                                   n_parms, n_parms_dim, &
                                   L_GP_print, GP_print_unit )


        !--------------------------------------------------------------------------

        ! don't replace original child SSE if lmdif SSE indicates a bad result

        if( abs( temp_SSE ) < 1.0D12 )then
            
            GP_child_individual_SSE(i_2_individual) = temp_SSE

        endif ! abs(temp_SSE)...

        !--------------------------------------------------------------------------

        !if( L_GP_print .and. i_2_individual == 3 )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:7 723 myid, n_parms, i_2_individual  AFTER LMDIF ', &
        !                 myid, n_parms, i_2_individual
        !    write(GP_print_unit,'(A/(5(1x,E15.7)))') &
        !     'gplp:7 child_parameters(1:n_parms,i_2_individual)', &
        !             child_parameters(1:n_parms,i_2_individual)
        !    write(GP_print_unit,'(A,2(1x,I6),1x,E15.7)') &
        !     'gplp:7 myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)', &
        !             myid, i_2_individual, GP_Child_Individual_SSE(i_2_individual)
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !     'gplp:7 myid, i_2_individual, individual_quality(i_2_individual)', &
        !             myid, i_2_individual, individual_quality(i_2_individual)
        !endif ! L_GP_print

        !if( L_GP_print )then
        !    write(GP_print_unit,'(A,3(1x,I6))') &
        !      'gplp:3 AFTER call setup_run_para_lmdif  myid, i_2_individual', &
        !                                               myid, i_2_individual
        !endif ! L_GP_print

        !-------------------------------------------------------------------------

        n_parms = GP_n_parms( i_2_individual )

        buffer2(1:n_maximum_number_parameters) = &
                 child_parameters(1:n_maximum_number_parameters,i_2_individual)
        buffer2(n_maximum_number_parameters+1) = GP_Child_Individual_SSE(i_2_individual)
        buffer2(n_maximum_number_parameters+2) = &
            real( individual_quality(i_2_individual), kind=8 )


        ! send the R-K integration results for individual i_2_individual to processor 0

        call MPI_SEND( buffer2, n_maximum_number_parameters+2, &
                       MPI_DOUBLE_PRECISION, 0, i_2_individual, MPI_COMM_WORLD, ierr )


        !---------------------------------------------------------------

        ! code to ensure that an error does not allow this loop to run forever

        nsafe = nsafe + 1

        if( nsafe > 100 * n_GP_individuals ) then
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

        n_parms = 0

        do  i_CODE_equation=1,n_CODE_equations
            n_parms = n_parms + 1
            GP_Population_Initial_Conditions(i_CODE_Equation, i_GP_Individual) = &
                                   child_parameters( n_parms, i_GP_individual)
        enddo  ! i_CODE_equation

        do  i_tree=1,n_trees
            do  i_node=1,n_nodes

                if( GP_Adult_population_Node_Type(i_Node,i_Tree,i_GP_individual ) == 0 ) then
                    n_parms = n_parms + 1
                    GP_population_node_parameters(i_node,i_tree,i_GP_individual) = &
                                child_parameters( n_parms, i_GP_individual)
                endif ! GP_Individual_Node_Type(i_Node,i_Tree) .eq. 0

            enddo ! i_node
        enddo  ! i_tree


    enddo ! i_GP_individual


    !-----------------------------------------------------------------------------
    ! moved to myid=0 section 20131022 jjm

    GP_Adult_Population_Parameter_Solution(1:n_maximum_number_parameters,1:n_GP_individuals) = &
                      child_parameters(1:n_maximum_number_parameters,1:n_GP_individuals)

    !-----------------------------------------------------------------------------

endif ! myid == 0

!if( i_GP_generation > 1 )then
!    if( L_GP_print )then
!        write(GP_print_unit,'(A,2(1x,I6))') &
!          'gplp: at stop  i_GP_generation, myid =', &
!                          i_GP_generation, myid
!        call MPI_FINALIZE(ierr) ! debug_only
!        stop ! debug_only
!    endif ! L_GP_print
!endif ! i_GP_generation > 1


!----------------------------------------------------------------------


! wait until all processors have finished the generation loop

call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

!if( L_GP_print )then
!    write(GP_print_unit,'(A,1x,I6)') &
!          'gplp: after barrier 3 myid =', myid
!endif ! L_GP_print

!------------------------------------------------------------------------

! broadcast individual_quality


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast individual_quality myid =', myid
!endif ! L_GP_print

message_len = n_GP_individuals
call MPI_BCAST( individual_quality, message_len,    &
                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast individual_quality  ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_Child_Individual_SSE


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast Individual_SSE_best_parent myid =', myid
!endif ! L_GP_print

message_len = n_GP_individuals
call MPI_BCAST( GP_Child_Individual_SSE, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast Individual_SSE_best_parent  ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_population_node_parameters


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast GP_population_node_parameters myid =', myid
!endif ! L_GP_print

message_len = n_trees * n_nodes * n_GP_individuals

call MPI_BCAST( GP_population_node_parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast GP_population_node_parameters  ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_Population_Initial_Conditions


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast GP_Population_Initial_Conditions   myid =', myid
!endif ! L_GP_print

message_len = n_CODE_equations  * n_GP_individuals

call MPI_BCAST( GP_Population_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast GP_Population_Initial_Conditions ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------

! broadcast GP_Adult_Population_Parameter_Solution


!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: broadcast GP_Adult_Population_Parameter_Solution  myid =', myid
!endif ! L_GP_print

message_len = n_maximum_number_parameters * n_GP_individuals

call MPI_BCAST( GP_Adult_Population_Parameter_Solution, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_GP_print )then
!    write(GP_print_unit,'(/A,1x,I6)') &
!     'gplp: aft broadcast GP_Adult_Population_Parameter_Solution ierr =', ierr
!endif ! L_GP_print


!------------------------------------------------------------------------


call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

!------------------------------------------------------------------------


if( myid == 0  )then

    if( L_GP_print )then

        write(GP_print_unit,'(//A)') 'gplp:  final child SSE'
        write(GP_print_unit,'(/A)') &
              'gplp:  indiv   GP_Child_Individual_SSE '
        do  i_GP_individual = 1, n_GP_individuals
            write(GP_print_unit,'(4x,I6,5x, E20.10)') &
                  i_GP_individual, &
                 GP_Child_Individual_SSE(i_GP_individual)

        enddo ! i_GP_individual

        write(GP_print_unit,'(//A)') 'gplp:  final child parameters  '

        do  i_GP_individual = 1, n_GP_individuals

            write(GP_print_unit,'(/A)') &
              'gplp:  indiv   equation  GP_Population_Initial_Conditions'

            do  ii=1,n_CODE_equations
                write(GP_print_unit,'(4x,I6,1x,I6,6x, E20.10)') &
                 i_GP_individual,  ii, &
                 GP_Population_Initial_Conditions(ii, i_GP_individual)
            enddo ! ii

            write(GP_print_unit,'(/A)') &
              'gplp:  indiv   parameter  child_parameters '

            n_parms = GP_n_parms( i_GP_individual )

            do  nn= n_code_equations + 1, n_parms
                write(GP_print_unit,'(4x,I6,1x,I6,6x,E20.10)') &
                 i_GP_individual, nn, &
                 child_parameters(nn,i_GP_individual)
            enddo ! nn

        enddo !  i_GP_individual

    endif ! L_GP_print

endif ! myid == 0


return


end subroutine GP_para_lmdif_process
