subroutine GPCODE_GA_lmdif_Parameter_Optimization( )

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

integer(kind=4) :: individual_quality(n_GA_individuals)

real(kind=8), external :: indiv_fitness

logical :: L_stop_run


integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: jj
integer(kind=4) :: i_parameter

!----------------------------------------------------------------------

buffer(1:n_maximum_number_parameters+2)      = 0.0D0
buffer_recv(1:n_maximum_number_parameters+2) = 0.0D0



n_parameters = n_GP_parameters



if( myid == 0 )then
    write(GA_print_unit,'(//A)') 'GP_GA_opt: at entry  '
    write(GA_print_unit,'(A,1x,E15.7)') 'GP_GA_opt: dt ', dt
    write(GA_print_unit,'(A,1x,I10)') &
          'GP_GA_opt: n_parameters    =   ', n_parameters
    write(GA_print_unit,'(A,1x,I10)') &
          'GP_GA_opt: n_GP_parameters =   ', n_GP_parameters

endif ! myid == 0


! jjm 20130417 >>>>>>>>>>>>>>>
if( myid == 0) then

    write(GA_print_unit,'(A)')' '
    do  i_tree=1,n_trees
        do  i_node=1,n_nodes
            if( abs( GP_Individual_Node_Parameters(i_node,i_tree) ) > 1.0e-20 )then
                write(GA_print_unit,'(A,2(1x,I6),1x,E15.7)') &
                  'GP_GA_opt: i_node, i_tree, GP_Indiv_Node_Param', &
                              i_node, i_tree, GP_Individual_Node_Parameters(i_node,i_tree)
            endif ! abs( GP_Indiv_Node_Param(i_node,i_tree) ) > 1.0e-20
        enddo ! i_node
    enddo  ! i_tree

    write(GA_print_unit,'(A)')' '


    do  i_tree=1,n_trees
        do  i_node=1,n_nodes
            if( GP_Individual_Node_Type(i_node,i_tree) > -9999 )then
                write(GA_print_unit,'(A,3(1x,I6))') &
                 'GP_GA_opt: i_node, i_tree, GP_Indiv_Node_Type', &
                             i_node, i_tree, GP_Individual_Node_Type(i_node,i_tree)
            endif ! GP_Indiv_Node_Type(i_node,i_tree) > -9999
        enddo ! i_node
    enddo  ! i_tree

    write(GA_print_unit,'(A)')' '

endif ! myid == 0
! jjm 20130417 <<<<<<<<<<<<<<<




if( n_parameters .le. 0) then
    write(GA_print_unit,'(A)')        &
          'GP_GA_opt: ERROR: n_parameters </= 0'
    write(GA_print_unit,'(A,1x,I10)') &
          'GP_GA_opt: n_parameters =   ', n_parameters
    stop 'n_par<=0'
endif

if( n_time_steps .lt. n_parameters) then
    write(GA_print_unit,'(A)') &
          'GP_GA_opt: ERROR: n_time_steps < n_parameters'
    write(GA_print_unit,'(A,1x,I10)') &
          'GP_GA_opt:  n_time_steps = ', n_time_steps
    write(GA_print_unit,'(A,1x,I10)') &
          'GP_GA_opt: n_parameters =   ', n_parameters
    stop 'n_time < n_par'
endif

if( tol .lt. zero) then
    write(GA_print_unit,'(A)') 'GP_GA_opt: ERROR: tol < 0'
    stop 'tol < 0.0'
endif


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!CALL RANDOM_SEED(size = n_seed)
!
!if( .not. allocated( seed ) )then
!
!    ALLOCATE(seed(n_seed))
!
!endif ! .not. allocated( seed )
!
!CALL SYSTEM_CLOCK(COUNT=clock)
!seed = clock + 37 * (/ (i_seed - 1, i_seed = 1, n_seed) /)
!
!CALL RANDOM_SEED(PUT = seed)
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


if( GA_Crossover_Probability+GA_Mutation_Probability .gt. 1.0d0 ) then
    if( myid == 0 )then
        write(GA_print_unit,'(A)') &
         'GP_GA_opt: Sum of Crossover and Mutation Probabilities are too high'
    endif ! myid == 0
endif !   GA_Crossover_Probability+GA_Mutation_Probability .gt. 1.0d0


! calculate the number of GA Crossovers
n_GA_Crossovers = nint(GA_Crossover_Probability * n_GA_individuals)

! calculate the number of GA Mutations
n_GA_Mutations  = nint(GA_Mutation_Probability  * n_GA_individuals)


! calculate the number of GA elites
n_GA_save_elites = nint(GA_save_elites_Probability  * n_GA_individuals)


! calculate the number of GA random replacements 
n_GA_rand_replaces  = nint(GA_rand_replace_Probability  * n_GA_individuals)



if( myid == 0 )then
    write(GA_print_unit,'(A,1x,I6)') &
          'GP_GA_opt: n_GA_Crossovers  ', n_GA_Crossovers
    write(GA_print_unit,'(A,1x,I6)') &
          'GP_GA_opt: n_GA_Mutations   ', n_GA_Mutations
    write(GA_print_unit,'(A,1x,I6)') &
          'GP_GA_opt: n_GA_save_elites ', n_GA_save_elites
    write(GA_print_unit,'(A,1x,I6)') &
          'GP_GA_opt: n_GA_rand_replaces   ', n_GA_rand_replaces
    write(GA_print_unit,'(A,1x,E12.5)') &
          'GP_GA_opt: GA_rand_replace_Probability ', GA_rand_replace_Probability
endif ! myid == 0


child_parameters( 1:n_maximum_number_parameters, 1:n_GA_individuals) = 0.0d0


!-----------------------------------------------------------------------------

! set up MPI process 


if( myid == 0 )then

    write(GA_print_unit,'(/A,2(1x,I6))')&
          'GP_GA_opt: n_GA_individuals, numprocs ', &
                      n_GA_individuals, numprocs
    write(GA_print_unit,'(A,1x,I6/)')&
          'GP_GA_opt: n_parameters ', n_parameters

endif ! myid == 0

!-----------------------------------------------------------------------------

L_stop_run  = .FALSE.

do  i_GA_generation=1,n_GA_Generations

    ! Run_GA_lmdif determines if the new child
    ! has to be sent to lmdif for 'local' optimization

    Run_GA_lmdif=.false.

    if( myid == 0 )then
        write(GA_print_unit,'(/A,1x,I6,1x,A/)') &
              'GA Generation ',i_GA_generation,' is underway'
    endif ! myid == 0



    if( myid == 0 )then

        if( i_GA_generation .eq. 1 ) then

            ! on the first generation,
            ! randomly create the initial individual parameter arrays
            ! for each individual

            call Initialize_GA_Child_Parameters(Child_Parameters)

            Run_GA_lmdif=.true.

            ! print child parameters at start of the generation

            write(GA_print_unit,'(/A,1x,I6)') &
            'GP_GA_opt: child parameters at start of generation: ', &
                                                i_GA_generation
            do  i_GA_individual = 1, n_GA_Individuals
                write(GA_print_unit,'(I6,1x,12(1x,E15.7))') &
                      i_GA_individual, &
                      ( child_parameters(jj,i_GA_individual), jj = 1,n_parameters )
            enddo ! i_GA_individual



        else  ! i_GA_generation > 1


            ! create the second 'generation' of parameter estimates using either:

            !    i) replace individuals with bad RK process results
            !   ii) 'Fitness-Proportionate Reproduction;
            !  iii) GA Crossover;
            !   iv) GA Mutation
            !    v) save elites from last generation from being changed


            !-------------------------------------------------------------------------------

            !   save the most fit individuals for the next generation

            write(GA_print_unit,'(/A)')'GP_GA_opt: call GA_save_elites '

            call GA_save_elites( ) 

            !-------------------------------------------------------------------------------

            !  replace the parameters of any individual with quality < 0 with new
            !  random numbers

            write(GA_print_unit,'(/A)')'GP_GA_opt: call GA_replace_bad_individuals  '
            call GA_replace_bad_individuals(Child_Parameters, individual_quality )

            !-------------------------------------------------------------------------------

            !   do initial "GA Fitness-Proportionate Reproduction"
            !   to create a new population of children for all n_GA_individual

            write(GA_print_unit,'(/A)')&
                  'GP_GA_opt: call GA_Fitness_Proportionate_Reproduction '
            !write(GA_print_unit,'(/A)')&
            !     'GP_GA_opt: SKIP GA_Fitness_Proportionate_Reproduction '


            !call system_clock( count=clock1, count_rate=ratec, count_max= maxclk)

            call GA_Fitness_Proportionate_Reproduction(&
                            Parent_Parameters,Child_Parameters, &
                                                individual_quality )

            !call system_clock( count=clock2, count_rate=ratec, count_max= maxclk)
            !!write(GA_print_unit,*) 'GP_GA_opt: clock1,clock2,ratec,maxclk ', &
            !!                                   clock1,clock2,ratec,maxclk
            !write(GA_print_unit,*) &
            !  'GP_GA_opt: GA_Fitness_Proportionate_Reproduction time = ', &
            !            real(clock2-clock1,kind=4)/real(ratec,kind=4) , ' seconds'

            !-------------------------------------------------------------------------------

            !   do "GA Parameter Crossover" Operations Using Tournament-Style Selection
            !   and randomly use it to replace the parents with the new children

            if( n_GA_Crossovers .gt. 0) then

                write(GA_print_unit,'(/A,1x,I6)')&
                      'GP_GA_opt: call GA_Tournament_Style_Sexual_Repro  n_GA_Crossovers',  &
                                                                         n_GA_Crossovers
                call GA_Tournament_Style_Sexual_Reproduction( &
                            Parent_Parameters, Child_Parameters, individual_quality )

            endif !   n_GA_Crossovers .gt. 0


            !-------------------------------------------------------------------------------

            !   do "GA Parameter Mutation" Operations

            !   select a random individual and put a new random number into one of
            !   its parameters

            if( n_GA_Mutations .gt. 0) then

                write(GA_print_unit,'(/A,1x,I6)')&
                      'GP_GA_opt: call GA_Mutations  n_GA_Mutations',  n_GA_Mutations

                call GA_Mutations( Child_Parameters, individual_quality )

            endif !   n_GA_Mutations .gt. 0


            !-------------------------------------------------------------------------------

            !   do "GA Parameter rand_replace" Operations

            !   select a random, non-elite individual and put new random numbers into 
            !   its parameters

            if( n_GA_rand_replaces > 0) then
                write(GA_print_unit,'(/A,1x,I6)')&
                      'GP_GA_opt: call GA_rand_replace  n_GA_rand_replaces',  &
                                                        n_GA_rand_replaces
                call GA_random_replace( Child_Parameters, individual_quality )
            endif !   n_GA_rand_replaces .gt. 0

            !-------------------------------------------------------------------------------


            ! print child parameters at start of the generation

            !if( i_GA_generation == n_GA_generations               .or.  &
            !    mod( i_GA_generation, GA_child_print_interval ) == 0      )then

                write(GA_print_unit,'(/A,1x,I6)') &
                'GP_GA_opt: child parameters at start of generation: ', &
                                                    i_GA_generation

                do  i_GA_individual = 1, n_GA_Individuals
                    write(GA_print_unit,'(I6,1x,12(1x,E15.7))') &
                          i_GA_individual, &
                         (child_parameters(jj, i_GA_individual),jj = 1,n_parameters )
                enddo ! i_GA_individual

            !endif ! i_GA_generation == n_GA_generations ...



        endif ! i_GA_generation .eq. 1

    endif ! myid == 0


    !------------------------------------------------------------------------

    !  broadcast child parameters


    !write(GA_print_unit,'(/A,2(1x,I6))') &
    ! 'GP_GA_opt:  broadcast child parameters myid, i_GA_generation ', myid, i_GA_generation

    child_number =  n_GA_Individuals * n_maximum_number_parameters

    !write(GA_print_unit,'(A,4(1x,I6)/)') &
    ! 'GP_GA_opt:  myid, n_GA_Individuals, n_maximum_number_parameters, child_number = ', &
    !              myid, n_GA_Individuals, n_maximum_number_parameters, child_number

    call MPI_BCAST( Child_Parameters,  child_number,    &
                    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    !write(GA_print_unit,'(/A,2(1x,I10)/)') &
    ! 'GP_GA_opt: child  broadcast myid, ierr = ', myid, ierr

    !write(GA_print_unit,'(/A,2(1x,I10)/)') &
    ! 'GP_GA_opt: myid, n_GA_Individuals = ', myid, n_GA_Individuals
    !do  i_GA_individual = 1, n_GA_Individuals
    !    write(GA_print_unit,'(I6,1x,I6,1x,12(1x,E15.7))') &
    !          myid, i_GA_individual, &
    !          ( child_parameters(jj,i_GA_individual), jj = 1,n_parameters )
    !enddo ! i_GA_individual


    !------------------------------------------------------------------------

    ! broadcast Run_GA_lmdif

    !write(GA_print_unit,'(/A,1x,I6)') &
    !      'GP_GA_opt:  broadcast Run_GA_lmdif i_GA_generation ', i_GA_generation

    call MPI_BCAST( Run_GA_lmdif,  n_GA_Individuals,    &
                        MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )

    !write(GA_print_unit,'(A,1x,I10/)') &
    ! 'GP_GA_opt: Run_GA_lmdif  broadcast ierr = ', ierr


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !   begin RK fcn integration segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    if( myid == 0 )then
        write(GA_print_unit,'(/A,1x,I6/)') &
              'GP_GA_opt: begin RK fcn integration segment i_GA_generation ', &
                                                           i_GA_generation
        write(GA_print_unit,'(A,1x,I6/)') 'GP_GA_opt: n_GA_individuals ', &
                                                      n_GA_individuals
    !else
    !    write(GA_print_unit,'(A,2(1x,I6)/)') 'GP_GA_opt: myid, n_GA_individuals ', &
    !                                                     myid, n_GA_individuals

    endif !  myid == 0


    individual_quality( 1: n_GA_individuals ) = 1

    nsafe = 0


    if( myid == 0  )then

        ! processor 0 sends job assignments to N  processors 
        ! where N is the smaller of (number of processors -1) 
        !                    and the number of individuals 

        ! numsent is the number of messages sent up to now

        numsent = 0 
        i_GA_individual = 0

        do  isource = 1, min( numprocs-1, n_GA_individuals ) 

            
            i_GA_individual = i_GA_individual + 1

            !write(GA_print_unit,'(A,1x,I6, 4x, L1)') &
            ! 'GP_GA_opt:1 494 i_GA_individual, Run_GA_lmdif(i_GA_individual)  ', &
            !                  i_GA_individual, Run_GA_lmdif(i_GA_individual)


            call MPI_SEND( i_dummy,  1, MPI_INTEGER,    &
                           isource, isource,  MPI_COMM_WORLD, ierr )                          
            numsent = numsent + 1 

            !write(GA_print_unit,'(A,4(1x,I6))') &
            !     'GP_GA_opt:1 504 myid, isource, i_GA_individual, numsent ', &
            !                      myid, isource, i_GA_individual, numsent
                
        enddo ! isource


        ! at this point i_GA_individual = numsent 


        !write(GA_print_unit,'(A,4(1x,I6))') &
        !         'GP_GA_opt: aft source loop 1 myid, i_GA_individual, numsent ', &
        !                                       myid, i_GA_individual, numsent

        !-------------------------------------------------------------------------------------

        ! processor 0 loops over the number of individuals and waits for a message
        ! from the other processors 

        do  isource = 1, n_GA_individuals 

            !write(GA_print_unit,'(A,3(1x,I6))') &
            ! 'GP_GA_opt:2 521 myid, isource, numsent ', &
            !                  myid, isource, numsent


            buffer_recv = 0.0d0
            call MPI_RECV( buffer_recv, n_parameters+2, &
                           MPI_DOUBLE_PRECISION, &
                           MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                           MPI_COMM_WORLD, MPI_STAT,  ierr )

            sender       = MPI_STAT( MPI_SOURCE ) 
            i_individual = MPI_STAT( MPI_TAG ) 


            ! received a message from processor "sender" which processed 
            ! individual "i_individual"


            !write(GA_print_unit,'(A,5(1x,I6))') &
            ! 'GP_GA_opt:2 529 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )', &
            !                  myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )
            !write(GA_print_unit,'(A,5(1x,I6))') &
            ! 'GP_GA_opt:2 529 myid, isource, numsent, sender, i_individual ', &
            !                  myid, isource, numsent, sender, i_individual
            !write(GA_print_unit,'(A,2(1x,I6),4x,L1)') &
            !  'GP_GA_opt:2 myid, i_individual, Run_GA_lmdif(i_individual)', &
            !               myid, i_individual, Run_GA_lmdif(i_individual)


            ! store the information received in the above message

            if( Run_GA_lmdif(i_individual) ) then

                child_parameters(1:n_parameters,i_individual) =  &
                                     buffer_recv(1:n_parameters)

                individual_SSE(i_individual)     =       buffer_recv( n_parameters+1) 
                individual_quality(i_individual) = nint( buffer_recv( n_parameters+2) )

                !write(GA_print_unit,'(A,3(1x,I6))') &
                ! 'GP_GA_opt:2 554 myid, n_parameters, i_individual ', &
                !                  myid, n_parameters, i_individual
                !write(GA_print_unit,'(A/(5(1x,E15.7)))') &
                ! 'GP_GA_opt:2 child_parameters(1:n_parameters,i_individual) ', &
                !              child_parameters(1:n_parameters,i_individual) 
                !write(GA_print_unit,'(A,2(1x,I6),1x,E15.7)') &
                ! 'GP_GA_opt:2 myid, i_individual, individual_SSE(i_individual)', &
                !              myid, i_individual, individual_SSE(i_individual)
                !write(GA_print_unit,'(A,3(1x,I6))') &
                ! 'GP_GA_opt:2 myid, i_individual, individual_quality(i_individual) ', &
                !              myid, i_individual,  individual_quality(i_individual)
               
            endif ! Run_GA_lmdif(i_individual) 

            !-------------------------------------------------------------------------------------

            ! check to see if all individuals have been processed

            !write(GA_print_unit,'(A,1x,I6, 4x,L1)') &
            ! 'GP_GA_opt:2 542 myid, numsent < n_GA_individuals ', &
            !                  myid, numsent < n_GA_individuals

            if( numsent <  n_GA_individuals )then 

                ! numsent <  n_GA_individuals    means not all individuals have been processed

                ! send a message to the processor "sender" which just sent a message saying it
                ! completed an individual, and tell it to process
                ! the individual "i_GA_individual" as the  "numsent+1"  task

                i_GA_individual = i_GA_individual + 1

                call MPI_SEND( i_GA_individual, 1, MPI_INTEGER,    &
                               sender, numsent+1,  MPI_COMM_WORLD, ierr )                          

                ! just sent a new task, so increment the number sent

                numsent = numsent + 1

                !write(GA_print_unit,'(A,4(1x,I6))') &
                !      'GP_GA_opt:2 556  myid, sender, numsent, i_GA_individual ', &
                !                        myid, sender, numsent, i_GA_individual


            else

                ! DONE!

                ! number of tasks sent out is >= number of individuals, so 
                ! all the work has been completed

                ! tell the "sender" processor that it is done and  
                ! send it a message to stop
                                                                                                
                !write(GA_print_unit,'(A,3(1x,I6))') &
                !      'GP_GA_opt:2 send msg to stop  myid, numsent, i_GA_individual ', &
                !                                     myid, numsent, i_GA_individual

                call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &  
                               sender, 0,  MPI_COMM_WORLD, ierr )              


            endif ! numsent

        enddo ! isource


        !----------------------------------------------------------------------

        ! this section takes care of the case where there are fewer GA individuals
        ! than (number of procs) - 1

        ! without the code below,  the program hangs because the processors
        ! with numbers  (n_GA_individuals+1)  to (numprocs-1)  
        ! are waiting for a signal to stop
        ! and that is never going to be sent from the loop above.

        ! so when the above loop is finished, send a stop signal to the unused
        ! processors so the program can continue
        
        if( n_GA_individuals < numprocs -1 )then

            do  i = n_GA_individuals+1, numprocs-1
    
                call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &  
                               i , 0,  MPI_COMM_WORLD, ierr )   
   
            enddo ! i 

        endif ! n_GA_individuals < numprocs -1

        !----------------------------------------------------------------------



    else  ! not myid == 0


        ! code for processors 1 - ( numprocs - 1 ) 


        ! these processors wait until a message is received from 
        ! processor 0 telling it to process
        ! the individual named in the message tag = MPI_STAT( MPI_TAG )

        !write(GA_print_unit,'(A,1x,I6)') &
        !          'GP_GA_opt:3  myid, before recv_loop    ', myid

        recv_loop:&
        do 

            
            call MPI_RECV( i_dummy, 1, MPI_INTEGER,    &                                      
                           0, MPI_ANY_TAG,  MPI_COMM_WORLD, MPI_STAT, ierr )      

            !---------------------------------------------------------------
            !write(GA_print_unit,'(A,3(1x,I6))') &
            !      'GP_GA_opt:3  myid, MPI_STAT( MPI_TAG ) ', &
            !                    myid, MPI_STAT( MPI_TAG ) 

            ! was a stop signal received ?

            ! if the tag is <= 0, this is a stop signal

            if(  MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop

            !---------------------------------------------------------------

            ! process the individual named in the message tag

            i_2_individual = MPI_STAT( MPI_TAG ) 

            !write(GA_print_unit,'(A,4(1x,I6))') &
            !      'GP_GA_opt:3  myid, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual ', &
            !                    myid, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual

            buffer = 0.0D0

            !write(GA_print_unit,'(A,2(1x,I6),4x,L1)') &
            !  'GP_GA_opt:3 myid, i_2_individual, Run_GA_lmdif(i_2_individual)', &
            !               myid, i_2_individual, Run_GA_lmdif(i_2_individual)

            if( Run_GA_lmdif(i_2_individual)) then

                !write(GA_print_unit,'(A,3(1x,I6))') &
                !  'GP_GA_opt:3 call setup_run_fcn  myid, i_2_individual   ', &
                !                                   myid, i_2_individual

                ! do the Runge-Kutta integration for individual i_2_individual

                call setup_run_fcn( i_2_individual, child_parameters, individual_quality )

                !write(GA_print_unit,'(A,3(1x,I6))') &
                !  'GP_GA_opt:3 AFTER call setup_run_fcn  myid, i_2_individual   ', &
                !                                         myid, i_2_individual

                !-------------------------------------------------------------------------

                buffer(1:n_parameters) = child_parameters(1:n_parameters,i_2_individual)
                buffer(n_parameters+1) = individual_SSE(i_2_individual)
                buffer(n_parameters+2) = real( individual_quality(i_2_individual), kind=8 )

            endif !  Run_GA_lmdif(i_2_individual)

            !write(GA_print_unit,'(A,3(1x,I6))') &
            !     'GP_GA_opt:3 705 myid, n_parameters, i_2_individual ', &
            !                      myid, n_parameters, i_2_individual
            !write(GA_print_unit,'(A/(5(1x,E15.7)))') &
            !     'GP_GA_opt:3 child_parameters(1:n_parameters,i_2_individual) ', &
            !                  child_parameters(1:n_parameters,i_2_individual) 
            !write(GA_print_unit,'(A,3(1x,I6))') &
            !     'GP_GA_opt:3 myid, i_2_individual, individual_quality(i_2_individual) ', &
            !                  myid, i_2_individual,  individual_quality(i_2_individual)
            !write(GA_print_unit,'(A,2(1x,I6), 1x, E15.7)') &
            !      'GP_GA_opt:3 send results myid, i_2_individual, individual_SSE(i_2_individual)   ', &
            !                                myid, i_2_individual, individual_SSE(i_2_individual)



            ! send the R-K integration results for individual i_2_individual to processor 0

            call MPI_SEND( buffer, n_parameters+2, &
                           MPI_DOUBLE_PRECISION, 0, i_2_individual, MPI_COMM_WORLD, ierr )


            !---------------------------------------------------------------

            ! code to ensure that an error does not allow this loop to run forever

            nsafe = nsafe + 1

            !write(GA_print_unit,'(A,2(1x,I10))') &
            !      'GP_GA_opt: myid,   nsafe = ', myid, nsafe 

            if( nsafe > 10 * n_GA_individuals ) then
                write(GA_print_unit,'(A,1x,I10)') &
                  'GP_GA_opt: too many iterations  nsafe = ', nsafe 
                call MPI_FINALIZE(ierr) 
                stop 'bad nsafe'
            endif ! nsafe

            !---------------------------------------------------------------

         enddo  recv_loop


    endif ! myid == 0

    !write(GA_print_unit,'(A,2(1x,I6))') &
    !      'GP_GA_opt: after recv_loop  myid = ', myid


    !-------------------------------------------------------------------

    ! wait until all n_GA_individuals have been processed

    call MPI_BARRIER( MPI_COMM_WORLD, ierr )

    !write(GA_print_unit,'(A,2(1x,I6))') &
    !      'GP_GA_opt: after barrier 2 i_GA_generation, myid = ', &
    !                                  i_GA_generation, myid

    !-------------------------------------------------------------------


    !  calculate the fitness for this generation


    if( myid == 0  )then

        !write(GA_print_unit,'(A)')  'GP_GA_opt: individual_SSE  '
        !write(GA_print_unit,'(5(1x,E15.7))')  individual_SSE(1:n_GA_individuals)  

        write(GA_print_unit,'(/A,1x,I6)') &
              'GP_GA_opt: call calc_fitness i_GA_generation ', &
                                            i_GA_generation

        call calc_fitness( child_parameters, individual_quality, &
                           i_GA_Best_Parent, Parent_Parameters, L_stop_run )
 

        write(GA_print_unit,'(/A,1x,I6/)') &
              'GP_GA_opt: aft call calc_fitness i_GA_generation ', &
                                                i_GA_generation
                           

        !write(GA_print_unit,'(/A)') &
        !  'GP_GA_opt: call select_best_RK_lmdif_result '

        ! this subroutine, if called, will run lmdif on the best individual
        ! for this generation, instead of waiting to run lmdif after the
        ! generation loop

        !call select_best_RK_lmdif_result( &
        !            i_GA_best_parent, parent_parameters, &
        !            child_parameters, &
        !            individual_quality, L_stop_run  )     

        !write(GA_print_unit,'(A/)') &
        !  'GP_GA_opt: aft call select_best_RK_lmdif_result '


    endif ! myid == 0

    !if( i_GA_generation > 1 )then
    !    write(GA_print_unit,'(A,2(1x,I6))') &
    !      'GP_GA_opt: at stop  i_GA_generation, myid = ', &
    !                           i_GA_generation, myid
    !    !call MPI_FINALIZE(ierr) ! debug only
    !    !stop ! debug only
    !endif

    !-------------------------------------------------------------------


    call MPI_BCAST( L_stop_run,  1,    &
                    MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )

    call MPI_BARRIER( MPI_COMM_WORLD, ierr )   ! necessary ?

    if( L_stop_run )then

        !write(GA_print_unit,'(/A,1x,I6/)') &
        !'GP_GA_opt: L_stop_fitness is true &
        ! &so last generation is i_GA_generation ', &
        !                        i_GA_generation

        i_GA_generation_last = i_GA_generation

        exit

    endif ! L_stop_run

enddo  ! i_generation


!write(GA_print_unit,'(//A,1x,I6/)') &
!  'GP_GA_opt:  finished the loop on i_GA_generation  myid =  ', myid

!----------------------------------------------------------------------

! wait until all processors have finished the generation loop

call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

!write(GA_print_unit,'(A,1x,I6)') &
!      'GP_GA_opt: after barrier 3 myid = ', myid

!----------------------------------------------------------------------


! finished all generations,

! now call lmdif on the best individual of the last generation
! and determine if lmdif has improved the fitness of this individual
! then save the parameters of the fitter of the two results, the RK result
! and the lmdif result

                           
if( myid == 0  )then

    write(GA_print_unit,'(/A)') &
      'GP_GA_opt: call select_best_RK_lmdif_result '

    call select_best_RK_lmdif_result( &
                i_GA_best_parent, parent_parameters, &
                child_parameters, &
                individual_quality, L_stop_run  )     

    write(GA_print_unit,'(A/)') &
      'GP_GA_opt: aft call select_best_RK_lmdif_result '


endif ! myid == 0



!------------------------------------------------------------------------

! broadcast individual_fitness


!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: broadcast individual_fitness myid = ', myid

message_len = 1
call MPI_BCAST( individual_fitness, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: aft broadcast individual_fitness  ierr = ', ierr


!------------------------------------------------------------------------

! broadcast Individual_SSE_best_parent


!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: broadcast Individual_SSE_best_parent myid = ', myid

message_len = 1
call MPI_BCAST( Individual_SSE_best_parent, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: aft broadcast Individual_SSE_best_parent  ierr = ', ierr


!------------------------------------------------------------------------

! broadcast GP_Individual_Node_Parameters


!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: broadcast GP_Individual_Node_Parameters  myid = ', myid

message_len = n_trees * n_nodes

call MPI_BCAST( GP_Individual_Node_Parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: aft broadcast GP_Individual_Node_Parameters  ierr = ', ierr


!------------------------------------------------------------------------

! broadcast GP_Individual_Initial_Conditions


!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: broadcast GP_Individual_Initial_Conditions myid = ', myid

message_len = n_CODE_equations

call MPI_BCAST( GP_Individual_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!write(GA_print_unit,'(/A,1x,I6)') &
! 'GP_GA_opt: aft broadcast GP_Individual_Initial_Conditions ierr = ', ierr


!------------------------------------------------------------------------


if( myid == 0  )then
    write(GA_print_unit,'(//A/)') 'GP_GA_opt:  final parent parameters  '
    write(GA_print_unit,'(A)') 'i_GA_individual                  parent_parameters '
    do  i_GA_individual = 1, n_GA_individuals
        write(GA_print_unit,'(I6,12(1x,E15.7 ))') &
          i_GA_individual, parent_parameters(1:n_parameters,i_GA_individual)
    enddo !  i_GA_individual
endif ! myid == 0


return


end subroutine GPCODE_GA_lmdif_Parameter_Optimization
