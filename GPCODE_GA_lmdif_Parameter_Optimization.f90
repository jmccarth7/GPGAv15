subroutine GPCODE_GA_lmdif_Parameter_Optimization( &
                  i_GP_Generation,i_GP_individual, &
                  new_group, new_comm, n_GP_inds  )

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

integer(kind=4),intent(in) :: new_group
integer(kind=4),intent(in) :: new_comm 
integer(kind=4),intent(in) :: n_GP_inds            

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


!  divider is the number of cpus in the current partition

real(kind=8),&
 dimension(n_GP_parameters,divider) ::  parent_parameters
real(kind=8),&
 dimension(n_GP_parameters,divider) ::  child_parameters


real(kind=8), dimension(n_GP_parameters + 2)  :: buffer
real(kind=8), dimension(n_GP_parameters + 2)  :: buffer_recv


integer(kind=4) ::      i
integer(kind=4) :: i_GA_Best_Parent

integer(kind=4) :: i_GA_generation_last


real(kind=8),parameter :: zero = 0.0d0


! individual_quality contains information on the result of lmdif
! if lmdif encounters an error, set individual_quality to -1
! if < 0 , reject this individual  ! jjm

integer(kind=4) :: individual_quality(divider)

real(kind=8), external :: indiv_fitness

logical :: L_stop_run

logical :: L_too_many_iters


integer(kind=4) :: jj
integer(kind=4) :: i_ga_ind
integer(kind=4) :: i_code_equation


integer(kind=4) :: new_rank 

integer(kind=4) :: n_procs   
!----------------------------------------------------------------------

!write(6,'(A,1x,I5,4x,L1)') 'GP_GA_opt: myid, L_GA_print   ', myid, L_GA_print
!write(6,'(A,1x,I5,1x,I5)') 'GP_GA_opt: myid, GA_print_unit', myid, GA_print_unit



call mpi_comm_rank( new_comm, new_rank, ierr ) 
call MPI_COMM_SIZE( new_comm, n_procs, ierr)  


write(6,'(A,1x,I5,1x,I5)') 'GP_GA_opt:1  myid, new_rank ', myid, new_rank


L_too_many_iters = .FALSE.

i_dummy = 0


do  jj = 1, n_GP_parameters+2
    buffer(jj)      = 0.0D0
    buffer_recv(jj) = 0.0D0
enddo ! jj


n_parameters = n_GP_parameters



!!if( myid == 0 )then
!!    !if( L_ga_print )then
!!    !    write(GA_print_unit,'(//A)') 'GP_GA_opt: at entry  '
!!    !    write(GA_print_unit,'(A,1x,E24.16)') 'GP_GA_opt: dt ', dt
!!    !    write(GA_print_unit,'(/A,1x,I10)') &
!!    !          'GP_GA_opt: n_procs        =   ', n_procs
!!    !    write(GA_print_unit,'(A,1x,I10)') &
!!    !          'GP_GA_opt: n_parameters    =   ', n_parameters
!!    !    write(GA_print_unit,'(A,1x,I10)') &
!!    !          'GP_GA_opt: n_GP_parameters =   ', n_GP_parameters
!!        write(6,'(//A)') 'GP_GA_opt: at entry  '
!!        write(6,'(A,1x,E24.16)') 'GP_GA_opt: dt ', dt
!!        write(6,'(A,1x,I10)') &
!!              'GP_GA_opt: n_parameters    =   ', n_parameters
!!        write(6,'(A,1x,I10)') &
!!              'GP_GA_opt: n_GP_parameters =   ', n_GP_parameters
!!        flush(6)
!!    !endif ! L_ga_print
!!endif ! myid == 0
!
!
!! jjm 20130417 >>>>>>>>>>>>>>>
!!if( myid == 0) then
!!
!!    if( L_ga_print )then
!!
!!        write(GA_print_unit,'(/A)') &
!!                      'GP_GA_opt: i_tree, i_node, GP_Indiv_Node_Param'
!!        do  i_tree=1,n_trees
!!            do  i_node=1,n_nodes
!!                !if( abs( GP_Individual_Node_Parameters(i_node,i_tree) ) > 1.0e-20 )then
!!                    !write(GA_print_unit,'(A,2(1x,I6),1x,E15.7)') &
!!                    !  'GP_GA_opt: i_tree, i_node, GP_Indiv_Node_Param', &
!!                    write(GA_print_unit,'(8x,2(1x,I6),1x,E15.7)') &
!!                                  i_tree, i_node, GP_Individual_Node_Parameters(i_node,i_tree)
!!                !endif ! abs( GP_Indiv_Node_Param(i_node,i_tree) ) > 1.0e-20
!!            enddo ! i_node
!!        enddo  ! i_tree
!!
!!        write(GA_print_unit,'(/A)') &
!!                     'GP_GA_opt: i_tree, i_node, GP_Indiv_Node_Type'
!!        do  i_tree=1,n_trees
!!            do  i_node=1,n_nodes
!!                if( GP_Individual_Node_Type(i_node,i_tree) > -9999 )then
!!                    write(GA_print_unit,'(8x,3(1x,I6))') &
!!                         i_tree, i_node, GP_Individual_Node_Type(i_node,i_tree)
!!                endif ! GP_Indiv_Node_Type(i_node,i_tree) > -9999
!!            enddo ! i_node
!!        enddo  ! i_tree
!!
!!        write(GA_print_unit,'(A)')' '
!!    endif ! L_ga_print
!!
!!endif ! myid == 0
!! jjm 20130417 <<<<<<<<<<<<<<<
!
!

!-----------------------------------------------------------------------------

if( n_parameters .le. 0) then

    if( L_GA_print )then
        write(GA_print_unit,'(A)')        &
              'GP_GA_opt: ERROR: n_parameters </= 0'
        write(GA_print_unit,'(A,1x,I10)') &
              'GP_GA_opt: n_parameters =   ', n_parameters
    endif ! L_GA_print 

    write(6,'(A)')        &
          'GP_GA_opt: ERROR: n_parameters </= 0'
    write(6,'(A,1x,I10)') &
          'GP_GA_opt: n_parameters =   ', n_parameters

    return ! 'n_par<=0'

endif !   n_parameters .le. 0


!if( n_time_steps .lt. n_parameters) then
!    write(GA_print_unit,'(A)') &
!          'GP_GA_opt: ERROR: n_time_steps < n_parameters'
!    write(GA_print_unit,'(A,1x,I10)') &
!          'GP_GA_opt:  n_time_steps = ', n_time_steps
!    write(GA_print_unit,'(A,1x,I10)') &
!          'GP_GA_opt: n_parameters =   ', n_parameters
!    write(6,'(A)') &
!          'GP_GA_opt: ERROR: n_time_steps < n_parameters'
!    write(6,'(A,1x,I10)') &
!          'GP_GA_opt:  n_time_steps = ', n_time_steps
!    write(6,'(A,1x,I10)') &
!          'GP_GA_opt: n_parameters =   ', n_parameters
!    stop 'n_time < n_par'
!endif ! n_time_steps .lt. n_parameters



child_parameters( 1:n_GP_parameters, 1:divider) = 0.0d0


!-----------------------------------------------------------------------------

! set up MPI process


!if( new_rank == 0 )then
!
!    if( L_ga_print )then
!        write(GA_print_unit,'(/A,1x,I10, 2(1x,I6)/)')&
!          'GP_GA_opt: divider, n_parameters, n_procs ', &
!                      divider, n_parameters, n_procs
!    endif ! L_ga_print
!
!endif ! new_rank == 0

!-----------------------------------------------------------------------------


L_stop_run  = .FALSE.
!L_stop_run  = .TRUE.

do  i_GA_generation = 1, n_GA_Generations

    ! Run_GA_lmdif determines if the new child
    ! has to be sent to lmdif for 'local' optimization

    Run_GA_lmdif=.false.

    !if( new_rank == 0 )then
    !    if( L_ga_print )then
    !        write(GA_print_unit,'(/A,1x,I6,1x,A/)') &
    !              'GA Generation ',i_GA_generation,' is underway'
    !        flush(GA_print_unit)
    !    endif ! L_ga_print
    !    write(6,'(/A,1x,I6,1x,A/)') &
    !              'GA Generation ',i_GA_generation,' is underway'
    !    flush(6)
    !endif ! new_rank == 0



    if( new_rank == 0 )then

        if( i_GA_generation .eq. 1 ) then

            ! on the first generation,
            ! randomly create the initial individual parameter arrays
            ! for each individual


            ! sets:
            !  child_parameters

            call Initialize_GA_Child_Parameters( Child_Parameters )


            Run_GA_lmdif=.true.

            ! print child parameters at start of the generation

            !if( L_ga_print )then
            !    !write(GA_print_unit,'(/A,1x,I6)') &
                write(6,'(/A,1x,I6)') &
                'GP_GA_opt:1 child parameters at start of generation: ', &
                                                    i_GA_generation
            !    do  i_ga_ind = 1, divider
            !        !write(GA_print_unit,'(I6,4(1x,E24.16)/(4(1x,E24.16)))') &
            !        write(6,'(I6,4(1x,E24.16)/(4(1x,E24.16)))') &
            !              i_ga_ind, &
            !              ( child_parameters(jj,i_ga_ind), jj = 1,n_parameters )
            !    enddo ! i_ga_ind
            !endif ! L_ga_print

            !flush(6)


        else  ! i_GA_generation > 1


            ! create the second 'generation' of parameter estimates using either:

            !    i) save elites from last generation from being changed
            !   ii) 'Fitness-Proportionate Reproduction;
            !  iii) GA Crossover;
            !   iv) GA Mutation

            !-------------------------------------------------------------------------------

            !   save the most fit individuals for the next generation

            !if( L_ga_print )then
            !    write(GA_print_unit,'(/A)')'GP_GA_opt: call GA_save_elites '
            !endif ! L_ga_print

            ! uses:
            ! individual_ranked_fitness

            ! sets:
            ! ga_individual_elites

            call GA_save_elites( )

            !-------------------------------------------------------------------------------

            !!  replace the parameters of any individual with quality < 0 with new
            !!  random numbers

            !if( L_ga_print )then
            !    write(GA_print_unit,'(/A)')'GP_GA_opt: call GA_replace_bad_individuals  '
            !endif ! L_ga_print
            !    write(6,'(/A)')'GP_GA_opt: call GA_replace_bad_individuals  '

            !!call GA_replace_bad_individuals(Child_Parameters, individual_quality )

            !-------------------------------------------------------------------------------

            !   do initial "GA Fitness-Proportionate Reproduction"
            !   to create a new population of children for all n_GA_individual

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A)')&
            !          'GP_GA_opt: call GA_Fitness_Proportionate_Reproduction '
            !endif ! L_ga_print


            ! uses:
            !  individual_quality
            !  Individual_Ranked_Fitness
            !  Parent_Parameters

            ! sets:
            !  Run_GA_lmdif
            !  individual_quality
            !  Individual_Ranked_Fitness
            !  Child_Parameters


            call GA_Fitness_Proportionate_Reproduction(&
                            Parent_Parameters,Child_Parameters, &
                                                individual_quality )

            !-------------------------------------------------------------------------------

            !   do "GA Parameter Crossover" Operations Using Tournament-Style Selection
            !   and randomly use it to replace the parents with the new children

            if( n_GA_Crossovers .gt. 0) then

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,1x,I6)')&
                !      'GP_GA_opt: call GA_Tournament_Style_Sexual_Repro  &
                !      &n_GA_Crossovers',  n_GA_Crossovers
                !endif ! L_ga_print

                ! uses:
                !  Individual_Ranked_Fitness

                ! sets:
                !  Child_Parameters
                !  Run_GA_lmdif
                !  individual_quality

                call GA_Tournament_Style_Sexual_Reproduction( &
                            Parent_Parameters, Child_Parameters, individual_quality )

            endif !   n_GA_Crossovers .gt. 0


            !-------------------------------------------------------------------------------

            !   do "GA Parameter Mutation" Operations

            !   select a random individual and put a new random number into one of
            !   its parameters


            if( n_GA_Mutations .gt. 0) then

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,1x,I6)')&
                !      'GP_GA_opt: call GA_Mutations  n_GA_Mutations',  n_GA_Mutations
                !endif ! L_ga_print

                ! uses:

                ! sets:
                !  child_parameters
                !  Run_GA_lmdif
                !  individual_quality

                call GA_Mutations( Child_Parameters, individual_quality )

            endif !   n_GA_Mutations .gt. 0


            !-------------------------------------------------------------------------------

            !   do "GA Parameter rand_replace" Operations

            !   select a random, non-elite individual and put new random numbers into
            !   its parameters

            if( n_GA_rand_replaces > 0) then

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,1x,I6)')&
                !      'GP_GA_opt: call GA_rand_replace  n_GA_rand_replaces',  &
                !                                        n_GA_rand_replaces
                !endif ! L_ga_print


                ! uses:

                ! sets:
                !  child_parameters
                !  Run_GA_lmdif
                !  individual_quality


                call GA_random_replace( Child_Parameters, individual_quality )


            endif !   n_GA_rand_replaces .gt. 0

            !-------------------------------------------------------------------------------


            ! print child parameters at start of the generation

            !if( i_GA_generation == n_GA_generations               .or.  &
            !    mod( i_GA_generation, GA_child_print_interval ) == 0      )then

            !if( i_GA_generation == n_GA_generations  )then

                !if( L_ga_print )then
                !    write(GA_print_unit,'(/A,1x,I6)') &
                !    'GP_GA_opt:2 child parameters at start of GA generation:', &
                !                                            i_GA_generation
                !    do  i_ga_ind = 1, divider
                !        write(GA_print_unit,'(I6,1x,10(1x,E15.7)/(10(1x,E15.7)))') &
                !              i_ga_ind, &
                !             (child_parameters(jj, i_ga_ind),&
                !                               jj = 1,n_parameters )
                !    enddo ! i_ga_ind
                !endif ! L_ga_print

            !endif ! i_GA_generation == n_GA_generations ...


        endif ! i_GA_generation .eq. 1




    endif ! new_rank == 0  ! new_rank == 0


    !------------------------------------------------------------------------

    !  broadcast child parameters


    !if( L_ga_print )then
    !    write(GA_print_unit,'(/A,2(1x,I6))') &
    !    'GP_GA_opt:  broadcast child parameters new_rank, i_GA_generation ', &
    !                                            new_rank, i_GA_generation
        write(6,'(/A,2(1x,I6))') &
        'GP_GA_opt:  broadcast child parameters new_rank, i_GA_generation ', &
                                                new_rank, i_GA_generation
    !endif ! L_ga_print

    child_number =  divider * n_GP_parameters

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,4(1x,I6)/)') &
    !    'GP_GA_opt: new_rank, divider, n_GP_parameters, child_number =', &
    !                new_rank, divider, n_GP_parameters, child_number
        write(6,'(A,4(1x,I6)/)') &
        'GP_GA_opt: new_rank, divider, n_GP_parameters, child_number =', &
                    new_rank, divider, n_GP_parameters, child_number
    !endif ! L_ga_print


    call MPI_BCAST( Child_Parameters,  child_number,    &
                    MPI_DOUBLE_PRECISION, 0, new_comm, ierr )


    !if( L_ga_print )then
    !    write(GA_print_unit,'(/A,2(1x,I10)/)') &
    !     'GP_GA_opt: child  broadcast new_rank, ierr = ', new_rank, ierr
    !    write(GA_print_unit,'(/A,2(1x,I10)/)') &
    !     'GP_GA_opt: new_rank, divider = ', new_rank, divider
    !    do  i_ga_ind = 1, divider
    !        write(GA_print_unit,'(I3,1x,I3,1x,12(1x,E15.7))') &
    !              new_rank, i_ga_ind, &
    !              ( child_parameters(jj,i_ga_ind), jj = 1,n_parameters )
    !    enddo ! i_ga_ind
    !endif ! L_ga_print




    !------------------------------------------------------------------------

    ! broadcast Run_GA_lmdif

    !if( L_ga_print )then
    !    write(GA_print_unit,'(/A,1x,I6)') &
    !      'GP_GA_opt:  broadcast Run_GA_lmdif i_GA_generation ', i_GA_generation
    !endif ! L_ga_print



    call MPI_BCAST( Run_GA_lmdif,  divider,    &
                        MPI_LOGICAL, 0, new_comm, ierr )

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,1x,I10/)') &
    !     'GP_GA_opt: Run_GA_lmdif  broadcast ierr = ', ierr
    !    write(GA_print_unit,'(/A,1x,I6,/(10(2x,L1)))') &
    !     'GP_GA_opt: new_rank, Run_GA_lmdif  ', new_rank, Run_GA_lmdif
        write(6,'(/A,1x,I6,/(10(2x,L1)))') &
         'GP_GA_opt: new_rank, Run_GA_lmdif  ', new_rank, Run_GA_lmdif
    !endif ! L_ga_print


    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !   begin RK fcn integration segment
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    !if( new_rank == 0 )then
    !    if( L_ga_print )then
    !        write(GA_print_unit,'(/A,1x,I4,1x,A,1x,I6/)') &
    !          'GP_GA_opt: begin RK fcn integration   i_GA_generation =', &
    !                                                 i_GA_generation,    &
    !                                         '       divider =', &
    !                                                 divider
            write(6,'(/A,1x,I4,1x,A,1x,I6/)') &
              'GP_GA_opt: begin RK fcn integration   i_GA_generation =', &
                                                     i_GA_generation,    &
                                             '       divider =', &
                                                     divider
    !    endif ! L_ga_print
    !endif !  new_rank == 0


    individual_quality( 1: divider ) = 1

    nsafe = 0

    if( new_rank == 0  )then

        ! processor 0 sends job assignments to N  processors
        ! where N is the smaller of (number of processors -1)
        !                    and the number of individuals

        ! numsent is the number of messages sent up to now

        numsent = 0
        i_ga_ind = 0

        do  isource = 1, min( n_procs-1, divider )


            i_ga_ind = i_ga_ind + 1

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,1x,I6, 4x, L1)') &
            !     'GP_GA_opt:1 494 i_ga_ind, Run_GA_lmdif(i_ga_ind)  ', &
            !                      i_ga_ind, Run_GA_lmdif(i_ga_ind)
            !endif ! L_ga_print


            call MPI_SEND( i_dummy,  1, MPI_INTEGER,    &
                           isource, isource,  new_comm, ierr )

            numsent = numsent + 1

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,4(1x,I6))') &
            !     'GP_GA_opt:1 504 new_rank, isource, i_ga_ind, numsent ', &
            !                      new_rank, isource, i_ga_ind, numsent
                write(6,'(A,4(1x,I6))') &
                 'GP_GA_opt:1 504 new_rank, isource, i_ga_ind, numsent ', &
                                  new_rank, isource, i_ga_ind, numsent
            !endif ! L_ga_print

        enddo ! isource


        ! at this point i_ga_ind = numsent


        !if( L_ga_print )then
        !    write(GA_print_unit,'(A,4(1x,I6))') &
        !         'GP_GA_opt: aft source loop 1 new_rank, i_ga_ind, numsent ', &
        !                                       new_rank, i_ga_ind, numsent
            write(6,'(A,4(1x,I6))') &
                 'GP_GA_opt: aft source loop 1 new_rank, i_ga_ind, numsent ', &
                                               new_rank, i_ga_ind, numsent
        !endif ! L_ga_print

        !-------------------------------------------------------------------------------------

        ! processor 0 loops over the number of individuals and waits for a message
        ! from the other processors

        do  isource = 1, divider

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,3(1x,I6))') &
            !     'GP_GA_opt:2 521 new_rank, isource, numsent ', &
            !                      new_rank, isource, numsent
            !endif ! L_ga_print


            buffer_recv = 0.0d0
            call MPI_RECV( buffer_recv, n_GP_parameters+2, &
                           MPI_DOUBLE_PRECISION, &
                           MPI_ANY_SOURCE, MPI_ANY_TAG,  &
                           new_comm, MPI_STAT,  ierr )

            sender       = MPI_STAT( MPI_SOURCE )
            i_individual = MPI_STAT( MPI_TAG )


            ! received a message from processor "sender" which processed
            ! individual "i_individual"


            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,5(1x,I6))') &
            !     'GP_GA_opt:2 529 new_rank, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )', &
            !                      new_rank, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )
            !    write(GA_print_unit,'(A,5(1x,I6))') &
            !     'GP_GA_opt:2 529 new_rank, isource, numsent, sender, i_individual ', &
            !                      new_rank, isource, numsent, sender, i_individual
            !    write(GA_print_unit,'(A,2(1x,I6),4x,L1)') &
            !      'GP_GA_opt:2 new_rank, i_individual, Run_GA_lmdif(i_individual)', &
            !                   new_rank, i_individual, Run_GA_lmdif(i_individual)
            !endif ! L_ga_print


            ! store the information received in the above message

            if( Run_GA_lmdif(i_individual) ) then


                do  jj = 1, n_GP_parameters
                    child_parameters(jj,i_individual) =  buffer_recv(jj)
                enddo ! jj

                individual_SSE(i_individual)     =       buffer_recv( n_GP_parameters+1)
                individual_quality(i_individual) = nint( buffer_recv( n_GP_parameters+2) )

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,3(1x,I6))') &
                !     'GP_GA_opt:2 554 new_rank, n_parameters, i_individual', &
                !                      new_rank, n_parameters, i_individual
                !    !write(GA_print_unit,'(A,1x,I6,1x, A/(5(1x,E15.7)))') &
                !    ! 'GP_GA_opt:2 new_rank = ',new_rank, &
                !    ! ' child_parameters(1:n_parameters,i_individual)', &
                !    !   child_parameters(1:n_parameters,i_individual)
                !    write(GA_print_unit,'(A,2(1x,I6),1x,E15.7)') &
                !     'GP_GA_opt:2 new_rank, i_individual, individual_SSE(i_individual)', &
                !                  new_rank, i_individual, individual_SSE(i_individual)
                !    write(GA_print_unit,'(A,3(1x,I6))') &
                !     'GP_GA_opt:2 new_rank, i_individual, individual_quality(i_individual)', &
                !                  new_rank, i_individual, individual_quality(i_individual)
                    write(6,'(A,2(1x,I6),1x,E15.7)') &
                     'GP_GA_opt:2 new_rank, i_individual, individual_SSE(i_individual)', &
                                  new_rank, i_individual, individual_SSE(i_individual)
                    write(6,'(A,3(1x,I6))') &
                     'GP_GA_opt:2 new_rank, i_individual, individual_quality(i_individual)', &
                                  new_rank, i_individual, individual_quality(i_individual)
                !endif ! L_ga_print

            endif ! Run_GA_lmdif(i_individual)

            !--------------------------------------------------------------------------------


            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,1x,I6, 4x,L1)') &
            !     'GP_GA_opt:2 542 new_rank, numsent < divider ', &
            !                      new_rank, numsent < divider
            !endif ! L_ga_print


            ! check to see if all individuals have been processed

            if( numsent <  divider )then

                ! numsent <  divider
                ! means not all individuals have been processed

                ! send a message to the processor "sender"
                ! which just sent a message saying it has
                ! completed an individual, and tell it to process
                ! the individual "i_ga_ind" as the  "numsent+1"  task

                i_ga_ind = i_ga_ind + 1

                call MPI_SEND( i_ga_ind, 1, MPI_INTEGER,    &
                               sender, numsent+1,  new_comm, ierr )

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,4(1x,I6))') &
                !      'GP_GA_opt:2 554 send new_rank, sender, numsent, i_ga_ind', &
                !                            new_rank, sender, numsent, i_ga_ind
                !endif ! L_ga_print


                ! just sent a new task, so increment the number sent

                numsent = numsent + 1

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,4(1x,I6))') &
                !     'GP_GA_opt:2 556  new_rank, sender, numsent, i_ga_ind ', &
                !                       new_rank, sender, numsent, i_ga_ind
                !endif ! L_ga_print


            else

                ! DONE !

                ! number of tasks sent out is >= number of individuals, so
                ! all the work has been completed

                ! tell the "sender" processor that it is done and
                ! send it a message to stop

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,3(1x,I6))') &
                !      'GP_GA_opt:2 send msg to stop  new_rank, numsent, i_ga_ind ', &
                !                                     new_rank, numsent, i_ga_ind
                !endif ! L_ga_print

                call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                               sender, 0,  new_comm, ierr )


            endif ! numsent

        enddo ! isource


        !----------------------------------------------------------------------

        ! this section takes care of the case where there are fewer GA individuals
        ! than (number of procs) - 1

        ! without the code below,  the program hangs because the processors
        ! with numbers  (divider+1)  to (n_procs-1)
        ! are waiting for a signal to stop
        ! and that is never going to be sent from the loop above.

        ! so when the above loop is finished, send a stop signal to the unused
        ! processors so the program can continue

        if( divider < n_procs -1 )then

            do  i = divider+1, n_procs-1

                call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
                               i , 0, new_comm , ierr )
            enddo ! i

        endif ! divider < n_procs -1

        !----------------------------------------------------------------------



    else  ! not new_rank == 0   ! new_rank == 0


        ! code for processors 1 - ( divider )  ! - 1 )


        ! these processors wait until a message is received from
        ! processor 0 telling it to process
        ! the individual named in the message tag = MPI_STAT( MPI_TAG )

        recv_loop:&
        do

            if( new_rank == 1 )then
                !write(GA_print_unit,'(A,1x,I6)') &
                write(6,'(A,1x,I6)') &
                      'GP_GA_opt: new_rank  proc1 in recv_loop', new_rank
            endif


            call MPI_RECV( i_dummy, 1, MPI_INTEGER,    &
                           0, MPI_ANY_TAG,  new_comm , MPI_STAT, ierr )

            !---------------------------------------------------------------

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,2(1x,I6))') &
            !      'GP_GA_opt:3  new_rank, MPI_STAT( MPI_TAG ) ', &
            !                    new_rank, MPI_STAT( MPI_TAG )
            !endif ! L_ga_print


            ! was a stop signal received ?

            ! if the tag is <= 0, this is a stop signal

            if( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop

            !---------------------------------------------------------------

            ! process the individual named in the message tag

            i_2_individual = MPI_STAT( MPI_TAG )

            buffer = 0.0D0

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,4(1x,I6))') &
            !    'GP_GA_opt:3 new_rank, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual', &
            !                 new_rank, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual
                write(6,'(A,4(1x,I6))') &
                'GP_GA_opt:3 new_rank, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual', &
                             new_rank, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual
            !    write(GA_print_unit,'(A,2(1x,I6),4x,L1)') &
            !    'GP_GA_opt:3 new_rank, i_2_individual, Run_GA_lmdif(i_2_individual)', &
            !                 new_rank, i_2_individual, Run_GA_lmdif(i_2_individual)
            !endif ! L_ga_print

            if( Run_GA_lmdif(i_2_individual)) then

                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,2(1x,I6))') &
                !      'GP_GA_opt:3 call setup_run_fcn  new_rank, i_2_individual', &
                !                                       new_rank, i_2_individual
                !endif ! L_ga_print
                write(6,'(A,2(1x,I6))') &
                      'GP_GA_opt:3 call setup_run_fcn  new_rank, i_2_individual', &
                                                       new_rank, i_2_individual
                !flush(6)



                ! do the Runge-Kutta integration for individual i_2_individual

                ! uses:
                !  child_parameters

                ! sets:
                !  individual_quality
                !  individual_SSE
                !  child_parameters


                call setup_run_fcn( i_2_individual, &
                                    child_parameters,individual_quality, &
                                    new_group, new_comm )



                !if( L_ga_print )then
                !    write(GA_print_unit,'(A,3(1x,I6))') &
                !      'GP_GA_opt:3 AFTER call setup_run_fcn new_rank, i_2_individual', &
                !                                            new_rank, i_2_individual
                !    flush(GA_print_unit)
                !endif ! L_ga_print

                write(6,'(A,3(1x,I6))') &
                  'GP_GA_opt:3 AFTER call setup_run_fcn new_rank, i_2_individual', &
                                                        new_rank, i_2_individual
                !flush(6)

                !-------------------------------------------------------------------------

                do  jj = 1, n_GP_parameters
                    buffer(jj) =  child_parameters(jj, i_2_individual)
                enddo ! jj

                buffer(n_GP_parameters+1) = &
                      individual_SSE(i_2_individual)
                buffer(n_GP_parameters+2) = &
                      real( individual_quality(i_2_individual), kind=8 )

            endif !  Run_GA_lmdif(i_2_individual)

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,3(1x,I6))') &
            !    'GP_GA_opt:3 705 new_rank, n_GP_parameters, i_2_individual', &
            !                     new_rank, n_GP_parameters, i_2_individual
            !    write(GA_print_unit,'(A/(5(1x,E24.16)))') &
            !    'GP_GA_opt:3 child_parameters(1:n_parameters,i_2_individual)', &
            !                 child_parameters(1:n_parameters,i_2_individual)
            !    write(GA_print_unit,'(A,3(1x,I6))') &
            !    'GP_GA_opt:3 new_rank, i_2_individual, individual_quality(i_2_individual)', &
            !                 new_rank, i_2_individual, individual_quality(i_2_individual)
            !    write(GA_print_unit,'(A,2(1x,I6), 1x, E15.7)') &
            !    'GP_GA_opt:3 send results &
            !      &new_rank, i_2_individual, individual_SSE(i_2_individual)', &
            !       new_rank, i_2_individual, individual_SSE(i_2_individual)
            !endif ! L_ga_print



            ! send the R-K integration results
            ! for individual i_2_individual to processor 0

            call MPI_SEND( buffer, n_GP_parameters+2,  &
                           MPI_DOUBLE_PRECISION, 0, &
                           i_2_individual, new_comm, ierr )


            !---------------------------------------------------------------

            ! code to ensure that an error does not allow this loop to run forever

            nsafe = nsafe + 1

            !if( L_ga_print )then
            !    write(GA_print_unit,'(A,2(1x,I10))') &
            !      'GP_GA_opt: new_rank, nsafe =', new_rank, nsafe
            !endif ! L_ga_print

            if( nsafe > 100 * divider ) then

                if( L_GA_print )then
                    write(GA_print_unit,'(A,1x,I10)') &
                      'GP_GA_opt: too many iterations  nsafe =', nsafe
                    flush(GA_print_unit) 
                endif ! L_GA_print 

                write(6,'(A,1x,I10)') &
                  'GP_GA_opt: too many iterations  nsafe =', nsafe
                flush(6) 

                L_too_many_iters = .TRUE.  
                exit recv_loop

            endif ! nsafe

            !---------------------------------------------------------------


         enddo  recv_loop


    endif ! new_rank == 0

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,2(1x,I6))') &
    !      'GP_GA_opt: after recv_loop  new_rank = ', new_rank
        write(6,'(A,2(1x,I6))') &
          'GP_GA_opt: after recv_loop  new_rank = ', new_rank
    !endif ! L_ga_print

    if( L_too_many_iters )then
        call MPI_FINALIZE(ierr)
        stop 'bad nsafe'
    endif !  L_too_many_iters 

    !-------------------------------------------------------------------

    ! wait until all divider have been processed

    call MPI_BARRIER( new_comm, ierr )

    !if( L_ga_print )then
    !    write(GA_print_unit,'(A,2(1x,I6))') &
    !      'GP_GA_opt: after barrier 2 i_GA_generation, new_rank = ', &
    !                                  i_GA_generation, new_rank
        write(6,'(A,2(1x,I6))') &
          'GP_GA_opt: after barrier 2 i_GA_generation, new_rank = ', &
                                      i_GA_generation, new_rank
    !endif ! L_ga_print

    !-------------------------------------------------------------------


    !  calculate the fitness for this generation


    if( new_rank == 0  )then

        !if( L_ga_print )then
        !    write(GA_print_unit,'(A)')  'GP_GA_opt: individual_SSE  '
        !    write(GA_print_unit,'(5(1x,E12.5))')  individual_SSE(1:divider)
        !    write(GA_print_unit,'(/A,1x,I6)') &
        !          'GP_GA_opt: call calc_fitness i_GA_generation ', &
        !                                        i_GA_generation
        !    write(GA_print_unit,'(/A,6(1x,E15.7))') &
        !          'GP_GA_opt: GP_Individual_Initial_Conditions ', &
        !                      GP_Individual_Initial_Conditions(1:n_code_equations)
        !endif ! L_ga_print
        write(6,'(/A,1x,I6)') &
              'GP_GA_opt: call calc_fitness i_GA_generation ', &
                                            i_GA_generation
        !flush(6)

        ! uses:
        !  child_parameters
        !  individual_quality
        !  individual_SSE

        ! sets:
        !  individual_quality
        !  individual_ranked_fitness
        !  integrated_SSE
        !  integrated_ranked_fitness


        call calc_fitness( child_parameters, individual_quality, &
                           i_GA_Best_Parent, Parent_Parameters, L_stop_run, &
                           i_GP_Generation, i_GP_individual, &
                           new_group, new_comm  )


        !if( L_ga_print )then
        !    write(GA_print_unit,'(/A,1x,I6/)') &
        !          'GP_GA_opt: aft call calc_fitness i_GA_generation ', &
        !                                            i_GA_generation
        !endif ! L_ga_print
        write(6,'(/A,1x,I6/)') &
               'GP_GA_opt: aft call calc_fitness i_GA_generation ', &
                                                 i_GA_generation
        !flush(6)

        !---------------------------------------------------------------------

        ! this subroutine, if called, will run lmdif on the best individual
        ! for this generation, instead of waiting to run lmdif after the
        ! generation loop

        !    call select_best_RK_lmdif_result( &
        !                i_GP_Generation,i_GP_individual, &
        !                i_GA_best_parent, parent_parameters, &
        !                child_parameters, &
        !                L_stop_run  )
        !                !individual_quality, L_stop_run  )

        !---------------------------------------------------------------------



    endif ! new_rank == 0

    !if( i_GA_generation > 1 )then
    !    if( L_ga_print )then
    !        write(GA_print_unit,'(A,2(1x,I6))') &
    !          'GP_GA_opt: at stop  i_GA_generation, new_rank = ', &
    !                               i_GA_generation, new_rank
    !        call MPI_FINALIZE(ierr) ! debug only
    !        stop ! debug only
    !    endif ! L_ga_print
    !endif ! i_GA_generation > 1

    !-------------------------------------------------------------------


    call MPI_BCAST( L_stop_run,  1,    &
                    MPI_LOGICAL, 0, new_comm, ierr )

    if( L_stop_run )then

        !if( L_ga_print )then
        !    write(GA_print_unit,'(/A,1x,I6/)') &
        !    'GP_GA_opt: L_stop_fitness is true &
        !     &so last generation is i_GA_generation ', &
        !                            i_GA_generation
        !endif ! L_ga_print

        i_GA_generation_last = i_GA_generation

        exit

    endif ! L_stop_run


enddo  ! i_generation






!if( L_ga_print )then
!    write(GA_print_unit,'(//A,1x,I6/)') &
!      'GP_GA_opt:  finished the loop on i_GA_generation  new_rank =', &
!                                                         new_rank
    write(6,'(//A,1x,I6/)') &
      'GP_GA_opt:  finished the loop on i_GA_generation  new_rank =', &
                                                         new_rank
!endif ! L_ga_print

!----------------------------------------------------------------------

! wait until all processors have finished the generation loop

call MPI_BARRIER( new_comm, ierr )    ! necessary?

!if( L_ga_print )then
!    write(GA_print_unit,'(A,1x,I6)') &
!          'GP_GA_opt: after barrier 3 new_rank = ', new_rank
    write(6,'(A,1x,I6)') &
          'GP_GA_opt: after barrier 3 new_rank = ', new_rank
!endif ! L_ga_print

!----------------------------------------------------------------------


! finished all generations,

! now call lmdif on the best individual of the last generation
! and determine if lmdif has improved the fitness of this individual
! then save the parameters of the fitter of the two results, the RK result
! and the lmdif result

! GP_Individual_Node_Parameters loaded from the best parent parameters (RK or lmdif)
! after the RK process and lmdif


if( new_rank == 0  )then

    !if( L_ga_print )then
        write(6,'(//A/A)') &
          'GP_GA_opt: finished all generations', &
          'GP_GA_opt: call select_best_RK_lmdif_result'
    !    flush(6)
    !endif ! L_ga_print

    ! uses:

    ! sets:

    call select_best_RK_lmdif_result( &
                i_GP_Generation,i_GP_individual, &
                i_GA_best_parent, parent_parameters, &
                child_parameters, &
                L_stop_run, new_group, new_comm  )

    !if( L_ga_print )then
    !    write(GA_print_unit,'(/A//)') &
    !      'GP_GA_opt: aft call select_best_RK_lmdif_result '
    !endif ! L_ga_print
    write(6,'(/A//)') &
         'GP_GA_opt: aft call select_best_RK_lmdif_result '
    !flush(6)


endif ! new_rank == 0


!------------------------------------------------------------------------

! broadcast individual_fitness


!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: broadcast individual_fitness myid = ', myid
!endif ! L_ga_print
!write(6,'(A,1x,I5,1x,I5)') 'GP_GA_opt:2  myid, new_rank ', myid, new_rank

!individual_fitness = real( i_GP_individual * new_rank, kind=8)  ! debug only
!individual_fitness = real( new_rank, kind=8)  ! debug only
!individual_fitness = real( myid, kind=8)  ! debug only

write(6,'(/A,2(1x,I6), 1x,F10.2/)') &
     'GP_GA_opt: broadcast individual_fitness  myid, new_rank, individual_fitness', &
                                               myid, new_rank, individual_fitness

!message_len = 1
!call MPI_REDUCE( individual_fitness, sum_if, message_len,    &
!                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, new_comm, ierr )



message_len = 1
call MPI_BCAST( individual_fitness, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
                !MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: aft broadcast individual_fitness  ierr = ', ierr
!if( new_rank == 0 )then
!write(6,'(/A,2(1x,I6), 1x,F10.2/)') &
!     'GP_GA_opt: aft REDUCE myid, new_rank, sum_if', &
!                            myid, new_rank, sum_if
!endif 
!endif ! L_ga_print


!------------------------------------------------------------------------

! broadcast Individual_SSE_best_parent


!if( L_ga_print )then
    write(6,'(/A,1x,I6)') &
     'GP_GA_opt: broadcast Individual_SSE_best_parent new_rank = ', new_rank
!endif ! L_ga_print

message_len = 1
call MPI_BCAST( Individual_SSE_best_parent, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
                !MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: aft broadcast Individual_SSE_best_parent  ierr = ', ierr
!endif ! L_ga_print


!------------------------------------------------------------------------

! broadcast GP_Individual_Node_Parameters


!if( L_ga_print )then
    write(6,'(/A,1x,I6)') &
     'GP_GA_opt: broadcast GP_Individual_Node_Parameters  new_rank = ', new_rank
!endif ! L_ga_print

message_len = n_trees * n_nodes

call MPI_BCAST( GP_Individual_Node_Parameters, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
                !MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: aft broadcast GP_Individual_Node_Parameters  ierr = ', ierr
!endif ! L_ga_print


!------------------------------------------------------------------------

! broadcast GP_Individual_Initial_Conditions


!if( L_ga_print )then
    write(6,'(/A,1x,I6)') &
     'GP_GA_opt: broadcast GP_Individual_Initial_Conditions new_rank = ', new_rank
!endif ! L_ga_print

message_len = n_CODE_equations

call MPI_BCAST( GP_Individual_Initial_Conditions, message_len,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
                !MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: aft broadcast GP_Individual_Initial_Conditions ierr = ', ierr
!endif ! L_ga_print


!------------------------------------------------------------------------


!if( myid == 0  )then
!    if( L_ga_print )then
!        write(GA_print_unit,'(//A/)') &
!        'GP_GA_opt:  final initial_conditions '
!        write(GA_print_unit,'(A)') &
!        'i_CODE_equation                  GP_Individual_Initial_Conditions '
!        do  i_code_equation = 1, n_code_equations
!            write(GA_print_unit,'(I6,1x,E15.7 )') &
!              i_code_equation, GP_Individual_Initial_Conditions(i_code_equation)
!        enddo !  i_code_equation
!
!        write(GA_print_unit,'(//A/)') &
!        'GP_GA_opt:  final parent parameters  '
!        write(GA_print_unit,'(A)') &
!        'i_ga_ind                  parent_parameters '
!        do  i_ga_ind = 1, divider
!            write(GA_print_unit,'(I6,12(1x,E15.7 ))') &
!              i_ga_ind, parent_parameters(1:n_parameters,i_ga_ind)
!        enddo !  i_ga_ind
!    endif ! L_ga_print
!endif ! myid == 0


return


end subroutine GPCODE_GA_lmdif_Parameter_Optimization
