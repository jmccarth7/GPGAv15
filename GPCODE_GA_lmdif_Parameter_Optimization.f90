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



!call MPI_GROUP_RANK( new_group, new_rank, ierr )
call mpi_comm_rank( new_comm, new_rank, ierr ) 
call MPI_COMM_SIZE( new_comm, n_procs, ierr)  

write(6,'(A,1x,I5,1x,I5)') 'GP_GA_opt:1  myid, new_rank ', myid, new_rank


!! debug only  !L_too_many_iters = .FALSE.
!! debug only  !
!! debug only  !i_dummy = 0
!! debug only  !
!! debug only  !
!! debug only  !do  jj = 1, n_GP_parameters+2
!! debug only  !    buffer(jj)      = 0.0D0
!! debug only  !    buffer_recv(jj) = 0.0D0
!! debug only  !enddo ! jj
!! debug only  !
!! debug only  !
!! debug only  !n_parameters = n_GP_parameters
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !!if( myid == 0 )then
!! debug only  !!    !if( L_ga_print )then
!! debug only  !!    !    write(GA_print_unit,'(//A)') 'GP_GA_opt: at entry  '
!! debug only  !!    !    write(GA_print_unit,'(A,1x,E24.16)') 'GP_GA_opt: dt ', dt
!! debug only  !!    !    write(GA_print_unit,'(/A,1x,I10)') &
!! debug only  !!    !          'GP_GA_opt: n_procs        =   ', n_procs
!! debug only  !!    !    write(GA_print_unit,'(A,1x,I10)') &
!! debug only  !!    !          'GP_GA_opt: n_parameters    =   ', n_parameters
!! debug only  !!    !    write(GA_print_unit,'(A,1x,I10)') &
!! debug only  !!    !          'GP_GA_opt: n_GP_parameters =   ', n_GP_parameters
!! debug only  !!        write(6,'(//A)') 'GP_GA_opt: at entry  '
!! debug only  !!        write(6,'(A,1x,E24.16)') 'GP_GA_opt: dt ', dt
!! debug only  !!        write(6,'(A,1x,I10)') &
!! debug only  !!              'GP_GA_opt: n_parameters    =   ', n_parameters
!! debug only  !!        write(6,'(A,1x,I10)') &
!! debug only  !!              'GP_GA_opt: n_GP_parameters =   ', n_GP_parameters
!! debug only  !!        flush(6)
!! debug only  !!    !endif ! L_ga_print
!! debug only  !!endif ! myid == 0
!! debug only  !
!! debug only  !
!! debug only  !! jjm 20130417 >>>>>>>>>>>>>>>
!! debug only  !!if( myid == 0) then
!! debug only  !!
!! debug only  !!    if( L_ga_print )then
!! debug only  !!
!! debug only  !!        write(GA_print_unit,'(/A)') &
!! debug only  !!                      'GP_GA_opt: i_tree, i_node, GP_Indiv_Node_Param'
!! debug only  !!        do  i_tree=1,n_trees
!! debug only  !!            do  i_node=1,n_nodes
!! debug only  !!                !if( abs( GP_Individual_Node_Parameters(i_node,i_tree) ) > 1.0e-20 )then
!! debug only  !!                    !write(GA_print_unit,'(A,2(1x,I6),1x,E15.7)') &
!! debug only  !!                    !  'GP_GA_opt: i_tree, i_node, GP_Indiv_Node_Param', &
!! debug only  !!                    write(GA_print_unit,'(8x,2(1x,I6),1x,E15.7)') &
!! debug only  !!                                  i_tree, i_node, GP_Individual_Node_Parameters(i_node,i_tree)
!! debug only  !!                !endif ! abs( GP_Indiv_Node_Param(i_node,i_tree) ) > 1.0e-20
!! debug only  !!            enddo ! i_node
!! debug only  !!        enddo  ! i_tree
!! debug only  !!
!! debug only  !!        write(GA_print_unit,'(/A)') &
!! debug only  !!                     'GP_GA_opt: i_tree, i_node, GP_Indiv_Node_Type'
!! debug only  !!        do  i_tree=1,n_trees
!! debug only  !!            do  i_node=1,n_nodes
!! debug only  !!                if( GP_Individual_Node_Type(i_node,i_tree) > -9999 )then
!! debug only  !!                    write(GA_print_unit,'(8x,3(1x,I6))') &
!! debug only  !!                         i_tree, i_node, GP_Individual_Node_Type(i_node,i_tree)
!! debug only  !!                endif ! GP_Indiv_Node_Type(i_node,i_tree) > -9999
!! debug only  !!            enddo ! i_node
!! debug only  !!        enddo  ! i_tree
!! debug only  !!
!! debug only  !!        write(GA_print_unit,'(A)')' '
!! debug only  !!    endif ! L_ga_print
!! debug only  !!
!! debug only  !!endif ! myid == 0
!! debug only  !! jjm 20130417 <<<<<<<<<<<<<<<
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !!-----------------------------------------------------------------------------
!! debug only  !
!! debug only  !if( n_parameters .le. 0) then
!! debug only  !
!! debug only  !    if( L_GA_print )then
!! debug only  !        write(GA_print_unit,'(A)')        &
!! debug only  !              'GP_GA_opt: ERROR: n_parameters </= 0'
!! debug only  !        write(GA_print_unit,'(A,1x,I10)') &
!! debug only  !              'GP_GA_opt: n_parameters =   ', n_parameters
!! debug only  !    endif ! L_GA_print 
!! debug only  !
!! debug only  !    write(6,'(A)')        &
!! debug only  !          'GP_GA_opt: ERROR: n_parameters </= 0'
!! debug only  !    write(6,'(A,1x,I10)') &
!! debug only  !          'GP_GA_opt: n_parameters =   ', n_parameters
!! debug only  !
!! debug only  !    return ! 'n_par<=0'
!! debug only  !
!! debug only  !endif !   n_parameters .le. 0
!! debug only  !
!! debug only  !
!! debug only  !!if( n_time_steps .lt. n_parameters) then
!! debug only  !!    write(GA_print_unit,'(A)') &
!! debug only  !!          'GP_GA_opt: ERROR: n_time_steps < n_parameters'
!! debug only  !!    write(GA_print_unit,'(A,1x,I10)') &
!! debug only  !!          'GP_GA_opt:  n_time_steps = ', n_time_steps
!! debug only  !!    write(GA_print_unit,'(A,1x,I10)') &
!! debug only  !!          'GP_GA_opt: n_parameters =   ', n_parameters
!! debug only  !!    write(6,'(A)') &
!! debug only  !!          'GP_GA_opt: ERROR: n_time_steps < n_parameters'
!! debug only  !!    write(6,'(A,1x,I10)') &
!! debug only  !!          'GP_GA_opt:  n_time_steps = ', n_time_steps
!! debug only  !!    write(6,'(A,1x,I10)') &
!! debug only  !!          'GP_GA_opt: n_parameters =   ', n_parameters
!! debug only  !!    stop 'n_time < n_par'
!! debug only  !!endif ! n_time_steps .lt. n_parameters
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !child_parameters( 1:n_GP_parameters, 1:divider) = 0.0d0
!! debug only  !
!! debug only  !
!! debug only  !!-----------------------------------------------------------------------------
!! debug only  !
!! debug only  !! set up MPI process
!! debug only  !
!! debug only  !
!! debug only  !!if( myid == 0 )then
!! debug only  !!
!! debug only  !!    if( L_ga_print )then
!! debug only  !!        write(GA_print_unit,'(/A,1x,I10, 2(1x,I6)/)')&
!! debug only  !!          'GP_GA_opt: divider, n_parameters, n_procs ', &
!! debug only  !!                      divider, n_parameters, n_procs
!! debug only  !!    endif ! L_ga_print
!! debug only  !!
!! debug only  !!endif ! myid == 0
!! debug only  !
!! debug only  !!-----------------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !L_stop_run  = .FALSE.
!! debug only  !!L_stop_run  = .TRUE.
!! debug only  !
!! debug only  !do  i_GA_generation = 1, n_GA_Generations
!! debug only  !
!! debug only  !    ! Run_GA_lmdif determines if the new child
!! debug only  !    ! has to be sent to lmdif for 'local' optimization
!! debug only  !
!! debug only  !    Run_GA_lmdif=.false.
!! debug only  !
!! debug only  !    !if( myid == 0 )then
!! debug only  !    !    if( L_ga_print )then
!! debug only  !    !        write(GA_print_unit,'(/A,1x,I6,1x,A/)') &
!! debug only  !    !              'GA Generation ',i_GA_generation,' is underway'
!! debug only  !    !        flush(GA_print_unit)
!! debug only  !    !    endif ! L_ga_print
!! debug only  !    !    write(6,'(/A,1x,I6,1x,A/)') &
!! debug only  !    !              'GA Generation ',i_GA_generation,' is underway'
!! debug only  !    !    flush(6)
!! debug only  !    !endif ! myid == 0
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !    !if( myid == 0 )then
!! debug only  !    if( new_rank == 0 )then
!! debug only  !
!! debug only  !        if( i_GA_generation .eq. 1 ) then
!! debug only  !
!! debug only  !            ! on the first generation,
!! debug only  !            ! randomly create the initial individual parameter arrays
!! debug only  !            ! for each individual
!! debug only  !
!! debug only  !
!! debug only  !            ! sets:
!! debug only  !            !  child_parameters
!! debug only  !
!! debug only  !            call Initialize_GA_Child_Parameters( Child_Parameters )
!! debug only  !
!! debug only  !
!! debug only  !            Run_GA_lmdif=.true.
!! debug only  !
!! debug only  !            ! print child parameters at start of the generation
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    !write(GA_print_unit,'(/A,1x,I6)') &
!! debug only  !            !    write(6,'(/A,1x,I6)') &
!! debug only  !            !    'GP_GA_opt:1 child parameters at start of generation: ', &
!! debug only  !            !                                        i_GA_generation
!! debug only  !            !    do  i_ga_ind = 1, divider
!! debug only  !            !        !write(GA_print_unit,'(I6,4(1x,E24.16)/(4(1x,E24.16)))') &
!! debug only  !            !        write(6,'(I6,4(1x,E24.16)/(4(1x,E24.16)))') &
!! debug only  !            !              i_ga_ind, &
!! debug only  !            !              ( child_parameters(jj,i_ga_ind), jj = 1,n_parameters )
!! debug only  !            !    enddo ! i_ga_ind
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !            !flush(6)
!! debug only  !
!! debug only  !
!! debug only  !        else  ! i_GA_generation > 1
!! debug only  !
!! debug only  !
!! debug only  !            ! create the second 'generation' of parameter estimates using either:
!! debug only  !
!! debug only  !            !    i) save elites from last generation from being changed
!! debug only  !            !   ii) 'Fitness-Proportionate Reproduction;
!! debug only  !            !  iii) GA Crossover;
!! debug only  !            !   iv) GA Mutation
!! debug only  !
!! debug only  !            !-------------------------------------------------------------------------------
!! debug only  !
!! debug only  !            !   save the most fit individuals for the next generation
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(/A)')'GP_GA_opt: call GA_save_elites '
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !            ! uses:
!! debug only  !            ! individual_ranked_fitness
!! debug only  !
!! debug only  !            ! sets:
!! debug only  !            ! ga_individual_elites
!! debug only  !
!! debug only  !            call GA_save_elites( )
!! debug only  !
!! debug only  !            !-------------------------------------------------------------------------------
!! debug only  !
!! debug only  !            !!  replace the parameters of any individual with quality < 0 with new
!! debug only  !            !!  random numbers
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(/A)')'GP_GA_opt: call GA_replace_bad_individuals  '
!! debug only  !            !endif ! L_ga_print
!! debug only  !            !    write(6,'(/A)')'GP_GA_opt: call GA_replace_bad_individuals  '
!! debug only  !
!! debug only  !            !!call GA_replace_bad_individuals(Child_Parameters, individual_quality )
!! debug only  !
!! debug only  !            !-------------------------------------------------------------------------------
!! debug only  !
!! debug only  !            !   do initial "GA Fitness-Proportionate Reproduction"
!! debug only  !            !   to create a new population of children for all n_GA_individual
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A)')&
!! debug only  !            !          'GP_GA_opt: call GA_Fitness_Proportionate_Reproduction '
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !            ! uses:
!! debug only  !            !  individual_quality
!! debug only  !            !  Individual_Ranked_Fitness
!! debug only  !            !  Parent_Parameters
!! debug only  !
!! debug only  !            ! sets:
!! debug only  !            !  Run_GA_lmdif
!! debug only  !            !  individual_quality
!! debug only  !            !  Individual_Ranked_Fitness
!! debug only  !            !  Child_Parameters
!! debug only  !
!! debug only  !
!! debug only  !            call GA_Fitness_Proportionate_Reproduction(&
!! debug only  !                            Parent_Parameters,Child_Parameters, &
!! debug only  !                                                individual_quality )
!! debug only  !
!! debug only  !            !-------------------------------------------------------------------------------
!! debug only  !
!! debug only  !            !   do "GA Parameter Crossover" Operations Using Tournament-Style Selection
!! debug only  !            !   and randomly use it to replace the parents with the new children
!! debug only  !
!! debug only  !            if( n_GA_Crossovers .gt. 0) then
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,1x,I6)')&
!! debug only  !                !      'GP_GA_opt: call GA_Tournament_Style_Sexual_Repro  &
!! debug only  !                !      &n_GA_Crossovers',  n_GA_Crossovers
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !                ! uses:
!! debug only  !                !  Individual_Ranked_Fitness
!! debug only  !
!! debug only  !                ! sets:
!! debug only  !                !  Child_Parameters
!! debug only  !                !  Run_GA_lmdif
!! debug only  !                !  individual_quality
!! debug only  !
!! debug only  !                call GA_Tournament_Style_Sexual_Reproduction( &
!! debug only  !                            Parent_Parameters, Child_Parameters, individual_quality )
!! debug only  !
!! debug only  !            endif !   n_GA_Crossovers .gt. 0
!! debug only  !
!! debug only  !
!! debug only  !            !-------------------------------------------------------------------------------
!! debug only  !
!! debug only  !            !   do "GA Parameter Mutation" Operations
!! debug only  !
!! debug only  !            !   select a random individual and put a new random number into one of
!! debug only  !            !   its parameters
!! debug only  !
!! debug only  !
!! debug only  !            if( n_GA_Mutations .gt. 0) then
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,1x,I6)')&
!! debug only  !                !      'GP_GA_opt: call GA_Mutations  n_GA_Mutations',  n_GA_Mutations
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !                ! uses:
!! debug only  !
!! debug only  !                ! sets:
!! debug only  !                !  child_parameters
!! debug only  !                !  Run_GA_lmdif
!! debug only  !                !  individual_quality
!! debug only  !
!! debug only  !                call GA_Mutations( Child_Parameters, individual_quality )
!! debug only  !
!! debug only  !            endif !   n_GA_Mutations .gt. 0
!! debug only  !
!! debug only  !
!! debug only  !            !-------------------------------------------------------------------------------
!! debug only  !
!! debug only  !            !   do "GA Parameter rand_replace" Operations
!! debug only  !
!! debug only  !            !   select a random, non-elite individual and put new random numbers into
!! debug only  !            !   its parameters
!! debug only  !
!! debug only  !            if( n_GA_rand_replaces > 0) then
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,1x,I6)')&
!! debug only  !                !      'GP_GA_opt: call GA_rand_replace  n_GA_rand_replaces',  &
!! debug only  !                !                                        n_GA_rand_replaces
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !                ! uses:
!! debug only  !
!! debug only  !                ! sets:
!! debug only  !                !  child_parameters
!! debug only  !                !  Run_GA_lmdif
!! debug only  !                !  individual_quality
!! debug only  !
!! debug only  !
!! debug only  !                call GA_random_replace( Child_Parameters, individual_quality )
!! debug only  !
!! debug only  !
!! debug only  !            endif !   n_GA_rand_replaces .gt. 0
!! debug only  !
!! debug only  !            !-------------------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !            ! print child parameters at start of the generation
!! debug only  !
!! debug only  !            !if( i_GA_generation == n_GA_generations               .or.  &
!! debug only  !            !    mod( i_GA_generation, GA_child_print_interval ) == 0      )then
!! debug only  !
!! debug only  !            !if( i_GA_generation == n_GA_generations  )then
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(/A,1x,I6)') &
!! debug only  !                !    'GP_GA_opt:2 child parameters at start of GA generation:', &
!! debug only  !                !                                            i_GA_generation
!! debug only  !                !    do  i_ga_ind = 1, divider
!! debug only  !                !        write(GA_print_unit,'(I6,1x,10(1x,E15.7)/(10(1x,E15.7)))') &
!! debug only  !                !              i_ga_ind, &
!! debug only  !                !             (child_parameters(jj, i_ga_ind),&
!! debug only  !                !                               jj = 1,n_parameters )
!! debug only  !                !    enddo ! i_ga_ind
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !            !endif ! i_GA_generation == n_GA_generations ...
!! debug only  !
!! debug only  !
!! debug only  !        endif ! i_GA_generation .eq. 1
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !    endif ! new_rank == 0  ! myid == 0
!! debug only  !
!! debug only  !
!! debug only  !    !------------------------------------------------------------------------
!! debug only  !
!! debug only  !    !  broadcast child parameters
!! debug only  !
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(/A,2(1x,I6))') &
!! debug only  !    !    'GP_GA_opt:  broadcast child parameters myid, i_GA_generation ', &
!! debug only  !    !                                            myid, i_GA_generation
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !    child_number =  divider * n_GP_parameters
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(A,4(1x,I6)/)') &
!! debug only  !    !    'GP_GA_opt: myid, divider, n_GP_parameters, child_number =', &
!! debug only  !    !                myid, divider, n_GP_parameters, child_number
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !    call MPI_BCAST( Child_Parameters,  child_number,    &
!! debug only  !                    MPI_DOUBLE_PRECISION, 0, new_comm, ierr )
!! debug only  !
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(/A,2(1x,I10)/)') &
!! debug only  !    !     'GP_GA_opt: child  broadcast myid, ierr = ', myid, ierr
!! debug only  !    !    write(GA_print_unit,'(/A,2(1x,I10)/)') &
!! debug only  !    !     'GP_GA_opt: myid, divider = ', myid, divider
!! debug only  !    !    do  i_ga_ind = 1, divider
!! debug only  !    !        write(GA_print_unit,'(I3,1x,I3,1x,12(1x,E15.7))') &
!! debug only  !    !              myid, i_ga_ind, &
!! debug only  !    !              ( child_parameters(jj,i_ga_ind), jj = 1,n_parameters )
!! debug only  !    !    enddo ! i_ga_ind
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !    !------------------------------------------------------------------------
!! debug only  !
!! debug only  !    ! broadcast Run_GA_lmdif
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(/A,1x,I6)') &
!! debug only  !    !      'GP_GA_opt:  broadcast Run_GA_lmdif i_GA_generation ', i_GA_generation
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !    call MPI_BCAST( Run_GA_lmdif,  divider,    &
!! debug only  !                        MPI_LOGICAL, 0, new_comm, ierr )
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(A,1x,I10/)') &
!! debug only  !    !     'GP_GA_opt: Run_GA_lmdif  broadcast ierr = ', ierr
!! debug only  !    !    write(GA_print_unit,'(/A,1x,I6,/(10(2x,L1)))') &
!! debug only  !    !     'GP_GA_opt: myid, Run_GA_lmdif  ', myid, Run_GA_lmdif
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!! debug only  !    !   begin RK fcn integration segment
!! debug only  !    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!! debug only  !
!! debug only  !
!! debug only  !    !if( myid == 0 )then
!! debug only  !    !    if( L_ga_print )then
!! debug only  !    !        write(GA_print_unit,'(/A,1x,I4,1x,A,1x,I6/)') &
!! debug only  !    !          'GP_GA_opt: begin RK fcn integration   i_GA_generation =', &
!! debug only  !    !                                                 i_GA_generation,    &
!! debug only  !    !                                         '       divider =', &
!! debug only  !    !                                                 divider
!! debug only  !    !    endif ! L_ga_print
!! debug only  !    !endif !  myid == 0
!! debug only  !
!! debug only  !
!! debug only  !    individual_quality( 1: divider ) = 1
!! debug only  !
!! debug only  !    nsafe = 0
!! debug only  !
!! debug only  !    !if( myid == 0  )then
!! debug only  !    if( new_rank == 0  )then
!! debug only  !
!! debug only  !        ! processor 0 sends job assignments to N  processors
!! debug only  !        ! where N is the smaller of (number of processors -1)
!! debug only  !        !                    and the number of individuals
!! debug only  !
!! debug only  !        ! numsent is the number of messages sent up to now
!! debug only  !
!! debug only  !        numsent = 0
!! debug only  !        i_ga_ind = 0
!! debug only  !
!! debug only  !        do  isource = 1, min( n_procs-1, divider )
!! debug only  !
!! debug only  !
!! debug only  !            i_ga_ind = i_ga_ind + 1
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,1x,I6, 4x, L1)') &
!! debug only  !            !     'GP_GA_opt:1 494 i_ga_ind, Run_GA_lmdif(i_ga_ind)  ', &
!! debug only  !            !                      i_ga_ind, Run_GA_lmdif(i_ga_ind)
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !            call MPI_SEND( i_dummy,  1, MPI_INTEGER,    &
!! debug only  !                           isource, isource,  new_comm, ierr )
!! debug only  !
!! debug only  !            numsent = numsent + 1
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,4(1x,I6))') &
!! debug only  !            !     'GP_GA_opt:1 504 myid, isource, i_ga_ind, numsent ', &
!! debug only  !            !                      myid, isource, i_ga_ind, numsent
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !        enddo ! isource
!! debug only  !
!! debug only  !
!! debug only  !        ! at this point i_ga_ind = numsent
!! debug only  !
!! debug only  !
!! debug only  !        !if( L_ga_print )then
!! debug only  !        !    write(GA_print_unit,'(A,4(1x,I6))') &
!! debug only  !        !         'GP_GA_opt: aft source loop 1 myid, i_ga_ind, numsent ', &
!! debug only  !        !                                       myid, i_ga_ind, numsent
!! debug only  !        !endif ! L_ga_print
!! debug only  !
!! debug only  !        !-------------------------------------------------------------------------------------
!! debug only  !
!! debug only  !        ! processor 0 loops over the number of individuals and waits for a message
!! debug only  !        ! from the other processors
!! debug only  !
!! debug only  !        do  isource = 1, divider
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,3(1x,I6))') &
!! debug only  !            !     'GP_GA_opt:2 521 myid, isource, numsent ', &
!! debug only  !            !                      myid, isource, numsent
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !            buffer_recv = 0.0d0
!! debug only  !            call MPI_RECV( buffer_recv, n_GP_parameters+2, &
!! debug only  !                           MPI_DOUBLE_PRECISION, &
!! debug only  !                           MPI_ANY_SOURCE, MPI_ANY_TAG,  &
!! debug only  !                           new_comm, MPI_STAT,  ierr )
!! debug only  !
!! debug only  !            sender       = MPI_STAT( MPI_SOURCE )
!! debug only  !            i_individual = MPI_STAT( MPI_TAG )
!! debug only  !
!! debug only  !
!! debug only  !            ! received a message from processor "sender" which processed
!! debug only  !            ! individual "i_individual"
!! debug only  !
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,5(1x,I6))') &
!! debug only  !            !     'GP_GA_opt:2 529 myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )', &
!! debug only  !            !                      myid, MPI_STAT( MPI_SOURCE ), MPI_STAT( MPI_TAG )
!! debug only  !            !    write(GA_print_unit,'(A,5(1x,I6))') &
!! debug only  !            !     'GP_GA_opt:2 529 myid, isource, numsent, sender, i_individual ', &
!! debug only  !            !                      myid, isource, numsent, sender, i_individual
!! debug only  !            !    write(GA_print_unit,'(A,2(1x,I6),4x,L1)') &
!! debug only  !            !      'GP_GA_opt:2 myid, i_individual, Run_GA_lmdif(i_individual)', &
!! debug only  !            !                   myid, i_individual, Run_GA_lmdif(i_individual)
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !            ! store the information received in the above message
!! debug only  !
!! debug only  !            if( Run_GA_lmdif(i_individual) ) then
!! debug only  !
!! debug only  !
!! debug only  !                do  jj = 1, n_GP_parameters
!! debug only  !                    child_parameters(jj,i_individual) =  buffer_recv(jj)
!! debug only  !                enddo ! jj
!! debug only  !
!! debug only  !                individual_SSE(i_individual)     =       buffer_recv( n_GP_parameters+1)
!! debug only  !                individual_quality(i_individual) = nint( buffer_recv( n_GP_parameters+2) )
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,3(1x,I6))') &
!! debug only  !                !     'GP_GA_opt:2 554 myid, n_parameters, i_individual', &
!! debug only  !                !                      myid, n_parameters, i_individual
!! debug only  !                !    !write(GA_print_unit,'(A,1x,I6,1x, A/(5(1x,E15.7)))') &
!! debug only  !                !    ! 'GP_GA_opt:2 myid = ',myid, &
!! debug only  !                !    ! ' child_parameters(1:n_parameters,i_individual)', &
!! debug only  !                !    !   child_parameters(1:n_parameters,i_individual)
!! debug only  !                !    write(GA_print_unit,'(A,2(1x,I6),1x,E15.7)') &
!! debug only  !                !     'GP_GA_opt:2 myid, i_individual, individual_SSE(i_individual)', &
!! debug only  !                !                  myid, i_individual, individual_SSE(i_individual)
!! debug only  !                !    write(GA_print_unit,'(A,3(1x,I6))') &
!! debug only  !                !     'GP_GA_opt:2 myid, i_individual, individual_quality(i_individual)', &
!! debug only  !                !                  myid, i_individual, individual_quality(i_individual)
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !            endif ! Run_GA_lmdif(i_individual)
!! debug only  !
!! debug only  !            !--------------------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,1x,I6, 4x,L1)') &
!! debug only  !            !     'GP_GA_opt:2 542 myid, numsent < divider ', &
!! debug only  !            !                      myid, numsent < divider
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !            ! check to see if all individuals have been processed
!! debug only  !
!! debug only  !            if( numsent <  divider )then
!! debug only  !
!! debug only  !                ! numsent <  divider
!! debug only  !                ! means not all individuals have been processed
!! debug only  !
!! debug only  !                ! send a message to the processor "sender"
!! debug only  !                ! which just sent a message saying it has
!! debug only  !                ! completed an individual, and tell it to process
!! debug only  !                ! the individual "i_ga_ind" as the  "numsent+1"  task
!! debug only  !
!! debug only  !                i_ga_ind = i_ga_ind + 1
!! debug only  !
!! debug only  !                call MPI_SEND( i_ga_ind, 1, MPI_INTEGER,    &
!! debug only  !                               sender, numsent+1,  new_comm, ierr )
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,4(1x,I6))') &
!! debug only  !                !      'GP_GA_opt:2 554 send myid, sender, numsent, i_ga_ind', &
!! debug only  !                !                            myid, sender, numsent, i_ga_ind
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !                ! just sent a new task, so increment the number sent
!! debug only  !
!! debug only  !                numsent = numsent + 1
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,4(1x,I6))') &
!! debug only  !                !     'GP_GA_opt:2 556  myid, sender, numsent, i_ga_ind ', &
!! debug only  !                !                       myid, sender, numsent, i_ga_ind
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !            else
!! debug only  !
!! debug only  !                ! DONE !
!! debug only  !
!! debug only  !                ! number of tasks sent out is >= number of individuals, so
!! debug only  !                ! all the work has been completed
!! debug only  !
!! debug only  !                ! tell the "sender" processor that it is done and
!! debug only  !                ! send it a message to stop
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,3(1x,I6))') &
!! debug only  !                !      'GP_GA_opt:2 send msg to stop  myid, numsent, i_ga_ind ', &
!! debug only  !                !                                     myid, numsent, i_ga_ind
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !                call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
!! debug only  !                               sender, 0,  new_comm, ierr )
!! debug only  !
!! debug only  !
!! debug only  !            endif ! numsent
!! debug only  !
!! debug only  !        enddo ! isource
!! debug only  !
!! debug only  !
!! debug only  !        !----------------------------------------------------------------------
!! debug only  !
!! debug only  !        ! this section takes care of the case where there are fewer GA individuals
!! debug only  !        ! than (number of procs) - 1
!! debug only  !
!! debug only  !        ! without the code below,  the program hangs because the processors
!! debug only  !        ! with numbers  (divider+1)  to (n_procs-1)
!! debug only  !        ! are waiting for a signal to stop
!! debug only  !        ! and that is never going to be sent from the loop above.
!! debug only  !
!! debug only  !        ! so when the above loop is finished, send a stop signal to the unused
!! debug only  !        ! processors so the program can continue
!! debug only  !
!! debug only  !        if( divider < n_procs -1 )then
!! debug only  !
!! debug only  !            do  i = divider+1, n_procs-1
!! debug only  !
!! debug only  !                call MPI_SEND( MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,    &
!! debug only  !                               i , 0, new_comm , ierr )
!! debug only  !            enddo ! i
!! debug only  !
!! debug only  !        endif ! divider < n_procs -1
!! debug only  !
!! debug only  !        !----------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !    else  ! not new_rank == 0   ! myid == 0
!! debug only  !
!! debug only  !
!! debug only  !        ! code for processors 1 - ( divider )  ! - 1 )
!! debug only  !
!! debug only  !
!! debug only  !        ! these processors wait until a message is received from
!! debug only  !        ! processor 0 telling it to process
!! debug only  !        ! the individual named in the message tag = MPI_STAT( MPI_TAG )
!! debug only  !
!! debug only  !        recv_loop:&
!! debug only  !        do
!! debug only  !
!! debug only  !            !if( myid == 1 )then
!! debug only  !            !    write(GA_print_unit,'(A,1x,I6)') &
!! debug only  !            !          'GP_GA_opt: myid  proc1 in recv_loop', myid
!! debug only  !            !endif
!! debug only  !
!! debug only  !
!! debug only  !            call MPI_RECV( i_dummy, 1, MPI_INTEGER,    &
!! debug only  !                           0, MPI_ANY_TAG,  new_comm , MPI_STAT, ierr )
!! debug only  !
!! debug only  !            !---------------------------------------------------------------
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,2(1x,I6))') &
!! debug only  !            !      'GP_GA_opt:3  myid, MPI_STAT( MPI_TAG ) ', &
!! debug only  !            !                    myid, MPI_STAT( MPI_TAG )
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !            ! was a stop signal received ?
!! debug only  !
!! debug only  !            ! if the tag is <= 0, this is a stop signal
!! debug only  !
!! debug only  !            if( MPI_STAT( MPI_TAG ) <= 0 ) exit recv_loop
!! debug only  !
!! debug only  !            !---------------------------------------------------------------
!! debug only  !
!! debug only  !            ! process the individual named in the message tag
!! debug only  !
!! debug only  !            i_2_individual = MPI_STAT( MPI_TAG )
!! debug only  !
!! debug only  !            buffer = 0.0D0
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,4(1x,I6))') &
!! debug only  !            !    'GP_GA_opt:3 myid, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual', &
!! debug only  !            !                 myid, i_dummy, MPI_STAT( MPI_TAG ), i_2_individual
!! debug only  !            !    write(GA_print_unit,'(A,2(1x,I6),4x,L1)') &
!! debug only  !            !    'GP_GA_opt:3 myid, i_2_individual, Run_GA_lmdif(i_2_individual)', &
!! debug only  !            !                 myid, i_2_individual, Run_GA_lmdif(i_2_individual)
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !            if( Run_GA_lmdif(i_2_individual)) then
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,2(1x,I6))') &
!! debug only  !                !      'GP_GA_opt:3 call setup_run_fcn  myid, i_2_individual', &
!! debug only  !                !                                       myid, i_2_individual
!! debug only  !                !endif ! L_ga_print
!! debug only  !                !write(6,'(A,2(1x,I6))') &
!! debug only  !                !      'GP_GA_opt:3 call setup_run_fcn  myid, i_2_individual', &
!! debug only  !                !                                       myid, i_2_individual
!! debug only  !                !flush(6)
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !                ! do the Runge-Kutta integration for individual i_2_individual
!! debug only  !
!! debug only  !                ! uses:
!! debug only  !                !  child_parameters
!! debug only  !
!! debug only  !                ! sets:
!! debug only  !                !  individual_quality
!! debug only  !                !  individual_SSE
!! debug only  !                !  child_parameters
!! debug only  !
!! debug only  !
!! debug only  !                call setup_run_fcn( i_2_individual, &
!! debug only  !                                    child_parameters,individual_quality, &
!! debug only  !                                    new_group, new_comm )
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !                !if( L_ga_print )then
!! debug only  !                !    write(GA_print_unit,'(A,3(1x,I6))') &
!! debug only  !                !      'GP_GA_opt:3 AFTER call setup_run_fcn myid, i_2_individual', &
!! debug only  !                !                                            myid, i_2_individual
!! debug only  !                !    flush(GA_print_unit)
!! debug only  !                !endif ! L_ga_print
!! debug only  !
!! debug only  !                !write(6,'(A,3(1x,I6))') &
!! debug only  !                !  'GP_GA_opt:3 AFTER call setup_run_fcn myid, i_2_individual', &
!! debug only  !                !                                        myid, i_2_individual
!! debug only  !                !flush(6)
!! debug only  !
!! debug only  !                !-------------------------------------------------------------------------
!! debug only  !
!! debug only  !                do  jj = 1, n_GP_parameters
!! debug only  !                    buffer(jj) =  child_parameters(jj, i_2_individual)
!! debug only  !                enddo ! jj
!! debug only  !
!! debug only  !                buffer(n_GP_parameters+1) = &
!! debug only  !                      individual_SSE(i_2_individual)
!! debug only  !                buffer(n_GP_parameters+2) = &
!! debug only  !                      real( individual_quality(i_2_individual), kind=8 )
!! debug only  !
!! debug only  !            endif !  Run_GA_lmdif(i_2_individual)
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,3(1x,I6))') &
!! debug only  !            !    'GP_GA_opt:3 705 myid, n_GP_parameters, i_2_individual', &
!! debug only  !            !                     myid, n_GP_parameters, i_2_individual
!! debug only  !            !    write(GA_print_unit,'(A/(5(1x,E24.16)))') &
!! debug only  !            !    'GP_GA_opt:3 child_parameters(1:n_parameters,i_2_individual)', &
!! debug only  !            !                 child_parameters(1:n_parameters,i_2_individual)
!! debug only  !            !    write(GA_print_unit,'(A,3(1x,I6))') &
!! debug only  !            !    'GP_GA_opt:3 myid, i_2_individual, individual_quality(i_2_individual)', &
!! debug only  !            !                 myid, i_2_individual, individual_quality(i_2_individual)
!! debug only  !            !    write(GA_print_unit,'(A,2(1x,I6), 1x, E15.7)') &
!! debug only  !            !    'GP_GA_opt:3 send results &
!! debug only  !            !      &myid, i_2_individual, individual_SSE(i_2_individual)', &
!! debug only  !            !       myid, i_2_individual, individual_SSE(i_2_individual)
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !            ! send the R-K integration results
!! debug only  !            ! for individual i_2_individual to processor 0
!! debug only  !
!! debug only  !            call MPI_SEND( buffer, n_GP_parameters+2,  &
!! debug only  !                           MPI_DOUBLE_PRECISION, 0, &
!! debug only  !                           i_2_individual, new_comm, ierr )
!! debug only  !
!! debug only  !
!! debug only  !            !---------------------------------------------------------------
!! debug only  !
!! debug only  !            ! code to ensure that an error does not allow this loop to run forever
!! debug only  !
!! debug only  !            nsafe = nsafe + 1
!! debug only  !
!! debug only  !            !if( L_ga_print )then
!! debug only  !            !    write(GA_print_unit,'(A,2(1x,I10))') &
!! debug only  !            !      'GP_GA_opt: myid, nsafe =', myid, nsafe
!! debug only  !            !endif ! L_ga_print
!! debug only  !
!! debug only  !            if( nsafe > 100 * divider ) then
!! debug only  !
!! debug only  !                if( L_GA_print )then
!! debug only  !                    write(GA_print_unit,'(A,1x,I10)') &
!! debug only  !                      'GP_GA_opt: too many iterations  nsafe =', nsafe
!! debug only  !                    flush(GA_print_unit) 
!! debug only  !                endif ! L_GA_print 
!! debug only  !
!! debug only  !                write(6,'(A,1x,I10)') &
!! debug only  !                  'GP_GA_opt: too many iterations  nsafe =', nsafe
!! debug only  !                flush(6) 
!! debug only  !
!! debug only  !                L_too_many_iters = .TRUE.  
!! debug only  !                exit recv_loop
!! debug only  !
!! debug only  !            endif ! nsafe
!! debug only  !
!! debug only  !            !---------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !         enddo  recv_loop
!! debug only  !
!! debug only  !
!! debug only  !    endif ! myid == 0
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(A,2(1x,I6))') &
!! debug only  !    !      'GP_GA_opt: after recv_loop  myid = ', myid
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !    if( L_too_many_iters )then
!! debug only  !        call MPI_FINALIZE(ierr)
!! debug only  !        stop 'bad nsafe'
!! debug only  !    endif !  L_too_many_iters 
!! debug only  !
!! debug only  !    !-------------------------------------------------------------------
!! debug only  !
!! debug only  !    ! wait until all divider have been processed
!! debug only  !
!! debug only  !    call MPI_BARRIER( new_comm, ierr )
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(A,2(1x,I6))') &
!! debug only  !    !      'GP_GA_opt: after barrier 2 i_GA_generation, myid = ', &
!! debug only  !    !                                  i_GA_generation, myid
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !    !-------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !    !  calculate the fitness for this generation
!! debug only  !
!! debug only  !
!! debug only  !    if( new_rank == 0  )then
!! debug only  !
!! debug only  !        !if( L_ga_print )then
!! debug only  !        !    write(GA_print_unit,'(A)')  'GP_GA_opt: individual_SSE  '
!! debug only  !        !    write(GA_print_unit,'(5(1x,E12.5))')  individual_SSE(1:divider)
!! debug only  !        !    write(GA_print_unit,'(/A,1x,I6)') &
!! debug only  !        !          'GP_GA_opt: call calc_fitness i_GA_generation ', &
!! debug only  !        !                                        i_GA_generation
!! debug only  !        !    write(GA_print_unit,'(/A,6(1x,E15.7))') &
!! debug only  !        !          'GP_GA_opt: GP_Individual_Initial_Conditions ', &
!! debug only  !        !                      GP_Individual_Initial_Conditions(1:n_code_equations)
!! debug only  !        !endif ! L_ga_print
!! debug only  !        !write(6,'(/A,1x,I6)') &
!! debug only  !        !      'GP_GA_opt: call calc_fitness i_GA_generation ', &
!! debug only  !        !                                    i_GA_generation
!! debug only  !        !flush(6)
!! debug only  !
!! debug only  !        ! uses:
!! debug only  !        !  child_parameters
!! debug only  !        !  individual_quality
!! debug only  !        !  individual_SSE
!! debug only  !
!! debug only  !        ! sets:
!! debug only  !        !  individual_quality
!! debug only  !        !  individual_ranked_fitness
!! debug only  !        !  integrated_SSE
!! debug only  !        !  integrated_ranked_fitness
!! debug only  !
!! debug only  !
!! debug only  !        call calc_fitness( child_parameters, individual_quality, &
!! debug only  !                           i_GA_Best_Parent, Parent_Parameters, L_stop_run, &
!! debug only  !                           i_GP_Generation, i_GP_individual, &
!! debug only  !                           new_group, new_comm  )
!! debug only  !
!! debug only  !
!! debug only  !        !if( L_ga_print )then
!! debug only  !        !    write(GA_print_unit,'(/A,1x,I6/)') &
!! debug only  !        !          'GP_GA_opt: aft call calc_fitness i_GA_generation ', &
!! debug only  !        !                                            i_GA_generation
!! debug only  !        !endif ! L_ga_print
!! debug only  !        !write(6,'(/A,1x,I6/)') &
!! debug only  !        !       'GP_GA_opt: aft call calc_fitness i_GA_generation ', &
!! debug only  !        !                                         i_GA_generation
!! debug only  !        !flush(6)
!! debug only  !
!! debug only  !        !---------------------------------------------------------------------
!! debug only  !
!! debug only  !        ! this subroutine, if called, will run lmdif on the best individual
!! debug only  !        ! for this generation, instead of waiting to run lmdif after the
!! debug only  !        ! generation loop
!! debug only  !
!! debug only  !        !    call select_best_RK_lmdif_result( &
!! debug only  !        !                i_GP_Generation,i_GP_individual, &
!! debug only  !        !                i_GA_best_parent, parent_parameters, &
!! debug only  !        !                child_parameters, &
!! debug only  !        !                L_stop_run  )
!! debug only  !        !                !individual_quality, L_stop_run  )
!! debug only  !
!! debug only  !        !---------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !    endif ! new_rank == 0
!! debug only  !
!! debug only  !    !if( i_GA_generation > 1 )then
!! debug only  !    !    if( L_ga_print )then
!! debug only  !    !        write(GA_print_unit,'(A,2(1x,I6))') &
!! debug only  !    !          'GP_GA_opt: at stop  i_GA_generation, myid = ', &
!! debug only  !    !                               i_GA_generation, myid
!! debug only  !    !        call MPI_FINALIZE(ierr) ! debug only
!! debug only  !    !        stop ! debug only
!! debug only  !    !    endif ! L_ga_print
!! debug only  !    !endif ! i_GA_generation > 1
!! debug only  !
!! debug only  !    !-------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !    call MPI_BCAST( L_stop_run,  1,    &
!! debug only  !                    MPI_LOGICAL, 0, new_comm, ierr )
!! debug only  !
!! debug only  !    if( L_stop_run )then
!! debug only  !
!! debug only  !        !if( L_ga_print )then
!! debug only  !        !    write(GA_print_unit,'(/A,1x,I6/)') &
!! debug only  !        !    'GP_GA_opt: L_stop_fitness is true &
!! debug only  !        !     &so last generation is i_GA_generation ', &
!! debug only  !        !                            i_GA_generation
!! debug only  !        !endif ! L_ga_print
!! debug only  !
!! debug only  !        i_GA_generation_last = i_GA_generation
!! debug only  !
!! debug only  !        exit
!! debug only  !
!! debug only  !    endif ! L_stop_run
!! debug only  !
!! debug only  !
!! debug only  !enddo  ! i_generation
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !
!! debug only  !!if( L_ga_print )then
!! debug only  !!    write(GA_print_unit,'(//A,1x,I6/)') &
!! debug only  !!      'GP_GA_opt:  finished the loop on i_GA_generation  new_rank =', &
!! debug only  !!                                                         new_rank
!! debug only  !!endif ! L_ga_print
!! debug only  !
!! debug only  !!----------------------------------------------------------------------
!! debug only  !
!! debug only  !! wait until all processors have finished the generation loop
!! debug only  !
!! debug only  !call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?
!! debug only  !
!! debug only  !!if( L_ga_print )then
!! debug only  !!    write(GA_print_unit,'(A,1x,I6)') &
!! debug only  !!          'GP_GA_opt: after barrier 3 myid = ', myid
!! debug only  !!endif ! L_ga_print
!! debug only  !
!! debug only  !!----------------------------------------------------------------------
!! debug only  !
!! debug only  !
!! debug only  !! finished all generations,
!! debug only  !
!! debug only  !! now call lmdif on the best individual of the last generation
!! debug only  !! and determine if lmdif has improved the fitness of this individual
!! debug only  !! then save the parameters of the fitter of the two results, the RK result
!! debug only  !! and the lmdif result
!! debug only  !
!! debug only  !! GP_Individual_Node_Parameters loaded from the best parent parameters (RK or lmdif)
!! debug only  !! after the RK process and lmdif
!! debug only  !
!! debug only  !
!! debug only  !if( new_rank == 0  )then
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(6,'(//A/A)') &
!! debug only  !    !      'GP_GA_opt: finished all generations', &
!! debug only  !    !      'GP_GA_opt: call select_best_RK_lmdif_result'
!! debug only  !    !    flush(6)
!! debug only  !    !endif ! L_ga_print
!! debug only  !
!! debug only  !    ! uses:
!! debug only  !
!! debug only  !    ! sets:
!! debug only  !
!! debug only  !    call select_best_RK_lmdif_result( &
!! debug only  !                i_GP_Generation,i_GP_individual, &
!! debug only  !                i_GA_best_parent, parent_parameters, &
!! debug only  !                child_parameters, &
!! debug only  !                L_stop_run, new_group, new_comm  )
!! debug only  !
!! debug only  !    !if( L_ga_print )then
!! debug only  !    !    write(GA_print_unit,'(/A//)') &
!! debug only  !    !      'GP_GA_opt: aft call select_best_RK_lmdif_result '
!! debug only  !    !endif ! L_ga_print
!! debug only  !    !write(6,'(/A//)') &
!! debug only  !    !     'GP_GA_opt: aft call select_best_RK_lmdif_result '
!! debug only  !    !flush(6)
!! debug only  !
!! debug only  !
!! debug only  !endif ! myid == 0


!------------------------------------------------------------------------

! broadcast individual_fitness


!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: broadcast individual_fitness myid = ', myid
!endif ! L_ga_print
!write(6,'(A,1x,I5,1x,I5)') 'GP_GA_opt:2  myid, new_rank ', myid, new_rank

!individual_fitness = real( i_GP_individual * new_rank, kind=8)  ! debug only
!individual_fitness = real( new_rank, kind=8)  ! debug only
individual_fitness = real( myid, kind=8)  ! debug only

write(6,'(/A,2(1x,I6), 1x,F10.2/)') &
     'GP_GA_opt: broadcast individual_fitness  myid, new_rank, individual_fitness', &
                                               myid, new_rank, individual_fitness

message_len = 1
call MPI_REDUCE( individual_fitness, sum_if, message_len,    &
                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, new_comm, ierr )



!call MPI_BCAST( individual_fitness, message_len,    &
!                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )
!                !MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: aft broadcast individual_fitness  ierr = ', ierr
if( new_rank == 0 )then
write(6,'(/A,2(1x,I6), 1x,F10.2/)') &
     'GP_GA_opt: aft REDUCE myid, new_rank, sum_if', &
                            myid, new_rank, sum_if
endif 
!endif ! L_ga_print


!------------------------------------------------------------------------

! broadcast Individual_SSE_best_parent


!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: broadcast Individual_SSE_best_parent myid = ', myid
!endif ! L_ga_print

!! debug only  !message_len = 1
!! debug only  !call MPI_BCAST( Individual_SSE_best_parent, message_len,    &
!! debug only  !                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: aft broadcast Individual_SSE_best_parent  ierr = ', ierr
!endif ! L_ga_print


!------------------------------------------------------------------------

! broadcast GP_Individual_Node_Parameters


!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: broadcast GP_Individual_Node_Parameters  myid = ', myid
!endif ! L_ga_print

!! debug only  !message_len = n_trees * n_nodes
!! debug only  !
!! debug only  !call MPI_BCAST( GP_Individual_Node_Parameters, message_len,    &
!! debug only  !                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: aft broadcast GP_Individual_Node_Parameters  ierr = ', ierr
!endif ! L_ga_print


!------------------------------------------------------------------------

! broadcast GP_Individual_Initial_Conditions


!if( L_ga_print )then
!    write(GA_print_unit,'(/A,1x,I6)') &
!     'GP_GA_opt: broadcast GP_Individual_Initial_Conditions myid = ', myid
!endif ! L_ga_print

!! debug only  !message_len = n_CODE_equations
!! debug only  !
!! debug only  !call MPI_BCAST( GP_Individual_Initial_Conditions, message_len,    &
!! debug only  !                MPI_DOUBLE_PRECISION, 0, new_comm, ierr )

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
