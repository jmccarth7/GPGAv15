subroutine bcast1()


use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use GP_variables_module


implicit none

integer :: array_len
!integer :: i

!----------------------------------------------------------------------------------------

call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary?

! broadcast the values read in by cpu 0 to others

call MPI_BCAST( GA_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_save_elites_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_rand_replace_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( dt, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Tree_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Elitist_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Asexual_Reproduction_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GP_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )


call MPI_BCAST( sse_min_time, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( sse_max_time, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( sse_low_wt, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------

call MPI_BCAST( model, 30,     &
                MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------


call MPI_BCAST( n_GA_Generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_GA_Individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_time_steps, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_gp_individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_gp_generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_node_functions, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( user_input_random_seed, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( number_GA_child_prints, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( number_GP_child_prints, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_levels, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_partitions, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
!-----------------------------------------------------------------

!write(6,'(A,2(1x,I6))') 'bc1: myid, n_seed = ', myid, n_seed

!if( myid == 0 )then
!    write(6,'(A,5x,L1)') 'bc1: L_restart', L_restart
!endif ! myid == 0 


call MPI_BCAST( temp_seed, n_seed,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( L_restart, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------

call MPI_BCAST( L_GA_print, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

!call MPI_BCAST( L_GP_print, 1,    &
!                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( L_GP_log, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( L_unit50_output, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( L_GP_output_parameters, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )

!-----------------------------------------------------------------
!-----------------------------------------------------------------

return

end subroutine bcast1
