subroutine bcast1()


use mpi
use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module
use Runge_Kutta_Variables_module


implicit none




!----------------------------------------------------------------------------------------

! broadcast the values read in by cpu 0 to others

call MPI_BCAST( GA_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_save_elites_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_GA_Generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_GA_Individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_time_steps, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( dt, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( model, 30,     &
                MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
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
call MPI_BCAST( n_gp_individuals, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( n_gp_generations, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

call MPI_BCAST( n_node_functions, 1,    &
                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )

return

end subroutine bcast1