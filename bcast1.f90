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

integer,dimension( 7 )      :: integer_buffer
real(kind=8),dimension( 9 ) :: real_buffer


!----------------------------------------------------------------------------------------

call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary? 

! broadcast the values read in by cpu 0 to others

call MPI_BCAST( GA_Crossover_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_Mutation_Probability, 1,    &
                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( GA_save_elites_Probability, 1,    &
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

!-----------------------------------------------------------------
!if( myid == 0 )then
!
!    real_buffer(1) = GA_Crossover_Probability
!    real_buffer(2) = GA_Mutation_Probability
!    real_buffer(3) = GA_save_elites_Probability
!    real_buffer(4) = dt
!    real_buffer(5) = GP_Tree_Probability
!    real_buffer(6) = GP_Elitist_Probability
!    real_buffer(7) = GP_Asexual_Reproduction_Probability
!    real_buffer(8) = GP_Crossover_Probability
!    real_buffer(9) = GP_Mutation_Probability
!
!endif ! myid == 0 
!                
!
!call MPI_BCAST( real_buffer, 9,                &
!                MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
!
!
!GA_Crossover_Probability             = real_buffer(1) 
!GA_Mutation_Probability              = real_buffer(2) 
!GA_save_elites_Probability           = real_buffer(3) 
!dt                                   = real_buffer(4)  
!GP_Tree_Probability                  = real_buffer(5)
!GP_Elitist_Probability               = real_buffer(6)  
!GP_Asexual_Reproduction_Probability  = real_buffer(7)
!GP_Crossover_Probability             = real_buffer(8)
!GP_Mutation_Probability              = real_buffer(9)


!write(6,'(A,1x,I6,1x,E15.7)') 'bc1: myid, GA_Crossover_Probability ', &
!                                    myid, GA_Crossover_Probability 
!write(6,'(A,1x,I6,1x,E15.7)') 'bc1: myid, GP_Mutation_Probability ', &
!                                    myid, GP_Mutation_Probability 




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

!-----------------------------------------------------------------
!
!if( myid == 0 )then
!    integer_buffer(1) = n_GA_Generations
!    integer_buffer(2) = n_GA_Individuals
!    integer_buffer(3) = n_time_steps
!    integer_buffer(4) = n_gp_individuals
!    integer_buffer(5) = n_gp_generations
!    integer_buffer(6) = n_node_functions
!    integer_buffer(7) = user_input_random_seed
!endif ! myid == 0 
!
!call MPI_BCAST( integer_buffer, 7,    &
!                MPI_INTEGER,  0, MPI_COMM_WORLD, ierr )
!
!n_GA_Generations       = integer_buffer(1) 
!n_GA_Individuals       = integer_buffer(2) 
!n_time_steps           = integer_buffer(3) 
!n_gp_individuals       = integer_buffer(4) 
!n_gp_generations       = integer_buffer(5) 
!n_node_functions       = integer_buffer(6) 
!user_input_random_seed = integer_buffer(7) 

!write(6,'(A,2(1x,I6))') 'bc1: myid, n_gp_individuals         ', &
!                              myid, n_gp_individuals         
!write(6,'(A,2(1x,I6))') 'bc1: myid, n_node_functions         ', &
!                              myid, n_node_functions         
!write(6,'(A,1x,I6,1x,I12)') 'bc1: myid, user_input_random_seed ', &
!                                  myid, user_input_random_seed 
!-----------------------------------------------------------------


call MPI_BCAST( L_GA_print, 1,    &
                MPI_LOGICAL,  0, MPI_COMM_WORLD, ierr )


call MPI_BARRIER( MPI_COMM_WORLD, ierr )    ! necessary? 

return

end subroutine bcast1
