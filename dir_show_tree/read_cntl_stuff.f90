subroutine read_cntl_stuff( )


!use mpi
!use mpi_module

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module
use GP_Data_module


IMPLICIT NONE


integer(kind=4) :: istat

integer(kind=4), parameter :: cntl_unitnum  = 501
integer(kind=4), parameter :: line_length   = 150

CHARACTER(line_length) :: Aline


real(kind=8) :: dt_min

!----------------------------------------------------------------------


! START OF EXECUTABLE CODE

! open the control input file

open( unit = cntl_unitnum, file = 'GPGACODE_cntl', form = 'formatted',&
      status = 'old' )


rewind(cntl_unitnum)


!---------------------------------------------------------------------

! echo control input


WRITE(6,'(//A)' ) &
'Input Echo Listing------------------------------------------------'

echoloop: &
do

    Aline(1:) = ' '

    istat  = 0
    READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline
    if( istat > 0 ) then
        WRITE(6,*) &
        'rcntl: ERROR *** Problem reading GPGACODE_cntl &
                             &in subroutine read_cntl_stuff'
        STOP 'bad read 1 cntl file'
   endif
   if( istat < 0 ) then
       EXIT echoloop
   endif

   WRITE(6,'(A)') trim( Aline )

enddo echoloop


WRITE(6,'(A//)' )&
 'End of Input Echo Listing-----------------------------------------'



!---------------------------------------------------------------------

! defaults
n_Node_Functions = 7
n_GP_individuals = 9
n_GP_generations = 1

GP_Tree_Probability=0.5d0

! Note: The next 4 parameters must add up to 1.0
GP_Elitist_Probability = 0.1d0
GP_Asexual_Reproduction_Probability = 0.4d0
GP_Crossover_Probability = 0.4d0
GP_Mutation_Probability = 0.1d0



!---------------------------------------------------------------------

rewind(cntl_unitnum)


cntlloop: &
do

    Aline = ' '

    istat  = 0
    READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline

    if( istat > 0 ) then
        WRITE(6,'(/A/)') &
         'rcntl: ERROR *** Problem reading GPGACODE_cntl.'
        STOP 'bad read 2 cntl file'
    endif
    if( istat < 0 ) then
        EXIT cntlloop
    endif


!------------------------------------------------------------------------------


!GA_Crossover_Probability = 0.3d0   ! probability of sexual crossing of parameter strings in GA_lmdif


    if( Aline(1:len('GA_Crossover_Probability')) == "GA_Crossover_Probability" .or.     &
        Aline(1:len('GA_Crossover_Probability')) == "ga_crossover_probability" ) then

        READ(Aline(len('GA_Crossover_Probability')+1:), * ) GA_Crossover_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GA_Crossover_Probability   = ', &
                                        GA_Crossover_Probability


!--------------------------------------------------------------------


!GA_Mutation_Probability  = 0.1d0   ! probability of mutation in parameter string of GA_lmdif

    elseif( Aline(1:len('GA_Mutation_Probability')) == "GA_Mutation_Probability" .or.     &
            Aline(1:len('GA_Mutation_Probability')) == "ga_mutation_probability" ) then

        READ(Aline(len('GA_Mutation_Probability')+1:), * ) GA_Mutation_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GA_Mutation_Probability    = ', &
                                        GA_Mutation_Probability


!--------------------------------------------------------------------


!GA_save_elites_Probability  = 0.005d0   ! probability of saving an individual as an elite individual

    elseif( Aline(1:len('GA_save_elites_Probability')) == "GA_save_elites_Probability" .or.     &
            Aline(1:len('GA_save_elites_Probability')) == "ga_save_elites_probability" ) then

        READ(Aline(len('GA_save_elites_Probability')+1:), * ) GA_save_elites_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GA_save_elites_Probability = ', &
                                        GA_save_elites_Probability


!--------------------------------------------------------------------


!GP_Tree_Probability  = 0.005d0   ! Estimated from previous work by Joel Cohen

    elseif( Aline(1:len('GP_Tree_Probability')) == "GP_Tree_Probability" .or.     &
            Aline(1:len('GP_Tree_Probability')) == "gp_tree_probability" ) then

        READ(Aline(len('GP_Tree_Probability')+1:), * ) GP_Tree_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GP_Tree_Probability = ', &
                                        GP_Tree_Probability


!--------------------------------------------------------------------


!GP_Elitist_Probability  = 0.005d0   ! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation

    elseif( Aline(1:len('GP_Elitist_Probability')) == "GP_Elitist_Probability" .or.     &
            Aline(1:len('GP_Elitist_Probability')) == "gp_elitist_probability" ) then

        READ(Aline(len('GP_Elitist_Probability')+1:), * ) GP_Elitist_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GP_Elitist_Probability = ', &
                                        GP_Elitist_Probability


!--------------------------------------------------------------------

!GP_Asexual_Reproduction_Probability  = 0.005d0   ! probability of asexual reproduction

    elseif( Aline(1:len('GP_Asexual_Reproduction_Probability')) == "GP_Asexual_Reproduction_Probability" .or.     &
            Aline(1:len('GP_Asexual_Reproduction_Probability')) == "gp_asexual_reproduction_probability" ) then

        READ(Aline(len('GP_Asexual_Reproduction_Probability')+1:), * ) GP_Asexual_Reproduction_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GP_Asexual_Reproduction_Probability = ', &
                                        GP_Asexual_Reproduction_Probability

!------------------------------------------------------------------------------

!GP_Crossover_Probability  = 0.005d0   !  probability of sexual crossing of binary string

    elseif( Aline(1:len('GP_Crossover_Probability')) == "GP_Crossover_Probability" .or.     &
            Aline(1:len('GP_Crossover_Probability')) == "gp_crossover_probability" ) then

        READ(Aline(len('GP_Crossover_Probability')+1:), * ) GP_Crossover_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GP_Crossover_Probability = ', &
                                        GP_Crossover_Probability



!--------------------------------------------------------------------

!GP_Mutation_Probability  = 0.005d0   ! probability of mutation in binary string

    elseif( Aline(1:len('GP_Mutation_Probability')) == "GP_Mutation_Probability" .or.     &
            Aline(1:len('GP_Mutation_Probability')) == "gp_mutation_probability" ) then

        READ(Aline(len('GP_Mutation_Probability')+1:), * ) GP_Mutation_Probability

        write(6,'(A,1x,F10.4)') 'rcntl: GP_Mutation_Probability = ', &
                                        GP_Mutation_Probability



!--------------------------------------------------------------------


!n_GA_Generations =  1000 !  5000 ! 1000

    elseif( Aline(1:len('n_GA_Generations')) == "n_GA_Generations" .or.     &
            Aline(1:len('n_GA_Generations')) == "n_ga_generations" ) then

        READ(Aline(len('n_GA_Generations')+1:), * ) n_GA_Generations

        write(6,'(A,1x,I6)') 'rcntl: n_GA_Generations = ', &
                                     n_GA_Generations


!--------------------------------------------------------------------


!n_GA_Individuals =  2000 ! 3  ! 100

    elseif( Aline(1:len('n_GA_Individuals')) == "n_GA_Individuals" .or.     &
            Aline(1:len('n_GA_Individuals')) == "n_GA_Individuals" ) then

        READ(Aline(len('n_GA_Individuals')+1:), * ) n_GA_Individuals

        write(6,'(A,1x,I6)') 'rcntl: n_GA_Individuals = ', &
                                     n_GA_Individuals


!--------------------------------------------------------------------

!n_time_steps= 2000 ! 8 ! 10     ! debug

    elseif( Aline(1:len('n_time_steps')) == "N_Time_Steps" .or.     &
            Aline(1:len('n_time_steps')) == "n_time_steps" ) then

        READ(Aline(len('n_time_steps')+1:), * ) n_time_steps

        write(6,'(A,1x,I6)') 'rcntl: n_time_steps     = ', &
                                     n_time_steps


!--------------------------------------------------------------------


!dt = 1.0D+1/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

    elseif( Aline(1:len('DT')) == "DT" .or.     &
            Aline(1:len('DT')) == "dt" ) then

        READ(Aline(len('DT')+1:), * )  dt_min

        write(6,'(A,1x,F10.4)') 'rcntl: dt (minutes) = ', dt_min

        dt = dt_min / 1440.0d0

        write(6,'(A,1x,F10.4)') 'rcntl: dt (days)    = ', dt



!--------------------------------------------------------------------


!dt = 1.0D+1/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

    elseif( Aline(1:len('model')) == "MODEL" .or.     &
            Aline(1:len('model')) == "model" ) then

        READ(Aline(len('model')+1:), * )  model

        write(6,'(A,1x,A)') 'rcntl: model = ', trim( model )



!--------------------------------------------------------------------


!N_GP_individuals

    elseif( Aline(1:len('n_gp_individuals')) == "N_GP_INDIVIDUALS" .or.     &
            Aline(1:len('n_gp_individuals')) == "n_gp_individuals" ) then

        READ(Aline(len('n_gp_individuals')+1:), * )  n_gp_individuals

        write(6,'(A,1x,I6)') 'rcntl: n_gp_individuals = ', n_gp_individuals



!--------------------------------------------------------------------


!N_GP_generations

    elseif( Aline(1:len('n_gp_generations')) == "N_GP_GENERATIONS" .or.     &
            Aline(1:len('n_gp_generations')) == "n_gp_generations" ) then

        READ(Aline(len('n_gp_generations')+1:), * )  n_gp_generations

        write(6,'(A,1x,I6)') 'rcntl: n_gp_generations = ', n_gp_generations




!--------------------------------------------------------------------


!n_Node_Functions

    elseif( Aline(1:len('n_Node_Functions')) == "n_Node_Functions" .or.     &
            Aline(1:len('n_Node_Functions')) == "n_node_functions" ) then

        READ(Aline(len('n_Node_Functions')+1:), * )  n_Node_Functions

        write(6,'(A,1x,I6)') 'rcntl: n_Node_Functions = ', n_Node_Functions






!--------------------------------------------------------------------

    else


    endif !   Aline(1:6) == ???


enddo cntlloop


return

!---------------------------------------------------------------------
END subroutine read_cntl_stuff
