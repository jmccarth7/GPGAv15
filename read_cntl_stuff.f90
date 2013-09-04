subroutine read_cntl_stuff( )


use mpi
use mpi_module

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


write(GP_print_unit,'(//A)' ) &
'Input Echo Listing------------------------------------------------'

echoloop: &
do

    Aline(1:) = ' '

    istat  = 0
    READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline
    if( istat > 0 ) then
        write(GP_print_unit,*) &
        'rcntl: ERROR *** Problem reading GPGACODE_cntl &
                             &in subroutine read_cntl_stuff'
        STOP 'bad read 1 cntl file'
   endif
   if( istat < 0 ) then
       EXIT echoloop
   endif

   write(GP_print_unit,'(A)') trim( Aline )

enddo echoloop


write(GP_print_unit,'(A//)' )&
 'End of Input Echo Listing-----------------------------------------'



!---------------------------------------------------------------------

! defaults

user_input_random_seed = 0

random_scale_large    = 50.0d0
random_scale_small    =  1.0d0
random_scale_fraction =  0.6d0

n_Node_Functions = 7
n_GP_individuals = 1  !  9
n_GP_generations = 1

GA_Crossover_Probability     = 0.4d0
GA_Mutation_Probability      = 0.2d0
GA_rand_replace_Probability  = 0.01d0
GA_save_elites_Probability   = 0.0d0

GP_Tree_Probability=0.5d0

! Note: The next 4 parameters must add up to 1.0
GP_Elitist_Probability              = 0.1d0
GP_Asexual_Reproduction_Probability = 0.4d0
GP_Crossover_Probability            = 0.4d0
GP_Mutation_Probability             = 0.1d0

ga_tournament_style = 0

n_time_steps = 2500

dt = 10.0d0 / 1440.0d0  ! 10 minutes

model = 'LV'


GA_print_flag = 0
L_GA_print = .FALSE.

GA_output_parameters_flag  = 0
L_GA_output_parameters = .FALSE.

GP_output_parameters_flag  = 0
L_GP_output_parameters = .FALSE.

fort333_output_flag  = 0
L_fort333_output = .FALSE.

fort444_output_flag  = 0
L_fort444_output = .FALSE.

GA_log_flag  = 0
L_GA_log = .FALSE.

GP_log_flag  = 0
L_GP_log = .FALSE.

unit50_output_flag  = 0
L_unit50_output = .FALSE.

print_equations_flag = 0 
L_print_equations = .FALSE.

!---------------------------------------------------------------------

rewind(cntl_unitnum)


cntlloop: &
do

    Aline = ' '

    istat  = 0
    READ( cntl_unitnum, '(A)', IOSTAT = istat ) Aline

    if( istat > 0 ) then
        write(GP_print_unit,'(/A/)') &
         'rcntl: ERROR *** Problem reading GPGACODE_cntl.'
        STOP 'bad read 2 cntl file'
    endif
    if( istat < 0 ) then
        EXIT cntlloop
    endif


!------------------------------------------------------------------------------


!GA_Crossover_Probability = 0.3d0   
! probability of sexual crossing of parameter strings in GA_lmdif


    if( Aline(1:len('GA_Crossover_Probability')) == "GA_Crossover_Probability" .or.     &
        Aline(1:len('GA_Crossover_Probability')) == "ga_crossover_probability" ) then

        READ(Aline(len('GA_Crossover_Probability')+1:), * ) GA_Crossover_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Crossover_Probability   = ', &
                                                    GA_Crossover_Probability


!--------------------------------------------------------------------


!GA_Mutation_Probability  = 0.1d0   
! probability of mutation in parameter string of GA_lmdif

    elseif( Aline(1:len('GA_Mutation_Probability')) == "GA_Mutation_Probability" .or.     &
            Aline(1:len('GA_Mutation_Probability')) == "ga_mutation_probability" ) then

        READ(Aline(len('GA_Mutation_Probability')+1:), * ) GA_Mutation_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_Mutation_Probability    = ', &
                                                    GA_Mutation_Probability



!--------------------------------------------------------------------

!GA_rand_replace_Probability  = 0.005d0   ! probability of rand_replace in binary string

    elseif( Aline(1:len('GA_rand_replace_Probability')) == "GA_Rand_Replace_Probability" .or. &
            Aline(1:len('GA_rand_replace_Probability')) == "ga_rand_replace_probability" ) then

        READ(Aline(len('GA_rand_replace_Probability')+1:), * ) GA_rand_replace_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_rand_replace_Probability = ', &
                                                    GA_rand_replace_Probability



!--------------------------------------------------------------------


!GA_save_elites_Probability  = 0.005d0   
! probability of saving an individual as an elite individual

    elseif( Aline(1:len('GA_save_elites_Probability')) == "GA_save_elites_Probability" .or.     &
            Aline(1:len('GA_save_elites_Probability')) == "ga_save_elites_probability" ) then

        READ(Aline(len('GA_save_elites_Probability')+1:), * ) GA_save_elites_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GA_save_elites_Probability = ', &
                                                    GA_save_elites_Probability


!--------------------------------------------------------------------


!GP_Tree_Probability  = 0.005d0   ! Estimated from previous work by Joel Cohen

    elseif( Aline(1:len('GP_Tree_Probability')) == "GP_Tree_Probability" .or.     &
            Aline(1:len('GP_Tree_Probability')) == "gp_tree_probability" ) then

        READ(Aline(len('GP_Tree_Probability')+1:), * ) GP_Tree_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Tree_Probability = ', &
                                                    GP_Tree_Probability


!--------------------------------------------------------------------


!GP_Elitist_Probability  = 0.005d0
! Keeps the top n_GP_Elitists of the Best Fit Individuals from Generation to Generation

    elseif( Aline(1:len('GP_Elitist_Probability')) == "GP_Elitist_Probability" .or.     &
            Aline(1:len('GP_Elitist_Probability')) == "gp_elitist_probability" ) then

        READ(Aline(len('GP_Elitist_Probability')+1:), * ) GP_Elitist_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Elitist_Probability = ', &
                                                    GP_Elitist_Probability


!--------------------------------------------------------------------

!GP_Asexual_Reproduction_Probability  = 0.005d0   ! probability of asexual reproduction

    elseif( Aline(1:len('GP_Asexual_Reproduction_Probability')) ==               &
                                  "GP_Asexual_Reproduction_Probability" .or.     &
            Aline(1:len('GP_Asexual_Reproduction_Probability')) ==               &
                                  "gp_asexual_reproduction_probability" ) then

        READ(Aline(len('GP_Asexual_Reproduction_Probability')+1:), * ) &
                        GP_Asexual_Reproduction_Probability

        write(GP_print_unit,'(A,1x,F10.4)') &
              'rcntl: GP_Asexual_Reproduction_Probability = ', &
                      GP_Asexual_Reproduction_Probability

!------------------------------------------------------------------------------

!GP_Crossover_Probability  = 0.005d0   !  probability of sexual crossing of binary string

    elseif( Aline(1:len('GP_Crossover_Probability')) == "GP_Crossover_Probability" .or.     &
            Aline(1:len('GP_Crossover_Probability')) == "gp_crossover_probability" ) then

        READ(Aline(len('GP_Crossover_Probability')+1:), * ) GP_Crossover_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Crossover_Probability = ', &
                                                    GP_Crossover_Probability



!--------------------------------------------------------------------

!GP_Mutation_Probability  = 0.005d0   ! probability of mutation in binary string

    elseif( Aline(1:len('GP_Mutation_Probability')) == "GP_Mutation_Probability" .or.     &
            Aline(1:len('GP_Mutation_Probability')) == "gp_mutation_probability" ) then

        READ(Aline(len('GP_Mutation_Probability')+1:), * ) GP_Mutation_Probability

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: GP_Mutation_Probability = ', &
                                                    GP_Mutation_Probability



!--------------------------------------------------------------------

!n_GA_Generations 

    elseif( Aline(1:len('n_GA_Generations')) == "n_GA_Generations" .or.     &
            Aline(1:len('n_GA_Generations')) == "n_ga_generations" ) then

        READ(Aline(len('n_GA_Generations')+1:), * ) n_GA_Generations

        write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Generations = ', &
                                                 n_GA_Generations


!--------------------------------------------------------------------


!n_GA_Individuals 

    elseif( Aline(1:len('n_GA_Individuals')) == "n_GA_Individuals" .or.     &
            Aline(1:len('n_GA_Individuals')) == "n_ga_individuals" ) then

        READ(Aline(len('n_GA_Individuals')+1:), * ) n_GA_Individuals

        write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_GA_Individuals = ', &
                                                 n_GA_Individuals


!--------------------------------------------------------------------

!n_time_steps

    elseif( Aline(1:len('n_time_steps')) == "N_Time_Steps" .or.     &
            Aline(1:len('n_time_steps')) == "n_time_steps" ) then

        READ(Aline(len('n_time_steps')+1:), * ) n_time_steps

        write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_time_steps     = ', &
                                                 n_time_steps


!--------------------------------------------------------------------


!dt = 1.0D+1/(24.0D+0*60.0D+0)   ! [d^-1; 10 minute time step]

    elseif( Aline(1:len('DT')) == "DT" .or.     &
            Aline(1:len('DT')) == "dt" ) then

        READ(Aline(len('DT')+1:), * )  dt_min

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (minutes) = ', dt_min

        dt = dt_min / 1440.0d0

        write(GP_print_unit,'(A,1x,F10.4)') 'rcntl: dt (days)    = ', dt



!--------------------------------------------------------------------


!model = LV  or  NPZ

    elseif( Aline(1:len('model')) == "MODEL" .or.     &
            Aline(1:len('model')) == "model" ) then

        READ(Aline(len('model')+1:), * )  model

        write(GP_print_unit,'(A,1x,A)') 'rcntl: model = ', trim( model )



!--------------------------------------------------------------------


!N_GP_individuals

    elseif( Aline(1:len('n_gp_individuals')) == "N_GP_INDIVIDUALS" .or.     &
            Aline(1:len('n_gp_individuals')) == "n_gp_individuals" ) then

        READ(Aline(len('n_gp_individuals')+1:), * )  n_gp_individuals

        write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_individuals = ', n_gp_individuals



!--------------------------------------------------------------------


!N_GP_generations

    elseif( Aline(1:len('n_gp_generations')) == "N_GP_GENERATIONS" .or.     &
            Aline(1:len('n_gp_generations')) == "n_gp_generations" ) then

        READ(Aline(len('n_gp_generations')+1:), * )  n_gp_generations

        write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_gp_generations = ', n_gp_generations




!--------------------------------------------------------------------


!n_Node_Functions

    elseif( Aline(1:len('n_Node_Functions')) == "n_Node_Functions" .or.     &
            Aline(1:len('n_Node_Functions')) == "n_node_functions" ) then

        READ(Aline(len('n_Node_Functions')+1:), * )  n_Node_Functions

        write(GP_print_unit,'(A,1x,I6)') 'rcntl: n_Node_Functions = ', n_Node_Functions


!--------------------------------------------------------------------

! random_scale_small 

! in random_real, random_scale_small is the smaller of the two scales
! used to scale the random number

    elseif( Aline(1:len('random_scale_small')) == "RANDOM_SCALE_SMALL" .or.  &
            Aline(1:len('random_scale_small')) == "random_scale_small" ) then

        READ(Aline(len('random_scale_small')+1:), * )  random_scale_small

        write(GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_small = ', &
                                                    random_scale_small


!--------------------------------------------------------------------

! random_scale_large 

! in random_real, random_scale_large is the larger of the two scales
! used to scale the random number

    elseif( Aline(1:len('random_scale_large')) == "RANDOM_SCALE_LARGE" .or.  &
            Aline(1:len('random_scale_large')) == "random_scale_large" ) then

        READ(Aline(len('random_scale_large')+1:), * )  random_scale_large

        write(GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_large = ', &
                                                    random_scale_large




!--------------------------------------------------------------------


! random scale fraction

! if a random number is less than the random scale fraction, then
! the small scale is chosen to scale the random number in random_real

    elseif( Aline(1:len('random_scale_fraction')) == &
                        "RANDOM_SCALE_FRACTION"        .or.     &
            Aline(1:len('random_scale_fraction')) == &
                        "random_scale_fraction"           ) then

        READ(Aline(len('random_scale_fraction')+1:), * )  random_scale_fraction

        write(GP_print_unit,'(A,1x,E15.7)') 'rcntl: random_scale_fraction = ', &
                                                    random_scale_fraction




!--------------------------------------------------------------------


! ga_tournament_style     

! = 0  - swap unmodified segments of parents 
! = 1  - swap segments of parents and randomly reset node at segment boundaries
! = 2  - swap segments of parents and reset node at segment boundaries using JM 
!        formula involving the mean and std. dev


    elseif( Aline(1:len('ga_tournament_style')) == "ga_tournament_style" .or. &
            Aline(1:len('ga_tournament_style')) == "GA_tournament_style" .or.     &      
            Aline(1:len('ga_tournament_style')) == "GA_TOURNAMENT_STYLE"           ) then

        READ(Aline(len('ga_tournament_style')+1:), * )  ga_tournament_style

        write(GP_print_unit,'(A,1x,I6)') 'rcntl: ga_tournament_style = ', &
                                                 ga_tournament_style


!--------------------------------------------------------------------


!  user_input_random_seed

! user_input_random_seed = 0  -  use system clock value for random number seed

! user_input_random_seed > 0  -  use this value for random number seed

! user_input_random_seed is used for debugging since it allows multiple
! runs to be made which have the same set of random numbers



    elseif( Aline(1:len('user_input_random_seed')) == "user_input_random_seed"  .or.     &
            Aline(1:len('user_input_random_seed')) == "USER_INPUT_RANDOM_SEED"           ) then

        READ(Aline(len('user_input_random_seed')+1:), * )  user_input_random_seed


        user_input_random_seed = abs( user_input_random_seed  ) 

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: user_input_random_seed =', &
                                                  user_input_random_seed


!--------------------------------------------------------------------

!  GA_print   

! if GA_print_flag >  0 - write printout to GA_print_unit
! if GA_print_flag <= 0 - do not write printout to GA_print_unit

! DEFAULT =   GA_print_flag =  0 -  do not write printout to GA_print_unit



    elseif( Aline(1:len('GA_print')) == "GA_print"  .or.     &
            Aline(1:len('GA_print')) == "ga_print"           ) then

        READ(Aline(len('GA_print')+1:), * )  GA_print_flag


        if( GA_print_flag > 0 )then
            L_GA_print = .TRUE.
        else
            L_GA_print = .FALSE. 
        endif ! GA_print_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: GA_print_flag =', &
                                                  GA_print_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_print =', &
                                                  L_GA_print


!--------------------------------------------------------------------

! GA output parameters  - formerly the file was called "output_parameters"


! if GA_output_parameters_flag >  0 - write printout to GA_output_parameters_unit
! if GA_output_parameters_flag <= 0 - do not write printout to GA_output_parameters_unit

!  DEFAULT =   GA_output_parameters_flag == 0 
!              - do not write printout to GA_output_parameters_unit



    elseif( Aline(1:len('GA_output_parameters')) == "GA_output_parameters"  .or.     &
            Aline(1:len('GA_output_parameters')) == "ga_output_parameters"           ) then

        READ(Aline(len('GA_output_parameters')+1:), * )  GA_output_parameters_flag


        if( GA_output_parameters_flag > 0 )then
            L_GA_output_parameters = .TRUE.
        else
            L_GA_output_parameters = .FALSE. 
        endif ! GA_output_parameters_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: GA_output_parameters_flag =', &
                                                  GA_output_parameters_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_output_parameters =', &
                                                  L_GA_output_parameters

!--------------------------------------------------------------------

! GP_output_parameters

! if GP_output_parameters_flag >  0 - write printout to GP_output_parameters_unit
! if GP_output_parameters_flag <= 0 - do not write printout to GP_output_parameters_unit

!  DEFAULT =   GP_output_parameters_flag == 0 
!              - do not write printout to GA_print_unit



    elseif( Aline(1:len('GP_output_parameters')) == "GP_output_parameters"  .or.     &
            Aline(1:len('GP_output_parameters')) == "gp_output_parameters"           ) then

        READ(Aline(len('GP_output_parameters')+1:), * )  GP_output_parameters_flag


        if( GP_output_parameters_flag > 0 )then
            L_GP_output_parameters = .TRUE.
        else
            L_GP_output_parameters = .FALSE. 
        endif ! GP_output_parameters_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: GP_output_parameters_flag =', &
                                                  GP_output_parameters_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_output_parameters =', &
                                                  L_GP_output_parameters

!--------------------------------------------------------------------

! fort333_output   

! if fort333_output_flag >  0 - write printout to fort333_output_unit
! if fort333_output_flag <= 0 - do not write printout to fort333_output_unit

!  DEFAULT =   fort333_output_flag == 0 
!              - do not write printout to fort333_output_unit



    elseif( Aline(1:len('fort333_output')) == "fort333_output"  ) then


        READ(Aline(len('fort333_output')+1:), * )  fort333_output_flag

        if( fort333_output_flag > 0 )then
            L_fort333_output = .TRUE.
        else
            L_fort333_output = .FALSE. 
        endif ! fort333_output_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: fort333_output_flag =', &
                                                  fort333_output_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort333_output =', &
                                                  L_fort333_output


!--------------------------------------------------------------------

! fort444_output   

! if fort444_output_flag >  0 - write printout to fort444_output_unit
! if fort444_output_flag <= 0 - do not write printout to fort444_output_unit

!  DEFAULT =   fort444_output_flag == 0 
!              - do not write printout to fort444_output_unit



    elseif( Aline(1:len('fort444_output')) == "fort444_output"  ) then


        READ(Aline(len('fort444_output')+1:), * )  fort444_output_flag

        if( fort444_output_flag > 0 )then
            L_fort444_output = .TRUE.
        else
            L_fort444_output = .FALSE. 
        endif ! fort444_output_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: fort444_output_flag =', &
                                                  fort444_output_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_fort444_output =', &
                                                  L_fort444_output


!--------------------------------------------------------------------

! GA log

!  if GA_log_flag >  0 - write printout to GA_log_unit
!  if GA_log_flag <= 0 - do not write printout to GA_log_unit

!  DEFAULT =   GA_log_flag == 0 
!              - do not write printout to GA_log_unit



    elseif( Aline(1:len('GA_log')) == "GA_log"  .or.     &
            Aline(1:len('GA_log')) == "ga_log"           ) then


        READ(Aline(len('GA_log')+1:), * )  GA_log_flag

        if( GA_log_flag > 0 )then
            L_GA_log = .TRUE.
        else
            L_GA_log = .FALSE. 
        endif ! GA_log_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: GA_log_flag =', &
                                                  GA_log_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GA_log =', &
                                                  L_GA_log


!--------------------------------------------------------------------

! GP log

!  if GP_log_flag >  0 - write printout to GP_log_unit
!  if GP_log_flag <= 0 - do not write printout to GP_log_unit

!  DEFAULT =   GP_log_flag == 0 
!             - do not write printout to GP_log_unit



    elseif( Aline(1:len('GP_log')) == "GP_log"  .or.     &
            Aline(1:len('GP_log')) == "gp_log"           ) then


        READ(Aline(len('GP_log')+1:), * )  GP_log_flag

        if( GP_log_flag > 0 )then
            L_GP_log = .TRUE.
        else
            L_GP_log = .FALSE. 
        endif ! GP_log_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: GP_log_flag =', &
                                                  GP_log_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_GP_log =', &
                                                  L_GP_log

!--------------------------------------------------------------------


! unit50_output   


! if unit50_output_flag >  0 - write printout to unit50_output_unit
! if unit50_output_flag <= 0 - do not write printout to unit50_output_unit

!  DEFAULT =   unit50_output_flag ==  0 
!              - do not write printout to unit50_output_unit



    elseif( Aline(1:len('unit50_output')) == "unit50_output" ) then


        READ(Aline(len('unit50_output')+1:), * )  unit50_output_flag

        if( unit50_output_flag > 0 )then
            L_unit50_output = .TRUE.
        else
            L_unit50_output = .FALSE. 
        endif ! unit50_output_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: unit50_output_flag =', &
                                                  unit50_output_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_unit50_output =', &
                                                  L_unit50_output



!--------------------------------------------------------------------

! print equations 

! if print_equations_flag >  0 - write equations together with tree 
!                                structures in subroutine print_trees
! if print_equations_flag <= 0 - do not write equations

!  DEFAULT =   print_equations_flag ==  0 
!                - do not write quations



    elseif( Aline(1:len('print_equations')) == "print_equations" ) then


        READ(Aline(len('print_equations')+1:), * )  print_equations_flag

        if( print_equations_flag > 0 )then
            L_print_equations = .TRUE.
        else
            L_print_equations = .FALSE. 
        endif ! print_equations_flag > 0

        write(GP_print_unit,'(A,1x,I12)') 'rcntl: print_equations_flag =', &
                                                  print_equations_flag
        write(GP_print_unit,'(A,4x,L1 )') 'rcntl: L_print_equations =', &
                                                  L_print_equations


!--------------------------------------------------------------------

    else

        write(GP_print_unit,'(/A)') 'rcntl: WARNING: UNRECOGNIZED OPTION '

        write(GP_print_unit,'(A,1x,A)') 'rcntl: Aline =', trim( Aline )           
        write(GP_print_unit,'(A/)') 'rcntl: WARNING: UNRECOGNIZED OPTION '

    endif !   Aline(1:6) == ???


enddo cntlloop


return

!---------------------------------------------------------------------
END subroutine read_cntl_stuff
