Description of the unit 5 input control file                        20130724
                                                                    J. McCarthy



Some of the keywords can be input with some letters capitalized. These variants
are shown in parentheses in the keyword definition sections below.



!------------------------------------------------------------------------------

GA_Crossover_Probability  =   probability of sexual crossing of
                              parameter strings in GA_lmdif



col 1-24	GA_Crossover_Probability  ( ga_crossover_probability ) 
col 26-80	value of the GA crossover probability
Example:

GA_Crossover_Probability  0.2

Note:  default value = 0.4


!--------------------------------------------------------------------


GA_Mutation_Probability  = probability of mutation in a parameter
                           string of GA_lmdif



col 1-23	GA_Mutation_Probability  (ga_mutation_probability) 
col 25-80	value of the GA mutation probability
Example:

GA_Mutation_Probability  0.2

Note:  default value = 0.2



!--------------------------------------------------------------------


GA_Rand_Replace_Probability = probability of replacing all the parameters
                              in a non-elite individual with random 
                              numbers 


col 1-27	GA_Rand_Replace_Probability ( ga_rand_replace_probability ) 
col 29-80	value of the GA random replace probability
Example:

GA_Rand_Replace_Probability 0.2

Note:  default value = 0.01 

!--------------------------------------------------------------------


GA_save_elites_Probability  = probability of saving an individual
                              as an elite individual.  "Saving" means
                              that elite individuals will not be modified
                              by mutation, etc. in the current generation



col 1-26	GA_save_elites_Probability ( ga_save_elites_probability )
col 28-80	value of the GA save_elites probability

Example:

GA_save_elites_Probability  0.001

Note:  default value = 0.0


!--------------------------------------------------------------------


GP_Tree_Probability  = ??
                       ! Estimated from previous work by Joel Cohen



col 1-19	GP_Tree_Probability ( gp_tree_probability )
col 21-80	value of the GP tree probability

Example:

GP_Tree_Probability   0.001

Note:  default value = 0.5


!--------------------------------------------------------------------


GP_Elitist_Probability  -  Keeps the top  n_GP_Elitists  of the
                           Best Fit Individuals from
                           Generation to Generation


col 1-22	GP_Elitist_Probability  ( gp_elitist_probability )
col 24-80	value of the GP elitist probability

Example:

GP_Elitist_Probability  0.001

Note:  default value = 0.1


!--------------------------------------------------------------------

GP_Asexual_Reproduction_Probability  = probability of asexual reproduction



col 1-35	GP_Asexual_Reproduction_Probability
                        ( gp_asexual_reproduction_probability )
col 37-80	value of the GP elitist probability

Example:

GP_Asexual_Reproduction_Probability 0.001

Note:  default value = 0.4

!------------------------------------------------------------------------------

GP_Crossover_Probability  = probability of sexual crossing of binary string



col 1-24	GP_Crossover_Probability  ( gp_crossover_probability )
col 26-80	value of the GP crossover probability

Example:

GP_Crossover_Probability  0.001

Note:  default value = 0.4



!--------------------------------------------------------------------

GP_Mutation_Probability  =  probability of mutation in a GP string



col 1-23	GP_Mutation_Probability  ( gp_mutation_probability ) 
col 25-80	value of the GP mutation probability

Example:

GP_Mutation_Probability 0.001

Note:  default value = 0.1



!--------------------------------------------------------------------


n_GA_Generations  = number of generations of the GA_lmdif process



col 1-16	n_GA_Generations  ( n_ga_generations ) 
col 18-80	number of generations of the GA_lmdif process

Example:

n_GA_Generations  10

Note:  default value = 0


!--------------------------------------------------------------------


n_GA_Individuals =  number of individuals used in the GA_lmdif process


col 1-16	n_GA_Individuals ( n_ga_individuals )
col 18-80	number of individuals in the GA_lmdif process

Example:

n_GA_Individuals  100

Note:  default value = 0

!--------------------------------------------------------------------

N_Time_Steps = number of time steps of length "dt" used to integrate
               the equations with the Runge-Kutta process



col 1-12	N_Time_Steps  ( n_time_steps ) 
col 14-80	number of n_time_steps in the Runge-Kutta integration

Example:

n_time_steps  2500

Note:  default value = 2500



!--------------------------------------------------------------------


DT = length of the time step for the Runge-Kutta integration (unit = minute**-1)



col 1-2		DT ( dt ) 
col 4-80	length of the Runge-Kutta time step in MINUTES

Example:

dt  10.0

Note:	dt is converted internally to units of days**-1

	default value = 10 minutes


!--------------------------------------------------------------------


MODEL = name of the model used to generate the truth data
        Currently, only the values "LV"  or  "NPZ"  are allowed


col 1-5		MODEL ( model ) 
col 7-80	name of the model used to generate the truth data
                Allowed values:  "LV"  and "NPZ"

Example:

model  LV

Note:	default value = LV 


!--------------------------------------------------------------------


N_GP_INDIVIDUALS  = number of GP individuals used in the GP process


col 1-16	N_GP_INDIVIDUALS ( n_gp_individuals ) 
col 18-80	number of GP individuals used in the GP process

Example:

N_GP_individuals  100

Note:	default value = 1


!--------------------------------------------------------------------


N_GP_GENERATIONS  = number of GP generations used in the GP process


col 1-16	N_GP_GENERATIONS ( n_gp_generations ) 
col 18-80	number of GP generations used in the GP process

Example:

n_gp_generations  100

Note:	default value = 1



!--------------------------------------------------------------------


n_Node_Functions


col 1-16	n_Node_Functions ( n_node_functions )
col 18-80	number of node functions

Example:

n_Node_Functions  7

Note:	default value = 7


!--------------------------------------------------------------------

RANDOM_SCALE_SMALL = random_scale_small is the smaller of the two scales
                     used to scale the random number in subroutine random_real
                     ( see RANDOM_SCALE_FRACTION for a detailed explanation )


col 1-18	RANDOM_SCALE_SMALL ( random_scale_small )
col 20-80	size of the small scale  used in random_real

Example:

RANDOM_SCALE_SMALL 3.0

Note:	default value = 1.0  


!--------------------------------------------------------------------


RANDOM_SCALE_LARGE = random_scale_large is the larger of the two scales
                     used to scale the random number in subroutine random_real
                     ( see RANDOM_SCALE_FRACTION for a detailed explanation )



col 1-18	RANDOM_SCALE_LARGE  ( random_scale_large )
col 20-80	size of the large scale  used in random_real

Example:

RANDOM_SCALE_LARGE 10.0

Note:	default value = 50.0



!--------------------------------------------------------------------



RANDOM_SCALE_FRACTION = if a random number is less than the random scale fraction,
                        then the small scale is chosen to scale the random number.
                        Otherwise the large scale is used to scale the random number
                        in subroutine random_real



col 1-21	RANDOM_SCALE_FRACTION ( random_scale_fraction )
col 23-80	size of the large scale  used in random_real

Example:

RANDOM_SCALE_FRACTION 0.5

Note:	default value = 0.6


 The output of the subroutine random_real is a randomly chosen value for one of the
model parameters.  Since the intrinsic random_number function returns uniformly
distributed random numbers in the range 0 to 1,  the number output by random_number
must be scaled to provide a larger range of parameter values.

 The current code operates as follows:

    1)  random_number is called to output a random number, R1, in the [0., 1.] range
    2)  R1 is compared to the value of RANDOM_SCALE_FRACTION
    3)  if R1 is less than the value of RANDOM_SCALE_FRACTION, then a new call
        is made to random_number, producing R2, and the number output by
        random_real  is R2 * RANDOM_SCALE_SMALL
    4)  if R1 is greater than the value of RANDOM_SCALE_FRACTION, then a new call
        is made to random_number, producing R3, and the number output by
        random_real  is R3 * RANDOM_SCALE_LARGE

 The reason for the above procedure is to get a reasonable set of numbers which
cover the possible range of the parameters.  Some parameters have large nominal values,
and others have very small values (all greater than zero).  If a single scale is
used to scale the random numbers from "random_number", then, if the scale is large
enough to produce numbers greater than the largest nominal parameter value, then
very few outputs of random_real will have values on the order of the smallest nominal
parameter values. The method described above allows the user to get random numbers
scaled to the size of the largest parameter values, but also can provide a substantial
fraction of small values for those parameters with small nominal values.
For example, the current Lotka-Volterra model has 7 parameters, one with nominal
value 30.0, one with nominal value 2.0, and the other 5 with values less than 1.0.  If
a single scale is used, it has to be greater than 30, otherwise the nominal value of
parameter 1 will not be in the range of the random numbers.  If you use a single scale
of 50, chosen to be larger than the parameter 1 nominal value, then  the fraction of
values returned by random real with values 2 or less will be 0.04,  and with values
equal to 1 or less, 0.02.  This will produce a large number of individuals with almost
no chance of having values close to the values in the nominal vector.

!--------------------------------------------------------------------

GA_TOURNAMENT_STYLE  -  chooses the method used in the
                        GA_Tournament_Style_Sexual_Reproduction subroutine

                        = 0  - swap unmodified segments of parents
                        = 1  - swap segments of parents and randomly reset
                               node value at the segment boundaries
                        = 2  - swap segments of parents and reset node at
                               segment boundaries using the JM formula
                               involving the mean and std. dev


col 1-19	GA_TOURNAMENT_STYLE  ( GA_tournament_style or ga_tournament_style ) 
col 21-80	flag to select the method used in the
                GA_Tournament_Style_Sexual_Reproduction subroutine


Example:

GA_TOURNAMENT_STYLE  0

Note:	default value = 0


!--------------------------------------------------------------------

USER_INPUT_RANDOM_SEED  - user_input_random_seed is used for debugging
                          since it allows multiple runs to be made which
                          have the same set of random numbers


col 1-22        USER_INPUT_RANDOM_SEED ( user_input_random_seed  )
col 24-80	integer value to select the method used in the
                pseudo-random number generation

		user_input_random_seed = 0  -  use system clock value for random number seed
		user_input_random_seed > 0  -  use this value for random number seed


Example:

USER_INPUT_RANDOM_SEED  345538

Note:	default value = 0 -- i.e. system clock value is used



