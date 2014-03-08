subroutine init_values_fasham( icall  )


!     written by John R. Moisan [14 November 2012]
!     for GPCODE testing/developing

! Lotka_Volterra_Example_Set_Up
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This is the tree representation of the CODE System
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! this examples is a simple Lotka-Volterra model.

! dP/dt = (grow * P)  - (graze * P * Z)
! dZ/dt = (graze * P * Z) - ((1 - efficiency) * graze * P *  Z) - (amort * Z)
! [Note: In this example, the (1-efficiency) parameter is set to a new unique parameter, 'effic']

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! 1: [ 3, 0,-1, N, N, N, N, N, N, N, N, N, N, N, N]
! 2: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 3: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 4: [ 3, 0, 3, N, N,-1,-2, N, N, N, N, N, N, N, N]
! 5: [ 1, 3, 3, 0,-2, 3, 3, N, N, N, N, 0,-2, 0,-1]
! 6: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!    [01,02,03,04,05,06,07,08,09,10,11,12,13,14,15]

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use mpi
use mpi_module

use GP_parameters_module
use GP_variables_module
use fasham_variables_module

implicit none


integer,intent(in)  :: icall

!logical :: LV_model1 = .TRUE.

integer(kind=4) :: i_Tree
integer(kind=4) :: i_Node

integer(kind=4) :: i
real(kind=8) :: increment

!-------------------------------------------------------------------------

if(  icall  == 0  )then




    n_CODE_equations =   7

    n_variables = 7

    !!! n_code_forcing = 4   ! parameter

    !!! n_Tracked_Resources = 1 ! parameter  ! number of different resources being tracked


    n_trees=  ((n_CODE_equations+1)**2)-(n_CODE_equations+1)

    n_nodes = pow2_table( n_levels )  ! n_nodes = int(2**n_levels)-1


    !orig n_maximum_number_parameters = n_CODE_equations +  n_nodes
    !n_maximum_number_parameters       = n_CODE_equations *  n_nodes

    !n_maximum_number_parameters       = n_CODE_equations +  n_nodes
    n_maximum_number_parameters       =  n_trees * 2    +  n_nodes   ! jjm 20140307


    if( myid == 0 )then
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_levels           ', n_levels
        write(GP_print_unit,'(A,2(1x,I6))')&
              'ivFA: int(2**n_levels)-1 , pow2_table( n_levels )', &
                     int(2**n_levels)-1 , pow2_table( n_levels )
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_CODE_equations   ', n_CODE_equations
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_input_vars       ', n_input_vars
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_trees            ', n_trees
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_variables        ', n_variables
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_code_forcing     ', n_code_forcing
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_Tracked_Resources', n_Tracked_Resources
        write(GP_print_unit,'(A,1x,I6)') 'ivFA: n_nodes            ', n_nodes
        write(GP_print_unit,'(A,1x,I6/)')'ivFA: n_maximum_number_parameters  ', &
                                                n_maximum_number_parameters
    endif ! myid == 0

    return

endif ! icall == 0




! load the arrays



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     set the GPCODE tree control values
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

do  i_tree = 1,n_trees

    do  i_node = 1,n_nodes
        GP_Individual_Node_Parameters(i_node,i_tree) = 0.0d0
        tree_evaluation(i_node,i_tree) = 0.0d0
        GP_Individual_Node_Type(i_node,i_tree)       = -9999
    enddo ! i_node

enddo ! i_tree

!--------------------------------------------------------------------------------------

    ! Set Variables

    alpha   =  0.025D+0    ! initial slope of the P-I curve [(W m-2)-1 d-1]
    aK1     =  0.5D+0      ! half-saturation for phytoplankton NO3 uptake [mMol N m-3]
    aK2     =  0.5D+0      ! half-saturation for phytoplankton NH4 uptake [mMol N m-3]
    amu1    =  0.045D+0    ! phytoplankton specific mortality rate [d-1]
    akc     =  0.03D+0     ! light attenuation by phytoplankton [m^2 mMol N)-1]
    gamma1  =  0.05D+0     ! fraction of total primary production that is exuded [n.d.]
    phi     =  1.5D+0      ! phytoplankton ammonium inhibition parameter [(mMol N)-1]
    g       =  1.0D+0      ! maximum zooplankton growth rate [d-1]
    beta1   =  0.75D+0     ! zooplankton assimilation efficiency of zooplankton [n.d.]
    beta2   =  0.75D+0     ! zooplankton assimilation efficiency of phytoplankton [n.d.]
    beta3   =  0.75D+0     ! zooplankton assimilation efficiency of bacteria [n.d.]
    amu2    =  0.1D+0      ! zooplankton specific excretion rate [d-1]
    amu5    =  0.05D+0     ! zooplankton specific mortality rate [d-1]
    aK3     =  1.0D+0      ! zooplankton half-saturation conts. for ingestion [d-1]
    omega   =  0.33D+0     ! detrital fraction of zooplankton mortality [n.d.]
    epsilon =  0.75D+0     ! ammonium fraction of zooplankton excretion [n.d.]
    Vb      =  2.0D+0      ! bacteria maximum growth rate [d-1]
    Vp      =  2.9D+0      ! phyto maximum growth rate [d-1]
    amu3    =  0.05D+0     ! bacteria specific excretion rate [d-1]
    aK4     =  0.5D+0      ! bacteria half-saturation rate for uptake [(mMol N) m-3]
    eta     =  0.6D+0      ! ammonium/DON uptake ratio [n.d.]
    amu4    =  0.05D+0     ! detrital breakdown rate [d-1]
    V       =  1.0D+0      ! detrital sinking rate [m d-1]
    p1      =  1.0D+0      ! zooplankton preference for phytoplankton [n.d.]
    p2      =  1.0D+0      ! zooplankton preference for bacteria [n.d.]
    p3      =  1.0D+0      ! zooplankton preference for detritus [n.d.]
    aN0     =  2.0D+0      ! concentration of NO3 below the mixed-layer [(mMol N) m-3]


    !-----------------------------------------------------------------------------------------

    ! parameters as in Table 1; Fasham et al. [JMR, 48, 591-639, 1990]

    akw = 0.04D+0     ! light attenuation due to sea water [m-1]
    am  = 0.1D+0      ! cross-thermocline mixing rate [m d-1]

    alatd=50.0 !Latitude

    !-----------------------------------------------------------------------------------------


    ! initialize the biological data fields

    aNO3    =   0.2D+0 ! Nitrate           [mmol N m-3]
    aNH4    =   0.1D+0 ! Ammonium          [mmol N m-3]
    DON     =   0.1D+0 ! DON               [mmol N m-3]
    DET     =   0.1D+0 ! DET [Detritus]    [mmol N m-3]
    bact    =   0.1D+0 ! Bacteria          [mmol N m-3]
    phyto   =   0.1D+0 ! Phytoplankton     [mmol N m-3]
    zoo     =   0.1D+0 ! Zooplankton       [mmol N m-3]

    !-----------------------------------------------------------------------------------------

    ! Enumerations that represent model variables.
    ! These are used by the binary tree parsing algorithm to select the index of the
    ! species or forcing function variable's value

    SPECIES_NITRATE = -1
    SPECIES_AMMONIUM = -2
    SPECIES_DISSOLVED_ORGANIC_NITROGEN = -3
    SPECIES_DETRITUS = -4
    SPECIES_BACTERIA = -5
    SPECIES_PHYTOPLANKTON = -6
    SPECIES_ZOOPLANKTON = -7


!---------------------------------------------------------------------------------------------

    ! made parameters and init in fasham*mod*

    ! See comment in GP_Variables

    bioflo_map(:,1) = (/ SPECIES_NITRATE, &
                         SPECIES_AMMONIUM, &
                         SPECIES_DISSOLVED_ORGANIC_NITROGEN, &
                         SPECIES_DETRITUS, &
                         SPECIES_BACTERIA, &
                         SPECIES_PHYTOPLANKTON, &
                         SPECIES_ZOOPLANKTON /)

    ! Since indexes are all negative, take the absolute value
    bioflo_map = abs(bioflo_map)

    !-------------------------------------------------------------------------------

    ! initialized as parameters in fasham_variables_module

    !FORCING_MIXED_LAYER_DEPTH         = -5001
    !FORCING_MLD_CHANGE_NON_MOTILE     = -5002
    !FORCING_MLD_CHANGE_MOTILE         = -5003
    !FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004



    Numerical_CODE_Initial_Conditions = (/aNO3, aNH4, DON, DET, bact, phyto, zoo/)
    Numerical_CODE_Forcing_Functions = 0.0D+0



!--------------------------------------------------------------------------------------

! Initial Conditions


if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,I6)')   'ivFA: n_CODE_equations  ', n_CODE_equations
    write(GP_print_unit,'(A,1x,E15.7/)') &
          'ivFA: Numerical_CODE_Initial_Conditions(1:n_code_equations)', &
                 Numerical_CODE_Initial_Conditions(1:n_code_equations)
endif ! myid == 0



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

if( n_levels == 6 )then
!   n_levels = 6
    Node_Probability = (/0.8d0,0.7d0,6.d0, &
                         0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]


elseif( n_levels == 7 )then
    !!  n_levels = 7
    Node_Probability = (/0.8d0,0.7d0,6.d0, &
                         0.5d0,0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]

elseif( n_levels == 8 )then
    !   n_levels = 8
    Node_Probability = (/0.9d0,0.8d0,0.7d0,6.d0, &
                         0.5d0,0.4d0,0.3d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]
else

    increment = 1.0d0 / real( n_levels, kind=8 )

    do  i = 1, n_levels-1
        Node_Probability(i) = 1.0d0 - increment * real(i,kind=8)
    enddo
    Node_Probability(n_levels) = 0.0d0

endif ! n_levels == 6

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,I6)')   'ivFA: n_levels ', n_levels
    write(GP_print_unit,'(A/(10(1x,E12.5)))') 'ivFA: Node_Probability', &
                                                     Node_Probability
    write(GP_print_unit,'(A)') ' '
endif ! myid == 0


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




return

END subroutine init_values_fasham
