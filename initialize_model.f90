subroutine Initialize_Model( buildTrees, L_myprint, myprint_unit )


use mpi
use mpi_module

use fasham_variables_module


use GP_parameters_module
use GA_parameters_module
use GP_variables_module

implicit none

logical :: buildTrees

integer :: i

logical, intent(in)  ::  L_myprint
integer, intent(in)  ::  myprint_unit


!------------------------------------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(/A,6x,L1/)') 'inmod: entry  Initialize_Model buildtrees = ', buildtrees
endif ! myid == 0

!if( L_myprint  .and.  myid == 1 )then
!    write(myprint_unit,'(/A/)') 'inmod: entry  Initialize_Model '
!endif ! myid == 1


!-------------------------------------------------------------------------------
!!if(  model == 'fasham'  .or. &
!!     model == 'FASHAM'         )then


!!    !-----------------------------------------------------------------------------------------
!!
!!    ! Set Variables
!!
!!    alpha   =  0.025D+0    ! initial slope of the P-I curve [(W m-2)-1 d-1]
!!    aK1     =  0.5D+0      ! half-saturation for phytoplankton NO3 uptake [mMol N m-3]
!!    aK2     =  0.5D+0      ! half-saturation for phytoplankton NH4 uptake [mMol N m-3]
!!    amu1    =  0.045D+0    ! phytoplankton specific mortality rate [d-1]
!!    akc     =  0.03D+0     ! light attenuation by phytoplankton [m^2 mMol N)-1]
!!    gamma1  =  0.05D+0     ! fraction of total primary production that is exuded [n.d.]
!!    phi     =  1.5D+0      ! phytoplankton ammonium inhibition parameter [(mMol N)-1]
!!    g       =  1.0D+0      ! maximum zooplankton growth rate [d-1]
!!    beta1   =  0.75D+0     ! zooplankton assimilation efficiency of zooplankton [n.d.]
!!    beta2   =  0.75D+0     ! zooplankton assimilation efficiency of phytoplankton [n.d.]
!!    beta3   =  0.75D+0     ! zooplankton assimilation efficiency of bacteria [n.d.]
!!    amu2    =  0.1D+0      ! zooplankton specific excretion rate [d-1]
!!    amu5    =  0.05D+0     ! zooplankton specific mortality rate [d-1]
!!    aK3     =  1.0D+0      ! zooplankton half-saturation conts. for ingestion [d-1]
!!    omega   =  0.33D+0     ! detrital fraction of zooplankton mortality [n.d.]
!!    epsilon =  0.75D+0     ! ammonium fraction of zooplankton excretion [n.d.]
!!    Vb      =  2.0D+0      ! bacteria maximum growth rate [d-1]
!!    Vp      =  2.9D+0      ! phyto maximum growth rate [d-1]
!!    amu3    =  0.05D+0     ! bacteria specific excretion rate [d-1]
!!    aK4     =  0.5D+0      ! bacteria half-saturation rate for uptake [(mMol N) m-3]
!!    eta     =  0.6D+0      ! ammonium/DON uptake ratio [n.d.]
!!    amu4    =  0.05D+0     ! detrital breakdown rate [d-1]
!!    V       =  1.0D+0      ! detrital sinking rate [m d-1]
!!    p1      =  1.0D+0      ! zooplankton preference for phytoplankton [n.d.]
!!    p2      =  1.0D+0      ! zooplankton preference for bacteria [n.d.]
!!    p3      =  1.0D+0      ! zooplankton preference for detritus [n.d.]
!!    aN0     =  2.0D+0      ! concentration of NO3 below the mixed-layer [(mMol N) m-3]
!!
!!
!!    !-----------------------------------------------------------------------------------------
!!
!!    ! parameters as in Table 1; Fasham et al. [JMR, 48, 591-639, 1990]
!!    akw = 0.04D+0     ! light attenuation due to sea water [m-1]
!!    am  = 0.1D+0      ! cross-thermocline mixing rate [m d-1]
!!
!!    alatd=50.0 !Latitude
!!
!!    !-----------------------------------------------------------------------------------------
!!
!!    ! initialize the biological data fields
!!    aNO3    =   0.2D+0 ! Nitrate           [mmol N m-3]
!!    aNH4    =   0.1D+0 ! Ammonium          [mmol N m-3]
!!    DON     =   0.1D+0 ! DON               [mmol N m-3]
!!    DET     =   0.1D+0 ! DET [Detritus]    [mmol N m-3]
!!    bact    =   0.1D+0 ! Bacteria          [mmol N m-3]
!!    phyto   =   0.1D+0 ! Phytoplankton     [mmol N m-3]
!!    zoo     =   0.1D+0 ! Zooplankton       [mmol N m-3]
!!
!!    !-----------------------------------------------------------------------------------------
!!
!!    ! Enumerations that represent model variables. 
!!    ! These are used by the binary tree parsing algorithm to select the index of the
!!    ! species or forcing function variable's value
!!
!!    SPECIES_NITRATE = -1
!!    SPECIES_AMMONIUM = -2
!!    SPECIES_DISSOLVED_ORGANIC_NITROGEN = -3
!!    SPECIES_DETRITUS = -4
!!    SPECIES_BACTERIA = -5
!!    SPECIES_PHYTOPLANKTON = -6
!!    SPECIES_ZOOPLANKTON = -7
!!
!!
!!!---------------------------------------------------------------------------------------------
!!
!!    ! made parameters and init in fasham*mod*
!!
!!    ! See comment in GP_Variables
!!
!!    bioflo_map(:,1) = (/ SPECIES_NITRATE, &
!!                         SPECIES_AMMONIUM, &
!!                         SPECIES_DISSOLVED_ORGANIC_NITROGEN, &
!!                         SPECIES_DETRITUS, &
!!                         SPECIES_BACTERIA, &
!!                         SPECIES_PHYTOPLANKTON, &
!!
!!                         SPECIES_ZOOPLANKTON /)
!!    ! Since indexes are all negative, take the absolute value
!!
!!    bioflo_map = abs(bioflo_map)
!!
!!    !-------------------------------------------------------------------------------
!!
!!    FORCING_MIXED_LAYER_DEPTH  = -5001
!!    FORCING_MLD_CHANGE_NON_MOTILE = -5002
!!    FORCING_MLD_CHANGE_MOTILE = -5003
!!    FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004
!!
!!    Numerical_CODE_Initial_Conditions = (/aNO3, aNH4, DON, DET, bact, phyto, zoo/)
!!    Numerical_CODE_Forcing_Functions = 0.0D+0

!!else



!-------------------------------------------------------------------------------

!if( myid == 0 )then
!    write(6,'(A,1x,I5)') 'inmod: n_CODE_equations ', n_CODE_equations
!endif ! myid == 0

!if( L_myprint  .and. myid == 1 )then
!    write(myprint_unit,'(A,1x,I5)') 'inmod: n_CODE_equations ', n_CODE_equations
!endif ! L_myprint  .and. myid == 1


! See comment in GP_Variables

do  i = 1, n_CODE_equations
    bioflo_map(i,1) = -i
enddo ! i


!if( L_myprint  .and. myid == 1 )then
!    write(myprint_unit,'(A,10(1x,I5))') 'inmod: bioflo_map(:,1) ', bioflo_map(:,1)
!endif ! L_myprint  .and. myid == 1


! Since indexes are all negative, take the absolute value

bioflo_map = abs(bioflo_map)


!if( L_myprint  .and. myid == 1 )then
!    write(myprint_unit,'(A,10(1x,I5))') 'inmod:2 abs bioflo_map(:,1) ', bioflo_map(:,1)
!endif ! L_myprint  .and. myid == 1


!-----------------------------------------------------------------------------------

!if( L_myprint  .and. myid == 1 )then
!    write(myprint_unit,'(A,2(1x,E15.7))') 'inmod: Numerical_CODE_Initial_Conditions ', &
!                                                  Numerical_CODE_Initial_Conditions
!endif ! L_myprint  .and. myid == 1


Numerical_CODE_Forcing_Functions = 0.0D+0

!write(6,'(A,4(1x,E15.7))') 'inmod: Numerical_CODE_Forcing_Functions ', &
!                                   Numerical_CODE_Forcing_Functions

btmp(1:n_code_equations) = 0.0D0


!!endif  ! model == 'fasham' ....



!-------------------------------------------------------------------------------


!if( myid == 0 )then
!    write(6,'(//A)') 'inmod: call Build_Trees  '
!endif ! myid == 0

!if( L_myprint  .and. myid == 1 )then
!    write(myprint_unit,'(//A)') 'inmod: call Build_Trees  '
!endif ! L_myprint  .and. myid == 1



call Build_Trees( GP_Trees(:, 1) ,  buildTrees )


!if( myid == 0 )then
!    write(6,'(A//)') 'inmod: aft call Build_Trees  '
!endif ! myid == 0

!if( L_myprint  .and. myid == 1 )then
!    write(myprint_unit,'(A//)') 'inmod: aft call Build_Trees  '
!endif ! L_myprint  .and. myid == 1



!call Deserialize_Trees( GP_Trees(:,:,:), &
!                        n_Trees, n_Tracked_Resources, output_dir )


!-------------------------------------------------------------------------------


end subroutine Initialize_Model







subroutine DoForcing(b_tmp_local, time_step_fraction, i_Time_Step)

use fasham_variables_module
use GP_variables_module

implicit none

real (kind=8) :: b_tmp_local(n_CODE_Equations)
real (kind=8) :: time_step_fraction, day, h, hplus, aMLD, aJ
integer (kind=4) :: i_Time_Step

!------------------------------------------------------------------------------------

!!!date=(i_Time_Step+time_step_fraction)*Delta_Time_in_Days/(365.D+0)  ! number of years
!!!thour=mod(((i_Time_Step+time_step_fraction)*Delta_Time_in_Days*24),24.D+0) ! time of day in hours
!!!dayn=(i_Time_Step+time_step_fraction)*Delta_Time_in_Days ! day number

date=(i_Time_Step+time_step_fraction)* dt /(365.D+0)  ! number of years
thour=mod(((i_Time_Step+time_step_fraction)* dt *24),24.D+0) ! time of day in hours
dayn=(i_Time_Step+time_step_fraction)* dt ! day number
day=mod(dayn,365.D+0) ! year day [0.D+0 to 365.D+0]

!write(6,'(A,1x,I6,4(1x,E15.7))') 'dof: i_time_step, dt, date, thour, day ', &
!                                       i_time_step, dt, date, thour, day

call mldforce(day, h, aMLD)

call JQforce(b_tmp_local, day, aMLD, aJ)

if( h .ge. 0.D+0) then
    hplus=h
else
    hplus=0.D+0
endif

Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MLD_CHANGE_MOTILE))         = h
Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MLD_CHANGE_NON_MOTILE))     = hplus
Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MIXED_LAYER_DEPTH))         = aMLD
Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_LIGHT_LIMITED_GROWTH_RATE)) = aJ

!write(6,'(A,1x,I6,4(1x,E15.7))') 'dof: i_time_step, h, hplus, aMLD, aJ ', &
!                                       i_time_step, h, hplus, aMLD, aJ 

return

end subroutine



subroutine SecondaryForcing()
    ! Do nothing - no secondary forcing
end subroutine



subroutine Model_Diagnostics()
    !   TODO: Create this routine if need be
end subroutine
