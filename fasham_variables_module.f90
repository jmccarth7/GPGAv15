module Fasham_Variables_module

implicit none


real(kind=8) :: aNO3, aNH4, DON, DET, bact
real(kind=8) :: Q1, Q2, Q

real(kind=8) :: nitro
real(kind=8) :: phyto
real(kind=8) :: zoo

real(kind=8) :: par0, cloudy, cloud, z, delz, oktas, tday, ratio, trans, albedo, saltd, cosz
real(kind=8) :: zenith, alatd, alatr, dec, th, E_0, value, sigma, parz, tmp, beta, tau
real(kind=8) :: U2, U1, S, enot, date, temp
real(kind=8) :: thour,  dayn


real(kind=8) :: prey_growth_rate
real(kind=8) :: predator_feeding_rate
real(kind=8) :: predator_assim
real(kind=8) :: predator_biomass_feeding_rate
real(kind=8) :: predator_mortality_rate

!-----------------------------------------------------------------------------------------

real(kind=8),parameter ::     alpha   =  0.025D+0    ! initial slope of the P-I curve [(W m-2)-1 d-1]
real(kind=8),parameter ::     aK1     =  0.5D+0      ! half-saturation for phytoplankton NO3 uptake [mMol N m-3]
real(kind=8),parameter ::     aK2     =  0.5D+0      ! half-saturation for phytoplankton NH4 uptake [mMol N m-3]
real(kind=8),parameter ::     amu1    =  0.045D+0    ! phytoplankton specific mortality rate [d-1]
real(kind=8),parameter ::     akc     =  0.03D+0     ! light attenuation by phytoplankton [m^2 mMol N)-1]
real(kind=8),parameter ::     gamma1  =  0.05D+0     ! fraction of total primary production that is exuded [n.d.]
real(kind=8),parameter ::     phi     =  1.5D+0      ! phytoplankton ammonium inhibition parameter [(mMol N)-1]
real(kind=8),parameter ::     g       =  1.0D+0      ! maximum zooplankton growth rate [d-1]
real(kind=8),parameter ::     beta1   =  0.75D+0     ! zooplankton assimilation efficiency of zooplankton [n.d.]
real(kind=8),parameter ::     beta2   =  0.75D+0     ! zooplankton assimilation efficiency of phytoplankton [n.d.]
real(kind=8),parameter ::     beta3   =  0.75D+0     ! zooplankton assimilation efficiency of bacteria [n.d.]
real(kind=8),parameter ::     amu2    =  0.1D+0      ! zooplankton specific excretion rate [d-1]
real(kind=8),parameter ::     amu5    =  0.05D+0     ! zooplankton specific mortality rate [d-1]
real(kind=8),parameter ::     aK3     =  1.0D+0      ! zooplankton half-saturation conts. for ingestion [d-1]
real(kind=8),parameter ::     omega   =  0.33D+0     ! detrital fraction of zooplankton mortality [n.d.]
real(kind=8),parameter ::     epsilon =  0.75D+0     ! ammonium fraction of zooplankton excretion [n.d.]
real(kind=8),parameter ::     Vb      =  2.0D+0      ! bacteria maximum growth rate [d-1]
real(kind=8),parameter ::     Vp      =  2.9D+0      ! phyto maximum growth rate [d-1]
real(kind=8),parameter ::     amu3    =  0.05D+0     ! bacteria specific excretion rate [d-1]
real(kind=8),parameter ::     aK4     =  0.5D+0      ! bacteria half-saturation rate for uptake [(mMol N) m-3]
real(kind=8),parameter ::     eta     =  0.6D+0      ! ammonium/DON uptake ratio [n.d.]
real(kind=8),parameter ::     amu4    =  0.05D+0     ! detrital breakdown rate [d-1]
real(kind=8),parameter ::     V       =  1.0D+0      ! detrital sinking rate [m d-1]
real(kind=8),parameter ::     p1      =  1.0D+0      ! zooplankton preference for phytoplankton [n.d.]
real(kind=8),parameter ::     p2      =  1.0D+0      ! zooplankton preference for bacteria [n.d.]
real(kind=8),parameter ::     p3      =  1.0D+0      ! zooplankton preference for detritus [n.d.]
real(kind=8),parameter ::     aN0     =  2.0D+0      ! concentration of NO3 below the mixed-layer [(mMol N) m-3]


!-----------------------------------------------------------------------------------------

! parameters as in Table 1; Fasham et al. [JMR, 48, 591-639, 1990]

real(kind=8),parameter ::     akw = 0.04D+0     ! light attenuation due to sea water [m-1]
real(kind=8),parameter ::     am  = 0.1D+0      ! cross-thermocline mixing rate [m d-1]


!-----------------------------------------------------------------------------------------

! Variables to define indexes of species in b_tmp and in Numerical Code Solutions

integer(kind=4),parameter :: SPECIES_NITRATE                      = -1
integer(kind=4),parameter :: SPECIES_AMMONIUM                     = -2
integer(kind=4),parameter :: SPECIES_DISSOLVED_ORGANIC_NITROGEN   = -3
integer(kind=4),parameter :: SPECIES_DETRITUS                     = -4
integer(kind=4),parameter :: SPECIES_BACTERIA                     = -5
integer(kind=4),parameter :: SPECIES_PHYTOPLANKTON                = -6
integer(kind=4),parameter :: SPECIES_ZOOPLANKTON                  = -7

integer(kind=4) :: SPECIES_Zoo
integer(kind=4) :: SPECIES_Phyto
integer(kind=4) :: SPECIES_Nitro



! Variables to define indexes of forcing functions in Numerical_CODE_Forcing_Functions

integer(kind=4),parameter :: FORCING_MIXED_LAYER_DEPTH         = -5001
integer(kind=4),parameter :: FORCING_MLD_CHANGE_NON_MOTILE     = -5002
integer(kind=4),parameter :: FORCING_MLD_CHANGE_MOTILE         = -5003
integer(kind=4),parameter :: FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004



end module Fasham_Variables_module
