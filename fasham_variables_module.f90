module Fasham_Variables_module    

implicit none

real(kind=8) :: akw,am,alpha,aK1,aK2,amu1,akc,gamma1
real(kind=8) :: phi,g,beta1,beta2,beta3,amu2,amu5,aK3,omega,epsilon,Vb,amu3,aK4,eta

real(kind=8) :: amu4,V,p1,p2,p3,aN0,aNO3,aNH4,DON,DET,bact
real(kind=8) :: Q1,Q2,Q

real(kind=8) :: nitro 
real(kind=8) :: phyto
real(kind=8) :: zoo

real(kind=8) :: par0,cloudy,cloud,z,delz,oktas,tday,ratio,trans,albedo,saltd,cosz
real(kind=8) :: zenith,alatd,alatr,dec,th,E_0,value,sigma,parz,tmp,beta,tau
real(kind=8) :: U2,U1,S,enot,date,temp,Vp
real(kind=8) :: thour, dayn


real(kind=8) :: prey_growth_rate
real(kind=8) :: predator_feeding_rate
real(kind=8) :: predator_assim
real(kind=8) :: predator_biomass_feeding_rate
real(kind=8) :: predator_mortality_rate 


! Variables to define indexes of species in b_tmp and in Numerical Code Solutions

integer(kind=4) :: SPECIES_PHYTOPLANKTON
integer(kind=4) :: SPECIES_BACTERIA
integer(kind=4) :: SPECIES_ZOOPLANKTON
integer(kind=4) :: SPECIES_NITRATE
integer(kind=4) :: SPECIES_AMMONIUM
integer(kind=4) :: SPECIES_DETRITUS
integer(kind=4) :: SPECIES_DISSOLVED_ORGANIC_NITROGEN

integer(kind=4) :: SPECIES_Zoo                         
integer(kind=4) :: SPECIES_Phyto                       
integer(kind=4) :: SPECIES_Nitro                       



! Variables to define indexes of forcing functions in Numerical_CODE_Forcing_Functions

integer(kind=4) :: FORCING_LIGHT_LIMITED_GROWTH_RATE
integer(kind=4) :: FORCING_MLD_CHANGE_NON_MOTILE
integer(kind=4) :: FORCING_MLD_CHANGE_MOTILE
integer(kind=4) :: FORCING_MIXED_LAYER_DEPTH


end module Fasham_Variables_module
