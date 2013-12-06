subroutine Initialize_Model(buildTrees)

                                                                                                                                          
use mpi
use mpi_module

use fasham_variables_module
use GP_model_parameters_module
use GP_variables_module
    
implicit none
    
logical :: buildTrees

integer :: i
    
!----------------------------------------------------------------------------------------------------

if( myid == 0 )then
    write(6,'(/A/)') 'inmod: entry  Initialize_Model '
endif ! myid == 0 


! Set Variables


! initialize the biological data fields


nitro   = 0.0D0 
phyto   = 0.0D0 ! 30.D0 ! Phytoplankton     [mmol N m-3]
zoo     = 0.0D0 ! 2.0D0 ! Zooplankton       [mmol N m-3]

!----------------------------------------------------------------------------------------------------

! Enumerations that represent model variables. 
! These are used by the binary tree parsing algorithm 
! to select the index of the species or forcing function variable's value

!SPECIES_Nitro    = -1
!SPECIES_Phyto    = -2
!SPECIES_Zoo      = -3

!write(6,'(A,1x,I5)') 'inmod: SPECIES_Phyto ', SPECIES_Phyto
!write(6,'(A,1x,I5)') 'inmod: SPECIES_Zoo   ', SPECIES_Zoo

if( myid == 0 )then
    write(6,'(A,1x,I5)') 'inmod: n_CODE_equations ', n_CODE_equations 
endif ! myid == 0 


! See comment in GP_Variables

do  i = 1, n_CODE_equations

    bioflo_map(i,1) = -i

enddo ! i 

if( myid == 0 )then
    write(6,'(A,10(1x,I5))') 'inmod: bioflo_map(:,1) ', bioflo_map(:,1)
endif ! myid == 0 


! Since indexes are all negative, take the absolute value

bioflo_map = abs(bioflo_map)


if( myid == 0 )then
    write(6,'(A,10(1x,I5))') 'inmod:2 abs bioflo_map(:,1) ', bioflo_map(:,1)
endif ! myid == 0 


!----------------------------------------------------------------------------------------------------

!    FORCING_MIXED_LAYER_DEPTH  = -5001
!    FORCING_MLD_CHANGE_NON_MOTILE = -5002
!    FORCING_MLD_CHANGE_MOTILE = -5003
!    FORCING_LIGHT_LIMITED_GROWTH_RATE = -5004

!Numerical_CODE_Initial_Conditions = (/aNO3, aNH4, DON, DET, bact, phyto, zoo/)

!Numerical_CODE_Initial_Conditions = (/phyto, zoo/)


if( myid == 0 )then
    write(6,'(A,2(1x,E15.7))') 'inmod: Numerical_CODE_Initial_Conditions ', &
                                       Numerical_CODE_Initial_Conditions
endif ! myid == 0 


Numerical_CODE_Forcing_Functions = 0.0D+0 

!write(6,'(A,2(1x,E15.7))') 'inmod: Numerical_CODE_Forcing_Functions ', &
!                                   Numerical_CODE_Forcing_Functions



if( buildTrees ) then

    if( myid == 0 )then
        write(6,'(//A)') 'inmod: call Build_Trees  '
    endif ! myid == 0 


    call Build_Trees( GP_Trees(:,1) )


    if( myid == 0 )then
        write(6,'(A//)') 'inmod: aft call Build_Trees  '
    endif ! myid == 0 

else

    call Deserialize_Trees( GP_Trees(:,:), n_Trees, n_Tracked_Resources, output_dir )

endif ! buildTrees



end subroutine Initialize_Model




!subroutine DoForcing(b_tmp_local, time_step_fraction, i_Time_Step)
!
!use Fasham_Variables
!!use Model_Parameters
!use GP_model_parameters_module
!!use GP_Variables
!use GP_variables_module
!    
!implicit none
!    
!real (kind=8) :: b_tmp_local(n_CODE_Equations)
!real (kind=8) :: time_step_fraction, day, h, hplus, aMLD, aJ
!integer (kind=4) :: i_Time_Step
!
!!------------------------------------------------------------------------------------
!    
!date=(i_Time_Step+time_step_fraction)*Delta_Time_in_Days/(365.D+0)  ! number of years
!thour=mod(((i_Time_Step+time_step_fraction)*Delta_Time_in_Days*24),24.D+0) ! time of day in hours
!dayn=(i_Time_Step+time_step_fraction)*Delta_Time_in_Days ! day number
!day=mod(dayn,365.D+0) ! year day [0.D+0 to 365.D+0]
!    
!call mldforce(day, h, aMLD)
!
!call JQforce(b_tmp_local, day, aMLD, aJ)
!    
!if( h .ge. 0.D+0) then
!    hplus=h
!else
!    hplus=0.D+0
!endif
!    
!Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MLD_CHANGE_MOTILE)) = h
!Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MLD_CHANGE_NON_MOTILE)) = hplus
!Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_MIXED_LAYER_DEPTH)) = aMLD
!Numerical_CODE_Forcing_Functions(abs(5000 + FORCING_LIGHT_LIMITED_GROWTH_RATE)) = aJ
!    
!end subroutine
!
!
!
!subroutine SecondaryForcing()
!    ! Do nothing - no secondary forcing
!end subroutine
!
!
!
!subroutine Model_Diagnostics()
!    !   TODO: Create this routine if need be
!end subroutine
