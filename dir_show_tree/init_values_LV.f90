subroutine init_values_LV( icall  )


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


use GP_parameters_module
use GP_variables_module
use Runge_Kutta_Variables_module

implicit none


integer,intent(in)  :: icall

logical :: LV_model1 = .TRUE.

!-------------------------------------------------------------------------

if(  icall  == 0  )then


    n_levels    = 4
    n_functions = 4

    n_CODE_equations = 2

    n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)

    n_nodes = int(2**n_levels)-1

    n_maximum_number_parameters = n_CODE_equations +  n_nodes


    n_linked_parms = 0
    n_linked_parms_dim = max( n_linked_parms, 1 )

    write(6,'(A,1x,I6)') 'ivLV: n_levels           ', n_levels
    write(6,'(A,1x,I6)') 'ivLV: n_functions        ', n_functions
    write(6,'(A,1x,I6)') 'ivLV: n_CODE_equations   ', n_CODE_equations
    write(6,'(A,1x,I6)') 'ivLV: n_trees            ', n_trees
    write(6,'(A,1x,I6)') 'ivLV: n_nodes            ', n_nodes
    write(6,'(A,1x,I6)') 'ivLV: n_linked_parms     ', n_linked_parms
    write(6,'(A,1x,I6)') 'ivLV: n_linked_parms_dim ', n_linked_parms_dim
    write(6,'(/A,1x,I6/)') 'ivLV: n_maximum_number_parameters  ', &
                                  n_maximum_number_parameters




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

! Lotka_Volterra_Example_Set_Up


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! This is the tree representation of the CODE System
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! this example is a simple Lotka-Volterra model.

! dP/dt = (grow * P)  - (graze * P * Z)
! dZ/dt = (graze * P * Z) - ((1 - efficiency) * graze * P *  Z) - (amort * Z)
! [Note: In this example, the (1-efficiency) parameter is set to a new unique paramater, 'effic']

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! 1: [ 3, 0,-1, N, N, N, N, N, N, N, N, N, N, N, N]
! 2: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 3: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 4: [ 3, 0, 3, N, N,-1,-2, N, N, N, N, N, N, N, N]
! 5: [ 1, 3, 3, 0,-2, 3, 3, N, N, N, N, 0,-2, 0,-1]
! 6: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!    [01,02,03,04,05,06,07,08,09,10,11,12,13,14,15]

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

write(6,'(A, 4x,L1)') 'ivLV: LV_model1 ', &
                             LV_model1

! Initial Conditions

if( LV_model1 )then
    Runge_Kutta_Initial_Conditions(1) = 30.0D+0  ! [prey]         [mmol N m-3]
    Runge_Kutta_Initial_Conditions(2) = 2.0D+0   ! [predator]     [mmol N m-3]
else
    Runge_Kutta_Initial_Conditions(1) = 19.66561   ! 30.0D+0  ! [prey]         [mmol N m-3]
    Runge_Kutta_Initial_Conditions(2) = 0.3960451  ! 2.0D+0   ! [predator]     [mmol N m-3]
endif  ! LV_model1

write(6,'(A, 2(1x,F10.2))') 'ivLV: Runge_Kutta_Initial_Conditions(1:2) ', &
                                   Runge_Kutta_Initial_Conditions(1:2)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Node_Probability = (/0.8d0,0.6d0,0.4d0,0.d0/)  ! NOTE: Last value MUST BE 0.0!!!]


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! This term is calculated in units of [mg N d-1]
GP_Individual_Node_Type(1,1) = 3   ! '*'
GP_Individual_Node_Type(2,1) = 0            ! prey growth rate

if( LV_model1 )then
    GP_Individual_Node_Parameters(2,1) = 0.4d0  ! [0.04, 0.4; prey growth rate [d-1]
else
    GP_Individual_Node_Parameters(2,1) = 5.599795d0  ! 0.4    ! [0.04, 0.4; prey growth rate [d-1]
endif  ! LV_model1


GP_Individual_Node_Type(3,1) = -1  ! Phyto
GP_Individual_Node_Type(1,4) = 3   ! '*'
GP_Individual_Node_Type(2,4) = 0            ! predator biomass-specific feeding rate [d-1]

if( LV_model1 )then
    GP_Individual_Node_Parameters(2,4) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
else
    GP_Individual_Node_Parameters(2,4) = 1.56521d0 !0.02! predator biomass-specific feeding rate [d-1]
endif  ! LV_model1


GP_Individual_Node_Type(3,4) = 3   ! '*'
GP_Individual_Node_Type(6,4) = -1  ! Phyto
GP_Individual_Node_Type(7,4) = -2  ! Zoo

GP_Individual_Node_Type(1,5) = 1   ! '+'
GP_Individual_Node_Type(2,5) = 3   ! '*'
GP_Individual_Node_Type(3,5) = 3   ! '*'
GP_Individual_Node_Type(4,5) = 0           ! predator biomass-specific mortality rate [d-1]

if( LV_model1 )then
    GP_Individual_Node_Parameters(4,5) = 0.6d0 ! [0.1, 0.6; predator biomass-specific mortality rate [d-1]
else
    GP_Individual_Node_Parameters(4,5) = 0.8346865d-06 !0.6![ predator biomass-specific mortality rate [d-1]
endif  ! LV_model1


GP_Individual_Node_Type(5,5) = -2  ! Zoo
GP_Individual_Node_Type(6,5) = 3   ! '*'
GP_Individual_Node_Type(7,5) = 3   ! '*'
GP_Individual_Node_Type(12,5) = 0           ! predator assimilation efficiency [fraction 0<==>1]

if( LV_model1 )then
    GP_Individual_Node_Parameters(12,5) = 0.5d0 ! [0.2, 0.5; predator assimilation efficiency [fraction 0<==>1]
else
    GP_Individual_Node_Parameters(12,5) = 0.2416847d+01 ! 0.5!  predator assimilation efficiency [fraction 0<==>1]
endif  ! LV_model1


GP_Individual_Node_Type(13,5) = -2 ! Zoo
GP_Individual_Node_Type(14,5) = 0            ! predator biomass-specific feeding rate [d-1]
if( LV_model1 )then
    GP_Individual_Node_Parameters(14,5) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
else
    GP_Individual_Node_Parameters(14,5) = 0.2585400E+00  ! 0.02  ! predator biomass-specific feeding rate [d-1]
endif  ! LV_model1

GP_Individual_Node_Type(15,5) = -1 ! Phytoplankton

!-------------------------------------------------------------------------------                                            
                                                                                                                            
! strings                                                                                                                   
node_type_string = ' '  
node_parameters_string   = ' '    ! [0.04, 0.4; prey growth rate [d-1] 
                                                                                                                            
! This term is calculated in units of [mg N d-1]
!GP_Individual_Node_Type(1,1) = 3   ! '*'
node_type_string(1,1)= ' *  '      ! '*'                                                                               
!GP_Individual_Node_Type(2,1) = 0            ! prey growth rate
node_type_string(2,1)= ' pgr  '    ! prey growth rate

!if( LV_model1 )then
!    GP_Individual_Node_Parameters(2,1) = 0.4d0  ! [0.04, 0.4; prey growth rate [d-1]
!else
!    GP_Individual_Node_Parameters(2,1) = 5.599795d0  ! 0.4    ! [0.04, 0.4; prey growth rate [d-1]
!endif  ! LV_model1
node_parameters_string(2,1)=' pgr'    ! [0.04, 0.4; prey growth rate [d-1] 


!GP_Individual_Node_Type(3,1) = -1  ! Phyto
node_type_string(3,1)= ' P  '      ! Phyto               

!---------------------

!GP_Individual_Node_Type(1,4) = 3   ! '*'
node_type_string(1,4)= ' *  '      ! '*'                                                                               
!GP_Individual_Node_Type(2,4) = 0   ! predator biomass-specific feeding rate [d-1]
node_type_string(2,4)= ' PFR  '    ! predator biomass-specific feeding rate [d-1]

!if( LV_model1 )then
!    GP_Individual_Node_Parameters(2,4) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
!else
!    GP_Individual_Node_Parameters(2,4) = 1.56521d0 !0.02! predator biomass-specific feeding rate [d-1]
!endif  ! LV_model1

node_parameters_string(2,4)=' PFR'    ! 0.02! predator biomass-specific feeding rate [d-1] 


!GP_Individual_Node_Type(3,4) = 3   ! '*'
node_type_string(3,4)= ' *  '      ! '*'                                                                               
!GP_Individual_Node_Type(6,4) = -1  ! Phyto
node_type_string(6,4)= ' P  '      ! Phyto
!GP_Individual_Node_Type(7,4) = -2  ! Zoo
node_type_string(7,4)= ' Z  '      ! Zoo

!---------------------
!GP_Individual_Node_Type(1,5) = 1   ! '+'
node_type_string(1,5)= ' +  '      ! '+'                                                                               
!GP_Individual_Node_Type(2,5) = 3   ! '*'
node_type_string(2,5)= ' *  '          ! '*'                                                                               
!GP_Individual_Node_Type(3,5) = 3   ! '*'
node_type_string(3,5)= ' *  '          ! '*'                                                                               
!GP_Individual_Node_Type(4,5) = 0           ! predator biomass-specific mortality rate [d-1]
node_type_string(4,5)= 'pmor'          ! '*'                                                                               

!if( LV_model1 )then
!    GP_Individual_Node_Parameters(4,5) = 0.6d0 ! [0.1, 0.6; predator biomass-specific mortality rate [d-1]
!else
!    GP_Individual_Node_Parameters(4,5) = 0.8346865d-06 !0.6![ predator biomass-specific mortality rate [d-1]
!endif  ! LV_model1
node_parameters_string(4,5)='pmor'    ! 0.6![ predator biomass-specific mortality rate [d-1]


!GP_Individual_Node_Type(5,5) = -2  ! Zoo
node_type_string(5,5)= ' Z  '      ! Zoo
!GP_Individual_Node_Type(6,5) = 3   ! '*'
node_type_string(6,5)= ' *  '      ! '*'                                                                               
!GP_Individual_Node_Type(7,5) = 3   ! '*'
node_type_string(7,5)= ' *  '      ! '*'                                                                               
!GP_Individual_Node_Type(12,5) = 0           ! predator assimilation efficiency [fraction 0<==>1]
node_type_string(12,5)= 'assm'     ! predator assimilation efficiency [fraction 0<==>1]

!if( LV_model1 )then
!    GP_Individual_Node_Parameters(12,5) = 0.5d0 ! [0.2, 0.5; predator assimilation efficiency [fraction 0<==>1]
!else
!    GP_Individual_Node_Parameters(12,5) = 0.2416847d+01 ! 0.5!  predator assimilation efficiency [fraction 0<==>1]
!endif  ! LV_model1
node_parameters_string(12,5)='assm'    ! 0.5!  predator assimilation efficiency [fraction 0<==>1] 


!GP_Individual_Node_Type(13,5) = -2 ! Zoo
node_type_string(13,5)= ' Z  '     ! Zoo                                                                               
!GP_Individual_Node_Type(14,5) = 0  ! predator biomass-specific feeding rate [d-1]
node_type_string(14,5)= 'PFR '     ! predator biomass-specific feeding rate [d-1]                                      
!if( LV_model1 )then
!    GP_Individual_Node_Parameters(14,5) = 0.02d0 ! [0.0005, 0.02; predator biomass-specific feeding rate [d-1]
!else
!    GP_Individual_Node_Parameters(14,5) = 0.2585400E+00  ! 0.02  ! predator biomass-specific feeding rate [d-1]
!endif  ! LV_model1
node_parameters_string(14,5)='PFR '    ! Phytoplankton Maximum Growth Rate, Ranges between 0.20 <==> 3.0 [d-1]              

!GP_Individual_Node_Type(15,5) = -1 ! Phytoplankton
node_type_string(15,5)= ' P  '     ! Phytoplankton                                                                          


                                                                                                                            
!--------------------------------------------------------------------------------------       

write(6,'(A,1x,F10.2)') 'ivLV: GP_Individual_Node_Parameters(12,5) ', &
                               GP_Individual_Node_Parameters(12,5)
write(6,'(A,1(1x,I6))') 'ivLV: GP_Individual_Node_Type(1,1)        ', &
                               GP_Individual_Node_Type(1,1)
write(6,'(A, 1x,F10.2)') 'ivLV: GP_Individual_Node_Parameters(2,1)  ', &
                                GP_Individual_Node_Parameters(2,1)
write(6,'(A, 1x,F10.2)') 'ivLV: GP_Individual_Node_Parameters(2,4)  ', &
                                GP_Individual_Node_Parameters(2,4)
write(6,'(A, 1x,F10.2)') 'ivLV: GP_Individual_Node_Parameters(4,5)  ', &
                                GP_Individual_Node_Parameters(4,5)
write(6,'(A, 1x,F10.2)') 'ivLV: GP_Individual_Node_Parameters(14,5) ', &
                                GP_Individual_Node_Parameters(14,5)


if( n_linked_parms > 0 )then
    linked_parms(2, 1:n_linked_parms) = 0
endif ! n_linked_parms > 0

!--------------------------------------------------------------------------------------


return

END subroutine init_values_LV
