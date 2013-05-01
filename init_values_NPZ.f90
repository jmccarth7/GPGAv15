subroutine init_values_NPZ( icall  )


!     written by John R. Moisan [14 November 2012]
!     for GPCODE testing/developing
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! this example is a simple NPZ (Nitrogen-Phytoplankton-Zooplankton) model
! developed by Peter Franks to demonstrate the impact of various grazing
! formulations on a model's performance

! The paper formulates two seperate grazing formulations:
! FORMULATION A:  Ivlev Grazing Formulation
!   dN/dt = -(V_m * P*( N/(K_N+N) ) )  + (mort * P) + (death * Z) + (effic * beta * (1-e^(-lambda*P)) * Z)
!   dP/dt =  (V_m * P*( N/(K_N+N) ) )  - (mort * P) - (beta * (1-e^(-lambda*P)) * Z)
!   dZ/dt = ((1 - effic) * beta * (1-e^(-lambda*P)) * Z) - (death * Z)

! FORMULATION B: Mayzaud-Poulet Grazing Formulation
!   dN/dt = -(V_m * P*( N/(K_N+N) ) )  + (mort * P) + (death * Z) + (effic * beta * Lambda * P * (1-e^(-lambda*P)) * Z)
!   dP/dt =  (V_m * P*( N/(K_N+N) ) )  - (mort * P) - (beta * lambda * P * (1-e^(-lambda*P)) * Z)
!   dZ/dt = ((1 - effic) * beta * lambda * P * (1-e^(-lambda*P)) * Z) - (death * Z)


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


use mpi
use mpi_module
use GP_parameters_module
use GP_variables_module
use Runge_Kutta_Variables_module

implicit none


integer,intent(in)  :: icall

integer, parameter :: nbio   = 3
integer, parameter :: iNO3   = 1
integer, parameter :: iphyto = 2
integer, parameter :: izoo   = 3

!     set the model coefficients

real(kind=8), parameter :: phyto_growth_maximum = 2.0d0 ! Ranges between 0.20 < =  = > 3.0 [d-1]
real(kind=8), parameter :: zoo_grazing_maximum  = 1.5d0 ! Ranges between 0.16 < =  = > 1.5 [d-1]
real(kind=8), parameter :: grazing_Control      = 1.0d0 ! Ranges between 0.10 < =  = > 2.0 [d-1]
real(kind=8), parameter :: K_NO3                = 1.0d0 ! [ug-at N l-1]
real(kind=8), parameter :: phyto_mortality      = 0.1d0 ! [d-1]
real(kind=8), parameter :: zoo_death_rate       = 0.2d0 ! [d-1]
real(kind=8), parameter :: assim                = 0.3d0 ! [d-1]


logical Ivlev

!integer(kind=4) :: i_GP_individual                                                                           
!integer(kind=4) :: i_GP_Best_Parent                                                                          
!integer(kind=4) :: i_GP_Generation                                                                           
!integer(kind=4) :: j_GP_Individual                                                                           
integer(kind=4) :: i_Tree                                                                                    
integer(kind=4) :: i_Node            

!-------------------------------------------------------------------------

if(  icall  == 0  )then


    n_levels    = 5
    n_functions = 5

    n_CODE_equations=3

    n_trees=((n_CODE_equations+1)**2)-(n_CODE_equations+1)

    n_nodes = int(2**n_levels)-1

    n_maximum_number_parameters = n_CODE_equations +  n_nodes

    n_linked_parms = 2

    n_linked_parms_dim = max( n_linked_parms, 1 )

    if( myid == 0 )then
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_levels           ', n_levels
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_functions        ', n_functions
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_CODE_equations   ', n_CODE_equations
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_trees            ', n_trees
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_nodes            ', n_nodes
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_linked_parms     ', n_linked_parms
        write(GP_print_unit,'(A,1x,I6)') 'ivNPZ: n_linked_parms_dim ', n_linked_parms_dim
        write(GP_print_unit,'(/A,1x,I6/)') 'ivNPZ: n_maximum_number_parameters ', &
                                       n_maximum_number_parameters
    endif ! myid == 0



    return

endif ! icall == 0




! load the arrays




Ivlev = .true.



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

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     This is the tree representation of the NPZ equation set
!     This terms are calculated in
!          units of [d-1] relative to the source term
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


! need 12 integer strings to represent the full CODE set
!  1: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  2: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  3: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  4: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  5: [ 3, 3, 0, 6,-2, N, N, 0,-1, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  6: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  7: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  8: [ 3,-2, 0, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
!  9: [ 3,-3, 3, N, N, 0, 5, N, N, N, N, N, N, 0,-2, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 10: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 11: [ 1, 3, 3, 0, 3,-3, 0, N, N, 3, 5, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N]
! 12: [ N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, 0,-3, 0,-2, N, N, N, N, N, N, N, N]
!     [01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]

!--------------------------------------------------------------------------------

! numbers

!     initialize the biological data fields
!     RENM: In the GPCODE these initial conditions will need
!           to be pooled into the parameter data set for
!           optimization in the mixed GA/LMDIF.f algorithm

!     the total must equal 2.0



Runge_Kutta_Initial_Conditions( iNO3 )   = 1.6D0   ! [NO3]            [mmol N m-3]
Runge_Kutta_Initial_Conditions( iphyto ) = 0.3D0   ! [Phytoplankton]  [mmol N m-3]
Runge_Kutta_Initial_Conditions( izoo )   = 0.1D0   ! [Zooplankton]    [mmol N m-3]




! numbers


GP_Individual_Node_Type(1,5)=3               ! '*'
GP_Individual_Node_Type(2,5)=3               ! '*'
GP_Individual_Node_Type(3,5)=0               ! Phytoplankton Maximum Growth Rate
GP_Individual_Node_Parameters(3,5)=2.0D+0    ! Phytoplankton Maximum Growth Rate, Ranges between 0.20 <==> 3.0 [d-1]
GP_Individual_Node_Type(4,5)=6               ! Michealis-Menton Term'
GP_Individual_Node_Type(5,5)=-2              ! Phytoplankton
GP_Individual_Node_Type(8,5)=0               ! K_NO3, Half-Saturation Term for Michaelis-Menton Term  [ug-at N l-1]
GP_Individual_Node_Parameters(8,5)=1.0D+0    ! K_NO3, Half-Saturation Term for Michaelis-Menton Term  [ug-at N l-1]
GP_Individual_Node_Type(9,5)=-1              ! [NO3]

GP_Individual_Node_Type(1,8)=3               ! '*'
GP_Individual_Node_Type(2,8)=-2              ! Phytoplankton
GP_Individual_Node_Type(3,8)=0               ! Phytoplankton Mortality Rate; [d-1]
GP_Individual_Node_Parameters(1,8)=0.1D+0    ! Phytoplankton Mortality Rate; [d-1]

GP_Individual_Node_Type(1,9)=3               ! '*'
GP_Individual_Node_Type(2,9)=-3              ! Zooplankton
GP_Individual_Node_Type(3,9)=3               ! '*'
GP_Individual_Node_Type(6,9)=0               ! Zooplankton Maximum Grazing Rate
GP_Individual_Node_Parameters(6,9)=1.5D+0    ! Zooplankton Maximum Grazing Rate; Ranges between 0.16 <==> 1.5 [d-1]
GP_Individual_Node_Type(7,9)=5               ! Ivlev Exponential Function (1 - e^-abs(left*right))
!off GA_Individual_Node_Type(7,4)=7          ! Mayzaud-Poulet Exponential Function abs(left*right)*(1 - e^-abs(left*right))
GP_Individual_Node_Type(14,9)=0              ! Grazing_Control
GP_Individual_Node_Parameters(14,9)=1.0D+0   ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
GP_Individual_Node_Type(15,9)=-2             ! Phytoplankton

GP_Individual_Node_Type(1,11)=1              ! '+'
GP_Individual_Node_Type(2,11)=3              ! '*'
GP_Individual_Node_Type(3,11)=3              ! '*'
GP_Individual_Node_Type(4,11)=0              ! Zooplankton Assimilation Rate; [d-1]
GP_Individual_Node_Parameters(4,11)=0.3D+0   ! Zooplankton Assimilation Rate; [d-1]
GP_Individual_Node_Type(5,11)=3              ! '*'
GP_Individual_Node_Type(6,11)=-3             ! Zooplankton
GP_Individual_Node_Type(7,11)=0              ! Zooplnakton Mortality Rate; [d-1]
GP_Individual_Node_Parameters(7,11)=0.2D+0   ! Zooplnakton Mortality Rate; [d-1]
GP_Individual_Node_Type(10,11)=3             ! '*'
GP_Individual_Node_Type(11,11)=5             ! Ivlev Exponential Function (1 - e^-abs(left*right))
!off GA_Individual_Node_Type(11,11)=7        ! Mayzaud-Poulet Exponential Function abs(left*right)*(1 - e^-abs(left*right))
GP_Individual_Node_Type(20,11)=0             ! Zooplankton Maximum Grazing Rate
GP_Individual_Node_Parameters(20,11)=1.5D+0  ! Zooplankton Maximum Grazing Rate; Ranges between 0.16 <==> 1.5 [d-1]
GP_Individual_Node_Type(21,11)=-3            ! Zooplankton
GP_Individual_Node_Type(22,11)=0             ! Grazing_control
GP_Individual_Node_Parameters(22,11)=1.0D+0  ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
GP_Individual_Node_Type(23,11)=-2            ! Phytoplankton

!-------------------------------------------------------------------------------

! strings

!node_type_string(1,5)= ' *  '          ! '*'
!node_type_string(2,5)= ' *  '               ! '*'
!node_type_string(3,5)= ' V  '               ! Phytoplankton Maximum Growth Rate
!node_parameters_string(3,5)=' V  '    ! Phytoplankton Maximum Growth Rate, Ranges between 0.20 <==> 3.0 [d-1]
!node_type_string(4,5)= ' I2 '         ! Michealis-Menton Term'
!node_type_string(5,5)= ' P  '         ! Phytoplankton
!node_type_string(8,5)= 'KNO3'         ! K_NO3, Half-Saturation Term for Michaelis-Menton Term  [ug-at N l-1]
!node_parameters_string(8,5)= 'KNO3'   ! K_NO3, Half-Saturation Term for Michaelis-Menton Term  [ug-at N l-1]
!node_type_string(9,5)= ' N  '         ! [NO3]
!
!node_type_string(1,8)= ' *  '         ! '*'
!node_type_string(2,8)= ' P  '         ! Phytoplankton
!node_type_string(3,8)= ' m  '         ! Phytoplankton Mortality Rate; [d-1]
!node_parameters_string(1,8)='pmor'    ! Phytoplankton Mortality Rate; [d-1]
!
!node_type_string(1,9)= ' *  '         ! '*'
!node_type_string(2,9)= ' Z  '         ! Zooplankton
!node_type_string(3,9)= ' *  '         ! '*'
!node_type_string(6,9)= ' b  '         ! Zooplankton Maximum Grazing Rate
!node_parameters_string(6,9)='zgm '    ! Zooplankton Maximum Grazing Rate; Ranges between 0.16 <==> 1.5 [d-1]
!node_type_string(7,9)= ' I1 '         ! Ivlev Exponential Function (1 - e^-abs(left*right))
!!off GA_Individual_Node_Type(7,4)=7          ! Mayzaud-Poulet Exponential Function abs(left*right)*(1 - e^-abs(left*right))
!node_type_string(14,9)= 'Rm '         ! Grazing_Control
!node_parameters_string(14,9)='GrCt'   ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
!node_type_string(15,9)= ' P  '        ! Phytoplankton
!
!node_type_string(1,11)= ' +  '        ! '+'
!node_type_string(2,11)= ' *  '        ! '*'
!node_type_string(3,11)= ' *  '        ! '*'
!node_type_string(4,11)= 'assm'        ! Zooplankton Assimilation Rate; [d-1]
!node_parameters_string(4,11)='assp'   ! Zooplankton Assimilation Rate; [d-1]
!node_type_string(5,11)= ' *  '        ! '*'
!node_type_string(6,11)= ' Z  '        ! Zooplankton
!node_type_string(7,11)= ' omg'        ! Zooplnakton Mortality Rate; [d-1]
!node_parameters_string(7,11)= 'zdr '  ! Zooplnakton Mortality Rate; [d-1]
!node_type_string(10,11)= ' *  '       ! '*'
!node_type_string(11,11)= ' I1 '       ! Ivlev Exponential Function (1 - e^-abs(left*right))
!!off GA_Individual_Node_Type(11,11)=7        ! Mayzaud-Poulet Exponential Function abs(left*right)*(1 - e^-abs(left*right))
!node_type_string(20,11)= ' b  '       ! Zooplankton Maximum Grazing Rate
!node_parameters_string(20,11)='zgm '  ! Zooplankton Maximum Grazing Rate; Ranges between 0.16 <==> 1.5 [d-1]
!node_type_string(21,11)= ' Z  '       ! Zooplankton
!node_type_string(22,11)= 'Rm  '       ! Grazing_control
!node_parameters_string(22,11)='GrCt'  ! Grazing Control; Ranges between 0.10 <==> 2.0 [d-1]
!node_type_string(23,11)= ' P  '       ! Phytoplankton

!--------------------------------------------------------------------------------------

! the parameter array will be filled with the same value for parms 7 and 11 ,
! and for parms 8 and 12

! the linked array will be used when random values are put into parameter values to
! ensure that the same random value will be in 7 and 11   and in 8 and 12


linked_parms(1,1) = 7
linked_parms(2,1) = 11

linked_parms(1,2) = 8
linked_parms(2,2) = 12

return

END subroutine init_values_NPZ
