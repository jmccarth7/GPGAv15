subroutine init_values( icall  )


use GP_parameters_module
use GP_variables_module
use Runge_Kutta_Variables_module

implicit none


integer,intent(in)  :: icall




!-------------------------------------------------------------------------

write(6,'(/A,1x,A)')  'iv: model ', trim(model)
write(6,'(A,1x,I6/)') 'iv: icall ', icall


if( trim(model) == 'NPZ' ) then

    call init_values_NPZ( icall )
    if( icall == 0 ) return

elseif( trim(model) == 'LV' )then

    call init_values_LV( icall )
    if( icall == 0 ) return

endif ! trim(model) == 'NPZ'



return

END subroutine init_values
