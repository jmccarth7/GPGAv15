subroutine  count_parens( test_string, paren_type, n_parens )


use GP_Parameters_module
use GP_Variables_module
use Runge_Kutta_Variables_module

implicit none



character(1),intent(in)  ::  paren_type 

character(1),parameter  ::  left_type  = 'L' 
character(1),parameter  ::  right_type = 'R' 

character(1),parameter  ::  left_paren  = '(' 
character(1),parameter  ::  right_paren = ')' 
character(1)  ::  search_char

integer :: n_parens 

integer :: i    
integer :: len_work
!integer :: j1    
!integer :: j2    
!integer :: indx 
!integer :: ksafe

!integer,parameter :: ksafe_max = 50

character(*), intent(in) ::  test_string

!---------------------------------------------------------------------

search_char = ' '
if( paren_type == left_type  ) search_char = left_paren
if( paren_type == right_type ) search_char = right_paren

if( search_char == ' ' ) stop 'bad search char'


len_work = len( trim( test_string ) )   

n_parens = 0 

do  i = 1, len_work                   

 
    if( test_string(i:i) == search_char )then
        n_parens = n_parens + 1 
    endif 

enddo  ! i  


!write(6,'(//A,1x,A)') &
!          'cpr:1 test_string :   ',  trim(test_string)
!write(6,'(A,1x,I6)') 'cpr:  n_parens = ', n_parens

return


end subroutine  count_parens
