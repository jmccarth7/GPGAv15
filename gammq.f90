FUNCTION gammq(a,x) 

use mpi                                                                                                   
use mpi_module

implicit none

!U    USES gcf,gser                                                     
REAL(kind=8) ::  a,gammq,x 
REAL(kind=8) ::  gammcf,gamser,gln 

!--------------------------------------------------------------

if( x.lt.0.0d0.or.a.le.0.0d0 )then
    call MPI_FINALIZE(ierr)
    stop 'bad arguments in gammq' 
endif

if( x.lt.a+1.0d0 )then 
    call gser(gamser,a,x,gln) 
    gammq=1.0d0-gamser 
else 
    call gcf(gammcf,a,x,gln) 
    gammq=gammcf 
endif 

return 

END                                           
