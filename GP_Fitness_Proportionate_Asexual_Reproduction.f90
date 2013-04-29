subroutine GP_Fitness_Proportionate_Asexual_Reproduction
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use mpi
use mpi_module
use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff
!real(kind=8) :: dff

integer(kind=4) :: icff
integer(kind=4) :: i_GP_Asexual_Reproduction


!-----------------------------------------------------------------------------

icff = -1

i_GP_Individual = n_GP_Elitists 

if( myid == 0 )then
    write(GP_print_unit,'(/A,1x,I6)' ) &
          'gpfpar: n_GP_Asexual_Reproductions ', n_GP_Asexual_Reproductions
    write(GP_print_unit,'(A,1x,I6)' ) &
          'gpfpar: n_GP_Elitists         ', n_GP_Elitists
    write(GP_print_unit,'(A,1x,I6)' ) &
          'gpfpar: start i_GP_individual ', n_GP_Elitists  + 1
endif ! myid == 0 


do  i_GP_Asexual_Reproduction=1,n_GP_Asexual_Reproductions

    i_GP_Individual=i_GP_Individual+1
  
    call Random_Number(cff) ! uniform random number generator
  
    ! normalize to the maximum values so that the range is from [0. to 1.]
  
    do  j_GP_Individual=1,n_GP_Individuals
  
        if( myid == 0 )then
            write(GP_print_unit,'(A,1x,I6,2(1x,E15.7))' ) &
              'gpfpar: j_GP_Individual, cff, &
              &GP_Integrated_Population_Ranked_Fitness(j_GP_Individual) ', &
                       j_GP_Individual, cff, &
               GP_Integrated_Population_Ranked_Fitness(j_GP_Individual) 
        endif ! myid == 0 
  
        if( cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)) then
  
            icff=j_GP_Individual
            exit
  
        endif !   cff .le. GP_Integrated_Population_Ranked_Fitness(j_GP_Individual)
  
    enddo ! j_GP_Individual
  
    ! index to move over both the parent parameters and the individual fitness levels
  
    j_GP_Individual=icff
  
    if( myid == 0 )then
        write(GP_print_unit,'(/A,2(1x,I6)/)' ) 'gpfpar: j_GP_Individual, icff ', &
                                                        j_GP_Individual, icff
    endif ! myid == 0 
  
    GP_Child_Population_Node_Type(i_GP_Individual,1:n_Nodes,1:n_Trees) = &
       GP_Adult_Population_Node_Type(j_GP_Individual,1:n_Nodes,1:n_Trees)
  
    GP_Child_Population_Parameter_Solution(i_GP_Individual,1:n_Maximum_Number_Parameters) = &
       GP_Adult_Population_Parameter_Solution(j_GP_Individual,1:n_Maximum_Number_Parameters)
  
    ! give the child the adult's SSE value
  
  
    GP_Child_Population_SSE(i_GP_Individual) = GP_Adult_Population_SSE(j_GP_Individual)
  
  
    if( myid == 0 )then
        !write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !      'gpfpar: j_GP_individual, GP_Adult_Individual_SSE(j_GP_Individual)',&
        !               j_GP_individual, GP_Adult_Individual_SSE(j_GP_Individual)
        write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
              'gpfpar: j_GP_individual, GP_Adult_Population_SSE(j_GP_Individual)',&
                       j_GP_individual, GP_Adult_Population_SSE(j_GP_Individual)
      
        !write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
        !      'gpfpar: i_GP_individual, GP_Child_Individual_SSE(i_GP_Individual)',&
        !               i_GP_individual, GP_Child_Individual_SSE(i_GP_Individual)
        write(GP_print_unit,'(A,1x,I6,1x,E15.7)' ) &
              'gpfpar: i_GP_individual, GP_Child_Population_SSE(i_GP_Individual)',&
                       i_GP_individual, GP_Child_Population_SSE(i_GP_Individual)
  
    endif ! myid == 0
  
  
  
    !!!Run_GP_Calculate_Fitness(i_GP_Individual)=.false.
  

enddo ! i_GP_Asexual_Reproduction


return

end subroutine GP_Fitness_Proportionate_Asexual_Reproduction
