subroutine GP_Tree_Build
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use GP_Parameters_module
use GA_Parameters_module
use GP_Variables_module
use GA_Variables_module

implicit none

real(kind=4) :: cff
real(kind=8) :: dff

integer(kind=4) :: n_Nodes_at_Level
integer(kind=4) :: i_Level_Node
integer(kind=4) :: Node_Function
integer(kind=4) :: Node_Variable

!-----------------------------------------------------------------------------


GP_Adult_Population_Node_Type=-9999 ! set all to null [-9999]



do  i_GP_Individual=1,n_GP_Individuals       ! for each GP individual

    do  i_Tree=1,n_Trees                       ! for each GPCODE tree

        call random_number(cff) ! uniform random number generator

        if( cff .le. GP_Tree_Probability ) then  ! go ahead and put in an equation

            GP_Adult_Population_Node_Type(i_GP_Individual,1,i_Tree)=0 ! set the first node to zero

            i_Node=0
            do  i_Level=1,n_Levels-1

                n_Nodes_at_Level=int(2**(i_Level-1))
                do  i_Level_Node=1,n_Nodes_at_Level

                    i_Node=i_Node+1

                    if( GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) &
                                                                          .eq. 0 ) then   

                        ! randomly decide function or terminal

                        call random_number(cff) ! uniform random number generator


                        if( cff .lt. Node_Probability(i_Level) ) then  ! set as a terminal

                            call random_number(cff) ! uniform random number generator
                            node_function=int( float(n_Node_Functions)*cff  + 1.0 )

                            ! should not need to do this but compilers may vary w.r.t. 
                            ! possible rand = 1.

                            !if( Node_Function .gt. n_Node_Functions ) then
                            !    Node_Function  =   n_Node_Functions
                            !endif !  Node_Function .gt. n_Node_Functions

                            Node_Function = min( Node_Function, n_Node_Functions )
                            GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) = &
                                                                         Node_Function

                            !  set the node vs terminal selection capability 
                            !  for the node inputs at the next level

                            if( i_Level .lt. N_Levels-1 ) then
     
                                ! set the node lowel level inputs to open

                                GP_Adult_Population_Node_Type(i_GP_Individual,2*i_Node   ,i_Tree)=0          
                                GP_Adult_Population_Node_Type(i_GP_Individual,2*i_Node +1,i_Tree)=0      

                            else

                                ! complete setting the node lowest level nodes with terminals

                                GP_Adult_Population_Node_Type(i_GP_Individual,2*i_Node   ,i_Tree)=-1         
                                GP_Adult_Population_Node_Type(i_GP_Individual,2*i_Node +1,i_Tree)=-1     

                            endif !   i_Level .lt. N_Levels-1


                        else ! set it as a Parameter or Variable at a later point in the code

                            ! set a parameter or variable later
    
                            GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)=-1               


                        endif !   cff .lt. Node_Probability(i_Level)

                    endif !  GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) .eq. 0

                enddo !  i_Level_Node

            enddo !  i_Level

        endif !   cff .le. GP_Tree_Probability

    enddo !  i_Tree

enddo !  i_GP_Individual

!--------------------------------------------------------------------------------------------------

! randomly fill the terminals of the GP_Adult_Population_Node_Type array 
! with parameter or variable 'types'


do  i_GP_Individual=1,n_GP_Individuals

    do  i_Tree=1,n_Trees

        i_Node=0
        do  i_Level=1,n_Levels

            n_Nodes_at_Level=int(2**(i_Level-1))
            do  i_Level_Node=1,n_Nodes_at_Level

                i_Node=i_Node+1

                if( GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) .eq. -1) then

                    call random_number(cff)   ! uniform random number generator

                    if( cff .le. GP_Set_Terminal_to_Parameter_Probability ) then    

                        ! Set the Terminal to a Parameter

                        call random_number(cff) ! uniform random number generator

                        ! One of the OBSERVATIONS, one for each equations N, P, Z, etc.

                        Node_Variable=int( float(n_CODE_Equations)*cff + 1.0 )  


                        !  should not need to do this but compilers may vary w.r.t. 
                        !  possible rand = 1.
                        !if( Node_Variable .gt. n_CODE_Equations) then
                        !    Node_Variable=n_CODE_Equations
                        !endif !   Node_Variable .gt. n_CODE_Equations

                        Node_Variable = min( Node_Variable, n_CODE_Equations )

                        GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)=-Node_Variable


                    else  

                        ! set as a random parameter

                        ! The setting to zero allows the parameters to be set in GA_lmdif

                        GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree)=0  


                    endif !   cff .le. GP_Set_Terminal_to_Parameter_Probability

                endif !   GP_Adult_Population_Node_Type(i_GP_Individual,i_Node,i_Tree) .eq. -1

            enddo !  i_Level_Node

        enddo !  i_Level

    enddo !  i_Tree

enddo !  i_GP_Individual

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return
end subroutine GP_Tree_Build
!23456789012345678901234567890123456789012345678901234567890123456789012
