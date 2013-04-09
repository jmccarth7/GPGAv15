program prog_read50


implicit none

integer(kind=4), parameter :: iunit = 50
integer(kind=4), parameter :: n_Levels=4
integer(kind=4), parameter :: n_CODE_Equations=2
integer(kind=4), parameter :: n_Trees=((n_CODE_Equations+1)**2)-(n_CODE_Equations+1)
integer(kind=4), parameter :: n_Nodes=int(2**n_Levels)-1

integer(kind=4) :: GP_Node_Type_for_Plotting(9,n_Nodes,n_Trees)

integer(kind=4) :: i
integer(kind=4) :: nrec
integer(kind=4) :: istat
integer(kind=4) :: itree
integer(kind=4) :: inode

!--------------------------------------------------------------------------

open( iunit, file = 'fort.50', form = 'unformatted', &
      access='sequential', status = 'old' ) 

nrec = 0
readloop:&
do

    read(iunit, iostat = istat ) GP_Node_Type_for_Plotting
    if( istat /= 0 ) exit readloop

    nrec = nrec + 1

    write(6,'(/A,1x,I10)') 'nrec = ', nrec                           
    do  i = 1, 9
        do  itree = 1, n_trees

            !do  inode = 1, n_nodes
            !    write(6,'(A,4(1x,I10))') &
            !         'i, itree, inode, GP_Node_Type(i, inode, itree) ', &
            !          i, itree, inode, GP_Node_Type_for_Plotting(i, inode, itree)
            !enddo ! inode = 1, n_nodes


            write(6,'(A,17(1x,I3))') &
                     'i, itree, GP_Node_Type(i, 1:n_nodes, itree) ', &
                      i, itree, GP_Node_Type_for_Plotting(i, 1:n_nodes, itree)


        enddo  ! itree

    enddo ! i 


enddo readloop


close( iunit )

stop 


end program prog_read50
