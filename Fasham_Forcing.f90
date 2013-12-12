!     
! File:   Fasham_Forcing.f90
! Author: Dave
!
! Created on June 26, 2013, 1:37 PM
!
subroutine JQforce(species, day, aMLD, aJ)
    
use fasham_variables_module
use GP_variables_module

    implicit none
    
    real (kind=8), intent(in) :: species(n_CODE_Equations), day, aMLD
    real (kind=8), intent(out) :: aJ
    real (kind=8) :: daym(14),coktas(14),pi
    real (kind=8), parameter :: solar=1353.D+0 ! the solar max is from Brock, 1981
    integer (kind=4) :: iz,i
    
    ! Copy phytoplankton to 
    phyto = species(abs(SPECIES_PHYTOPLANKTON))
    
    data daym /0.D+0,16.D+0,46.D+0,75.D+0,105.D+0,136.D+0,166.D+0,&
               197.D+0,228.D+0,258.D+0,289.D+0, 319.D+0, 350.D+0, 365.D+0/

    ! OWS Papa CLOUD CLIMATOLOGY
    !offdata coktas /6.29D+0,6.26D+0,6.31D+0,6.31D+0,6.32D+0,6.70D+0,&
    !           7.12D+0,7.26D+0,6.93D+0,6.25D+0,6.19D+0,6.23D+0,6.31D+0,6.29D+0/
    ! BATS CLOUD CLIMATOLOGY
    data coktas /4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,&
                 4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0,4.00D+0/
    ! BATS COADS CLOUD CLIMATOLOGY
    !off data coktas /5.25D+0,5.36D+0,5.38D+0,5.18D+0,4.84D+0,4.65D+0,&
    !                 4.58D+0,3.87D+0,3.81D+0,4.15D+0,4.70D+0,4.81D+0,5.14D+0,5.25D+0/

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    pi=acos(-1.D+0)

    dec=-0.406D+0*cos(2.D+0*pi*day/365.D+0) ! in radians and from Evans and Parlsow, 1985
    alatr=alatd*2.D+0*pi/360.D+0
!   correction of solar constant due to the ellipse if the earth orbit
    enot=1.D+0+(0.033D+0*cos(2.D+0*pi*day/365.D+0))              ! E_0 from [Duffie and Beckman, 1980]

    th=(abs(thour-12.D+0)*15.D+0)*(2.D+0*pi/360.D+0)     ! hour angle [15.D+0 ==> degrees per hour]
    cosz=((sin(dec)*sin(alatr))+(cos(dec)*cos(alatr)*cos(th)))
    zenith=acos(cosz)  ! zenith ==> solar zenith angle in radians
    saltd=90.D+0-((360.D+0*zenith)/(2.D+0*pi))  ! solar altitude in degrees
    albedo=0.04D+0
    trans=1.0D+0

    do i=1,13
      if (day .ge. daym(i) .and. day .le. daym(i+1)) then
        ratio=(day-daym(i))/(daym(i+1)-daym(i))
        oktas=coktas(i)+(ratio*(coktas(i+1)-coktas(i)))
      endif
    enddo
    cloud=oktas/8.D+0
!   cloud corr from Smith and Dobson, [Ecological Modeling, 14, 1-19, 1981]
    if (zenith .le. 1.55D+0) then
      cloudy=0.0375D+0+(cosz*exp(-0.24D+0/cosz)*((cloud*exp(-0.07D+0/cosz))+1.D+0-cloud))
    else
      cloudy=0.0375D+0
    endif

!   calculate the light value only during the day [night = zenith .ge. pi/2.D+0]
    if (zenith .lt. pi/2.D+0) then
      par0=solar*trans*enot*cosz*0.43D+0*(1.D+0-albedo)*cloudy
    else
      par0=0.D+0
    endif

!   calculate aJ by integrating (simple trapezoidal) PvsI terms over the aMLD
    aJ=0.D+0
    if (par0 .gt. 0.D+0) then
      delz=aMLD/100.D+0
      do iz=0,100
        z=float(iz)*delz
        parz=par0*exp(-z*(akw+(akc*phyto)))
        tmp=(Vp*alpha*parz)/(sqrt((Vp*Vp)+(alpha*alpha*parz*parz)))
        if (iz .eq. 0 .or. iz .eq. 100) then
          aJ=aJ+(0.5D+0*tmp)
        elseif (iz .gt. 0 .and. iz .lt. 100) then
          aJ=aJ+tmp
        endif
      enddo
      aJ=aJ*delz/aMLD
      
    endif
    
end subroutine JQforce

subroutine mldforce(day, h, aMLD)
    ! subroutine to determine the aMLD and h(t) terms
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

use fasham_variables_module
use GP_variables_module

    implicit none
    
    real (kind=8), intent(in) :: day
    real (kind=8), intent(out) :: h, aMLD
    integer (kind=4), parameter :: n=14
    real (kind=8) :: daym(14),cmld(14)
    integer (kind=4) nloop,iloop,i
    

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    data daym /0.D+0,16.D+0,46.D+0,75.D+0,105.D+0,136.D+0,&
               166.D+0,197.D+0,228.D+0,258.D+0,289.D+0,319.D+0,350.D+0,365.D+0/

    ! forcing data for SUPER [a.k.a. OWS Papa] site
    ! forcing data for OWS Papa site [data digitized from Matear, 1995]

    data cmld /87.0D+0,92.5D+0,100.0D+0,93.5D+0,67.0D+0,45.0D+0,&
               30.0D+0,19.0D+0,22.5D+0,30.0D+0,45.0D+0,65.0D+0,82.0D+0,87.0D+0/

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! calculate the aMLD and h values
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    ! find a [nloop+1] day moving average for h
    nloop=30     ! make the number EVEN to keep it symmetrical
    th=0.D+0
    do iloop=1,nloop+1
      tday=day-(nloop/2.D+0)+iloop
      if (tday .lt. 0.D+0) tday=tday+365.D+0
      if (tday .ge. 365.D+0) tday=tday-365.D+0
      do i=1,n-1
        if (tday .ge. daym(i) .and. tday .lt. daym(i+1)) then
          th=th+((cmld(i+1)-cmld(i))/(daym(i+1)-daym(i)))
        endif
      enddo
    enddo
    h=th/(nloop+1)

    do i=1,n-1
      if (day .ge. daym(i) .and. day .le. daym(i+1)) then
        ratio=(day-daym(i))/(daym(i+1)-daym(i))
        aMLD=cmld(i)+(ratio*(cmld(i+1)-cmld(i)))
      endif
    enddo
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    return
end subroutine mldforce
