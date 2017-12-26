!************************************************************************
module molecule_data
!************************************************************************
  implicit none


  !cut sphere / oblate sphero params
  real(8)    ::     LD
  real(8)    ::     D
  real(8)    ::     LDsq
  real(8)    ::     dsq
  real(8)    ::     LDh
  real(8)    ::     LDhsq
  real(8)    ::     rad
  real(8)    ::     radsq
  real(8)    ::    rfsq
  real(8)    ::    rf
  real(8)    ::    rt
  real(8)    ::    rtsq
  real(8)    ::      rh
  real(8)    ::      rhsq

  real(8)    :: e, asp, bsp                          !spheroids

  real(8)    :: Xii(3), Xjj(3),side_1(3), side_2(3)  !boards

  real(8)    ::   invR, sigma

  real(8)    :: XLD2                                 ! spherocylinder
  
  real(8) :: rmax, rmaxsq                     ! maximum separation (and squared)
  real(8) :: vo                               ! volume of particle


!************************************************************************
end module molecule_data
!************************************************************************



!************************************************************************
module control_data
!************************************************************************  

implicit none

integer          :: ntrial,nbatch                      !number of trials and number of batches
real(8)          :: S2, alpha, sh_alpha
real(8)          :: pi, twopi
!real(8)          :: caux(10)
!character(40)    :: caux
!************************************************************************
end module control_data
!************************************************************************
