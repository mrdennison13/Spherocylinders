subroutine calc_virial(level, shape, param_in, n_loop_in, S2_in, n,  B_out, error_out)

  use B2
  use B3
  use B4
  use B5
  use overlaps
  use molecule_data
  use control_data


  implicit none

  integer :: n
  real(8), dimension(n) :: param_in
  integer  :: level, n_loop_in
  real(8) :: S2_in, B_out, error_out!, stmt(10)
  character(40) :: shape!, stmt

  !f2py intent(in) :: level, shape, param_in, n_loop_in, n_points_in, S2_in
  !f2py intent(hide), depend(param_in) :: n = shape(param_in, 0)
  !f2py intent(out) B_out, error_out


  pi=4.0d0*atan(1.0d0)
  twopi=2.0d0*pi

  shape = trim(shape)

  if (shape .eq. 'sphere') then

     rmax   = param_in(1)
     rmaxsq = rmax**2.0d0
     vo = 4.0d0*pi*param_in(1)*param_in(1)*param_in(1)/3.0d0
     
     overlap => overlap_sphere

     
  elseif (shape .eq. 'spherocyl') then

     
     XLD2 = 0.5d0*param_in(1)/param_in(2)
     rmax   = param_in(1)+param_in(2)
     rmaxsq = rmax**2.0d0
     vo = pi*param_in(2)*param_in(2)*param_in(1)/4.0d0 + pi*param_in(2)*param_in(2)*param_in(2)/6.0d0
     
     overlap => overlap_spherocylinder



  elseif (shape .eq. 'oblate_spherocyl') then

     D = param_in(2)
     LD = param_in(1)

     dsq = D*D
     lDsq = LD*LD
     
     rad= 0.5d0*(1.0d0-LD)
     radsq= rad*rad
     invR = 1.0d0/rad
     
     sigma = D-LD

     rmax   = D
     rmaxsq = rmax**2.0d0
     vo = (pi/6.0d0)*LD*LD*LD + 0.125d0*pi*pi*sigma*LD*LD + 0.25d0*pi*sigma*sigma*LD
     
     overlap => overlap_oblate_spherocylinder

  elseif (shape .eq. 'cut_sphere') then

     LD = param_in(1)*param_in(2)
     d = param_in(2)
     LDsq=LD*LD
     dsq=d*d
     LDh=LD/2.0d0
     LDhsq=LDsq/4.0d0
     rad=d/2.0d0
     radsq=rad*rad
     rfsq=radsq-LDhsq !!*CORRECTED, HOPEFULLY
     rf=sqrt(rfsq)
     rt=LD/d
     rtsq=rt*rt
     rh=(LD+d)/2.0d0
     rhsq=rh*rh

     rmax   = D
     rmaxsq = rmax**2.0d0

     overlap => overlap_cut_sphere

     
  elseif (shape .eq. 'spheroid') then
     asp = param_in(1)
     bsp = param_in(2)
     e = asp/bsp

     rmax   = max(asp,bsp)
     rmaxsq = rmax**2.0d0

     overlap => overlap_spheroid

  endif

  if (level .eq. 2) then
     call onsager_B2(n_loop_in, S2_in, B_out, error_out)
  elseif (level .eq. 3) then
     call onsager_B3(n_loop_in, S2_in, B_out, error_out)
  elseif (level .eq. 4) then
     call onsager_B4(n_loop_in, S2_in, B_out, error_out)
  elseif (level .eq. 5) then
     call onsager_B5(n_loop_in, S2_in, B_out, error_out)
  end if

!  stmt = caux
  !stmt =trim(caux)

end subroutine calc_virial
