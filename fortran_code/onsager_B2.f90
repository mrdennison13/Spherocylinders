module B2
  contains
  !************************************************************************
    subroutine onsager_B2(n_loop, S2_o, B2_o, error_o)

      use molecule_data
      use control_data
      use overlaps
      use functions

      implicit none

      real(8) :: r1(3),u1(3),r2(3),u2(3)
      real(8) :: sep(3)
      real(8) :: vexc, vexcsq, diff
      real(8) :: vexc_av, vexcsq_av, error, vexc_sub

      integer :: n_loop
      real(8) :: S2_o, B2_o, error_o

      integer  :: i, j, n_over
      logical  :: over

      !seed random number generator
      CALL RANDOM_SEED()

      ntrial = 100
      nbatch = n_loop/ntrial


      !set particle 1 to be in the centre of the box (0,0,0)
      r1(:)=0.0d0

      !Set the nematic order parameter S2
      S2 = S2_o

      !and find the corresponding alpha value (from the Onsager trial function)
      call find_alpha

      !set counters to zero
      vexc=0.0d0
      vexcsq=0.0d0

      do  i=1,nbatch

         !number of overlaps = 0
         n_over=0

         do j=1,ntrial
            call onsager(u1)            !generate random director
            call place(r1,r2,u2)        !place particle 2



            sep(:) = r2(:)-r1(:)            !separation between particle 1 and 2
            call overlap(sep, u1,u2, over)  !check for overlap

 
            !if over, add 1
            if (over) then
               n_over=n_over+1
            end if

         END DO

         vexc_sub = dble(n_over)/dble(ntrial)
         vexc=vexc+vexc_sub
         vexcsq=vexcsq+vexc_sub**2.0d0

      END DO



      vexc_av=vexc/dble(nbatch)
      vexcsq_av=vexcsq/dble(nbatch)
      diff = vexcsq_av-vexc_av**2.0d0

      if (diff .gt. 0.0d0) then
         error=sqrt(diff/dble(nbatch-1))
      else
         error=sqrt(-diff/dble(nbatch-1))
      end if

      B2_o = 4.0d0*vexc_av*rmax*rmax*rmax
      error_o = 4.0d0*error*rmax*rmax*rmax

      return

    end subroutine onsager_B2
end module B2
