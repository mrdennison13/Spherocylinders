module B3
contains
  subroutine onsager_B3(n_loop, S2_o, B3_o, error_o)

    use molecule_data
    use control_data
    use overlaps
    use functions


    implicit none

    REAL(8) :: r1(3),u1(3)
    REAL(8) :: r2(3),u2(3)
    REAL(8) :: r3(3),u3(3)
    real(8) :: sep(3)
    real(8) :: vexc, vexcsq, diff
    real(8) :: vexc_av, vexcsq_av, error, vexc_sub

    integer :: n_loop
    real(8) :: S2_o, B3_o, error_o

    integer  :: i, j, n_over
    logical  :: over

    !seed random number generator
    CALL RANDOM_SEED()


    ntrial = 100
    nbatch = n_loop/ntrial

    !set particle 1 to be in the centre of the box (0,0,0)
    r1(:)=0.0d0


    !loop through number of alpha values
    !do p = 1,n_points

    !Set the nematic order parameter S2
    S2 = S2_o
    !if (n_points .eq. 1) then
    !   S2 = 0.0d0
    !else
    !   S2 = dble(p-1)/dble(n_points-1)
    !end if

    !and find the corresponding alpha value (from the Onsager trial function)
    call find_alpha

    !set counters to zero
    vexc=0.0d0
    vexcsq=0.0d0

    do  i=1,nbatch

       n_over=0

       do j=1,ntrial

          call onsager(u1)

          call place(r1,r2,u2)
          sep(:) = r2(:)-r1(:)
          call overlap(sep, u1,u2, over)
100       continue
          if (.not. over) then
             call place(r1,r2,u2)
             sep(:) = r2(:)-r1(:)
             call overlap(sep, u1,u2, over)
             goto 100
          end if


          call place(r2,r3,u3)
          sep(:) = r3(:)-r2(:)
          call overlap(sep, u2,u3, over)
200       continue          
          if (.not. over) then
             call place(r2,r3,u3)
             sep(:) = r3(:)-r2(:)
             call overlap(sep, u2,u3, over)
             goto 200
          end if

          sep(:) = r3(:)-r1(:)
          call overlap(sep, u1,u3, over)

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

       error=sqrt((vexcsq_av-vexc_av**2.0d0)/dble(nbatch-1))

    else

       error=sqrt((vexc_av**2.0d0-vexcsq_av)/dble(nbatch-1))
       !   stop

    end if

    ! S2_o(p) = S2
    B3_o = (4.0d0/3.0d0)*vexc_av
    error_o = (4.0d0/3.0d0)*error

    ! end do



  end subroutine onsager_B3

end module B3
