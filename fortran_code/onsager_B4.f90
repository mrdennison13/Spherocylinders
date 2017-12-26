module B4
contains
  subroutine onsager_B4(n_loop, S2_o, B4_o, error_o)

    use molecule_data
    use control_data
    use overlaps
    use functions

    implicit none

    integer :: n_loop
    real(8) :: S2_o, B4_o, error_o

    REAL(8) :: r1(3),u1(3)
    REAL(8) :: r2(3),u2(3)
    REAL(8) :: r3(3),u3(3)
    REAL(8) :: r4(3),u4(3)
    REAL(8) :: sep(3)
    REAL(8) :: vexc, vexcsq, diff
    REAL(8) :: vexc_av, vexcsq_av, error, vexc_sub
    INTEGER ::  i, j, n_over
    LOGICAL ::    over
    LOGICAL ::    over13
    LOGICAL ::    over24


    !seed random number generator
    CALL RANDOM_SEED()

    ntrial = 100
    nbatch = n_loop/ntrial

    r1(:)=0.0d0

    !do p = 1,n_points

    S2 = S2_o
    !  if (n_points .eq. 1) then
    !     S2 = 0.0d0
    !  else
    !     S2 = dble(p-1)/dble(n_points-1)
    !  end if


    call find_alpha

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


          call place(r3,r4,u4)
          sep(:) = r4(:)-r3(:)
          call overlap(sep, u3,u4, over)
300       continue          
          if (.not. over) then
             call place(r3,r4,u4)
             sep(:) = r4(:)-r3(:)
             call overlap(sep, u3,u4, over)
             goto 300
          end if


          sep(:) = r4(:)-r1(:)
          call overlap(sep, u1,u4, over)

          if (over) then


             sep(:) = r3(:)-r1(:)
             call overlap(sep, u1,u3, over13)
             sep(:) = r4(:)-r2(:)
             call overlap(sep, u2,u4, over24)

             if (over13 .and. over24) then
                ! n_over=n_over-3
                n_over=n_over+2
             elseif (.not. over13 .and. .not. over24) then
                ! n_over=n_over+2
                n_over=n_over-3
             end if

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
       stop
    end if


    ! S2_o(p) = S2
    B4_o = vexc_av
    error_o = error

    !  write(*,'(2x,a1,3(f10.6,3x))') '#', S2,vexc_av,error

    !end do


  end subroutine onsager_B4
end module B4
