module B5
contains
!************************************************************************
  subroutine onsager_B5(n_loop, S2_o, B5_o, error_o)

    use molecule_data
    use control_data
    use overlaps
    use functions


    implicit none


    ! real(8) :: L, D
    integer :: n_loop
    real(8) :: S2_o, B5_o, error_o

    REAL(8) :: r1(3),u1(3)
    REAL(8) :: r2(3),u2(3)
    REAL(8) :: r3(3),u3(3)
    REAL(8) :: r4(3),u4(3)
    REAL(8) :: r5(3),u5(3)
    REAL(8) :: sep(3)
    real(8) :: vexc, vexcsq, diff, n_over
    REAL(8) :: vexc_av, vexcsq_av, error, vexc_sub
    INTEGER ::  i, j
    LOGICAL ::    over
    LOGICAL ::    over23
    LOGICAL ::    over24
    LOGICAL ::    over25
    LOGICAL ::    over14
    LOGICAL ::    over15
    LOGICAL ::    over35

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

       n_over=0.0d0

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


          call place(r1,r3,u3)
          sep(:) = r3(:)-r1(:)
          call overlap(sep, u1,u3, over)
200       continue          
          if (.not. over) then
             call place(r1,r3,u3)
             sep(:) = r3(:)-r1(:)
             call overlap(sep, u1,u3, over)
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


          call place(r4,r5,u5)
          sep(:) = r5(:)-r4(:)
          call overlap(sep, u4,u5, over)
400       continue          
          if (.not. over) then
             call place(r4,r5,u5)
             sep(:) = r5(:)-r4(:)
             call overlap(sep, u4,u5, over)
             goto 400
          end if

          sep(:) = r3(:)-r2(:)
          call overlap(sep, u2,u3, over23)
          sep(:) = r4(:)-r2(:)
          call overlap(sep, u2,u4, over24)
          sep(:) = r5(:)-r2(:)
          call overlap(sep, u2,u5, over25)

          sep(:) = r4(:)-r1(:)
          call overlap(sep, u1,u4, over14)
          sep(:) = r5(:)-r1(:)
          call overlap(sep, u1,u5, over15)

          sep(:) = r5(:)-r3(:)
          call overlap(sep, u3,u5, over35)



          if (over25 .and..not. over14 .and..not. over15 .and..not. over23 .and..not. over24 .and..not. over35) then
             n_over=n_over-12.0d0


          elseif(over24 .and..not. over14 .and. over15.and..not. over23 .and..not. over25 .and..not. over35) then
             n_over=n_over+10.0d0


          elseif(.not. over24 .and..not. over15 .and..not. over23 .and.over14.and.over35.and.over25) then
             n_over=n_over+(60.0d0/7.0d0)              
          elseif(.not. over25 .and..not. over23 .and..not. over14 .and.over24.and.over15.and.over35) then
             n_over=n_over+(60.0d0/7.0d0)
          elseif(.not. over23 .and..not. over14 .and..not. over15 .and.over24.and.over25.and.over35) then
             n_over=n_over+(60.0d0/7.0d0)
          elseif(.not. over14 .and..not. over23 .and..not. over35 .and.over15.and.over24.and.over25) then
             n_over=n_over+(60.0d0/7.0d0)
          elseif(.not. over24 .and..not. over14 .and..not. over35 .and.over23.and.over25.and.over15) then
             n_over=n_over+(60.0d0/7.0d0)                
          elseif(.not. over24 .and..not. over15 .and..not. over35 .and.over14.and.over23.and.over25) then
             n_over=n_over+(60.0d0/7.0d0)
          elseif(.not. over14 .and..not. over25 .and..not. over35 .and.over15.and.over24.and.over23) then
             n_over=n_over+(60.0d0/7.0d0)


          elseif(.not. over23 .and..not. over14 .and.over15.and.over24.and.over25.and.over35) then
             n_over=n_over+7.5d0
          elseif(.not. over23 .and..not. over15 .and.over14.and.over24.and.over25.and.over35) then
             n_over=n_over+7.5d0
          elseif(.not. over24 .and..not. over15 .and.over35.and.over25.and.over14.and.over23) then
             n_over=n_over+7.5d0
          elseif(.not. over14 .and..not. over25 .and.over35.and.over23.and.over24.and.over15) then
             n_over=n_over+7.5d0
          elseif(.not. over14 .and..not. over35 .and.over15.and.over24.and.over23.and.over25) then
             n_over=n_over+7.5d0
          elseif(.not. over24 .and..not. over35 .and.over14.and.over15.and.over25.and.over23) then
             n_over=n_over+7.5d0


          elseif(over24.and.over35.and.over14.and.over15.and.over25.and.over23) then
             n_over=n_over-6.0d0

          end if


       END DO

       vexc_sub = n_over/dble(ntrial)

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

    !S2_o(p) = S2
    B5_o = -(30.0d0/56.0d0)*vexc_av
    error_o = (30.0d0/56.0d0)*error

    ! end do


  end subroutine  onsager_B5
end module B5
