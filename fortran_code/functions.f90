module functions



contains


  function random_num()
    !function to generate random numbers  
    real(8) :: random_num

    CALL RANDOM_NUMBER(random_num)

    return 
  end function random_num



  subroutine find_alpha
    !subroutine to find alpha for a given nematic order parameter
    use control_data
    use molecule_data

    implicit none

    integer :: n
    real(8) :: x1,x2,f1,f2

    !find corresponding alpha value
    if (S2 .lt. 1e-6) then   !S2 = 0, alpha = 0
       alpha    = 0.0d0
       sh_alpha = 0.0d0
    elseif (abs(S2-1.0d0) .lt. 1e-6) then !S2 = 1, alpha -> infinity
       alpha = 1e14
       sh_alpha = 1e14
    else          !0 < S2 < 1 
       x1 = 10.0d0*S2
       n = 0
       DO
          f1 = 1.0d0 - S2 + (3.0d0 / (x1**2.0d0)) - (3.0d0 / (x1*tanh(x1)))
          f2 = (3.0d0/((x1**2.0d0)*tanh(x1)))-(6.0d0/X1**3.0d0)+(3.0d0/(x1*(TANH(x1)*COSH(X1))**2.0d0))
          x2 = x1 - f1/f2
          IF ((abs(x2 - x1)/x2) < 1e-6) THEN 
             !success, exit the loop
             alpha = x2
             goto 232
          END IF
          n = n + 1
          if(n.gt.10000000) then
             alpha = x2
             goto 232
             !give up trying to find alpha
             !print*, '# error finding alpha', S2, x1,x2
             !stop         
          end if
          x1 = x2
       END DO
232    continue
       alpha = abs(alpha)
       if (alpha .lt. 20.0d0) then
          sh_alpha=sinh(alpha)
       end if
    end if
    return
  end subroutine find_alpha




  subroutine onsager(u)
    !subroutine to generate a (random) unit vector
    use control_data

    implicit none

    real(8) :: u(3), x, ux, uy, uz, fac

    if(S2.lt.1e-6) then
       uz=random_num()
    elseif(abs(S2-1.0d0).lt.1e-6) then
       uz=1.0d0
       ux = 0.0d0
       uy = 0.0d0
       goto 234
    else
       if (alpha .lt. 20.0d0) then
          x=random_num()*sh_alpha
          uz=log(x+sqrt(1.0d0+x*x) )/alpha
       else
          x=random_num()
          uz = (log(x) + alpha)/alpha
       end if
    end if

    x = twopi*random_num()
    fac=sqrt( (1.0d0-uz*uz))
    ux = fac * cos(x)
    uy = fac * sin(x)

234 continue

    u(1)=ux
    u(2)=uy
    u(3)=uz

    return
  end subroutine onsager



  subroutine place(r1,r2,u2)
    !subroutine to place particle 2, relative to particle 1
    use molecule_data
    implicit none

    REAL(8) :: r1(3),r2(3),u2(3)
    integer :: i

    do  i=1,3
       r2(i)=r1(i)+rmax*(2.0d0*random_num()-1.0d0)
    end do
    !and generate its orientation
    call onsager(u2)
    return
  end subroutine place




end module functions
