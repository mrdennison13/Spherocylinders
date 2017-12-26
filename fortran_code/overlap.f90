module overlaps

  interface

     subroutine sub_interface(r, u1, u2, over)
       real(8) :: r(3), u1(3), u2(3)
       logical :: over
     end subroutine sub_interface
  end interface


  procedure(sub_interface),pointer :: overlap
contains





  !  interface
  !     function fun(n,x) result(f)
  !       integer, intent(in) :: n
  !       double precision, intent(in) :: x(n)
  !       double precision :: f(n)
  !     end function fun
  !     subroutine sub(r, u1, u2, over)
  !       real(8), intent(in) :: r(3), u1(3), u2(3)
  !       logical :: over
  !     end subroutine sub
  !  end interface




  subroutine overlap_spherocylinder(r12, u1, u2, over)

    use molecule_data
    use control_data

    implicit none

    real(8) :: r12(3)
    real(8) :: u1(3), u2(3)
    logical :: over

    real(8) :: RUI, RUJ
    real(8) :: UIJ
    real(8) :: SINSQ
    real(8) :: CI, CJ
    real(8) :: DI, DJ
    real(8) :: RIJSQ

    over = .false.

    RUI = sum(r12(:)*u1(:))
    RUJ = sum(r12(:)*u2(:))

    UIJ = sum(u1(:)*u2(:))

    SINSQ=1.0d0 - UIJ*UIJ

    if (SINSQ .lt. 1e-14) then
       CI = -RUI*0.5d0
       CJ = RUJ*0.5d0
    else
       CI = (-RUI+UIJ*RUJ)/SINSQ
       CJ = (RUJ-UIJ*RUI)/SINSQ
    endif

    if (DMAX1(DABS(CI),DABS(CJ)).lt.XLD2) goto 641

    if (DABS(CI).ge.DABS(CJ)) then
       CI=DSIGN(XLD2,CI)
       CJ=RUJ+CI*UIJ
    else
       CJ=DSIGN(XLD2,CJ)
       CI=CJ*UIJ-RUI
    endif

    if (DABS(CI).GT.XLD2) CI=DSIGN(XLD2,CI)
    if (DABS(CJ).GT.XLD2) CJ=DSIGN(XLD2,CJ)

641 continue

    DI=2.0d0*RUI+CI-CJ*UIJ
    DJ=-2.0d0*RUJ+CJ-CI*UIJ
    RIJSQ=sum(r12(:)**2)
    if(RIJSQ+CI*DI+CJ*DJ < 1.0d0) then
       over = .true.
    end if
    
    return

  end subroutine overlap_spherocylinder









  SUBROUTINE overlap_cut_sphere(rij, ui, uj, over)

    use molecule_data
    use control_data

    implicit none

    real(8) :: rij(3)
    real(8) :: ui(3), uj(3)
    logical :: over


    INTEGER :: i
    REAL(8):: ul(3),rijsq, xnorm, uij, rui, ruj, rijul
    REAL(8):: ruisq, rujsq, rumaxsq, csij, snij, rmin, xhuij
    REAL(8):: xla, xlasq, xmu, xmusq, di, dj, xh, dl, dlsq, rc, ds, dssq
    REAL(8):: dll, snijsq


    over=.FALSE.


    rijsq=sum(rij(:)*rij(:))

    ! 	sphere test

    IF( rijsq .LE. LDsq ) THEN

       over=.TRUE.
       RETURN

    END IF

    IF( rijsq .GT. dsq) RETURN
    !
    ! 	rim-rim test
    !

    ! 	cross product calculations

    ul(1)=ui(2)*uj(3)-ui(3)*uj(2)
    ul(2)=ui(3)*uj(1)-ui(1)*uj(3)
    ul(3)=ui(1)*uj(2)-ui(2)*uj(1)

    xnorm=0.0d0
    DO   i=1,3
       xnorm=xnorm+ul(i)*ul(i)
    END DO
    xnorm=1.0d0/sqrt(xnorm)
    DO  i=1,3
       ul(i)=ul(i)*xnorm
    END DO

    !       inproduct calculation

    uij=0.0d0
    rui=0.0d0
    ruj=0.0d0
    rijul=0.0d0

    DO  i=1,3
       uij=uij+ui(i)*uj(i)
       rui=rui+rij(i)*ui(i)
       ruj=ruj+rij(i)*uj(i)
       rijul=rijul+rij(i)*ul(i)
    END DO

    ruisq=rui*rui
    rujsq=ruj*ruj
    rumaxsq=MAX(ruisq,rujsq) 
    IF(rumaxsq .LT. (rijsq*rtsq) ) THEN

       over=.TRUE.
       RETURN

    END IF

    csij = abs(uij)
    snijsq = 1.0d0-csij*csij
    snij = sqrt(snijsq)

    ! too close or too far test

    IF(csij.LT.rt)THEN
       IF(rijsq.LE.rhsq) THEN

          over=.TRUE.
          RETURN
       END IF

       IF(rumaxsq .GT. rhsq) RETURN

    ELSE

       rmin=LDh*(1.0d0+csij)+rf*snij
       rmin=rmin*rmin
       IF(rijsq.LE.rmin) THEN

          over=.TRUE.
          RETURN
       END IF

       IF(rumaxsq .GT. rmin) RETURN

    END IF

    ! faces overlap

    rijul=abs(rijul)
    xhuij=LDh*uij

    !	check overlap upper face particle j with lower face particle i

    xla=ruj+xhuij+LDh
    xlasq=xla*xla/snijsq
    if(xlasq.gt.rfsq)goto 200
    xmu=rui+LDh+xhuij
    xmusq=xmu*xmu/snijsq

    if(xmusq.gt.rfsq)goto 200
    di=sqrt(rfsq-xlasq)
    dj=sqrt(rfsq-xmusq)

    IF((di+dj).GT.rijul)THEN

       over=.TRUE.
       RETURN
    END IF

200 CONTINUE

    ! 	check overlap lower face particle j with lower face particle i

    xla=ruj+xhuij-LDh
    xlasq=xla*xla/snijsq

    IF(xlasq.GT.rfsq)GOTO 300
    xmu=rui+LDh-xhuij
    xmusq=xmu*xmu/snijsq

    IF(xmusq.GT.rfsq)GOTO 300
    di=sqrt(rfsq-xlasq)
    dj=sqrt(rfsq-xmusq)

    IF((di+dj).GT.rijul)THEN

       over=.TRUE.
       RETURN
    END IF

300 CONTINUE

    ! 	check overlap upper face particle j with upper face particle i

    xla=ruj-xhuij+LDh
    xlasq=xla*xla/snijsq

    IF (xlasq.GT.rfsq) GOTO 400
    xmu=rui-LDh+xhuij
    xmusq=xmu*xmu/snijsq

    IF (xmusq.GT.rfsq)GOTO 400
    di=sqrt(rfsq-xlasq)
    dj=sqrt(rfsq-xmusq)

    IF ((di+dj).GT.rijul) THEN

       over=.TRUE.
       RETURN
    END IF
400 CONTINUE

    !	check overlap lower face particle j with upper face particle i

    xla=ruj-xhuij-LDh
    xlasq=xla*xla/snijsq

    IF (xlasq.GT.rfsq) GOTO 600
    xmu=rui-LDh-xhuij
    xmusq=xmu*xmu/snijsq

    IF (xmusq.GT.rfsq) GOTO 600
    di=sqrt(rfsq-xlasq)
    dj=sqrt(rfsq-xmusq)

    IF ((di+dj).GT.rijul) THEN

       over=.TRUE.
       RETURN
    end if

600 CONTINUE

    !	overlap between face of one particle and rim of other particle

    IF (ruj.LE.0.0d0)THEN
       xh=LDh
    ELSE
       xh=-LDh
    END IF
    dl=abs(ruj+xh)

    IF (dl.LE.rad)THEN
       dlsq=dl*dl
       rc=sqrt(radsq-dlsq)
       dssq=rijsq-rujsq
       ds=sqrt(dssq)
       IF (ds.LE.(rf+rc))THEN
          IF (ds.LE.rf)THEN
             IF ((dl*csij).LE.LDh)THEN
                over=.TRUE.
                RETURN
             END IF
          ELSE 
             dll=rui+xh*uij-rf*(rui-ruj*uij)/ds
             dll=abs(dll)
             IF (dll.LE.LDh)THEN

                over=.TRUE.
                RETURN
             END IF
          END IF
       END IF
    END IF

    IF (rui.LE.0.0d0)THEN
       xh=-LDh
    ELSE
       xh=LDh
    END IF

    dl=abs(rui-xh)
    IF (dl.LE.rad)THEN
       dlsq=dl*dl
       rc=sqrt(radsq-dlsq)
       dssq=rijsq-ruisq
       ds=sqrt(dssq)
       IF (ds.LE.(rf+rc))THEN
          IF (ds.LE.rf)THEN
             IF ((dl*csij).LE.LDh)THEN
                over=.TRUE.
                RETURN
             END IF
          ELSE
             dll=-ruj+xh*uij-rf*(-ruj+rui*uij)/ds
             dll=abs(dll)
             IF (dll.LE.LDh)THEN

                over=.TRUE.
                RETURN
             END IF
          END IF
       END IF
    END IF

    RETURN 

  END SUBROUTINE overlap_cut_sphere









  subroutine overlap_sphere(r12, u1, u2, over)

    use molecule_data
    use control_data

    implicit none

    real(8) :: r12(3)
    real(8) :: u1(3), u2(3)
    logical :: over


    over = .false.

    if(sum(r12(:)*r12(:)) .le. rmaxsq) then
       over = .true.
    end if

    return

  end subroutine overlap_sphere











  subroutine overlap_oblate_spherocylinder(r12,u1,u2,over)
    use molecule_data
    use control_data

    implicit none

    REAL(8) :: r12(3),u1(3),u2(3)
    logical :: over

    REAL(8) :: rsq,rx,ry,rz,u(3),r12dotu
    real(8) :: dot,drperp_sqr,dr,norm,invnorm
    real(8) :: res

    INTEGER :: success,two_pass

    over = .false.

    rx = r12(1)
    ry = r12(2)
    rz = r12(3)

    rsq=rx*rx + ry*ry + rz*rz

    if(rsq .gt. Dsq) then
       return
    end if

    dot=u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3)

    if (1.0d0-dot .lt. 1.0e-11) then

       u(:) = 0.5d0*(u1(:)+u2(:))

       dot=rx*u(1)+ry*u(2)+rz*u(3)

       drperp_sqr = rsq - dot*dot

       if(drperp_sqr .lt. 4.0d0*RADSQ) then

          if (dot*dot .lt. LDsq) then

             over = .true.

             return

          else

             return

          end if

       else

          dr=sqrt(drperp_sqr)-2.0d0*RAD

          if (dot*dot+dr*dr .lt. LDsq) then 

             over = .true.

             return

          else

             return

          end if


       end if

       print*, 'here'

    end if

    u(1)=u1(2)*u2(3) - u1(3)*u2(2)
    u(2)=u1(3)*u2(1) - u1(1)*u2(3)
    u(3)=u1(1)*u2(2) - u1(2)*u2(1)

    norm=sqrt(1.0d0-dot*dot)

    invnorm = 1.0d0/norm

    u(:) = u(:)*invnorm

    r12dotu=rx*u(1)+ry*u(2)+rz*u(3)

    res=try_analytic(-r12(:),u1,u2,rsq,u,dot,norm,invnorm,r12dotu,success)

    if (success .ne. 0) then

       if (res .lt. LDsq) then

          over = .true.
          return

       else

          return

       end if
    end if

    u(:) = -u(:)

    res=try_analytic(r12(:),u2,u1,rsq,u,dot,norm,invnorm,r12dotu,success)

    if (success .ne. 0) then

       if (res .lt. LDsq) then

          over = .true.
          return

       else

          return

       end if
    end if


    res=calc_numeric(-r12(:),u1,u2,LD)

    if (res .lt. LDsq) then

       over = .true.
       return

    else

       return

    end if

    return

  end subroutine overlap_oblate_spherocylinder



  ! analytic part of the overlap check, see JChemPhys_129 p. 214706
  function try_analytic(rr,u1,u2,rsqr,u,dot,norm, invnorm,r12dotu,success)
    use molecule_data
    use control_data

    implicit none

    REAL(8) :: try_analytic
    REAL(8) :: rr(3),u1(3),u2(3),u(3),rsqr
    logical :: over

    REAL(8) :: rsq,rx,ry,rz,r12dotu
    real(8) :: dot,drperp_sqr,dr,norm,invnorm
    real(8) :: res,rp_sqr,qf,qt,qp,qm!,rp_sqr
    real(8) :: r12dot1,test1,bdot1,xnp,xnm,r12dotb

    INTEGER :: success

    !print*, 'func'

    success=1

    r12dot1=rr(1)*u1(1) + rr(2)*u1(2) + rr(3)*u1(3)!dcx*dix+dcy*diy+dcz*diz

    test1 = -(r12dot1*invnorm*invR)

    bdot1=norm

    xnp=0.0d0
    xnm=0.0d0

    !r12dotb=(r12dot1-dot*(djx*dcx+djy*dcy+djz*dcz))*invnorm
    r12dotb=(r12dot1-dot*(u2(1)*rr(1)+u2(2)*rr(2)+u2(3)*rr(3)))*invnorm

    if(abs(test1) .lt. 1.0d0) then
       qf=rsqr+RADSQ+2.0d0*RAD*test1*(r12dotb)
       qt=2.0d0*RAD*sqrt(1.0d0-test1*test1)*r12dotu

       qp=qf+qt
       qm=qf-qt

       if(qm.lt.0.0d0 .or. qp.lt.0.0d0)then

          !  erador_exit(0,"OOOOPS!!!")
          print*, 'erador'
          stop

       end if
       if(qp.lt.RADSQ .or. qm.lt.RADSQ) then

          try_analytic = 0.0d0

          return 

       end if

    else

       if(LD .eq. 0.0) then

          try_analytic = 0.01d0

          return

       end if

       xnp=(r12dot1+RAD*bdot1)
       xnm=(r12dot1-RAD*bdot1)

       !  if(optim_for_overlap) then
       if(abs(xnp) .gt. LD+1.0e-4 .and. abs(xnm).gt.LD+1.0e-4) then

          try_analytic = (LD+0.01d0)*(LD+0.01d0)

          return 
       end if
       !  end if

       if(abs(xnp) .lt. abs(xnm)) then
          rp_sqr=rsqr+RADSQ+xnp*xnp-2.0*xnp*r12dot1+2.0*RAD*(r12dotb-xnp*bdot1)
          if(rp_sqr .lt. RADSQ) then

             try_analytic = xnp*xnp

             return 
          end if

       else
          rp_sqr=rsqr+RADSQ+xnm*xnm-2.0*xnm*r12dot1-2.0*RAD*(r12dotb-xnm*bdot1)
          if(rp_sqr .lt. RADSQ) then

             try_analytic = xnm*xnm

             return 
          end if
       end if

    end if

    success=0
    return! 0.0/0.0

  end function try_analytic







  function calc_numeric(rr,u1,u2,dist_max)
    use molecule_data
    use control_data

    implicit none

    real(8) :: calc_numeric
    real(8) :: rr(3),u1(3),u2(3),ri(3), rj(3)
    real(8) :: dist_max,res,oldres
    integer :: i	
    res = dist_max*dist_max

    ri(1)=rr(1)!dcx
    ri(2)=rr(2)!dcy
    ri(3)=rr(3)!dcz

    do i = 1,1000000

       oldres=res
       !   if(i%2){
       if (mod(i,2) .ne. 0) then
          ! res=update_point(-dcx,-dcy,-dcz,dix,diy,diz,ri,rj)
          res=update_point(-rr,u1,ri,rj)
       else
          ! res=update_point(dcx,dcy,dcz,djx,djy,djz,rj,ri)
          res=update_point(rr,u2,rj,ri)
       end if
       if(res .lt. dist_max*dist_max) then
          calc_numeric = res
          return 
       end if
       if(abs(oldres-res) .lt. 1.0e-9*res) then
          calc_numeric = res
          return 
       end if
       ! erador_exit("calc_numeric()","Could not find the minimum distance in 1000000 iterations!");
    end do
    calc_numeric = res
    return
  end function calc_numeric



  function update_point(rr,u,r_in,r_out)
    use molecule_data
    use control_data

    implicit none

    real(8) :: update_point
    real(8) :: rr(3),u(3),r_in(3),r_out(3)
    real(8) :: rperp(3), dot,rsqr,fact,tdr(3)

    dot=r_in(1)*u(1)+r_in(2)*u(2)+r_in(3)*u(3)

    rperp(1)=r_in(1)-dot*u(1)
    rperp(2)=r_in(2)-dot*u(2)
    rperp(3)=r_in(3)-dot*u(3)

    rsqr=rperp(1)*rperp(1)+rperp(2)*rperp(2)+rperp(3)*rperp(3)

    if(rsqr .lt. RADSQ) then
       r_out(1)=rr(1)+rperp(1)
       r_out(2)=rr(2)+rperp(2)
       r_out(3)=rr(3)+rperp(3)

       update_point = dot*dot

       return 

    else

       fact=RAD/sqrt(rsqr)

       r_out(1)=rr(1)+fact*rperp(1)
       r_out(2)=rr(2)+fact*rperp(2)
       r_out(3)=rr(3)+fact*rperp(3)

       tdr(:) = r_in(:)

       tdr(:) =  tdr(:) - fact*rperp(:)

       update_point = tdr(1)*tdr(1) + tdr(2)*tdr(2) + tdr(3)*tdr(3)

       return 

    end if

  end function update_point
























    !! ###################################################################
  subroutine overlap_spheroid(ra,ua,ub,over)
    use molecule_data
    use control_data

    implicit none



    !    ** No include files for this routine **
    !    ===================================================================

    !    ===================================================================
    !    |  Routine to calculate the reduced Perram-Wertheim f function    |
    !    |  we specialize to the case where the r vector lies along z      |
    !    |  We return f = F/rij**2                                         |
    !    ===================================================================

    !    ** Arguments **
    REAL(8) ::  a, b, n
    REAL(8) :: rxi,ryi,rzi,uxi,uyi,uzi
    REAL(8) :: rxj,ryj,rzj,uxj,uyj,uzj
    logical :: over

    !    ** Local variables **
    REAL(8) :: ui, uj, uiuj
    REAL(8) :: rxij,ryij,rzij,r2,rr
    REAL(8) :: l
    REAL(8) :: lmin, lmax, lmean, ldiff, lv, lw, lx
    REAL(8) :: f
    REAL(8) :: fv, fw, fx, ldif, p, q, r, dd
    REAL(8) :: tol, tol2
    REAL(8) :: aalpha, beta, gamma
    REAL(8) :: x, y, denom, bsq

    real(8) :: ra(3),ua(3),rb(3),ub(3)

    INTEGER:: i

    !    ** Parameters **
    REAL(8) :: golden, eps, toler

    golden = 0.3819966d0 
    eps = 0.0000000001d0
    toler = 0.000000000001d0



    over =.false.

    rxij = ra(1)!-rb(1)
    ryij = ra(2)!-rb(2)
    rzij = ra(3)!-rb(3)

    r2 = rxij**2.0d0 + ryij**2.0d0 + rzij**2.0d0
    rr=sqrt(r2)

    if (r2.gt.rmaxsq) then 

       return

    end if

    uxi=ua(1)
    uyi=ua(2)
    uzi=ua(3)

    uxj=ub(1)
    uyj=ub(2)
    uzj=ub(3)

    uiuj = uxi*uxj + uyi*uyj + uzi*uzj

    ui = uxi*rxij + uyi*ryij + uzi*rzij

    uj = uxj*rxij + uyj*ryij + uzj*rzij

    ui=ui/rr
    uj=uj/rr

    l     = golden

    n=0.0d0

    i =0

    !    ** initial function evaluation **
    x     = (1.0d0-l)*(1.0d0-e**2.0d0)
    y     =    l *(1.0d0-e**2.0d0)
    denom = (1.0d0-x)*(1.0d0-y)-x*y*((uiuj)**2.0d0)
    aalpha = x*(1.0d0-y)
    beta  = y*(1.0d0-x)
    gamma = x*y*uiuj  
    f = -l*(1.0d0-l)*(1.0d0+(aalpha*(ui**2.0d0)+beta*(uj**2.0d0)+2.0d0*gamma*ui*uj)/denom)

    lmin = 0.0d0
    lmax = 1.0d0
    lv   = l
    lw   = l
    lx   = l
    fv   = f
    fw   = f
    fx   = f
    ldif   = 0.0d0

    !    ** Begin iteration **

    !     ten: do i = 1,10000 







10  continue


    lmean = 0.5d0 * ( lmin + lmax )
    tol   = eps * abs ( lx ) + toler
    tol2  = 2.0d0 * tol
    ldiff = lmax - lmin

    !    ** This is the test to stop iteration **

    !        print*, abs ( lx - lmean ), ( tol2 - 0.5d0 * ldiff )

    if ( abs ( lx - lmean ) .gt. ( tol2 - 0.5d0 * ldiff ) ) then

       p = 0.0d0
       q = 0.0d0
       r = 0.0d0

       !          print*, i

       if ( abs ( ldif ) .gt. tol ) then

          !          ** Parabolic fit **

          q = ( lx - lv ) * ( fx - fw )
          r = ( lx - lw ) * ( fx - fv )
          p = ( lx - lv ) * q - ( lx - lw ) * r
          q = 2.0d0 * ( q - r )

          if ( q .gt. 0.0d0 ) then
             p = -p
          else
             q = -q
          endif

          r  = ldif
          ldif = dd

       endif

       if ( ( abs ( p ) .lt. abs ( 0.5d0 * q * r ) ) .and.&
            &          ( p .gt. ( q * ( lmin - lx ) )       )   .and.&
            &         ( p .lt. ( q * ( lmax - lx ) )       ) ) then

          !          ** Parabolic interpolation step **
          dd = p / q
          l = lx + dd

          if ( ( ( l - lmin ) .lt. tol2 ) .or.&
               &             ( ( lmax - l ) .lt. tol2 ) ) then
             if ( lx .lt. lmean ) then
                dd = tol
             else
                dd = -tol
             endif
          endif

       else

          !          ** Golden section step **
          if ( lx .lt. lmean ) then
             ldif = lmax - lx
          else
             ldif = lmin - lx
          endif

          dd = golden * ldif

       endif

       if ( abs ( dd ) .ge. tol ) then
          l = lx + dd
       elseif ( dd .gt. 0.0d0 ) then
          l = lx + tol
       else
          l = lx - tol
       endif

       !       ** Function evaluation for this step **
       x     = (1.0d0-l)*(1.0d0-e**2.0d0)
       y     =    l *(1.0d0-e**2.0d0)
       denom = (1.0d0-x)*(1.0d0-y)-x*y*uiuj**2.0d0
       aalpha = x*(1.0d0-y)
       beta  = y*(1.0d0-x)
       gamma = x*y*uiuj  
       f = -l*(1.0d0-l)*(1.0d0+(aalpha*ui**2.0d0+beta*uj**2.0d0+2.0d0*gamma*ui*uj)/denom)

       if ( f .le. fx ) then

          if ( l .lt. lx ) then
             lmax = lx
          else
             lmin = lx
          endif

          lv = lw
          lw = lx
          lx = l
          fv = fw
          fw = fx
          fx = f

       else

          if ( l .lt. lx ) then
             lmin = l
          else
             lmax = l
          endif

          if ( ( f .le. fw ) .or. ( lw .eq. lx ) ) then

             lv = lw
             lw = l
             fv = fw
             fw = f

          elseif ( ( f .le. fv ) .or.&
               &                 ( lv .eq. lx ) .or. ( lv .eq. lw ) ) then

             lv = l
             fv = f

          endif

       endif

       i = i+1

       goto 10

    endif

    !    ** Minimum of -f has been reached **    

    bsq=0.25d0*bsp**2.0d0

    f=-f*r2/bsq

    if (f .lt. 1.0d0) then
       over = .true.

    END IF



    return

  end subroutine overlap_spheroid






end module overlaps
