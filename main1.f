      program example
          use dHwang1
          implicit none
      external fobj, jacobj
      double precision :: E, nu, sigy, N
      integer :: m
      double precision, allocatable, dimension(:) :: deps
      
      double precision :: t0, t1, dt
      double precision, dimension(:), allocatable :: y0, y1
      double precision, allocatable, dimension(:) :: ydot1, ydot2
      double precision, allocatable, dimension(:,:) :: pd

      integer :: neq, lrw, liw, mf, ml, mu, nrpd

      double precision, allocatable, dimension(:) :: rtol, atol
      integer :: itol, itask, istate, iopt, iout
      integer :: niters, num

      double precision, allocatable, dimension(:) :: rwork 
      integer, allocatable, dimension(:) :: iwork
      double precision, allocatable, dimension(:) :: rpar
      integer, allocatable, dimension(:) :: ipar

      double precision, allocatable, dimension(:) :: sig, eps_e, eps_p
      double precision, allocatable, dimension(:) :: Cpl, eps
      double precision :: deps_pl_p, deps_pl_p0, deps_pl_p1
      double precision :: epsp
      allocate(deps(6))
      allocate(sig(6))
      allocate(Cpl(6), eps(6)) 
      allocate(eps_e(6))
      allocate(eps_p(6)) 

      neq = 13
      nrpd = neq
      allocate(y0(neq))
      allocate(y1(neq))
      allocate(atol(neq))
      allocate(rtol(neq))
      
      allocate(ydot1(neq))
      allocate(ydot2(neq))

      allocate(pd(nrpd, neq))

      ydot1 = 0 
      ydot2 = 0
      pd = 0

      itol = 2
      rtol = 1.d-6

      atol(1:6) = 1.d-8
      atol(7:12) = 1.d-8
      atol(13) = 1.d-8

      itask = 1
      istate = 1
      iopt = 0

      !lrw = 373
      !liw = 43
      ml = 2
      mu = 5
      lrw = 22 + 11*neq + (3*ml + 2*mu)*neq
      liw = 30 + neq

      lrw = 10*lrw
      liw = 10*liw
      allocate(rwork(lrw))
      allocate(iwork(liw))

      E = 1
      nu = 0.3
      sigy = 0.2e-2*E
      N = 0.2
      m = 20
      
      num = 1
      deps = 0
      deps(1) = 1e-3

      allocate(rpar(11))
      allocate(ipar(2))
      rpar(1:6) = deps
      rpar(7) = E  
      rpar(8) = nu

      rpar(9) = sigy 
      rpar(10) = N

      ipar(1) = m

      t0 = 0
      t1 = 4
      dt = 1e-4
      niters = int((t1-t0)/dt)+1
      ipar(2) = niters

      y0 = 0
      y1 = 0
      t1 = dt
      mf = 22
      y0(1) = E*deps(1)*dt
      open(1, FILE="Hwang1d.dat")
      do iout = 1, niters
         call dvode(fobj, neq, y0, t0, t1, itol, rtol, atol, itask,
     1 istate, iopt, rwork, lrw, iwork, liw, jacobj, mf, rpar, ipar)
         write(6, 20) t0, y0(1), y0(2), y0(3), y0(13)
 20      format(' At t = ', d12.4,'   y=', 4d14.6)
         eps_e = y0(1:6)
         eps_p = y0(7:12)
         eps = eps_e + eps_p
         epsp = y0(13)
         call hooke1(sig, eps_e, E, nu)
         write(1, 10) eps(1), sig(1)/sigy
 10      format(e14.6, e14.6)

         if (istate .lt. 0) then
            write(6,90) istate
 90      format(///' Error halt: istate =', i3)
            stop
         end if
         t1 = t1 + dt
      end do
         write(6,60) iwork(11), iwork(12), iwork(13), iwork(19),
     1        iwork(20), iwork(21), iwork(22)
 60      format(/' NO. steps =', i4,'   No. f-s =',i4,
     1        '   No. J-s =',i4,'   No. lu-s =', i4/
     2        '  No. nonlinear iterations =', i4/
     3        '  No. nonlinear convergence failures =', i4/
     4        '  No. error test failure =',i4/)


C         call svode(fobj, neq, y0, t0, t1, itol, rtol, atol, itask,
C     1 istate, iopt, rwork, lrw, iwork, liw, jacobj, mf, rpar, ipar)
c      eps_e = deps*dt
c      epsp = dt
c      sig = eps_e
c      call hooke0(sig, E, nu)
c      deps_pl_p = deps_pl(sig, deps, epsp, sigy, E, N, m)
c      y0(1:6) = eps_e
c      call jac(neq, t0, y0, ml, mu, pd, nrpd, rpar, ipar)
c      !deps(4) = deps(4) - dt
c      eps_e(num) = eps_e(num) - dt
c      sig = eps_e
c      call hooke0(sig, E, nu)
c      deps_pl_p0 = deps_pl(sig, deps, epsp, sigy, E, N, m)
c      y1(num) = eps_e(num)
c      call f (neq, t0-dt, y1, ydot1, rpar, ipar) 
c      !deps(4) = deps(4) + 2*dt
c      eps_e(num) = eps_e(num) + 2*dt
c      sig = eps_e
c      call hooke0(sig, E, nu)
c      deps_pl_p1 = deps_pl(sig, deps, epsp, sigy, E, N, m)
c      y1(num) = eps_e(num)
c      call f (neq, t0+2*dt, y1, ydot2, rpar, ipar) 
c
c      write(6,*) (deps_pl_p1-deps_pl_p0)/(2*dt)
c      write(6,*) (ydot2 - ydot1)/(2*dt)
c
c      eps_e(num) = eps_e(num) - dt 
c      sig = eps_e
c      call hooke0(sig, E, nu)
c      call deps_pl_deps_e(Cpl, sig, deps, epsp, sigy, E, nu,
c     1 N, m)
c      write(6,*) Cpl
c      write(6,*) pd(num,:)
c
c       
c      call hooke1(sig, eps_e, E, nu)
c      deps_pl_p = deps_pl(sig, deps, epsp, sigy, E, N, m) 
c      epsp = epsp - dt 
c      deps_pl_p0 = deps_pl(sig, deps, epsp, sigy, E, N, m) 
c      epsp = epsp + 2*dt
c      deps_pl_p1 = deps_pl(sig, deps, epsp, sigy, E, N, m)
c      write(6,*) (deps_pl_p1-deps_pl_p0)/(2*dt)
c
c      epsp = epsp - dt
c      deps_pl_p = deps_pl_depsp(sig, deps, epsp, sigy, E, N, m)
c      write(6,*) deps_pl_p 
      end program example

      subroutine fobj(neq, t, y, ydot, rpar, ipar) 
          use dHwang1
      !    implicit none
      !integer, intent(in) :: neq 
      !double precision, dimension(:), intent(in) :: y
      !double precision, dimension(:) :: ydot
      !double precision, intent(in) :: t
      !double precision, dimension(:), intent(in) :: rpar
      !integer, dimension(:), intent(in) :: ipar
      integer neq, ipar
      double precision t, y, ydot, rpar
      dimension y(neq), ydot(neq), rpar(*), ipar(*)

      double precision, allocatable, dimension(:) :: eps_e
      double precision, allocatable, dimension(:) :: eps_p
      double precision :: epsp

      double precision, allocatable, dimension(:):: sig, deps
      double precision :: E, nu, sigy, N
      integer :: m, i, j

      double precision, allocatable, dimension(:) :: deps_e_dt
      double precision, allocatable, dimension(:) :: deps_p_dt
      double precision :: depsp_dt
      allocate(eps_e(6), eps_p(6), sig(6), deps(6))
      allocate(deps_e_dt(6), deps_p_dt(6))
      ! eps_e = 0
      ! write(6, *) "allocate", eps_e
      ! write(6, *) "rpar", rpar
      ! write(6, *) "ipar", ipar
      ! write(6, *) y
      do i = 1, 6
      eps_e(i) = y(i)
      end do
      do j = 1, 6
      eps_p(j) = y(6+j)
      end do
      epsp = y(13)
      
      deps = rpar(1:6)
      E = rpar(7)
      nu = rpar(8)

      sigy = rpar(9)
      N = rpar(10)

      m = ipar(1)

      call hooke1(sig, eps_e, E, nu)

      call hooke2(deps_e_dt, deps, E, nu)
      
! deps_p
      call deps_p(deps_p_dt, sig, deps, epsp, sigy, E, N, m) 
      ydot(7:12) = deps_p_dt
! deps_e
      !ydot(1:6) = deps - deps_p_dt
      ydot(1:6) = deps_e_dt
! depsp
      depsp_dt = deps_pl(sig, deps, epsp, sigy, E, N, m)
      ydot(13) = depsp_dt

      deallocate(eps_e, eps_p, sig, deps)
      deallocate(deps_e_dt, deps_p_dt)
      contains
          subroutine hooke2(eps, sig, E, nu)
              double precision, dimension(:), intent(in) :: sig
              double precision, intent(in) :: E, nu
              double precision, dimension(:) :: eps
              
              integer :: i,j
              double precision, dimension(:,:), allocatable :: S
              allocate(S(6,6))
              S = 0
              S(1,1) = 1
              S(2,2) = 1
              S(3,3) = 1
              S(1,2) = -nu
              S(1,3) = -nu
              S(2,3) = -nu
              S(4,4) = 2.0 + 2.0*nu
              S(5,5) = 2.0 + 2.0*nu
              S(6,6) = 2.0 + 2.0*nu

              do i = 1, 6
                do j = 1, i
                  S(i,j) = S(j,i)
                end do
              end do

              S = 1.0/E*S

              eps = matmul(S, sig)

              eps(4:6) = eps(4:6)/2.0
 
          end subroutine


 
      end subroutine 


      subroutine jacobj(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
          use dHwang1
      !    implicit none
      !double precision, dimension(:,:) :: pd
      !double precision, dimension(:) :: y
      !integer :: neq, ml, mu, nrpd
      !double precision :: t
      !double precision, dimension(:) :: rpar
      !integer, dimension(:) ::ipar
      integer neq, ml, mu, nrpd, ipar
      double precision t, y, pd, rpar
      dimension y(neq), pd(nrpd, neq), rpar(*), ipar(*)

      double precision, dimension(:), allocatable :: eps_e
      double precision, dimension(:), allocatable :: eps_p
      double precision :: epsp

      double precision, dimension(:), allocatable :: sig, deps
      double precision :: E, nu, sigy, N
      integer :: m

      double precision, allocatable, dimension(:, :) :: deps_e_dt_deps_e
      double precision, allocatable, dimension(:, :) :: deps_p_dt_deps_e
      double precision, allocatable, dimension(:) :: depsp_dt_deps_e

      double precision, allocatable, dimension(:) :: deps_e_dt_depsp
      double precision, allocatable, dimension(:) :: deps_p_dt_depsp
      double precision :: depsp_dt_depsp

      allocate(eps_e(6), eps_p(6), sig(6), deps(6))
      allocate(deps_e_dt_deps_e(6,6), deps_p_dt_deps_e(6,6))
      allocate(depsp_dt_deps_e(6), deps_e_dt_depsp(6),
     1  deps_p_dt_depsp(6))

      eps_e = y(1:6)
      eps_p = y(7:12)
      epsp = y(13)

      deps = rpar(1:6)
      E = rpar(7)
      nu = rpar(8)
      sigy = rpar(9)
      N = rpar(10)

      m = ipar(1)
      call hooke1(sig, eps_e, E, nu)
      
      call deps_p_deps_e(deps_p_dt_deps_e, sig, deps, epsp, sigy, E, nu,
     1 N, m) 
      call deps_p_depsp(deps_p_dt_depsp, sig, deps, epsp, sigy, E,
     1 N, m)
      call deps_pl_deps_e(depsp_dt_deps_e, sig, deps, epsp, sigy, E, nu,
     1 N, m) 
      depsp_dt_depsp = deps_pl_depsp(sig, deps, epsp, sigy, E,
     1 N, m)

      deps_e_dt_deps_e = -deps_p_dt_deps_e
      deps_e_dt_depsp = -deps_p_dt_depsp

      pd(1:6,1:6) = deps_e_dt_deps_e
      pd(1:6,13) = deps_e_dt_depsp
      pd(7:12,1:6) = deps_p_dt_deps_e
      pd(7:12,13) = deps_p_dt_depsp
      pd(13, 1:6) = depsp_dt_deps_e
      pd(13, 13) = depsp_dt_depsp

      deallocate(eps_e, eps_p, sig, deps)
      deallocate(deps_e_dt_deps_e, deps_e_dt_depsp)
      deallocate(deps_p_dt_deps_e, deps_p_dt_depsp, depsp_dt_deps_e)
 
      end subroutine

