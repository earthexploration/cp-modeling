      module dHwang1
      implicit none
            
      contains

      subroutine hooke0(eps_e, E, nu) 
      double precision, dimension(:,:), allocatable :: C
      double precision, dimension(:) :: eps_e
      double precision :: E, nu

      integer :: i,j
      allocate(C(6,6))
      C = 0
      C(1,1) = 1 - nu
      C(2,2) = 1 - nu
      C(3,3) = 1 - nu
      C(1,2) = nu
      C(1,3) = nu
      C(2,3) = nu

      C(4,4) = (1-2*nu)/2.0
      C(5,5) = (1-2*nu)/2.0
      C(6,6) = (1-2*nu)/2.0

      do i = 1, 6
         do j = 1, i
            C(i,j) = C(j,i)
         end do
      end do

      C = E/(1+nu)/(1-2*nu)*C

      eps_e(4:6) = 2.0*eps_e(4:6)

      eps_e = matmul(C, eps_e)
      deallocate(C)
      end subroutine
       
      subroutine hooke1(sig, eps_e, E, nu) 
      double precision, dimension(:,:), allocatable :: C
      double precision, dimension(:), intent(in) :: eps_e
      double precision, dimension(:) :: sig
      double precision :: E, nu
      double precision, dimension(:), allocatable :: eps_e_v 
      integer :: i,j
      allocate(C(6,6), eps_e_v(6))
      C = 0
      C(1,1) = 1 - nu
      C(2,2) = 1 - nu
      C(3,3) = 1 - nu
      C(1,2) = nu
      C(1,3) = nu
      C(2,3) = nu

      C(4,4) = (1-2*nu)/2.0
      C(5,5) = (1-2*nu)/2.0
      C(6,6) = (1-2*nu)/2.0

      do i = 1, 6
         do j = 1, i
            C(i,j) = C(j,i)
         end do
      end do

      C = E/(1+nu)/(1-2*nu)*C

      eps_e_v = eps_e 
      eps_e_v(4:6) = 2.0*eps_e(4:6)

      sig = matmul(C, eps_e_v)
      deallocate(C, eps_e_v)
      end subroutine

      subroutine hooke(eps, E, nu) 
      double precision, dimension(6,6) :: C
      double precision, dimension(13) :: eps
      double precision :: E, nu

      integer :: i,j
      C = 0
      C(1,1) = 1 - nu
      C(2,2) = 1 - nu
      C(3,3) = 1 - nu
      C(1,2) = nu
      C(1,3) = nu
      C(2,3) = nu

      C(4,4) = (1-2*nu)/2.0
      C(5,5) = (1-2*nu)/2.0
      C(6,6) = (1-2*nu)/2.0

      do i = 1, 6
         do j = 1, i
            C(i,j) = C(j,i)
         end do
      end do

      C = E/(1+nu)/(1-2*nu)*C

      ! eps = [eps_e, eps_p, epsp]
      eps(7:12) = eps(7:12) + eps(1:6)

      eps(10:12) = 2.0*eps(10:12)
      eps(4:6) = 2.0*eps(4:6)

      ! eps = [sig, eps, epsp]
      eps(1:6) = matmul(C, eps(1:6))

      end subroutine
      
      function hard1(epsp, sigy, E, N)
      double precision :: epsp, sigy, E, N
      double precision :: hard1

      hard1 = (1.0 + E/sigy*epsp)**N
      
      return
      end

      function dhard1_depsp(epsp, sigy, E, N) result(dhard1_dp) 
      double precision :: epsp, sigy, E, N
      double precision :: dhard1_dp

      dhard1_dp = E/sigy*N*(1.0 + E/sigy*epsp)**(N-1)
      end

      subroutine dsig_dev_deps_e(Cklpq, E, nu) 
      double precision, intent(in) :: E, nu
      double precision, dimension(:,:) :: Cklpq
      double precision, dimension(:,:), allocatable :: C
      integer :: i, j
      allocate(C(6,6))
      C = 0
      C(1,1) = 1 - nu
      C(2,2) = 1 - nu
      C(3,3) = 1 - nu
      C(1,2) = nu
      C(1,3) = nu
      C(2,3) = nu

      C(4,4) = (1-2*nu)/2.0
      C(5,5) = (1-2*nu)/2.0
      C(6,6) = (1-2*nu)/2.0

      do i = 1, 6
         do j = 1, i
            C(i,j) = C(j,i)
         end do
      end do

      C = E/(1+nu)/(1-2*nu)*C

      Cklpq = 0
      Cklpq(1, :) = (C(1, :) + C(2, :) + C(3, :))/3.0
      Cklpq(2, :) = (C(1, :) + C(2, :) + C(3, :))/3.0
      Cklpq(3, :) = (C(1, :) + C(2, :) + C(3, :))/3.0

      Cklpq = C - Cklpq
      deallocate(C)
      end 

      function deps_pl(sig, deps, epsp, sigy, E, N, m) 
      double precision, dimension(:), intent(in) :: sig, deps
      double precision :: epsp, sigy, E, N
      integer :: m 
      double precision :: deps_pl
      double precision, dimension(:), allocatable :: sig_dev, deps_dev
      double precision :: sJ1, sJ2, deJ1, deJ2
      allocate(sig_dev(6), deps_dev(6))
      sig_dev = sig
      sJ1 = sig(1) + sig(2) + sig(3)
      sig_dev(1:3) = sig(1:3) - sJ1/3.0
      sJ2 = sqrt(3.0/2.0*dot_product(sig_dev, sig_dev))
      deps_dev = deps
      deJ1 = deps(1) + deps(2) + deps(3)
      deps_dev(1:3) = deps(1:3) - deJ1/3.0
      deJ2 = sqrt(2.0/3.0*dot_product(deps_dev, deps_dev))

      deps_pl = deJ2*(sJ2/sigy/hard1(epsp, sigy, E, N))**m
      deallocate(sig_dev, deps_dev)
      return
      end

      subroutine deps_pl_deps_e(Cpl, sig, deps, epsp, sigy, E, nu, N, 
     1 m) 
      double precision, dimension(:) :: Cpl
      double precision, dimension(:), intent(in) :: sig, deps
      double precision, dimension(:), allocatable :: sig_dev, deps_dev
      double precision :: epsp, sigy, E, nu, N
      integer :: m 
      double precision, dimension(:,:), allocatable :: Cklpq
      double precision :: sJ1, sJ2, deJ1, deJ2
      double precision :: deps_pl
      allocate(sig_dev(6), deps_dev(6))
      allocate(Cklpq(6,6)) 

      sig_dev = sig
      sJ1 = sig(1) + sig(2) + sig(3)
      sig_dev(1:3) = sig(1:3) - sJ1/3.0
      sJ2 = sqrt(3.0/2.0*dot_product(sig_dev, sig_dev))

      deps_dev = deps
      deJ1 = deps(1) + deps(2) + deps(3)
      deps_dev(1:3) = deps(1:3) - deJ1/3.0
      deJ2 = sqrt(2.0/3.0*dot_product(deps_dev, deps_dev))

      deps_pl = deJ2*(sJ2/sigy/hard1(epsp, sigy, E, N))**m

      Cpl = deps_pl*m/sJ2*3.0/2.0/sJ2

      call dsig_dev_deps_e(Cklpq, E, nu)
      Cpl = Cpl*matmul(transpose(Cklpq), sig_dev)

      deallocate(sig_dev, deps_dev)
      deallocate(Cklpq)
      end 

      function deps_pl_depsp(sig, deps, epsp, sigy, E, N, m)  
      double precision, dimension(:), intent(in) :: sig, deps
      double precision :: epsp, sigy, E, N
      integer :: m 
      double precision :: deps_pl
      double precision, dimension(:), allocatable :: sig_dev, deps_dev
      double precision :: sJ1, sJ2, deJ1, deJ2
      double precision :: h1, dh1, deps_pl_depsp
      allocate(sig_dev(6), deps_dev(6))
      sig_dev = sig
      sJ1 = sig(1) + sig(2) + sig(3)
      sig_dev(1:3) = sig(1:3) - sJ1/3.0
      sJ2 = sqrt(3.0/2.0*dot_product(sig_dev, sig_dev))
      deps_dev = deps
      deJ1 = deps(1) + deps(2) + deps(3)
      deps_dev(1:3) = deps(1:3) - deJ1/3.0
      deJ2 = sqrt(2.0/3.0*dot_product(deps_dev, deps_dev))

      deps_pl = deJ2*(sJ2/sigy/hard1(epsp, sigy, E, N))**m

      h1 = hard1(epsp, sigy, E, N)
      dh1 = dhard1_depsp(epsp, sigy, E, N)
      deps_pl_depsp = deps_pl*(-m)/h1*dh1
      
      deallocate(sig_dev, deps_dev)
      end 

      subroutine deps_p(deps_pe, sig, deps, epsp, sigy, E, N, m) 
      double precision, dimension(:), intent(in) :: sig, deps
      double precision :: epsp, sigy, E, N
      integer :: m
      double precision, dimension(:) :: deps_pe
      
      double precision, dimension(:), allocatable :: sig_dev, deps_dev
      double precision :: depsp 
      double precision :: sJ1, sJ2, deJ1, deJ2
      allocate(sig_dev(6), deps_dev(6))
      sig_dev = sig
      sJ1 = sig(1) + sig(2) + sig(3)
      sig_dev(1:3) = sig(1:3) - sJ1/3.0
      sJ2 = sqrt(3.0/2.0*dot_product(sig_dev, sig_dev))
      deps_dev = deps
      deJ1 = deps(1) + deps(2) + deps(3)
      deps_dev(1:3) = deps(1:3) - deJ1/3.0
      deJ2 = sqrt(2.0/3.0*dot_product(deps_dev, deps_dev))

      depsp = deJ2*(sJ2/sigy/hard1(epsp, sigy, E, N))**m

      if (sJ2 > 0) then
         deps_pe = 3.0*depsp/2.0/sJ2*sig_dev
      else
         deps_pe = 0
      end if 
      
      deallocate(sig_dev, deps_dev)
      end subroutine

      subroutine deps_p_deps_e(Cijpq, sig, deps, epsp, sigy, E, nu, N, 
     1 m)  
      double precision, dimension(:), intent(in) :: sig, deps
      double precision :: epsp, sigy, E, nu, N
      integer :: m, i, j
      double precision, dimension(:,:) :: Cijpq
      double precision, dimension(:,:), allocatable :: Cklpq
      double precision, dimension(:,:), allocatable :: Cijpq_es, 
     1 Cijpq_et, Cijpq_p 
      double precision, dimension(:), allocatable :: sig_dev, deps_dev,
     1 Cpl
      double precision, dimension(:), allocatable :: v1, v2
      double precision :: depsp 
      double precision :: sJ1, sJ2, deJ1, deJ2
      allocate(sig_dev(6), deps_dev(6), Cpl(6))
      allocate(Cklpq(6,6), Cijpq_es(6,6), Cijpq_et(6,6), Cijpq_p(6,6))
      allocate(v1(6), v2(6))
      sig_dev = sig
      sJ1 = sig(1) + sig(2) + sig(3)
      sig_dev(1:3) = sig(1:3) - sJ1/3.0
      sJ2 = sqrt(3.0/2.0*dot_product(sig_dev, sig_dev))
      deps_dev = deps
      deJ1 = deps(1) + deps(2) + deps(3)
      deps_dev(1:3) = deps(1:3) - deJ1/3.0
      deJ2 = sqrt(2.0/3.0*dot_product(deps_dev, deps_dev))

      depsp = deJ2*(sJ2/sigy/hard1(epsp, sigy, E, N))**m
      call dsig_dev_deps_e(Cklpq, E, nu)

      v2 = matmul(sig_dev, Cklpq)
      Cijpq_es = -9.0/4.0*depsp/sJ2**3
c      ! matmul(sig_dev, matmul(transpose(sig_dev), Cklpq))
      do i = 1, 6
        do j = 1, 6
          Cijpq_es(i, j) = Cijpq_es(i, j)*sig_dev(i)*v2(j)
        end do
      end do 

      Cijpq_et = 3.0*depsp/2.0/sJ2*Cklpq

      call deps_pl_deps_e(Cpl, sig, deps, epsp, sigy, E, nu, N, m)

      v1 = 3.0*sig_dev/2.0/sJ2 
c     ! Cijpq_p = matmul(3.0*sig_dev/2.0/sJ2, transpose(Cpl))
      do i = 1, 6
        do j = 1, 6
          Cijpq_p(i, j) = v1(i)*Cpl(j)
        end do
      end do

      Cijpq = Cijpq_es + Cijpq_et + Cijpq_p

      deallocate(sig_dev, deps_dev)
      deallocate(Cklpq, Cijpq_es, Cijpq_et, Cijpq_p)
      deallocate(v1, v2)
      end
      
      subroutine deps_p_depsp(Cijpl, sig, deps, epsp, sigy, E, N, m) 
      double precision, dimension(:), intent(in) :: sig, deps
      double precision :: epsp, sigy, E, N
      integer :: m
      double precision, dimension(:) :: Cijpl
      double precision, dimension(:), allocatable :: sig_dev, deps_dev
      double precision :: depsp, deps_pl_dp
      double precision :: sJ1, sJ2, deJ1, deJ2
      allocate(sig_dev(6), deps_dev(6))

      sig_dev = sig
      sJ1 = sig(1) + sig(2) + sig(3)
      sig_dev(1:3) = sig(1:3) - sJ1/3.0
      sJ2 = sqrt(3.0/2.0*dot_product(sig_dev, sig_dev))
      deps_dev = deps
      deJ1 = deps(1) + deps(2) + deps(3)
      deps_dev(1:3) = deps(1:3) - deJ1/3.0
      deJ2 = sqrt(2.0/3.0*dot_product(deps_dev, deps_dev))

      depsp = deJ2*(sJ2/sigy/hard1(epsp, sigy, E, N))**m
      deps_pl_dp = deps_pl_depsp(sig, deps, epsp, sigy, E, N, m)
      
      Cijpl = 3.0*sig_dev/2.0/sJ2*deps_pl_dp

      deallocate(sig_dev, deps_dev)
      end subroutine

      subroutine F(neq, t, y, ydot, rpar, ipar) 
      integer, intent(in) :: neq 
      double precision, dimension(:), intent(in) :: y
      double precision, dimension(:) :: ydot
      double precision, intent(in) :: t
      double precision, dimension(:), intent(in) :: rpar
      integer, dimension(:), intent(in) :: ipar

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
      eps_e = 0
      ! write(6, *) "allocate", eps_e
      ! write(6, *) "rpar", rpar
      ! write(6, *) "ipar", ipar
      ! write(6, *) y

      do i = 1, 6
      eps_e(i) = y(i)
      end do
      do j = 1, 6
      eps_p(j) = y(j)
      end do
      epsp = y(13)
      
      deps = rpar(1:6)
      E = rpar(7)
      nu = rpar(8)

      sigy = rpar(9)
      N = rpar(10)

      m = ipar(1)

      call hooke1(sig, eps_e, E, nu)
      
! deps_p
      call deps_p(deps_p_dt, sig, deps, epsp, sigy, E, N, m) 
      ydot(7:12) = deps_p_dt
! deps_e
      ydot(1:6) = deps - deps_p_dt
! depsp
      depsp_dt = deps_pl(sig, deps, epsp, sigy, E, N, m)
      ydot(13) = depsp_dt

      deallocate(eps_e, eps_p, sig, deps)
      deallocate(deps_e_dt, deps_p_dt)
      end subroutine

      subroutine JAC(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
      double precision, dimension(:,:) :: pd
      double precision, dimension(:) :: y
      integer :: neq, ml, mu, nrpd
      double precision :: t
      double precision, dimension(:) :: rpar
      integer, dimension(:) ::ipar

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
      return
      end
 

      end module dHwang1

