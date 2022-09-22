      module bassani1
          implicit none
      
          integer, dimension(:,:), allocatable :: fcc_slipnor 
          integer, dimension(:,:), allocatable :: fcc_primary
          integer, dimension(:,:), allocatable :: fcc_conjugate
          integer, dimension(:,:), allocatable :: fcc_crossglide
          integer, dimension(:,:), allocatable :: fcc_critical
          
          real, dimension(:,:), allocatable :: fcc_C

          real :: fcc_c11, fcc_c12, fcc_c44
          real :: bassani_h0, bassani_tau0, bassani_taus, bassani_hs
          real :: bassani_q

          integer :: fcc_a1, fcc_a2, fcc_a3, fcc_a4, fcc_a5
         
          public :: fcc_slipnor
          public :: fcc_primary, fcc_conjugate
          public :: fcc_crossglide, fcc_critical

          private :: fcc_c11, fcc_c12, fcc_c44
          private :: bassani_h0, bassani_tau0, bassani_taus, bassani_hs
          private :: bassani_q
          private :: fcc_a1, fcc_a2, fcc_a3, fcc_a4, fcc_a5
      contains

          subroutine init(c11, c12, c44, h0, tau0, taus, hs, q)
          real, intent(in) :: c11, c12, c44, h0, tau0, taus, hs
          real, intent(in) :: q
              allocate(fcc_slipnor(3, 4))
              allocate(fcc_primary(3, 3))
              allocate(fcc_conjugate(3, 3))
              allocate(fcc_crossglide(3, 3))
              allocate(fcc_critical(3, 3))

              !allocate(fcc_C(6,6))

              fcc_slipnor = reshape((/
     1                  1, 1, -1, 
     2                  1, -1, 1,
     3                  1, -1, -1,
     4                  1,  1,  1 
     5                 /), shape(fcc_slipnor))
                  
              fcc_primary = reshape((/
     1                   1, 0, 1,
     2                  0, -1, -1,
     3                  -1, 1, 0
     4                  /), shape(fcc_primary))

              fcc_conjugate = reshape((/
     1                  -1, -1, 0,
     2                  0, 1, 1,
     3                  1, 0, -1
     4                  /), shape(fcc_conjugate))

              fcc_crossglide = reshape((/
     1                  -1, 0, -1,
     2                  1, 1, 0,
     3                  0, -1, 1
     4                  /), shape(fcc_crossglide))

              fcc_critical = reshape((/
     1                  -1, 0, 1,
     2                  0, 1, -1,
     3                  1, -1, 0
     4                  /), shape(fcc_critical))
 

              fcc_c11 = c11
              fcc_c12 = c12
              fcc_c44 = c44

              bassani_h0 = h0
              bassani_tau0 = tau0
              bassani_taus = taus
              bassani_hs = hs

              bassani_q = q
              !fcc_C = 0

              fcc_a1 = 8
              fcc_a2 = 8
              fcc_a3 = 8
              fcc_a4 = 15
              fcc_a5 = 20
          end subroutine 

          subroutine elastic_moduli(C, c11, c12, c44)
          real, intent(in) :: c11, c12, c44
          real, dimension(:,:), intent(out) :: C

          integer :: i,j

          C(1,1) = c11
          C(2,2) = c11
          C(3,3) = c11

          C(4,4) = c44
          C(5,5) = c44
          C(6,6) = c44

          C(1,2) = c12
          C(1,3) = c12
          C(2,3) = c12

          ! isotropic
          ! c44 = (c11-c12)/2.0

          do i = 1, 6
            do j = 1, i
              C(i, j) = C(j, i)
            end do
          end do

          end subroutine

          ! axb = c
          subroutine cross_product(c, b, a)
          integer, dimension(:), intent(in) :: b, a
          integer, dimension(:), intent(out) :: c

          c(1) = a(2)*b(3) - a(3)*b(2)
          c(2) = a(3)*b(1) - a(1)*b(3)
          c(3) = a(1)*b(2) - a(2)*b(1)

          end subroutine

          subroutine slip_system(slipsys, slipdir, slipnor)
          integer, dimension(:), intent(in) :: slipdir, slipnor
          real, dimension(:,:), intent(out) :: slipsys

          integer, dimension(3) :: sliplin
          real :: er, es, et
          er = 0
          es = 0
          et = 0

          call cross_product(sliplin, slipdir, slipnor)

          er = sqrt(real(dot_product(slipnor, slipnor)))
          es = sqrt(real(dot_product(slipdir, slipdir)))
          et = sqrt(real(dot_product(sliplin, sliplin)))

          slipsys(:,1) = slipnor/er
          slipsys(:,2) = slipdir/es
          slipsys(:,3) = sliplin/et 

          end subroutine

          ! (x, y, z) -> (r, s, t)
          subroutine rotation_r(R, ex, ey, ez, er, es, et)
          real, dimension(:), intent(in) :: ex, ey, ez
          real, dimension(:), intent(in) :: er, es, et
          real, dimension(:,:), intent(out) :: R
          
          R(1,1) = dot_product(ex, er)
          R(2,1) = dot_product(ex, es)
          R(3,1) = dot_product(ex, et)

          R(1,2) = dot_product(ey, er)
          R(2,2) = dot_product(ey, es)
          R(3,2) = dot_product(ey, et)

          R(1,3) = dot_product(ez, er)
          R(2,3) = dot_product(ez, es)
          R(3,3) = dot_product(ez, et)

          end subroutine

          subroutine rotation_q(Q, ex, ey, ez, er, es, et)
          real, dimension(:), intent(in) :: ex, ey, ez
          real, dimension(:), intent(in) :: er, es, et
          real, dimension(:,:), intent(out) :: Q

          real :: l1, l2, l3
          real :: m1, m2, m3
          real :: n1, n2, n3
          real, dimension(:,:), allocatable :: R
          l1 = 0
          l2 = 0
          l3 = 0
          m1 = 0
          m2 = 0
          m3 = 0
          n1 = 0
          n2 = 0
          n3 = 0

          allocate(R(3,3))
          R = 0
          
          l1 = dot_product(ex, er)
          l2 = dot_product(ex, es)
          l3 = dot_product(ex, et)

          m1 = dot_product(ey, er)
          m2 = dot_product(ey, es)
          m3 = dot_product(ey, et)

          n1 = dot_product(ez, er)
          n2 = dot_product(ez, es)
          n3 = dot_product(ez, et)

          R(1,1) = l1
          R(2,1) = l2
          R(3,1) = l3
          R(1,2) = m1
          R(2,2) = m2
          R(3,2) = m3
          R(1,3) = n1
          R(2,3) = n2
          R(3,3) = n3

          Q(1:3, 1:3) = R*R
          ! voigt
          Q(4, 1) = l2*l3
          Q(5, 1) = l1*l3
          Q(6, 1) = l1*l2

          Q(4, 2) = m2*m3
          Q(5, 2) = m1*m3
          Q(6, 2) = m1*m2

          Q(4, 3) = n2*n3
          Q(5, 3) = n1*n3
          Q(6, 3) = n1*n2

          Q(1, 4) = 2.0*m1*n1
          Q(2, 4) = 2.0*m2*n2
          Q(3, 4) = 2.0*m3*n3

          Q(1, 5) = 2.0*l1*n1
          Q(2, 5) = 2.0*l2*n2
          Q(3, 5) = 2.0*l3*n3

          Q(1, 6) = 2.0*l1*m1
          Q(2, 6) = 2.0*l2*m2
          Q(3, 6) = 2.0*l3*m3

          Q(4, 4) = m2*n3 + n2*m3
          Q(5, 4) = m1*n3 + n1*m3
          Q(6, 4) = m1*n2 + n1*m2

          Q(4, 5) = l2*n3 + n2*l3
          Q(5, 5) = l1*n3 + n1*l3
          Q(6, 5) = l1*n2 + n1*l2
          
          Q(4, 6) = l2*m3 + m2*l3
          Q(5, 6) = l1*m3 + m1*l3
          Q(6, 6) = l1*m2 + m1*l2

          deallocate(R)
          end subroutine

          subroutine slip_system_r(R, slipdir, slipnor)
          integer, dimension(:), intent(in) :: slipdir, slipnor
          real, dimension(:,:), intent(out) :: R
          
          real, dimension(:), allocatable :: ex, ey, ez
          real, dimension(:,:), allocatable :: slipsys
          allocate(ex(3), ey(3), ez(3))
          allocate(slipsys(3,3))

          ex = (/1, 0, 0/)
          ey = (/0, 1, 0/)
          ez = (/0, 0, 1/)
          slipsys = 0

          call slip_system(slipsys, slipdir, slipnor)
            
          call rotation_r(R, ex, ey, ez, 
     1 slipsys(:,1), slipsys(:,2), slipsys(:,3))


          deallocate(slipsys)
          deallocate(ex, ey, ez)
          end subroutine

          subroutine slip_system_q(Q, slipdir, slipnor)
          integer, dimension(:), intent(in) :: slipdir, slipnor
          real, dimension(:,:), intent(out) :: Q
          
          real, dimension(:), allocatable :: ex, ey, ez
          real, dimension(:,:), allocatable :: slipsys
          allocate(ex(3), ey(3), ez(3))
          allocate(slipsys(3,3))

          ex = (/1, 0, 0/)
          ey = (/0, 1, 0/)
          ez = (/0, 0, 1/)
          slipsys = 0

          call slip_system(slipsys, slipdir, slipnor)
            
          call rotation_q(Q, ex, ey, ez, 
     1 slipsys(:,1), slipsys(:,2), slipsys(:,3))


          deallocate(slipsys)
          deallocate(ex, ey, ez)
          end subroutine

          subroutine slipsys_Schmid(schmid, slipdir, slipnor)
          integer, dimension(:), intent(in) :: slipdir, slipnor
          real, dimension(:,:), intent(out) :: schmid

          real, dimension(:,:), allocatable :: slipsys
          real, dimension(:), allocatable :: en, em
          integer :: i, j
          allocate(slipsys(3,3))
          allocate(en(3), em(3))
          slipsys = 0
          en = 0
          em = 0
          call slip_system(slipsys, slipdir, slipnor)
          en = slipsys(:,1)
          em = slipsys(:,2)
          do i = 1, 3
            do j = 1, 3
              schmid(i,j) = em(i)*en(j)
            end do
          end do

          deallocate(slipsys)
          deallocate(en, em)
          end subroutine

          subroutine slipsys_symSchmid(symSchmid, slipdir, slipnor)
          integer, dimension(:), intent(in) :: slipdir, slipnor
          real, dimension(:,:), intent(out) :: symSchmid

          real, dimension(:,:), allocatable :: schmid
          allocate(schmid(3,3))
          schmid = 0
          call slipsys_Schmid(schmid, slipdir, slipnor)
          symSchmid = 0.5*(schmid + transpose(schmid))
          deallocate(schmid)
          end subroutine

          function schmid_shear(stress, slipdir, slipnor)
          integer, dimension(:), intent(in) :: slipdir, slipnor
          real, dimension(:,:), intent(in) :: stress
          real :: schmid_shear

          real :: shr1, shr2, shr3
          real, dimension(:,:), allocatable :: symSchmid
          allocate(symSchmid(3,3))
          symSchmid = 0
          call slipsys_symSchmid(symSchmid, slipdir, slipnor)
          shr1 = 0
          shr2 = 0
          shr3 = 0
          shr1= dot_product(stress(:,1), symSchmid(:,1))
          shr2= dot_product(stress(:,2), symSchmid(:,2))
          shr3= dot_product(stress(:,3), symSchmid(:,3))

          schmid_shear = shr1 + shr2 + shr3
          deallocate(symSchmid)
          end function

          function hard(p, pgamma, h0, tau0, taus, hs, iflags)
          real, intent(in) :: p
          real, intent(in) :: h0, tau0, taus, hs
          real, dimension(:), intent(in) :: pgamma
          integer :: n
          integer, dimension(:), intent(in) :: iflags
          real :: hard

          integer, dimension(:,:), allocatable :: hcoeff, icoeff
          real :: p0
          real :: F, G

          integer, dimension(:), allocatable :: nslipsys
          integer :: i, j

          n = 0
          do i = 1, 12
            if (iflags(i) >= 0) then
                n = n + 1
            end if
          end do
          !print *, n
          allocate(nslipsys(n))

          n = 0
          do i = 1, 12 
          if (iflags(i) >= 0) then
                n = n + 1
                ! index from 1
                nslipsys(n) = i
            end if
          end do
 

          allocate(hcoeff(12,12))
          allocate(icoeff(n, n))
          hcoeff = 0
          F = 0
          G = 0
          p0 = 0
          icoeff = 0

          if (taus < tau0) then
              print *, ' tau0 < taus '
              call abort
          end if
          if (h0 < hs) then
              print *, ' hs < h0 '
              call abort
          end if
          if (n > 12) then
              print *, " maximum slip systems: 12"
              call abort
          end if

          !hard = h0*(1.0/cosh(h0*p/(taus-tau0)))**2
          p0 = (taus-tau0)/h0
          F = (h0-hs)*(1.0/cosh((h0-hs)*p/(taus-tau0)))**2
     &      + hs
          
          call hardcoeff_matrix(hcoeff, 
     1      fcc_a1, fcc_a2, fcc_a3, fcc_a4, fcc_a5)

          do i = 1, n
          do j = 1, n 
            icoeff(i,j) = hcoeff(nslipsys(i), nslipsys(j))
          end do
          end do

          G = 1+sum(matmul(icoeff, tanh(pgamma/p0)))

          hard = F*G
          return
          end function


          subroutine hard_moduli(H, N, p, pgamma, 
     1      h0, tau0, taus, hs, q, iflags)
          real, intent(in) :: p
          real, dimension(:), intent(in) :: pgamma
          real, intent(in) :: h0, tau0, taus, hs, q
          integer, intent(in) :: N
          integer, dimension(:), intent(in) :: iflags

          real, dimension(:,:), intent(out) :: H

          integer :: i, j

          do i = 1, N
            do j = 1, N
              if (i .eq. j) then
                  H(i, j) = hard(p, pgamma, h0, tau0, taus, hs, 
     1                          iflags)
              else 
                  H(i, j) = q*hard(p, pgamma, h0, tau0, taus, hs, 
     1                          iflags)
              end if 
            end do
          end do   

          end subroutine

          function hard_coeff(nslipdir1, nslipnor1, 
     1      nslipdir2, nslipnor2, a1, a2, a3, a4, a5)
          integer :: nslipdir1, nslipnor1
          integer :: nslipdir2, nslipnor2
          ! character :: a1, a2, a3, a4, a5
          integer :: a1, a2, a3, a4, a5
          integer :: a0 
          ! character :: hard_coeff
          integer :: hard_coeff

          integer, dimension(:), allocatable :: slipdir1, slipdir2
          allocate(slipdir1(3), slipdir2(3))
          slipdir1 = 0
          slipdir2 = 0

          call slipsys_dir(slipdir1, nslipdir1, nslipnor1)
          call slipsys_dir(slipdir2, nslipdir2, nslipnor2)

          a0 = 0
          hard_coeff = 0

          a0 = dot_product(slipdir1, slipdir2)

          if (nslipnor1 .ne. nslipnor2) then
              if (a0 .eq. -2) then
                  ! 'N : no junction'
                  hard_coeff = a1
                  return 
              end if
              if (a0 .eq. 0) then
                  ! 'H: Hirth lock'
                  hard_coeff = a2
                  return
              end if
              if (a0 .eq. 1) then
                  ! 'G: glissile junction'
                  hard_coeff = a4
                  return
              end if
              if (a0 .eq. -1) then
                  ! 'S: sessile junction'
                  hard_coeff = a5
                  return
              end if
          end if


          if (nslipnor1 .eq. nslipnor2) then
              if (a0 .eq. -1) then
                ! 'C: coplanar junction'
                hard_coeff = a3
                return
              end if
          end if

          deallocate(slipdir1, slipdir2)
          return
          end function

          subroutine hardcoeff_matrix(hcoeff, a1, a2, a3, a4, a5)
          !character, dimension(:,:), intent(out) :: hcoeff
          integer, dimension(:,:), intent(out) :: hcoeff
          !character, intent(in) :: a1, a2, a3, a4, a5
          integer, intent(in) :: a1, a2, a3, a4, a5

          integer :: i, j, k, l
          integer :: m, n
          !character, dimension(:,:,:,:), allocatable :: hcoeff1
          integer, dimension(:,:,:,:), allocatable :: hcoeff1
          allocate(hcoeff1(4, 3, 4, 3))
          hcoeff1 = 1

          do i = 1, 4
            do j = 1, 3
              do k = 1, 4
                do l = 1, 3
                  hcoeff1(i, j, k, l) = hard_coeff(j-1, i-1, l-1, k-1,
     &                                a1, a2, a3, a4, a5)
                  ! write(6, "(a3)") hcoeff1(i, j, k, l)
                  m = 3*(i-1) + (j-1)
                  n = 3*(k-1) + (l-1)
                  hcoeff(m+1,n+1) = hcoeff1(i, j, k, l)
                end do 
              end do
              ! print *
            end do
          end do

          deallocate(hcoeff1)
          end subroutine

          subroutine slipsys_index(nslipdir, nslipnor, n)
          integer, intent(in) :: n
          integer, intent(out) :: nslipdir, nslipnor
          nslipdir = 0
          nslipnor = 0

          if ( (n .lt. 0) .or. (n .ge. 12)) then
              print *, " slipsys index fault "
              call abort
          end if

          nslipnor = n/3
          nslipdir = mod(n, 3)

          end subroutine 

          subroutine slipsys_dir(slipdir, nslipdir, nslipnor)
          integer, intent(in) :: nslipdir, nslipnor
          integer, dimension(:), intent(out) :: slipdir
          integer :: ndir
          slipdir = 0
          ndir = 0 
          if ( (nslipnor .lt. 0) .or. (nslipnor .ge. 4)) then
              print *, " slipnor index fault "
              call abort
          end if
          ndir = nslipdir + 1
          if ( (nslipdir .lt. 0) .or. (nslipdir .ge. 3)) then
              print *, " slipdir index fault "
              call abort
          end if

          if (nslipnor .eq. 0) then
              slipdir = fcc_primary(:, ndir)
          else if (nslipnor .eq. 1) then
              slipdir = fcc_conjugate(:, ndir)
          else if (nslipnor .eq. 2) then
              slipdir = fcc_crossglide(:, ndir)
          else if (nslipnor .eq. 3) then
              slipdir = fcc_critical(:, ndir)
          end if

          end subroutine

          subroutine slipsys_nor(slipnor, nslipnor)
          integer, intent(in) :: nslipnor
          integer, dimension(:), intent(out) :: slipnor
          integer :: nnor
          slipnor = 0
          nnor = 0
          if ( (nslipnor .lt. 0) .or. (nslipnor .ge. 4)) then
              print *, " slipnor index fault "
              call abort
          end if
          nnor = nslipnor + 1 
          slipnor = fcc_slipnor(:, nnor)
          end subroutine

          subroutine voigt(v, symm)
          real, dimension(:,:), intent(in) :: symm
          real, dimension(:), intent(out) :: v

          integer :: i, j

          do i = 1, 3
            do j = 1, i
              if (j .eq. i) then
                  v(i) = symm(i, j)
              else
                  v(9-i-j) = symm(i, j)
              end if
            end do
          end do

          end subroutine

          subroutine voigt2(symm, v)
          real, dimension(:), intent(in) :: v
          real, dimension(:,:), intent(out) :: symm

          integer :: i, j
          
          do i = 1, 3
            do j = 1, 3
              if (i .eq. j) then
                  symm(i, j) = v(i)
              else
                  symm(i, j) = v(9-i-j)
              end if
            end do
          end do

          end subroutine

          subroutine invsym(a, n)
          real, dimension(:,:), allocatable :: a
          integer :: n

          character :: uplo
          integer :: nrhs
          integer :: lda
          integer, dimension(:), allocatable :: ipiv
          real, dimension(:,:), allocatable :: b
          integer :: ldb
          real, dimension(:), allocatable :: work
          integer :: lwork
          integer :: info

          real, dimension(:,:), allocatable :: x,y
          integer :: i, j

          uplo = 'U'
          info = -1
          lda = n

          ldb = n
          nrhs = n

          lwork = lda*n*ldb*nrhs

          allocate(b(ldb, nrhs))
          allocate(ipiv(n))
          allocate(work(lwork))
          ipiv = 0
          work = 0
          
          b = 0
          do i = 1, n
            do j = 1, n
              if (i==j) then
                b(i, j) = 1
              end if
            end do
          end do

          allocate(x(lda, n))
          allocate(y(ldb, nrhs))
          x = a
          y = b 

          call ssysv(uplo, n, nrhs, a, lda, ipiv,
     1      b, ldb, work, lwork, info)
              
          if (info == 0) then
              a = b
              return
          end if

          if (info > 0) then
              print *, " singular matrix "
              a = x
              return
          end if

          if (info < 0) then
              print *,  -info, " illegal value "
              a = y
              return
          end if

          deallocate(b, ipiv, work)
          deallocate(x, y)
          end subroutine 
          

          subroutine run_forward(neq, strains, t0, dt,
     1      rparams, iparams)
          integer, intent(in) :: neq
          real, dimension(:) :: strains
          real, intent(in) :: t0
          real, intent(in) :: dt
          real, dimension(:), intent(in) :: rparams
          integer, dimension(:), intent(in) :: iparams
          
          integer :: i, j
          real :: t, h
          real, dimension(:), allocatable :: y0, ydot
          real, dimension(:), allocatable :: y1, y2
          real, dimension(:), allocatable :: y3, y4

          allocate(y0(neq), ydot(neq))
          allocate(y1(neq), y2(neq))
          allocate(y3(neq), y4(neq))
          y0 = 0
          ydot = 0
          y1 = 0
          y2 = 0
          y3 = 0
          y4 = 0

          y0 = strains(1:neq)
          y1 = y0
          y2 = y0
          y3 = y0
          y4 = y0

          ! rk4
            t = t0          
            call shear_rate(neq, t, y0, ydot, 
     1          rparams, iparams)
            y1 = ydot
            
            t = t0 + dt/2.0
            y2 = y0 + dt/2.0*y1
            call shear_rate(neq, t, y2, ydot,
     2          rparams, iparams)
            y2 = ydot
            
            t = t0 + dt/2.0
            y3 = y0 + dt/2.0*y2
            call shear_rate(neq, t, y3, ydot,
     3          rparams, iparams)
            y3 = ydot
            
            t = t0 + dt
            y4 = y0 + dt*y3
            call shear_rate(neq, t, y4, ydot,
     4          rparams, iparams)
            y4 = ydot
     
         
          
          strains(1:neq) = y0 + dt/6.0*(y1 + 2.0*y2 + 2.0*y3 + y4)

          deallocate(y0, ydot)
          deallocate(y1, y2)
          deallocate(y3, y4)
          end subroutine


          subroutine shear_rate(neq, t, strains, dstrains, 
     1        rparams, iparams)
          integer, intent(in) :: neq
          real, intent(in) :: t
          real, dimension(:), intent(in) :: strains
          real, dimension(:), intent(in) :: rparams
          integer, dimension(:), intent(in) :: iparams
          real, dimension(:), intent(out) :: dstrains

          real, dimension(:), allocatable :: dstress
          integer, dimension(:), allocatable :: iflags
          real, dimension(:), allocatable :: nslipsys

          real, dimension(:), allocatable :: pgamma
          real, dimension(:,:), allocatable :: H
          real, dimension(:), allocatable :: dtau
          real, dimension(:), allocatable :: dgamma

          integer :: i, j, k, l
          integer :: n

          real, dimension(:), allocatable :: pstrain
          real :: p
          integer, dimension(:), allocatable :: slipdir, slipnor

          real, dimension(:,:), allocatable :: v, Lp
          real, dimension(:), allocatable :: dpstrain

          allocate(dstress(6))
          allocate(iflags(12))

          dstress = rparams(1:6)

          iflags = iparams(1:12)

          n = 0
          do i = 1, 12
            if (iflags(i) >= 0) then
                n = n + 1
            end if
          end do
          !print *, n
          allocate(nslipsys(n))
          allocate(pgamma(n))
          allocate(H(n,n))
          allocate(dtau(n))
          allocate(dgamma(n))

          n = 0
          do i = 1, 12 
          if (iflags(i) >= 0) then
                n = n + 1
                ! index from 0 
                nslipsys(n) = i-1
            end if
          end do
          

          allocate(pstrain(6))
          pstrain = 0
          p = 0

          pstrain = strains(1:6)
          p = strains(7)

          allocate(slipdir(3), slipnor(3))
          slipdir = 0
          slipnor = 0

          allocate(v(3,3), Lp(3,3))
          allocate(dpstrain(6))
          v = 0
          Lp = 0
          dpstrain = 0

        
          pgamma = 0
          dtau = 0
          do j = 1, n
            i = nslipsys(j)
            call slipsys_index(k, l, i)
            !slipnor = fcc_slipnor(:, l+1)
            call slipsys_nor(slipnor, l)
            call slipsys_dir(slipdir, k, l)

            call voigt2(v, pstrain)
            ! shear strain
            pgamma(j) = schmid_shear(v, slipdir, slipnor)
            ! shear stress
            call voigt2(v, dstress)
            dtau(j) = schmid_shear(v, slipdir, slipnor)

          end do 
            
          H = 0
          call hard_moduli(H, n, p, pgamma,
     1    bassani_h0, bassani_tau0, bassani_taus, bassani_hs, bassani_q,
     2    iflags)
          call invsym(H, n)
         
          dgamma = 0
          dgamma = matmul(H, dtau)

          Lp = 0
          do j = 1, n
            i = nslipsys(j)
            call slipsys_index(k, l, i)
            !slipnor = fcc_slipnor(:, l+1)
            call slipsys_nor(slipnor, l)
            call slipsys_dir(slipdir, k, l)
            
            call slipsys_symSchmid(v, slipdir, slipnor) 
            Lp = Lp + dgamma(j)*v

          end do
          
          call voigt(dpstrain, Lp)

          dstrains(1:6) = dpstrain
          dstrains(7) = sum(abs(dgamma))

          deallocate(dstress, iflags, nslipsys)
          deallocate(pgamma, H, dtau, dgamma) 
          deallocate(pstrain, slipdir, slipnor)
          deallocate(v, Lp)
          deallocate(dpstrain)
          end subroutine



C----------------------------------------------------------------------


      end module bassani1
