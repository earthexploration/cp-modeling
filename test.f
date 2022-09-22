      program main
          ! use bassani1
          use asaro2
          implicit none
          
      
          real, dimension(:), allocatable :: strains, stress
          real, dimension(:), allocatable :: rparams
          integer, dimension(:), allocatable :: iparams
          real :: t, dt
          integer :: i, j
          
          allocate(strains(7))
          allocate(rparams(6), iparams(12))
          allocate(stress(6))
         

          ! call init(1.0, 0.1, 0.5, 1.0, 1.0, 2.0, 0.01, 1.4)
          call init(1.0, 0.1, 0.5, 1.0, 1.0, 2.0, 1.4)

          strains = 0.0
          rparams = (/0.0, 0.0, 0.0, 1.0, 0.0, 0.0/)
          iparams = (/1, 1, -1, 
     1                -1, -1,-1,
     2                -1, -1,-1, 
     3                1, -1,-1/)
        
          stress = 0.0

          open(1, FILE="test_asaro2.dat")
          dt = 1e-2
          t = 0.0
          write(1, "(F10.4, F10.4)") strains(4), stress(4)
          do i = 1, 75
          write(6, "(a5, i3)") "step ", i
          call run_forward(7, strains, 0.0, dt, 
     1            rparams, iparams)
          t = t+dt
          stress = stress + rparams*dt
          write(6, "(9F10.4)") strains
          write(1, "(F10.4, F10.4)") strains(4), stress(4)
          end do
          

          close(1)

      end program
