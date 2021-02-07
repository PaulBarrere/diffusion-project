!> @file iso_diff_abs.f90
!! @author Paul Barrère
!! @brief Le programme iso_diff_abs.f90 calcule la réflectance de la diffusion 
!! isotrope en fonction de l'albedo.
!! @param Na_max Nombre de valeurs pour l'albedo
!! @param N_abs Nombre de photons absorbés
!! @param N_diff Nombre de phtons diffusés
!! @param N_retro Nombre de photons rétrodiffusés
!! @param R Réflectance
!! @param albedo Tableau stockant les albedos
program iso_diff_abs
  use param
  use random

  implicit none
  integer, parameter :: Na_max=1e2
  integer :: Na, i 
  integer :: N, N_abs, N_retro, N_diff
  real(kind=kd) :: b
  real(kind=kd) :: mu, s, t
  real(kind=kd) :: R
  real(kind=kd) :: U0, U1, U2
  real(kind=kd), allocatable :: albedo(:)

  ! Seed
  !call seed_file('seed.txt')
  call put_seed('seed.txt')

  ! Albedos
  allocate(albedo(Na_max))
  albedo=(/ (i, i=0,Na_max-1, 1) /)
  albedo=real(albedo,kind=kd)/real(Na_max,kind=kd)

  ! Mise en place des paramètres
  b=1._kd-exp(-t_max)

  open(11, file='iso_diff_abs.dat',status='unknown')
    
  do Na=1,Na_max
    N_abs=0; N_diff=0; N_retro=0
    do N=1,N_part
      ! 1ère diffusion
      mu=1._kd 
      call random_number(U0)      
      t=-log(1._kd-b*U0)

      ! Diffusions suivantes
      do while (t>0._kd .and. t<t_max)
        call random_number(U0)  
        if (U0<=albedo(Na)) then
          call random_number(U1)
          call random_number(U2)
          mu=2._kd*U1-1._kd
          s=-log(U2)
          t=t+s*mu
        else
          N_abs=N_abs+1
          exit
        end if
      end do
      
      ! Sortie
      if (U0<=albedo(Na)) then
        if (t<=0) then
          N_retro=N_retro+1
        else
          N_diff=N_diff+1
        end if
      end if
    end do

    ! Sauvegarde
    R=b*real(N_retro,kind=kd)/real(N_part,kind=kd)

    write(11,*) albedo(Na), R
  end do
  close(11)
end program iso_diff_abs
