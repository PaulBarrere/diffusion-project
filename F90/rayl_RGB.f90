!> @file rayl_RGB.f90
!! @author Paul Barrère
!! @brief Le programme rayl_RGB.f90 simule la diffusion de Rayleigh pour les couleurs RGB.
!! @param tau_max Profondeur optique
!! @param l Libre parcours moyen
!! @param b Paramètre lié à la profondeur optique
!! @param phase Tableau stockant la fonction de phase en fonction
!! de l'angle de diffusion
!! @param d_mu Pas du cosinus de l'angle de diffusion alpha
program rayl_RGB
  use param
  use random
  use sortie_diff
  use all_diff

  implicit none
  character(len=*), parameter :: files(3) = (/character(1024) :: "Data/rayl_diff_R.dat" &
  , "Data/rayl_diff_G.dat", "Data/rayl_diff_B.dat"/)
  real(kind=kd), allocatable :: l(:), tau_max(:), b(:)
  real(kind=kd), allocatable :: phase(:)
  real(kind=kd) :: d_mu
  integer :: i

  ! Mise en place système RBG
  allocate(l(3)); allocate(b(3)); allocate(tau_max(3))
  l=(/(lambda_R/lambda_B)**4, (lambda_G/lambda_B)**4, 1._kd/)
  tau_max=D/l

  ! Seed
  !call seed_file('seed.txt')
  call put_seed('seed.txt')

  ! Affichage de quelques paramètres
  print *, "Nombre de phases: ", N_phase
  print *, "Nombre de photons: ", N_part
  print *, "Paramètre D: ", D
  print *, " "
  print *, "Libres parcours moyens :"
  print *, "Rouge :", l(1)
  print *, "Vert :", l(2)
  print *, "Bleu :", l(3)
  print *, " "

  do i=1,3 
    ! Ouverture fichier stockant les résultats
    open(11+i, file=files(i),status='unknown')

    ! Calcul de la fonction de phase
    call get_phase('rayl',0._kd,0._kd,tau_max(i),b(i),d_mu,phase)

    ! Sauvegarde des résulats
    call save_res((11+i), d_mu, phase)
    deallocate(phase)
    close(11+i)
  end do
  open(16,file="Data/nb_photons_RGB_rayl.dat")
  write(16,*) real(N_part,kind=kd)*b(1)/b(3), real(N_part,kind=kd)*b(2)/b(3), real(N_part,kind=kd)
  close(16)

  print *, "Nombre de photons rouges diffusés : ", real(N_part,kind=kd)*b(1)/b(3)
  print *, "Nombre de photons verts diffusés : ", real(N_part,kind=kd)*b(2)/b(3)
  print *, "Nombre de photons bleus diffusés : ", real(N_part,kind=kd)

end program rayl_RGB
