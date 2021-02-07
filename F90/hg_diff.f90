!> @file hg_diff.f90
!! @author Paul Barrère
!! @brief Le programme hg_diff.f90 simule la diffusion avec une fonction de
!! phase de Henyey-Greenstein.
!! @param d_mu Pas du cosinus de l'angle de diffusion
!! @param b Paramètre lié à la profondeur optique
!! @param phase Tableau stockant la fonction de phase
program hg_diff
  use param
  use random
  use sortie_diff
  use all_diff

  implicit none
  real(kind=kd) :: d_mu, b
  real(kind=kd), allocatable :: phase(:)

  ! Seed
  !call seed_file('seed.txt')
  call put_seed('seed.txt')

  ! Affichage de quelques paramètres
  print *, "Nombre de phases: ", N_phase
  print *, "Nombre de particules :", N_part
  print *, "Nombre de diffusions max :" , N_diff_max
  print *, "Chemin optique max :", t_max

  ! Ouverture fichier stockant les résultats
  open(12, file='Data/hg_diff.dat',status='unknown')

  ! Calcul de la fonction de phase
  call get_phase('hg',g,0._kd,t_max,b,d_mu,phase)

  !Sauvegarde des phases et rho
  call save_res(12, d_mu, phase)

  deallocate(phase)
  close(12)
end program hg_diff

