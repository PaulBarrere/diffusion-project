!> @file iso_diff.f90
!! @author Paul Barrère
!! @brief Le programme iso_diff.f90 simule la diffusion isotrope.
!! @param d_mu Pas du cosinus de l'angle de diffusion
!! @param b Paramètre lié à la profondeur optique
!! @param phase Tableau stockant la fonction de phase
program iso_diff
  use random
  use param
  use init_diff
  use sortie_diff
  use iso_mod

  implicit none
  real(kind=kd) :: d_mu, b
  real(kind=kd), allocatable :: phase(:)

  integer :: n_diff, N, k
  real(kind=kd) :: mu, t
  
  ! Seed
  !call seed_file('seed.txt')
  call put_seed('seed.txt')

  ! Mise en place fonction de phase
  call init(d_mu, b, t_max, phase)

  ! Affichage de quelques paramètres
  print *, "Nombre de phases: ", N_phase
  print *, "Nombre de particules :", N_part
  print *, "Nombre de diffusions max :" , N_diff_max
  print *, "Chemin optique max :", t_max

  !Ouverture fichier stockant les résultats
  open(11, file='iso_diff.dat',status='unknown')

  do N=1,N_part
    ! 1ère diffusion
    call first_diff(b, n_diff, mu, t) 

    ! Diffusions suivantes
    do while (t>0._kd .and. t<t_max)
      ! Vérification du nombre de diffusion
      if (n_diff>N_diff_max) then
        print *, 'Erreur: photon ne sort pas'
        exit
      end if
      
      ! Mise à jour du tau
      call update_mu_iso(mu)
      call update_path(n_diff, mu, t)
    end do 
    
    ! Sortie
    call update_phase(d_mu, mu, k, phase) 
  end do

  ! Sauvegarde de la fonction de phase
  call save_res(11, d_mu, phase)
  
  deallocate(phase)
  close(11)
end program iso_diff
