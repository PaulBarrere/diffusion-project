!> @file mie_diff_CloudySky.f90
!! @author Paul Barrère
!! @brief Le programme mie_diff_CloudySky.f90 simule la diffusion de Mie à
!! travers un ciel nuageux.
!! @param d_alpha Pas de l'angle de diffusion pour l'enveloppe et bhmie
!! @param d_mu Pas du cosinus de l'angle de diffusion alpha
!! @param d_theta Pas de l'angle de diffusion pour la diffusion de Mie
!! @param d_phi Pas de l'angle azimuthal
!! @param S Vecteur de Stokes initial
!! @param Ss Vecteur de Stokes en sortie 
!! @param b Paramètre lié à la profondeur optique
!! @param u Axe de rotation
!! @param p Vecteur d'onde
!! @param phase_theta Tableau stockant la fonction de phase en fonction 
!! de l'angle de diffusion
!! @param phase_phi Tableau stockant la fonction de phase en fonction 
!! de l'azimuth
!! @param deg_pol Tableau stockant le degré de polarisation en fonction
!! de l'angle de diffusion
!! @param ang_pol Tableau stockant l'angle de polarisation en fonction
!! de l'angle de diffusion
!! @param phase_bhmie_elem Tableaux stockant la fonction de phase de bhmie
!! avec l'élément 'elem'
!! @param phase_env_elem Tableaux stockant l'enveloppe avec l'élément 'elem'
!! @param phase_norm_elem Rapport phase_bhmie_elem/phase_env_elem
!! @param ref_elem Rapport entre l'indice de réfraction de la particule 
!! diffusante 'elem' et du milieu environnant (ici l'air)
!! @param x_1,x_2 Facteur de forme du ciel clair et du nuage
!! @param abd_elem Albedo calculé pour l'élément 'elem'
!! @param S11_elem,S12_elem,S33_elem,S34_elem Composantes de la matrice de Mueller
!! pour l'élément 'elem'
!! @param eff_elem Efficacité de l'enveloppe la plus adaptée pour l'élément 'elem'
!! @param diff_type_elem Type de l'enveloppe: 'rayl' pour Rayleigh, 
!! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN pour l'élément 'elem'
!! @param para_elem Paramètre associé à diff_type_elem
!! @param N_abs Nombre de photons absorbés
!! @param N_diff Nombre de phtons diffusés
!! @param N_retro Nombre de photons rétrodiffusés

program mie_diff_CloudySky
  use param
  use random
  use init_diff
  use sortie_diff
  use simu_diff
  use init_mie_diff
  use sortie_mie_diff
  use simu_mie_diff
  use env_mod
  use fzero
  use quick_build_stereo_map
  use pixmap

  implicit none
  real(kind=kd) :: d_alpha, d_theta, d_phi, d_mu, b
  real(kind=kd), dimension(4) :: S, Ss
  real(kind=kd), dimension(3) :: u, p
  real(kind=kd), allocatable :: phase_theta(:,:), phase_phi(:,:), deg_pol(:), ang_pol(:)
  real(kind=kd), allocatable :: phase_bhmie_H2O(:), phase_env_H2O(:), phase_norm_H2O(:) 
  real(kind=kd), allocatable :: phase_bhmie_N2(:), phase_env_N2(:), phase_norm_N2(:)
  real(kind=kd), allocatable :: phase_bhmie_O2(:), phase_env_O2(:), phase_norm_O2(:)

  complex(kind=kd) :: ref_H2O, ref_N2, ref_O2
  real(kind=kd) :: x_1, x_2
  real(kind=kd) :: abd, abd_H2O, abd_N2, abd_O2
  real(kind=kd), allocatable :: S11_H2O(:), S12_H2O(:), S33_H2O(:), S34_H2O(:)
  real(kind=kd), allocatable :: S11_N2(:), S12_N2(:), S33_N2(:), S34_N2(:)
  real(kind=kd), allocatable :: S11_O2(:), S12_O2(:), S33_O2(:), S34_O2(:)

  real(kind=kd) :: eff_H2O, eff_N2, eff_O2
  real(kind=kd) :: para_H2O, para_N2, para_O2
  character(len=32) :: diff_type_H2O, diff_type_N2, diff_type_O2

  real(kind=kd) :: rd_abs, rd_O2N2
  real(kind=kd) :: t
  real(kind=kd) :: theta, phi
  real(kind=kd) :: R

  integer :: N_abs, N_retro, N_diff
  integer :: nb_diff, i 

  !integer, allocatable :: MapF(:,:), MapB(:,:)

  !> Seed
  !call seed_file('seed.txt')
  call put_seed('seed.txt')

  ! Affichage de quelques paramètres
  print *, "Nombres de phases: ", N_phase, N_theta, N_phi
  print *, "Nombre de particules :", N_part
  print *, "Nombre de diffusions max :" , N_diff_max
  print *, "Chemin optique max :", t_max
  print *, "Longueur d'onde:", lambda_B
  print *, " "

  ! Initialisation des variables diffusion mie
  call init_var_mie(d_alpha,d_theta,d_phi,b,t_max,N_abs,N_diff,N_retro)

  ! Initialisation des facteurs de forme
  x_1=2._kd*pi*0.01_kd/lambda_R
  x_2=2._kd*pi*10._kd/lambda_R
  print *, x_1, x_2
  ! Initialisation des arrays diffusion mie
  ! Pour l'eau
  call init_array_mie(phase_env_H2O,phase_norm_H2O,phase_theta,phase_phi,deg_pol,ang_pol)
  ! Pour le diazote
  call init_array_mie(phase_env_N2,phase_norm_N2,phase_theta,phase_phi,deg_pol,ang_pol)
  ! Pour le dioxygène
  call init_array_mie(phase_env_O2,phase_norm_O2,phase_theta,phase_phi,deg_pol,ang_pol)
  
  ! Initialisation des variables de bhmie et de la matrice de Mueller 
  ! Pour l'eau
  call init_bhmie('blue','H2O',x_2,abd_H2O,ref_H2O,S11_H2O,S12_H2O,S33_H2O,S34_H2O,phase_bhmie_H2O)
  ! Pour le diazote
  call init_bhmie('blue','N2',x_1,abd_N2,ref_N2,S11_N2,S12_N2,S33_N2,S34_N2,phase_bhmie_N2)
  ! Pour le dioxygène
  call init_bhmie('blue','O2',x_1,abd_O2,ref_O2,S11_O2,S12_O2,S33_O2,S34_O2,phase_bhmie_O2)
  
  print *, "Initialisation finie"

  ! Calcul de la meilleure enveloppe
  ! Pour l'eau
  call best_env(x_2,200,'fast',ref_H2O,eff_H2O,para_H2O,diff_type_H2O)
  call get_phase_th(diff_type_H2O,para_H2O,d_mu,phase_env_H2O)
  phase_norm_H2O=eff_H2O*phase_bhmie_H2O/phase_env_H2O
  
  ! Pour le diazote
  call best_env(x_1,200,'fast',ref_N2,eff_N2,para_N2,diff_type_N2)
  call get_phase_th(diff_type_N2,para_N2,d_mu,phase_env_N2)
  phase_norm_N2=eff_N2*phase_bhmie_N2/phase_env_N2

  ! Pour le dioxygène
  call best_env(x_1,200,'fast',ref_O2,eff_O2,para_O2,diff_type_O2)
  call get_phase_th(diff_type_O2,para_O2,d_mu,phase_env_O2)
  phase_norm_O2=eff_O2*phase_bhmie_O2/phase_env_O2


  ! Ouverture des fichiers sauvegardant les vecteurs de Stokes
  open(11,file='Data/theta_mie_diff_CloudySky.dat',status='unknown')
  open(12,file='Data/phi_mie_diff_CloudySky.dat',status='unknown')

  ! Au cas où polar=false, on a une lumière non polarisée
  S=S0

  do i=1,N_part
    ! Génération d'un vecteur de Stokes initial
    if (polar) then
      call generate_photon(S)
    end if
    ! 1ère Diffusion
    call first_diff_mie(b,S,Ss,u,p,nb_diff,t)
    
    ! Diffusions suivantes
    do while (exp(-t)<=1._kd .and. exp(-t)>=(1._kd-b))
      if (t<=0.3_kd) then
        call random_number(rd_O2N2)
        if (rd_O2N2<=0.78_kd) then
          call new_path(diff_type_N2,abd_N2,para_N2,phase_norm_N2,d_alpha,S11_N2,S12_N2,S33_N2,S34_N2,Ss,u,p,nb_diff,t,rd_abs)
          abd=abd_N2
        else
          call new_path(diff_type_O2,abd_O2,para_O2,phase_norm_O2,d_alpha,S11_O2,S12_O2,S33_O2,S34_N2,Ss,u,p,nb_diff,t,rd_abs)
          abd=abd_O2
        end if
      else
        call new_path(diff_type_H2O,abd_H2O,para_H2O,phase_norm_H2O,d_alpha,S11_H2O,S12_H2O,S33_H2O,S34_H2O,Ss,u,p,nb_diff,t,rd_abs)
        abd=abd_H2O
        !call new_path(diff_type_N2,abd_N2,para_N2,phase_norm_N2,d_alpha,S11_N2,S12_N2,S33_N2,S34_N2,Ss,u,p,nb_diff,t,rd_abs)
        !abd=abd_N2
      end if

      ! On regarde si le photon a été absorbé
      if ((rd_abs>abd) .or. (nb_diff>N_diff_max)) then
        exit
      end if
    end do
    
    ! Mise à jour du nombre de diffusion/rétrodiffusions
    if ((rd_abs>abd)  .or. (nb_diff>N_diff_max)) then
      N_abs=N_abs+1
    else if ((t<=0._kd) .and. (rd_abs<=abd)) then
      N_retro=N_retro+1
    else if ((t>0._kd) .and. (rd_abs<=abd)) then
      N_diff=N_diff+1
    end if

    ! Calcul de theta, de phi, et du vecteur de Stokes
    call final_stokes(u,p,theta,phi,Ss)
    
    ! Mise à jour du stockage des vecteurs de Stokes
    call update_phases_mie(d_theta,d_phi,theta,phi,Ss,phase_theta,phase_phi)
  end do
  print *, "Simulations terminées"

  ! Sauvegarde des résultats
  call save_res_mie(11,d_theta,N_theta,phase_theta)
  call save_res_mie(12,d_phi,N_phi,phase_phi)

  close(11)
  close(12) 
  print *, "Sauvegarde faite"

  ! Calcul de la réflectance, du degré et l'angle de polarisation
  R=b*N_retro/N_part
  deg_pol=dsqrt(phase_theta(:,2)**2+phase_theta(:,3)**2+phase_theta(:,4)**2)/phase_theta(:,1)
  ang_pol=0.5_kd*datan2(phase_theta(:,3),phase_theta(:,2))
  ! Projection stéréographique
  !call stereo_map(N_theta,deg_pol,MapF,MapB) 
  !call ppm(13,'deg_pol_mie_diff.ppm',255,MapF)
  !call ppm(14,'deg_pol_mie_rdiff.ppm',255,MapB)
  print *, "Projection faite"
  print *, " "

  print *, 'Nombre de photons diffusés: ',N_diff
  print *, 'Nombre de photons rétrodiffusés: ',N_retro 
  print *, 'Nombre de photons absorbés: ',N_part-N_diff-N_retro
  print *, 'Réflectance du milieu: ',R

end program mie_diff_CloudySky
