!> @file mie_diff.f90
!! @author Paul Barrère
!! @brief Le programme mie_diff.f90 simule la diffusion de Mie.
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
!! @param phase_bhmie Tableaux stockant la fonction de phase de bhmie
!! @param phase_env Tableaux stockant l'enveloppe
!! @param phase_norm Rapport phase_bhmie/phase_env
!! @param refrel Rapport entre l'indice de réfraction de la particule 
!! diffusante et du milieu environnant (ici l'air)
!! @param albedo Albedo calculé 
!! @param S11,S12,S33,S34 Composantes de la matrice de Mueller
!! @param eff Efficacité de l'enveloppe la plus adaptée
!! @param diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
!! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN
!! @param para Paramètre associé à diff_type
!! @param N_abs Nombre de photons absorbés
!! @param N_diff Nombre de phtons diffusés
!! @param N_retro Nombre de photons rétrodiffusés
program mie_diff
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
  real(kind=kd), allocatable :: phase_bhmie(:), phase_env(:), phase_norm(:), deg_pol(:), ang_pol(:) 
  real(kind=kd), allocatable :: phase_theta(:,:), phase_phi(:,:)

  complex(kind=kd) :: refrel
  real(kind=kd), dimension(3) :: u, p
  real(kind=kd) :: albedo
  real(kind=kd), allocatable :: S11(:), S12(:), S33(:), S34(:)

  real(kind=kd) :: eff, para
  character(len=32) :: diff_type

  real(kind=kd) :: t, rd
  real(kind=kd) :: theta, phi
  real(kind=kd) :: R

  integer :: N_abs, N_retro, N_diff
  integer :: nb_diff, i

  integer, allocatable :: MapF(:,:), MapB(:,:)

  ! Seed
  !call seed_file('seed.txt')
  call put_seed('seed.txt')

  ! Affichage de quelques paramètres
  print *, "Nombres de phases: ", N_phase, N_theta, N_phi
  print *, "Nombre de particules :", N_part
  print *, "Nombre de diffusions max :" , N_diff_max
  print *, "Chemin optique max :", t_max
  print *, "Facteur de forme :", x_form
  print *, " "
  
  ! Initialisation des variables et arrays diffusion mie
  call init_array_mie(phase_env,phase_norm,phase_theta,phase_phi,deg_pol,ang_pol)
  call init_var_mie(d_alpha,d_theta,d_phi,b,t_max,N_abs,N_diff,N_retro)
  ! Initialisation des variables de bhmie et de la matrice de Mueller
  call init_bhmie('red','H2O',x_form,albedo,refrel,S11,S12,S33,S34,phase_bhmie)
  print *, "Initialisation finie"

  ! Calcul de la meilleure enveloppe
  call best_env(x_form,200,'fast',refrel,eff,para,diff_type)
  call get_phase_th(diff_type,para,d_mu,phase_env)
  phase_norm=eff*phase_bhmie/phase_env
  print *, "La meilleure enveloppe est: ", trim(diff_type)
  
  ! Ouverture des fichiers sauvegardant les vecteurs de Stokes
  open(11,file='Data/theta_mie_diff.dat',status='unknown')
  open(12,file='Data/phi_mie_diff.dat',status='unknown')
  print *, "Albedo: ",albedo
  
  ! Au cas où polar=false, on a une lumière non polarisée
  S=S0

  do i=1,N_part
    ! Génération d'un vecteur de Stokes initial
    if (polar) then
      call generate_photon(S)
    end if  

    ! Mise en place de la 1ère Diffusion
    call first_diff_mie(b,S,Ss,u,p,nb_diff,t)
        
    ! Diffusions suivantes
    do while (exp(-t)<=1._kd .and. exp(-t)>=(1._kd-b))
      
      ! Mise à jour de la trajectoire du photon
      call new_path(diff_type,albedo,para,phase_norm,d_alpha,S11,S12,S33,S34,Ss,u,p,nb_diff,t,rd)

      ! On regarde si le photon a été absorbé
      if ((rd>albedo) .or. (nb_diff>N_diff_max)) then
        exit
      end if
    end do
    
    ! Mise à jour du nombre de diffusion/rétrodiffusions
    if ((rd>albedo) .or. (nb_diff>N_diff_max)) then
      N_abs=N_abs+1
    else if (t<=0._kd .and. rd<=albedo) then
      N_retro=N_retro+1
    else if (t>0._kd .and. rd<=albedo) then
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

  ! Calcul de la réflectance, du degré et de l'angle de polarisation
  R=b*N_retro/N_part
  deg_pol=dsqrt(phase_theta(:,2)**2+phase_theta(:,3)**2+phase_theta(:,4)**2)/phase_theta(:,1)
  ang_pol=0.5_kd*datan2(phase_theta(:,3),phase_theta(:,2)) 
  open(132,file='Data/degpol_nonpol_1_1.dat')
  open(133,file='Data/angpol_nonpol_1_1.dat')
  call save_res(132,d_theta,deg_pol)
  call save_res(133,d_theta,ang_pol)
  close(132)
  close(133)

  ! Projection stéréographique du degré de polarisation
  !call stereo_map(N_theta,deg_pol,MapF,MapB) 
  !call ppm(13,'deg_pol_mie_diff.ppm',255,MapF)
  !call ppm(14,'deg_pol_mie_rdiff.ppm',255,MapB)

  ! Projection stéréographique de I,Q,U et V
  call stereo_map(N_theta,phase_theta(:,1),MapF,MapB)
  call ppm(15,'stereo_mie_I_diff_nonpol_1_1.ppm',255,MapF)
  call ppm(16,'stereo_mie_I_rdiff_nonpol_1_1.ppm',255,MapB)
  call stereo_map(N_theta,phase_theta(:,2),MapF,MapB)
  call ppm(17,'stereo_mie_Q_diff_nonpol_1_1.ppm',255,MapF)
  call ppm(18,'stereo_mie_Q_rdiff_nonpol_1_1.ppm',255,MapB)
  call stereo_map(N_theta,phase_theta(:,3),MapF,MapB)
  call ppm(19,'stereo_mie_U_diff_nonpol_1_1.ppm',255,MapF)
  call ppm(20,'stereo_mie_U_rdiff_nonpol_1_1.ppm',255,MapB)
  call stereo_map(N_theta,phase_theta(:,4),MapF,MapB)
  call ppm(21,'stereo_mie_V_diff_nonpol_1_1.ppm',255,MapF)
  call ppm(22,'stereo_mie_V_rdiff_nonpol_1_1.ppm',255,MapB)

  print *, "Projection faite"
  print *, " "

  print *, 'Nombre de photons diffusés: ',N_diff
  print *, 'Nombre de photons rétrodiffusés: ',N_retro 
  print *, 'Nombre de photons absorbés: ', N_abs
  print *, 'Réflectance du milieu: ',R

end program mie_diff
