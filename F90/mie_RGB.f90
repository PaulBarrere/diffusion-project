!> @file mie_RGB.f90
!! @author Paul Barrère
!! @brief Le programme mie_RGB.f90 simule la diffusion de Mie pour les couleurs RGB.
!! @param d_alpha Pas de l'angle de diffusion pour l'enveloppe et bhmie
!! @param d_mu Pas du cosinus de l'angle de diffusion alpha
!! @param d_theta Pas de l'angle de diffusion pour la diffusion de Mie
!! @param d_phi Pas de l'angle azimuthal
!! @param S Vecteur de Stokes initial
!! @param Ss Vecteur de Stokes en sortie 
!! @param tau_max Profondeur optique
!! @param l Libre parcours moyen
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

program mie_RGB
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
  use proj_stereo

  implicit none
  character(len=*), parameter :: colors(3) = (/character(1024) :: "red", "green", "blue"/)
  character(len=*), parameter :: files_theta(3) = (/character(1024) :: "Data/theta_mie_diff_R.dat" &
  , "Data/theta_mie_diff_G.dat", "Data/theta_mie_diff_B.dat"/)
  character(len=*), parameter :: files_phi(3) = (/character(1024) :: "Data/phi_mie_diff_R.dat" &
  , "Data/phi_mie_diff_G.dat", "Data/phi_mie_diff_B.dat"/)
  real(kind=kd), allocatable :: l(:), tau_max(:), b(:)
  
  real(kind=kd) :: d_alpha, d_theta, d_phi, d_mu
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
  integer :: i, nb_diff, N

  ! Mise en place système RBG
  allocate(l(3)); allocate(tau_max(3)); allocate(b(3))
  l=(/(lambda_R/lambda_B)**4, (lambda_G/lambda_B)**4, 1._kd/)
  tau_max=D/l

  ! Seed
  !call seed_file('seed.txt')
  call put_seed('seed.txt')

  ! Affichage de quelques paramètres
  print *, "Nombres de phases: ", N_phase, N_theta, N_phi
  print *, "Nombre de particules :", N_part
  print *, "Nombre de diffusions max :" , N_diff_max
  print *, "Chemin optique max du rouge :", tau_max(1)
  print *, "Chemin optique max du vert :", tau_max(2)
  print *, "Chemin optique max du bleu :", tau_max(3)
  print *, "Facteur de forme :", x_form
  print *, " "
  
  do i=1,3
    ! Initialisation des variables et arrays diffusion mie
    call init_array_mie(phase_env,phase_norm,phase_theta,phase_phi,deg_pol,ang_pol)
    call init_var_mie(d_alpha,d_theta,d_phi,b(i),tau_max(i),N_abs,N_diff,N_retro)
  
    ! Initialisation des variables de bhmie et de la matrice de Mueller
    call init_bhmie(colors(i),'N2',x_form,albedo,refrel,S11,S12,S33,S34,phase_bhmie)
    print *, "Initialisation finie"

    ! Calcul de la meilleure enveloppe
    call best_env(x_form,200,'fast',refrel,eff,para,diff_type)
    call get_phase_th(diff_type,para,d_mu,phase_env)
    phase_norm=eff*phase_bhmie/phase_env
    print *, "La meilleure enveloppe est: ", trim(diff_type)

    ! Ouverture des fichiers sauvegardant les vecteurs de Stokes
    open(20+i,file=files_theta(i),status='unknown')
    open(23+i,file=files_phi(i),status='unknown')
    print *, "Albedo: ",albedo
  
    ! Au cas où polar=false, on a une lumière non polarisée
    S=S0

    do N=1,N_part
      ! Génération d'un vecteur de Stokes initial
      if (polar) then
        call generate_photon(S)
      end if  

      ! Mise en place de la 1ère Diffusion
      call first_diff_mie(b(i),S,Ss,u,p,nb_diff,t)
        
      ! Diffusions suivantes
      do while (exp(-t)<=1._kd .and. exp(-t)>=(1._kd-b(i)))
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
    call save_res_mie(20+i,d_theta,N_theta,phase_theta)
    call save_res_mie(23+i,d_phi,N_phi,phase_phi)

    close(20+i)
    close(23+i) 
    print *, "Sauvegarde faite"

    ! Calcul de la réflectance
    R=b(i)*real(N_retro,kind=kd)/real(N_part,kind=kd)
    deg_pol=dsqrt(phase_theta(:,2)**2+phase_theta(:,3)**2+phase_theta(:,4)**2)/phase_theta(:,1)
    
    print *, 'Nombre de photons diffusés: ',N_diff
    print *, 'Nombre de photons rétrodiffusés: ',N_retro 
    print *, 'Nombre de photons absorbés: ', N_abs
    print *, 'Réflectance du milieu: ',R
    print *, " "
    print *, " "
  end do

  ! Calcul et sauvegarde du nombre de photons pour chaque couleur
  print *, "Nombre de photons rouges diffusés : ", real(N_part,kind=kd)*b(1)/b(3)
  print *, "Nombre de photons verts diffusés : ", real(N_part,kind=kd)*b(2)/b(3)
  print *, "Nombre de photons bleus diffusés : ", real(N_part,kind=kd)
  
  open(33, file="Data/nb_photons_RGB_mie.dat")
  write(33,*) real(N_part,kind=kd)*b(1)/b(3), real(N_part,kind=kd)*b(2)/b(3), real(N_part,kind=kd)
  close(33)
end program mie_RGB
