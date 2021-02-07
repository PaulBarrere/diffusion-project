!> @file mie_stereo_RGB.f90 
!! @author Paul Barrère
!! @brief Le programme mie_stereo_RGB.f90 fait la projection stéréographique
!! de la couleur de la couche diffusante pour la diffusion de Mie.
!! @param N_R,N_G,N_B Nombres de photons des couleurs RGB
!! @param phase_R,phase_G,phase_B Tableaux stockant la fonction de phase 
!! des couleurs RGB
!! @param theta Tableau stockant les angles de diffusion
program mie_stereo_RGB
  use param
  use read_phase
  use build_stereo_map

  implicit none
  real(kind=kd) :: N_R, N_G, N_B
  real(kind=kd), dimension(N_theta,4) :: phase_R, phase_G, phase_B
  real(kind=kd), dimension(N_theta) :: theta

  open(33, file="Data/nb_photons_RGB_mie.dat",status='old')
  read(33,*) N_R, N_G, N_B
  close(33)
  print *, N_R, N_G, N_B

  call read_phase_mie(N_theta,'Data/theta_mie_diff_R.dat',phase_R,theta)
  call read_phase_mie(N_theta,'Data/theta_mie_diff_G.dat',phase_G,theta)
  call read_phase_mie(N_theta,'Data/theta_mie_diff_B.dat',phase_B,theta)

  phase_R=N_R*phase_R
  phase_G=N_G*phase_G
  phase_B=N_B*phase_B

  call stereo_map_RGB(phase_R(:,1),phase_G(:,1),phase_B(:,1),theta,"Fig/mie_stereo_diff_10_2.ppm","Fig/mie_stereo_rdiff_10_2.ppm")

end program mie_stereo_RGB
