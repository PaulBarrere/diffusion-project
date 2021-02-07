!> @file rayl_stereo_RGB.f90 
!! @author Paul Barrère
!! @brief Le programme rayl_stereo_RGB.f90 fait la projection stéréographique
!! de la couleur de la couche diffusante pour la diffusion de Rayleigh.
!! @param N_R,N_G,N_B Nombres de photons des couleurs RGB
!! @param phase_R,phase_G,phase_B Tableaux stockant la fonction de phase 
!! des couleurs RGB
!! @param theta Tableau stockant les angles de diffusion
program rayl_stereo_RGB
  use param
  use read_phase
  use build_stereo_map

  implicit none
  real(kind=kd) :: N_R, N_G, N_B
  real(kind=kd), dimension(N_theta,4) :: phase_R, phase_G, phase_B
  real(kind=kd), dimension(N_theta) :: theta

  open(33, file="Data/nb_photons_RGB_rayl.dat",status='old')
  read(33,*) N_R, N_G, N_B
  close(33)
  print *, N_R, N_G, N_B

  call read_phase_rayl(N_theta,'Data/rayl_diff_R.dat',phase_R,theta)
  call read_phase_rayl(N_theta,'Data/rayl_diff_G.dat',phase_G,theta)
  call read_phase_rayl(N_theta,'Data/rayl_diff_B.dat',phase_B,theta)

  phase_R=N_R*phase_R
  phase_G=N_G*phase_G
  phase_B=N_B*phase_B

  call stereo_map_RGB(phase_R,phase_G,phase_B,theta,"rayl_stereo_diff.ppm","rayl_stereo_rdiff.ppm")

end program rayl_stereo_RGB
