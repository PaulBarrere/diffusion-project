!> @file mod_param.f90
!! @author Paul Barrère
!! @brief Le fichier mod_param.f90 répertorie toutes les paramètres utilisés dans les subroutines et programmes principaux.

!> @author Paul Barrère
!! @brief Le module param répertorie toutes les paramètres utilisés dans les subroutines et programmes principaux.
module param
  implicit none
  
  integer, parameter :: ks=selected_real_kind(8) !< Simple précision
  integer, parameter :: kd=selected_real_kind(15) !< Double precision

  real(kind=kd), parameter :: pi=atan(1._kd)*4._kd !< Valeur de pi
  real(kind=kd), parameter :: eps=epsilon(1._kd) !< Valeur du plus petit réel

  real(kind=kd), parameter :: lambda_R=0.71_kd, lambda_G=0.51_kd ,lambda_B=0.44_kd !< Longeurs d'onde RGB
  real(kind=kd), parameter :: neau_R=1.3296_kd, neau_G=1.33535_kd, neau_B=1.34096 !< Partie réel de l'indice de réfraction de l'eau RGB
  real(kind=kd), parameter :: keau_R=5.5_kd*10._kd**(-8), keau_G=2._kd*10._kd**(-9), keau_B=1._kd*10._kd**(-9) !< Partie réel de l'indice de réfraction de l'eau RGB
  real(kind=kd), parameter :: nair_R=1.00027537_kd, nair_G=1.00027814_kd, nair_B=1.00028124_kd !< Indice de réfraction de l'air RGB
  real(kind=kd), parameter :: nN2_R=1.00025039_kd, nN2_G=1.00025383_kd, nN2_B=1.00025638_kd !< Indice de réfraction du diazote RGB
  real(kind=kd), parameter :: nO2_R=1.00029307_kd, nO2_G=1.00029604_kd, nO2_B=1.00029824_kd !< Indice de réfraction du dioxygène RGB

  integer, parameter :: N_phase=1997 ,N_part=5e6, N_diff_max=1e3 !< Nombres de phases, de photons, et de diffusion maximum pour un photon
  integer, parameter :: N_theta=1997, N_phi=1997 !< Nombres de phases suivant théta, et phi (utilisés pour la diffusion de Mie)
  integer, parameter :: N_pix=2048 !< Nombre de pixels pour l'image ppm
  real(kind=kd), parameter :: t_max=1._kd !< Profondeur optique du milieu par défaut
  real(kind=kd), parameter :: D=0.4_kd !< Paramètre caractérisant le milieu (utilisé pour normalisé le t_max du rouge et du vert avec le bleu
  real(kind=kd), parameter :: x_form=1._kd !< Facteur de forme par défaut
  real(kind=kd), parameter :: g=0.1_kd !< Paramètre associé à l'enveloppe Henyey-Greenstein par défaut
  real(kind=kd), parameter :: w=0.1_kd, n=3._kd !< Paramètres associés à l'enveloppe PolyN par défaut
  real(kind=kd), dimension(4), parameter :: S0=(/1._kd,0._kd,0._kd,0._kd/) !< Vecteur de Stokes initial par défaut
  logical, parameter :: polar=.false. !< Photon incident polarisé aléatoirement si .true., photon incident polarisé ou non avec comme
  !! vecteur de Stokes S0 
  

end module param
