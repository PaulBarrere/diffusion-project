!> @file mod_iso.f90
!! @author Paul Barrère
!! @brief Le fichier mod_iso.f90 contient un module utile seulement 
!! pour la diffusion isotrope.

!> @author Paul Barrère
!! @brief Le module iso_mod contient un module utile seulement
!! pour la diffusion isotrope.
module iso_mod
  use param
  use random

  implicit none

  contains

    !> @brief La subroutine update_mu_iso met à jour l'angle de diffusion. 
    !! @param [out] mu Cosinus de l'angle de diffusion
    subroutine update_mu_iso(mu)
      implicit none
      real(kind=kd), intent(out) :: mu
      real(kind=kd) :: U

      call random_number(U)
      mu=2._kd*U-1._kd
    end subroutine update_mu_iso

end module iso_mod
