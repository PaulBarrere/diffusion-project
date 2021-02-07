!> @file mod_enveloppe.f90
!! @author Paul Barrère
!! @brief Le fichier mod_enveloppe.f90 contient 1 routine, 
!! qui rassemble les subroutines utiles pour déterminer 
!! la meilleure enveloppe et pour effectuer la méthode par réjection de 
!! von Neumann.

!> @author Paul Barrère
!! @brief Le module env_mod rassemble les subroutines utiles pour déterminer 
!! la meilleure enveloppe et pour effectuer la méthode par réjection de 
!! von Neumann.
module env_mod
  use param
  use random
  use mie
  use all_diff
  use interpolation 

  implicit none

  contains

    !> @brief La subroutine efficacite calcule l'efficacité de la méthode.
    !! @param [in] x Facteur de forme
    !! @param [in] refrel Rapport entre l'indice de réfraction 
    !! de la particule diffusante et du milieu environnant
    !! @param [in] speed Vitesse pour calculer l'enveloppe: 'fast' 
    !! via la formule théorique, et 'slow' via la méthode de Monte Carlo
    !! @param [in] diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
    !! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.
    !! @param [in] param Paramètre associé à diff_type
    !! @param [in] n Degré du polynôme dans le cas où diff_type='polyn'
    !! @param [in] tau_max Profondeur optique du milieu
    !! @param [out] eff Efficacité
    subroutine efficacite(x,refrel,speed,diff_type,param,n,tau_max,eff)
      implicit none
      character(len=*), intent(in) :: speed, diff_type
      real(kind=kd), intent(in) :: x, tau_max, n, param
      complex(kind=kd), intent(in) :: refrel
      real(kind=kd), intent(out) :: eff

      ! Variables liées à BHMIE
      complex(kind=kd), dimension(N_phase) :: S1, S2
      real(kind=kd), dimension(N_phase) :: I_s
      real(kind=kd) :: Qext, Qsca, Qback, Gsca

      ! Variables liées à l'enveloppe
      real(kind=kd) :: d_mu, b
      real(kind=kd), allocatable :: phase_env(:)
      
      ! Variables liées au lissage de la fonction de phase
      integer, parameter :: deg=21
      integer :: i
      real(kind=kd), dimension(deg+1) :: coeff
      real(kind=kd), dimension(N_phase) :: theta
      
      ! Calcul de la fonction de phase à laquelle on doit trouver une enveloppe
      call bhmie(x, refrel, (N_phase+1)/2, S1, S2, Qext, Qsca, Qback, Gsca)
      I_s=2._kd*real(0.5_kd*(S2*conjg(S2)+S1*conjg(S1)),kind=kd)/(x*x*Qsca)

      ! Calcul de l'enveloppe par méthode de Monte Carlo
      if (speed .eq. 'slow') then
        call get_phase(diff_type,param,n,tau_max,b,d_mu,phase_env)
        phase_env=real(N_phase,kind=kd)*phase_env/(2._kd*real(N_part,kind=kd))
        theta=dacos(d_mu*(real((/(i,i=1,N_phase)/),kind=kd)-0.5_kd)-1._kd)
        call lisse(deg,theta,phase_env,coeff)
        do i=1,N_phase
          call interpol(deg,coeff,dcos(theta(i)),phase_env(i))
        end do

      ! Calcul de l'enveloppe via la formule théorique de la fonction
      ! de phase
      else if (speed .eq. 'fast') then
        call get_phase_th(diff_type,param,d_mu,phase_env) 
      else
        print *, "Speed ne peut qu'être 'fast' ou 'slow'."
      end if 
      ! Calcul de l'efficacité
      eff=1._kd/maxval(I_s/phase_env)
    end subroutine efficacite

    !> @brief La subroutine best_param_g détermine la valeur du paramètre
    !! g (-1<g<1) maximisant l'efficacité de l'enveloppe de Henyey-Greenstein.
    !! @param [in] x Facteur de forme
    !! @param [in] step Pas entre 2 valeurs de g à tester
    !! @param [in] speed Vitesse pour calculer l'enveloppe: 'fast' 
    !! via la formule théorique, et 'slow' via la méthode de Monte Carlo
    !! @param [in] refrel Rapport entre l'indice de réfraction 
    !! de la particule diffusante et du milieu environnant
    !! @param [out] best_g Paramètre g maximisant l'efficacité
    !! @param [out] best_eff Valeur de l'efficacité associée à best_g
    subroutine best_param_g(x,step,speed,refrel,best_g,best_eff)
      implicit none
      real(kind=kd), intent(in) :: x
      integer, intent(in) :: step
      character(len=*), intent(in) :: speed
      complex(kind=kd), intent(in) :: refrel
      real(kind=kd), intent(out) :: best_g, best_eff
      real(kind=kd) :: g, eff
      integer :: j

      best_g=1._kd/real(step/2,kind=kd)-1._kd
      call efficacite(x,refrel,speed,'hg',best_g,0._kd,0._kd,best_eff)

      do j=2,step-1
        g=real(j,kind=kd)/real(step/2,kind=kd)-1._kd
        call efficacite(x,refrel,speed,'hg',g,0._kd,0._kd,eff)
        if (best_eff<eff) then
           best_eff=eff
           best_g=g
        end if
      end do
    end subroutine best_param_g

    !> @brief La subroutine best_param_w détermine la valeur du paramètre
    !! w (0<w<1) maximisant l'efficacité de l'enveloppe PolyN.
    !! @param [in] x Facteur de forme
    !! @param [in] n Degré du polynôme  
    !! @param [in] step Pas entre 2 valeurs de w à tester
    !! @param [in] speed Vitesse pour calculer l'enveloppe: 'fast' 
    !! via la formule théorique, et 'slow' via la méthode de Monte Carlo
    !! @param [in] refrel Rapport entre l'indice de réfraction 
    !! de la particule diffusante et du milieu environnant
    !! @param [out] best_w Paramètre w maximisant l'efficacité
    !! @param [out] best_eff Valeur de l'efficacité associée à best_w
    subroutine best_param_w(x,n,step,speed,refrel,best_w,best_eff)
      implicit none
      real(kind=kd), intent(in) :: x, n
      integer, intent(in) :: step
      character(len=*), intent(in) :: speed
      complex(kind=kd), intent(in) :: refrel
      real(kind=kd), intent(out) :: best_w, best_eff
      real(kind=kd) :: w, eff
      integer ::j

      best_w=0._kd
      call efficacite(x,refrel,speed,'polyn',best_w,n,0._kd,best_eff)

      do j=1,step
        w=real(j,kind=kd)/real(step,kind=kd)
        call efficacite(x,refrel,speed,'polyn',w,n,0._kd,eff)
        if (best_eff<eff) then
           best_eff=eff
           best_w=w
        end if
      end do
    end subroutine best_param_w

    !> @brief La subroutine best_param_w détermine la valeur du paramètre
    !! w (0<w<1) maximisant l'efficacité de l'enveloppe PolyN.
    !! @param [in] x Facteur de forme
    !! @param [in] step Pas entre 2 valeurs de g et de w à tester
    !! @param [in] speed Vitesse pour calculer l'enveloppe: 'fast' 
    !! via la formule théorique, et 'slow' via la méthode de Monte Carlo
    !! @param [in] refrel Rapport entre l'indice de réfraction 
    !! de la particule diffusante et du milieu environnant
    !! @param [out] eff Efficacité de la meilleure enveloppe
    !! @param [out] param Paramètre associé à la meilleure enveloppe
    !! @param [out] diff_type Type de la meilleure enveloppe: 
    !! 'rayl' pour Rayleigh, 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.  
    subroutine best_env(x,step,speed,refrel,eff,param,diff_type)
      implicit none
      real(kind=kd), intent(in) :: x
      integer, intent(in) :: step
      character(len=*), intent(in) :: speed
      complex(kind=kd), intent(in) :: refrel
      real(kind=kd), intent(out) :: eff, param
      character(len=*), intent(out) :: diff_type
      real(kind=kd) :: eff_rayl, eff_hg, eff_polyn, best_g, best_w
      
      call efficacite(x,refrel,speed,'rayl',0._kd,0._kd,0._kd,eff_rayl)
      call best_param_g(x,step,speed,refrel,best_g,eff_hg)
      call best_param_w(x,3._kd,step,speed,refrel,best_w,eff_polyn)
      print *,"L'efficacité de Rayleigh", eff_rayl
      print *,"L'efficacité de Henyey-Greenstein", eff_hg
      print *,"L'efficacité du PolyN", eff_polyn
      if (eff_rayl>eff_hg .and. eff_rayl>=eff_polyn) then
        param=0._kd
        diff_type='rayl'
        eff=eff_rayl
      else if (eff_hg>eff_rayl .and. eff_hg>eff_polyn) then
        param=best_g
        diff_type='hg'
        eff=eff_hg
      else
        param=best_w
        diff_type='polyn'
        eff=eff_polyn
      end if
    end subroutine best_env

end module env_mod

