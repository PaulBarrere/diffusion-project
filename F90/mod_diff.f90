!> @file mod_diff.f90
!! @author Paul Barrère
!! @brief Le fichier mod_diff.f90 rassemble les modules utiles à simuler 
!! la diffusion de Rayleigh et de sortir les fonctions de phase 
!! servant d'enveloppe (Henyey-Greenstein et PolyN).

!> @author Paul Barrère
!! @brief Le module init_diff permet d'initialiser plusieurs quantités 
!! pour obtenir la fonction de phase.
module init_diff
  use param
  use random

  implicit none

  contains
    
    !> @brief La subroutine init initialise le tableau de phase,
    !! et les paramètres b et d_mu.
    !! @param [in] tau_max Profondeur optique du milieu
    !! @param [out] b Paramètre lié à la pronfondeur optique
    !! @param [out] d_mu Pas pour le cosinus de l'angle de diffusion
    !! @param [out] phase Tableau de la fonction de phase
    subroutine init(d_mu, b, tau_max, phase)
      implicit none
      real(kind=kd), intent(in) :: tau_max
      real(kind=kd), intent(out) :: d_mu, b
      real(kind=kd), allocatable, intent(out) :: phase(:)

      d_mu=2._kd*(1._kd+eps)/real(N_phase,kind=kd)
      b=1._kd-exp(-tau_max)

      allocate(phase(N_phase))
      phase=0._kd
    end subroutine init

    !> @brief La subroutine first_diff simule une première diffusion
    !! du photon, donc sort un premier angle de diffusion
    !! et une première profondeur optique conditionnelle.
    !! @param [in] b Paramètre lié à la pronfondeur optique
    !! @param [out] n_diff Nombre de diffusions
    !! @param [out] mu Cosinus de l'angle de diffusion
    !! @param [out] tau Profondeur optique conditionnelle
    subroutine first_diff(b, n_diff, mu, tau)
      implicit none
      real(kind=kd), intent(in) :: b
      integer, intent(out) :: n_diff
      real(kind=kd), intent(out) :: mu, tau
      real(kind=kd) :: U

      n_diff=1
      mu=1._kd

      call random_number(U)
      tau=-log(1._kd-b*U)
    end subroutine first_diff

end module init_diff
!******************************************************************************!

!******************************************************************************!

!> @author Paul Barrère
!> @brief Le module sortie_diff rassemble les subroutines utilisées 
!! en sortie de boucle à la fin d'une diffusion et à la fin des N_part 
!! simulations.
module sortie_diff
  use param

  implicit none

  contains

    !> @brief La subroutine update_path met à jour le chemin du photon
    !! après avoir été diffusé, donc sort une nouvelle profondeur 
    !! optique conditionnelle.
    !! @param [in] n_diff Nombre de diffusions déjà subies par le photon
    !! @param [inout] mu Cosinus de l'angle de diffusion
    !! @param [inout] tau Profondeur optique conditionnelle
    subroutine update_path(n_diff, mu, tau)
      implicit none
      real(kind=kd), intent(in) :: mu
      real(kind=kd) :: s
      real(kind=kd), intent(inout) :: tau
      integer, intent(inout) :: n_diff
      real(kind=kd) :: U

      call random_number(U)
      s=-log(U)
      tau=tau+s*mu
      n_diff=n_diff+1
    end subroutine update_path

    !> @brief La subroutine update_phase met à jour le tableau de la fonction
    !! de phase après que le photon soit sorti de la couche diffusante.
    !! @param [in] d_mu Pas pour le cosinus de l'angle de diffusion
    !! @param [in] mu Cosinus de l'angle de diffusion en sortie de la couche
    !! diffusante
    !! @param [inout] phase Tableau de la fonction de phase à mettre à jour
    !! @param [out] k Indice du tableau de la fonction de phase correspond 
    !! à l'angle de diffusion de sortie
    subroutine update_phase(d_mu, mu, k, phase)
      implicit none
      real(kind=kd), intent(in) :: d_mu, mu
      real(kind=kd), dimension(:), intent(inout) :: phase
      integer, intent(out) :: k

      k=floor((mu+1._kd)/d_mu)+1
      phase(k)=phase(k)+1
    end subroutine update_phase

    !> @brief Le subroutine save_res sauvegarde les résultats dans un fichier
    !! .dat en 3 colonnes: indice, phase, angle de diffusion.
    !! @param [in] num Indice
    !! @param [in] d_mu Pas pour le cosinus de l'angle de diffusion
    !! @param [inout] phase Tableau de la fonction de phase à mettre à jour
    subroutine save_res(num, d_mu, phase)
      implicit none
      integer, intent(in) :: num
      real(kind=kd), intent(in) :: d_mu
      real(kind=kd), dimension(:), intent(inout) :: phase
      real(kind=kd) :: rho, theta
      integer :: i

      do i=1,N_phase
        rho=(phase(i)*real(N_phase,kind=kd))/(2._kd*real(N_part,kind=kd))
        theta=dacos(d_mu*(real(i,kind=kd)-0.5_kd)-1._kd)
        write(num,*) i, rho, theta
      end do
    end subroutine save_res

    !> @brief Le subroutine save_res_bhmie sauvegarde la fonction de phase, 
    !! déduite de S1 et S2 de bhmie, dans un fichier .dat en 3 colonnes:
    !! indice, phase, angle de diffusion.
    !! @param [in] num Indice
    !! @param [inout] phase Tableau de la fonction de phase à mettre à jour
    subroutine save_res_bhmie(num, phase)
      implicit none
      integer, intent(in) :: num
      real(kind=kd), dimension(:), intent(inout) :: phase
      real(kind=kd) :: rho, theta
      integer :: i

      do i=1,N_phase
        rho=phase(i)
        theta=real((i-1),kind=kd)*pi/real(N_phase,kind=kd)
        write(num,*) i, rho, theta
      end do
    end subroutine save_res_bhmie


end module sortie_diff
!******************************************************************************!

!******************************************************************************!

!> @author Paul Barrère
!! @brief Le module simu_diff rassemble les subroutines utiles pour simuler 
!! la diffusion elle-même, c'est dire choisir un nouvel angle de diffusion.
module simu_diff
  use param
  use random
  use init_diff
  use sortie_diff
  implicit none

  contains
    
    !> @brief La subroutine choose_m_rayl calcule un angle polaire pour
    !! la diffusion de Rayleigh.
    !! @param [out] m Cosinus de l'angle polaire
    subroutine choose_m_rayl(m)
      implicit none
      real(kind=kd), intent(out) :: m
      real(kind=kd) :: U1, U2, U3

      call random_number(U1)
      call random_number(U2)
      call random_number(U3)

      if ((U2-U1)*(U3-U1)<0) then
        m=2._kd*U3-1._kd
      else 
        m=2._kd*U1-1._kd
      end if
    end subroutine choose_m_rayl
    
    !> @brief La subroutine choose_m_hg calcule un angle polaire pour
    !! l'enveloppe de Henyey-Greenstein.
    !! @param [in] g Paramètre de l'enveloppe de Henyey-Greenstein
    !! @param [out] m Cosinus de l'angle polaire
    subroutine choose_m_hg(g,m)
      implicit none
      real(kind=kd), intent(in) :: g
      real(kind=kd), intent(out) :: m
      real(kind=kd) :: c0, c1, c2, U, Z

      c0=0.5_kd*g*(3._kd-g*g)
      c1=1._kd+g*g
      c2=0.5_kd*g*(1._kd+g*g)

      call random_number(U)
      Z=2._kd*U-1._kd

      m=(c0+Z*(c1+c2*Z))/((1._kd+g*Z)*(1._kd+g*Z))
    end subroutine choose_m_hg

    !> @brief La subroutine choose_m_polyn1 calcule un angle polaire pour
    !! l'enveloppe PolyN.
    !! @param [in] n Degré du polynôme 
    !! @param [in] w Paramètre de l'enveloppe PolyN
    !! @param [out] m Cosinus de l'angle polaire
    subroutine choose_m_polyn1(n,w,m)
      implicit none
      real(kind=kd), intent(in) :: w
      real(kind=kd), intent(out) :: m
      integer, intent(in) :: n
      real(kind=kd), dimension(n+1) :: U
      integer :: i

      do i=1,n+1
        call random_number(U(i))
      end do

      if (U(n+1)<=w) then 
        m=2._kd*maxval(U)-1._kd
      else
        m=1._kd-2._kd*maxval(U)
      end if
    end subroutine choose_m_polyn1

    !> @brief La subroutine choose_m_polyn1 calcule un angle polaire pour
    !! l'enveloppe PolyN en incluant les degrés non-entiers.
    !! @param [in] n Degré du polynôme qui peut être non-entier
    !! @param [in] w Paramètre de l'enveloppe PolyN
    !! @param [out] m Cosinus de l'angle polaire
    subroutine choose_m_polyn2(n,w,m)
      implicit none
      real(kind=kd), intent(in) :: w, n
      real(kind=kd), intent(out) :: m
      real(kind=kd) :: U1, U2

      call random_number(U1)
      call random_number(U2)

      if (U2<=w) then
        m=2._kd*U1**(1._kd/n)-1._kd
      else
        m=1._kd-2._kd*U1**(1._kd/n)
      end if
    end subroutine choose_m_polyn2
    
    !> @brief La subroutine update_mu calcule le nouvel angle de diffusion
    !! depuis les angles polaire et azimutal selon l'enveloppe que l'on choisit.
    !! @param [in] diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
    !! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.
    !! @param [in] n Degré du polynôme dans le cas où diff_type='polyn'
    !! @param [in] param Paramètre associé à diff_type
    subroutine update_mu(diff_type,param,n,mu)
      implicit none
      character(len=*), intent(in) :: diff_type
      real(kind=kd), intent(in) :: n, param
      real(kind=kd), intent(out) :: mu
      real(kind=kd) :: m, beta
      real(kind=kd) :: U0

      ! Choix du m
      if (diff_type .eq. 'rayl') then
        call choose_m_rayl(m)
      elseif (diff_type .eq. 'hg') then
        call choose_m_hg(param,m)
      elseif (diff_type .eq. 'polyn') then
        call choose_m_polyn2(n,param,m)
      end if

      ! Calcul de mu
      call random_number(U0)
      beta=2._kd*pi*U0
      mu=mu*m+dsqrt((1._kd-mu*mu)*(1._kd-m*m))*dcos(beta)
    end subroutine update_mu

end module simu_diff 
!******************************************************************************!

!> @author Paul Barrère
!! @brief Le module all_diff rassemble 2 subroutines calculant 
!! la fonction de phase pour les 3 types d'enveloppes: une théorique et 
!! une via la méthode de Monte Carlo.
module all_diff
  use param
  use random
  use init_diff
  use sortie_diff
  use simu_diff
  implicit none

  contains

    subroutine get_phase_th(diff_type,param,d_mu,phase)
      implicit none
      character(len=*), intent(in) :: diff_type
      real(kind=kd), intent(in) :: param 
      real(kind=kd), intent(out) :: d_mu
      real(kind=kd), allocatable, intent(out) :: phase(:)
      real(kind=kd) :: mu
      integer :: i

      allocate(phase(N_phase))
      d_mu=pi/real(N_phase-1,kind=kd)
      do i=1,N_phase
        mu=dcos(real((i-1),kind=kd)*pi/real(N_phase-1,kind=kd))
        if (diff_type .eq. 'rayl') then
          phase(i)=3._kd*(1._kd+mu*mu)/8._kd
        else if (diff_type .eq. 'hg') then
          phase(i)=0.5_kd*(1._kd-param*param)/((1._kd-2._kd*param*mu+param*param)**(1.5_kd))
        else if (diff_type .eq. 'polyn') then
          phase(i)=3._kd*(param*(1._kd+mu)*(1._kd+mu)+(1._kd-param)*(1._kd-mu)*(1._kd-mu))/8._kd
        else 
          print *, diff_type
        end if
      end do
    end subroutine get_phase_th

    !> @brief La subroutine get_phase calcule la fonction de phase via 
    !! la méthode de Monte Carlo.
    !! @param [in] diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
    !! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.
    !! @param [in] param Paramètre associé à diff_type
    !! @param [in] n Degré du polynôme dans le cas où diff_type='polyn'
    !! @param [in] tau_max Profondeur optique du milieu
    !! @param [out] b Paramètre lié à tau_max
    !! @param [out] d_mu Pas pour le cosinus de l'angle de diffusion
    !! @param [out] phase Tableau de la fonction de phase associée à diff_type
    subroutine get_phase(diff_type,param,n,tau_max,b,d_mu,phase)
      implicit none
      character(len=*), intent(in) :: diff_type
      real(kind=kd), intent(in) :: tau_max, n, param
      integer :: n_diff, i, k
      real(kind=kd), intent(out) :: d_mu, b
      real(kind=kd) :: mu, t
      real(kind=kd), allocatable, intent(out) :: phase(:)

      ! Initialisation
      call init(d_mu, b, tau_max, phase)

      do i=1,N_part
        ! 1ère diffusion
        call first_diff(b, n_diff, mu, t)

        ! Diffusions suivantes
        do while (exp(-t)<=1._kd .and. exp(-t)>=(1._kd-b))
          ! Vérification du nombre de diffusion
          if (n_diff>N_diff_max) then
            print *, 'Erreur: photon ne sort pas'
            exit
          end if

          ! Calcul de mu
          call update_mu(diff_type,param,n,mu)

          ! Maj de tau
          call update_path(n_diff, mu, t)
        end do

        ! Sortie
        call update_phase(d_mu, mu, k, phase)
      end do
    end subroutine get_phase

end module all_diff
