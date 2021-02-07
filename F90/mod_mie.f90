!> @file mod_mie.f90
!! @author Paul Barrère
!! @brief Le fichier mod_mie.f90 rassemble les modules utiles à simuler 
!! la diffusion de Mie.

!> @author Paul Barrère
!! @brief Le module init_mie_diff permet d'initialiser plusieurs quantités 
!! pour obtenir la fonction de phase.
module init_mie_diff
  use param
  use random
  use mie

  implicit none

  contains

    !> @brief La subroutine generate_photon génère de manière aléatoire
    !! un photon polarisé.
    !! @param [out] S Vecteur de Stokes du photon polarisé généré. 
    subroutine generate_photon(S)
      implicit none
      real(kind=kd), dimension(4), intent(out) :: S
      real(kind=kd) :: U, Q, V
      real(kind=kd) :: sgn1, sgn2, sgn3
      
      ! Tirage aléatoire de paramètres de Stokes jusqu'à 
      ! satisfaction de la propriété U**2+Q**2+V**2<=I**2
      call random_number(U)
      call random_number(Q)
      call random_number(V)
      do while (U*U+Q*Q+V*V>1._kd)
        call random_number(U)
        call random_number(Q)
        call random_number(V)
      end do

      ! Tirage aléatoire du signe
      call random_number(sgn1)
      call random_number(sgn2)
      call random_number(sgn3)
      if (sgn1>0.5_kd) then
        U=-U
      end if
      if (sgn2>0.5_kd) then
        Q=-Q
      end if
      if (sgn3>0.5_kd) then
        V=-V
      end if

      S=(/1._kd,U,Q,V/)
    end subroutine generate_photon
    
    !> @brief La subroutine calcul_refrel calcule le paramètre refrel 
    !! utile pour la subroutine bhmie en fonction de la couleur et de
    !! l'élément diffusant.
    !! @param [in] color Couleur: 'red', 'green', ou 'blue'.
    !! @param [in] elem_diff Elément diffusant: 'O2', 'N2', ou 'H2O'
    !! @param [out] refrel Rapport entre l'indice de réfraction 
    !! de la particule diffusante et du milieu environnant
    subroutine calcul_refrel(color,elem_diff,refrel)
      implicit none
      character(len=*), intent(in) :: color, elem_diff
      complex(kind=kd), intent(out) :: refrel

      ! Molécule d'eau
      if (elem_diff .eq. 'H2O') then
        if (color .eq. 'red') then
          refrel=cmplx(neau_R,keau_R,kind=kd)/cmplx(nair_R,0._kd,kind=kd)
        else if (color .eq. 'green') then
          refrel=cmplx(neau_G,keau_G,kind=kd)/cmplx(nair_G,0._kd,kind=kd)
        else if (color .eq. 'blue') then
          refrel=cmplx(neau_B,keau_B,kind=kd)/cmplx(nair_B,0._kd,kind=kd)
        end if 
     
      ! Molécule de diazote
      else if (elem_diff .eq. 'N2') then
        if (color .eq. 'red') then
          refrel=cmplx(nN2_R,0._kd,kind=kd)/cmplx(nair_R,0._kd,kind=kd)
        else if (color .eq. 'green') then
          refrel=cmplx(nN2_G,0._kd,kind=kd)/cmplx(nair_G,0._kd,kind=kd)
        else if (color .eq. 'blue') then
          refrel=cmplx(nN2_B,0._kd,kind=kd)/cmplx(nair_B,0._kd,kind=kd)
        end if
 
      ! Molécule de dioxygène
      else if (elem_diff .eq. 'O2') then
        if (color .eq. 'red') then
          refrel=cmplx(nO2_R,0._kd,kind=kd)/cmplx(nair_R,0._kd,kind=kd)
        else if (color .eq. 'green') then
          refrel=cmplx(nO2_G,0._kd,kind=kd)/cmplx(nair_G,0._kd,kind=kd)
        else if (color .eq. 'blue') then
          refrel=cmplx(nO2_B,0._kd,kind=kd)/cmplx(nair_B,0._kd,kind=kd)
        end if
      end if 
    end subroutine calcul_refrel

    !> @brief La subroutine init_array_mie initialise les tableaux. 
    !! @param [out] phase_env Tableau stockant l'enveloppe 
    !! @param [out] phase_norm Tableau stockant le rapport entre 
    !! la fonction de phase de bhmie et celle de l'enveloppe
    !! @param [out] phase_theta Tableau stockant la fonction de phase finale
    !! en fonction de l'angle de diffusion theta
    !! @param [out] phase_phi Tableau stockant la fonction de phase finale
    !! en fonction de l'angle de diffusion phi
    !! @param [out] deg_pol Tableau stockant le degré de polarisation 
    !! en fonction de l'angle de diffusion theta
    !! @param [out] ang_pol Tableau stockant l'angle de polarisation 
    !! en fonction de l'angle de diffusion theta
    subroutine init_array_mie(phase_env,phase_norm,phase_theta,phase_phi,deg_pol,ang_pol)    
      implicit none
      real(kind=kd), allocatable, intent(out) :: phase_env(:), phase_norm(:), phase_theta(:,:), phase_phi(:,:), & 
      deg_pol(:), ang_pol(:)

      allocate(phase_env(N_phase)); allocate(phase_norm(N_phase)); allocate(deg_pol(N_theta)); allocate(ang_pol(N_theta))
      allocate(phase_theta(N_theta,4)); allocate(phase_phi(N_phi,4))
      phase_theta=0._kd; phase_phi=0._kd
    end subroutine init_array_mie

    !> @brief La subroutine init_var_mie initialise des quantités utiles.
    !! @param [in] tau_max Profondeur optique du milieu
    !! @param [out] d_alpha Pas pour l'angle alpha
    !! @param [out] d_theta Pas pour l'angle de diffusion theta
    !! @param [out] d_phi Pas pour l'angle de diffusion phi
    !! @param [out] b Paramètre lié à la pronfondeur optique
    !! @param [out] N_abs Nombre de photons absorbés
    !! @param [out] N_diff Nombre de photons diffusés
    !! @param [out] N_retro Nombre de photons rétrodiffusés
    subroutine init_var_mie(d_alpha,d_theta,d_phi,b,tau_max,N_abs,N_diff,N_retro)
      implicit none
      real(kind=kd), intent(in) :: tau_max
      real(kind=kd), intent(out) :: d_alpha, d_theta, d_phi, b
      integer, intent(out) :: N_abs, N_diff, N_retro

      N_abs=0; N_diff=0; N_retro=0

      d_alpha=pi/real(N_phase-1,kind=kd)
      d_theta=2._kd*(1._kd+eps)/N_theta
      d_phi=2._kd*pi/real(N_phi-1,kind=kd)

      b=1._kd-dexp(-tau_max)
    end subroutine init_var_mie

    !> @brief La subroutine init_bhmie initialise les quantités utiles pour 
    !! la subroutine bhmie.
    !! @param [in] color Couleur: 'red', 'green', ou 'blue'.
    !! @param [in] elem_diff Elément diffusant: 'O2', 'N2', ou 'H2O'
    !! @param [in] x Facteur de forme
    !! @param [out] albedo Probabilité que le photon soit diffusé 
    !! @param [out] refrel Rapport entre l'indice de réfraction 
    !! de la particule diffusante et du milieu environnant
    !! @param [out] S11,S12,S33,S34 Composantes de la matrice de Mueller
    !! @param [out] phase Tableau stockant la fonction de phase
    subroutine init_bhmie(color,elem_diff,x,albedo,refrel,S11,S12,S33,S34,phase)
      implicit none
      character(len=*), intent(in) :: color, elem_diff
      real(kind=kd), intent(in) :: x
      real(kind=kd), intent(out) :: albedo
      real(kind=kd), allocatable, intent(out) :: S11(:), S12(:), S33(:), S34(:), phase(:)
      complex(kind=kd), intent(out) :: refrel
      complex(kind=kd), dimension(N_phase) :: S1, S2
      real(kind=kd) :: Qext, Qsca, Qback, Gsca

      allocate(S11(N_phase)); allocate(S12(N_phase)); allocate(S33(N_phase)); allocate(S34(N_phase))
      allocate(phase(N_phase))

      call calcul_refrel(color,elem_diff,refrel)
      call bhmie(x, refrel, (N_phase+1)/2, S1, S2, Qext, Qsca, Qback, Gsca)

      S11=real(0.5_kd*(S2*conjg(S2)+S1*conjg(S1)),kind=kd)
      S12=real(0.5_kd*(S2*conjg(S2)-S1*conjg(S1)),kind=kd)
      S33=real(S2*conjg(S1),kind=kd)
      S34=aimag(S2*conjg(S1))

      albedo=Qsca/Qext

      phase=2._kd*S11/(x*x*Qsca)
      print *, "L'albedo est:", albedo
    end subroutine init_bhmie

    !> @brief La subroutine first_diff_mie initialise le vecteur d'onde 
    !! et l'axe de rotation du photon, et calcule une première diffusion.
    !! @param [in] b Paramètre lié à la pronfondeur optique
    !! @param [in] S Vecteur de Stokes initial
    !! @param [out] Ss Vecteur de Stokes en sortie
    !! @param [out] u Vecteur d'onde
    !! @param [out] p Axe de rotation
    !! @param [out] n_diff Nombre de diffusions déjà subies par le photon
    !! @param [out] tau Profondeur optique conditionnelle 
    subroutine first_diff_mie(b,S,Ss,u,p,n_diff,tau)
      implicit none
      real(kind=kd), intent(in) :: b
      real(kind=kd), dimension(4), intent(in) :: S
      integer, intent(out) :: n_diff
      real(kind=kd), intent(out) :: tau
      real(kind=kd), dimension(3), intent(out) :: u, p
      real(kind=kd), dimension(4), intent(out) :: Ss
      real(kind=kd) :: U0

      n_diff=1      
      u=(/0,1,0/)
      p=(/0,0,1/)
      Ss=S

      call random_number(U0)
      tau=-log(1._kd-b*U0)
    end subroutine first_diff_mie

end module init_mie_diff
!******************************************************************************!

!******************************************************************************!

!> @author Paul Barrère
!! @brief Le module sortie_mie_diff rassemble les subroutines utilisées 
!! en sortie de boucle à la fin d'une diffusion et à la fin des N_part 
!! simulations.
module sortie_mie_diff
  use param
  use random

  implicit none

  contains

    !> @brief La subroutine final_stokes sort un nouveau vecteur de Stokes.
    !! @param [in] u Axe de rotation
    !! @param [in] p Vecteur d'onde
    !! @param [out] theta, phi Angles de diffusion
    !! @param [out] Ss Vecteur de Stokes en sortie
    subroutine final_stokes(u,p,theta,phi,Ss)
      implicit none 
      real(kind=kd), dimension(3), intent(in) :: u, p
      real(kind=kd), intent(out) :: theta, phi
      real(kind=kd), dimension(4), intent(inout) :: Ss
      real(kind=kd) :: beta
      real(kind=kd), dimension(4,4) :: R_beta

      theta=dacos(p(3))
      phi=datan2(p(2),p(1))

      beta=datan2(u(3),p(1)*u(2)-p(2)*u(1))
      R_beta=0._kd
      R_beta(1,1)=1._kd; R_beta(4,4)=1._kd
      R_beta(2,2)=dcos(2._kd*beta); R_beta(3,3)=dcos(2._kd*beta)
      R_beta(2,3)=-dsin(2._kd*beta); R_beta(3,2)=dsin(2._kd*beta)

      Ss=matmul(R_beta,Ss)
    end subroutine final_stokes
      
    !> @brief La subroutine update_phases_mie met à jour les tableaux stockant 
    !! les fonction de phases.
    !! @param [in] d_theta Pas pour l'angle de diffusion theta
    !! @param [in] d_phi Pas pour l'angle de diffusion phi
    !! @param [in] theta,phi Angles de diffusion
    !! @param [in] Ss Vecteur de Stokes en sortie
    !! @param [out] phase_theta Tableau stockant la fonction de phase finale
    !! en fonction de l'angle de diffusion theta
    !! @param [out] phase_phi Tableau stockant la fonction de phase finale
    !! en fonction de l'angle de diffusion phi
    subroutine update_phases_mie(d_theta,d_phi,theta,phi,Ss,phase_theta,phase_phi)
      implicit none
      real(kind=kd), intent(in) :: d_theta, d_phi, theta, phi
      real(kind=kd), dimension(4),intent(in) :: Ss
      real(kind=kd), dimension(N_phase,4), intent(out) :: phase_theta, phase_phi
      integer :: k
      
      k=floor((dcos(theta)+1._kd)/d_theta)+1
      phase_theta(k,:)=phase_theta(k,:)+Ss
      phase_phi(floor(phi/d_phi)+floor(0.5_kd*(real(N_phi,kind=kd))+1._kd),:) &
      =phase_phi(floor(phi/d_phi)+floor(0.5_kd*(real(N_phi,kind=kd))+1._kd),:)+Ss
    end subroutine update_phases_mie

    !> @brief La subroutine save_res_mie sauvegarde la fonction de phase 
    !! finale dans un fichier .dat.
    !! @param [in] num Indice associé au fichier à ouvrir
    !! @param [in] d_mu Pas pour le cosinus de l'angle de diffusion
    !! @param [in] N Longueur du tableau
    !! @param [inout] phase Tableau stockant la fonction de phase
    subroutine save_res_mie(num,d_mu,N,phase)
      implicit none
      integer, intent(in) :: num, N
      real(kind=kd), intent(in) :: d_mu
      real(kind=kd), dimension(:,:), intent(inout) :: phase
      real(kind=kd), dimension (4) :: rho
      real(kind=kd) :: theta
      integer :: i

      do i=1,N
        rho=(phase(i,:)*real(N_phase,kind=kd))/(2._kd*real(N_part,kind=kd))
        theta=dacos(d_mu*(real(i,kind=kd)-0.5_kd)-1._kd)
        write(num,*) i, rho(1), rho(2), rho(3), rho(4), theta
      end do
    end subroutine save_res_mie

end module sortie_mie_diff
!******************************************************************************!

!******************************************************************************!

!> @author Paul Barrère
!! @brief Le module simu_mie_diff rassemble les subroutines utiles pour simuler 
!! la diffusion elle-même, c'est dire choisir les nouveaux angles de diffusion
!! et un nouveau vecteur de Stokes.
module simu_mie_diff
  use param
  use random
  use simu_diff
  use fzero
  use interpolation

  implicit none

  contains

    !> @brief La subroutine calculate_alpha_mie calcule un angle 
    !! de diffusion. Il est à noter la partie interpollant 'phase' 
    !! en utilisant le module 'interpolation' ralentit (beaucoup) 
    !! les calculs donc on peut utiliser une interpolation linéaire 
    !! pour aller plus vite.
    !! @param [in] diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
    !! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.
    !! @param [in] param Paramètre associé à diff_type
    !! @param [inout] phase Tableau stockant le rapport entre 
    !! la fonction de phase de bhmie et celle de l'enveloppe
    !! @param [in] d_alpha Pas pour l'angle de diffusion theta
    !! @param [out] alpha Angle de diffusion
    !! @parma [out] P_alpha Evaluation de la fonction de phase en alpha
    !! @param [out] Réel entre 0 et 1 tiré aléatoirement
    subroutine calculate_alpha_mie(diff_type,param,phase,d_alpha,alpha,P_alpha,U)
      implicit none
      character(len=*), intent(in) :: diff_type
      real(kind=kd), intent(in) :: param, d_alpha
      real(kind=kd), dimension(*), intent(inout) :: phase
      real(kind=kd), intent(out) :: alpha, P_alpha, U
      real(kind=kd) :: mu
      
      ! Variables liées au lissage de la fonction de phase
      integer, parameter :: deg=21
      integer :: i
      real(kind=kd), dimension(deg+1) :: coeff
      real(kind=kd), dimension(N_phase) :: tab_alpha

      ! Choix du mu
      if (diff_type .eq. 'rayl') then
        call choose_m_rayl(mu)
      elseif (diff_type .eq. 'hg') then
        call choose_m_hg(param,mu)
      elseif (diff_type .eq. 'polyn') then
        call choose_m_polyn2(3._kd,param,mu)
      end if

      alpha=dacos(mu)

      ! Interpolation de phase via le module 'interpolation'
      !tab_alpha=real((/(i,i=0,N_phase-1)/),kind=kd)*d_alpha
      !call lisse(deg,tab_alpha,phase,coeff)
      !call interpol(deg,coeff,alpha,P_alpha)

      ! Interpolation linéaire de 'phase'
      P_alpha=(phase(floor(alpha/d_alpha)+1)+phase(floor(alpha/d_alpha)+2))/2._kd
      
      call random_number(U)
    end subroutine calculate_alpha_mie

    !> @brief La subroutine choose_alpha choisit un angle de diffusion 
    !! via la méthode de réjection de von Neumann.
    !! @param [in] diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
    !! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.
    !! @param [in] param Paramètre associé à diff_type
    !! @param [inout] phase Tableau stockant le rapport entre 
    !! la fonction de phase de bhmie et celle de l'enveloppe
    !! @param [in] d_alpha Pas pour l'angle de diffusion theta
    !! @param [out] alpha Angle de diffusion
    subroutine choose_alpha_mie(diff_type,param,phase,d_alpha,alpha)
      implicit none
      character(len=*), intent(in) :: diff_type
      real(kind=kd), intent(in) :: param, d_alpha
      real(kind=kd), dimension(*), intent(inout) :: phase
      real(kind=kd), intent(out) :: alpha
      real(kind=kd) :: P_alpha, U

      call calculate_alpha_mie(diff_type,param,phase,d_alpha,alpha,P_alpha,U)

      do while (U>P_alpha)
        call calculate_alpha_mie(diff_type,param,phase,d_alpha,alpha,P_alpha,U)
      end do
    end subroutine choose_alpha_mie

    !> @brief La subroutine choose_beta_mie calcule l'angle azimuthal 
    !! et évalue les composante S11 et S12 en alpha.
    !! Il est à noter la partie interpollant S11 et S12 en utilisant le 
    !! module 'interpolation' ralentit (beaucoup) les calculs donc on peut
    !! utiliser une interpolation linéaire pour aller plus vite.
    !! @param [in] alpha Anlge de diffusion
    !! @param [in] d_alpha Pas pour l'angle de diffusion
    !! @param [in] Si Vecteur de Stokes précédent à Ss
    !! @param [out] S11,S12 Composantes de la matrice de Mueller
    !! @param [out] S11_alpha,S12_alpha Evaluation de S11 et S12 en alpha
    !! @param [out] beta Angle azimuthal
    subroutine choose_beta_mie(alpha,d_alpha,Si,S11,S12,S11_alpha,S12_alpha,beta)
      implicit none
      real(kind=kd), intent(in) :: alpha, d_alpha
      real(kind=kd), dimension(4), intent(in) :: Si
      real(kind=kd), dimension(*), intent(inout) :: S11, S12
      real(kind=kd), intent(out) :: S11_alpha, S12_alpha, beta
      real(kind=kd) :: P_lin, A, B, chi, U

      ! Variables liées au lissage de la fonction de phase
      integer, parameter :: deg=21
      integer :: i
      real(kind=kd), dimension(deg+1) :: coeff_S11, coeff_S12
      real(kind=kd), dimension(N_phase) :: tab_alpha

      ! Interpolation de S11 et S12 via le module 'interpolation'
      !tab_alpha=real((/(i,i=0,N_phase-1)/),kind=kd)*d_alpha
      !call lisse(deg,tab_alpha,S11,coeff_S11)
      !call lisse(deg,tab_alpha,S12,coeff_S12)
      !call interpol(deg,coeff_S11,alpha,S11_alpha)
      !call interpol(deg,coeff_S12,alpha,S12_alpha)

      ! Interpolation linéaire de S11 et S12
      S11_alpha=(S11(floor(alpha/d_alpha)+1)+S11(floor(alpha/d_alpha)+2))/2._kd
      S12_alpha=(S12(floor(alpha/d_alpha)+1)+S12(floor(alpha/d_alpha)+2))/2._kd
      P_lin=dsqrt(Si(2)*Si(2)+Si(3)*Si(3))

      if (P_lin .eq. 0._kd) then
        chi=0._kd 
      else
        chi=0.5_kd*atan2(Si(3),Si(2))
      end if

      call random_number(U)
      A=0.5_kd*P_lin*S12_alpha/S11_alpha
      B=A*dsin(2._kd*chi)-2._kd*pi*U
      if (U<=0.5_kd) then
        call dicho(A,B,chi,0._kd,pi,10._kd**(-8),beta)
      else 
        call dicho(A,B,chi,pi,2._kd*pi,10._kd**(-8),beta)
      end if
    end subroutine choose_beta_mie

    !> @brief La subroutine update_stokes met à jour le vecteur de Stokes.
    !! Il est à noter la partie interpollant S33 et S34 en utilisant le 
    !! module 'interpolation' ralentit (beaucoup) les calculs donc on peut
    !! utiliser une interpolation linéaire pour aller plus vite.
    !! @param [in] diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
    !! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.
    !! @param [in] param Paramètre associé à diff_type
    !! @param [inout] phase Tableau stockant le rapport entre 
    !! la fonction de phase de bhmie et celle de l'enveloppe
    !! @param [in] d_alpha Pas pour l'angle de diffusion theta
    !! @param [inout] S11,S12,S33,S34 Composantes de la matrice de Mueller
    !! @param [in] Si Vecteur de Stokes précédent à Ss
    !! @param [out] alpha Angle de diffusion
    !! @param [out] beta Angle azimuthal
    !! @param [out] Vecteur de Stokes en sortie
    subroutine update_stokes(diff_type,param,phase,d_alpha,S11,S12,S33,S34,Si,alpha,beta,Ss)
      implicit none
      character(len=*), intent(in) :: diff_type
      real(kind=kd), intent(in) :: param, d_alpha
      real(kind=kd), dimension(*), intent(inout) :: S11, S12, S33, S34, phase
      real(kind=kd), dimension(4), intent(in) :: Si
      real(kind=kd), dimension(4), intent(out) :: Ss
      real(kind=kd) :: S11_alpha, S12_alpha, S33_alpha, S34_alpha
      real(kind=kd), intent(out) :: alpha, beta
      real(kind=kd), dimension(4,4) :: M_alpha, R_beta

      ! Variables liées au lissage de la fonction de phase
      integer, parameter :: deg=21
      integer :: i
      real(kind=kd), dimension(deg+1) :: coeff_S33, coeff_S34
      real(kind=kd), dimension(N_phase) :: tab_alpha
      
      ! Interpolation de S33 et S34 via le module 'interpolation'
      !tab_alpha=real((/(i,i=0,N_phase-1)/),kind=kd)*d_alpha
      !call lisse(deg,tab_alpha,S33,coeff_S33)
      !call lisse(deg,tab_alpha,S34,coeff_S34)
      !call interpol(deg,coeff_S33,alpha,S33_alpha)
      !call interpol(deg,coeff_S34,alpha,S34_alpha)

      ! Interpolation linéaire de S33 et S34
      S33_alpha=(S33(floor(alpha/d_alpha)+1)+S33(floor(alpha/d_alpha)+2))/2._kd
      S34_alpha=(S34(floor(alpha/d_alpha)+1)+S34(floor(alpha/d_alpha)+2))/2._kd

      ! Calcul des angles alpha et beta
      call choose_alpha_mie(diff_type,param,phase,d_alpha,alpha)
      call choose_beta_mie(alpha,d_alpha,Si,S11,S12,S11_alpha,S12_alpha,beta)

      ! Construction de la matrice de Mueller
      M_alpha=0._kd
      M_alpha(1,1)=S11_alpha; M_alpha(2,2)=S11_alpha
      M_alpha(1,2)=S12_alpha; M_alpha(2,1)=S12_alpha
      M_alpha(3,3)=S33_alpha; M_alpha(4,4)=S33_alpha
      M_alpha(3,4)=S34_alpha; M_alpha(4,3)=-S34_alpha

      ! Construction de la matrice de rotation
      R_beta=0._kd
      R_beta(1,1)=1._kd; R_beta(4,4)=1._kd
      R_beta(2,2)=dcos(2._kd*beta); R_beta(3,3)=dcos(2._kd*beta)
      R_beta(2,3)=dsin(2._kd*beta); R_beta(3,2)=-dsin(2._kd*beta)

      ! Reconstruction du vecteur de Stokes
      Ss=matmul(M_alpha,matmul(R_beta,Si))
      Ss=Ss/Ss(1)
    end subroutine update_stokes

    !> @brief La subroutine update_p_u met à jour le vecteur d'onde et
    !! l'axe de rotation.
    !! @param [in] alpha Angle de diffusion
    !! @param [in] beta Angle azimuthal
    !! @parma [in] u0,p0 Vecteurs permettant la mise à jour de p et u
    !! @param [out] u Axe de rotation
    !! @parma [out] p Vecteur d'onde
    subroutine update_p_u(alpha,beta,u0,p0,u,p)
      implicit none
      real(kind=kd), intent(in) :: alpha, beta
      real(kind=kd), dimension(3), intent(in) :: u0, p0
      real(kind=kd), dimension(3), intent(out) :: u, p

      u(1)=u0(1)*dcos(beta)+(p0(2)*u0(3)-p0(3)*u0(2))*dsin(beta)
      u(2)=u0(2)*dcos(beta)+(p0(3)*u0(1)-p0(1)*u0(3))*dsin(beta)
      u(3)=u0(3)*dcos(beta)+(p0(1)*u0(2)-p0(2)*u0(1))*dsin(beta)

      u=u/dsqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))

      p(1)=p0(1)*dcos(alpha)+(u(2)*p0(3)-u(3)*p0(2))*dsin(alpha)
      p(2)=p0(2)*dcos(alpha)+(u(3)*p0(1)-u(1)*p0(3))*dsin(alpha)
      p(3)=p0(3)*dcos(alpha)+(u(1)*p0(2)-u(2)*p0(1))*dsin(alpha)

      p=p/dsqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
    end subroutine update_p_u

    !> @brief La subroutine new_path met à jour le nouveau chemin du photon diffusé.
    !! @param [in] diff_type Type de l'enveloppe: 'rayl' pour Rayleigh, 
    !! 'hg' pour Henyey-Greenstein, et 'polyn' pour PolyN.
    !! @param [inout] albedo Probabilité que le photon soit diffusé
    !! @param [in] para Paramètre associé à diff_type
    !! @param [inout] phase_norm Tableau stockant le rapport entre 
    !! la fonction de phase de bhmie et celle de l'enveloppe
    !! @param [in] d_alpha Pas pour l'angle de diffusion
    !! @param [inout] S11,S12,S33,S34 Composantes de la matrice de Mueller
    !! @param [inout] Ss Vecteur de Stokes en sortie
    !! @param [inout] u Axe de rotation
    !! @param [inout] p Vecteur d'onde
    !! @param [inout] nb_diff Nombre de diffusions déjè subies par le photon
    !! @param [inout] t Profondeur optique conditionnelle
    !! @param [inout] rd Réel entre 0 et 1 tiré aléatoirement pour vérifier
    !! s'il y a diffusion
    subroutine new_path(diff_type,albedo,para,phase_norm,d_alpha,S11,S12,S33,S34,Ss,u,p,nb_diff,t,rd)
      implicit none
      character(len=*), intent(in) :: diff_type
      real(kind=kd), dimension(*), intent(inout) :: phase_norm, S11, S12, S33, S34
      real(kind=kd), intent(in) :: d_alpha, para
      real(kind=kd), dimension(4), intent(inout) :: Ss
      real(kind=kd), dimension(3), intent(inout) :: u, p
      real(kind=kd), intent(inout) :: t, rd, albedo
      integer, intent(inout) :: nb_diff

      real(kind=kd), dimension(4) :: Si
      real(kind=kd), dimension(3) :: u0, p0
      real(kind=kd) :: alpha, beta

      ! Absorption ou diffusion
      call random_number(rd)
      if ((rd<=albedo) .and. (nb_diff<=N_diff_max)) then
        p0=p; u0=u
        Si=Ss

        ! Calcul du nouveau vecteur de Stokes
        call update_stokes(diff_type,para,phase_norm,d_alpha,S11,S12,S33,S34,Si,alpha,beta,Ss)

        ! Calcul du nouveau vecteur d'onde
        call update_p_u(alpha,beta,u0,p0,u,p)
      
        ! Mise à jour de tau
        call update_path(nb_diff, p(3), t)
      end if
      
    end subroutine new_path

end module simu_mie_diff

