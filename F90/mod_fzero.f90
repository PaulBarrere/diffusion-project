!> @file mod_fzero.f90
!! @author Paul Barrère
!! @brief Le fichier mod_fzero.f90 rassemble 2 modules pour calculer 
!! un zéro d'une fonction.

!> @author Paul Barrère
!! @brief Le module func évalue la valeur de la fonction 
!! f: beta |--> beta+A*sin(2*(beta-chi))+B et de sa dérivée.
module func
  use param

  implicit none

  contains

    !> @brief La subroutine f évalue la fonction
    !! f: beta |--> beta+A*sin(2*(beta-chi))+B en un beta donné.
    !! @param [in] A Paramètre lié au degré de polarisation linéaire
    !! @param [in] B Paramètre dépendant de A et de l'angle de polarisation chi
    !! @param [in] chi Angle de polarisation
    !! @param [in] beta Angle azimuthal de diffusion
    !! @param [out] y y=f(beta)
    subroutine f(A,B,chi,beta,y) 
      implicit none
      real(kind=kd), intent(in) :: A, B, chi, beta
      real(kind=kd), intent(out) :: y

      y=beta+A*dsin(2._kd*(beta-chi))+B
    end subroutine f

    !> @brief La subroutine f_prime évalue la dérivée de la fonction
    !! f: beta |--> beta+A*sin(2*(beta-chi))+B en un beta donné.
    !! @param [in] A Paramètre lié au degré de polarisation linéaire
    !! @param [in] chi Angle de polarisation
    !! @param [in] beta Angle azimuthal de diffusion
    !! @param [out] y y=f'(beta)
    subroutine f_prime(A,chi,beta,y) 
      implicit none
      real(kind=kd), intent(in) :: A, chi, beta
      real(kind=kd), intent(out) :: y

      y=1._kd+2._kd*A*dcos(2._kd*(beta-chi))
    end subroutine f_prime

end module func

!> @author Paul Barrère
!! @brief Le module fzero rassemble différentes pour calculer un zéro d'une
!! fonction.
module fzero
  use param
  use func

  implicit none

  contains

    !> @brief  La subroutine dicho cherche un zéro par dichotomie.
    !! @param [in] A Paramètre lié au degré de polarisation linéaire
    !! @param [in] B Paramètre dépendant de A et de l'angle de polarisation chi
    !! @param [in] chi Angle de polarisation
    !! @param [in] x0 Valeur de la borne inférieure de l'intervalle 
    !! @param [in] x1 Valeur de la borne supérieure de l'intervalle
    !! @param [in] delta Ecart au zéro maximal toléré
    !! @param [out] zero Zéro de la fonction 
    subroutine dicho(A,B,chi,x0,x1,delta,zero)
      implicit none
      real(kind=kd), intent(in) :: delta, A, B, chi
      real(kind=kd), intent(in) :: x0, x1
      real(kind=kd), intent(out) :: zero
      real(kind=kd) :: x_max, x_min, x_mid, y_min, y_mid
      integer :: i=0

      x_min=x0; x_max=x1
      call f(A,B,chi,x_min,y_min)
      
      do while (abs(y_min)>delta) 
        call f(A,B,chi,x_min,y_min)
        i=i+1
        x_mid=0.5_kd*(x_max+x_min)
        call f(A,B,chi,x_mid,y_mid)
        if (y_mid*y_min>0._kd) then
          x_min=x_mid
        else
          x_max=x_mid
        end if
      end do
      
      zero=x_mid
    end subroutine dicho

    !> @brief  La subroutine sec cherche un zéro par méthode de la sécante.
    !! @param [in] A Paramètre lié au degré de polarisation linéaire
    !! @param [in] B Paramètre dépendant de A et de l'angle de polarisation chi
    !! @param [in] chi Angle de polarisation
    !! @param [in] x0 Valeur de la borne inférieure de l'intervalle
    !! @param [in] x1 Valeur de la borne supérieure de l'intervalle
    !! @param [in] delta Ecart au zéro maximal toléré
    !! @param [out] zero Zéro de la fonction
    subroutine sec(A,B,chi,x0,x1,delta,zero)
      implicit none
      real(kind=kd), intent(in) :: delta, A, B, chi
      real(kind=kd), intent(in) :: x0, x1
      real(kind=kd), intent(out) :: zero
      real(kind=kd) :: x, x_min, x_max, y_max, y_min
      integer :: i=0
      
      x_min=x0; x_max=x1
      call f(A,B,chi,x_max,y_max)
      
      do while (y_max>delta)
        call f(A,B,chi,x_max,y_max)
        call f(A,B,chi,x_min,y_min)
        i=i+1
        x=x_max-y_max*(x_max-x_min)/(y_max-y_min)
        x_min=x_max
        x_max=x
      end do
      
      zero=x_max
    end subroutine sec

    !> @brief  La subroutine dicho_sec cherche un zéro par la méthode
    !! de la fausse position.
    !! @param [in] A Paramètre lié au degré de polarisation linéaire
    !! @param [in] B Paramètre dépendant de A et de l'angle de polarisation chi
    !! @param [in] chi Angle de polarisation
    !! @param [in] x0 Valeur de la borne inférieure de l'intervalle
    !! @param [in] x1 Valeur de la borne supérieure de l'intervalle
    !! @param [in] delta Ecart au zéro maximal toléré
    !! @param [out] zero Zéro de la fonction
    subroutine dicho_sec(A,B,chi,x0,x1,delta,zero)
      implicit none
      real(kind=kd), intent(in) :: delta, A, B, chi
      real(kind=kd), intent(in) :: x0, x1
      real(kind=kd), intent(out) :: zero
      real(kind=kd) :: x_min, x_max, x_mid, y_min, y_max, y_mid
      integer :: i=0

      x_min=x0; x_max=x1
      do while (abs(x_max-x_min)>delta)
        call f(A,B,chi,x_max,y_max)
        call f(A,B,chi,x_min,y_min)
        i=i+1
        x_mid=x_max-y_max*(x_min-x_max)/(y_min-y_max)
        call f(A,B,chi,x_mid,y_mid)
        if (y_mid*y_min>0._kd) then
          x_min=x_mid
        else
          x_max=x_mid
        end if
      end do

      zero=x_mid
    end subroutine dicho_sec

    !> @brief  La subroutine newton cherche un zéro par la méthode de Newton
    !! quand on est pas trop loin du zéro.
    !! @param [in] A Paramètre lié au degré de polarisation linéaire
    !! @param [in] B Paramètre dépendant de A et de l'angle de polarisation chi
    !! @param [in] chi Angle de polarisation
    !! @param [in] x_min Valeur de la borne inférieure de l'intervalle
    !! @param [in] x_max Valeur de la borne supérieure de l'intervalle
    !! @param [in] delta Ecart au zéro maximal toléré
    !! @param [out] zero Zéro de la fonction
    subroutine newton(A,B,chi,x_min,x_max,delta,zero)
      implicit none
      real(kind=kd), intent(in) :: delta, A, B, chi
      real(kind=kd), intent(in) :: x_min, x_max
      real(kind=kd), intent(out) :: zero
      real(kind=kd) :: x, y, y_prime
      integer :: i=0

      ! On se rapproche du zéro par dichotomie pour assurer la convergence.
      call dicho(A,B,chi,x_min,x_max,10._kd**(-2),x)

      ! On applique la méthode de Newton.
      call f(A,B,chi,x,y)
      do while (abs(y)>delta)
        call f(A,B,chi,x,y)
        call f_prime(A,chi,x,y_prime)
        i=i+1
        x=x-y/y_prime
      end do

      zero=x
    end subroutine newton

end module fzero
