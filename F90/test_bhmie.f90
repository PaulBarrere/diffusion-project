!> @file test_bhmie.f90
!! @author Paul Barrère
!! @brief Le programme test_bhmie.f90 effectue quelques tests de la 
!! subroutine bhmie.
!! @param test Test que l'on veut effectuer: 'x_inf','f_phase','lambda_dpd','a_dpd'
program test_bhmie
  use param
  use sortie_diff
  use mie

  implicit none
  character(len=32), parameter :: test='f_phase'

  real(kind=kd), parameter :: a=250 !En micromètres
  real(kind=kd) :: a
  real(kind=kd), parameter :: lambda = 0.5 !En micromètres
  real(kind=kd) :: lambda

  real(kind=kd) :: x, d_mu
  complex(kind=kd) :: refrel
  integer, parameter :: nang=999
  complex(kind=kd), dimension(2*nang-1) :: S1, S2
  real(kind=kd) :: Qext, Qsca, Qback, Gsca

  real(kind=kd), dimension(2*nang-1) :: I_s
  real(kind=kd) :: C, k

  integer :: i

  !!! Facteur de forme infini
  if (test .eq. 'x_inf') then
    x=125000._kd
    print *, Qext

  !!! Fonction de phase
  else if (test .eq. 'f_phase') then
    x=0.01_kd
    refrel=cmplx(neau_B,keau_B,kind=kd)/cmplx(nair_B,0._kd,kind=kd)

    call bhmie(x, refrel, nang, S1, S2, Qext, Qsca, Qback, Gsca)

    I_s=0.5_kd*real((S2*conjg(S2)+S1*conjg(S1)),kind=kd)
    C=2._kd/(x*x*Qsca)
    I_s=C*I_s

    open(unit=15, file="Data/bhmie_diff.dat", status="unknown")
    d_mu=2._kd*(1._kd+eps)/real(N_phase,kind=kd)
    call save_res_bhmie(15, I_s)
    close(15)

  !!! Dépendance en lambda dans la limite x<=10
  else if (test .eq. 'lambda_dpd') then
    do i=1,15
      lambda=1.25_kd*2._kd**i
      x=2._kd*pi*a/lambda

      call bhmie(x, refrel, nang, S1, S2, Qext, Qsca, Qback, Gsca)
      k=lambda/(2._kd*pi)
      I_s=k*k*real(0.5_kd*(S2*conjg(S2)+S1*conjg(S1)),kind=kd)

      open(12, file='Data/Ivlambda_rayl_test_bhmie.dat', status='unknown')
      write(12,*) I_s
    end do
    close(12)

  !!! Dépendance en a dans la limite x<=10
  else if (test .eq. 'a_dpd') then
    do i=1,20
      a=7.95_kd*10._kd**(-6)*2**i
      x=2._kd*pi*a/lambda

      call bhmie(x, refrel, nang, S1, S2, Qext, Qsca, Qback, Gsca)
      k=lambda/(2._kd*pi)
      I_s=k*k*real(0.5_kd*(S2*conjg(S2)+S1*conjg(S1)),kind=kd)

      open(12, file='Data/Iva_rayl_test_bhmie.dat', status='unknown')
      write(12,*) I_s
    end do
    close(12)
  end if
end program test_bhmie 
