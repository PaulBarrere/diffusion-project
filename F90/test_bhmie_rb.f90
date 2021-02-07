!> @file test_bhmie_rb.f90
!! @author Paul Barrère
!! @brief Le programme test_bhmie_rb.f90 teste la subroutine bhmie en mettant
!! évidence le phénomène d'arc-en-ciel.
program test_bhmie_rb
  use param
  use quick_build_stereo_map 
  use mie

  implicit none
  real(kind=kd), parameter :: a=30 !En micromètres
  
  real(kind=kd) :: x_R, x_G, x_B
  complex(kind=kd) :: refrel_R, refrel_G, refrel_B
  integer, parameter :: nang=999
  complex(kind=kd), dimension(2*nang-1) :: S1_R, S1_G, S1_B
  complex(kind=kd), dimension(2*nang-1) :: S2_R, S2_G, S2_B

  real(kind=kd) :: Qext_R, Qsca_R, Qback_R, Gsca_R
  real(kind=kd) :: Qext_G, Qsca_G, Qback_G, Gsca_G
  real(kind=kd) :: Qext_B, Qsca_B, Qback_B, Gsca_B  

  real(kind=kd), dimension(2*nang-1) :: I_sR, I_sG, I_sB

  
  !!! Arc-en-ciel
  refrel_R=cmplx(neau_R,keau_R,kind=kd)/cmplx(nair_R,0._kd,kind=kd)
  refrel_G=cmplx(neau_G,keau_G,kind=kd)/cmplx(nair_G,0._kd,kind=kd)
  refrel_B=cmplx(neau_B,keau_B,kind=kd)/cmplx(nair_B,0._kd,kind=kd)

  x_R=2._kd*pi*a/lambda_R
  x_G=2._kd*pi*a/lambda_G
  x_B=2._kd*pi*a/lambda_B

  call bhmie(x_R, refrel_R, nang, S1_R, S2_R, Qext_R, Qsca_R, Qback_R, Gsca_R)
  call bhmie(x_G, refrel_G, nang, S1_G, S2_G, Qext_G, Qsca_G, Qback_G, Gsca_G)
  call bhmie(x_B, refrel_B, nang, S1_B, S2_B, Qext_B, Qsca_B, Qback_B, Gsca_B)

  I_sR=real(0.5_kd*(S2_R*conjg(S2_R)+S1_R*conjg(S1_R)),kind=kd)
  I_sG=real(0.5_kd*(S2_G*conjg(S2_G)+S1_G*conjg(S1_G)),kind=kd)
  I_sB=real(0.5_kd*(S2_B*conjg(S2_B)+S1_B*conjg(S1_B)),kind=kd)

  I_sR=2._kd*I_sR/(x_R*x_R*Qsca_R)
  I_sG=2._kd*I_sG/(x_G*x_G*Qsca_G)
  I_sB=2._kd*I_sB/(x_B*x_B*Qsca_B)

  open(12, file='Data/Ivtheta_test_bhmie_rb.dat', status='unknown')
  write(12,*) I_sR
  write(12,*) I_sG
  write(12,*) I_sB
  close(12)

  call stereo_map_RGB_rb(nang,I_sR,I_sG,I_sB,"Fig/stereo_wb_bhmie_diff.ppm","Fig/stereo_wb_bhmie_rdiff_contraste.ppm")

end program test_bhmie_rb
