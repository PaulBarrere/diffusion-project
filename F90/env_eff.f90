!> @file env_eff.f90
!! @author Paul Barrère
!! @brief Le programme env_eff.f90 calcule  l'efficacité des enveloppes
!! de Rayleigh, de Henyey-Greenstein, et de Poly3 en fonction du facteur de
!! forme et en choisissant le paramètre maximisant l'efficacité.
!! @param x Facteur de forme
!! @param g Paramètre de l'enveloppe de Henyey-Greenstein
!! @param w Paramètre de l'enveloppe PolyN
!! @param refrel Rapport entre l'indice de réfraction 
!! de la particule diffusante et du milieu environnant
!! @param best_eff Efficacité de l'enveloppe la plus adaptée
!! @param best_g Paramètre maximisant l'efficacité de l'enveloppe 
!! Henyey-Greenstein
!! @param best_w Paramètre maximisant l'efficacité de l'enveloppe
!! PolyN

program env_eff
  use param
  use random
  use init_diff
  use init_mie_diff
  use sortie_diff
  use simu_diff
  use env_mod

  implicit none
  real(kind=kd) :: x 
  real(kind=kd) :: g
  real(kind=kd) :: w 
  real(kind=kd) :: best_eff, best_g, best_w
  
  complex(kind=kd) :: refrel
  integer :: i, j

  open(11, file='Data/env_eff_rayl.dat', status='unknown')
  open(12, file='Data/env_eff_hg.dat', status='unknown')
  open(13, file='Data/env_eff_polyn.dat', status='unknown')

  call calcul_refrel('blue','H2O',refrel)

  ! Meilleures enveloppes avec l'approximation de Rayleigh
  do i=0,468
    x1=1.025_kd**(i-187)
    call efficacite(x1,refrel,'fast','rayl',0._kd,0._kd,0._kd,eff)
    write(11,*) x1, eff
  end do
  print *, "Test pour une enveloppe de type Rayleigh fini"

  ! Meilleures enveloppes avec les fonctions HG
  do i=0,468
    x1=1.025_kd**(i-187)
    call efficacite(x1,refrel,'fast','hg',-0.99_kd,0._kd,0._kd,best_eff)
    best_g=-0.9_kd
    do j=2,19
      g=real(j,kind=kd)/10._kd-1._kd
      call efficacite(x1,refrel,'fast','hg',g,0._kd,0._kd,eff) 
      if (best_eff<eff) then
         best_eff=eff
         best_g=g
      end if
    end do
    write(12,*) x1, best_eff, best_g
  end do
  print *, "Test pour une enveloppe de type HG fini"

  ! Meilleures enveloppes avec les fonctions PolyN
  do i=0,468
    x1=1.025_kd**(i-187)
    call efficacite(x1,refrel,'fast','polyn',0._kd,3._kd,0._kd,best_eff)
    best_w=0._kd
    do j=1,10
      w=real(j,kind=kd)/10._kd
      call efficacite(x1,refrel,'fast','polyn',w,3._kd,0._kd,eff)
      if (best_eff<eff) then
         best_eff=eff
         best_w=w
      end if
    end do
    write(13,*) x1, best_eff, best_w
  end do
  print *, "Test pour une enveloppe de type Polyn fini"
end program env_eff
