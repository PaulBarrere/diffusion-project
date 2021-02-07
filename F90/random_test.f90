!> @file random_test.f90
!! @author Paul Barrère
!! @brief Le programme random_test.f90 teste les fonctions liées au tirage 
!! aléatoire.
program random_test
  use param

  implicit none
  integer :: k
  integer, allocatable :: tab_seed(:)
  real(kind=kd), dimension(10) :: tab

  call random_seed

  ! Get the seed
  call random_seed(size=k)
  allocate(tab_seed(k))
  call random_seed(get=tab_seed(1:k))

  ! Write the seed in a file
  open(11, file='seed.txt', form='formatted', status='unknown')
  write(11,*) tab_seed
  close(11)

  call random_number(harvest=tab)  
  print *, tab
end program random_test
