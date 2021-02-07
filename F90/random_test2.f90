!> @file random_test2.f90
!! @author Paul Barrère
!! @brief Le programme random_test2.f90 teste les fonctions liées au tirage 
!! aléatoire.
program random_test2
  use param

  implicit none
  integer :: k
  integer, allocatable :: tab_seed(:)
  real(kind=kd), dimension(10) :: tab

  !call random_seed

  ! Write the seed in a file
  call random_seed(size=k)
  allocate(tab_seed(k))

  open(11, file='seed.txt', form='formatted', status='unknown')
  read(11,*) tab_seed(1:k)
  close(11)
  
  call random_seed(put=tab_seed(1:k))

  call random_number(harvest=tab)  
  print *, tab
end program random_test2
