!> @file mod_random.f90
!! @author Paul Barrère
!! @brief Le fichier mod_random.f90 contient les subroutines pour gérer le seed
!! du générateur de nombre aléatoire.

!> @author Paul Barrère
!! @brief Le module random contient les subroutines pour gérer le seed
!! du générateur de nombre aléatoire.
module random
  use param
  implicit none

  integer :: seed_size !< Taille du tableau stockant le seed
  integer, allocatable :: tab_seed(:) !< Tableau stockant le seed

  contains

    !> @brief La subroutine seed_file sauvegarde le seed.
    !! @param [in] file_name Nom du fichier stockant le seed
    subroutine seed_file(file_name)
      implicit none
      character(len=1024), intent(in) :: file_name

      call random_seed

      ! Get the seed
      call random_seed(size=seed_size)
      allocate(tab_seed(seed_size))
      call random_seed(get=tab_seed(1:seed_size))

      ! Write the seed in a file
      open(11, file=trim(file_name), form='formatted', status='unknown')
      write(11,*) tab_seed
      close(11)

      deallocate(tab_seed)
    end subroutine seed_file

    !> @brief La subroutine put_seed permet d'utiliser un seed en particulier
    !! pour avoir la même suite de nombres.
    !! @param [in] Nom du fichier stockant le seed
    subroutine put_seed(seed_name)
      implicit none
      character(len=*), intent(in) :: seed_name

      call random_seed(size=seed_size)
      allocate(tab_seed(seed_size))
                                                                                                                              
      open(11, file=trim(seed_name), form='formatted', status='unknown')
      read(11,*) tab_seed(1:seed_size)
      close(11)

      call random_seed(put=tab_seed(1:seed_size))

      deallocate(tab_seed)
    end subroutine put_seed

end module random

