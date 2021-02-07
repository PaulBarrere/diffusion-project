!> @file mod_proj_stereo.f90
!! @author Paul Barrère
!! @brief Le fichier mod_proj_stereo.f90 rassemble les modules pour effectuer  
!! une projection stéréographiques et créer une image en format pixmap (.ppm)

!> @author Paul Barrère
!! @brief Le module read_phase contient 2 subroutines permettant de lire les 
!! fichiers .dat en distinguant la diffusion de Rayleigh (applicable aussi aux 
!! diffusions isotropes, Henyey_Greenstein et polyN) et la diffusion de Mie qui
!! ne stockent pas les mêmes données.
module read_phase
  use param
  implicit none  

  contains

    !> @brief La subroutine read_phase_rayl lit les fichiers .dat de 3
    !! colonnes: indice, intensité, angle
    !! @param [in] N Taille du tableau phase
    !! @param [in] file_name Nom du fichier .dat
    !! @param [out] phase Tableau stockant la fonction de phase
    !! @param [out] theta Tableau stockant les angles de diffusion
    subroutine read_phase_rayl(N,file_name,phase,theta)
      implicit none
      integer, intent(in) :: N
      character(len=*), intent(in) :: file_name
      real(kind=kd), dimension(*), intent(out) :: phase, theta   
      integer :: num, i  

      open(61,file=trim(file_name))
      do i=1,N
        read(61,*) num, phase(i), theta(i)
      end do
      close(61)
    end subroutine read_phase_rayl

    !> @brief La subroutine read_phase_mie lit les fichiers .dat de 6
    !! colonnes: indice, 4 composantes du vecteur de Stokes, angle
    !! @param [in] N Taille du tableau phase
    !! @param [in] file_name Nom du fichier .dat
    !! @param [out] phase Tableau stockant la fonction de phase
    !! @param [out] theta Tableau stockant les angles de diffusion
    subroutine read_phase_mie(N,file_name,phase,theta)
      implicit none
      integer, intent(in) :: N
      character(len=*), intent(in) :: file_name
      real(kind=kd), dimension(*), intent(out) :: theta
      real(kind=kd), dimension(N,4), intent(out) :: phase 
      integer :: num, i

      open(21,file=trim(file_name))
      do i=1,N
        read(21,*) num,phase(i,1),phase(i,2),phase(i,3),phase(i,4),theta(i)
      end do
      close(21)
    end subroutine read_phase_mie

end module read_phase
!******************************************************************************!

!******************************************************************************!
!> @author Paul Barrère
!! @brief Le module interpolation contient une subroutine lissant une fonction
!! et une subroutine permettant d'évaluer la fonction en un point.
module interpolation
  use param
  implicit none

  contains
    
    !> @brief La subroutine lisse permet de lisser une fonction en polynômes
    !! cos(theta).
    !! @param [in] deg Degré du polynôme
    !! @param [inout] theta Tableau stockant les angles de diffusion
    !! @param [inout] phase Tableau stockant la fonction de phase
    !! @param [out] coeff Tableau stockant les coefficients du polynôme dans 
    !! l'ordre croissant des puissances
    subroutine lisse(deg,theta,phase,coeff)
      implicit none
      integer, intent(in) :: deg
      real(kind=kd), dimension(deg+1), intent(out) :: coeff
      real(kind=kd), dimension(N_phase), intent(inout) :: theta, phase
      
      real(kind=kd), allocatable :: A(:,:) 
      real(kind=kd), allocatable :: B(:) 
      real(kind=kd), allocatable :: work(:)
      integer :: info, lwork, i  
        
      ! Construction de la matrice A
      allocate(A(N_phase,deg+1))
      do i=0,deg
        A(:,i+1)=dcos(theta)**i
      end do
      
      ! Initialisation de la matrice avec les coeff à calculer
      allocate(B(N_phase))
      B=phase

      ! Autres arguments pour DGELS
      lwork=10*N_phase
      allocate(work(lwork))
  
      call dgels('N', N_phase, deg+1, 1, A, N_phase, B, N_phase, work, lwork, info)
      coeff=B(1:deg+1)
      open(35,file="Data/coefficients")
      write(35,*) coeff
      close(35)
      deallocate(A, B, work) 
    end subroutine lisse

    !> @brief La subroutine interpol évalue la fonction en un point à partir 
    !! du polynôme approximant la fonction de phase.
    !! @param [in] deg Degré du polynôme
    !! @param [in] coeff Tableau stockant les coefficients du polynôme dans 
    !! l'ordre croissant des puissances
    !! @param [in] Cosinus de l'angle de diffusion
    !! @param [out] fit_mu Valeur de la fonction au point mu
    subroutine interpol(deg,coeff,mu,fit_mu)
      implicit none
      integer, intent(in) :: deg
      real(kind=kd), dimension(deg+1), intent(in) :: coeff
      real(kind=kd), intent(in) :: mu
      real(kind=kd), intent(out) :: fit_mu
      integer :: k

      fit_mu=0._kd

      do k=0,deg
         fit_mu=fit_mu+coeff(k+1)*mu**k
      end do
    end subroutine interpol

end module interpolation
!******************************************************************************!

!******************************************************************************!
!> @author Paul Barrère
!! @brief Le module pixmap contient une subroutines créant une image ppm en noir
!!  et blanc, et une image ppm en couleur RGB.
module pixmap
  use param
  implicit none

  contains

    !> @brief La subroutine ppm construit une image ppm en noir et blanc.
    !! @param [in] num Numéro associé au fichier .ppm
    !! @param [in] ppm_name Nom du fichier .ppm
    !! @param [in] maxi Valeur maximale du pixel
    !! @param [in] Map Tableau stockant la valeur des pixels
    subroutine ppm(num,ppm_name,maxi,Map)
      implicit none
      character(len=*), intent(in) :: ppm_name
      integer, intent(in) :: num, maxi
      integer :: n_x, n_y
      integer, dimension(N_pix,N_pix), intent(in) :: Map

      n_x=size(Map,dim=1)
      n_y=size(Map,dim=2)

      open(unit=num, file=trim(ppm_name), position='Rewind')
      write(num,'(a)') "P2"
      write(num,*) N_pix, N_pix
      write(num,*) maxi
      write(num,*) Map
      close(num)
    end subroutine ppm

    !> @brief La subroutine ppm construit une image ppm en couleurs RGB.
    !! @param [in] ppm_name Nom du fichier .ppm
    !! @param [in] maxi Valeur maximale du pixel
    !! @param [in] MapR,MapG,MapB Tableaux stockant la valeur des pixels
    !! pour chaque couleur
    subroutine ppm_RGB(ppm_name,maxi,MapR,MapG,MapB)
      use param

      implicit none
      character(len=*), intent(in) :: ppm_name
      integer, intent(in) :: maxi
      integer :: n_x, n_y
      integer, dimension(N_pix,N_pix), intent(in) :: MapR, MapG, MapB
      integer, dimension(3,N_pix,N_pix) :: tabRGB
     
      n_x=size(tabRGB,dim=2)
      n_y=size(tabRGB,dim=3)

      tabRGB(1,:,:)=MapR
      tabRGB(2,:,:)=MapG
      tabRGB(3,:,:)=MapB

      open(unit=12, file=trim(ppm_name), position='Rewind')
      write(12,'(a)') "P3"
      write(12,*) n_x, n_y
      write(12,*) maxi
      write(12,*) tabRGB
      close(12)       
    end subroutine ppm_RGB

end module pixmap
!******************************************************************************!

!******************************************************************************!
!> @author Paul Barrère
!! @brief Le module build_stereo_map contient les subroutines pour faire
!! la projection stéréographique.
module build_stereo_map
  use param
  use pixmap
  use interpolation
  implicit none

  contains
    
    !> @brief La subroutine contraste_lin permet de gérer linéairement 
    !! le contraste de l'image.
    !! @param [inout] Map Tableau stockant la valeur des pixels
    !! @param [in] inf Borne inférieure en dessous de laquelle le pixel devient 0
    !! @param [in] sup Borne supérieure en dessous de laquelle le pixel devient 
    !! 255
    subroutine contraste_lin(Map,inf,sup)
      use param

      implicit none
      integer, dimension(N_pix,N_pix), intent(inout) :: Map
      integer, intent(in) :: inf, sup
      integer :: i, j

      do i=1,N_pix
        do j=1,N_pix
          if (Map(i,j)<inf) then
            Map(i,j)=0
          elseif (Map(i,j)>sup) then
            Map(i,j)=255
          else
            Map(i,j)=floor(real(255,kind=kd)*real(Map(i,j)-inf,kind=kd)/real(sup-inf,kind=kd))
          end if
        end do
      end do

    end subroutine contraste_lin

    !> @brief La subroutine stereo_map_liss crée deux projections 
    !! stéréographiques (diffusion et rétrodiffusion) via le lissage 
    !! de la fonction de phase.
    !! @param [in] deg Degré du polynôme
    !! @param [in] N Taille du tableau I
    !! @param [inout] I Tableau stockant la fonction de phase
    !! @param [inout] theta Tableau stckant les angles de diffusion
    !! @param [out] MapF Tableau stockant la valeur des pixels pour la diffusion
    !! @param [out] MapB Tableau stockant la valeur des pixels pour la 
    !! rétrodiffusion
    subroutine stereo_map_liss(deg,N,I,theta,MapF,MapB)
      implicit none
      integer, intent(in) :: deg, N
      real(kind=kd), dimension(N), intent(inout) :: I, theta
      real(kind=kd) :: dt, x, y, t, mu
      integer :: j, k, l
      real(kind=kd), allocatable, intent(out) :: MapF(:,:), MapB(:,:)
      real(kind=kd), allocatable :: tScale(:), coeff(:)

      ! Initialisation des images
      allocate(MapF(N_pix,N_pix)); allocate(MapB(N_pix,N_pix))

      ! Lissage de la phase
      allocate(coeff(deg+1))
      call lisse(deg,theta,I,coeff)

      ! Construction du tableau d'échelles
      allocate(tScale(N_pix))
      dt=2._kd/real(N_pix,kind=kd)
      do l=1,N_pix
        tScale(l)=-1._kd+real(l,kind=kd)*dt-dt/2._kd
      end do
      
      ! Construction des images
      do j=1, N_pix
        x=tScale(j)
        do k=1, N_pix
          y=tScale(k)
          t=x*x+y*y
          mu=(1._kd-t)/(1._kd+t)
          if (t<=1._kd) then
          call interpol(deg,coeff,mu,MapF(j,k))
          call interpol(deg,coeff,-mu,MapB(j,k))
          else
            MapF(j,k)=0._kd
            MapB(j,k)=0._kd
          end if
        end do
      end do
    end subroutine stereo_map_liss

    !> @brief La subroutine stereo_map_RGB crée les projections stéréographiques
    !! pour les couleurs RGB en image .ppm.
    !! @param [inout] I_R,I_G,I_B Tableau stockant la fonction de phase 
    !! des couleurs RGB
    !! @param [inout] theta Tableau stckant les angles de diffusion
    !! @param [in] file_name_diff Nom de l'image .ppm pour la diffusion
    !! @param [in] file_name_rdiff Nom de l'image .ppm pour la rétrodiffusion
    subroutine stereo_map_RGB(I_R,I_G,I_B,theta,file_name_diff,file_name_rdiff)
      use param

      implicit none
      character(len=*), intent(in) :: file_name_diff, file_name_rdiff
      real(kind=kd), dimension(N_theta), intent(inout) :: theta
      real(kind=kd), dimension(N_theta), intent(inout) :: I_R, I_G, I_B
      real(kind=kd), allocatable :: MapFR(:,:), MapFG(:,:), MapFB(:,:)
      real(kind=kd), allocatable :: MapBR(:,:), MapBG(:,:), MapBB(:,:)
      real(kind=kd) :: val_max, val_max_F, val_max_B

      ! Projection stéréographique
      call stereo_map_liss(21,N_theta,I_R,theta,MapFR,MapBR)
      call stereo_map_liss(21,N_theta,I_G,theta,MapFG,MapBG)
      call stereo_map_liss(21,N_theta,I_B,theta,MapFB,MapBB)

      ! Normalisation
      val_max=max(maxval(MapFR),maxval(MapFG),maxval(MapFB),maxval(MapBR),maxval(MapBG),maxval(MapBB))
      val_max_F=max(maxval(MapFR),maxval(MapFG),maxval(MapFB))
      val_max_B=max(maxval(MapBR),maxval(MapBG),maxval(MapBB))

      MapFR=floor(MapFR*255._kd/val_max_F); MapFG=floor(MapFG*255._kd/val_max_F); MapFB=floor(MapFB*255._kd/val_max_F)
      MapBR=floor(MapBR*255._kd/val_max_B); MapBG=floor(MapBG*255._kd/val_max_B); MapBB=floor(MapBB*255._kd/val_max_B)

      ! Création de l'image pixmap
      call ppm_RGB(trim(file_name_diff),255,int(MapFR),int(MapFG),int(MapFB))
      call ppm_RGB(trim(file_name_rdiff),255,int(MapBR),int(MapBG),int(MapBB))
    end subroutine stereo_map_RGB

end module build_stereo_map
!******************************************************************************!

!******************************************************************************!
!> @author Paul Barrère
!! @brief Le module quick_build_stereo_map permet de faire rapidement 
!! une projection stéréographique. Il est juste utiliser pour faire des tests
!! (comme pour tester bhmie).
module quick_build_stereo_map
  use param
  use pixmap
  use interpolation
  use build_stereo_map
  implicit none

  contains

    !> @brief La subroutine stereo_map_liss crée deux projections 
    !! stéréographiques (diffusion et rétrodiffusion) via une interpolation
    !! linéaire.
    !! @param [in] N Taille du tableau I
    !! @param [inout] I Tableau stockant la fonction de phase
    !! @param [out] MapF_norm Tableau stockant la valeur des pixels 
    !! pour la diffusion
    !! @param [out] MapB_norm Tableau stockant la valeur des pixels 
    !! pour la rétrodiffusion
    subroutine stereo_map(N,I,MapF_norm,MapB_norm)
      implicit none
      integer, intent(in) :: N
      real(kind=kd), dimension(N), intent(inout) :: I
      real(kind=kd) :: x, y, t, mu
      integer :: j, k, l
      real(kind=kd), allocatable :: MapF(:,:), MapB(:,:)
      integer, allocatable, intent(out) :: MapF_norm(:,:), MapB_norm(:,:)
      real(kind=kd), allocatable :: tScale(:)

      ! Initialisation des images
      allocate(MapF(N_pix,N_pix)); allocate(MapB(N_pix,N_pix))
      allocate(MapF_norm(N_pix,N_pix)); allocate(MapB_norm(N_pix,N_pix))

      ! Construction du tableau d'échelles
      allocate(tScale(N_pix))
      tScale=-1._kd+2._kd*(/(l,l=1,N_pix)/)/real(N_pix, kind=kd)-1._kd/real(N_pix, kind=kd)

      do j=1, N_pix
        x=tScale(j)
        do k=1, N_pix
          y=tScale(k)
          t=x*x+y*y
          mu=(1._kd-t)/(1._kd+t)
          if (t<=1._kd) then
            MapF(j,k)=abs((I(floor(dacos(mu)*real(N-1,kind=kd)/pi)+1)+I(floor(dacos(mu)*real(N-1,kind=kd)/pi)+2))/2._kd)
            MapB(j,k)=abs((I(floor(dacos(-mu)*real(N-1,kind=kd)/pi)+1)+I(floor(dacos(-mu)*real(N-1,kind=kd)/pi)+2))/2._kd)
          else
            MapF(j,k)=0._kd
            MapB(j,k)=0._kd
          end if
        end do
      end do

      ! Normalisation
      MapF_norm=floor(MapF*255._kd/maxval(MapF))
      MapB_norm=floor(MapB*255._kd/maxval(MapB))
    end subroutine stereo_map

    !> @brief La subroutine stereo_map_RGB crée les projections stéréographiques
    !! pour les couleurs RGB en image .ppm juste pour faire le test de bhmie 
    !! consistant à mettre en avant le phénomène d'arc-en-ciel.
    !! @param [in] nang Nombre d'angles entre 
    !! @param [inout] I_R,I_G,I_B Tableau stockant la fonction de phase 
    !! des couleurs RGB
    !! @param [in] file_name_diff Nom de l'image .ppm pour la diffusion
    !! @param [in] file_name_rdiff Nom de l'image .ppm pour la rétrodiffusion
    subroutine stereo_map_RGB_rb(nang,I_R,I_G,I_B,file_name_diff,file_name_rdiff)
      implicit none
      character(len=*), intent(in) :: file_name_diff, file_name_rdiff
      integer, intent(in) :: nang
      real(kind=kd), dimension(2*nang-1), intent(inout) :: I_R, I_G, I_B
      integer, allocatable :: MapFR(:,:), MapFG(:,:), MapFB(:,:)
      integer, allocatable :: MapBR(:,:), MapBG(:,:), MapBB(:,:)
      
      ! Projection stéréographique
      call stereo_map(2*nang-1,I_R,MapFR,MapBR)
      call stereo_map(2*nang-1,I_G,MapFG,MapBG)
      call stereo_map(2*nang-1,I_B,MapFB,MapBB)

      ! Augmentation du contraste
      print *, "Attention: le contraste a été augmenté"
      call contraste_lin(MapBR,0,130)
      call contraste_lin(MapBG,0,130)
      call contraste_lin(MapBB,0,130)

      ! Création de l'image pixmap
      call ppm_RGB(trim(file_name_diff),255,MapFR,MapFG,MapFB)
      call ppm_RGB(trim(file_name_rdiff),255,MapBR,MapBG,MapBB)
    end subroutine stereo_map_RGB_rb

end module quick_build_stereo_map
