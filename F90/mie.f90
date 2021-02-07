Module Mie
Implicit None
Private
Public :: BHMIE
Contains

  !> @file mie.f90
  !! @author C.F. Bohren et D. Huffman
  !! @brief Le fichier mie.f90 contient 1 routine qui calcule les paramètres 
  !! de la diffusion de Mie d'une particule à travers une particule sphérique 
  !! de facteur de forme X.

  !> @author C.F. Bohren et D. Huffman
  !! @brief Le module Mie calcule les paramètres de la diffusion de Mie
  !! d'une particule à travers une particule sphérique de facteur de forme X.
  !! @param [in] X Facteur de forme
  !! @param [in] NANG Nombre d'angles entre 0 et pi/2
  !! @param [in] REFREL Rapport entre l'indice de réfraction 
  !! de la particule diffusante et du milieu environnant
  !! @param [out] S1,S2 Vecteurs pour construire la matrice de Mueller  
  !! @param [out] QEXT Facteur d'efficacité pour l'extinction
  !! @param [out] QSCA Facteur d'efficacité pour la diffusion
  !! @param [out] QBACK Facteur d'efficacité pour la rétrodiffusion
  !! @param [out] GSCA Valeur moyenne du cosinus de l'angle de diffusion
  !! @param wp Double précision
  Subroutine BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
    Implicit None
    Integer,  Parameter :: wp = Selected_real_kind(15) ! uziebfc
    Real(wp), Parameter :: one=1._wp, zero=0._wp, two=2._wp
    ! Declare parameters:
    Integer, Parameter::  MXNANG=1000,NMXX=150000
    ! Arguments:
    Integer  :: NANG
    Real(wp) :: GSCA,QBACK,QEXT,QSCA,X
    Complex(wp) :: REFREL
!!$    Complex(wp), Dimension(2*MXNANG-1) :: S1, S2
    Complex(wp), Dimension(:) :: S1, S2
    ! Local variables:
    Integer :: J,JJ,N,NSTOP,NMX,NN
    Real(wp) :: CHI,CHI0,CHI1,DANG,DX,EN,FN,P,PII,PSI,PSI0,PSI1,&
         THETA,XSTOP,YMOD
    Real(wp) ::  AMU(MXNANG),PI(MXNANG),PI0(MXNANG),PI1(MXNANG),&
         TAU(MXNANG)
    Complex(wp) ::  AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
    Complex(wp), Dimension(NMXX) :: D
    !***********************************************************************
    ! Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
    !    to calculate scattering and absorption by a homogenous isotropic
    !    sphere.
    ! Given:
    !    X = 2*pi*a/lambda
    !    REFREL = (complex refr. index of sphere)/(real index of medium)
    !    NANG = number of angles between 0 and 90 degrees
    !           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
    !           if called with NANG<2, will set NANG=2 and will compute
    !           scattering for theta=0,90,180.
    ! Returns:
    !    S1(1 : 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
    !                                scatt. E perp. to scatt. plane)
    !    S2(1 : 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
    !                                scatt. E parr. to scatt. plane)
    !    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
    !    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
    !    QBACK = (dC_sca/domega)/pi*a**2
    !          = backscattering efficiency [NB: this is (1/4*pi) smaller
    !            than the "radar backscattering efficiency"; see Bohren &
    !            Huffman 1983 pp. 120-123]
    !    GSCA = <cos(theta)> for scattering
    !
    ! Original program taken from Bohren and Huffman (1983), Appendix A
    ! Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
    ! in order to compute <cos(theta)>
    ! 91/05/07 (BTD): Modified to allow NANG=1
    ! 91/08/15 (BTD): Corrected error (failure to initialize P)
    ! 91/08/15 (BTD): Modified to enhance vectorizability.
    ! 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
    ! 91/08/15 (BTD): Changed definition of QBACK.
    ! 92/01/08 (BTD): Converted to full double precision and double complex
    !                 eliminated 2 unneed lines of code
    !                 eliminated redundant variables (e.g. APSI,APSI0)
    !                 renamed RN -> EN = double precision N
    !                 Note that DOUBLE COMPLEX and DCMPLX are not part
    !                 of f77 standard, so this version may not be fully
    !                 portable.  In event that portable version is
    !                 needed, use src/bhmie_f77.f
    ! 93/06/01 (BTD): Changed AMAX1 to generic function MAX
    !***********************************************************************
    !*** Safety checks
    If(NANG.Gt.MXNANG) Stop '***Error: NANG > MXNANG in bhmie'
    If(NANG.Lt.2)NANG=2
    !*** Obtain pi:
    PII = 4._wp*Atan(1._wp)
    DX = X
    DREFRL = REFREL
    Y = X*DREFRL
    YMOD = Abs(Y)
    !
    !*** Series expansion terminated after NSTOP terms
    !    Logarithmic derivatives calculated from NMX on down
    XSTOP = X + 4._wp*X**(one/3._wp) + two
    NMX = Max(XSTOP,YMOD) + 15
    ! BTD experiment 91/1/15: add one more term to series and compare results
    !      NMX=AMAX1(XSTOP,YMOD)+16
    ! test: compute 7001 wavelengths between .0001 and 1000 micron
    ! for a=1.0micron SiC grain.  When NMX increased by 1, only a single
    ! computed number changed (out of 4*7001) and it only changed by 1/8387
    ! conclusion: we are indeed retaining enough terms in series!
    NSTOP = XSTOP

    If(NMX.Gt.NMXX)Then
       Write(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
       Stop
    Endif
    !*** Require NANG.GE.1 in order to calculate scattering intensities
    DANG = zero
    If(NANG.Gt.1)DANG=.5_wp*PII/Real(NANG-1,wp)
    Do J=1,NANG
       THETA = Real(J-1,wp)*DANG
       AMU(J) = Cos(THETA)
    End Do

    PI0(1:NANG) = zero
    PI1(1:NANG) = one

    NN=2*NANG-1 
    S1(1:NN) = (zero,zero)
    S2(1:NN) = (zero,zero)

    !*** Logarithmic derivative D(J) calculated by downward recurrence
    !    beginning with initial value (0.,0.) at J=NMX

    D(NMX)=(zero,zero)
    NN=NMX-1
    Do N=1,NN
       EN = NMX-N+1
       D(NMX-N) = (EN/Y) - (one/(D(NMX-N+1)+EN/Y))
    End Do

    !*** Riccati-Bessel functions with real argument X
    !    calculated by upward recurrence

    PSI0 = Cos(DX)
    PSI1 = Sin(DX)
    CHI0 =-Sin(DX)
    CHI1 = Cos(DX)
    XI1 = Cmplx(PSI1,-CHI1,wp)
    QSCA = zero
    GSCA = zero
    P = -one
    Do N=1,NSTOP
       EN = N
       FN = (two*EN+one)/(EN*(EN+one))
       ! for given N, PSI  = psi_n        CHI  = chi_n
       !              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
       !              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
       ! Calculate psi_n and chi_n
       PSI = (two*EN-one)*PSI1/DX-PSI0
       CHI = (two*EN-one)*CHI1/DX-CHI0
       XI=DCMPLX(PSI,-CHI)

       !*** Store previous values of AN and BN for use
       !    in computation of g=<cos(theta)>
       If(N.Gt.1)Then
          AN1 = AN
          BN1 = BN
       End If

       !*** Compute AN and BN:
       AN = (D(N)/DREFRL+EN/DX)*PSI-PSI1
       AN = AN/((D(N)/DREFRL+EN/DX)*XI-XI1)
       BN = (DREFRL*D(N)+EN/DX)*PSI-PSI1
       BN = BN/((DREFRL*D(N)+EN/DX)*XI-XI1)

       !*** Augment sums for Qsca and g=<cos(theta)>
       QSCA = QSCA+(two*EN+one)*(Abs(AN)**2+Abs(BN)**2)
       GSCA = GSCA+((two*EN+one)/(EN*(EN+one)))* &
            (Real(AN,wp)*Real(BN,wp)+IMAG(AN)*IMAG(BN))
       If(N.Gt.1)Then
          GSCA=GSCA+((EN-1.)*(EN+1.)/EN)* &
               (Real(AN1,wp)*Real(AN,wp)+Aimag(AN1)*Aimag(AN)+ &
               Real(BN1,wp)*Real(BN,wp)+Aimag(BN1)*Aimag(BN))
       End If

       !*** Now calculate scattering intensity pattern
       !    First do angles from 0 to 90
       Do J=1,NANG
          JJ = 2*NANG-J
          PI(J) = PI1(J)
          TAU(J) = EN*AMU(J)*PI(J)-(EN+one)*PI0(J)
          S1(J) = S1(J) + FN*(AN*PI(J)+BN*TAU(J))
          S2(J) = S2(J) + FN*(AN*TAU(J)+BN*PI(J))
       End Do

       !*** Now do angles greater than 90 using PI and TAU from
       !    angles less than 90.
       !    P=1 for N=1,3,...; P=-1 for N=2,4,...
       P = -P
       Do J=1,NANG-1
          JJ = 2*NANG-J
          S1(JJ) = S1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
          S2(JJ) = S2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
       End Do
       PSI0 = PSI1
       PSI1 = PSI
       CHI0 = CHI1
       CHI1 = CHI
       XI1 = Cmplx(PSI1,-CHI1,wp)

       !*** Compute pi_n for next value of n
       !    For each angle J, compute pi_n+1
       !    from PI = pi_n , PI0 = pi_n-1
       Do J=1,NANG
          PI1(J) = ((two*EN+one)*AMU(J)*PI(J)-(EN+one)*PI0(J))/EN
          PI0(J) = PI(J)
       End Do
    End Do

    !*** Have summed sufficient terms.
    !    Now compute QSCA,QEXT,QBACK,and GSCA
    GSCA = two*GSCA/QSCA
    QSCA = (two/(DX*DX))*QSCA
    QEXT = (4._wp/(DX*DX))*Real(S1(1),wp)
    QBACK = (Abs(S1(2*NANG-1))/DX)**2/PII
  End Subroutine BHMIE
End Module Mie
