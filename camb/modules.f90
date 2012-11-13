! Modules used by cmbmain and other routines.

!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info) and Anthony Challinor
!     See readme.html for documentation.
!
!     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
!     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
!     Original CMBFAST copyright and disclaimer:
!
!     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
!     the Massachusetts Institute of Technology.  All rights reserved.
!
!     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
!     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
!     By way of example, but not limitation,
!     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
!     MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
!     THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
!     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
!
!     portions of this software are based on the COSMICS package of
!     E. Bertschinger.  See the LICENSE file of the COSMICS distribution
!     for restrictions on the modification and distribution of this software.


        module ModelParams
        use precision
        use Ranges
        use InitialPower
        use Reionization
        use Recombination
        use Errors

        implicit none
        public

        character(LEN=*), parameter :: version = 'Oct_12'

        integer :: FeedbackLevel = 0 !if >0 print out useful information about the model

        logical, parameter :: DebugMsgs=.false. !Set to true to view progress and timing

        logical, parameter :: DebugEvolution = .false. !Set to true to do all the evolution for all k

        logical ::  do_bispectrum  = .false.
        logical, parameter :: hard_bispectrum = .false. ! e.g. warm inflation where delicate cancellations

        logical, parameter :: full_bessel_integration = .false. !(go into the tails when calculating the sources)

        integer, parameter :: Nu_int = 0, Nu_trunc=1, Nu_approx = 2, Nu_best = 3
         !For CAMBparams%MassiveNuMethod
         !Nu_int: always integrate distribution function
         !Nu_trunc: switch to expansion in velocity once non-relativistic
         !Nu_approx: approximate scheme - good for CMB, but not formally correct and no good for matter power
         !Nu_best: automatically use mixture which is fastest and most accurate

        integer, parameter :: max_Nu = 5 !Maximum number of neutrino species
        integer, parameter :: max_transfer_redshifts = 150
        integer, parameter :: fileio_unit = 13 !Any number not used elsewhere will do
        integer, parameter :: outCOBE=0, outNone=1

        integer :: max_bessels_l_index  = 1000000
        real(dl) :: max_bessels_etak = 1000000*2


        real(dl), parameter ::  OutputDenominator =twopi
       !When using outNone the output is l(l+1)Cl/OutputDenominator


        Type(Regions) :: TimeSteps


        type TransferParams
            logical     ::  high_precision
            integer     ::  num_redshifts
            real(dl)    ::  kmax         !these are acutally q values, but same as k for flat
            integer     ::  k_per_logint ! ..
            real(dl)    ::  redshifts(max_transfer_redshifts)
        end type TransferParams

!other variables, options, derived variables, etc.

         integer, parameter :: NonLinear_none=0, NonLinear_Pk =1, NonLinear_Lens=2

! Main parameters type
        type CAMBparams

         logical   :: WantCls, WantTransfer
         logical   :: WantScalars, WantTensors, WantVectors
         logical   :: DoLensing
         logical   :: want_zstar, want_zdrag     !!JH for updated BAO likelihood.
         integer   :: NonLinear

         integer   :: Max_l, Max_l_tensor
         real(dl)  :: Max_eta_k, Max_eta_k_tensor
          ! _tensor settings only used in initialization,
          !Max_l and Max_eta_k are set to the tensor variables if only tensors requested

         real(dl)  :: omegab, omegac, omegav, omegan
         !Omega baryon, CDM, Lambda and massive neutrino
         real(dl)  :: H0,TCMB,yhe,Num_Nu_massless
         integer   :: Num_Nu_massive

         logical :: Nu_mass_splittings
         integer   :: Nu_mass_eigenstates  !1 for degenerate masses
         real(dl)  :: Nu_mass_degeneracies(max_nu)
         real(dl)  :: Nu_mass_fractions(max_nu)
             !The ratios of the masses

         integer   :: Scalar_initial_condition
         !must be one of the initial_xxx values defined in GaugeInterface

         integer   :: OutputNormalization
         !outNone, outCOBE, or C_OutputNormalization=1 if > 1

         logical   :: AccuratePolarization
           !Do you care about the accuracy of the polarization Cls?

         logical   :: AccurateBB
           !Do you care about BB accuracy (e.g. in lensing)

!Reionization settings - used if Reion%Reionization=.true.
         logical   :: AccurateReionization
           !Do you care about pecent level accuracy on EE signal from reionization?

         integer   :: MassiveNuMethod

         type(InitialPowerParams) :: InitPower  !see power_tilt.f90 - you can change this
         type(ReionizationParams) :: Reion
         type(RecombinationParams):: Recomb
         type(TransferParams)     :: Transfer

         real(dl) ::  InitialConditionVector(1:10) !Allow up to 10 for future extensions
          !ignored unless Scalar_initial_condition == initial_vector

         logical OnlyTransfers !Don't use initial power spectrum data, instead get Delta_q_l array
          !If true, sigma_8 is not calculated either

         logical DerivedParameters !calculate various derived parameters  (ThermoDerivedParams)

!Derived parameters, not set initially
         type(ReionizationHistory) :: ReionHist

         logical flat,closed,open
         real(dl) omegak
         real(dl) curv,r, Ksign !CP%r = 1/sqrt(|CP%curv|), CP%Ksign = 1,0 or -1
         real(dl) tau0,chi0 !time today and rofChi(CP%tau0/CP%r)

         end type CAMBparams

        type(CAMBparams) CP  !Global collection of parameters


       real(dl) scale !relative to CP%flat. e.g. for scaling lSamp%l sampling.

       logical ::call_again = .false.
          !if being called again with same parameters to get different thing


!     grhom =kappa*a^2*rho_m0
!     grhornomass=grhor*number of massless neutrino species
!     taurst,taurend - time at start/end of recombination
!     dtaurec - dtau during recombination
!     adotrad - a(tau) in radiation era

        real(dl) grhom,grhog,grhor,grhob,grhoc,grhov,grhornomass,grhok
        real(dl) taurst,dtaurec,taurend, tau_maxvis,adotrad

!Neutrinos
        real(dl) grhormass(max_nu)

!     nu_masses=m_nu*c**2/(k_B*T_nu0)
       real(dl) :: nu_masses(max_nu)

       real(dl) akthom !sigma_T * (number density of protons now)
       real(dl) fHe !n_He_tot / n_H_tot
       real(dl) Nnow


      integer :: ThreadNum = 0
       !If zero assigned automatically, obviously only used if parallelised

!Parameters for checking/changing overall accuracy
!If HighAccuracyDefault=.false., the other parameters equal to 1 corresponds to ~0.3% scalar C_l accuracy
!If HighAccuracyDefault=.true., the other parameters equal to 1 corresponds to ~0.1% scalar C_l accuracy (at L>600)
      logical :: HighAccuracyDefault = .false.

      real(dl) :: lSampleBoost=1._dl
          !Increase lSampleBoost to increase sampling in lSamp%l for Cl interpolation

      real(dl) :: AccuracyBoost =1._dl

          !Decrease step sizes, etc. by this parameter. Useful for checking accuracy.
          !Can also be used to improve speed significantly if less accuracy is required.
          !or improving accuracy for extreme models.
          !Note this does not increase lSamp%l sampling or massive neutrino q-sampling

      real(sp) :: lAccuracyBoost=1.
          !Boost number of multipoles integrated in Boltzman heirarchy

      integer, parameter :: lmin = 2
          !must be either 1 or 2

      real(dl), parameter :: OmegaKFlat = 5e-7_dl !Value at which to use flat code

      real(dl),parameter :: tol=1.0d-4 !Base tolerance for integrations

!     used as parameter for spline - tells it to use 'natural' end values
      real(dl), parameter :: spl_large=1.e40_dl

      integer, parameter:: l0max=4000

!     lmax is max possible number of l's evaluated
      integer, parameter :: lmax_arr = l0max

      character(LEN=1024) :: highL_unlensed_cl_template = 'HighLExtrapTemplate_lenspotentialCls.dat'
           !fiducial high-accuracy high-L C_L used for making small cosmology-independent numerical corrections
           !to lensing and C_L interpolation. Ideally close to models of interest, but dependence is weak.
      logical :: use_spline_template = .true.
      integer, parameter :: lmax_extrap_highl = 6000
      real(dl), allocatable :: highL_CL_template(:,:)

      integer, parameter :: derived_zstar=1, derived_rstar=2, derived_thetastar=3,derived_zdrag=4, &
         derived_rdrag=5,derived_kD=6,derived_thetaD=7 , derived_zEQ =8, derived_thetaEQ=9 !, derived_mnu=10
      integer, parameter :: nthermo_derived = 9
      real(dl) ThermoDerivedParams(nthermo_derived)

       contains


         subroutine CAMBParams_Set(P, error, DoReion)
           use constants
           type(CAMBparams), intent(in) :: P
           real(dl) GetOmegak, fractional_number, conv
           integer, optional :: error !Zero if OK
           logical, optional :: DoReion
           logical WantReion
           integer nu_i,actual_massless
           external GetOmegak
           real(dl), save :: last_tau0
           !Constants in SI units

            global_error_flag = 0

            if ((P%WantTensors .or. P%WantVectors).and. P%WantTransfer .and. .not. P%WantScalars) then
              call GlobalError( 'Cannot generate tensor C_l and transfer without scalar C_l',error_unsupported_params)
            end if

            if (present(error)) error = global_error_flag
            if (global_error_flag/=0) return

           if (present(DoReion)) then
            WantReion = DoReion
           else
            WantReion = .true.
           end if

           CP=P
           if (call_again) CP%DerivedParameters = .false.

           CP%Max_eta_k = max(CP%Max_eta_k,CP%Max_eta_k_tensor)

           if (CP%WantTransfer) then
              CP%WantScalars=.true.
              if (.not. CP%WantCls) then
                 CP%AccuratePolarization = .false.
                 CP%Reion%Reionization = .false.
              end if
           else
              CP%transfer%num_redshifts=0
           end if

           if (CP%Omegan == 0 .and. CP%Num_Nu_Massive /=0) then
              CP%Num_Nu_Massless = CP%Num_Nu_Massless + CP%Num_Nu_Massive
              CP%Num_Nu_Massive  = 0
           end if

           if (CP%Num_nu_massive > 0) then
               if (.not. CP%Nu_mass_splittings) then
                 !Default totally degenerate masses
                 CP%Nu_mass_eigenstates = 1
                 CP%Nu_mass_degeneracies(1) = CP%Num_Nu_Massive
                 CP%Nu_mass_fractions(1) = 1
               else
                 if (CP%Nu_mass_degeneracies(1)==0) CP%Nu_mass_degeneracies(1) = CP%Num_Nu_Massive
                 if (abs(sum(CP%Nu_mass_fractions(1:CP%Nu_mass_eigenstates))-1) > 1e-4) &
                   stop 'Nu_mass_fractions do not add up to 1'

                 if (abs(sum(CP%Nu_mass_degeneracies(1:CP%Nu_mass_eigenstates))-CP%Num_nu_massive) >1e-4 ) &
                    stop 'nu_mass_eigenstates do not add up to num_nu_massive'
                 if (CP%Nu_mass_eigenstates==0) stop 'Have Num_nu_massive>0 but no nu_mass_eigenstates'

               end if
           else
            CP%Nu_mass_eigenstates = 0
           end if

           if ((CP%WantTransfer).and. CP%MassiveNuMethod==Nu_approx) then
              CP%MassiveNuMethod = Nu_trunc
           end if

           CP%omegak = GetOmegak()

           CP%flat = (abs(CP%omegak) <= OmegaKFlat)
           CP%closed = CP%omegak < -OmegaKFlat

           CP%open = .not.CP%flat.and..not.CP%closed
           if (CP%flat) then
              CP%curv=0
              CP%Ksign=0
              CP%r=1._dl !so we can use tau/CP%r, etc, where CP%r's cancel
           else
           CP%curv=-CP%omegak/((c/1000)/CP%h0)**2
           CP%Ksign =sign(1._dl,CP%curv)
           CP%r=1._dl/sqrt(abs(CP%curv))
           end if
!  grho gives the contribution to the expansion rate from: (g) photons,
!  (r) one flavor of relativistic neutrino (2 degrees of freedom),
!  (m) nonrelativistic matter (for Omega=1).  grho is actually
!  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
!  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
!  (Used only to set the initial conformal time.)

           !H0 is in km/s/Mpc

           grhom = 3*CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)

          !grhom=3.3379d-11*h0*h0
           grhog = kappa/c**2*4*sigma_boltz/c**3*CP%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
          ! grhog=1.4952d-13*tcmb**4
           grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
          !grhor=3.3957d-14*tcmb**4
           !correction for fractional number of neutrinos, e.g. 3.04 to give slightly higher T_nu hence rhor
           !Num_Nu_massive is already integer, Num_Nu_massless can contain fraction
           !We assume all eigenstates affected the same way
           fractional_number  = CP%Num_Nu_massless + CP%Num_Nu_massive
           actual_massless = int(CP%Num_Nu_massless + 1e-6)
           if (actual_massless + CP%Num_Nu_massive /= 0) then
             grhor = grhor * fractional_number/(actual_massless + CP%Num_Nu_massive)
             grhornomass=grhor*actual_massless
           else
             !Prevent problems with n_eff < 1; thanks Zhen Hou
             grhornomass=grhor*CP%Num_Nu_massless
           end if
           grhormass=0
           do nu_i = 1, CP%Nu_mass_eigenstates
            grhormass(nu_i)=grhor*CP%Nu_mass_degeneracies(nu_i)
           end do
           grhoc=grhom*CP%omegac
           grhob=grhom*CP%omegab
           grhov=grhom*CP%omegav
           grhok=grhom*CP%omegak
!  adotrad gives the relation a(tau) in the radiation era:
           adotrad = sqrt((grhog+grhornomass+sum(grhormass(1:CP%Nu_mass_eigenstates)))/3)


           Nnow = CP%omegab*(1-CP%yhe)*grhom*c**2/kappa/m_H/Mpc**2

           akthom = sigma_thomson*Nnow*Mpc
              !sigma_T * (number density of protons now)

           fHe = CP%YHe/(mass_ratio_He_H*(1.d0-CP%YHe))  !n_He_tot / n_H_tot

           if (CP%omegan==0) then
              CP%Num_nu_massless = CP%Num_nu_massless + CP%Num_nu_massive
              CP%Num_nu_massive = 0
           end if

           if (.not.call_again) then

            call init_massive_nu(CP%omegan /=0)
            call init_background
            if (global_error_flag==0) then
             CP%tau0=TimeOfz(0._dl)
         ! print *, 'chi = ',  (CP%tau0 - TimeOfz(0.15_dl)) * CP%h0/100
             last_tau0=CP%tau0
             if (WantReion) call Reionization_Init(CP%Reion,CP%ReionHist, CP%YHe, akthom, CP%tau0, FeedbackLevel)
            end if
           else
              CP%tau0=last_tau0
           end if

           if ( CP%NonLinear==NonLinear_Lens) then
             CP%Transfer%kmax = max(CP%Transfer%kmax, CP%Max_eta_k/CP%tau0)
             if (FeedbackLevel > 0 .and. CP%Transfer%kmax== CP%Max_eta_k/CP%tau0) &
                  write (*,*) 'max_eta_k changed to ', CP%Max_eta_k
           end if

           if (CP%closed .and. CP%tau0/CP%r >3.14) then
             call GlobalError('chi >= pi in closed model not supported',error_unsupported_params)
           end if

           if (global_error_flag/=0) then
             if (present(error)) error = global_error_flag
             return
           end if

           if (present(error)) then
              error = 0
           else if (FeedbackLevel > 0 .and. .not. call_again) then
              write(*,'("Om_b h^2             = ",f9.6)') CP%omegab*(CP%H0/100)**2
              write(*,'("Om_c h^2             = ",f9.6)') CP%omegac*(CP%H0/100)**2
              write(*,'("Om_nu h^2            = ",f9.6)') CP%omegan*(CP%H0/100)**2
              write(*,'("Om_Lambda            = ",f9.6)') CP%omegav
              write(*,'("Om_K                 = ",f9.6)') CP%omegak
              write(*,'("Om_m (1-Om_K-Om_L)   = ",f9.6)') 1-CP%omegak-CP%omegav
              write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*CosmomcTheta()
              if (CP%Num_Nu_Massive > 0) then
                conv = k_B*(8*grhor/grhog/7)**0.25*CP%tcmb/eV !approx 1.68e-4
                do nu_i=1, CP%Nu_mass_eigenstates
                 write(*,'(f5.2, " nu, m_nu*c^2/k_B/T_nu0   = ",f8.2," (m_nu = ",f6.3," eV)")') &
                     CP%nu_mass_degeneracies(nu_i), nu_masses(nu_i),conv*nu_masses(nu_i)
                end do
              end if
           end if
           CP%chi0=rofChi(CP%tau0/CP%r)
           scale= CP%chi0*CP%r/CP%tau0  !e.g. changel sampling depending on approx peak spacing

         end subroutine CAMBParams_Set


         function GetTestTime()
           real(sp) GetTestTime
           real(sp) atime

!           GetTestTime = etime(tarray)
         !Can replace this if etime gives problems
         !Or just comment out - only used if DebugMsgs = .true.
           call cpu_time(atime)
           GetTestTime = atime

         end function GetTestTime


         function rofChi(Chi) !sinh(chi) for open, sin(chi) for closed.
         real(dl) Chi,rofChi

         if (CP%closed) then
            rofChi=sin(chi)
         else if (CP%open) then
            rofChi=sinh(chi)
         else
            rofChi=chi
         endif
         end function rofChi


         function cosfunc (Chi)
         real(dl) Chi,cosfunc

         if (CP%closed) then
            cosfunc= cos(chi)
         else if (CP%open) then
            cosfunc=cosh(chi)
         else
            cosfunc = 1._dl
         endif
         end function cosfunc

         function tanfunc(Chi)
         real(dl) Chi,tanfunc
         if (CP%closed) then
            tanfunc=tan(Chi)
         else if (CP%open) then
            tanfunc=tanh(Chi)
         else
            tanfunc=Chi
         end if

         end  function tanfunc

         function invsinfunc(x)
         real(dl) invsinfunc,x

         if (CP%closed) then
          invsinfunc=asin(x)
          else if (CP%open) then
          invsinfunc=log((x+sqrt(1._dl+x**2)))
          else
          invsinfunc = x
         endif
         end function invsinfunc

        function f_K(x)
         real(dl) :: f_K
         real(dl), intent(in) :: x
         f_K = CP%r*rofChi(x/CP%r)

        end function f_K


        function DeltaTime(a1,a2, in_tol)
        implicit none
        real(dl) DeltaTime, atol
        real(dl), intent(IN) :: a1,a2
        real(dl), optional, intent(in) :: in_tol
        real(dl) dtauda, rombint !diff of tau w.CP%r.t a and integration
        external dtauda, rombint

        if (present(in_tol)) then
         atol = in_tol
        else
         atol = tol/1000/exp(AccuracyBoost-1)
        end if
        DeltaTime=rombint(dtauda,a1,a2,atol)

        end function DeltaTime

        function TimeOfz(z)
        implicit none
        real(dl) TimeOfz
        real(dl), intent(IN) :: z

        TimeOfz=DeltaTime(0._dl,1._dl/(z+1._dl))
        end function TimeOfz

        function AngularDiameterDistance(z)
          real(dl) AngularDiameterDistance
          real(dl), intent(in) :: z

          AngularDiameterDistance = CP%r/(1+z)*rofchi(DeltaTime(1/(1+z),1._dl)/CP%r)

        end function AngularDiameterDistance

       function dsound_da(a)
          implicit none
          real(dl) dsound_da,dtauda,a,R,cs
          external dtauda

           R=3.0d4*a*CP%omegab*(CP%h0/100.0d0)**2
!          R = 3*grhob*a / (4*grhog) //above is mostly within 0.2% and used for previous consistency
           cs=1.0d0/sqrt(3*(1+R))
           dsound_da=dtauda(a)*cs

       end function dsound_da


       function CosmomcTheta()
         real(dl) zstar, astar, atol, rs, DA
         real(dl) CosmomcTheta
         real(dl) ombh2, omdmh2
         real(dl) rombint
         external rombint

         ombh2 = CP%omegab*(CP%h0/100.0d0)**2
         omdmh2 = (CP%omegac+CP%omegan)*(CP%h0/100.0d0)**2

    !!From Hu & Sugiyama
           zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
            (0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
               (omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))

           astar = 1/(1+zstar)
           atol = 1e-6
           rs = rombint(dsound_da,1d-8,astar,atol)
           DA = AngularDiameterDistance(zstar)/astar
           CosmomcTheta = rs/DA
    !       print *,'z* = ',zstar, 'r_s = ',rs, 'DA = ',DA, rs/DA

       end function CosmomcTheta

   end module ModelParams



!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module lvalues
        use precision
        use ModelParams
        implicit none
        public

        Type lSamples
            integer l0
            integer l(lmax_arr)
        end Type lSamples

        Type(lSamples) :: lSamp

       contains


        subroutine initlval(lSet,max_l)

! This subroutines initializes lSet%l arrays. Other values will be interpolated.

        implicit none
        type(lSamples) :: lSet

        integer, intent(IN) :: max_l
        integer lind, lvar, step,top,bot,ls(lmax_arr)
        real(dl) AScale

        Ascale=scale/lSampleBoost

        if (lSampleBoost >=50) then
         !just do all of them
         lind=0
         do lvar=lmin, max_l
           lind=lind+1
           ls(lind)=lvar
         end do
         lSet%l0=lind
         lSet%l(1:lind) = ls(1:lind)
         return
        end if

        lind=0
        do lvar=lmin, 10
           lind=lind+1
           ls(lind)=lvar
        end do

        if (CP%AccurateReionization) then
            if (lSampleBoost > 1) then
             do lvar=11, 37,1
               lind=lind+1
               ls(lind)=lvar
             end do
            else
             do lvar=11, 37,2
               lind=lind+1
               ls(lind)=lvar
             end do
            end if

            step = max(nint(5*Ascale),2)
            bot=40
            top=bot + step*10
        else

            if (lSampleBoost >1) then
             do lvar=11, 15
               lind=lind+1
               ls(lind)=lvar
             end do
            else
             lind=lind+1
             ls(lind)=12
             lind=lind+1
             ls(lind)=15
            end if
            step = max(nint(10*Ascale),3)
            bot=15+max(step/2,2)
            top=bot + step*7
        end if

        do lvar=bot, top, step
           lind=lind+1
           ls(lind)=lvar
        end do

        step=max(nint(20*Ascale),4)
        bot=ls(lind)+step
        top=bot+step*2

        do lvar = bot,top,step
          lind=lind+1
          ls(lind)=lvar
        end do

        if (ls(lind)>=max_l) then
           do lvar=lind,1,-1
            if (ls(lvar)<=max_l) exit
           end do
           lind=lvar
           if (ls(lind)<max_l) then
              lind=lind+1
              ls(lind)=max_l
           end if
        else

        step=max(nint(25*Ascale),4)
!Get EE right around l=200 by putting extra point at 175
        bot=ls(lind)+step
        top=bot+step

        do lvar = bot,top,step
          lind=lind+1
          ls(lind)=lvar
        end do


        if (ls(lind)>=max_l) then
           do lvar=lind,1,-1
            if (ls(lvar)<=max_l) exit
           end do
           lind=lvar
           if (ls(lind)<max_l) then
              lind=lind+1
              ls(lind)=max_l
           end if
        else

        if (HighAccuracyDefault .and. .not. use_spline_template) then
         step=max(nint(42*Ascale),7)
        else
         step=max(nint(50*Ascale),7)
        end if
        bot=ls(lind)+step
        top=min(5000,max_l)

         do lvar = bot,top,step
          lind=lind+1
          ls(lind)=lvar
         end do

         if (max_l > 5000) then
             !Should be pretty smooth or tiny out here
             step=max(nint(400*Ascale),50)
             lvar = ls(lind)

             do
              lvar = lvar + step
              if (lvar > max_l) exit
              lind=lind+1
              ls(lind)=lvar
              step = nint(step*1.5) !log spacing
             end do

         end if

         if (ls(lind) /=max_l) then
           lind=lind+1
           ls(lind)=max_l
         end if
        if (.not. CP%flat) ls(lind-1)=int(max_l+ls(lind-2))/2
        !Not in CP%flat case so interpolation table is the same when using lower l_max
        end if
        end if
        lSet%l0=lind
        lSet%l(1:lind) = ls(1:lind)

      end subroutine initlval

      subroutine InterpolateClArr(lSet,iCl, all_Cl, max_ind)
      type (lSamples), intent(in) :: lSet
      real(dl), intent(in) :: iCl(*)
      real(dl), intent(out):: all_Cl(lmin:*)
      integer, intent(in) :: max_ind
      integer il,llo,lhi, xi
      real(dl) ddCl(lSet%l0)
      real(dl) xl(lSet%l0)

      real(dl) a0,b0,ho
      real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

      if (max_ind > lSet%l0) stop 'Wrong max_ind in InterpolateClArr'

      xl = real(lSet%l(1:lSet%l0),dl)
      call spline(xl,iCL(1),max_ind,cllo,clhi,ddCl(1))

            llo=1
            do il=lmin,lSet%l(max_ind)
               xi=il
               if ((xi > lSet%l(llo+1)).and.(llo < max_ind)) then
                  llo=llo+1
               end if
               lhi=llo+1
               ho=lSet%l(lhi)-lSet%l(llo)
               a0=(lSet%l(lhi)-xi)/ho
               b0=(xi-lSet%l(llo))/ho

               all_Cl(il) = a0*iCl(llo)+ b0*iCl(lhi)+((a0**3-a0)* ddCl(llo) &
                       +(b0**3-b0)*ddCl(lhi))*ho**2/6

            end do

      end subroutine InterpolateClArr

       subroutine InterpolateClArrTemplated(lSet,iCl, all_Cl, max_ind, template_index)
      type (lSamples), intent(in) :: lSet
      real(dl), intent(in) :: iCl(*)
      real(dl), intent(out):: all_Cl(lmin:*)
      integer, intent(in) :: max_ind
      integer, intent(in), optional :: template_index
      integer maxdelta, il
      real(dl) DeltaCL(lSet%l0)
      real(dl), allocatable :: tmpall(:)

      if (max_ind > lSet%l0) stop 'Wrong max_ind in InterpolateClArrTemplated'

      if (use_spline_template .and. present(template_index)) then
        if (template_index<=3) then
          !interpolate only the difference between the C_l and an accurately interpolated template. Temp only for the mo.
          !Using unlensed for template, seems to be good enough
           maxdelta=max_ind
           do while (lSet%l(maxdelta) > lmax_extrap_highl)
            maxdelta=maxdelta-1
           end do
           DeltaCL(1:maxdelta)=iCL(1:maxdelta)- highL_CL_template(lSet%l(1:maxdelta), template_index)

           call InterpolateClArr(lSet,DeltaCl, all_Cl, maxdelta)

           do il=lmin,lSet%l(maxdelta)
                all_Cl(il) = all_Cl(il) +  highL_CL_template(il,template_index)
           end do

           if (maxdelta < max_ind) then
           !directly interpolate high L where no template (doesn't effect lensing spectrum much anyway)
            allocate(tmpall(lmin:lSet%l(max_ind)))
            call InterpolateClArr(lSet,iCl, tmpall, max_ind)
            !overlap to reduce interpolation artefacts
            all_cl(lSet%l(maxdelta-2):lSet%l(max_ind) ) = tmpall(lSet%l(maxdelta-2):lSet%l(max_ind))
            deallocate(tmpall)
           end if
          return
        end if
       end if

       call InterpolateClArr(lSet,iCl, all_Cl, max_ind)


      end subroutine InterpolateClArrTemplated





!ccccccccccccccccccccccccccc

        end module lvalues



!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module ModelData
        use precision
        use ModelParams
        use InitialPower
        use lValues
        use Ranges
        use AMlUtils
        implicit none
        public

         Type ClTransferData
      !Cl transfer function variables
       !values of q for integration over q to get C_ls
          Type (lSamples) :: ls ! scalar and tensor l that are computed
          integer :: NumSources
          !Changes -scalars:  2 for just CMB, 3 for lensing
          !- tensors: T and E and phi (for lensing), and T, E, B respectively

          Type (Regions) :: q

          real(dl), dimension(:,:,:), pointer :: Delta_p_l_k => NULL()

         end Type ClTransferData

         Type(ClTransferData), save, target :: CTransScal, CTransTens, CTransVec

        !Computed output power spectra data

        integer, parameter :: C_Temp = 1, C_E = 2, C_Cross =3, C_Phi = 4, C_PhiTemp = 5, C_PhiE=6
        integer :: C_last = C_PhiE
        integer, parameter :: CT_Temp =1, CT_E = 2, CT_B = 3, CT_Cross=  4


        real(dl), dimension (:,:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
        !Indices are Cl_xxx( l , intial_power_index, Cl_type)
        !where Cl_type is one of the above constants

        !The following are set only if doing lensing
        integer lmax_lensed !Only accurate to rather less than this
        real(dl) , dimension (:,:,:), allocatable :: Cl_lensed
          !Cl_lensed(l, power_index, Cl_type) are the interpolated Cls

        real(dl), dimension (:), allocatable ::  COBElikelihoods,COBE_scales
        !Set by COBEnormalize if using outCOBE
        contains


        subroutine Init_ClTransfer(CTrans)
        !Need to set the Ranges array q before calling this
          Type(ClTransferData) :: CTrans
          integer st

          deallocate(CTrans%Delta_p_l_k, STAT = st)
          call Ranges_getArray(CTrans%q, .true.)

          allocate(CTrans%Delta_p_l_k(CTrans%NumSources,min(max_bessels_l_index,CTrans%ls%l0), CTrans%q%npoints))
          CTrans%Delta_p_l_k = 0


         end subroutine Init_ClTransfer


        subroutine Free_ClTransfer(CTrans)
          Type(ClTransferData) :: CTrans
          integer st

           deallocate(CTrans%Delta_p_l_k, STAT = st)
           nullify(CTrans%Delta_p_l_k)
           call Ranges_Free(CTrans%q)

        end subroutine Free_ClTransfer

        subroutine CheckLoadedHighLTemplate
          integer L
          real(dl) array(7)

         if (.not. allocated(highL_CL_template)) then

             highL_CL_template = reshape((/ &
     2.00000E+00 3.00000E+00 4.00000E+00 5.00000E+00 6.00000E+00 7.00000E+00 8.00000E+00 9.00000E+00 &
     1.00000E+01 1.10000E+01 1.20000E+01 1.30000E+01 1.40000E+01 1.50000E+01 1.60000E+01 1.70000E+01 &
     1.80000E+01 1.90000E+01 2.00000E+01 2.10000E+01 2.20000E+01 2.30000E+01 2.40000E+01 2.50000E+01 &
     2.60000E+01 2.70000E+01 2.80000E+01 2.90000E+01 3.00000E+01 3.10000E+01 3.20000E+01 3.30000E+01 &
     3.40000E+01 3.50000E+01 3.60000E+01 3.70000E+01 3.80000E+01 3.90000E+01 4.00000E+01 4.10000E+01 &
     4.20000E+01 4.30000E+01 4.40000E+01 4.50000E+01 4.60000E+01 4.70000E+01 4.80000E+01 4.90000E+01 &
     5.00000E+01 5.10000E+01 5.20000E+01 5.30000E+01 5.40000E+01 5.50000E+01 5.60000E+01 5.70000E+01 &
     5.80000E+01 5.90000E+01 6.00000E+01 6.10000E+01 6.20000E+01 6.30000E+01 6.40000E+01 6.50000E+01 &
     6.60000E+01 6.70000E+01 6.80000E+01 6.90000E+01 7.00000E+01 7.10000E+01 7.20000E+01 7.30000E+01 &
     7.40000E+01 7.50000E+01 7.60000E+01 7.70000E+01 7.80000E+01 7.90000E+01 8.00000E+01 8.10000E+01 &
     8.20000E+01 8.30000E+01 8.40000E+01 8.50000E+01 8.60000E+01 8.70000E+01 8.80000E+01 8.90000E+01 &
     9.00000E+01 9.10000E+01 9.20000E+01 9.30000E+01 9.40000E+01 9.50000E+01 9.60000E+01 9.70000E+01 &
     9.80000E+01 9.90000E+01 1.00000E+02 1.01000E+02 1.02000E+02 1.03000E+02 1.04000E+02 1.05000E+02 &
     1.06000E+02 1.07000E+02 1.08000E+02 1.09000E+02 1.10000E+02 1.11000E+02 1.12000E+02 1.13000E+02 &
     1.14000E+02 1.15000E+02 1.16000E+02 1.17000E+02 1.18000E+02 1.19000E+02 1.20000E+02 1.21000E+02 &
     1.22000E+02 1.23000E+02 1.24000E+02 1.25000E+02 1.26000E+02 1.27000E+02 1.28000E+02 1.29000E+02 &
     1.30000E+02 1.31000E+02 1.32000E+02 1.33000E+02 1.34000E+02 1.35000E+02 1.36000E+02 1.37000E+02 &
     1.38000E+02 1.39000E+02 1.40000E+02 1.41000E+02 1.42000E+02 1.43000E+02 1.44000E+02 1.45000E+02 &
     1.46000E+02 1.47000E+02 1.48000E+02 1.49000E+02 1.50000E+02 1.51000E+02 1.52000E+02 1.53000E+02 &
     1.54000E+02 1.55000E+02 1.56000E+02 1.57000E+02 1.58000E+02 1.59000E+02 1.60000E+02 1.61000E+02 &
     1.62000E+02 1.63000E+02 1.64000E+02 1.65000E+02 1.66000E+02 1.67000E+02 1.68000E+02 1.69000E+02 &
     1.70000E+02 1.71000E+02 1.72000E+02 1.73000E+02 1.74000E+02 1.75000E+02 1.76000E+02 1.77000E+02 &
     1.78000E+02 1.79000E+02 1.80000E+02 1.81000E+02 1.82000E+02 1.83000E+02 1.84000E+02 1.85000E+02 &
     1.86000E+02 1.87000E+02 1.88000E+02 1.89000E+02 1.90000E+02 1.91000E+02 1.92000E+02 1.93000E+02 &
     1.94000E+02 1.95000E+02 1.96000E+02 1.97000E+02 1.98000E+02 1.99000E+02 2.00000E+02 2.01000E+02 &
     2.02000E+02 2.03000E+02 2.04000E+02 2.05000E+02 2.06000E+02 2.07000E+02 2.08000E+02 2.09000E+02 &
     2.10000E+02 2.11000E+02 2.12000E+02 2.13000E+02 2.14000E+02 2.15000E+02 2.16000E+02 2.17000E+02 &
     2.18000E+02 2.19000E+02 2.20000E+02 2.21000E+02 2.22000E+02 2.23000E+02 2.24000E+02 2.25000E+02 &
     2.26000E+02 2.27000E+02 2.28000E+02 2.29000E+02 2.30000E+02 2.31000E+02 2.32000E+02 2.33000E+02 &
     2.34000E+02 2.35000E+02 2.36000E+02 2.37000E+02 2.38000E+02 2.39000E+02 2.40000E+02 2.41000E+02 &
     2.42000E+02 2.43000E+02 2.44000E+02 2.45000E+02 2.46000E+02 2.47000E+02 2.48000E+02 2.49000E+02 &
     2.50000E+02 2.51000E+02 2.52000E+02 2.53000E+02 2.54000E+02 2.55000E+02 2.56000E+02 2.57000E+02 &
     2.58000E+02 2.59000E+02 2.60000E+02 2.61000E+02 2.62000E+02 2.63000E+02 2.64000E+02 2.65000E+02 &
     2.66000E+02 2.67000E+02 2.68000E+02 2.69000E+02 2.70000E+02 2.71000E+02 2.72000E+02 2.73000E+02 &
     2.74000E+02 2.75000E+02 2.76000E+02 2.77000E+02 2.78000E+02 2.79000E+02 2.80000E+02 2.81000E+02 &
     2.82000E+02 2.83000E+02 2.84000E+02 2.85000E+02 2.86000E+02 2.87000E+02 2.88000E+02 2.89000E+02 &
     2.90000E+02 2.91000E+02 2.92000E+02 2.93000E+02 2.94000E+02 2.95000E+02 2.96000E+02 2.97000E+02 &
     2.98000E+02 2.99000E+02 3.00000E+02 3.01000E+02 3.02000E+02 3.03000E+02 3.04000E+02 3.05000E+02 &
     3.06000E+02 3.07000E+02 3.08000E+02 3.09000E+02 3.10000E+02 3.11000E+02 3.12000E+02 3.13000E+02 &
     3.14000E+02 3.15000E+02 3.16000E+02 3.17000E+02 3.18000E+02 3.19000E+02 3.20000E+02 3.21000E+02 &
     3.22000E+02 3.23000E+02 3.24000E+02 3.25000E+02 3.26000E+02 3.27000E+02 3.28000E+02 3.29000E+02 &
     3.30000E+02 3.31000E+02 3.32000E+02 3.33000E+02 3.34000E+02 3.35000E+02 3.36000E+02 3.37000E+02 &
     3.38000E+02 3.39000E+02 3.40000E+02 3.41000E+02 3.42000E+02 3.43000E+02 3.44000E+02 3.45000E+02 &
     3.46000E+02 3.47000E+02 3.48000E+02 3.49000E+02 3.50000E+02 3.51000E+02 3.52000E+02 3.53000E+02 &
     3.54000E+02 3.55000E+02 3.56000E+02 3.57000E+02 3.58000E+02 3.59000E+02 3.60000E+02 3.61000E+02 &
     3.62000E+02 3.63000E+02 3.64000E+02 3.65000E+02 3.66000E+02 3.67000E+02 3.68000E+02 3.69000E+02 &
     3.70000E+02 3.71000E+02 3.72000E+02 3.73000E+02 3.74000E+02 3.75000E+02 3.76000E+02 3.77000E+02 &
     3.78000E+02 3.79000E+02 3.80000E+02 3.81000E+02 3.82000E+02 3.83000E+02 3.84000E+02 3.85000E+02 &
     3.86000E+02 3.87000E+02 3.88000E+02 3.89000E+02 3.90000E+02 3.91000E+02 3.92000E+02 3.93000E+02 &
     3.94000E+02 3.95000E+02 3.96000E+02 3.97000E+02 3.98000E+02 3.99000E+02 4.00000E+02 4.01000E+02 &
     4.02000E+02 4.03000E+02 4.04000E+02 4.05000E+02 4.06000E+02 4.07000E+02 4.08000E+02 4.09000E+02 &
     4.10000E+02 4.11000E+02 4.12000E+02 4.13000E+02 4.14000E+02 4.15000E+02 4.16000E+02 4.17000E+02 &
     4.18000E+02 4.19000E+02 4.20000E+02 4.21000E+02 4.22000E+02 4.23000E+02 4.24000E+02 4.25000E+02 &
     4.26000E+02 4.27000E+02 4.28000E+02 4.29000E+02 4.30000E+02 4.31000E+02 4.32000E+02 4.33000E+02 &
     4.34000E+02 4.35000E+02 4.36000E+02 4.37000E+02 4.38000E+02 4.39000E+02 4.40000E+02 4.41000E+02 &
     4.42000E+02 4.43000E+02 4.44000E+02 4.45000E+02 4.46000E+02 4.47000E+02 4.48000E+02 4.49000E+02 &
     4.50000E+02 4.51000E+02 4.52000E+02 4.53000E+02 4.54000E+02 4.55000E+02 4.56000E+02 4.57000E+02 &
     4.58000E+02 4.59000E+02 4.60000E+02 4.61000E+02 4.62000E+02 4.63000E+02 4.64000E+02 4.65000E+02 &
     4.66000E+02 4.67000E+02 4.68000E+02 4.69000E+02 4.70000E+02 4.71000E+02 4.72000E+02 4.73000E+02 &
     4.74000E+02 4.75000E+02 4.76000E+02 4.77000E+02 4.78000E+02 4.79000E+02 4.80000E+02 4.81000E+02 &
     4.82000E+02 4.83000E+02 4.84000E+02 4.85000E+02 4.86000E+02 4.87000E+02 4.88000E+02 4.89000E+02 &
     4.90000E+02 4.91000E+02 4.92000E+02 4.93000E+02 4.94000E+02 4.95000E+02 4.96000E+02 4.97000E+02 &
     4.98000E+02 4.99000E+02 5.00000E+02 5.01000E+02 5.02000E+02 5.03000E+02 5.04000E+02 5.05000E+02 &
     5.06000E+02 5.07000E+02 5.08000E+02 5.09000E+02 5.10000E+02 5.11000E+02 5.12000E+02 5.13000E+02 &
     5.14000E+02 5.15000E+02 5.16000E+02 5.17000E+02 5.18000E+02 5.19000E+02 5.20000E+02 5.21000E+02 &
     5.22000E+02 5.23000E+02 5.24000E+02 5.25000E+02 5.26000E+02 5.27000E+02 5.28000E+02 5.29000E+02 &
     5.30000E+02 5.31000E+02 5.32000E+02 5.33000E+02 5.34000E+02 5.35000E+02 5.36000E+02 5.37000E+02 &
     5.38000E+02 5.39000E+02 5.40000E+02 5.41000E+02 5.42000E+02 5.43000E+02 5.44000E+02 5.45000E+02 &
     5.46000E+02 5.47000E+02 5.48000E+02 5.49000E+02 5.50000E+02 5.51000E+02 5.52000E+02 5.53000E+02 &
     5.54000E+02 5.55000E+02 5.56000E+02 5.57000E+02 5.58000E+02 5.59000E+02 5.60000E+02 5.61000E+02 &
     5.62000E+02 5.63000E+02 5.64000E+02 5.65000E+02 5.66000E+02 5.67000E+02 5.68000E+02 5.69000E+02 &
     5.70000E+02 5.71000E+02 5.72000E+02 5.73000E+02 5.74000E+02 5.75000E+02 5.76000E+02 5.77000E+02 &
     5.78000E+02 5.79000E+02 5.80000E+02 5.81000E+02 5.82000E+02 5.83000E+02 5.84000E+02 5.85000E+02 &
     5.86000E+02 5.87000E+02 5.88000E+02 5.89000E+02 5.90000E+02 5.91000E+02 5.92000E+02 5.93000E+02 &
     5.94000E+02 5.95000E+02 5.96000E+02 5.97000E+02 5.98000E+02 5.99000E+02 6.00000E+02 6.01000E+02 &
     6.02000E+02 6.03000E+02 6.04000E+02 6.05000E+02 6.06000E+02 6.07000E+02 6.08000E+02 6.09000E+02 &
     6.10000E+02 6.11000E+02 6.12000E+02 6.13000E+02 6.14000E+02 6.15000E+02 6.16000E+02 6.17000E+02 &
     6.18000E+02 6.19000E+02 6.20000E+02 6.21000E+02 6.22000E+02 6.23000E+02 6.24000E+02 6.25000E+02 &
     6.26000E+02 6.27000E+02 6.28000E+02 6.29000E+02 6.30000E+02 6.31000E+02 6.32000E+02 6.33000E+02 &
     6.34000E+02 6.35000E+02 6.36000E+02 6.37000E+02 6.38000E+02 6.39000E+02 6.40000E+02 6.41000E+02 &
     6.42000E+02 6.43000E+02 6.44000E+02 6.45000E+02 6.46000E+02 6.47000E+02 6.48000E+02 6.49000E+02 &
     6.50000E+02 6.51000E+02 6.52000E+02 6.53000E+02 6.54000E+02 6.55000E+02 6.56000E+02 6.57000E+02 &
     6.58000E+02 6.59000E+02 6.60000E+02 6.61000E+02 6.62000E+02 6.63000E+02 6.64000E+02 6.65000E+02 &
     6.66000E+02 6.67000E+02 6.68000E+02 6.69000E+02 6.70000E+02 6.71000E+02 6.72000E+02 6.73000E+02 &
     6.74000E+02 6.75000E+02 6.76000E+02 6.77000E+02 6.78000E+02 6.79000E+02 6.80000E+02 6.81000E+02 &
     6.82000E+02 6.83000E+02 6.84000E+02 6.85000E+02 6.86000E+02 6.87000E+02 6.88000E+02 6.89000E+02 &
     6.90000E+02 6.91000E+02 6.92000E+02 6.93000E+02 6.94000E+02 6.95000E+02 6.96000E+02 6.97000E+02 &
     6.98000E+02 6.99000E+02 7.00000E+02 7.01000E+02 7.02000E+02 7.03000E+02 7.04000E+02 7.05000E+02 &
     7.06000E+02 7.07000E+02 7.08000E+02 7.09000E+02 7.10000E+02 7.11000E+02 7.12000E+02 7.13000E+02 &
     7.14000E+02 7.15000E+02 7.16000E+02 7.17000E+02 7.18000E+02 7.19000E+02 7.20000E+02 7.21000E+02 &
     7.22000E+02 7.23000E+02 7.24000E+02 7.25000E+02 7.26000E+02 7.27000E+02 7.28000E+02 7.29000E+02 &
     7.30000E+02 7.31000E+02 7.32000E+02 7.33000E+02 7.34000E+02 7.35000E+02 7.36000E+02 7.37000E+02 &
     7.38000E+02 7.39000E+02 7.40000E+02 7.41000E+02 7.42000E+02 7.43000E+02 7.44000E+02 7.45000E+02 &
     7.46000E+02 7.47000E+02 7.48000E+02 7.49000E+02 7.50000E+02 7.51000E+02 7.52000E+02 7.53000E+02 &
     7.54000E+02 7.55000E+02 7.56000E+02 7.57000E+02 7.58000E+02 7.59000E+02 7.60000E+02 7.61000E+02 &
     7.62000E+02 7.63000E+02 7.64000E+02 7.65000E+02 7.66000E+02 7.67000E+02 7.68000E+02 7.69000E+02 &
     7.70000E+02 7.71000E+02 7.72000E+02 7.73000E+02 7.74000E+02 7.75000E+02 7.76000E+02 7.77000E+02 &
     7.78000E+02 7.79000E+02 7.80000E+02 7.81000E+02 7.82000E+02 7.83000E+02 7.84000E+02 7.85000E+02 &
     7.86000E+02 7.87000E+02 7.88000E+02 7.89000E+02 7.90000E+02 7.91000E+02 7.92000E+02 7.93000E+02 &
     7.94000E+02 7.95000E+02 7.96000E+02 7.97000E+02 7.98000E+02 7.99000E+02 8.00000E+02 8.01000E+02 &
     8.02000E+02 8.03000E+02 8.04000E+02 8.05000E+02 8.06000E+02 8.07000E+02 8.08000E+02 8.09000E+02 &
     8.10000E+02 8.11000E+02 8.12000E+02 8.13000E+02 8.14000E+02 8.15000E+02 8.16000E+02 8.17000E+02 &
     8.18000E+02 8.19000E+02 8.20000E+02 8.21000E+02 8.22000E+02 8.23000E+02 8.24000E+02 8.25000E+02 &
     8.26000E+02 8.27000E+02 8.28000E+02 8.29000E+02 8.30000E+02 8.31000E+02 8.32000E+02 8.33000E+02 &
     8.34000E+02 8.35000E+02 8.36000E+02 8.37000E+02 8.38000E+02 8.39000E+02 8.40000E+02 8.41000E+02 &
     8.42000E+02 8.43000E+02 8.44000E+02 8.45000E+02 8.46000E+02 8.47000E+02 8.48000E+02 8.49000E+02 &
     8.50000E+02 8.51000E+02 8.52000E+02 8.53000E+02 8.54000E+02 8.55000E+02 8.56000E+02 8.57000E+02 &
     8.58000E+02 8.59000E+02 8.60000E+02 8.61000E+02 8.62000E+02 8.63000E+02 8.64000E+02 8.65000E+02 &
     8.66000E+02 8.67000E+02 8.68000E+02 8.69000E+02 8.70000E+02 8.71000E+02 8.72000E+02 8.73000E+02 &
     8.74000E+02 8.75000E+02 8.76000E+02 8.77000E+02 8.78000E+02 8.79000E+02 8.80000E+02 8.81000E+02 &
     8.82000E+02 8.83000E+02 8.84000E+02 8.85000E+02 8.86000E+02 8.87000E+02 8.88000E+02 8.89000E+02 &
     8.90000E+02 8.91000E+02 8.92000E+02 8.93000E+02 8.94000E+02 8.95000E+02 8.96000E+02 8.97000E+02 &
     8.98000E+02 8.99000E+02 9.00000E+02 9.01000E+02 9.02000E+02 9.03000E+02 9.04000E+02 9.05000E+02 &
     9.06000E+02 9.07000E+02 9.08000E+02 9.09000E+02 9.10000E+02 9.11000E+02 9.12000E+02 9.13000E+02 &
     9.14000E+02 9.15000E+02 9.16000E+02 9.17000E+02 9.18000E+02 9.19000E+02 9.20000E+02 9.21000E+02 &
     9.22000E+02 9.23000E+02 9.24000E+02 9.25000E+02 9.26000E+02 9.27000E+02 9.28000E+02 9.29000E+02 &
     9.30000E+02 9.31000E+02 9.32000E+02 9.33000E+02 9.34000E+02 9.35000E+02 9.36000E+02 9.37000E+02 &
     9.38000E+02 9.39000E+02 9.40000E+02 9.41000E+02 9.42000E+02 9.43000E+02 9.44000E+02 9.45000E+02 &
     9.46000E+02 9.47000E+02 9.48000E+02 9.49000E+02 9.50000E+02 9.51000E+02 9.52000E+02 9.53000E+02 &
     9.54000E+02 9.55000E+02 9.56000E+02 9.57000E+02 9.58000E+02 9.59000E+02 9.60000E+02 9.61000E+02 &
     9.62000E+02 9.63000E+02 9.64000E+02 9.65000E+02 9.66000E+02 9.67000E+02 9.68000E+02 9.69000E+02 &
     9.70000E+02 9.71000E+02 9.72000E+02 9.73000E+02 9.74000E+02 9.75000E+02 9.76000E+02 9.77000E+02 &
     9.78000E+02 9.79000E+02 9.80000E+02 9.81000E+02 9.82000E+02 9.83000E+02 9.84000E+02 9.85000E+02 &
     9.86000E+02 9.87000E+02 9.88000E+02 9.89000E+02 9.90000E+02 9.91000E+02 9.92000E+02 9.93000E+02 &
     9.94000E+02 9.95000E+02 9.96000E+02 9.97000E+02 9.98000E+02 9.99000E+02 1.00000E+03 1.00100E+03 &
     1.00200E+03 1.00300E+03 1.00400E+03 1.00500E+03 1.00600E+03 1.00700E+03 1.00800E+03 1.00900E+03 &
     1.01000E+03 1.01100E+03 1.01200E+03 1.01300E+03 1.01400E+03 1.01500E+03 1.01600E+03 1.01700E+03 &
     1.01800E+03 1.01900E+03 1.02000E+03 1.02100E+03 1.02200E+03 1.02300E+03 1.02400E+03 1.02500E+03 &
     1.02600E+03 1.02700E+03 1.02800E+03 1.02900E+03 1.03000E+03 1.03100E+03 1.03200E+03 1.03300E+03 &
     1.03400E+03 1.03500E+03 1.03600E+03 1.03700E+03 1.03800E+03 1.03900E+03 1.04000E+03 1.04100E+03 &
     1.04200E+03 1.04300E+03 1.04400E+03 1.04500E+03 1.04600E+03 1.04700E+03 1.04800E+03 1.04900E+03 &
     1.05000E+03 1.05100E+03 1.05200E+03 1.05300E+03 1.05400E+03 1.05500E+03 1.05600E+03 1.05700E+03 &
     1.05800E+03 1.05900E+03 1.06000E+03 1.06100E+03 1.06200E+03 1.06300E+03 1.06400E+03 1.06500E+03 &
     1.06600E+03 1.06700E+03 1.06800E+03 1.06900E+03 1.07000E+03 1.07100E+03 1.07200E+03 1.07300E+03 &
     1.07400E+03 1.07500E+03 1.07600E+03 1.07700E+03 1.07800E+03 1.07900E+03 1.08000E+03 1.08100E+03 &
     1.08200E+03 1.08300E+03 1.08400E+03 1.08500E+03 1.08600E+03 1.08700E+03 1.08800E+03 1.08900E+03 &
     1.09000E+03 1.09100E+03 1.09200E+03 1.09300E+03 1.09400E+03 1.09500E+03 1.09600E+03 1.09700E+03 &
     1.09800E+03 1.09900E+03 1.10000E+03 1.10100E+03 1.10200E+03 1.10300E+03 1.10400E+03 1.10500E+03 &
     1.10600E+03 1.10700E+03 1.10800E+03 1.10900E+03 1.11000E+03 1.11100E+03 1.11200E+03 1.11300E+03 &
     1.11400E+03 1.11500E+03 1.11600E+03 1.11700E+03 1.11800E+03 1.11900E+03 1.12000E+03 1.12100E+03 &
     1.12200E+03 1.12300E+03 1.12400E+03 1.12500E+03 1.12600E+03 1.12700E+03 1.12800E+03 1.12900E+03 &
     1.13000E+03 1.13100E+03 1.13200E+03 1.13300E+03 1.13400E+03 1.13500E+03 1.13600E+03 1.13700E+03 &
     1.13800E+03 1.13900E+03 1.14000E+03 1.14100E+03 1.14200E+03 1.14300E+03 1.14400E+03 1.14500E+03 &
     1.14600E+03 1.14700E+03 1.14800E+03 1.14900E+03 1.15000E+03 1.15100E+03 1.15200E+03 1.15300E+03 &
     1.15400E+03 1.15500E+03 1.15600E+03 1.15700E+03 1.15800E+03 1.15900E+03 1.16000E+03 1.16100E+03 &
     1.16200E+03 1.16300E+03 1.16400E+03 1.16500E+03 1.16600E+03 1.16700E+03 1.16800E+03 1.16900E+03 &
     1.17000E+03 1.17100E+03 1.17200E+03 1.17300E+03 1.17400E+03 1.17500E+03 1.17600E+03 1.17700E+03 &
     1.17800E+03 1.17900E+03 1.18000E+03 1.18100E+03 1.18200E+03 1.18300E+03 1.18400E+03 1.18500E+03 &
     1.18600E+03 1.18700E+03 1.18800E+03 1.18900E+03 1.19000E+03 1.19100E+03 1.19200E+03 1.19300E+03 &
     1.19400E+03 1.19500E+03 1.19600E+03 1.19700E+03 1.19800E+03 1.19900E+03 1.20000E+03 1.20100E+03 &
     1.20200E+03 1.20300E+03 1.20400E+03 1.20500E+03 1.20600E+03 1.20700E+03 1.20800E+03 1.20900E+03 &
     1.21000E+03 1.21100E+03 1.21200E+03 1.21300E+03 1.21400E+03 1.21500E+03 1.21600E+03 1.21700E+03 &
     1.21800E+03 1.21900E+03 1.22000E+03 1.22100E+03 1.22200E+03 1.22300E+03 1.22400E+03 1.22500E+03 &
     1.22600E+03 1.22700E+03 1.22800E+03 1.22900E+03 1.23000E+03 1.23100E+03 1.23200E+03 1.23300E+03 &
     1.23400E+03 1.23500E+03 1.23600E+03 1.23700E+03 1.23800E+03 1.23900E+03 1.24000E+03 1.24100E+03 &
     1.24200E+03 1.24300E+03 1.24400E+03 1.24500E+03 1.24600E+03 1.24700E+03 1.24800E+03 1.24900E+03 &
     1.25000E+03 1.25100E+03 1.25200E+03 1.25300E+03 1.25400E+03 1.25500E+03 1.25600E+03 1.25700E+03 &
     1.25800E+03 1.25900E+03 1.26000E+03 1.26100E+03 1.26200E+03 1.26300E+03 1.26400E+03 1.26500E+03 &
     1.26600E+03 1.26700E+03 1.26800E+03 1.26900E+03 1.27000E+03 1.27100E+03 1.27200E+03 1.27300E+03 &
     1.27400E+03 1.27500E+03 1.27600E+03 1.27700E+03 1.27800E+03 1.27900E+03 1.28000E+03 1.28100E+03 &
     1.28200E+03 1.28300E+03 1.28400E+03 1.28500E+03 1.28600E+03 1.28700E+03 1.28800E+03 1.28900E+03 &
     1.29000E+03 1.29100E+03 1.29200E+03 1.29300E+03 1.29400E+03 1.29500E+03 1.29600E+03 1.29700E+03 &
     1.29800E+03 1.29900E+03 1.30000E+03 1.30100E+03 1.30200E+03 1.30300E+03 1.30400E+03 1.30500E+03 &
     1.30600E+03 1.30700E+03 1.30800E+03 1.30900E+03 1.31000E+03 1.31100E+03 1.31200E+03 1.31300E+03 &
     1.31400E+03 1.31500E+03 1.31600E+03 1.31700E+03 1.31800E+03 1.31900E+03 1.32000E+03 1.32100E+03 &
     1.32200E+03 1.32300E+03 1.32400E+03 1.32500E+03 1.32600E+03 1.32700E+03 1.32800E+03 1.32900E+03 &
     1.33000E+03 1.33100E+03 1.33200E+03 1.33300E+03 1.33400E+03 1.33500E+03 1.33600E+03 1.33700E+03 &
     1.33800E+03 1.33900E+03 1.34000E+03 1.34100E+03 1.34200E+03 1.34300E+03 1.34400E+03 1.34500E+03 &
     1.34600E+03 1.34700E+03 1.34800E+03 1.34900E+03 1.35000E+03 1.35100E+03 1.35200E+03 1.35300E+03 &
     1.35400E+03 1.35500E+03 1.35600E+03 1.35700E+03 1.35800E+03 1.35900E+03 1.36000E+03 1.36100E+03 &
     1.36200E+03 1.36300E+03 1.36400E+03 1.36500E+03 1.36600E+03 1.36700E+03 1.36800E+03 1.36900E+03 &
     1.37000E+03 1.37100E+03 1.37200E+03 1.37300E+03 1.37400E+03 1.37500E+03 1.37600E+03 1.37700E+03 &
     1.37800E+03 1.37900E+03 1.38000E+03 1.38100E+03 1.38200E+03 1.38300E+03 1.38400E+03 1.38500E+03 &
     1.38600E+03 1.38700E+03 1.38800E+03 1.38900E+03 1.39000E+03 1.39100E+03 1.39200E+03 1.39300E+03 &
     1.39400E+03 1.39500E+03 1.39600E+03 1.39700E+03 1.39800E+03 1.39900E+03 1.40000E+03 1.40100E+03 &
     1.40200E+03 1.40300E+03 1.40400E+03 1.40500E+03 1.40600E+03 1.40700E+03 1.40800E+03 1.40900E+03 &
     1.41000E+03 1.41100E+03 1.41200E+03 1.41300E+03 1.41400E+03 1.41500E+03 1.41600E+03 1.41700E+03 &
     1.41800E+03 1.41900E+03 1.42000E+03 1.42100E+03 1.42200E+03 1.42300E+03 1.42400E+03 1.42500E+03 &
     1.42600E+03 1.42700E+03 1.42800E+03 1.42900E+03 1.43000E+03 1.43100E+03 1.43200E+03 1.43300E+03 &
     1.43400E+03 1.43500E+03 1.43600E+03 1.43700E+03 1.43800E+03 1.43900E+03 1.44000E+03 1.44100E+03 &
     1.44200E+03 1.44300E+03 1.44400E+03 1.44500E+03 1.44600E+03 1.44700E+03 1.44800E+03 1.44900E+03 &
     1.45000E+03 1.45100E+03 1.45200E+03 1.45300E+03 1.45400E+03 1.45500E+03 1.45600E+03 1.45700E+03 &
     1.45800E+03 1.45900E+03 1.46000E+03 1.46100E+03 1.46200E+03 1.46300E+03 1.46400E+03 1.46500E+03 &
     1.46600E+03 1.46700E+03 1.46800E+03 1.46900E+03 1.47000E+03 1.47100E+03 1.47200E+03 1.47300E+03 &
     1.47400E+03 1.47500E+03 1.47600E+03 1.47700E+03 1.47800E+03 1.47900E+03 1.48000E+03 1.48100E+03 &
     1.48200E+03 1.48300E+03 1.48400E+03 1.48500E+03 1.48600E+03 1.48700E+03 1.48800E+03 1.48900E+03 &
     1.49000E+03 1.49100E+03 1.49200E+03 1.49300E+03 1.49400E+03 1.49500E+03 1.49600E+03 1.49700E+03 &
     1.49800E+03 1.49900E+03 1.50000E+03 1.50100E+03 1.50200E+03 1.50300E+03 1.50400E+03 1.50500E+03 &
     1.50600E+03 1.50700E+03 1.50800E+03 1.50900E+03 1.51000E+03 1.51100E+03 1.51200E+03 1.51300E+03 &
     1.51400E+03 1.51500E+03 1.51600E+03 1.51700E+03 1.51800E+03 1.51900E+03 1.52000E+03 1.52100E+03 &
     1.52200E+03 1.52300E+03 1.52400E+03 1.52500E+03 1.52600E+03 1.52700E+03 1.52800E+03 1.52900E+03 &
     1.53000E+03 1.53100E+03 1.53200E+03 1.53300E+03 1.53400E+03 1.53500E+03 1.53600E+03 1.53700E+03 &
     1.53800E+03 1.53900E+03 1.54000E+03 1.54100E+03 1.54200E+03 1.54300E+03 1.54400E+03 1.54500E+03 &
     1.54600E+03 1.54700E+03 1.54800E+03 1.54900E+03 1.55000E+03 1.55100E+03 1.55200E+03 1.55300E+03 &
     1.55400E+03 1.55500E+03 1.55600E+03 1.55700E+03 1.55800E+03 1.55900E+03 1.56000E+03 1.56100E+03 &
     1.56200E+03 1.56300E+03 1.56400E+03 1.56500E+03 1.56600E+03 1.56700E+03 1.56800E+03 1.56900E+03 &
     1.57000E+03 1.57100E+03 1.57200E+03 1.57300E+03 1.57400E+03 1.57500E+03 1.57600E+03 1.57700E+03 &
     1.57800E+03 1.57900E+03 1.58000E+03 1.58100E+03 1.58200E+03 1.58300E+03 1.58400E+03 1.58500E+03 &
     1.58600E+03 1.58700E+03 1.58800E+03 1.58900E+03 1.59000E+03 1.59100E+03 1.59200E+03 1.59300E+03 &
     1.59400E+03 1.59500E+03 1.59600E+03 1.59700E+03 1.59800E+03 1.59900E+03 1.60000E+03 1.60100E+03 &
     1.60200E+03 1.60300E+03 1.60400E+03 1.60500E+03 1.60600E+03 1.60700E+03 1.60800E+03 1.60900E+03 &
     1.61000E+03 1.61100E+03 1.61200E+03 1.61300E+03 1.61400E+03 1.61500E+03 1.61600E+03 1.61700E+03 &
     1.61800E+03 1.61900E+03 1.62000E+03 1.62100E+03 1.62200E+03 1.62300E+03 1.62400E+03 1.62500E+03 &
     1.62600E+03 1.62700E+03 1.62800E+03 1.62900E+03 1.63000E+03 1.63100E+03 1.63200E+03 1.63300E+03 &
     1.63400E+03 1.63500E+03 1.63600E+03 1.63700E+03 1.63800E+03 1.63900E+03 1.64000E+03 1.64100E+03 &
     1.64200E+03 1.64300E+03 1.64400E+03 1.64500E+03 1.64600E+03 1.64700E+03 1.64800E+03 1.64900E+03 &
     1.65000E+03 1.65100E+03 1.65200E+03 1.65300E+03 1.65400E+03 1.65500E+03 1.65600E+03 1.65700E+03 &
     1.65800E+03 1.65900E+03 1.66000E+03 1.66100E+03 1.66200E+03 1.66300E+03 1.66400E+03 1.66500E+03 &
     1.66600E+03 1.66700E+03 1.66800E+03 1.66900E+03 1.67000E+03 1.67100E+03 1.67200E+03 1.67300E+03 &
     1.67400E+03 1.67500E+03 1.67600E+03 1.67700E+03 1.67800E+03 1.67900E+03 1.68000E+03 1.68100E+03 &
     1.68200E+03 1.68300E+03 1.68400E+03 1.68500E+03 1.68600E+03 1.68700E+03 1.68800E+03 1.68900E+03 &
     1.69000E+03 1.69100E+03 1.69200E+03 1.69300E+03 1.69400E+03 1.69500E+03 1.69600E+03 1.69700E+03 &
     1.69800E+03 1.69900E+03 1.70000E+03 1.70100E+03 1.70200E+03 1.70300E+03 1.70400E+03 1.70500E+03 &
     1.70600E+03 1.70700E+03 1.70800E+03 1.70900E+03 1.71000E+03 1.71100E+03 1.71200E+03 1.71300E+03 &
     1.71400E+03 1.71500E+03 1.71600E+03 1.71700E+03 1.71800E+03 1.71900E+03 1.72000E+03 1.72100E+03 &
     1.72200E+03 1.72300E+03 1.72400E+03 1.72500E+03 1.72600E+03 1.72700E+03 1.72800E+03 1.72900E+03 &
     1.73000E+03 1.73100E+03 1.73200E+03 1.73300E+03 1.73400E+03 1.73500E+03 1.73600E+03 1.73700E+03 &
     1.73800E+03 1.73900E+03 1.74000E+03 1.74100E+03 1.74200E+03 1.74300E+03 1.74400E+03 1.74500E+03 &
     1.74600E+03 1.74700E+03 1.74800E+03 1.74900E+03 1.75000E+03 1.75100E+03 1.75200E+03 1.75300E+03 &
     1.75400E+03 1.75500E+03 1.75600E+03 1.75700E+03 1.75800E+03 1.75900E+03 1.76000E+03 1.76100E+03 &
     1.76200E+03 1.76300E+03 1.76400E+03 1.76500E+03 1.76600E+03 1.76700E+03 1.76800E+03 1.76900E+03 &
     1.77000E+03 1.77100E+03 1.77200E+03 1.77300E+03 1.77400E+03 1.77500E+03 1.77600E+03 1.77700E+03 &
     1.77800E+03 1.77900E+03 1.78000E+03 1.78100E+03 1.78200E+03 1.78300E+03 1.78400E+03 1.78500E+03 &
     1.78600E+03 1.78700E+03 1.78800E+03 1.78900E+03 1.79000E+03 1.79100E+03 1.79200E+03 1.79300E+03 &
     1.79400E+03 1.79500E+03 1.79600E+03 1.79700E+03 1.79800E+03 1.79900E+03 1.80000E+03 1.80100E+03 &
     1.80200E+03 1.80300E+03 1.80400E+03 1.80500E+03 1.80600E+03 1.80700E+03 1.80800E+03 1.80900E+03 &
     1.81000E+03 1.81100E+03 1.81200E+03 1.81300E+03 1.81400E+03 1.81500E+03 1.81600E+03 1.81700E+03 &
     1.81800E+03 1.81900E+03 1.82000E+03 1.82100E+03 1.82200E+03 1.82300E+03 1.82400E+03 1.82500E+03 &
     1.82600E+03 1.82700E+03 1.82800E+03 1.82900E+03 1.83000E+03 1.83100E+03 1.83200E+03 1.83300E+03 &
     1.83400E+03 1.83500E+03 1.83600E+03 1.83700E+03 1.83800E+03 1.83900E+03 1.84000E+03 1.84100E+03 &
     1.84200E+03 1.84300E+03 1.84400E+03 1.84500E+03 1.84600E+03 1.84700E+03 1.84800E+03 1.84900E+03 &
     1.85000E+03 1.85100E+03 1.85200E+03 1.85300E+03 1.85400E+03 1.85500E+03 1.85600E+03 1.85700E+03 &
     1.85800E+03 1.85900E+03 1.86000E+03 1.86100E+03 1.86200E+03 1.86300E+03 1.86400E+03 1.86500E+03 &
     1.86600E+03 1.86700E+03 1.86800E+03 1.86900E+03 1.87000E+03 1.87100E+03 1.87200E+03 1.87300E+03 &
     1.87400E+03 1.87500E+03 1.87600E+03 1.87700E+03 1.87800E+03 1.87900E+03 1.88000E+03 1.88100E+03 &
     1.88200E+03 1.88300E+03 1.88400E+03 1.88500E+03 1.88600E+03 1.88700E+03 1.88800E+03 1.88900E+03 &
     1.89000E+03 1.89100E+03 1.89200E+03 1.89300E+03 1.89400E+03 1.89500E+03 1.89600E+03 1.89700E+03 &
     1.89800E+03 1.89900E+03 1.90000E+03 1.90100E+03 1.90200E+03 1.90300E+03 1.90400E+03 1.90500E+03 &
     1.90600E+03 1.90700E+03 1.90800E+03 1.90900E+03 1.91000E+03 1.91100E+03 1.91200E+03 1.91300E+03 &
     1.91400E+03 1.91500E+03 1.91600E+03 1.91700E+03 1.91800E+03 1.91900E+03 1.92000E+03 1.92100E+03 &
     1.92200E+03 1.92300E+03 1.92400E+03 1.92500E+03 1.92600E+03 1.92700E+03 1.92800E+03 1.92900E+03 &
     1.93000E+03 1.93100E+03 1.93200E+03 1.93300E+03 1.93400E+03 1.93500E+03 1.93600E+03 1.93700E+03 &
     1.93800E+03 1.93900E+03 1.94000E+03 1.94100E+03 1.94200E+03 1.94300E+03 1.94400E+03 1.94500E+03 &
     1.94600E+03 1.94700E+03 1.94800E+03 1.94900E+03 1.95000E+03 1.95100E+03 1.95200E+03 1.95300E+03 &
     1.95400E+03 1.95500E+03 1.95600E+03 1.95700E+03 1.95800E+03 1.95900E+03 1.96000E+03 1.96100E+03 &
     1.96200E+03 1.96300E+03 1.96400E+03 1.96500E+03 1.96600E+03 1.96700E+03 1.96800E+03 1.96900E+03 &
     1.97000E+03 1.97100E+03 1.97200E+03 1.97300E+03 1.97400E+03 1.97500E+03 1.97600E+03 1.97700E+03 &
     1.97800E+03 1.97900E+03 1.98000E+03 1.98100E+03 1.98200E+03 1.98300E+03 1.98400E+03 1.98500E+03 &
     1.98600E+03 1.98700E+03 1.98800E+03 1.98900E+03 1.99000E+03 1.99100E+03 1.99200E+03 1.99300E+03 &
     1.99400E+03 1.99500E+03 1.99600E+03 1.99700E+03 1.99800E+03 1.99900E+03 2.00000E+03 2.00100E+03 &
     2.00200E+03 2.00300E+03 2.00400E+03 2.00500E+03 2.00600E+03 2.00700E+03 2.00800E+03 2.00900E+03 &
     2.01000E+03 2.01100E+03 2.01200E+03 2.01300E+03 2.01400E+03 2.01500E+03 2.01600E+03 2.01700E+03 &
     2.01800E+03 2.01900E+03 2.02000E+03 2.02100E+03 2.02200E+03 2.02300E+03 2.02400E+03 2.02500E+03 &
     2.02600E+03 2.02700E+03 2.02800E+03 2.02900E+03 2.03000E+03 2.03100E+03 2.03200E+03 2.03300E+03 &
     2.03400E+03 2.03500E+03 2.03600E+03 2.03700E+03 2.03800E+03 2.03900E+03 2.04000E+03 2.04100E+03 &
     2.04200E+03 2.04300E+03 2.04400E+03 2.04500E+03 2.04600E+03 2.04700E+03 2.04800E+03 2.04900E+03 &
     2.05000E+03 2.05100E+03 2.05200E+03 2.05300E+03 2.05400E+03 2.05500E+03 2.05600E+03 2.05700E+03 &
     2.05800E+03 2.05900E+03 2.06000E+03 2.06100E+03 2.06200E+03 2.06300E+03 2.06400E+03 2.06500E+03 &
     2.06600E+03 2.06700E+03 2.06800E+03 2.06900E+03 2.07000E+03 2.07100E+03 2.07200E+03 2.07300E+03 &
     2.07400E+03 2.07500E+03 2.07600E+03 2.07700E+03 2.07800E+03 2.07900E+03 2.08000E+03 2.08100E+03 &
     2.08200E+03 2.08300E+03 2.08400E+03 2.08500E+03 2.08600E+03 2.08700E+03 2.08800E+03 2.08900E+03 &
     2.09000E+03 2.09100E+03 2.09200E+03 2.09300E+03 2.09400E+03 2.09500E+03 2.09600E+03 2.09700E+03 &
     2.09800E+03 2.09900E+03 2.10000E+03 2.10100E+03 2.10200E+03 2.10300E+03 2.10400E+03 2.10500E+03 &
     2.10600E+03 2.10700E+03 2.10800E+03 2.10900E+03 2.11000E+03 2.11100E+03 2.11200E+03 2.11300E+03 &
     2.11400E+03 2.11500E+03 2.11600E+03 2.11700E+03 2.11800E+03 2.11900E+03 2.12000E+03 2.12100E+03 &
     2.12200E+03 2.12300E+03 2.12400E+03 2.12500E+03 2.12600E+03 2.12700E+03 2.12800E+03 2.12900E+03 &
     2.13000E+03 2.13100E+03 2.13200E+03 2.13300E+03 2.13400E+03 2.13500E+03 2.13600E+03 2.13700E+03 &
     2.13800E+03 2.13900E+03 2.14000E+03 2.14100E+03 2.14200E+03 2.14300E+03 2.14400E+03 2.14500E+03 &
     2.14600E+03 2.14700E+03 2.14800E+03 2.14900E+03 2.15000E+03 2.15100E+03 2.15200E+03 2.15300E+03 &
     2.15400E+03 2.15500E+03 2.15600E+03 2.15700E+03 2.15800E+03 2.15900E+03 2.16000E+03 2.16100E+03 &
     2.16200E+03 2.16300E+03 2.16400E+03 2.16500E+03 2.16600E+03 2.16700E+03 2.16800E+03 2.16900E+03 &
     2.17000E+03 2.17100E+03 2.17200E+03 2.17300E+03 2.17400E+03 2.17500E+03 2.17600E+03 2.17700E+03 &
     2.17800E+03 2.17900E+03 2.18000E+03 2.18100E+03 2.18200E+03 2.18300E+03 2.18400E+03 2.18500E+03 &
     2.18600E+03 2.18700E+03 2.18800E+03 2.18900E+03 2.19000E+03 2.19100E+03 2.19200E+03 2.19300E+03 &
     2.19400E+03 2.19500E+03 2.19600E+03 2.19700E+03 2.19800E+03 2.19900E+03 2.20000E+03 2.20100E+03 &
     2.20200E+03 2.20300E+03 2.20400E+03 2.20500E+03 2.20600E+03 2.20700E+03 2.20800E+03 2.20900E+03 &
     2.21000E+03 2.21100E+03 2.21200E+03 2.21300E+03 2.21400E+03 2.21500E+03 2.21600E+03 2.21700E+03 &
     2.21800E+03 2.21900E+03 2.22000E+03 2.22100E+03 2.22200E+03 2.22300E+03 2.22400E+03 2.22500E+03 &
     2.22600E+03 2.22700E+03 2.22800E+03 2.22900E+03 2.23000E+03 2.23100E+03 2.23200E+03 2.23300E+03 &
     2.23400E+03 2.23500E+03 2.23600E+03 2.23700E+03 2.23800E+03 2.23900E+03 2.24000E+03 2.24100E+03 &
     2.24200E+03 2.24300E+03 2.24400E+03 2.24500E+03 2.24600E+03 2.24700E+03 2.24800E+03 2.24900E+03 &
     2.25000E+03 2.25100E+03 2.25200E+03 2.25300E+03 2.25400E+03 2.25500E+03 2.25600E+03 2.25700E+03 &
     2.25800E+03 2.25900E+03 2.26000E+03 2.26100E+03 2.26200E+03 2.26300E+03 2.26400E+03 2.26500E+03 &
     2.26600E+03 2.26700E+03 2.26800E+03 2.26900E+03 2.27000E+03 2.27100E+03 2.27200E+03 2.27300E+03 &
     2.27400E+03 2.27500E+03 2.27600E+03 2.27700E+03 2.27800E+03 2.27900E+03 2.28000E+03 2.28100E+03 &
     2.28200E+03 2.28300E+03 2.28400E+03 2.28500E+03 2.28600E+03 2.28700E+03 2.28800E+03 2.28900E+03 &
     2.29000E+03 2.29100E+03 2.29200E+03 2.29300E+03 2.29400E+03 2.29500E+03 2.29600E+03 2.29700E+03 &
     2.29800E+03 2.29900E+03 2.30000E+03 2.30100E+03 2.30200E+03 2.30300E+03 2.30400E+03 2.30500E+03 &
     2.30600E+03 2.30700E+03 2.30800E+03 2.30900E+03 2.31000E+03 2.31100E+03 2.31200E+03 2.31300E+03 &
     2.31400E+03 2.31500E+03 2.31600E+03 2.31700E+03 2.31800E+03 2.31900E+03 2.32000E+03 2.32100E+03 &
     2.32200E+03 2.32300E+03 2.32400E+03 2.32500E+03 2.32600E+03 2.32700E+03 2.32800E+03 2.32900E+03 &
     2.33000E+03 2.33100E+03 2.33200E+03 2.33300E+03 2.33400E+03 2.33500E+03 2.33600E+03 2.33700E+03 &
     2.33800E+03 2.33900E+03 2.34000E+03 2.34100E+03 2.34200E+03 2.34300E+03 2.34400E+03 2.34500E+03 &
     2.34600E+03 2.34700E+03 2.34800E+03 2.34900E+03 2.35000E+03 2.35100E+03 2.35200E+03 2.35300E+03 &
     2.35400E+03 2.35500E+03 2.35600E+03 2.35700E+03 2.35800E+03 2.35900E+03 2.36000E+03 2.36100E+03 &
     2.36200E+03 2.36300E+03 2.36400E+03 2.36500E+03 2.36600E+03 2.36700E+03 2.36800E+03 2.36900E+03 &
     2.37000E+03 2.37100E+03 2.37200E+03 2.37300E+03 2.37400E+03 2.37500E+03 2.37600E+03 2.37700E+03 &
     2.37800E+03 2.37900E+03 2.38000E+03 2.38100E+03 2.38200E+03 2.38300E+03 2.38400E+03 2.38500E+03 &
     2.38600E+03 2.38700E+03 2.38800E+03 2.38900E+03 2.39000E+03 2.39100E+03 2.39200E+03 2.39300E+03 &
     2.39400E+03 2.39500E+03 2.39600E+03 2.39700E+03 2.39800E+03 2.39900E+03 2.40000E+03 2.40100E+03 &
     2.40200E+03 2.40300E+03 2.40400E+03 2.40500E+03 2.40600E+03 2.40700E+03 2.40800E+03 2.40900E+03 &
     2.41000E+03 2.41100E+03 2.41200E+03 2.41300E+03 2.41400E+03 2.41500E+03 2.41600E+03 2.41700E+03 &
     2.41800E+03 2.41900E+03 2.42000E+03 2.42100E+03 2.42200E+03 2.42300E+03 2.42400E+03 2.42500E+03 &
     2.42600E+03 2.42700E+03 2.42800E+03 2.42900E+03 2.43000E+03 2.43100E+03 2.43200E+03 2.43300E+03 &
     2.43400E+03 2.43500E+03 2.43600E+03 2.43700E+03 2.43800E+03 2.43900E+03 2.44000E+03 2.44100E+03 &
     2.44200E+03 2.44300E+03 2.44400E+03 2.44500E+03 2.44600E+03 2.44700E+03 2.44800E+03 2.44900E+03 &
     2.45000E+03 2.45100E+03 2.45200E+03 2.45300E+03 2.45400E+03 2.45500E+03 2.45600E+03 2.45700E+03 &
     2.45800E+03 2.45900E+03 2.46000E+03 2.46100E+03 2.46200E+03 2.46300E+03 2.46400E+03 2.46500E+03 &
     2.46600E+03 2.46700E+03 2.46800E+03 2.46900E+03 2.47000E+03 2.47100E+03 2.47200E+03 2.47300E+03 &
     2.47400E+03 2.47500E+03 2.47600E+03 2.47700E+03 2.47800E+03 2.47900E+03 2.48000E+03 2.48100E+03 &
     2.48200E+03 2.48300E+03 2.48400E+03 2.48500E+03 2.48600E+03 2.48700E+03 2.48800E+03 2.48900E+03 &
     2.49000E+03 2.49100E+03 2.49200E+03 2.49300E+03 2.49400E+03 2.49500E+03 2.49600E+03 2.49700E+03 &
     2.49800E+03 2.49900E+03 2.50000E+03 2.50100E+03 2.50200E+03 2.50300E+03 2.50400E+03 2.50500E+03 &
     2.50600E+03 2.50700E+03 2.50800E+03 2.50900E+03 2.51000E+03 2.51100E+03 2.51200E+03 2.51300E+03 &
     2.51400E+03 2.51500E+03 2.51600E+03 2.51700E+03 2.51800E+03 2.51900E+03 2.52000E+03 2.52100E+03 &
     2.52200E+03 2.52300E+03 2.52400E+03 2.52500E+03 2.52600E+03 2.52700E+03 2.52800E+03 2.52900E+03 &
     2.53000E+03 2.53100E+03 2.53200E+03 2.53300E+03 2.53400E+03 2.53500E+03 2.53600E+03 2.53700E+03 &
     2.53800E+03 2.53900E+03 2.54000E+03 2.54100E+03 2.54200E+03 2.54300E+03 2.54400E+03 2.54500E+03 &
     2.54600E+03 2.54700E+03 2.54800E+03 2.54900E+03 2.55000E+03 2.55100E+03 2.55200E+03 2.55300E+03 &
     2.55400E+03 2.55500E+03 2.55600E+03 2.55700E+03 2.55800E+03 2.55900E+03 2.56000E+03 2.56100E+03 &
     2.56200E+03 2.56300E+03 2.56400E+03 2.56500E+03 2.56600E+03 2.56700E+03 2.56800E+03 2.56900E+03 &
     2.57000E+03 2.57100E+03 2.57200E+03 2.57300E+03 2.57400E+03 2.57500E+03 2.57600E+03 2.57700E+03 &
     2.57800E+03 2.57900E+03 2.58000E+03 2.58100E+03 2.58200E+03 2.58300E+03 2.58400E+03 2.58500E+03 &
     2.58600E+03 2.58700E+03 2.58800E+03 2.58900E+03 2.59000E+03 2.59100E+03 2.59200E+03 2.59300E+03 &
     2.59400E+03 2.59500E+03 2.59600E+03 2.59700E+03 2.59800E+03 2.59900E+03 2.60000E+03 2.60100E+03 &
     2.60200E+03 2.60300E+03 2.60400E+03 2.60500E+03 2.60600E+03 2.60700E+03 2.60800E+03 2.60900E+03 &
     2.61000E+03 2.61100E+03 2.61200E+03 2.61300E+03 2.61400E+03 2.61500E+03 2.61600E+03 2.61700E+03 &
     2.61800E+03 2.61900E+03 2.62000E+03 2.62100E+03 2.62200E+03 2.62300E+03 2.62400E+03 2.62500E+03 &
     2.62600E+03 2.62700E+03 2.62800E+03 2.62900E+03 2.63000E+03 2.63100E+03 2.63200E+03 2.63300E+03 &
     2.63400E+03 2.63500E+03 2.63600E+03 2.63700E+03 2.63800E+03 2.63900E+03 2.64000E+03 2.64100E+03 &
     2.64200E+03 2.64300E+03 2.64400E+03 2.64500E+03 2.64600E+03 2.64700E+03 2.64800E+03 2.64900E+03 &
     2.65000E+03 2.65100E+03 2.65200E+03 2.65300E+03 2.65400E+03 2.65500E+03 2.65600E+03 2.65700E+03 &
     2.65800E+03 2.65900E+03 2.66000E+03 2.66100E+03 2.66200E+03 2.66300E+03 2.66400E+03 2.66500E+03 &
     2.66600E+03 2.66700E+03 2.66800E+03 2.66900E+03 2.67000E+03 2.67100E+03 2.67200E+03 2.67300E+03 &
     2.67400E+03 2.67500E+03 2.67600E+03 2.67700E+03 2.67800E+03 2.67900E+03 2.68000E+03 2.68100E+03 &
     2.68200E+03 2.68300E+03 2.68400E+03 2.68500E+03 2.68600E+03 2.68700E+03 2.68800E+03 2.68900E+03 &
     2.69000E+03 2.69100E+03 2.69200E+03 2.69300E+03 2.69400E+03 2.69500E+03 2.69600E+03 2.69700E+03 &
     2.69800E+03 2.69900E+03 2.70000E+03 2.70100E+03 2.70200E+03 2.70300E+03 2.70400E+03 2.70500E+03 &
     2.70600E+03 2.70700E+03 2.70800E+03 2.70900E+03 2.71000E+03 2.71100E+03 2.71200E+03 2.71300E+03 &
     2.71400E+03 2.71500E+03 2.71600E+03 2.71700E+03 2.71800E+03 2.71900E+03 2.72000E+03 2.72100E+03 &
     2.72200E+03 2.72300E+03 2.72400E+03 2.72500E+03 2.72600E+03 2.72700E+03 2.72800E+03 2.72900E+03 &
     2.73000E+03 2.73100E+03 2.73200E+03 2.73300E+03 2.73400E+03 2.73500E+03 2.73600E+03 2.73700E+03 &
     2.73800E+03 2.73900E+03 2.74000E+03 2.74100E+03 2.74200E+03 2.74300E+03 2.74400E+03 2.74500E+03 &
     2.74600E+03 2.74700E+03 2.74800E+03 2.74900E+03 2.75000E+03 2.75100E+03 2.75200E+03 2.75300E+03 &
     2.75400E+03 2.75500E+03 2.75600E+03 2.75700E+03 2.75800E+03 2.75900E+03 2.76000E+03 2.76100E+03 &
     2.76200E+03 2.76300E+03 2.76400E+03 2.76500E+03 2.76600E+03 2.76700E+03 2.76800E+03 2.76900E+03 &
     2.77000E+03 2.77100E+03 2.77200E+03 2.77300E+03 2.77400E+03 2.77500E+03 2.77600E+03 2.77700E+03 &
     2.77800E+03 2.77900E+03 2.78000E+03 2.78100E+03 2.78200E+03 2.78300E+03 2.78400E+03 2.78500E+03 &
     2.78600E+03 2.78700E+03 2.78800E+03 2.78900E+03 2.79000E+03 2.79100E+03 2.79200E+03 2.79300E+03 &
     2.79400E+03 2.79500E+03 2.79600E+03 2.79700E+03 2.79800E+03 2.79900E+03 2.80000E+03 2.80100E+03 &
     2.80200E+03 2.80300E+03 2.80400E+03 2.80500E+03 2.80600E+03 2.80700E+03 2.80800E+03 2.80900E+03 &
     2.81000E+03 2.81100E+03 2.81200E+03 2.81300E+03 2.81400E+03 2.81500E+03 2.81600E+03 2.81700E+03 &
     2.81800E+03 2.81900E+03 2.82000E+03 2.82100E+03 2.82200E+03 2.82300E+03 2.82400E+03 2.82500E+03 &
     2.82600E+03 2.82700E+03 2.82800E+03 2.82900E+03 2.83000E+03 2.83100E+03 2.83200E+03 2.83300E+03 &
     2.83400E+03 2.83500E+03 2.83600E+03 2.83700E+03 2.83800E+03 2.83900E+03 2.84000E+03 2.84100E+03 &
     2.84200E+03 2.84300E+03 2.84400E+03 2.84500E+03 2.84600E+03 2.84700E+03 2.84800E+03 2.84900E+03 &
     2.85000E+03 2.85100E+03 2.85200E+03 2.85300E+03 2.85400E+03 2.85500E+03 2.85600E+03 2.85700E+03 &
     2.85800E+03 2.85900E+03 2.86000E+03 2.86100E+03 2.86200E+03 2.86300E+03 2.86400E+03 2.86500E+03 &
     2.86600E+03 2.86700E+03 2.86800E+03 2.86900E+03 2.87000E+03 2.87100E+03 2.87200E+03 2.87300E+03 &
     2.87400E+03 2.87500E+03 2.87600E+03 2.87700E+03 2.87800E+03 2.87900E+03 2.88000E+03 2.88100E+03 &
     2.88200E+03 2.88300E+03 2.88400E+03 2.88500E+03 2.88600E+03 2.88700E+03 2.88800E+03 2.88900E+03 &
     2.89000E+03 2.89100E+03 2.89200E+03 2.89300E+03 2.89400E+03 2.89500E+03 2.89600E+03 2.89700E+03 &
     2.89800E+03 2.89900E+03 2.90000E+03 2.90100E+03 2.90200E+03 2.90300E+03 2.90400E+03 2.90500E+03 &
     2.90600E+03 2.90700E+03 2.90800E+03 2.90900E+03 2.91000E+03 2.91100E+03 2.91200E+03 2.91300E+03 &
     2.91400E+03 2.91500E+03 2.91600E+03 2.91700E+03 2.91800E+03 2.91900E+03 2.92000E+03 2.92100E+03 &
     2.92200E+03 2.92300E+03 2.92400E+03 2.92500E+03 2.92600E+03 2.92700E+03 2.92800E+03 2.92900E+03 &
     2.93000E+03 2.93100E+03 2.93200E+03 2.93300E+03 2.93400E+03 2.93500E+03 2.93600E+03 2.93700E+03 &
     2.93800E+03 2.93900E+03 2.94000E+03 2.94100E+03 2.94200E+03 2.94300E+03 2.94400E+03 2.94500E+03 &
     2.94600E+03 2.94700E+03 2.94800E+03 2.94900E+03 2.95000E+03 2.95100E+03 2.95200E+03 2.95300E+03 &
     2.95400E+03 2.95500E+03 2.95600E+03 2.95700E+03 2.95800E+03 2.95900E+03 2.96000E+03 2.96100E+03 &
     2.96200E+03 2.96300E+03 2.96400E+03 2.96500E+03 2.96600E+03 2.96700E+03 2.96800E+03 2.96900E+03 &
     2.97000E+03 2.97100E+03 2.97200E+03 2.97300E+03 2.97400E+03 2.97500E+03 2.97600E+03 2.97700E+03 &
     2.97800E+03 2.97900E+03 2.98000E+03 2.98100E+03 2.98200E+03 2.98300E+03 2.98400E+03 2.98500E+03 &
     2.98600E+03 2.98700E+03 2.98800E+03 2.98900E+03 2.99000E+03 2.99100E+03 2.99200E+03 2.99300E+03 &
     2.99400E+03 2.99500E+03 2.99600E+03 2.99700E+03 2.99800E+03 2.99900E+03 3.00000E+03 3.00100E+03 &
     3.00200E+03 3.00300E+03 3.00400E+03 3.00500E+03 3.00600E+03 3.00700E+03 3.00800E+03 3.00900E+03 &
     3.01000E+03 3.01100E+03 3.01200E+03 3.01300E+03 3.01400E+03 3.01500E+03 3.01600E+03 3.01700E+03 &
     3.01800E+03 3.01900E+03 3.02000E+03 3.02100E+03 3.02200E+03 3.02300E+03 3.02400E+03 3.02500E+03 &
     3.02600E+03 3.02700E+03 3.02800E+03 3.02900E+03 3.03000E+03 3.03100E+03 3.03200E+03 3.03300E+03 &
     3.03400E+03 3.03500E+03 3.03600E+03 3.03700E+03 3.03800E+03 3.03900E+03 3.04000E+03 3.04100E+03 &
     3.04200E+03 3.04300E+03 3.04400E+03 3.04500E+03 3.04600E+03 3.04700E+03 3.04800E+03 3.04900E+03 &
     3.05000E+03 3.05100E+03 3.05200E+03 3.05300E+03 3.05400E+03 3.05500E+03 3.05600E+03 3.05700E+03 &
     3.05800E+03 3.05900E+03 3.06000E+03 3.06100E+03 3.06200E+03 3.06300E+03 3.06400E+03 3.06500E+03 &
     3.06600E+03 3.06700E+03 3.06800E+03 3.06900E+03 3.07000E+03 3.07100E+03 3.07200E+03 3.07300E+03 &
     3.07400E+03 3.07500E+03 3.07600E+03 3.07700E+03 3.07800E+03 3.07900E+03 3.08000E+03 3.08100E+03 &
     3.08200E+03 3.08300E+03 3.08400E+03 3.08500E+03 3.08600E+03 3.08700E+03 3.08800E+03 3.08900E+03 &
     3.09000E+03 3.09100E+03 3.09200E+03 3.09300E+03 3.09400E+03 3.09500E+03 3.09600E+03 3.09700E+03 &
     3.09800E+03 3.09900E+03 3.10000E+03 3.10100E+03 3.10200E+03 3.10300E+03 3.10400E+03 3.10500E+03 &
     3.10600E+03 3.10700E+03 3.10800E+03 3.10900E+03 3.11000E+03 3.11100E+03 3.11200E+03 3.11300E+03 &
     3.11400E+03 3.11500E+03 3.11600E+03 3.11700E+03 3.11800E+03 3.11900E+03 3.12000E+03 3.12100E+03 &
     3.12200E+03 3.12300E+03 3.12400E+03 3.12500E+03 3.12600E+03 3.12700E+03 3.12800E+03 3.12900E+03 &
     3.13000E+03 3.13100E+03 3.13200E+03 3.13300E+03 3.13400E+03 3.13500E+03 3.13600E+03 3.13700E+03 &
     3.13800E+03 3.13900E+03 3.14000E+03 3.14100E+03 3.14200E+03 3.14300E+03 3.14400E+03 3.14500E+03 &
     3.14600E+03 3.14700E+03 3.14800E+03 3.14900E+03 3.15000E+03 3.15100E+03 3.15200E+03 3.15300E+03 &
     3.15400E+03 3.15500E+03 3.15600E+03 3.15700E+03 3.15800E+03 3.15900E+03 3.16000E+03 3.16100E+03 &
     3.16200E+03 3.16300E+03 3.16400E+03 3.16500E+03 3.16600E+03 3.16700E+03 3.16800E+03 3.16900E+03 &
     3.17000E+03 3.17100E+03 3.17200E+03 3.17300E+03 3.17400E+03 3.17500E+03 3.17600E+03 3.17700E+03 &
     3.17800E+03 3.17900E+03 3.18000E+03 3.18100E+03 3.18200E+03 3.18300E+03 3.18400E+03 3.18500E+03 &
     3.18600E+03 3.18700E+03 3.18800E+03 3.18900E+03 3.19000E+03 3.19100E+03 3.19200E+03 3.19300E+03 &
     3.19400E+03 3.19500E+03 3.19600E+03 3.19700E+03 3.19800E+03 3.19900E+03 3.20000E+03 3.20100E+03 &
     3.20200E+03 3.20300E+03 3.20400E+03 3.20500E+03 3.20600E+03 3.20700E+03 3.20800E+03 3.20900E+03 &
     3.21000E+03 3.21100E+03 3.21200E+03 3.21300E+03 3.21400E+03 3.21500E+03 3.21600E+03 3.21700E+03 &
     3.21800E+03 3.21900E+03 3.22000E+03 3.22100E+03 3.22200E+03 3.22300E+03 3.22400E+03 3.22500E+03 &
     3.22600E+03 3.22700E+03 3.22800E+03 3.22900E+03 3.23000E+03 3.23100E+03 3.23200E+03 3.23300E+03 &
     3.23400E+03 3.23500E+03 3.23600E+03 3.23700E+03 3.23800E+03 3.23900E+03 3.24000E+03 3.24100E+03 &
     3.24200E+03 3.24300E+03 3.24400E+03 3.24500E+03 3.24600E+03 3.24700E+03 3.24800E+03 3.24900E+03 &
     3.25000E+03 3.25100E+03 3.25200E+03 3.25300E+03 3.25400E+03 3.25500E+03 3.25600E+03 3.25700E+03 &
     3.25800E+03 3.25900E+03 3.26000E+03 3.26100E+03 3.26200E+03 3.26300E+03 3.26400E+03 3.26500E+03 &
     3.26600E+03 3.26700E+03 3.26800E+03 3.26900E+03 3.27000E+03 3.27100E+03 3.27200E+03 3.27300E+03 &
     3.27400E+03 3.27500E+03 3.27600E+03 3.27700E+03 3.27800E+03 3.27900E+03 3.28000E+03 3.28100E+03 &
     3.28200E+03 3.28300E+03 3.28400E+03 3.28500E+03 3.28600E+03 3.28700E+03 3.28800E+03 3.28900E+03 &
     3.29000E+03 3.29100E+03 3.29200E+03 3.29300E+03 3.29400E+03 3.29500E+03 3.29600E+03 3.29700E+03 &
     3.29800E+03 3.29900E+03 3.30000E+03 3.30100E+03 3.30200E+03 3.30300E+03 3.30400E+03 3.30500E+03 &
     3.30600E+03 3.30700E+03 3.30800E+03 3.30900E+03 3.31000E+03 3.31100E+03 3.31200E+03 3.31300E+03 &
     3.31400E+03 3.31500E+03 3.31600E+03 3.31700E+03 3.31800E+03 3.31900E+03 3.32000E+03 3.32100E+03 &
     3.32200E+03 3.32300E+03 3.32400E+03 3.32500E+03 3.32600E+03 3.32700E+03 3.32800E+03 3.32900E+03 &
     3.33000E+03 3.33100E+03 3.33200E+03 3.33300E+03 3.33400E+03 3.33500E+03 3.33600E+03 3.33700E+03 &
     3.33800E+03 3.33900E+03 3.34000E+03 3.34100E+03 3.34200E+03 3.34300E+03 3.34400E+03 3.34500E+03 &
     3.34600E+03 3.34700E+03 3.34800E+03 3.34900E+03 3.35000E+03 3.35100E+03 3.35200E+03 3.35300E+03 &
     3.35400E+03 3.35500E+03 3.35600E+03 3.35700E+03 3.35800E+03 3.35900E+03 3.36000E+03 3.36100E+03 &
     3.36200E+03 3.36300E+03 3.36400E+03 3.36500E+03 3.36600E+03 3.36700E+03 3.36800E+03 3.36900E+03 &
     3.37000E+03 3.37100E+03 3.37200E+03 3.37300E+03 3.37400E+03 3.37500E+03 3.37600E+03 3.37700E+03 &
     3.37800E+03 3.37900E+03 3.38000E+03 3.38100E+03 3.38200E+03 3.38300E+03 3.38400E+03 3.38500E+03 &
     3.38600E+03 3.38700E+03 3.38800E+03 3.38900E+03 3.39000E+03 3.39100E+03 3.39200E+03 3.39300E+03 &
     3.39400E+03 3.39500E+03 3.39600E+03 3.39700E+03 3.39800E+03 3.39900E+03 3.40000E+03 3.40100E+03 &
     3.40200E+03 3.40300E+03 3.40400E+03 3.40500E+03 3.40600E+03 3.40700E+03 3.40800E+03 3.40900E+03 &
     3.41000E+03 3.41100E+03 3.41200E+03 3.41300E+03 3.41400E+03 3.41500E+03 3.41600E+03 3.41700E+03 &
     3.41800E+03 3.41900E+03 3.42000E+03 3.42100E+03 3.42200E+03 3.42300E+03 3.42400E+03 3.42500E+03 &
     3.42600E+03 3.42700E+03 3.42800E+03 3.42900E+03 3.43000E+03 3.43100E+03 3.43200E+03 3.43300E+03 &
     3.43400E+03 3.43500E+03 3.43600E+03 3.43700E+03 3.43800E+03 3.43900E+03 3.44000E+03 3.44100E+03 &
     3.44200E+03 3.44300E+03 3.44400E+03 3.44500E+03 3.44600E+03 3.44700E+03 3.44800E+03 3.44900E+03 &
     3.45000E+03 3.45100E+03 3.45200E+03 3.45300E+03 3.45400E+03 3.45500E+03 3.45600E+03 3.45700E+03 &
     3.45800E+03 3.45900E+03 3.46000E+03 3.46100E+03 3.46200E+03 3.46300E+03 3.46400E+03 3.46500E+03 &
     3.46600E+03 3.46700E+03 3.46800E+03 3.46900E+03 3.47000E+03 3.47100E+03 3.47200E+03 3.47300E+03 &
     3.47400E+03 3.47500E+03 3.47600E+03 3.47700E+03 3.47800E+03 3.47900E+03 3.48000E+03 3.48100E+03 &
     3.48200E+03 3.48300E+03 3.48400E+03 3.48500E+03 3.48600E+03 3.48700E+03 3.48800E+03 3.48900E+03 &
     3.49000E+03 3.49100E+03 3.49200E+03 3.49300E+03 3.49400E+03 3.49500E+03 3.49600E+03 3.49700E+03 &
     3.49800E+03 3.49900E+03 3.50000E+03 3.50100E+03 3.50200E+03 3.50300E+03 3.50400E+03 3.50500E+03 &
     3.50600E+03 3.50700E+03 3.50800E+03 3.50900E+03 3.51000E+03 3.51100E+03 3.51200E+03 3.51300E+03 &
     3.51400E+03 3.51500E+03 3.51600E+03 3.51700E+03 3.51800E+03 3.51900E+03 3.52000E+03 3.52100E+03 &
     3.52200E+03 3.52300E+03 3.52400E+03 3.52500E+03 3.52600E+03 3.52700E+03 3.52800E+03 3.52900E+03 &
     3.53000E+03 3.53100E+03 3.53200E+03 3.53300E+03 3.53400E+03 3.53500E+03 3.53600E+03 3.53700E+03 &
     3.53800E+03 3.53900E+03 3.54000E+03 3.54100E+03 3.54200E+03 3.54300E+03 3.54400E+03 3.54500E+03 &
     3.54600E+03 3.54700E+03 3.54800E+03 3.54900E+03 3.55000E+03 3.55100E+03 3.55200E+03 3.55300E+03 &
     3.55400E+03 3.55500E+03 3.55600E+03 3.55700E+03 3.55800E+03 3.55900E+03 3.56000E+03 3.56100E+03 &
     3.56200E+03 3.56300E+03 3.56400E+03 3.56500E+03 3.56600E+03 3.56700E+03 3.56800E+03 3.56900E+03 &
     3.57000E+03 3.57100E+03 3.57200E+03 3.57300E+03 3.57400E+03 3.57500E+03 3.57600E+03 3.57700E+03 &
     3.57800E+03 3.57900E+03 3.58000E+03 3.58100E+03 3.58200E+03 3.58300E+03 3.58400E+03 3.58500E+03 &
     3.58600E+03 3.58700E+03 3.58800E+03 3.58900E+03 3.59000E+03 3.59100E+03 3.59200E+03 3.59300E+03 &
     3.59400E+03 3.59500E+03 3.59600E+03 3.59700E+03 3.59800E+03 3.59900E+03 3.60000E+03 3.60100E+03 &
     3.60200E+03 3.60300E+03 3.60400E+03 3.60500E+03 3.60600E+03 3.60700E+03 3.60800E+03 3.60900E+03 &
     3.61000E+03 3.61100E+03 3.61200E+03 3.61300E+03 3.61400E+03 3.61500E+03 3.61600E+03 3.61700E+03 &
     3.61800E+03 3.61900E+03 3.62000E+03 3.62100E+03 3.62200E+03 3.62300E+03 3.62400E+03 3.62500E+03 &
     3.62600E+03 3.62700E+03 3.62800E+03 3.62900E+03 3.63000E+03 3.63100E+03 3.63200E+03 3.63300E+03 &
     3.63400E+03 3.63500E+03 3.63600E+03 3.63700E+03 3.63800E+03 3.63900E+03 3.64000E+03 3.64100E+03 &
     3.64200E+03 3.64300E+03 3.64400E+03 3.64500E+03 3.64600E+03 3.64700E+03 3.64800E+03 3.64900E+03 &
     3.65000E+03 3.65100E+03 3.65200E+03 3.65300E+03 3.65400E+03 3.65500E+03 3.65600E+03 3.65700E+03 &
     3.65800E+03 3.65900E+03 3.66000E+03 3.66100E+03 3.66200E+03 3.66300E+03 3.66400E+03 3.66500E+03 &
     3.66600E+03 3.66700E+03 3.66800E+03 3.66900E+03 3.67000E+03 3.67100E+03 3.67200E+03 3.67300E+03 &
     3.67400E+03 3.67500E+03 3.67600E+03 3.67700E+03 3.67800E+03 3.67900E+03 3.68000E+03 3.68100E+03 &
     3.68200E+03 3.68300E+03 3.68400E+03 3.68500E+03 3.68600E+03 3.68700E+03 3.68800E+03 3.68900E+03 &
     3.69000E+03 3.69100E+03 3.69200E+03 3.69300E+03 3.69400E+03 3.69500E+03 3.69600E+03 3.69700E+03 &
     3.69800E+03 3.69900E+03 3.70000E+03 3.70100E+03 3.70200E+03 3.70300E+03 3.70400E+03 3.70500E+03 &
     3.70600E+03 3.70700E+03 3.70800E+03 3.70900E+03 3.71000E+03 3.71100E+03 3.71200E+03 3.71300E+03 &
     3.71400E+03 3.71500E+03 3.71600E+03 3.71700E+03 3.71800E+03 3.71900E+03 3.72000E+03 3.72100E+03 &
     3.72200E+03 3.72300E+03 3.72400E+03 3.72500E+03 3.72600E+03 3.72700E+03 3.72800E+03 3.72900E+03 &
     3.73000E+03 3.73100E+03 3.73200E+03 3.73300E+03 3.73400E+03 3.73500E+03 3.73600E+03 3.73700E+03 &
     3.73800E+03 3.73900E+03 3.74000E+03 3.74100E+03 3.74200E+03 3.74300E+03 3.74400E+03 3.74500E+03 &
     3.74600E+03 3.74700E+03 3.74800E+03 3.74900E+03 3.75000E+03 3.75100E+03 3.75200E+03 3.75300E+03 &
     3.75400E+03 3.75500E+03 3.75600E+03 3.75700E+03 3.75800E+03 3.75900E+03 3.76000E+03 3.76100E+03 &
     3.76200E+03 3.76300E+03 3.76400E+03 3.76500E+03 3.76600E+03 3.76700E+03 3.76800E+03 3.76900E+03 &
     3.77000E+03 3.77100E+03 3.77200E+03 3.77300E+03 3.77400E+03 3.77500E+03 3.77600E+03 3.77700E+03 &
     3.77800E+03 3.77900E+03 3.78000E+03 3.78100E+03 3.78200E+03 3.78300E+03 3.78400E+03 3.78500E+03 &
     3.78600E+03 3.78700E+03 3.78800E+03 3.78900E+03 3.79000E+03 3.79100E+03 3.79200E+03 3.79300E+03 &
     3.79400E+03 3.79500E+03 3.79600E+03 3.79700E+03 3.79800E+03 3.79900E+03 3.80000E+03 3.80100E+03 &
     3.80200E+03 3.80300E+03 3.80400E+03 3.80500E+03 3.80600E+03 3.80700E+03 3.80800E+03 3.80900E+03 &
     3.81000E+03 3.81100E+03 3.81200E+03 3.81300E+03 3.81400E+03 3.81500E+03 3.81600E+03 3.81700E+03 &
     3.81800E+03 3.81900E+03 3.82000E+03 3.82100E+03 3.82200E+03 3.82300E+03 3.82400E+03 3.82500E+03 &
     3.82600E+03 3.82700E+03 3.82800E+03 3.82900E+03 3.83000E+03 3.83100E+03 3.83200E+03 3.83300E+03 &
     3.83400E+03 3.83500E+03 3.83600E+03 3.83700E+03 3.83800E+03 3.83900E+03 3.84000E+03 3.84100E+03 &
     3.84200E+03 3.84300E+03 3.84400E+03 3.84500E+03 3.84600E+03 3.84700E+03 3.84800E+03 3.84900E+03 &
     3.85000E+03 3.85100E+03 3.85200E+03 3.85300E+03 3.85400E+03 3.85500E+03 3.85600E+03 3.85700E+03 &
     3.85800E+03 3.85900E+03 3.86000E+03 3.86100E+03 3.86200E+03 3.86300E+03 3.86400E+03 3.86500E+03 &
     3.86600E+03 3.86700E+03 3.86800E+03 3.86900E+03 3.87000E+03 3.87100E+03 3.87200E+03 3.87300E+03 &
     3.87400E+03 3.87500E+03 3.87600E+03 3.87700E+03 3.87800E+03 3.87900E+03 3.88000E+03 3.88100E+03 &
     3.88200E+03 3.88300E+03 3.88400E+03 3.88500E+03 3.88600E+03 3.88700E+03 3.88800E+03 3.88900E+03 &
     3.89000E+03 3.89100E+03 3.89200E+03 3.89300E+03 3.89400E+03 3.89500E+03 3.89600E+03 3.89700E+03 &
     3.89800E+03 3.89900E+03 3.90000E+03 3.90100E+03 3.90200E+03 3.90300E+03 3.90400E+03 3.90500E+03 &
     3.90600E+03 3.90700E+03 3.90800E+03 3.90900E+03 3.91000E+03 3.91100E+03 3.91200E+03 3.91300E+03 &
     3.91400E+03 3.91500E+03 3.91600E+03 3.91700E+03 3.91800E+03 3.91900E+03 3.92000E+03 3.92100E+03 &
     3.92200E+03 3.92300E+03 3.92400E+03 3.92500E+03 3.92600E+03 3.92700E+03 3.92800E+03 3.92900E+03 &
     3.93000E+03 3.93100E+03 3.93200E+03 3.93300E+03 3.93400E+03 3.93500E+03 3.93600E+03 3.93700E+03 &
     3.93800E+03 3.93900E+03 3.94000E+03 3.94100E+03 3.94200E+03 3.94300E+03 3.94400E+03 3.94500E+03 &
     3.94600E+03 3.94700E+03 3.94800E+03 3.94900E+03 3.95000E+03 3.95100E+03 3.95200E+03 3.95300E+03 &
     3.95400E+03 3.95500E+03 3.95600E+03 3.95700E+03 3.95800E+03 3.95900E+03 3.96000E+03 3.96100E+03 &
     3.96200E+03 3.96300E+03 3.96400E+03 3.96500E+03 3.96600E+03 3.96700E+03 3.96800E+03 3.96900E+03 &
     3.97000E+03 3.97100E+03 3.97200E+03 3.97300E+03 3.97400E+03 3.97500E+03 3.97600E+03 3.97700E+03 &
     3.97800E+03 3.97900E+03 3.98000E+03 3.98100E+03 3.98200E+03 3.98300E+03 3.98400E+03 3.98500E+03 &
     3.98600E+03 3.98700E+03 3.98800E+03 3.98900E+03 3.99000E+03 3.99100E+03 3.99200E+03 3.99300E+03 &
     3.99400E+03 3.99500E+03 3.99600E+03 3.99700E+03 3.99800E+03 3.99900E+03 4.00000E+03 4.00100E+03 &
     4.00200E+03 4.00300E+03 4.00400E+03 4.00500E+03 4.00600E+03 4.00700E+03 4.00800E+03 4.00900E+03 &
     4.01000E+03 4.01100E+03 4.01200E+03 4.01300E+03 4.01400E+03 4.01500E+03 4.01600E+03 4.01700E+03 &
     4.01800E+03 4.01900E+03 4.02000E+03 4.02100E+03 4.02200E+03 4.02300E+03 4.02400E+03 4.02500E+03 &
     4.02600E+03 4.02700E+03 4.02800E+03 4.02900E+03 4.03000E+03 4.03100E+03 4.03200E+03 4.03300E+03 &
     4.03400E+03 4.03500E+03 4.03600E+03 4.03700E+03 4.03800E+03 4.03900E+03 4.04000E+03 4.04100E+03 &
     4.04200E+03 4.04300E+03 4.04400E+03 4.04500E+03 4.04600E+03 4.04700E+03 4.04800E+03 4.04900E+03 &
     4.05000E+03 4.05100E+03 4.05200E+03 4.05300E+03 4.05400E+03 4.05500E+03 4.05600E+03 4.05700E+03 &
     4.05800E+03 4.05900E+03 4.06000E+03 4.06100E+03 4.06200E+03 4.06300E+03 4.06400E+03 4.06500E+03 &
     4.06600E+03 4.06700E+03 4.06800E+03 4.06900E+03 4.07000E+03 4.07100E+03 4.07200E+03 4.07300E+03 &
     4.07400E+03 4.07500E+03 4.07600E+03 4.07700E+03 4.07800E+03 4.07900E+03 4.08000E+03 4.08100E+03 &
     4.08200E+03 4.08300E+03 4.08400E+03 4.08500E+03 4.08600E+03 4.08700E+03 4.08800E+03 4.08900E+03 &
     4.09000E+03 4.09100E+03 4.09200E+03 4.09300E+03 4.09400E+03 4.09500E+03 4.09600E+03 4.09700E+03 &
     4.09800E+03 4.09900E+03 4.10000E+03 4.10100E+03 4.10200E+03 4.10300E+03 4.10400E+03 4.10500E+03 &
     4.10600E+03 4.10700E+03 4.10800E+03 4.10900E+03 4.11000E+03 4.11100E+03 4.11200E+03 4.11300E+03 &
     4.11400E+03 4.11500E+03 4.11600E+03 4.11700E+03 4.11800E+03 4.11900E+03 4.12000E+03 4.12100E+03 &
     4.12200E+03 4.12300E+03 4.12400E+03 4.12500E+03 4.12600E+03 4.12700E+03 4.12800E+03 4.12900E+03 &
     4.13000E+03 4.13100E+03 4.13200E+03 4.13300E+03 4.13400E+03 4.13500E+03 4.13600E+03 4.13700E+03 &
     4.13800E+03 4.13900E+03 4.14000E+03 4.14100E+03 4.14200E+03 4.14300E+03 4.14400E+03 4.14500E+03 &
     4.14600E+03 4.14700E+03 4.14800E+03 4.14900E+03 4.15000E+03 4.15100E+03 4.15200E+03 4.15300E+03 &
     4.15400E+03 4.15500E+03 4.15600E+03 4.15700E+03 4.15800E+03 4.15900E+03 4.16000E+03 4.16100E+03 &
     4.16200E+03 4.16300E+03 4.16400E+03 4.16500E+03 4.16600E+03 4.16700E+03 4.16800E+03 4.16900E+03 &
     4.17000E+03 4.17100E+03 4.17200E+03 4.17300E+03 4.17400E+03 4.17500E+03 4.17600E+03 4.17700E+03 &
     4.17800E+03 4.17900E+03 4.18000E+03 4.18100E+03 4.18200E+03 4.18300E+03 4.18400E+03 4.18500E+03 &
     4.18600E+03 4.18700E+03 4.18800E+03 4.18900E+03 4.19000E+03 4.19100E+03 4.19200E+03 4.19300E+03 &
     4.19400E+03 4.19500E+03 4.19600E+03 4.19700E+03 4.19800E+03 4.19900E+03 4.20000E+03 4.20100E+03 &
     4.20200E+03 4.20300E+03 4.20400E+03 4.20500E+03 4.20600E+03 4.20700E+03 4.20800E+03 4.20900E+03 &
     4.21000E+03 4.21100E+03 4.21200E+03 4.21300E+03 4.21400E+03 4.21500E+03 4.21600E+03 4.21700E+03 &
     4.21800E+03 4.21900E+03 4.22000E+03 4.22100E+03 4.22200E+03 4.22300E+03 4.22400E+03 4.22500E+03 &
     4.22600E+03 4.22700E+03 4.22800E+03 4.22900E+03 4.23000E+03 4.23100E+03 4.23200E+03 4.23300E+03 &
     4.23400E+03 4.23500E+03 4.23600E+03 4.23700E+03 4.23800E+03 4.23900E+03 4.24000E+03 4.24100E+03 &
     4.24200E+03 4.24300E+03 4.24400E+03 4.24500E+03 4.24600E+03 4.24700E+03 4.24800E+03 4.24900E+03 &
     4.25000E+03 4.25100E+03 4.25200E+03 4.25300E+03 4.25400E+03 4.25500E+03 4.25600E+03 4.25700E+03 &
     4.25800E+03 4.25900E+03 4.26000E+03 4.26100E+03 4.26200E+03 4.26300E+03 4.26400E+03 4.26500E+03 &
     4.26600E+03 4.26700E+03 4.26800E+03 4.26900E+03 4.27000E+03 4.27100E+03 4.27200E+03 4.27300E+03 &
     4.27400E+03 4.27500E+03 4.27600E+03 4.27700E+03 4.27800E+03 4.27900E+03 4.28000E+03 4.28100E+03 &
     4.28200E+03 4.28300E+03 4.28400E+03 4.28500E+03 4.28600E+03 4.28700E+03 4.28800E+03 4.28900E+03 &
     4.29000E+03 4.29100E+03 4.29200E+03 4.29300E+03 4.29400E+03 4.29500E+03 4.29600E+03 4.29700E+03 &
     4.29800E+03 4.29900E+03 4.30000E+03 4.30100E+03 4.30200E+03 4.30300E+03 4.30400E+03 4.30500E+03 &
     4.30600E+03 4.30700E+03 4.30800E+03 4.30900E+03 4.31000E+03 4.31100E+03 4.31200E+03 4.31300E+03 &
     4.31400E+03 4.31500E+03 4.31600E+03 4.31700E+03 4.31800E+03 4.31900E+03 4.32000E+03 4.32100E+03 &
     4.32200E+03 4.32300E+03 4.32400E+03 4.32500E+03 4.32600E+03 4.32700E+03 4.32800E+03 4.32900E+03 &
     4.33000E+03 4.33100E+03 4.33200E+03 4.33300E+03 4.33400E+03 4.33500E+03 4.33600E+03 4.33700E+03 &
     4.33800E+03 4.33900E+03 4.34000E+03 4.34100E+03 4.34200E+03 4.34300E+03 4.34400E+03 4.34500E+03 &
     4.34600E+03 4.34700E+03 4.34800E+03 4.34900E+03 4.35000E+03 4.35100E+03 4.35200E+03 4.35300E+03 &
     4.35400E+03 4.35500E+03 4.35600E+03 4.35700E+03 4.35800E+03 4.35900E+03 4.36000E+03 4.36100E+03 &
     4.36200E+03 4.36300E+03 4.36400E+03 4.36500E+03 4.36600E+03 4.36700E+03 4.36800E+03 4.36900E+03 &
     4.37000E+03 4.37100E+03 4.37200E+03 4.37300E+03 4.37400E+03 4.37500E+03 4.37600E+03 4.37700E+03 &
     4.37800E+03 4.37900E+03 4.38000E+03 4.38100E+03 4.38200E+03 4.38300E+03 4.38400E+03 4.38500E+03 &
     4.38600E+03 4.38700E+03 4.38800E+03 4.38900E+03 4.39000E+03 4.39100E+03 4.39200E+03 4.39300E+03 &
     4.39400E+03 4.39500E+03 4.39600E+03 4.39700E+03 4.39800E+03 4.39900E+03 4.40000E+03 4.40100E+03 &
     4.40200E+03 4.40300E+03 4.40400E+03 4.40500E+03 4.40600E+03 4.40700E+03 4.40800E+03 4.40900E+03 &
     4.41000E+03 4.41100E+03 4.41200E+03 4.41300E+03 4.41400E+03 4.41500E+03 4.41600E+03 4.41700E+03 &
     4.41800E+03 4.41900E+03 4.42000E+03 4.42100E+03 4.42200E+03 4.42300E+03 4.42400E+03 4.42500E+03 &
     4.42600E+03 4.42700E+03 4.42800E+03 4.42900E+03 4.43000E+03 4.43100E+03 4.43200E+03 4.43300E+03 &
     4.43400E+03 4.43500E+03 4.43600E+03 4.43700E+03 4.43800E+03 4.43900E+03 4.44000E+03 4.44100E+03 &
     4.44200E+03 4.44300E+03 4.44400E+03 4.44500E+03 4.44600E+03 4.44700E+03 4.44800E+03 4.44900E+03 &
     4.45000E+03 4.45100E+03 4.45200E+03 4.45300E+03 4.45400E+03 4.45500E+03 4.45600E+03 4.45700E+03 &
     4.45800E+03 4.45900E+03 4.46000E+03 4.46100E+03 4.46200E+03 4.46300E+03 4.46400E+03 4.46500E+03 &
     4.46600E+03 4.46700E+03 4.46800E+03 4.46900E+03 4.47000E+03 4.47100E+03 4.47200E+03 4.47300E+03 &
     4.47400E+03 4.47500E+03 4.47600E+03 4.47700E+03 4.47800E+03 4.47900E+03 4.48000E+03 4.48100E+03 &
     4.48200E+03 4.48300E+03 4.48400E+03 4.48500E+03 4.48600E+03 4.48700E+03 4.48800E+03 4.48900E+03 &
     4.49000E+03 4.49100E+03 4.49200E+03 4.49300E+03 4.49400E+03 4.49500E+03 4.49600E+03 4.49700E+03 &
     4.49800E+03 4.49900E+03 4.50000E+03 4.50100E+03 4.50200E+03 4.50300E+03 4.50400E+03 4.50500E+03 &
     4.50600E+03 4.50700E+03 4.50800E+03 4.50900E+03 4.51000E+03 4.51100E+03 4.51200E+03 4.51300E+03 &
     4.51400E+03 4.51500E+03 4.51600E+03 4.51700E+03 4.51800E+03 4.51900E+03 4.52000E+03 4.52100E+03 &
     4.52200E+03 4.52300E+03 4.52400E+03 4.52500E+03 4.52600E+03 4.52700E+03 4.52800E+03 4.52900E+03 &
     4.53000E+03 4.53100E+03 4.53200E+03 4.53300E+03 4.53400E+03 4.53500E+03 4.53600E+03 4.53700E+03 &
     4.53800E+03 4.53900E+03 4.54000E+03 4.54100E+03 4.54200E+03 4.54300E+03 4.54400E+03 4.54500E+03 &
     4.54600E+03 4.54700E+03 4.54800E+03 4.54900E+03 4.55000E+03 4.55100E+03 4.55200E+03 4.55300E+03 &
     4.55400E+03 4.55500E+03 4.55600E+03 4.55700E+03 4.55800E+03 4.55900E+03 4.56000E+03 4.56100E+03 &
     4.56200E+03 4.56300E+03 4.56400E+03 4.56500E+03 4.56600E+03 4.56700E+03 4.56800E+03 4.56900E+03 &
     4.57000E+03 4.57100E+03 4.57200E+03 4.57300E+03 4.57400E+03 4.57500E+03 4.57600E+03 4.57700E+03 &
     4.57800E+03 4.57900E+03 4.58000E+03 4.58100E+03 4.58200E+03 4.58300E+03 4.58400E+03 4.58500E+03 &
     4.58600E+03 4.58700E+03 4.58800E+03 4.58900E+03 4.59000E+03 4.59100E+03 4.59200E+03 4.59300E+03 &
     4.59400E+03 4.59500E+03 4.59600E+03 4.59700E+03 4.59800E+03 4.59900E+03 4.60000E+03 4.60100E+03 &
     4.60200E+03 4.60300E+03 4.60400E+03 4.60500E+03 4.60600E+03 4.60700E+03 4.60800E+03 4.60900E+03 &
     4.61000E+03 4.61100E+03 4.61200E+03 4.61300E+03 4.61400E+03 4.61500E+03 4.61600E+03 4.61700E+03 &
     4.61800E+03 4.61900E+03 4.62000E+03 4.62100E+03 4.62200E+03 4.62300E+03 4.62400E+03 4.62500E+03 &
     4.62600E+03 4.62700E+03 4.62800E+03 4.62900E+03 4.63000E+03 4.63100E+03 4.63200E+03 4.63300E+03 &
     4.63400E+03 4.63500E+03 4.63600E+03 4.63700E+03 4.63800E+03 4.63900E+03 4.64000E+03 4.64100E+03 &
     4.64200E+03 4.64300E+03 4.64400E+03 4.64500E+03 4.64600E+03 4.64700E+03 4.64800E+03 4.64900E+03 &
     4.65000E+03 4.65100E+03 4.65200E+03 4.65300E+03 4.65400E+03 4.65500E+03 4.65600E+03 4.65700E+03 &
     4.65800E+03 4.65900E+03 4.66000E+03 4.66100E+03 4.66200E+03 4.66300E+03 4.66400E+03 4.66500E+03 &
     4.66600E+03 4.66700E+03 4.66800E+03 4.66900E+03 4.67000E+03 4.67100E+03 4.67200E+03 4.67300E+03 &
     4.67400E+03 4.67500E+03 4.67600E+03 4.67700E+03 4.67800E+03 4.67900E+03 4.68000E+03 4.68100E+03 &
     4.68200E+03 4.68300E+03 4.68400E+03 4.68500E+03 4.68600E+03 4.68700E+03 4.68800E+03 4.68900E+03 &
     4.69000E+03 4.69100E+03 4.69200E+03 4.69300E+03 4.69400E+03 4.69500E+03 4.69600E+03 4.69700E+03 &
     4.69800E+03 4.69900E+03 4.70000E+03 4.70100E+03 4.70200E+03 4.70300E+03 4.70400E+03 4.70500E+03 &
     4.70600E+03 4.70700E+03 4.70800E+03 4.70900E+03 4.71000E+03 4.71100E+03 4.71200E+03 4.71300E+03 &
     4.71400E+03 4.71500E+03 4.71600E+03 4.71700E+03 4.71800E+03 4.71900E+03 4.72000E+03 4.72100E+03 &
     4.72200E+03 4.72300E+03 4.72400E+03 4.72500E+03 4.72600E+03 4.72700E+03 4.72800E+03 4.72900E+03 &
     4.73000E+03 4.73100E+03 4.73200E+03 4.73300E+03 4.73400E+03 4.73500E+03 4.73600E+03 4.73700E+03 &
     4.73800E+03 4.73900E+03 4.74000E+03 4.74100E+03 4.74200E+03 4.74300E+03 4.74400E+03 4.74500E+03 &
     4.74600E+03 4.74700E+03 4.74800E+03 4.74900E+03 4.75000E+03 4.75100E+03 4.75200E+03 4.75300E+03 &
     4.75400E+03 4.75500E+03 4.75600E+03 4.75700E+03 4.75800E+03 4.75900E+03 4.76000E+03 4.76100E+03 &
     4.76200E+03 4.76300E+03 4.76400E+03 4.76500E+03 4.76600E+03 4.76700E+03 4.76800E+03 4.76900E+03 &
     4.77000E+03 4.77100E+03 4.77200E+03 4.77300E+03 4.77400E+03 4.77500E+03 4.77600E+03 4.77700E+03 &
     4.77800E+03 4.77900E+03 4.78000E+03 4.78100E+03 4.78200E+03 4.78300E+03 4.78400E+03 4.78500E+03 &
     4.78600E+03 4.78700E+03 4.78800E+03 4.78900E+03 4.79000E+03 4.79100E+03 4.79200E+03 4.79300E+03 &
     4.79400E+03 4.79500E+03 4.79600E+03 4.79700E+03 4.79800E+03 4.79900E+03 4.80000E+03 4.80100E+03 &
     4.80200E+03 4.80300E+03 4.80400E+03 4.80500E+03 4.80600E+03 4.80700E+03 4.80800E+03 4.80900E+03 &
     4.81000E+03 4.81100E+03 4.81200E+03 4.81300E+03 4.81400E+03 4.81500E+03 4.81600E+03 4.81700E+03 &
     4.81800E+03 4.81900E+03 4.82000E+03 4.82100E+03 4.82200E+03 4.82300E+03 4.82400E+03 4.82500E+03 &
     4.82600E+03 4.82700E+03 4.82800E+03 4.82900E+03 4.83000E+03 4.83100E+03 4.83200E+03 4.83300E+03 &
     4.83400E+03 4.83500E+03 4.83600E+03 4.83700E+03 4.83800E+03 4.83900E+03 4.84000E+03 4.84100E+03 &
     4.84200E+03 4.84300E+03 4.84400E+03 4.84500E+03 4.84600E+03 4.84700E+03 4.84800E+03 4.84900E+03 &
     4.85000E+03 4.85100E+03 4.85200E+03 4.85300E+03 4.85400E+03 4.85500E+03 4.85600E+03 4.85700E+03 &
     4.85800E+03 4.85900E+03 4.86000E+03 4.86100E+03 4.86200E+03 4.86300E+03 4.86400E+03 4.86500E+03 &
     4.86600E+03 4.86700E+03 4.86800E+03 4.86900E+03 4.87000E+03 4.87100E+03 4.87200E+03 4.87300E+03 &
     4.87400E+03 4.87500E+03 4.87600E+03 4.87700E+03 4.87800E+03 4.87900E+03 4.88000E+03 4.88100E+03 &
     4.88200E+03 4.88300E+03 4.88400E+03 4.88500E+03 4.88600E+03 4.88700E+03 4.88800E+03 4.88900E+03 &
     4.89000E+03 4.89100E+03 4.89200E+03 4.89300E+03 4.89400E+03 4.89500E+03 4.89600E+03 4.89700E+03 &
     4.89800E+03 4.89900E+03 4.90000E+03 4.90100E+03 4.90200E+03 4.90300E+03 4.90400E+03 4.90500E+03 &
     4.90600E+03 4.90700E+03 4.90800E+03 4.90900E+03 4.91000E+03 4.91100E+03 4.91200E+03 4.91300E+03 &
     4.91400E+03 4.91500E+03 4.91600E+03 4.91700E+03 4.91800E+03 4.91900E+03 4.92000E+03 4.92100E+03 &
     4.92200E+03 4.92300E+03 4.92400E+03 4.92500E+03 4.92600E+03 4.92700E+03 4.92800E+03 4.92900E+03 &
     4.93000E+03 4.93100E+03 4.93200E+03 4.93300E+03 4.93400E+03 4.93500E+03 4.93600E+03 4.93700E+03 &
     4.93800E+03 4.93900E+03 4.94000E+03 4.94100E+03 4.94200E+03 4.94300E+03 4.94400E+03 4.94500E+03 &
     4.94600E+03 4.94700E+03 4.94800E+03 4.94900E+03 4.95000E+03 4.95100E+03 4.95200E+03 4.95300E+03 &
     4.95400E+03 4.95500E+03 4.95600E+03 4.95700E+03 4.95800E+03 4.95900E+03 4.96000E+03 4.96100E+03 &
     4.96200E+03 4.96300E+03 4.96400E+03 4.96500E+03 4.96600E+03 4.96700E+03 4.96800E+03 4.96900E+03 &
     4.97000E+03 4.97100E+03 4.97200E+03 4.97300E+03 4.97400E+03 4.97500E+03 4.97600E+03 4.97700E+03 &
     4.97800E+03 4.97900E+03 4.98000E+03 4.98100E+03 4.98200E+03 4.98300E+03 4.98400E+03 4.98500E+03 &
     4.98600E+03 4.98700E+03 4.98800E+03 4.98900E+03 4.99000E+03 4.99100E+03 4.99200E+03 4.99300E+03 &
     4.99400E+03 4.99500E+03 4.99600E+03 4.99700E+03 4.99800E+03 4.99900E+03 5.00000E+03 5.00100E+03 &
     5.00200E+03 5.00300E+03 5.00400E+03 5.00500E+03 5.00600E+03 5.00700E+03 5.00800E+03 5.00900E+03 &
     5.01000E+03 5.01100E+03 5.01200E+03 5.01300E+03 5.01400E+03 5.01500E+03 5.01600E+03 5.01700E+03 &
     5.01800E+03 5.01900E+03 5.02000E+03 5.02100E+03 5.02200E+03 5.02300E+03 5.02400E+03 5.02500E+03 &
     5.02600E+03 5.02700E+03 5.02800E+03 5.02900E+03 5.03000E+03 5.03100E+03 5.03200E+03 5.03300E+03 &
     5.03400E+03 5.03500E+03 5.03600E+03 5.03700E+03 5.03800E+03 5.03900E+03 5.04000E+03 5.04100E+03 &
     5.04200E+03 5.04300E+03 5.04400E+03 5.04500E+03 5.04600E+03 5.04700E+03 5.04800E+03 5.04900E+03 &
     5.05000E+03 5.05100E+03 5.05200E+03 5.05300E+03 5.05400E+03 5.05500E+03 5.05600E+03 5.05700E+03 &
     5.05800E+03 5.05900E+03 5.06000E+03 5.06100E+03 5.06200E+03 5.06300E+03 5.06400E+03 5.06500E+03 &
     5.06600E+03 5.06700E+03 5.06800E+03 5.06900E+03 5.07000E+03 5.07100E+03 5.07200E+03 5.07300E+03 &
     5.07400E+03 5.07500E+03 5.07600E+03 5.07700E+03 5.07800E+03 5.07900E+03 5.08000E+03 5.08100E+03 &
     5.08200E+03 5.08300E+03 5.08400E+03 5.08500E+03 5.08600E+03 5.08700E+03 5.08800E+03 5.08900E+03 &
     5.09000E+03 5.09100E+03 5.09200E+03 5.09300E+03 5.09400E+03 5.09500E+03 5.09600E+03 5.09700E+03 &
     5.09800E+03 5.09900E+03 5.10000E+03 5.10100E+03 5.10200E+03 5.10300E+03 5.10400E+03 5.10500E+03 &
     5.10600E+03 5.10700E+03 5.10800E+03 5.10900E+03 5.11000E+03 5.11100E+03 5.11200E+03 5.11300E+03 &
     5.11400E+03 5.11500E+03 5.11600E+03 5.11700E+03 5.11800E+03 5.11900E+03 5.12000E+03 5.12100E+03 &
     5.12200E+03 5.12300E+03 5.12400E+03 5.12500E+03 5.12600E+03 5.12700E+03 5.12800E+03 5.12900E+03 &
     5.13000E+03 5.13100E+03 5.13200E+03 5.13300E+03 5.13400E+03 5.13500E+03 5.13600E+03 5.13700E+03 &
     5.13800E+03 5.13900E+03 5.14000E+03 5.14100E+03 5.14200E+03 5.14300E+03 5.14400E+03 5.14500E+03 &
     5.14600E+03 5.14700E+03 5.14800E+03 5.14900E+03 5.15000E+03 5.15100E+03 5.15200E+03 5.15300E+03 &
     5.15400E+03 5.15500E+03 5.15600E+03 5.15700E+03 5.15800E+03 5.15900E+03 5.16000E+03 5.16100E+03 &
     5.16200E+03 5.16300E+03 5.16400E+03 5.16500E+03 5.16600E+03 5.16700E+03 5.16800E+03 5.16900E+03 &
     5.17000E+03 5.17100E+03 5.17200E+03 5.17300E+03 5.17400E+03 5.17500E+03 5.17600E+03 5.17700E+03 &
     5.17800E+03 5.17900E+03 5.18000E+03 5.18100E+03 5.18200E+03 5.18300E+03 5.18400E+03 5.18500E+03 &
     5.18600E+03 5.18700E+03 5.18800E+03 5.18900E+03 5.19000E+03 5.19100E+03 5.19200E+03 5.19300E+03 &
     5.19400E+03 5.19500E+03 5.19600E+03 5.19700E+03 5.19800E+03 5.19900E+03 5.20000E+03 5.20100E+03 &
     5.20200E+03 5.20300E+03 5.20400E+03 5.20500E+03 5.20600E+03 5.20700E+03 5.20800E+03 5.20900E+03 &
     5.21000E+03 5.21100E+03 5.21200E+03 5.21300E+03 5.21400E+03 5.21500E+03 5.21600E+03 5.21700E+03 &
     5.21800E+03 5.21900E+03 5.22000E+03 5.22100E+03 5.22200E+03 5.22300E+03 5.22400E+03 5.22500E+03 &
     5.22600E+03 5.22700E+03 5.22800E+03 5.22900E+03 5.23000E+03 5.23100E+03 5.23200E+03 5.23300E+03 &
     5.23400E+03 5.23500E+03 5.23600E+03 5.23700E+03 5.23800E+03 5.23900E+03 5.24000E+03 5.24100E+03 &
     5.24200E+03 5.24300E+03 5.24400E+03 5.24500E+03 5.24600E+03 5.24700E+03 5.24800E+03 5.24900E+03 &
     5.25000E+03 5.25100E+03 5.25200E+03 5.25300E+03 5.25400E+03 5.25500E+03 5.25600E+03 5.25700E+03 &
     5.25800E+03 5.25900E+03 5.26000E+03 5.26100E+03 5.26200E+03 5.26300E+03 5.26400E+03 5.26500E+03 &
     5.26600E+03 5.26700E+03 5.26800E+03 5.26900E+03 5.27000E+03 5.27100E+03 5.27200E+03 5.27300E+03 &
     5.27400E+03 5.27500E+03 5.27600E+03 5.27700E+03 5.27800E+03 5.27900E+03 5.28000E+03 5.28100E+03 &
     5.28200E+03 5.28300E+03 5.28400E+03 5.28500E+03 5.28600E+03 5.28700E+03 5.28800E+03 5.28900E+03 &
     5.29000E+03 5.29100E+03 5.29200E+03 5.29300E+03 5.29400E+03 5.29500E+03 5.29600E+03 5.29700E+03 &
     5.29800E+03 5.29900E+03 5.30000E+03 5.30100E+03 5.30200E+03 5.30300E+03 5.30400E+03 5.30500E+03 &
     5.30600E+03 5.30700E+03 5.30800E+03 5.30900E+03 5.31000E+03 5.31100E+03 5.31200E+03 5.31300E+03 &
     5.31400E+03 5.31500E+03 5.31600E+03 5.31700E+03 5.31800E+03 5.31900E+03 5.32000E+03 5.32100E+03 &
     5.32200E+03 5.32300E+03 5.32400E+03 5.32500E+03 5.32600E+03 5.32700E+03 5.32800E+03 5.32900E+03 &
     5.33000E+03 5.33100E+03 5.33200E+03 5.33300E+03 5.33400E+03 5.33500E+03 5.33600E+03 5.33700E+03 &
     5.33800E+03 5.33900E+03 5.34000E+03 5.34100E+03 5.34200E+03 5.34300E+03 5.34400E+03 5.34500E+03 &
     5.34600E+03 5.34700E+03 5.34800E+03 5.34900E+03 5.35000E+03 5.35100E+03 5.35200E+03 5.35300E+03 &
     5.35400E+03 5.35500E+03 5.35600E+03 5.35700E+03 5.35800E+03 5.35900E+03 5.36000E+03 5.36100E+03 &
     5.36200E+03 5.36300E+03 5.36400E+03 5.36500E+03 5.36600E+03 5.36700E+03 5.36800E+03 5.36900E+03 &
     5.37000E+03 5.37100E+03 5.37200E+03 5.37300E+03 5.37400E+03 5.37500E+03 5.37600E+03 5.37700E+03 &
     5.37800E+03 5.37900E+03 5.38000E+03 5.38100E+03 5.38200E+03 5.38300E+03 5.38400E+03 5.38500E+03 &
     5.38600E+03 5.38700E+03 5.38800E+03 5.38900E+03 5.39000E+03 5.39100E+03 5.39200E+03 5.39300E+03 &
     5.39400E+03 5.39500E+03 5.39600E+03 5.39700E+03 5.39800E+03 5.39900E+03 5.40000E+03 5.40100E+03 &
     5.40200E+03 5.40300E+03 5.40400E+03 5.40500E+03 5.40600E+03 5.40700E+03 5.40800E+03 5.40900E+03 &
     5.41000E+03 5.41100E+03 5.41200E+03 5.41300E+03 5.41400E+03 5.41500E+03 5.41600E+03 5.41700E+03 &
     5.41800E+03 5.41900E+03 5.42000E+03 5.42100E+03 5.42200E+03 5.42300E+03 5.42400E+03 5.42500E+03 &
     5.42600E+03 5.42700E+03 5.42800E+03 5.42900E+03 5.43000E+03 5.43100E+03 5.43200E+03 5.43300E+03 &
     5.43400E+03 5.43500E+03 5.43600E+03 5.43700E+03 5.43800E+03 5.43900E+03 5.44000E+03 5.44100E+03 &
     5.44200E+03 5.44300E+03 5.44400E+03 5.44500E+03 5.44600E+03 5.44700E+03 5.44800E+03 5.44900E+03 &
     5.45000E+03 5.45100E+03 5.45200E+03 5.45300E+03 5.45400E+03 5.45500E+03 5.45600E+03 5.45700E+03 &
     5.45800E+03 5.45900E+03 5.46000E+03 5.46100E+03 5.46200E+03 5.46300E+03 5.46400E+03 5.46500E+03 &
     5.46600E+03 5.46700E+03 5.46800E+03 5.46900E+03 5.47000E+03 5.47100E+03 5.47200E+03 5.47300E+03 &
     5.47400E+03 5.47500E+03 5.47600E+03 5.47700E+03 5.47800E+03 5.47900E+03 5.48000E+03 5.48100E+03 &
     5.48200E+03 5.48300E+03 5.48400E+03 5.48500E+03 5.48600E+03 5.48700E+03 5.48800E+03 5.48900E+03 &
     5.49000E+03 5.49100E+03 5.49200E+03 5.49300E+03 5.49400E+03 5.49500E+03 5.49600E+03 5.49700E+03 &
     5.49800E+03 5.49900E+03 5.50000E+03 5.50100E+03 5.50200E+03 5.50300E+03 5.50400E+03 5.50500E+03 &
     5.50600E+03 5.50700E+03 5.50800E+03 5.50900E+03 5.51000E+03 5.51100E+03 5.51200E+03 5.51300E+03 &
     5.51400E+03 5.51500E+03 5.51600E+03 5.51700E+03 5.51800E+03 5.51900E+03 5.52000E+03 5.52100E+03 &
     5.52200E+03 5.52300E+03 5.52400E+03 5.52500E+03 5.52600E+03 5.52700E+03 5.52800E+03 5.52900E+03 &
     5.53000E+03 5.53100E+03 5.53200E+03 5.53300E+03 5.53400E+03 5.53500E+03 5.53600E+03 5.53700E+03 &
     5.53800E+03 5.53900E+03 5.54000E+03 5.54100E+03 5.54200E+03 5.54300E+03 5.54400E+03 5.54500E+03 &
     5.54600E+03 5.54700E+03 5.54800E+03 5.54900E+03 5.55000E+03 5.55100E+03 5.55200E+03 5.55300E+03 &
     5.55400E+03 5.55500E+03 5.55600E+03 5.55700E+03 5.55800E+03 5.55900E+03 5.56000E+03 5.56100E+03 &
     5.56200E+03 5.56300E+03 5.56400E+03 5.56500E+03 5.56600E+03 5.56700E+03 5.56800E+03 5.56900E+03 &
     5.57000E+03 5.57100E+03 5.57200E+03 5.57300E+03 5.57400E+03 5.57500E+03 5.57600E+03 5.57700E+03 &
     5.57800E+03 5.57900E+03 5.58000E+03 5.58100E+03 5.58200E+03 5.58300E+03 5.58400E+03 5.58500E+03 &
     5.58600E+03 5.58700E+03 5.58800E+03 5.58900E+03 5.59000E+03 5.59100E+03 5.59200E+03 5.59300E+03 &
     5.59400E+03 5.59500E+03 5.59600E+03 5.59700E+03 5.59800E+03 5.59900E+03 5.60000E+03 5.60100E+03 &
     5.60200E+03 5.60300E+03 5.60400E+03 5.60500E+03 5.60600E+03 5.60700E+03 5.60800E+03 5.60900E+03 &
     5.61000E+03 5.61100E+03 5.61200E+03 5.61300E+03 5.61400E+03 5.61500E+03 5.61600E+03 5.61700E+03 &
     5.61800E+03 5.61900E+03 5.62000E+03 5.62100E+03 5.62200E+03 5.62300E+03 5.62400E+03 5.62500E+03 &
     5.62600E+03 5.62700E+03 5.62800E+03 5.62900E+03 5.63000E+03 5.63100E+03 5.63200E+03 5.63300E+03 &
     5.63400E+03 5.63500E+03 5.63600E+03 5.63700E+03 5.63800E+03 5.63900E+03 5.64000E+03 5.64100E+03 &
     5.64200E+03 5.64300E+03 5.64400E+03 5.64500E+03 5.64600E+03 5.64700E+03 5.64800E+03 5.64900E+03 &
     5.65000E+03 5.65100E+03 5.65200E+03 5.65300E+03 5.65400E+03 5.65500E+03 5.65600E+03 5.65700E+03 &
     5.65800E+03 5.65900E+03 5.66000E+03 5.66100E+03 5.66200E+03 5.66300E+03 5.66400E+03 5.66500E+03 &
     5.66600E+03 5.66700E+03 5.66800E+03 5.66900E+03 5.67000E+03 5.67100E+03 5.67200E+03 5.67300E+03 &
     5.67400E+03 5.67500E+03 5.67600E+03 5.67700E+03 5.67800E+03 5.67900E+03 5.68000E+03 5.68100E+03 &
     5.68200E+03 5.68300E+03 5.68400E+03 5.68500E+03 5.68600E+03 5.68700E+03 5.68800E+03 5.68900E+03 &
     5.69000E+03 5.69100E+03 5.69200E+03 5.69300E+03 5.69400E+03 5.69500E+03 5.69600E+03 5.69700E+03 &
     5.69800E+03 5.69900E+03 5.70000E+03 5.70100E+03 5.70200E+03 5.70300E+03 5.70400E+03 5.70500E+03 &
     5.70600E+03 5.70700E+03 5.70800E+03 5.70900E+03 5.71000E+03 5.71100E+03 5.71200E+03 5.71300E+03 &
     5.71400E+03 5.71500E+03 5.71600E+03 5.71700E+03 5.71800E+03 5.71900E+03 5.72000E+03 5.72100E+03 &
     5.72200E+03 5.72300E+03 5.72400E+03 5.72500E+03 5.72600E+03 5.72700E+03 5.72800E+03 5.72900E+03 &
     5.73000E+03 5.73100E+03 5.73200E+03 5.73300E+03 5.73400E+03 5.73500E+03 5.73600E+03 5.73700E+03 &
     5.73800E+03 5.73900E+03 5.74000E+03 5.74100E+03 5.74200E+03 5.74300E+03 5.74400E+03 5.74500E+03 &
     5.74600E+03 5.74700E+03 5.74800E+03 5.74900E+03 5.75000E+03 5.75100E+03 5.75200E+03 5.75300E+03 &
     5.75400E+03 5.75500E+03 5.75600E+03 5.75700E+03 5.75800E+03 5.75900E+03 5.76000E+03 5.76100E+03 &
     5.76200E+03 5.76300E+03 5.76400E+03 5.76500E+03 5.76600E+03 5.76700E+03 5.76800E+03 5.76900E+03 &
     5.77000E+03 5.77100E+03 5.77200E+03 5.77300E+03 5.77400E+03 5.77500E+03 5.77600E+03 5.77700E+03 &
     5.77800E+03 5.77900E+03 5.78000E+03 5.78100E+03 5.78200E+03 5.78300E+03 5.78400E+03 5.78500E+03 &
     5.78600E+03 5.78700E+03 5.78800E+03 5.78900E+03 5.79000E+03 5.79100E+03 5.79200E+03 5.79300E+03 &
     5.79400E+03 5.79500E+03 5.79600E+03 5.79700E+03 5.79800E+03 5.79900E+03 5.80000E+03 5.80100E+03 &
     5.80200E+03 5.80300E+03 5.80400E+03 5.80500E+03 5.80600E+03 5.80700E+03 5.80800E+03 5.80900E+03 &
     5.81000E+03 5.81100E+03 5.81200E+03 5.81300E+03 5.81400E+03 5.81500E+03 5.81600E+03 5.81700E+03 &
     5.81800E+03 5.81900E+03 5.82000E+03 5.82100E+03 5.82200E+03 5.82300E+03 5.82400E+03 5.82500E+03 &
     5.82600E+03 5.82700E+03 5.82800E+03 5.82900E+03 5.83000E+03 5.83100E+03 5.83200E+03 5.83300E+03 &
     5.83400E+03 5.83500E+03 5.83600E+03 5.83700E+03 5.83800E+03 5.83900E+03 5.84000E+03 5.84100E+03 &
     5.84200E+03 5.84300E+03 5.84400E+03 5.84500E+03 5.84600E+03 5.84700E+03 5.84800E+03 5.84900E+03 &
     5.85000E+03 5.85100E+03 5.85200E+03 5.85300E+03 5.85400E+03 5.85500E+03 5.85600E+03 5.85700E+03 &
     5.85800E+03 5.85900E+03 5.86000E+03 5.86100E+03 5.86200E+03 5.86300E+03 5.86400E+03 5.86500E+03 &
     5.86600E+03 5.86700E+03 5.86800E+03 5.86900E+03 5.87000E+03 5.87100E+03 5.87200E+03 5.87300E+03 &
     5.87400E+03 5.87500E+03 5.87600E+03 5.87700E+03 5.87800E+03 5.87900E+03 5.88000E+03 5.88100E+03 &
     5.88200E+03 5.88300E+03 5.88400E+03 5.88500E+03 5.88600E+03 5.88700E+03 5.88800E+03 5.88900E+03 &
     5.89000E+03 5.89100E+03 5.89200E+03 5.89300E+03 5.89400E+03 5.89500E+03 5.89600E+03 5.89700E+03 &
     5.89800E+03 5.89900E+03 5.90000E+03 5.90100E+03 5.90200E+03 5.90300E+03 5.90400E+03 5.90500E+03 &
     5.90600E+03 5.90700E+03 5.90800E+03 5.90900E+03 5.91000E+03 5.91100E+03 5.91200E+03 5.91300E+03 &
     5.91400E+03 5.91500E+03 5.91600E+03 5.91700E+03 5.91800E+03 5.91900E+03 5.92000E+03 5.92100E+03 &
     5.92200E+03 5.92300E+03 5.92400E+03 5.92500E+03 5.92600E+03 5.92700E+03 5.92800E+03 5.92900E+03 &
     5.93000E+03 5.93100E+03 5.93200E+03 5.93300E+03 5.93400E+03 5.93500E+03 5.93600E+03 5.93700E+03 &
     5.93800E+03 5.93900E+03 5.94000E+03 5.94100E+03 5.94200E+03 5.94300E+03 5.94400E+03 5.94500E+03 &
     5.94600E+03 5.94700E+03 5.94800E+03 5.94900E+03 5.95000E+03 5.95100E+03 5.95200E+03 5.95300E+03 &
     5.95400E+03 5.95500E+03 5.95600E+03 5.95700E+03 5.95800E+03 5.95900E+03 5.96000E+03 5.96100E+03 &
     5.96200E+03 5.96300E+03 5.96400E+03 5.96500E+03 5.96600E+03 5.96700E+03 5.96800E+03 5.96900E+03 &
     5.97000E+03 5.97100E+03 5.97200E+03 5.97300E+03 5.97400E+03 5.97500E+03 5.97600E+03 5.97700E+03 &
     5.97800E+03 5.97900E+03 5.98000E+03 5.98100E+03 5.98200E+03 5.98300E+03 5.98400E+03 5.98500E+03 &
     5.98600E+03 5.98700E+03 5.98800E+03 5.98900E+03 5.99000E+03 5.99100E+03 5.99200E+03 5.99300E+03 &
     5.99400E+03 5.99500E+03 5.99600E+03 5.99700E+03 5.99800E+03 5.99900E+03 6.00000E+03 1.59390E-10 &
     1.47540E-10 1.36810E-10 1.28590E-10 1.22630E-10 1.18460E-10 1.15610E-10 1.13770E-10 1.12720E-10 &
     1.12220E-10 1.12290E-10 1.12650E-10 1.13270E-10 1.13990E-10 1.15010E-10 1.16310E-10 1.17720E-10 &
     1.19350E-10 1.20960E-10 1.22660E-10 1.24440E-10 1.26240E-10 1.28110E-10 1.30060E-10 1.32020E-10 &
     1.34050E-10 1.36130E-10 1.38230E-10 1.40400E-10 1.42550E-10 1.44780E-10 1.47030E-10 1.49300E-10 &
     1.51650E-10 1.53970E-10 1.56340E-10 1.58750E-10 1.61210E-10 1.63680E-10 1.66160E-10 1.68650E-10 &
     1.71160E-10 1.73690E-10 1.76270E-10 1.78890E-10 1.81550E-10 1.84230E-10 1.86940E-10 1.89670E-10 &
     1.92420E-10 1.95200E-10 1.97990E-10 2.00800E-10 2.03620E-10 2.06480E-10 2.09390E-10 2.12340E-10 &
     2.15330E-10 2.18340E-10 2.21390E-10 2.24460E-10 2.27550E-10 2.30660E-10 2.33790E-10 2.36950E-10 &
     2.40130E-10 2.43340E-10 2.46590E-10 2.49880E-10 2.53210E-10 2.56600E-10 2.60030E-10 2.63500E-10 &
     2.67000E-10 2.70500E-10 2.74000E-10 2.77510E-10 2.81040E-10 2.84600E-10 2.88220E-10 2.91880E-10 &
     2.95590E-10 2.99340E-10 3.03120E-10 3.06940E-10 3.10790E-10 3.14660E-10 3.18560E-10 3.22490E-10 &
     3.26440E-10 3.30430E-10 3.34440E-10 3.38480E-10 3.42550E-10 3.46640E-10 3.50770E-10 3.54930E-10 &
     3.59110E-10 3.63330E-10 3.67570E-10 3.71840E-10 3.76140E-10 3.80460E-10 3.84810E-10 3.89180E-10 &
     3.93570E-10 3.97990E-10 4.02430E-10 4.06890E-10 4.11360E-10 4.15850E-10 4.20370E-10 4.24890E-10 &
     4.29440E-10 4.33990E-10 4.38560E-10 4.43150E-10 4.47740E-10 4.52350E-10 4.56960E-10 4.61580E-10 &
     4.66210E-10 4.70850E-10 4.75490E-10 4.80140E-10 4.84790E-10 4.89440E-10 4.94090E-10 4.98740E-10 &
     5.03390E-10 5.08040E-10 5.12680E-10 5.17320E-10 5.21960E-10 5.26580E-10 5.31190E-10 5.35800E-10 &
     5.40380E-10 5.44960E-10 5.49520E-10 5.54060E-10 5.58590E-10 5.63090E-10 5.67570E-10 5.72030E-10 &
     5.76460E-10 5.80860E-10 5.85240E-10 5.89590E-10 5.93910E-10 5.98190E-10 6.02450E-10 6.06680E-10 &
     6.10880E-10 6.15040E-10 6.19180E-10 6.23290E-10 6.27370E-10 6.31420E-10 6.35430E-10 6.39410E-10 &
     6.43360E-10 6.47260E-10 6.51130E-10 6.54950E-10 6.58740E-10 6.62470E-10 6.66170E-10 6.69810E-10 &
     6.73410E-10 6.76950E-10 6.80440E-10 6.83880E-10 6.87270E-10 6.90590E-10 6.93860E-10 6.97070E-10 &
     7.00210E-10 7.03300E-10 7.06320E-10 7.09280E-10 7.12170E-10 7.15000E-10 7.17770E-10 7.20460E-10 &
     7.23090E-10 7.25660E-10 7.28150E-10 7.30580E-10 7.32930E-10 7.35220E-10 7.37430E-10 7.39580E-10 &
     7.41650E-10 7.43650E-10 7.45580E-10 7.47440E-10 7.49230E-10 7.50940E-10 7.52580E-10 7.54140E-10 &
     7.55620E-10 7.57020E-10 7.58340E-10 7.59580E-10 7.60740E-10 7.61810E-10 7.62790E-10 7.63700E-10 &
     7.64520E-10 7.65260E-10 7.65920E-10 7.66490E-10 7.66980E-10 7.67390E-10 7.67720E-10 7.67960E-10 &
     7.68120E-10 7.68190E-10 7.68180E-10 7.68080E-10 7.67890E-10 7.67620E-10 7.67260E-10 7.66820E-10 &
     7.66290E-10 7.65680E-10 7.64990E-10 7.64210E-10 7.63340E-10 7.62400E-10 7.61370E-10 7.60250E-10 &
     7.59060E-10 7.57780E-10 7.56430E-10 7.54990E-10 7.53470E-10 7.51870E-10 7.50200E-10 7.48440E-10 &
     7.46610E-10 7.44700E-10 7.42710E-10 7.40650E-10 7.38510E-10 7.36300E-10 7.34020E-10 7.31670E-10 &
     7.29240E-10 7.26740E-10 7.24180E-10 7.21540E-10 7.18830E-10 7.16050E-10 7.13210E-10 7.10300E-10 &
     7.07320E-10 7.04290E-10 7.01190E-10 6.98030E-10 6.94820E-10 6.91550E-10 6.88220E-10 6.84840E-10 &
     6.81400E-10 6.77910E-10 6.74370E-10 6.70770E-10 6.67120E-10 6.63430E-10 6.59690E-10 6.55900E-10 &
     6.52060E-10 6.48190E-10 6.44260E-10 6.40300E-10 6.36300E-10 6.32260E-10 6.28190E-10 6.24070E-10 &
     6.19930E-10 6.15750E-10 6.11540E-10 6.07290E-10 6.03020E-10 5.98720E-10 5.94400E-10 5.90050E-10 &
     5.85680E-10 5.81280E-10 5.76870E-10 5.72430E-10 5.67980E-10 5.63510E-10 5.59020E-10 5.54510E-10 &
     5.50000E-10 5.45470E-10 5.40930E-10 5.36380E-10 5.31830E-10 5.27270E-10 5.22710E-10 5.18140E-10 &
     5.13570E-10 5.09000E-10 5.04440E-10 4.99870E-10 4.95310E-10 4.90760E-10 4.86210E-10 4.81670E-10 &
     4.77130E-10 4.72610E-10 4.68100E-10 4.63610E-10 4.59130E-10 4.54660E-10 4.50210E-10 4.45780E-10 &
     4.41370E-10 4.36980E-10 4.32620E-10 4.28270E-10 4.23950E-10 4.19660E-10 4.15390E-10 4.11150E-10 &
     4.06940E-10 4.02760E-10 3.98610E-10 3.94500E-10 3.90420E-10 3.86370E-10 3.82360E-10 3.78380E-10 &
     3.74450E-10 3.70550E-10 3.66690E-10 3.62870E-10 3.59100E-10 3.55370E-10 3.51670E-10 3.48030E-10 &
     3.44420E-10 3.40860E-10 3.37350E-10 3.33880E-10 3.30470E-10 3.27100E-10 3.23780E-10 3.20510E-10 &
     3.17290E-10 3.14130E-10 3.11020E-10 3.07960E-10 3.04950E-10 3.02000E-10 2.99100E-10 2.96260E-10 &
     2.93470E-10 2.90740E-10 2.88060E-10 2.85440E-10 2.82880E-10 2.80370E-10 2.77930E-10 2.75540E-10 &
     2.73210E-10 2.70930E-10 2.68720E-10 2.66570E-10 2.64470E-10 2.62430E-10 2.60460E-10 2.58540E-10 &
     2.56680E-10 2.54890E-10 2.53150E-10 2.51470E-10 2.49850E-10 2.48300E-10 2.46800E-10 2.45360E-10 &
     2.43980E-10 2.42660E-10 2.41400E-10 2.40200E-10 2.39060E-10 2.37970E-10 2.36950E-10 2.35980E-10 &
     2.35070E-10 2.34220E-10 2.33420E-10 2.32680E-10 2.32000E-10 2.31370E-10 2.30800E-10 2.30290E-10 &
     2.29830E-10 2.29420E-10 2.29070E-10 2.28770E-10 2.28530E-10 2.28340E-10 2.28200E-10 2.28110E-10 &
     2.28070E-10 2.28080E-10 2.28140E-10 2.28240E-10 2.28400E-10 2.28600E-10 2.28850E-10 2.29140E-10 &
     2.29480E-10 2.29870E-10 2.30300E-10 2.30770E-10 2.31280E-10 2.31830E-10 2.32430E-10 2.33070E-10 &
     2.33740E-10 2.34450E-10 2.35200E-10 2.35990E-10 2.36810E-10 2.37660E-10 2.38550E-10 2.39470E-10 &
     2.40420E-10 2.41400E-10 2.42420E-10 2.43460E-10 2.44540E-10 2.45640E-10 2.46770E-10 2.47920E-10 &
     2.49100E-10 2.50300E-10 2.51530E-10 2.52770E-10 2.54040E-10 2.55320E-10 2.56620E-10 2.57950E-10 &
     2.59290E-10 2.60640E-10 2.62020E-10 2.63400E-10 2.64800E-10 2.66220E-10 2.67640E-10 2.69070E-10 &
     2.70520E-10 2.71970E-10 2.73430E-10 2.74890E-10 2.76360E-10 2.77830E-10 2.79310E-10 2.80790E-10 &
     2.82280E-10 2.83760E-10 2.85250E-10 2.86740E-10 2.88220E-10 2.89700E-10 2.91180E-10 2.92650E-10 &
     2.94120E-10 2.95580E-10 2.97030E-10 2.98480E-10 2.99920E-10 3.01350E-10 3.02770E-10 3.04190E-10 &
     3.05590E-10 3.06970E-10 3.08350E-10 3.09710E-10 3.11050E-10 3.12380E-10 3.13690E-10 3.14990E-10 &
     3.16260E-10 3.17520E-10 3.18760E-10 3.19990E-10 3.21190E-10 3.22370E-10 3.23530E-10 3.24670E-10 &
     3.25790E-10 3.26890E-10 3.27950E-10 3.29000E-10 3.30020E-10 3.31010E-10 3.31980E-10 3.32920E-10 &
     3.33840E-10 3.34730E-10 3.35590E-10 3.36420E-10 3.37230E-10 3.38010E-10 3.38760E-10 3.39490E-10 &
     3.40180E-10 3.40850E-10 3.41480E-10 3.42090E-10 3.42670E-10 3.43210E-10 3.43720E-10 3.44200E-10 &
     3.44650E-10 3.45070E-10 3.45460E-10 3.45810E-10 3.46130E-10 3.46420E-10 3.46680E-10 3.46910E-10 &
     3.47110E-10 3.47270E-10 3.47400E-10 3.47500E-10 3.47570E-10 3.47610E-10 3.47610E-10 3.47580E-10 &
     3.47520E-10 3.47430E-10 3.47310E-10 3.47160E-10 3.46970E-10 3.46760E-10 3.46510E-10 3.46240E-10 &
     3.45930E-10 3.45590E-10 3.45230E-10 3.44830E-10 3.44410E-10 3.43950E-10 3.43470E-10 3.42960E-10 &
     3.42420E-10 3.41850E-10 3.41260E-10 3.40640E-10 3.39990E-10 3.39320E-10 3.38620E-10 3.37900E-10 &
     3.37150E-10 3.36380E-10 3.35590E-10 3.34770E-10 3.33920E-10 3.33060E-10 3.32170E-10 3.31260E-10 &
     3.30330E-10 3.29370E-10 3.28400E-10 3.27410E-10 3.26390E-10 3.25360E-10 3.24310E-10 3.23250E-10 &
     3.22160E-10 3.21070E-10 3.19950E-10 3.18820E-10 3.17680E-10 3.16530E-10 3.15360E-10 3.14180E-10 &
     3.12990E-10 3.11780E-10 3.10570E-10 3.09340E-10 3.08110E-10 3.06870E-10 3.05620E-10 3.04370E-10 &
     3.03100E-10 3.01830E-10 3.00560E-10 2.99280E-10 2.98000E-10 2.96710E-10 2.95420E-10 2.94130E-10 &
     2.92830E-10 2.91540E-10 2.90240E-10 2.88950E-10 2.87650E-10 2.86360E-10 2.85070E-10 2.83780E-10 &
     2.82500E-10 2.81220E-10 2.79940E-10 2.78670E-10 2.77400E-10 2.76140E-10 2.74890E-10 2.73650E-10 &
     2.72410E-10 2.71190E-10 2.69970E-10 2.68760E-10 2.67570E-10 2.66380E-10 2.65210E-10 2.64050E-10 &
     2.62910E-10 2.61770E-10 2.60660E-10 2.59550E-10 2.58460E-10 2.57390E-10 2.56330E-10 2.55290E-10 &
     2.54270E-10 2.53260E-10 2.52270E-10 2.51300E-10 2.50350E-10 2.49420E-10 2.48510E-10 2.47610E-10 &
     2.46740E-10 2.45890E-10 2.45060E-10 2.44250E-10 2.43470E-10 2.42700E-10 2.41960E-10 2.41240E-10 &
     2.40550E-10 2.39870E-10 2.39230E-10 2.38600E-10 2.38000E-10 2.37430E-10 2.36880E-10 2.36350E-10 &
     2.35850E-10 2.35380E-10 2.34930E-10 2.34510E-10 2.34110E-10 2.33740E-10 2.33390E-10 2.33070E-10 &
     2.32780E-10 2.32510E-10 2.32280E-10 2.32060E-10 2.31880E-10 2.31720E-10 2.31590E-10 2.31480E-10 &
     2.31400E-10 2.31350E-10 2.31320E-10 2.31320E-10 2.31350E-10 2.31410E-10 2.31490E-10 2.31590E-10 &
     2.31730E-10 2.31890E-10 2.32070E-10 2.32280E-10 2.32520E-10 2.32780E-10 2.33070E-10 2.33380E-10 &
     2.33720E-10 2.34080E-10 2.34470E-10 2.34880E-10 2.35320E-10 2.35780E-10 2.36270E-10 2.36770E-10 &
     2.37300E-10 2.37860E-10 2.38430E-10 2.39030E-10 2.39650E-10 2.40290E-10 2.40950E-10 2.41630E-10 &
     2.42330E-10 2.43060E-10 2.43800E-10 2.44560E-10 2.45340E-10 2.46140E-10 2.46960E-10 2.47790E-10 &
     2.48640E-10 2.49510E-10 2.50390E-10 2.51290E-10 2.52210E-10 2.53130E-10 2.54080E-10 2.55030E-10 &
     2.56010E-10 2.56990E-10 2.57990E-10 2.59000E-10 2.60020E-10 2.61050E-10 2.62090E-10 2.63140E-10 &
     2.64210E-10 2.65280E-10 2.66360E-10 2.67450E-10 2.68550E-10 2.69650E-10 2.70760E-10 2.71880E-10 &
     2.73000E-10 2.74120E-10 2.75250E-10 2.76380E-10 2.77520E-10 2.78660E-10 2.79800E-10 2.80930E-10 &
     2.82070E-10 2.83210E-10 2.84350E-10 2.85490E-10 2.86630E-10 2.87760E-10 2.88890E-10 2.90010E-10 &
     2.91140E-10 2.92250E-10 2.93370E-10 2.94470E-10 2.95580E-10 2.96670E-10 2.97760E-10 2.98850E-10 &
     2.99920E-10 3.00990E-10 3.02040E-10 3.03090E-10 3.04130E-10 3.05150E-10 3.06170E-10 3.07170E-10 &
     3.08160E-10 3.09130E-10 3.10090E-10 3.11040E-10 3.11970E-10 3.12880E-10 3.13780E-10 3.14670E-10 &
     3.15540E-10 3.16390E-10 3.17230E-10 3.18050E-10 3.18860E-10 3.19650E-10 3.20410E-10 3.21160E-10 &
     3.21890E-10 3.22600E-10 3.23290E-10 3.23960E-10 3.24600E-10 3.25230E-10 3.25830E-10 3.26410E-10 &
     3.26960E-10 3.27500E-10 3.28010E-10 3.28500E-10 3.28960E-10 3.29410E-10 3.29830E-10 3.30220E-10 &
     3.30600E-10 3.30940E-10 3.31270E-10 3.31570E-10 3.31840E-10 3.32080E-10 3.32300E-10 3.32490E-10 &
     3.32650E-10 3.32790E-10 3.32900E-10 3.32990E-10 3.33050E-10 3.33080E-10 3.33090E-10 3.33080E-10 &
     3.33040E-10 3.32970E-10 3.32870E-10 3.32740E-10 3.32590E-10 3.32400E-10 3.32190E-10 3.31950E-10 &
     3.31680E-10 3.31390E-10 3.31070E-10 3.30730E-10 3.30360E-10 3.29960E-10 3.29540E-10 3.29090E-10 &
     3.28620E-10 3.28110E-10 3.27580E-10 3.27020E-10 3.26430E-10 3.25820E-10 3.25180E-10 3.24520E-10 &
     3.23830E-10 3.23110E-10 3.22370E-10 3.21610E-10 3.20830E-10 3.20020E-10 3.19180E-10 3.18320E-10 &
     3.17430E-10 3.16520E-10 3.15590E-10 3.14630E-10 3.13640E-10 3.12640E-10 3.11610E-10 3.10560E-10 &
     3.09490E-10 3.08400E-10 3.07290E-10 3.06160E-10 3.05010E-10 3.03830E-10 3.02640E-10 3.01440E-10 &
     3.00210E-10 2.98960E-10 2.97700E-10 2.96420E-10 2.95120E-10 2.93800E-10 2.92460E-10 2.91110E-10 &
     2.89740E-10 2.88350E-10 2.86960E-10 2.85540E-10 2.84120E-10 2.82680E-10 2.81230E-10 2.79760E-10 &
     2.78290E-10 2.76800E-10 2.75300E-10 2.73790E-10 2.72270E-10 2.70740E-10 2.69190E-10 2.67640E-10 &
     2.66080E-10 2.64510E-10 2.62930E-10 2.61350E-10 2.59750E-10 2.58150E-10 2.56540E-10 2.54930E-10 &
     2.53310E-10 2.51690E-10 2.50070E-10 2.48440E-10 2.46800E-10 2.45170E-10 2.43530E-10 2.41890E-10 &
     2.40250E-10 2.38610E-10 2.36970E-10 2.35330E-10 2.33690E-10 2.32050E-10 2.30410E-10 2.28780E-10 &
     2.27140E-10 2.25510E-10 2.23890E-10 2.22270E-10 2.20650E-10 2.19030E-10 2.17430E-10 2.15820E-10 &
     2.14220E-10 2.12630E-10 2.11050E-10 2.09470E-10 2.07900E-10 2.06340E-10 2.04780E-10 2.03240E-10 &
     2.01710E-10 2.00180E-10 1.98670E-10 1.97160E-10 1.95670E-10 1.94190E-10 1.92710E-10 1.91250E-10 &
     1.89810E-10 1.88370E-10 1.86950E-10 1.85540E-10 1.84150E-10 1.82770E-10 1.81400E-10 1.80050E-10 &
     1.78710E-10 1.77380E-10 1.76080E-10 1.74790E-10 1.73510E-10 1.72250E-10 1.71010E-10 1.69780E-10 &
     1.68570E-10 1.67370E-10 1.66200E-10 1.65040E-10 1.63900E-10 1.62780E-10 1.61670E-10 1.60580E-10 &
     1.59520E-10 1.58470E-10 1.57440E-10 1.56430E-10 1.55440E-10 1.54460E-10 1.53510E-10 1.52580E-10 &
     1.51660E-10 1.50770E-10 1.49890E-10 1.49040E-10 1.48210E-10 1.47390E-10 1.46600E-10 1.45820E-10 &
     1.45070E-10 1.44340E-10 1.43620E-10 1.42930E-10 1.42260E-10 1.41610E-10 1.40980E-10 1.40360E-10 &
     1.39770E-10 1.39200E-10 1.38650E-10 1.38120E-10 1.37610E-10 1.37120E-10 1.36650E-10 1.36200E-10 &
     1.35780E-10 1.35370E-10 1.34980E-10 1.34610E-10 1.34260E-10 1.33930E-10 1.33620E-10 1.33330E-10 &
     1.33060E-10 1.32810E-10 1.32570E-10 1.32360E-10 1.32160E-10 1.31980E-10 1.31820E-10 1.31680E-10 &
     1.31550E-10 1.31440E-10 1.31350E-10 1.31280E-10 1.31220E-10 1.31180E-10 1.31160E-10 1.31150E-10 &
     1.31160E-10 1.31190E-10 1.31230E-10 1.31290E-10 1.31360E-10 1.31450E-10 1.31550E-10 1.31670E-10 &
     1.31800E-10 1.31940E-10 1.32100E-10 1.32270E-10 1.32450E-10 1.32650E-10 1.32850E-10 1.33080E-10 &
     1.33310E-10 1.33560E-10 1.33810E-10 1.34080E-10 1.34360E-10 1.34640E-10 1.34940E-10 1.35250E-10 &
     1.35560E-10 1.35890E-10 1.36220E-10 1.36570E-10 1.36920E-10 1.37270E-10 1.37640E-10 1.38010E-10 &
     1.38390E-10 1.38780E-10 1.39170E-10 1.39570E-10 1.39980E-10 1.40380E-10 1.40800E-10 1.41220E-10 &
     1.41640E-10 1.42060E-10 1.42490E-10 1.42930E-10 1.43360E-10 1.43800E-10 1.44240E-10 1.44680E-10 &
     1.45130E-10 1.45570E-10 1.46020E-10 1.46470E-10 1.46910E-10 1.47360E-10 1.47810E-10 1.48250E-10 &
     1.48700E-10 1.49140E-10 1.49590E-10 1.50030E-10 1.50460E-10 1.50900E-10 1.51330E-10 1.51760E-10 &
     1.52190E-10 1.52610E-10 1.53030E-10 1.53450E-10 1.53860E-10 1.54270E-10 1.54670E-10 1.55060E-10 &
     1.55450E-10 1.55840E-10 1.56210E-10 1.56580E-10 1.56950E-10 1.57310E-10 1.57660E-10 1.58000E-10 &
     1.58340E-10 1.58670E-10 1.58990E-10 1.59310E-10 1.59610E-10 1.59910E-10 1.60200E-10 1.60480E-10 &
     1.60750E-10 1.61010E-10 1.61260E-10 1.61500E-10 1.61730E-10 1.61960E-10 1.62170E-10 1.62370E-10 &
     1.62560E-10 1.62750E-10 1.62920E-10 1.63080E-10 1.63230E-10 1.63370E-10 1.63500E-10 1.63620E-10 &
     1.63720E-10 1.63820E-10 1.63900E-10 1.63970E-10 1.64030E-10 1.64080E-10 1.64120E-10 1.64140E-10 &
     1.64160E-10 1.64160E-10 1.64150E-10 1.64130E-10 1.64100E-10 1.64050E-10 1.64000E-10 1.63930E-10 &
     1.63850E-10 1.63750E-10 1.63650E-10 1.63530E-10 1.63400E-10 1.63260E-10 1.63110E-10 1.62950E-10 &
     1.62770E-10 1.62580E-10 1.62380E-10 1.62170E-10 1.61950E-10 1.61710E-10 1.61470E-10 1.61210E-10 &
     1.60940E-10 1.60660E-10 1.60360E-10 1.60060E-10 1.59750E-10 1.59420E-10 1.59090E-10 1.58740E-10 &
     1.58380E-10 1.58020E-10 1.57640E-10 1.57250E-10 1.56860E-10 1.56450E-10 1.56030E-10 1.55610E-10 &
     1.55170E-10 1.54730E-10 1.54270E-10 1.53810E-10 1.53340E-10 1.52860E-10 1.52370E-10 1.51870E-10 &
     1.51370E-10 1.50860E-10 1.50340E-10 1.49810E-10 1.49270E-10 1.48730E-10 1.48180E-10 1.47620E-10 &
     1.47060E-10 1.46490E-10 1.45910E-10 1.45330E-10 1.44740E-10 1.44150E-10 1.43550E-10 1.42950E-10 &
     1.42340E-10 1.41720E-10 1.41100E-10 1.40480E-10 1.39850E-10 1.39220E-10 1.38590E-10 1.37950E-10 &
     1.37310E-10 1.36660E-10 1.36010E-10 1.35360E-10 1.34710E-10 1.34060E-10 1.33400E-10 1.32740E-10 &
     1.32080E-10 1.31420E-10 1.30760E-10 1.30090E-10 1.29430E-10 1.28770E-10 1.28100E-10 1.27440E-10 &
     1.26770E-10 1.26110E-10 1.25450E-10 1.24780E-10 1.24120E-10 1.23460E-10 1.22810E-10 1.22150E-10 &
     1.21500E-10 1.20840E-10 1.20190E-10 1.19550E-10 1.18900E-10 1.18260E-10 1.17630E-10 1.16990E-10 &
     1.16360E-10 1.15740E-10 1.15110E-10 1.14500E-10 1.13880E-10 1.13270E-10 1.12670E-10 1.12070E-10 &
     1.11470E-10 1.10880E-10 1.10300E-10 1.09720E-10 1.09150E-10 1.08580E-10 1.08020E-10 1.07470E-10 &
     1.06920E-10 1.06380E-10 1.05840E-10 1.05310E-10 1.04790E-10 1.04280E-10 1.03770E-10 1.03270E-10 &
     1.02770E-10 1.02290E-10 1.01810E-10 1.01340E-10 1.00880E-10 1.00430E-10 9.99820E-11 9.95440E-11 &
     9.91150E-11 9.86940E-11 9.82820E-11 9.78780E-11 9.74830E-11 9.70970E-11 9.67190E-11 9.63500E-11 &
     9.59900E-11 9.56400E-11 9.52980E-11 9.49650E-11 9.46410E-11 9.43270E-11 9.40210E-11 9.37250E-11 &
     9.34380E-11 9.31610E-11 9.28920E-11 9.26330E-11 9.23830E-11 9.21430E-11 9.19120E-11 9.16900E-11 &
     9.14780E-11 9.12750E-11 9.10810E-11 9.08970E-11 9.07220E-11 9.05560E-11 9.04000E-11 9.02520E-11 &
     9.01140E-11 8.99850E-11 8.98660E-11 8.97550E-11 8.96530E-11 8.95600E-11 8.94760E-11 8.94010E-11 &
     8.93350E-11 8.92780E-11 8.92290E-11 8.91890E-11 8.91570E-11 8.91340E-11 8.91190E-11 8.91120E-11 &
     8.91130E-11 8.91230E-11 8.91410E-11 8.91660E-11 8.91990E-11 8.92400E-11 8.92880E-11 8.93440E-11 &
     8.94070E-11 8.94780E-11 8.95550E-11 8.96400E-11 8.97310E-11 8.98300E-11 8.99350E-11 9.00460E-11 &
     9.01640E-11 9.02880E-11 9.04180E-11 9.05540E-11 9.06960E-11 9.08430E-11 9.09970E-11 9.11550E-11 &
     9.13190E-11 9.14880E-11 9.16610E-11 9.18400E-11 9.20230E-11 9.22110E-11 9.24020E-11 9.25980E-11 &
     9.27980E-11 9.30020E-11 9.32100E-11 9.34200E-11 9.36350E-11 9.38520E-11 9.40730E-11 9.42960E-11 &
     9.45220E-11 9.47510E-11 9.49820E-11 9.52150E-11 9.54500E-11 9.56870E-11 9.59260E-11 9.61670E-11 &
     9.64080E-11 9.66510E-11 9.68950E-11 9.71400E-11 9.73860E-11 9.76320E-11 9.78790E-11 9.81250E-11 &
     9.83720E-11 9.86190E-11 9.88650E-11 9.91110E-11 9.93570E-11 9.96010E-11 9.98450E-11 1.00090E-10 &
     1.00330E-10 1.00570E-10 1.00810E-10 1.01040E-10 1.01280E-10 1.01510E-10 1.01740E-10 1.01970E-10 &
     1.02190E-10 1.02420E-10 1.02640E-10 1.02850E-10 1.03070E-10 1.03280E-10 1.03490E-10 1.03690E-10 &
     1.03890E-10 1.04090E-10 1.04280E-10 1.04470E-10 1.04650E-10 1.04830E-10 1.05010E-10 1.05180E-10 &
     1.05350E-10 1.05510E-10 1.05670E-10 1.05820E-10 1.05960E-10 1.06110E-10 1.06240E-10 1.06370E-10 &
     1.06500E-10 1.06620E-10 1.06730E-10 1.06840E-10 1.06940E-10 1.07040E-10 1.07130E-10 1.07210E-10 &
     1.07290E-10 1.07360E-10 1.07430E-10 1.07490E-10 1.07540E-10 1.07590E-10 1.07630E-10 1.07660E-10 &
     1.07690E-10 1.07710E-10 1.07720E-10 1.07730E-10 1.07730E-10 1.07720E-10 1.07710E-10 1.07690E-10 &
     1.07660E-10 1.07630E-10 1.07580E-10 1.07530E-10 1.07480E-10 1.07420E-10 1.07350E-10 1.07270E-10 &
     1.07180E-10 1.07090E-10 1.06990E-10 1.06890E-10 1.06770E-10 1.06650E-10 1.06520E-10 1.06390E-10 &
     1.06250E-10 1.06100E-10 1.05940E-10 1.05780E-10 1.05610E-10 1.05430E-10 1.05240E-10 1.05050E-10 &
     1.04850E-10 1.04650E-10 1.04440E-10 1.04220E-10 1.03990E-10 1.03760E-10 1.03520E-10 1.03280E-10 &
     1.03030E-10 1.02770E-10 1.02510E-10 1.02230E-10 1.01960E-10 1.01680E-10 1.01390E-10 1.01090E-10 &
     1.00790E-10 1.00480E-10 1.00170E-10 9.98540E-11 9.95300E-11 9.92010E-11 9.88660E-11 9.85270E-11 &
     9.81820E-11 9.78320E-11 9.74770E-11 9.71170E-11 9.67530E-11 9.63830E-11 9.60090E-11 9.56300E-11 &
     9.52470E-11 9.48600E-11 9.44680E-11 9.40720E-11 9.36720E-11 9.32680E-11 9.28600E-11 9.24490E-11 &
     9.20340E-11 9.16150E-11 9.11930E-11 9.07680E-11 9.03400E-11 8.99080E-11 8.94740E-11 8.90370E-11 &
     8.85970E-11 8.81550E-11 8.77100E-11 8.72630E-11 8.68130E-11 8.63610E-11 8.59080E-11 8.54530E-11 &
     8.49960E-11 8.45380E-11 8.40780E-11 8.36170E-11 8.31550E-11 8.26920E-11 8.22280E-11 8.17630E-11 &
     8.12980E-11 8.08320E-11 8.03650E-11 7.98990E-11 7.94320E-11 7.89650E-11 7.84980E-11 7.80320E-11 &
     7.75650E-11 7.71000E-11 7.66340E-11 7.61700E-11 7.57060E-11 7.52430E-11 7.47810E-11 7.43210E-11 &
     7.38610E-11 7.34030E-11 7.29470E-11 7.24920E-11 7.20380E-11 7.15860E-11 7.11360E-11 7.06880E-11 &
     7.02420E-11 6.97990E-11 6.93570E-11 6.89190E-11 6.84820E-11 6.80490E-11 6.76180E-11 6.71900E-11 &
     6.67640E-11 6.63420E-11 6.59230E-11 6.55070E-11 6.50940E-11 6.46850E-11 6.42790E-11 6.38760E-11 &
     6.34770E-11 6.30820E-11 6.26910E-11 6.23030E-11 6.19190E-11 6.15390E-11 6.11640E-11 6.07920E-11 &
     6.04250E-11 6.00610E-11 5.97020E-11 5.93480E-11 5.89980E-11 5.86520E-11 5.83110E-11 5.79740E-11 &
     5.76430E-11 5.73160E-11 5.69930E-11 5.66760E-11 5.63630E-11 5.60560E-11 5.57530E-11 5.54550E-11 &
     5.51620E-11 5.48750E-11 5.45920E-11 5.43140E-11 5.40420E-11 5.37740E-11 5.35120E-11 5.32550E-11 &
     5.30030E-11 5.27570E-11 5.25150E-11 5.22790E-11 5.20480E-11 5.18230E-11 5.16020E-11 5.13870E-11 &
     5.11780E-11 5.09730E-11 5.07740E-11 5.05800E-11 5.03910E-11 5.02070E-11 5.00290E-11 4.98560E-11 &
     4.96880E-11 4.95260E-11 4.93680E-11 4.92160E-11 4.90690E-11 4.89270E-11 4.87900E-11 4.86580E-11 &
     4.85310E-11 4.84100E-11 4.82930E-11 4.81810E-11 4.80740E-11 4.79720E-11 4.78750E-11 4.77830E-11 &
     4.76950E-11 4.76120E-11 4.75340E-11 4.74600E-11 4.73910E-11 4.73270E-11 4.72670E-11 4.72110E-11 &
     4.71600E-11 4.71130E-11 4.70700E-11 4.70310E-11 4.69970E-11 4.69670E-11 4.69400E-11 4.69180E-11 &
     4.69000E-11 4.68850E-11 4.68740E-11 4.68670E-11 4.68640E-11 4.68640E-11 4.68680E-11 4.68750E-11 &
     4.68850E-11 4.68990E-11 4.69150E-11 4.69350E-11 4.69580E-11 4.69840E-11 4.70130E-11 4.70450E-11 &
     4.70800E-11 4.71170E-11 4.71570E-11 4.71990E-11 4.72440E-11 4.72910E-11 4.73400E-11 4.73920E-11 &
     4.74450E-11 4.75010E-11 4.75590E-11 4.76180E-11 4.76800E-11 4.77430E-11 4.78070E-11 4.78740E-11 &
     4.79410E-11 4.80100E-11 4.80810E-11 4.81520E-11 4.82250E-11 4.82990E-11 4.83740E-11 4.84490E-11 &
     4.85260E-11 4.86030E-11 4.86810E-11 4.87600E-11 4.88380E-11 4.89180E-11 4.89980E-11 4.90770E-11 &
     4.91580E-11 4.92380E-11 4.93180E-11 4.93980E-11 4.94780E-11 4.95580E-11 4.96370E-11 4.97160E-11 &
     4.97950E-11 4.98730E-11 4.99500E-11 5.00270E-11 5.01030E-11 5.01780E-11 5.02530E-11 5.03260E-11 &
     5.03990E-11 5.04700E-11 5.05400E-11 5.06090E-11 5.06770E-11 5.07440E-11 5.08090E-11 5.08720E-11 &
     5.09340E-11 5.09950E-11 5.10540E-11 5.11110E-11 5.11670E-11 5.12210E-11 5.12730E-11 5.13230E-11 &
     5.13710E-11 5.14170E-11 5.14610E-11 5.15030E-11 5.15430E-11 5.15810E-11 5.16170E-11 5.16500E-11 &
     5.16810E-11 5.17100E-11 5.17360E-11 5.17600E-11 5.17810E-11 5.18000E-11 5.18160E-11 5.18300E-11 &
     5.18420E-11 5.18500E-11 5.18570E-11 5.18600E-11 5.18610E-11 5.18590E-11 5.18550E-11 5.18470E-11 &
     5.18370E-11 5.18240E-11 5.18080E-11 5.17900E-11 5.17680E-11 5.17440E-11 5.17170E-11 5.16870E-11 &
     5.16540E-11 5.16190E-11 5.15800E-11 5.15380E-11 5.14940E-11 5.14470E-11 5.13960E-11 5.13430E-11 &
     5.12870E-11 5.12280E-11 5.11660E-11 5.11010E-11 5.10340E-11 5.09630E-11 5.08900E-11 5.08140E-11 &
     5.07350E-11 5.06540E-11 5.05690E-11 5.04820E-11 5.03920E-11 5.03000E-11 5.02040E-11 5.01070E-11 &
     5.00060E-11 4.99030E-11 4.97970E-11 4.96880E-11 4.95770E-11 4.94640E-11 4.93480E-11 4.92290E-11 &
     4.91080E-11 4.89850E-11 4.88590E-11 4.87310E-11 4.86010E-11 4.84680E-11 4.83330E-11 4.81960E-11 &
     4.80570E-11 4.79150E-11 4.77720E-11 4.76260E-11 4.74780E-11 4.73280E-11 4.71770E-11 4.70230E-11 &
     4.68680E-11 4.67100E-11 4.65510E-11 4.63910E-11 4.62280E-11 4.60640E-11 4.58990E-11 4.57320E-11 &
     4.55630E-11 4.53930E-11 4.52220E-11 4.50490E-11 4.48750E-11 4.47000E-11 4.45230E-11 4.43460E-11 &
     4.41670E-11 4.39870E-11 4.38070E-11 4.36250E-11 4.34430E-11 4.32590E-11 4.30750E-11 4.28900E-11 &
     4.27050E-11 4.25190E-11 4.23320E-11 4.21450E-11 4.19570E-11 4.17690E-11 4.15810E-11 4.13930E-11 &
     4.12040E-11 4.10150E-11 4.08260E-11 4.06370E-11 4.04470E-11 4.02580E-11 4.00690E-11 3.98800E-11 &
     3.96910E-11 3.95030E-11 3.93150E-11 3.91270E-11 3.89390E-11 3.87520E-11 3.85650E-11 3.83790E-11 &
     3.81930E-11 3.80090E-11 3.78240E-11 3.76410E-11 3.74580E-11 3.72760E-11 3.70950E-11 3.69140E-11 &
     3.67350E-11 3.65560E-11 3.63790E-11 3.62030E-11 3.60270E-11 3.58530E-11 3.56800E-11 3.55090E-11 &
     3.53380E-11 3.51690E-11 3.50020E-11 3.48350E-11 3.46700E-11 3.45070E-11 3.43450E-11 3.41850E-11 &
     3.40260E-11 3.38690E-11 3.37130E-11 3.35590E-11 3.34070E-11 3.32570E-11 3.31080E-11 3.29610E-11 &
     3.28160E-11 3.26730E-11 3.25320E-11 3.23920E-11 3.22540E-11 3.21190E-11 3.19850E-11 3.18530E-11 &
     3.17240E-11 3.15960E-11 3.14700E-11 3.13470E-11 3.12250E-11 3.11060E-11 3.09890E-11 3.08740E-11 &
     3.07610E-11 3.06500E-11 3.05410E-11 3.04350E-11 3.03300E-11 3.02280E-11 3.01280E-11 3.00300E-11 &
     2.99350E-11 2.98410E-11 2.97500E-11 2.96610E-11 2.95740E-11 2.94900E-11 2.94080E-11 2.93280E-11 &
     2.92500E-11 2.91740E-11 2.91010E-11 2.90300E-11 2.89610E-11 2.88940E-11 2.88290E-11 2.87670E-11 &
     2.87070E-11 2.86480E-11 2.85930E-11 2.85390E-11 2.84870E-11 2.84380E-11 2.83900E-11 2.83450E-11 &
     2.83020E-11 2.82610E-11 2.82220E-11 2.81850E-11 2.81500E-11 2.81170E-11 2.80850E-11 2.80560E-11 &
     2.80290E-11 2.80040E-11 2.79810E-11 2.79590E-11 2.79400E-11 2.79220E-11 2.79060E-11 2.78920E-11 &
     2.78790E-11 2.78690E-11 2.78590E-11 2.78520E-11 2.78460E-11 2.78420E-11 2.78400E-11 2.78390E-11 &
     2.78390E-11 2.78410E-11 2.78450E-11 2.78490E-11 2.78560E-11 2.78630E-11 2.78720E-11 2.78830E-11 &
     2.78940E-11 2.79070E-11 2.79210E-11 2.79360E-11 2.79530E-11 2.79700E-11 2.79890E-11 2.80080E-11 &
     2.80290E-11 2.80500E-11 2.80720E-11 2.80960E-11 2.81200E-11 2.81450E-11 2.81710E-11 2.81970E-11 &
     2.82240E-11 2.82520E-11 2.82800E-11 2.83090E-11 2.83390E-11 2.83690E-11 2.84000E-11 2.84310E-11 &
     2.84620E-11 2.84940E-11 2.85260E-11 2.85590E-11 2.85910E-11 2.86240E-11 2.86570E-11 2.86900E-11 &
     2.87240E-11 2.87570E-11 2.87900E-11 2.88240E-11 2.88570E-11 2.88900E-11 2.89230E-11 2.89560E-11 &
     2.89890E-11 2.90220E-11 2.90540E-11 2.90860E-11 2.91170E-11 2.91490E-11 2.91800E-11 2.92100E-11 &
     2.92400E-11 2.92700E-11 2.92980E-11 2.93270E-11 2.93550E-11 2.93820E-11 2.94080E-11 2.94340E-11 &
     2.94590E-11 2.94840E-11 2.95080E-11 2.95310E-11 2.95530E-11 2.95740E-11 2.95940E-11 2.96140E-11 &
     2.96330E-11 2.96500E-11 2.96670E-11 2.96830E-11 2.96970E-11 2.97110E-11 2.97230E-11 2.97350E-11 &
     2.97450E-11 2.97540E-11 2.97630E-11 2.97690E-11 2.97750E-11 2.97800E-11 2.97830E-11 2.97850E-11 &
     2.97850E-11 2.97850E-11 2.97830E-11 2.97800E-11 2.97750E-11 2.97690E-11 2.97620E-11 2.97530E-11 &
     2.97430E-11 2.97320E-11 2.97190E-11 2.97050E-11 2.96890E-11 2.96720E-11 2.96540E-11 2.96340E-11 &
     2.96120E-11 2.95890E-11 2.95650E-11 2.95390E-11 2.95120E-11 2.94830E-11 2.94520E-11 2.94210E-11 &
     2.93870E-11 2.93520E-11 2.93160E-11 2.92780E-11 2.92390E-11 2.91980E-11 2.91560E-11 2.91120E-11 &
     2.90670E-11 2.90200E-11 2.89720E-11 2.89230E-11 2.88710E-11 2.88190E-11 2.87650E-11 2.87090E-11 &
     2.86530E-11 2.85940E-11 2.85350E-11 2.84740E-11 2.84110E-11 2.83470E-11 2.82820E-11 2.82150E-11 &
     2.81480E-11 2.80780E-11 2.80080E-11 2.79360E-11 2.78630E-11 2.77880E-11 2.77120E-11 2.76350E-11 &
     2.75570E-11 2.74780E-11 2.73970E-11 2.73150E-11 2.72320E-11 2.71480E-11 2.70620E-11 2.69760E-11 &
     2.68880E-11 2.67990E-11 2.67100E-11 2.66190E-11 2.65270E-11 2.64340E-11 2.63400E-11 2.62460E-11 &
     2.61500E-11 2.60530E-11 2.59560E-11 2.58570E-11 2.57580E-11 2.56580E-11 2.55570E-11 2.54550E-11 &
     2.53530E-11 2.52500E-11 2.51460E-11 2.50420E-11 2.49360E-11 2.48310E-11 2.47240E-11 2.46170E-11 &
     2.45100E-11 2.44020E-11 2.42940E-11 2.41850E-11 2.40750E-11 2.39650E-11 2.38550E-11 2.37440E-11 &
     2.36330E-11 2.35220E-11 2.34110E-11 2.32990E-11 2.31870E-11 2.30740E-11 2.29620E-11 2.28490E-11 &
     2.27360E-11 2.26240E-11 2.25110E-11 2.23980E-11 2.22850E-11 2.21720E-11 2.20590E-11 2.19460E-11 &
     2.18330E-11 2.17200E-11 2.16070E-11 2.14950E-11 2.13820E-11 2.12700E-11 2.11580E-11 2.10460E-11 &
     2.09350E-11 2.08230E-11 2.07120E-11 2.06020E-11 2.04910E-11 2.03810E-11 2.02720E-11 2.01630E-11 &
     2.00540E-11 1.99460E-11 1.98380E-11 1.97310E-11 1.96240E-11 1.95180E-11 1.94130E-11 1.93080E-11 &
     1.92040E-11 1.91000E-11 1.89970E-11 1.88950E-11 1.87930E-11 1.86920E-11 1.85920E-11 1.84920E-11 &
     1.83940E-11 1.82960E-11 1.81990E-11 1.81030E-11 1.80070E-11 1.79130E-11 1.78190E-11 1.77260E-11 &
     1.76340E-11 1.75430E-11 1.74530E-11 1.73640E-11 1.72750E-11 1.71880E-11 1.71020E-11 1.70160E-11 &
     1.69320E-11 1.68480E-11 1.67660E-11 1.66850E-11 1.66040E-11 1.65250E-11 1.64470E-11 1.63700E-11 &
     1.62940E-11 1.62190E-11 1.61450E-11 1.60720E-11 1.60000E-11 1.59290E-11 1.58600E-11 1.57910E-11 &
     1.57240E-11 1.56580E-11 1.55920E-11 1.55280E-11 1.54660E-11 1.54040E-11 1.53430E-11 1.52840E-11 &
     1.52250E-11 1.51680E-11 1.51120E-11 1.50570E-11 1.50030E-11 1.49500E-11 1.48980E-11 1.48480E-11 &
     1.47980E-11 1.47500E-11 1.47030E-11 1.46570E-11 1.46120E-11 1.45680E-11 1.45250E-11 1.44840E-11 &
     1.44430E-11 1.44040E-11 1.43650E-11 1.43280E-11 1.42920E-11 1.42570E-11 1.42230E-11 1.41900E-11 &
     1.41570E-11 1.41260E-11 1.40960E-11 1.40670E-11 1.40390E-11 1.40120E-11 1.39860E-11 1.39610E-11 &
     1.39370E-11 1.39140E-11 1.38920E-11 1.38710E-11 1.38500E-11 1.38310E-11 1.38120E-11 1.37940E-11 &
     1.37780E-11 1.37620E-11 1.37460E-11 1.37320E-11 1.37190E-11 1.37060E-11 1.36940E-11 1.36830E-11 &
     1.36720E-11 1.36630E-11 1.36540E-11 1.36450E-11 1.36380E-11 1.36310E-11 1.36250E-11 1.36190E-11 &
     1.36140E-11 1.36100E-11 1.36060E-11 1.36030E-11 1.36000E-11 1.35990E-11 1.35970E-11 1.35960E-11 &
     1.35960E-11 1.35960E-11 1.35960E-11 1.35970E-11 1.35990E-11 1.36000E-11 1.36030E-11 1.36050E-11 &
     1.36080E-11 1.36120E-11 1.36150E-11 1.36190E-11 1.36240E-11 1.36280E-11 1.36330E-11 1.36380E-11 &
     1.36430E-11 1.36490E-11 1.36550E-11 1.36610E-11 1.36670E-11 1.36730E-11 1.36790E-11 1.36860E-11 &
     1.36930E-11 1.36990E-11 1.37060E-11 1.37130E-11 1.37200E-11 1.37270E-11 1.37340E-11 1.37410E-11 &
     1.37470E-11 1.37540E-11 1.37610E-11 1.37680E-11 1.37750E-11 1.37810E-11 1.37880E-11 1.37940E-11 &
     1.38000E-11 1.38060E-11 1.38120E-11 1.38180E-11 1.38230E-11 1.38290E-11 1.38340E-11 1.38390E-11 &
     1.38430E-11 1.38470E-11 1.38520E-11 1.38550E-11 1.38590E-11 1.38620E-11 1.38650E-11 1.38680E-11 &
     1.38700E-11 1.38720E-11 1.38730E-11 1.38740E-11 1.38750E-11 1.38750E-11 1.38750E-11 1.38750E-11 &
     1.38740E-11 1.38730E-11 1.38710E-11 1.38690E-11 1.38660E-11 1.38630E-11 1.38600E-11 1.38560E-11 &
     1.38510E-11 1.38460E-11 1.38410E-11 1.38350E-11 1.38290E-11 1.38220E-11 1.38140E-11 1.38060E-11 &
     1.37980E-11 1.37890E-11 1.37790E-11 1.37690E-11 1.37580E-11 1.37470E-11 1.37350E-11 1.37230E-11 &
     1.37100E-11 1.36970E-11 1.36830E-11 1.36690E-11 1.36530E-11 1.36380E-11 1.36220E-11 1.36050E-11 &
     1.35880E-11 1.35700E-11 1.35520E-11 1.35330E-11 1.35130E-11 1.34930E-11 1.34730E-11 1.34510E-11 &
     1.34300E-11 1.34080E-11 1.33850E-11 1.33610E-11 1.33380E-11 1.33130E-11 1.32880E-11 1.32630E-11 &
     1.32370E-11 1.32100E-11 1.31830E-11 1.31560E-11 1.31280E-11 1.30990E-11 1.30700E-11 1.30410E-11 &
     1.30110E-11 1.29800E-11 1.29490E-11 1.29170E-11 1.28850E-11 1.28530E-11 1.28200E-11 1.27870E-11 &
     1.27530E-11 1.27190E-11 1.26840E-11 1.26490E-11 1.26130E-11 1.25770E-11 1.25410E-11 1.25040E-11 &
     1.24670E-11 1.24300E-11 1.23920E-11 1.23540E-11 1.23150E-11 1.22770E-11 1.22370E-11 1.21980E-11 &
     1.21580E-11 1.21180E-11 1.20780E-11 1.20370E-11 1.19960E-11 1.19550E-11 1.19130E-11 1.18720E-11 &
     1.18300E-11 1.17880E-11 1.17450E-11 1.17030E-11 1.16600E-11 1.16170E-11 1.15740E-11 1.15300E-11 &
     1.14870E-11 1.14430E-11 1.14000E-11 1.13560E-11 1.13120E-11 1.12680E-11 1.12240E-11 1.11790E-11 &
     1.11350E-11 1.10910E-11 1.10460E-11 1.10020E-11 1.09570E-11 1.09130E-11 1.08680E-11 1.08230E-11 &
     1.07790E-11 1.07340E-11 1.06900E-11 1.06450E-11 1.06010E-11 1.05570E-11 1.05120E-11 1.04680E-11 &
     1.04240E-11 1.03800E-11 1.03360E-11 1.02920E-11 1.02480E-11 1.02050E-11 1.01610E-11 1.01180E-11 &
     1.00750E-11 1.00320E-11 9.98940E-12 9.94690E-12 9.90450E-12 9.86240E-12 9.82040E-12 9.77870E-12 &
     9.73730E-12 9.69600E-12 9.65500E-12 9.61430E-12 9.57380E-12 9.53360E-12 9.49370E-12 9.45400E-12 &
     9.41470E-12 9.37560E-12 9.33690E-12 9.29840E-12 9.26030E-12 9.22250E-12 9.18510E-12 9.14790E-12 &
     9.11110E-12 9.07470E-12 9.03860E-12 9.00290E-12 8.96760E-12 8.93260E-12 8.89800E-12 8.86380E-12 &
     8.82990E-12 8.79650E-12 8.76340E-12 8.73070E-12 8.69850E-12 8.66660E-12 8.63510E-12 8.60410E-12 &
     8.57350E-12 8.54320E-12 8.51340E-12 8.48410E-12 8.45510E-12 8.42660E-12 8.39850E-12 8.37080E-12 &
     8.34360E-12 8.31680E-12 8.29050E-12 8.26460E-12 8.23910E-12 8.21400E-12 8.18940E-12 8.16530E-12 &
     8.14160E-12 8.11830E-12 8.09550E-12 8.07310E-12 8.05110E-12 8.02960E-12 8.00860E-12 7.98800E-12 &
     7.96780E-12 7.94800E-12 7.92870E-12 7.90990E-12 7.89140E-12 7.87340E-12 7.85590E-12 7.83870E-12 &
     7.82200E-12 7.80570E-12 7.78980E-12 7.77430E-12 7.75930E-12 7.74470E-12 7.73040E-12 7.71660E-12 &
     7.70320E-12 7.69020E-12 7.67760E-12 7.66530E-12 7.65350E-12 7.64200E-12 7.63090E-12 7.62020E-12 &
     7.60990E-12 7.59990E-12 7.59030E-12 7.58110E-12 7.57220E-12 7.56370E-12 7.55550E-12 7.54760E-12 &
     7.54010E-12 7.53290E-12 7.52600E-12 7.51940E-12 7.51310E-12 7.50720E-12 7.50150E-12 7.49620E-12 &
     7.49110E-12 7.48630E-12 7.48170E-12 7.47750E-12 7.47350E-12 7.46970E-12 7.46620E-12 7.46300E-12 &
     7.45990E-12 7.45710E-12 7.45460E-12 7.45220E-12 7.45000E-12 7.44810E-12 7.44630E-12 7.44470E-12 &
     7.44330E-12 7.44210E-12 7.44100E-12 7.44010E-12 7.43930E-12 7.43870E-12 7.43820E-12 7.43790E-12 &
     7.43770E-12 7.43760E-12 7.43760E-12 7.43770E-12 7.43790E-12 7.43820E-12 7.43860E-12 7.43900E-12 &
     7.43950E-12 7.44010E-12 7.44070E-12 7.44140E-12 7.44210E-12 7.44290E-12 7.44360E-12 7.44440E-12 &
     7.44520E-12 7.44600E-12 7.44680E-12 7.44760E-12 7.44840E-12 7.44920E-12 7.44990E-12 7.45060E-12 &
     7.45120E-12 7.45180E-12 7.45240E-12 7.45290E-12 7.45330E-12 7.45370E-12 7.45400E-12 7.45420E-12 &
     7.45430E-12 7.45430E-12 7.45420E-12 7.45400E-12 7.45370E-12 7.45330E-12 7.45270E-12 7.45200E-12 &
     7.45120E-12 7.45020E-12 7.44910E-12 7.44780E-12 7.44640E-12 7.44480E-12 7.44300E-12 7.44110E-12 &
     7.43900E-12 7.43670E-12 7.43420E-12 7.43150E-12 7.42860E-12 7.42560E-12 7.42230E-12 7.41880E-12 &
     7.41510E-12 7.41120E-12 7.40710E-12 7.40280E-12 7.39820E-12 7.39340E-12 7.38840E-12 7.38320E-12 &
     7.37770E-12 7.37200E-12 7.36600E-12 7.35980E-12 7.35340E-12 7.34670E-12 7.33980E-12 7.33260E-12 &
     7.32520E-12 7.31750E-12 7.30950E-12 7.30130E-12 7.29290E-12 7.28420E-12 7.27520E-12 7.26600E-12 &
     7.25650E-12 7.24670E-12 7.23660E-12 7.22630E-12 7.21570E-12 7.20490E-12 7.19380E-12 7.18240E-12 &
     7.17080E-12 7.15890E-12 7.14670E-12 7.13430E-12 7.12160E-12 7.10870E-12 7.09550E-12 7.08200E-12 &
     7.06820E-12 7.05420E-12 7.03990E-12 7.02540E-12 7.01050E-12 6.99550E-12 6.98020E-12 6.96460E-12 &
     6.94880E-12 6.93270E-12 6.91630E-12 6.89980E-12 6.88290E-12 6.86580E-12 6.84850E-12 6.83090E-12 &
     6.81310E-12 6.79510E-12 6.77680E-12 6.75830E-12 6.73960E-12 6.72060E-12 6.70150E-12 6.68210E-12 &
     6.66250E-12 6.64270E-12 6.62270E-12 6.60250E-12 6.58210E-12 6.56160E-12 6.54080E-12 6.51980E-12 &
     6.49860E-12 6.47730E-12 6.45580E-12 6.43410E-12 6.41230E-12 6.39030E-12 6.36810E-12 6.34570E-12 &
     6.32330E-12 6.30060E-12 6.27780E-12 6.25490E-12 6.23180E-12 6.20860E-12 6.18530E-12 6.16190E-12 &
     6.13830E-12 6.11460E-12 6.09090E-12 6.06700E-12 6.04300E-12 6.01890E-12 5.99470E-12 5.97050E-12 &
     5.94610E-12 5.92170E-12 5.89720E-12 5.87270E-12 5.84800E-12 5.82330E-12 5.79860E-12 5.77380E-12 &
     5.74900E-12 5.72410E-12 5.69920E-12 5.67420E-12 5.64930E-12 5.62430E-12 5.59930E-12 5.57430E-12 &
     5.54920E-12 5.52420E-12 5.49920E-12 5.47420E-12 5.44920E-12 5.42420E-12 5.39920E-12 5.37430E-12 &
     5.34940E-12 5.32450E-12 5.29960E-12 5.27490E-12 5.25010E-12 5.22550E-12 5.20090E-12 5.17630E-12 &
     5.15180E-12 5.12740E-12 5.10310E-12 5.07880E-12 5.05470E-12 5.03060E-12 5.00660E-12 4.98270E-12 &
     4.95890E-12 4.93520E-12 4.91170E-12 4.88820E-12 4.86490E-12 4.84160E-12 4.81850E-12 4.79550E-12 &
     4.77270E-12 4.74990E-12 4.72740E-12 4.70490E-12 4.68260E-12 4.66040E-12 4.63840E-12 4.61660E-12 &
     4.59490E-12 4.57330E-12 4.55200E-12 4.53070E-12 4.50970E-12 4.48880E-12 4.46810E-12 4.44760E-12 &
     4.42720E-12 4.40710E-12 4.38710E-12 4.36730E-12 4.34760E-12 4.32820E-12 4.30900E-12 4.28990E-12 &
     4.27110E-12 4.25240E-12 4.23390E-12 4.21570E-12 4.19760E-12 4.17970E-12 4.16210E-12 4.14460E-12 &
     4.12740E-12 4.11040E-12 4.09350E-12 4.07690E-12 4.06050E-12 4.04430E-12 4.02830E-12 4.01260E-12 &
     3.99700E-12 3.98170E-12 3.96660E-12 3.95170E-12 3.93700E-12 3.92250E-12 3.90830E-12 3.89420E-12 &
     3.88040E-12 3.86680E-12 3.85340E-12 3.84030E-12 3.82730E-12 3.81460E-12 3.80210E-12 3.78980E-12 &
     3.77770E-12 3.76580E-12 3.75410E-12 3.74270E-12 3.73140E-12 3.72040E-12 3.70960E-12 3.69900E-12 &
     3.68860E-12 3.67840E-12 3.66840E-12 3.65860E-12 3.64900E-12 3.63960E-12 3.63040E-12 3.62150E-12 &
     3.61270E-12 3.60410E-12 3.59570E-12 3.58750E-12 3.57940E-12 3.57160E-12 3.56400E-12 3.55650E-12 &
     3.54920E-12 3.54210E-12 3.53520E-12 3.52850E-12 3.52190E-12 3.51550E-12 3.50930E-12 3.50320E-12 &
     3.49730E-12 3.49160E-12 3.48600E-12 3.48060E-12 3.47530E-12 3.47020E-12 3.46520E-12 3.46040E-12 &
     3.45580E-12 3.45120E-12 3.44680E-12 3.44260E-12 3.43850E-12 3.43450E-12 3.43060E-12 3.42690E-12 &
     3.42330E-12 3.41980E-12 3.41640E-12 3.41310E-12 3.41000E-12 3.40690E-12 3.40400E-12 3.40120E-12 &
     3.39840E-12 3.39580E-12 3.39330E-12 3.39080E-12 3.38840E-12 3.38610E-12 3.38390E-12 3.38180E-12 &
     3.37980E-12 3.37780E-12 3.37590E-12 3.37400E-12 3.37220E-12 3.37050E-12 3.36880E-12 3.36720E-12 &
     3.36570E-12 3.36420E-12 3.36270E-12 3.36130E-12 3.35990E-12 3.35850E-12 3.35720E-12 3.35590E-12 &
     3.35470E-12 3.35340E-12 3.35220E-12 3.35100E-12 3.34980E-12 3.34870E-12 3.34750E-12 3.34630E-12 &
     3.34520E-12 3.34410E-12 3.34290E-12 3.34180E-12 3.34060E-12 3.33940E-12 3.33830E-12 3.33710E-12 &
     3.33590E-12 3.33460E-12 3.33340E-12 3.33210E-12 3.33080E-12 3.32950E-12 3.32810E-12 3.32670E-12 &
     3.32530E-12 3.32380E-12 3.32230E-12 3.32080E-12 3.31920E-12 3.31750E-12 3.31590E-12 3.31410E-12 &
     3.31230E-12 3.31050E-12 3.30860E-12 3.30660E-12 3.30460E-12 3.30260E-12 3.30040E-12 3.29820E-12 &
     3.29600E-12 3.29360E-12 3.29120E-12 3.28880E-12 3.28620E-12 3.28360E-12 3.28090E-12 3.27810E-12 &
     3.27530E-12 3.27240E-12 3.26940E-12 3.26630E-12 3.26310E-12 3.25990E-12 3.25660E-12 3.25320E-12 &
     3.24970E-12 3.24610E-12 3.24240E-12 3.23870E-12 3.23480E-12 3.23090E-12 3.22690E-12 3.22280E-12 &
     3.21860E-12 3.21430E-12 3.20990E-12 3.20550E-12 3.20090E-12 3.19630E-12 3.19150E-12 3.18670E-12 &
     3.18180E-12 3.17680E-12 3.17160E-12 3.16640E-12 3.16110E-12 3.15570E-12 3.15030E-12 3.14470E-12 &
     3.13900E-12 3.13330E-12 3.12740E-12 3.12150E-12 3.11550E-12 3.10940E-12 3.10320E-12 3.09690E-12 &
     3.09050E-12 3.08410E-12 3.07750E-12 3.07090E-12 3.06420E-12 3.05740E-12 3.05050E-12 3.04360E-12 &
     3.03650E-12 3.02940E-12 3.02220E-12 3.01490E-12 3.00760E-12 3.00020E-12 2.99270E-12 2.98510E-12 &
     2.97740E-12 2.96970E-12 2.96190E-12 2.95410E-12 2.94620E-12 2.93820E-12 2.93010E-12 2.92200E-12 &
     2.91380E-12 2.90550E-12 2.89720E-12 2.88890E-12 2.88050E-12 2.87200E-12 2.86350E-12 2.85490E-12 &
     2.84630E-12 2.83760E-12 2.82890E-12 2.82010E-12 2.81130E-12 2.80240E-12 2.79350E-12 2.78460E-12 &
     2.77560E-12 2.76660E-12 2.75750E-12 2.74850E-12 2.73930E-12 2.73020E-12 2.72100E-12 2.71180E-12 &
     2.70260E-12 2.69340E-12 2.68410E-12 2.67480E-12 2.66550E-12 2.65620E-12 2.64690E-12 2.63750E-12 &
     2.62820E-12 2.61880E-12 2.60950E-12 2.60010E-12 2.59070E-12 2.58130E-12 2.57200E-12 2.56260E-12 &
     2.55320E-12 2.54390E-12 2.53450E-12 2.52520E-12 2.51580E-12 2.50650E-12 2.49720E-12 2.48790E-12 &
     2.47860E-12 2.46940E-12 2.46010E-12 2.45090E-12 2.44170E-12 2.43260E-12 2.42340E-12 2.41430E-12 &
     2.40520E-12 2.39620E-12 2.38710E-12 2.37820E-12 2.36920E-12 2.36030E-12 2.35140E-12 2.34260E-12 &
     2.33380E-12 2.32500E-12 2.31630E-12 2.30770E-12 2.29910E-12 2.29050E-12 2.28200E-12 2.27350E-12 &
     2.26510E-12 2.25670E-12 2.24840E-12 2.24020E-12 2.23200E-12 2.22380E-12 2.21580E-12 2.20770E-12 &
     2.19980E-12 2.19190E-12 2.18400E-12 2.17630E-12 2.16860E-12 2.16090E-12 2.15330E-12 2.14580E-12 &
     2.13840E-12 2.13100E-12 2.12370E-12 2.11640E-12 2.10930E-12 2.10220E-12 2.09510E-12 2.08820E-12 &
     2.08130E-12 2.07450E-12 2.06780E-12 2.06110E-12 2.05450E-12 2.04800E-12 2.04160E-12 2.03520E-12 &
     2.02890E-12 2.02270E-12 2.01660E-12 2.01050E-12 2.00460E-12 1.99870E-12 1.99290E-12 1.98710E-12 &
     1.98140E-12 1.97590E-12 1.97030E-12 1.96490E-12 1.95960E-12 1.95430E-12 1.94910E-12 1.94400E-12 &
     1.93890E-12 1.93390E-12 1.92900E-12 1.92420E-12 1.91950E-12 1.91480E-12 1.91020E-12 1.90570E-12 &
     1.90130E-12 1.89690E-12 1.89260E-12 1.88840E-12 1.88430E-12 1.88020E-12 1.87620E-12 1.87230E-12 &
     1.86850E-12 1.86470E-12 1.86100E-12 1.85730E-12 1.85380E-12 1.85030E-12 1.84680E-12 1.84340E-12 &
     1.84010E-12 1.83690E-12 1.83370E-12 1.83060E-12 1.82760E-12 1.82460E-12 1.82160E-12 1.81880E-12 &
     1.81600E-12 1.81320E-12 1.81050E-12 1.80790E-12 1.80530E-12 1.80280E-12 1.80030E-12 1.79790E-12 &
     1.79560E-12 1.79320E-12 1.79100E-12 1.78880E-12 1.78660E-12 1.78450E-12 1.78240E-12 1.78040E-12 &
     1.77840E-12 1.77640E-12 1.77450E-12 1.77270E-12 1.77080E-12 1.76910E-12 1.76730E-12 1.76560E-12 &
     1.76390E-12 1.76230E-12 1.76070E-12 1.75910E-12 1.75750E-12 1.75600E-12 1.75450E-12 1.75300E-12 &
     1.75160E-12 1.75010E-12 1.74870E-12 1.74740E-12 1.74600E-12 1.74470E-12 1.74330E-12 1.74200E-12 &
     1.74080E-12 1.73950E-12 1.73820E-12 1.73700E-12 1.73570E-12 1.73450E-12 1.73330E-12 1.73210E-12 &
     1.73080E-12 1.72960E-12 1.72840E-12 1.72720E-12 1.72600E-12 1.72480E-12 1.72360E-12 1.72240E-12 &
     1.72120E-12 1.72000E-12 1.71880E-12 1.71760E-12 1.71640E-12 1.71510E-12 1.71390E-12 1.71260E-12 &
     1.71130E-12 1.71010E-12 1.70880E-12 1.70750E-12 1.70610E-12 1.70480E-12 1.70340E-12 1.70210E-12 &
     1.70070E-12 1.69930E-12 1.69780E-12 1.69640E-12 1.69490E-12 1.69340E-12 1.69190E-12 1.69030E-12 &
     1.68880E-12 1.68710E-12 1.68550E-12 1.68390E-12 1.68220E-12 1.68050E-12 1.67870E-12 1.67690E-12 &
     1.67510E-12 1.67330E-12 1.67140E-12 1.66950E-12 1.66760E-12 1.66560E-12 1.66360E-12 1.66160E-12 &
     1.65950E-12 1.65740E-12 1.65520E-12 1.65300E-12 1.65080E-12 1.64860E-12 1.64630E-12 1.64390E-12 &
     1.64150E-12 1.63910E-12 1.63660E-12 1.63410E-12 1.63160E-12 1.62900E-12 1.62640E-12 1.62370E-12 &
     1.62100E-12 1.61830E-12 1.61550E-12 1.61270E-12 1.60980E-12 1.60690E-12 1.60390E-12 1.60090E-12 &
     1.59780E-12 1.59480E-12 1.59160E-12 1.58850E-12 1.58520E-12 1.58200E-12 1.57870E-12 1.57530E-12 &
     1.57200E-12 1.56860E-12 1.56510E-12 1.56160E-12 1.55800E-12 1.55450E-12 1.55080E-12 1.54720E-12 &
     1.54350E-12 1.53970E-12 1.53590E-12 1.53210E-12 1.52830E-12 1.52440E-12 1.52040E-12 1.51650E-12 &
     1.51240E-12 1.50840E-12 1.50430E-12 1.50020E-12 1.49600E-12 1.49180E-12 1.48760E-12 1.48330E-12 &
     1.47900E-12 1.47470E-12 1.47030E-12 1.46590E-12 1.46150E-12 1.45710E-12 1.45260E-12 1.44800E-12 &
     1.44350E-12 1.43890E-12 1.43430E-12 1.42970E-12 1.42500E-12 1.42030E-12 1.41560E-12 1.41090E-12 &
     1.40610E-12 1.40130E-12 1.39650E-12 1.39170E-12 1.38680E-12 1.38190E-12 1.37700E-12 1.37210E-12 &
     1.36720E-12 1.36220E-12 1.35720E-12 1.35220E-12 1.34720E-12 1.34220E-12 1.33720E-12 1.33210E-12 &
     1.32710E-12 1.32200E-12 1.31690E-12 1.31180E-12 1.30670E-12 1.30160E-12 1.29650E-12 1.29140E-12 &
     1.28620E-12 1.28110E-12 1.27590E-12 1.27080E-12 1.26570E-12 1.26050E-12 1.25530E-12 1.25020E-12 &
     1.24500E-12 1.23990E-12 1.23470E-12 1.22960E-12 1.22440E-12 1.21930E-12 1.21410E-12 1.20900E-12 &
     1.20390E-12 1.19870E-12 1.19360E-12 1.18850E-12 1.18340E-12 1.17830E-12 1.17330E-12 1.16820E-12 &
     1.16310E-12 1.15810E-12 1.15310E-12 1.14800E-12 1.14300E-12 1.13810E-12 1.13310E-12 1.12810E-12 &
     1.12320E-12 1.11830E-12 1.11340E-12 1.10850E-12 1.10360E-12 1.09880E-12 1.09400E-12 1.08920E-12 &
     1.08440E-12 1.07970E-12 1.07490E-12 1.07020E-12 1.06550E-12 1.06090E-12 1.05630E-12 1.05170E-12 &
     1.04710E-12 1.04250E-12 1.03800E-12 1.03350E-12 1.02910E-12 1.02460E-12 1.02020E-12 1.01590E-12 &
     1.01150E-12 1.00720E-12 1.00300E-12 9.98710E-13 9.94500E-13 9.90320E-13 9.86170E-13 9.82060E-13 &
     9.77990E-13 9.73940E-13 9.69940E-13 9.65960E-13 9.62020E-13 9.58120E-13 9.54250E-13 9.50420E-13 &
     9.46630E-13 9.42870E-13 9.39140E-13 9.35460E-13 9.31810E-13 9.28200E-13 9.24620E-13 9.21080E-13 &
     9.17580E-13 9.14120E-13 9.10700E-13 9.07310E-13 9.03970E-13 9.00660E-13 8.97390E-13 8.94150E-13 &
     8.90960E-13 8.87800E-13 8.84680E-13 8.81600E-13 8.78560E-13 8.75560E-13 8.72590E-13 8.69670E-13 &
     8.66780E-13 8.63930E-13 8.61120E-13 8.58340E-13 8.55610E-13 8.52910E-13 8.50250E-13 8.47630E-13 &
     8.45040E-13 8.42490E-13 8.39980E-13 8.37510E-13 8.35070E-13 8.32670E-13 8.30300E-13 8.27980E-13 &
     8.25680E-13 8.23430E-13 8.21210E-13 8.19020E-13 8.16870E-13 8.14760E-13 8.12680E-13 8.10630E-13 &
     8.08610E-13 8.06630E-13 8.04690E-13 8.02770E-13 8.00890E-13 7.99040E-13 7.97230E-13 7.95440E-13 &
     7.93690E-13 7.91960E-13 7.90270E-13 7.88610E-13 7.86970E-13 7.85370E-13 7.83790E-13 7.82250E-13 &
     7.80730E-13 7.79230E-13 7.77770E-13 7.76330E-13 7.74920E-13 7.73530E-13 7.72170E-13 7.70830E-13 &
     7.69520E-13 7.68230E-13 7.66970E-13 7.65730E-13 7.64510E-13 7.63320E-13 7.62140E-13 7.60990E-13 &
     7.59860E-13 7.58750E-13 7.57660E-13 7.56580E-13 7.55530E-13 7.54490E-13 7.53480E-13 7.52480E-13 &
     7.51490E-13 7.50530E-13 7.49580E-13 7.48640E-13 7.47720E-13 7.46810E-13 7.45920E-13 7.45040E-13 &
     7.44170E-13 7.43320E-13 7.42480E-13 7.41650E-13 7.40830E-13 7.40020E-13 7.39220E-13 7.38430E-13 &
     7.37650E-13 7.36880E-13 7.36120E-13 7.35360E-13 7.34610E-13 7.33870E-13 7.33130E-13 7.32400E-13 &
     7.31680E-13 7.30960E-13 7.30240E-13 7.29530E-13 7.28820E-13 7.28120E-13 7.27410E-13 7.26710E-13 &
     7.26010E-13 7.25310E-13 7.24620E-13 7.23920E-13 7.23220E-13 7.22530E-13 7.21830E-13 7.21130E-13 &
     7.20430E-13 7.19720E-13 7.19020E-13 7.18310E-13 7.17600E-13 7.16890E-13 7.16170E-13 7.15450E-13 &
     7.14720E-13 7.13990E-13 7.13250E-13 7.12510E-13 7.11760E-13 7.11000E-13 7.10240E-13 7.09470E-13 &
     7.08700E-13 7.07920E-13 7.07130E-13 7.06330E-13 7.05520E-13 7.04710E-13 7.03880E-13 7.03050E-13 &
     7.02210E-13 7.01360E-13 7.00500E-13 6.99630E-13 6.98740E-13 6.97850E-13 6.96950E-13 6.96040E-13 &
     6.95110E-13 6.94180E-13 6.93230E-13 6.92280E-13 6.91310E-13 6.90330E-13 6.89340E-13 6.88330E-13 &
     6.87320E-13 6.86290E-13 6.85250E-13 6.84200E-13 6.83140E-13 6.82060E-13 6.80970E-13 6.79870E-13 &
     6.78760E-13 6.77630E-13 6.76490E-13 6.75340E-13 6.74170E-13 6.72990E-13 6.71800E-13 6.70600E-13 &
     6.69380E-13 6.68150E-13 6.66910E-13 6.65660E-13 6.64390E-13 6.63110E-13 6.61810E-13 6.60510E-13 &
     6.59190E-13 6.57860E-13 6.56510E-13 6.55160E-13 6.53790E-13 6.52410E-13 6.51010E-13 6.49610E-13 &
     6.48190E-13 6.46760E-13 6.45320E-13 6.43870E-13 6.42400E-13 6.40930E-13 6.39440E-13 6.37940E-13 &
     6.36430E-13 6.34910E-13 6.33380E-13 6.31840E-13 6.30280E-13 6.28720E-13 6.27150E-13 6.25560E-13 &
     6.23970E-13 6.22370E-13 6.20760E-13 6.19140E-13 6.17510E-13 6.15870E-13 6.14220E-13 6.12570E-13 &
     6.10910E-13 6.09230E-13 6.07560E-13 6.05870E-13 6.04180E-13 6.02480E-13 6.00770E-13 5.99060E-13 &
     5.97340E-13 5.95610E-13 5.93880E-13 5.92150E-13 5.90400E-13 5.88660E-13 5.86910E-13 5.85150E-13 &
     5.83390E-13 5.81620E-13 5.79850E-13 5.78080E-13 5.76300E-13 5.74530E-13 5.72740E-13 5.70960E-13 &
     5.69170E-13 5.67390E-13 5.65600E-13 5.63800E-13 5.62010E-13 5.60220E-13 5.58420E-13 5.56630E-13 &
     5.54830E-13 5.53040E-13 5.51240E-13 5.49450E-13 5.47650E-13 5.45860E-13 5.44070E-13 5.42280E-13 &
     5.40490E-13 5.38710E-13 5.36920E-13 5.35140E-13 5.33360E-13 5.31590E-13 5.29820E-13 5.28050E-13 &
     5.26280E-13 5.24520E-13 5.22770E-13 5.21010E-13 5.19270E-13 5.17530E-13 5.15790E-13 5.14060E-13 &
     5.12330E-13 5.10610E-13 5.08900E-13 5.07190E-13 5.05490E-13 5.03790E-13 5.02110E-13 5.00430E-13 &
     4.98750E-13 4.97090E-13 4.95430E-13 4.93780E-13 4.92140E-13 4.90510E-13 4.88880E-13 4.87270E-13 &
     4.85660E-13 4.84060E-13 4.82470E-13 4.80890E-13 4.79320E-13 4.77760E-13 4.76210E-13 4.74670E-13 &
     4.73130E-13 4.71610E-13 4.70100E-13 4.68600E-13 4.67110E-13 4.65630E-13 4.64160E-13 4.62700E-13 &
     4.61260E-13 4.59820E-13 4.58400E-13 4.56980E-13 4.55580E-13 4.54190E-13 4.52810E-13 4.51440E-13 &
     4.50090E-13 4.48740E-13 4.47410E-13 4.46090E-13 4.44780E-13 4.43480E-13 4.42200E-13 4.40930E-13 &
     4.39660E-13 4.38420E-13 4.37180E-13 4.35950E-13 4.34740E-13 4.33540E-13 4.32350E-13 4.31170E-13 &
     4.30010E-13 4.28860E-13 4.27720E-13 4.26590E-13 4.25470E-13 4.24370E-13 4.23280E-13 4.22200E-13 &
     4.21130E-13 4.20070E-13 4.19030E-13 4.18000E-13 4.16980E-13 4.15970E-13 4.14970E-13 4.13990E-13 &
     4.13010E-13 4.12050E-13 4.11100E-13 4.10160E-13 4.09240E-13 4.08320E-13 4.07410E-13 4.06520E-13 &
     4.05640E-13 4.04760E-13 4.03900E-13 4.03050E-13 4.02210E-13 4.01380E-13 4.00560E-13 3.99750E-13 &
     3.98950E-13 3.98160E-13 3.97380E-13 3.96610E-13 3.95850E-13 3.95100E-13 3.94360E-13 3.93630E-13 &
     3.92910E-13 3.92190E-13 3.91490E-13 3.90790E-13 3.90100E-13 3.89420E-13 3.88750E-13 3.88080E-13 &
     3.87430E-13 3.86780E-13 3.86140E-13 3.85510E-13 3.84880E-13 3.84260E-13 3.83650E-13 3.83040E-13 &
     3.82440E-13 3.81850E-13 3.81260E-13 3.80680E-13 3.80110E-13 3.79540E-13 3.78970E-13 3.78410E-13 &
     3.77860E-13 3.77310E-13 3.76770E-13 3.76230E-13 3.75690E-13 3.75160E-13 3.74640E-13 3.74110E-13 &
     3.73590E-13 3.73080E-13 3.72570E-13 3.72060E-13 3.71550E-13 3.71050E-13 3.70550E-13 3.70060E-13 &
     3.69560E-13 3.69070E-13 3.68580E-13 3.68090E-13 3.67600E-13 3.67120E-13 3.66630E-13 3.66150E-13 &
     3.65660E-13 3.65180E-13 3.64700E-13 3.64220E-13 3.63740E-13 3.63260E-13 3.62780E-13 3.62300E-13 &
     3.61820E-13 3.61340E-13 3.60860E-13 3.60380E-13 3.59900E-13 3.59410E-13 3.58930E-13 3.58440E-13 &
     3.57950E-13 3.57460E-13 3.56970E-13 3.56470E-13 3.55980E-13 3.55480E-13 3.54980E-13 3.54470E-13 &
     3.53970E-13 3.53460E-13 3.52940E-13 3.52430E-13 3.51910E-13 3.51380E-13 3.50860E-13 3.50330E-13 &
     3.49790E-13 3.49260E-13 3.48710E-13 3.48170E-13 3.47620E-13 3.47070E-13 3.46510E-13 3.45950E-13 &
     3.45380E-13 3.44810E-13 3.44240E-13 3.43660E-13 3.43070E-13 3.42480E-13 3.41890E-13 3.41290E-13 &
     3.40690E-13 3.40080E-13 3.39460E-13 3.38840E-13 3.38220E-13 3.37590E-13 3.36950E-13 3.36310E-13 &
     3.35670E-13 3.35020E-13 3.34360E-13 3.33700E-13 3.33030E-13 3.32350E-13 3.31670E-13 3.30990E-13 &
     3.30300E-13 3.29600E-13 3.28900E-13 3.28190E-13 3.27480E-13 3.26760E-13 3.26030E-13 3.25300E-13 &
     3.24560E-13 3.23820E-13 3.23070E-13 3.22310E-13 3.21550E-13 3.20790E-13 3.20010E-13 3.19240E-13 &
     3.18450E-13 3.17660E-13 3.16870E-13 3.16060E-13 3.15260E-13 3.14450E-13 3.13630E-13 3.12810E-13 &
     3.11980E-13 3.11140E-13 3.10300E-13 3.09460E-13 3.08610E-13 3.07750E-13 3.06890E-13 3.06030E-13 &
     3.05160E-13 3.04280E-13 3.03400E-13 3.02520E-13 3.01630E-13 3.00740E-13 2.99840E-13 2.98940E-13 &
     2.98030E-13 2.97120E-13 2.96210E-13 2.95290E-13 2.94360E-13 2.93440E-13 2.92500E-13 2.91570E-13 &
     2.90630E-13 2.89690E-13 2.88740E-13 2.87800E-13 2.86840E-13 2.85890E-13 2.84930E-13 2.83970E-13 &
     2.83000E-13 2.82030E-13 2.81060E-13 2.80090E-13 2.79110E-13 2.78140E-13 2.77160E-13 2.76180E-13 &
     2.75190E-13 2.74210E-13 2.73220E-13 2.72230E-13 2.71240E-13 2.70240E-13 2.69250E-13 2.68250E-13 &
     2.67260E-13 2.66260E-13 2.65260E-13 2.64260E-13 2.63260E-13 2.62260E-13 2.61260E-13 2.60260E-13 &
     2.59260E-13 2.58260E-13 2.57250E-13 2.56250E-13 2.55250E-13 2.54250E-13 2.53250E-13 2.52250E-13 &
     2.51250E-13 2.50250E-13 2.49250E-13 2.48260E-13 2.47260E-13 2.46270E-13 2.45280E-13 2.44290E-13 &
     2.43300E-13 2.42310E-13 2.41320E-13 2.40340E-13 2.39360E-13 2.38380E-13 2.37410E-13 2.36430E-13 &
     2.35460E-13 2.34490E-13 2.33530E-13 2.32570E-13 2.31610E-13 2.30650E-13 2.29700E-13 2.28750E-13 &
     2.27800E-13 2.26860E-13 2.25920E-13 2.24980E-13 2.24050E-13 2.23120E-13 2.22200E-13 2.21280E-13 &
     2.20360E-13 2.19450E-13 2.18540E-13 2.17640E-13 2.16740E-13 2.15850E-13 2.14960E-13 2.14070E-13 &
     2.13190E-13 2.12310E-13 2.11440E-13 2.10580E-13 2.09720E-13 2.08860E-13 2.08010E-13 2.07170E-13 &
     2.06330E-13 2.05490E-13 2.04660E-13 2.03840E-13 2.03020E-13 2.02210E-13 2.01410E-13 2.00600E-13 &
     1.99810E-13 1.99020E-13 1.98240E-13 1.97460E-13 1.96690E-13 1.95930E-13 1.95170E-13 1.94410E-13 &
     1.93670E-13 1.92930E-13 1.92190E-13 1.91470E-13 1.90740E-13 1.90030E-13 1.89320E-13 1.88620E-13 &
     1.87920E-13 1.87230E-13 1.86550E-13 1.85870E-13 1.85200E-13 1.84540E-13 1.83880E-13 1.83230E-13 &
     1.82580E-13 1.81940E-13 1.81310E-13 1.80680E-13 1.80060E-13 1.79450E-13 1.78850E-13 1.78240E-13 &
     1.77650E-13 1.77060E-13 1.76480E-13 1.75910E-13 1.75340E-13 1.74780E-13 1.74220E-13 1.73670E-13 &
     1.73130E-13 1.72590E-13 1.72060E-13 1.71540E-13 1.71020E-13 1.70510E-13 1.70000E-13 1.69500E-13 &
     1.69010E-13 1.68520E-13 1.68040E-13 1.67560E-13 1.67090E-13 1.66630E-13 1.66170E-13 1.65720E-13 &
     1.65270E-13 1.64830E-13 1.64400E-13 1.63970E-13 1.63540E-13 1.63120E-13 1.62710E-13 1.62300E-13 &
     1.61900E-13 1.61500E-13 1.61110E-13 1.60730E-13 1.60340E-13 1.59970E-13 1.59600E-13 1.59230E-13 &
     1.58870E-13 1.58510E-13 1.58160E-13 1.57810E-13 1.57470E-13 1.57130E-13 1.56790E-13 1.56460E-13 &
     1.56140E-13 1.55820E-13 1.55500E-13 1.55190E-13 1.54880E-13 1.54580E-13 1.54280E-13 1.53980E-13 &
     1.53690E-13 1.53400E-13 1.53110E-13 1.52830E-13 1.52550E-13 1.52280E-13 1.52000E-13 1.51730E-13 &
     1.51470E-13 1.51210E-13 1.50950E-13 1.50690E-13 1.50440E-13 1.50190E-13 1.49940E-13 1.49700E-13 &
     1.49450E-13 1.49210E-13 1.48970E-13 1.48740E-13 1.48510E-13 1.48280E-13 1.48050E-13 1.47820E-13 &
     1.47600E-13 1.47370E-13 1.47150E-13 1.46940E-13 1.46720E-13 1.46500E-13 1.46290E-13 1.46080E-13 &
     1.45870E-13 1.45660E-13 1.45450E-13 1.45240E-13 1.45040E-13 1.44830E-13 1.44630E-13 1.44420E-13 &
     1.44220E-13 1.44020E-13 1.43820E-13 1.43620E-13 1.43420E-13 1.43230E-13 1.43030E-13 1.42830E-13 &
     1.42630E-13 1.42440E-13 1.42240E-13 1.42050E-13 1.41850E-13 1.41650E-13 1.41460E-13 1.41260E-13 &
     1.41070E-13 1.40870E-13 1.40670E-13 1.40480E-13 1.40280E-13 1.40080E-13 1.39890E-13 1.39690E-13 &
     1.39490E-13 1.39290E-13 1.39090E-13 1.38890E-13 1.38690E-13 1.38490E-13 1.38290E-13 1.38090E-13 &
     1.37880E-13 1.37680E-13 1.37470E-13 1.37270E-13 1.37060E-13 1.36850E-13 1.36640E-13 1.36430E-13 &
     1.36220E-13 1.36010E-13 1.35800E-13 1.35580E-13 1.35370E-13 1.35150E-13 1.34930E-13 1.34710E-13 &
     1.34490E-13 1.34270E-13 1.34050E-13 1.33820E-13 1.33600E-13 1.33370E-13 1.33140E-13 1.32910E-13 &
     1.32680E-13 1.32440E-13 1.32210E-13 1.31970E-13 1.31740E-13 1.31500E-13 1.31260E-13 1.31020E-13 &
     1.30770E-13 1.30530E-13 1.30280E-13 1.30030E-13 1.29780E-13 1.29530E-13 1.29280E-13 1.29030E-13 &
     1.28770E-13 1.28510E-13 1.28250E-13 1.27990E-13 1.27730E-13 1.27470E-13 1.27210E-13 1.26940E-13 &
     1.26670E-13 1.26400E-13 1.26130E-13 1.25860E-13 1.25590E-13 1.25310E-13 1.25040E-13 1.24760E-13 &
     1.24480E-13 1.24200E-13 1.23920E-13 1.23640E-13 1.23350E-13 1.23070E-13 1.22780E-13 1.22500E-13 &
     1.22210E-13 1.21920E-13 1.21630E-13 1.21330E-13 1.21040E-13 1.20750E-13 1.20450E-13 1.20150E-13 &
     1.19860E-13 1.19560E-13 1.19260E-13 1.18960E-13 1.18660E-13 1.18350E-13 1.18050E-13 1.17750E-13 &
     1.17440E-13 1.17140E-13 1.16830E-13 1.16530E-13 1.16220E-13 1.15910E-13 1.15600E-13 1.15290E-13 &
     1.14980E-13 1.14670E-13 1.14360E-13 1.14050E-13 1.13740E-13 1.13430E-13 1.13110E-13 1.12800E-13 &
     1.12490E-13 1.12170E-13 1.11860E-13 1.11550E-13 1.11230E-13 1.10920E-13 1.10610E-13 1.10290E-13 &
     1.09980E-13 1.09660E-13 1.09350E-13 1.09040E-13 1.08720E-13 1.08410E-13 1.08090E-13 1.07780E-13 &
     1.07470E-13 1.07160E-13 1.06840E-13 1.06530E-13 1.06220E-13 1.05910E-13 1.05600E-13 1.05290E-13 &
     1.04980E-13 1.04670E-13 1.04360E-13 1.04050E-13 1.03750E-13 1.03440E-13 1.03140E-13 1.02830E-13 &
     1.02530E-13 1.02220E-13 1.01920E-13 1.01620E-13 1.01320E-13 1.01020E-13 1.00720E-13 1.00420E-13 &
     1.00120E-13 9.98270E-14 9.95310E-14 9.92370E-14 9.89440E-14 9.86520E-14 9.83620E-14 9.80720E-14 &
     9.77840E-14 9.74960E-14 9.72100E-14 9.69260E-14 9.66420E-14 9.63600E-14 9.60790E-14 9.58000E-14 &
     9.55210E-14 9.52450E-14 9.49690E-14 9.46950E-14 9.44230E-14 9.41520E-14 9.38820E-14 9.36140E-14 &
     9.33470E-14 9.30820E-14 9.28180E-14 9.25560E-14 9.22960E-14 9.20370E-14 9.17790E-14 9.15240E-14 &
     9.12690E-14 9.10170E-14 9.07660E-14 9.05160E-14 9.02680E-14 9.00220E-14 8.97770E-14 8.95340E-14 &
     8.92930E-14 8.90530E-14 8.88150E-14 8.85790E-14 8.83440E-14 8.81110E-14 8.78800E-14 8.76500E-14 &
     8.74220E-14 8.71960E-14 8.69710E-14 8.67480E-14 8.65260E-14 8.63070E-14 8.60890E-14 8.58720E-14 &
     8.56570E-14 8.54440E-14 8.52330E-14 8.50230E-14 8.48150E-14 8.46090E-14 8.44040E-14 8.42010E-14 &
     8.40000E-14 8.38000E-14 8.36020E-14 8.34050E-14 8.32090E-14 8.30160E-14 8.28230E-14 8.26330E-14 &
     8.24430E-14 8.22560E-14 8.20700E-14 8.18850E-14 8.17020E-14 8.15200E-14 8.13390E-14 8.11610E-14 &
     8.09830E-14 8.08070E-14 8.06320E-14 8.04590E-14 8.02870E-14 8.01160E-14 7.99470E-14 7.97780E-14 &
     7.96110E-14 7.94460E-14 7.92810E-14 7.91180E-14 7.89560E-14 7.87950E-14 7.86360E-14 7.84770E-14 &
     7.83200E-14 7.81640E-14 7.80090E-14 7.78550E-14 7.77030E-14 7.75510E-14 7.74010E-14 7.72510E-14 &
     7.71020E-14 7.69550E-14 7.68080E-14 7.66620E-14 7.65170E-14 7.63730E-14 7.62300E-14 7.60870E-14 &
     7.59460E-14 7.58050E-14 7.56640E-14 7.55250E-14 7.53850E-14 7.52470E-14 7.51090E-14 7.49710E-14 &
     7.48350E-14 7.46980E-14 7.45630E-14 7.44270E-14 7.42930E-14 7.41590E-14 7.40250E-14 7.38920E-14 &
     7.37590E-14 7.36260E-14 7.34940E-14 7.33630E-14 7.32310E-14 7.31000E-14 7.29690E-14 7.28380E-14 &
     7.27080E-14 7.25780E-14 7.24480E-14 7.23180E-14 7.21880E-14 7.20590E-14 7.19290E-14 7.18000E-14 &
     7.16710E-14 7.15420E-14 7.14120E-14 7.12830E-14 7.11540E-14 7.10250E-14 7.08960E-14 7.07660E-14 &
     7.06370E-14 7.05070E-14 7.03770E-14 7.02470E-14 7.01170E-14 6.99870E-14 6.98560E-14 6.97250E-14 &
     6.95940E-14 6.94630E-14 6.93310E-14 6.91990E-14 6.90670E-14 6.89340E-14 6.88010E-14 6.86670E-14 &
     6.85330E-14 6.83990E-14 6.82640E-14 6.81290E-14 6.79940E-14 6.78580E-14 6.77210E-14 6.75840E-14 &
     6.74460E-14 6.73080E-14 6.71700E-14 6.70310E-14 6.68910E-14 6.67510E-14 6.66110E-14 6.64700E-14 &
     6.63280E-14 6.61860E-14 6.60430E-14 6.58990E-14 6.57550E-14 6.56100E-14 6.54650E-14 6.53180E-14 &
     6.51720E-14 6.50240E-14 6.48760E-14 6.47270E-14 6.45780E-14 6.44280E-14 6.42770E-14 6.41260E-14 &
     6.39740E-14 6.38210E-14 6.36680E-14 6.35140E-14 6.33590E-14 6.32040E-14 6.30470E-14 6.28910E-14 &
     6.27330E-14 6.25750E-14 6.24160E-14 6.22560E-14 6.20960E-14 6.19350E-14 6.17730E-14 6.16110E-14 &
     6.14480E-14 6.12850E-14 6.11210E-14 6.09560E-14 6.07910E-14 6.06250E-14 6.04590E-14 6.02910E-14 &
     6.01230E-14 5.99550E-14 5.97860E-14 5.96160E-14 5.94460E-14 5.92750E-14 5.91040E-14 5.89320E-14 &
     5.87590E-14 5.85860E-14 5.84120E-14 5.82370E-14 5.80620E-14 5.78860E-14 5.77100E-14 5.75330E-14 &
     5.73560E-14 5.71780E-14 5.70000E-14 5.68210E-14 5.66420E-14 5.64620E-14 5.62830E-14 5.61020E-14 &
     5.59220E-14 5.57410E-14 5.55590E-14 5.53770E-14 5.51950E-14 5.50130E-14 5.48300E-14 5.46470E-14 &
     5.44640E-14 5.42800E-14 5.40960E-14 5.39120E-14 5.37280E-14 5.35430E-14 5.33590E-14 5.31740E-14 &
     5.29890E-14 5.28040E-14 5.26190E-14 5.24340E-14 5.22480E-14 5.20630E-14 5.18770E-14 5.16920E-14 &
     5.15060E-14 5.13200E-14 5.11340E-14 5.09480E-14 5.07630E-14 5.05770E-14 5.03910E-14 5.02060E-14 &
     5.00200E-14 4.98340E-14 4.96490E-14 4.94640E-14 4.92790E-14 4.90940E-14 4.89090E-14 4.87240E-14 &
     4.85390E-14 4.83550E-14 4.81710E-14 4.79870E-14 4.78040E-14 4.76210E-14 4.74380E-14 4.72550E-14 &
     4.70730E-14 4.68920E-14 4.67100E-14 4.65300E-14 4.63490E-14 4.61690E-14 4.59900E-14 4.58110E-14 &
     4.56320E-14 4.54540E-14 4.52760E-14 4.50990E-14 4.49220E-14 4.47460E-14 4.45710E-14 4.43960E-14 &
     4.42210E-14 4.40470E-14 4.38740E-14 4.37020E-14 4.35300E-14 4.33580E-14 4.31880E-14 4.30180E-14 &
     4.28490E-14 4.26800E-14 4.25130E-14 4.23450E-14 4.21790E-14 4.20140E-14 4.18490E-14 4.16850E-14 &
     4.15220E-14 4.13590E-14 4.11970E-14 4.10370E-14 4.08770E-14 4.07180E-14 4.05600E-14 4.04020E-14 &
     4.02460E-14 4.00900E-14 3.99360E-14 3.97820E-14 3.96290E-14 3.94770E-14 3.93260E-14 3.91760E-14 &
     3.90270E-14 3.88790E-14 3.87310E-14 3.85850E-14 3.84400E-14 3.82950E-14 3.81520E-14 3.80100E-14 &
     3.78680E-14 3.77280E-14 3.75880E-14 3.74500E-14 3.73120E-14 3.71750E-14 3.70400E-14 3.69050E-14 &
     3.67720E-14 3.66400E-14 3.65080E-14 3.63780E-14 3.62490E-14 3.61210E-14 3.59940E-14 3.58670E-14 &
     3.57420E-14 3.56180E-14 3.54950E-14 3.53740E-14 3.52530E-14 3.51330E-14 3.50140E-14 3.48970E-14 &
     3.47800E-14 3.46650E-14 3.45500E-14 3.44370E-14 3.43250E-14 3.42130E-14 3.41030E-14 3.39940E-14 &
     3.38860E-14 3.37790E-14 3.36730E-14 3.35680E-14 3.34640E-14 3.33610E-14 3.32590E-14 3.31570E-14 &
     3.30570E-14 3.29580E-14 3.28600E-14 3.27630E-14 3.26670E-14 3.25720E-14 3.24780E-14 3.23850E-14 &
     3.22930E-14 3.22020E-14 3.21120E-14 3.20220E-14 3.19340E-14 3.18470E-14 3.17600E-14 3.16750E-14 &
     3.15900E-14 3.15070E-14 3.14240E-14 3.13420E-14 3.12610E-14 3.11810E-14 3.11020E-14 3.10240E-14 &
     3.09460E-14 3.08700E-14 3.07940E-14 3.07190E-14 3.06450E-14 3.05720E-14 3.04990E-14 3.04280E-14 &
     3.03570E-14 3.02870E-14 3.02170E-14 3.01490E-14 3.00810E-14 3.00140E-14 2.99480E-14 2.98820E-14 &
     2.98170E-14 2.97530E-14 2.96900E-14 2.96270E-14 2.95650E-14 2.95030E-14 2.94420E-14 2.93820E-14 &
     2.93230E-14 2.92640E-14 2.92060E-14 2.91480E-14 2.90910E-14 2.90340E-14 2.89780E-14 2.89230E-14 &
     2.88680E-14 2.88140E-14 2.87600E-14 2.87060E-14 2.86540E-14 2.86010E-14 2.85490E-14 2.84980E-14 &
     2.84470E-14 2.83970E-14 2.83470E-14 2.82970E-14 2.82480E-14 2.81990E-14 2.81510E-14 2.81030E-14 &
     2.80550E-14 2.80080E-14 2.79610E-14 2.79150E-14 2.78690E-14 2.78230E-14 2.77770E-14 2.77320E-14 &
     2.76870E-14 2.76420E-14 2.75980E-14 2.75540E-14 2.75100E-14 2.74660E-14 2.74230E-14 2.73800E-14 &
     2.73370E-14 2.72940E-14 2.72520E-14 2.72090E-14 2.71670E-14 2.71250E-14 2.70840E-14 2.70420E-14 &
     2.70010E-14 2.69590E-14 2.69180E-14 2.68770E-14 2.68360E-14 2.67950E-14 2.67550E-14 2.67140E-14 &
     2.66730E-14 2.66330E-14 2.65920E-14 2.65520E-14 2.65120E-14 2.64710E-14 2.64310E-14 2.63910E-14 &
     2.63510E-14 2.63110E-14 2.62700E-14 2.62300E-14 2.61900E-14 2.61500E-14 2.61100E-14 2.60690E-14 &
     2.60290E-14 2.59890E-14 2.59480E-14 2.59080E-14 2.58680E-14 2.58270E-14 2.57870E-14 2.57460E-14 &
     2.57050E-14 2.56650E-14 2.56240E-14 2.55830E-14 2.55420E-14 2.55010E-14 2.54600E-14 2.54180E-14 &
     2.53770E-14 2.53360E-14 2.52940E-14 2.52520E-14 2.52110E-14 2.51690E-14 2.51270E-14 2.50850E-14 &
     2.50430E-14 2.50000E-14 2.49580E-14 2.49150E-14 2.48720E-14 2.48290E-14 2.47860E-14 2.47430E-14 &
     2.47000E-14 2.46570E-14 2.46130E-14 2.45690E-14 2.45260E-14 2.44820E-14 2.44370E-14 2.43930E-14 &
     2.43490E-14 2.43040E-14 2.42590E-14 2.42150E-14 2.41700E-14 2.41240E-14 2.40790E-14 2.40340E-14 &
     2.39880E-14 2.39420E-14 2.38960E-14 2.38500E-14 2.38040E-14 2.37580E-14 2.37110E-14 2.36640E-14 &
     2.36180E-14 2.35710E-14 2.35240E-14 2.34760E-14 2.34290E-14 2.33820E-14 2.33340E-14 2.32860E-14 &
     2.32390E-14 2.31910E-14 2.31430E-14 2.30950E-14 2.30460E-14 2.29980E-14 2.29500E-14 2.29010E-14 &
     2.28520E-14 2.28040E-14 2.27550E-14 2.27060E-14 2.26570E-14 2.26080E-14 2.25580E-14 2.25090E-14 &
     2.24600E-14 2.24100E-14 2.23610E-14 2.23110E-14 2.22610E-14 2.22110E-14 2.21610E-14 2.21110E-14 &
     2.20610E-14 2.20110E-14 2.19610E-14 2.19110E-14 2.18610E-14 2.18100E-14 2.17600E-14 2.17100E-14 &
     2.16590E-14 2.16090E-14 2.15590E-14 2.15080E-14 2.14580E-14 2.14070E-14 2.13570E-14 2.13060E-14 &
     2.12560E-14 2.12060E-14 2.11550E-14 2.11050E-14 2.10540E-14 2.10040E-14 2.09540E-14 2.09030E-14 &
     2.08530E-14 2.08030E-14 2.07530E-14 2.07030E-14 2.06530E-14 2.06030E-14 2.05530E-14 2.05030E-14 &
     2.04530E-14 2.04040E-14 2.03540E-14 2.03050E-14 2.02550E-14 2.02060E-14 2.01570E-14 2.01070E-14 &
     2.00580E-14 2.00100E-14 1.99610E-14 1.99120E-14 1.98640E-14 1.98150E-14 1.97670E-14 1.97190E-14 &
     1.96710E-14 1.96230E-14 1.95750E-14 1.95270E-14 1.94800E-14 1.94320E-14 1.93850E-14 1.93380E-14 &
     1.92910E-14 1.92440E-14 1.91970E-14 1.91500E-14 1.91040E-14 1.90570E-14 1.90110E-14 1.89650E-14 &
     1.89190E-14 1.88730E-14 1.88270E-14 1.87820E-14 1.87370E-14 1.86910E-14 1.86460E-14 1.86010E-14 &
     1.85560E-14 1.85120E-14 1.84670E-14 1.84230E-14 1.83790E-14 1.83350E-14 1.82910E-14 1.82470E-14 &
     1.82030E-14 1.81600E-14 1.81170E-14 1.80730E-14 1.80300E-14 1.79880E-14 1.79450E-14 1.79030E-14 &
     1.78600E-14 1.78180E-14 1.77760E-14 1.77340E-14 1.76930E-14 1.76510E-14 1.76100E-14 1.75690E-14 &
     1.75280E-14 1.74870E-14 1.74460E-14 1.74060E-14 1.73650E-14 1.73250E-14 1.72850E-14 1.72460E-14 &
     1.72060E-14 1.71670E-14 1.71280E-14 1.70890E-14 1.70500E-14 1.70110E-14 1.69730E-14 1.69340E-14 &
     1.68960E-14 1.68580E-14 1.68210E-14 1.67830E-14 1.67460E-14 1.67090E-14 1.66720E-14 1.66350E-14 &
     1.65980E-14 1.65620E-14 1.65260E-14 1.64900E-14 1.64540E-14 1.64180E-14 1.63830E-14 1.63470E-14 &
     1.63120E-14 1.62770E-14 1.62420E-14 1.62080E-14 1.61730E-14 1.61390E-14 1.61040E-14 1.60700E-14 &
     1.60360E-14 1.60030E-14 1.59690E-14 1.59350E-14 1.59020E-14 1.58690E-14 1.58360E-14 1.58030E-14 &
     1.57700E-14 1.57370E-14 1.57040E-14 1.56720E-14 1.56400E-14 1.56070E-14 1.55750E-14 1.55430E-14 &
     1.55110E-14 1.54790E-14 1.54470E-14 1.54160E-14 1.53840E-14 1.53530E-14 1.53210E-14 1.52900E-14 &
     1.52590E-14 1.52280E-14 1.51970E-14 1.51660E-14 1.51350E-14 1.51040E-14 1.50730E-14 1.50420E-14 &
     1.50120E-14 1.49810E-14 1.49510E-14 1.49200E-14 1.48900E-14 1.48600E-14 1.48290E-14 1.47990E-14 &
     1.47690E-14 1.47390E-14 1.47080E-14 1.46780E-14 1.46480E-14 1.46180E-14 1.45880E-14 1.45580E-14 &
     1.45280E-14 1.44980E-14 1.44690E-14 1.44390E-14 1.44090E-14 1.43790E-14 1.43490E-14 1.43190E-14 &
     1.42890E-14 1.42600E-14 1.42300E-14 1.42000E-14 1.41700E-14 1.41400E-14 1.41100E-14 1.40800E-14 &
     1.40500E-14 1.40200E-14 1.39910E-14 1.39610E-14 1.39310E-14 1.39010E-14 1.38700E-14 1.38400E-14 &
     1.38100E-14 1.37800E-14 1.37500E-14 1.37200E-14 1.36890E-14 1.36590E-14 1.36290E-14 1.35980E-14 &
     1.35680E-14 1.35370E-14 1.35060E-14 1.34760E-14 1.34450E-14 1.34140E-14 1.33830E-14 1.33520E-14 &
     1.33210E-14 1.32900E-14 1.32590E-14 1.32270E-14 1.31960E-14 1.31640E-14 1.31330E-14 1.31010E-14 &
     1.30690E-14 1.30380E-14 1.30060E-14 1.29740E-14 1.29420E-14 1.29090E-14 1.28770E-14 1.28450E-14 &
     1.28120E-14 1.27800E-14 1.27470E-14 1.27150E-14 1.26820E-14 1.26490E-14 1.26160E-14 1.25830E-14 &
     1.25500E-14 1.25170E-14 1.24840E-14 1.24510E-14 1.24180E-14 1.23840E-14 1.23510E-14 1.23180E-14 &
     1.22840E-14 1.22510E-14 1.22170E-14 1.21840E-14 1.21500E-14 1.21160E-14 1.20830E-14 1.20490E-14 &
     1.20150E-14 1.19810E-14 1.19470E-14 1.19130E-14 1.18790E-14 1.18450E-14 1.18110E-14 1.17770E-14 &
     1.17430E-14 1.17090E-14 1.16750E-14 1.16410E-14 1.16070E-14 1.15730E-14 1.15380E-14 1.15040E-14 &
     1.14700E-14 1.14360E-14 1.14010E-14 1.13670E-14 1.13330E-14 1.12990E-14 1.12640E-14 1.12300E-14 &
     1.11960E-14 1.11610E-14 1.11270E-14 1.10930E-14 1.10590E-14 1.10240E-14 1.09900E-14 1.09560E-14 &
     1.09210E-14 1.08870E-14 1.08530E-14 1.08190E-14 1.07840E-14 1.07500E-14 1.07160E-14 1.06820E-14 &
     1.06480E-14 1.06140E-14 1.05800E-14 1.05450E-14 1.05110E-14 1.04770E-14 1.04430E-14 1.04100E-14 &
     1.03760E-14 1.03420E-14 1.03080E-14 1.02740E-14 1.02400E-14 1.02070E-14 1.01730E-14 1.01390E-14 &
     1.01060E-14 1.00720E-14 1.00390E-14 1.00050E-14 9.97200E-15 9.93860E-15 9.90540E-15 9.87220E-15 &
     9.83900E-15 9.80590E-15 9.77290E-15 9.73990E-15 9.70700E-15 9.67420E-15 9.64140E-15 9.60880E-15 &
     9.57610E-15 9.54360E-15 9.51110E-15 9.47870E-15 9.44640E-15 9.41410E-15 9.38200E-15 9.34990E-15 &
     9.31790E-15 9.28600E-15 9.25420E-15 9.22240E-15 9.19080E-15 9.15920E-15 9.12780E-15 9.09640E-15 &
     9.06510E-15 9.03400E-15 9.00290E-15 8.97190E-15 8.94100E-15 8.91030E-15 8.87960E-15 8.84910E-15 &
     8.81860E-15 8.78830E-15 8.75810E-15 8.72790E-15 8.69790E-15 8.66810E-15 8.63830E-15 8.60870E-15 &
     8.57910E-15 8.54970E-15 8.52050E-15 8.49130E-15 8.46230E-15 8.43340E-15 8.40460E-15 8.37600E-15 &
     8.34750E-15 8.31910E-15 8.29090E-15 8.26280E-15 8.23480E-15 8.20700E-15 8.17930E-15 8.15180E-15 &
     8.12440E-15 8.09720E-15 8.07010E-15 8.04320E-15 8.01640E-15 7.98980E-15 7.96330E-15 7.93700E-15 &
     7.91080E-15 7.88480E-15 7.85900E-15 7.83330E-15 7.80780E-15 7.78240E-15 7.75720E-15 7.73220E-15 &
     7.70740E-15 7.68270E-15 7.65820E-15 7.63390E-15 7.60970E-15 7.58570E-15 7.56190E-15 7.53820E-15 &
     7.51470E-15 7.49130E-15 7.46810E-15 7.44510E-15 7.42230E-15 7.39960E-15 7.37700E-15 7.35470E-15 &
     7.33240E-15 7.31040E-15 7.28850E-15 7.26670E-15 7.24510E-15 7.22370E-15 7.20240E-15 7.18130E-15 &
     7.16030E-15 7.13950E-15 7.11880E-15 7.09830E-15 7.07790E-15 7.05770E-15 7.03760E-15 7.01770E-15 &
     6.99790E-15 6.97830E-15 6.95880E-15 6.93940E-15 6.92020E-15 6.90120E-15 6.88220E-15 6.86350E-15 &
     6.84480E-15 6.82630E-15 6.80800E-15 6.78980E-15 6.77170E-15 6.75380E-15 6.73600E-15 6.71830E-15 &
     6.70080E-15 6.68340E-15 6.66610E-15 6.64900E-15 6.63200E-15 6.61510E-15 6.59840E-15 6.58180E-15 &
     6.56530E-15 6.54890E-15 6.53270E-15 6.51660E-15 6.50070E-15 6.48480E-15 6.46910E-15 6.45350E-15 &
     6.43800E-15 6.42270E-15 6.40750E-15 6.39240E-15 6.37740E-15 6.36250E-15 6.34780E-15 6.33320E-15 &
     6.31870E-15 6.30430E-15 6.29000E-15 6.27580E-15 6.26180E-15 6.24790E-15 6.23400E-15 6.22030E-15 &
     6.20670E-15 6.19330E-15 6.17990E-15 6.16660E-15 6.15350E-15 6.14040E-15 6.12750E-15 6.11470E-15 &
     6.10200E-15 6.08930E-15 6.07680E-15 6.06440E-15 6.05210E-15 6.03990E-15 6.02780E-15 6.01580E-15 &
     6.00390E-15 5.99210E-15 5.98040E-15 5.96880E-15 5.95730E-15 5.94590E-15 5.93450E-15 5.92330E-15 &
     5.91220E-15 5.90120E-15 5.89020E-15 5.87940E-15 5.86860E-15 5.85790E-15 5.84740E-15 5.83690E-15 &
     5.82650E-15 5.81620E-15 5.80590E-15 5.79580E-15 5.78570E-15 5.77580E-15 5.76590E-15 5.75610E-15 &
     5.74640E-15 5.73670E-15 5.72720E-15 5.71770E-15 5.70830E-15 5.69900E-15 5.68970E-15 5.68060E-15 &
     5.67150E-15 5.66250E-15 5.65360E-15 5.64470E-15 5.63590E-15 5.62720E-15 5.61860E-15 5.61000E-15 &
     5.60150E-15 5.59310E-15 5.58470E-15 5.57640E-15 5.56820E-15 5.56010E-15 5.55200E-15 5.54390E-15 &
     5.53600E-15 5.52810E-15 5.52030E-15 5.51250E-15 5.50480E-15 5.49720E-15 5.48960E-15 5.48210E-15 &
     5.47460E-15 5.46720E-15 5.45980E-15 5.45260E-15 5.44530E-15 5.43820E-15 5.43100E-15 5.42400E-15 &
     5.41690E-15 5.41000E-15 5.40310E-15 5.39620E-15 5.38940E-15 5.38260E-15 5.37590E-15 5.36930E-15 &
     5.36270E-15 5.35610E-15 5.34960E-15 5.34310E-15 5.33670E-15 5.33030E-15 5.32390E-15 5.31760E-15 &
     5.31140E-15 5.30520E-15 5.29900E-15 5.29280E-15 5.28670E-15 5.28070E-15 5.27470E-15 5.26870E-15 &
     5.26270E-15 5.25680E-15 5.25090E-15 5.24510E-15 5.23920E-15 5.23350E-15 5.22770E-15 5.22200E-15 &
     5.21630E-15 5.21060E-15 5.20500E-15 5.19940E-15 5.19380E-15 5.18830E-15 5.18270E-15 5.17720E-15 &
     5.17180E-15 5.16630E-15 5.16090E-15 5.15550E-15 5.15010E-15 5.14470E-15 5.13930E-15 5.13400E-15 &
     5.12870E-15 5.12340E-15 5.11810E-15 5.11290E-15 5.10760E-15 5.10240E-15 5.09720E-15 5.09200E-15 &
     5.08680E-15 5.08160E-15 5.07640E-15 5.07130E-15 5.06610E-15 5.06100E-15 5.05590E-15 5.05070E-15 &
     5.04560E-15 5.04050E-15 5.03540E-15 5.03030E-15 5.02520E-15 5.02010E-15 5.01500E-15 5.00990E-15 &
     5.00480E-15 4.99970E-15 4.99460E-15 4.98950E-15 4.98450E-15 4.97940E-15 4.97430E-15 4.96910E-15 &
     4.96400E-15 4.95890E-15 4.95380E-15 4.94870E-15 4.94350E-15 4.93840E-15 4.93320E-15 4.92810E-15 &
     4.92290E-15 4.91770E-15 4.91250E-15 4.90730E-15 4.90200E-15 4.89680E-15 4.89150E-15 4.88630E-15 &
     4.88100E-15 4.87570E-15 4.87030E-15 4.86500E-15 4.85970E-15 4.85430E-15 4.84890E-15 4.84350E-15 &
     4.83810E-15 4.83260E-15 4.82720E-15 4.82170E-15 4.81630E-15 4.81080E-15 4.80530E-15 4.79970E-15 &
     4.79420E-15 4.78860E-15 4.78310E-15 4.77750E-15 4.77190E-15 4.76630E-15 4.76060E-15 4.75500E-15 &
     4.74930E-15 4.74370E-15 4.73800E-15 4.73230E-15 4.72650E-15 4.72080E-15 4.71510E-15 4.70930E-15 &
     4.70350E-15 4.69770E-15 4.69190E-15 4.68610E-15 4.68030E-15 4.67440E-15 4.66860E-15 4.66270E-15 &
     4.65680E-15 4.65090E-15 4.64500E-15 4.63900E-15 4.63310E-15 4.62710E-15 4.62120E-15 4.61520E-15 &
     4.60920E-15 4.60320E-15 4.59710E-15 4.59110E-15 4.58500E-15 4.57900E-15 4.57290E-15 4.56680E-15 &
     4.56070E-15 4.55450E-15 4.54840E-15 4.54230E-15 4.53610E-15 4.52990E-15 4.52370E-15 4.51750E-15 &
     4.51130E-15 4.50510E-15 4.49880E-15 4.49260E-15 4.48630E-15 4.48000E-15 4.47370E-15 4.46740E-15 &
     4.46110E-15 4.45480E-15 4.44850E-15 4.44210E-15 4.43570E-15 4.42940E-15 4.42300E-15 4.41660E-15 &
     4.41010E-15 4.40370E-15 4.39730E-15 4.39080E-15 4.38440E-15 4.37790E-15 4.37140E-15 4.36490E-15 &
     4.35840E-15 4.35180E-15 4.34530E-15 4.33880E-15 4.33220E-15 4.32560E-15 4.31910E-15 4.31250E-15 &
     4.30590E-15 4.29920E-15 4.29260E-15 4.28600E-15 4.27930E-15 4.27270E-15 4.26600E-15 4.25930E-15 &
     4.25260E-15 4.24590E-15 4.23920E-15 4.23240E-15 4.22570E-15 4.21900E-15 4.21220E-15 4.20540E-15 &
     4.19860E-15 4.19180E-15 4.18500E-15 4.17820E-15 4.17140E-15 4.16460E-15 4.15770E-15 4.15090E-15 &
     4.14400E-15 4.13710E-15 4.13020E-15 4.12330E-15 4.11640E-15 4.10950E-15 4.10260E-15 4.09560E-15 &
     4.08870E-15 4.08170E-15 4.07480E-15 4.06780E-15 4.06080E-15 4.05380E-15 4.04680E-15 4.03980E-15 &
     4.03270E-15 4.02570E-15 4.01870E-15 4.01160E-15 4.00450E-15 3.99750E-15 3.99040E-15 3.98330E-15 &
     3.97620E-15 3.96910E-15 3.96190E-15 3.95480E-15 3.94770E-15 3.94050E-15 3.93340E-15 3.92620E-15 &
     3.91900E-15 3.91190E-15 3.90470E-15 3.89750E-15 3.89020E-15 3.88300E-15 3.87580E-15 3.86860E-15 &
     3.86130E-15 3.85410E-15 3.84680E-15 3.83950E-15 3.83230E-15 3.82500E-15 3.81770E-15 3.81040E-15 &
     3.80310E-15 3.79580E-15 3.78840E-15 3.78110E-15 3.77380E-15 3.76640E-15 3.75900E-15 3.75170E-15 &
     3.74430E-15 3.73690E-15 3.72950E-15 3.72210E-15 3.71470E-15 3.70730E-15 3.69990E-15 3.69250E-15 &
     3.68500E-15 3.67760E-15 3.67020E-15 3.66270E-15 3.65520E-15 3.64780E-15 3.64030E-15 3.63280E-15 &
     3.62530E-15 3.61780E-15 3.61030E-15 3.60280E-15 3.59530E-15 3.58780E-15 3.58020E-15 3.57270E-15 &
     3.56510E-15 3.55760E-15 3.55000E-15 3.54250E-15 3.53490E-15 3.52730E-15 3.51970E-15 3.51210E-15 &
     3.50450E-15 3.49690E-15 3.48930E-15 3.48170E-15 3.47410E-15 3.46650E-15 3.45880E-15 3.45120E-15 &
     3.44360E-15 3.43590E-15 3.42820E-15 3.42060E-15 3.41290E-15 3.40520E-15 3.39760E-15 3.38990E-15 &
     3.38220E-15 3.37450E-15 3.36680E-15 3.35910E-15 3.35140E-15 3.34370E-15 3.33590E-15 3.32820E-15 &
     3.32050E-15 3.31270E-15 3.30500E-15 3.29720E-15 3.28950E-15 3.28170E-15 3.27400E-15 3.26620E-15 &
     3.25840E-15 3.25070E-15 3.24290E-15 3.23510E-15 3.22730E-15 3.21950E-15 3.21170E-15 3.20390E-15 &
     3.19610E-15 3.18830E-15 3.18050E-15 3.17260E-15 3.16480E-15 3.15700E-15 3.14910E-15 3.14130E-15 &
     3.13350E-15 3.12560E-15 3.11780E-15 3.10990E-15 3.10210E-15 3.09420E-15 3.08630E-15 3.07850E-15 &
     3.07060E-15 3.06270E-15 3.05480E-15 3.04690E-15 3.03910E-15 3.03120E-15 3.02330E-15 3.01540E-15 &
     3.00750E-15 2.99960E-15 2.99170E-15 2.98380E-15 2.97580E-15 2.96790E-15 2.96000E-15 2.95210E-15 &
     2.94420E-15 2.93620E-15 2.92830E-15 2.92040E-15 2.91240E-15 2.90450E-15 2.89660E-15 2.88860E-15 &
     2.88070E-15 2.87270E-15 2.86480E-15 2.85680E-15 2.84890E-15 2.84090E-15 2.83290E-15 2.82500E-15 &
     2.81700E-15 2.80910E-15 2.80110E-15 2.79310E-15 2.78510E-15 2.77720E-15 2.76920E-15 2.76120E-15 &
     2.75320E-15 2.74530E-15 2.73730E-15 2.72930E-15 2.72130E-15 2.71330E-15 2.70530E-15 2.69730E-15 &
     2.68930E-15 2.68130E-15 2.67330E-15 2.66540E-15 2.65740E-15 2.64940E-15 2.64140E-15 2.63340E-15 &
     2.62540E-15 2.61740E-15 2.60930E-15 2.60130E-15 2.59330E-15 2.58530E-15 2.57730E-15 2.56930E-15 &
     2.56130E-15 2.55330E-15 2.54530E-15 2.53730E-15 2.52930E-15 2.52130E-15 2.51320E-15 2.50520E-15 &
     2.49720E-15 2.48920E-15 2.48120E-15 2.47320E-15 2.46520E-15 2.45720E-15 7.56760E-15 1.21470E-14 &
     1.33940E-14 1.16410E-14 8.42030E-15 5.25410E-15 2.99730E-15 1.77920E-15 1.29520E-15 1.14810E-15 &
     1.06650E-15 9.52020E-16 8.22860E-16 7.29170E-16 6.95920E-16 7.13630E-16 7.57410E-16 8.06810E-16 &
     8.57520E-16 9.20910E-16 1.01150E-15 1.13540E-15 1.29430E-15 1.48720E-15 1.71070E-15 1.96000E-15 &
     2.23120E-15 2.52730E-15 2.85030E-15 3.19620E-15 3.56520E-15 3.95390E-15 4.36990E-15 4.81610E-15 &
     5.29460E-15 5.80570E-15 6.34640E-15 6.91940E-15 7.52910E-15 8.17850E-15 8.86660E-15 9.59060E-15 &
     1.03480E-14 1.11360E-14 1.19560E-14 1.28090E-14 1.36990E-14 1.46290E-14 1.56020E-14 1.66150E-14 &
     1.76670E-14 1.87580E-14 1.98860E-14 2.10490E-14 2.22490E-14 2.34850E-14 2.47580E-14 2.60700E-14 &
     2.74190E-14 2.88060E-14 3.02310E-14 3.16910E-14 3.31860E-14 3.47160E-14 3.62780E-14 3.78710E-14 &
     3.94950E-14 4.11470E-14 4.28290E-14 4.45380E-14 4.62770E-14 4.80440E-14 4.98380E-14 5.16570E-14 &
     5.34980E-14 5.53600E-14 5.72420E-14 5.91420E-14 6.10590E-14 6.29900E-14 6.49360E-14 6.68940E-14 &
     6.88630E-14 7.08400E-14 7.28260E-14 7.48170E-14 7.68130E-14 7.88110E-14 8.08100E-14 8.28080E-14 &
     8.48030E-14 8.67940E-14 8.87780E-14 9.07540E-14 9.27190E-14 9.46730E-14 9.66120E-14 9.85350E-14 &
     1.00440E-13 1.02330E-13 1.04190E-13 1.06030E-13 1.07840E-13 1.09630E-13 1.11390E-13 1.13110E-13 &
     1.14810E-13 1.16470E-13 1.18090E-13 1.19680E-13 1.21220E-13 1.22730E-13 1.24190E-13 1.25610E-13 &
     1.26980E-13 1.28300E-13 1.29580E-13 1.30810E-13 1.31980E-13 1.33110E-13 1.34180E-13 1.35190E-13 &
     1.36150E-13 1.37050E-13 1.37900E-13 1.38680E-13 1.39410E-13 1.40070E-13 1.40670E-13 1.41210E-13 &
     1.41690E-13 1.42110E-13 1.42460E-13 1.42750E-13 1.42970E-13 1.43130E-13 1.43230E-13 1.43270E-13 &
     1.43240E-13 1.43140E-13 1.42990E-13 1.42770E-13 1.42490E-13 1.42150E-13 1.41750E-13 1.41290E-13 &
     1.40770E-13 1.40190E-13 1.39560E-13 1.38870E-13 1.38120E-13 1.37330E-13 1.36480E-13 1.35580E-13 &
     1.34630E-13 1.33640E-13 1.32600E-13 1.31520E-13 1.30390E-13 1.29230E-13 1.28030E-13 1.26800E-13 &
     1.25530E-13 1.24240E-13 1.22920E-13 1.21570E-13 1.20210E-13 1.18820E-13 1.17430E-13 1.16020E-13 &
     1.14600E-13 1.13170E-13 1.11750E-13 1.10330E-13 1.08910E-13 1.07500E-13 1.06110E-13 1.04730E-13 &
     1.03370E-13 1.02040E-13 1.00740E-13 9.94640E-14 9.82310E-14 9.70390E-14 9.58920E-14 9.47940E-14 &
     9.37510E-14 9.27670E-14 9.18460E-14 9.09920E-14 9.02100E-14 8.95050E-14 8.88820E-14 8.83440E-14 &
     8.78980E-14 8.75470E-14 8.72970E-14 8.71520E-14 8.71170E-14 8.71970E-14 8.73970E-14 8.77220E-14 &
     8.81770E-14 8.87660E-14 8.94950E-14 9.03680E-14 9.13910E-14 9.25670E-14 9.39030E-14 9.54020E-14 &
     9.70710E-14 9.89130E-14 1.00930E-13 1.03140E-13 1.05530E-13 1.08110E-13 1.10890E-13 1.13870E-13 &
     1.17060E-13 1.20460E-13 1.24070E-13 1.27910E-13 1.31960E-13 1.36250E-13 1.40760E-13 1.45520E-13 &
     1.50510E-13 1.55750E-13 1.61230E-13 1.66960E-13 1.72950E-13 1.79200E-13 1.85710E-13 1.92490E-13 &
     1.99530E-13 2.06850E-13 2.14440E-13 2.22310E-13 2.30450E-13 2.38880E-13 2.47600E-13 2.56600E-13 &
     2.65890E-13 2.75480E-13 2.85360E-13 2.95530E-13 3.06000E-13 3.16770E-13 3.27850E-13 3.39220E-13 &
     3.50900E-13 3.62870E-13 3.75160E-13 3.87740E-13 4.00640E-13 4.13830E-13 4.27340E-13 4.41140E-13 &
     4.55260E-13 4.69670E-13 4.84400E-13 4.99430E-13 5.14760E-13 5.30390E-13 5.46330E-13 5.62560E-13 &
     5.79090E-13 5.95920E-13 6.13040E-13 6.30450E-13 6.48160E-13 6.66140E-13 6.84420E-13 7.02970E-13 &
     7.21800E-13 7.40910E-13 7.60300E-13 7.79950E-13 7.99860E-13 8.20040E-13 8.40470E-13 8.61150E-13 &
     8.82080E-13 9.03260E-13 9.24670E-13 9.46320E-13 9.68200E-13 9.90300E-13 1.01260E-12 1.03520E-12 &
     1.05790E-12 1.08080E-12 1.10400E-12 1.12730E-12 1.15080E-12 1.17440E-12 1.19830E-12 1.22230E-12 &
     1.24640E-12 1.27070E-12 1.29510E-12 1.31970E-12 1.34440E-12 1.36920E-12 1.39410E-12 1.41910E-12 &
     1.44430E-12 1.46950E-12 1.49470E-12 1.52010E-12 1.54550E-12 1.57100E-12 1.59650E-12 1.62200E-12 &
     1.64760E-12 1.67310E-12 1.69870E-12 1.72430E-12 1.74990E-12 1.77540E-12 1.80090E-12 1.82640E-12 &
     1.85190E-12 1.87720E-12 1.90250E-12 1.92780E-12 1.95290E-12 1.97800E-12 2.00290E-12 2.02770E-12 &
     2.05240E-12 2.07700E-12 2.10140E-12 2.12570E-12 2.14970E-12 2.17370E-12 2.19740E-12 2.22100E-12 &
     2.24430E-12 2.26740E-12 2.29040E-12 2.31300E-12 2.33550E-12 2.35760E-12 2.37960E-12 2.40120E-12 &
     2.42250E-12 2.44360E-12 2.46430E-12 2.48480E-12 2.50490E-12 2.52470E-12 2.54420E-12 2.56340E-12 &
     2.58220E-12 2.60070E-12 2.61880E-12 2.63650E-12 2.65380E-12 2.67080E-12 2.68730E-12 2.70350E-12 &
     2.71920E-12 2.73450E-12 2.74940E-12 2.76380E-12 2.77780E-12 2.79140E-12 2.80450E-12 2.81710E-12 &
     2.82930E-12 2.84100E-12 2.85230E-12 2.86300E-12 2.87330E-12 2.88310E-12 2.89230E-12 2.90110E-12 &
     2.90930E-12 2.91710E-12 2.92430E-12 2.93100E-12 2.93710E-12 2.94280E-12 2.94790E-12 2.95250E-12 &
     2.95650E-12 2.96000E-12 2.96300E-12 2.96550E-12 2.96740E-12 2.96880E-12 2.96960E-12 2.96990E-12 &
     2.96970E-12 2.96890E-12 2.96760E-12 2.96570E-12 2.96330E-12 2.96030E-12 2.95680E-12 2.95280E-12 &
     2.94820E-12 2.94300E-12 2.93740E-12 2.93120E-12 2.92450E-12 2.91730E-12 2.90950E-12 2.90120E-12 &
     2.89240E-12 2.88310E-12 2.87330E-12 2.86300E-12 2.85220E-12 2.84090E-12 2.82910E-12 2.81680E-12 &
     2.80410E-12 2.79090E-12 2.77720E-12 2.76300E-12 2.74840E-12 2.73330E-12 2.71780E-12 2.70190E-12 &
     2.68550E-12 2.66870E-12 2.65160E-12 2.63400E-12 2.61600E-12 2.59770E-12 2.57900E-12 2.56000E-12 &
     2.54060E-12 2.52090E-12 2.50080E-12 2.48040E-12 2.45970E-12 2.43870E-12 2.41740E-12 2.39580E-12 &
     2.37400E-12 2.35180E-12 2.32950E-12 2.30690E-12 2.28410E-12 2.26100E-12 2.23780E-12 2.21440E-12 &
     2.19080E-12 2.16700E-12 2.14310E-12 2.11900E-12 2.09480E-12 2.07050E-12 2.04610E-12 2.02170E-12 &
     1.99710E-12 1.97250E-12 1.94780E-12 1.92310E-12 1.89840E-12 1.87370E-12 1.84890E-12 1.82420E-12 &
     1.79950E-12 1.77480E-12 1.75020E-12 1.72570E-12 1.70130E-12 1.67700E-12 1.65280E-12 1.62870E-12 &
     1.60480E-12 1.58110E-12 1.55750E-12 1.53410E-12 1.51090E-12 1.48790E-12 1.46510E-12 1.44260E-12 &
     1.42030E-12 1.39830E-12 1.37660E-12 1.35520E-12 1.33400E-12 1.31320E-12 1.29270E-12 1.27260E-12 &
     1.25280E-12 1.23330E-12 1.21430E-12 1.19560E-12 1.17740E-12 1.15960E-12 1.14210E-12 1.12520E-12 &
     1.10860E-12 1.09260E-12 1.07700E-12 1.06190E-12 1.04720E-12 1.03310E-12 1.01950E-12 1.00640E-12 &
     9.93910E-13 9.81920E-13 9.70490E-13 9.59610E-13 9.49310E-13 9.39590E-13 9.30460E-13 9.21930E-13 &
     9.14000E-13 9.06690E-13 8.99990E-13 8.93930E-13 8.88490E-13 8.83700E-13 8.79560E-13 8.76080E-13 &
     8.73250E-13 8.71090E-13 8.69600E-13 8.68780E-13 8.68630E-13 8.69180E-13 8.70410E-13 8.72320E-13 &
     8.74940E-13 8.78240E-13 8.82250E-13 8.86960E-13 8.92370E-13 8.98490E-13 9.05310E-13 9.12840E-13 &
     9.21090E-13 9.30040E-13 9.39700E-13 9.50070E-13 9.61160E-13 9.72950E-13 9.85440E-13 9.98650E-13 &
     1.01260E-12 1.02720E-12 1.04250E-12 1.05850E-12 1.07520E-12 1.09250E-12 1.11060E-12 1.12940E-12 &
     1.14880E-12 1.16890E-12 1.18970E-12 1.21110E-12 1.23320E-12 1.25590E-12 1.27930E-12 1.30330E-12 &
     1.32790E-12 1.35320E-12 1.37910E-12 1.40550E-12 1.43260E-12 1.46030E-12 1.48850E-12 1.51740E-12 &
     1.54680E-12 1.57670E-12 1.60720E-12 1.63820E-12 1.66970E-12 1.70180E-12 1.73430E-12 1.76740E-12 &
     1.80090E-12 1.83490E-12 1.86930E-12 1.90420E-12 1.93950E-12 1.97520E-12 2.01130E-12 2.04780E-12 &
     2.08470E-12 2.12200E-12 2.15960E-12 2.19760E-12 2.23580E-12 2.27440E-12 2.31330E-12 2.35240E-12 &
     2.39190E-12 2.43150E-12 2.47140E-12 2.51150E-12 2.55190E-12 2.59240E-12 2.63300E-12 2.67390E-12 &
     2.71490E-12 2.75600E-12 2.79720E-12 2.83850E-12 2.87990E-12 2.92130E-12 2.96280E-12 3.00430E-12 &
     3.04580E-12 3.08730E-12 3.12880E-12 3.17020E-12 3.21160E-12 3.25290E-12 3.29420E-12 3.33530E-12 &
     3.37630E-12 3.41720E-12 3.45800E-12 3.49850E-12 3.53890E-12 3.57910E-12 3.61910E-12 3.65880E-12 &
     3.69830E-12 3.73760E-12 3.77650E-12 3.81520E-12 3.85350E-12 3.89160E-12 3.92920E-12 3.96650E-12 &
     4.00350E-12 4.04000E-12 4.07620E-12 4.11190E-12 4.14720E-12 4.18210E-12 4.21640E-12 4.25040E-12 &
     4.28380E-12 4.31670E-12 4.34910E-12 4.38090E-12 4.41220E-12 4.44290E-12 4.47300E-12 4.50260E-12 &
     4.53160E-12 4.55990E-12 4.58770E-12 4.61480E-12 4.64130E-12 4.66720E-12 4.69240E-12 4.71700E-12 &
     4.74080E-12 4.76400E-12 4.78650E-12 4.80830E-12 4.82930E-12 4.84970E-12 4.86930E-12 4.88810E-12 &
     4.90630E-12 4.92360E-12 4.94020E-12 4.95610E-12 4.97110E-12 4.98540E-12 4.99890E-12 5.01160E-12 &
     5.02350E-12 5.03460E-12 5.04490E-12 5.05440E-12 5.06310E-12 5.07090E-12 5.07790E-12 5.08410E-12 &
     5.08940E-12 5.09390E-12 5.09760E-12 5.10040E-12 5.10240E-12 5.10350E-12 5.10390E-12 5.10340E-12 &
     5.10210E-12 5.10000E-12 5.09700E-12 5.09330E-12 5.08870E-12 5.08330E-12 5.07700E-12 5.07000E-12 &
     5.06210E-12 5.05340E-12 5.04390E-12 5.03360E-12 5.02260E-12 5.01070E-12 4.99800E-12 4.98460E-12 &
     4.97040E-12 4.95550E-12 4.93980E-12 4.92330E-12 4.90610E-12 4.88820E-12 4.86960E-12 4.85020E-12 &
     4.83010E-12 4.80930E-12 4.78780E-12 4.76560E-12 4.74270E-12 4.71920E-12 4.69500E-12 4.67020E-12 &
     4.64480E-12 4.61870E-12 4.59210E-12 4.56490E-12 4.53710E-12 4.50870E-12 4.47980E-12 4.45030E-12 &
     4.42040E-12 4.38990E-12 4.35890E-12 4.32740E-12 4.29540E-12 4.26300E-12 4.23020E-12 4.19690E-12 &
     4.16320E-12 4.12910E-12 4.09470E-12 4.05990E-12 4.02470E-12 3.98920E-12 3.95340E-12 3.91730E-12 &
     3.88100E-12 3.84430E-12 3.80740E-12 3.77030E-12 3.73300E-12 3.69540E-12 3.65770E-12 3.61980E-12 &
     3.58180E-12 3.54370E-12 3.50540E-12 3.46700E-12 3.42860E-12 3.39010E-12 3.35160E-12 3.31310E-12 &
     3.27460E-12 3.23600E-12 3.19750E-12 3.15910E-12 3.12070E-12 3.08230E-12 3.04410E-12 3.00600E-12 &
     2.96800E-12 2.93020E-12 2.89260E-12 2.85510E-12 2.81790E-12 2.78090E-12 2.74410E-12 2.70760E-12 &
     2.67130E-12 2.63540E-12 2.59970E-12 2.56440E-12 2.52940E-12 2.49470E-12 2.46040E-12 2.42650E-12 &
     2.39290E-12 2.35980E-12 2.32710E-12 2.29480E-12 2.26300E-12 2.23170E-12 2.20090E-12 2.17060E-12 &
     2.14090E-12 2.11160E-12 2.08300E-12 2.05490E-12 2.02730E-12 2.00040E-12 1.97400E-12 1.94820E-12 &
     1.92300E-12 1.89850E-12 1.87450E-12 1.85130E-12 1.82870E-12 1.80690E-12 1.78570E-12 1.76520E-12 &
     1.74550E-12 1.72640E-12 1.70810E-12 1.69050E-12 1.67370E-12 1.65750E-12 1.64220E-12 1.62760E-12 &
     1.61380E-12 1.60070E-12 1.58850E-12 1.57700E-12 1.56640E-12 1.55650E-12 1.54750E-12 1.53930E-12 &
     1.53190E-12 1.52530E-12 1.51960E-12 1.51460E-12 1.51050E-12 1.50730E-12 1.50480E-12 1.50330E-12 &
     1.50250E-12 1.50270E-12 1.50370E-12 1.50550E-12 1.50820E-12 1.51170E-12 1.51600E-12 1.52120E-12 &
     1.52720E-12 1.53400E-12 1.54170E-12 1.55020E-12 1.55950E-12 1.56970E-12 1.58070E-12 1.59260E-12 &
     1.60530E-12 1.61880E-12 1.63310E-12 1.64820E-12 1.66410E-12 1.68080E-12 1.69830E-12 1.71640E-12 &
     1.73540E-12 1.75510E-12 1.77560E-12 1.79680E-12 1.81870E-12 1.84140E-12 1.86490E-12 1.88900E-12 &
     1.91390E-12 1.93950E-12 1.96580E-12 1.99270E-12 2.02030E-12 2.04860E-12 2.07750E-12 2.10700E-12 &
     2.13720E-12 2.16790E-12 2.19920E-12 2.23110E-12 2.26350E-12 2.29640E-12 2.32990E-12 2.36400E-12 &
     2.39850E-12 2.43350E-12 2.46910E-12 2.50510E-12 2.54160E-12 2.57850E-12 2.61590E-12 2.65370E-12 &
     2.69180E-12 2.73040E-12 2.76930E-12 2.80860E-12 2.84820E-12 2.88820E-12 2.92840E-12 2.96890E-12 &
     3.00970E-12 3.05070E-12 3.09190E-12 3.13330E-12 3.17490E-12 3.21670E-12 3.25870E-12 3.30080E-12 &
     3.34300E-12 3.38530E-12 3.42780E-12 3.47030E-12 3.51290E-12 3.55550E-12 3.59820E-12 3.64080E-12 &
     3.68350E-12 3.72620E-12 3.76880E-12 3.81140E-12 3.85390E-12 3.89630E-12 3.93860E-12 3.98070E-12 &
     4.02280E-12 4.06460E-12 4.10640E-12 4.14790E-12 4.18920E-12 4.23030E-12 4.27110E-12 4.31170E-12 &
     4.35200E-12 4.39200E-12 4.43180E-12 4.47110E-12 4.51020E-12 4.54890E-12 4.58720E-12 4.62510E-12 &
     4.66260E-12 4.69970E-12 4.73640E-12 4.77260E-12 4.80830E-12 4.84350E-12 4.87830E-12 4.91250E-12 &
     4.94630E-12 4.97940E-12 5.01210E-12 5.04420E-12 5.07570E-12 5.10660E-12 5.13690E-12 5.16660E-12 &
     5.19560E-12 5.22400E-12 5.25180E-12 5.27890E-12 5.30540E-12 5.33120E-12 5.35630E-12 5.38070E-12 &
     5.40440E-12 5.42740E-12 5.44970E-12 5.47130E-12 5.49210E-12 5.51220E-12 5.53150E-12 5.55010E-12 &
     5.56790E-12 5.58490E-12 5.60120E-12 5.61670E-12 5.63130E-12 5.64520E-12 5.65830E-12 5.67060E-12 &
     5.68210E-12 5.69270E-12 5.70260E-12 5.71160E-12 5.71990E-12 5.72730E-12 5.73390E-12 5.73960E-12 &
     5.74450E-12 5.74860E-12 5.75190E-12 5.75430E-12 5.75580E-12 5.75650E-12 5.75640E-12 5.75540E-12 &
     5.75350E-12 5.75090E-12 5.74730E-12 5.74300E-12 5.73780E-12 5.73170E-12 5.72480E-12 5.71710E-12 &
     5.70850E-12 5.69910E-12 5.68880E-12 5.67770E-12 5.66580E-12 5.65320E-12 5.63980E-12 5.62570E-12 &
     5.61080E-12 5.59520E-12 5.57890E-12 5.56180E-12 5.54400E-12 5.52550E-12 5.50610E-12 5.48610E-12 &
     5.46520E-12 5.44370E-12 5.42140E-12 5.39840E-12 5.37470E-12 5.35040E-12 5.32540E-12 5.29980E-12 &
     5.27350E-12 5.24660E-12 5.21910E-12 5.19100E-12 5.16230E-12 5.13300E-12 5.10320E-12 5.07280E-12 &
     5.04190E-12 5.01050E-12 4.97860E-12 4.94620E-12 4.91330E-12 4.88000E-12 4.84620E-12 4.81190E-12 &
     4.77710E-12 4.74200E-12 4.70630E-12 4.67020E-12 4.63370E-12 4.59680E-12 4.55950E-12 4.52190E-12 &
     4.48390E-12 4.44570E-12 4.40710E-12 4.36840E-12 4.32940E-12 4.29010E-12 4.25070E-12 4.21100E-12 &
     4.17110E-12 4.13110E-12 4.09090E-12 4.05050E-12 4.01000E-12 3.96930E-12 3.92860E-12 3.88780E-12 &
     3.84680E-12 3.80590E-12 3.76480E-12 3.72380E-12 3.68270E-12 3.64160E-12 3.60060E-12 3.55960E-12 &
     3.51860E-12 3.47770E-12 3.43690E-12 3.39610E-12 3.35550E-12 3.31510E-12 3.27470E-12 3.23450E-12 &
     3.19450E-12 3.15470E-12 3.11510E-12 3.07570E-12 3.03650E-12 2.99760E-12 2.95890E-12 2.92040E-12 &
     2.88230E-12 2.84450E-12 2.80690E-12 2.76980E-12 2.73290E-12 2.69640E-12 2.66030E-12 2.62460E-12 &
     2.58930E-12 2.55430E-12 2.51980E-12 2.48570E-12 2.45210E-12 2.41890E-12 2.38610E-12 2.35390E-12 &
     2.32210E-12 2.29080E-12 2.26000E-12 2.22980E-12 2.20010E-12 2.17090E-12 2.14220E-12 2.11410E-12 &
     2.08660E-12 2.05960E-12 2.03330E-12 2.00750E-12 1.98230E-12 1.95770E-12 1.93370E-12 1.91040E-12 &
     1.88770E-12 1.86560E-12 1.84410E-12 1.82330E-12 1.80320E-12 1.78370E-12 1.76490E-12 1.74670E-12 &
     1.72920E-12 1.71250E-12 1.69630E-12 1.68090E-12 1.66620E-12 1.65210E-12 1.63880E-12 1.62610E-12 &
     1.61420E-12 1.60290E-12 1.59240E-12 1.58260E-12 1.57340E-12 1.56500E-12 1.55730E-12 1.55030E-12 &
     1.54400E-12 1.53840E-12 1.53360E-12 1.52940E-12 1.52590E-12 1.52320E-12 1.52110E-12 1.51980E-12 &
     1.51910E-12 1.51910E-12 1.51990E-12 1.52130E-12 1.52340E-12 1.52620E-12 1.52960E-12 1.53380E-12 &
     1.53860E-12 1.54410E-12 1.55020E-12 1.55700E-12 1.56440E-12 1.57250E-12 1.58120E-12 1.59050E-12 &
     1.60050E-12 1.61110E-12 1.62220E-12 1.63400E-12 1.64640E-12 1.65940E-12 1.67290E-12 1.68710E-12 &
     1.70180E-12 1.71700E-12 1.73280E-12 1.74910E-12 1.76590E-12 1.78330E-12 1.80120E-12 1.81960E-12 &
     1.83840E-12 1.85770E-12 1.87750E-12 1.89780E-12 1.91850E-12 1.93960E-12 1.96120E-12 1.98320E-12 &
     2.00560E-12 2.02830E-12 2.05150E-12 2.07500E-12 2.09890E-12 2.12310E-12 2.14760E-12 2.17250E-12 &
     2.19770E-12 2.22310E-12 2.24890E-12 2.27490E-12 2.30120E-12 2.32770E-12 2.35440E-12 2.38140E-12 &
     2.40860E-12 2.43590E-12 2.46350E-12 2.49120E-12 2.51910E-12 2.54710E-12 2.57530E-12 2.60350E-12 &
     2.63190E-12 2.66030E-12 2.68890E-12 2.71750E-12 2.74610E-12 2.77480E-12 2.80350E-12 2.83220E-12 &
     2.86090E-12 2.88960E-12 2.91830E-12 2.94690E-12 2.97550E-12 3.00410E-12 3.03250E-12 3.06090E-12 &
     3.08910E-12 3.11730E-12 3.14530E-12 3.17310E-12 3.20080E-12 3.22840E-12 3.25580E-12 3.28290E-12 &
     3.30990E-12 3.33670E-12 3.36320E-12 3.38950E-12 3.41550E-12 3.44130E-12 3.46680E-12 3.49200E-12 &
     3.51690E-12 3.54150E-12 3.56580E-12 3.58980E-12 3.61340E-12 3.63670E-12 3.65960E-12 3.68220E-12 &
     3.70440E-12 3.72630E-12 3.74770E-12 3.76870E-12 3.78940E-12 3.80960E-12 3.82940E-12 3.84870E-12 &
     3.86770E-12 3.88610E-12 3.90410E-12 3.92170E-12 3.93880E-12 3.95540E-12 3.97150E-12 3.98720E-12 &
     4.00230E-12 4.01700E-12 4.03110E-12 4.04480E-12 4.05790E-12 4.07050E-12 4.08250E-12 4.09400E-12 &
     4.10500E-12 4.11550E-12 4.12530E-12 4.13470E-12 4.14340E-12 4.15170E-12 4.15930E-12 4.16640E-12 &
     4.17290E-12 4.17890E-12 4.18430E-12 4.18910E-12 4.19340E-12 4.19710E-12 4.20020E-12 4.20270E-12 &
     4.20460E-12 4.20600E-12 4.20680E-12 4.20700E-12 4.20660E-12 4.20560E-12 4.20410E-12 4.20200E-12 &
     4.19930E-12 4.19600E-12 4.19220E-12 4.18780E-12 4.18280E-12 4.17730E-12 4.17120E-12 4.16460E-12 &
     4.15740E-12 4.14970E-12 4.14140E-12 4.13250E-12 4.12310E-12 4.11320E-12 4.10270E-12 4.09170E-12 &
     4.08020E-12 4.06820E-12 4.05570E-12 4.04270E-12 4.02920E-12 4.01520E-12 4.00070E-12 3.98570E-12 &
     3.97030E-12 3.95440E-12 3.93810E-12 3.92130E-12 3.90400E-12 3.88630E-12 3.86820E-12 3.84960E-12 &
     3.83060E-12 3.81130E-12 3.79150E-12 3.77130E-12 3.75070E-12 3.72980E-12 3.70850E-12 3.68680E-12 &
     3.66470E-12 3.64230E-12 3.61960E-12 3.59650E-12 3.57310E-12 3.54950E-12 3.52550E-12 3.50120E-12 &
     3.47670E-12 3.45190E-12 3.42680E-12 3.40150E-12 3.37600E-12 3.35020E-12 3.32420E-12 3.29800E-12 &
     3.27160E-12 3.24500E-12 3.21830E-12 3.19130E-12 3.16430E-12 3.13700E-12 3.10970E-12 3.08220E-12 &
     3.05460E-12 3.02690E-12 2.99910E-12 2.97120E-12 2.94330E-12 2.91530E-12 2.88720E-12 2.85920E-12 &
     2.83100E-12 2.80290E-12 2.77480E-12 2.74660E-12 2.71850E-12 2.69030E-12 2.66220E-12 2.63420E-12 &
     2.60620E-12 2.57830E-12 2.55050E-12 2.52270E-12 2.49500E-12 2.46750E-12 2.44010E-12 2.41270E-12 &
     2.38560E-12 2.35850E-12 2.33160E-12 2.30490E-12 2.27840E-12 2.25200E-12 2.22580E-12 2.19980E-12 &
     2.17410E-12 2.14850E-12 2.12320E-12 2.09810E-12 2.07320E-12 2.04860E-12 2.02420E-12 2.00010E-12 &
     1.97630E-12 1.95280E-12 1.92960E-12 1.90660E-12 1.88400E-12 1.86160E-12 1.83960E-12 1.81790E-12 &
     1.79660E-12 1.77560E-12 1.75490E-12 1.73460E-12 1.71470E-12 1.69510E-12 1.67590E-12 1.65710E-12 &
     1.63860E-12 1.62060E-12 1.60290E-12 1.58560E-12 1.56880E-12 1.55230E-12 1.53630E-12 1.52060E-12 &
     1.50540E-12 1.49070E-12 1.47630E-12 1.46240E-12 1.44890E-12 1.43580E-12 1.42320E-12 1.41110E-12 &
     1.39930E-12 1.38810E-12 1.37720E-12 1.36690E-12 1.35700E-12 1.34750E-12 1.33850E-12 1.33000E-12 &
     1.32190E-12 1.31430E-12 1.30710E-12 1.30040E-12 1.29420E-12 1.28840E-12 1.28310E-12 1.27820E-12 &
     1.27390E-12 1.26990E-12 1.26650E-12 1.26350E-12 1.26090E-12 1.25880E-12 1.25720E-12 1.25600E-12 &
     1.25530E-12 1.25500E-12 1.25510E-12 1.25570E-12 1.25680E-12 1.25830E-12 1.26020E-12 1.26250E-12 &
     1.26530E-12 1.26850E-12 1.27210E-12 1.27620E-12 1.28070E-12 1.28550E-12 1.29080E-12 1.29650E-12 &
     1.30260E-12 1.30900E-12 1.31590E-12 1.32310E-12 1.33080E-12 1.33880E-12 1.34710E-12 1.35580E-12 &
     1.36490E-12 1.37430E-12 1.38410E-12 1.39420E-12 1.40470E-12 1.41540E-12 1.42650E-12 1.43790E-12 &
     1.44960E-12 1.46160E-12 1.47390E-12 1.48650E-12 1.49940E-12 1.51250E-12 1.52590E-12 1.53960E-12 &
     1.55350E-12 1.56760E-12 1.58200E-12 1.59660E-12 1.61140E-12 1.62650E-12 1.64170E-12 1.65720E-12 &
     1.67280E-12 1.68860E-12 1.70460E-12 1.72070E-12 1.73700E-12 1.75350E-12 1.77010E-12 1.78680E-12 &
     1.80370E-12 1.82070E-12 1.83780E-12 1.85490E-12 1.87220E-12 1.88960E-12 1.90700E-12 1.92450E-12 &
     1.94210E-12 1.95970E-12 1.97740E-12 1.99510E-12 2.01280E-12 2.03050E-12 2.04830E-12 2.06600E-12 &
     2.08370E-12 2.10150E-12 2.11920E-12 2.13680E-12 2.15450E-12 2.17200E-12 2.18950E-12 2.20700E-12 &
     2.22440E-12 2.24170E-12 2.25890E-12 2.27600E-12 2.29300E-12 2.30990E-12 2.32670E-12 2.34340E-12 &
     2.35990E-12 2.37630E-12 2.39260E-12 2.40860E-12 2.42460E-12 2.44030E-12 2.45590E-12 2.47130E-12 &
     2.48660E-12 2.50160E-12 2.51640E-12 2.53100E-12 2.54550E-12 2.55960E-12 2.57360E-12 2.58730E-12 &
     2.60080E-12 2.61410E-12 2.62710E-12 2.63990E-12 2.65230E-12 2.66460E-12 2.67660E-12 2.68820E-12 &
     2.69970E-12 2.71080E-12 2.72170E-12 2.73220E-12 2.74250E-12 2.75250E-12 2.76210E-12 2.77150E-12 &
     2.78060E-12 2.78930E-12 2.79770E-12 2.80580E-12 2.81360E-12 2.82100E-12 2.82810E-12 2.83490E-12 &
     2.84130E-12 2.84740E-12 2.85310E-12 2.85850E-12 2.86360E-12 2.86830E-12 2.87260E-12 2.87660E-12 &
     2.88020E-12 2.88350E-12 2.88640E-12 2.88900E-12 2.89120E-12 2.89300E-12 2.89450E-12 2.89560E-12 &
     2.89630E-12 2.89670E-12 2.89670E-12 2.89640E-12 2.89570E-12 2.89470E-12 2.89330E-12 2.89150E-12 &
     2.88940E-12 2.88690E-12 2.88410E-12 2.88090E-12 2.87740E-12 2.87350E-12 2.86930E-12 2.86470E-12 &
     2.85980E-12 2.85450E-12 2.84890E-12 2.84300E-12 2.83670E-12 2.83010E-12 2.82310E-12 2.81590E-12 &
     2.80820E-12 2.80030E-12 2.79210E-12 2.78350E-12 2.77460E-12 2.76540E-12 2.75590E-12 2.74620E-12 &
     2.73610E-12 2.72570E-12 2.71500E-12 2.70410E-12 2.69280E-12 2.68130E-12 2.66950E-12 2.65750E-12 &
     2.64510E-12 2.63250E-12 2.61970E-12 2.60660E-12 2.59320E-12 2.57970E-12 2.56580E-12 2.55180E-12 &
     2.53750E-12 2.52300E-12 2.50830E-12 2.49340E-12 2.47820E-12 2.46290E-12 2.44740E-12 2.43170E-12 &
     2.41580E-12 2.39970E-12 2.38350E-12 2.36710E-12 2.35050E-12 2.33380E-12 2.31700E-12 2.30000E-12 &
     2.28280E-12 2.26560E-12 2.24820E-12 2.23080E-12 2.21320E-12 2.19550E-12 2.17770E-12 2.15990E-12 &
     2.14190E-12 2.12390E-12 2.10580E-12 2.08760E-12 2.06940E-12 2.05120E-12 2.03280E-12 2.01450E-12 &
     1.99610E-12 1.97770E-12 1.95930E-12 1.94080E-12 1.92240E-12 1.90390E-12 1.88540E-12 1.86700E-12 &
     1.84850E-12 1.83010E-12 1.81170E-12 1.79340E-12 1.77500E-12 1.75680E-12 1.73860E-12 1.72040E-12 &
     1.70230E-12 1.68430E-12 1.66630E-12 1.64850E-12 1.63070E-12 1.61300E-12 1.59540E-12 1.57790E-12 &
     1.56050E-12 1.54320E-12 1.52610E-12 1.50900E-12 1.49210E-12 1.47540E-12 1.45870E-12 1.44220E-12 &
     1.42590E-12 1.40970E-12 1.39370E-12 1.37780E-12 1.36210E-12 1.34660E-12 1.33120E-12 1.31600E-12 &
     1.30100E-12 1.28620E-12 1.27160E-12 1.25720E-12 1.24290E-12 1.22890E-12 1.21510E-12 1.20150E-12 &
     1.18810E-12 1.17490E-12 1.16200E-12 1.14930E-12 1.13680E-12 1.12450E-12 1.11250E-12 1.10070E-12 &
     1.08910E-12 1.07780E-12 1.06670E-12 1.05580E-12 1.04520E-12 1.03490E-12 1.02480E-12 1.01500E-12 &
     1.00540E-12 9.96050E-13 9.86970E-13 9.78150E-13 9.69590E-13 9.61300E-13 9.53260E-13 9.45490E-13 &
     9.37980E-13 9.30730E-13 9.23750E-13 9.17030E-13 9.10590E-13 9.04410E-13 8.98490E-13 8.92850E-13 &
     8.87470E-13 8.82350E-13 8.77510E-13 8.72930E-13 8.68610E-13 8.64570E-13 8.60780E-13 8.57260E-13 &
     8.54010E-13 8.51010E-13 8.48280E-13 8.45810E-13 8.43600E-13 8.41650E-13 8.39950E-13 8.38510E-13 &
     8.37320E-13 8.36390E-13 8.35710E-13 8.35270E-13 8.35080E-13 8.35140E-13 8.35430E-13 8.35970E-13 &
     8.36750E-13 8.37760E-13 8.39000E-13 8.40480E-13 8.42180E-13 8.44110E-13 8.46260E-13 8.48630E-13 &
     8.51220E-13 8.54020E-13 8.57040E-13 8.60260E-13 8.63680E-13 8.67310E-13 8.71130E-13 8.75150E-13 &
     8.79360E-13 8.83760E-13 8.88340E-13 8.93100E-13 8.98030E-13 9.03140E-13 9.08420E-13 9.13870E-13 &
     9.19470E-13 9.25240E-13 9.31150E-13 9.37220E-13 9.43430E-13 9.49780E-13 9.56270E-13 9.62890E-13 &
     9.69640E-13 9.76510E-13 9.83510E-13 9.90620E-13 9.97840E-13 1.00520E-12 1.01260E-12 1.02010E-12 &
     1.02780E-12 1.03550E-12 1.04330E-12 1.05120E-12 1.05910E-12 1.06720E-12 1.07530E-12 1.08340E-12 &
     1.09160E-12 1.09990E-12 1.10820E-12 1.11660E-12 1.12500E-12 1.13350E-12 1.14190E-12 1.15040E-12 &
     1.15900E-12 1.16750E-12 1.17610E-12 1.18460E-12 1.19320E-12 1.20180E-12 1.21030E-12 1.21890E-12 &
     1.22740E-12 1.23590E-12 1.24440E-12 1.25290E-12 1.26140E-12 1.26980E-12 1.27820E-12 1.28650E-12 &
     1.29480E-12 1.30300E-12 1.31120E-12 1.31930E-12 1.32740E-12 1.33540E-12 1.34330E-12 1.35120E-12 &
     1.35900E-12 1.36670E-12 1.37430E-12 1.38180E-12 1.38930E-12 1.39660E-12 1.40390E-12 1.41100E-12 &
     1.41810E-12 1.42500E-12 1.43190E-12 1.43860E-12 1.44520E-12 1.45160E-12 1.45800E-12 1.46420E-12 &
     1.47030E-12 1.47630E-12 1.48210E-12 1.48780E-12 1.49340E-12 1.49880E-12 1.50400E-12 1.50920E-12 &
     1.51410E-12 1.51900E-12 1.52360E-12 1.52810E-12 1.53250E-12 1.53670E-12 1.54070E-12 1.54460E-12 &
     1.54830E-12 1.55180E-12 1.55520E-12 1.55840E-12 1.56140E-12 1.56430E-12 1.56690E-12 1.56940E-12 &
     1.57180E-12 1.57390E-12 1.57590E-12 1.57770E-12 1.57930E-12 1.58080E-12 1.58200E-12 1.58310E-12 &
     1.58400E-12 1.58470E-12 1.58530E-12 1.58560E-12 1.58580E-12 1.58580E-12 1.58560E-12 1.58520E-12 &
     1.58470E-12 1.58390E-12 1.58300E-12 1.58190E-12 1.58070E-12 1.57920E-12 1.57760E-12 1.57580E-12 &
     1.57380E-12 1.57170E-12 1.56930E-12 1.56680E-12 1.56410E-12 1.56130E-12 1.55830E-12 1.55510E-12 &
     1.55170E-12 1.54810E-12 1.54440E-12 1.54060E-12 1.53650E-12 1.53230E-12 1.52790E-12 1.52340E-12 &
     1.51870E-12 1.51390E-12 1.50890E-12 1.50370E-12 1.49840E-12 1.49300E-12 1.48740E-12 1.48170E-12 &
     1.47580E-12 1.46980E-12 1.46360E-12 1.45730E-12 1.45090E-12 1.44430E-12 1.43770E-12 1.43090E-12 &
     1.42390E-12 1.41690E-12 1.40970E-12 1.40240E-12 1.39500E-12 1.38750E-12 1.37990E-12 1.37220E-12 &
     1.36440E-12 1.35650E-12 1.34850E-12 1.34040E-12 1.33220E-12 1.32390E-12 1.31550E-12 1.30710E-12 &
     1.29860E-12 1.29000E-12 1.28140E-12 1.27260E-12 1.26380E-12 1.25500E-12 1.24610E-12 1.23710E-12 &
     1.22810E-12 1.21900E-12 1.20990E-12 1.20070E-12 1.19150E-12 1.18230E-12 1.17300E-12 1.16370E-12 &
     1.15440E-12 1.14500E-12 1.13570E-12 1.12630E-12 1.11680E-12 1.10740E-12 1.09800E-12 1.08850E-12 &
     1.07910E-12 1.06960E-12 1.06020E-12 1.05070E-12 1.04130E-12 1.03190E-12 1.02250E-12 1.01310E-12 &
     1.00370E-12 9.94390E-13 9.85080E-13 9.75800E-13 9.66550E-13 9.57340E-13 9.48160E-13 9.39030E-13 &
     9.29950E-13 9.20900E-13 9.11910E-13 9.02980E-13 8.94090E-13 8.85270E-13 8.76510E-13 8.67810E-13 &
     8.59170E-13 8.50610E-13 8.42110E-13 8.33690E-13 8.25340E-13 8.17070E-13 8.08880E-13 8.00770E-13 &
     7.92740E-13 7.84800E-13 7.76950E-13 7.69190E-13 7.61530E-13 7.53960E-13 7.46480E-13 7.39110E-13 &
     7.31840E-13 7.24670E-13 7.17600E-13 7.10640E-13 7.03780E-13 6.97040E-13 6.90400E-13 6.83880E-13 &
     6.77470E-13 6.71170E-13 6.64990E-13 6.58930E-13 6.52990E-13 6.47170E-13 6.41470E-13 6.35890E-13 &
     6.30430E-13 6.25100E-13 6.19890E-13 6.14810E-13 6.09860E-13 6.05030E-13 6.00330E-13 5.95760E-13 &
     5.91320E-13 5.87010E-13 5.82830E-13 5.78780E-13 5.74860E-13 5.71080E-13 5.67420E-13 5.63900E-13 &
     5.60520E-13 5.57260E-13 5.54140E-13 5.51150E-13 5.48290E-13 5.45570E-13 5.42980E-13 5.40530E-13 &
     5.38200E-13 5.36010E-13 5.33950E-13 5.32020E-13 5.30220E-13 5.28550E-13 5.27010E-13 5.25600E-13 &
     5.24330E-13 5.23170E-13 5.22150E-13 5.21250E-13 5.20480E-13 5.19830E-13 5.19310E-13 5.18910E-13 &
     5.18640E-13 5.18480E-13 5.18450E-13 5.18530E-13 5.18730E-13 5.19050E-13 5.19480E-13 5.20030E-13 &
     5.20690E-13 5.21470E-13 5.22350E-13 5.23340E-13 5.24440E-13 5.25640E-13 5.26950E-13 5.28360E-13 &
     5.29870E-13 5.31490E-13 5.33190E-13 5.35000E-13 5.36900E-13 5.38890E-13 5.40970E-13 5.43150E-13 &
     5.45400E-13 5.47750E-13 5.50180E-13 5.52690E-13 5.55270E-13 5.57940E-13 5.60680E-13 5.63500E-13 &
     5.66390E-13 5.69340E-13 5.72370E-13 5.75460E-13 5.78610E-13 5.81830E-13 5.85100E-13 5.88430E-13 &
     5.91820E-13 5.95260E-13 5.98750E-13 6.02290E-13 6.05870E-13 6.09500E-13 6.13170E-13 6.16880E-13 &
     6.20630E-13 6.24410E-13 6.28230E-13 6.32080E-13 6.35950E-13 6.39860E-13 6.43780E-13 6.47730E-13 &
     6.51700E-13 6.55690E-13 6.59690E-13 6.63700E-13 6.67730E-13 6.71760E-13 6.75810E-13 6.79850E-13 &
     6.83900E-13 6.87950E-13 6.92000E-13 6.96050E-13 7.00080E-13 7.04120E-13 7.08140E-13 7.12150E-13 &
     7.16150E-13 7.20130E-13 7.24090E-13 7.28030E-13 7.31950E-13 7.35840E-13 7.39710E-13 7.43550E-13 &
     7.47370E-13 7.51150E-13 7.54890E-13 7.58610E-13 7.62280E-13 7.65910E-13 7.69510E-13 7.73060E-13 &
     7.76560E-13 7.80020E-13 7.83430E-13 7.86790E-13 7.90100E-13 7.93360E-13 7.96560E-13 7.99710E-13 &
     8.02800E-13 8.05830E-13 8.08800E-13 8.11710E-13 8.14550E-13 8.17330E-13 8.20050E-13 8.22700E-13 &
     8.25270E-13 8.27780E-13 8.30220E-13 8.32590E-13 8.34880E-13 8.37100E-13 8.39240E-13 8.41310E-13 &
     8.43300E-13 8.45210E-13 8.47050E-13 8.48800E-13 8.50480E-13 8.52080E-13 8.53590E-13 8.55020E-13 &
     8.56370E-13 8.57630E-13 8.58800E-13 8.59890E-13 8.60900E-13 8.61820E-13 8.62650E-13 8.63400E-13 &
     8.64060E-13 8.64630E-13 8.65120E-13 8.65520E-13 8.65830E-13 8.66050E-13 8.66180E-13 8.66230E-13 &
     8.66180E-13 8.66050E-13 8.65830E-13 8.65520E-13 8.65120E-13 8.64630E-13 8.64050E-13 8.63380E-13 &
     8.62630E-13 8.61780E-13 8.60850E-13 8.59830E-13 8.58720E-13 8.57520E-13 8.56230E-13 8.54860E-13 &
     8.53400E-13 8.51850E-13 8.50220E-13 8.48500E-13 8.46700E-13 8.44810E-13 8.42850E-13 8.40800E-13 &
     8.38670E-13 8.36460E-13 8.34170E-13 8.31800E-13 8.29360E-13 8.26830E-13 8.24230E-13 8.21560E-13 &
     8.18810E-13 8.15980E-13 8.13090E-13 8.10120E-13 8.07090E-13 8.03980E-13 8.00810E-13 7.97560E-13 &
     7.94260E-13 7.90890E-13 7.87450E-13 7.83950E-13 7.80400E-13 7.76780E-13 7.73100E-13 7.69360E-13 &
     7.65570E-13 7.61720E-13 7.57810E-13 7.53860E-13 7.49850E-13 7.45790E-13 7.41690E-13 7.37530E-13 &
     7.33330E-13 7.29090E-13 7.24810E-13 7.20480E-13 7.16120E-13 7.11710E-13 7.07270E-13 7.02790E-13 &
     6.98280E-13 6.93730E-13 6.89150E-13 6.84540E-13 6.79900E-13 6.75240E-13 6.70550E-13 6.65830E-13 &
     6.61090E-13 6.56330E-13 6.51550E-13 6.46750E-13 6.41940E-13 6.37100E-13 6.32260E-13 6.27400E-13 &
     6.22530E-13 6.17650E-13 6.12770E-13 6.07870E-13 6.02980E-13 5.98080E-13 5.93180E-13 5.88280E-13 &
     5.83380E-13 5.78490E-13 5.73600E-13 5.68710E-13 5.63840E-13 5.58970E-13 5.54110E-13 5.49260E-13 &
     5.44420E-13 5.39600E-13 5.34790E-13 5.30010E-13 5.25230E-13 5.20480E-13 5.15750E-13 5.11040E-13 &
     5.06360E-13 5.01690E-13 4.97060E-13 4.92450E-13 4.87870E-13 4.83320E-13 4.78800E-13 4.74310E-13 &
     4.69860E-13 4.65440E-13 4.61050E-13 4.56700E-13 4.52390E-13 4.48120E-13 4.43890E-13 4.39700E-13 &
     4.35540E-13 4.31440E-13 4.27370E-13 4.23350E-13 4.19370E-13 4.15450E-13 4.11560E-13 4.07730E-13 &
     4.03940E-13 4.00210E-13 3.96520E-13 3.92890E-13 3.89310E-13 3.85780E-13 3.82300E-13 3.78880E-13 &
     3.75510E-13 3.72190E-13 3.68940E-13 3.65730E-13 3.62590E-13 3.59500E-13 3.56470E-13 3.53500E-13 &
     3.50590E-13 3.47730E-13 3.44940E-13 3.42210E-13 3.39540E-13 3.36920E-13 3.34370E-13 3.31880E-13 &
     3.29460E-13 3.27090E-13 3.24790E-13 3.22540E-13 3.20370E-13 3.18250E-13 3.16200E-13 3.14200E-13 &
     3.12280E-13 3.10410E-13 3.08610E-13 3.06870E-13 3.05190E-13 3.03570E-13 3.02020E-13 3.00530E-13 &
     2.99100E-13 2.97740E-13 2.96430E-13 2.95190E-13 2.94010E-13 2.92890E-13 2.91830E-13 2.90830E-13 &
     2.89900E-13 2.89020E-13 2.88200E-13 2.87440E-13 2.86740E-13 2.86100E-13 2.85520E-13 2.84990E-13 &
     2.84530E-13 2.84110E-13 2.83760E-13 2.83460E-13 2.83210E-13 2.83020E-13 2.82880E-13 2.82800E-13 &
     2.82770E-13 2.82790E-13 2.82860E-13 2.82990E-13 2.83160E-13 2.83380E-13 2.83650E-13 2.83970E-13 &
     2.84340E-13 2.84750E-13 2.85210E-13 2.85710E-13 2.86250E-13 2.86840E-13 2.87480E-13 2.88150E-13 &
     2.88860E-13 2.89620E-13 2.90410E-13 2.91240E-13 2.92110E-13 2.93010E-13 2.93950E-13 2.94930E-13 &
     2.95930E-13 2.96970E-13 2.98050E-13 2.99150E-13 3.00280E-13 3.01440E-13 3.02630E-13 3.03850E-13 &
     3.05090E-13 3.06360E-13 3.07650E-13 3.08960E-13 3.10300E-13 3.11660E-13 3.13040E-13 3.14430E-13 &
     3.15850E-13 3.17280E-13 3.18730E-13 3.20190E-13 3.21670E-13 3.23160E-13 3.24660E-13 3.26180E-13 &
     3.27700E-13 3.29240E-13 3.30780E-13 3.32330E-13 3.33880E-13 3.35450E-13 3.37010E-13 3.38580E-13 &
     3.40150E-13 3.41730E-13 3.43300E-13 3.44870E-13 3.46450E-13 3.48020E-13 3.49580E-13 3.51150E-13 &
     3.52700E-13 3.54250E-13 3.55800E-13 3.57330E-13 3.58860E-13 3.60380E-13 3.61890E-13 3.63380E-13 &
     3.64870E-13 3.66340E-13 3.67790E-13 3.69240E-13 3.70660E-13 3.72070E-13 3.73470E-13 3.74840E-13 &
     3.76200E-13 3.77540E-13 3.78860E-13 3.80160E-13 3.81430E-13 3.82690E-13 3.83920E-13 3.85130E-13 &
     3.86310E-13 3.87470E-13 3.88600E-13 3.89710E-13 3.90790E-13 3.91850E-13 3.92880E-13 3.93880E-13 &
     3.94850E-13 3.95790E-13 3.96700E-13 3.97580E-13 3.98430E-13 3.99250E-13 4.00040E-13 4.00800E-13 &
     4.01520E-13 4.02210E-13 4.02870E-13 4.03490E-13 4.04080E-13 4.04640E-13 4.05160E-13 4.05650E-13 &
     4.06100E-13 4.06520E-13 4.06900E-13 4.07250E-13 4.07560E-13 4.07830E-13 4.08070E-13 4.08270E-13 &
     4.08430E-13 4.08550E-13 4.08640E-13 4.08690E-13 4.08710E-13 4.08680E-13 4.08620E-13 4.08520E-13 &
     4.08380E-13 4.08210E-13 4.07990E-13 4.07740E-13 4.07450E-13 4.07120E-13 4.06750E-13 4.06350E-13 &
     4.05910E-13 4.05430E-13 4.04920E-13 4.04370E-13 4.03780E-13 4.03160E-13 4.02500E-13 4.01800E-13 &
     4.01070E-13 4.00300E-13 3.99500E-13 3.98660E-13 3.97780E-13 3.96880E-13 3.95930E-13 3.94960E-13 &
     3.93950E-13 3.92900E-13 3.91820E-13 3.90720E-13 3.89570E-13 3.88400E-13 3.87200E-13 3.85960E-13 &
     3.84690E-13 3.83390E-13 3.82070E-13 3.80710E-13 3.79320E-13 3.77910E-13 3.76460E-13 3.74990E-13 &
     3.73490E-13 3.71960E-13 3.70410E-13 3.68830E-13 3.67230E-13 3.65600E-13 3.63950E-13 3.62270E-13 &
     3.60570E-13 3.58850E-13 3.57110E-13 3.55340E-13 3.53560E-13 3.51750E-13 3.49920E-13 3.48080E-13 &
     3.46210E-13 3.44330E-13 3.42430E-13 3.40510E-13 3.38580E-13 3.36630E-13 3.34670E-13 3.32690E-13 &
     3.30700E-13 3.28690E-13 3.26680E-13 3.24650E-13 3.22610E-13 3.20560E-13 3.18500E-13 3.16430E-13 &
     3.14350E-13 3.12260E-13 3.10170E-13 3.08070E-13 3.05970E-13 3.03860E-13 3.01750E-13 2.99630E-13 &
     2.97510E-13 2.95380E-13 2.93260E-13 2.91130E-13 2.89000E-13 2.86870E-13 2.84740E-13 2.82620E-13 &
     2.80490E-13 2.78370E-13 2.76250E-13 2.74130E-13 2.72020E-13 2.69910E-13 2.67800E-13 2.65710E-13 &
     2.63610E-13 2.61530E-13 2.59450E-13 2.57380E-13 2.55320E-13 2.53270E-13 2.51230E-13 2.49200E-13 &
     2.47180E-13 2.45170E-13 2.43180E-13 2.41190E-13 2.39220E-13 2.37260E-13 2.35320E-13 2.33390E-13 &
     2.31470E-13 2.29570E-13 2.27690E-13 2.25820E-13 2.23970E-13 2.22130E-13 2.20320E-13 2.18520E-13 &
     2.16740E-13 2.14980E-13 2.13230E-13 2.11510E-13 2.09810E-13 2.08130E-13 2.06460E-13 2.04820E-13 &
     2.03200E-13 2.01600E-13 2.00030E-13 1.98470E-13 1.96940E-13 1.95430E-13 1.93950E-13 1.92490E-13 &
     1.91050E-13 1.89630E-13 1.88240E-13 1.86880E-13 1.85540E-13 1.84220E-13 1.82930E-13 1.81660E-13 &
     1.80420E-13 1.79210E-13 1.78020E-13 1.76860E-13 1.75720E-13 1.74610E-13 1.73520E-13 1.72460E-13 &
     1.71430E-13 1.70420E-13 1.69440E-13 1.68490E-13 1.67560E-13 1.66660E-13 1.65780E-13 1.64940E-13 &
     1.64120E-13 1.63320E-13 1.62550E-13 1.61810E-13 1.61100E-13 1.60410E-13 1.59750E-13 1.59110E-13 &
     1.58500E-13 1.57920E-13 1.57360E-13 1.56830E-13 1.56330E-13 1.55850E-13 1.55400E-13 1.54970E-13 &
     1.54570E-13 1.54190E-13 1.53840E-13 1.53510E-13 1.53200E-13 1.52930E-13 1.52670E-13 1.52440E-13 &
     1.52240E-13 1.52050E-13 1.51890E-13 1.51760E-13 1.51640E-13 1.51550E-13 1.51490E-13 1.51440E-13 &
     1.51420E-13 1.51410E-13 1.51430E-13 1.51470E-13 1.51530E-13 1.51610E-13 1.51710E-13 1.51840E-13 &
     1.51980E-13 1.52130E-13 1.52310E-13 1.52510E-13 1.52720E-13 1.52960E-13 1.53200E-13 1.53470E-13 &
     1.53750E-13 1.54050E-13 1.54370E-13 1.54700E-13 1.55050E-13 1.55410E-13 1.55780E-13 1.56170E-13 &
     1.56570E-13 1.56990E-13 1.57420E-13 1.57860E-13 1.58310E-13 1.58780E-13 1.59250E-13 1.59740E-13 &
     1.60230E-13 1.60740E-13 1.61260E-13 1.61780E-13 1.62310E-13 1.62860E-13 1.63410E-13 1.63960E-13 &
     1.64530E-13 1.65100E-13 1.65680E-13 1.66260E-13 1.66850E-13 1.67440E-13 1.68040E-13 1.68640E-13 &
     1.69240E-13 1.69850E-13 1.70460E-13 1.71080E-13 1.71690E-13 1.72310E-13 1.72930E-13 1.73550E-13 &
     1.74160E-13 1.74780E-13 1.75400E-13 1.76020E-13 1.76630E-13 1.77250E-13 1.77860E-13 1.78470E-13 &
     1.79080E-13 1.79680E-13 1.80280E-13 1.80880E-13 1.81470E-13 1.82050E-13 1.82640E-13 1.83210E-13 &
     1.83780E-13 1.84350E-13 1.84910E-13 1.85460E-13 1.86010E-13 1.86540E-13 1.87070E-13 1.87590E-13 &
     1.88110E-13 1.88610E-13 1.89110E-13 1.89600E-13 1.90070E-13 1.90540E-13 1.91000E-13 1.91440E-13 &
     1.91880E-13 1.92310E-13 1.92720E-13 1.93120E-13 1.93520E-13 1.93900E-13 1.94260E-13 1.94620E-13 &
     1.94960E-13 1.95290E-13 1.95610E-13 1.95910E-13 1.96200E-13 1.96480E-13 1.96740E-13 1.96990E-13 &
     1.97230E-13 1.97450E-13 1.97660E-13 1.97850E-13 1.98020E-13 1.98190E-13 1.98330E-13 1.98470E-13 &
     1.98580E-13 1.98680E-13 1.98770E-13 1.98840E-13 1.98890E-13 1.98930E-13 1.98960E-13 1.98960E-13 &
     1.98950E-13 1.98930E-13 1.98890E-13 1.98830E-13 1.98760E-13 1.98670E-13 1.98560E-13 1.98440E-13 &
     1.98310E-13 1.98160E-13 1.97990E-13 1.97800E-13 1.97600E-13 1.97390E-13 1.97150E-13 1.96910E-13 &
     1.96640E-13 1.96360E-13 1.96070E-13 1.95760E-13 1.95430E-13 1.95090E-13 1.94730E-13 1.94360E-13 &
     1.93970E-13 1.93570E-13 1.93150E-13 1.92710E-13 1.92270E-13 1.91800E-13 1.91320E-13 1.90830E-13 &
     1.90330E-13 1.89810E-13 1.89270E-13 1.88720E-13 1.88160E-13 1.87580E-13 1.87000E-13 1.86390E-13 &
     1.85780E-13 1.85150E-13 1.84510E-13 1.83850E-13 1.83190E-13 1.82510E-13 1.81820E-13 1.81110E-13 &
     1.80400E-13 1.79680E-13 1.78940E-13 1.78190E-13 1.77430E-13 1.76670E-13 1.75890E-13 1.75100E-13 &
     1.74300E-13 1.73490E-13 1.72670E-13 1.71850E-13 1.71010E-13 1.70170E-13 1.69310E-13 1.68450E-13 &
     1.67590E-13 1.66710E-13 1.65830E-13 1.64940E-13 1.64040E-13 1.63140E-13 1.62230E-13 1.61320E-13 &
     1.60400E-13 1.59470E-13 1.58540E-13 1.57600E-13 1.56660E-13 1.55720E-13 1.54770E-13 1.53820E-13 &
     1.52860E-13 1.51900E-13 1.50940E-13 1.49970E-13 1.49010E-13 1.48040E-13 1.47060E-13 1.46090E-13 &
     1.45110E-13 1.44140E-13 1.43160E-13 1.42180E-13 1.41200E-13 1.40220E-13 1.39240E-13 1.38260E-13 &
     1.37290E-13 1.36310E-13 1.35330E-13 1.34360E-13 1.33380E-13 1.32410E-13 1.31440E-13 1.30480E-13 &
     1.29510E-13 1.28550E-13 1.27590E-13 1.26640E-13 1.25690E-13 1.24740E-13 1.23800E-13 1.22860E-13 &
     1.21920E-13 1.21000E-13 1.20070E-13 1.19150E-13 1.18240E-13 1.17330E-13 1.16420E-13 1.15530E-13 &
     1.14640E-13 1.13750E-13 1.12870E-13 1.12000E-13 1.11140E-13 1.10280E-13 1.09430E-13 1.08590E-13 &
     1.07750E-13 1.06930E-13 1.06110E-13 1.05300E-13 1.04490E-13 1.03700E-13 1.02920E-13 1.02140E-13 &
     1.01370E-13 1.00610E-13 9.98640E-14 9.91230E-14 9.83930E-14 9.76720E-14 9.69600E-14 9.62580E-14 &
     9.55660E-14 9.48840E-14 9.42120E-14 9.35500E-14 9.28980E-14 9.22570E-14 9.16260E-14 9.10050E-14 &
     9.03950E-14 8.97960E-14 8.92080E-14 8.86300E-14 8.80630E-14 8.75060E-14 8.69610E-14 8.64270E-14 &
     8.59030E-14 8.53910E-14 8.48900E-14 8.43990E-14 8.39200E-14 8.34520E-14 8.29950E-14 8.25490E-14 &
     8.21150E-14 8.16910E-14 8.12790E-14 8.08780E-14 8.04870E-14 8.01080E-14 7.97410E-14 7.93840E-14 &
     7.90380E-14 7.87030E-14 7.83790E-14 7.80660E-14 7.77640E-14 7.74730E-14 7.71930E-14 7.69230E-14 &
     7.66640E-14 7.64160E-14 7.61780E-14 7.59510E-14 7.57340E-14 7.55270E-14 7.53310E-14 7.51450E-14 &
     7.49700E-14 7.48040E-14 7.46480E-14 7.45020E-14 7.43660E-14 7.42390E-14 7.41220E-14 7.40140E-14 &
     7.39160E-14 7.38260E-14 7.37460E-14 7.36750E-14 7.36130E-14 7.35590E-14 7.35140E-14 7.34780E-14 &
     7.34500E-14 7.34300E-14 7.34180E-14 7.34140E-14 7.34180E-14 7.34300E-14 7.34490E-14 7.34750E-14 &
     7.35090E-14 7.35500E-14 7.35980E-14 7.36520E-14 7.37130E-14 7.37810E-14 7.38550E-14 7.39350E-14 &
     7.40220E-14 7.41140E-14 7.42120E-14 7.43150E-14 7.44240E-14 7.45380E-14 7.46570E-14 7.47810E-14 &
     7.49100E-14 7.50430E-14 7.51800E-14 7.53220E-14 7.54680E-14 7.56180E-14 7.57710E-14 7.59280E-14 &
     7.60880E-14 7.62520E-14 7.64180E-14 7.65880E-14 7.67600E-14 7.69350E-14 7.71120E-14 7.72910E-14 &
     7.74720E-14 7.76560E-14 7.78410E-14 7.80270E-14 7.82150E-14 7.84040E-14 7.85940E-14 7.87850E-14 &
     7.89770E-14 7.91700E-14 7.93620E-14 7.95550E-14 7.97490E-14 7.99420E-14 8.01350E-14 8.03270E-14 &
     8.05190E-14 8.07110E-14 8.09010E-14 8.10910E-14 8.12800E-14 8.14670E-14 8.16530E-14 8.18370E-14 &
     8.20200E-14 8.22000E-14 8.23790E-14 8.25560E-14 8.27300E-14 8.29020E-14 8.30720E-14 8.32390E-14 &
     8.34030E-14 8.35650E-14 8.37230E-14 8.38780E-14 8.40300E-14 8.41790E-14 8.43240E-14 8.44660E-14 &
     8.46040E-14 8.47380E-14 8.48680E-14 8.49940E-14 8.51160E-14 8.52340E-14 8.53470E-14 8.54560E-14 &
     8.55610E-14 8.56610E-14 8.57560E-14 8.58460E-14 8.59320E-14 8.60120E-14 8.60880E-14 8.61580E-14 &
     8.62240E-14 8.62840E-14 8.63390E-14 8.63890E-14 8.64340E-14 8.64730E-14 8.65070E-14 8.65350E-14 &
     8.65580E-14 8.65750E-14 8.65870E-14 8.65930E-14 8.65930E-14 8.65870E-14 8.65760E-14 8.65590E-14 &
     8.65360E-14 8.65070E-14 8.64730E-14 8.64320E-14 8.63860E-14 8.63340E-14 8.62750E-14 8.62110E-14 &
     8.61410E-14 8.60650E-14 8.59820E-14 8.58940E-14 8.58000E-14 8.57000E-14 8.55940E-14 8.54830E-14 &
     8.53650E-14 8.52410E-14 8.51120E-14 8.49760E-14 8.48350E-14 8.46880E-14 8.45350E-14 8.43770E-14 &
     8.42120E-14 8.40420E-14 8.38660E-14 8.36850E-14 8.34980E-14 8.33050E-14 8.31070E-14 8.29040E-14 &
     8.26950E-14 8.24810E-14 8.22610E-14 8.20370E-14 8.18070E-14 8.15720E-14 8.13320E-14 8.10870E-14 &
     8.08370E-14 8.05820E-14 8.03220E-14 8.00580E-14 7.97890E-14 7.95150E-14 7.92370E-14 7.89550E-14 &
     7.86680E-14 7.83780E-14 7.80830E-14 7.77840E-14 7.74810E-14 7.71740E-14 7.68630E-14 7.65490E-14 &
     7.62310E-14 7.59090E-14 7.55840E-14 7.52560E-14 7.49240E-14 7.45890E-14 7.42510E-14 7.39110E-14 &
     7.35670E-14 7.32200E-14 7.28710E-14 7.25190E-14 7.21640E-14 7.18070E-14 7.14470E-14 7.10850E-14 &
     7.07210E-14 7.03550E-14 6.99870E-14 6.96170E-14 6.92460E-14 6.88730E-14 6.84980E-14 6.81220E-14 &
     6.77440E-14 6.73660E-14 6.69860E-14 6.66050E-14 6.62230E-14 6.58400E-14 6.54570E-14 6.50730E-14 &
     6.46880E-14 6.43030E-14 6.39180E-14 6.35320E-14 6.31470E-14 6.27610E-14 6.23750E-14 6.19900E-14 &
     6.16050E-14 6.12200E-14 6.08350E-14 6.04520E-14 6.00680E-14 5.96860E-14 5.93040E-14 5.89240E-14 &
     5.85440E-14 5.81660E-14 5.77880E-14 5.74130E-14 5.70380E-14 5.66650E-14 5.62940E-14 5.59240E-14 &
     5.55560E-14 5.51900E-14 5.48260E-14 5.44640E-14 5.41040E-14 5.37460E-14 5.33900E-14 5.30370E-14 &
     5.26860E-14 5.23370E-14 5.19910E-14 5.16480E-14 5.13070E-14 5.09690E-14 5.06340E-14 5.03010E-14 &
     4.99720E-14 4.96450E-14 4.93220E-14 4.90020E-14 4.86850E-14 4.83710E-14 4.80600E-14 4.77530E-14 &
     4.74490E-14 4.71490E-14 4.68520E-14 4.65580E-14 4.62680E-14 4.59820E-14 4.57000E-14 4.54210E-14 &
     4.51460E-14 4.48750E-14 4.46080E-14 4.43450E-14 4.40850E-14 4.38300E-14 4.35780E-14 4.33310E-14 &
     4.30870E-14 4.28480E-14 4.26130E-14 4.23810E-14 4.21540E-14 4.19320E-14 4.17130E-14 4.14980E-14 &
     4.12880E-14 4.10820E-14 4.08800E-14 4.06820E-14 4.04890E-14 4.02990E-14 4.01140E-14 3.99340E-14 &
     3.97570E-14 3.95850E-14 3.94170E-14 3.92540E-14 3.90940E-14 3.89390E-14 3.87880E-14 3.86410E-14 &
     3.84990E-14 3.83610E-14 3.82270E-14 3.80970E-14 3.79710E-14 3.78490E-14 3.77320E-14 3.76190E-14 &
     3.75090E-14 3.74040E-14 3.73030E-14 3.72060E-14 3.71120E-14 3.70230E-14 3.69380E-14 3.68560E-14 &
     3.67780E-14 3.67040E-14 3.66340E-14 3.65680E-14 3.65050E-14 3.64460E-14 3.63910E-14 3.63390E-14 &
     3.62910E-14 3.62460E-14 3.62050E-14 3.61670E-14 3.61320E-14 3.61010E-14 3.60730E-14 3.60480E-14 &
     3.60270E-14 3.60080E-14 3.59930E-14 3.59800E-14 3.59710E-14 3.59640E-14 3.59610E-14 3.59600E-14 &
     3.59620E-14 3.59660E-14 3.59730E-14 3.59830E-14 3.59950E-14 3.60100E-14 3.60270E-14 3.60470E-14 &
     3.60680E-14 3.60920E-14 3.61180E-14 3.61460E-14 3.61760E-14 3.62090E-14 3.62430E-14 3.62780E-14 &
     3.63160E-14 3.63550E-14 3.63960E-14 3.64390E-14 3.64830E-14 3.65290E-14 3.65760E-14 3.66240E-14 &
     3.66730E-14 3.67240E-14 3.67760E-14 3.68290E-14 3.68830E-14 3.69380E-14 3.69940E-14 3.70500E-14 &
     3.71080E-14 3.71660E-14 3.72240E-14 3.72830E-14 3.73430E-14 3.74030E-14 3.74640E-14 3.75250E-14 &
     3.75860E-14 3.76470E-14 3.77080E-14 3.77700E-14 3.78310E-14 3.78920E-14 3.79540E-14 3.80150E-14 &
     3.80760E-14 3.81360E-14 3.81960E-14 3.82560E-14 3.83150E-14 3.83740E-14 3.84320E-14 3.84900E-14 &
     3.85460E-14 3.86020E-14 3.86580E-14 3.87120E-14 3.87660E-14 3.88180E-14 3.88700E-14 3.89200E-14 &
     3.89690E-14 3.90170E-14 3.90640E-14 3.91100E-14 3.91540E-14 3.91970E-14 3.92390E-14 3.92790E-14 &
     3.93180E-14 3.93560E-14 3.93910E-14 3.94260E-14 3.94580E-14 3.94890E-14 3.95190E-14 3.95460E-14 &
     3.95720E-14 3.95960E-14 3.96190E-14 3.96390E-14 3.96580E-14 3.96740E-14 3.96890E-14 3.97020E-14 &
     3.97120E-14 3.97210E-14 3.97280E-14 3.97330E-14 3.97350E-14 3.97350E-14 3.97340E-14 3.97300E-14 &
     3.97240E-14 3.97160E-14 3.97050E-14 3.96930E-14 3.96780E-14 3.96610E-14 3.96410E-14 3.96200E-14 &
     3.95960E-14 3.95690E-14 3.95410E-14 3.95100E-14 3.94770E-14 3.94410E-14 3.94030E-14 3.93630E-14 &
     3.93210E-14 3.92760E-14 3.92290E-14 3.91790E-14 3.91270E-14 3.90730E-14 3.90160E-14 3.89570E-14 &
     3.88960E-14 3.88330E-14 3.87670E-14 3.86990E-14 3.86280E-14 3.85560E-14 3.84810E-14 3.84030E-14 &
     3.83240E-14 3.82420E-14 3.81580E-14 3.80720E-14 3.79840E-14 3.78930E-14 3.78010E-14 3.77060E-14 &
     3.76090E-14 3.75100E-14 3.74100E-14 3.73070E-14 3.72020E-14 3.70950E-14 3.69860E-14 3.68750E-14 &
     3.67620E-14 3.66480E-14 3.65310E-14 3.64130E-14 3.62930E-14 3.61710E-14 3.60470E-14 3.59220E-14 &
     3.57950E-14 3.56660E-14 3.55350E-14 3.54030E-14 3.52690E-14 3.51340E-14 3.49970E-14 3.48590E-14 &
     3.47190E-14 3.45780E-14 3.44350E-14 3.42910E-14 3.41460E-14 3.39990E-14 3.38520E-14 3.37030E-14 &
     3.35520E-14 3.34010E-14 3.32490E-14 3.30950E-14 3.29410E-14 3.27850E-14 3.26290E-14 3.24720E-14 &
     3.23130E-14 3.21540E-14 3.19950E-14 3.18340E-14 3.16730E-14 3.15110E-14 3.13480E-14 3.11850E-14 &
     3.10210E-14 3.08570E-14 3.06920E-14 3.05270E-14 3.03610E-14 3.01950E-14 3.00290E-14 2.98620E-14 &
     2.96950E-14 2.95280E-14 2.93610E-14 2.91940E-14 2.90270E-14 2.88590E-14 2.86920E-14 2.85250E-14 &
     2.83570E-14 2.81900E-14 2.80230E-14 2.78560E-14 2.76900E-14 2.75230E-14 2.73570E-14 2.71920E-14 &
     2.70260E-14 2.68610E-14 2.66970E-14 2.65330E-14 2.63690E-14 2.62060E-14 2.60440E-14 2.58820E-14 &
     2.57200E-14 2.55600E-14 2.54000E-14 2.52410E-14 2.50820E-14 2.49240E-14 2.47680E-14 2.46120E-14 &
     2.44570E-14 2.43020E-14 2.41490E-14 2.39970E-14 2.38460E-14 2.36960E-14 2.35460E-14 2.33980E-14 &
     2.32510E-14 2.31060E-14 2.29610E-14 2.28170E-14 2.26750E-14 2.25340E-14 2.23940E-14 2.22560E-14 &
     2.21190E-14 2.19830E-14 2.18480E-14 2.17150E-14 2.15830E-14 2.14530E-14 2.13240E-14 2.11960E-14 &
     2.10700E-14 2.09450E-14 2.08220E-14 2.07000E-14 2.05800E-14 2.04610E-14 2.03440E-14 2.02290E-14 &
     2.01150E-14 2.00020E-14 1.98910E-14 1.97820E-14 1.96740E-14 1.95680E-14 1.94640E-14 1.93610E-14 &
     1.92600E-14 1.91600E-14 1.90620E-14 1.89660E-14 1.88710E-14 1.87780E-14 1.86870E-14 1.85970E-14 &
     1.85100E-14 1.84230E-14 1.83380E-14 1.82550E-14 1.81740E-14 1.80940E-14 1.80160E-14 1.79400E-14 &
     1.78650E-14 1.77920E-14 1.77210E-14 1.76510E-14 1.75830E-14 1.75170E-14 1.74520E-14 1.73880E-14 &
     1.73270E-14 1.72670E-14 1.72080E-14 1.71520E-14 1.70960E-14 1.70430E-14 1.69900E-14 1.69400E-14 &
     1.68910E-14 1.68430E-14 1.67970E-14 1.67530E-14 1.67100E-14 1.66680E-14 1.66280E-14 1.65900E-14 &
     1.65530E-14 1.65170E-14 1.64820E-14 1.64490E-14 1.64180E-14 1.63880E-14 1.63590E-14 1.63310E-14 &
     1.63050E-14 1.62800E-14 1.62560E-14 1.62340E-14 1.62120E-14 1.61920E-14 1.61740E-14 1.61560E-14 &
     1.61390E-14 1.61240E-14 1.61100E-14 1.60960E-14 1.60840E-14 1.60730E-14 1.60630E-14 1.60540E-14 &
     1.60460E-14 1.60390E-14 1.60330E-14 1.60270E-14 1.60230E-14 1.60200E-14 1.60170E-14 1.60150E-14 &
     1.60140E-14 1.60140E-14 1.60140E-14 1.60150E-14 1.60170E-14 1.60200E-14 1.60230E-14 1.60270E-14 &
     1.60310E-14 1.60360E-14 1.60420E-14 1.60480E-14 1.60550E-14 1.60620E-14 1.60700E-14 1.60780E-14 &
     1.60860E-14 1.60950E-14 1.61040E-14 1.61140E-14 1.61240E-14 1.61340E-14 1.61450E-14 1.61550E-14 &
     1.61660E-14 1.61770E-14 1.61890E-14 1.62000E-14 1.62120E-14 1.62240E-14 1.62360E-14 1.62470E-14 &
     1.62590E-14 1.62710E-14 1.62830E-14 1.62950E-14 1.63070E-14 1.63190E-14 1.63310E-14 1.63430E-14 &
     1.63540E-14 1.63660E-14 1.63770E-14 1.63880E-14 1.63990E-14 1.64100E-14 1.64200E-14 1.64300E-14 &
     1.64400E-14 1.64500E-14 1.64590E-14 1.64680E-14 1.64770E-14 1.64850E-14 1.64930E-14 1.65000E-14 &
     1.65070E-14 1.65140E-14 1.65200E-14 1.65260E-14 1.65310E-14 1.65360E-14 1.65400E-14 1.65440E-14 &
     1.65470E-14 1.65500E-14 1.65520E-14 1.65540E-14 1.65550E-14 1.65550E-14 1.65550E-14 1.65540E-14 &
     1.65520E-14 1.65500E-14 1.65480E-14 1.65440E-14 1.65400E-14 1.65350E-14 1.65300E-14 1.65230E-14 &
     1.65170E-14 1.65090E-14 1.65010E-14 1.64920E-14 1.64820E-14 1.64710E-14 1.64600E-14 1.64480E-14 &
     1.64350E-14 1.64220E-14 1.64070E-14 1.63920E-14 1.63760E-14 1.63600E-14 1.63420E-14 1.63240E-14 &
     1.63050E-14 1.62860E-14 1.62650E-14 1.62440E-14 1.62220E-14 1.61990E-14 1.61760E-14 1.61510E-14 &
     1.61260E-14 1.61000E-14 1.60730E-14 1.60460E-14 1.60180E-14 1.59890E-14 1.59590E-14 1.59280E-14 &
     1.58970E-14 1.58650E-14 1.58320E-14 1.57990E-14 1.57640E-14 1.57290E-14 1.56940E-14 1.56570E-14 &
     1.56200E-14 1.55820E-14 1.55430E-14 1.55040E-14 1.54640E-14 1.54230E-14 1.53820E-14 1.53400E-14 &
     1.52970E-14 1.52530E-14 1.52090E-14 1.51650E-14 1.51190E-14 1.50730E-14 1.50270E-14 1.49790E-14 &
     1.49320E-14 1.48830E-14 1.48340E-14 1.47850E-14 1.47350E-14 1.46840E-14 1.46330E-14 1.45810E-14 &
     1.45290E-14 1.44760E-14 1.44230E-14 1.43690E-14 1.43150E-14 1.42610E-14 1.42060E-14 1.41500E-14 &
     1.40940E-14 1.40380E-14 1.39810E-14 1.39240E-14 1.38670E-14 1.38090E-14 1.37510E-14 1.36930E-14 &
     1.36340E-14 1.35750E-14 1.35160E-14 1.34570E-14 1.33970E-14 1.33370E-14 1.32770E-14 1.32160E-14 &
     1.31560E-14 1.30950E-14 1.30340E-14 1.29730E-14 1.29120E-14 1.28500E-14 1.27890E-14 1.27270E-14 &
     1.26650E-14 1.26040E-14 1.25420E-14 1.24800E-14 1.24180E-14 1.23560E-14 1.22940E-14 1.22320E-14 &
     1.21700E-14 1.21080E-14 1.20460E-14 1.19840E-14 1.19220E-14 1.18610E-14 1.17990E-14 1.17370E-14 &
     1.16760E-14 1.16150E-14 1.15540E-14 1.14930E-14 1.14320E-14 1.13710E-14 1.13110E-14 1.12510E-14 &
     1.11910E-14 1.11310E-14 1.10720E-14 1.10120E-14 1.09530E-14 1.08950E-14 1.08360E-14 1.07780E-14 &
     1.07200E-14 1.06630E-14 1.06060E-14 1.05490E-14 1.04920E-14 1.04360E-14 1.03800E-14 1.03250E-14 &
     1.02700E-14 1.02150E-14 1.01610E-14 1.01080E-14 1.00540E-14 1.00010E-14 9.94880E-15 9.89680E-15 &
     9.84520E-15 9.79410E-15 9.74350E-15 9.69330E-15 9.64360E-15 9.59440E-15 9.54570E-15 9.49750E-15 &
     9.44990E-15 9.40270E-15 9.35600E-15 9.30990E-15 9.26430E-15 9.21920E-15 9.17470E-15 9.13070E-15 &
     9.08730E-15 9.04440E-15 9.00200E-15 8.96020E-15 8.91900E-15 8.87840E-15 8.83830E-15 8.79880E-15 &
     8.75990E-15 8.72150E-15 8.68380E-15 8.64660E-15 8.61000E-15 8.57390E-15 8.53850E-15 8.50370E-15 &
     8.46940E-15 8.43570E-15 8.40260E-15 8.37010E-15 8.33820E-15 8.30690E-15 8.27620E-15 8.24600E-15 &
     8.21650E-15 8.18750E-15 8.15910E-15 8.13130E-15 8.10410E-15 8.07740E-15 8.05130E-15 8.02580E-15 &
     8.00090E-15 7.97650E-15 7.95270E-15 7.92950E-15 7.90680E-15 7.88470E-15 7.86320E-15 7.84210E-15 &
     7.82170E-15 7.80180E-15 7.78240E-15 7.76350E-15 7.74520E-15 7.72740E-15 7.71020E-15 7.69340E-15 &
     7.67720E-15 7.66140E-15 7.64620E-15 7.63150E-15 7.61720E-15 7.60340E-15 7.59010E-15 7.57730E-15 &
     7.56490E-15 7.55300E-15 7.54150E-15 7.53050E-15 7.51990E-15 7.50980E-15 7.50010E-15 7.49080E-15 &
     7.48190E-15 7.47340E-15 7.46530E-15 7.45750E-15 7.45020E-15 7.44320E-15 7.43660E-15 7.43040E-15 &
     7.42450E-15 7.41890E-15 7.41370E-15 7.40880E-15 7.40420E-15 7.39990E-15 7.39590E-15 7.39230E-15 &
     7.38890E-15 7.38570E-15 7.38290E-15 7.38020E-15 7.37790E-15 7.37580E-15 7.37390E-15 7.37220E-15 &
     7.37080E-15 7.36950E-15 7.36850E-15 7.36760E-15 7.36690E-15 7.36640E-15 7.36610E-15 7.36590E-15 &
     7.36590E-15 7.36600E-15 7.36620E-15 7.36660E-15 7.36700E-15 7.36760E-15 7.36830E-15 7.36910E-15 &
     7.36990E-15 7.37090E-15 7.37180E-15 7.37290E-15 7.37400E-15 7.37510E-15 7.37630E-15 7.37750E-15 &
     7.37870E-15 7.37990E-15 7.38110E-15 7.38230E-15 7.38350E-15 7.38460E-15 7.38570E-15 7.38680E-15 &
     7.38790E-15 7.38880E-15 7.38970E-15 7.39060E-15 7.39140E-15 7.39210E-15 7.39270E-15 7.39320E-15 &
     7.39360E-15 7.39390E-15 7.39410E-15 7.39420E-15 7.39410E-15 7.39390E-15 7.39360E-15 7.39310E-15 &
     7.39240E-15 7.39160E-15 7.39070E-15 7.38950E-15 7.38820E-15 7.38670E-15 7.38500E-15 7.38310E-15 &
     7.38100E-15 7.37880E-15 7.37630E-15 7.37360E-15 7.37070E-15 7.36750E-15 7.36420E-15 7.36060E-15 &
     7.35680E-15 7.35270E-15 7.34840E-15 7.34390E-15 7.33910E-15 7.33410E-15 7.32880E-15 7.32320E-15 &
     7.31740E-15 7.31140E-15 7.30500E-15 7.29850E-15 7.29160E-15 7.28450E-15 7.27710E-15 7.26940E-15 &
     7.26140E-15 7.25320E-15 7.24470E-15 7.23590E-15 7.22680E-15 7.21750E-15 7.20780E-15 7.19790E-15 &
     7.18770E-15 7.17720E-15 7.16640E-15 7.15530E-15 7.14400E-15 7.13230E-15 7.12040E-15 7.10810E-15 &
     7.09560E-15 7.08280E-15 7.06970E-15 7.05630E-15 7.04270E-15 7.02870E-15 7.01450E-15 6.99990E-15 &
     6.98510E-15 6.97000E-15 6.95460E-15 6.93890E-15 6.92300E-15 6.90680E-15 6.89030E-15 6.87350E-15 &
     6.85650E-15 6.83920E-15 6.82160E-15 6.80380E-15 6.78570E-15 6.76730E-15 6.74870E-15 6.72980E-15 &
     6.71070E-15 6.69130E-15 6.67170E-15 6.65180E-15 6.63170E-15 6.61130E-15 6.59080E-15 6.57000E-15 &
     6.54890E-15 6.52770E-15 6.50630E-15 6.48460E-15 6.46280E-15 6.44070E-15 6.41850E-15 6.39600E-15 &
     6.37340E-15 6.35050E-15 6.32750E-15 6.30430E-15 6.28100E-15 6.25740E-15 6.23370E-15 6.20980E-15 &
     6.18580E-15 6.16160E-15 6.13730E-15 6.11280E-15 6.08820E-15 6.06350E-15 6.03860E-15 6.01360E-15 &
     5.98850E-15 5.96320E-15 5.93780E-15 5.91240E-15 5.88680E-15 5.86110E-15 5.83540E-15 5.80950E-15 &
     5.78360E-15 5.75760E-15 5.73150E-15 5.70540E-15 5.67920E-15 5.65290E-15 5.62660E-15 5.60020E-15 &
     5.57380E-15 5.54740E-15 5.52090E-15 5.49440E-15 5.46790E-15 5.44140E-15 5.41480E-15 5.38830E-15 &
     5.36180E-15 5.33520E-15 5.30870E-15 5.28220E-15 5.25570E-15 5.22920E-15 5.20270E-15 5.17630E-15 &
     5.14990E-15 5.12360E-15 5.09730E-15 5.07110E-15 5.04490E-15 5.01880E-15 4.99270E-15 4.96680E-15 &
     4.94090E-15 4.91500E-15 4.88930E-15 4.86370E-15 4.83810E-15 4.81260E-15 4.78730E-15 4.76200E-15 &
     4.73690E-15 4.71180E-15 4.68690E-15 4.66210E-15 4.63750E-15 4.61290E-15 4.58850E-15 4.56420E-15 &
     4.54000E-15 4.51600E-15 4.49220E-15 4.46850E-15 4.44490E-15 4.42150E-15 4.39830E-15 4.37520E-15 &
     4.35230E-15 4.32950E-15 4.30700E-15 4.28460E-15 4.26230E-15 4.24030E-15 4.21840E-15 4.19670E-15 &
     4.17530E-15 4.15390E-15 4.13280E-15 4.11190E-15 4.09120E-15 4.07060E-15 4.05030E-15 4.03020E-15 &
     4.01030E-15 3.99050E-15 3.97100E-15 3.95170E-15 3.93260E-15 3.91370E-15 3.89510E-15 3.87660E-15 &
     3.85840E-15 3.84030E-15 3.82250E-15 3.80490E-15 3.78760E-15 3.77040E-15 3.75350E-15 3.73680E-15 &
     3.72030E-15 3.70400E-15 3.68800E-15 3.67220E-15 3.65660E-15 3.64120E-15 3.62610E-15 3.61120E-15 &
     3.59650E-15 3.58200E-15 3.56780E-15 3.55380E-15 3.54000E-15 3.52640E-15 3.51310E-15 3.49990E-15 &
     3.48700E-15 3.47430E-15 3.46190E-15 3.44960E-15 3.43760E-15 3.42580E-15 3.41420E-15 3.40280E-15 &
     3.39170E-15 3.38070E-15 3.37000E-15 3.35940E-15 3.34910E-15 3.33900E-15 3.32910E-15 3.31940E-15 &
     3.30990E-15 3.30060E-15 3.29150E-15 3.28260E-15 3.27390E-15 3.26540E-15 3.25710E-15 3.24890E-15 &
     3.24100E-15 3.23320E-15 3.22570E-15 3.21830E-15 3.21110E-15 3.20410E-15 3.19720E-15 3.19050E-15 &
     3.18400E-15 3.17760E-15 3.17140E-15 3.16540E-15 3.15950E-15 3.15380E-15 3.14830E-15 3.14290E-15 &
     3.13760E-15 3.13250E-15 3.12760E-15 3.12280E-15 3.11810E-15 3.11360E-15 3.10920E-15 3.10490E-15 &
     3.10070E-15 3.09670E-15 3.09280E-15 3.08910E-15 3.08540E-15 3.08190E-15 3.07840E-15 3.07510E-15 &
     3.07190E-15 3.06880E-15 3.06570E-15 3.06280E-15 3.06000E-15 3.05720E-15 3.05460E-15 3.05200E-15 &
     3.04950E-15 3.04710E-15 3.04470E-15 3.04250E-15 3.04030E-15 3.03810E-15 3.03600E-15 3.03400E-15 &
     3.03210E-15 3.03020E-15 3.02830E-15 3.02650E-15 3.02470E-15 3.02300E-15 3.02130E-15 3.01970E-15 &
     3.01810E-15 3.01650E-15 3.01490E-15 3.01340E-15 3.01190E-15 3.01040E-15 3.00890E-15 3.00750E-15 &
     3.00600E-15 3.00460E-15 3.00310E-15 3.00170E-15 3.00030E-15 2.99880E-15 2.99740E-15 2.99590E-15 &
     2.99450E-15 2.99300E-15 2.99150E-15 2.99000E-15 2.98850E-15 2.98700E-15 2.98540E-15 2.98380E-15 &
     2.98220E-15 2.98050E-15 2.97880E-15 2.97710E-15 2.97530E-15 2.97350E-15 2.97170E-15 2.96980E-15 &
     2.96790E-15 2.96590E-15 2.96390E-15 2.96180E-15 2.95960E-15 2.95740E-15 2.95520E-15 2.95290E-15 &
     2.95050E-15 2.94810E-15 2.94560E-15 2.94310E-15 2.94050E-15 2.93780E-15 2.93510E-15 2.93220E-15 &
     2.92940E-15 2.92640E-15 2.92340E-15 2.92030E-15 2.91710E-15 2.91380E-15 2.91050E-15 2.90710E-15 &
     2.90360E-15 2.90010E-15 2.89640E-15 2.89270E-15 2.88890E-15 2.88500E-15 2.88100E-15 2.87690E-15 &
     2.87280E-15 2.86860E-15 2.86430E-15 2.85990E-15 2.85540E-15 2.85080E-15 2.84610E-15 2.84140E-15 &
     2.83650E-15 2.83160E-15 2.82650E-15 2.82140E-15 2.81620E-15 2.81090E-15 2.80560E-15 2.80010E-15 &
     2.79460E-15 2.78890E-15 2.78320E-15 2.77740E-15 2.77160E-15 2.76560E-15 2.75950E-15 2.75340E-15 &
     2.74720E-15 2.74090E-15 2.73450E-15 2.72800E-15 2.72150E-15 2.71480E-15 2.70810E-15 2.70130E-15 &
     2.69450E-15 2.68750E-15 2.68050E-15 2.67340E-15 2.66620E-15 2.65900E-15 2.65160E-15 2.64420E-15 &
     2.63680E-15 2.62920E-15 2.62160E-15 2.61400E-15 2.60620E-15 2.59840E-15 2.59050E-15 2.58260E-15 &
     2.57460E-15 2.56650E-15 2.55840E-15 2.55020E-15 2.54200E-15 2.53370E-15 2.52530E-15 2.51690E-15 &
     2.50840E-15 2.49990E-15 2.49140E-15 2.48280E-15 2.47410E-15 2.46540E-15 2.45670E-15 2.44790E-15 &
     2.43900E-15 2.43020E-15 2.42130E-15 2.41230E-15 2.40340E-15 2.39430E-15 2.38530E-15 2.37620E-15 &
     2.36710E-15 2.35800E-15 2.34890E-15 2.33970E-15 2.33050E-15 2.32130E-15 2.31210E-15 2.30280E-15 &
     2.29360E-15 2.28430E-15 2.27500E-15 2.26570E-15 2.25640E-15 2.24710E-15 2.23780E-15 2.22840E-15 &
     2.21910E-15 2.20980E-15 2.20050E-15 2.19110E-15 2.18180E-15 2.17250E-15 2.16320E-15 2.15390E-15 &
     2.14460E-15 2.13540E-15 2.12610E-15 2.11690E-15 2.10760E-15 2.09840E-15 2.08920E-15 2.08000E-15 &
     2.07090E-15 2.06180E-15 2.05270E-15 2.04360E-15 2.03450E-15 2.02550E-15 2.01650E-15 2.00760E-15 &
     1.99870E-15 1.98980E-15 1.98090E-15 1.97210E-15 1.96340E-15 1.95460E-15 1.94590E-15 1.93730E-15 &
     1.92870E-15 1.92010E-15 1.91160E-15 1.90310E-15 1.89470E-15 1.88640E-15 1.87800E-15 1.86980E-15 &
     1.86160E-15 1.85340E-15 1.84530E-15 1.83720E-15 1.82920E-15 1.82130E-15 1.81340E-15 1.80560E-15 &
     1.79790E-15 1.79020E-15 1.78250E-15 1.77490E-15 1.76740E-15 1.76000E-15 1.75260E-15 1.74530E-15 &
     1.73800E-15 1.73080E-15 1.72370E-15 1.71670E-15 1.70970E-15 1.70280E-15 1.69590E-15 1.68920E-15 &
     1.68250E-15 1.67580E-15 1.66930E-15 1.66280E-15 1.65630E-15 1.65000E-15 1.64370E-15 1.63750E-15 &
     1.63140E-15 1.62530E-15 1.61930E-15 1.61340E-15 1.60760E-15 1.60180E-15 1.59610E-15 1.59050E-15 &
     1.58500E-15 1.57950E-15 1.57410E-15 1.56880E-15 1.56350E-15 1.55830E-15 1.55320E-15 1.54820E-15 &
     1.54320E-15 1.53830E-15 1.53350E-15 1.52880E-15 1.52410E-15 1.51950E-15 1.51500E-15 1.51050E-15 &
     1.50610E-15 1.50180E-15 1.49760E-15 1.49340E-15 1.48930E-15 1.48520E-15 1.48130E-15 1.47730E-15 &
     1.47350E-15 1.46970E-15 1.46600E-15 1.46240E-15 1.45880E-15 1.45530E-15 1.45180E-15 1.44840E-15 &
     1.44510E-15 1.44180E-15 1.43860E-15 1.43550E-15 1.43240E-15 1.42940E-15 1.42640E-15 1.42350E-15 &
     1.42060E-15 1.41790E-15 1.41510E-15 1.41240E-15 1.40980E-15 1.40720E-15 1.40470E-15 1.40220E-15 &
     1.39980E-15 1.39740E-15 1.39500E-15 1.39280E-15 1.39050E-15 1.38830E-15 1.38620E-15 1.38400E-15 &
     1.38200E-15 1.38000E-15 1.37800E-15 1.37600E-15 1.37410E-15 1.37220E-15 1.37040E-15 1.36860E-15 &
     1.36680E-15 1.36510E-15 1.36340E-15 1.36170E-15 1.36010E-15 1.35850E-15 1.35690E-15 1.35540E-15 &
     1.35380E-15 1.35230E-15 1.35090E-15 1.34940E-15 1.34800E-15 1.34660E-15 1.34520E-15 1.34380E-15 &
     1.34250E-15 1.34110E-15 1.33980E-15 1.33850E-15 1.33720E-15 1.33600E-15 1.33470E-15 1.33350E-15 &
     1.33220E-15 1.33100E-15 1.32980E-15 1.32860E-15 1.32730E-15 1.32610E-15 1.32490E-15 1.32370E-15 &
     1.32260E-15 1.32140E-15 1.32020E-15 1.31900E-15 1.31780E-15 1.31660E-15 1.31540E-15 1.31420E-15 &
     1.31300E-15 1.31180E-15 1.31060E-15 1.30930E-15 1.30810E-15 1.30690E-15 1.30560E-15 1.30440E-15 &
     1.30310E-15 1.30180E-15 1.30050E-15 1.29920E-15 1.29790E-15 1.29660E-15 1.29530E-15 1.29390E-15 &
     1.29250E-15 1.29110E-15 1.28970E-15 1.28830E-15 1.28680E-15 1.28540E-15 1.28390E-15 1.28240E-15 &
     1.28080E-15 1.27930E-15 1.27770E-15 1.27610E-15 1.27450E-15 1.27290E-15 1.27120E-15 1.26950E-15 &
     1.26780E-15 1.26600E-15 1.26430E-15 1.26240E-15 1.26060E-15 1.25880E-15 1.25690E-15 1.25490E-15 &
     1.25300E-15 1.25100E-15 1.24900E-15 1.24700E-15 1.24490E-15 1.24280E-15 1.24070E-15 1.23860E-15 &
     1.23640E-15 1.23420E-15 1.23200E-15 1.22970E-15 1.22740E-15 1.22510E-15 1.22270E-15 1.22030E-15 &
     1.21790E-15 1.21550E-15 1.21300E-15 1.21050E-15 1.20790E-15 1.20540E-15 1.20280E-15 1.20010E-15 &
     1.19750E-15 1.19480E-15 1.19200E-15 1.18930E-15 1.18650E-15 1.18370E-15 1.18080E-15 1.17790E-15 &
     1.17500E-15 1.17210E-15 1.16910E-15 1.16610E-15 1.16310E-15 1.16000E-15 1.15690E-15 1.15380E-15 &
     1.15060E-15 1.14750E-15 1.14430E-15 1.14100E-15 1.13780E-15 1.13450E-15 1.13120E-15 1.12780E-15 &
     1.12450E-15 1.12110E-15 1.11770E-15 1.11420E-15 1.11080E-15 1.10730E-15 1.10380E-15 1.10020E-15 &
     1.09670E-15 1.09310E-15 1.08950E-15 1.08590E-15 1.08220E-15 1.07860E-15 1.07490E-15 1.07120E-15 &
     1.06750E-15 1.06370E-15 1.06000E-15 1.05620E-15 1.05240E-15 1.04860E-15 1.04470E-15 1.04090E-15 &
     1.03700E-15 1.03320E-15 1.02930E-15 1.02540E-15 1.02150E-15 1.01750E-15 1.01360E-15 1.00970E-15 &
     1.00570E-15 1.00170E-15 9.97770E-16 9.93780E-16 9.89790E-16 9.85790E-16 9.81780E-16 9.77770E-16 &
     9.73750E-16 9.69730E-16 9.65700E-16 9.61670E-16 9.57640E-16 9.53600E-16 9.49560E-16 9.45510E-16 &
     9.41470E-16 9.37420E-16 9.33380E-16 9.29330E-16 9.25290E-16 9.21250E-16 9.17200E-16 9.13170E-16 &
     9.09130E-16 9.05100E-16 9.01070E-16 8.97050E-16 8.93040E-16 8.89020E-16 8.85020E-16 8.81020E-16 &
     8.77030E-16 8.73050E-16 8.69080E-16 8.65110E-16 8.61160E-16 8.57210E-16 8.53280E-16 8.49360E-16 &
     8.45450E-16 8.41550E-16 8.37670E-16 8.33800E-16 8.29950E-16 8.26110E-16 8.22290E-16 8.18480E-16 &
     8.14690E-16 8.10910E-16 8.07150E-16 8.03410E-16 7.99680E-16 7.95980E-16 7.92290E-16 7.88620E-16 &
     7.84970E-16 7.81340E-16 7.77740E-16 7.74150E-16 7.70580E-16 7.67030E-16 7.63500E-16 7.60000E-16 &
     7.56520E-16 7.53050E-16 7.49610E-16 7.46200E-16 7.42810E-16 7.39440E-16 7.36090E-16 7.32770E-16 &
     7.29480E-16 7.26210E-16 7.22960E-16 7.19750E-16 7.16550E-16 7.13380E-16 7.10240E-16 7.07130E-16 &
     7.04040E-16 7.00980E-16 6.97940E-16 6.94930E-16 6.91950E-16 6.89000E-16 6.86070E-16 6.83170E-16 &
     6.80300E-16 6.77460E-16 6.74640E-16 6.71850E-16 6.69090E-16 6.66360E-16 6.63650E-16 6.60980E-16 &
     6.58330E-16 6.55710E-16 6.53120E-16 6.50560E-16 6.48020E-16 6.45520E-16 6.43040E-16 6.40590E-16 &
     6.38170E-16 6.35780E-16 6.33420E-16 6.31080E-16 6.28770E-16 6.26490E-16 6.24240E-16 6.22020E-16 &
     6.19830E-16 6.17660E-16 6.15520E-16 6.13420E-16 6.11330E-16 6.09280E-16 6.07260E-16 6.05260E-16 &
     6.03290E-16 6.01350E-16 5.99440E-16 5.97550E-16 5.95690E-16 5.93860E-16 5.92050E-16 5.90270E-16 &
     5.88510E-16 5.86790E-16 5.85080E-16 5.83400E-16 5.81750E-16 5.80120E-16 5.78510E-16 5.76930E-16 &
     5.75370E-16 5.73830E-16 5.72320E-16 5.70830E-16 5.69370E-16 5.67920E-16 5.66500E-16 5.65100E-16 &
     5.63710E-16 5.62360E-16 5.61020E-16 5.59700E-16 5.58400E-16 5.57120E-16 5.55860E-16 5.54630E-16 &
     5.53410E-16 5.52200E-16 5.51020E-16 5.49860E-16 5.48710E-16 5.47580E-16 5.46470E-16 5.45380E-16 &
     5.44300E-16 5.43240E-16 5.42200E-16 5.41170E-16 5.40160E-16 5.39160E-16 5.38180E-16 5.37210E-16 &
     5.36260E-16 5.35320E-16 5.34400E-16 5.33490E-16 5.32590E-16 5.31710E-16 5.30840E-16 5.29980E-16 &
     5.29130E-16 5.28300E-16 5.27480E-16 5.26660E-16 5.25860E-16 5.25070E-16 5.24290E-16 5.23510E-16 &
     5.22750E-16 5.22000E-16 5.21250E-16 5.20510E-16 5.19780E-16 5.19060E-16 5.18350E-16 5.17640E-16 &
     5.16940E-16 5.16240E-16 5.15550E-16 5.14870E-16 5.14190E-16 5.13510E-16 5.12840E-16 5.12170E-16 &
     5.11510E-16 5.10850E-16 5.10200E-16 5.09540E-16 5.08890E-16 5.08240E-16 5.07590E-16 5.06950E-16 &
     5.06300E-16 5.05660E-16 5.05010E-16 5.04370E-16 5.03720E-16 5.03080E-16 5.02430E-16 5.01780E-16 &
     5.01130E-16 5.00480E-16 4.99820E-16 4.99160E-16 4.98500E-16 4.97840E-16 4.97170E-16 4.96500E-16 &
     4.95820E-16 4.95140E-16 4.94450E-16 4.93760E-16 4.93060E-16 4.92350E-16 4.91640E-16 4.90920E-16 &
     4.90190E-16 4.89460E-16 4.88720E-16 4.87960E-16 4.87210E-16 4.86440E-16 4.85660E-16 4.84870E-16 &
     4.84070E-16 4.83270E-16 4.82450E-16 4.81620E-16 4.80770E-16 4.79920E-16 4.79060E-16 4.78180E-16 &
     4.77300E-16 4.76400E-16 4.75490E-16 4.74580E-16 4.73650E-16 4.72710E-16 4.71760E-16 4.70800E-16 &
     4.69840E-16 4.68860E-16 4.67880E-16 4.66880E-16 4.65880E-16 4.64860E-16 4.63840E-16 4.62810E-16 &
     4.61770E-16 4.60730E-16 4.59670E-16 4.58610E-16 4.57540E-16 4.56460E-16 4.55370E-16 4.54280E-16 &
     4.53180E-16 4.52080E-16 4.50960E-16 4.49840E-16 4.48710E-16 4.47580E-16 4.46440E-16 4.45300E-16 &
     4.44150E-16 4.42990E-16 4.41830E-16 4.40660E-16 4.39480E-16 4.38310E-16 4.37120E-16 4.35930E-16 &
     4.34740E-16 4.33540E-16 4.32340E-16 4.31140E-16 4.29930E-16 4.28710E-16 4.27500E-16 4.26270E-16 &
     4.25050E-16 4.23820E-16 4.22590E-16 4.21360E-16 4.20120E-16 4.18880E-16 4.17640E-16 4.16400E-16 &
     4.15150E-16 4.13900E-16 4.12650E-16 4.11400E-16 4.10150E-16 4.08900E-16 4.07640E-16 4.06380E-16 &
     4.05130E-16 4.03870E-16 4.02610E-16 4.01350E-16 4.00090E-16 3.98830E-16 3.97570E-16 3.96310E-16 &
     3.95050E-16 3.93790E-16 3.92530E-16 3.91280E-16 3.90020E-16 3.88760E-16 3.87510E-16 3.86260E-16 &
     3.85010E-16 3.83760E-16 3.82510E-16 3.81270E-16 3.80020E-16 3.78780E-16 3.77540E-16 3.76310E-16 &
     3.75070E-16 3.73840E-16 3.72620E-16 3.71390E-16 3.70170E-16 3.68960E-16 3.67750E-16 3.66540E-16 &
     3.65330E-16 3.64130E-16 3.62940E-16 3.61740E-16 3.60560E-16 3.59380E-16 3.58200E-16 3.57030E-16 &
     3.55860E-16 3.54700E-16 3.53540E-16 3.52390E-16 3.51250E-16 3.50110E-16 3.48980E-16 3.47850E-16 &
     3.46730E-16 3.45620E-16 3.44510E-16 3.43410E-16 3.42320E-16 3.41220E-16 3.40140E-16 3.39060E-16 &
     3.37990E-16 3.36920E-16 3.35860E-16 3.34800E-16 3.33750E-16 3.32710E-16 3.31670E-16 3.30630E-16 &
     3.29610E-16 3.28580E-16 3.27570E-16 3.26550E-16 3.25550E-16 3.24540E-16 3.23550E-16 3.22560E-16 &
     3.21570E-16 3.20590E-16 3.19610E-16 3.18640E-16 3.17680E-16 3.16720E-16 3.15760E-16 3.14810E-16 &
     3.13870E-16 3.12920E-16 3.11990E-16 3.11060E-16 3.10130E-16 3.09210E-16 3.08290E-16 3.07380E-16 &
     3.06470E-16 3.05570E-16 3.04670E-16 3.03780E-16 3.02890E-16 3.02000E-16 3.01120E-16 3.00250E-16 &
     2.99380E-16 2.98510E-16 2.97650E-16 2.96790E-16 2.95930E-16 2.95080E-16 2.94240E-16 2.93400E-16 &
     2.92560E-16 2.91720E-16 2.90890E-16 2.90070E-16 2.89250E-16 2.88430E-16 2.87620E-16 2.86810E-16 &
     2.86000E-16 2.85200E-16 2.84400E-16 2.83610E-16 2.82820E-16 2.82030E-16 2.81250E-16 2.80470E-16 &
     2.79690E-16 2.78920E-16 2.78150E-16 2.77380E-16 2.76620E-16 2.75860E-16 2.75110E-16 2.74360E-16 &
     2.73610E-16 2.72860E-16 2.72120E-16 2.71380E-16 2.70650E-16 2.69910E-16 2.69180E-16 2.68460E-16 &
     2.67730E-16 2.67010E-16 2.66300E-16 2.65580E-16 2.64870E-16 2.64160E-16 2.63460E-16 2.62760E-16 &
     2.62060E-16 2.61360E-16 2.60660E-16 2.59970E-16 2.59280E-16 2.58600E-16 2.57910E-16 2.57230E-16 &
     2.56550E-16 2.55880E-16 2.55200E-16 2.54530E-16 2.53860E-16 2.53200E-16 2.52530E-16 2.51870E-16 &
     2.51210E-16 2.50560E-16 2.49900E-16 2.49250E-16 2.48600E-16 2.47950E-16 2.47300E-16 2.46660E-16 &
     2.46010E-16 2.45370E-16 2.44730E-16 2.44100E-16 2.43460E-16 2.42830E-16 2.42200E-16 2.41570E-16 &
     2.40940E-16 2.40310E-16 2.39690E-16 2.39070E-16 2.38450E-16 2.37830E-16 2.37210E-16 2.36590E-16 &
     2.35980E-16 2.35360E-16 2.34750E-16 2.34140E-16 2.33530E-16 2.32920E-16 2.32310E-16 2.31710E-16 &
     2.31100E-16 2.30500E-16 2.29900E-16 2.29300E-16 2.28700E-16 2.28100E-16 2.27500E-16 2.26900E-16 &
     2.26310E-16 2.25710E-16 2.25120E-16 2.24520E-16 2.23930E-16 2.23340E-16 2.22750E-16 2.22160E-16 &
     2.21570E-16 2.20980E-16 2.20390E-16 2.19800E-16 2.19220E-16 2.18630E-16 2.18040E-16 2.17460E-16 &
     2.16870E-16 2.16280E-16 2.15700E-16 2.15120E-16 2.14530E-16 2.13950E-16 2.13370E-16 2.12780E-16 &
     2.12200E-16 2.11620E-16 2.11040E-16 2.10460E-16 2.09880E-16 2.09300E-16 2.08720E-16 2.08140E-16 &
     2.07560E-16 2.06990E-16 2.06410E-16 2.05830E-16 2.05260E-16 2.04680E-16 2.04110E-16 2.03530E-16 &
     2.02960E-16 2.02390E-16 2.01810E-16 2.01240E-16 2.00670E-16 2.00100E-16 1.99530E-16 1.98960E-16 &
     1.98390E-16 1.97820E-16 1.97260E-16 1.96690E-16 1.96120E-16 1.95550E-16 1.94990E-16 1.94420E-16 &
     1.93860E-16 1.93300E-16 1.92730E-16 1.92170E-16 1.91610E-16 1.91050E-16 1.90490E-16 1.89930E-16 &
     1.89370E-16 1.88810E-16 1.88250E-16 1.87690E-16 1.87140E-16 1.86580E-16 1.86030E-16 1.85470E-16 &
     1.84920E-16 1.84360E-16 1.83810E-16 1.83260E-16 1.82710E-16 1.82160E-16 1.81610E-16 1.81060E-16 &
     1.80510E-16 1.79960E-16 1.79420E-16 1.78870E-16 1.78330E-16 1.77780E-16 1.77240E-16 1.76690E-16 &
     1.76150E-16 1.75610E-16 1.75070E-16 1.74530E-16 1.73990E-16 1.73450E-16 1.72910E-16 1.72370E-16 &
     1.71840E-16 1.71300E-16 1.70770E-16 1.70230E-16 1.69700E-16 1.69170E-16 1.68640E-16 1.68110E-16 &
     1.67580E-16 1.67050E-16 1.66520E-16 1.65990E-16 1.65470E-16 1.64940E-16 1.64410E-16 1.63890E-16 &
     1.63370E-16 1.62840E-16 1.62320E-16 1.61800E-16 1.61280E-16 1.60760E-16 1.60250E-16 1.59730E-16 &
     1.59210E-16 1.58700E-16 1.58180E-16 1.57670E-16 1.57160E-16 1.56640E-16 1.56130E-16 1.55620E-16 &
     1.55110E-16 1.54610E-16 1.54100E-16 1.53590E-16 1.53090E-16 1.52580E-16 1.52080E-16 1.51580E-16 &
     1.51080E-16 1.50570E-16 1.50080E-16 1.49580E-16 1.49080E-16 1.48580E-16 1.48090E-16 1.47590E-16 &
     1.47100E-16 1.46600E-16 1.46110E-16 1.45620E-16 1.45130E-16 1.44640E-16 1.44160E-16 1.43670E-16 &
     1.43180E-16 1.42700E-16 1.42210E-16 1.41730E-16 1.41250E-16 1.40770E-16 1.40290E-16 1.39810E-16 &
     1.39330E-16 1.38860E-16 1.38380E-16 1.37910E-16 1.37430E-16 1.36960E-16 1.36490E-16 1.36020E-16 &
     1.35550E-16 1.35080E-16 1.34620E-16 1.34150E-16 1.33690E-16 1.33220E-16 1.32760E-16 1.32300E-16 &
     1.31840E-16 1.31380E-16 1.30920E-16 1.30470E-16 1.30010E-16 1.29560E-16 1.29100E-16 1.28650E-16 &
     1.28200E-16 1.27750E-16 1.27300E-16 1.26850E-16 1.26410E-16 1.25960E-16 1.25520E-16 1.25080E-16 &
     1.24630E-16 1.24190E-16 1.23760E-16 1.23320E-16 1.22880E-16 1.22440E-16 1.22010E-16 1.21580E-16 &
     1.21150E-16 1.20710E-16 1.20290E-16 1.19860E-16 1.19430E-16 1.19000E-16 1.18580E-16 1.18160E-16 &
     1.17730E-16 1.17310E-16 1.16890E-16 1.16480E-16 1.16060E-16 1.15640E-16 1.15230E-16 1.14820E-16 &
     1.14400E-16 1.13990E-16 1.13580E-16 1.13180E-16 1.12770E-16 1.12360E-16 1.11960E-16 1.11560E-16 &
     1.11160E-16 1.10760E-16 1.10360E-16 1.09960E-16 1.09560E-16 1.09170E-16 1.08780E-16 1.08380E-16 &
     1.07990E-16 1.07600E-16 1.07220E-16 1.06830E-16 1.06440E-16 1.06060E-16 1.05680E-16 1.05300E-16 &
     1.04920E-16 1.04540E-16 1.04160E-16 1.03790E-16 1.03410E-16 1.03040E-16 1.02670E-16 1.02300E-16 &
     1.01930E-16 1.01560E-16 1.01200E-16 1.00830E-16 1.00470E-16 1.00110E-16 9.97500E-17 9.93920E-17 &
     9.90350E-17 9.86800E-17 9.83260E-17 9.79740E-17 9.76230E-17 9.72730E-17 9.69250E-17 9.65790E-17 &
     9.62340E-17 9.58900E-17 9.55480E-17 9.52070E-17 9.48680E-17 9.45310E-17 9.41950E-17 9.38600E-17 &
     9.35270E-17 9.31960E-17 9.28660E-17 9.25370E-17 9.22100E-17 9.18850E-17 9.15610E-17 9.12380E-17 &
     9.09170E-17 9.05980E-17 9.02800E-17 8.99630E-17 8.96480E-17 8.93340E-17 8.90210E-17 8.87110E-17 &
     8.84010E-17 8.80930E-17 8.77860E-17 8.74810E-17 8.71780E-17 8.68750E-17 8.65740E-17 8.62750E-17 &
     8.59770E-17 8.56800E-17 8.53850E-17 8.50910E-17 8.47980E-17 8.45070E-17 8.42170E-17 8.39290E-17 &
     8.36410E-17 8.33560E-17 8.30710E-17 8.27880E-17 8.25070E-17 8.22260E-17 8.19470E-17 8.16700E-17 &
     8.13930E-17 8.11180E-17 8.08450E-17 8.05720E-17 8.03010E-17 8.00320E-17 7.97630E-17 7.94960E-17 &
     7.92300E-17 7.89660E-17 7.87020E-17 7.84410E-17 7.81800E-17 7.79200E-17 7.76620E-17 7.74050E-17 &
     7.71500E-17 7.68950E-17 7.66420E-17 7.63900E-17 7.61400E-17 7.58900E-17 7.56420E-17 7.53950E-17 &
     7.51500E-17 7.49050E-17 7.46620E-17 7.44200E-17 7.41790E-17 7.39390E-17 7.37010E-17 7.34630E-17 &
     7.32270E-17 7.29930E-17 7.27590E-17 7.25260E-17 7.22950E-17 7.20650E-17 7.18360E-17 7.16080E-17 &
     7.13810E-17 7.11560E-17 7.09310E-17 7.07080E-17 7.04860E-17 7.02650E-17 7.00450E-17 6.98260E-17 &
     6.96080E-17 6.93920E-17 6.91760E-17 6.89620E-17 6.87490E-17 6.85370E-17 6.83260E-17 6.81160E-17 &
     6.79070E-17 6.76990E-17 6.74920E-17 6.72870E-17 6.70820E-17 6.68790E-17 6.66760E-17 6.64750E-17 &
     6.62740E-17 6.60750E-17 6.58770E-17 6.56800E-17 6.54830E-17 6.52880E-17 6.50940E-17 6.49010E-17 &
     6.47090E-17 6.45180E-17 6.43270E-17 6.41380E-17 6.39500E-17 6.37630E-17 6.35770E-17 6.33920E-17 &
     6.32080E-17 6.30240E-17 6.28420E-17 6.26610E-17 6.24800E-17 6.23010E-17 6.21230E-17 6.19450E-17 &
     6.17690E-17 6.15930E-17 6.14180E-17 6.12450E-17 6.10720E-17 6.09000E-17 6.07290E-17 6.05590E-17 &
     6.03900E-17 6.02210E-17 6.00540E-17 5.98880E-17 5.97220E-17 5.95570E-17 5.93930E-17 5.92300E-17 &
     5.90680E-17 5.89070E-17 5.87470E-17 5.85870E-17 5.84290E-17 5.82710E-17 5.81140E-17 5.79580E-17 &
     5.78030E-17 5.76480E-17 5.74940E-17 5.73420E-17 5.71900E-17 5.70390E-17 5.68880E-17 5.67390E-17 &
     5.65900E-17 5.64420E-17 5.62950E-17 5.61480E-17 5.60030E-17 5.58580E-17 5.57140E-17 5.55710E-17 &
     5.54280E-17 5.52860E-17 5.51450E-17 5.50050E-17 5.48660E-17 5.47270E-17 5.45890E-17 5.44520E-17 &
     5.43150E-17 5.41790E-17 5.40440E-17 5.39100E-17 5.37760E-17 5.36430E-17 5.35110E-17 5.33790E-17 &
     5.32480E-17 5.31180E-17 5.29890E-17 5.28600E-17 5.27320E-17 5.26040E-17 5.24780E-17 5.23510E-17 &
     5.22260E-17 5.21010E-17 5.19770E-17 5.18530E-17 5.17310E-17 5.16080E-17 5.14870E-17 5.13660E-17 &
     5.12450E-17 5.11260E-17 5.10060E-17 5.08880E-17 5.07700E-17 5.06530E-17 5.05360E-17 5.04200E-17 &
     5.03040E-17 5.01890E-17 5.00750E-17 4.99610E-17 4.98480E-17 4.97350E-17 4.96230E-17 4.95120E-17 &
     4.94010E-17 4.92900E-17 4.91800E-17 4.90710E-17 4.89620E-17 4.88540E-17 4.87460E-17 4.86390E-17 &
     4.85320E-17 4.84260E-17 4.83200E-17 4.82150E-17 4.81100E-17 4.80060E-17 4.79020E-17 4.77990E-17 &
     4.76960E-17 4.75940E-17 4.74920E-17 4.73910E-17 4.72900E-17 4.71890E-17 4.70890E-17 4.69900E-17 &
     4.68910E-17 4.67920E-17 4.66940E-17 4.65960E-17 4.64990E-17 4.64020E-17 4.63060E-17 4.62090E-17 &
     4.61140E-17 4.60190E-17 4.59240E-17 4.58290E-17 4.57350E-17 4.56420E-17 4.55480E-17 4.54550E-17 &
     4.53630E-17 4.52710E-17 4.51790E-17 4.50880E-17 4.49960E-17 4.49060E-17 4.48150E-17 4.47250E-17 &
     4.46360E-17 4.45460E-17 4.44570E-17 4.43690E-17 4.42800E-17 4.41920E-17 4.41050E-17 4.40170E-17 &
     4.39300E-17 4.38430E-17 4.37570E-17 4.36700E-17 4.35840E-17 4.34990E-17 4.34130E-17 4.33280E-17 &
     4.32440E-17 4.31590E-17 4.30750E-17 4.29910E-17 4.29070E-17 4.28230E-17 4.27400E-17 4.26570E-17 &
     4.25740E-17 4.24910E-17 4.24090E-17 4.23270E-17 4.22450E-17 4.21630E-17 4.20820E-17 4.20000E-17 &
     4.19190E-17 4.18380E-17 4.17580E-17 4.16770E-17 4.15970E-17 4.15170E-17 4.14360E-17 4.13570E-17 &
     4.12770E-17 4.11970E-17 4.11180E-17 4.10390E-17 4.09600E-17 4.08810E-17 4.08020E-17 4.07230E-17 &
     4.06450E-17 4.05660E-17 4.04880E-17 4.04100E-17 4.03320E-17 4.02540E-17 4.01760E-17 4.00980E-17 &
     4.00210E-17 3.99430E-17 3.98660E-17 3.97880E-17 3.97110E-17 3.96330E-17 3.95560E-17 3.94790E-17 &
     3.94020E-17 3.93250E-17 3.92480E-17 3.91710E-17 3.90940E-17 3.90170E-17 3.89400E-17 3.88640E-17 &
     3.87870E-17 3.87100E-17 3.86330E-17 3.85570E-17 3.84800E-17 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 &
     0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 4.59740E-13 5.77720E-13 5.93110E-13 5.48060E-13 &
     4.75080E-13 3.94620E-13 3.19100E-13 2.54700E-13 2.04740E-13 1.70300E-13 1.50920E-13 1.43540E-13 &
     1.43970E-13 1.49280E-13 1.57710E-13 1.67280E-13 1.76300E-13 1.84090E-13 1.89920E-13 1.94780E-13 &
     1.99410E-13 2.04130E-13 2.09250E-13 2.14350E-13 2.20180E-13 2.26640E-13 2.32710E-13 2.38780E-13 &
     2.43530E-13 2.47860E-13 2.50570E-13 2.51320E-13 2.50610E-13 2.47850E-13 2.44610E-13 2.39960E-13 &
     2.33780E-13 2.26180E-13 2.17170E-13 2.06840E-13 1.95570E-13 1.83630E-13 1.70490E-13 1.55530E-13 &
     1.38450E-13 1.19090E-13 9.76260E-14 7.43140E-14 4.94640E-14 2.33240E-14 -4.12010E-15 -3.29560E-14 &
     -6.32960E-14 -9.52900E-14 -1.29220E-13 -1.65320E-13 -2.03490E-13 -2.43570E-13 -2.85470E-13 -3.29150E-13 &
     -3.74570E-13 -4.21700E-13 -4.70520E-13 -5.21000E-13 -5.73100E-13 -6.26790E-13 -6.82040E-13 -7.38830E-13 &
     -7.97180E-13 -8.57140E-13 -9.18750E-13 -9.82010E-13 -1.04680E-12 -1.11310E-12 -1.18070E-12 -1.24940E-12 &
     -1.31950E-12 -1.39080E-12 -1.46340E-12 -1.53740E-12 -1.61270E-12 -1.68930E-12 -1.76690E-12 -1.84550E-12 &
     -1.92510E-12 -2.00550E-12 -2.08660E-12 -2.16850E-12 -2.25100E-12 -2.33420E-12 -2.41800E-12 -2.50240E-12 &
     -2.58740E-12 -2.67280E-12 -2.75880E-12 -2.84500E-12 -2.93160E-12 -3.01850E-12 -3.10540E-12 -3.19250E-12 &
     -3.27950E-12 -3.36650E-12 -3.45330E-12 -3.53990E-12 -3.62620E-12 -3.71210E-12 -3.79760E-12 -3.88250E-12 &
     -3.96680E-12 -4.05040E-12 -4.13330E-12 -4.21540E-12 -4.29660E-12 -4.37680E-12 -4.45600E-12 -4.53410E-12 &
     -4.61110E-12 -4.68680E-12 -4.76120E-12 -4.83420E-12 -4.90570E-12 -4.97580E-12 -5.04420E-12 -5.11090E-12 &
     -5.17590E-12 -5.23920E-12 -5.30050E-12 -5.35990E-12 -5.41720E-12 -5.47250E-12 -5.52560E-12 -5.57660E-12 &
     -5.62520E-12 -5.67150E-12 -5.71530E-12 -5.75670E-12 -5.79550E-12 -5.83170E-12 -5.86520E-12 -5.89600E-12 &
     -5.92400E-12 -5.94910E-12 -5.97140E-12 -5.99070E-12 -6.00690E-12 -6.02010E-12 -6.03020E-12 -6.03720E-12 &
     -6.04090E-12 -6.04140E-12 -6.03870E-12 -6.03260E-12 -6.02320E-12 -6.01050E-12 -5.99440E-12 -5.97500E-12 &
     -5.95210E-12 -5.92580E-12 -5.89600E-12 -5.86280E-12 -5.82610E-12 -5.78580E-12 -5.74200E-12 -5.69460E-12 &
     -5.64360E-12 -5.58910E-12 -5.53090E-12 -5.46910E-12 -5.40370E-12 -5.33460E-12 -5.26180E-12 -5.18540E-12 &
     -5.10540E-12 -5.02170E-12 -4.93430E-12 -4.84330E-12 -4.74860E-12 -4.65030E-12 -4.54830E-12 -4.44280E-12 &
     -4.33360E-12 -4.22090E-12 -4.10460E-12 -3.98480E-12 -3.86150E-12 -3.73470E-12 -3.60450E-12 -3.47090E-12 &
     -3.33390E-12 -3.19350E-12 -3.04980E-12 -2.90290E-12 -2.75270E-12 -2.59940E-12 -2.44290E-12 -2.28340E-12 &
     -2.12080E-12 -1.95520E-12 -1.78670E-12 -1.61530E-12 -1.44110E-12 -1.26400E-12 -1.08430E-12 -9.01910E-13 &
     -7.16890E-13 -5.29320E-13 -3.39270E-13 -1.46800E-13 4.80080E-14 2.45080E-13 4.44340E-13 6.45730E-13 &
     8.49160E-13 1.05450E-12 1.26180E-12 1.47090E-12 1.68170E-12 1.89420E-12 2.10820E-12 2.32370E-12 &
     2.54060E-12 2.75880E-12 2.97820E-12 3.19880E-12 3.42040E-12 3.64300E-12 3.86640E-12 4.09050E-12 &
     4.31540E-12 4.54080E-12 4.76680E-12 4.99310E-12 5.21970E-12 5.44650E-12 5.67340E-12 5.90030E-12 &
     6.12710E-12 6.35380E-12 6.58010E-12 6.80610E-12 7.03150E-12 7.25640E-12 7.48060E-12 7.70400E-12 &
     7.92650E-12 8.14810E-12 8.36860E-12 8.58790E-12 8.80590E-12 9.02260E-12 9.23780E-12 9.45140E-12 &
     9.66330E-12 9.87340E-12 1.00820E-11 1.02880E-11 1.04920E-11 1.06940E-11 1.08940E-11 1.10910E-11 &
     1.12860E-11 1.14790E-11 1.16690E-11 1.18560E-11 1.20410E-11 1.22220E-11 1.24010E-11 1.25770E-11 &
     1.27500E-11 1.29190E-11 1.30850E-11 1.32480E-11 1.34070E-11 1.35630E-11 1.37150E-11 1.38640E-11 &
     1.40080E-11 1.41490E-11 1.42870E-11 1.44200E-11 1.45490E-11 1.46740E-11 1.47950E-11 1.49120E-11 &
     1.50240E-11 1.51320E-11 1.52360E-11 1.53360E-11 1.54310E-11 1.55210E-11 1.56070E-11 1.56880E-11 &
     1.57650E-11 1.58360E-11 1.59030E-11 1.59660E-11 1.60230E-11 1.60750E-11 1.61220E-11 1.61650E-11 &
     1.62020E-11 1.62350E-11 1.62620E-11 1.62850E-11 1.63020E-11 1.63140E-11 1.63220E-11 1.63240E-11 &
     1.63210E-11 1.63130E-11 1.63000E-11 1.62810E-11 1.62580E-11 1.62290E-11 1.61950E-11 1.61560E-11 &
     1.61120E-11 1.60630E-11 1.60090E-11 1.59490E-11 1.58850E-11 1.58160E-11 1.57410E-11 1.56620E-11 &
     1.55770E-11 1.54870E-11 1.53930E-11 1.52940E-11 1.51890E-11 1.50800E-11 1.49670E-11 1.48480E-11 &
     1.47250E-11 1.45970E-11 1.44640E-11 1.43270E-11 1.41850E-11 1.40390E-11 1.38890E-11 1.37340E-11 &
     1.35750E-11 1.34120E-11 1.32440E-11 1.30720E-11 1.28970E-11 1.27170E-11 1.25330E-11 1.23450E-11 &
     1.21540E-11 1.19590E-11 1.17600E-11 1.15580E-11 1.13520E-11 1.11430E-11 1.09310E-11 1.07160E-11 &
     1.04980E-11 1.02760E-11 1.00520E-11 9.82460E-12 9.59450E-12 9.36170E-12 9.12640E-12 8.88860E-12 &
     8.64840E-12 8.40590E-12 8.16130E-12 7.91450E-12 7.66590E-12 7.41530E-12 7.16300E-12 6.90900E-12 &
     6.65340E-12 6.39640E-12 6.13800E-12 5.87840E-12 5.61750E-12 5.35570E-12 5.09280E-12 4.82910E-12 &
     4.56470E-12 4.29950E-12 4.03390E-12 3.76780E-12 3.50140E-12 3.23480E-12 2.96810E-12 2.70140E-12 &
     2.43490E-12 2.16860E-12 1.90260E-12 1.63710E-12 1.37210E-12 1.10780E-12 8.44280E-13 5.81640E-13 &
     3.19990E-13 5.94170E-14 -1.99960E-13 -4.58050E-13 -7.14730E-13 -9.69920E-13 -1.22350E-12 -1.47540E-12 &
     -1.72540E-12 -1.97360E-12 -2.21970E-12 -2.46370E-12 -2.70550E-12 -2.94490E-12 -3.18200E-12 -3.41660E-12 &
     -3.64860E-12 -3.87790E-12 -4.10450E-12 -4.32830E-12 -4.54910E-12 -4.76690E-12 -4.98160E-12 -5.19310E-12 &
     -5.40140E-12 -5.60620E-12 -5.80760E-12 -6.00550E-12 -6.19970E-12 -6.39040E-12 -6.57730E-12 -6.76050E-12 &
     -6.93990E-12 -7.11540E-12 -7.28700E-12 -7.45450E-12 -7.61810E-12 -7.77750E-12 -7.93280E-12 -8.08380E-12 &
     -8.23050E-12 -8.37290E-12 -8.51090E-12 -8.64450E-12 -8.77360E-12 -8.89820E-12 -9.01820E-12 -9.13360E-12 &
     -9.24440E-12 -9.35060E-12 -9.45210E-12 -9.54890E-12 -9.64090E-12 -9.72820E-12 -9.81080E-12 -9.88860E-12 &
     -9.96160E-12 -1.00300E-11 -1.00930E-11 -1.01520E-11 -1.02050E-11 -1.02540E-11 -1.02980E-11 -1.03370E-11 &
     -1.03710E-11 -1.04010E-11 -1.04250E-11 -1.04450E-11 -1.04590E-11 -1.04690E-11 -1.04740E-11 -1.04750E-11 &
     -1.04700E-11 -1.04610E-11 -1.04470E-11 -1.04290E-11 -1.04060E-11 -1.03780E-11 -1.03460E-11 -1.03100E-11 &
     -1.02680E-11 -1.02230E-11 -1.01730E-11 -1.01180E-11 -1.00590E-11 -9.99620E-12 -9.92880E-12 -9.85720E-12 &
     -9.78150E-12 -9.70170E-12 -9.61790E-12 -9.53010E-12 -9.43830E-12 -9.34280E-12 -9.24350E-12 -9.14050E-12 &
     -9.03380E-12 -8.92360E-12 -8.80980E-12 -8.69270E-12 -8.57220E-12 -8.44830E-12 -8.32130E-12 -8.19110E-12 &
     -8.05790E-12 -7.92170E-12 -7.78260E-12 -7.64080E-12 -7.49630E-12 -7.34910E-12 -7.19940E-12 -7.04740E-12 &
     -6.89300E-12 -6.73630E-12 -6.57750E-12 -6.41670E-12 -6.25390E-12 -6.08920E-12 -5.92270E-12 -5.75460E-12 &
     -5.58490E-12 -5.41370E-12 -5.24110E-12 -5.06730E-12 -4.89220E-12 -4.71610E-12 -4.53900E-12 -4.36090E-12 &
     -4.18210E-12 -4.00260E-12 -3.82250E-12 -3.64190E-12 -3.46090E-12 -3.27960E-12 -3.09810E-12 -2.91660E-12 &
     -2.73510E-12 -2.55370E-12 -2.37260E-12 -2.19180E-12 -2.01150E-12 -1.83170E-12 -1.65260E-12 -1.47420E-12 &
     -1.29670E-12 -1.12010E-12 -9.44590E-13 -7.70210E-13 -5.97060E-13 -4.25230E-13 -2.54830E-13 -8.59300E-14 &
     8.13590E-14 2.46950E-13 4.10740E-13 5.72640E-13 7.32550E-13 8.90370E-13 1.04600E-12 1.19940E-12 &
     1.35040E-12 1.49890E-12 1.64490E-12 1.78830E-12 1.92900E-12 2.06700E-12 2.20200E-12 2.33410E-12 &
     2.46320E-12 2.58920E-12 2.71200E-12 2.83170E-12 2.94800E-12 3.06090E-12 3.17040E-12 3.27630E-12 &
     3.37880E-12 3.47760E-12 3.57280E-12 3.66420E-12 3.75190E-12 3.83580E-12 3.91570E-12 3.99180E-12 &
     4.06380E-12 4.13180E-12 4.19570E-12 4.25550E-12 4.31110E-12 4.36250E-12 4.40960E-12 4.45250E-12 &
     4.49100E-12 4.52520E-12 4.55490E-12 4.58030E-12 4.60120E-12 4.61770E-12 4.62960E-12 4.63710E-12 &
     4.64000E-12 4.63840E-12 4.63230E-12 4.62160E-12 4.60630E-12 4.58650E-12 4.56210E-12 4.53310E-12 &
     4.49960E-12 4.46150E-12 4.41880E-12 4.37160E-12 4.31980E-12 4.26350E-12 4.20270E-12 4.13740E-12 &
     4.06750E-12 3.99320E-12 3.91450E-12 3.83140E-12 3.74390E-12 3.65210E-12 3.55610E-12 3.45570E-12 &
     3.35110E-12 3.24240E-12 3.12950E-12 3.01240E-12 2.89130E-12 2.76620E-12 2.63710E-12 2.50400E-12 &
     2.36700E-12 2.22620E-12 2.08160E-12 1.93330E-12 1.78130E-12 1.62570E-12 1.46660E-12 1.30390E-12 &
     1.13790E-12 9.68440E-13 7.95730E-13 6.19790E-13 4.40700E-13 2.58540E-13 7.33700E-14 -1.14720E-13 &
     -3.05660E-13 -4.99370E-13 -6.95780E-13 -8.94800E-13 -1.09630E-12 -1.30030E-12 -1.50670E-12 -1.71530E-12 &
     -1.92600E-12 -2.13890E-12 -2.35380E-12 -2.57060E-12 -2.78920E-12 -3.00950E-12 -3.23140E-12 -3.45490E-12 &
     -3.67980E-12 -3.90610E-12 -4.13360E-12 -4.36230E-12 -4.59200E-12 -4.82260E-12 -5.05420E-12 -5.28650E-12 &
     -5.51940E-12 -5.75290E-12 -5.98690E-12 -6.22120E-12 -6.45580E-12 -6.69060E-12 -6.92540E-12 -7.16020E-12 &
     -7.39490E-12 -7.62930E-12 -7.86340E-12 -8.09690E-12 -8.33000E-12 -8.56230E-12 -8.79390E-12 -9.02460E-12 &
     -9.25440E-12 -9.48320E-12 -9.71080E-12 -9.93720E-12 -1.01620E-11 -1.03860E-11 -1.06080E-11 -1.08290E-11 &
     -1.10470E-11 -1.12640E-11 -1.14800E-11 -1.16930E-11 -1.19040E-11 -1.21130E-11 -1.23190E-11 -1.25240E-11 &
     -1.27250E-11 -1.29250E-11 -1.31210E-11 -1.33150E-11 -1.35070E-11 -1.36950E-11 -1.38800E-11 -1.40630E-11 &
     -1.42420E-11 -1.44180E-11 -1.45910E-11 -1.47600E-11 -1.49260E-11 -1.50890E-11 -1.52480E-11 -1.54030E-11 &
     -1.55550E-11 -1.57030E-11 -1.58470E-11 -1.59870E-11 -1.61230E-11 -1.62550E-11 -1.63840E-11 -1.65080E-11 &
     -1.66280E-11 -1.67430E-11 -1.68550E-11 -1.69620E-11 -1.70650E-11 -1.71630E-11 -1.72570E-11 -1.73470E-11 &
     -1.74310E-11 -1.75120E-11 -1.75870E-11 -1.76580E-11 -1.77250E-11 -1.77860E-11 -1.78430E-11 -1.78950E-11 &
     -1.79420E-11 -1.79840E-11 -1.80210E-11 -1.80540E-11 -1.80810E-11 -1.81040E-11 -1.81210E-11 -1.81340E-11 &
     -1.81420E-11 -1.81450E-11 -1.81420E-11 -1.81350E-11 -1.81230E-11 -1.81060E-11 -1.80840E-11 -1.80570E-11 &
     -1.80250E-11 -1.79880E-11 -1.79470E-11 -1.79000E-11 -1.78490E-11 -1.77930E-11 -1.77320E-11 -1.76660E-11 &
     -1.75950E-11 -1.75190E-11 -1.74380E-11 -1.73530E-11 -1.72630E-11 -1.71680E-11 -1.70680E-11 -1.69640E-11 &
     -1.68550E-11 -1.67420E-11 -1.66250E-11 -1.65030E-11 -1.63760E-11 -1.62460E-11 -1.61110E-11 -1.59720E-11 &
     -1.58290E-11 -1.56820E-11 -1.55310E-11 -1.53760E-11 -1.52160E-11 -1.50530E-11 -1.48860E-11 -1.47150E-11 &
     -1.45400E-11 -1.43610E-11 -1.41790E-11 -1.39930E-11 -1.38040E-11 -1.36120E-11 -1.34170E-11 -1.32180E-11 &
     -1.30170E-11 -1.28130E-11 -1.26050E-11 -1.23950E-11 -1.21820E-11 -1.19660E-11 -1.17470E-11 -1.15260E-11 &
     -1.13020E-11 -1.10760E-11 -1.08470E-11 -1.06170E-11 -1.03840E-11 -1.01490E-11 -9.91190E-12 -9.67320E-12 &
     -9.43290E-12 -9.19090E-12 -8.94760E-12 -8.70280E-12 -8.45670E-12 -8.20930E-12 -7.96060E-12 -7.71080E-12 &
     -7.45980E-12 -7.20770E-12 -6.95460E-12 -6.70080E-12 -6.44630E-12 -6.19130E-12 -5.93590E-12 -5.68040E-12 &
     -5.42480E-12 -5.16920E-12 -4.91370E-12 -4.65820E-12 -4.40290E-12 -4.14770E-12 -3.89280E-12 -3.63820E-12 &
     -3.38400E-12 -3.13030E-12 -2.87730E-12 -2.62510E-12 -2.37380E-12 -2.12360E-12 -1.87450E-12 -1.62670E-12 &
     -1.38030E-12 -1.13520E-12 -8.91690E-13 -6.49700E-13 -4.09360E-13 -1.70730E-13 6.61160E-14 3.01110E-13 &
     5.34190E-13 7.65280E-13 9.94310E-13 1.22120E-12 1.44590E-12 1.66840E-12 1.88830E-12 2.10580E-12 &
     2.32060E-12 2.53260E-12 2.74180E-12 2.94790E-12 3.15110E-12 3.35120E-12 3.54820E-12 3.74220E-12 &
     3.93300E-12 4.12070E-12 4.30520E-12 4.48650E-12 4.66440E-12 4.83900E-12 5.01010E-12 5.17770E-12 &
     5.34160E-12 5.50190E-12 5.65840E-12 5.81110E-12 5.95980E-12 6.10460E-12 6.24540E-12 6.38210E-12 &
     6.51460E-12 6.64300E-12 6.76720E-12 6.88720E-12 7.00290E-12 7.11450E-12 7.22180E-12 7.32490E-12 &
     7.42370E-12 7.51830E-12 7.60850E-12 7.69440E-12 7.77600E-12 7.85320E-12 7.92600E-12 7.99440E-12 &
     8.05840E-12 8.11790E-12 8.17290E-12 8.22350E-12 8.26950E-12 8.31100E-12 8.34800E-12 8.38050E-12 &
     8.40850E-12 8.43200E-12 8.45100E-12 8.46540E-12 8.47540E-12 8.48090E-12 8.48190E-12 8.47850E-12 &
     8.47060E-12 8.45830E-12 8.44150E-12 8.42040E-12 8.39490E-12 8.36510E-12 8.33100E-12 8.29260E-12 &
     8.25000E-12 8.20330E-12 8.15240E-12 8.09740E-12 8.03830E-12 7.97520E-12 7.90810E-12 7.83710E-12 &
     7.76210E-12 7.68320E-12 7.60040E-12 7.51380E-12 7.42350E-12 7.32940E-12 7.23170E-12 7.13030E-12 &
     7.02530E-12 6.91680E-12 6.80480E-12 6.68950E-12 6.57070E-12 6.44870E-12 6.32340E-12 6.19500E-12 &
     6.06340E-12 5.92870E-12 5.79110E-12 5.65050E-12 5.50710E-12 5.36080E-12 5.21170E-12 5.06000E-12 &
     4.90570E-12 4.74880E-12 4.58960E-12 4.42790E-12 4.26410E-12 4.09800E-12 3.92980E-12 3.75950E-12 &
     3.58730E-12 3.41320E-12 3.23720E-12 3.05950E-12 2.88010E-12 2.69910E-12 2.51660E-12 2.33270E-12 &
     2.14730E-12 1.96070E-12 1.77290E-12 1.58400E-12 1.39410E-12 1.20320E-12 1.01140E-12 8.18750E-13 &
     6.25410E-13 4.31410E-13 2.36830E-13 4.17480E-14 -1.53770E-13 -3.49670E-13 -5.45860E-13 -7.42280E-13 &
     -9.38860E-13 -1.13550E-12 -1.33220E-12 -1.52880E-12 -1.72520E-12 -1.92140E-12 -2.11720E-12 -2.31270E-12 &
     -2.50760E-12 -2.70190E-12 -2.89560E-12 -3.08840E-12 -3.28040E-12 -3.47140E-12 -3.66130E-12 -3.85000E-12 &
     -4.03760E-12 -4.22390E-12 -4.40890E-12 -4.59250E-12 -4.77480E-12 -4.95560E-12 -5.13500E-12 -5.31280E-12 &
     -5.48910E-12 -5.66370E-12 -5.83660E-12 -6.00770E-12 -6.17700E-12 -6.34440E-12 -6.50980E-12 -6.67310E-12 &
     -6.83430E-12 -6.99310E-12 -7.14960E-12 -7.30360E-12 -7.45510E-12 -7.60390E-12 -7.75020E-12 -7.89400E-12 &
     -8.03520E-12 -8.17390E-12 -8.31010E-12 -8.44370E-12 -8.57470E-12 -8.70310E-12 -8.82870E-12 -8.95140E-12 &
     -9.07110E-12 -9.18790E-12 -9.30150E-12 -9.41190E-12 -9.51920E-12 -9.62340E-12 -9.72440E-12 -9.82220E-12 &
     -9.91690E-12 -1.00080E-11 -1.00970E-11 -1.01820E-11 -1.02630E-11 -1.03420E-11 -1.04170E-11 -1.04890E-11 &
     -1.05580E-11 -1.06230E-11 -1.06850E-11 -1.07440E-11 -1.07990E-11 -1.08510E-11 -1.08990E-11 -1.09440E-11 &
     -1.09860E-11 -1.10240E-11 -1.10590E-11 -1.10910E-11 -1.11190E-11 -1.11430E-11 -1.11640E-11 -1.11820E-11 &
     -1.11970E-11 -1.12080E-11 -1.12160E-11 -1.12200E-11 -1.12210E-11 -1.12190E-11 -1.12130E-11 -1.12050E-11 &
     -1.11920E-11 -1.11770E-11 -1.11590E-11 -1.11370E-11 -1.11120E-11 -1.10840E-11 -1.10530E-11 -1.10180E-11 &
     -1.09810E-11 -1.09400E-11 -1.08970E-11 -1.08500E-11 -1.08010E-11 -1.07480E-11 -1.06930E-11 -1.06340E-11 &
     -1.05730E-11 -1.05100E-11 -1.04430E-11 -1.03740E-11 -1.03020E-11 -1.02280E-11 -1.01510E-11 -1.00710E-11 &
     -9.98900E-12 -9.90460E-12 -9.81770E-12 -9.72860E-12 -9.63710E-12 -9.54340E-12 -9.44750E-12 -9.34950E-12 &
     -9.24940E-12 -9.14720E-12 -9.04310E-12 -8.93710E-12 -8.82920E-12 -8.71950E-12 -8.60810E-12 -8.49500E-12 &
     -8.38020E-12 -8.26380E-12 -8.14590E-12 -8.02650E-12 -7.90570E-12 -7.78360E-12 -7.66010E-12 -7.53540E-12 &
     -7.40950E-12 -7.28240E-12 -7.15440E-12 -7.02530E-12 -6.89530E-12 -6.76450E-12 -6.63290E-12 -6.50050E-12 &
     -6.36740E-12 -6.23370E-12 -6.09950E-12 -5.96470E-12 -5.82950E-12 -5.69390E-12 -5.55800E-12 -5.42180E-12 &
     -5.28540E-12 -5.14890E-12 -5.01230E-12 -4.87570E-12 -4.73910E-12 -4.60260E-12 -4.46630E-12 -4.33020E-12 &
     -4.19440E-12 -4.05890E-12 -3.92380E-12 -3.78920E-12 -3.65510E-12 -3.52150E-12 -3.38860E-12 -3.25630E-12 &
     -3.12480E-12 -2.99410E-12 -2.86420E-12 -2.73530E-12 -2.60730E-12 -2.48030E-12 -2.35440E-12 -2.22960E-12 &
     -2.10590E-12 -1.98350E-12 -1.86230E-12 -1.74250E-12 -1.62400E-12 -1.50690E-12 -1.39130E-12 -1.27710E-12 &
     -1.16450E-12 -1.05350E-12 -9.44170E-13 -8.36480E-13 -7.30510E-13 -6.26290E-13 -5.23870E-13 -4.23290E-13 &
     -3.24580E-13 -2.27790E-13 -1.32950E-13 -4.00960E-14 5.07320E-14 1.39500E-13 2.26170E-13 3.10720E-13 &
     3.93120E-13 4.73350E-13 5.51360E-13 6.27140E-13 7.00650E-13 7.71880E-13 8.40790E-13 9.07370E-13 &
     9.71600E-13 1.03340E-12 1.09290E-12 1.14990E-12 1.20450E-12 1.25660E-12 1.30630E-12 1.35350E-12 &
     1.39820E-12 1.44040E-12 1.48010E-12 1.51730E-12 1.55190E-12 1.58400E-12 1.61360E-12 1.64070E-12 &
     1.66520E-12 1.68710E-12 1.70650E-12 1.72330E-12 1.73760E-12 1.74940E-12 1.75860E-12 1.76520E-12 &
     1.76930E-12 1.77090E-12 1.76990E-12 1.76640E-12 1.76040E-12 1.75190E-12 1.74100E-12 1.72750E-12 &
     1.71160E-12 1.69330E-12 1.67260E-12 1.64940E-12 1.62390E-12 1.59600E-12 1.56570E-12 1.53310E-12 &
     1.49810E-12 1.46080E-12 1.42130E-12 1.37950E-12 1.33550E-12 1.28930E-12 1.24090E-12 1.19030E-12 &
     1.13770E-12 1.08290E-12 1.02600E-12 9.67110E-13 9.06180E-13 8.43270E-13 7.78410E-13 7.11640E-13 &
     6.43010E-13 5.72530E-13 5.00270E-13 4.26240E-13 3.50490E-13 2.73060E-13 1.93990E-13 1.13310E-13 &
     3.10590E-14 -5.27200E-14 -1.37990E-13 -2.24710E-13 -3.12840E-13 -4.02320E-13 -4.93130E-13 -5.85200E-13 &
     -6.78510E-13 -7.73010E-13 -8.68650E-13 -9.65390E-13 -1.06320E-12 -1.16200E-12 -1.26170E-12 -1.36240E-12 &
     -1.46390E-12 -1.56630E-12 -1.66950E-12 -1.77330E-12 -1.87790E-12 -1.98310E-12 -2.08890E-12 -2.19530E-12 &
     -2.30220E-12 -2.40950E-12 -2.51720E-12 -2.62520E-12 -2.73360E-12 -2.84230E-12 -2.95110E-12 -3.06010E-12 &
     -3.16920E-12 -3.27840E-12 -3.38760E-12 -3.49680E-12 -3.60590E-12 -3.71500E-12 -3.82380E-12 -3.93250E-12 &
     -4.04080E-12 -4.14890E-12 -4.25670E-12 -4.36410E-12 -4.47100E-12 -4.57740E-12 -4.68340E-12 -4.78870E-12 &
     -4.89350E-12 -4.99760E-12 -5.10100E-12 -5.20370E-12 -5.30550E-12 -5.40660E-12 -5.50680E-12 -5.60610E-12 &
     -5.70440E-12 -5.80180E-12 -5.89810E-12 -5.99330E-12 -6.08740E-12 -6.18040E-12 -6.27230E-12 -6.36290E-12 &
     -6.45220E-12 -6.54030E-12 -6.62700E-12 -6.71240E-12 -6.79650E-12 -6.87900E-12 -6.96020E-12 -7.03990E-12 &
     -7.11800E-12 -7.19460E-12 -7.26960E-12 -7.34300E-12 -7.41480E-12 -7.48490E-12 -7.55340E-12 -7.62010E-12 &
     -7.68510E-12 -7.74830E-12 -7.80980E-12 -7.86940E-12 -7.92720E-12 -7.98310E-12 -8.03720E-12 -8.08930E-12 &
     -8.13950E-12 -8.18780E-12 -8.23420E-12 -8.27860E-12 -8.32110E-12 -8.36150E-12 -8.40000E-12 -8.43650E-12 &
     -8.47100E-12 -8.50340E-12 -8.53380E-12 -8.56220E-12 -8.58840E-12 -8.61270E-12 -8.63490E-12 -8.65490E-12 &
     -8.67300E-12 -8.68890E-12 -8.70280E-12 -8.71450E-12 -8.72420E-12 -8.73180E-12 -8.73730E-12 -8.74080E-12 &
     -8.74210E-12 -8.74140E-12 -8.73860E-12 -8.73380E-12 -8.72690E-12 -8.71790E-12 -8.70690E-12 -8.69380E-12 &
     -8.67860E-12 -8.66140E-12 -8.64220E-12 -8.62100E-12 -8.59780E-12 -8.57260E-12 -8.54540E-12 -8.51630E-12 &
     -8.48530E-12 -8.45230E-12 -8.41740E-12 -8.38070E-12 -8.34210E-12 -8.30160E-12 -8.25930E-12 -8.21520E-12 &
     -8.16920E-12 -8.12150E-12 -8.07200E-12 -8.02080E-12 -7.96790E-12 -7.91330E-12 -7.85700E-12 -7.79910E-12 &
     -7.73960E-12 -7.67850E-12 -7.61590E-12 -7.55170E-12 -7.48600E-12 -7.41880E-12 -7.35010E-12 -7.28010E-12 &
     -7.20860E-12 -7.13580E-12 -7.06170E-12 -6.98630E-12 -6.90960E-12 -6.83170E-12 -6.75260E-12 -6.67240E-12 &
     -6.59100E-12 -6.50850E-12 -6.42500E-12 -6.34050E-12 -6.25490E-12 -6.16840E-12 -6.08100E-12 -5.99260E-12 &
     -5.90350E-12 -5.81350E-12 -5.72270E-12 -5.63110E-12 -5.53890E-12 -5.44590E-12 -5.35240E-12 -5.25820E-12 &
     -5.16350E-12 -5.06820E-12 -4.97240E-12 -4.87620E-12 -4.77960E-12 -4.68250E-12 -4.58510E-12 -4.48740E-12 &
     -4.38940E-12 -4.29120E-12 -4.19280E-12 -4.09420E-12 -3.99540E-12 -3.89660E-12 -3.79780E-12 -3.69890E-12 &
     -3.60000E-12 -3.50120E-12 -3.40260E-12 -3.30400E-12 -3.20560E-12 -3.10740E-12 -3.00940E-12 -2.91170E-12 &
     -2.81430E-12 -2.71730E-12 -2.62060E-12 -2.52430E-12 -2.42840E-12 -2.33310E-12 -2.23820E-12 -2.14390E-12 &
     -2.05010E-12 -1.95700E-12 -1.86450E-12 -1.77270E-12 -1.68150E-12 -1.59110E-12 -1.50150E-12 -1.41260E-12 &
     -1.32460E-12 -1.23740E-12 -1.15110E-12 -1.06570E-12 -9.81230E-13 -8.97740E-13 -8.15240E-13 -7.33750E-13 &
     -6.53300E-13 -5.73920E-13 -4.95630E-13 -4.18460E-13 -3.42430E-13 -2.67570E-13 -1.93890E-13 -1.21420E-13 &
     -5.01890E-14 1.97890E-14 8.84900E-14 1.55900E-13 2.21990E-13 2.86750E-13 3.50160E-13 4.12210E-13 &
     4.72860E-13 5.32110E-13 5.89950E-13 6.46350E-13 7.01310E-13 7.54800E-13 8.06810E-13 8.57340E-13 &
     9.06360E-13 9.53880E-13 9.99870E-13 1.04430E-12 1.08730E-12 1.12870E-12 1.16850E-12 1.20680E-12 &
     1.24360E-12 1.27870E-12 1.31230E-12 1.34440E-12 1.37480E-12 1.40370E-12 1.43090E-12 1.45660E-12 &
     1.48070E-12 1.50310E-12 1.52400E-12 1.54330E-12 1.56100E-12 1.57710E-12 1.59160E-12 1.60460E-12 &
     1.61590E-12 1.62570E-12 1.63390E-12 1.64050E-12 1.64560E-12 1.64910E-12 1.65100E-12 1.65140E-12 &
     1.65030E-12 1.64760E-12 1.64340E-12 1.63760E-12 1.63030E-12 1.62160E-12 1.61140E-12 1.59970E-12 &
     1.58650E-12 1.57190E-12 1.55590E-12 1.53850E-12 1.51960E-12 1.49940E-12 1.47780E-12 1.45490E-12 &
     1.43060E-12 1.40500E-12 1.37810E-12 1.34990E-12 1.32040E-12 1.28970E-12 1.25780E-12 1.22460E-12 &
     1.19030E-12 1.15480E-12 1.11810E-12 1.08030E-12 1.04140E-12 1.00150E-12 9.60420E-13 9.18340E-13 &
     8.75230E-13 8.31140E-13 7.86080E-13 7.40080E-13 6.93160E-13 6.45350E-13 5.96680E-13 5.47160E-13 &
     4.96820E-13 4.45690E-13 3.93780E-13 3.41130E-13 2.87760E-13 2.33680E-13 1.78930E-13 1.23540E-13 &
     6.75280E-14 1.09210E-14 -4.62550E-14 -1.03970E-13 -1.62210E-13 -2.20930E-13 -2.80110E-13 -3.39730E-13 &
     -3.99750E-13 -4.60160E-13 -5.20920E-13 -5.82020E-13 -6.43410E-13 -7.05080E-13 -7.67000E-13 -8.29150E-13 &
     -8.91510E-13 -9.54030E-13 -1.01670E-12 -1.07950E-12 -1.14240E-12 -1.20530E-12 -1.26830E-12 -1.33130E-12 &
     -1.39430E-12 -1.45720E-12 -1.52000E-12 -1.58280E-12 -1.64540E-12 -1.70790E-12 -1.77020E-12 -1.83230E-12 &
     -1.89420E-12 -1.95590E-12 -2.01730E-12 -2.07840E-12 -2.13920E-12 -2.19960E-12 -2.25970E-12 -2.31940E-12 &
     -2.37870E-12 -2.43750E-12 -2.49590E-12 -2.55380E-12 -2.61120E-12 -2.66810E-12 -2.72440E-12 -2.78020E-12 &
     -2.83540E-12 -2.89000E-12 -2.94400E-12 -2.99730E-12 -3.04990E-12 -3.10190E-12 -3.15320E-12 -3.20380E-12 &
     -3.25360E-12 -3.30260E-12 -3.35090E-12 -3.39850E-12 -3.44520E-12 -3.49110E-12 -3.53610E-12 -3.58030E-12 &
     -3.62370E-12 -3.66610E-12 -3.70770E-12 -3.74840E-12 -3.78820E-12 -3.82700E-12 -3.86490E-12 -3.90190E-12 &
     -3.93790E-12 -3.97290E-12 -4.00690E-12 -4.04000E-12 -4.07210E-12 -4.10310E-12 -4.13320E-12 -4.16220E-12 &
     -4.19010E-12 -4.21710E-12 -4.24300E-12 -4.26780E-12 -4.29160E-12 -4.31430E-12 -4.33600E-12 -4.35660E-12 &
     -4.37610E-12 -4.39460E-12 -4.41200E-12 -4.42830E-12 -4.44360E-12 -4.45770E-12 -4.47080E-12 -4.48280E-12 &
     -4.49370E-12 -4.50350E-12 -4.51220E-12 -4.51980E-12 -4.52640E-12 -4.53190E-12 -4.53630E-12 -4.53960E-12 &
     -4.54190E-12 -4.54310E-12 -4.54320E-12 -4.54230E-12 -4.54030E-12 -4.53720E-12 -4.53310E-12 -4.52800E-12 &
     -4.52180E-12 -4.51460E-12 -4.50640E-12 -4.49720E-12 -4.48700E-12 -4.47580E-12 -4.46360E-12 -4.45050E-12 &
     -4.43640E-12 -4.42130E-12 -4.40530E-12 -4.38840E-12 -4.37050E-12 -4.35180E-12 -4.33210E-12 -4.31160E-12 &
     -4.29020E-12 -4.26790E-12 -4.24480E-12 -4.22080E-12 -4.19610E-12 -4.17050E-12 -4.14420E-12 -4.11700E-12 &
     -4.08920E-12 -4.06050E-12 -4.03120E-12 -4.00110E-12 -3.97030E-12 -3.93890E-12 -3.90670E-12 -3.87390E-12 &
     -3.84050E-12 -3.80640E-12 -3.77180E-12 -3.73650E-12 -3.70070E-12 -3.66430E-12 -3.62740E-12 -3.59000E-12 &
     -3.55210E-12 -3.51370E-12 -3.47480E-12 -3.43550E-12 -3.39580E-12 -3.35560E-12 -3.31510E-12 -3.27410E-12 &
     -3.23280E-12 -3.19120E-12 -3.14920E-12 -3.10700E-12 -3.06440E-12 -3.02160E-12 -2.97860E-12 -2.93530E-12 &
     -2.89180E-12 -2.84810E-12 -2.80420E-12 -2.76020E-12 -2.71600E-12 -2.67170E-12 -2.62730E-12 -2.58280E-12 &
     -2.53830E-12 -2.49370E-12 -2.44910E-12 -2.40450E-12 -2.35990E-12 -2.31530E-12 -2.27080E-12 -2.22630E-12 &
     -2.18200E-12 -2.13770E-12 -2.09350E-12 -2.04950E-12 -2.00560E-12 -1.96190E-12 -1.91830E-12 -1.87500E-12 &
     -1.83190E-12 -1.78900E-12 -1.74640E-12 -1.70400E-12 -1.66190E-12 -1.62010E-12 -1.57870E-12 -1.53750E-12 &
     -1.49670E-12 -1.45630E-12 -1.41620E-12 -1.37650E-12 -1.33720E-12 -1.29830E-12 -1.25990E-12 -1.22190E-12 &
     -1.18430E-12 -1.14720E-12 -1.11060E-12 -1.07450E-12 -1.03890E-12 -1.00390E-12 -9.69300E-13 -9.35290E-13 &
     -9.01830E-13 -8.68940E-13 -8.36610E-13 -8.04870E-13 -7.73720E-13 -7.43170E-13 -7.13240E-13 -6.83920E-13 &
     -6.55230E-13 -6.27170E-13 -5.99770E-13 -5.73010E-13 -5.46910E-13 -5.21480E-13 -4.96720E-13 -4.72650E-13 &
     -4.49260E-13 -4.26560E-13 -4.04560E-13 -3.83270E-13 -3.62680E-13 -3.42800E-13 -3.23650E-13 -3.05210E-13 &
     -2.87490E-13 -2.70500E-13 -2.54250E-13 -2.38720E-13 -2.23940E-13 -2.09880E-13 -1.96570E-13 -1.84000E-13 &
     -1.72170E-13 -1.61090E-13 -1.50750E-13 -1.41150E-13 -1.32300E-13 -1.24190E-13 -1.16830E-13 -1.10210E-13 &
     -1.04330E-13 -9.91960E-14 -9.47980E-14 -9.11370E-14 -8.82100E-14 -8.60150E-14 -8.45510E-14 -8.38130E-14 &
     -8.37990E-14 -8.45030E-14 -8.59240E-14 -8.80560E-14 -9.08970E-14 -9.44410E-14 -9.86830E-14 -1.03620E-13 &
     -1.09230E-13 -1.15520E-13 -1.22480E-13 -1.30100E-13 -1.38380E-13 -1.47300E-13 -1.56870E-13 -1.67060E-13 &
     -1.77880E-13 -1.89320E-13 -2.01370E-13 -2.14020E-13 -2.27260E-13 -2.41080E-13 -2.55480E-13 -2.70440E-13 &
     -2.85950E-13 -3.02010E-13 -3.18610E-13 -3.35730E-13 -3.53360E-13 -3.71500E-13 -3.90130E-13 -4.09250E-13 &
     -4.28840E-13 -4.48900E-13 -4.69400E-13 -4.90350E-13 -5.11730E-13 -5.33520E-13 -5.55720E-13 -5.78300E-13 &
     -6.01270E-13 -6.24610E-13 -6.48300E-13 -6.72330E-13 -6.96700E-13 -7.21390E-13 -7.46380E-13 -7.71680E-13 &
     -7.97260E-13 -8.23110E-13 -8.49210E-13 -8.75570E-13 -9.02160E-13 -9.28970E-13 -9.55980E-13 -9.83200E-13 &
     -1.01060E-12 -1.03820E-12 -1.06590E-12 -1.09370E-12 -1.12170E-12 -1.14980E-12 -1.17800E-12 -1.20630E-12 &
     -1.23470E-12 -1.26320E-12 -1.29170E-12 -1.32020E-12 -1.34880E-12 -1.37730E-12 -1.40590E-12 -1.43450E-12 &
     -1.46300E-12 -1.49150E-12 -1.52000E-12 -1.54840E-12 -1.57670E-12 -1.60490E-12 -1.63300E-12 -1.66100E-12 &
     -1.68890E-12 -1.71670E-12 -1.74430E-12 -1.77170E-12 -1.79900E-12 -1.82610E-12 -1.85300E-12 -1.87970E-12 &
     -1.90620E-12 -1.93250E-12 -1.95850E-12 -1.98420E-12 -2.00980E-12 -2.03500E-12 -2.06000E-12 -2.08460E-12 &
     -2.10900E-12 -2.13300E-12 -2.15680E-12 -2.18020E-12 -2.20320E-12 -2.22600E-12 -2.24830E-12 -2.27030E-12 &
     -2.29190E-12 -2.31320E-12 -2.33410E-12 -2.35450E-12 -2.37460E-12 -2.39420E-12 -2.41340E-12 -2.43220E-12 &
     -2.45060E-12 -2.46850E-12 -2.48600E-12 -2.50300E-12 -2.51960E-12 -2.53570E-12 -2.55130E-12 -2.56650E-12 &
     -2.58120E-12 -2.59540E-12 -2.60910E-12 -2.62230E-12 -2.63500E-12 -2.64720E-12 -2.65890E-12 -2.67000E-12 &
     -2.68070E-12 -2.69090E-12 -2.70050E-12 -2.70960E-12 -2.71820E-12 -2.72620E-12 -2.73370E-12 -2.74070E-12 &
     -2.74720E-12 -2.75310E-12 -2.75850E-12 -2.76330E-12 -2.76760E-12 -2.77140E-12 -2.77460E-12 -2.77730E-12 &
     -2.77940E-12 -2.78100E-12 -2.78200E-12 -2.78250E-12 -2.78250E-12 -2.78190E-12 -2.78080E-12 -2.77910E-12 &
     -2.77700E-12 -2.77420E-12 -2.77100E-12 -2.76720E-12 -2.76290E-12 -2.75800E-12 -2.75260E-12 -2.74680E-12 &
     -2.74030E-12 -2.73340E-12 -2.72600E-12 -2.71810E-12 -2.70970E-12 -2.70070E-12 -2.69130E-12 -2.68150E-12 &
     -2.67110E-12 -2.66030E-12 -2.64900E-12 -2.63720E-12 -2.62500E-12 -2.61230E-12 -2.59920E-12 -2.58560E-12 &
     -2.57160E-12 -2.55720E-12 -2.54240E-12 -2.52710E-12 -2.51140E-12 -2.49540E-12 -2.47890E-12 -2.46200E-12 &
     -2.44480E-12 -2.42720E-12 -2.40920E-12 -2.39090E-12 -2.37220E-12 -2.35310E-12 -2.33370E-12 -2.31400E-12 &
     -2.29400E-12 -2.27370E-12 -2.25300E-12 -2.23210E-12 -2.21090E-12 -2.18940E-12 -2.16760E-12 -2.14550E-12 &
     -2.12320E-12 -2.10070E-12 -2.07790E-12 -2.05490E-12 -2.03170E-12 -2.00830E-12 -1.98460E-12 -1.96080E-12 &
     -1.93680E-12 -1.91260E-12 -1.88820E-12 -1.86370E-12 -1.83910E-12 -1.81430E-12 -1.78940E-12 -1.76430E-12 &
     -1.73920E-12 -1.71400E-12 -1.68860E-12 -1.66320E-12 -1.63770E-12 -1.61210E-12 -1.58650E-12 -1.56090E-12 &
     -1.53520E-12 -1.50950E-12 -1.48370E-12 -1.45800E-12 -1.43220E-12 -1.40650E-12 -1.38080E-12 -1.35510E-12 &
     -1.32940E-12 -1.30380E-12 -1.27830E-12 -1.25280E-12 -1.22740E-12 -1.20200E-12 -1.17680E-12 -1.15160E-12 &
     -1.12650E-12 -1.10160E-12 -1.07680E-12 -1.05210E-12 -1.02760E-12 -1.00310E-12 -9.78900E-13 -9.54810E-13 &
     -9.30880E-13 -9.07130E-13 -8.83570E-13 -8.60200E-13 -8.37020E-13 -8.14050E-13 -7.91290E-13 -7.68740E-13 &
     -7.46430E-13 -7.24350E-13 -7.02500E-13 -6.80910E-13 -6.59560E-13 -6.38470E-13 -6.17640E-13 -5.97090E-13 &
     -5.76800E-13 -5.56800E-13 -5.37080E-13 -5.17650E-13 -4.98520E-13 -4.79690E-13 -4.61160E-13 -4.42940E-13 &
     -4.25040E-13 -4.07460E-13 -3.90200E-13 -3.73260E-13 -3.56660E-13 -3.40400E-13 -3.24470E-13 -3.08890E-13 &
     -2.93650E-13 -2.78760E-13 -2.64230E-13 -2.50050E-13 -2.36230E-13 -2.22780E-13 -2.09690E-13 -1.96970E-13 &
     -1.84610E-13 -1.72630E-13 -1.61020E-13 -1.49790E-13 -1.38930E-13 -1.28450E-13 -1.18350E-13 -1.08630E-13 &
     -9.92900E-14 -9.03350E-14 -8.17640E-14 -7.35800E-14 -6.57810E-14 -5.83680E-14 -5.13400E-14 -4.46980E-14 &
     -3.84410E-14 -3.25700E-14 -2.70830E-14 -2.19800E-14 -1.72610E-14 -1.29240E-14 -8.96820E-15 -5.39300E-15 &
     -2.19750E-15 6.19650E-16 3.06060E-15 5.12810E-15 6.82460E-15 8.15280E-15 9.11520E-15 9.71440E-15 &
     9.95310E-15 9.83390E-15 9.35950E-15 8.53270E-15 7.35620E-15 5.83260E-15 3.96470E-15 1.75540E-15 &
     -7.91650E-16 -3.67250E-15 -6.88330E-15 -1.04200E-14 -1.42790E-14 -1.84570E-14 -2.29480E-14 -2.77500E-14 &
     -3.28580E-14 -3.82680E-14 -4.39760E-14 -4.99780E-14 -5.62690E-14 -6.28460E-14 -6.97020E-14 -7.68340E-14 &
     -8.42350E-14 -9.19010E-14 -9.98250E-14 -1.08000E-13 -1.16430E-13 -1.25100E-13 -1.34010E-13 -1.43150E-13 &
     -1.52520E-13 -1.62110E-13 -1.71920E-13 -1.81940E-13 -1.92170E-13 -2.02600E-13 -2.13220E-13 -2.24030E-13 &
     -2.35020E-13 -2.46190E-13 -2.57520E-13 -2.69020E-13 -2.80680E-13 -2.92490E-13 -3.04450E-13 -3.16550E-13 &
     -3.28770E-13 -3.41130E-13 -3.53610E-13 -3.66200E-13 -3.78900E-13 -3.91700E-13 -4.04600E-13 -4.17580E-13 &
     -4.30660E-13 -4.43800E-13 -4.57030E-13 -4.70310E-13 -4.83650E-13 -4.97050E-13 -5.10490E-13 -5.23970E-13 &
     -5.37490E-13 -5.51030E-13 -5.64590E-13 -5.78170E-13 -5.91760E-13 -6.05350E-13 -6.18930E-13 -6.32510E-13 &
     -6.46080E-13 -6.59620E-13 -6.73140E-13 -6.86620E-13 -7.00070E-13 -7.13470E-13 -7.26830E-13 -7.40130E-13 &
     -7.53370E-13 -7.66540E-13 -7.79650E-13 -7.92680E-13 -8.05620E-13 -8.18480E-13 -8.31250E-13 -8.43920E-13 &
     -8.56490E-13 -8.68950E-13 -8.81300E-13 -8.93530E-13 -9.05650E-13 -9.17640E-13 -9.29490E-13 -9.41220E-13 &
     -9.52800E-13 -9.64240E-13 -9.75530E-13 -9.86670E-13 -9.97660E-13 -1.00850E-12 -1.01910E-12 -1.02960E-12 &
     -1.04000E-12 -1.05010E-12 -1.06010E-12 -1.06990E-12 -1.07950E-12 -1.08890E-12 -1.09810E-12 -1.10710E-12 &
     -1.11600E-12 -1.12460E-12 -1.13300E-12 -1.14120E-12 -1.14920E-12 -1.15700E-12 -1.16460E-12 -1.17200E-12 &
     -1.17910E-12 -1.18600E-12 -1.19270E-12 -1.19920E-12 -1.20550E-12 -1.21150E-12 -1.21730E-12 -1.22280E-12 &
     -1.22820E-12 -1.23330E-12 -1.23810E-12 -1.24270E-12 -1.24710E-12 -1.25130E-12 -1.25520E-12 -1.25880E-12 &
     -1.26220E-12 -1.26540E-12 -1.26840E-12 -1.27110E-12 -1.27350E-12 -1.27570E-12 -1.27770E-12 -1.27940E-12 &
     -1.28090E-12 -1.28220E-12 -1.28320E-12 -1.28390E-12 -1.28440E-12 -1.28470E-12 -1.28480E-12 -1.28460E-12 &
     -1.28420E-12 -1.28350E-12 -1.28260E-12 -1.28150E-12 -1.28010E-12 -1.27850E-12 -1.27670E-12 -1.27470E-12 &
     -1.27240E-12 -1.27000E-12 -1.26730E-12 -1.26430E-12 -1.26120E-12 -1.25780E-12 -1.25430E-12 -1.25050E-12 &
     -1.24650E-12 -1.24230E-12 -1.23790E-12 -1.23330E-12 -1.22850E-12 -1.22350E-12 -1.21840E-12 -1.21300E-12 &
     -1.20740E-12 -1.20170E-12 -1.19570E-12 -1.18960E-12 -1.18330E-12 -1.17690E-12 -1.17020E-12 -1.16340E-12 &
     -1.15650E-12 -1.14940E-12 -1.14210E-12 -1.13470E-12 -1.12710E-12 -1.11930E-12 -1.11150E-12 -1.10350E-12 &
     -1.09530E-12 -1.08700E-12 -1.07860E-12 -1.07010E-12 -1.06150E-12 -1.05270E-12 -1.04380E-12 -1.03480E-12 &
     -1.02570E-12 -1.01650E-12 -1.00730E-12 -9.97870E-13 -9.88400E-13 -9.78840E-13 -9.69210E-13 -9.59490E-13 &
     -9.49700E-13 -9.39840E-13 -9.29910E-13 -9.19930E-13 -9.09880E-13 -8.99780E-13 -8.89630E-13 -8.79440E-13 &
     -8.69210E-13 -8.58940E-13 -8.48640E-13 -8.38310E-13 -8.27950E-13 -8.17580E-13 -8.07190E-13 -7.96780E-13 &
     -7.86370E-13 -7.75950E-13 -7.65530E-13 -7.55110E-13 -7.44700E-13 -7.34300E-13 -7.23910E-13 -7.13540E-13 &
     -7.03190E-13 -6.92870E-13 -6.82570E-13 -6.72310E-13 -6.62080E-13 -6.51900E-13 -6.41750E-13 -6.31650E-13 &
     -6.21610E-13 -6.11610E-13 -6.01670E-13 -5.91790E-13 -5.81980E-13 -5.72240E-13 -5.62560E-13 -5.52960E-13 &
     -5.43430E-13 -5.33990E-13 -5.24630E-13 -5.15350E-13 -5.06160E-13 -4.97060E-13 -4.88060E-13 -4.79160E-13 &
     -4.70350E-13 -4.61650E-13 -4.53060E-13 -4.44570E-13 -4.36190E-13 -4.27930E-13 -4.19780E-13 -4.11750E-13 &
     -4.03830E-13 -3.96050E-13 -3.88380E-13 -3.80840E-13 -3.73430E-13 -3.66150E-13 -3.59000E-13 -3.51990E-13 &
     -3.45110E-13 -3.38370E-13 -3.31770E-13 -3.25310E-13 -3.18990E-13 -3.12820E-13 -3.06790E-13 -3.00910E-13 &
     -2.95170E-13 -2.89590E-13 -2.84150E-13 -2.78870E-13 -2.73740E-13 -2.68770E-13 -2.63940E-13 -2.59280E-13 &
     -2.54770E-13 -2.50410E-13 -2.46220E-13 -2.42180E-13 -2.38300E-13 -2.34570E-13 -2.31010E-13 -2.27600E-13 &
     -2.24360E-13 -2.21270E-13 -2.18340E-13 -2.15570E-13 -2.12960E-13 -2.10510E-13 -2.08220E-13 -2.06090E-13 &
     -2.04110E-13 -2.02300E-13 -2.00630E-13 -1.99130E-13 -1.97780E-13 -1.96590E-13 -1.95550E-13 -1.94660E-13 &
     -1.93930E-13 -1.93350E-13 -1.92920E-13 -1.92640E-13 -1.92510E-13 -1.92530E-13 -1.92700E-13 -1.93010E-13 &
     -1.93460E-13 -1.94060E-13 -1.94800E-13 -1.95680E-13 -1.96700E-13 -1.97860E-13 -1.99150E-13 -2.00580E-13 &
     -2.02140E-13 -2.03830E-13 -2.05650E-13 -2.07600E-13 -2.09670E-13 -2.11870E-13 -2.14190E-13 -2.16640E-13 &
     -2.19200E-13 -2.21870E-13 -2.24660E-13 -2.27570E-13 -2.30580E-13 -2.33700E-13 -2.36920E-13 -2.40250E-13 &
     -2.43680E-13 -2.47210E-13 -2.50840E-13 -2.54560E-13 -2.58370E-13 -2.62280E-13 -2.66270E-13 -2.70340E-13 &
     -2.74490E-13 -2.78730E-13 -2.83040E-13 -2.87430E-13 -2.91890E-13 -2.96410E-13 -3.01010E-13 -3.05670E-13 &
     -3.10390E-13 -3.15170E-13 -3.20000E-13 -3.24890E-13 -3.29830E-13 -3.34820E-13 -3.39860E-13 -3.44930E-13 &
     -3.50050E-13 -3.55210E-13 -3.60400E-13 -3.65630E-13 -3.70880E-13 -3.76160E-13 -3.81460E-13 -3.86790E-13 &
     -3.92140E-13 -3.97500E-13 -4.02870E-13 -4.08260E-13 -4.13650E-13 -4.19050E-13 -4.24460E-13 -4.29860E-13 &
     -4.35270E-13 -4.40660E-13 -4.46050E-13 -4.51440E-13 -4.56810E-13 -4.62160E-13 -4.67500E-13 -4.72820E-13 &
     -4.78120E-13 -4.83390E-13 -4.88630E-13 -4.93850E-13 -4.99040E-13 -5.04190E-13 -5.09310E-13 -5.14380E-13 &
     -5.19420E-13 -5.24420E-13 -5.29370E-13 -5.34280E-13 -5.39130E-13 -5.43940E-13 -5.48690E-13 -5.53390E-13 &
     -5.58030E-13 -5.62610E-13 -5.67130E-13 -5.71580E-13 -5.75980E-13 -5.80300E-13 -5.84560E-13 -5.88740E-13 &
     -5.92860E-13 -5.96900E-13 -6.00860E-13 -6.04750E-13 -6.08560E-13 -6.12290E-13 -6.15940E-13 -6.19510E-13 &
     -6.22990E-13 -6.26380E-13 -6.29690E-13 -6.32900E-13 -6.36030E-13 -6.39060E-13 -6.42000E-13 -6.44850E-13 &
     -6.47600E-13 -6.50260E-13 -6.52820E-13 -6.55290E-13 -6.57650E-13 -6.59920E-13 -6.62090E-13 -6.64150E-13 &
     -6.66120E-13 -6.67980E-13 -6.69730E-13 -6.71390E-13 -6.72940E-13 -6.74380E-13 -6.75720E-13 -6.76960E-13 &
     -6.78090E-13 -6.79110E-13 -6.80030E-13 -6.80840E-13 -6.81540E-13 -6.82140E-13 -6.82630E-13 -6.83010E-13 &
     -6.83280E-13 -6.83450E-13 -6.83510E-13 -6.83460E-13 -6.83300E-13 -6.83050E-13 -6.82680E-13 -6.82210E-13 &
     -6.81630E-13 -6.80950E-13 -6.80170E-13 -6.79280E-13 -6.78290E-13 -6.77190E-13 -6.75990E-13 -6.74690E-13 &
     -6.73290E-13 -6.71780E-13 -6.70180E-13 -6.68480E-13 -6.66670E-13 -6.64770E-13 -6.62770E-13 -6.60680E-13 &
     -6.58490E-13 -6.56200E-13 -6.53820E-13 -6.51350E-13 -6.48790E-13 -6.46130E-13 -6.43380E-13 -6.40550E-13 &
     -6.37630E-13 -6.34620E-13 -6.31530E-13 -6.28350E-13 -6.25100E-13 -6.21760E-13 -6.18340E-13 -6.14840E-13 &
     -6.11270E-13 -6.07620E-13 -6.03900E-13 -6.00100E-13 -5.96240E-13 -5.92300E-13 -5.88300E-13 -5.84230E-13 &
     -5.80090E-13 -5.75890E-13 -5.71630E-13 -5.67310E-13 -5.62930E-13 -5.58490E-13 -5.54000E-13 -5.49450E-13 &
     -5.44860E-13 -5.40210E-13 -5.35520E-13 -5.30780E-13 -5.25990E-13 -5.21160E-13 -5.16280E-13 -5.11370E-13 &
     -5.06420E-13 -5.01430E-13 -4.96410E-13 -4.91350E-13 -4.86270E-13 -4.81150E-13 -4.76010E-13 -4.70850E-13 &
     -4.65660E-13 -4.60450E-13 -4.55210E-13 -4.49960E-13 -4.44700E-13 -4.39410E-13 -4.34120E-13 -4.28810E-13 &
     -4.23490E-13 -4.18170E-13 -4.12840E-13 -4.07500E-13 -4.02160E-13 -3.96830E-13 -3.91490E-13 -3.86150E-13 &
     -3.80820E-13 -3.75490E-13 -3.70170E-13 -3.64860E-13 -3.59560E-13 -3.54280E-13 -3.49000E-13 -3.43750E-13 &
     -3.38510E-13 -3.33290E-13 -3.28090E-13 -3.22910E-13 -3.17760E-13 -3.12630E-13 -3.07530E-13 -3.02460E-13 &
     -2.97420E-13 -2.92410E-13 -2.87430E-13 -2.82490E-13 -2.77580E-13 -2.72710E-13 -2.67870E-13 -2.63080E-13 &
     -2.58330E-13 -2.53620E-13 -2.48950E-13 -2.44330E-13 -2.39760E-13 -2.35230E-13 -2.30750E-13 -2.26320E-13 &
     -2.21940E-13 -2.17620E-13 -2.13350E-13 -2.09130E-13 -2.04960E-13 -2.00860E-13 -1.96810E-13 -1.92820E-13 &
     -1.88880E-13 -1.85010E-13 -1.81200E-13 -1.77450E-13 -1.73770E-13 -1.70140E-13 -1.66580E-13 -1.63090E-13 &
     -1.59660E-13 -1.56300E-13 -1.53010E-13 -1.49780E-13 -1.46620E-13 -1.43530E-13 -1.40510E-13 -1.37560E-13 &
     -1.34680E-13 -1.31860E-13 -1.29120E-13 -1.26460E-13 -1.23860E-13 -1.21330E-13 -1.18880E-13 -1.16500E-13 &
     -1.14190E-13 -1.11960E-13 -1.09800E-13 -1.07710E-13 -1.05690E-13 -1.03750E-13 -1.01880E-13 -1.00090E-13 &
     -9.83690E-14 -9.67220E-14 -9.51480E-14 -9.36480E-14 -9.22210E-14 -9.08670E-14 -8.95860E-14 -8.83790E-14 &
     -8.72440E-14 -8.61820E-14 -8.51910E-14 -8.42730E-14 -8.34270E-14 -8.26520E-14 -8.19470E-14 -8.13140E-14 &
     -8.07510E-14 -8.02580E-14 -7.98340E-14 -7.94800E-14 -7.91940E-14 -7.89770E-14 -7.88280E-14 -7.87460E-14 &
     -7.87300E-14 -7.87800E-14 -7.88950E-14 -7.90740E-14 -7.93170E-14 -7.96230E-14 -7.99910E-14 -8.04200E-14 &
     -8.09110E-14 -8.14610E-14 -8.20710E-14 -8.27390E-14 -8.34640E-14 -8.42470E-14 -8.50850E-14 -8.59780E-14 &
     -8.69250E-14 -8.79250E-14 -8.89770E-14 -9.00810E-14 -9.12350E-14 -9.24380E-14 -9.36890E-14 -9.49870E-14 &
     -9.63320E-14 -9.77230E-14 -9.91580E-14 -1.00640E-13 -1.02160E-13 -1.03720E-13 -1.05320E-13 -1.06960E-13 &
     -1.08640E-13 -1.10350E-13 -1.12100E-13 -1.13890E-13 -1.15700E-13 -1.17550E-13 -1.19430E-13 -1.21350E-13 &
     -1.23290E-13 -1.25250E-13 -1.27250E-13 -1.29270E-13 -1.31320E-13 -1.33380E-13 -1.35470E-13 -1.37590E-13 &
     -1.39720E-13 -1.41870E-13 -1.44040E-13 -1.46220E-13 -1.48420E-13 -1.50640E-13 -1.52860E-13 -1.55100E-13 &
     -1.57350E-13 -1.59610E-13 -1.61880E-13 -1.64150E-13 -1.66440E-13 -1.68720E-13 -1.71020E-13 -1.73310E-13 &
     -1.75610E-13 -1.77900E-13 -1.80200E-13 -1.82500E-13 -1.84790E-13 -1.87080E-13 -1.89360E-13 -1.91640E-13 &
     -1.93910E-13 -1.96170E-13 -1.98430E-13 -2.00670E-13 -2.02910E-13 -2.05130E-13 -2.07340E-13 -2.09540E-13 &
     -2.11720E-13 -2.13880E-13 -2.16030E-13 -2.18160E-13 -2.20280E-13 -2.22370E-13 -2.24440E-13 -2.26490E-13 &
     -2.28520E-13 -2.30530E-13 -2.32510E-13 -2.34470E-13 -2.36400E-13 -2.38310E-13 -2.40190E-13 -2.42040E-13 &
     -2.43870E-13 -2.45660E-13 -2.47430E-13 -2.49160E-13 -2.50860E-13 -2.52540E-13 -2.54180E-13 -2.55780E-13 &
     -2.57360E-13 -2.58890E-13 -2.60400E-13 -2.61870E-13 -2.63300E-13 -2.64700E-13 -2.66060E-13 -2.67380E-13 &
     -2.68660E-13 -2.69910E-13 -2.71120E-13 -2.72290E-13 -2.73420E-13 -2.74510E-13 -2.75560E-13 -2.76570E-13 &
     -2.77540E-13 -2.78470E-13 -2.79360E-13 -2.80210E-13 -2.81010E-13 -2.81770E-13 -2.82500E-13 -2.83180E-13 &
     -2.83810E-13 -2.84410E-13 -2.84960E-13 -2.85470E-13 -2.85930E-13 -2.86360E-13 -2.86740E-13 -2.87070E-13 &
     -2.87370E-13 -2.87620E-13 -2.87830E-13 -2.87990E-13 -2.88120E-13 -2.88190E-13 -2.88230E-13 -2.88230E-13 &
     -2.88180E-13 -2.88090E-13 -2.87950E-13 -2.87770E-13 -2.87560E-13 -2.87300E-13 -2.86990E-13 -2.86650E-13 &
     -2.86260E-13 -2.85840E-13 -2.85370E-13 -2.84860E-13 -2.84310E-13 -2.83730E-13 -2.83100E-13 -2.82440E-13 &
     -2.81730E-13 -2.80990E-13 -2.80210E-13 -2.79390E-13 -2.78540E-13 -2.77650E-13 -2.76720E-13 -2.75760E-13 &
     -2.74760E-13 -2.73730E-13 -2.72660E-13 -2.71560E-13 -2.70430E-13 -2.69270E-13 -2.68070E-13 -2.66840E-13 &
     -2.65580E-13 -2.64290E-13 -2.62970E-13 -2.61620E-13 -2.60240E-13 -2.58830E-13 -2.57400E-13 -2.55930E-13 &
     -2.54450E-13 -2.52930E-13 -2.51400E-13 -2.49830E-13 -2.48250E-13 -2.46640E-13 -2.45010E-13 -2.43350E-13 &
     -2.41680E-13 -2.39990E-13 -2.38270E-13 -2.36540E-13 -2.34790E-13 -2.33020E-13 -2.31230E-13 -2.29430E-13 &
     -2.27610E-13 -2.25770E-13 -2.23930E-13 -2.22060E-13 -2.20190E-13 -2.18300E-13 -2.16400E-13 -2.14500E-13 &
     -2.12580E-13 -2.10650E-13 -2.08710E-13 -2.06760E-13 -2.04810E-13 -2.02850E-13 -2.00880E-13 -1.98910E-13 &
     -1.96940E-13 -1.94960E-13 -1.92980E-13 -1.91000E-13 -1.89010E-13 -1.87030E-13 -1.85040E-13 -1.83050E-13 &
     -1.81070E-13 -1.79090E-13 -1.77110E-13 -1.75130E-13 -1.73160E-13 -1.71190E-13 -1.69220E-13 -1.67270E-13 &
     -1.65310E-13 -1.63370E-13 -1.61430E-13 -1.59500E-13 -1.57580E-13 -1.55670E-13 -1.53770E-13 -1.51880E-13 &
     -1.50010E-13 -1.48140E-13 -1.46280E-13 -1.44440E-13 -1.42620E-13 -1.40800E-13 -1.39000E-13 -1.37220E-13 &
     -1.35450E-13 -1.33700E-13 -1.31960E-13 -1.30240E-13 -1.28540E-13 -1.26860E-13 -1.25190E-13 -1.23550E-13 &
     -1.21920E-13 -1.20310E-13 -1.18730E-13 -1.17160E-13 -1.15620E-13 -1.14090E-13 -1.12590E-13 -1.11110E-13 &
     -1.09650E-13 -1.08220E-13 -1.06810E-13 -1.05420E-13 -1.04060E-13 -1.02720E-13 -1.01400E-13 -1.00110E-13 &
     -9.88420E-14 -9.76010E-14 -9.63850E-14 -9.51940E-14 -9.40290E-14 -9.28910E-14 -9.17790E-14 -9.06930E-14 &
     -8.96330E-14 -8.86010E-14 -8.75950E-14 -8.66160E-14 -8.56640E-14 -8.47400E-14 -8.38430E-14 -8.29730E-14 &
     -8.21310E-14 -8.13160E-14 -8.05290E-14 -7.97700E-14 -7.90380E-14 -7.83340E-14 -7.76580E-14 -7.70100E-14 &
     -7.63890E-14 -7.57960E-14 -7.52310E-14 -7.46940E-14 -7.41850E-14 -7.37030E-14 -7.32480E-14 -7.28210E-14 &
     -7.24210E-14 -7.20490E-14 -7.17040E-14 -7.13850E-14 -7.10940E-14 -7.08290E-14 -7.05910E-14 -7.03790E-14 &
     -7.01930E-14 -7.00340E-14 -6.99000E-14 -6.97920E-14 -6.97100E-14 -6.96520E-14 -6.96200E-14 -6.96120E-14 &
     -6.96280E-14 -6.96690E-14 -6.97340E-14 -6.98220E-14 -6.99340E-14 -7.00680E-14 -7.02260E-14 -7.04060E-14 &
     -7.06080E-14 -7.08320E-14 -7.10780E-14 -7.13450E-14 -7.16320E-14 -7.19410E-14 -7.22690E-14 -7.26180E-14 &
     -7.29850E-14 -7.33720E-14 -7.37780E-14 -7.42020E-14 -7.46430E-14 -7.51020E-14 -7.55790E-14 -7.60720E-14 &
     -7.65810E-14 -7.71060E-14 -7.76470E-14 -7.82030E-14 -7.87740E-14 -7.93590E-14 -7.99580E-14 -8.05700E-14 &
     -8.11960E-14 -8.18340E-14 -8.24840E-14 -8.31460E-14 -8.38200E-14 -8.45040E-14 -8.51990E-14 -8.59040E-14 &
     -8.66180E-14 -8.73410E-14 -8.80730E-14 -8.88140E-14 -8.95620E-14 -9.03180E-14 -9.10800E-14 -9.18490E-14 &
     -9.26240E-14 -9.34040E-14 -9.41890E-14 -9.49790E-14 -9.57740E-14 -9.65710E-14 -9.73730E-14 -9.81770E-14 &
     -9.89830E-14 -9.97910E-14 -1.00600E-13 -1.01410E-13 -1.02220E-13 -1.03040E-13 -1.03850E-13 -1.04660E-13 &
     -1.05470E-13 -1.06280E-13 -1.07080E-13 -1.07890E-13 -1.08690E-13 -1.09480E-13 -1.10280E-13 -1.11070E-13 &
     -1.11850E-13 -1.12640E-13 -1.13410E-13 -1.14180E-13 -1.14940E-13 -1.15700E-13 -1.16450E-13 -1.17200E-13 &
     -1.17930E-13 -1.18660E-13 -1.19380E-13 -1.20090E-13 -1.20790E-13 -1.21490E-13 -1.22170E-13 -1.22840E-13 &
     -1.23510E-13 -1.24160E-13 -1.24800E-13 -1.25430E-13 -1.26050E-13 -1.26660E-13 -1.27260E-13 -1.27840E-13 &
     -1.28410E-13 -1.28970E-13 -1.29520E-13 -1.30050E-13 -1.30570E-13 -1.31070E-13 -1.31560E-13 -1.32040E-13 &
     -1.32500E-13 -1.32950E-13 -1.33380E-13 -1.33800E-13 -1.34210E-13 -1.34590E-13 -1.34960E-13 -1.35320E-13 &
     -1.35660E-13 -1.35990E-13 -1.36290E-13 -1.36590E-13 -1.36860E-13 -1.37120E-13 -1.37360E-13 -1.37590E-13 &
     -1.37800E-13 -1.37990E-13 -1.38170E-13 -1.38320E-13 -1.38460E-13 -1.38590E-13 -1.38690E-13 -1.38780E-13 &
     -1.38850E-13 -1.38910E-13 -1.38940E-13 -1.38960E-13 -1.38960E-13 -1.38940E-13 -1.38910E-13 -1.38860E-13 &
     -1.38790E-13 -1.38700E-13 -1.38590E-13 -1.38470E-13 -1.38330E-13 -1.38170E-13 -1.38000E-13 -1.37810E-13 &
     -1.37600E-13 -1.37370E-13 -1.37130E-13 -1.36870E-13 -1.36590E-13 -1.36300E-13 -1.35990E-13 -1.35660E-13 &
     -1.35310E-13 -1.34960E-13 -1.34580E-13 -1.34190E-13 -1.33780E-13 -1.33360E-13 -1.32920E-13 -1.32470E-13 &
     -1.32000E-13 -1.31510E-13 -1.31010E-13 -1.30500E-13 -1.29970E-13 -1.29430E-13 -1.28870E-13 -1.28300E-13 &
     -1.27720E-13 -1.27120E-13 -1.26510E-13 -1.25890E-13 -1.25250E-13 -1.24600E-13 -1.23940E-13 -1.23270E-13 &
     -1.22580E-13 -1.21880E-13 -1.21170E-13 -1.20450E-13 -1.19720E-13 -1.18980E-13 -1.18230E-13 -1.17470E-13 &
     -1.16690E-13 -1.15910E-13 -1.15120E-13 -1.14320E-13 -1.13510E-13 -1.12690E-13 -1.11870E-13 -1.11030E-13 &
     -1.10190E-13 -1.09340E-13 -1.08480E-13 -1.07620E-13 -1.06740E-13 -1.05870E-13 -1.04980E-13 -1.04090E-13 &
     -1.03200E-13 -1.02300E-13 -1.01390E-13 -1.00480E-13 -9.95670E-14 -9.86490E-14 -9.77270E-14 -9.68020E-14 &
     -9.58730E-14 -9.49410E-14 -9.40070E-14 -9.30710E-14 -9.21330E-14 -9.11920E-14 -9.02510E-14 -8.93080E-14 &
     -8.83650E-14 -8.74200E-14 -8.64760E-14 -8.55320E-14 -8.45880E-14 -8.36440E-14 -8.27010E-14 -8.17590E-14 &
     -8.08190E-14 -7.98800E-14 -7.89430E-14 -7.80090E-14 -7.70770E-14 -7.61470E-14 -7.52210E-14 -7.42980E-14 &
     -7.33790E-14 -7.24630E-14 -7.15510E-14 -7.06440E-14 -6.97410E-14 -6.88430E-14 -6.79500E-14 -6.70620E-14 &
     -6.61800E-14 -6.53030E-14 -6.44320E-14 -6.35670E-14 -6.27090E-14 -6.18570E-14 -6.10120E-14 -6.01740E-14 &
     -5.93430E-14 -5.85190E-14 -5.77030E-14 -5.68950E-14 -5.60950E-14 -5.53020E-14 -5.45190E-14 -5.37430E-14 &
     -5.29770E-14 -5.22190E-14 -5.14700E-14 -5.07310E-14 -5.00010E-14 -4.92800E-14 -4.85690E-14 -4.78680E-14 &
     -4.71770E-14 -4.64960E-14 -4.58250E-14 -4.51640E-14 -4.45140E-14 -4.38740E-14 -4.32450E-14 -4.26270E-14 &
     -4.20200E-14 -4.14240E-14 -4.08390E-14 -4.02650E-14 -3.97030E-14 -3.91520E-14 -3.86120E-14 -3.80840E-14 &
     -3.75680E-14 -3.70630E-14 -3.65700E-14 -3.60890E-14 -3.56200E-14 -3.51620E-14 -3.47170E-14 -3.42830E-14 &
     -3.38620E-14 -3.34520E-14 -3.30550E-14 -3.26700E-14 -3.22960E-14 -3.19350E-14 -3.15860E-14 -3.12490E-14 &
     -3.09240E-14 -3.06110E-14 -3.03100E-14 -3.00210E-14 -2.97440E-14 -2.94790E-14 -2.92260E-14 -2.89850E-14 &
     -2.87550E-14 -2.85380E-14 -2.83320E-14 -2.81380E-14 -2.79550E-14 -2.77840E-14 -2.76240E-14 -2.74760E-14 &
     -2.73400E-14 -2.72140E-14 -2.71000E-14 -2.69970E-14 -2.69050E-14 -2.68240E-14 -2.67540E-14 -2.66950E-14 &
     -2.66460E-14 -2.66080E-14 -2.65800E-14 -2.65620E-14 -2.65550E-14 -2.65580E-14 -2.65710E-14 -2.65930E-14 &
     -2.66260E-14 -2.66680E-14 -2.67200E-14 -2.67810E-14 -2.68510E-14 -2.69310E-14 -2.70190E-14 -2.71160E-14 &
     -2.72220E-14 -2.73370E-14 -2.74600E-14 -2.75910E-14 -2.77300E-14 -2.78780E-14 -2.80330E-14 -2.81950E-14 &
     -2.83660E-14 -2.85430E-14 -2.87280E-14 -2.89200E-14 -2.91180E-14 -2.93230E-14 -2.95340E-14 -2.97520E-14 &
     -2.99760E-14 -3.02050E-14 -3.04410E-14 -3.06820E-14 -3.09280E-14 -3.11800E-14 -3.14370E-14 -3.16980E-14 &
     -3.19640E-14 -3.22350E-14 -3.25100E-14 -3.27890E-14 -3.30720E-14 -3.33590E-14 -3.36490E-14 -3.39430E-14 &
     -3.42400E-14 -3.45400E-14 -3.48430E-14 -3.51480E-14 -3.54560E-14 -3.57660E-14 -3.60790E-14 -3.63930E-14 &
     -3.67090E-14 -3.70270E-14 -3.73460E-14 -3.76670E-14 -3.79880E-14 -3.83100E-14 -3.86330E-14 -3.89570E-14 &
     -3.92800E-14 -3.96050E-14 -3.99290E-14 -4.02530E-14 -4.05760E-14 -4.08990E-14 -4.12220E-14 -4.15440E-14 &
     -4.18640E-14 -4.21840E-14 -4.25020E-14 -4.28190E-14 -4.31340E-14 -4.34480E-14 -4.37600E-14 -4.40690E-14 &
     -4.43770E-14 -4.46820E-14 -4.49840E-14 -4.52840E-14 -4.55820E-14 -4.58760E-14 -4.61670E-14 -4.64550E-14 &
     -4.67400E-14 -4.70210E-14 -4.72990E-14 -4.75730E-14 -4.78440E-14 -4.81100E-14 -4.83720E-14 -4.86300E-14 &
     -4.88840E-14 -4.91340E-14 -4.93790E-14 -4.96190E-14 -4.98550E-14 -5.00860E-14 -5.03110E-14 -5.05320E-14 &
     -5.07480E-14 -5.09590E-14 -5.11640E-14 -5.13640E-14 -5.15590E-14 -5.17480E-14 -5.19310E-14 -5.21090E-14 &
     -5.22810E-14 -5.24470E-14 -5.26070E-14 -5.27610E-14 -5.29100E-14 -5.30520E-14 -5.31880E-14 -5.33180E-14 &
     -5.34420E-14 -5.35590E-14 -5.36700E-14 -5.37750E-14 -5.38740E-14 -5.39660E-14 -5.40520E-14 -5.41310E-14 &
     -5.42040E-14 -5.42700E-14 -5.43300E-14 -5.43830E-14 -5.44300E-14 -5.44700E-14 -5.45040E-14 -5.45310E-14 &
     -5.45520E-14 -5.45660E-14 -5.45740E-14 -5.45750E-14 -5.45690E-14 -5.45570E-14 -5.45380E-14 -5.45130E-14 &
     -5.44820E-14 -5.44430E-14 -5.43990E-14 -5.43480E-14 -5.42910E-14 -5.42280E-14 -5.41580E-14 -5.40820E-14 &
     -5.40000E-14 -5.39110E-14 -5.38170E-14 -5.37160E-14 -5.36090E-14 -5.34960E-14 -5.33780E-14 -5.32530E-14 &
     -5.31230E-14 -5.29860E-14 -5.28440E-14 -5.26970E-14 -5.25430E-14 -5.23840E-14 -5.22200E-14 -5.20500E-14 &
     -5.18750E-14 -5.16940E-14 -5.15090E-14 -5.13180E-14 -5.11220E-14 -5.09210E-14 -5.07160E-14 -5.05050E-14 &
     -5.02900E-14 -5.00690E-14 -4.98450E-14 -4.96160E-14 -4.93820E-14 -4.91450E-14 -4.89030E-14 -4.86570E-14 &
     -4.84070E-14 -4.81530E-14 -4.78950E-14 -4.76340E-14 -4.73690E-14 -4.71010E-14 -4.68290E-14 -4.65540E-14 &
     -4.62760E-14 -4.59940E-14 -4.57100E-14 -4.54230E-14 -4.51330E-14 -4.48400E-14 -4.45450E-14 -4.42470E-14 &
     -4.39470E-14 -4.36440E-14 -4.33390E-14 -4.30330E-14 -4.27240E-14 -4.24140E-14 -4.21010E-14 -4.17870E-14 &
     -4.14720E-14 -4.11550E-14 -4.08370E-14 -4.05170E-14 -4.01970E-14 -3.98760E-14 -3.95530E-14 -3.92300E-14 &
     -3.89070E-14 -3.85820E-14 -3.82570E-14 -3.79320E-14 -3.76070E-14 -3.72810E-14 -3.69560E-14 -3.66300E-14 &
     -3.63040E-14 -3.59790E-14 -3.56540E-14 -3.53300E-14 -3.50060E-14 -3.46830E-14 -3.43600E-14 -3.40380E-14 &
     -3.37180E-14 -3.33980E-14 -3.30790E-14 -3.27610E-14 -3.24450E-14 -3.21300E-14 -3.18170E-14 -3.15050E-14 &
     -3.11950E-14 -3.08870E-14 -3.05810E-14 -3.02760E-14 -2.99740E-14 -2.96730E-14 -2.93750E-14 -2.90790E-14 &
     -2.87850E-14 -2.84940E-14 -2.82050E-14 -2.79190E-14 -2.76350E-14 -2.73540E-14 -2.70760E-14 -2.68000E-14 &
     -2.65280E-14 -2.62580E-14 -2.59920E-14 -2.57280E-14 -2.54680E-14 -2.52100E-14 -2.49560E-14 -2.47060E-14 &
     -2.44580E-14 -2.42140E-14 -2.39740E-14 -2.37370E-14 -2.35030E-14 -2.32740E-14 -2.30470E-14 -2.28250E-14 &
     -2.26060E-14 -2.23910E-14 -2.21800E-14 -2.19730E-14 -2.17690E-14 -2.15700E-14 -2.13740E-14 -2.11820E-14 &
     -2.09940E-14 -2.08110E-14 -2.06310E-14 -2.04550E-14 -2.02840E-14 -2.01160E-14 -1.99530E-14 -1.97930E-14 &
     -1.96380E-14 -1.94870E-14 -1.93400E-14 -1.91970E-14 -1.90590E-14 -1.89240E-14 -1.87940E-14 -1.86680E-14 &
     -1.85460E-14 -1.84280E-14 -1.83140E-14 -1.82050E-14 -1.80990E-14 -1.79980E-14 -1.79010E-14 -1.78080E-14 &
     -1.77190E-14 -1.76340E-14 -1.75530E-14 -1.74770E-14 -1.74040E-14 -1.73350E-14 -1.72710E-14 -1.72100E-14 &
     -1.71530E-14 -1.71000E-14 -1.70510E-14 -1.70060E-14 -1.69640E-14 -1.69270E-14 -1.68930E-14 -1.68630E-14 &
     -1.68370E-14 -1.68140E-14 -1.67950E-14 -1.67790E-14 -1.67670E-14 -1.67590E-14 -1.67540E-14 -1.67520E-14 &
     -1.67540E-14 -1.67590E-14 -1.67680E-14 -1.67790E-14 -1.67940E-14 -1.68120E-14 -1.68330E-14 -1.68570E-14 &
     -1.68840E-14 -1.69140E-14 -1.69470E-14 -1.69830E-14 -1.70210E-14 -1.70620E-14 -1.71060E-14 -1.71520E-14 &
     -1.72010E-14 -1.72520E-14 -1.73060E-14 -1.73620E-14 -1.74210E-14 -1.74810E-14 -1.75440E-14 -1.76090E-14 &
     -1.76760E-14 -1.77450E-14 -1.78150E-14 -1.78880E-14 -1.79620E-14 -1.80390E-14 -1.81160E-14 -1.81960E-14 &
     -1.82770E-14 -1.83590E-14 -1.84430E-14 -1.85280E-14 -1.86140E-14 -1.87010E-14 -1.87900E-14 -1.88790E-14 &
     -1.89700E-14 -1.90620E-14 -1.91540E-14 -1.92470E-14 -1.93410E-14 -1.94360E-14 -1.95310E-14 -1.96270E-14 &
     -1.97230E-14 -1.98200E-14 -1.99170E-14 -2.00140E-14 -2.01120E-14 -2.02090E-14 -2.03070E-14 -2.04050E-14 &
     -2.05030E-14 -2.06010E-14 -2.06980E-14 -2.07960E-14 -2.08930E-14 -2.09890E-14 -2.10860E-14 -2.11820E-14 &
     -2.12770E-14 -2.13720E-14 -2.14660E-14 -2.15600E-14 -2.16530E-14 -2.17450E-14 -2.18360E-14 -2.19260E-14 &
     -2.20160E-14 -2.21040E-14 -2.21910E-14 -2.22770E-14 -2.23620E-14 -2.24460E-14 -2.25290E-14 -2.26100E-14 &
     -2.26900E-14 -2.27690E-14 -2.28460E-14 -2.29210E-14 -2.29950E-14 -2.30670E-14 -2.31380E-14 -2.32070E-14 &
     -2.32750E-14 -2.33400E-14 -2.34040E-14 -2.34660E-14 -2.35260E-14 -2.35840E-14 -2.36410E-14 -2.36950E-14 &
     -2.37480E-14 -2.37980E-14 -2.38460E-14 -2.38920E-14 -2.39360E-14 -2.39780E-14 -2.40180E-14 -2.40560E-14 &
     -2.40910E-14 -2.41240E-14 -2.41550E-14 -2.41840E-14 -2.42100E-14 -2.42340E-14 -2.42560E-14 -2.42750E-14 &
     -2.42920E-14 -2.43060E-14 -2.43190E-14 -2.43280E-14 -2.43360E-14 -2.43410E-14 -2.43430E-14 -2.43430E-14 &
     -2.43410E-14 -2.43360E-14 -2.43290E-14 -2.43190E-14 -2.43070E-14 -2.42920E-14 -2.42750E-14 -2.42560E-14 &
     -2.42340E-14 -2.42090E-14 -2.41830E-14 -2.41530E-14 -2.41220E-14 -2.40870E-14 -2.40510E-14 -2.40120E-14 &
     -2.39700E-14 -2.39260E-14 -2.38800E-14 -2.38310E-14 -2.37800E-14 -2.37260E-14 -2.36700E-14 -2.36120E-14 &
     -2.35510E-14 -2.34880E-14 -2.34230E-14 -2.33550E-14 -2.32860E-14 -2.32140E-14 -2.31390E-14 -2.30630E-14 &
     -2.29840E-14 -2.29030E-14 -2.28210E-14 -2.27360E-14 -2.26480E-14 -2.25590E-14 -2.24680E-14 -2.23750E-14 &
     -2.22800E-14 -2.21830E-14 -2.20840E-14 -2.19830E-14 -2.18810E-14 -2.17760E-14 -2.16700E-14 -2.15620E-14 &
     -2.14520E-14 -2.13410E-14 -2.12280E-14 -2.11130E-14 -2.09970E-14 -2.08800E-14 -2.07610E-14 -2.06400E-14 &
     -2.05180E-14 -2.03940E-14 -2.02700E-14 -2.01440E-14 -2.00160E-14 -1.98870E-14 -1.97580E-14 -1.96270E-14 &
     -1.94950E-14 -1.93610E-14 -1.92270E-14 -1.90920E-14 -1.89560E-14 -1.88190E-14 -1.86810E-14 -1.85420E-14 &
     -1.84020E-14 -1.82620E-14 -1.81210E-14 -1.79790E-14 -1.78360E-14 -1.76930E-14 -1.75490E-14 -1.74050E-14 &
     -1.72610E-14 -1.71150E-14 -1.69700E-14 -1.68240E-14 -1.66780E-14 -1.65310E-14 -1.63850E-14 -1.62380E-14 &
     -1.60910E-14 -1.59430E-14 -1.57960E-14 -1.56490E-14 -1.55020E-14 -1.53550E-14 -1.52070E-14 -1.50600E-14 &
     -1.49140E-14 -1.47670E-14 -1.46210E-14 -1.44750E-14 -1.43290E-14 -1.41840E-14 -1.40390E-14 -1.38940E-14 &
     -1.37500E-14 -1.36070E-14 -1.34640E-14 -1.33210E-14 -1.31790E-14 -1.30380E-14 -1.28980E-14 -1.27580E-14 &
     -1.26190E-14 -1.24810E-14 -1.23440E-14 -1.22070E-14 -1.20720E-14 -1.19370E-14 -1.18030E-14 -1.16700E-14 &
     -1.15390E-14 -1.14080E-14 -1.12780E-14 -1.11500E-14 -1.10230E-14 -1.08960E-14 -1.07710E-14 -1.06480E-14 &
     -1.05250E-14 -1.04040E-14 -1.02840E-14 -1.01650E-14 -1.00480E-14 -9.93180E-15 -9.81720E-15 -9.70390E-15 &
     -9.59210E-15 -9.48170E-15 -9.37280E-15 -9.26530E-15 -9.15930E-15 -9.05490E-15 -8.95190E-15 -8.85060E-15 &
     -8.75070E-15 -8.65240E-15 -8.55570E-15 -8.46060E-15 -8.36720E-15 -8.27530E-15 -8.18510E-15 -8.09650E-15 &
     -8.00950E-15 -7.92420E-15 -7.84060E-15 -7.75870E-15 -7.67850E-15 -7.60000E-15 -7.52320E-15 -7.44810E-15 &
     -7.37470E-15 -7.30300E-15 -7.23310E-15 -7.16500E-15 -7.09850E-15 -7.03380E-15 -6.97090E-15 -6.90970E-15 &
     -6.85030E-15 -6.79250E-15 -6.73660E-15 -6.68240E-15 -6.62990E-15 -6.57920E-15 -6.53030E-15 -6.48310E-15 &
     -6.43760E-15 -6.39390E-15 -6.35190E-15 -6.31160E-15 -6.27310E-15 -6.23620E-15 -6.20100E-15 -6.16750E-15 &
     -6.13560E-15 -6.10540E-15 -6.07690E-15 -6.05000E-15 -6.02470E-15 -6.00100E-15 -5.97900E-15 -5.95850E-15 &
     -5.93960E-15 -5.92230E-15 -5.90650E-15 -5.89220E-15 -5.87950E-15 -5.86820E-15 -5.85850E-15 -5.85020E-15 &
     -5.84340E-15 -5.83800E-15 -5.83400E-15 -5.83140E-15 -5.83020E-15 -5.83040E-15 -5.83190E-15 -5.83470E-15 &
     -5.83890E-15 -5.84430E-15 -5.85100E-15 -5.85890E-15 -5.86810E-15 -5.87850E-15 -5.89000E-15 -5.90270E-15 &
     -5.91660E-15 -5.93160E-15 -5.94770E-15 -5.96480E-15 -5.98300E-15 -6.00230E-15 -6.02250E-15 -6.04370E-15 &
     -6.06590E-15 -6.08900E-15 -6.11310E-15 -6.13800E-15 -6.16380E-15 -6.19040E-15 -6.21780E-15 -6.24610E-15 &
     -6.27510E-15 -6.30480E-15 -6.33530E-15 -6.36650E-15 -6.39830E-15 -6.43070E-15 -6.46380E-15 -6.49740E-15 &
     -6.53160E-15 -6.56630E-15 -6.60150E-15 -6.63720E-15 -6.67340E-15 -6.71000E-15 -6.74690E-15 -6.78430E-15 &
     -6.82200E-15 -6.86010E-15 -6.89840E-15 -6.93710E-15 -6.97600E-15 -7.01510E-15 -7.05440E-15 -7.09400E-15 &
     -7.13360E-15 -7.17350E-15 -7.21340E-15 -7.25340E-15 -7.29350E-15 -7.33360E-15 -7.37370E-15 -7.41380E-15 &
     -7.45390E-15 -7.49390E-15 -7.53390E-15 -7.57370E-15 -7.61340E-15 -7.65290E-15 -7.69230E-15 -7.73150E-15 &
     -7.77040E-15 -7.80910E-15 -7.84760E-15 -7.88580E-15 -7.92370E-15 -7.96130E-15 -7.99850E-15 -8.03540E-15 &
     -8.07190E-15 -8.10800E-15 -8.14370E-15 -8.17900E-15 -8.21380E-15 -8.24820E-15 -8.28210E-15 -8.31550E-15 &
     -8.34840E-15 -8.38070E-15 -8.41250E-15 -8.44380E-15 -8.47440E-15 -8.50450E-15 -8.53390E-15 -8.56270E-15 &
     -8.59090E-15 -8.61840E-15 -8.64520E-15 -8.67140E-15 -8.69690E-15 -8.72160E-15 -8.74570E-15 -8.76910E-15 &
     -8.79180E-15 -8.81370E-15 -8.83480E-15 -8.85520E-15 -8.87480E-15 -8.89370E-15 -8.91170E-15 -8.92890E-15 &
     -8.94540E-15 -8.96100E-15 -8.97580E-15 -8.98980E-15 -9.00290E-15 -9.01520E-15 -9.02660E-15 -9.03710E-15 &
     -9.04680E-15 -9.05570E-15 -9.06360E-15 -9.07070E-15 -9.07680E-15 -9.08210E-15 -9.08650E-15 -9.09000E-15 &
     -9.09270E-15 -9.09440E-15 -9.09530E-15 -9.09530E-15 -9.09440E-15 -9.09260E-15 -9.09000E-15 -9.08640E-15 &
     -9.08200E-15 -9.07670E-15 -9.07050E-15 -9.06340E-15 -9.05540E-15 -9.04660E-15 -9.03690E-15 -9.02640E-15 &
     -9.01500E-15 -9.00270E-15 -8.98960E-15 -8.97560E-15 -8.96080E-15 -8.94510E-15 -8.92870E-15 -8.91140E-15 &
     -8.89330E-15 -8.87440E-15 -8.85460E-15 -8.83410E-15 -8.81280E-15 -8.79070E-15 -8.76790E-15 -8.74420E-15 &
     -8.71980E-15 -8.69470E-15 -8.66880E-15 -8.64220E-15 -8.61480E-15 -8.58680E-15 -8.55800E-15 -8.52860E-15 &
     -8.49850E-15 -8.46770E-15 -8.43620E-15 -8.40410E-15 -8.37140E-15 -8.33800E-15 -8.30400E-15 -8.26940E-15 &
     -8.23430E-15 -8.19850E-15 -8.16220E-15 -8.12530E-15 -8.08790E-15 -8.04990E-15 -8.01140E-15 -7.97240E-15 &
     -7.93300E-15 -7.89300E-15 -7.85260E-15 -7.81180E-15 -7.77050E-15 -7.72880E-15 -7.68670E-15 -7.64420E-15 &
     -7.60130E-15 -7.55810E-15 -7.51450E-15 -7.47050E-15 -7.42620E-15 -7.38160E-15 -7.33680E-15 -7.29160E-15 &
     -7.24620E-15 -7.20050E-15 -7.15460E-15 -7.10840E-15 -7.06210E-15 -7.01550E-15 -6.96880E-15 -6.92190E-15 &
     -6.87480E-15 -6.82760E-15 -6.78030E-15 -6.73290E-15 -6.68530E-15 -6.63770E-15 -6.59000E-15 -6.54220E-15 &
     -6.49440E-15 -6.44660E-15 -6.39880E-15 -6.35090E-15 -6.30310E-15 -6.25530E-15 -6.20760E-15 -6.15990E-15 &
     -6.11220E-15 -6.06470E-15 -6.01720E-15 -5.96990E-15 -5.92270E-15 -5.87560E-15 -5.82860E-15 -5.78180E-15 &
     -5.73520E-15 -5.68870E-15 -5.64240E-15 -5.59630E-15 -5.55050E-15 -5.50480E-15 -5.45940E-15 -5.41430E-15 &
     -5.36940E-15 -5.32480E-15 -5.28040E-15 -5.23640E-15 -5.19260E-15 -5.14920E-15 -5.10610E-15 -5.06330E-15 &
     -5.02090E-15 -4.97880E-15 -4.93710E-15 -4.89570E-15 -4.85480E-15 -4.81420E-15 -4.77400E-15 -4.73420E-15 &
     -4.69480E-15 -4.65580E-15 -4.61730E-15 -4.57920E-15 -4.54150E-15 -4.50430E-15 -4.46750E-15 -4.43120E-15 &
     -4.39530E-15 -4.35990E-15 -4.32500E-15 -4.29060E-15 -4.25670E-15 -4.22320E-15 -4.19020E-15 -4.15780E-15 &
     -4.12580E-15 -4.09440E-15 -4.06340E-15 -4.03300E-15 -4.00310E-15 -3.97380E-15 -3.94490E-15 -3.91670E-15 &
     -3.88890E-15 -3.86170E-15 -3.83500E-15 -3.80890E-15 -3.78330E-15 -3.75830E-15 -3.73380E-15 -3.70990E-15 &
     -3.68650E-15 -3.66370E-15 -3.64140E-15 -3.61970E-15 -3.59850E-15 -3.57790E-15 -3.55790E-15 -3.53830E-15 &
     -3.51940E-15 -3.50090E-15 -3.48300E-15 -3.46570E-15 -3.44890E-15 -3.43260E-15 -3.41690E-15 -3.40170E-15 &
     -3.38700E-15 -3.37290E-15 -3.35930E-15 -3.34620E-15 -3.33370E-15 -3.32170E-15 -3.31020E-15 -3.29920E-15 &
     -3.28870E-15 -3.27870E-15 -3.26920E-15 -3.26020E-15 -3.25170E-15 -3.24370E-15 -3.23620E-15 -3.22910E-15 &
     -3.22250E-15 -3.21640E-15 -3.21070E-15 -3.20550E-15 -3.20080E-15 -3.19650E-15 -3.19260E-15 -3.18910E-15 &
     -3.18610E-15 -3.18350E-15 -3.18130E-15 -3.17960E-15 -3.17820E-15 -3.17730E-15 -3.17670E-15 -3.17660E-15 &
     -3.17680E-15 -3.17740E-15 -3.17830E-15 -3.17960E-15 -3.18120E-15 -3.18310E-15 -3.18530E-15 -3.18790E-15 &
     -3.19080E-15 -3.19400E-15 -3.19750E-15 -3.20120E-15 -3.20530E-15 -3.20960E-15 -3.21420E-15 -3.21900E-15 &
     -3.22410E-15 -3.22950E-15 -3.23500E-15 -3.24080E-15 -3.24680E-15 -3.25300E-15 -3.25940E-15 -3.26610E-15 &
     -3.27290E-15 -3.27990E-15 -3.28700E-15 -3.29430E-15 -3.30180E-15 -3.30950E-15 -3.31720E-15 -3.32510E-15 &
     -3.33310E-15 -3.34110E-15 -3.34930E-15 -3.35760E-15 -3.36590E-15 -3.37430E-15 -3.38280E-15 -3.39130E-15 &
     -3.39990E-15 -3.40850E-15 -3.41710E-15 -3.42580E-15 -3.43440E-15 -3.44310E-15 -3.45180E-15 -3.46050E-15 &
     -3.46910E-15 -3.47770E-15 -3.48630E-15 -3.49490E-15 -3.50330E-15 -3.51170E-15 -3.52010E-15 -3.52830E-15 &
     -3.53650E-15 -3.54460E-15 -3.55260E-15 -3.56050E-15 -3.56840E-15 -3.57610E-15 -3.58370E-15 -3.59120E-15 &
     -3.59860E-15 -3.60590E-15 -3.61300E-15 -3.61990E-15 -3.62680E-15 -3.63340E-15 -3.63990E-15 -3.64620E-15 &
     -3.65230E-15 -3.65830E-15 -3.66410E-15 -3.66960E-15 -3.67500E-15 -3.68020E-15 -3.68520E-15 -3.69000E-15 &
     -3.69450E-15 -3.69880E-15 -3.70290E-15 -3.70680E-15 -3.71040E-15 -3.71370E-15 -3.71680E-15 -3.71970E-15 &
     -3.72230E-15 -3.72460E-15 -3.72670E-15 -3.72850E-15 -3.73010E-15 -3.73140E-15 -3.73240E-15 -3.73310E-15 &
     -3.73360E-15 -3.73370E-15 -3.73360E-15 -3.73320E-15 -3.73250E-15 -3.73160E-15 -3.73030E-15 -3.72870E-15 &
     -3.72690E-15 -3.72470E-15 -3.72230E-15 -3.71950E-15 -3.71650E-15 -3.71320E-15 -3.70950E-15 -3.70560E-15 &
     -3.70140E-15 -3.69680E-15 -3.69200E-15 -3.68690E-15 -3.68150E-15 -3.67580E-15 -3.66980E-15 -3.66360E-15 &
     -3.65700E-15 -3.65010E-15 -3.64290E-15 -3.63540E-15 -3.62760E-15 -3.61960E-15 -3.61120E-15 -3.60250E-15 &
     -3.59350E-15 -3.58430E-15 -3.57470E-15 -3.56490E-15 -3.55480E-15 -3.54440E-15 -3.53370E-15 -3.52280E-15 &
     -3.51150E-15 -3.50000E-15 -3.48830E-15 -3.47620E-15 -3.46390E-15 -3.45140E-15 -3.43850E-15 -3.42550E-15 &
     -3.41210E-15 -3.39850E-15 -3.38470E-15 -3.37060E-15 -3.35630E-15 -3.34180E-15 -3.32700E-15 -3.31200E-15 &
     -3.29680E-15 -3.28130E-15 -3.26570E-15 -3.24980E-15 -3.23370E-15 -3.21740E-15 -3.20100E-15 -3.18430E-15 &
     -3.16740E-15 -3.15030E-15 -3.13310E-15 -3.11570E-15 -3.09810E-15 -3.08030E-15 -3.06240E-15 -3.04430E-15 &
     -3.02610E-15 -3.00770E-15 -2.98920E-15 -2.97050E-15 -2.95160E-15 -2.93260E-15 -2.91350E-15 -2.89430E-15 &
     -2.87490E-15 -2.85550E-15 -2.83590E-15 -2.81630E-15 -2.79660E-15 -2.77680E-15 -2.75690E-15 -2.73690E-15 &
     -2.71690E-15 -2.69680E-15 -2.67660E-15 -2.65640E-15 -2.63620E-15 -2.61580E-15 -2.59550E-15 -2.57510E-15 &
     -2.55470E-15 -2.53420E-15 -2.51370E-15 -2.49320E-15 -2.47280E-15 -2.45220E-15 -2.43170E-15 -2.41120E-15 &
     -2.39080E-15 -2.37030E-15 -2.34980E-15 -2.32940E-15 -2.30900E-15 -2.28860E-15 -2.26820E-15 -2.24790E-15 &
     -2.22770E-15 -2.20750E-15 -2.18730E-15 -2.16720E-15 -2.14720E-15 -2.12720E-15 -2.10730E-15 -2.08750E-15 &
     -2.06780E-15 -2.04820E-15 -2.02860E-15 -2.00920E-15 -1.98990E-15 -1.97060E-15 -1.95150E-15 -1.93250E-15 &
     -1.91360E-15 -1.89480E-15 -1.87610E-15 -1.85760E-15 -1.83920E-15 -1.82100E-15 -1.80290E-15 -1.78490E-15 &
     -1.76710E-15 -1.74950E-15 -1.73200E-15 -1.71470E-15 -1.69760E-15 -1.68060E-15 -1.66380E-15 -1.64720E-15 &
     -1.63070E-15 -1.61450E-15 -1.59840E-15 -1.58240E-15 -1.56670E-15 -1.55120E-15 -1.53580E-15 -1.52060E-15 &
     -1.50570E-15 -1.49090E-15 -1.47630E-15 -1.46200E-15 -1.44780E-15 -1.43390E-15 -1.42010E-15 -1.40660E-15 &
     -1.39320E-15 -1.38010E-15 -1.36720E-15 -1.35450E-15 -1.34210E-15 -1.32980E-15 -1.31780E-15 -1.30600E-15 &
     -1.29440E-15 -1.28300E-15 -1.27180E-15 -1.26090E-15 -1.25020E-15 -1.23970E-15 -1.22950E-15 -1.21950E-15 &
     -1.20970E-15 -1.20010E-15 -1.19080E-15 -1.18170E-15 -1.17280E-15 -1.16420E-15 -1.15580E-15 -1.14760E-15 &
     -1.13960E-15 -1.13190E-15 -1.12440E-15 -1.11710E-15 -1.11010E-15 -1.10330E-15 -1.09670E-15 -1.09030E-15 &
     -1.08420E-15 -1.07820E-15 -1.07250E-15 -1.06700E-15 -1.06180E-15 -1.05670E-15 -1.05190E-15 -1.04730E-15 &
     -1.04290E-15 -1.03870E-15 -1.03470E-15 -1.03100E-15 -1.02740E-15 -1.02410E-15 -1.02090E-15 -1.01800E-15 &
     -1.01530E-15 -1.01280E-15 -1.01040E-15 -1.00830E-15 -1.00640E-15 -1.00460E-15 -1.00310E-15 -1.00170E-15 &
     -1.00050E-15 -9.99490E-16 -9.98670E-16 -9.98020E-16 -9.97550E-16 -9.97260E-16 -9.97140E-16 -9.97190E-16 &
     -9.97410E-16 -9.97790E-16 -9.98330E-16 -9.99030E-16 -9.99880E-16 -1.00090E-15 -1.00200E-15 -1.00330E-15 &
     -1.00470E-15 -1.00630E-15 -1.00800E-15 -1.00980E-15 -1.01180E-15 -1.01390E-15 -1.01610E-15 -1.01840E-15 &
     -1.02090E-15 -1.02350E-15 -1.02610E-15 -1.02890E-15 -1.03180E-15 -1.03480E-15 -1.03790E-15 -1.04110E-15 &
     -1.04440E-15 -1.04780E-15 -1.05130E-15 -1.05480E-15 -1.05840E-15 -1.06220E-15 -1.06590E-15 -1.06980E-15 &
     -1.07370E-15 -1.07770E-15 -1.08180E-15 -1.08590E-15 -1.09010E-15 -1.09440E-15 -1.09860E-15 -1.10300E-15 &
     -1.10730E-15 -1.11180E-15 -1.11620E-15 -1.12070E-15 -1.12520E-15 -1.12970E-15 -1.13430E-15 -1.13880E-15 &
     -1.14340E-15 -1.14800E-15 -1.15260E-15 -1.15730E-15 -1.16190E-15 -1.16650E-15 -1.17120E-15 -1.17580E-15 &
     -1.18040E-15 -1.18500E-15 -1.18960E-15 -1.19420E-15 -1.19870E-15 -1.20330E-15 -1.20780E-15 -1.21220E-15 &
     -1.21670E-15 -1.22110E-15 -1.22550E-15 -1.22980E-15 -1.23410E-15 -1.23840E-15 -1.24260E-15 -1.24680E-15 &
     -1.25090E-15 -1.25500E-15 -1.25900E-15 -1.26300E-15 -1.26690E-15 -1.27080E-15 -1.27460E-15 -1.27830E-15 &
     -1.28200E-15 -1.28560E-15 -1.28910E-15 -1.29250E-15 -1.29590E-15 -1.29920E-15 -1.30240E-15 -1.30560E-15 &
     -1.30870E-15 -1.31160E-15 -1.31450E-15 -1.31740E-15 -1.32010E-15 -1.32270E-15 -1.32530E-15 -1.32770E-15 &
     -1.33010E-15 -1.33230E-15 -1.33450E-15 -1.33660E-15 -1.33850E-15 -1.34040E-15 -1.34220E-15 -1.34390E-15 &
     -1.34540E-15 -1.34690E-15 -1.34830E-15 -1.34950E-15 -1.35070E-15 -1.35170E-15 -1.35270E-15 -1.35350E-15 &
     -1.35420E-15 -1.35480E-15 -1.35530E-15 -1.35570E-15 -1.35600E-15 -1.35620E-15 -1.35630E-15 -1.35620E-15 &
     -1.35610E-15 -1.35580E-15 -1.35550E-15 -1.35500E-15 -1.35440E-15 -1.35370E-15 -1.35290E-15 -1.35200E-15 &
     -1.35090E-15 -1.34980E-15 -1.34850E-15 -1.34710E-15 -1.34570E-15 -1.34410E-15 -1.34240E-15 -1.34060E-15 &
     -1.33860E-15 -1.33660E-15 -1.33450E-15 -1.33220E-15 -1.32990E-15 -1.32750E-15 -1.32490E-15 -1.32220E-15 &
     -1.31950E-15 -1.31660E-15 -1.31370E-15 -1.31060E-15 -1.30750E-15 -1.30420E-15 -1.30090E-15 -1.29740E-15 &
     -1.29390E-15 -1.29020E-15 -1.28650E-15 -1.28270E-15 -1.27880E-15 -1.27480E-15 -1.27070E-15 -1.26650E-15 &
     -1.26230E-15 -1.25800E-15 -1.25350E-15 -1.24900E-15 -1.24450E-15 -1.23980E-15 -1.23510E-15 -1.23020E-15 &
     -1.22540E-15 -1.22040E-15 -1.21540E-15 -1.21030E-15 -1.20510E-15 -1.19990E-15 -1.19460E-15 -1.18920E-15 &
     -1.18380E-15 -1.17830E-15 -1.17270E-15 -1.16710E-15 -1.16150E-15 -1.15570E-15 -1.15000E-15 -1.14410E-15 &
     -1.13830E-15 -1.13240E-15 -1.12640E-15 -1.12040E-15 -1.11430E-15 -1.10830E-15 -1.10210E-15 -1.09600E-15 &
     -1.08980E-15 -1.08360E-15 -1.07730E-15 -1.07100E-15 -1.06470E-15 -1.05840E-15 -1.05200E-15 -1.04560E-15 &
     -1.03920E-15 -1.03280E-15 -1.02640E-15 -1.01990E-15 -1.01350E-15 -1.00700E-15 -1.00050E-15 -9.93970E-16 &
     -9.87460E-16 -9.80930E-16 -9.74410E-16 -9.67870E-16 -9.61340E-16 -9.54810E-16 -9.48280E-16 -9.41750E-16 &
     -9.35230E-16 -9.28720E-16 -9.22220E-16 -9.15730E-16 -9.09250E-16 -9.02790E-16 -8.96340E-16 -8.89900E-16 &
     -8.83490E-16 -8.77090E-16 -8.70720E-16 -8.64360E-16 -8.58030E-16 -8.51730E-16 -8.45450E-16 -8.39200E-16 &
     -8.32980E-16 -8.26790E-16 -8.20630E-16 -8.14510E-16 -8.08420E-16 -8.02370E-16 -7.96360E-16 -7.90380E-16 &
     -7.84450E-16 -7.78550E-16 -7.72700E-16 -7.66880E-16 -7.61110E-16 -7.55390E-16 -7.49710E-16 -7.44070E-16 &
     -7.38480E-16 -7.32930E-16 -7.27440E-16 -7.21990E-16 -7.16590E-16 -7.11240E-16 -7.05940E-16 -7.00700E-16 &
     -6.95500E-16 -6.90360E-16 -6.85270E-16 -6.80240E-16 -6.75260E-16 -6.70350E-16 -6.65480E-16 -6.60680E-16 &
     -6.55930E-16 -6.51250E-16 -6.46620E-16 -6.42060E-16 -6.37560E-16 -6.33120E-16 -6.28750E-16 -6.24440E-16 &
     -6.20190E-16 -6.16010E-16 -6.11900E-16 -6.07860E-16 -6.03890E-16 -5.99980E-16 -5.96150E-16 -5.92380E-16 &
     -5.88690E-16 -5.85070E-16 -5.81530E-16 -5.78060E-16 -5.74670E-16 -5.71350E-16 -5.68110E-16 -5.64940E-16 &
     -5.61850E-16 -5.58840E-16 -5.55900E-16 -5.53040E-16 -5.50250E-16 -5.47530E-16 -5.44880E-16 -5.42300E-16 &
     -5.39790E-16 -5.37350E-16 -5.34970E-16 -5.32660E-16 -5.30410E-16 -5.28230E-16 -5.26110E-16 -5.24050E-16 &
     -5.22050E-16 -5.20110E-16 -5.18230E-16 -5.16410E-16 -5.14640E-16 -5.12930E-16 -5.11270E-16 -5.09670E-16 &
     -5.08110E-16 -5.06610E-16 -5.05160E-16 -5.03760E-16 -5.02410E-16 -5.01100E-16 -4.99840E-16 -4.98630E-16 &
     -4.97450E-16 -4.96330E-16 -4.95240E-16 -4.94200E-16 -4.93190E-16 -4.92230E-16 -4.91300E-16 -4.90410E-16 &
     -4.89550E-16 -4.88730E-16 -4.87950E-16 -4.87200E-16 -4.86480E-16 -4.85790E-16 -4.85130E-16 -4.84490E-16 &
     -4.83890E-16 -4.83310E-16 -4.82760E-16 -4.82240E-16 -4.81730E-16 -4.81250E-16 -4.80800E-16 -4.80360E-16 &
     -4.79940E-16 -4.79540E-16 -4.79160E-16 -4.78800E-16 -4.78450E-16 -4.78110E-16 -4.77790E-16 -4.77480E-16 &
     -4.77190E-16 -4.76900E-16 -4.76620E-16 -4.76360E-16 -4.76100E-16 -4.75840E-16 -4.75590E-16 -4.75350E-16 &
     -4.75110E-16 -4.74870E-16 -4.74630E-16 -4.74390E-16 -4.74160E-16 -4.73920E-16 -4.73690E-16 -4.73460E-16 &
     -4.73220E-16 -4.72990E-16 -4.72760E-16 -4.72530E-16 -4.72290E-16 -4.72060E-16 -4.71820E-16 -4.71590E-16 &
     -4.71350E-16 -4.71120E-16 -4.70880E-16 -4.70640E-16 -4.70400E-16 -4.70150E-16 -4.69910E-16 -4.69660E-16 &
     -4.69410E-16 -4.69160E-16 -4.68900E-16 -4.68640E-16 -4.68380E-16 -4.68120E-16 -4.67850E-16 -4.67580E-16 &
     -4.67300E-16 -4.67030E-16 -4.66740E-16 -4.66460E-16 -4.66170E-16 -4.65870E-16 -4.65570E-16 -4.65270E-16 &
     -4.64960E-16 -4.64640E-16 -4.64320E-16 -4.63990E-16 -4.63660E-16 -4.63330E-16 -4.62980E-16 -4.62630E-16 &
     -4.62280E-16 -4.61910E-16 -4.61540E-16 -4.61170E-16 -4.60790E-16 -4.60400E-16 -4.60000E-16 -4.59590E-16 &
     -4.59180E-16 -4.58760E-16 -4.58330E-16 -4.57900E-16 -4.57450E-16 -4.57000E-16 -4.56540E-16 -4.56060E-16 &
     -4.55580E-16 -4.55100E-16 -4.54600E-16 -4.54090E-16 -4.53570E-16 -4.53040E-16 -4.52510E-16 -4.51960E-16 &
     -4.51400E-16 -4.50830E-16 -4.50250E-16 -4.49660E-16 -4.49060E-16 -4.48450E-16 -4.47830E-16 -4.47190E-16 &
     -4.46550E-16 -4.45890E-16 -4.45220E-16 -4.44540E-16 -4.43840E-16 -4.43130E-16 -4.42410E-16 -4.41680E-16 &
     -4.40930E-16 -4.40170E-16 -4.39400E-16 -4.38620E-16 -4.37820E-16 -4.37000E-16 -4.36180E-16 -4.35330E-16 &
     -4.34480E-16 -4.33610E-16 -4.32720E-16 -4.31820E-16 -4.30910E-16 -4.29980E-16 -4.29030E-16 -4.28070E-16 &
     -4.27100E-16 -4.26100E-16 -4.25090E-16 -4.24070E-16 -4.23030E-16 -4.21970E-16 -4.20900E-16 -4.19810E-16 &
     -4.18700E-16 -4.17570E-16 -4.16430E-16 -4.15270E-16 -4.14090E-16 -4.12900E-16 -4.11690E-16 -4.10460E-16 &
     -4.09210E-16 -4.07950E-16 -4.06680E-16 -4.05380E-16 -4.04070E-16 -4.02750E-16 -4.01410E-16 -4.00060E-16 &
     -3.98690E-16 -3.97300E-16 -3.95900E-16 -3.94490E-16 -3.93060E-16 -3.91620E-16 -3.90160E-16 -3.88690E-16 &
     -3.87210E-16 -3.85720E-16 -3.84210E-16 -3.82680E-16 -3.81150E-16 -3.79600E-16 -3.78040E-16 -3.76470E-16 &
     -3.74890E-16 -3.73300E-16 -3.71690E-16 -3.70070E-16 -3.68440E-16 -3.66800E-16 -3.65150E-16 -3.63490E-16 &
     -3.61820E-16 -3.60140E-16 -3.58450E-16 -3.56750E-16 -3.55040E-16 -3.53320E-16 -3.51590E-16 -3.49850E-16 &
     -3.48110E-16 -3.46350E-16 -3.44590E-16 -3.42820E-16 -3.41040E-16 -3.39250E-16 -3.37460E-16 -3.35650E-16 &
     -3.33840E-16 -3.32030E-16 -3.30210E-16 -3.28380E-16 -3.26540E-16 -3.24700E-16 -3.22850E-16 -3.20990E-16 &
     -3.19130E-16 -3.17270E-16 -3.15400E-16 -3.13520E-16 -3.11640E-16 -3.09760E-16 -3.07870E-16 -3.05970E-16 &
     -3.04080E-16 -3.02170E-16 -3.00270E-16 -2.98360E-16 -2.96450E-16 -2.94530E-16 -2.92610E-16 -2.90690E-16 &
     -2.88770E-16 -2.86840E-16 -2.84920E-16 -2.82990E-16 -2.81050E-16 -2.79120E-16 -2.77190E-16 -2.75250E-16 &
     -2.73310E-16 -2.71380E-16 -2.69440E-16 -2.67500E-16 -2.65560E-16 -2.63620E-16 -2.61680E-16 -2.59750E-16 &
     -2.57810E-16 -2.55870E-16 -2.53940E-16 -2.52000E-16 -2.50070E-16 -2.48140E-16 -2.46210E-16 -2.44280E-16 &
     -2.42350E-16 -2.40430E-16 -2.38510E-16 -2.36590E-16 -2.34680E-16 -2.32760E-16 -2.30860E-16 -2.28950E-16 &
     -2.27050E-16 -2.25150E-16 -2.23260E-16 -2.21370E-16 -2.19480E-16 -2.17600E-16 -2.15730E-16 -2.13860E-16 &
     -2.11990E-16 -2.10130E-16 -2.08280E-16 -2.06430E-16 -2.04590E-16 -2.02760E-16 -2.00930E-16 -1.99100E-16 &
     -1.97290E-16 -1.95480E-16 -1.93680E-16 -1.91880E-16 -1.90100E-16 -1.88320E-16 -1.86550E-16 -1.84790E-16 &
     -1.83030E-16 -1.81290E-16 -1.79550E-16 -1.77820E-16 -1.76100E-16 -1.74390E-16 -1.72690E-16 -1.71000E-16 &
     -1.69320E-16 -1.67650E-16 -1.65990E-16 -1.64350E-16 -1.62710E-16 -1.61080E-16 -1.59460E-16 -1.57860E-16 &
     -1.56260E-16 -1.54680E-16 -1.53110E-16 -1.51550E-16 -1.50010E-16 -1.48480E-16 -1.46950E-16 -1.45450E-16 &
     -1.43950E-16 -1.42470E-16 -1.41000E-16 -1.39550E-16 -1.38110E-16 -1.36680E-16 -1.35270E-16 -1.33870E-16 &
     -1.32490E-16 -1.31120E-16 -1.29770E-16 -1.28430E-16 -1.27110E-16 -1.25800E-16 -1.24510E-16 -1.23230E-16 &
     -1.21970E-16 -1.20720E-16 -1.19490E-16 -1.18280E-16 -1.17080E-16 -1.15890E-16 -1.14720E-16 -1.13560E-16 &
     -1.12420E-16 -1.11300E-16 -1.10190E-16 -1.09090E-16 -1.08010E-16 -1.06940E-16 -1.05880E-16 -1.04840E-16 &
     -1.03820E-16 -1.02810E-16 -1.01810E-16 -1.00820E-16 -9.98540E-17 -9.88970E-17 -9.79550E-17 -9.70260E-17 &
     -9.61100E-17 -9.52080E-17 -9.43190E-17 -9.34440E-17 -9.25820E-17 -9.17330E-17 -9.08970E-17 -9.00740E-17 &
     -8.92640E-17 -8.84670E-17 -8.76820E-17 -8.69110E-17 -8.61520E-17 -8.54050E-17 -8.46710E-17 -8.39490E-17 &
     -8.32400E-17 -8.25430E-17 -8.18580E-17 -8.11850E-17 -8.05240E-17 -7.98750E-17 -7.92370E-17 -7.86120E-17 &
     -7.79980E-17 -7.73960E-17 -7.68060E-17 -7.62270E-17 -7.56590E-17 -7.51030E-17 -7.45570E-17 -7.40230E-17 &
     -7.35010E-17 -7.29890E-17 -7.24880E-17 -7.19980E-17 -7.15190E-17 -7.10500E-17 -7.05920E-17 -7.01450E-17 &
     -6.97080E-17 -6.92810E-17 -6.88650E-17 -6.84590E-17 -6.80640E-17 -6.76780E-17 -6.73030E-17 -6.69370E-17 &
     -6.65810E-17 -6.62350E-17 -6.58990E-17 -6.55720E-17 -6.52550E-17 -6.49480E-17 -6.46500E-17 -6.43610E-17 &
     -6.40810E-17 -6.38110E-17 -6.35500E-17 -6.32980E-17 -6.30550E-17 -6.28200E-17 -6.25950E-17 -6.23780E-17 &
     -6.21700E-17 -6.19700E-17 -6.17790E-17 -6.15970E-17 -6.14220E-17 -6.12560E-17 -6.10990E-17 -6.09490E-17 &
     -6.08070E-17 -6.06740E-17 -6.05480E-17 -6.04300E-17 -6.03200E-17 -6.02170E-17 -6.01220E-17 -6.00350E-17 &
     -5.99550E-17 -5.98820E-17 -5.98170E-17 -5.97590E-17 -5.97080E-17 -5.96640E-17 -5.96270E-17 -5.95970E-17 &
     -5.95740E-17 -5.95570E-17 -5.95470E-17 -5.95440E-17 -5.95480E-17 -5.95570E-17 -5.95730E-17 -5.95960E-17 &
     -5.96250E-17 -5.96590E-17 -5.97000E-17 -5.97470E-17 -5.98000E-17 -5.98590E-17 -5.99230E-17 -5.99930E-17 &
     -6.00690E-17 -6.01500E-17 -6.02370E-17 -6.03290E-17 -6.04260E-17 -6.05290E-17 -6.06370E-17 -6.07500E-17 &
     -6.08680E-17 -6.09900E-17 -6.11180E-17 -6.12510E-17 -6.13880E-17 -6.15290E-17 -6.16760E-17 -6.18270E-17 &
     -6.19820E-17 -6.21410E-17 -6.23050E-17 -6.24730E-17 -6.26450E-17 -6.28210E-17 -6.30010E-17 -6.31850E-17 &
     -6.33720E-17 -6.35640E-17 -6.37590E-17 -6.39570E-17 -6.41590E-17 -6.43640E-17 -6.45730E-17 -6.47850E-17 &
     -6.50000E-17 -6.52190E-17 -6.54400E-17 -6.56640E-17 -6.58910E-17 -6.61210E-17 -6.63540E-17 -6.65890E-17 &
     -6.68270E-17 -6.70670E-17 -6.73100E-17 -6.75550E-17 -6.78020E-17 -6.80520E-17 -6.83040E-17 -6.85570E-17 &
     -6.88130E-17 -6.90700E-17 -6.93300E-17 -6.95910E-17 -6.98540E-17 -7.01180E-17 -7.03840E-17 -7.06510E-17 &
     -7.09200E-17 -7.11900E-17 -7.14610E-17 -7.17330E-17 -7.20060E-17 -7.22800E-17 -7.25550E-17 -7.28310E-17 &
     -7.31080E-17 -7.33850E-17 -7.36630E-17 -7.39420E-17 -7.42210E-17 -7.45000E-17 -7.47800E-17 -7.50590E-17 &
     -7.53390E-17 -7.56190E-17 -7.58990E-17 -7.61790E-17 -7.64590E-17 -7.67380E-17 -7.70170E-17 -7.72960E-17 &
     -7.75740E-17 -7.78520E-17 -7.81290E-17 -7.84050E-17 -7.86810E-17 -7.89560E-17 -7.92300E-17 -7.95020E-17 &
     -7.97740E-17 -8.00450E-17 -8.03140E-17 -8.05820E-17 -8.08490E-17 -8.11140E-17 -8.13770E-17 -8.16390E-17 &
     -8.19000E-17 -8.21580E-17 -8.24150E-17 -8.26700E-17 -8.29230E-17 -8.31730E-17 -8.34220E-17 -8.36680E-17 &
     -8.39120E-17 -8.41540E-17 -8.43930E-17 -8.46300E-17 -8.48640E-17 -8.50950E-17 -8.53240E-17 -8.55500E-17 &
     -8.57730E-17 -8.59930E-17 -8.62100E-17 -8.64230E-17 -8.66340E-17 -8.68410E-17 -8.70450E-17 -8.72450E-17 &
     -8.74420E-17 -8.76360E-17 -8.78250E-17 -8.80110E-17 -8.81930E-17 -8.83720E-17 -8.85460E-17 -8.87170E-17 &
     -8.88840E-17 -8.90470E-17 -8.92060E-17 -8.93620E-17 -8.95140E-17 -8.96620E-17 -8.98070E-17 -8.99480E-17 &
     -9.00850E-17 -9.02180E-17 -9.03480E-17 -9.04740E-17 -9.05960E-17 -9.07150E-17 -9.08310E-17 -9.09420E-17 &
     -9.10500E-17 -9.11550E-17 -9.12560E-17 -9.13530E-17 -9.14470E-17 -9.15370E-17 -9.16230E-17 -9.17070E-17 &
     -9.17860E-17 -9.18620E-17 -9.19350E-17 -9.20040E-17 -9.20700E-17 -9.21320E-17 -9.21910E-17 -9.22460E-17 &
     -9.22980E-17 -9.23470E-17 -9.23920E-17 -9.24340E-17 -9.24720E-17 -9.25070E-17 -9.25390E-17 -9.25670E-17 &
     -9.25920E-17 -9.26140E-17 -9.26320E-17 -9.26470E-17 -9.26590E-17 -9.26680E-17 -9.26730E-17 -9.26750E-17 &
     -9.26740E-17 -9.26690E-17 -9.26620E-17 -9.26510E-17 -9.26370E-17 -9.26190E-17 -9.25990E-17 -9.25750E-17 &
     -9.25490E-17 -9.25190E-17 -9.24860E-17 -9.24490E-17 -9.24100E-17 -9.23680E-17 -9.23220E-17 -9.22740E-17 &
     -9.22220E-17 -9.21680E-17 -9.21100E-17 -9.20490E-17 -9.19860E-17 -9.19190E-17 -9.18490E-17 -9.17760E-17 &
     -9.17010E-17 -9.16220E-17 -9.15400E-17 -9.14560E-17 -9.13680E-17 -9.12780E-17 -9.11850E-17 -9.10880E-17 &
     -9.09890E-17 -9.08870E-17 -9.07830E-17 -9.06750E-17 -9.05640E-17 -9.04510E-17 -9.03350E-17 -9.02160E-17 &
     -9.00940E-17 -8.99700E-17 -8.98420E-17 -8.97120E-17 -8.95790E-17 -8.94440E-17 -8.93050E-17 -8.91640E-17 &
     -8.90210E-17 -8.88740E-17 -8.87250E-17 -8.85730E-17 -8.84190E-17 -8.82620E-17 -8.81020E-17 -8.79400E-17 &
     -8.77750E-17 -8.76070E-17 -8.74370E-17 -8.72640E-17 -8.70890E-17 -8.69110E-17 -8.67300E-17 -8.65470E-17 &
     -8.63620E-17 -8.61740E-17 -8.59830E-17 -8.57900E-17 -8.55950E-17 -8.53970E-17 -8.51960E-17 -8.49930E-17 &
     -8.47880E-17 -8.45800E-17 -8.43700E-17 -8.41570E-17 -8.39420E-17 -8.37250E-17 -8.35050E-17 -8.32830E-17 &
     -8.30580E-17 -8.28320E-17 -8.26030E-17 -8.23710E-17 -8.21370E-17 -8.19010E-17 -8.16630E-17 -8.14220E-17 &
     -8.11790E-17 -8.09340E-17 -8.06870E-17 -8.04370E-17 -8.01850E-17 -7.99310E-17 -7.96750E-17 -7.94170E-17 &
     -7.91560E-17 -7.88940E-17 -7.86290E-17 -7.83620E-17 -7.80920E-17 -7.78210E-17 -7.75480E-17 -7.72720E-17 &
     -7.69950E-17 -7.67150E-17 -7.64340E-17 -7.61500E-17 -7.58640E-17 -7.55760E-17 -7.52860E-17 -7.49950E-17 &
     -7.47010E-17 -7.44050E-17 -7.41070E-17 -7.38070E-17 -7.35060E-17 -7.32020E-17 -7.28970E-17 -7.25890E-17 &
     -7.22800E-17 -7.19690E-17 -7.16550E-17 -7.13400E-17 -7.10240E-17 -7.07050E-17 -7.03840E-17 -7.00620E-17 &
     -6.97380E-17 -6.94120E-17 -6.90840E-17 -6.87540E-17 -6.84230E-17 -6.80900E-17 -6.77550E-17 -6.74180E-17 &
     -6.70800E-17 -6.67400E-17 -6.63980E-17 -6.60550E-17 -6.57100E-17 -6.53630E-17 -6.50140E-17 -6.46640E-17 &
     -6.43130E-17 -6.39590E-17 -6.36040E-17 -6.32480E-17 -6.28890E-17 -6.25300E-17 -6.21680E-17 -6.18050E-17 &
     -6.14410E-17 -6.10750E-17 -6.07080E-17 -6.03390E-17 -5.99680E-17 -5.95960E-17 -5.92220E-17 -5.88480E-17 &
     -5.84710E-17 -5.80930E-17 -5.77140E-17 -5.73330E-17 -5.69510E-17 -5.65670E-17 -5.61830E-17 -5.57960E-17 &
     -5.54080E-17 -5.50190E-17 -5.46290E-17 -5.42370E-17 -5.38440E-17 -5.34500E-17 -5.30540E-17 -5.26570E-17 &
     -5.22590E-17 -5.18590E-17 -5.14590E-17 -5.10570E-17 -5.06530E-17 -5.02490E-17 -4.98430E-17 -4.94360E-17 &
     -4.90280E-17 -4.86190E-17 -4.82080E-17 -4.77970E-17 -4.73840E-17 -4.69700E-17 -4.65550E-17 -4.61390E-17 &
     -4.57220E-17 -4.53030E-17 -4.48840E-17 -4.44630E-17 -4.40420E-17 -4.36190E-17 -4.31960E-17 -4.27710E-17 &
     -4.23450E-17 -4.19190E-17 -4.14910E-17 -4.10620E-17 -4.06330E-17 -4.02020E-17 -3.97710E-17 -3.93380E-17 &
     -3.89050E-17 -3.84710E-17 -3.80350E-17 -3.75990E-17 -3.71620E-17 -3.67240E-17 -3.62860E-17 -3.58460E-17 &
     -3.54060E-17 -3.49650E-17 -3.45230E-17 -3.40800E-17 -3.36360E-17 -3.31920E-17 -3.27470E-17 -3.23010E-17 &
     -3.18540E-17 -3.14070E-17 -3.09590E-17 -3.05100E-17 -3.00600E-17 -2.96100E-17 -2.91590E-17 -2.87080E-17 &
     -2.82550E-17 -2.78030E-17 -2.73490E-17 -2.68950E-17 -2.64400E-17 -2.59850E-17 -2.55290E-17 -2.50720E-17 &
     -2.46150E-17 -2.41580E-17 -2.37000E-17 -2.32410E-17 -2.27820E-17 -2.23220E-17 -2.18620E-17 -2.14010E-17 &
     -2.09400E-17 -2.04780E-17 -2.00160E-17 -1.95530E-17 -1.90900E-17 -1.86270E-17 -1.81630E-17 -1.76990E-17 &
     -1.72340E-17 -1.67690E-17 -1.63030E-17 -1.58380E-17 -1.53720E-17 -1.49050E-17 -1.44380E-17 -1.39710E-17 &
     -1.35040E-17 -1.30360E-17 -1.25680E-17 -1.21000E-17 -1.16320E-17 -1.11630E-17 -1.06940E-17 -1.02250E-17 &
     -9.75550E-18 -9.28590E-18 -8.81620E-18 -8.34620E-18 -7.87610E-18 -7.40580E-18 -6.93540E-18 -6.46490E-18 &
     -5.99420E-18 -5.52340E-18 -5.05250E-18 -4.58160E-18 -4.11050E-18 -3.63940E-18 -3.16820E-18 -2.69700E-18 &
     -2.22580E-18 -1.75450E-18 -1.28320E-18 4.89670E-08 5.96180E-08 6.87050E-08 7.65790E-08 8.34390E-08 &
     8.94670E-08 9.47930E-08 9.95240E-08 1.03730E-07 1.07500E-07 1.10840E-07 1.13860E-07 1.16530E-07 &
     1.18940E-07 1.21090E-07 1.22990E-07 1.24700E-07 1.26220E-07 1.27550E-07 1.28720E-07 1.29750E-07 &
     1.30650E-07 1.31430E-07 1.32090E-07 1.32650E-07 1.33110E-07 1.33490E-07 1.33780E-07 1.34000E-07 &
     1.34130E-07 1.34200E-07 1.34210E-07 1.34160E-07 1.34070E-07 1.33930E-07 1.33740E-07 1.33510E-07 &
     1.33240E-07 1.32930E-07 1.32590E-07 1.32220E-07 1.31820E-07 1.31390E-07 1.30920E-07 1.30440E-07 &
     1.29940E-07 1.29410E-07 1.28870E-07 1.28320E-07 1.27740E-07 1.27160E-07 1.26560E-07 1.25940E-07 &
     1.25320E-07 1.24690E-07 1.24040E-07 1.23390E-07 1.22730E-07 1.22060E-07 1.21380E-07 1.20700E-07 &
     1.20010E-07 1.19320E-07 1.18630E-07 1.17930E-07 1.17230E-07 1.16530E-07 1.15820E-07 1.15120E-07 &
     1.14410E-07 1.13700E-07 1.13000E-07 1.12290E-07 1.11580E-07 1.10880E-07 1.10170E-07 1.09470E-07 &
     1.08760E-07 1.08060E-07 1.07360E-07 1.06660E-07 1.05960E-07 1.05270E-07 1.04570E-07 1.03880E-07 &
     1.03200E-07 1.02510E-07 1.01830E-07 1.01160E-07 1.00480E-07 9.98130E-08 9.91470E-08 9.84850E-08 &
     9.78270E-08 9.71730E-08 9.65220E-08 9.58760E-08 9.52330E-08 9.45940E-08 9.39590E-08 9.33280E-08 &
     9.27010E-08 9.20780E-08 9.14590E-08 9.08450E-08 9.02340E-08 8.96280E-08 8.90270E-08 8.84300E-08 &
     8.78370E-08 8.72500E-08 8.66660E-08 8.60880E-08 8.55140E-08 8.49450E-08 8.43800E-08 8.38200E-08 &
     8.32640E-08 8.27130E-08 8.21660E-08 8.16230E-08 8.10850E-08 8.05500E-08 8.00200E-08 7.94940E-08 &
     7.89720E-08 7.84540E-08 7.79410E-08 7.74320E-08 7.69280E-08 7.64280E-08 7.59330E-08 7.54410E-08 &
     7.49550E-08 7.44720E-08 7.39930E-08 7.35190E-08 7.30490E-08 7.25830E-08 7.21210E-08 7.16630E-08 &
     7.12100E-08 7.07600E-08 7.03140E-08 6.98720E-08 6.94340E-08 6.90000E-08 6.85680E-08 6.81410E-08 &
     6.77160E-08 6.72950E-08 6.68780E-08 6.64640E-08 6.60530E-08 6.56460E-08 6.52440E-08 6.48450E-08 &
     6.44490E-08 6.40580E-08 6.36700E-08 6.32860E-08 6.29060E-08 6.25290E-08 6.21560E-08 6.17860E-08 &
     6.14190E-08 6.10550E-08 6.06950E-08 6.03380E-08 5.99830E-08 5.96320E-08 5.92840E-08 5.89380E-08 &
     5.85960E-08 5.82560E-08 5.79200E-08 5.75860E-08 5.72560E-08 5.69290E-08 5.66040E-08 5.62830E-08 &
     5.59640E-08 5.56490E-08 5.53360E-08 5.50260E-08 5.47190E-08 5.44150E-08 5.41130E-08 5.38140E-08 &
     5.35170E-08 5.32230E-08 5.29310E-08 5.26410E-08 5.23530E-08 5.20680E-08 5.17850E-08 5.15050E-08 &
     5.12270E-08 5.09520E-08 5.06800E-08 5.04100E-08 5.01420E-08 4.98770E-08 4.96130E-08 4.93520E-08 &
     4.90930E-08 4.88350E-08 4.85790E-08 4.83250E-08 4.80740E-08 4.78240E-08 4.75760E-08 4.73300E-08 &
     4.70870E-08 4.68450E-08 4.66050E-08 4.63680E-08 4.61320E-08 4.58990E-08 4.56680E-08 4.54380E-08 &
     4.52110E-08 4.49850E-08 4.47610E-08 4.45390E-08 4.43190E-08 4.41010E-08 4.38840E-08 4.36690E-08 &
     4.34550E-08 4.32430E-08 4.30330E-08 4.28240E-08 4.26170E-08 4.24120E-08 4.22080E-08 4.20060E-08 &
     4.18050E-08 4.16060E-08 4.14080E-08 4.12120E-08 4.10180E-08 4.08250E-08 4.06330E-08 4.04430E-08 &
     4.02550E-08 4.00670E-08 3.98820E-08 3.96970E-08 3.95140E-08 3.93320E-08 3.91520E-08 3.89720E-08 &
     3.87940E-08 3.86160E-08 3.84400E-08 3.82640E-08 3.80900E-08 3.79170E-08 3.77450E-08 3.75750E-08 &
     3.74070E-08 3.72390E-08 3.70740E-08 3.69090E-08 3.67460E-08 3.65830E-08 3.64220E-08 3.62620E-08 &
     3.61030E-08 3.59450E-08 3.57880E-08 3.56320E-08 3.54770E-08 3.53230E-08 3.51700E-08 3.50180E-08 &
     3.48670E-08 3.47180E-08 3.45690E-08 3.44210E-08 3.42750E-08 3.41290E-08 3.39840E-08 3.38410E-08 &
     3.36980E-08 3.35560E-08 3.34150E-08 3.32750E-08 3.31360E-08 3.29980E-08 3.28610E-08 3.27240E-08 &
     3.25880E-08 3.24530E-08 3.23190E-08 3.21850E-08 3.20520E-08 3.19200E-08 3.17880E-08 3.16570E-08 &
     3.15270E-08 3.13990E-08 3.12710E-08 3.11440E-08 3.10170E-08 3.08920E-08 3.07680E-08 3.06440E-08 &
     3.05210E-08 3.03980E-08 3.02770E-08 3.01560E-08 3.00350E-08 2.99150E-08 2.97960E-08 2.96780E-08 &
     2.95610E-08 2.94440E-08 2.93280E-08 2.92120E-08 2.90970E-08 2.89830E-08 2.88700E-08 2.87570E-08 &
     2.86450E-08 2.85340E-08 2.84230E-08 2.83130E-08 2.82040E-08 2.80950E-08 2.79860E-08 2.78780E-08 &
     2.77710E-08 2.76640E-08 2.75580E-08 2.74530E-08 2.73480E-08 2.72440E-08 2.71400E-08 2.70370E-08 &
     2.69350E-08 2.68330E-08 2.67310E-08 2.66300E-08 2.65290E-08 2.64280E-08 2.63280E-08 2.62280E-08 &
     2.61290E-08 2.60300E-08 2.59320E-08 2.58350E-08 2.57390E-08 2.56430E-08 2.55480E-08 2.54530E-08 &
     2.53590E-08 2.52650E-08 2.51720E-08 2.50790E-08 2.49860E-08 2.48940E-08 2.48020E-08 2.47110E-08 &
     2.46200E-08 2.45290E-08 2.44400E-08 2.43500E-08 2.42620E-08 2.41740E-08 2.40860E-08 2.39990E-08 &
     2.39120E-08 2.38260E-08 2.37400E-08 2.36550E-08 2.35690E-08 2.34850E-08 2.34000E-08 2.33150E-08 &
     2.32310E-08 2.31470E-08 2.30630E-08 2.29800E-08 2.28970E-08 2.28150E-08 2.27330E-08 2.26520E-08 &
     2.25720E-08 2.24910E-08 2.24120E-08 2.23320E-08 2.22540E-08 2.21750E-08 2.20970E-08 2.20190E-08 &
     2.19420E-08 2.18650E-08 2.17880E-08 2.17110E-08 2.16350E-08 2.15590E-08 2.14840E-08 2.14090E-08 &
     2.13340E-08 2.12590E-08 2.11850E-08 2.11120E-08 2.10380E-08 2.09650E-08 2.08930E-08 2.08210E-08 &
     2.07490E-08 2.06780E-08 2.06070E-08 2.05360E-08 2.04660E-08 2.03960E-08 2.03260E-08 2.02560E-08 &
     2.01870E-08 2.01170E-08 2.00480E-08 1.99790E-08 1.99100E-08 1.98420E-08 1.97740E-08 1.97060E-08 &
     1.96390E-08 1.95720E-08 1.95060E-08 1.94400E-08 1.93740E-08 1.93090E-08 1.92440E-08 1.91790E-08 &
     1.91150E-08 1.90500E-08 1.89870E-08 1.89230E-08 1.88600E-08 1.87970E-08 1.87340E-08 1.86710E-08 &
     1.86090E-08 1.85470E-08 1.84850E-08 1.84240E-08 1.83620E-08 1.83010E-08 1.82410E-08 1.81810E-08 &
     1.81200E-08 1.80610E-08 1.80010E-08 1.79420E-08 1.78830E-08 1.78250E-08 1.77660E-08 1.77080E-08 &
     1.76500E-08 1.75930E-08 1.75350E-08 1.74780E-08 1.74210E-08 1.73630E-08 1.73070E-08 1.72500E-08 &
     1.71930E-08 1.71370E-08 1.70810E-08 1.70260E-08 1.69710E-08 1.69160E-08 1.68610E-08 1.68070E-08 &
     1.67520E-08 1.66980E-08 1.66450E-08 1.65910E-08 1.65380E-08 1.64860E-08 1.64330E-08 1.63810E-08 &
     1.63280E-08 1.62760E-08 1.62250E-08 1.61730E-08 1.61220E-08 1.60710E-08 1.60200E-08 1.59700E-08 &
     1.59190E-08 1.58690E-08 1.58190E-08 1.57690E-08 1.57190E-08 1.56700E-08 1.56200E-08 1.55710E-08 &
     1.55220E-08 1.54730E-08 1.54240E-08 1.53760E-08 1.53280E-08 1.52800E-08 1.52320E-08 1.51850E-08 &
     1.51380E-08 1.50910E-08 1.50450E-08 1.49980E-08 1.49520E-08 1.49060E-08 1.48600E-08 1.48140E-08 &
     1.47680E-08 1.47230E-08 1.46770E-08 1.46320E-08 1.45880E-08 1.45430E-08 1.44990E-08 1.44550E-08 &
     1.44110E-08 1.43670E-08 1.43240E-08 1.42800E-08 1.42370E-08 1.41940E-08 1.41520E-08 1.41090E-08 &
     1.40660E-08 1.40240E-08 1.39810E-08 1.39390E-08 1.38970E-08 1.38550E-08 1.38130E-08 1.37720E-08 &
     1.37300E-08 1.36890E-08 1.36480E-08 1.36080E-08 1.35670E-08 1.35270E-08 1.34870E-08 1.34470E-08 &
     1.34080E-08 1.33680E-08 1.33290E-08 1.32890E-08 1.32500E-08 1.32110E-08 1.31730E-08 1.31340E-08 &
     1.30960E-08 1.30580E-08 1.30200E-08 1.29820E-08 1.29440E-08 1.29070E-08 1.28690E-08 1.28320E-08 &
     1.27950E-08 1.27580E-08 1.27210E-08 1.26840E-08 1.26470E-08 1.26100E-08 1.25740E-08 1.25380E-08 &
     1.25020E-08 1.24660E-08 1.24300E-08 1.23940E-08 1.23590E-08 1.23240E-08 1.22890E-08 1.22540E-08 &
     1.22190E-08 1.21850E-08 1.21500E-08 1.21160E-08 1.20820E-08 1.20480E-08 1.20140E-08 1.19800E-08 &
     1.19460E-08 1.19130E-08 1.18790E-08 1.18460E-08 1.18130E-08 1.17800E-08 1.17470E-08 1.17150E-08 &
     1.16820E-08 1.16500E-08 1.16180E-08 1.15860E-08 1.15540E-08 1.15220E-08 1.14900E-08 1.14590E-08 &
     1.14270E-08 1.13960E-08 1.13650E-08 1.13340E-08 1.13030E-08 1.12720E-08 1.12410E-08 1.12100E-08 &
     1.11790E-08 1.11490E-08 1.11180E-08 1.10880E-08 1.10580E-08 1.10280E-08 1.09980E-08 1.09680E-08 &
     1.09390E-08 1.09100E-08 1.08800E-08 1.08510E-08 1.08220E-08 1.07940E-08 1.07650E-08 1.07360E-08 &
     1.07080E-08 1.06790E-08 1.06510E-08 1.06230E-08 1.05940E-08 1.05660E-08 1.05380E-08 1.05110E-08 &
     1.04830E-08 1.04550E-08 1.04280E-08 1.04010E-08 1.03740E-08 1.03470E-08 1.03200E-08 1.02930E-08 &
     1.02660E-08 1.02390E-08 1.02130E-08 1.01860E-08 1.01590E-08 1.01330E-08 1.01060E-08 1.00800E-08 &
     1.00540E-08 1.00280E-08 1.00020E-08 9.97650E-09 9.95080E-09 9.92530E-09 9.90000E-09 9.87470E-09 &
     9.84950E-09 9.82450E-09 9.79950E-09 9.77470E-09 9.74990E-09 9.72530E-09 9.70080E-09 9.67640E-09 &
     9.65200E-09 9.62780E-09 9.60360E-09 9.57950E-09 9.55550E-09 9.53150E-09 9.50760E-09 9.48380E-09 &
     9.46010E-09 9.43660E-09 9.41310E-09 9.38970E-09 9.36650E-09 9.34340E-09 9.32030E-09 9.29720E-09 &
     9.27400E-09 9.25060E-09 9.22710E-09 9.20330E-09 9.17960E-09 9.15610E-09 9.13310E-09 9.11100E-09 &
     9.08980E-09 9.07000E-09 9.05160E-09 9.03430E-09 9.01780E-09 9.00150E-09 8.98500E-09 8.96800E-09 &
     8.95000E-09 8.93060E-09 8.91010E-09 8.88870E-09 8.86670E-09 8.84430E-09 8.82190E-09 8.79960E-09 &
     8.77790E-09 8.75650E-09 8.73540E-09 8.71460E-09 8.69400E-09 8.67350E-09 8.65310E-09 8.63260E-09 &
     8.61220E-09 8.59170E-09 8.57120E-09 8.55080E-09 8.53050E-09 8.51020E-09 8.49000E-09 8.47000E-09 &
     8.45000E-09 8.43010E-09 8.41020E-09 8.39050E-09 8.37080E-09 8.35120E-09 8.33160E-09 8.31210E-09 &
     8.29270E-09 8.27330E-09 8.25400E-09 8.23480E-09 8.21560E-09 8.19650E-09 8.17750E-09 8.15850E-09 &
     8.13970E-09 8.12080E-09 8.10210E-09 8.08340E-09 8.06470E-09 8.04620E-09 8.02760E-09 8.00920E-09 &
     7.99080E-09 7.97250E-09 7.95420E-09 7.93600E-09 7.91790E-09 7.89980E-09 7.88180E-09 7.86390E-09 &
     7.84600E-09 7.82820E-09 7.81040E-09 7.79270E-09 7.77500E-09 7.75740E-09 7.73990E-09 7.72240E-09 &
     7.70500E-09 7.68770E-09 7.67030E-09 7.65310E-09 7.63590E-09 7.61880E-09 7.60170E-09 7.58470E-09 &
     7.56770E-09 7.55080E-09 7.53400E-09 7.51720E-09 7.50040E-09 7.48370E-09 7.46710E-09 7.45050E-09 &
     7.43400E-09 7.41750E-09 7.40110E-09 7.38470E-09 7.36840E-09 7.35220E-09 7.33590E-09 7.31980E-09 &
     7.30370E-09 7.28760E-09 7.27160E-09 7.25570E-09 7.23980E-09 7.22390E-09 7.20810E-09 7.19230E-09 &
     7.17660E-09 7.16100E-09 7.14540E-09 7.12980E-09 7.11430E-09 7.09880E-09 7.08340E-09 7.06810E-09 &
     7.05270E-09 7.03750E-09 7.02230E-09 7.00710E-09 6.99190E-09 6.97690E-09 6.96180E-09 6.94680E-09 &
     6.93190E-09 6.91700E-09 6.90210E-09 6.88730E-09 6.87260E-09 6.85790E-09 6.84320E-09 6.82860E-09 &
     6.81400E-09 6.79940E-09 6.78490E-09 6.77050E-09 6.75600E-09 6.74170E-09 6.72740E-09 6.71310E-09 &
     6.69880E-09 6.68460E-09 6.67050E-09 6.65640E-09 6.64230E-09 6.62820E-09 6.61430E-09 6.60030E-09 &
     6.58640E-09 6.57250E-09 6.55870E-09 6.54490E-09 6.53120E-09 6.51750E-09 6.50380E-09 6.49020E-09 &
     6.47660E-09 6.46310E-09 6.44960E-09 6.43610E-09 6.42270E-09 6.40930E-09 6.39590E-09 6.38260E-09 &
     6.36930E-09 6.35610E-09 6.34290E-09 6.32970E-09 6.31660E-09 6.30350E-09 6.29050E-09 6.27750E-09 &
     6.26450E-09 6.25160E-09 6.23870E-09 6.22580E-09 6.21300E-09 6.20020E-09 6.18750E-09 6.17470E-09 &
     6.16210E-09 6.14940E-09 6.13680E-09 6.12420E-09 6.11170E-09 6.09920E-09 6.08670E-09 6.07430E-09 &
     6.06190E-09 6.04950E-09 6.03720E-09 6.02490E-09 6.01260E-09 6.00040E-09 5.98820E-09 5.97610E-09 &
     5.96390E-09 5.95190E-09 5.93980E-09 5.92780E-09 5.91580E-09 5.90380E-09 5.89190E-09 5.88000E-09 &
     5.86810E-09 5.85630E-09 5.84450E-09 5.83280E-09 5.82100E-09 5.80930E-09 5.79770E-09 5.78600E-09 &
     5.77440E-09 5.76280E-09 5.75130E-09 5.73980E-09 5.72830E-09 5.71690E-09 5.70550E-09 5.69410E-09 &
     5.68270E-09 5.67140E-09 5.66010E-09 5.64880E-09 5.63760E-09 5.62640E-09 5.61520E-09 5.60410E-09 &
     5.59300E-09 5.58190E-09 5.57080E-09 5.55980E-09 5.54880E-09 5.53790E-09 5.52690E-09 5.51600E-09 &
     5.50510E-09 5.49430E-09 5.48350E-09 5.47270E-09 5.46190E-09 5.45120E-09 5.44050E-09 5.42980E-09 &
     5.41910E-09 5.40850E-09 5.39790E-09 5.38740E-09 5.37680E-09 5.36630E-09 5.35580E-09 5.34540E-09 &
     5.33500E-09 5.32460E-09 5.31420E-09 5.30380E-09 5.29350E-09 5.28320E-09 5.27300E-09 5.26270E-09 &
     5.25250E-09 5.24230E-09 5.23220E-09 5.22210E-09 5.21190E-09 5.20190E-09 5.19180E-09 5.18180E-09 &
     5.17180E-09 5.16180E-09 5.15190E-09 5.14200E-09 5.13210E-09 5.12220E-09 5.11240E-09 5.10250E-09 &
     5.09270E-09 5.08300E-09 5.07320E-09 5.06350E-09 5.05380E-09 5.04420E-09 5.03450E-09 5.02490E-09 &
     5.01530E-09 5.00580E-09 4.99620E-09 4.98670E-09 4.97720E-09 4.96770E-09 4.95830E-09 4.94890E-09 &
     4.93950E-09 4.93010E-09 4.92080E-09 4.91140E-09 4.90220E-09 4.89290E-09 4.88360E-09 4.87440E-09 &
     4.86520E-09 4.85600E-09 4.84690E-09 4.83770E-09 4.82860E-09 4.81950E-09 4.81050E-09 4.80140E-09 &
     4.79240E-09 4.78340E-09 4.77450E-09 4.76550E-09 4.75660E-09 4.74770E-09 4.73880E-09 4.73000E-09 &
     4.72110E-09 4.71230E-09 4.70350E-09 4.69480E-09 4.68600E-09 4.67730E-09 4.66860E-09 4.65990E-09 &
     4.65130E-09 4.64260E-09 4.63400E-09 4.62540E-09 4.61690E-09 4.60830E-09 4.59980E-09 4.59130E-09 &
     4.58280E-09 4.57440E-09 4.56590E-09 4.55750E-09 4.54910E-09 4.54070E-09 4.53240E-09 4.52410E-09 &
     4.51570E-09 4.50750E-09 4.49920E-09 4.49090E-09 4.48270E-09 4.47450E-09 4.46630E-09 4.45820E-09 &
     4.45000E-09 4.44190E-09 4.43380E-09 4.42570E-09 4.41770E-09 4.40960E-09 4.40160E-09 4.39360E-09 &
     4.38560E-09 4.37770E-09 4.36970E-09 4.36180E-09 4.35390E-09 4.34600E-09 4.33820E-09 4.33030E-09 &
     4.32250E-09 4.31470E-09 4.30690E-09 4.29910E-09 4.29140E-09 4.28370E-09 4.27600E-09 4.26830E-09 &
     4.26060E-09 4.25300E-09 4.24540E-09 4.23770E-09 4.23020E-09 4.22260E-09 4.21500E-09 4.20750E-09 &
     4.20000E-09 4.19250E-09 4.18500E-09 4.17760E-09 4.17010E-09 4.16270E-09 4.15530E-09 4.14790E-09 &
     4.14060E-09 4.13320E-09 4.12590E-09 4.11860E-09 4.11130E-09 4.10400E-09 4.09680E-09 4.08950E-09 &
     4.08230E-09 4.07510E-09 4.06790E-09 4.06080E-09 4.05360E-09 4.04650E-09 4.03940E-09 4.03230E-09 &
     4.02520E-09 4.01820E-09 4.01110E-09 4.00410E-09 3.99710E-09 3.99010E-09 3.98320E-09 3.97620E-09 &
     3.96930E-09 3.96240E-09 3.95550E-09 3.94860E-09 3.94170E-09 3.93490E-09 3.92800E-09 3.92120E-09 &
     3.91440E-09 3.90760E-09 3.90090E-09 3.89410E-09 3.88740E-09 3.88070E-09 3.87400E-09 3.86730E-09 &
     3.86060E-09 3.85400E-09 3.84730E-09 3.84070E-09 3.83410E-09 3.82750E-09 3.82100E-09 3.81440E-09 &
     3.80790E-09 3.80140E-09 3.79490E-09 3.78840E-09 3.78190E-09 3.77540E-09 3.76900E-09 3.76260E-09 &
     3.75620E-09 3.74980E-09 3.74340E-09 3.73700E-09 3.73070E-09 3.72440E-09 3.71800E-09 3.71170E-09 &
     3.70550E-09 3.69920E-09 3.69290E-09 3.68670E-09 3.68050E-09 3.67430E-09 3.66810E-09 3.66190E-09 &
     3.65570E-09 3.64960E-09 3.64340E-09 3.63730E-09 3.63120E-09 3.62510E-09 3.61910E-09 3.61300E-09 &
     3.60700E-09 3.60090E-09 3.59490E-09 3.58890E-09 3.58290E-09 3.57700E-09 3.57100E-09 3.56510E-09 &
     3.55910E-09 3.55320E-09 3.54730E-09 3.54140E-09 3.53560E-09 3.52970E-09 3.52390E-09 3.51800E-09 &
     3.51220E-09 3.50640E-09 3.50060E-09 3.49490E-09 3.48910E-09 3.48340E-09 3.47760E-09 3.47190E-09 &
     3.46620E-09 3.46050E-09 3.45480E-09 3.44920E-09 3.44350E-09 3.43790E-09 3.43230E-09 3.42670E-09 &
     3.42110E-09 3.41550E-09 3.40990E-09 3.40440E-09 3.39880E-09 3.39330E-09 3.38780E-09 3.38230E-09 &
     3.37680E-09 3.37130E-09 3.36580E-09 3.36040E-09 3.35490E-09 3.34950E-09 3.34410E-09 3.33870E-09 &
     3.33330E-09 3.32790E-09 3.32260E-09 3.31720E-09 3.31190E-09 3.30650E-09 3.30120E-09 3.29590E-09 &
     3.29060E-09 3.28540E-09 3.28010E-09 3.27480E-09 3.26960E-09 3.26440E-09 3.25920E-09 3.25400E-09 &
     3.24880E-09 3.24360E-09 3.23840E-09 3.23330E-09 3.22810E-09 3.22300E-09 3.21790E-09 3.21280E-09 &
     3.20770E-09 3.20260E-09 3.19750E-09 3.19250E-09 3.18740E-09 3.18240E-09 3.17730E-09 3.17230E-09 &
     3.16730E-09 3.16230E-09 3.15740E-09 3.15240E-09 3.14740E-09 3.14250E-09 3.13750E-09 3.13260E-09 &
     3.12770E-09 3.12280E-09 3.11790E-09 3.11300E-09 3.10820E-09 3.10330E-09 3.09850E-09 3.09360E-09 &
     3.08880E-09 3.08400E-09 3.07920E-09 3.07440E-09 3.06960E-09 3.06490E-09 3.06010E-09 3.05540E-09 &
     3.05060E-09 3.04590E-09 3.04120E-09 3.03650E-09 3.03180E-09 3.02710E-09 3.02240E-09 3.01780E-09 &
     3.01310E-09 3.00850E-09 3.00380E-09 2.99920E-09 2.99460E-09 2.99000E-09 2.98540E-09 2.98080E-09 &
     2.97620E-09 2.97170E-09 2.96710E-09 2.96260E-09 2.95810E-09 2.95350E-09 2.94900E-09 2.94450E-09 &
     2.94000E-09 2.93550E-09 2.93110E-09 2.92660E-09 2.92220E-09 2.91770E-09 2.91330E-09 2.90890E-09 &
     2.90440E-09 2.90000E-09 2.89560E-09 2.89130E-09 2.88690E-09 2.88250E-09 2.87820E-09 2.87380E-09 &
     2.86950E-09 2.86520E-09 2.86080E-09 2.85650E-09 2.85220E-09 2.84790E-09 2.84360E-09 2.83940E-09 &
     2.83510E-09 2.83090E-09 2.82660E-09 2.82240E-09 2.81810E-09 2.81390E-09 2.80970E-09 2.80550E-09 &
     2.80130E-09 2.79710E-09 2.79300E-09 2.78880E-09 2.78470E-09 2.78050E-09 2.77640E-09 2.77220E-09 &
     2.76810E-09 2.76400E-09 2.75990E-09 2.75580E-09 2.75170E-09 2.74760E-09 2.74360E-09 2.73950E-09 &
     2.73550E-09 2.73140E-09 2.72740E-09 2.72340E-09 2.71930E-09 2.71530E-09 2.71130E-09 2.70730E-09 &
     2.70340E-09 2.69940E-09 2.69540E-09 2.69140E-09 2.68750E-09 2.68350E-09 2.67960E-09 2.67570E-09 &
     2.67180E-09 2.66780E-09 2.66390E-09 2.66000E-09 2.65620E-09 2.65230E-09 2.64840E-09 2.64450E-09 &
     2.64070E-09 2.63680E-09 2.63300E-09 2.62920E-09 2.62530E-09 2.62150E-09 2.61770E-09 2.61390E-09 &
     2.61010E-09 2.60630E-09 2.60260E-09 2.59880E-09 2.59500E-09 2.59130E-09 2.58750E-09 2.58380E-09 &
     2.58000E-09 2.57630E-09 2.57260E-09 2.56890E-09 2.56520E-09 2.56150E-09 2.55780E-09 2.55410E-09 &
     2.55040E-09 2.54680E-09 2.54310E-09 2.53950E-09 2.53580E-09 2.53220E-09 2.52860E-09 2.52490E-09 &
     2.52130E-09 2.51770E-09 2.51410E-09 2.51050E-09 2.50690E-09 2.50340E-09 2.49980E-09 2.49620E-09 &
     2.49270E-09 2.48910E-09 2.48560E-09 2.48200E-09 2.47850E-09 2.47500E-09 2.47150E-09 2.46790E-09 &
     2.46440E-09 2.46090E-09 2.45750E-09 2.45400E-09 2.45050E-09 2.44700E-09 2.44360E-09 2.44010E-09 &
     2.43670E-09 2.43320E-09 2.42980E-09 2.42640E-09 2.42290E-09 2.41950E-09 2.41610E-09 2.41270E-09 &
     2.40930E-09 2.40590E-09 2.40250E-09 2.39920E-09 2.39580E-09 2.39240E-09 2.38910E-09 2.38570E-09 &
     2.38240E-09 2.37900E-09 2.37570E-09 2.37240E-09 2.36910E-09 2.36580E-09 2.36240E-09 2.35910E-09 &
     2.35590E-09 2.35260E-09 2.34930E-09 2.34600E-09 2.34270E-09 2.33950E-09 2.33620E-09 2.33300E-09 &
     2.32970E-09 2.32650E-09 2.32330E-09 2.32000E-09 2.31680E-09 2.31360E-09 2.31040E-09 2.30720E-09 &
     2.30400E-09 2.30080E-09 2.29760E-09 2.29450E-09 2.29130E-09 2.28810E-09 2.28500E-09 2.28180E-09 &
     2.27870E-09 2.27550E-09 2.27240E-09 2.26930E-09 2.26610E-09 2.26300E-09 2.25990E-09 2.25680E-09 &
     2.25370E-09 2.25060E-09 2.24750E-09 2.24440E-09 2.24140E-09 2.23830E-09 2.23520E-09 2.23220E-09 &
     2.22910E-09 2.22610E-09 2.22300E-09 2.22000E-09 2.21690E-09 2.21390E-09 2.21090E-09 2.20790E-09 &
     2.20490E-09 2.20190E-09 2.19890E-09 2.19590E-09 2.19290E-09 2.18990E-09 2.18690E-09 2.18400E-09 &
     2.18100E-09 2.17800E-09 2.17510E-09 2.17210E-09 2.16920E-09 2.16620E-09 2.16330E-09 2.16040E-09 &
     2.15750E-09 2.15450E-09 2.15160E-09 2.14870E-09 2.14580E-09 2.14290E-09 2.14000E-09 2.13720E-09 &
     2.13430E-09 2.13140E-09 2.12850E-09 2.12570E-09 2.12280E-09 2.12000E-09 2.11710E-09 2.11430E-09 &
     2.11140E-09 2.10860E-09 2.10580E-09 2.10290E-09 2.10010E-09 2.09730E-09 2.09450E-09 2.09170E-09 &
     2.08890E-09 2.08610E-09 2.08330E-09 2.08050E-09 2.07780E-09 2.07500E-09 2.07220E-09 2.06950E-09 &
     2.06670E-09 2.06390E-09 2.06120E-09 2.05850E-09 2.05570E-09 2.05300E-09 2.05030E-09 2.04750E-09 &
     2.04480E-09 2.04210E-09 2.03940E-09 2.03670E-09 2.03400E-09 2.03130E-09 2.02860E-09 2.02590E-09 &
     2.02330E-09 2.02060E-09 2.01790E-09 2.01530E-09 2.01260E-09 2.00990E-09 2.00730E-09 2.00460E-09 &
     2.00200E-09 1.99940E-09 1.99670E-09 1.99410E-09 1.99150E-09 1.98890E-09 1.98630E-09 1.98370E-09 &
     1.98100E-09 1.97850E-09 1.97590E-09 1.97330E-09 1.97070E-09 1.96810E-09 1.96550E-09 1.96300E-09 &
     1.96040E-09 1.95780E-09 1.95530E-09 1.95270E-09 1.95020E-09 1.94760E-09 1.94510E-09 1.94260E-09 &
     1.94000E-09 1.93750E-09 1.93500E-09 1.93250E-09 1.92990E-09 1.92740E-09 1.92490E-09 1.92240E-09 &
     1.91990E-09 1.91740E-09 1.91500E-09 1.91250E-09 1.91000E-09 1.90750E-09 1.90510E-09 1.90260E-09 &
     1.90010E-09 1.89770E-09 1.89520E-09 1.89280E-09 1.89030E-09 1.88790E-09 1.88550E-09 1.88300E-09 &
     1.88060E-09 1.87820E-09 1.87580E-09 1.87340E-09 1.87100E-09 1.86860E-09 1.86620E-09 1.86380E-09 &
     1.86140E-09 1.85900E-09 1.85660E-09 1.85420E-09 1.85190E-09 1.84950E-09 1.84710E-09 1.84480E-09 &
     1.84240E-09 1.84000E-09 1.83770E-09 1.83540E-09 1.83300E-09 1.83070E-09 1.82830E-09 1.82600E-09 &
     1.82370E-09 1.82140E-09 1.81910E-09 1.81670E-09 1.81440E-09 1.81210E-09 1.80980E-09 1.80750E-09 &
     1.80520E-09 1.80290E-09 1.80070E-09 1.79840E-09 1.79610E-09 1.79380E-09 1.79160E-09 1.78930E-09 &
     1.78700E-09 1.78480E-09 1.78250E-09 1.78030E-09 1.77800E-09 1.77580E-09 1.77360E-09 1.77130E-09 &
     1.76910E-09 1.76690E-09 1.76460E-09 1.76240E-09 1.76020E-09 1.75800E-09 1.75580E-09 1.75360E-09 &
     1.75140E-09 1.74920E-09 1.74700E-09 1.74480E-09 1.74260E-09 1.74050E-09 1.73830E-09 1.73610E-09 &
     1.73390E-09 1.73180E-09 1.72960E-09 1.72750E-09 1.72530E-09 1.72320E-09 1.72100E-09 1.71890E-09 &
     1.71670E-09 1.71460E-09 1.71250E-09 1.71030E-09 1.70820E-09 1.70610E-09 1.70400E-09 1.70190E-09 &
     1.69970E-09 1.69760E-09 1.69550E-09 1.69340E-09 1.69130E-09 1.68930E-09 1.68720E-09 1.68510E-09 &
     1.68300E-09 1.68090E-09 1.67890E-09 1.67680E-09 1.67470E-09 1.67270E-09 1.67060E-09 1.66850E-09 &
     1.66650E-09 1.66440E-09 1.66240E-09 1.66040E-09 1.65830E-09 1.65630E-09 1.65420E-09 1.65220E-09 &
     1.65020E-09 1.64820E-09 1.64620E-09 1.64410E-09 1.64210E-09 1.64010E-09 1.63810E-09 1.63610E-09 &
     1.63410E-09 1.63210E-09 1.63010E-09 1.62820E-09 1.62620E-09 1.62420E-09 1.62220E-09 1.62020E-09 &
     1.61830E-09 1.61630E-09 1.61430E-09 1.61240E-09 1.61040E-09 1.60850E-09 1.60650E-09 1.60460E-09 &
     1.60260E-09 1.60070E-09 1.59880E-09 1.59680E-09 1.59490E-09 1.59300E-09 1.59100E-09 1.58910E-09 &
     1.58720E-09 1.58530E-09 1.58340E-09 1.58150E-09 1.57960E-09 1.57770E-09 1.57580E-09 1.57390E-09 &
     1.57200E-09 1.57010E-09 1.56820E-09 1.56630E-09 1.56440E-09 1.56260E-09 1.56070E-09 1.55880E-09 &
     1.55700E-09 1.55510E-09 1.55320E-09 1.55140E-09 1.54950E-09 1.54770E-09 1.54580E-09 1.54400E-09 &
     1.54210E-09 1.54030E-09 1.53850E-09 1.53660E-09 1.53480E-09 1.53300E-09 1.53110E-09 1.52930E-09 &
     1.52750E-09 1.52570E-09 1.52390E-09 1.52210E-09 1.52030E-09 1.51850E-09 1.51670E-09 1.51490E-09 &
     1.51310E-09 1.51130E-09 1.50950E-09 1.50770E-09 1.50590E-09 1.50420E-09 1.50240E-09 1.50060E-09 &
     1.49880E-09 1.49710E-09 1.49530E-09 1.49350E-09 1.49180E-09 1.49000E-09 1.48830E-09 1.48650E-09 &
     1.48480E-09 1.48300E-09 1.48130E-09 1.47960E-09 1.47780E-09 1.47610E-09 1.47440E-09 1.47260E-09 &
     1.47090E-09 1.46920E-09 1.46750E-09 1.46580E-09 1.46400E-09 1.46230E-09 1.46060E-09 1.45890E-09 &
     1.45720E-09 1.45550E-09 1.45380E-09 1.45210E-09 1.45040E-09 1.44880E-09 1.44710E-09 1.44540E-09 &
     1.44370E-09 1.44200E-09 1.44040E-09 1.43870E-09 1.43700E-09 1.43530E-09 1.43370E-09 1.43200E-09 &
     1.43040E-09 1.42870E-09 1.42710E-09 1.42540E-09 1.42380E-09 1.42210E-09 1.42050E-09 1.41880E-09 &
     1.41720E-09 1.41560E-09 1.41390E-09 1.41230E-09 1.41070E-09 1.40900E-09 1.40740E-09 1.40580E-09 &
     1.40420E-09 1.40260E-09 1.40100E-09 1.39940E-09 1.39780E-09 1.39620E-09 1.39460E-09 1.39300E-09 &
     1.39140E-09 1.38980E-09 1.38820E-09 1.38660E-09 1.38500E-09 1.38340E-09 1.38180E-09 1.38030E-09 &
     1.37870E-09 1.37710E-09 1.37550E-09 1.37400E-09 1.37240E-09 1.37080E-09 1.36930E-09 1.36770E-09 &
     1.36620E-09 1.36460E-09 1.36310E-09 1.36150E-09 1.36000E-09 1.35840E-09 1.35690E-09 1.35530E-09 &
     1.35380E-09 1.35230E-09 1.35070E-09 1.34920E-09 1.34770E-09 1.34620E-09 1.34460E-09 1.34310E-09 &
     1.34160E-09 1.34010E-09 1.33860E-09 1.33710E-09 1.33560E-09 1.33410E-09 1.33260E-09 1.33100E-09 &
     1.32960E-09 1.32810E-09 1.32660E-09 1.32510E-09 1.32360E-09 1.32210E-09 1.32060E-09 1.31910E-09 &
     1.31760E-09 1.31620E-09 1.31470E-09 1.31320E-09 1.31170E-09 1.31030E-09 1.30880E-09 1.30730E-09 &
     1.30590E-09 1.30440E-09 1.30300E-09 1.30150E-09 1.30010E-09 1.29860E-09 1.29720E-09 1.29570E-09 &
     1.29430E-09 1.29280E-09 1.29140E-09 1.29000E-09 1.28850E-09 1.28710E-09 1.28570E-09 1.28420E-09 &
     1.28280E-09 1.28140E-09 1.28000E-09 1.27850E-09 1.27710E-09 1.27570E-09 1.27430E-09 1.27290E-09 &
     1.27150E-09 1.27010E-09 1.26870E-09 1.26720E-09 1.26580E-09 1.26440E-09 1.26310E-09 1.26170E-09 &
     1.26030E-09 1.25890E-09 1.25750E-09 1.25610E-09 1.25470E-09 1.25330E-09 1.25200E-09 1.25060E-09 &
     1.24920E-09 1.24780E-09 1.24650E-09 1.24510E-09 1.24370E-09 1.24230E-09 1.24100E-09 1.23960E-09 &
     1.23830E-09 1.23690E-09 1.23560E-09 1.23420E-09 1.23280E-09 1.23150E-09 1.23010E-09 1.22880E-09 &
     1.22750E-09 1.22610E-09 1.22480E-09 1.22340E-09 1.22210E-09 1.22080E-09 1.21940E-09 1.21810E-09 &
     1.21680E-09 1.21550E-09 1.21410E-09 1.21280E-09 1.21150E-09 1.21020E-09 1.20890E-09 1.20750E-09 &
     1.20620E-09 1.20490E-09 1.20360E-09 1.20230E-09 1.20100E-09 1.19970E-09 1.19840E-09 1.19710E-09 &
     1.19580E-09 1.19450E-09 1.19320E-09 1.19190E-09 1.19060E-09 1.18930E-09 1.18810E-09 1.18680E-09 &
     1.18550E-09 1.18420E-09 1.18290E-09 1.18170E-09 1.18040E-09 1.17910E-09 1.17780E-09 1.17660E-09 &
     1.17530E-09 1.17400E-09 1.17280E-09 1.17150E-09 1.17030E-09 1.16900E-09 1.16770E-09 1.16650E-09 &
     1.16520E-09 1.16400E-09 1.16270E-09 1.16150E-09 1.16020E-09 1.15900E-09 1.15780E-09 1.15650E-09 &
     1.15530E-09 1.15400E-09 1.15280E-09 1.15160E-09 1.15040E-09 1.14910E-09 1.14790E-09 1.14670E-09 &
     1.14540E-09 1.14420E-09 1.14300E-09 1.14180E-09 1.14060E-09 1.13940E-09 1.13810E-09 1.13690E-09 &
     1.13570E-09 1.13450E-09 1.13330E-09 1.13210E-09 1.13090E-09 1.12970E-09 1.12850E-09 1.12730E-09 &
     1.12610E-09 1.12490E-09 1.12370E-09 1.12250E-09 1.12130E-09 1.12010E-09 1.11900E-09 1.11780E-09 &
     1.11660E-09 1.11540E-09 1.11420E-09 1.11300E-09 1.11190E-09 1.11070E-09 1.10950E-09 1.10840E-09 &
     1.10720E-09 1.10600E-09 1.10480E-09 1.10370E-09 1.10250E-09 1.10140E-09 1.10020E-09 1.09900E-09 &
     1.09790E-09 1.09670E-09 1.09560E-09 1.09440E-09 1.09330E-09 1.09210E-09 1.09100E-09 1.08980E-09 &
     1.08870E-09 1.08760E-09 1.08640E-09 1.08530E-09 1.08410E-09 1.08300E-09 1.08190E-09 1.08070E-09 &
     1.07960E-09 1.07850E-09 1.07730E-09 1.07620E-09 1.07510E-09 1.07400E-09 1.07290E-09 1.07170E-09 &
     1.07060E-09 1.06950E-09 1.06840E-09 1.06730E-09 1.06620E-09 1.06500E-09 1.06390E-09 1.06280E-09 &
     1.06170E-09 1.06060E-09 1.05950E-09 1.05840E-09 1.05730E-09 1.05620E-09 1.05510E-09 1.05400E-09 &
     1.05290E-09 1.05180E-09 1.05070E-09 1.04970E-09 1.04860E-09 1.04750E-09 1.04640E-09 1.04530E-09 &
     1.04420E-09 1.04320E-09 1.04210E-09 1.04100E-09 1.03990E-09 1.03880E-09 1.03780E-09 1.03670E-09 &
     1.03560E-09 1.03460E-09 1.03350E-09 1.03240E-09 1.03140E-09 1.03030E-09 1.02920E-09 1.02820E-09 &
     1.02710E-09 1.02610E-09 1.02500E-09 1.02400E-09 1.02290E-09 1.02190E-09 1.02080E-09 1.01980E-09 &
     1.01870E-09 1.01770E-09 1.01660E-09 1.01560E-09 1.01450E-09 1.01350E-09 1.01250E-09 1.01140E-09 &
     1.01040E-09 1.00940E-09 1.00830E-09 1.00730E-09 1.00630E-09 1.00520E-09 1.00420E-09 1.00320E-09 &
     1.00220E-09 1.00110E-09 1.00010E-09 9.99090E-10 9.98070E-10 9.97060E-10 9.96040E-10 9.95030E-10 &
     9.94010E-10 9.93000E-10 9.91990E-10 9.90980E-10 9.89980E-10 9.88970E-10 9.87960E-10 9.86960E-10 &
     9.85960E-10 9.84960E-10 9.83960E-10 9.82960E-10 9.81960E-10 9.80970E-10 9.79970E-10 9.78980E-10 &
     9.77990E-10 9.77000E-10 9.76010E-10 9.75020E-10 9.74040E-10 9.73050E-10 9.72070E-10 9.71090E-10 &
     9.70110E-10 9.69130E-10 9.68150E-10 9.67170E-10 9.66200E-10 9.65220E-10 9.64250E-10 9.63280E-10 &
     9.62310E-10 9.61340E-10 9.60370E-10 9.59400E-10 9.58440E-10 9.57470E-10 9.56510E-10 9.55550E-10 &
     9.54590E-10 9.53630E-10 9.52680E-10 9.51720E-10 9.50760E-10 9.49810E-10 9.48860E-10 9.47910E-10 &
     9.46960E-10 9.46010E-10 9.45060E-10 9.44120E-10 9.43170E-10 9.42230E-10 9.41290E-10 9.40350E-10 &
     9.39410E-10 9.38470E-10 9.37530E-10 9.36590E-10 9.35660E-10 9.34730E-10 9.33790E-10 9.32860E-10 &
     9.31930E-10 9.31010E-10 9.30080E-10 9.29150E-10 9.28230E-10 9.27300E-10 9.26380E-10 9.25460E-10 &
     9.24540E-10 9.23620E-10 9.22700E-10 9.21790E-10 9.20870E-10 9.19960E-10 9.19050E-10 9.18140E-10 &
     9.17230E-10 9.16320E-10 9.15410E-10 9.14500E-10 9.13600E-10 9.12690E-10 9.11790E-10 9.10890E-10 &
     9.09990E-10 9.09090E-10 9.08190E-10 9.07290E-10 9.06400E-10 9.05500E-10 9.04610E-10 9.03720E-10 &
     9.02830E-10 9.01940E-10 9.01050E-10 9.00160E-10 8.99270E-10 8.98390E-10 8.97500E-10 8.96620E-10 &
     8.95740E-10 8.94860E-10 8.93980E-10 8.93100E-10 8.92220E-10 8.91350E-10 8.90470E-10 8.89600E-10 &
     8.88730E-10 8.87850E-10 8.86980E-10 8.86120E-10 8.85250E-10 8.84380E-10 8.83510E-10 8.82650E-10 &
     8.81790E-10 8.80920E-10 8.80060E-10 8.79200E-10 8.78340E-10 8.77490E-10 8.76630E-10 8.75770E-10 &
     8.74920E-10 8.74070E-10 8.73210E-10 8.72360E-10 8.71510E-10 8.70660E-10 8.69820E-10 8.68970E-10 &
     8.68120E-10 8.67280E-10 8.66440E-10 8.65590E-10 8.64750E-10 8.63910E-10 8.63070E-10 8.62230E-10 &
     8.61400E-10 8.60560E-10 8.59730E-10 8.58890E-10 8.58060E-10 8.57230E-10 8.56400E-10 8.55570E-10 &
     8.54740E-10 8.53910E-10 8.53090E-10 8.52260E-10 8.51440E-10 8.50620E-10 8.49800E-10 8.48970E-10 &
     8.48150E-10 8.47340E-10 8.46520E-10 8.45700E-10 8.44890E-10 8.44070E-10 8.43260E-10 8.42450E-10 &
     8.41630E-10 8.40820E-10 8.40010E-10 8.39210E-10 8.38400E-10 8.37590E-10 8.36790E-10 8.35980E-10 &
     8.35180E-10 8.34380E-10 8.33580E-10 8.32780E-10 8.31980E-10 8.31180E-10 8.30380E-10 8.29590E-10 &
     8.28790E-10 8.28000E-10 8.27210E-10 8.26420E-10 8.25620E-10 8.24830E-10 8.24050E-10 8.23260E-10 &
     8.22470E-10 8.21690E-10 8.20900E-10 8.20120E-10 8.19330E-10 8.18550E-10 8.17770E-10 8.16990E-10 &
     8.16210E-10 8.15440E-10 8.14660E-10 8.13880E-10 8.13110E-10 8.12330E-10 8.11560E-10 8.10790E-10 &
     8.10020E-10 8.09250E-10 8.08480E-10 8.07710E-10 8.06950E-10 8.06180E-10 8.05410E-10 8.04650E-10 &
     8.03890E-10 8.03120E-10 8.02360E-10 8.01600E-10 8.00840E-10 8.00090E-10 7.99330E-10 7.98570E-10 &
     7.97820E-10 7.97060E-10 7.96310E-10 7.95560E-10 7.94810E-10 7.94050E-10 7.93300E-10 7.92560E-10 &
     7.91810E-10 7.91060E-10 7.90320E-10 7.89570E-10 7.88830E-10 7.88080E-10 7.87340E-10 7.86600E-10 &
     7.85860E-10 7.85120E-10 7.84380E-10 7.83640E-10 7.82910E-10 7.82170E-10 7.81440E-10 7.80700E-10 &
     7.79970E-10 7.79240E-10 7.78510E-10 7.77780E-10 7.77050E-10 7.76320E-10 7.75590E-10 7.74860E-10 &
     7.74140E-10 7.73410E-10 7.72690E-10 7.71970E-10 7.71240E-10 7.70520E-10 7.69800E-10 7.69080E-10 &
     7.68360E-10 7.67650E-10 7.66930E-10 7.66210E-10 7.65500E-10 7.64780E-10 7.64070E-10 7.63360E-10 &
     7.62650E-10 7.61940E-10 7.61230E-10 7.60520E-10 7.59810E-10 7.59100E-10 7.58400E-10 7.57690E-10 &
     7.56990E-10 7.56280E-10 7.55580E-10 7.54880E-10 7.54180E-10 7.53480E-10 7.52780E-10 7.52080E-10 &
     7.51380E-10 7.50690E-10 7.49990E-10 7.49300E-10 7.48600E-10 7.47910E-10 7.47220E-10 7.46520E-10 &
     7.45830E-10 7.45140E-10 7.44460E-10 7.43770E-10 7.43080E-10 7.42390E-10 7.41710E-10 7.41020E-10 &
     7.40340E-10 7.39660E-10 7.38970E-10 7.38290E-10 7.37610E-10 7.36930E-10 7.36250E-10 7.35570E-10 &
     7.34900E-10 7.34220E-10 7.33540E-10 7.32870E-10 7.32200E-10 7.31520E-10 7.30850E-10 7.30180E-10 &
     7.29510E-10 7.28840E-10 7.28170E-10 7.27500E-10 7.26830E-10 7.26160E-10 7.25500E-10 7.24830E-10 &
     7.24170E-10 7.23510E-10 7.22840E-10 7.22180E-10 7.21520E-10 7.20860E-10 7.20200E-10 7.19540E-10 &
     7.18880E-10 7.18230E-10 7.17570E-10 7.16910E-10 7.16260E-10 7.15600E-10 7.14950E-10 7.14300E-10 &
     7.13650E-10 7.13000E-10 7.12350E-10 7.11700E-10 7.11050E-10 7.10400E-10 7.09750E-10 7.09110E-10 &
     7.08460E-10 7.07820E-10 7.07170E-10 7.06530E-10 7.05890E-10 7.05250E-10 7.04610E-10 7.03970E-10 &
     7.03330E-10 7.02690E-10 7.02050E-10 7.01410E-10 7.00780E-10 7.00140E-10 6.99510E-10 6.98870E-10 &
     6.98240E-10 6.97610E-10 6.96970E-10 6.96340E-10 6.95710E-10 6.95080E-10 6.94450E-10 6.93830E-10 &
     6.93200E-10 6.92570E-10 6.91950E-10 6.91320E-10 6.90700E-10 6.90070E-10 6.89450E-10 6.88830E-10 &
     6.88210E-10 6.87590E-10 6.86970E-10 6.86350E-10 6.85730E-10 6.85110E-10 6.84490E-10 6.83880E-10 &
     6.83260E-10 6.82650E-10 6.82030E-10 6.81420E-10 6.80810E-10 6.80200E-10 6.79580E-10 6.78970E-10 &
     6.78360E-10 6.77750E-10 6.77150E-10 6.76540E-10 6.75930E-10 6.75330E-10 6.74720E-10 6.74110E-10 &
     6.73510E-10 6.72910E-10 6.72300E-10 6.71700E-10 6.71100E-10 6.70500E-10 6.69900E-10 6.69300E-10 &
     6.68700E-10 6.68100E-10 6.67510E-10 6.66910E-10 6.66320E-10 6.65720E-10 6.65130E-10 6.64530E-10 &
     6.63940E-10 6.63350E-10 6.62750E-10 6.62160E-10 6.61570E-10 6.60980E-10 6.60400E-10 6.59810E-10 &
     6.59220E-10 6.58630E-10 6.58050E-10 6.57460E-10 6.56880E-10 6.56290E-10 6.55710E-10 6.55130E-10 &
     6.54540E-10 6.53960E-10 6.53380E-10 6.52800E-10 6.52220E-10 6.51640E-10 6.51060E-10 6.50490E-10 &
     6.49910E-10 6.49330E-10 6.48760E-10 6.48180E-10 6.47610E-10 6.47040E-10 6.46460E-10 6.45890E-10 &
     6.45320E-10 6.44750E-10 6.44180E-10 6.43610E-10 6.43040E-10 6.42470E-10 6.41910E-10 6.41340E-10 &
     6.40770E-10 6.40210E-10 6.39640E-10 6.39080E-10 6.38510E-10 6.37950E-10 6.37390E-10 6.36830E-10 &
     6.36260E-10 6.35700E-10 6.35140E-10 6.34580E-10 6.34030E-10 6.33470E-10 6.32910E-10 6.32350E-10 &
     6.31800E-10 6.31240E-10 6.30690E-10 6.30130E-10 6.29580E-10 6.29030E-10 6.28470E-10 6.27920E-10 &
     6.27370E-10 6.26820E-10 6.26270E-10 6.25720E-10 6.25170E-10 6.24620E-10 6.24080E-10 6.23530E-10 &
     6.22980E-10 6.22440E-10 6.21890E-10 6.21350E-10 6.20810E-10 6.20260E-10 6.19720E-10 6.19180E-10 &
     6.18640E-10 6.18100E-10 6.17560E-10 6.17020E-10 6.16480E-10 6.15940E-10 6.15400E-10 6.14870E-10 &
     6.14330E-10 6.13790E-10 6.13260E-10 6.12720E-10 6.12190E-10 6.11660E-10 6.11120E-10 6.10590E-10 &
     6.10060E-10 6.09530E-10 6.09000E-10 6.08470E-10 6.07940E-10 6.07410E-10 6.06880E-10 6.06360E-10 &
     6.05830E-10 6.05300E-10 6.04780E-10 6.04250E-10 6.03730E-10 6.03200E-10 6.02680E-10 6.02160E-10 &
     6.01640E-10 6.01110E-10 6.00590E-10 6.00070E-10 5.99550E-10 5.99030E-10 5.98510E-10 5.98000E-10 &
     5.97480E-10 5.96960E-10 5.96450E-10 5.95930E-10 5.95410E-10 5.94900E-10 5.94390E-10 5.93870E-10 &
     5.93360E-10 5.92850E-10 5.92330E-10 5.91820E-10 5.91310E-10 5.90800E-10 5.90290E-10 5.89780E-10 &
     5.89280E-10 5.88770E-10 5.88260E-10 5.87750E-10 5.87250E-10 5.86740E-10 5.86240E-10 5.85730E-10 &
     5.85230E-10 5.84720E-10 5.84220E-10 5.83720E-10 5.83220E-10 5.82720E-10 5.82220E-10 5.81720E-10 &
     5.81220E-10 5.80720E-10 5.80220E-10 5.79720E-10 5.79220E-10 5.78730E-10 5.78230E-10 5.77730E-10 &
     5.77240E-10 5.76740E-10 5.76250E-10 5.75760E-10 5.75260E-10 5.74770E-10 5.74280E-10 5.73790E-10 &
     5.73300E-10 5.72810E-10 5.72320E-10 5.71830E-10 5.71340E-10 5.70850E-10 5.70360E-10 5.69880E-10 &
     5.69390E-10 5.68900E-10 5.68420E-10 5.67930E-10 5.67450E-10 5.66960E-10 5.66480E-10 5.66000E-10 &
     5.65520E-10 5.65030E-10 5.64550E-10 5.64070E-10 5.63590E-10 5.63110E-10 5.62630E-10 5.62150E-10 &
     5.61680E-10 5.61200E-10 5.60720E-10 5.60240E-10 5.59770E-10 5.59290E-10 5.58820E-10 5.58340E-10 &
     5.57870E-10 5.57400E-10 5.56920E-10 5.56450E-10 5.55980E-10 5.55510E-10 5.55040E-10 5.54560E-10 &
     5.54090E-10 5.53630E-10 5.53160E-10 5.52690E-10 5.52220E-10 5.51750E-10 5.51290E-10 5.50820E-10 &
     5.50350E-10 5.49890E-10 5.49420E-10 5.48960E-10 5.48490E-10 5.48030E-10 5.47570E-10 5.47110E-10 &
     5.46640E-10 5.46180E-10 5.45720E-10 5.45260E-10 5.44800E-10 5.44340E-10 5.43880E-10 5.43420E-10 &
     5.42970E-10 5.42510E-10 5.42050E-10 5.41600E-10 5.41140E-10 5.40680E-10 5.40230E-10 5.39770E-10 &
     5.39320E-10 5.38870E-10 5.38410E-10 5.37960E-10 5.37510E-10 5.37060E-10 5.36610E-10 5.36160E-10 &
     5.35710E-10 5.35260E-10 5.34810E-10 5.34360E-10 5.33910E-10 5.33460E-10 5.33010E-10 5.32570E-10 &
     5.32120E-10 5.31680E-10 5.31230E-10 5.30790E-10 5.30340E-10 5.29900E-10 5.29450E-10 5.29010E-10 &
     5.28570E-10 5.28130E-10 5.27680E-10 5.27240E-10 5.26800E-10 5.26360E-10 5.25920E-10 5.25480E-10 &
     5.25050E-10 5.24610E-10 5.24170E-10 5.23730E-10 5.23290E-10 5.22860E-10 5.22420E-10 5.21990E-10 &
     5.21550E-10 5.21120E-10 5.20680E-10 5.20250E-10 5.19820E-10 5.19380E-10 5.18950E-10 5.18520E-10 &
     5.18090E-10 5.17660E-10 5.17230E-10 5.16800E-10 5.16370E-10 5.15940E-10 5.15510E-10 5.15080E-10 &
     5.14650E-10 5.14230E-10 5.13800E-10 5.13370E-10 5.12950E-10 5.12520E-10 5.12100E-10 5.11670E-10 &
     5.11250E-10 5.10820E-10 5.10400E-10 5.09980E-10 5.09550E-10 5.09130E-10 5.08710E-10 5.08290E-10 &
     5.07870E-10 5.07450E-10 5.07030E-10 5.06610E-10 5.06190E-10 5.05770E-10 5.05350E-10 5.04940E-10 &
     5.04520E-10 5.04100E-10 5.03690E-10 5.03270E-10 5.02860E-10 5.02440E-10 5.02030E-10 5.01610E-10 &
     5.01200E-10 5.00790E-10 5.00370E-10 4.99960E-10 4.99550E-10 4.99140E-10 4.98730E-10 4.98320E-10 &
     4.97910E-10 4.97500E-10 4.97090E-10 4.96680E-10 4.96270E-10 4.95860E-10 4.95450E-10 4.95050E-10 &
     4.94640E-10 4.94230E-10 4.93830E-10 4.93420E-10 4.93020E-10 4.92610E-10 4.92210E-10 4.91800E-10 &
     4.91400E-10 4.91000E-10 4.90600E-10 4.90190E-10 4.89790E-10 4.89390E-10 4.88990E-10 4.88590E-10 &
     4.88190E-10 4.87790E-10 4.87390E-10 4.86990E-10 4.86590E-10 4.86200E-10 4.85800E-10 4.85400E-10 &
     4.85000E-10 4.84610E-10 4.84210E-10 4.83820E-10 4.83420E-10 4.83030E-10 4.82630E-10 4.82240E-10 &
     4.81840E-10 4.81450E-10 4.81060E-10 4.80670E-10 4.80270E-10 4.79880E-10 4.79490E-10 4.79100E-10 &
     4.78710E-10 4.78320E-10 4.77930E-10 4.77540E-10 4.77150E-10 4.76770E-10 4.76380E-10 4.75990E-10 &
     4.75600E-10 4.75220E-10 4.74830E-10 4.74450E-10 4.74060E-10 4.73670E-10 4.73290E-10 4.72910E-10 &
     4.72520E-10 4.72140E-10 4.71760E-10 4.71370E-10 4.70990E-10 4.70610E-10 4.70230E-10 4.69850E-10 &
     4.69470E-10 4.69090E-10 4.68710E-10 4.68330E-10 4.67950E-10 4.67570E-10 4.67190E-10 4.66810E-10 &
     4.66430E-10 4.66060E-10 4.65680E-10 4.65300E-10 4.64930E-10 4.64550E-10 4.64180E-10 4.63800E-10 &
     4.63430E-10 4.63050E-10 4.62680E-10 4.62310E-10 4.61930E-10 4.61560E-10 4.61190E-10 4.60820E-10 &
     4.60440E-10 4.60070E-10 4.59700E-10 4.59330E-10 4.58960E-10 4.58590E-10 4.58220E-10 4.57850E-10 &
     4.57480E-10 4.57120E-10 4.56750E-10 4.56380E-10 4.56010E-10 4.55650E-10 4.55280E-10 4.54920E-10 &
     4.54550E-10 4.54180E-10 4.53820E-10 4.53460E-10 4.53090E-10 4.52730E-10 4.52360E-10 4.52000E-10 &
     4.51640E-10 4.51280E-10 4.50910E-10 4.50550E-10 4.50190E-10 4.49830E-10 4.49470E-10 4.49110E-10 &
     4.48750E-10 4.48390E-10 4.48030E-10 4.47670E-10 4.47320E-10 4.46960E-10 4.46600E-10 4.46240E-10 &
     4.45890E-10 4.45530E-10 4.45170E-10 4.44820E-10 4.44460E-10 4.44110E-10 4.43750E-10 4.43400E-10 &
     4.43040E-10 4.42690E-10 4.42340E-10 4.41980E-10 4.41630E-10 4.41280E-10 4.40930E-10 4.40580E-10 &
     4.40220E-10 4.39870E-10 4.39520E-10 4.39170E-10 4.38820E-10 4.38470E-10 4.38120E-10 4.37780E-10 &
     4.37430E-10 4.37080E-10 4.36730E-10 4.36380E-10 4.36040E-10 4.35690E-10 4.35340E-10 4.35000E-10 &
     4.34650E-10 4.34310E-10 4.33960E-10 4.33620E-10 4.33270E-10 4.32930E-10 4.32590E-10 4.32240E-10 &
     4.31900E-10 4.31560E-10 4.31220E-10 4.30870E-10 4.30530E-10 4.30190E-10 4.29850E-10 4.29510E-10 &
     4.29170E-10 4.28830E-10 4.28490E-10 4.28150E-10 4.27810E-10 4.27470E-10 4.27130E-10 4.26800E-10 &
     4.26460E-10 4.26120E-10 4.25790E-10 4.25450E-10 4.25110E-10 4.24780E-10 4.24440E-10 4.24110E-10 &
     4.23770E-10 4.23440E-10 4.23100E-10 4.22770E-10 4.22430E-10 4.22100E-10 4.21770E-10 4.21440E-10 &
     4.21100E-10 4.20770E-10 4.20440E-10 4.20110E-10 4.19780E-10 4.19450E-10 4.19120E-10 4.18790E-10 &
     4.18460E-10 4.18130E-10 4.17800E-10 4.17470E-10 4.17140E-10 4.16810E-10 4.16490E-10 4.16160E-10 &
     4.15830E-10 4.15500E-10 4.15180E-10 4.14850E-10 4.14530E-10 4.14200E-10 4.13880E-10 4.13550E-10 &
     4.13230E-10 4.12900E-10 4.12580E-10 4.12250E-10 4.11930E-10 4.11610E-10 4.11290E-10 4.10960E-10 &
     4.10640E-10 4.10320E-10 4.10000E-10 4.09680E-10 4.09360E-10 4.09040E-10 4.08720E-10 4.08400E-10 &
     4.08080E-10 4.07760E-10 4.07440E-10 4.07120E-10 4.06800E-10 4.06480E-10 4.06160E-10 4.05850E-10 &
     4.05530E-10 4.05210E-10 4.04900E-10 4.04580E-10 4.04260E-10 4.03950E-10 4.03630E-10 4.03320E-10 &
     4.03000E-10 4.02690E-10 4.02380E-10 4.02060E-10 4.01750E-10 4.01440E-10 4.01120E-10 4.00810E-10 &
     4.00500E-10 4.00190E-10 3.99870E-10 3.99560E-10 3.99250E-10 3.98940E-10 3.98630E-10 3.98320E-10 &
     3.98010E-10 3.97700E-10 3.97390E-10 3.97080E-10 3.96770E-10 3.96470E-10 3.96160E-10 3.95850E-10 &
     3.95540E-10 3.95230E-10 3.94930E-10 3.94620E-10 3.94310E-10 3.94010E-10 3.93700E-10 3.93400E-10 &
     3.93090E-10 3.92790E-10 3.92480E-10 3.92180E-10 3.91870E-10 3.91570E-10 3.91270E-10 3.90960E-10 &
     3.90660E-10 3.90360E-10 3.90060E-10 3.89750E-10 3.89450E-10 3.89150E-10 3.88850E-10 3.88550E-10 &
     3.88250E-10 3.87950E-10 3.87650E-10 3.87350E-10 3.87050E-10 3.86750E-10 3.86450E-10 3.86150E-10 &
     3.85850E-10 3.85560E-10 3.85260E-10 3.84960E-10 3.84660E-10 3.84370E-10 3.84070E-10 3.83770E-10 &
     3.83480E-10 3.83180E-10 3.82890E-10 3.82590E-10 3.82300E-10 3.82000E-10 3.81710E-10 3.81410E-10 &
     3.81120E-10 3.80830E-10 3.80530E-10 3.80240E-10 3.79950E-10 3.79660E-10 3.79360E-10 3.79070E-10 &
     3.78780E-10 3.78490E-10 3.78200E-10 3.77910E-10 3.77620E-10 3.77330E-10 3.77040E-10 3.76750E-10 &
     3.76460E-10 3.76170E-10 3.75880E-10 3.75590E-10 3.75300E-10 3.75010E-10 3.74730E-10 3.74440E-10 &
     3.74150E-10 3.73860E-10 3.73580E-10 3.73290E-10 3.73000E-10 3.72720E-10 3.72430E-10 3.72150E-10 &
     3.71860E-10 3.71580E-10 3.71290E-10 3.71010E-10 3.70720E-10 3.70440E-10 3.70160E-10 3.69870E-10 &
     3.69590E-10 3.69310E-10 3.69030E-10 3.68740E-10 3.68460E-10 3.68180E-10 3.67900E-10 3.67620E-10 &
     3.67340E-10 3.67060E-10 3.66780E-10 3.66500E-10 3.66220E-10 3.65940E-10 3.65660E-10 3.65380E-10 &
     3.65100E-10 3.64820E-10 3.64540E-10 3.64260E-10 3.63990E-10 3.63710E-10 3.63430E-10 3.63150E-10 &
     3.62880E-10 3.62600E-10 3.62320E-10 3.62050E-10 3.61770E-10 3.61500E-10 3.61220E-10 3.60950E-10 &
     3.60670E-10 3.60400E-10 3.60120E-10 3.59850E-10 3.59580E-10 3.59300E-10 3.59030E-10 3.58760E-10 &
     3.58480E-10 3.58210E-10 3.57940E-10 3.57670E-10 3.57400E-10 3.57120E-10 3.56850E-10 3.56580E-10 &
     3.56310E-10 3.56040E-10 3.55770E-10 3.55500E-10 3.55230E-10 3.54960E-10 3.54690E-10 3.54420E-10 &
     3.54150E-10 3.53890E-10 3.53620E-10 3.53350E-10 3.53080E-10 3.52810E-10 3.52550E-10 3.52280E-10 &
     3.52010E-10 3.51750E-10 3.51480E-10 3.51210E-10 3.50950E-10 3.50680E-10 3.50420E-10 3.50150E-10 &
     3.49890E-10 3.49620E-10 3.49360E-10 3.49100E-10 3.48830E-10 3.48570E-10 3.48300E-10 3.48040E-10 &
     3.47780E-10 3.47520E-10 3.47250E-10 3.46990E-10 3.46730E-10 3.46470E-10 3.46210E-10 3.45950E-10 &
     3.45680E-10 3.45420E-10 3.45160E-10 3.44900E-10 3.44640E-10 3.44380E-10 3.44120E-10 3.43860E-10 &
     3.43610E-10 3.43350E-10 3.43090E-10 3.42830E-10 3.42570E-10 3.42310E-10 3.42060E-10 3.41800E-10 &
     3.41540E-10 3.41280E-10 3.41030E-10 3.40770E-10 3.40520E-10 3.40260E-10 3.40000E-10 3.39750E-10 &
     3.39490E-10 3.39240E-10 3.38980E-10 3.38730E-10 3.38470E-10 3.38220E-10 3.37970E-10 3.37710E-10 &
     3.37460E-10 3.37200E-10 3.36950E-10 3.36700E-10 3.36450E-10 3.36190E-10 3.35940E-10 3.35690E-10 &
     3.35440E-10 3.35190E-10 3.34940E-10 3.34680E-10 3.34430E-10 3.34180E-10 3.33930E-10 3.33680E-10 &
     3.33430E-10 3.33180E-10 3.32930E-10 3.32690E-10 3.32440E-10 3.32190E-10 3.31940E-10 3.31690E-10 &
     3.31440E-10 3.31200E-10 3.30950E-10 3.30700E-10 3.30450E-10 3.30210E-10 3.29960E-10 3.29710E-10 &
     3.29470E-10 3.29220E-10 3.28970E-10 3.28730E-10 3.28480E-10 3.28240E-10 3.27990E-10 3.27750E-10 &
     3.27500E-10 3.27260E-10 3.27020E-10 3.26770E-10 3.26530E-10 3.26290E-10 3.26040E-10 3.25800E-10 &
     3.25560E-10 3.25310E-10 3.25070E-10 3.24830E-10 3.24590E-10 3.24350E-10 3.24100E-10 3.23860E-10 &
     3.23620E-10 3.23380E-10 3.23140E-10 3.22900E-10 3.22660E-10 3.22420E-10 3.22180E-10 3.21940E-10 &
     3.21700E-10 3.21460E-10 3.21220E-10 3.20980E-10 3.20750E-10 3.20510E-10 3.20270E-10 3.20030E-10 &
     3.19790E-10 3.19560E-10 3.19320E-10 3.19080E-10 3.18850E-10 3.18610E-10 3.18370E-10 3.18140E-10 &
     3.17900E-10 3.17660E-10 3.17430E-10 3.17190E-10 3.16960E-10 3.16720E-10 3.16490E-10 3.16250E-10 &
     3.16020E-10 3.15790E-10 3.15550E-10 3.15320E-10 3.15090E-10 3.14850E-10 3.14620E-10 3.14390E-10 &
     3.14150E-10 3.13920E-10 3.13690E-10 3.13460E-10 3.13230E-10 3.12990E-10 3.12760E-10 3.12530E-10 &
     3.12300E-10 3.12070E-10 3.11840E-10 3.11610E-10 3.11380E-10 3.11150E-10 3.10920E-10 3.10690E-10 &
     3.10460E-10 3.10230E-10 3.10000E-10 3.09770E-10 3.09540E-10 3.09320E-10 3.09090E-10 3.08860E-10 &
     3.08630E-10 3.08400E-10 3.08180E-10 3.07950E-10 3.07720E-10 3.07500E-10 3.07270E-10 3.07040E-10 &
     3.06820E-10 3.06590E-10 3.06370E-10 3.06140E-10 3.05910E-10 3.05690E-10 3.05460E-10 3.05240E-10 &
     3.05020E-10 3.04790E-10 3.04570E-10 3.04340E-10 3.04120E-10 3.03900E-10 3.03670E-10 3.03450E-10 &
     3.03230E-10 3.03000E-10 3.02780E-10 3.02560E-10 3.02340E-10 3.02110E-10 3.01890E-10 3.01670E-10 &
     3.01450E-10 3.01230E-10 3.01010E-10 3.00790E-10 3.00560E-10 3.00340E-10 3.00120E-10 2.99900E-10 &
     2.99680E-10 2.99460E-10 2.99240E-10 2.99030E-10 2.98810E-10 2.98590E-10 2.98370E-10 2.98150E-10 &
     2.97930E-10 2.97710E-10 2.97500E-10 2.97280E-10 2.97060E-10 2.96840E-10 2.96630E-10 2.96410E-10 &
     2.96190E-10 2.95970E-10 2.95760E-10 2.95540E-10 2.95330E-10 2.95110E-10 2.94890E-10 2.94680E-10 &
     2.94460E-10 2.94250E-10 2.94030E-10 2.93820E-10 2.93600E-10 2.93390E-10 2.93170E-10 2.92960E-10 &
     2.92750E-10 2.92530E-10 2.92320E-10 2.92110E-10 2.91890E-10 2.91680E-10 2.91470E-10 2.91250E-10 &
     2.91040E-10 2.90830E-10 2.90620E-10 2.90410E-10 2.90190E-10 2.89980E-10 2.89770E-10 2.89560E-10 &
     2.89350E-10 2.89140E-10 2.88930E-10 2.88720E-10 2.88510E-10 2.88300E-10 2.88090E-10 2.87880E-10 &
     2.87670E-10 2.87460E-10 2.87250E-10 2.87040E-10 2.86830E-10 2.86620E-10 2.86410E-10 2.86210E-10 &
     2.86000E-10 2.85790E-10 2.85580E-10 2.85380E-10 2.85170E-10 2.84960E-10 2.84750E-10 2.84550E-10 &
     2.84340E-10 2.84130E-10 2.83930E-10 2.83720E-10 2.83520E-10 2.83310E-10 2.83100E-10 2.82900E-10 &
     2.82690E-10 2.82490E-10 2.82280E-10 2.82080E-10 2.81870E-10 2.81670E-10 2.81470E-10 2.81260E-10 &
     2.81060E-10 2.80850E-10 2.80650E-10 2.80450E-10 2.80240E-10 2.80040E-10 2.79840E-10 2.79640E-10 &
     2.79430E-10 2.79230E-10 2.79030E-10 2.78830E-10 2.78620E-10 2.78420E-10 2.78220E-10 2.78020E-10 &
     2.77820E-10 2.77620E-10 2.77420E-10 2.77220E-10 2.77020E-10 2.76820E-10 2.76620E-10 2.76420E-10 &
     2.76220E-10 2.76020E-10 2.75820E-10 2.75620E-10 2.75420E-10 2.75220E-10 2.75020E-10 2.74820E-10 &
     2.74620E-10 2.74430E-10 2.74230E-10 2.74030E-10 2.73830E-10 2.73640E-10 2.73440E-10 2.73240E-10 &
     2.73040E-10 2.72850E-10 2.72650E-10 2.72450E-10 2.72260E-10 2.72060E-10 2.71860E-10 2.71670E-10 &
     2.71470E-10 2.71280E-10 2.71080E-10 2.70890E-10 2.70690E-10 2.70500E-10 2.70300E-10 2.70110E-10 &
     2.69910E-10 2.69720E-10 2.69520E-10 2.69330E-10 2.69140E-10 2.68940E-10 2.68750E-10 2.68560E-10 &
     2.68360E-10 2.68170E-10 2.67980E-10 2.67790E-10 2.67590E-10 2.67400E-10 2.67210E-10 2.67020E-10 &
     2.66820E-10 2.66630E-10 2.66440E-10 2.66250E-10 2.66060E-10 2.65870E-10 2.65680E-10 2.65490E-10 &
     2.65300E-10 2.65110E-10 2.64920E-10 2.64730E-10 2.64540E-10 2.64350E-10 2.64160E-10 2.63970E-10 &
     2.63780E-10 2.63590E-10 2.63400E-10 2.63210E-10 2.63020E-10 2.62830E-10 2.62650E-10 2.62460E-10 &
     2.62270E-10 2.62080E-10 2.61890E-10 2.61710E-10 2.61520E-10 2.61330E-10 2.61140E-10 2.60960E-10 &
     2.60770E-10 2.60580E-10 2.60400E-10 2.60210E-10 2.60030E-10 2.59840E-10 2.59650E-10 2.59470E-10 &
     2.59280E-10 2.59100E-10 2.58910E-10 2.58730E-10 2.58540E-10 2.58360E-10 2.58170E-10 2.57990E-10 &
     2.57800E-10 2.57620E-10 2.57440E-10 2.57250E-10 2.57070E-10 2.56880E-10 2.56700E-10 2.56520E-10 &
     2.56330E-10 2.56150E-10 2.55970E-10 2.55790E-10 2.55600E-10 2.55420E-10 2.55240E-10 2.55060E-10 &
     2.54880E-10 2.54690E-10 2.54510E-10 2.54330E-10 2.54150E-10 2.53970E-10 2.53790E-10 2.53610E-10 &
     2.53430E-10 2.53240E-10 2.53060E-10 2.52880E-10 2.52700E-10 2.52520E-10 2.52340E-10 2.52160E-10 &
     2.51990E-10 2.51810E-10 2.51630E-10 2.51450E-10 2.51270E-10 2.51090E-10 2.50910E-10 2.50730E-10 &
     2.50550E-10 2.50380E-10 2.50200E-10 2.50020E-10 2.49840E-10 2.49670E-10 2.49490E-10 2.49310E-10 &
     2.49130E-10 2.48960E-10 2.48780E-10 2.48600E-10 2.48430E-10 2.48250E-10 2.48070E-10 2.47900E-10 &
     2.47720E-10 2.47550E-10 2.47370E-10 2.47190E-10 2.47020E-10 2.46840E-10 2.46670E-10 2.46490E-10 &
     2.46320E-10 2.46140E-10 2.45970E-10 2.45790E-10 2.45620E-10 2.45450E-10 2.45270E-10 2.45100E-10 &
     2.44920E-10 2.44750E-10 2.44580E-10 2.44400E-10 2.44230E-10 2.44060E-10 2.43890E-10 2.43710E-10 &
     2.43540E-10 2.43370E-10 2.43200E-10 2.43020E-10 2.42850E-10 2.42680E-10 2.42510E-10 2.42340E-10 &
     2.42160E-10 2.41990E-10 2.41820E-10 2.41650E-10 2.41480E-10 2.41310E-10 2.41140E-10 2.40970E-10 &
     2.40800E-10 2.40630E-10 2.40460E-10 2.40290E-10 2.40120E-10 2.39950E-10 2.39780E-10 2.39610E-10 &
     2.39440E-10 2.39270E-10 2.39100E-10 2.38930E-10 2.38760E-10 2.38600E-10 2.38430E-10 2.38260E-10 &
     2.38090E-10 2.37920E-10 2.37750E-10 2.37590E-10 2.37420E-10 2.37250E-10 2.37080E-10 2.36920E-10 &
     2.36750E-10 2.36580E-10 2.36420E-10 2.36250E-10 2.36080E-10 2.35920E-10 2.35750E-10 2.35580E-10 &
     2.35420E-10 2.35250E-10 2.35090E-10 2.34920E-10 2.34760E-10 2.34590E-10 2.34430E-10 2.34260E-10 &
     2.34100E-10 2.33930E-10 2.33770E-10 2.33600E-10 2.33440E-10 2.33270E-10 2.33110E-10 2.32940E-10 &
     2.32780E-10 2.32620E-10 2.32450E-10 2.32290E-10 2.32130E-10 2.31960E-10 2.31800E-10 2.31640E-10 &
     2.31470E-10 2.31310E-10 2.31150E-10 2.30990E-10 2.30820E-10 2.30660E-10 2.30500E-10 2.30340E-10 &
     2.30180E-10 2.30010E-10 2.29850E-10 2.29690E-10 2.29530E-10 2.29370E-10 2.29210E-10 2.29050E-10 &
     2.28890E-10 2.28730E-10 2.28560E-10 2.28400E-10 2.28240E-10 2.28080E-10 2.27920E-10 2.27760E-10 &
     2.27600E-10 2.27440E-10 2.27290E-10 2.27130E-10 2.26970E-10 2.26810E-10 2.26650E-10 2.26490E-10 &
     2.26330E-10 2.26170E-10 2.26010E-10 2.25860E-10 2.25700E-10 2.25540E-10 2.25380E-10 2.25220E-10 &
     2.25070E-10 2.24910E-10 2.24750E-10 2.24590E-10 2.24440E-10 2.24280E-10 2.24120E-10 2.23970E-10 &
     2.23810E-10 2.23650E-10 2.23500E-10 2.23340E-10 2.23180E-10 2.23030E-10 2.22870E-10 2.22710E-10 &
     2.22560E-10 2.22400E-10 2.22250E-10 2.22090E-10 2.21940E-10 2.21780E-10 2.21630E-10 2.21470E-10 &
     2.21320E-10 2.21160E-10 2.21010E-10 2.20850E-10 2.20700E-10 2.20550E-10 2.20390E-10 2.20240E-10 &
     2.20080E-10 2.19930E-10 2.19780E-10 2.19620E-10 2.19470E-10 2.19320E-10 2.19160E-10 2.19010E-10 &
     2.18860E-10 2.18710E-10 2.18550E-10 2.18400E-10 2.18250E-10 2.18100E-10 2.17940E-10 2.17790E-10 &
     2.17640E-10 2.17490E-10 2.17340E-10 2.17180E-10 2.17030E-10 2.16880E-10 2.16730E-10 2.16580E-10 &
     2.16430E-10 2.16280E-10 2.16130E-10 2.15980E-10 2.15830E-10 2.15680E-10 2.15530E-10 2.15380E-10 &
     2.15230E-10 2.15080E-10 2.14930E-10 2.14780E-10 2.14630E-10 2.14480E-10 2.14330E-10 2.14180E-10 &
     2.14030E-10 2.13880E-10 2.13730E-10 2.13580E-10 2.13440E-10 2.13290E-10 2.13140E-10 2.12990E-10 &
     2.12840E-10 2.12690E-10 2.12550E-10 2.12400E-10 2.12250E-10 2.12100E-10 2.11960E-10 2.11810E-10 &
     2.11660E-10 2.11520E-10 2.11370E-10 2.11220E-10 2.11070E-10 2.10930E-10 2.10780E-10 2.10640E-10 &
     2.10490E-10 2.10340E-10 2.10200E-10 2.10050E-10 2.09910E-10 2.09760E-10 2.09610E-10 2.09470E-10 &
     2.09320E-10 2.09180E-10 2.09030E-10 2.08890E-10 2.08740E-10 2.08600E-10 2.08450E-10 2.08310E-10 &
     2.08160E-10 2.08020E-10 2.07880E-10 2.07730E-10 2.07590E-10 2.07440E-10 2.07300E-10 2.07160E-10 &
     2.07010E-10 2.06870E-10 2.06730E-10 2.06580E-10 2.06440E-10 2.06300E-10 2.06150E-10 2.06010E-10 &
     2.05870E-10 2.05730E-10 2.05580E-10 2.05440E-10 2.05300E-10 2.05160E-10 2.05010E-10 2.04870E-10 &
     2.04730E-10 2.04590E-10 2.04450E-10 2.04310E-10 2.04170E-10 2.04020E-10 2.03880E-10 2.03740E-10 &
     2.03600E-10 2.03460E-10 2.03320E-10 2.03180E-10 2.03040E-10 2.02900E-10 2.02760E-10 2.02620E-10 &
     2.02480E-10 2.02340E-10 2.02200E-10 2.02060E-10 2.01920E-10 2.01780E-10 2.01640E-10 2.01500E-10 &
     2.01360E-10 2.01220E-10 2.01080E-10 2.00940E-10 2.00810E-10 2.00670E-10 2.00530E-10 2.00390E-10 &
     2.00250E-10 2.00110E-10 1.99980E-10 1.99840E-10 1.99700E-10 1.99560E-10 1.99420E-10 1.99290E-10 &
     1.99150E-10 1.99010E-10 1.98870E-10 1.98740E-10 1.98600E-10 1.98460E-10 1.98330E-10 1.98190E-10 &
     1.98050E-10 1.97920E-10 1.97780E-10 1.97640E-10 1.97510E-10 1.97370E-10 1.97240E-10 1.97100E-10 &
     1.96960E-10 1.96830E-10 1.96690E-10 1.96560E-10 1.96420E-10 1.96290E-10 1.96150E-10 1.96020E-10 &
     1.95880E-10 1.95750E-10 1.95610E-10 1.95480E-10 1.95340E-10 1.95210E-10 1.95070E-10 1.94940E-10 &
     1.94800E-10 1.94670E-10 1.94540E-10 1.94400E-10 1.94270E-10 1.94140E-10 1.94000E-10 1.93870E-10 &
     1.93740E-10 1.93600E-10 1.93470E-10 1.93340E-10 1.93200E-10 1.93070E-10 1.92940E-10 1.92800E-10 &
     1.92670E-10 1.92540E-10 1.92410E-10 1.92280E-10 1.92140E-10 1.92010E-10 1.91880E-10 1.91750E-10 &
     1.91620E-10 1.91480E-10 1.91350E-10 1.91220E-10 1.91090E-10 1.90960E-10 1.90830E-10 1.90700E-10 &
     1.90570E-10 1.90430E-10 1.90300E-10 1.90170E-10 1.90040E-10 1.89910E-10 1.89780E-10 1.89650E-10 &
     1.89520E-10 1.89390E-10 1.89260E-10 1.89130E-10 1.89000E-10 1.88870E-10 1.88740E-10 1.88610E-10 &
     1.88480E-10 1.88360E-10 1.88230E-10 1.88100E-10 1.87970E-10 1.87840E-10 1.87710E-10 1.87580E-10 &
     1.87450E-10 1.87320E-10 1.87200E-10 1.87070E-10 1.86940E-10 1.86810E-10 1.86680E-10 1.86560E-10 &
     1.86430E-10 1.86300E-10 1.86170E-10 1.86050E-10 1.85920E-10 1.85790E-10 1.85660E-10 1.85540E-10 &
     1.85410E-10 1.85280E-10 1.85160E-10 1.85030E-10 1.84900E-10 1.84780E-10 1.84650E-10 1.84520E-10 &
     1.84400E-10 1.84270E-10 1.84140E-10 1.84020E-10 1.83890E-10 1.83770E-10 1.83640E-10 1.83520E-10 &
     1.83390E-10 1.83260E-10 1.83140E-10 1.83010E-10 1.82890E-10 1.82760E-10 1.82640E-10 1.82510E-10 &
     1.82390E-10 1.82260E-10 1.82140E-10 1.82020E-10 1.81890E-10 1.81770E-10 1.81640E-10 1.81520E-10 &
     1.81390E-10 1.81270E-10 1.81150E-10 1.81020E-10 1.80900E-10 1.80780E-10 1.80650E-10 1.80530E-10 &
     1.80410E-10 1.80280E-10 1.80160E-10 1.80040E-10 1.79910E-10 1.79790E-10 1.79670E-10 1.79550E-10 &
     1.79420E-10 1.79300E-10 1.79180E-10 1.79060E-10 1.78930E-10 1.78810E-10 1.78690E-10 1.78570E-10 &
     1.78450E-10 1.78320E-10 1.78200E-10 1.78080E-10 1.77960E-10 1.77840E-10 1.77720E-10 1.77600E-10 &
     1.77470E-10 1.77350E-10 1.77230E-10 1.77110E-10 1.76990E-10 1.76870E-10 1.76750E-10 1.76630E-10 &
     1.76510E-10 1.76390E-10 1.76270E-10 1.76150E-10 1.76030E-10 1.75910E-10 1.75790E-10 1.75670E-10 &
     1.75550E-10 1.75430E-10 1.75310E-10 1.75190E-10 1.75070E-10 1.74950E-10 1.74830E-10 1.74710E-10 &
     1.74600E-10 1.74480E-10 1.74360E-10 1.74240E-10 1.74120E-10 1.74000E-10 1.73880E-10 1.73770E-10 &
     1.73650E-10 1.73530E-10 1.73410E-10 1.73290E-10 1.73180E-10 1.73060E-10 1.72940E-10 1.72820E-10 &
     1.72700E-10 1.72590E-10 1.72470E-10 1.72350E-10 1.72240E-10 1.72120E-10 1.72000E-10 1.71880E-10 &
     1.71770E-10 1.71650E-10 1.71530E-10 1.71420E-10 1.71300E-10 1.71180E-10 1.71070E-10 1.70950E-10 &
     1.70840E-10 1.70720E-10 1.70600E-10 1.70490E-10 1.70370E-10 1.70260E-10 1.70140E-10 1.70030E-10 &
     1.69910E-10 1.69790E-10 1.69680E-10 1.69560E-10 1.69450E-10 1.69330E-10 1.69220E-10 1.69100E-10 &
     1.68990E-10 1.68870E-10 1.68760E-10 1.68650E-10 1.68530E-10 1.68420E-10 1.68300E-10 1.68190E-10 &
     1.68070E-10 1.67960E-10 1.67850E-10 1.67730E-10 1.67620E-10 1.67500E-10 1.67390E-10 1.67280E-10 &
     1.67160E-10 1.67050E-10 1.66940E-10 1.66820E-10 1.66710E-10 1.66600E-10 1.66480E-10 1.66370E-10 &
     1.66260E-10 1.66150E-10 1.66030E-10 1.65920E-10 1.65810E-10 1.65700E-10 1.65580E-10 1.65470E-10 &
     1.65360E-10 1.65250E-10 1.65140E-10 1.65020E-10 1.64910E-10 1.64800E-10 1.64690E-10 1.64580E-10 &
     1.64470E-10 1.64350E-10 1.64240E-10 1.64130E-10 1.64020E-10 1.63910E-10 1.63800E-10 1.63690E-10 &
     1.63580E-10 1.63470E-10 1.63360E-10 1.63240E-10 1.63130E-10 1.63020E-10 1.62910E-10 1.62800E-10 &
     1.62690E-10 1.62580E-10 1.62470E-10 1.62360E-10 1.62250E-10 1.62140E-10 1.62030E-10 1.61920E-10 &
     1.61810E-10 1.61700E-10 1.61600E-10 1.61490E-10 1.61380E-10 1.61270E-10 1.61160E-10 1.61050E-10 &
     1.60940E-10 1.60830E-10 1.60720E-10 1.60610E-10 1.60510E-10 1.60400E-10 1.60290E-10 1.60180E-10 &
     1.60070E-10 1.59960E-10 1.59860E-10 1.59750E-10 1.59640E-10 1.59530E-10 1.59420E-10 1.59320E-10 &
     1.59210E-10 1.59100E-10 1.58990E-10 1.58890E-10 1.58780E-10 1.58670E-10 1.58560E-10 1.58460E-10 &
     1.58350E-10 1.58240E-10 1.58140E-10 1.58030E-10 1.57920E-10 1.57820E-10 1.57710E-10 1.57600E-10 &
     1.57500E-10 1.57390E-10 1.57280E-10 1.57180E-10 1.57070E-10 1.56970E-10 1.56860E-10 1.56750E-10 &
     1.56650E-10 1.56540E-10 1.56440E-10 1.56330E-10 1.56230E-10 1.56120E-10 1.56010E-10 1.55910E-10 &
     1.55800E-10 1.55700E-10 1.55590E-10 1.55490E-10 1.55380E-10 1.55280E-10 1.55170E-10 1.55070E-10 &
     1.54970E-10 1.54860E-10 1.54760E-10 1.54650E-10 1.54550E-10 1.54440E-10 1.54340E-10 1.54240E-10 &
     1.54130E-10 1.54030E-10 1.53920E-10 1.53820E-10 1.53720E-10 1.53610E-10 1.53510E-10 1.53410E-10 &
     1.53300E-10 1.53200E-10 1.53100E-10 1.52990E-10 1.52890E-10 1.52790E-10 1.52680E-10 1.52580E-10 &
     1.52480E-10 1.52370E-10 1.52270E-10 1.52170E-10 1.52070E-10 1.51960E-10 1.51860E-10 1.51760E-10 &
     1.51660E-10 1.51560E-10 1.51450E-10 1.51350E-10 1.51250E-10 1.51150E-10 1.51050E-10 1.50940E-10 &
     1.50840E-10 1.50740E-10 1.50640E-10 1.50540E-10 1.50440E-10 1.50340E-10 1.50230E-10 1.50130E-10 &
     1.50030E-10 1.49930E-10 1.49830E-10 1.49730E-10 1.49630E-10 1.49530E-10 1.49430E-10 1.49330E-10 &
     1.49230E-10 1.49130E-10 1.49030E-10 1.48920E-10 1.48820E-10 1.48720E-10 1.48620E-10 1.48520E-10 &
     1.48420E-10 1.48320E-10 1.48220E-10 1.48130E-10 1.48030E-10 1.47930E-10 1.47830E-10 1.47730E-10 &
     1.47630E-10 1.47530E-10 1.47430E-10 1.47330E-10 1.47230E-10 1.47130E-10 1.47030E-10 1.46930E-10 &
     1.46840E-10 1.46740E-10 1.46640E-10 1.46540E-10 1.46440E-10 1.46340E-10 1.46240E-10 1.46150E-10 &
     1.46050E-10 1.45950E-10 1.45850E-10 1.45750E-10 1.45660E-10 1.45560E-10 1.45460E-10 1.45360E-10 &
     1.45260E-10 1.45170E-10 1.45070E-10 1.44970E-10 1.44870E-10 1.44780E-10 1.44680E-10 1.44580E-10 &
     1.44490E-10 1.44390E-10 1.44290E-10 1.44190E-10 1.44100E-10 1.44000E-10 1.43900E-10 1.43810E-10 &
     1.43710E-10 1.43610E-10 1.43520E-10 1.43420E-10 1.43330E-10 1.43230E-10 1.43130E-10 1.43040E-10 &
     1.42940E-10 1.42850E-10 1.42750E-10 1.42650E-10 1.42560E-10 1.42460E-10 1.42370E-10 1.42270E-10 &
     1.42180E-10 1.42080E-10 1.41980E-10 1.41890E-10 1.41790E-10 1.41700E-10 1.41600E-10 1.41510E-10 &
     1.41410E-10 1.41320E-10 1.41220E-10 1.41130E-10 1.41030E-10 1.40940E-10 1.40850E-10 1.40750E-10 &
     1.40660E-10 1.40560E-10 1.40470E-10 1.40370E-10 1.40280E-10 1.40190E-10 1.40090E-10 1.40000E-10 &
     1.39900E-10 1.39810E-10 1.39720E-10 1.39620E-10 1.39530E-10 1.39440E-10 1.39340E-10 1.39250E-10 &
     1.39160E-10 1.39060E-10 1.38970E-10 1.38880E-10 1.38780E-10 1.38690E-10 1.38600E-10 1.38500E-10 &
     1.38410E-10 1.38320E-10 1.38230E-10 1.38130E-10 1.38040E-10 1.37950E-10 1.37850E-10 1.37760E-10 &
     1.37670E-10 1.37580E-10 1.37490E-10 1.37390E-10 1.37300E-10 1.37210E-10 1.37120E-10 1.37030E-10 &
     1.36930E-10 1.36840E-10 1.36750E-10 1.36660E-10 1.36570E-10 1.36480E-10 1.36380E-10 1.36290E-10 &
     1.36200E-10 1.36110E-10 1.36020E-10 1.35930E-10 1.35840E-10 1.35750E-10 1.35660E-10 1.35560E-10 &
     1.35470E-10 1.35380E-10 1.35290E-10 1.35200E-10 1.35110E-10 1.35020E-10 1.34930E-10 1.34840E-10 &
     1.34750E-10 1.34660E-10 1.34570E-10 1.34480E-10 1.34390E-10 1.34300E-10 1.34210E-10 1.34120E-10 &
     1.34030E-10 1.33940E-10 1.33850E-10 1.33760E-10 1.33670E-10 1.33580E-10 1.33490E-10 1.33400E-10 &
     1.33310E-10 1.33220E-10 1.33140E-10 1.33050E-10 1.32960E-10 1.32870E-10 1.32780E-10 1.32690E-10 &
     1.32600E-10 1.32510E-10 1.32420E-10 1.32340E-10 1.32250E-10 1.32160E-10 1.32070E-10 1.31980E-10 &
     1.31890E-10 1.31810E-10 1.31720E-10 1.31630E-10 1.31540E-10 1.31450E-10 1.31370E-10 1.31280E-10 &
     1.31190E-10 1.31100E-10 1.31010E-10 1.30930E-10 1.30840E-10 1.30750E-10 1.30660E-10 1.30580E-10 &
     1.30490E-10 1.30400E-10 1.30320E-10 1.30230E-10 1.30140E-10 1.30050E-10 1.29970E-10 1.29880E-10 &
     1.29790E-10 1.29710E-10 1.29620E-10 1.29530E-10 1.29450E-10 1.29360E-10 1.29270E-10 1.29190E-10 &
     1.29100E-10 1.29010E-10 1.28930E-10 1.28840E-10 1.28760E-10 1.28670E-10 1.28580E-10 1.28500E-10 &
     1.28410E-10 1.28330E-10 1.28240E-10 1.28150E-10 1.28070E-10 1.27980E-10 1.27900E-10 1.27810E-10 &
     1.27730E-10 1.27640E-10 1.27560E-10 1.27470E-10 1.27390E-10 1.27300E-10 1.27220E-10 1.27130E-10 &
     1.27050E-10 1.26960E-10 1.26880E-10 1.26790E-10 1.26710E-10 1.26620E-10 1.26540E-10 1.26450E-10 &
     1.26370E-10 1.26280E-10 1.26200E-10 1.26120E-10 1.26030E-10 1.25950E-10 1.25860E-10 1.25780E-10 &
     1.25700E-10 1.25610E-10 1.25530E-10 1.25440E-10 1.25360E-10 1.25280E-10 1.25190E-10 1.25110E-10 &
     1.25030E-10 1.24940E-10 1.24860E-10 1.24770E-10 1.24690E-10 1.24610E-10 1.24530E-10 1.24440E-10 &
     1.24360E-10 1.24280E-10 1.24190E-10 1.24110E-10 1.24030E-10 1.23940E-10 1.23860E-10 1.23780E-10 &
     1.23700E-10 1.23610E-10 1.23530E-10 1.23450E-10 1.23370E-10 1.23280E-10 1.23200E-10 1.23120E-10 &
     1.23040E-10 1.22950E-10 1.22870E-10 1.22790E-10 1.22710E-10 1.22630E-10 1.22550E-10 1.22460E-10 &
     1.22380E-10 1.22300E-10 1.22220E-10 1.22140E-10 1.22060E-10 1.21970E-10 1.21890E-10 1.21810E-10 &
     1.21730E-10 1.21650E-10 1.21570E-10 1.21490E-10 1.21410E-10 1.21320E-10 1.21240E-10 1.21160E-10 &
     1.21080E-10 1.21000E-10 1.20920E-10 1.20840E-10 1.20760E-10 1.20680E-10 1.20600E-10 1.20520E-10 &
     1.20440E-10 1.20360E-10 1.20280E-10 1.20200E-10 1.20120E-10 1.20040E-10 1.19960E-10 1.19880E-10 &
     1.19800E-10 1.19720E-10 1.19640E-10 1.19560E-10 1.19480E-10 1.19400E-10 1.19320E-10 1.19240E-10 &
     1.19160E-10 1.19080E-10 1.19000E-10 1.18920E-10 1.18840E-10 1.18760E-10 1.18680E-10 1.18600E-10 &
     1.18520E-10 1.18440E-10 1.18370E-10 1.18290E-10 1.18210E-10 1.18130E-10 1.18050E-10 1.17970E-10 &
     1.17890E-10 1.17810E-10 1.17740E-10 1.17660E-10 1.17580E-10 1.17500E-10 1.17420E-10 1.17340E-10 &
     1.17270E-10 1.17190E-10 1.17110E-10 1.17030E-10 1.16950E-10 1.16870E-10 1.16800E-10 1.16720E-10 &
     1.16640E-10 1.16560E-10 1.16490E-10 1.16410E-10 1.16330E-10 1.16250E-10 1.16180E-10 1.16100E-10 &
     1.16020E-10 1.15940E-10 1.15870E-10 1.15790E-10 1.15710E-10 1.15630E-10 1.15560E-10 1.15480E-10 &
     1.15400E-10 1.15330E-10 1.15250E-10 1.15170E-10 1.15100E-10 1.15020E-10 1.14940E-10 1.14870E-10 &
     1.14790E-10 1.14710E-10 1.14640E-10 1.14560E-10 1.14480E-10 1.14410E-10 1.14330E-10 1.14250E-10 &
     1.14180E-10 1.14100E-10 1.14030E-10 1.13950E-10 1.13870E-10 1.13800E-10 1.13720E-10 1.13650E-10 &
     1.13570E-10 1.13500E-10 1.13420E-10 1.13340E-10 1.13270E-10 1.13190E-10 1.13120E-10 1.13040E-10 &
     1.12970E-10 1.12890E-10 1.12820E-10 1.12740E-10 1.12670E-10 1.12590E-10 1.12520E-10 1.12440E-10 &
     1.12370E-10 1.12290E-10 1.12220E-10 1.12140E-10 1.12070E-10 1.11990E-10 1.11920E-10 1.11840E-10 &
     1.11770E-10 1.11690E-10 1.11620E-10 1.11540E-10 1.11470E-10 1.11400E-10 1.11320E-10 1.11250E-10 &
     1.11170E-10 1.11100E-10 1.11030E-10 1.10950E-10 1.10880E-10 1.10800E-10 1.10730E-10 1.10660E-10 &
     1.10580E-10 1.10510E-10 1.10430E-10 1.10360E-10 1.10290E-10 1.10210E-10 1.10140E-10 1.10070E-10 &
     1.09990E-10 1.09920E-10 1.09850E-10 1.09770E-10 1.09700E-10 1.09630E-10 1.09560E-10 1.09480E-10 &
     1.09410E-10 1.09340E-10 1.09260E-10 1.09190E-10 1.09120E-10 1.09050E-10 1.08970E-10 1.08900E-10 &
     1.08830E-10 1.08750E-10 1.08680E-10 1.08610E-10 1.08540E-10 1.08470E-10 1.08390E-10 1.08320E-10 &
     1.08250E-10 1.08180E-10 1.08100E-10 1.08030E-10 1.07960E-10 1.07890E-10 1.07820E-10 1.07740E-10 &
     1.07670E-10 1.07600E-10 1.07530E-10 1.07460E-10 1.07390E-10 1.07310E-10 1.07240E-10 1.07170E-10 &
     1.07100E-10 1.07030E-10 1.06960E-10 1.06890E-10 1.06820E-10 1.06740E-10 1.06670E-10 1.06600E-10 &
     1.06530E-10 1.06460E-10 1.06390E-10 1.06320E-10 1.06250E-10 1.06180E-10 1.06110E-10 1.06030E-10 &
     1.05960E-10 1.05890E-10 1.05820E-10 1.05750E-10 1.05680E-10 1.05610E-10 1.05540E-10 1.05470E-10 &
     1.05400E-10 1.05330E-10 1.05260E-10 1.05190E-10 1.05120E-10 1.05050E-10 1.04980E-10 1.04910E-10 &
     1.04840E-10 1.04770E-10 1.04700E-10 1.04630E-10 1.04560E-10 1.04490E-10 1.04420E-10 1.04350E-10 &
     1.04280E-10 1.04210E-10 1.04140E-10 1.04070E-10 1.04000E-10 1.03940E-10 1.03870E-10 1.03800E-10 &
     1.03730E-10 1.03660E-10 1.03590E-10 1.03520E-10 1.03450E-10 1.03380E-10 1.03310E-10 1.03250E-10 &
     1.03180E-10 1.03110E-10 1.03040E-10 1.02970E-10 1.02900E-10 1.02830E-10 1.02760E-10 1.02700E-10 &
     1.02630E-10 1.02560E-10 1.02490E-10 1.02420E-10 1.02350E-10 1.02290E-10 1.02220E-10 1.02150E-10 &
     1.02080E-10 1.02010E-10 1.01950E-10 1.01880E-10 1.01810E-10 1.01740E-10 1.01670E-10 1.01610E-10 &
     1.01540E-10 1.01470E-10 1.01400E-10 1.01340E-10 1.01270E-10 1.01200E-10 1.01130E-10 1.01070E-10 &
     1.01000E-10 1.00930E-10 1.00860E-10 1.00800E-10 1.00730E-10 1.00660E-10 1.00600E-10 1.00530E-10 &
     1.00460E-10 1.00390E-10 1.00330E-10 1.00260E-10 1.00190E-10 1.00130E-10 1.00060E-10 9.99930E-11 &
     9.99270E-11 9.98600E-11 9.97940E-11 9.97270E-11 9.96610E-11 9.95940E-11 9.95280E-11 9.94610E-11 &
     9.93950E-11 9.93290E-11 9.92630E-11 9.91960E-11 9.91300E-11 9.90640E-11 9.89980E-11 9.89320E-11 &
     9.88660E-11 9.88000E-11 9.87350E-11 9.86690E-11 9.86030E-11 9.85370E-11 9.84720E-11 9.84060E-11 &
     9.83400E-11 9.82750E-11 9.82090E-11 9.81440E-11 9.80780E-11 9.80130E-11 9.79480E-11 9.78820E-11 &
     9.78170E-11 9.77520E-11 9.76870E-11 9.76220E-11 9.75570E-11 9.74920E-11 9.74270E-11 9.73620E-11 &
     9.72970E-11 9.72320E-11 9.71670E-11 9.71020E-11 9.70380E-11 9.69730E-11 9.69080E-11 9.68440E-11 &
     9.67790E-11 9.67150E-11 9.66500E-11 9.65860E-11 9.65210E-11 9.64570E-11 9.63930E-11 9.63280E-11 &
     9.62640E-11 9.62000E-11 9.61360E-11 9.60720E-11 9.60080E-11 9.59440E-11 9.58800E-11 9.58160E-11 &
     9.57520E-11 9.56880E-11 9.56240E-11 9.55610E-11 9.54970E-11 9.54330E-11 9.53700E-11 9.53060E-11 &
     9.52430E-11 9.51790E-11 9.51160E-11 9.50520E-11 9.49890E-11 9.49260E-11 9.48620E-11 9.47990E-11 &
     9.47360E-11 9.46730E-11 9.46100E-11 9.45470E-11 9.44840E-11 9.44210E-11 9.43580E-11 9.42950E-11 &
     9.42320E-11 9.41690E-11 9.41060E-11 9.40430E-11 9.39810E-11 9.39180E-11 9.38560E-11 9.37930E-11 &
     9.37300E-11 9.36680E-11 9.36050E-11 9.35430E-11 9.34810E-11 9.34180E-11 9.33560E-11 9.32940E-11 &
     9.32320E-11 9.31690E-11 9.31070E-11 9.30450E-11 9.29830E-11 9.29210E-11 9.28590E-11 9.27970E-11 &
     9.27350E-11 9.26730E-11 9.26120E-11 9.25500E-11 9.24880E-11 9.24260E-11 9.23650E-11 9.23030E-11 &
     9.22420E-11 9.21800E-11 9.21190E-11 9.20570E-11 9.19960E-11 9.19340E-11 9.18730E-11 9.18120E-11 &
     9.17500E-11 9.16890E-11 9.16280E-11 9.15670E-11 9.15060E-11 9.14450E-11 9.13840E-11 9.13230E-11 &
     9.12620E-11 9.12010E-11 9.11400E-11 9.10790E-11 9.10180E-11 9.09580E-11 9.08970E-11 9.08360E-11 &
     9.07760E-11 9.07150E-11 9.06540E-11 9.05940E-11 9.05330E-11 9.04730E-11 9.04130E-11 9.03520E-11 &
     9.02920E-11 9.02320E-11 9.01710E-11 9.01110E-11 9.00510E-11 8.99910E-11 8.99310E-11 8.98710E-11 &
     8.98110E-11 8.97510E-11 8.96910E-11 8.96310E-11 8.95710E-11 8.95110E-11 8.94510E-11 8.93910E-11 &
     8.93320E-11 8.92720E-11 8.92120E-11 8.91530E-11 8.90930E-11 8.90340E-11 8.89740E-11 8.89150E-11 &
     8.88550E-11 8.87960E-11 8.87370E-11 8.86770E-11 8.86180E-11 8.85590E-11 8.84990E-11 8.84400E-11 &
     8.83810E-11 8.83220E-11 8.82630E-11 8.82040E-11 8.81450E-11 8.80860E-11 8.80270E-11 8.79680E-11 &
     8.79090E-11 8.78510E-11 8.77920E-11 8.77330E-11 8.76740E-11 8.76160E-11 8.75570E-11 8.74980E-11 &
     8.74400E-11 8.73810E-11 8.73230E-11 8.72640E-11 8.72060E-11 8.71480E-11 8.70890E-11 8.70310E-11 &
     8.69730E-11 8.69140E-11 8.68560E-11 8.67980E-11 8.67400E-11 8.66820E-11 8.66240E-11 8.65660E-11 &
     8.65080E-11 8.64500E-11 8.63920E-11 8.63340E-11 8.62760E-11 8.62180E-11 8.61610E-11 8.61030E-11 &
     8.60450E-11 8.59870E-11 8.59300E-11 8.58720E-11 8.58150E-11 8.57570E-11 8.57000E-11 8.56420E-11 &
     8.55850E-11 8.55270E-11 8.54700E-11 8.54130E-11 8.53550E-11 8.52980E-11 8.52410E-11 8.51840E-11 &
     8.51260E-11 8.50690E-11 8.50120E-11 8.49550E-11 8.48980E-11 8.48410E-11 8.47840E-11 8.47270E-11 &
     8.46710E-11 8.46140E-11 8.45570E-11 8.45000E-11 8.44430E-11 8.43870E-11 8.43300E-11 8.42730E-11 &
     8.42170E-11 8.41600E-11 8.41040E-11 8.40470E-11 8.39910E-11 8.39340E-11 8.38780E-11 8.38220E-11 &
     8.37650E-11 8.37090E-11 8.36530E-11 8.35970E-11 8.35400E-11 8.34840E-11 8.34280E-11 8.33720E-11 &
     8.33160E-11 8.32600E-11 8.32040E-11 8.31480E-11 8.30920E-11 8.30370E-11 8.29810E-11 8.29250E-11 &
     8.28690E-11 8.28140E-11 8.27580E-11 8.27020E-11 8.26470E-11 8.25910E-11 8.25360E-11 8.24800E-11 &
     8.24250E-11 8.23690E-11 8.23140E-11 8.22580E-11 8.22030E-11 8.21480E-11 8.20930E-11 8.20370E-11 &
     8.19820E-11 8.19270E-11 8.18720E-11 8.18170E-11 8.17620E-11 8.17070E-11 8.16520E-11 8.15970E-11 &
     8.15420E-11 8.14870E-11 8.14330E-11 8.13780E-11 8.13230E-11 8.12680E-11 8.12140E-11 8.11590E-11 &
     8.11050E-11 8.10500E-11 8.09950E-11 8.09410E-11 8.08870E-11 8.08320E-11 8.07780E-11 8.07230E-11 &
     8.06690E-11 8.06150E-11 8.05610E-11 8.05060E-11 8.04520E-11 8.03980E-11 8.03440E-11 8.02900E-11 &
     8.02360E-11 8.01820E-11 8.01280E-11 8.00740E-11 8.00200E-11 7.99670E-11 7.99130E-11 7.98590E-11 &
     7.98050E-11 7.97520E-11 7.96980E-11 7.96440E-11 7.95910E-11 7.95370E-11 7.94840E-11 7.94300E-11 &
     7.93770E-11 7.93230E-11 7.92700E-11 7.92170E-11 7.91630E-11 7.91100E-11 7.90570E-11 7.90040E-11 &
     7.89510E-11 7.88980E-11 7.88450E-11 7.87920E-11 7.87390E-11 7.86860E-11 7.86330E-11 7.85800E-11 &
     7.85270E-11 7.84740E-11 7.84210E-11 7.83690E-11 7.83160E-11 7.82630E-11 7.82110E-11 7.81580E-11 &
     7.81060E-11 7.80530E-11 7.80010E-11 7.79480E-11 7.78960E-11 7.78430E-11 7.77910E-11 7.77390E-11 &
     7.76870E-11 7.76340E-11 7.75820E-11 7.75300E-11 7.74780E-11 7.74260E-11 7.73740E-11 7.73220E-11 &
     7.72700E-11 7.72180E-11 7.71660E-11 7.71140E-11 7.70620E-11 7.70110E-11 7.69590E-11 7.69070E-11 &
     7.68560E-11 7.68040E-11 7.67520E-11 7.67010E-11 7.66490E-11 7.65980E-11 7.65460E-11 7.64950E-11 &
     7.64440E-11 7.63920E-11 7.63410E-11 7.62900E-11 7.62390E-11 7.61880E-11 7.61360E-11 7.60850E-11 &
     7.60340E-11 7.59830E-11 7.59320E-11 7.58810E-11 7.58300E-11 7.57800E-11 7.57290E-11 7.56780E-11 &
     7.56270E-11 7.55760E-11 7.55260E-11 7.54750E-11 7.54250E-11 7.53740E-11 7.53230E-11 7.52730E-11 &
     7.52220E-11 7.51720E-11 7.51220E-11 7.50710E-11 7.50210E-11 7.49710E-11 7.49210E-11 7.48700E-11 &
     7.48200E-11 7.47700E-11 7.47200E-11 7.46700E-11 7.46200E-11 7.45700E-11 7.45200E-11 7.44700E-11 &
     7.44200E-11 7.43710E-11 7.43210E-11 7.42710E-11 7.42210E-11 7.41720E-11 7.41220E-11 7.40730E-11 &
     7.40230E-11 7.39730E-11 7.39240E-11 7.38750E-11 7.38250E-11 7.37760E-11 7.37260E-11 7.36770E-11 &
     7.36280E-11 7.35790E-11 7.35300E-11 7.34810E-11 7.34310E-11 7.33820E-11 7.33330E-11 7.32840E-11 &
     7.32350E-11 7.31870E-11 7.31380E-11 7.30890E-11 7.30400E-11 7.29910E-11 7.29430E-11 7.28940E-11 &
     7.28450E-11 7.27970E-11 7.27480E-11 7.27000E-11 7.26510E-11 7.26030E-11 7.25540E-11 7.25060E-11 &
     7.24580E-11 7.24090E-11 7.23610E-11 7.23130E-11 7.22650E-11 7.22170E-11 7.21690E-11 7.21210E-11 &
     7.20730E-11 7.20250E-11 7.19770E-11 7.19290E-11 7.18810E-11 7.18330E-11 7.17850E-11 7.17370E-11 &
     7.16900E-11 7.16420E-11 7.15940E-11 7.15470E-11 7.14990E-11 7.14520E-11 7.14040E-11 7.13570E-11 &
     7.13090E-11 7.12620E-11 7.12140E-11 7.11670E-11 7.11200E-11 7.10720E-11 7.10250E-11 7.09780E-11 &
     7.09310E-11 7.08840E-11 7.08370E-11 7.07900E-11 7.07430E-11 7.06960E-11 7.06490E-11 7.06020E-11 &
     7.05550E-11 7.05080E-11 7.04610E-11 7.04140E-11 7.03680E-11 7.03210E-11 7.02740E-11 7.02270E-11 &
     7.01810E-11 7.01340E-11 7.00880E-11 7.00410E-11 6.99950E-11 6.99480E-11 6.99020E-11 6.98550E-11 &
     6.98090E-11 6.97630E-11 6.97160E-11 6.96700E-11 6.96240E-11 6.95770E-11 6.95310E-11 6.94850E-11 &
     6.94390E-11 6.93930E-11 6.93470E-11 6.93010E-11 6.92550E-11 6.92090E-11 6.91630E-11 6.91170E-11 &
     6.90710E-11 6.90250E-11 6.89790E-11 6.89340E-11 6.88880E-11 6.88420E-11 6.87960E-11 6.87510E-11 &
     6.87050E-11 6.86590E-11 6.86140E-11 6.85680E-11 6.85230E-11 6.84770E-11 6.84320E-11 6.83860E-11 &
     6.83410E-11 6.82950E-11 6.82500E-11 6.82050E-11 6.81590E-11 6.81140E-11 6.80690E-11 6.80240E-11 &
     6.79790E-11 6.79330E-11 6.78880E-11 6.78430E-11 6.77980E-11 6.77530E-11 6.77080E-11 6.76630E-11 &
     6.76180E-11 6.75730E-11 6.75280E-11 6.74830E-11 6.74380E-11 6.73930E-11 6.73480E-11 6.73040E-11 &
     6.72590E-11 6.72140E-11 6.71690E-11 6.71250E-11 6.70800E-11 6.70350E-11 6.69910E-11 6.69460E-11 &
     6.69010E-11 6.68570E-11 6.68120E-11 6.67680E-11 6.67230E-11 6.66790E-11 6.66340E-11 6.65900E-11 &
     6.65460E-11 6.65010E-11 6.64570E-11 6.64130E-11 6.63680E-11 6.63240E-11 6.62800E-11 6.62360E-11 &
     6.61910E-11 6.61470E-11 6.61030E-11 6.60590E-11 6.60150E-11 6.59710E-11 6.59270E-11 6.58830E-11 &
     6.58380E-11 6.57940E-11 6.57500E-11 6.57070E-11 6.56630E-11 6.56190E-11 6.55750E-11 6.55310E-11 &
     6.54870E-11 6.54430E-11 6.53990E-11 6.53560E-11 6.53120E-11 6.52680E-11 6.52240E-11 6.51810E-11 &
     6.51370E-11 6.50930E-11 6.50500E-11 6.50060E-11 6.49620E-11 6.49190E-11 6.48750E-11 6.48320E-11 &
     6.47880E-11 6.47440E-11 6.47010E-11 6.46570E-11 6.46140E-11 6.45710E-11 6.45270E-11 6.44840E-11 &
     6.44400E-11 6.43970E-11 6.43540E-11 6.43100E-11 6.42670E-11 6.42240E-11 6.41800E-11 6.41370E-11 &
     6.40940E-11 6.40510E-11 6.40070E-11 6.39640E-11 6.39210E-11 6.38780E-11 6.38350E-11 6.37910E-11 &
     6.37480E-11 6.37050E-11 6.36620E-11 6.36190E-11 6.35760E-11 6.35330E-11 6.34900E-11 6.34470E-11 &
     6.34040E-11 6.33610E-11 6.33180E-11 6.32750E-11 6.32320E-11 6.31890E-11 6.31460E-11 6.31030E-11 &
     6.30610E-11 6.30180E-11 6.29750E-11 6.29320E-11 6.28890E-11 6.28460E-11 6.28040E-11 6.27610E-11 &
     6.27180E-11 6.26750E-11 6.26330E-11 6.25900E-11 6.25470E-11 6.25040E-11 6.24620E-11 6.24190E-11 &
     6.23760E-11 6.23340E-11 6.22910E-11 6.22490E-11 6.22060E-11 6.21630E-11 6.21210E-11 6.20780E-11 &
     6.20360E-11 6.19930E-11 6.19510E-11 6.19080E-11 6.18660E-11 6.18230E-11 6.17810E-11 6.17380E-11 &
     6.16960E-11 6.16530E-11 6.16110E-11 6.15680E-11 6.15260E-11 6.14830E-11 6.14410E-11 6.13990E-11 &
     6.13560E-11 6.13140E-11 6.12720E-11 6.12290E-11 6.11870E-11 6.11450E-11 6.11020E-11 6.10600E-11 &
     6.10180E-11 6.09750E-11 6.09330E-11 6.08910E-11 6.08480E-11 6.08060E-11 6.07640E-11 6.07220E-11 &
     6.06790E-11 6.06370E-11 6.05950E-11 6.05530E-11 6.05110E-11 6.04680E-11 6.04260E-11 6.03840E-11 &
     6.03420E-11 6.03000E-11 6.02580E-11 6.02150E-11 6.01730E-11 6.01310E-11 6.00890E-11 6.00470E-11 &
     6.00050E-11 5.99630E-11 5.99210E-11 5.98780E-11 5.98360E-11 5.97940E-11 5.97520E-11 5.97100E-11 &
     5.96680E-11 5.96260E-11 5.95840E-11 5.95420E-11 5.95000E-11 5.94580E-11 5.94160E-11 5.93740E-11 &
     5.93320E-11 5.92900E-11 5.92480E-11 5.92060E-11 5.91640E-11 5.91220E-11 5.90800E-11 5.90380E-11 &
     5.89960E-11 5.89540E-11 5.89120E-11 5.88700E-11 5.88280E-11 5.87860E-11 5.87440E-11 5.87020E-11 &
     5.86600E-11 5.86180E-11 5.85760E-11 5.85340E-11 5.84920E-11 5.84500E-11 5.84080E-11 5.83660E-11 &
     5.83240E-11 5.82820E-11 5.82400E-11 5.81980E-11 5.81560E-11 5.81140E-11 5.80720E-11 5.80310E-11 &
     5.79890E-11 5.79470E-11 5.79050E-11 5.78630E-11 5.78210E-11 5.77790E-11 5.77370E-11 5.76950E-11 &
     5.76530E-11 5.76110E-11 1.33820E-09 1.32760E-09 1.29090E-09 1.24460E-09 1.19540E-09 1.14700E-09 &
     1.09960E-09 1.05300E-09 1.00790E-09 9.66160E-10 9.27750E-10 8.93940E-10 8.62100E-10 8.32270E-10 &
     8.03160E-10 7.75210E-10 7.49280E-10 7.24550E-10 7.01300E-10 6.78710E-10 6.56980E-10 6.36390E-10 &
     6.16620E-10 5.98260E-10 5.80800E-10 5.64010E-10 5.48190E-10 5.32680E-10 5.18440E-10 5.04420E-10 &
     4.91130E-10 4.78320E-10 4.65630E-10 4.53570E-10 4.41680E-10 4.30400E-10 4.19840E-10 4.09840E-10 &
     4.00260E-10 3.90930E-10 3.81800E-10 3.72840E-10 3.64130E-10 3.55770E-10 3.47720E-10 3.39950E-10 &
     3.32470E-10 3.25270E-10 3.18210E-10 3.11170E-10 3.04190E-10 2.97380E-10 2.90800E-10 2.84530E-10 &
     2.78550E-10 2.72860E-10 2.67490E-10 2.62470E-10 2.57650E-10 2.52870E-10 2.48100E-10 2.43340E-10 &
     2.38570E-10 2.33810E-10 2.29100E-10 2.24480E-10 2.19990E-10 2.15680E-10 2.11550E-10 2.07600E-10 &
     2.03830E-10 2.00220E-10 1.96740E-10 1.93370E-10 1.90070E-10 1.86800E-10 1.83560E-10 1.80350E-10 &
     1.77160E-10 1.73980E-10 1.70840E-10 1.67740E-10 1.64700E-10 1.61730E-10 1.58830E-10 1.56030E-10 &
     1.53340E-10 1.50740E-10 1.48230E-10 1.45800E-10 1.43450E-10 1.41170E-10 1.38950E-10 1.36780E-10 &
     1.34660E-10 1.32590E-10 1.30550E-10 1.28550E-10 1.26580E-10 1.24630E-10 1.22710E-10 1.20810E-10 &
     1.18930E-10 1.17080E-10 1.15250E-10 1.13450E-10 1.11680E-10 1.09930E-10 1.08210E-10 1.06530E-10 &
     1.04870E-10 1.03250E-10 1.01660E-10 1.00100E-10 9.85840E-11 9.71010E-11 9.56560E-11 9.42480E-11 &
     9.28780E-11 9.15450E-11 9.02520E-11 8.89960E-11 8.77770E-11 8.65890E-11 8.54300E-11 8.42970E-11 &
     8.31860E-11 8.20930E-11 8.10170E-11 7.99540E-11 7.89070E-11 7.78730E-11 7.68530E-11 7.58450E-11 &
     7.48500E-11 7.38680E-11 7.28980E-11 7.19420E-11 7.10010E-11 7.00760E-11 6.91680E-11 6.82770E-11 &
     6.74050E-11 6.65510E-11 6.57150E-11 6.48960E-11 6.40960E-11 6.33120E-11 6.25460E-11 6.17960E-11 &
     6.10620E-11 6.03420E-11 5.96360E-11 5.89420E-11 5.82590E-11 5.75860E-11 5.69210E-11 5.62660E-11 &
     5.56210E-11 5.49850E-11 5.43600E-11 5.37450E-11 5.31410E-11 5.25480E-11 5.19660E-11 5.13960E-11 &
     5.08360E-11 5.02880E-11 4.97500E-11 4.92230E-11 4.87080E-11 4.82030E-11 4.77090E-11 4.72250E-11 &
     4.67520E-11 4.62890E-11 4.58360E-11 4.53930E-11 4.49610E-11 4.45390E-11 4.41270E-11 4.37270E-11 &
     4.33380E-11 4.29600E-11 4.25930E-11 4.22370E-11 4.18900E-11 4.15520E-11 4.12200E-11 4.08950E-11 &
     4.05750E-11 4.02590E-11 3.99470E-11 3.96410E-11 3.93410E-11 3.90480E-11 3.87620E-11 3.84850E-11 &
     3.82160E-11 3.79540E-11 3.77000E-11 3.74510E-11 3.72060E-11 3.69650E-11 3.67260E-11 3.64880E-11 &
     3.62530E-11 3.60210E-11 3.57940E-11 3.55720E-11 3.53580E-11 3.51520E-11 3.49550E-11 3.47660E-11 &
     3.45860E-11 3.44110E-11 3.42420E-11 3.40780E-11 3.39170E-11 3.37590E-11 3.36040E-11 3.34530E-11 &
     3.33070E-11 3.31660E-11 3.30310E-11 3.29030E-11 3.27830E-11 3.26690E-11 3.25610E-11 3.24570E-11 &
     3.23560E-11 3.22580E-11 3.21610E-11 3.20650E-11 3.19710E-11 3.18780E-11 3.17880E-11 3.17030E-11 &
     3.16220E-11 3.15470E-11 3.14780E-11 3.14150E-11 3.13570E-11 3.13030E-11 3.12530E-11 3.12050E-11 &
     3.11600E-11 3.11170E-11 3.10760E-11 3.10400E-11 3.10100E-11 3.09880E-11 3.09750E-11 3.09730E-11 &
     3.09830E-11 3.10040E-11 3.10360E-11 3.10770E-11 3.11280E-11 3.11870E-11 3.12520E-11 3.13250E-11 &
     3.14050E-11 3.14930E-11 3.15900E-11 3.16980E-11 3.18180E-11 3.19500E-11 3.20950E-11 3.22510E-11 &
     3.24170E-11 3.25900E-11 3.27690E-11 3.29500E-11 3.31330E-11 3.33150E-11 3.34960E-11 3.36750E-11 &
     3.38520E-11 3.40280E-11 3.42010E-11 3.43710E-11 3.45390E-11 3.47020E-11 3.48590E-11 3.50070E-11 &
     3.51450E-11 3.52720E-11 3.53850E-11 3.54830E-11 3.55660E-11 3.56330E-11 3.56830E-11 3.57180E-11 &
     3.57350E-11 3.57360E-11 3.57200E-11 3.56870E-11 3.56370E-11 3.55720E-11 3.54900E-11 3.53930E-11 &
     3.52800E-11 3.51520E-11 3.50090E-11 3.48520E-11 3.46800E-11 3.44940E-11 3.42940E-11 3.40800E-11 &
     3.38530E-11 3.36150E-11 3.33670E-11 3.31110E-11 3.28490E-11 3.25840E-11 3.23170E-11 3.20500E-11 &
     3.17820E-11 3.15130E-11 3.12440E-11 3.09730E-11 3.07000E-11 3.04250E-11 3.01470E-11 2.98680E-11 &
     2.95890E-11 2.93110E-11 2.90370E-11 2.87670E-11 2.85020E-11 2.82440E-11 2.79930E-11 2.77490E-11 &
     2.75100E-11 2.72770E-11 2.70500E-11 2.68270E-11 2.66090E-11 2.63950E-11 2.61840E-11 2.59760E-11 &
     2.57700E-11 2.55650E-11 2.53620E-11 2.51590E-11 2.49570E-11 2.47570E-11 2.45610E-11 2.43670E-11 &
     2.41790E-11 2.39960E-11 2.38190E-11 2.36470E-11 2.34800E-11 2.33150E-11 2.31520E-11 2.29900E-11 &
     2.28270E-11 2.26630E-11 2.24980E-11 2.23320E-11 2.21660E-11 2.19990E-11 2.18330E-11 2.16680E-11 &
     2.15030E-11 2.13410E-11 2.11800E-11 2.10220E-11 2.08660E-11 2.07140E-11 2.05650E-11 2.04200E-11 &
     2.02780E-11 2.01390E-11 2.00010E-11 1.98630E-11 1.97260E-11 1.95870E-11 1.94470E-11 1.93060E-11 &
     1.91640E-11 1.90210E-11 1.88790E-11 1.87380E-11 1.85980E-11 1.84590E-11 1.83230E-11 1.81880E-11 &
     1.80570E-11 1.79270E-11 1.78010E-11 1.76780E-11 1.75590E-11 1.74420E-11 1.73270E-11 1.72140E-11 &
     1.71020E-11 1.69910E-11 1.68790E-11 1.67680E-11 1.66550E-11 1.65420E-11 1.64280E-11 1.63140E-11 &
     1.61990E-11 1.60830E-11 1.59660E-11 1.58500E-11 1.57340E-11 1.56190E-11 1.55070E-11 1.53970E-11 &
     1.52910E-11 1.51890E-11 1.50900E-11 1.49940E-11 1.49000E-11 1.48080E-11 1.47160E-11 1.46240E-11 &
     1.45300E-11 1.44360E-11 1.43410E-11 1.42450E-11 1.41470E-11 1.40470E-11 1.39460E-11 1.38440E-11 &
     1.37410E-11 1.36370E-11 1.35330E-11 1.34310E-11 1.33300E-11 1.32310E-11 1.31350E-11 1.30410E-11 &
     1.29500E-11 1.28610E-11 1.27750E-11 1.26900E-11 1.26060E-11 1.25240E-11 1.24440E-11 1.23640E-11 &
     1.22840E-11 1.22040E-11 1.21240E-11 1.20430E-11 1.19620E-11 1.18790E-11 1.17960E-11 1.17120E-11 &
     1.16290E-11 1.15450E-11 1.14620E-11 1.13800E-11 1.12990E-11 1.12180E-11 1.11370E-11 1.10570E-11 &
     1.09770E-11 1.08970E-11 1.08180E-11 1.07390E-11 1.06610E-11 1.05850E-11 1.05110E-11 1.04400E-11 &
     1.03720E-11 1.03080E-11 1.02470E-11 1.01890E-11 1.01310E-11 1.00730E-11 1.00150E-11 9.95560E-12 &
     9.89360E-12 9.82940E-12 9.76360E-12 9.69660E-12 9.62890E-12 9.56080E-12 9.49290E-12 9.42550E-12 &
     9.35870E-12 9.29250E-12 9.22690E-12 9.16210E-12 9.09780E-12 9.03430E-12 8.97150E-12 8.90950E-12 &
     8.84840E-12 8.78840E-12 8.72970E-12 8.67220E-12 8.61620E-12 8.56180E-12 8.50880E-12 8.45710E-12 &
     8.40670E-12 8.35730E-12 8.30890E-12 8.26140E-12 8.21460E-12 8.16830E-12 8.12220E-12 8.07630E-12 &
     8.03010E-12 7.98360E-12 7.93650E-12 7.88860E-12 7.84010E-12 7.79120E-12 7.74220E-12 7.69330E-12 &
     7.64480E-12 7.59690E-12 7.54970E-12 7.50300E-12 7.45630E-12 7.40920E-12 7.36140E-12 7.31240E-12 &
     7.26190E-12 7.20970E-12 7.15640E-12 7.10290E-12 7.05020E-12 6.99910E-12 6.95060E-12 6.90540E-12 &
     6.86420E-12 6.82640E-12 6.79110E-12 6.75730E-12 6.72420E-12 6.69080E-12 6.65630E-12 6.62000E-12 &
     6.58220E-12 6.54370E-12 6.50510E-12 6.46710E-12 6.43040E-12 6.39560E-12 6.36310E-12 6.33260E-12 &
     6.30330E-12 6.27460E-12 6.24590E-12 6.21640E-12 6.18550E-12 6.15270E-12 6.11840E-12 6.08320E-12 &
     6.04760E-12 6.01230E-12 5.97780E-12 5.94460E-12 5.91320E-12 5.88330E-12 5.85420E-12 5.82550E-12 &
     5.79670E-12 5.76720E-12 5.73640E-12 5.70410E-12 5.67040E-12 5.63580E-12 5.60070E-12 5.56560E-12 &
     5.53070E-12 5.49660E-12 5.46350E-12 5.43150E-12 5.40030E-12 5.36980E-12 5.33990E-12 5.31050E-12 &
     5.28130E-12 5.25220E-12 5.22330E-12 5.19450E-12 5.16580E-12 5.13710E-12 5.10850E-12 5.07990E-12 &
     5.05140E-12 5.02310E-12 4.99510E-12 4.96770E-12 4.94090E-12 4.91510E-12 4.89030E-12 4.86670E-12 &
     4.84410E-12 4.82220E-12 4.80060E-12 4.77920E-12 4.75770E-12 4.73570E-12 4.71320E-12 4.69030E-12 &
     4.66720E-12 4.64450E-12 4.62220E-12 4.60080E-12 4.58050E-12 4.56150E-12 4.54350E-12 4.52610E-12 &
     4.50880E-12 4.49130E-12 4.47310E-12 4.45380E-12 4.43300E-12 4.41130E-12 4.38900E-12 4.36670E-12 &
     4.34490E-12 4.32400E-12 4.30470E-12 4.28720E-12 4.27120E-12 4.25600E-12 4.24100E-12 4.22570E-12 &
     4.20950E-12 4.19180E-12 4.17220E-12 4.15120E-12 4.12940E-12 4.10750E-12 4.08620E-12 4.06600E-12 &
     4.04770E-12 4.03170E-12 4.01750E-12 4.00440E-12 3.99170E-12 3.97870E-12 3.96460E-12 3.94880E-12 &
     3.93080E-12 3.91110E-12 3.89040E-12 3.86940E-12 3.84890E-12 3.82970E-12 3.81230E-12 3.79740E-12 &
     3.78440E-12 3.77250E-12 3.76100E-12 3.74920E-12 3.73620E-12 3.72130E-12 3.70400E-12 3.68480E-12 &
     3.66450E-12 3.64370E-12 3.62330E-12 3.60400E-12 3.58660E-12 3.57150E-12 3.55840E-12 3.54680E-12 &
     3.53590E-12 3.52520E-12 3.51400E-12 3.50190E-12 3.48830E-12 3.47290E-12 3.45590E-12 3.43720E-12 &
     3.41660E-12 3.39430E-12 3.37000E-12 3.34430E-12 3.31910E-12 3.29690E-12 3.27990E-12 3.27080E-12 &
     3.27170E-12 3.28530E-12 3.31160E-12 3.34170E-12 3.36430E-12 3.36820E-12 3.34220E-12 3.27510E-12 &
     3.15560E-12 2.97710E-12 2.75200E-12 2.49700E-12 2.22910E-12 1.96510E-12 1.72200E-12 1.51650E-12 &
     1.36120E-12 1.25120E-12 1.17720E-12 1.12980E-12 1.09980E-12 1.07780E-12 1.05460E-12 1.02260E-12 &
     9.81950E-13 9.34350E-13 8.81640E-13 8.25650E-13 7.68200E-13 7.11090E-13 6.55820E-13 6.02460E-13 &
     5.50770E-13 5.00470E-13 4.51310E-13 4.03040E-13 3.55400E-13 3.08170E-13 2.61260E-13 2.14640E-13 &
     1.68270E-13 1.22100E-13 7.60950E-14 3.02120E-14 -1.55550E-14 -6.10650E-14 -1.06140E-13 -1.50600E-13 &
     -1.94280E-13 -2.36980E-13 -2.78550E-13 -3.18840E-13 -3.57940E-13 -3.95990E-13 -4.33100E-13 -4.69410E-13 &
     -5.05050E-13 -5.40140E-13 -5.74830E-13 -6.09220E-13 -6.43440E-13 -6.77610E-13 -7.11850E-13 -7.46280E-13 &
     -7.81020E-13 -8.16120E-13 -8.51320E-13 -8.86320E-13 -9.20780E-13 -9.54390E-13 -9.86840E-13 -1.01780E-12 &
     -1.04710E-12 -1.07500E-12 -1.10210E-12 -1.12890E-12 -1.15580E-12 -1.18340E-12 -1.21210E-12 -1.24230E-12 &
     -1.27320E-12 -1.30410E-12 -1.33410E-12 -1.36240E-12 -1.38820E-12 -1.41060E-12 -1.42920E-12 -1.44450E-12 &
     -1.45750E-12 -1.46930E-12 -1.48070E-12 -1.49270E-12 -1.50620E-12 -1.52200E-12 -1.53920E-12 -1.55700E-12 &
     -1.57430E-12 -1.59010E-12 -1.60340E-12 -1.61330E-12 -1.61900E-12 -1.62150E-12 -1.62210E-12 -1.62190E-12 &
     -1.62230E-12 -1.62440E-12 -1.62970E-12 -1.63880E-12 -1.65120E-12 -1.66550E-12 -1.68060E-12 -1.69550E-12 &
     -1.70890E-12 -1.71970E-12 -1.72700E-12 -1.73130E-12 -1.73340E-12 -1.73400E-12 -1.73400E-12 -1.73420E-12 &
     -1.73520E-12 -1.73780E-12 -1.74130E-12 -1.74530E-12 -1.74900E-12 -1.75180E-12 -1.75300E-12 -1.75200E-12 &
     -1.74820E-12 -1.74190E-12 -1.73350E-12 -1.72330E-12 -1.71170E-12 -1.69900E-12 -1.68560E-12 -1.67190E-12 &
     -1.65830E-12 -1.64490E-12 -1.63220E-12 -1.62060E-12 -1.61020E-12 -1.60140E-12 -1.59450E-12 -1.58860E-12 &
     -1.58320E-12 -1.57720E-12 -1.57000E-12 -1.56070E-12 -1.54850E-12 -1.53290E-12 -1.51440E-12 -1.49380E-12 &
     -1.47200E-12 -1.44980E-12 -1.42800E-12 -1.40730E-12 -1.38860E-12 -1.37130E-12 -1.35500E-12 -1.33930E-12 &
     -1.32340E-12 -1.30710E-12 -1.28960E-12 -1.27070E-12 -1.25050E-12 -1.22930E-12 -1.20740E-12 -1.18510E-12 &
     -1.16270E-12 -1.14040E-12 -1.11850E-12 -1.09710E-12 -1.07620E-12 -1.05560E-12 -1.03540E-12 -1.01560E-12 &
     -9.96150E-13 -9.76960E-13 -9.57780E-13 -9.38280E-13 -9.18130E-13 -8.97000E-13 -8.74570E-13 -8.50520E-13 &
     -8.24680E-13 -7.97600E-13 -7.69960E-13 -7.42470E-13 -7.15830E-13 -6.90750E-13 -6.67920E-13 -6.47780E-13 &
     -6.29780E-13 -6.13070E-13 -5.96830E-13 -5.80230E-13 -5.62460E-13 -5.42670E-13 -5.20260E-13 -4.95520E-13 &
     -4.68940E-13 -4.41030E-13 -4.12280E-13 -3.83190E-13 -3.54270E-13 -3.25960E-13 -2.98460E-13 -2.71920E-13 &
     -2.46470E-13 -2.22270E-13 -1.99460E-13 -1.78170E-13 -1.58430E-13 -1.39760E-13 -1.21550E-13 -1.03190E-13 &
     -8.40620E-14 -6.35700E-14 -4.11020E-14 -1.62690E-14 1.04320E-14 3.82860E-14 6.65770E-14 9.45870E-14 &
     1.21600E-13 1.46900E-13 1.70000E-13 1.91280E-13 2.11370E-13 2.30870E-13 2.50400E-13 2.70580E-13 &
     2.92030E-13 3.15210E-13 3.40040E-13 3.66290E-13 3.93740E-13 4.22160E-13 4.51320E-13 4.81000E-13 &
     5.10950E-13 5.40850E-13 5.70390E-13 5.99240E-13 6.27060E-13 6.53540E-13 6.78340E-13 7.01250E-13 &
     7.22480E-13 7.42380E-13 7.61270E-13 7.79490E-13 7.97370E-13 8.15240E-13 8.33380E-13 8.51810E-13 &
     8.70520E-13 8.89470E-13 9.08640E-13 9.28000E-13 9.47520E-13 9.67150E-13 9.86680E-13 1.00590E-12 &
     1.02450E-12 1.04240E-12 1.05920E-12 1.07470E-12 1.08880E-12 1.10190E-12 1.11440E-12 1.12660E-12 &
     1.13910E-12 1.15230E-12 1.16660E-12 1.18220E-12 1.19900E-12 1.21630E-12 1.23390E-12 1.25110E-12 &
     1.26760E-12 1.28290E-12 1.29670E-12 1.30910E-12 1.32040E-12 1.33100E-12 1.34100E-12 1.35090E-12 &
     1.36080E-12 1.37110E-12 1.38160E-12 1.39240E-12 1.40330E-12 1.41420E-12 1.42510E-12 1.43590E-12 &
     1.44660E-12 1.45690E-12 1.46690E-12 1.47640E-12 1.48540E-12 1.49380E-12 1.50150E-12 1.50850E-12 &
     1.51490E-12 1.52090E-12 1.52660E-12 1.53240E-12 1.53820E-12 1.54440E-12 1.55100E-12 1.55800E-12 &
     1.56500E-12 1.57200E-12 1.57880E-12 1.58520E-12 1.59090E-12 1.59580E-12 1.60010E-12 1.60380E-12 &
     1.60700E-12 1.60970E-12 1.61230E-12 1.61460E-12 1.61680E-12 1.61910E-12 1.62130E-12 1.62360E-12 &
     1.62610E-12 1.62880E-12 1.63170E-12 1.63480E-12 1.63800E-12 1.64090E-12 1.64330E-12 1.64500E-12 &
     1.64560E-12 1.64500E-12 1.64300E-12 1.63970E-12 1.63560E-12 1.63110E-12 1.62650E-12 1.62220E-12 &
     1.61860E-12 1.61580E-12 1.61380E-12 1.61240E-12 1.61130E-12 1.61040E-12 1.60930E-12 1.60780E-12 &
     1.60580E-12 1.60320E-12 1.60000E-12 1.59600E-12 1.59130E-12 1.58580E-12 1.57950E-12 1.57230E-12 &
     1.56450E-12 1.55620E-12 1.54780E-12 1.53930E-12 1.53100E-12 1.52320E-12 1.51590E-12 1.50910E-12 &
     1.50250E-12 1.49580E-12 1.48900E-12 1.48160E-12 1.47360E-12 1.46480E-12 1.45510E-12 1.44490E-12 &
     1.43420E-12 1.42310E-12 1.41190E-12 1.40060E-12 1.38930E-12 1.37820E-12 1.36710E-12 1.35600E-12 &
     1.34490E-12 1.33380E-12 1.32270E-12 1.31150E-12 1.30020E-12 1.28880E-12 1.27710E-12 1.26510E-12 &
     1.25270E-12 1.24000E-12 1.22680E-12 1.21310E-12 1.19910E-12 1.18480E-12 1.17020E-12 1.15530E-12 &
     1.14030E-12 1.12520E-12 1.11020E-12 1.09530E-12 1.08090E-12 1.06710E-12 1.05410E-12 1.04200E-12 &
     1.03100E-12 1.02090E-12 1.01140E-12 1.00230E-12 9.93360E-13 9.84260E-13 9.74780E-13 9.64750E-13 &
     9.54300E-13 9.43650E-13 9.32990E-13 9.22550E-13 9.12530E-13 9.03140E-13 8.94490E-13 8.86290E-13 &
     8.78160E-13 8.69690E-13 8.60510E-13 8.50210E-13 8.38420E-13 8.24850E-13 8.09780E-13 7.93580E-13 &
     7.76620E-13 7.59290E-13 7.41970E-13 7.25030E-13 7.08780E-13 6.93170E-13 6.78050E-13 6.63290E-13 &
     6.48770E-13 6.34340E-13 6.19870E-13 6.05250E-13 5.90470E-13 5.75530E-13 5.60440E-13 5.45210E-13 &
     5.29850E-13 5.14370E-13 4.98780E-13 4.83150E-13 4.67520E-13 4.51960E-13 4.36540E-13 4.21320E-13 &
     4.06350E-13 3.91690E-13 3.77300E-13 3.63130E-13 3.49120E-13 3.35240E-13 3.21430E-13 3.07630E-13 &
     2.93810E-13 2.79910E-13 2.65900E-13 2.51740E-13 2.37390E-13 2.22800E-13 2.07950E-13 1.92850E-13 &
     1.77740E-13 1.62930E-13 1.48720E-13 1.35430E-13 1.23340E-13 1.12780E-13 1.03920E-13 9.65240E-14 &
     9.02260E-14 8.46640E-14 7.94760E-14 7.43010E-14 6.87760E-14 6.26280E-14 5.59440E-14 4.88980E-14 &
     4.16640E-14 3.44180E-14 2.73340E-14 2.05880E-14 1.43050E-14 8.41920E-15 2.81540E-15 -2.62120E-15 &
     -8.00550E-15 -1.34530E-14 -1.90770E-14 -2.49750E-14 -3.11660E-14 -3.76490E-14 -4.44250E-14 -5.14940E-14 &
     -5.88560E-14 -6.65110E-14 -7.44410E-14 -8.25590E-14 -9.07590E-14 -9.89380E-14 -1.06990E-13 -1.14810E-13 &
     -1.22290E-13 -1.29350E-13 -1.36050E-13 -1.42430E-13 -1.48580E-13 -1.54560E-13 -1.60440E-13 -1.66280E-13 &
     -1.72140E-13 -1.77980E-13 -1.83750E-13 -1.89410E-13 -1.94900E-13 -2.00170E-13 -2.05170E-13 -2.09860E-13 &
     -2.14230E-13 -2.18300E-13 -2.22080E-13 -2.25560E-13 -2.28760E-13 -2.31690E-13 -2.34360E-13 -2.36820E-13 &
     -2.39120E-13 -2.41340E-13 -2.43510E-13 -2.45700E-13 -2.47970E-13 -2.50340E-13 -2.52670E-13 -2.54800E-13 &
     -2.56580E-13 -2.57830E-13 -2.58390E-13 -2.58090E-13 -2.56830E-13 -2.54750E-13 -2.52050E-13 -2.48920E-13 &
     -2.45560E-13 -2.42180E-13 -2.38960E-13 -2.36050E-13 -2.33440E-13 -2.31020E-13 -2.28710E-13 -2.26450E-13 &
     -2.24130E-13 -2.21680E-13 -2.19040E-13 -2.16200E-13 -2.13180E-13 -2.10010E-13 -2.06700E-13 -2.03260E-13 &
     -1.99730E-13 -1.96110E-13 -1.92440E-13 -1.88730E-13 -1.85010E-13 -1.81310E-13 -1.77640E-13 -1.74050E-13 &
     -1.70520E-13 -1.67010E-13 -1.63430E-13 -1.59720E-13 -1.55790E-13 -1.51580E-13 -1.46990E-13 -1.41990E-13 &
     -1.36600E-13 -1.30870E-13 -1.24870E-13 -1.18650E-13 -1.12250E-13 -1.05750E-13 -9.91750E-14 -9.25380E-14 &
     -8.58310E-14 -7.90490E-14 -7.21870E-14 -6.52390E-14 -5.82010E-14 -5.10630E-14 -4.38090E-14 -3.64190E-14 &
     -2.88720E-14 -2.11500E-14 -1.32310E-14 -5.09660E-15 3.26490E-15 1.18300E-14 2.05680E-14 2.94460E-14 &
     3.84320E-14 4.74950E-14 5.66030E-14 6.57250E-14 7.48410E-14 8.39290E-14 9.29700E-14 1.01940E-13 &
     1.10830E-13 1.19610E-13 1.28260E-13 1.36810E-13 1.45260E-13 1.53640E-13 1.61970E-13 1.70250E-13 &
     1.78510E-13 1.86760E-13 1.95000E-13 2.03210E-13 2.11390E-13 2.19510E-13 2.27580E-13 2.35570E-13 &
     2.43490E-13 2.51360E-13 2.59210E-13 2.67070E-13 2.74990E-13 2.82990E-13 2.91110E-13 2.99360E-13 &
     3.07710E-13 3.16130E-13 3.24560E-13 3.32970E-13 3.41300E-13 3.49510E-13 3.57580E-13 3.65550E-13 &
     3.73480E-13 3.81430E-13 3.89460E-13 3.97640E-13 4.06020E-13 4.14660E-13 4.23550E-13 4.32670E-13 &
     4.42020E-13 4.51580E-13 4.61320E-13 4.71240E-13 4.81320E-13 4.91520E-13 5.01780E-13 5.12070E-13 &
     5.22340E-13 5.32550E-13 5.42640E-13 5.52590E-13 5.62410E-13 5.72130E-13 5.81780E-13 5.91380E-13 &
     6.00970E-13 6.10580E-13 6.20210E-13 6.29880E-13 6.39560E-13 6.49250E-13 6.58940E-13 6.68610E-13 &
     6.78250E-13 6.87830E-13 6.97290E-13 7.06530E-13 7.15460E-13 7.24000E-13 7.32050E-13 7.39530E-13 &
     7.46380E-13 7.52680E-13 7.58550E-13 7.64110E-13 7.69460E-13 7.74720E-13 7.80020E-13 7.85430E-13 &
     7.90940E-13 7.96500E-13 8.02070E-13 8.07590E-13 8.13020E-13 8.18320E-13 8.23450E-13 8.28420E-13 &
     8.33270E-13 8.38030E-13 8.42710E-13 8.47360E-13 8.52000E-13 8.56650E-13 8.61310E-13 8.65970E-13 &
     8.70610E-13 8.75240E-13 8.79830E-13 8.84380E-13 8.88880E-13 8.93320E-13 8.97690E-13 9.01980E-13 &
     9.06180E-13 9.10280E-13 9.14270E-13 9.18150E-13 9.21950E-13 9.25690E-13 9.29410E-13 9.33140E-13 &
     9.36920E-13 9.40780E-13 9.44730E-13 9.48740E-13 9.52780E-13 9.56790E-13 9.60730E-13 9.64550E-13 &
     9.68210E-13 9.71680E-13 9.74960E-13 9.78060E-13 9.80980E-13 9.83750E-13 9.86370E-13 9.88850E-13 &
     9.91210E-13 9.93450E-13 9.95590E-13 9.97630E-13 9.99570E-13 1.00140E-12 1.00320E-12 1.00490E-12 &
     1.00650E-12 1.00790E-12 1.00900E-12 1.00980E-12 1.01020E-12 1.01020E-12 1.00960E-12 1.00860E-12 &
     1.00720E-12 1.00560E-12 1.00390E-12 1.00200E-12 1.00020E-12 9.98490E-13 9.96840E-13 9.95230E-13 &
     9.93640E-13 9.92030E-13 9.90360E-13 9.88610E-13 9.86750E-13 9.84770E-13 9.82700E-13 9.80540E-13 &
     9.78300E-13 9.75990E-13 9.73620E-13 9.71210E-13 9.68750E-13 9.66250E-13 9.63710E-13 9.61140E-13 &
     9.58530E-13 9.55900E-13 9.53230E-13 9.50530E-13 9.47770E-13 9.44950E-13 9.42060E-13 9.39070E-13 &
     9.35980E-13 9.32780E-13 9.29480E-13 9.26110E-13 9.22700E-13 9.19260E-13 9.15810E-13 9.12390E-13 &
     9.08990E-13 9.05600E-13 9.02150E-13 8.98610E-13 8.94930E-13 8.91070E-13 8.86980E-13 8.82620E-13 &
     8.78020E-13 8.73200E-13 8.68180E-13 8.62980E-13 8.57650E-13 8.52190E-13 8.46640E-13 8.41000E-13 &
     8.35270E-13 8.29470E-13 8.23590E-13 8.17650E-13 8.11630E-13 8.05560E-13 7.99420E-13 7.93220E-13 &
     7.86950E-13 7.80600E-13 7.74190E-13 7.67700E-13 7.61130E-13 7.54500E-13 7.47810E-13 7.41100E-13 &
     7.34350E-13 7.27600E-13 7.20860E-13 7.14130E-13 7.07420E-13 7.00720E-13 6.94030E-13 6.87350E-13 &
     6.80670E-13 6.73990E-13 6.67320E-13 6.60630E-13 6.53940E-13 6.47230E-13 6.40510E-13 6.33750E-13 &
     6.26970E-13 6.20160E-13 6.13330E-13 6.06480E-13 5.99620E-13 5.92760E-13 5.85920E-13 5.79090E-13 &
     5.72280E-13 5.65500E-13 5.58720E-13 5.51950E-13 5.45160E-13 5.38360E-13 5.31540E-13 5.24680E-13 &
     5.17790E-13 5.10870E-13 5.03930E-13 4.96970E-13 4.89990E-13 4.82990E-13 4.75990E-13 4.68980E-13 &
     4.61970E-13 4.54990E-13 4.48030E-13 4.41100E-13 4.34230E-13 4.27390E-13 4.20570E-13 4.13720E-13 &
     4.06800E-13 3.99750E-13 3.92540E-13 3.85130E-13 3.77480E-13 3.69630E-13 3.61640E-13 3.53540E-13 &
     3.45410E-13 3.37270E-13 3.29190E-13 3.21210E-13 3.13350E-13 3.05610E-13 2.98000E-13 2.90540E-13 &
     2.83240E-13 2.76100E-13 2.69120E-13 2.62300E-13 2.55600E-13 2.49000E-13 2.42470E-13 2.36000E-13 &
     2.29560E-13 2.23130E-13 2.16710E-13 2.10310E-13 2.03950E-13 1.97630E-13 1.91360E-13 1.85150E-13 &
     1.79010E-13 1.72930E-13 1.66920E-13 1.60950E-13 1.55030E-13 1.49150E-13 1.43290E-13 1.37470E-13 &
     1.31670E-13 1.25900E-13 1.20180E-13 1.14490E-13 1.08850E-13 1.03270E-13 9.77340E-14 9.22560E-14 &
     8.68310E-14 8.14590E-14 7.61390E-14 7.08710E-14 6.56530E-14 6.04850E-14 5.53630E-14 5.02850E-14 &
     4.52470E-14 4.02470E-14 3.52810E-14 3.03470E-14 2.54420E-14 2.05710E-14 1.57410E-14 1.09590E-14 &
     6.22900E-15 1.55900E-15 -3.04510E-15 -7.57930E-15 -1.20480E-14 -1.64600E-14 -2.08210E-14 -2.51380E-14 &
     -2.94190E-14 -3.36710E-14 -3.78990E-14 -4.20930E-14 -4.62440E-14 -5.03380E-14 -5.43670E-14 -5.83180E-14 &
     -6.21790E-14 -6.59410E-14 -6.95900E-14 -7.31160E-14 -7.65060E-14 -7.97480E-14 -8.28300E-14 -8.57410E-14 &
     -8.84750E-14 -9.10470E-14 -9.34830E-14 -9.58040E-14 -9.80330E-14 -1.00190E-13 -1.02310E-13 -1.04400E-13 &
     -1.06450E-13 -1.08470E-13 -1.10440E-13 -1.12360E-13 -1.14200E-13 -1.15980E-13 -1.17680E-13 -1.19300E-13 &
     -1.20850E-13 -1.22320E-13 -1.23730E-13 -1.25080E-13 -1.26370E-13 -1.27610E-13 -1.28790E-13 -1.29910E-13 &
     -1.30990E-13 -1.32010E-13 -1.32970E-13 -1.33890E-13 -1.34750E-13 -1.35550E-13 -1.36300E-13 -1.36980E-13 &
     -1.37600E-13 -1.38160E-13 -1.38640E-13 -1.39060E-13 -1.39410E-13 -1.39690E-13 -1.39920E-13 -1.40100E-13 &
     -1.40220E-13 -1.40300E-13 -1.40350E-13 -1.40340E-13 -1.40290E-13 -1.40190E-13 -1.40030E-13 -1.39800E-13 &
     -1.39510E-13 -1.39150E-13 -1.38720E-13 -1.38230E-13 -1.37680E-13 -1.37080E-13 -1.36430E-13 -1.35730E-13 &
     -1.35000E-13 -1.34210E-13 -1.33370E-13 -1.32450E-13 -1.31460E-13 -1.30380E-13 -1.29200E-13 -1.27910E-13 &
     -1.26520E-13 -1.25010E-13 -1.23400E-13 -1.21690E-13 -1.19870E-13 -1.17950E-13 -1.15930E-13 -1.13810E-13 &
     -1.11630E-13 -1.09370E-13 -1.07070E-13 -1.04730E-13 -1.02370E-13 -9.99930E-14 -9.76010E-14 -9.51920E-14 &
     -9.27610E-14 -9.03070E-14 -8.78250E-14 -8.53130E-14 -8.27680E-14 -8.01890E-14 -7.75780E-14 -7.49360E-14 &
     -7.22620E-14 -6.95590E-14 -6.68250E-14 -6.40630E-14 -6.12740E-14 -5.84590E-14 -5.56210E-14 -5.27610E-14 &
     -4.98800E-14 -4.69820E-14 -4.40660E-14 -4.11310E-14 -3.81770E-14 -3.52020E-14 -3.22040E-14 -2.91830E-14 &
     -2.61350E-14 -2.30610E-14 -1.99620E-14 -1.68380E-14 -1.36930E-14 -1.05260E-14 -7.34120E-15 -4.13860E-15 &
     -9.19800E-16 2.31500E-15 5.56620E-15 8.83400E-15 1.21190E-14 1.54210E-14 1.87400E-14 2.20770E-14 &
     2.54330E-14 2.88100E-14 3.22100E-14 3.56330E-14 3.90810E-14 4.25550E-14 4.60560E-14 4.95750E-14 &
     5.31000E-14 5.66230E-14 6.01320E-14 6.36170E-14 6.70670E-14 7.04750E-14 7.38410E-14 7.71680E-14 &
     8.04580E-14 8.37150E-14 8.69420E-14 9.01410E-14 9.33170E-14 9.64800E-14 9.96430E-14 1.02820E-13 &
     1.06020E-13 1.09250E-13 1.12540E-13 1.15880E-13 1.19260E-13 1.22680E-13 1.26110E-13 1.29540E-13 &
     1.32960E-13 1.36350E-13 1.39690E-13 1.42990E-13 1.46260E-13 1.49490E-13 1.52700E-13 1.55890E-13 &
     1.59060E-13 1.62220E-13 1.65370E-13 1.68500E-13 1.71610E-13 1.74690E-13 1.77750E-13 1.80780E-13 &
     1.83780E-13 1.86760E-13 1.89700E-13 1.92620E-13 1.95510E-13 1.98380E-13 2.01230E-13 2.04070E-13 &
     2.06880E-13 2.09670E-13 2.12440E-13 2.15180E-13 2.17900E-13 2.20590E-13 2.23250E-13 2.25890E-13 &
     2.28490E-13 2.31060E-13 2.33600E-13 2.36100E-13 2.38580E-13 2.41020E-13 2.43430E-13 2.45810E-13 &
     2.48160E-13 2.50490E-13 2.52790E-13 2.55070E-13 2.57320E-13 2.59550E-13 2.61740E-13 2.63890E-13 &
     2.66000E-13 2.68050E-13 2.70030E-13 2.71950E-13 2.73790E-13 2.75550E-13 2.77220E-13 2.78800E-13 &
     2.80280E-13 2.81650E-13 2.82910E-13 2.84080E-13 2.85170E-13 2.86210E-13 2.87210E-13 2.88200E-13 &
     2.89190E-13 2.90190E-13 2.91210E-13 2.92230E-13 2.93240E-13 2.94230E-13 2.95180E-13 2.96090E-13 &
     2.96950E-13 2.97760E-13 2.98510E-13 2.99220E-13 2.99890E-13 3.00510E-13 3.01100E-13 3.01650E-13 &
     3.02170E-13 3.02650E-13 3.03090E-13 3.03490E-13 3.03840E-13 3.04160E-13 3.04430E-13 3.04650E-13 &
     3.04840E-13 3.04980E-13 3.05080E-13 3.05140E-13 3.05160E-13 3.05140E-13 3.05080E-13 3.04990E-13 &
     3.04860E-13 3.04690E-13 3.04480E-13 3.04240E-13 3.03960E-13 3.03650E-13 3.03300E-13 3.02900E-13 &
     3.02470E-13 3.02000E-13 3.01480E-13 3.00930E-13 3.00320E-13 2.99680E-13 2.99000E-13 2.98270E-13 &
     2.97510E-13 2.96710E-13 2.95870E-13 2.95000E-13 2.94080E-13 2.93130E-13 2.92140E-13 2.91110E-13 &
     2.90040E-13 2.88930E-13 2.87780E-13 2.86590E-13 2.85360E-13 2.84110E-13 2.82820E-13 2.81500E-13 &
     2.80150E-13 2.78790E-13 2.77400E-13 2.76010E-13 2.74620E-13 2.73240E-13 2.71860E-13 2.70500E-13 &
     2.69130E-13 2.67750E-13 2.66330E-13 2.64850E-13 2.63290E-13 2.61640E-13 2.59880E-13 2.58020E-13 &
     2.56070E-13 2.54060E-13 2.52000E-13 2.49910E-13 2.47800E-13 2.45680E-13 2.43550E-13 2.41420E-13 &
     2.39270E-13 2.37110E-13 2.34930E-13 2.32730E-13 2.30510E-13 2.28270E-13 2.26000E-13 2.23710E-13 &
     2.21410E-13 2.19080E-13 2.16730E-13 2.14370E-13 2.11990E-13 2.09590E-13 2.07180E-13 2.04750E-13 &
     2.02310E-13 1.99860E-13 1.97390E-13 1.94910E-13 1.92430E-13 1.89920E-13 1.87410E-13 1.84880E-13 &
     1.82340E-13 1.79790E-13 1.77220E-13 1.74640E-13 1.72040E-13 1.69420E-13 1.66790E-13 1.64140E-13 &
     1.61480E-13 1.58800E-13 1.56120E-13 1.53440E-13 1.50770E-13 1.48110E-13 1.45470E-13 1.42860E-13 &
     1.40270E-13 1.37720E-13 1.35210E-13 1.32740E-13 1.30320E-13 1.27950E-13 1.25630E-13 1.23360E-13 &
     1.21130E-13 1.18930E-13 1.16740E-13 1.14560E-13 1.12370E-13 1.10180E-13 1.07970E-13 1.05770E-13 &
     1.03560E-13 1.01360E-13 9.91730E-14 9.70030E-14 9.48500E-14 9.27050E-14 9.05510E-14 8.83760E-14 &
     8.61640E-14 8.39000E-14 8.15710E-14 7.91660E-14 7.66980E-14 7.41830E-14 7.16380E-14 6.90810E-14 &
     6.65280E-14 6.39970E-14 6.15000E-14 5.90370E-14 5.66040E-14 5.41970E-14 5.18100E-14 4.94410E-14 &
     4.70840E-14 4.47360E-14 4.23990E-14 4.00730E-14 3.77600E-14 3.54620E-14 3.31790E-14 3.09150E-14 &
     2.86700E-14 2.64450E-14 2.42400E-14 2.20550E-14 1.98910E-14 1.77480E-14 1.56270E-14 1.35270E-14 &
     1.14480E-14 9.38540E-15 7.33880E-15 5.30570E-15 3.28390E-15 1.27140E-15 -7.32050E-16 -2.71910E-15 &
     -4.68070E-15 -6.60740E-15 -8.49020E-15 -1.03200E-14 -1.20870E-14 -1.37840E-14 -1.54110E-14 -1.69670E-14 &
     -1.84540E-14 -1.98730E-14 -2.12240E-14 -2.25070E-14 -2.37260E-14 -2.48880E-14 -2.60020E-14 -2.70770E-14 &
     -2.81220E-14 -2.91470E-14 -3.01600E-14 -3.11680E-14 -3.21710E-14 -3.31660E-14 -3.41510E-14 -3.51220E-14 &
     -3.60790E-14 -3.70170E-14 -3.79360E-14 -3.88360E-14 -3.97180E-14 -4.05850E-14 -4.14380E-14 -4.22780E-14 &
     -4.31070E-14 -4.39280E-14 -4.47370E-14 -4.55350E-14 -4.63190E-14 -4.70880E-14 -4.78400E-14 -4.85750E-14 &
     -4.92900E-14 -4.99840E-14 -5.06560E-14 -5.13030E-14 -5.19230E-14 -5.25150E-14 -5.30780E-14 -5.36090E-14 &
     -5.41090E-14 -5.45790E-14 -5.50180E-14 -5.54280E-14 -5.58090E-14 -5.61620E-14 -5.64860E-14 -5.67820E-14 &
     -5.70500E-14 -5.72880E-14 -5.74980E-14 -5.76780E-14 -5.78280E-14 -5.79490E-14 -5.80400E-14 -5.81040E-14 &
     -5.81420E-14 -5.81540E-14 -5.81420E-14 -5.81080E-14 -5.80510E-14 -5.79670E-14 -5.78530E-14 -5.77050E-14 &
     -5.75160E-14 -5.72840E-14 -5.70040E-14 -5.66720E-14 -5.62900E-14 -5.58610E-14 -5.53880E-14 -5.48740E-14 &
     -5.43230E-14 -5.37360E-14 -5.31170E-14 -5.24690E-14 -5.17970E-14 -5.11040E-14 -5.03930E-14 -4.96690E-14 &
     -4.89350E-14 -4.81940E-14 -4.74460E-14 -4.66890E-14 -4.59220E-14 -4.51430E-14 -4.43510E-14 -4.35450E-14 &
     -4.27230E-14 -4.18880E-14 -4.10380E-14 -4.01770E-14 -3.93040E-14 -3.84220E-14 -3.75300E-14 -3.66300E-14 &
     -3.57210E-14 -3.47980E-14 -3.38610E-14 -3.29060E-14 -3.19310E-14 -3.09340E-14 -2.99130E-14 -2.88640E-14 &
     -2.77870E-14 -2.66800E-14 -2.55400E-14 -2.43660E-14 -2.31550E-14 -2.19080E-14 -2.06260E-14 -1.93140E-14 &
     -1.79770E-14 -1.66180E-14 -1.52420E-14 -1.38520E-14 -1.24530E-14 -1.10440E-14 -9.62470E-15 -8.19500E-15 &
     -6.75440E-15 -5.30230E-15 -3.83830E-15 -2.36200E-15 -8.73100E-16 6.28700E-16 2.14360E-15 3.67200E-15 &
     5.21390E-15 6.76980E-15 8.33940E-15 9.92090E-15 1.15120E-14 1.31110E-14 1.47140E-14 1.63210E-14 &
     1.79280E-14 1.95330E-14 2.11360E-14 2.27360E-14 2.43310E-14 2.59220E-14 2.75070E-14 2.90850E-14 &
     3.06560E-14 3.22200E-14 3.37790E-14 3.53320E-14 3.68810E-14 3.84270E-14 3.99700E-14 4.15100E-14 &
     4.30490E-14 4.45860E-14 4.61210E-14 4.76530E-14 4.91840E-14 5.07120E-14 5.22380E-14 5.37620E-14 &
     5.52850E-14 5.68060E-14 5.83260E-14 5.98460E-14 6.13650E-14 6.28830E-14 6.44010E-14 6.59160E-14 &
     6.74290E-14 6.89380E-14 7.04430E-14 7.19420E-14 7.34350E-14 7.49250E-14 7.64140E-14 7.79060E-14 &
     7.94030E-14 8.09090E-14 8.24250E-14 8.39540E-14 8.54950E-14 8.70460E-14 8.86040E-14 9.01680E-14 &
     9.17340E-14 9.33010E-14 9.48670E-14 9.64290E-14 9.79860E-14 9.95360E-14 1.01080E-13 1.02610E-13 &
     1.04120E-13 1.05630E-13 1.07120E-13 1.08600E-13 1.10060E-13 1.11530E-13 1.12980E-13 1.14430E-13 &
     1.15870E-13 1.17300E-13 1.18720E-13 1.20120E-13 1.21500E-13 1.22840E-13 1.24140E-13 1.25400E-13 &
     1.26620E-13 1.27810E-13 1.28960E-13 1.30070E-13 1.31150E-13 1.32200E-13 1.33220E-13 1.34220E-13 &
     1.35200E-13 1.36150E-13 1.37090E-13 1.38020E-13 1.38930E-13 1.39840E-13 1.40730E-13 1.41620E-13 &
     1.42490E-13 1.43360E-13 1.44210E-13 1.45040E-13 1.45870E-13 1.46670E-13 1.47470E-13 1.48250E-13 &
     1.49010E-13 1.49760E-13 1.50500E-13 1.51220E-13 1.51930E-13 1.52620E-13 1.53300E-13 1.53970E-13 &
     1.54620E-13 1.55260E-13 1.55890E-13 1.56500E-13 1.57100E-13 1.57690E-13 1.58270E-13 1.58840E-13 &
     1.59390E-13 1.59940E-13 1.60470E-13 1.60990E-13 1.61490E-13 1.61970E-13 1.62430E-13 1.62860E-13 &
     1.63260E-13 1.63640E-13 1.64000E-13 1.64330E-13 1.64630E-13 1.64910E-13 1.65160E-13 1.65390E-13 &
     1.65600E-13 1.65780E-13 1.65940E-13 1.66090E-13 1.66210E-13 1.66320E-13 1.66410E-13 1.66480E-13 &
     1.66530E-13 1.66550E-13 1.66540E-13 1.66500E-13 1.66430E-13 1.66310E-13 1.66160E-13 1.65980E-13 &
     1.65760E-13 1.65520E-13 1.65260E-13 1.64970E-13 1.64670E-13 1.64350E-13 1.64020E-13 1.63670E-13 &
     1.63300E-13 1.62930E-13 1.62540E-13 1.62140E-13 1.61730E-13 1.61310E-13 1.60870E-13 1.60430E-13 &
     1.59970E-13 1.59500E-13 1.59020E-13 1.58520E-13 1.58010E-13 1.57490E-13 1.56960E-13 1.56410E-13 &
     1.55850E-13 1.55290E-13 1.54710E-13 1.54120E-13 1.53520E-13 1.52900E-13 1.52280E-13 1.51650E-13 &
     1.51010E-13 1.50360E-13 1.49690E-13 1.49020E-13 1.48330E-13 1.47630E-13 1.46910E-13 1.46180E-13 &
     1.45440E-13 1.44690E-13 1.43930E-13 1.43170E-13 1.42400E-13 1.41630E-13 1.40860E-13 1.40080E-13 &
     1.39290E-13 1.38470E-13 1.37630E-13 1.36740E-13 1.35810E-13 1.34830E-13 1.33800E-13 1.32730E-13 &
     1.31630E-13 1.30500E-13 1.29360E-13 1.28200E-13 1.27030E-13 1.25860E-13 1.24690E-13 1.23510E-13 &
     1.22330E-13 1.21150E-13 1.19960E-13 1.18780E-13 1.17590E-13 1.16400E-13 1.15200E-13 1.14010E-13 &
     1.12820E-13 1.11630E-13 1.10450E-13 1.09260E-13 1.08080E-13 1.06900E-13 1.05720E-13 1.04540E-13 &
     1.03350E-13 1.02170E-13 1.00980E-13 9.97950E-14 9.86070E-14 9.74190E-14 9.62300E-14 9.50410E-14 &
     9.38520E-14 9.26630E-14 9.14740E-14 9.02850E-14 8.90950E-14 8.79060E-14 8.67160E-14 8.55260E-14 &
     8.43360E-14 8.31460E-14 8.19570E-14 8.07670E-14 7.95780E-14 7.83900E-14 7.72020E-14 7.60140E-14 &
     7.48280E-14 7.36420E-14 7.24560E-14 7.12710E-14 7.00870E-14 6.89040E-14 6.77220E-14 6.65400E-14 &
     6.53580E-14 6.41770E-14 6.29970E-14 6.18170E-14 6.06370E-14 5.94590E-14 5.82820E-14 5.71070E-14 &
     5.59360E-14 5.47690E-14 5.36060E-14 5.24470E-14 5.12890E-14 5.01270E-14 4.89570E-14 4.77750E-14 &
     4.65750E-14 4.53530E-14 4.41080E-14 4.28470E-14 4.15790E-14 4.03140E-14 3.90620E-14 3.78320E-14 &
     3.66350E-14 3.54770E-14 3.43570E-14 3.32700E-14 3.22110E-14 3.11760E-14 3.01610E-14 2.91620E-14 &
     2.81730E-14 2.71960E-14 2.62290E-14 2.52720E-14 2.43260E-14 2.33900E-14 2.24650E-14 2.15500E-14 &
     2.06440E-14 1.97490E-14 1.88630E-14 1.79860E-14 1.71180E-14 1.62580E-14 1.54070E-14 1.45650E-14 &
     1.37320E-14 1.29080E-14 1.20940E-14 1.12890E-14 1.04950E-14 9.71150E-15 8.93810E-15 8.17470E-15 &
     7.42130E-15 6.67770E-15 5.94350E-15 5.21870E-15 4.50310E-15 3.79670E-15 3.09980E-15 2.41250E-15 &
     1.73480E-15 1.06710E-15 4.09290E-16 -2.38400E-16 -8.76160E-16 -1.50420E-15 -2.12280E-15 -2.73210E-15 &
     -3.33230E-15 -3.92380E-15 -4.50650E-15 -5.08030E-15 -5.64480E-15 -6.19980E-15 -6.74470E-15 -7.27950E-15 &
     -7.80360E-15 -8.31700E-15 -8.82000E-15 -9.31330E-15 -9.79740E-15 -1.02730E-14 -1.07410E-14 -1.12010E-14 &
     -1.16540E-14 -1.20950E-14 -1.25200E-14 -1.29240E-14 -1.33040E-14 -1.36530E-14 -1.39680E-14 -1.42440E-14 &
     -1.44860E-14 -1.46970E-14 -1.48820E-14 -1.50440E-14 -1.51880E-14 -1.53180E-14 -1.54380E-14 -1.55480E-14 &
     -1.56480E-14 -1.57380E-14 -1.58180E-14 -1.58890E-14 -1.59490E-14 -1.59990E-14 -1.60390E-14 -1.60680E-14 &
     -1.60890E-14 -1.60990E-14 -1.60990E-14 -1.60900E-14 -1.60710E-14 -1.60420E-14 -1.60040E-14 -1.59560E-14 &
     -1.58990E-14 -1.58330E-14 -1.57570E-14 -1.56720E-14 -1.55770E-14 -1.54740E-14 -1.53610E-14 -1.52390E-14 &
     -1.51080E-14 -1.49680E-14 -1.48190E-14 -1.46610E-14 -1.44940E-14 -1.43180E-14 -1.41340E-14 -1.39400E-14 &
     -1.37380E-14 -1.35260E-14 -1.33060E-14 -1.30770E-14 -1.28400E-14 -1.25940E-14 -1.23390E-14 -1.20760E-14 &
     -1.18050E-14 -1.15250E-14 -1.12370E-14 -1.09410E-14 -1.06360E-14 -1.03220E-14 -1.00000E-14 -9.67010E-15 &
     -9.33140E-15 -8.98470E-15 -8.63010E-15 -8.26800E-15 -7.89860E-15 -7.52200E-15 -7.13850E-15 -6.74720E-15 &
     -6.34730E-15 -5.93790E-15 -5.51800E-15 -5.08690E-15 -4.64350E-15 -4.18710E-15 -3.71740E-15 -3.23420E-15 &
     -2.73720E-15 -2.22620E-15 -1.70100E-15 -1.16130E-15 -6.07310E-16 -4.09540E-17 5.35550E-16 1.11990E-15 &
     1.70990E-15 2.30320E-15 2.89760E-15 3.49130E-15 4.08440E-15 4.67740E-15 5.27110E-15 5.86590E-15 &
     6.46250E-15 7.06160E-15 7.66340E-15 8.26800E-15 8.87500E-15 9.48420E-15 1.00950E-14 1.07080E-14 &
     1.13220E-14 1.19380E-14 1.25540E-14 1.31720E-14 1.37910E-14 1.44120E-14 1.50340E-14 1.56580E-14 &
     1.62830E-14 1.69100E-14 1.75390E-14 1.81690E-14 1.88000E-14 1.94330E-14 2.00660E-14 2.07010E-14 &
     2.13360E-14 2.19720E-14 2.26080E-14 2.32460E-14 2.38840E-14 2.45220E-14 2.51620E-14 2.58010E-14 &
     2.64410E-14 2.70800E-14 2.77200E-14 2.83590E-14 2.89970E-14 2.96340E-14 3.02720E-14 3.09100E-14 &
     3.15490E-14 3.21910E-14 3.28360E-14 3.34850E-14 3.41370E-14 3.47910E-14 3.54410E-14 3.60830E-14 &
     3.67150E-14 3.73320E-14 3.79300E-14 3.85060E-14 3.90640E-14 3.96070E-14 4.01410E-14 4.06690E-14 &
     4.11960E-14 4.17260E-14 4.22620E-14 4.28050E-14 4.33520E-14 4.39020E-14 4.44550E-14 4.50070E-14 &
     4.55580E-14 4.61070E-14 4.66530E-14 4.71960E-14 4.77340E-14 4.82670E-14 4.87960E-14 4.93180E-14 &
     4.98350E-14 5.03450E-14 5.08480E-14 5.13460E-14 5.18380E-14 5.23230E-14 5.28030E-14 5.32760E-14 &
     5.37440E-14 5.42060E-14 5.46620E-14 5.51120E-14 5.55560E-14 5.59940E-14 5.64260E-14 5.68520E-14 &
     5.72720E-14 5.76870E-14 5.80950E-14 5.84960E-14 5.88920E-14 5.92810E-14 5.96640E-14 6.00410E-14 &
     6.04120E-14 6.07760E-14 6.11350E-14 6.14880E-14 6.18340E-14 6.21750E-14 6.25090E-14 6.28360E-14 &
     6.31570E-14 6.34700E-14 6.37760E-14 6.40740E-14 6.43660E-14 6.46510E-14 6.49310E-14 6.52070E-14 &
     6.54790E-14 6.57490E-14 6.60150E-14 6.62770E-14 6.65300E-14 6.67710E-14 6.69960E-14 6.72030E-14 &
     6.73880E-14 6.75490E-14 6.76870E-14 6.78070E-14 6.79120E-14 6.80060E-14 6.80920E-14 6.81730E-14 &
     6.82530E-14 6.83310E-14 6.84060E-14 6.84760E-14 6.85400E-14 6.85970E-14 6.86460E-14 6.86860E-14 &
     6.87160E-14 6.87380E-14 6.87510E-14 6.87550E-14 6.87520E-14 6.87410E-14 6.87220E-14 6.86960E-14 &
     6.86630E-14 6.86220E-14 6.85740E-14 6.85190E-14 6.84560E-14 6.83860E-14 6.83090E-14 6.82240E-14 &
     6.81330E-14 6.80340E-14 6.79290E-14 6.78170E-14 6.76980E-14 6.75720E-14 6.74390E-14 6.73000E-14 &
     6.71540E-14 6.70010E-14 6.68420E-14 6.66750E-14 6.65030E-14 6.63230E-14 6.61370E-14 6.59450E-14 &
     6.57460E-14 6.55420E-14 6.53310E-14 6.51130E-14 6.48900E-14 6.46600E-14 6.44240E-14 6.41820E-14 &
     6.39320E-14 6.36770E-14 6.34150E-14 6.31460E-14 6.28720E-14 6.25910E-14 6.23040E-14 6.20120E-14 &
     6.17130E-14 6.14090E-14 6.10990E-14 6.07820E-14 6.04600E-14 6.01300E-14 5.97940E-14 5.94520E-14 &
     5.91030E-14 5.87480E-14 5.83880E-14 5.80240E-14 5.76560E-14 5.72850E-14 5.69110E-14 5.65340E-14 &
     5.61550E-14 5.57750E-14 5.53940E-14 5.50110E-14 5.46280E-14 5.42440E-14 5.38580E-14 5.34670E-14 &
     5.30700E-14 5.26650E-14 5.22490E-14 5.18200E-14 5.13780E-14 5.09230E-14 5.04580E-14 4.99850E-14 &
     4.95060E-14 4.90240E-14 4.85410E-14 4.80590E-14 4.75770E-14 4.70940E-14 4.66110E-14 4.61270E-14 &
     4.56410E-14 4.51520E-14 4.46600E-14 4.41660E-14 4.36680E-14 4.31690E-14 4.26670E-14 4.21630E-14 &
     4.16580E-14 4.11510E-14 4.06440E-14 4.01350E-14 3.96240E-14 3.91130E-14 3.86010E-14 3.80880E-14 &
     3.75740E-14 3.70600E-14 3.65450E-14 3.60300E-14 3.55140E-14 3.49990E-14 3.44830E-14 3.39680E-14 &
     3.34530E-14 3.29360E-14 3.24180E-14 3.18980E-14 3.13750E-14 3.08490E-14 3.03190E-14 2.97880E-14 &
     2.92580E-14 2.87310E-14 2.82110E-14 2.76990E-14 2.71980E-14 2.67100E-14 2.62330E-14 2.57660E-14 &
     2.53050E-14 2.48500E-14 2.43980E-14 2.39470E-14 2.34950E-14 2.30430E-14 2.25900E-14 2.21380E-14 &
     2.16870E-14 2.12380E-14 2.07910E-14 2.03470E-14 1.99060E-14 1.94670E-14 1.90310E-14 1.85960E-14 &
     1.81640E-14 1.77330E-14 1.73030E-14 1.68740E-14 1.64460E-14 1.60180E-14 1.55900E-14 1.51610E-14 &
     1.47310E-14 1.43000E-14 1.38680E-14 1.34360E-14 1.30060E-14 1.25780E-14 1.21540E-14 1.17350E-14 &
     1.13210E-14 1.09120E-14 1.05080E-14 1.01080E-14 9.71370E-15 9.32320E-15 8.93660E-15 8.55390E-15 &
     8.17490E-15 7.79970E-15 7.42850E-15 7.06130E-15 6.69820E-15 6.33930E-15 5.98460E-15 5.63430E-15 &
     5.28850E-15 4.94740E-15 4.61110E-15 4.27980E-15 3.95370E-15 3.63270E-15 3.31640E-15 3.00420E-15 &
     2.69560E-15 2.38990E-15 2.08660E-15 1.78510E-15 1.48540E-15 1.18960E-15 9.00280E-16 6.20220E-16 &
     3.52060E-16 9.84630E-17 -1.37900E-16 -3.55280E-16 -5.55680E-16 -7.42020E-16 -9.17200E-16 -1.08420E-15 &
     -1.24580E-15 -1.40510E-15 -1.56420E-15 -1.72310E-15 -1.88070E-15 -2.03640E-15 -2.18910E-15 -2.33800E-15 &
     -2.48240E-15 -2.62150E-15 -2.75540E-15 -2.88450E-15 -3.00910E-15 -3.12960E-15 -3.24610E-15 -3.35910E-15 &
     -3.46890E-15 -3.57520E-15 -3.67790E-15 -3.77680E-15 -3.87180E-15 -3.96260E-15 -4.04910E-15 -4.13110E-15 &
     -4.20780E-15 -4.27860E-15 -4.34290E-15 -4.39990E-15 -4.44890E-15 -4.48930E-15 -4.52050E-15 -4.54320E-15 &
     -4.55800E-15 -4.56570E-15 -4.56720E-15 -4.56300E-15 -4.55410E-15 -4.54090E-15 -4.52350E-15 -4.50170E-15 &
     -4.47520E-15 -4.44400E-15 -4.40780E-15 -4.36650E-15 -4.31980E-15 -4.26790E-15 -4.21080E-15 -4.14860E-15 &
     -4.08140E-15 -4.00920E-15 -3.93220E-15 -3.85040E-15 -3.76390E-15 -3.67280E-15 -3.57720E-15 -3.47720E-15 &
     -3.37290E-15 -3.26440E-15 -3.15160E-15 -3.03410E-15 -2.91100E-15 -2.78160E-15 -2.64520E-15 -2.50100E-15 &
     -2.34850E-15 -2.18700E-15 -2.01760E-15 -1.84120E-15 -1.65910E-15 -1.47250E-15 -1.28240E-15 -1.09000E-15 &
     -8.96280E-16 -7.01250E-16 -5.04780E-16 -3.06710E-16 -1.06890E-16 9.48220E-17 2.98590E-16 5.04530E-16 &
     7.12680E-16 9.23060E-16 1.13570E-15 1.35050E-15 1.56760E-15 1.78700E-15 2.00860E-15 2.23240E-15 &
     2.45820E-15 2.68600E-15 2.91570E-15 3.14720E-15 3.38040E-15 3.61520E-15 3.85210E-15 4.09160E-15 &
     4.33410E-15 4.58020E-15 4.83030E-15 5.08490E-15 5.34440E-15 5.60910E-15 5.87910E-15 6.15460E-15 &
     6.43570E-15 6.72250E-15 7.01530E-15 7.31400E-15 7.61790E-15 7.92630E-15 8.23830E-15 8.55310E-15 &
     8.86990E-15 9.18790E-15 9.50640E-15 9.82540E-15 1.01450E-14 1.04660E-14 1.07870E-14 1.11090E-14 &
     1.14330E-14 1.17580E-14 1.20840E-14 1.24110E-14 1.27400E-14 1.30700E-14 1.34020E-14 1.37350E-14 &
     1.40690E-14 1.44040E-14 1.47390E-14 1.50730E-14 1.54050E-14 1.57340E-14 1.60590E-14 1.63800E-14 &
     1.66970E-14 1.70110E-14 1.73230E-14 1.76320E-14 1.79400E-14 1.82460E-14 1.85530E-14 1.88580E-14 &
     1.91640E-14 1.94680E-14 1.97720E-14 2.00750E-14 2.03770E-14 2.06780E-14 2.09770E-14 2.12760E-14 &
     2.15730E-14 2.18690E-14 2.21640E-14 2.24570E-14 2.27490E-14 2.30400E-14 2.33300E-14 2.36190E-14 &
     2.39060E-14 2.41920E-14 2.44780E-14 2.47620E-14 2.50440E-14 2.53260E-14 2.56070E-14 2.58860E-14 &
     2.61640E-14 2.64400E-14 2.67160E-14 2.69900E-14 2.72630E-14 2.75350E-14 2.78070E-14 2.80790E-14 &
     2.83520E-14 2.86240E-14 2.88970E-14 2.91690E-14 2.94410E-14 2.97110E-14 2.99810E-14 3.02490E-14 &
     3.05150E-14 3.07790E-14 3.10410E-14 3.13000E-14 3.15560E-14 3.18090E-14 3.20590E-14 3.23050E-14 &
     3.25490E-14 3.27890E-14 3.30260E-14 3.32610E-14 3.34940E-14 3.37240E-14 3.39530E-14 3.41780E-14 &
     3.43990E-14 3.46140E-14 3.48220E-14 3.50210E-14 3.52100E-14 3.53890E-14 3.55570E-14 3.57170E-14 &
     3.58710E-14 3.60190E-14 3.61630E-14 3.63050E-14 3.64460E-14 3.65860E-14 3.67240E-14 3.68610E-14 &
     3.69940E-14 3.71250E-14 3.72530E-14 3.73770E-14 3.74980E-14 3.76150E-14 3.77290E-14 3.78400E-14 &
     3.79480E-14 3.80540E-14 3.81570E-14 3.82570E-14 3.83540E-14 3.84490E-14 3.85410E-14 3.86310E-14 &
     3.87170E-14 3.88010E-14 3.88820E-14 3.89600E-14 3.90350E-14 3.91070E-14 3.91760E-14 3.92420E-14 &
     3.93050E-14 3.93650E-14 3.94230E-14 3.94780E-14 3.95300E-14 3.95810E-14 3.96290E-14 3.96760E-14 &
     3.97200E-14 3.97610E-14 3.97980E-14 3.98300E-14 3.98570E-14 3.98790E-14 3.98930E-14 3.99020E-14 &
     3.99050E-14 3.99020E-14 3.98950E-14 3.98820E-14 3.98650E-14 3.98440E-14 3.98200E-14 3.97910E-14 &
     3.97590E-14 3.97220E-14 3.96830E-14 3.96400E-14 3.95930E-14 3.95420E-14 3.94870E-14 3.94270E-14 &
     3.93610E-14 3.92900E-14 3.92110E-14 3.91260E-14 3.90350E-14 3.89380E-14 3.88370E-14 3.87330E-14 &
     3.86260E-14 3.85170E-14 3.84070E-14 3.82970E-14 3.81840E-14 3.80710E-14 3.79550E-14 3.78380E-14 &
     3.77180E-14 3.75950E-14 3.74700E-14 3.73430E-14 3.72130E-14 3.70820E-14 3.69490E-14 3.68140E-14 &
     3.66780E-14 3.65400E-14 3.64010E-14 3.62600E-14 3.61170E-14 3.59720E-14 3.58260E-14 3.56770E-14 &
     3.55260E-14 3.53730E-14 3.52180E-14 3.50610E-14 3.49030E-14 3.47440E-14 3.45820E-14 3.44200E-14 &
     3.42550E-14 3.40890E-14 3.39210E-14 3.37510E-14 3.35790E-14 3.34050E-14 3.32290E-14 3.30510E-14 &
     3.28720E-14 3.26910E-14 3.25100E-14 3.23290E-14 3.21470E-14 3.19640E-14 3.17790E-14 3.15910E-14 &
     3.13990E-14 3.12010E-14 3.09960E-14 3.07850E-14 3.05670E-14 3.03440E-14 3.01170E-14 2.98860E-14 &
     2.96530E-14 2.94190E-14 2.91850E-14 2.89500E-14 2.87160E-14 2.84820E-14 2.82490E-14 2.80170E-14 &
     2.77850E-14 2.75540E-14 2.73240E-14 2.70950E-14 2.68660E-14 2.66370E-14 2.64080E-14 2.61780E-14 &
     2.59480E-14 2.57170E-14 2.54860E-14 2.52540E-14 2.50230E-14 2.47920E-14 2.45620E-14 2.43320E-14 &
     2.41040E-14 2.38760E-14 2.36500E-14 2.34240E-14 2.32000E-14 2.29770E-14 2.27560E-14 2.25360E-14 &
     2.23170E-14 2.20980E-14 2.18790E-14 2.16600E-14 2.14400E-14 2.12200E-14 2.10000E-14 2.07780E-14 &
     2.05570E-14 2.03360E-14 2.01160E-14 1.98960E-14 1.96770E-14 1.94590E-14 1.92420E-14 1.90260E-14 &
     1.88110E-14 1.85980E-14 1.83850E-14 1.81740E-14 1.79640E-14 1.77550E-14 1.75460E-14 1.73390E-14 &
     1.71320E-14 1.69260E-14 1.67200E-14 1.65150E-14 1.63100E-14 1.61060E-14 1.59020E-14 1.56990E-14 &
     1.54960E-14 1.52940E-14 1.50920E-14 1.48910E-14 1.46900E-14 1.44890E-14 1.42890E-14 1.40890E-14 &
     1.38890E-14 1.36910E-14 1.34950E-14 1.33020E-14 1.31130E-14 1.29290E-14 1.27510E-14 1.25790E-14 &
     1.24140E-14 1.22540E-14 1.20980E-14 1.19450E-14 1.17940E-14 1.16440E-14 1.14950E-14 1.13460E-14 &
     1.11980E-14 1.10500E-14 1.09040E-14 1.07600E-14 1.06180E-14 1.04770E-14 1.03390E-14 1.02030E-14 &
     1.00700E-14 9.94000E-15 9.81250E-15 9.68800E-15 9.56660E-15 9.44810E-15 9.33220E-15 9.21860E-15 &
     9.10690E-15 8.99700E-15 8.88850E-15 8.78120E-15 8.67500E-15 8.57020E-15 8.46660E-15 8.36450E-15 &
     8.26390E-15 8.16480E-15 8.06740E-15 7.97190E-15 7.87840E-15 7.78710E-15 7.69840E-15 7.61230E-15 &
     7.52910E-15 7.44900E-15 7.37180E-15 7.29720E-15 7.22520E-15 7.15540E-15 7.08770E-15 7.02180E-15 &
     6.95770E-15 6.89510E-15 6.83420E-15 6.77480E-15 6.71700E-15 6.66070E-15 6.60580E-15 6.55250E-15 &
     6.50080E-15 6.45120E-15 6.40370E-15 6.35880E-15 6.31660E-15 6.27740E-15 6.24150E-15 6.20930E-15 &
     6.18130E-15 6.15790E-15 6.13960E-15 6.12680E-15 6.11990E-15 6.11920E-15 6.12410E-15 6.13370E-15 &
     6.14720E-15 6.16360E-15 6.18230E-15 6.20220E-15 6.22280E-15 6.24400E-15 6.26600E-15 6.28880E-15 &
     6.31270E-15 6.33760E-15 6.36390E-15 6.39160E-15 6.42080E-15 6.45190E-15 6.48490E-15 6.52010E-15 &
     6.55770E-15 6.59780E-15 6.64050E-15 6.68580E-15 6.73350E-15 6.78330E-15 6.83510E-15 6.88860E-15 &
     6.94380E-15 7.00030E-15 7.05820E-15 7.11710E-15 7.17700E-15 7.23770E-15 7.29910E-15 7.36110E-15 &
     7.42350E-15 7.48670E-15 7.55080E-15 7.61630E-15 7.68350E-15 7.75260E-15 7.82390E-15 7.89770E-15 &
     7.97400E-15 8.05240E-15 8.13280E-15 8.21510E-15 8.29900E-15 8.38430E-15 8.47090E-15 8.55870E-15 &
     8.64750E-15 8.73730E-15 8.82810E-15 8.91970E-15 9.01210E-15 9.10520E-15 9.19910E-15 9.29380E-15 &
     9.38950E-15 9.48630E-15 9.58420E-15 9.68340E-15 9.78400E-15 9.88590E-15 9.98930E-15 1.00940E-14 &
     1.02000E-14 1.03080E-14 1.04170E-14 1.05270E-14 1.06390E-14 1.07530E-14 1.08680E-14 1.09850E-14 &
     1.11030E-14 1.12230E-14 1.13440E-14 1.14670E-14 1.15900E-14 1.17140E-14 1.18390E-14 1.19630E-14 &
     1.20870E-14 1.22090E-14 1.23320E-14 1.24540E-14 1.25750E-14 1.26970E-14 1.28190E-14 1.29410E-14 &
     1.30640E-14 1.31880E-14 1.33120E-14 1.34360E-14 1.35610E-14 1.36870E-14 1.38120E-14 1.39380E-14 &
     1.40640E-14 1.41900E-14 1.43150E-14 1.44400E-14 1.45640E-14 1.46870E-14 1.48080E-14 1.49280E-14 &
     1.50470E-14 1.51660E-14 1.52850E-14 1.54030E-14 1.55220E-14 1.56410E-14 1.57610E-14 1.58820E-14 &
     1.60030E-14 1.61240E-14 1.62450E-14 1.63660E-14 1.64880E-14 1.66090E-14 1.67300E-14 1.68500E-14 &
     1.69700E-14 1.70880E-14 1.72060E-14 1.73230E-14 1.74380E-14 1.75520E-14 1.76650E-14 1.77780E-14 &
     1.78890E-14 1.80000E-14 1.81100E-14 1.82190E-14 1.83260E-14 1.84320E-14 1.85370E-14 1.86390E-14 &
     1.87390E-14 1.88360E-14 1.89310E-14 1.90240E-14 1.91160E-14 1.92060E-14 1.92950E-14 1.93840E-14 &
     1.94730E-14 1.95610E-14 1.96490E-14 1.97360E-14 1.98210E-14 1.99060E-14 1.99890E-14 2.00700E-14 &
     2.01490E-14 2.02270E-14 2.03030E-14 2.03780E-14 2.04510E-14 2.05230E-14 2.05940E-14 2.06640E-14 &
     2.07330E-14 2.08010E-14 2.08670E-14 2.09340E-14 2.09990E-14 2.10640E-14 2.11280E-14 2.11900E-14 &
     2.12510E-14 2.13110E-14 2.13690E-14 2.14240E-14 2.14780E-14 2.15290E-14 2.15780E-14 2.16250E-14 &
     2.16700E-14 2.17150E-14 2.17580E-14 2.18000E-14 2.18410E-14 2.18810E-14 2.19200E-14 2.19580E-14 &
     2.19960E-14 2.20340E-14 2.20700E-14 2.21060E-14 2.21410E-14 2.21750E-14 2.22070E-14 2.22370E-14 &
     2.22660E-14 2.22920E-14 2.23170E-14 2.23390E-14 2.23580E-14 2.23760E-14 2.23920E-14 2.24050E-14 &
     2.24160E-14 2.24240E-14 2.24310E-14 2.24350E-14 2.24370E-14 2.24370E-14 2.24350E-14 2.24300E-14 &
     2.24230E-14 2.24150E-14 2.24050E-14 2.23940E-14 2.23810E-14 2.23680E-14 2.23540E-14 2.23390E-14 &
     2.23230E-14 2.23050E-14 2.22850E-14 2.22620E-14 2.22370E-14 2.22080E-14 2.21760E-14 2.21420E-14 &
     2.21050E-14 2.20660E-14 2.20250E-14 2.19840E-14 2.19410E-14 2.18970E-14 2.18520E-14 2.18060E-14 &
     2.17590E-14 2.17110E-14 2.16630E-14 2.16140E-14 2.15640E-14 2.15130E-14 2.14600E-14 2.14070E-14 &
     2.13520E-14 2.12960E-14 2.12380E-14 2.11780E-14 2.11170E-14 2.10550E-14 2.09910E-14 2.09260E-14 &
     2.08590E-14 2.07910E-14 2.07210E-14 2.06510E-14 2.05800E-14 2.05080E-14 2.04360E-14 2.03640E-14 &
     2.02920E-14 2.02190E-14 2.01460E-14 2.00730E-14 1.99990E-14 1.99240E-14 1.98470E-14 1.97690E-14 &
     1.96900E-14 1.96080E-14 1.95260E-14 1.94410E-14 1.93550E-14 1.92670E-14 1.91770E-14 1.90860E-14 &
     1.89930E-14 1.89000E-14 1.88060E-14 1.87110E-14 1.86170E-14 1.85230E-14 1.84300E-14 1.83360E-14 &
     1.82420E-14 1.81480E-14 1.80540E-14 1.79590E-14 1.78630E-14 1.77660E-14 1.76680E-14 1.75710E-14 &
     1.74720E-14 1.73740E-14 1.72750E-14 1.71760E-14 1.70760E-14 1.69760E-14 1.68750E-14 1.67720E-14 &
     1.66680E-14 1.65610E-14 1.64530E-14 1.63420E-14 1.62300E-14 1.61170E-14 1.60030E-14 1.58890E-14 &
     1.57760E-14 1.56630E-14 1.55510E-14 1.54390E-14 1.53280E-14 1.52170E-14 1.51060E-14 1.49950E-14 &
     1.48830E-14 1.47720E-14 1.46600E-14 1.45470E-14 1.44350E-14 1.43220E-14 1.42080E-14 1.40940E-14 &
     1.39800E-14 1.38660E-14 1.37520E-14 1.36380E-14 1.35250E-14 1.34120E-14 1.33010E-14 1.31900E-14 &
     1.30810E-14 1.29710E-14 1.28620E-14 1.27540E-14 1.26450E-14 1.25360E-14 1.24270E-14 1.23180E-14 &
     1.22090E-14 1.21000E-14 1.19910E-14 1.18820E-14 1.17730E-14 1.16640E-14 1.15560E-14 1.14480E-14 &
     1.13420E-14 1.12380E-14 1.11350E-14 1.10350E-14 1.09360E-14 1.08400E-14 1.07440E-14 1.06500E-14 &
     1.05570E-14 1.04650E-14 1.03720E-14 1.02800E-14 1.01890E-14 1.00970E-14 1.00060E-14 9.91430E-15 &
     9.82300E-15 9.73160E-15 9.64040E-15 9.54950E-15 9.45890E-15 9.36900E-15 9.27970E-15 9.19130E-15 &
     9.10380E-15 9.01730E-15 8.93160E-15 8.84680E-15 8.76280E-15 8.67950E-15 8.59680E-15 8.51470E-15 &
     8.43340E-15 8.35270E-15 8.27280E-15 8.19370E-15 8.11550E-15 8.03830E-15 7.96200E-15 7.88670E-15 &
     7.81220E-15 7.73860E-15 7.66560E-15 7.59330E-15 7.52160E-15 7.45040E-15 7.37980E-15 7.31000E-15 &
     7.24110E-15 7.17320E-15 7.10640E-15 7.04090E-15 6.97670E-15 6.91390E-15 6.85240E-15 6.79210E-15 &
     6.73310E-15 6.67530E-15 6.61860E-15 6.56300E-15 6.50850E-15 6.45520E-15 6.40290E-15 6.35170E-15 &
     6.30150E-15 6.25240E-15 6.20440E-15 6.15760E-15 6.11220E-15 6.06830E-15 6.02620E-15 5.98600E-15 &
     5.94780E-15 5.91190E-15 5.87800E-15 5.84610E-15 5.81590E-15 5.78740E-15 5.76030E-15 5.73450E-15 &
     5.70990E-15 5.68630E-15 5.66350E-15 5.64140E-15 5.61970E-15 5.59840E-15 5.57720E-15 5.55600E-15 &
     5.53490E-15 5.51420E-15 5.49390E-15 5.47420E-15 5.45540E-15 5.43750E-15 5.42070E-15 5.40510E-15 &
     5.39060E-15 5.37730E-15 5.36530E-15 5.35450E-15 5.34500E-15 5.33680E-15 5.32990E-15 5.32440E-15 &
     5.32020E-15 5.31740E-15 5.31590E-15 5.31600E-15 5.31740E-15 5.32010E-15 5.32400E-15 5.32890E-15 &
     5.33470E-15 5.34130E-15 5.34850E-15 5.35620E-15 5.36460E-15 5.37360E-15 5.38340E-15 5.39400E-15 &
     5.40570E-15 5.41840E-15 5.43230E-15 5.44730E-15 5.46350E-15 5.48090E-15 5.49940E-15 5.51900E-15 &
     5.53980E-15 5.56160E-15 5.58450E-15 5.60830E-15 5.63300E-15 5.65860E-15 5.68490E-15 5.71190E-15 &
     5.73950E-15 5.76770E-15 5.79650E-15 5.82610E-15 5.85630E-15 5.88720E-15 5.91880E-15 5.95110E-15 &
     5.98420E-15 6.01820E-15 6.05320E-15 6.08910E-15 6.12610E-15 6.16420E-15 6.20350E-15 6.24370E-15 &
     6.28460E-15 6.32600E-15 6.36760E-15 6.40920E-15 6.45040E-15 6.49120E-15 6.53160E-15 6.57170E-15 &
     6.61160E-15 6.65160E-15 6.69160E-15 6.73180E-15 6.77230E-15 6.81330E-15 6.85480E-15 6.89690E-15 &
     6.93970E-15 6.98330E-15 7.02770E-15 7.07300E-15 7.11920E-15 7.16610E-15 7.21370E-15 7.26190E-15 &
     7.31050E-15 7.35960E-15 7.40890E-15 7.45860E-15 7.50860E-15 7.55900E-15 7.60980E-15 7.66100E-15 &
     7.71260E-15 7.76470E-15 7.81730E-15 7.87020E-15 7.92340E-15 7.97690E-15 8.03070E-15 8.08460E-15 &
     8.13870E-15 8.19300E-15 8.24750E-15 8.30220E-15 8.35720E-15 8.41250E-15 8.46820E-15 8.52420E-15 &
     8.58040E-15 8.63670E-15 8.69290E-15 8.74890E-15 8.80460E-15 8.85980E-15 8.91440E-15 8.96820E-15 &
     9.02140E-15 9.07390E-15 9.12570E-15 9.17660E-15 9.22680E-15 9.27620E-15 9.32500E-15 9.37330E-15 &
     9.42150E-15 9.46960E-15 9.51800E-15 9.56680E-15 9.61610E-15 9.66590E-15 9.71590E-15 9.76590E-15 &
     9.81570E-15 9.86520E-15 9.91410E-15 9.96230E-15 1.00100E-14 1.00560E-14 1.01020E-14 1.01470E-14 &
     1.01910E-14 1.02350E-14 1.02770E-14 1.03190E-14 1.03610E-14 1.04020E-14 1.04430E-14 1.04850E-14 &
     1.05270E-14 1.05690E-14 1.06120E-14 1.06550E-14 1.06980E-14 1.07410E-14 1.07840E-14 1.08260E-14 &
     1.08670E-14 1.09080E-14 1.09480E-14 1.09870E-14 1.10260E-14 1.10640E-14 1.11010E-14 1.11370E-14 &
     1.11730E-14 1.12080E-14 1.12430E-14 1.12770E-14 1.13100E-14 1.13430E-14 1.13760E-14 1.14080E-14 &
     1.14400E-14 1.14710E-14 1.15030E-14 1.15330E-14 1.15640E-14 1.15940E-14 1.16240E-14 1.16530E-14 &
     1.16810E-14 1.17080E-14 1.17340E-14 1.17590E-14 1.17820E-14 1.18040E-14 1.18250E-14 1.18440E-14 &
     1.18620E-14 1.18780E-14 1.18940E-14 1.19080E-14 1.19210E-14 1.19330E-14 1.19440E-14 1.19550E-14 &
     1.19660E-14 1.19760E-14 1.19860E-14 1.19970E-14 1.20070E-14 1.20170E-14 1.20270E-14 1.20360E-14 &
     1.20450E-14 1.20530E-14 1.20610E-14 1.20670E-14 1.20730E-14 1.20780E-14 1.20820E-14 1.20850E-14 &
     1.20860E-14 1.20860E-14 1.20850E-14 1.20840E-14 1.20820E-14 1.20800E-14 1.20780E-14 1.20760E-14 &
     1.20740E-14 1.20720E-14 1.20700E-14 1.20670E-14 1.20650E-14 1.20620E-14 1.20580E-14 1.20530E-14 &
     1.20480E-14 1.20420E-14 1.20350E-14 1.20280E-14 1.20190E-14 1.20100E-14 1.20000E-14 1.19880E-14 &
     1.19760E-14 1.19630E-14 1.19480E-14 1.19330E-14 1.19160E-14 1.18990E-14 1.18810E-14 1.18620E-14 &
     1.18420E-14 1.18210E-14 1.18010E-14 1.17790E-14 1.17580E-14 1.17350E-14 1.17130E-14 1.16890E-14 &
     1.16640E-14 1.16380E-14 1.16110E-14 1.15830E-14 1.15540E-14 1.15240E-14 1.14930E-14 1.14620E-14 &
     1.14290E-14 1.13960E-14 1.13630E-14 1.13290E-14 1.12950E-14 1.12600E-14 1.12250E-14 1.11900E-14 &
     1.11550E-14 1.11200E-14 1.10850E-14 1.10500E-14 1.10140E-14 1.09790E-14 1.09440E-14 1.09080E-14 &
     1.08730E-14 1.08370E-14 1.08000E-14 1.07630E-14 1.07250E-14 1.06870E-14 1.06470E-14 1.06070E-14 &
     1.05660E-14 1.05240E-14 1.04820E-14 1.04390E-14 1.03970E-14 1.03540E-14 1.03120E-14 1.02690E-14 &
     1.02270E-14 1.01850E-14 1.01440E-14 1.01020E-14 1.00610E-14 1.00200E-14 9.97820E-15 9.93650E-15 &
     9.89420E-15 9.85120E-15 9.80730E-15 9.76230E-15 9.71630E-15 9.66960E-15 9.62230E-15 9.57480E-15 &
     9.52710E-15 9.47950E-15 9.43220E-15 9.38500E-15 9.33790E-15 9.29050E-15 9.24290E-15 9.19480E-15 &
     9.14610E-15 9.09660E-15 9.04650E-15 8.99580E-15 8.94480E-15 8.89360E-15 8.84230E-15 8.79100E-15 &
     8.73980E-15 8.68880E-15 8.63800E-15 8.58730E-15 8.53670E-15 8.48620E-15 8.43580E-15 8.38550E-15 &
     8.33520E-15 8.28510E-15 8.23520E-15 8.18540E-15 8.13590E-15 8.08650E-15 8.03750E-15 7.98870E-15 &
     7.94030E-15 7.89220E-15 7.84450E-15 7.79720E-15 7.75040E-15 7.70400E-15 7.65800E-15 7.61220E-15 &
     7.56660E-15 7.52110E-15 7.47550E-15 7.42970E-15 7.38360E-15 7.33730E-15 7.29080E-15 7.24420E-15 &
     7.19760E-15 7.15090E-15 7.10430E-15 7.05780E-15 7.01150E-15 6.96550E-15 6.92000E-15 6.87510E-15 &
     6.83100E-15 6.78760E-15 6.74520E-15 6.70350E-15 6.66240E-15 6.62160E-15 6.58100E-15 6.54030E-15 &
     6.49930E-15 6.45790E-15 6.41610E-15 6.37400E-15 6.33180E-15 6.28940E-15 6.24700E-15 6.20460E-15 &
     6.16250E-15 6.12050E-15 6.07890E-15 6.03770E-15 5.99690E-15 5.95670E-15 5.91710E-15 5.87810E-15 &
     5.83980E-15 5.80220E-15 5.76530E-15 5.72920E-15 5.69390E-15 5.65940E-15 5.62570E-15 5.59270E-15 &
     5.56050E-15 5.52890E-15 5.49780E-15 5.46730E-15 5.43720E-15 5.40750E-15 5.37810E-15 5.34920E-15 &
     5.32070E-15 5.29280E-15 5.26530E-15 5.23840E-15 5.21210E-15 5.18630E-15 5.16110E-15 5.13650E-15 &
     5.11250E-15 5.08910E-15 5.06630E-15 5.04400E-15 5.02230E-15 5.00110E-15 4.98030E-15 4.95990E-15 &
     4.93990E-15 4.92020E-15 4.90070E-15 4.88150E-15 4.86250E-15 4.84370E-15 4.82500E-15 4.80640E-15 &
     4.78800E-15 4.76960E-15 4.75150E-15 4.73370E-15 4.71650E-15 4.69990E-15 4.68400E-15 4.66910E-15 &
     4.65530E-15 4.64240E-15 4.63020E-15 4.61870E-15 4.60770E-15 4.59710E-15 4.58670E-15 4.57640E-15 &
     4.56620E-15 4.55600E-15 4.54580E-15 4.53560E-15 4.52540E-15 4.51520E-15 4.50490E-15 4.49470E-15 &
     4.48480E-15 4.47520E-15 4.46620E-15 4.45780E-15 4.45030E-15 4.44370E-15 4.43800E-15 4.43330E-15 &
     4.42960E-15 4.42680E-15 4.42490E-15 4.42390E-15 4.42390E-15 4.42460E-15 4.42620E-15 4.42840E-15 &
     4.43110E-15 4.43440E-15 4.43810E-15 4.44220E-15 4.44660E-15 4.45140E-15 4.45660E-15 4.46230E-15 &
     4.46850E-15 4.47530E-15 4.48260E-15 4.49040E-15 4.49880E-15 4.50760E-15 4.51690E-15 4.52660E-15 &
     4.53660E-15 4.54700E-15 4.55770E-15 4.56860E-15 4.57980E-15 4.59130E-15 4.60290E-15 4.61470E-15 &
     4.62670E-15 4.63870E-15 4.65090E-15 4.66310E-15 4.67540E-15 4.68770E-15 4.70000E-15 4.71230E-15 &
     4.72460E-15 4.73710E-15 4.74990E-15 4.76300E-15 4.77660E-15 4.79070E-15 4.80540E-15 4.82070E-15 &
     4.83650E-15 4.85280E-15 4.86960E-15 4.88660E-15 4.90400E-15 4.92160E-15 4.93940E-15 4.95720E-15 &
     4.97500E-15 4.99250E-15 5.00990E-15 5.02680E-15 5.04340E-15 5.05960E-15 5.07580E-15 5.09200E-15 &
     5.10840E-15 5.12530E-15 5.14270E-15 5.16090E-15 5.17980E-15 5.19930E-15 5.21940E-15 5.24010E-15 &
     5.26130E-15 5.28280E-15 5.30480E-15 5.32700E-15 5.34930E-15 5.37160E-15 5.39380E-15 5.41580E-15 &
     5.43740E-15 5.45850E-15 5.47940E-15 5.50000E-15 5.52060E-15 5.54130E-15 5.56230E-15 5.58360E-15 &
     5.60550E-15 5.62780E-15 5.65060E-15 5.67360E-15 5.69690E-15 5.72040E-15 5.74400E-15 5.76750E-15 &
     5.79110E-15 5.81460E-15 5.83790E-15 5.86110E-15 5.88410E-15 5.90680E-15 5.92920E-15 5.95130E-15 &
     5.97320E-15 5.99480E-15 6.01610E-15 6.03720E-15 6.05800E-15 6.07870E-15 6.09910E-15 6.11950E-15 &
     6.13970E-15 6.15990E-15 6.18010E-15 6.20030E-15 6.22060E-15 6.24100E-15 6.26150E-15 6.28220E-15 &
     6.30290E-15 6.32380E-15 6.34480E-15 6.36600E-15 6.38720E-15 6.40840E-15 6.42930E-15 6.44990E-15 &
     6.47000E-15 6.48940E-15 6.50820E-15 6.52640E-15 6.54410E-15 6.56140E-15 6.57860E-15 6.59560E-15 &
     6.61270E-15 6.62990E-15 6.64730E-15 6.66480E-15 6.68240E-15 6.70010E-15 6.71780E-15 6.73550E-15 &
     6.75320E-15 6.77070E-15 6.78790E-15 6.80480E-15 6.82120E-15 6.83690E-15 6.85190E-15 6.86600E-15 &
     6.87940E-15 6.89210E-15 6.90420E-15 6.91600E-15 6.92740E-15 6.93850E-15 6.94960E-15 6.96060E-15 &
     6.97150E-15 6.98240E-15 6.99330E-15 7.00420E-15 7.01520E-15 7.02630E-15 7.03730E-15 7.04820E-15 &
     7.05890E-15 7.06920E-15 7.07910E-15 7.08850E-15 7.09720E-15 7.10540E-15 7.11310E-15 7.12020E-15 &
     7.12690E-15 7.13320E-15 7.13920E-15 7.14490E-15 7.15040E-15 7.15550E-15 7.16050E-15 7.16520E-15 &
     7.16960E-15 7.17390E-15 7.17790E-15 7.18180E-15 7.18550E-15 7.18910E-15 7.19250E-15 7.19590E-15 &
     7.19930E-15 7.20260E-15 7.20570E-15 7.20870E-15 7.21140E-15 7.21370E-15 7.21560E-15 7.21690E-15 &
     7.21760E-15 7.21770E-15 7.21730E-15 7.21640E-15 7.21510E-15 7.21340E-15 7.21140E-15 7.20910E-15 &
     7.20660E-15 7.20380E-15 7.20080E-15 7.19770E-15 7.19440E-15 7.19110E-15 7.18760E-15 7.18400E-15 &
     7.18020E-15 7.17600E-15 7.17150E-15 7.16640E-15 7.16070E-15 7.15450E-15 7.14760E-15 7.14010E-15 &
     7.13210E-15 7.12370E-15 7.11480E-15 7.10560E-15 7.09600E-15 7.08620E-15 7.07610E-15 7.06590E-15 &
     7.05570E-15 7.04540E-15 7.03510E-15 7.02490E-15 7.01460E-15 7.00430E-15 6.99370E-15 6.98290E-15 &
     6.97160E-15 6.95990E-15 6.94760E-15 6.93480E-15 6.92150E-15 6.90770E-15 6.89340E-15 6.87880E-15 &
     6.86380E-15 6.84850E-15 6.83300E-15 6.81720E-15 6.80130E-15 6.78520E-15 6.76910E-15 6.75300E-15 &
     6.73680E-15 6.72070E-15 6.70460E-15 6.68850E-15 6.67230E-15 6.65610E-15 6.63980E-15 6.62340E-15 &
     6.60690E-15 6.59020E-15 6.57330E-15 6.55620E-15 6.53890E-15 6.52120E-15 6.50330E-15 6.48510E-15 &
     6.46650E-15 6.44760E-15 6.42830E-15 6.40870E-15 6.38870E-15 6.36840E-15 6.34770E-15 6.32690E-15 &
     6.30610E-15 6.28520E-15 6.26450E-15 6.24400E-15 6.22390E-15 6.20390E-15 6.18420E-15 6.16460E-15 &
     6.14500E-15 6.12530E-15 6.10560E-15 6.08560E-15 6.06540E-15 6.04500E-15 6.02420E-15 6.00300E-15 &
     5.98130E-15 5.95920E-15 5.93660E-15 5.91370E-15 5.89050E-15 5.86730E-15 5.84430E-15 5.82150E-15 &
     5.79930E-15 5.77760E-15 5.75650E-15 5.73570E-15 5.71500E-15 5.69450E-15 5.67380E-15 5.65290E-15 &
     5.63160E-15 5.60990E-15 5.58790E-15 5.56550E-15 5.54280E-15 5.51970E-15 5.49640E-15 5.47280E-15 &
     5.44900E-15 5.42520E-15 5.40130E-15 5.37750E-15 5.35400E-15 5.33080E-15 5.30790E-15 5.28540E-15 &
     5.26310E-15 5.24110E-15 5.21940E-15 5.19770E-15 5.17620E-15 5.15470E-15 5.13320E-15 5.11180E-15 &
     5.09040E-15 5.06900E-15 5.04750E-15 5.02600E-15 5.00450E-15 4.98300E-15 4.96130E-15 4.93960E-15 &
     4.91770E-15 4.89580E-15 4.87360E-15 4.85140E-15 4.82910E-15 4.80700E-15 4.78500E-15 4.76340E-15 &
     4.74220E-15 4.72160E-15 4.70170E-15 4.68250E-15 4.66370E-15 4.64550E-15 4.62760E-15 4.61010E-15 &
     4.59270E-15 4.57550E-15 4.55840E-15 4.54130E-15 4.52410E-15 4.50670E-15 4.48900E-15 4.47100E-15 &
     4.45260E-15 4.43390E-15 4.41510E-15 4.39620E-15 4.37750E-15 4.35910E-15 4.34110E-15 4.32370E-15 &
     4.30680E-15 4.29040E-15 4.27440E-15 4.25890E-15 4.24380E-15 4.22910E-15 4.21470E-15 4.20050E-15 &
     4.18650E-15 4.17250E-15 4.15840E-15 4.14430E-15 4.12980E-15 4.11510E-15 4.10010E-15 4.08500E-15 &
     4.07000E-15 4.05500E-15 4.04040E-15 4.02610E-15 4.01240E-15 3.99910E-15 3.98620E-15 3.97390E-15 &
     3.96190E-15 3.95030E-15 3.93910E-15 3.92830E-15 3.91770E-15 3.90740E-15 3.89720E-15 3.88720E-15 &
     3.87720E-15 3.86730E-15 3.85740E-15 3.84740E-15 3.83740E-15 3.82740E-15 3.81750E-15 3.80760E-15 &
     3.79790E-15 3.78820E-15 3.77870E-15 3.76940E-15 3.76050E-15 3.75200E-15 3.74400E-15 3.73650E-15 &
     3.72970E-15 3.72360E-15 3.71790E-15 3.71280E-15 3.70810E-15 3.70380E-15 3.69980E-15 3.69610E-15 &
     3.69260E-15 3.68930E-15 3.68600E-15 3.68270E-15 3.67920E-15 3.67560E-15 3.67170E-15 3.66770E-15 &
     3.66350E-15 3.65930E-15 3.65510E-15 3.65110E-15 3.64720E-15 3.64360E-15 3.64030E-15 3.63740E-15 &
     3.63490E-15 3.63290E-15 3.63140E-15 3.63050E-15 3.63010E-15 3.63020E-15 3.63060E-15 3.63120E-15 &
     3.63180E-15 3.63240E-15 3.63260E-15 3.63250E-15 3.63220E-15 3.63170E-15 3.63120E-15 3.63070E-15 &
     3.63040E-15 3.63050E-15 3.63100E-15 3.63190E-15 3.63330E-15 3.63510E-15 3.63740E-15 3.64030E-15 &
     3.64350E-15 3.64730E-15 3.65150E-15 3.65600E-15 3.66070E-15 3.66540E-15 3.67020E-15 3.67480E-15 &
     3.67920E-15 3.68350E-15 3.68770E-15 3.69170E-15 3.69580E-15 3.69980E-15 3.70390E-15 3.70810E-15 &
     3.71240E-15 3.71690E-15 3.72170E-15 3.72670E-15 3.73190E-15 3.73750E-15 3.74350E-15 3.74980E-15 &
     3.75640E-15 3.76330E-15 3.77040E-15 3.77780E-15 3.78530E-15 3.79290E-15 3.80070E-15 3.80850E-15 &
     3.81630E-15 3.82410E-15 3.83170E-15 3.83930E-15 3.84670E-15 3.85400E-15 3.86110E-15 3.86810E-15 &
     3.87490E-15 3.88170E-15 3.88830E-15 3.89490E-15 3.90160E-15 3.90830E-15 3.91520E-15 3.92240E-15 &
     3.93000E-15 3.93800E-15 3.94660E-15 3.95560E-15 3.96500E-15 3.97440E-15 3.98400E-15 3.99340E-15 &
     4.00260E-15 4.01150E-15 4.02000E-15 4.02830E-15 4.03630E-15 4.04420E-15 4.05180E-15 4.05930E-15 &
     4.06670E-15 4.07420E-15 4.08170E-15 4.08940E-15 4.09740E-15 4.10590E-15 4.11480E-15 4.12440E-15 &
     4.13440E-15 4.14470E-15 4.15530E-15 4.16590E-15 4.17640E-15 4.18680E-15 4.19680E-15 4.20650E-15 &
     4.21590E-15 4.22510E-15 4.23390E-15 4.24260E-15 4.25110E-15 4.25940E-15 4.26760E-15 4.27570E-15 &
     4.28370E-15 4.29170E-15 4.29980E-15 4.30790E-15 4.31620E-15 4.32450E-15 4.33290E-15 4.34130E-15 &
     4.34980E-15 4.35820E-15 4.36660E-15 4.37500E-15 4.38330E-15 4.39150E-15 4.39960E-15 4.40760E-15 &
     4.41540E-15 4.42310E-15 4.43070E-15 4.43800E-15 4.44520E-15 4.45220E-15 4.45890E-15 4.46540E-15 &
     4.47170E-15 4.47780E-15 4.48360E-15 4.48940E-15 4.49520E-15 4.50100E-15 4.50690E-15 4.51310E-15 &
     4.51960E-15 4.52630E-15 4.53310E-15 4.54010E-15 4.54700E-15 4.55400E-15 4.56080E-15 4.56750E-15 &
     4.57390E-15 4.58010E-15 4.58590E-15 4.59140E-15 4.59650E-15 4.60110E-15 4.60520E-15 4.60900E-15 &
     4.61250E-15 4.61590E-15 4.61940E-15 4.62300E-15 4.62690E-15 4.63120E-15 4.63590E-15 4.64070E-15 &
     4.64570E-15 4.65070E-15 4.65570E-15 4.66040E-15 4.66480E-15 4.66900E-15 4.67290E-15 4.67660E-15 &
     4.68000E-15 4.68330E-15 4.68630E-15 4.68910E-15 4.69180E-15 4.69430E-15 4.69670E-15 4.69880E-15 &
     4.70080E-15 4.70270E-15 4.70430E-15 4.70580E-15 4.70710E-15 4.70830E-15 4.70930E-15 4.71010E-15 &
     4.71070E-15 4.71120E-15 4.71160E-15 4.71180E-15 4.71180E-15 4.71160E-15 4.71130E-15 4.71090E-15 &
     4.71030E-15 4.70950E-15 4.70860E-15 4.70750E-15 4.70630E-15 4.70490E-15 4.70340E-15 4.70170E-15 &
     4.69990E-15 4.69790E-15 4.69580E-15 4.69360E-15 4.69120E-15 4.68860E-15 4.68600E-15 4.68310E-15 &
     4.68020E-15 4.67710E-15 4.67380E-15 4.67040E-15 4.66690E-15 4.66330E-15 4.65950E-15 4.65560E-15 &
     4.65150E-15 4.64730E-15 4.64300E-15 4.63860E-15 4.63400E-15 4.62930E-15 4.62450E-15 4.61950E-15 &
     4.61440E-15 4.60920E-15 4.60390E-15 4.59850E-15 4.59290E-15 4.58730E-15 4.58150E-15 4.57560E-15 &
     4.56960E-15 4.56350E-15 4.55730E-15 4.55090E-15 4.54450E-15 4.53800E-15 4.53130E-15 4.52460E-15 &
     4.51770E-15 4.51080E-15 4.50380E-15 4.49660E-15 4.48940E-15 4.48210E-15 4.47470E-15 4.46720E-15 &
     4.45960E-15 4.45190E-15 4.44420E-15 4.43630E-15 4.42840E-15 4.42040E-15 4.41230E-15 4.40420E-15 &
     4.39590E-15 4.38760E-15 4.37920E-15 4.37080E-15 4.36220E-15 4.35370E-15 4.34500E-15 4.33630E-15 &
     4.32750E-15 4.31860E-15 4.30970E-15 4.30080E-15 4.29170E-15 4.28260E-15 4.27350E-15 4.26430E-15 &
     4.25500E-15 4.24570E-15 4.23640E-15 4.22700E-15 4.21760E-15 4.20810E-15 4.19850E-15 4.18900E-15 &
     4.17940E-15 4.16970E-15 4.16000E-15 4.15030E-15 4.14050E-15 4.13070E-15 4.12090E-15 4.11100E-15 &
     4.10120E-15 4.09120E-15 4.08130E-15 4.07130E-15 4.06140E-15 4.05140E-15 4.04130E-15 4.03130E-15 &
     4.02120E-15 4.01110E-15 4.00110E-15 3.99090E-15 3.98080E-15 3.97070E-15 3.96060E-15 3.95040E-15 &
     3.94030E-15 3.93010E-15 3.91990E-15 3.90980E-15 3.89960E-15 3.88940E-15 3.87920E-15 3.86910E-15 &
     3.85890E-15 3.84870E-15 3.83860E-15 3.82840E-15 3.81830E-15 3.80810E-15 3.79800E-15 3.78790E-15 &
     3.77780E-15 3.76770E-15 3.75760E-15 3.74750E-15 3.73740E-15 3.72740E-15 3.71740E-15 3.70740E-15 &
     3.69740E-15 3.68750E-15 3.67750E-15 3.66760E-15 3.65780E-15 3.64790E-15 3.63810E-15 3.62830E-15 &
     3.61850E-15 3.60880E-15 3.59910E-15 3.58940E-15 3.57980E-15 3.57020E-15 3.56060E-15 3.55110E-15 &
     3.54170E-15 3.53220E-15 3.52280E-15 3.51350E-15 3.50420E-15 3.49490E-15 3.48570E-15 3.47650E-15 &
     3.46740E-15 3.45840E-15 3.44940E-15 3.44040E-15 3.43150E-15 3.42270E-15 3.41390E-15 3.40510E-15 &
     3.39650E-15 3.38790E-15 3.37930E-15 3.37080E-15 3.36240E-15 3.35400E-15 3.34570E-15 3.33750E-15 &
     3.32940E-15 3.32130E-15 3.31330E-15 3.30530E-15 3.29750E-15 3.28970E-15 3.28190E-15 3.27430E-15 &
     3.26670E-15 3.25920E-15 3.25180E-15 3.24450E-15 3.23730E-15 3.23010E-15 3.22310E-15 3.21610E-15 &
     3.20920E-15 3.20240E-15 3.19570E-15 3.18900E-15 3.18250E-15 3.17610E-15 3.16970E-15 3.16350E-15 &
     3.15730E-15 3.15130E-15 3.14530E-15 3.13950E-15 3.13370E-15 3.12810E-15 3.12250E-15 3.11710E-15 &
     3.11170E-15 3.10650E-15 3.10140E-15 3.09640E-15 3.09150E-15 3.08670E-15 3.08210E-15 3.07750E-15 &
     3.07310E-15 3.06880E-15 3.06460E-15 3.06050E-15 3.05650E-15 3.05260E-15 3.04890E-15 3.04530E-15 &
     3.04170E-15 3.03830E-15 3.03500E-15 3.03180E-15 3.02870E-15 3.02570E-15 3.02290E-15 3.02010E-15 &
     3.01740E-15 3.01480E-15 3.01230E-15 3.01000E-15 3.00770E-15 3.00550E-15 3.00340E-15 3.00140E-15 &
     2.99950E-15 2.99770E-15 2.99600E-15 2.99440E-15 2.99280E-15 2.99140E-15 2.99000E-15 2.98880E-15 &
     2.98760E-15 2.98650E-15 2.98540E-15 2.98450E-15 2.98370E-15 2.98290E-15 2.98220E-15 2.98160E-15 &
     2.98100E-15 2.98060E-15 2.98020E-15 2.97990E-15 2.97960E-15 2.97950E-15 2.97940E-15 2.97930E-15 &
     2.97940E-15 2.97950E-15 2.97970E-15 2.97990E-15 2.98020E-15 2.98060E-15 2.98100E-15 2.98150E-15 &
     2.98210E-15 2.98270E-15 2.98340E-15 2.98410E-15 2.98490E-15 2.98580E-15 2.98670E-15 2.98770E-15 &
     2.98870E-15 2.98970E-15 2.99080E-15 2.99200E-15 2.99320E-15 2.99450E-15 2.99580E-15 2.99710E-15 &
     2.99850E-15 3.00000E-15 3.00150E-15 3.00300E-15 3.00450E-15 3.00610E-15 3.00780E-15 3.00950E-15 &
     3.01120E-15 3.01290E-15 3.01470E-15 3.01650E-15 3.01840E-15 3.02030E-15 3.02220E-15 3.02410E-15 &
     3.02610E-15 3.02810E-15 3.03010E-15 3.03220E-15 3.03420E-15 3.03630E-15 3.03840E-15 3.04060E-15 &
     3.04270E-15 3.04490E-15 3.04710E-15 3.04930E-15 3.05160E-15 3.05380E-15 3.05610E-15 3.05830E-15 &
     3.06060E-15 3.06290E-15 3.06520E-15 3.06750E-15 3.06980E-15 3.07210E-15 3.07450E-15 3.07680E-15 &
     3.07910E-15 3.08150E-15 3.08380E-15 3.08620E-15 3.08850E-15 3.09090E-15 3.09320E-15 3.09550E-15 &
     3.09790E-15 3.10020E-15 3.10250E-15 3.10480E-15 3.10710E-15 3.10940E-15 3.11170E-15 3.11390E-15 &
     3.11620E-15 3.11840E-15 3.12070E-15 3.12290E-15 3.12510E-15 3.12720E-15 3.12940E-15 3.13150E-15 &
     3.13360E-15 3.13570E-15 3.13780E-15 3.13980E-15 3.14190E-15 3.14390E-15 3.14580E-15 3.14780E-15 &
     3.14970E-15 3.15150E-15 3.15340E-15 3.15520E-15 3.15700E-15 3.15870E-15 3.16040E-15 3.16210E-15 &
     3.16380E-15 3.16540E-15 3.16690E-15 3.16840E-15 3.16990E-15 3.17130E-15 3.17270E-15 3.17410E-15 &
     3.17540E-15 3.17660E-15 3.17780E-15 3.17900E-15 3.18010E-15 3.18120E-15 3.18220E-15 3.18310E-15 &
     3.18400E-15 3.18490E-15 3.18560E-15 3.18640E-15 3.18700E-15 3.18770E-15 3.18820E-15 3.18870E-15 &
     3.18920E-15 3.18960E-15 3.18990E-15 3.19020E-15 3.19040E-15 3.19060E-15 3.19070E-15 3.19080E-15 &
     3.19080E-15 3.19080E-15 3.19070E-15 3.19050E-15 3.19030E-15 3.19010E-15 3.18980E-15 3.18950E-15 &
     3.18910E-15 3.18860E-15 3.18810E-15 3.18760E-15 3.18700E-15 3.18630E-15 3.18560E-15 3.18490E-15 &
     3.18410E-15 3.18330E-15 3.18240E-15 3.18150E-15 3.18050E-15 3.17950E-15 3.17840E-15 3.17730E-15 &
     3.17610E-15 3.17490E-15 3.17360E-15 3.17230E-15 3.17100E-15 3.16960E-15 3.16820E-15 3.16670E-15 &
     3.16520E-15 3.16360E-15 3.16200E-15 3.16040E-15 3.15870E-15 3.15700E-15 3.15520E-15 3.15340E-15 &
     3.15160E-15 3.14970E-15 3.14770E-15 3.14580E-15 3.14380E-15 3.14170E-15 3.13960E-15 3.13750E-15 &
     3.13530E-15 3.13310E-15 3.13090E-15 3.12860E-15 3.12630E-15 3.12400E-15 3.12160E-15 3.11920E-15 &
     3.11670E-15 3.11420E-15 3.11170E-15 3.10910E-15 3.10650E-15 3.10390E-15 3.10120E-15 3.09850E-15 &
     3.09580E-15 3.09300E-15 3.09020E-15 3.08740E-15 3.08450E-15 3.08170E-15 3.07870E-15 3.07580E-15 &
     3.07280E-15 3.06980E-15 3.06670E-15 3.06360E-15 3.06050E-15 3.05740E-15 3.05420E-15 3.05110E-15 &
     3.04780E-15 3.04460E-15 3.04130E-15 3.03800E-15 3.03470E-15 3.03130E-15 3.02790E-15 3.02450E-15 &
     3.02110E-15 3.01760E-15 3.01410E-15 3.01060E-15 3.00710E-15 3.00350E-15 3.00000E-15 2.99640E-15 &
     2.99270E-15 2.98910E-15 2.98540E-15 2.98170E-15 2.97800E-15 2.97420E-15 2.97050E-15 2.96670E-15 &
     2.96290E-15 2.95910E-15 2.95520E-15 2.95130E-15 2.94750E-15 2.94360E-15 2.93960E-15 2.93570E-15 &
     2.93170E-15 2.92780E-15 2.92380E-15 2.91970E-15 2.91570E-15 2.91170E-15 2.90760E-15 2.90350E-15 &
     2.89940E-15 2.89530E-15 2.89120E-15 2.88700E-15 2.88290E-15 2.87870E-15 2.87450E-15 2.87030E-15 &
     2.86610E-15 2.86190E-15 2.85770E-15 2.85340E-15 2.84910E-15 2.84490E-15 2.84060E-15 2.83630E-15 &
     2.83200E-15 2.82770E-15 2.82330E-15 2.81900E-15 2.81470E-15 2.81030E-15 2.80590E-15 2.80160E-15 &
     2.79720E-15 2.79280E-15 2.78840E-15 2.78400E-15 2.77960E-15 2.77510E-15 2.77070E-15 2.76630E-15 &
     2.76180E-15 2.75740E-15 2.75290E-15 2.74850E-15 2.74400E-15 2.73960E-15 2.73510E-15 2.73060E-15 &
     2.72610E-15 2.72170E-15 2.71720E-15 2.71270E-15 2.70820E-15 2.70370E-15 2.69920E-15 2.69470E-15 &
     2.69020E-15 2.68570E-15 2.68130E-15 2.67680E-15 2.67230E-15 2.66780E-15 2.66330E-15 2.65880E-15 &
     2.65430E-15 2.64980E-15 2.64530E-15 2.64080E-15 2.63640E-15 2.63190E-15 2.62740E-15 2.62290E-15 &
     2.61850E-15 2.61400E-15 2.60950E-15 2.60510E-15 2.60060E-15 2.59620E-15 2.59180E-15 2.58730E-15 &
     2.58290E-15 2.57850E-15 2.57410E-15 2.56970E-15 2.56530E-15 2.56090E-15 2.55650E-15 2.55210E-15 &
     2.54780E-15 2.54340E-15 2.53910E-15 2.53480E-15 2.53040E-15 2.52610E-15 2.52180E-15 2.51750E-15 &
     2.51330E-15 2.50900E-15 2.50470E-15 2.50050E-15 2.49630E-15 2.49200E-15 2.48780E-15 2.48360E-15 &
     2.47950E-15 2.47530E-15 2.47120E-15 2.46700E-15 2.46290E-15 2.45880E-15 2.45470E-15 2.45060E-15 &
     2.44660E-15 2.44250E-15 2.43850E-15 2.43450E-15 2.43050E-15 2.42660E-15 2.42260E-15 2.41870E-15 &
     2.41480E-15 2.41090E-15 2.40700E-15 2.40310E-15 2.39930E-15 2.39550E-15 2.39170E-15 2.38790E-15 &
     2.38420E-15 2.38040E-15 2.37670E-15 2.37300E-15 2.36940E-15 2.36570E-15 2.36210E-15 2.35850E-15 &
     2.35490E-15 2.35140E-15 2.34780E-15 2.34430E-15 2.34080E-15 2.33740E-15 2.33390E-15 2.33050E-15 &
     2.32710E-15 2.32370E-15 2.32030E-15 2.31700E-15 2.31370E-15 2.31030E-15 2.30710E-15 2.30380E-15 &
     2.30060E-15 2.29730E-15 2.29410E-15 2.29090E-15 2.28780E-15 2.28460E-15 2.28150E-15 2.27840E-15 &
     2.27530E-15 2.27230E-15 2.26920E-15 2.26620E-15 2.26320E-15 2.26020E-15 2.25730E-15 2.25430E-15 &
     2.25140E-15 2.24850E-15 2.24560E-15 2.24270E-15 2.23990E-15 2.23710E-15 2.23430E-15 2.23150E-15 &
     2.22870E-15 2.22600E-15 2.22320E-15 2.22050E-15 2.21780E-15 2.21510E-15 2.21250E-15 2.20980E-15 &
     2.20720E-15 2.20460E-15 2.20200E-15 2.19950E-15 2.19690E-15 2.19440E-15 2.19190E-15 2.18940E-15 &
     2.18690E-15 2.18450E-15 2.18200E-15 2.17960E-15 2.17720E-15 2.17480E-15 2.17240E-15 2.17010E-15 &
     2.16770E-15 2.16540E-15 2.16310E-15 2.16080E-15 2.15860E-15 2.15630E-15 2.15410E-15 2.15190E-15 &
     2.14970E-15 2.14750E-15 2.14530E-15 2.14320E-15 2.14110E-15 2.13900E-15 2.13690E-15 2.13480E-15 &
     2.13270E-15 2.13070E-15 2.12860E-15 2.12660E-15 2.12460E-15 2.12260E-15 2.12070E-15 2.11870E-15 &
     2.11680E-15 2.11490E-15 2.11300E-15 2.11110E-15 2.10920E-15 2.10730E-15 2.10550E-15 2.10370E-15 &
     2.10190E-15 2.10010E-15 2.09830E-15 2.09650E-15 2.09480E-15 2.09300E-15 2.09130E-15 2.08960E-15 &
     2.08790E-15 2.08630E-15 2.08460E-15 2.08300E-15 2.08130E-15 2.07970E-15 2.07810E-15 2.07650E-15 &
     2.07490E-15 2.07340E-15 2.07180E-15 2.07030E-15 2.06880E-15 2.06730E-15 2.06580E-15 2.06430E-15 &
     2.06290E-15 2.06140E-15 2.06000E-15 2.05860E-15 2.05720E-15 2.05580E-15 2.05440E-15 2.05300E-15 &
     2.05170E-15 2.05030E-15 2.04900E-15 2.04770E-15 2.04640E-15 2.04510E-15 2.04390E-15 2.04260E-15 &
     2.04130E-15 2.04010E-15 2.03890E-15 2.03770E-15 2.03650E-15 2.03530E-15 2.03410E-15 2.03300E-15 &
     2.03180E-15 2.03070E-15 2.02960E-15 2.02850E-15 2.02740E-15 2.02630E-15 2.02520E-15 2.02410E-15 &
     2.02310E-15 2.02200E-15 2.02100E-15 2.02000E-15 2.01900E-15 2.01800E-15 2.01700E-15 2.01610E-15 &
     2.01510E-15 2.01420E-15 2.01320E-15 2.01230E-15 2.01140E-15 2.01050E-15 2.00960E-15 2.00870E-15 &
     2.00780E-15 2.00700E-15 2.00610E-15 2.00530E-15 2.00450E-15 2.00360E-15 2.00280E-15 2.00200E-15 &
     2.00120E-15 2.00050E-15 1.99970E-15 1.99900E-15 1.99820E-15 1.99750E-15 1.99670E-15 1.99600E-15 &
     1.99530E-15 1.99460E-15 1.99390E-15 1.99330E-15 1.99260E-15 1.99190E-15 1.99130E-15 1.99060E-15 &
     1.99000E-15 1.98940E-15 1.98880E-15 1.98820E-15 1.98760E-15 1.98700E-15 1.98640E-15 1.98590E-15 &
     1.98530E-15 1.98470E-15 1.98420E-15 1.98370E-15 1.98320E-15 1.98260E-15 1.98210E-15 1.98160E-15 &
     1.98110E-15 1.98070E-15 1.98020E-15 1.97970E-15 1.97930E-15 1.97880E-15 1.97840E-15 1.97790E-15 &
     1.97750E-15 1.97710E-15 1.97670E-15 1.97630E-15 1.97590E-15 1.97550E-15 1.97510E-15 1.97470E-15 &
     1.97430E-15 1.97400E-15 1.97360E-15 1.97330E-15 1.97290E-15 1.97260E-15 1.97230E-15 1.97200E-15 &
     1.97170E-15 1.97140E-15 1.97110E-15 1.97080E-15 1.97050E-15 1.97020E-15 1.96990E-15 1.96970E-15 &
     1.96940E-15 1.96910E-15 1.96890E-15 1.96860E-15 1.96840E-15 1.96820E-15 1.96800E-15 1.96770E-15 &
     1.96750E-15 1.96730E-15 1.96710E-15 1.96690E-15 1.96670E-15 1.96650E-15 1.96630E-15 1.96620E-15 &
     1.96600E-15 1.96580E-15 1.96570E-15 1.96550E-15 1.96540E-15 1.96520E-15 1.96510E-15 1.96490E-15 &
     1.96480E-15 1.96470E-15 1.96460E-15 1.96440E-15 1.96430E-15 1.96420E-15 1.96410E-15 1.96400E-15 &
     1.96390E-15 1.96380E-15 1.96370E-15 1.96370E-15 1.96360E-15 1.96350E-15 1.96340E-15 1.96340E-15 &
     1.96330E-15 1.96320E-15 1.96320E-15 1.96310E-15 1.96310E-15 1.96300E-15 1.96300E-15 1.96290E-15 &
     1.96290E-15 1.96290E-15 1.96280E-15 1.96280E-15 1.96280E-15 1.96270E-15 1.96270E-15 1.96270E-15 &
     1.96270E-15 1.96270E-15 1.96270E-15 1.96260E-15 1.96260E-15 1.96260E-15 1.96260E-15 1.96260E-15 &
     1.96260E-15 1.96260E-15 1.96260E-15 1.96270E-15 1.96270E-15 1.96270E-15 1.96270E-15 1.96270E-15 &
     1.96270E-15 1.96270E-15 1.96280E-15 1.96280E-15 1.96280E-15 1.96280E-15 1.96290E-15 1.96290E-15 &
     1.96290E-15 1.96290E-15 1.96300E-15 1.96300E-15 1.96300E-15 1.96310E-15 1.96310E-15 1.96310E-15 &
     1.96320E-15 1.96320E-15 1.96320E-15 1.96330E-15 1.96330E-15 1.96330E-15 1.96340E-15 1.96340E-15 &
     1.96340E-15 -6.11880E-12 -7.42060E-12 -7.38310E-12 -6.40570E-12 -4.84990E-12 -3.05920E-12 -1.34210E-12 &
     6.12690E-14 1.02440E-12 1.53600E-12 1.64260E-12 1.40900E-12 9.26590E-13 3.35520E-13 -2.07710E-13 &
     -5.99380E-13 -8.21820E-13 -9.11090E-13 -8.98000E-13 -7.87410E-13 -5.89100E-13 -3.48020E-13 -1.39080E-13 &
     -2.33070E-14 -1.09350E-14 -6.55950E-14 -1.37400E-13 -2.00350E-13 -2.58400E-13 -3.24600E-13 -3.99970E-13 &
     -4.61140E-13 -4.78700E-13 -4.36090E-13 -3.43750E-13 -2.34880E-13 -1.52670E-13 -1.04410E-13 -8.73860E-14 &
     -9.41980E-14 -9.87580E-14 -7.76510E-14 -3.65660E-14 1.00450E-14 4.16480E-14 4.05710E-14 6.52280E-15 &
     -5.36720E-14 -1.22060E-13 -1.79820E-13 -2.15710E-13 -2.23300E-13 -2.07740E-13 -1.77600E-13 -1.43420E-13 &
     -1.14760E-13 -9.53370E-14 -8.57790E-14 -8.02440E-14 -7.25770E-14 -6.18130E-14 -4.82770E-14 -3.22990E-14 &
     -1.47840E-14 1.04160E-15 1.13780E-14 1.24280E-14 1.79380E-15 -1.73480E-14 -4.04340E-14 -6.29040E-14 &
     -8.07030E-14 -9.17630E-14 -9.45150E-14 -8.73920E-14 -6.99470E-14 -4.62100E-14 -2.13250E-14 -4.32310E-16 &
     1.24090E-14 1.74410E-14 1.59810E-14 9.34420E-15 -1.15620E-15 -1.42070E-14 -2.84970E-14 -4.27700E-14 &
     -5.59930E-14 -6.71890E-14 -7.53800E-14 -7.95910E-14 -7.88460E-14 -7.21700E-14 -5.90160E-14 -4.05520E-14 &
     -1.83700E-14 5.93790E-15 3.07790E-14 5.45640E-14 7.57020E-14 9.29020E-14 1.06050E-13 1.15320E-13 &
     1.20910E-13 1.23000E-13 1.21750E-13 1.17370E-13 1.10060E-13 1.00230E-13 8.83080E-14 7.47330E-14 &
     5.99390E-14 4.43620E-14 2.84370E-14 1.26180E-14 -2.55800E-15 -1.65340E-14 -2.87550E-14 -3.86640E-14 &
     -4.57070E-14 -4.93260E-14 -4.91490E-14 -4.55320E-14 -3.90090E-14 -3.01180E-14 -1.93960E-14 -7.37670E-15 &
     5.40280E-15 1.84250E-14 3.12460E-14 4.34390E-14 5.45790E-14 6.42370E-14 7.19900E-14 7.74110E-14 &
     8.02110E-14 8.06480E-14 7.91170E-14 7.60120E-14 7.17270E-14 6.66570E-14 6.11950E-14 5.56970E-14 &
     5.03620E-14 4.53520E-14 4.08280E-14 3.69490E-14 3.38790E-14 3.17760E-14 3.07430E-14 3.06480E-14 &
     3.12990E-14 3.25060E-14 3.40760E-14 3.58190E-14 3.75440E-14 3.90820E-14 4.03520E-14 4.12950E-14 &
     4.18520E-14 4.19660E-14 4.15760E-14 4.06250E-14 3.90780E-14 3.69960E-14 3.44650E-14 3.15690E-14 &
     2.83950E-14 2.50270E-14 2.15500E-14 1.80460E-14 1.45820E-14 1.12180E-14 8.01900E-15 5.04580E-15 &
     2.36090E-15 2.66830E-17 -1.90440E-15 -3.41930E-15 -4.51480E-15 -5.18750E-15 -5.43420E-15 -5.25160E-15 &
     -4.63660E-15 -3.59530E-15 -2.17290E-15 -4.23900E-16 1.59720E-15 3.83570E-15 6.23720E-15 8.74710E-15 &
     1.13110E-14 1.38770E-14 1.63930E-14 1.88060E-14 2.10650E-14 2.31170E-14 2.49110E-14 2.64040E-14 &
     2.75950E-14 2.84900E-14 2.90980E-14 2.94260E-14 2.94810E-14 2.92730E-14 2.88130E-14 2.81390E-14 &
     2.72910E-14 2.63120E-14 2.52440E-14 2.41270E-14 2.30050E-14 2.19170E-14 2.08990E-14 1.99880E-14 &
     1.92170E-14 1.86220E-14 1.82380E-14 1.81000E-14 1.82330E-14 1.86210E-14 1.92370E-14 2.00540E-14 &
     2.10470E-14 2.21890E-14 2.34530E-14 2.48130E-14 2.62420E-14 2.77110E-14 2.91930E-14 3.06600E-14 &
     3.20860E-14 3.34410E-14 3.46990E-14 3.58380E-14 3.68340E-14 3.76670E-14 3.83130E-14 3.87520E-14 &
     3.89600E-14 3.89230E-14 3.86560E-14 3.81820E-14 3.75220E-14 3.66990E-14 3.57360E-14 3.46540E-14 &
     3.34760E-14 3.22200E-14 3.09020E-14 2.95400E-14 2.81520E-14 2.67540E-14 2.53630E-14 2.39980E-14 &
     2.26770E-14 2.14180E-14 2.02390E-14 1.91600E-14 1.81990E-14 1.73750E-14 1.67000E-14 1.61660E-14 &
     1.57610E-14 1.54700E-14 1.52790E-14 1.51750E-14 1.51440E-14 1.51740E-14 1.52550E-14 1.53780E-14 &
     1.55340E-14 1.57140E-14 1.59100E-14 1.61130E-14 1.63130E-14 1.65020E-14 1.66690E-14 1.68050E-14 &
     1.68990E-14 1.69420E-14 1.69250E-14 1.68380E-14 1.66810E-14 1.64520E-14 1.61530E-14 1.57820E-14 &
     1.53400E-14 1.48270E-14 1.42450E-14 1.36050E-14 1.29230E-14 1.22140E-14 1.14910E-14 1.07680E-14 &
     1.00620E-14 9.38240E-15 8.72970E-15 8.10010E-15 7.49030E-15 6.89680E-15 6.31620E-15 5.74500E-15 &
     5.18250E-15 4.63890E-15 4.12680E-15 3.65930E-15 3.24900E-15 2.90880E-15 2.65150E-15 2.48540E-15 &
     2.40070E-15 2.38310E-15 2.41820E-15 2.49180E-15 2.58940E-15 2.69680E-15 2.80240E-15 2.90570E-15 &
     3.00880E-15 3.11390E-15 3.22330E-15 3.33900E-15 3.46340E-15 3.59800E-15 3.74270E-15 3.89690E-15 &
     4.05970E-15 4.23050E-15 4.40860E-15 4.59330E-15 4.78250E-15 4.96830E-15 5.14150E-15 5.29280E-15 &
     5.41300E-15 5.49270E-15 5.52290E-15 5.49760E-15 5.42550E-15 5.31830E-15 5.18820E-15 5.04710E-15 &
     4.90700E-15 4.77980E-15 4.67470E-15 4.59000E-15 4.52100E-15 4.46310E-15 4.41180E-15 4.36230E-15 &
     4.31030E-15 4.25130E-15 4.18270E-15 4.10210E-15 4.00710E-15 3.89520E-15 3.76420E-15 3.61170E-15 &
     3.43770E-15 3.25230E-15 3.06820E-15 2.89800E-15 2.75420E-15 2.64940E-15 2.59620E-15 2.60260E-15 &
     2.65790E-15 2.74650E-15 2.85310E-15 2.96240E-15 3.05880E-15 3.12700E-15 3.15480E-15 3.14290E-15 &
     3.09500E-15 3.01510E-15 2.90690E-15 2.77440E-15 2.62130E-15 2.45280E-15 2.27990E-15 2.11450E-15 &
     1.96890E-15 1.85520E-15 1.78550E-15 1.77210E-15 1.82190E-15 1.92230E-15 2.05520E-15 2.20270E-15 &
     2.34690E-15 2.46990E-15 2.55380E-15 2.58510E-15 2.56740E-15 2.50870E-15 2.41710E-15 2.30070E-15 &
     2.16730E-15 2.02520E-15 1.88230E-15 1.74720E-15 1.62830E-15 1.53430E-15 1.47370E-15 1.45500E-15 &
     1.48680E-15 1.57390E-15 1.70680E-15 1.87210E-15 2.05660E-15 2.24710E-15 2.43020E-15 2.59260E-15 &
     2.72310E-15 2.81800E-15 2.87590E-15 2.89500E-15 2.87380E-15 2.81070E-15 2.70400E-15 2.55490E-15 &
     2.37480E-15 2.17790E-15 1.97860E-15 1.79090E-15 1.62910E-15 1.50730E-15 1.43650E-15 1.41370E-15 &
     1.43300E-15 1.48810E-15 1.57290E-15 1.68120E-15 1.80690E-15 1.94350E-15 2.08260E-15 2.21580E-15 &
     2.33440E-15 2.42980E-15 2.49340E-15 2.51660E-15 2.49310E-15 2.42610E-15 2.32100E-15 2.18330E-15 &
     2.01860E-15 1.83230E-15 1.63000E-15 1.41770E-15 1.20440E-15 9.99660E-16 8.12980E-16 6.53940E-16 &
     5.32110E-16 4.57040E-16 4.35450E-16 4.62760E-16 5.31530E-16 6.34340E-16 7.63770E-16 9.12380E-16 &
     1.07270E-15 1.23700E-15 1.39550E-15 1.53810E-15 1.65470E-15 1.73520E-15 1.76950E-15 1.74750E-15 &
     1.66280E-15 1.52450E-15 1.34510E-15 1.13750E-15 9.14300E-16 6.88250E-16 4.72070E-16 2.76350E-16 &
     1.03300E-16 -4.70160E-17 -1.74510E-16 -2.79100E-16 -3.60700E-16 -4.19230E-16 -4.54860E-16 -4.68750E-16 &
     -4.62320E-16 -4.36980E-16 -3.94130E-16 -3.35200E-16 -2.61600E-16 -1.74840E-16 -7.68260E-17 3.04130E-17 &
     1.44860E-16 2.64510E-16 3.87340E-16 5.11340E-16 6.33600E-16 7.47700E-16 8.46320E-16 9.22150E-16 &
     9.67890E-16 9.76210E-16 9.39820E-16 8.54200E-16 7.26140E-16 5.65210E-16 3.81010E-16 1.83110E-16 &
     -1.88930E-17 -2.15420E-16 -3.98480E-16 -5.66420E-16 -7.19200E-16 -8.56760E-16 -9.79060E-16 -1.08600E-15 &
     -1.17760E-15 -1.25310E-15 -1.30870E-15 -1.34000E-15 -1.34250E-15 -1.31170E-15 -1.24340E-15 -1.13290E-15 &
     -9.78630E-16 -7.89930E-16 -5.78900E-16 -3.57670E-16 -1.38350E-16 6.69410E-17 2.46090E-16 3.90200E-16 &
     5.03240E-16 5.92400E-16 6.64880E-16 7.27860E-16 7.88540E-16 8.54100E-16 9.29340E-16 1.00950E-15 &
     1.08740E-15 1.15580E-15 1.20770E-15 1.23570E-15 1.23290E-15 1.19410E-15 1.12360E-15 1.02750E-15 &
     9.12180E-16 7.83940E-16 6.49050E-16 5.13820E-16 3.82990E-16 2.55190E-16 1.27490E-16 -3.04080E-18 &
     -1.39310E-16 -2.84250E-16 -4.40790E-16 -6.09920E-16 -7.85020E-16 -9.57550E-16 -1.11900E-15 -1.26070E-15 &
     -1.37420E-15 -1.45100E-15 -1.48500E-15 -1.48000E-15 -1.44240E-15 -1.37870E-15 -1.29510E-15 -1.19800E-15 &
     -1.09370E-15 -9.87560E-16 -8.80080E-16 -7.70740E-16 -6.58990E-16 -5.44290E-16 -4.26080E-16 -3.03820E-16 &
     -1.77240E-16 -4.72240E-17 8.50730E-17 2.18490E-16 3.51880E-16 4.84070E-16 6.13900E-16 7.40420E-16 &
     8.63420E-16 9.82900E-16 1.09890E-15 1.21130E-15 1.32010E-15 1.42550E-15 1.52710E-15 1.62400E-15 &
     1.71530E-15 1.79960E-15 1.87610E-15 1.94360E-15 2.00090E-15 2.04740E-15 2.08360E-15 2.11060E-15 &
     2.12920E-15 2.14050E-15 2.14540E-15 2.14490E-15 2.13980E-15 2.13060E-15 2.11750E-15 2.10100E-15 &
     2.08140E-15 2.05890E-15 2.03400E-15 2.00710E-15 1.97950E-15 1.95260E-15 1.92770E-15 1.90630E-15 &
     1.88980E-15 1.87960E-15 1.87620E-15 1.87720E-15 1.87930E-15 1.87910E-15 1.87330E-15 1.85870E-15 &
     1.83180E-15 1.79070E-15 1.73890E-15 1.68120E-15 1.62220E-15 1.56680E-15 1.51980E-15 1.48600E-15 &
     1.46820E-15 1.46190E-15 1.46070E-15 1.45810E-15 1.44760E-15 1.42300E-15 1.37760E-15 1.30760E-15 &
     1.21880E-15 1.11950E-15 1.01810E-15 9.22750E-16 8.41880E-16 7.83770E-16 7.53600E-16 7.44020E-16 &
     7.44550E-16 7.44700E-16 7.33990E-16 7.01940E-16 6.38060E-16 5.38470E-16 4.25680E-16 3.28800E-16 &
     2.76940E-16 2.99200E-16 4.24710E-16 6.82550E-16 1.07860E-15 1.52570E-15 1.91330E-15 2.13120E-15 &
     2.06870E-15 1.61560E-15 6.61450E-16 -8.15780E-16 -2.48410E-15 -3.92300E-15 -4.71200E-15 -4.43080E-15 &
     -2.65890E-15 1.02410E-15 6.70570E-15 1.31410E-14 1.87520E-14 2.19620E-14 2.11920E-14 1.48650E-14 &
     1.40340E-15 -1.95260E-14 -4.32770E-14 -6.39580E-14 -7.56780E-14 -7.25460E-14 -4.86720E-14 1.83530E-15 &
     8.23450E-14 1.86130E-13 3.03960E-13 4.26580E-13 5.44740E-13 6.49210E-13 7.30740E-13 7.82580E-13 &
     8.07950E-13 8.12560E-13 8.02110E-13 7.82310E-13 7.58890E-13 7.37530E-13 7.22770E-13 7.14370E-13 &
     7.10900E-13 7.10920E-13 7.13020E-13 7.15760E-13 7.17720E-13 7.17770E-13 7.15970E-13 7.12700E-13 &
     7.08310E-13 7.03180E-13 6.97660E-13 6.92130E-13 6.86880E-13 6.81910E-13 6.77140E-13 6.72490E-13 &
     6.67900E-13 6.63280E-13 6.58560E-13 6.53680E-13 6.48660E-13 6.43510E-13 6.38270E-13 6.32970E-13 &
     6.27620E-13 6.22260E-13 6.16920E-13 6.11570E-13 6.06230E-13 6.00880E-13 5.95510E-13 5.90110E-13 &
     5.84680E-13 5.79210E-13 5.73700E-13 5.68170E-13 5.62620E-13 5.57060E-13 5.51500E-13 5.45950E-13 &
     5.40410E-13 5.34880E-13 5.29330E-13 5.23770E-13 5.18170E-13 5.12520E-13 5.06820E-13 5.01050E-13 &
     4.95210E-13 4.89320E-13 4.83380E-13 4.77390E-13 4.71360E-13 4.65300E-13 4.59200E-13 4.53070E-13 &
     4.46890E-13 4.40650E-13 4.34340E-13 4.27950E-13 4.21460E-13 4.14870E-13 4.08190E-13 4.01450E-13 &
     3.94670E-13 3.87860E-13 3.81040E-13 3.74240E-13 3.67470E-13 3.60730E-13 3.54030E-13 3.47360E-13 &
     3.40720E-13 3.34100E-13 3.27510E-13 3.20950E-13 3.14430E-13 3.07950E-13 3.01520E-13 2.95160E-13 &
     2.88860E-13 2.82650E-13 2.76530E-13 2.70470E-13 2.64460E-13 2.58480E-13 2.52500E-13 2.46500E-13 &
     2.40470E-13 2.34370E-13 2.28220E-13 2.21990E-13 2.15690E-13 2.09320E-13 2.02870E-13 1.96340E-13 &
     1.89730E-13 1.83050E-13 1.76320E-13 1.69570E-13 1.62820E-13 1.56080E-13 1.49380E-13 1.42720E-13 &
     1.36130E-13 1.29620E-13 1.23180E-13 1.16820E-13 1.10570E-13 1.04420E-13 9.83730E-14 9.24360E-14 &
     8.65990E-14 8.08590E-14 7.52080E-14 6.96410E-14 6.41540E-14 5.87350E-14 5.33560E-14 4.79800E-14 &
     4.25730E-14 3.71010E-14 3.15290E-14 2.58220E-14 1.99580E-14 1.39690E-14 7.90090E-15 1.79850E-15 &
     -4.29310E-15 -1.03290E-14 -1.62630E-14 -2.20570E-14 -2.76990E-14 -3.31820E-14 -3.85000E-14 -4.36480E-14 &
     -4.86190E-14 -5.34080E-14 -5.80160E-14 -6.24740E-14 -6.68230E-14 -7.11010E-14 -7.53480E-14 -7.96050E-14 &
     -8.39090E-14 -8.82930E-14 -9.27570E-14 -9.72940E-14 -1.01900E-13 -1.06550E-13 -1.11260E-13 -1.16010E-13 &
     -1.20790E-13 -1.25560E-13 -1.30270E-13 -1.34870E-13 -1.39320E-13 -1.43580E-13 -1.47590E-13 -1.51320E-13 &
     -1.54820E-13 -1.58130E-13 -1.61300E-13 -1.64390E-13 -1.67440E-13 -1.70500E-13 -1.73630E-13 -1.76820E-13 &
     -1.80090E-13 -1.83430E-13 -1.86850E-13 -1.90350E-13 -1.93920E-13 -1.97570E-13 -2.01240E-13 -2.04860E-13 &
     -2.08370E-13 -2.11690E-13 -2.14780E-13 -2.17550E-13 -2.19970E-13 -2.22110E-13 -2.24050E-13 -2.25890E-13 &
     -2.27720E-13 -2.29620E-13 -2.31710E-13 -2.34020E-13 -2.36530E-13 -2.39180E-13 -2.41880E-13 -2.44570E-13 &
     -2.47200E-13 -2.49680E-13 -2.51970E-13 -2.54060E-13 -2.55990E-13 -2.57760E-13 -2.59390E-13 -2.60890E-13 &
     -2.62300E-13 -2.63610E-13 -2.64870E-13 -2.66080E-13 -2.67270E-13 -2.68470E-13 -2.69700E-13 -2.70970E-13 &
     -2.72310E-13 -2.73640E-13 -2.74880E-13 -2.75970E-13 -2.76820E-13 -2.77350E-13 -2.77490E-13 -2.77190E-13 &
     -2.76560E-13 -2.75720E-13 -2.74820E-13 -2.73980E-13 -2.73340E-13 -2.73030E-13 -2.73150E-13 -2.73600E-13 &
     -2.74250E-13 -2.74980E-13 -2.75650E-13 -2.76120E-13 -2.76270E-13 -2.76000E-13 -2.75380E-13 -2.74540E-13 &
     -2.73590E-13 -2.72630E-13 -2.71790E-13 -2.71180E-13 -2.70880E-13 -2.70810E-13 -2.70870E-13 -2.70960E-13 &
     -2.70950E-13 -2.70750E-13 -2.70250E-13 -2.69370E-13 -2.68190E-13 -2.66810E-13 -2.65340E-13 -2.63870E-13 &
     -2.62530E-13 -2.61400E-13 -2.60570E-13 -2.59970E-13 -2.59500E-13 -2.59060E-13 -2.58550E-13 -2.57880E-13 &
     -2.56940E-13 -2.55660E-13 -2.54110E-13 -2.52370E-13 -2.50530E-13 -2.48670E-13 -2.46890E-13 -2.45260E-13 &
     -2.43870E-13 -2.42650E-13 -2.41550E-13 -2.40510E-13 -2.39460E-13 -2.38340E-13 -2.37070E-13 -2.35620E-13 &
     -2.33990E-13 -2.32230E-13 -2.30370E-13 -2.28430E-13 -2.26450E-13 -2.24470E-13 -2.22510E-13 -2.20590E-13 &
     -2.18700E-13 -2.16850E-13 -2.15050E-13 -2.13290E-13 -2.11590E-13 -2.09940E-13 -2.08320E-13 -2.06690E-13 &
     -2.05030E-13 -2.03310E-13 -2.01490E-13 -1.99560E-13 -1.97470E-13 -1.95280E-13 -1.93000E-13 -1.90680E-13 &
     -1.88360E-13 -1.86070E-13 -1.83840E-13 -1.81720E-13 -1.79670E-13 -1.77670E-13 -1.75700E-13 -1.73730E-13 &
     -1.71730E-13 -1.69680E-13 -1.67540E-13 -1.65290E-13 -1.62900E-13 -1.60330E-13 -1.57560E-13 -1.54550E-13 &
     -1.51270E-13 -1.47730E-13 -1.44010E-13 -1.40250E-13 -1.36590E-13 -1.33140E-13 -1.30050E-13 -1.27430E-13 &
     -1.25380E-13 -1.23790E-13 -1.22520E-13 -1.21410E-13 -1.20320E-13 -1.19100E-13 -1.17600E-13 -1.15710E-13 &
     -1.13470E-13 -1.10950E-13 -1.08220E-13 -1.05360E-13 -1.02440E-13 -9.95260E-14 -9.66860E-14 -9.39190E-14 &
     -9.12130E-14 -8.85560E-14 -8.59360E-14 -8.33410E-14 -8.07590E-14 -7.81800E-14 -7.56040E-14 -7.30310E-14 &
     -7.04630E-14 -6.79010E-14 -6.53460E-14 -6.27990E-14 -6.02620E-14 -5.77370E-14 -5.52270E-14 -5.27340E-14 &
     -5.02600E-14 -4.78080E-14 -4.53810E-14 -4.29810E-14 -4.06060E-14 -3.82540E-14 -3.59230E-14 -3.36120E-14 &
     -3.13190E-14 -2.90420E-14 -2.67790E-14 -2.45340E-14 -2.23100E-14 -2.01100E-14 -1.79370E-14 -1.57950E-14 &
     -1.36880E-14 -1.16160E-14 -9.57400E-15 -7.55170E-15 -5.54080E-15 -3.53240E-15 -1.51790E-15 5.11720E-16 &
     2.55730E-15 4.58850E-15 6.56710E-15 8.45480E-15 1.02130E-14 1.18040E-14 1.31900E-14 1.43450E-14 &
     1.52980E-14 1.60900E-14 1.67650E-14 1.73620E-14 1.79250E-14 1.84950E-14 1.91050E-14 1.97520E-14 &
     2.04250E-14 2.11110E-14 2.17980E-14 2.24760E-14 2.31320E-14 2.37580E-14 2.43700E-14 2.49850E-14 &
     2.56210E-14 2.62980E-14 2.70330E-14 2.78460E-14 2.87470E-14 2.97180E-14 3.07360E-14 3.17750E-14 &
     3.28100E-14 3.38180E-14 3.47720E-14 3.56540E-14 3.64660E-14 3.72160E-14 3.79110E-14 3.85590E-14 &
     3.91670E-14 3.97430E-14 4.02940E-14 4.08180E-14 4.13120E-14 4.17750E-14 4.22040E-14 4.25970E-14 &
     4.29500E-14 4.32620E-14 4.35330E-14 4.37630E-14 4.39520E-14 4.41000E-14 4.42080E-14 4.42760E-14 &
     4.43030E-14 4.42920E-14 4.42430E-14 4.41570E-14 4.40350E-14 4.38790E-14 4.36900E-14 4.34680E-14 &
     4.32120E-14 4.29180E-14 4.25860E-14 4.22120E-14 4.17940E-14 4.13300E-14 4.08200E-14 4.02670E-14 &
     3.96780E-14 3.90590E-14 3.84160E-14 3.77560E-14 3.70850E-14 3.64070E-14 3.57250E-14 3.50390E-14 &
     3.43510E-14 3.36600E-14 3.29680E-14 3.22750E-14 3.15820E-14 3.08880E-14 3.01890E-14 2.94860E-14 &
     2.87750E-14 2.80550E-14 2.73240E-14 2.65800E-14 2.58250E-14 2.50590E-14 2.42840E-14 2.35020E-14 &
     2.27120E-14 2.19170E-14 2.11160E-14 2.02960E-14 1.94450E-14 1.85480E-14 1.75930E-14 1.65640E-14 &
     1.54490E-14 1.42380E-14 1.29410E-14 1.15710E-14 1.01420E-14 8.66710E-15 7.16120E-15 5.63740E-15 &
     4.10700E-15 2.57010E-15 1.02410E-15 -5.33490E-16 -2.10530E-15 -3.69380E-15 -5.30160E-15 -6.93080E-15 &
     -8.58190E-15 -1.02550E-14 -1.19500E-14 -1.36670E-14 -1.54070E-14 -1.71690E-14 -1.89520E-14 -2.07540E-14 &
     -2.25720E-14 -2.44010E-14 -2.62370E-14 -2.80780E-14 -2.99200E-14 -3.17610E-14 -3.36090E-14 -3.54740E-14 &
     -3.73660E-14 -3.92960E-14 -4.12720E-14 -4.33060E-14 -4.53990E-14 -4.75210E-14 -4.96320E-14 -5.16950E-14 &
     -5.36700E-14 -5.55180E-14 -5.72020E-14 -5.86950E-14 -6.00190E-14 -6.12120E-14 -6.23080E-14 -6.33430E-14 &
     -6.43540E-14 -6.53760E-14 -6.64370E-14 -6.75360E-14 -6.86660E-14 -6.98190E-14 -7.09870E-14 -7.21610E-14 &
     -7.33330E-14 -7.44970E-14 -7.56530E-14 -7.68000E-14 -7.79390E-14 -7.90710E-14 -8.01970E-14 -8.13160E-14 &
     -8.24310E-14 -8.35460E-14 -8.46670E-14 -8.58000E-14 -8.69500E-14 -8.81240E-14 -8.93270E-14 -9.05640E-14 &
     -9.18400E-14 -9.31590E-14 -9.45260E-14 -9.59460E-14 -9.74210E-14 -9.89570E-14 -1.00560E-13 -1.02200E-13 &
     -1.03890E-13 -1.05600E-13 -1.07310E-13 -1.09020E-13 -1.10700E-13 -1.12360E-13 -1.13980E-13 -1.15570E-13 &
     -1.17140E-13 -1.18700E-13 -1.20240E-13 -1.21760E-13 -1.23280E-13 -1.24800E-13 -1.26290E-13 -1.27780E-13 &
     -1.29240E-13 -1.30670E-13 -1.32080E-13 -1.33460E-13 -1.34820E-13 -1.36160E-13 -1.37490E-13 -1.38820E-13 &
     -1.40160E-13 -1.41510E-13 -1.42880E-13 -1.44230E-13 -1.45540E-13 -1.46770E-13 -1.47870E-13 -1.48830E-13 &
     -1.49600E-13 -1.50160E-13 -1.50530E-13 -1.50760E-13 -1.50880E-13 -1.50940E-13 -1.50960E-13 -1.50990E-13 &
     -1.51060E-13 -1.51160E-13 -1.51290E-13 -1.51440E-13 -1.51590E-13 -1.51730E-13 -1.51860E-13 -1.51970E-13 &
     -1.52060E-13 -1.52130E-13 -1.52190E-13 -1.52230E-13 -1.52270E-13 -1.52300E-13 -1.52330E-13 -1.52360E-13 &
     -1.52380E-13 -1.52390E-13 -1.52390E-13 -1.52380E-13 -1.52350E-13 -1.52300E-13 -1.52250E-13 -1.52180E-13 &
     -1.52110E-13 -1.52050E-13 -1.51990E-13 -1.51950E-13 -1.51920E-13 -1.51910E-13 -1.51910E-13 -1.51910E-13 &
     -1.51900E-13 -1.51900E-13 -1.51880E-13 -1.51840E-13 -1.51790E-13 -1.51720E-13 -1.51620E-13 -1.51500E-13 &
     -1.51350E-13 -1.51180E-13 -1.50970E-13 -1.50730E-13 -1.50460E-13 -1.50170E-13 -1.49840E-13 -1.49480E-13 &
     -1.49090E-13 -1.48670E-13 -1.48230E-13 -1.47760E-13 -1.47260E-13 -1.46750E-13 -1.46210E-13 -1.45660E-13 &
     -1.45090E-13 -1.44480E-13 -1.43840E-13 -1.43140E-13 -1.42380E-13 -1.41530E-13 -1.40580E-13 -1.39530E-13 &
     -1.38390E-13 -1.37170E-13 -1.35910E-13 -1.34610E-13 -1.33290E-13 -1.31980E-13 -1.30690E-13 -1.29410E-13 &
     -1.28150E-13 -1.26900E-13 -1.25640E-13 -1.24390E-13 -1.23120E-13 -1.21840E-13 -1.20550E-13 -1.19250E-13 &
     -1.17940E-13 -1.16630E-13 -1.15310E-13 -1.13980E-13 -1.12650E-13 -1.11320E-13 -1.09980E-13 -1.08630E-13 &
     -1.07290E-13 -1.05930E-13 -1.04570E-13 -1.03200E-13 -1.01820E-13 -1.00430E-13 -9.90420E-14 -9.76460E-14 &
     -9.62460E-14 -9.48420E-14 -9.34370E-14 -9.20250E-14 -9.06060E-14 -8.91740E-14 -8.77280E-14 -8.62640E-14 &
     -8.47790E-14 -8.32710E-14 -8.17390E-14 -8.01850E-14 -7.86100E-14 -7.70150E-14 -7.54000E-14 -7.37670E-14 &
     -7.21170E-14 -7.04510E-14 -6.87690E-14 -6.70730E-14 -6.53640E-14 -6.36430E-14 -6.19110E-14 -6.01700E-14 &
     -5.84180E-14 -5.66560E-14 -5.48840E-14 -5.31000E-14 -5.13060E-14 -4.95010E-14 -4.76850E-14 -4.58600E-14 &
     -4.40300E-14 -4.21970E-14 -4.03640E-14 -3.85350E-14 -3.67130E-14 -3.48990E-14 -3.30930E-14 -3.12960E-14 &
     -2.95070E-14 -2.77240E-14 -2.59480E-14 -2.41770E-14 -2.24110E-14 -2.06510E-14 -1.88960E-14 -1.71460E-14 &
     -1.54000E-14 -1.36600E-14 -1.19240E-14 -1.01930E-14 -8.46720E-15 -6.74610E-15 -5.02980E-15 -3.31860E-15 &
     -1.61240E-15 8.86900E-17 1.78470E-15 3.47550E-15 5.16130E-15 6.84180E-15 8.51720E-15 1.01880E-14 &
     1.18530E-14 1.35130E-14 1.51680E-14 1.68190E-14 1.84670E-14 2.01120E-14 2.17550E-14 2.33950E-14 &
     2.50340E-14 2.66710E-14 2.83040E-14 2.99320E-14 3.15540E-14 3.31690E-14 3.47760E-14 3.63730E-14 &
     3.79650E-14 3.95530E-14 4.11420E-14 4.27350E-14 4.43340E-14 4.59440E-14 4.75650E-14 4.91970E-14 &
     5.08330E-14 5.24690E-14 5.41030E-14 5.57280E-14 5.73410E-14 5.89380E-14 6.05190E-14 6.20860E-14 &
     6.36390E-14 6.51780E-14 6.67060E-14 6.82220E-14 6.97260E-14 7.12130E-14 7.26760E-14 7.41080E-14 &
     7.55030E-14 7.68520E-14 7.81510E-14 7.93930E-14 8.05850E-14 8.17350E-14 8.28510E-14 8.39410E-14 &
     8.50120E-14 8.60730E-14 8.71310E-14 8.81840E-14 8.92320E-14 9.02720E-14 9.13020E-14 9.23220E-14 &
     9.33280E-14 9.43190E-14 9.52960E-14 9.62590E-14 9.72090E-14 9.81460E-14 9.90710E-14 9.99840E-14 &
     1.00890E-13 1.01780E-13 1.02660E-13 1.03520E-13 1.04380E-13 1.05220E-13 1.06050E-13 1.06870E-13 &
     1.07680E-13 1.08470E-13 1.09250E-13 1.10020E-13 1.10780E-13 1.11530E-13 1.12270E-13 1.13000E-13 &
     1.13710E-13 1.14420E-13 1.15110E-13 1.15790E-13 1.16450E-13 1.17100E-13 1.17740E-13 1.18370E-13 &
     1.18980E-13 1.19590E-13 1.20190E-13 1.20780E-13 1.21360E-13 1.21920E-13 1.22480E-13 1.23000E-13 &
     1.23500E-13 1.23970E-13 1.24390E-13 1.24770E-13 1.25110E-13 1.25410E-13 1.25680E-13 1.25910E-13 &
     1.26120E-13 1.26310E-13 1.26480E-13 1.26630E-13 1.26760E-13 1.26850E-13 1.26920E-13 1.26960E-13 &
     1.26960E-13 1.26930E-13 1.26860E-13 1.26770E-13 1.26650E-13 1.26510E-13 1.26360E-13 1.26190E-13 &
     1.26030E-13 1.25850E-13 1.25670E-13 1.25480E-13 1.25280E-13 1.25080E-13 1.24860E-13 1.24620E-13 &
     1.24380E-13 1.24120E-13 1.23860E-13 1.23580E-13 1.23290E-13 1.22990E-13 1.22680E-13 1.22360E-13 &
     1.22040E-13 1.21700E-13 1.21350E-13 1.20990E-13 1.20620E-13 1.20230E-13 1.19840E-13 1.19440E-13 &
     1.19020E-13 1.18600E-13 1.18160E-13 1.17720E-13 1.17260E-13 1.16790E-13 1.16320E-13 1.15830E-13 &
     1.15330E-13 1.14820E-13 1.14300E-13 1.13770E-13 1.13230E-13 1.12680E-13 1.12110E-13 1.11540E-13 &
     1.10960E-13 1.10360E-13 1.09760E-13 1.09150E-13 1.08520E-13 1.07890E-13 1.07240E-13 1.06580E-13 &
     1.05910E-13 1.05230E-13 1.04520E-13 1.03800E-13 1.03040E-13 1.02250E-13 1.01410E-13 1.00530E-13 &
     9.95950E-14 9.86160E-14 9.76040E-14 9.65690E-14 9.55200E-14 9.44690E-14 9.34250E-14 9.23960E-14 &
     9.13800E-14 9.03760E-14 8.93780E-14 8.83830E-14 8.73880E-14 8.63890E-14 8.53840E-14 8.43730E-14 &
     8.33560E-14 8.23350E-14 8.13100E-14 8.02840E-14 7.92560E-14 7.82270E-14 7.71970E-14 7.61670E-14 &
     7.51360E-14 7.41030E-14 7.30690E-14 7.20340E-14 7.09960E-14 6.99560E-14 6.89150E-14 6.78730E-14 &
     6.68290E-14 6.57840E-14 6.47390E-14 6.36930E-14 6.26460E-14 6.15980E-14 6.05500E-14 5.95000E-14 &
     5.84500E-14 5.73980E-14 5.63450E-14 5.52910E-14 5.42360E-14 5.31800E-14 5.21240E-14 5.10680E-14 &
     5.00120E-14 4.89570E-14 4.79010E-14 4.68450E-14 4.57890E-14 4.47330E-14 4.36770E-14 4.26200E-14 &
     4.15630E-14 4.05050E-14 3.94480E-14 3.83920E-14 3.73370E-14 3.62840E-14 3.52340E-14 3.41860E-14 &
     3.31390E-14 3.20940E-14 3.10470E-14 2.99990E-14 2.89490E-14 2.78940E-14 2.68350E-14 2.57730E-14 &
     2.47100E-14 2.36490E-14 2.25910E-14 2.15400E-14 2.04980E-14 1.94650E-14 1.84430E-14 1.74310E-14 &
     1.64280E-14 1.54340E-14 1.44480E-14 1.34710E-14 1.25010E-14 1.15390E-14 1.05860E-14 9.63960E-15 &
     8.70170E-15 7.77190E-15 6.85020E-15 5.93680E-15 5.03170E-15 4.13480E-15 3.24610E-15 2.36560E-15 &
     1.49330E-15 6.29220E-16 -2.26660E-16 -1.07430E-15 -1.91370E-15 -2.74470E-15 -3.56730E-15 -4.38150E-15 &
     -5.18730E-15 -5.98450E-15 -6.77310E-15 -7.55320E-15 -8.32470E-15 -9.08760E-15 -9.84200E-15 -1.05880E-14 &
     -1.13250E-14 -1.20530E-14 -1.27720E-14 -1.34830E-14 -1.41840E-14 -1.48750E-14 -1.55570E-14 -1.62290E-14 &
     -1.68920E-14 -1.75450E-14 -1.81890E-14 -1.88230E-14 -1.94480E-14 -2.00640E-14 -2.06700E-14 -2.12680E-14 &
     -2.18560E-14 -2.24340E-14 -2.30030E-14 -2.35630E-14 -2.41130E-14 -2.46530E-14 -2.51830E-14 -2.57030E-14 &
     -2.62120E-14 -2.67110E-14 -2.71990E-14 -2.76760E-14 -2.81410E-14 -2.85940E-14 -2.90340E-14 -2.94610E-14 &
     -2.98720E-14 -3.02680E-14 -3.06480E-14 -3.10100E-14 -3.13560E-14 -3.16860E-14 -3.20010E-14 -3.23010E-14 &
     -3.25870E-14 -3.28590E-14 -3.31180E-14 -3.33650E-14 -3.36000E-14 -3.38240E-14 -3.40370E-14 -3.42390E-14 &
     -3.44310E-14 -3.46130E-14 -3.47860E-14 -3.49490E-14 -3.51010E-14 -3.52440E-14 -3.53760E-14 -3.54970E-14 &
     -3.56080E-14 -3.57080E-14 -3.57980E-14 -3.58780E-14 -3.59470E-14 -3.60060E-14 -3.60550E-14 -3.60950E-14 &
     -3.61240E-14 -3.61440E-14 -3.61550E-14 -3.61560E-14 -3.61480E-14 -3.61300E-14 -3.61030E-14 -3.60670E-14 &
     -3.60220E-14 -3.59670E-14 -3.59030E-14 -3.58290E-14 -3.57470E-14 -3.56540E-14 -3.55530E-14 -3.54420E-14 &
     -3.53220E-14 -3.51930E-14 -3.50540E-14 -3.49070E-14 -3.47510E-14 -3.45860E-14 -3.44120E-14 -3.42290E-14 &
     -3.40370E-14 -3.38360E-14 -3.36260E-14 -3.34070E-14 -3.31790E-14 -3.29430E-14 -3.26990E-14 -3.24480E-14 &
     -3.21900E-14 -3.19250E-14 -3.16550E-14 -3.13820E-14 -3.11100E-14 -3.08420E-14 -3.05810E-14 -3.03300E-14 &
     -3.00930E-14 -2.98720E-14 -2.96600E-14 -2.94510E-14 -2.92370E-14 -2.90110E-14 -2.87650E-14 -2.84930E-14 &
     -2.81880E-14 -2.78530E-14 -2.74920E-14 -2.71100E-14 -2.67100E-14 -2.62980E-14 -2.58760E-14 -2.54490E-14 &
     -2.50170E-14 -2.45810E-14 -2.41400E-14 -2.36960E-14 -2.32490E-14 -2.27970E-14 -2.23430E-14 -2.18860E-14 &
     -2.14250E-14 -2.09620E-14 -2.04950E-14 -2.00250E-14 -1.95530E-14 -1.90770E-14 -1.85990E-14 -1.81180E-14 &
     -1.76340E-14 -1.71480E-14 -1.66590E-14 -1.61680E-14 -1.56760E-14 -1.51810E-14 -1.46840E-14 -1.41850E-14 &
     -1.36840E-14 -1.31820E-14 -1.26770E-14 -1.21710E-14 -1.16630E-14 -1.11540E-14 -1.06430E-14 -1.01320E-14 &
     -9.61930E-15 -9.10610E-15 -8.59220E-15 -8.07740E-15 -7.56110E-15 -7.04300E-15 -6.52240E-15 -5.99910E-15 &
     -5.47250E-15 -4.94270E-15 -4.41120E-15 -3.88030E-15 -3.35200E-15 -2.82840E-15 -2.31160E-15 -1.80380E-15 &
     -1.30690E-15 -8.22350E-16 -3.51370E-16 1.04760E-16 5.44790E-16 9.67440E-16 1.37150E-15 1.75620E-15 &
     2.12360E-15 2.47620E-15 2.81640E-15 3.14680E-15 3.47000E-15 3.78850E-15 4.10490E-15 4.42270E-15 &
     4.74560E-15 5.07730E-15 5.42120E-15 5.78120E-15 6.16070E-15 6.56210E-15 6.98210E-15 7.41600E-15 &
     7.85900E-15 8.30660E-15 8.75390E-15 9.19640E-15 9.63020E-15 1.00550E-14 1.04730E-14 1.08830E-14 &
     1.12870E-14 1.16860E-14 1.20810E-14 1.24720E-14 1.28590E-14 1.32430E-14 1.36220E-14 1.39970E-14 &
     1.43670E-14 1.47320E-14 1.50910E-14 1.54460E-14 1.57950E-14 1.61400E-14 1.64790E-14 1.68140E-14 &
     1.71430E-14 1.74670E-14 1.77870E-14 1.81010E-14 1.84090E-14 1.87120E-14 1.90090E-14 1.93000E-14 &
     1.95850E-14 1.98650E-14 2.01390E-14 2.04080E-14 2.06720E-14 2.09330E-14 2.11900E-14 2.14430E-14 &
     2.16920E-14 2.19340E-14 2.21670E-14 2.23900E-14 2.26000E-14 2.27960E-14 2.29770E-14 2.31420E-14 &
     2.32900E-14 2.34230E-14 2.35390E-14 2.36400E-14 2.37240E-14 2.37920E-14 2.38470E-14 2.38910E-14 &
     2.39250E-14 2.39520E-14 2.39750E-14 2.39950E-14 2.40150E-14 2.40350E-14 2.40530E-14 2.40700E-14 &
     2.40840E-14 2.40960E-14 2.41040E-14 2.41080E-14 2.41090E-14 2.41060E-14 2.41000E-14 2.40910E-14 &
     2.40800E-14 2.40670E-14 2.40520E-14 2.40350E-14 2.40130E-14 2.39880E-14 2.39570E-14 2.39210E-14 &
     2.38770E-14 2.38260E-14 2.37680E-14 2.37020E-14 2.36300E-14 2.35500E-14 2.34650E-14 2.33730E-14 &
     2.32760E-14 2.31720E-14 2.30630E-14 2.29480E-14 2.28270E-14 2.26990E-14 2.25650E-14 2.24250E-14 &
     2.22790E-14 2.21270E-14 2.19690E-14 2.18050E-14 2.16350E-14 2.14590E-14 2.12780E-14 2.10910E-14 &
     2.08990E-14 2.07010E-14 2.04980E-14 2.02890E-14 2.00740E-14 1.98540E-14 1.96290E-14 1.93970E-14 &
     1.91590E-14 1.89140E-14 1.86630E-14 1.84050E-14 1.81390E-14 1.78670E-14 1.75880E-14 1.73050E-14 &
     1.70160E-14 1.67240E-14 1.64280E-14 1.61300E-14 1.58300E-14 1.55270E-14 1.52220E-14 1.49150E-14 &
     1.46060E-14 1.42960E-14 1.39840E-14 1.36710E-14 1.33560E-14 1.30400E-14 1.27220E-14 1.24030E-14 &
     1.20820E-14 1.17600E-14 1.14360E-14 1.11100E-14 1.07820E-14 1.04520E-14 1.01190E-14 9.78410E-15 &
     9.44620E-15 9.10470E-15 8.75850E-15 8.40680E-15 8.04860E-15 7.68290E-15 7.30870E-15 6.92560E-15 &
     6.53450E-15 6.13670E-15 5.73360E-15 5.32650E-15 4.91680E-15 4.50580E-15 4.09460E-15 3.68310E-15 &
     3.27110E-15 2.85800E-15 2.44380E-15 2.02780E-15 1.60990E-15 1.18990E-15 7.67740E-16 3.43870E-16 &
     -8.14940E-17 -5.08060E-16 -9.35570E-16 -1.36370E-15 -1.79240E-15 -2.22180E-15 -2.65220E-15 -3.08390E-15 &
     -3.51730E-15 -3.95280E-15 -4.39060E-15 -4.83070E-15 -5.27210E-15 -5.71350E-15 -6.15350E-15 -6.59060E-15 &
     -7.02370E-15 -7.45130E-15 -7.87220E-15 -8.28630E-15 -8.69340E-15 -9.09350E-15 -9.48650E-15 -9.87250E-15 &
     -1.02510E-14 -1.06230E-14 -1.09880E-14 -1.13490E-14 -1.17050E-14 -1.20580E-14 -1.24080E-14 -1.27580E-14 &
     -1.31070E-14 -1.34550E-14 -1.38030E-14 -1.41500E-14 -1.44960E-14 -1.48400E-14 -1.51830E-14 -1.55240E-14 &
     -1.58640E-14 -1.62010E-14 -1.65370E-14 -1.68720E-14 -1.72040E-14 -1.75350E-14 -1.78640E-14 -1.81920E-14 &
     -1.85190E-14 -1.88450E-14 -1.91710E-14 -1.94970E-14 -1.98240E-14 -2.01520E-14 -2.04800E-14 -2.08080E-14 &
     -2.11370E-14 -2.14670E-14 -2.17960E-14 -2.21250E-14 -2.24540E-14 -2.27820E-14 -2.31090E-14 -2.34340E-14 &
     -2.37560E-14 -2.40760E-14 -2.43910E-14 -2.47020E-14 -2.50090E-14 -2.53120E-14 -2.56110E-14 -2.59050E-14 &
     -2.61950E-14 -2.64810E-14 -2.67620E-14 -2.70400E-14 -2.73150E-14 -2.75860E-14 -2.78550E-14 -2.81210E-14 &
     -2.83850E-14 -2.86470E-14 -2.89060E-14 -2.91590E-14 -2.94050E-14 -2.96430E-14 -2.98690E-14 -3.00840E-14 &
     -3.02840E-14 -3.04720E-14 -3.06460E-14 -3.08090E-14 -3.09600E-14 -3.11010E-14 -3.12330E-14 -3.13550E-14 &
     -3.14700E-14 -3.15780E-14 -3.16800E-14 -3.17770E-14 -3.18700E-14 -3.19610E-14 -3.20500E-14 -3.21370E-14 &
     -3.22220E-14 -3.23040E-14 -3.23840E-14 -3.24610E-14 -3.25350E-14 -3.26050E-14 -3.26730E-14 -3.27370E-14 &
     -3.27970E-14 -3.28550E-14 -3.29100E-14 -3.29610E-14 -3.30100E-14 -3.30550E-14 -3.30980E-14 -3.31380E-14 &
     -3.31750E-14 -3.32090E-14 -3.32420E-14 -3.32710E-14 -3.32980E-14 -3.33220E-14 -3.33420E-14 -3.33580E-14 &
     -3.33690E-14 -3.33750E-14 -3.33760E-14 -3.33710E-14 -3.33610E-14 -3.33440E-14 -3.33210E-14 -3.32930E-14 &
     -3.32580E-14 -3.32160E-14 -3.31690E-14 -3.31160E-14 -3.30560E-14 -3.29920E-14 -3.29220E-14 -3.28470E-14 &
     -3.27670E-14 -3.26830E-14 -3.25940E-14 -3.25010E-14 -3.24040E-14 -3.23030E-14 -3.21990E-14 -3.20910E-14 &
     -3.19790E-14 -3.18610E-14 -3.17380E-14 -3.16070E-14 -3.14680E-14 -3.13200E-14 -3.11630E-14 -3.09960E-14 &
     -3.08220E-14 -3.06400E-14 -3.04510E-14 -3.02570E-14 -3.00580E-14 -2.98550E-14 -2.96480E-14 -2.94380E-14 &
     -2.92260E-14 -2.90100E-14 -2.87930E-14 -2.85740E-14 -2.83540E-14 -2.81320E-14 -2.79080E-14 -2.76830E-14 &
     -2.74570E-14 -2.72290E-14 -2.69980E-14 -2.67660E-14 -2.65320E-14 -2.62970E-14 -2.60590E-14 -2.58200E-14 &
     -2.55790E-14 -2.53360E-14 -2.50920E-14 -2.48460E-14 -2.45980E-14 -2.43480E-14 -2.40960E-14 -2.38420E-14 &
     -2.35850E-14 -2.33270E-14 -2.30660E-14 -2.28040E-14 -2.25400E-14 -2.22750E-14 -2.20100E-14 -2.17440E-14 &
     -2.14790E-14 -2.12120E-14 -2.09430E-14 -2.06700E-14 -2.03920E-14 -2.01090E-14 -1.98180E-14 -1.95190E-14 &
     -1.92130E-14 -1.89000E-14 -1.85820E-14 -1.82590E-14 -1.79330E-14 -1.76040E-14 -1.72730E-14 -1.69410E-14 &
     -1.66060E-14 -1.62710E-14 -1.59330E-14 -1.55950E-14 -1.52540E-14 -1.49120E-14 -1.45690E-14 -1.42250E-14 &
     -1.38800E-14 -1.35340E-14 -1.31890E-14 -1.28440E-14 -1.24990E-14 -1.21550E-14 -1.18120E-14 -1.14700E-14 &
     -1.11290E-14 -1.07890E-14 -1.04510E-14 -1.01150E-14 -9.78030E-15 -9.44700E-15 -9.11510E-15 -8.78420E-15 &
     -8.45440E-15 -8.12540E-15 -7.79710E-15 -7.46970E-15 -7.14300E-15 -6.81730E-15 -6.49250E-15 -6.16860E-15 &
     -5.84580E-15 -5.52410E-15 -5.20350E-15 -4.88400E-15 -4.56550E-15 -4.24800E-15 -3.93150E-15 -3.61600E-15 &
     -3.30140E-15 -2.98780E-15 -2.67530E-15 -2.36380E-15 -2.05340E-15 -1.74410E-15 -1.43600E-15 -1.12910E-15 &
     -8.23310E-16 -5.18640E-16 -2.15050E-16 8.75290E-17 3.89120E-16 6.89780E-16 9.89520E-16 1.28820E-15 &
     1.58590E-15 1.88220E-15 2.17720E-15 2.47080E-15 2.76280E-15 3.05310E-15 3.34210E-15 3.62990E-15 &
     3.91680E-15 4.20310E-15 4.48910E-15 4.77500E-15 5.06100E-15 5.34690E-15 5.63220E-15 5.91670E-15 &
     6.19980E-15 6.48130E-15 6.76090E-15 7.03790E-15 7.31190E-15 7.58220E-15 7.84810E-15 8.10880E-15 &
     8.36380E-15 8.61230E-15 8.85390E-15 9.08890E-15 9.31770E-15 9.54080E-15 9.75880E-15 9.97200E-15 &
     1.01810E-14 1.03860E-14 1.05880E-14 1.07860E-14 1.09810E-14 1.11730E-14 1.13610E-14 1.15470E-14 &
     1.17310E-14 1.19110E-14 1.20890E-14 1.22640E-14 1.24370E-14 1.26060E-14 1.27730E-14 1.29370E-14 &
     1.30990E-14 1.32570E-14 1.34130E-14 1.35670E-14 1.37170E-14 1.38650E-14 1.40100E-14 1.41530E-14 &
     1.42920E-14 1.44290E-14 1.45640E-14 1.46950E-14 1.48230E-14 1.49490E-14 1.50710E-14 1.51910E-14 &
     1.53070E-14 1.54210E-14 1.55330E-14 1.56410E-14 1.57470E-14 1.58510E-14 1.59510E-14 1.60490E-14 &
     1.61440E-14 1.62360E-14 1.63250E-14 1.64110E-14 1.64940E-14 1.65740E-14 1.66520E-14 1.67280E-14 &
     1.68010E-14 1.68720E-14 1.69400E-14 1.70060E-14 1.70690E-14 1.71260E-14 1.71780E-14 1.72230E-14 &
     1.72610E-14 1.72910E-14 1.73130E-14 1.73280E-14 1.73360E-14 1.73380E-14 1.73350E-14 1.73270E-14 &
     1.73150E-14 1.72990E-14 1.72800E-14 1.72570E-14 1.72310E-14 1.72020E-14 1.71700E-14 1.71360E-14 &
     1.70990E-14 1.70590E-14 1.70180E-14 1.69730E-14 1.69270E-14 1.68780E-14 1.68270E-14 1.67740E-14 &
     1.67180E-14 1.66610E-14 1.66010E-14 1.65390E-14 1.64750E-14 1.64090E-14 1.63410E-14 1.62710E-14 &
     1.61990E-14 1.61240E-14 1.60480E-14 1.59700E-14 1.58890E-14 1.58070E-14 1.57230E-14 1.56370E-14 &
     1.55490E-14 1.54590E-14 1.53670E-14 1.52730E-14 1.51770E-14 1.50800E-14 1.49810E-14 1.48800E-14 &
     1.47770E-14 1.46720E-14 1.45660E-14 1.44580E-14 1.43480E-14 1.42360E-14 1.41220E-14 1.40070E-14 &
     1.38890E-14 1.37700E-14 1.36490E-14 1.35260E-14 1.34020E-14 1.32770E-14 1.31500E-14 1.30210E-14 &
     1.28920E-14 1.27610E-14 1.26280E-14 1.24940E-14 1.23560E-14 1.22170E-14 1.20740E-14 1.19280E-14 &
     1.17790E-14 1.16270E-14 1.14710E-14 1.13110E-14 1.11470E-14 1.09800E-14 1.08080E-14 1.06340E-14 &
     1.04570E-14 1.02780E-14 1.00990E-14 9.91930E-15 9.74080E-15 9.56390E-15 9.38830E-15 9.21380E-15 &
     9.04010E-15 8.86680E-15 8.69370E-15 8.52040E-15 8.34680E-15 8.17280E-15 7.99860E-15 7.82430E-15 &
     7.64990E-15 7.47560E-15 7.30150E-15 7.12760E-15 6.95400E-15 6.78060E-15 6.60740E-15 6.43440E-15 &
     6.26150E-15 6.08880E-15 5.91620E-15 5.74370E-15 5.57140E-15 5.39920E-15 5.22730E-15 5.05570E-15 &
     4.88430E-15 4.71320E-15 4.54250E-15 4.37210E-15 4.20210E-15 4.03250E-15 3.86320E-15 3.69440E-15 &
     3.52590E-15 3.35790E-15 3.19020E-15 3.02300E-15 2.85610E-15 2.68960E-15 2.52350E-15 2.35770E-15 &
     2.19240E-15 2.02750E-15 1.86320E-15 1.69950E-15 1.53650E-15 1.37420E-15 1.21270E-15 1.05180E-15 &
     8.91370E-16 7.31260E-16 5.71290E-16 4.11300E-16 2.51130E-16 9.07540E-17 -6.93030E-17 -2.28370E-16 &
     -3.85780E-16 -5.40880E-16 -6.92980E-16 -8.41430E-16 -9.85770E-16 -1.12640E-15 -1.26400E-15 -1.39920E-15 &
     -1.53270E-15 -1.66500E-15 -1.79680E-15 -1.92860E-15 -2.06020E-15 -2.19130E-15 -2.32150E-15 -2.45050E-15 &
     -2.57800E-15 -2.70370E-15 -2.82720E-15 -2.94850E-15 -3.06790E-15 -3.18520E-15 -3.30070E-15 -3.41430E-15 &
     -3.52630E-15 -3.63650E-15 -3.74510E-15 -3.85200E-15 -3.95730E-15 -4.06080E-15 -4.16270E-15 -4.26280E-15 &
     -4.36120E-15 -4.45790E-15 -4.55290E-15 -4.64610E-15 -4.73770E-15 -4.82750E-15 -4.91570E-15 -5.00210E-15 &
     -5.08690E-15 -5.17000E-15 -5.25140E-15 -5.33110E-15 -5.40910E-15 -5.48550E-15 -5.56010E-15 -5.63310E-15 &
     -5.70440E-15 -5.77400E-15 -5.84190E-15 -5.90810E-15 -5.97270E-15 -6.03550E-15 -6.09660E-15 -6.15590E-15 &
     -6.21340E-15 -6.26910E-15 -6.32290E-15 -6.37480E-15 -6.42480E-15 -6.47300E-15 -6.51940E-15 -6.56430E-15 &
     -6.60770E-15 -6.64970E-15 -6.69050E-15 -6.73000E-15 -6.76810E-15 -6.80410E-15 -6.83770E-15 -6.86840E-15 &
     -6.89590E-15 -6.91960E-15 -6.93930E-15 -6.95510E-15 -6.96750E-15 -6.97670E-15 -6.98320E-15 -6.98720E-15 &
     -6.98910E-15 -6.98920E-15 -6.98750E-15 -6.98380E-15 -6.97810E-15 -6.97040E-15 -6.96050E-15 -6.94830E-15 &
     -6.93390E-15 -6.91730E-15 -6.89860E-15 -6.87800E-15 -6.85550E-15 -6.83140E-15 -6.80580E-15 -6.77870E-15 &
     -6.75020E-15 -6.72020E-15 -6.68870E-15 -6.65580E-15 -6.62140E-15 -6.58550E-15 -6.54810E-15 -6.50910E-15 &
     -6.46860E-15 -6.42660E-15 -6.38310E-15 -6.33820E-15 -6.29170E-15 -6.24380E-15 -6.19440E-15 -6.14370E-15 &
     -6.09150E-15 -6.03810E-15 -5.98330E-15 -5.92730E-15 -5.87000E-15 -5.81150E-15 -5.75170E-15 -5.69070E-15 &
     -5.62830E-15 -5.56460E-15 -5.49960E-15 -5.43320E-15 -5.36540E-15 -5.29640E-15 -5.22610E-15 -5.15450E-15 &
     -5.08170E-15 -5.00780E-15 -4.93270E-15 -4.85630E-15 -4.77870E-15 -4.69970E-15 -4.61920E-15 -4.53720E-15 &
     -4.45360E-15 -4.36840E-15 -4.28180E-15 -4.19420E-15 -4.10580E-15 -4.01710E-15 -3.92830E-15 -3.83970E-15 &
     -3.75160E-15 -3.66400E-15 -3.57690E-15 -3.49020E-15 -3.40390E-15 -3.31800E-15 -3.23230E-15 -3.14690E-15 &
     -3.06110E-15 -2.97470E-15 -2.88700E-15 -2.79760E-15 -2.70600E-15 -2.61160E-15 -2.51430E-15 -2.41420E-15 &
     -2.31180E-15 -2.20750E-15 -2.10170E-15 -1.99500E-15 -1.88770E-15 -1.78030E-15 -1.67270E-15 -1.56500E-15 &
     -1.45720E-15 -1.34930E-15 -1.24130E-15 -1.13320E-15 -1.02500E-15 -9.16690E-16 -8.08320E-16 -6.99870E-16 &
     -5.91330E-16 -4.82720E-16 -3.74040E-16 -2.65290E-16 -1.56480E-16 -4.76570E-17 6.11740E-17 1.69990E-16 &
     2.78760E-16 3.87460E-16 4.96090E-16 6.04640E-16 7.13140E-16 8.21580E-16 9.29990E-16 1.03840E-15 &
     1.14680E-15 1.25510E-15 1.36340E-15 1.47160E-15 1.57940E-15 1.68700E-15 1.79410E-15 1.90060E-15 &
     2.00670E-15 2.11230E-15 2.21780E-15 2.32350E-15 2.42950E-15 2.53600E-15 2.64340E-15 2.75170E-15 &
     2.86010E-15 2.96780E-15 3.07390E-15 3.17740E-15 3.27760E-15 3.37350E-15 3.46440E-15 3.55100E-15 &
     3.63400E-15 3.71420E-15 3.79240E-15 3.86930E-15 3.94590E-15 4.02260E-15 4.09960E-15 4.17640E-15 &
     4.25310E-15 4.32920E-15 4.40470E-15 4.47940E-15 4.55290E-15 4.62570E-15 4.69780E-15 4.76960E-15 &
     4.84130E-15 4.91320E-15 4.98550E-15 5.05830E-15 5.13150E-15 5.20480E-15 5.27780E-15 5.35020E-15 &
     5.42180E-15 5.49220E-15 5.56120E-15 5.62880E-15 5.69510E-15 5.76000E-15 5.82380E-15 5.88630E-15 &
     5.94780E-15 6.00820E-15 6.06750E-15 6.12580E-15 6.18290E-15 6.23900E-15 6.29400E-15 6.34780E-15 &
     6.40050E-15 6.45210E-15 6.50260E-15 6.55190E-15 6.60020E-15 6.64740E-15 6.69350E-15 6.73860E-15 &
     6.78260E-15 6.82530E-15 6.86680E-15 6.90700E-15 6.94580E-15 6.98320E-15 7.01910E-15 7.05360E-15 &
     7.08700E-15 7.11950E-15 7.15130E-15 7.18250E-15 7.21340E-15 7.24400E-15 7.27370E-15 7.30180E-15 &
     7.32750E-15 7.35020E-15 7.36900E-15 7.38320E-15 7.39230E-15 7.39700E-15 7.39790E-15 7.39580E-15 &
     7.39170E-15 7.38620E-15 7.38020E-15 7.37440E-15 7.36860E-15 7.36260E-15 7.35630E-15 7.34950E-15 &
     7.34180E-15 7.33320E-15 7.32340E-15 7.31260E-15 7.30070E-15 7.28790E-15 7.27430E-15 7.26000E-15 &
     7.24500E-15 7.22950E-15 7.21320E-15 7.19590E-15 7.17730E-15 7.15730E-15 7.13560E-15 7.11190E-15 &
     7.08610E-15 7.05820E-15 7.02860E-15 6.99740E-15 6.96480E-15 6.93100E-15 6.89620E-15 6.86050E-15 &
     6.82390E-15 6.78650E-15 6.74800E-15 6.70860E-15 6.66820E-15 6.62670E-15 6.58400E-15 6.54030E-15 &
     6.49540E-15 6.44950E-15 6.40250E-15 6.35460E-15 6.30560E-15 6.25560E-15 6.20470E-15 6.15280E-15 &
     6.09990E-15 6.04600E-15 5.99110E-15 5.93530E-15 5.87850E-15 5.82070E-15 5.76200E-15 5.70250E-15 &
     5.64210E-15 5.58100E-15 5.51910E-15 5.45660E-15 5.39320E-15 5.32890E-15 5.26340E-15 5.19670E-15 &
     5.12860E-15 5.05900E-15 4.98780E-15 4.91510E-15 4.84130E-15 4.76660E-15 4.69120E-15 4.61540E-15 &
     4.53950E-15 4.46370E-15 4.38790E-15 4.31200E-15 4.23610E-15 4.16000E-15 4.08370E-15 4.00700E-15 &
     3.93000E-15 3.85270E-15 3.77500E-15 3.69710E-15 3.61880E-15 3.54030E-15 3.46160E-15 3.38270E-15 &
     3.30360E-15 3.22420E-15 3.14450E-15 3.06440E-15 2.98390E-15 2.90290E-15 2.82140E-15 2.73940E-15 &
     2.65660E-15 2.57300E-15 2.48840E-15 2.40290E-15 2.31610E-15 2.22810E-15 2.13910E-15 2.04930E-15 &
     1.95890E-15 1.86810E-15 1.77710E-15 1.68630E-15 1.59570E-15 1.50530E-15 1.41510E-15 1.32500E-15 &
     1.23500E-15 1.14490E-15 1.05480E-15 9.64660E-16 8.74400E-16 7.84120E-16 6.93870E-16 6.03690E-16 &
     5.13640E-16 4.23760E-16 3.34090E-16 2.44590E-16 1.55200E-16 6.58610E-17 -2.34760E-17 -1.12870E-16 &
     -2.02380E-16 -2.91990E-16 -3.81470E-16 -4.70520E-16 -5.58830E-16 -6.46110E-16 -7.32040E-16 -8.16330E-16 &
     -8.98770E-16 -9.79530E-16 -1.05890E-15 -1.13700E-15 -1.21430E-15 -1.29100E-15 -1.36730E-15 -1.44340E-15 &
     -1.51930E-15 -1.59500E-15 -1.67040E-15 -1.74540E-15 -1.82010E-15 -1.89420E-15 -1.96780E-15 -2.04080E-15 &
     -2.11330E-15 -2.18530E-15 -2.25670E-15 -2.32770E-15 -2.39820E-15 -2.46820E-15 -2.53770E-15 -2.60670E-15 &
     -2.67530E-15 -2.74340E-15 -2.81100E-15 -2.87810E-15 -2.94480E-15 -3.01100E-15 -3.07680E-15 -3.14210E-15 &
     -3.20700E-15 -3.27140E-15 -3.33540E-15 -3.39900E-15 -3.46200E-15 -3.52450E-15 -3.58630E-15 -3.64750E-15 &
     -3.70780E-15 -3.76740E-15 -3.82600E-15 -3.88370E-15 -3.94050E-15 -3.99650E-15 -4.05160E-15 -4.10580E-15 &
     -4.15920E-15 -4.21180E-15 -4.26350E-15 -4.31440E-15 -4.36440E-15 -4.41350E-15 -4.46160E-15 -4.50880E-15 &
     -4.55500E-15 -4.60020E-15 -4.64460E-15 -4.68810E-15 -4.73100E-15 -4.77310E-15 -4.81450E-15 -4.85540E-15 &
     -4.89530E-15 -4.93400E-15 -4.97120E-15 -5.00650E-15 -5.03970E-15 -5.07030E-15 -5.09830E-15 -5.12380E-15 &
     -5.14730E-15 -5.16900E-15 -5.18950E-15 -5.20890E-15 -5.22780E-15 -5.24640E-15 -5.26470E-15 -5.28250E-15 &
     -5.29990E-15 -5.31680E-15 -5.33290E-15 -5.34840E-15 -5.36300E-15 -5.37690E-15 -5.38990E-15 -5.40220E-15 &
     -5.41380E-15 -5.42460E-15 -5.43470E-15 -5.44420E-15 -5.45290E-15 -5.46090E-15 -5.46820E-15 -5.47480E-15 &
     -5.48070E-15 -5.48580E-15 -5.49010E-15 -5.49370E-15 -5.49670E-15 -5.49900E-15 -5.50070E-15 -5.50190E-15 &
     -5.50250E-15 -5.50270E-15 -5.50210E-15 -5.50080E-15 -5.49860E-15 -5.49520E-15 -5.49060E-15 -5.48450E-15 &
     -5.47680E-15 -5.46760E-15 -5.45700E-15 -5.44480E-15 -5.43140E-15 -5.41660E-15 -5.40050E-15 -5.38320E-15 &
     -5.36480E-15 -5.34540E-15 -5.32490E-15 -5.30350E-15 -5.28120E-15 -5.25820E-15 -5.23440E-15 -5.20990E-15 &
     -5.18460E-15 -5.15870E-15 -5.13200E-15 -5.10460E-15 -5.07650E-15 -5.04760E-15 -5.01800E-15 -4.98750E-15 &
     -4.95600E-15 -4.92360E-15 -4.89000E-15 -4.85540E-15 -4.81950E-15 -4.78260E-15 -4.74480E-15 -4.70620E-15 &
     -4.66710E-15 -4.62750E-15 -4.58760E-15 -4.54760E-15 -4.50730E-15 -4.46690E-15 -4.42620E-15 -4.38530E-15 &
     -4.34400E-15 -4.30230E-15 -4.26020E-15 -4.21770E-15 -4.17480E-15 -4.13160E-15 -4.08810E-15 -4.04420E-15 &
     -4.00000E-15 -3.95550E-15 -3.91070E-15 -3.86560E-15 -3.82030E-15 -3.77460E-15 -3.72870E-15 -3.68260E-15 &
     -3.63620E-15 -3.58950E-15 -3.54250E-15 -3.49520E-15 -3.44750E-15 -3.39950E-15 -3.35110E-15 -3.30240E-15 &
     -3.25330E-15 -3.20390E-15 -3.15430E-15 -3.10470E-15 -3.05490E-15 -3.00520E-15 -2.95560E-15 -2.90570E-15 &
     -2.85560E-15 -2.80480E-15 -2.75330E-15 -2.70080E-15 -2.64700E-15 -2.59180E-15 -2.53540E-15 -2.47790E-15 &
     -2.41960E-15 -2.36070E-15 -2.30130E-15 -2.24180E-15 -2.18220E-15 -2.12260E-15 -2.06300E-15 -2.00340E-15 &
     -1.94360E-15 -1.88380E-15 -1.82390E-15 -1.76390E-15 -1.70390E-15 -1.64400E-15 -1.58440E-15 -1.52520E-15 &
     -1.46650E-15 -1.40850E-15 -1.35120E-15 -1.29450E-15 -1.23840E-15 -1.18290E-15 -1.12770E-15 -1.07280E-15 &
     -1.01820E-15 -9.63740E-16 -9.09420E-16 -8.55270E-16 -8.01310E-16 -7.47570E-16 -6.94050E-16 -6.40790E-16 &
     -5.87790E-16 -5.35060E-16 -4.82590E-16 -4.30390E-16 -3.78460E-16 -3.26780E-16 -2.75370E-16 -2.24220E-16 &
     -1.73340E-16 -1.22720E-16 -7.23660E-17 -2.22900E-17 2.75100E-17 7.70300E-17 1.26270E-16 1.75220E-16 &
     2.23880E-16 2.72240E-16 3.20320E-16 3.68090E-16 4.15560E-16 4.62730E-16 5.09600E-16 5.56170E-16 &
     6.02430E-16 6.48400E-16 6.94070E-16 7.39440E-16 7.84520E-16 8.29300E-16 8.73770E-16 9.17930E-16 &
     9.61770E-16 1.00530E-15 1.04850E-15 1.09130E-15 1.13390E-15 1.17610E-15 1.21800E-15 1.25970E-15 &
     1.30110E-15 1.34220E-15 1.38310E-15 1.42370E-15 1.46400E-15 1.50380E-15 1.54320E-15 1.58190E-15 &
     1.62010E-15 1.65760E-15 1.69430E-15 1.73010E-15 1.76490E-15 1.79860E-15 1.83110E-15 1.86220E-15 &
     1.89200E-15 1.92040E-15 1.94780E-15 1.97420E-15 1.99980E-15 2.02470E-15 2.04920E-15 2.07330E-15 &
     2.09710E-15 2.12040E-15 2.14330E-15 2.16570E-15 2.18760E-15 2.20900E-15 2.22970E-15 2.24980E-15 &
     2.26940E-15 2.28840E-15 2.30680E-15 2.32470E-15 2.34210E-15 2.35890E-15 2.37530E-15 2.39110E-15 &
     2.40640E-15 2.42120E-15 2.43560E-15 2.44940E-15 2.46270E-15 2.47560E-15 2.48790E-15 2.49980E-15 &
     2.51110E-15 2.52200E-15 2.53230E-15 2.54220E-15 2.55150E-15 2.56040E-15 2.56870E-15 2.57660E-15 &
     2.58390E-15 2.59080E-15 2.59710E-15 2.60300E-15 2.60840E-15 2.61320E-15 2.61760E-15 2.62140E-15 &
     2.62480E-15 2.62760E-15 2.63000E-15 2.63190E-15 2.63330E-15 2.63430E-15 2.63490E-15 2.63500E-15 &
     2.63480E-15 2.63400E-15 2.63250E-15 2.63020E-15 2.62690E-15 2.62250E-15 2.61680E-15 2.60970E-15 &
     2.60130E-15 2.59180E-15 2.58120E-15 2.56970E-15 2.55730E-15 2.54430E-15 2.53070E-15 2.51650E-15 &
     2.50190E-15 2.48680E-15 2.47130E-15 2.45540E-15 2.43920E-15 2.42260E-15 2.40580E-15 2.38860E-15 &
     2.37110E-15 2.35330E-15 2.33520E-15 2.31680E-15 2.29800E-15 2.27890E-15 2.25950E-15 2.23980E-15 &
     2.21980E-15 2.19940E-15 2.17870E-15 2.15770E-15 2.13640E-15 2.11480E-15 2.09290E-15 2.07070E-15 &
     2.04830E-15 2.02560E-15 2.00260E-15 1.97940E-15 1.95600E-15 1.93230E-15 1.90830E-15 1.88410E-15 &
     1.85960E-15 1.83490E-15 1.80990E-15 1.78470E-15 1.75920E-15 1.73350E-15 1.70750E-15 1.68130E-15 &
     1.65490E-15 1.62820E-15 1.60130E-15 1.57420E-15 1.54680E-15 1.51910E-15 1.49120E-15 1.46310E-15 &
     1.43480E-15 1.40620E-15 1.37740E-15 1.34830E-15 1.31910E-15 1.28960E-15 1.26000E-15 1.23010E-15 &
     1.20000E-15 1.16970E-15 1.13910E-15 1.10830E-15 1.07720E-15 1.04580E-15 1.01420E-15 9.82280E-16 &
     9.50150E-16 9.17790E-16 8.85220E-16 8.52430E-16 8.19450E-16 7.86330E-16 7.53120E-16 7.19860E-16 &
     6.86610E-16 6.53430E-16 6.20350E-16 5.87420E-16 5.54650E-16 5.22020E-16 4.89520E-16 4.57140E-16 &
     4.24880E-16 3.92710E-16 3.60640E-16 3.28670E-16 2.96790E-16 2.65020E-16 2.33370E-16 2.01830E-16 &
     1.70420E-16 1.39130E-16 1.07970E-16 7.69370E-17 4.60370E-17 1.52690E-17 -1.53670E-17 -4.58690E-17 &
     -7.62370E-17 -1.06470E-16 -1.36580E-16 -1.66560E-16 -1.96410E-16 -2.26130E-16 -2.55740E-16 -2.85220E-16 &
     -3.14570E-16 -3.43790E-16 -3.72870E-16 -4.01790E-16 -4.30560E-16 -4.59160E-16 -4.87590E-16 -5.15850E-16 &
     -5.43940E-16 -5.71860E-16 -5.99620E-16 -6.27220E-16 -6.54650E-16 -6.81930E-16 -7.09060E-16 -7.36030E-16 &
     -7.62860E-16 -7.89540E-16 -8.16080E-16 -8.42490E-16 -8.68760E-16 -8.94860E-16 -9.20740E-16 -9.46390E-16 &
     -9.71750E-16 -9.96780E-16 -1.02150E-15 -1.04570E-15 -1.06960E-15 -1.09290E-15 -1.11560E-15 -1.13780E-15 &
     -1.15920E-15 -1.18000E-15 -1.20000E-15 -1.21920E-15 -1.23790E-15 -1.25610E-15 -1.27380E-15 -1.29110E-15 &
     -1.30820E-15 -1.32510E-15 -1.34180E-15 -1.35820E-15 -1.37440E-15 -1.39030E-15 -1.40580E-15 -1.42100E-15 &
     -1.43580E-15 -1.45020E-15 -1.46430E-15 -1.47800E-15 -1.49130E-15 -1.50420E-15 -1.51680E-15 -1.52910E-15 &
     -1.54100E-15 -1.55250E-15 -1.56380E-15 -1.57460E-15 -1.58520E-15 -1.59540E-15 -1.60520E-15 -1.61470E-15 &
     -1.62390E-15 -1.63280E-15 -1.64130E-15 -1.64950E-15 -1.65740E-15 -1.66500E-15 -1.67220E-15 -1.67910E-15 &
     -1.68560E-15 -1.69180E-15 -1.69770E-15 -1.70330E-15 -1.70860E-15 -1.71350E-15 -1.71800E-15 -1.72230E-15 &
     -1.72620E-15 -1.72970E-15 -1.73300E-15 -1.73590E-15 -1.73840E-15 -1.74060E-15 -1.74250E-15 -1.74410E-15 &
     -1.74530E-15 -1.74630E-15 -1.74690E-15 -1.74720E-15 -1.74710E-15 -1.74670E-15 -1.74580E-15 -1.74450E-15 &
     -1.74260E-15 -1.74030E-15 -1.73750E-15 -1.73420E-15 -1.73030E-15 -1.72590E-15 -1.72110E-15 -1.71560E-15 &
     -1.70970E-15 -1.70330E-15 -1.69640E-15 -1.68900E-15 -1.68110E-15 -1.67280E-15 -1.66400E-15 -1.65490E-15 &
     -1.64530E-15 -1.63540E-15 -1.62500E-15 -1.61420E-15 -1.60310E-15 -1.59160E-15 -1.57970E-15 -1.56740E-15 &
     -1.55480E-15 -1.54200E-15 -1.52880E-15 -1.51540E-15 -1.50180E-15 -1.48790E-15 -1.47390E-15 -1.45970E-15 &
     -1.44520E-15 -1.43050E-15 -1.41560E-15 -1.40050E-15 -1.38520E-15 -1.36970E-15 -1.35390E-15 -1.33790E-15 &
     -1.32170E-15 -1.30530E-15 -1.28870E-15 -1.27180E-15 -1.25480E-15 -1.23750E-15 -1.22000E-15 -1.20240E-15 &
     -1.18450E-15 -1.16650E-15 -1.14830E-15 -1.12980E-15 -1.11130E-15 -1.09250E-15 -1.07360E-15 -1.05440E-15 &
     -1.03520E-15 -1.01570E-15 -9.96070E-16 -9.76270E-16 -9.56300E-16 -9.36150E-16 -9.15830E-16 -8.95340E-16 &
     -8.74660E-16 -8.53820E-16 -8.32810E-16 -8.11650E-16 -7.90330E-16 -7.68870E-16 -7.47270E-16 -7.25550E-16 &
     -7.03710E-16 -6.81790E-16 -6.59790E-16 -6.37750E-16 -6.15670E-16 -5.93580E-16 -5.71500E-16 -5.49430E-16 &
     -5.27390E-16 -5.05370E-16 -4.83380E-16 -4.61430E-16 -4.39520E-16 -4.17660E-16 -3.95790E-16 -3.73870E-16 &
     -3.51850E-16 -3.29680E-16 -3.07300E-16 -2.84660E-16 -2.61740E-16 -2.38560E-16 -2.15190E-16 -1.91670E-16 &
     -1.68070E-16 -1.44420E-16 -1.20790E-16 -9.72280E-17 -7.37340E-17 -5.03110E-17 -2.69580E-17 -3.67490E-18 &
     1.95380E-17 4.26810E-17 6.57540E-17 8.87570E-17 1.11690E-16 1.34550E-16 1.57340E-16 1.80050E-16 &
     2.02690E-16 2.25250E-16 2.47740E-16 2.70150E-16 2.92480E-16 3.14730E-16 3.36900E-16 3.58990E-16 &
     3.80990E-16 4.02900E-16 4.24720E-16 4.46440E-16 4.68050E-16 4.89550E-16 5.10940E-16 5.32200E-16 &
     5.53360E-16 5.74410E-16 5.95360E-16 6.16220E-16 6.36990E-16 6.57690E-16 6.78310E-16 6.98830E-16 &
     7.19210E-16 7.39400E-16 7.59370E-16 7.79090E-16 7.98510E-16 8.17610E-16 8.36360E-16 8.54750E-16 &
     8.72780E-16 8.90430E-16 9.07680E-16 9.24540E-16 9.40990E-16 9.57060E-16 9.72810E-16 9.88270E-16 &
     1.00350E-15 1.01850E-15 1.03330E-15 1.04810E-15 1.06260E-15 1.07710E-15 1.09150E-15 1.10570E-15 &
     1.11980E-15 1.13370E-15 1.14750E-15 1.16120E-15 1.17470E-15 1.18790E-15 1.20090E-15 1.21370E-15 &
     1.22610E-15 1.23820E-15 1.25010E-15 1.26160E-15 1.27290E-15 1.28390E-15 1.29460E-15 1.30510E-15 &
     1.31540E-15 1.32540E-15 1.33520E-15 1.34480E-15 1.35410E-15 1.36320E-15 1.37210E-15 1.38070E-15 &
     1.38900E-15 1.39710E-15 1.40500E-15 1.41260E-15 1.42000E-15 1.42710E-15 1.43400E-15 1.44060E-15 &
     1.44700E-15 1.45310E-15 1.45900E-15 1.46470E-15 1.47010E-15 1.47520E-15 1.48010E-15 1.48480E-15 &
     1.48930E-15 1.49360E-15 1.49770E-15 1.50150E-15 1.50520E-15 1.50870E-15 1.51190E-15 1.51480E-15 &
     1.51730E-15 1.51950E-15 1.52120E-15 1.52250E-15 1.52340E-15 1.52380E-15 1.52380E-15 1.52340E-15 &
     1.52260E-15 1.52140E-15 1.51980E-15 1.51790E-15 1.51570E-15 1.51320E-15 1.51050E-15 1.50760E-15 &
     1.50450E-15 1.50130E-15 1.49800E-15 1.49450E-15 1.49080E-15 1.48700E-15 1.48310E-15 1.47890E-15 &
     1.47460E-15 1.47010E-15 1.46540E-15 1.46060E-15 1.45560E-15 1.45040E-15 1.44500E-15 1.43960E-15 &
     1.43390E-15 1.42790E-15 1.42160E-15 1.41500E-15 1.40800E-15 1.40040E-15 1.39240E-15 1.38390E-15 &
     1.37500E-15 1.36570E-15 1.35620E-15 1.34650E-15 1.33670E-15 1.32670E-15 1.31670E-15 1.30650E-15 &
     1.29630E-15 1.28580E-15 1.27530E-15 1.26450E-15 1.25360E-15 1.24250E-15 1.23120E-15 1.21970E-15 &
     1.20810E-15 1.19630E-15 1.18440E-15 1.17230E-15 1.16010E-15 1.14780E-15 1.13530E-15 1.12270E-15 &
     1.11000E-15 1.09710E-15 1.08420E-15 1.07100E-15 1.05770E-15 1.04430E-15 1.03070E-15 1.01700E-15 &
     1.00310E-15 9.89020E-16 9.74780E-16 9.60390E-16 9.45870E-16 9.31220E-16 9.16470E-16 9.01630E-16 &
     8.86700E-16 8.71700E-16 8.56640E-16 8.41520E-16 8.26350E-16 8.11140E-16 7.95900E-16 7.80630E-16 &
     7.65330E-16 7.50010E-16 7.34650E-16 7.19270E-16 7.03870E-16 6.88430E-16 6.72960E-16 6.57460E-16 &
     6.41940E-16 6.26390E-16 6.10810E-16 5.95210E-16 5.79580E-16 5.63940E-16 5.48260E-16 5.32550E-16 &
     5.16790E-16 5.00990E-16 4.85130E-16 4.69200E-16 4.53210E-16 4.37130E-16 4.20970E-16 4.04690E-16 &
     3.88300E-16 3.71780E-16 3.55130E-16 3.38330E-16 3.21410E-16 3.04410E-16 2.87360E-16 2.70300E-16 &
     2.53270E-16 2.36300E-16 2.19410E-16 2.02610E-16 1.85890E-16 1.69240E-16 1.52650E-16 1.36130E-16 &
     1.19650E-16 1.03210E-16 8.68210E-17 7.04640E-17 5.41410E-17 3.78500E-17 2.15870E-17 5.35080E-18 &
     -1.08580E-17 -2.70170E-17 -4.31030E-17 -5.90890E-17 -7.49500E-17 -9.06600E-17 -1.06190E-16 -1.21530E-16 &
     -1.36670E-16 -1.51610E-16 -1.66360E-16 -1.80910E-16 -1.95280E-16 -2.09450E-16 -2.23440E-16 -2.37260E-16 &
     -2.50920E-16 -2.64430E-16 -2.77800E-16 -2.91060E-16 -3.04210E-16 -3.17260E-16 -3.30210E-16 -3.43060E-16 &
     -3.55810E-16 -3.68450E-16 -3.80990E-16 -3.93420E-16 -4.05730E-16 -4.17930E-16 -4.30020E-16 -4.41990E-16 &
     -4.53840E-16 -4.65570E-16 -4.77190E-16 -4.88680E-16 -5.00050E-16 -5.11300E-16 -5.22440E-16 -5.33460E-16 &
     -5.44360E-16 -5.55150E-16 -5.65830E-16 -5.76400E-16 -5.86840E-16 -5.97170E-16 -6.07360E-16 -6.17430E-16 &
     -6.27360E-16 -6.37150E-16 -6.46790E-16 -6.56280E-16 -6.65600E-16 -6.74740E-16 -6.83700E-16 -6.92460E-16 &
     -7.01010E-16 -7.09370E-16 -7.17520E-16 -7.25490E-16 -7.33260E-16 -7.40860E-16 -7.48280E-16 -7.55530E-16 &
     -7.62610E-16 -7.69540E-16 -7.76310E-16 -7.82940E-16 -7.89430E-16 -7.95790E-16 -8.02020E-16 -8.08100E-16 &
     -8.13980E-16 -8.19650E-16 -8.25060E-16 -8.30200E-16 -8.35020E-16 -8.39500E-16 -8.43660E-16 -8.47500E-16 &
     -8.51060E-16 -8.54340E-16 -8.57370E-16 -8.60160E-16 -8.62730E-16 -8.65100E-16 -8.67270E-16 -8.69260E-16 &
     -8.71080E-16 -8.72740E-16 -8.74260E-16 -8.75640E-16 -8.76890E-16 -8.78010E-16 -8.78990E-16 -8.79830E-16 &
     -8.80540E-16 -8.81110E-16 -8.81540E-16 -8.81840E-16 -8.82000E-16 -8.82020E-16 -8.81920E-16 -8.81690E-16 &
     -8.81320E-16 -8.80830E-16 -8.80220E-16 -8.79470E-16 -8.78590E-16 -8.77570E-16 -8.76410E-16 -8.75110E-16 &
     -8.73670E-16 -8.72100E-16 -8.70400E-16 -8.68580E-16 -8.66660E-16 -8.64630E-16 -8.62520E-16 -8.60330E-16 &
     -8.58020E-16 -8.55570E-16 -8.52950E-16 -8.50130E-16 -8.47080E-16 -8.43770E-16 -8.40170E-16 -8.36300E-16 &
     -8.32170E-16 -8.27800E-16 -8.23200E-16 -8.18390E-16 -8.13390E-16 -8.08210E-16 -8.02860E-16 -7.97370E-16 &
     -7.91730E-16 -7.85960E-16 -7.80070E-16 -7.74070E-16 -7.67980E-16 -7.61770E-16 -7.55440E-16 -7.48980E-16 &
     -7.42370E-16 -7.35600E-16 -7.28660E-16 -7.21540E-16 -7.14260E-16 -7.06840E-16 -6.99280E-16 -6.91600E-16 &
     -6.83840E-16 -6.75990E-16 -6.68090E-16 -6.60120E-16 -6.52090E-16 -6.44010E-16 -6.35870E-16 -6.27680E-16 &
     -6.19440E-16 -6.11140E-16 -6.02800E-16 -5.94400E-16 -5.85960E-16 -5.77480E-16 -5.68950E-16 -5.60380E-16 &
     -5.51770E-16 -5.43120E-16 -5.34430E-16 -5.25690E-16 -5.16920E-16 -5.08110E-16 -4.99260E-16 -4.90370E-16 &
     -4.81440E-16 -4.72480E-16 -4.63480E-16 -4.54450E-16 -4.45390E-16 -4.36310E-16 -4.27200E-16 -4.18060E-16 &
     -4.08890E-16 -3.99690E-16 -3.90450E-16 -3.81180E-16 -3.71860E-16 -3.62510E-16 -3.53120E-16 -3.43710E-16 &
     -3.34270E-16 -3.24820E-16 -3.15350E-16 -3.05890E-16 -2.96420E-16 -2.86940E-16 -2.77420E-16 -2.67830E-16 &
     -2.58160E-16 -2.48390E-16 -2.38490E-16 -2.28450E-16 -2.18280E-16 -2.08010E-16 -1.97650E-16 -1.87230E-16 &
     -1.76770E-16 -1.66290E-16 -1.55810E-16 -1.45350E-16 -1.34920E-16 -1.24540E-16 -1.14240E-16 -1.04010E-16 &
     -9.38870E-17 -8.38740E-17 -7.39730E-17 -6.41790E-17 -5.44860E-17 -4.48900E-17 -3.53850E-17 -2.59690E-17 &
     -1.66350E-17 -7.38100E-18 1.79620E-18 1.08990E-17 1.99310E-17 2.88950E-17 3.77930E-17 4.66270E-17 &
     5.53980E-17 6.41030E-17 7.27420E-17 8.13140E-17 8.98170E-17 9.82510E-17 1.06610E-16 1.14910E-16 &
     1.23130E-16 1.31290E-16 1.39390E-16 1.47430E-16 1.55400E-16 1.63320E-16 1.71170E-16 1.78960E-16 &
     1.86690E-16 1.94350E-16 2.01940E-16 2.09460E-16 2.16910E-16 2.24290E-16 2.31590E-16 2.38820E-16 &
     2.45990E-16 2.53090E-16 2.60120E-16 2.67090E-16 2.73990E-16 2.80830E-16 2.87600E-16 2.94290E-16 &
     3.00920E-16 3.07470E-16 3.13940E-16 3.20340E-16 3.26670E-16 3.32930E-16 3.39110E-16 3.45230E-16 &
     3.51270E-16 3.57250E-16 3.63160E-16 3.68990E-16 3.74730E-16 3.80390E-16 3.85960E-16 3.91430E-16 &
     3.96800E-16 4.02050E-16 4.07140E-16 4.12070E-16 4.16800E-16 4.21310E-16 4.25570E-16 4.29590E-16 &
     4.33350E-16 4.36900E-16 4.40240E-16 4.43400E-16 4.46410E-16 4.49280E-16 4.52030E-16 4.54670E-16 &
     4.57200E-16 4.59620E-16 4.61930E-16 4.64140E-16 4.66240E-16 4.68240E-16 4.70140E-16 4.71940E-16 &
     4.73640E-16 4.75240E-16 4.76750E-16 4.78150E-16 4.79460E-16 4.80680E-16 4.81800E-16 4.82820E-16 &
     4.83750E-16 4.84590E-16 4.85330E-16 4.85970E-16 4.86520E-16 4.86980E-16 4.87350E-16 4.87620E-16 &
     4.87810E-16 4.87900E-16 4.87910E-16 4.87830E-16 4.87660E-16 4.87400E-16 4.87060E-16 4.86620E-16 &
     4.86110E-16 4.85500E-16 4.84810E-16 4.84030E-16 4.83160E-16 4.82200E-16 4.81150E-16 4.80000E-16 &
     4.78760E-16 4.77430E-16 4.76000E-16 4.74480E-16 4.72870E-16 4.71160E-16 4.69360E-16 4.67470E-16 &
     4.65490E-16 4.63430E-16 4.61290E-16 4.59070E-16 4.56780E-16 4.54430E-16 4.52000E-16 4.49490E-16 &
     4.46840E-16 4.44040E-16 4.41050E-16 4.37840E-16 4.34380E-16 4.30660E-16 4.26680E-16 4.22500E-16 &
     4.18160E-16 4.13680E-16 4.09110E-16 4.04480E-16 3.99830E-16 3.95150E-16 3.90450E-16 3.85710E-16 &
     3.80930E-16 3.76110E-16 3.71230E-16 3.66300E-16 3.61310E-16 3.56280E-16 3.51190E-16 3.46070E-16 &
     3.40900E-16 3.35700E-16 3.30460E-16 3.25190E-16 3.19890E-16 3.14540E-16 3.09160E-16 3.03740E-16 &
     2.98280E-16 2.92780E-16 2.87230E-16 2.81650E-16 2.76030E-16 2.70370E-16 2.64670E-16 2.58940E-16 &
     2.53180E-16 2.47390E-16 2.41570E-16 2.35710E-16 2.29830E-16 2.23910E-16 2.17970E-16 2.11990E-16 &
     2.05990E-16 1.99960E-16 1.93900E-16 1.87820E-16 1.81700E-16 1.75560E-16 1.69390E-16 1.63200E-16 &
     1.56970E-16 1.50730E-16 1.44460E-16 1.38160E-16 1.31850E-16 1.25520E-16 1.19160E-16 1.12780E-16 &
     1.06380E-16 9.99450E-17 9.34800E-17 8.69820E-17 8.04510E-17 7.38950E-17 6.73260E-17 6.07530E-17 &
     5.41890E-17 4.76430E-17 4.11270E-17 3.46490E-17 2.82030E-17 2.17850E-17 1.53850E-17 8.99790E-18 &
     2.61600E-18 -3.76770E-18 -1.01580E-17 -1.65490E-17 -2.29340E-17 -2.93050E-17 -3.56540E-17 -4.19730E-17 &
     -4.82560E-17 -5.44950E-17 -6.06910E-17 -6.68440E-17 -7.29570E-17 -7.90300E-17 -8.50640E-17 -9.10610E-17 &
     -9.70220E-17 -1.02950E-16 -1.08830E-16 -1.14680E-16 -1.20500E-16 -1.26270E-16 -1.32010E-16 -1.37710E-16 &
     -1.43370E-16 -1.48990E-16 -1.54570E-16 -1.60110E-16 -1.65620E-16 -1.71080E-16 -1.76510E-16 -1.81890E-16 &
     -1.87230E-16 -1.92540E-16 -1.97800E-16 -2.03030E-16 -2.08210E-16 -2.13350E-16 -2.18460E-16 -2.23520E-16 &
     -2.28540E-16 -2.33520E-16 -2.38450E-16 -2.43350E-16 -2.48200E-16 -2.53010E-16 -2.57770E-16 -2.62480E-16 &
     -2.67130E-16 -2.71740E-16 -2.76280E-16 -2.80770E-16 -2.85200E-16 -2.89590E-16 -2.93950E-16 -2.98270E-16 &
     -3.02560E-16 -3.06850E-16 -3.11110E-16 -3.15340E-16 -3.19510E-16 -3.23600E-16 -3.27570E-16 -3.31400E-16 &
     -3.35060E-16 -3.38550E-16 -3.41860E-16 -3.45020E-16 -3.48040E-16 -3.50940E-16 -3.53750E-16 -3.56480E-16 &
     -3.59140E-16 -3.61730E-16 -3.64250E-16 -3.66690E-16 -3.69040E-16 -3.71310E-16 -3.73480E-16 -3.75550E-16 &
     -3.77530E-16 -3.79420E-16 -3.81230E-16 -3.82950E-16 -3.84610E-16 -3.86200E-16 -3.87730E-16 -3.89200E-16 &
     -3.90610E-16 -3.91960E-16 -3.93240E-16 -3.94460E-16 -3.95610E-16 -3.96700E-16 -3.97730E-16 -3.98680E-16 &
     -3.99580E-16 -4.00400E-16 -4.01160E-16 -4.01850E-16 -4.02470E-16 -4.03020E-16 -4.03510E-16 -4.03930E-16 &
     -4.04290E-16 -4.04590E-16 -4.04820E-16 -4.05000E-16 -4.05110E-16 -4.05170E-16 -4.05160E-16 -4.05100E-16 &
     -4.04970E-16 -4.04790E-16 -4.04540E-16 -4.04240E-16 -4.03870E-16 -4.03440E-16 -4.02940E-16 -4.02380E-16 &
     -4.01750E-16 -4.01050E-16 -4.00290E-16 -3.99470E-16 -3.98600E-16 -3.97680E-16 -3.96730E-16 -3.95730E-16 &
     -3.94700E-16 -3.93620E-16 -3.92460E-16 -3.91220E-16 -3.89870E-16 -3.88390E-16 -3.86760E-16 -3.84970E-16 &
     -3.83040E-16 -3.80990E-16 -3.78850E-16 -3.76630E-16 -3.74360E-16 -3.72080E-16 -3.69790E-16 -3.67480E-16 &
     -3.65140E-16 -3.62750E-16 -3.60290E-16 -3.57740E-16 -3.55100E-16 -3.52340E-16 -3.49470E-16 -3.46500E-16 &
     -3.43440E-16 -3.40300E-16 -3.37090E-16 -3.33820E-16 -3.30500E-16 -3.27130E-16 -3.23730E-16 -3.20280E-16 &
     -3.16790E-16 -3.13270E-16 -3.09710E-16 -3.06130E-16 -3.02520E-16 -2.98880E-16 -2.95210E-16 -2.91510E-16 &
     -2.87780E-16 -2.84020E-16 -2.80230E-16 -2.76410E-16 -2.72560E-16 -2.68680E-16 -2.64780E-16 -2.60840E-16 &
     -2.56890E-16 -2.52900E-16 -2.48890E-16 -2.44860E-16 -2.40800E-16 -2.36720E-16 -2.32610E-16 -2.28470E-16 &
     -2.24300E-16 -2.20110E-16 -2.15900E-16 -2.11670E-16 -2.07410E-16 -2.03140E-16 -1.98850E-16 -1.94540E-16 &
     -1.90210E-16 -1.85870E-16 -1.81510E-16 -1.77130E-16 -1.72730E-16 -1.68300E-16 -1.63850E-16 -1.59380E-16 &
     -1.54890E-16 -1.50390E-16 -1.45860E-16 -1.41330E-16 -1.36790E-16 -1.32240E-16 -1.27690E-16 -1.23130E-16 &
     -1.18570E-16 -1.14010E-16 -1.09460E-16 -1.04910E-16 -1.00360E-16 -9.58180E-17 -9.12840E-17 -8.67590E-17 &
     -8.22430E-17 -7.77370E-17 -7.32430E-17 -6.87590E-17 -6.42820E-17 -5.98060E-17 -5.53260E-17 -5.08360E-17 &
     -4.63320E-17 -4.18070E-17 -3.72590E-17 -3.26930E-17 -2.81150E-17 -2.35330E-17 -1.89530E-17 -1.43840E-17 &
     -9.83180E-18 -5.30270E-18 -7.98040E-19 3.68190E-18 8.13690E-18 1.25670E-17 1.69710E-17 2.13500E-17 &
     2.57020E-17 3.00290E-17 3.43300E-17 3.86040E-17 4.28530E-17 4.70750E-17 5.12700E-17 5.54400E-17 &
     5.95830E-17 6.37010E-17 6.77930E-17 7.18610E-17 7.59050E-17 7.99250E-17 8.39200E-17 8.78910E-17 &
     9.18340E-17 9.57490E-17 9.96330E-17 1.03490E-16 1.07300E-16 1.11090E-16 1.14840E-16 1.18560E-16 &
     1.22270E-16 1.25960E-16 1.29630E-16 1.33300E-16 1.36950E-16 1.40590E-16 1.44190E-16 1.47730E-16 &
     1.51200E-16 1.54580E-16 1.57850E-16 1.61000E-16 1.64040E-16 1.66980E-16 1.69850E-16 1.72660E-16 &
     1.75420E-16 1.78150E-16 1.80870E-16 1.83560E-16 1.86240E-16 1.88890E-16 1.91500E-16 1.94080E-16 &
     1.96610E-16 1.99090E-16 2.01530E-16 2.03930E-16 2.06290E-16 2.08610E-16 2.10910E-16 2.13170E-16 &
     2.15400E-16 2.17600E-16 2.19760E-16 2.21860E-16 2.23910E-16 2.25880E-16 2.27770E-16 2.29570E-16 &
     2.31290E-16 2.32940E-16 2.34520E-16 2.36040E-16 2.37510E-16 2.38940E-16 2.40320E-16 2.41670E-16 &
     2.42980E-16 2.44240E-16 2.45460E-16 2.46640E-16 2.47770E-16 2.48850E-16 2.49890E-16 2.50880E-16 &
     2.51820E-16 2.52710E-16 2.53560E-16 2.54370E-16 2.55130E-16 2.55850E-16 2.56520E-16 2.57140E-16 &
     2.57720E-16 2.58250E-16 2.58740E-16 2.59170E-16 2.59570E-16 2.59920E-16 2.60230E-16 2.60510E-16 &
     2.60760E-16 2.60980E-16 2.61170E-16 2.61330E-16 2.61440E-16 2.61470E-16 2.61430E-16 2.61280E-16 &
     2.61030E-16 2.60650E-16 2.60170E-16 2.59590E-16 2.58930E-16 2.58220E-16 2.57460E-16 2.56670E-16 &
     2.55870E-16 2.55060E-16 2.54220E-16 2.53370E-16 2.52480E-16 2.51570E-16 2.50620E-16 2.49630E-16 &
     2.48610E-16 2.47550E-16 2.46450E-16 2.45320E-16 2.44160E-16 2.42970E-16 2.41750E-16 2.40500E-16 &
     2.39210E-16 2.37900E-16 2.36560E-16 2.35190E-16 2.33790E-16 2.32360E-16 2.30890E-16 2.29380E-16 &
     2.27810E-16 2.26190E-16 2.24490E-16 2.22720E-16 2.20860E-16 2.18930E-16 2.16940E-16 2.14890E-16 &
     2.12810E-16 2.10690E-16 2.08560E-16 2.06420E-16 2.04270E-16 2.02100E-16 1.99920E-16 1.97720E-16 &
     1.95500E-16 1.93250E-16 1.90980E-16 1.88680E-16 1.86350E-16 1.84000E-16 1.81620E-16 1.79230E-16 &
     1.76810E-16 1.74360E-16 1.71900E-16 1.69420E-16 1.66920E-16 1.64400E-16 1.61870E-16 1.59310E-16 &
     1.56740E-16 1.54150E-16 1.51540E-16 1.48900E-16 1.46250E-16 1.43560E-16 1.40850E-16 1.38110E-16 &
     1.35350E-16 1.32570E-16 1.29770E-16 1.26970E-16 1.24160E-16 1.21340E-16 1.18540E-16 1.15730E-16 &
     1.12930E-16 1.10120E-16 1.07320E-16 1.04510E-16 1.01700E-16 9.88860E-17 9.60660E-17 9.32430E-17 &
     9.04190E-17 8.75920E-17 8.47660E-17 8.19400E-17 7.91160E-17 7.62930E-17 7.34730E-17 7.06540E-17 &
     6.78380E-17 6.50240E-17 6.22120E-17 5.94030E-17 5.65950E-17 5.37870E-17 5.09770E-17 4.81650E-17 &
     4.53490E-17 4.25280E-17 3.97000E-17 3.68670E-17 3.40280E-17 3.11840E-17 2.83370E-17 2.54860E-17 &
     2.26320E-17 1.97760E-17 1.69210E-17 1.40690E-17 1.12230E-17 8.38420E-18 5.55680E-18 2.74280E-18 &
     -5.55920E-20 -2.83780E-18 -5.60370E-18 -8.35330E-18 -1.10860E-17 -1.38030E-17 -1.65030E-17 -1.91860E-17 &
     -2.18530E-17 -2.45070E-17 -2.71470E-17 -2.97750E-17 -3.23930E-17 -3.50010E-17 -3.76000E-17 -4.01840E-17 &
     -4.27450E-17 -4.52760E-17 -4.77690E-17 -5.02180E-17 -5.26160E-17 -5.49570E-17 -5.72440E-17 -5.94840E-17 &
     -6.16830E-17 -6.38470E-17 -6.59820E-17 -6.80930E-17 -7.01850E-17 -7.22590E-17 -7.43140E-17 -7.63470E-17 &
     -7.83590E-17 -8.03480E-17 -8.23130E-17 -8.42530E-17 -8.61680E-17 -8.80580E-17 -8.99230E-17 -9.17630E-17 &
     -9.35790E-17 -9.53710E-17 -9.71380E-17 -9.88810E-17 -1.00600E-16 -1.02290E-16 -1.03970E-16 -1.05610E-16 &
     -1.07240E-16 -1.08840E-16 -1.10420E-16 -1.11970E-16 -1.13500E-16 -1.15010E-16 -1.16490E-16 -1.17960E-16 &
     -1.19400E-16 -1.20810E-16 -1.22200E-16 -1.23560E-16 -1.24890E-16 -1.26190E-16 -1.27460E-16 -1.28690E-16 &
     -1.29880E-16 -1.31040E-16 -1.32160E-16 -1.33240E-16 -1.34280E-16 -1.35290E-16 -1.36250E-16 -1.37170E-16 &
     -1.38060E-16 -1.38910E-16 -1.39720E-16 -1.40500E-16 -1.41250E-16 -1.41960E-16 -1.42650E-16 -1.43300E-16 &
     -1.43930E-16 -1.44530E-16 -1.45100E-16 -1.45640E-16 -1.46150E-16 -1.46630E-16 -1.47070E-16 -1.47470E-16 &
     -1.47810E-16 -1.48090E-16 -1.48310E-16 -1.48450E-16 -1.48530E-16 -1.48560E-16 -1.48530E-16 -1.48470E-16 &
     -1.48380E-16 -1.48270E-16 -1.48140E-16 -1.48000E-16 -1.47840E-16 -1.47660E-16 -1.47460E-16 -1.47230E-16 &
     -1.46980E-16 -1.46700E-16 -1.46390E-16 -1.46050E-16 -1.45690E-16 -1.45300E-16 -1.44890E-16 -1.44460E-16 &
     -1.44000E-16 -1.43520E-16 -1.43020E-16 -1.42490E-16 -1.41950E-16 -1.41380E-16 -1.40790E-16 -1.40190E-16 &
     -1.39560E-16 -1.38910E-16 -1.38230E-16 -1.37540E-16 -1.36810E-16 -1.36070E-16 -1.35300E-16 -1.34500E-16 &
     -1.33680E-16 -1.32850E-16 -1.31990E-16 -1.31110E-16 -1.30220E-16 -1.29320E-16 -1.28390E-16 -1.27440E-16 &
     -1.26460E-16 -1.25450E-16 -1.24400E-16 -1.23310E-16 -1.22170E-16 -1.20980E-16 -1.19760E-16 -1.18490E-16 &
     -1.17200E-16 -1.15870E-16 -1.14520E-16 -1.13150E-16 -1.11750E-16 -1.10340E-16 -1.08900E-16 -1.07450E-16 &
     -1.05980E-16 -1.04490E-16 -1.02990E-16 -1.01470E-16 -9.99310E-17 -9.83800E-17 -9.68140E-17 -9.52330E-17 &
     -9.36380E-17 -9.20300E-17 -9.04080E-17 -8.87730E-17 -8.71250E-17 -8.54650E-17 -8.37930E-17 -8.21090E-17 &
     -8.04140E-17 -7.87080E-17 -7.69920E-17 -7.52650E-17 -7.35290E-17 -7.17830E-17 -7.00280E-17 -6.82650E-17 &
     -6.64930E-17 -6.47130E-17 -6.29250E-17 -6.11310E-17 -5.93290E-17 -5.75200E-17 -5.57060E-17 -5.38860E-17 &
     -5.20600E-17 -5.02290E-17 -4.83940E-17 -4.65540E-17 -4.47100E-17 -4.28620E-17 -4.10120E-17 -3.91580E-17 &
     -3.73020E-17 -3.54430E-17 -3.35830E-17 -3.17210E-17 -2.98580E-17 -2.79950E-17 -2.61310E-17 -2.42670E-17 &
     -2.24030E-17 -2.05400E-17 -1.86780E-17 -1.68180E-17 -1.49600E-17 -1.31030E-17 -1.12490E-17 -9.39790E-18 &
     -7.54990E-18 -5.70530E-18 -3.86440E-18 -2.02770E-18 -1.95440E-19 1.63210E-18 3.45450E-18 5.27150E-18 &
     7.08270E-18 8.88780E-18 1.06870E-17 1.24790E-17 1.42630E-17 1.60410E-17 1.78110E-17 1.95720E-17 &
     2.13260E-17 2.30700E-17 2.48060E-17 2.65320E-17 2.82480E-17 2.99550E-17 3.16510E-17 3.33360E-17 &
     3.50110E-17 3.66740E-17 3.83250E-17 3.99640E-17 4.15910E-17 4.32050E-17 4.48060E-17 4.63940E-17 &
     4.79690E-17 4.95290E-17 5.10750E-17 5.26060E-17 5.41220E-17 5.56230E-17 5.71090E-17 5.85780E-17 &
     6.00310E-17 6.14680E-17 6.28870E-17 6.42890E-17 6.56740E-17 6.70410E-17 6.83890E-17 6.97190E-17 &
     7.10300E-17 7.23210E-17 7.35930E-17 7.48450E-17 7.60770E-17 7.72880E-17 7.84790E-17 7.96480E-17 &
     8.07950E-17 8.19210E-17 8.30240E-17 8.41050E-17 8.51630E-17 8.61980E-17 8.72090E-17 8.81960E-17 &
     8.91600E-17 9.00980E-17 9.10120E-17 9.19000E-17 9.27630E-17 9.36010E-17 9.44120E-17 9.51960E-17 &
     9.59540E-17 9.66870E-17 9.73930E-17 9.80740E-17 9.87290E-17 9.93580E-17 9.99630E-17 1.00540E-16 &
     1.01100E-16 1.01630E-16 1.02140E-16 1.02620E-16 1.03080E-16 1.03510E-16 1.03920E-16 1.04310E-16 &
     1.04680E-16 1.05020E-16 1.05340E-16 1.05640E-16 1.05910E-16 1.06170E-16 1.06400E-16 1.06610E-16 &
     1.06800E-16 1.06960E-16 1.07110E-16 1.07240E-16 1.07340E-16 1.07430E-16 1.07490E-16 1.07530E-16 &
     1.07560E-16 1.07570E-16 1.07550E-16 1.07520E-16 1.07470E-16 1.07400E-16 1.07310E-16 1.07200E-16 &
     1.07070E-16 1.06930E-16 1.06770E-16 1.06590E-16 1.06400E-16 1.06180E-16 1.05950E-16 1.05710E-16 &
     1.05440E-16 1.05160E-16 1.04870E-16 1.04560E-16 1.04230E-16 1.03890E-16 1.03530E-16 1.03160E-16 &
     1.02770E-16 1.02370E-16 1.01960E-16 1.01530E-16 1.01080E-16 1.00620E-16 1.00150E-16 9.96680E-17 &
     9.91700E-17 9.86590E-17 9.81350E-17 9.75990E-17 9.70510E-17 9.64900E-17 9.59170E-17 9.53330E-17 &
     9.47370E-17 9.41300E-17 9.35110E-17 9.28820E-17 9.22420E-17 9.15920E-17 9.09310E-17 9.02600E-17 &
     8.95800E-17 8.88890E-17 8.81900E-17 8.74810E-17 8.67630E-17 8.60360E-17 8.53000E-17 8.45560E-17 &
     8.38040E-17 8.30430E-17 8.22750E-17 8.15000E-17 8.07160E-17 7.99260E-17 7.91290E-17 7.83250E-17 &
     7.75140E-17 7.66970E-17 7.58730E-17 7.50440E-17 7.42090E-17 7.33690E-17 7.25230E-17 7.16720E-17 &
     7.08160E-17 6.99550E-17 6.90890E-17 6.82200E-17 6.73460E-17 6.64680E-17 6.55870E-17 6.47020E-17 &
     6.38140E-17 6.29220E-17 6.20280E-17 6.11300E-17 6.02290E-17 5.93260E-17 5.84190E-17 5.75100E-17 &
     5.65980E-17 5.56840E-17 5.47670E-17 5.38480E-17 5.29260E-17 5.20020E-17 5.10760E-17 5.01480E-17 &
     4.92180E-17 4.82860E-17 4.73520E-17 4.64160E-17 4.54790E-17 4.45400E-17 4.35990E-17 4.26570E-17 &
     4.17140E-17 4.07690E-17 3.98230E-17 3.88760E-17 3.79280E-17 3.69790E-17 3.60290E-17 3.50790E-17 &
     3.41270E-17 3.31750E-17 3.22230E-17 3.12700E-17 3.03160E-17 2.93620E-17 2.84080E-17 2.74540E-17 &
     2.65000E-17 2.55450E-17 2.45910E-17 2.36370E-17 2.26830E-17 2.17300E-17 2.07770E-17 1.98240E-17 &
     1.88720E-17 1.79210E-17 1.69700E-17 1.60200E-17 1.50710E-17 1.41230E-17 1.31760E-17 1.22300E-17 &
     1.12860E-17 1.03420E-17 9.40040E-18 8.45990E-18 7.52090E-18 6.58360E-18 5.64790E-18 4.71400E-18 &
     3.78200E-18 2.85190E-18 1.92380E-18 9.97870E-19 7.40730E-20 -8.47460E-19 -1.76670E-18 -2.68340E-18 &
     -3.59770E-18 -4.50940E-18 -5.41840E-18 -6.32470E-18 -7.22820E-18 -8.12880E-18 -9.02630E-18 -9.92090E-18 &
     -1.08120E-17 -1.17000E-17 -1.25850E-17 -1.34670E-17 -1.43450E-17 -1.52190E-17 -1.60900E-17 -1.69570E-17 &
     -1.78200E-17 -1.86800E-17 -1.95350E-17 -2.03870E-17 -2.12340E-17 -2.20770E-17 -2.29160E-17 -2.37500E-17 &
     -2.45810E-17 -2.54060E-17 -2.62270E-17 -2.70430E-17 -2.78550E-17 -2.86620E-17 -2.94630E-17 -3.02600E-17 &
     -3.10520E-17 -3.18380E-17 -3.26200E-17 -3.33960E-17 -3.41660E-17 -3.49310E-17 -3.56910E-17 -3.64450E-17 &
     -3.71930E-17 -3.79350E-17 -3.86710E-17 -3.94020E-17 -4.01260E-17 -4.08450E-17 -4.15570E-17 -4.22630E-17 &
     -4.29620E-17 -4.36550E-17 -4.43410E-17 -4.50210E-17 -4.56940E-17 -4.63610E-17 -4.70200E-17 -4.76730E-17 &
     -4.83180E-17 -4.89570E-17 -4.95880E-17 -5.02120E-17 -5.08290E-17 -5.14380E-17 -5.20400E-17 -5.26340E-17 &
     -5.32210E-17 -5.38000E-17 -5.43710E-17 -5.49340E-17 -5.54890E-17 -5.60360E-17 -5.65750E-17 -5.71050E-17 &
     -5.76280E-17 -5.81420E-17 -5.86470E-17 -5.91440E-17 -5.96320E-17 -6.01120E-17 -6.05830E-17 -6.10450E-17 &
     -6.14980E-17 -6.19420E-17 -6.23770E-17 -6.28020E-17 -6.32190E-17 -6.36260E-17 -6.40230E-17 -6.44110E-17 &
     -6.47900E-17 -6.51590E-17 -6.55180E-17 -6.58670E-17 -6.62060E-17 -6.65360E-17 -6.68550E-17 -6.71640E-17 &
     -6.74630E-17 -6.77510E-17 -6.80300E-17 -6.82970E-17 -6.85550E-17 -6.88020E-17 -6.90390E-17 -6.92660E-17 &
     -6.94830E-17 -6.96900E-17 -6.98870E-17 -7.00740E-17 -7.02520E-17 -7.04190E-17 -7.05770E-17 -7.07260E-17 &
     -7.08640E-17 -7.09930E-17 -7.11130E-17 -7.12230E-17 -7.13240E-17 -7.14160E-17 -7.14990E-17 -7.15720E-17 &
     -7.16360E-17 -7.16910E-17 -7.17370E-17 -7.17750E-17 -7.18030E-17 -7.18230E-17 -7.18330E-17 -7.18350E-17 &
     -7.18290E-17 -7.18140E-17 -7.17900E-17 -7.17580E-17 -7.17180E-17 -7.16690E-17 -7.16120E-17 -7.15460E-17 &
     -7.14730E-17 -7.13910E-17 -7.13020E-17 -7.12040E-17 -7.10990E-17 -7.09850E-17 -7.08640E-17 -7.07350E-17 &
     -7.05980E-17 -7.04540E-17 -7.03020E-17 -7.01430E-17 -6.99760E-17 -6.98020E-17 -6.96200E-17 -6.94310E-17 &
     -6.92350E-17 -6.90320E-17 -6.88210E-17 -6.86040E-17 -6.83800E-17 -6.81480E-17 -6.79100E-17 -6.76650E-17 &
     -6.74130E-17 -6.71550E-17 -6.68900E-17 -6.66180E-17 -6.63400E-17 -6.60550E-17 -6.57640E-17 -6.54670E-17 &
     -6.51630E-17 -6.48530E-17 -6.45370E-17 -6.42150E-17 -6.38870E-17 -6.35530E-17 -6.32130E-17 -6.28670E-17 &
     -6.25150E-17 -6.21580E-17 -6.17950E-17 -6.14260E-17 -6.10510E-17 -6.06710E-17 -6.02860E-17 -5.98950E-17 &
     -5.94990E-17 -5.90980E-17 -5.86910E-17 -5.82790E-17 -5.78630E-17 -5.74410E-17 -5.70140E-17 -5.65820E-17 &
     -5.61450E-17 -5.57040E-17 -5.52580E-17 -5.48070E-17 -5.43510E-17 -5.38910E-17 -5.34260E-17 -5.29570E-17 &
     -5.24840E-17 -5.20060E-17 -5.15240E-17 -5.10380E-17 -5.05470E-17 -5.00530E-17 -4.95540E-17 -4.90510E-17 &
     -4.85450E-17 -4.80340E-17 -4.75200E-17 -4.70020E-17 -4.64810E-17 -4.59550E-17 -4.54270E-17 -4.48940E-17 &
     -4.43580E-17 -4.38190E-17 -4.32770E-17 -4.27310E-17 -4.21820E-17 -4.16300E-17 -4.10740E-17 -4.05160E-17 &
     -3.99550E-17 -3.93910E-17 -3.88240E-17 -3.82540E-17 -3.76810E-17 -3.71060E-17 -3.65280E-17 -3.59470E-17 &
     -3.53640E-17 -3.47790E-17 -3.41910E-17 -3.36010E-17 -3.30080E-17 -3.24130E-17 -3.18170E-17 -3.12180E-17 &
     -3.06170E-17 -3.00140E-17 -2.94090E-17 -2.88020E-17 -2.81930E-17 -2.75830E-17 -2.69710E-17 -2.63580E-17 &
     -2.57420E-17 -2.51260E-17 -2.45080E-17 -2.38880E-17 -2.32670E-17 -2.26450E-17 -2.20220E-17 -2.13970E-17 &
     -2.07710E-17 -2.01450E-17 -1.95170E-17 -1.88890E-17 -1.82590E-17 -1.76290E-17 -1.69980E-17 -1.63660E-17 &
     -1.57340E-17 -1.51010E-17 -1.44680E-17 -1.38340E-17 -1.32000E-17 -1.25650E-17 -1.19310E-17 -1.12960E-17 &
     -1.06600E-17 -1.00250E-17 -9.38980E-18 -8.75450E-18 -8.11930E-18 -7.48430E-18 -6.84940E-18 -6.21480E-18 &
     -5.58050E-18 -4.94660E-18 -4.31320E-18 -3.68020E-18 -3.04770E-18 -2.41590E-18 -1.78480E-18 -1.15430E-18 &
     -5.24690E-19 1.04100E-19 7.31990E-19 1.35890E-18 1.98480E-18 2.60960E-18 3.23330E-18 3.85580E-18 &
     4.47700E-18 5.09690E-18 5.71550E-18 6.33250E-18 6.94810E-18 7.56220E-18 8.17460E-18 8.78540E-18 &
     9.39440E-18 1.00020E-17 1.06070E-17 1.12110E-17 1.18120E-17 1.24120E-17 1.30090E-17 1.36040E-17 &
     1.41970E-17 1.47880E-17 1.53770E-17 1.59630E-17 1.65470E-17 1.71280E-17 1.77060E-17 1.82820E-17 &
     1.88560E-17 1.94260E-17 1.99940E-17 2.05590E-17 2.11210E-17 2.16800E-17 2.22370E-17 2.27900E-17 &
     2.33390E-17 2.38860E-17 2.44300E-17 2.49700E-17 2.55060E-17 2.60400E-17 2.65690E-17 2.70960E-17 &
     2.76180E-17 2.81370E-17 2.86520E-17 2.91640E-17 2.96710E-17 3.01750E-17 3.06750E-17 3.11700E-17 &
     3.16620E-17 3.21490E-17 3.26330E-17 3.31120E-17 3.35860E-17 3.40560E-17 3.45220E-17 3.49840E-17 &
     3.54400E-17 3.58930E-17 3.63400E-17 3.67830E-17 3.72220E-17 3.76550E-17 3.80850E-17 3.85100E-17 &
     3.89300E-17 3.93460E-17 3.97570E-17 4.01640E-17 4.05660E-17 4.09640E-17 4.13570E-17 4.17460E-17 &
     4.21310E-17 4.25110E-17 4.28870E-17 4.32580E-17 4.36250E-17 4.39880E-17 4.43460E-17 4.47000E-17 &
     4.50500E-17 4.53950E-17 4.57370E-17 4.60740E-17 4.64060E-17 4.67350E-17 4.70590E-17 4.73790E-17 &
     4.76950E-17 4.80060E-17 4.83140E-17 4.86170E-17 4.89160E-17 4.92110E-17 4.95020E-17 4.97890E-17 &
     5.00720E-17 5.03500E-17 5.06250E-17 5.08950E-17 5.11620E-17 5.14240E-17 5.16830E-17 5.19370E-17 &
     5.21880E-17 5.24340E-17 5.26770E-17 5.29150E-17 5.31500E-17 5.33810E-17 5.36070E-17 5.38300E-17 &
     5.40500E-17 5.42650E-17 5.44760E-17 5.46840E-17 5.48870E-17 5.50870E-17 5.52830E-17 5.54760E-17 &
     5.56640E-17 5.58490E-17 5.60300E-17 5.62080E-17 5.63810E-17 5.65510E-17 5.67170E-17 5.68800E-17 &
     5.70390E-17 5.71940E-17 5.73460E-17 5.74940E-17 5.76380E-17 5.77790E-17 5.79160E-17 5.80500E-17 &
     5.81800E-17 5.83070E-17 5.84300E-17 5.85500E-17 5.86660E-17 5.87780E-17 5.88870E-17 5.89930E-17 &
     5.90950E-17 5.91940E-17 5.92890E-17 5.93810E-17 5.94700E-17 5.95550E-17 5.96370E-17 5.97150E-17 &
     5.97900E-17 5.98620E-17 5.99310E-17 5.99960E-17 6.00580E-17 6.01160E-17 6.01720E-17 6.02240E-17 &
     6.02730E-17 6.03190E-17 6.03610E-17 6.04010E-17 6.04370E-17 6.04700E-17 6.04990E-17 6.05260E-17 &
     6.05500E-17 6.05700E-17 6.05880E-17 6.06020E-17 6.06130E-17 6.06210E-17 6.06260E-17 6.06290E-17 &
     6.06280E-17 6.06240E-17 6.06170E-17 6.06070E-17 6.05940E-17 6.05780E-17 6.05600E-17 6.05380E-17 &
     6.05140E-17 6.04860E-17 6.04560E-17 6.04230E-17 6.03870E-17 6.03480E-17 6.03060E-17 6.02620E-17 &
     6.02140E-17 6.01640E-17 6.01110E-17 6.00560E-17 5.99970E-17 5.99360E-17 5.98730E-17 5.98060E-17 &
     5.97370E-17 5.96650E-17 5.95900E-17 5.95130E-17 5.94330E-17 5.93510E-17 5.92660E-17 5.91780E-17 &
     5.90880E-17 5.89950E-17 5.89000E-17 5.88020E-17 5.87020E-17 5.85990E-17 5.84930E-17 5.83850E-17 &
     5.82750E-17 5.81620E-17 5.80470E-17 5.79290E-17 5.78090E-17 5.76860E-17 5.75610E-17 5.74340E-17 &
     5.73040E-17 5.71720E-17 5.70370E-17 5.69010E-17 5.67620E-17 5.66200E-17 5.64770E-17 5.63310E-17 &
     5.61820E-17 5.60320E-17 5.58790E-17 5.57240E-17 5.55670E-17 5.54080E-17 5.52470E-17 5.50830E-17 &
     5.49170E-17 5.47490E-17 5.45790E-17 5.44070E-17 5.42330E-17 5.40560E-17 5.38780E-17 5.36970E-17 &
     5.35150E-17 5.33300E-17 5.31440E-17 5.29550E-17 5.27650E-17 5.25720E-17 5.23780E-17 5.21810E-17 &
     5.19830E-17 5.17820E-17 5.15800E-17 5.13760E-17 5.11700E-17 5.09620E-17 5.07520E-17 5.05410E-17 &
     5.03270E-17 5.01120E-17 4.98950E-17 4.96760E-17 4.94560E-17 4.92330E-17 4.90090E-17 4.87830E-17 &
     4.85560E-17 4.83260E-17 4.80950E-17 4.78630E-17 4.76280E-17 4.73930E-17 4.71550E-17 4.69160E-17 &
     4.66750E-17 4.64320E-17 4.61880E-17 4.59430E-17 4.56950E-17 4.54470E-17 4.51960E-17 4.49450E-17 &
     4.46910E-17 4.44370E-17 4.41800E-17 4.39230E-17 4.36630E-17 4.34030E-17 4.31410E-17 4.28770E-17 &
     4.26120E-17 4.23460E-17 4.20790E-17 4.18100E-17 4.15390E-17 4.12680E-17 4.09950E-17 4.07200E-17 &
     4.04450E-17 4.01680E-17 3.98900E-17 3.96100E-17 3.93300E-17 3.90480E-17 3.87650E-17 3.84810E-17 &
     3.81950E-17 3.79090E-17 3.76210E-17 3.73320E-17 3.70420E-17 3.67510E-17 3.64590E-17 3.61650E-17 &
     3.58710E-17 3.55750E-17 3.52790E-17 3.49810E-17 3.46830E-17 3.43830E-17 3.40830E-17 3.37810E-17 &
     3.34790E-17 3.31750E-17 3.28710E-17 3.25650E-17 3.22590E-17 3.19520E-17 3.16440E-17 3.13350E-17 &
     3.10250E-17 3.07140E-17 3.04030E-17 3.00910E-17 2.97770E-17 2.94640E-17 2.91490E-17 2.88330E-17 &
     2.85170E-17 2.82000E-17 2.78830E-17 2.75640E-17 2.72450E-17 2.69250E-17 2.66050E-17 2.62840E-17 &
     2.59620E-17 2.56400E-17 2.53170E-17 2.49930E-17 2.46690E-17 2.43440E-17 2.40190E-17 2.36930E-17 &
     2.33670E-17 2.30400E-17 2.27120E-17 2.23840E-17 2.20560E-17 2.17270E-17 2.13970E-17 2.10680E-17 &
     2.07370E-17 2.04070E-17 2.00760E-17 1.97440E-17 1.94120E-17 1.90800E-17 1.87480E-17 1.84150E-17 &
     1.80820E-17 1.77480E-17 1.74140E-17 1.70800E-17 1.67460E-17 1.64110E-17 1.60760E-17 1.57410E-17 &
     1.54060E-17 1.50700E-17 1.47350E-17 1.43990E-17 1.40630E-17 1.37270E-17 1.33900E-17 1.30540E-17 &
     1.27170E-17 1.23810E-17 1.20440E-17 1.17070E-17 1.13700E-17 1.10330E-17 1.06970E-17 1.03600E-17 &
             /) (/5999,8/))

         end if

        end subroutine CheckLoadedHighLTemplate



        subroutine Init_Cls

        call CheckLoadedHighLTemplate
        if (CP%WantScalars) then
         if (allocated(Cl_scalar)) deallocate(Cl_scalar)
         allocate(Cl_scalar(lmin:CP%Max_l, CP%InitPower%nn, C_Temp:C_last))
         Cl_scalar = 0
        end if

        if (CP%WantVectors) then
         if (allocated(Cl_vector)) deallocate(Cl_vector)
         allocate(Cl_vector(lmin:CP%Max_l, CP%InitPower%nn, CT_Temp:CT_Cross))
         Cl_vector = 0
        end if


        if (CP%WantTensors) then
          if (allocated(Cl_tensor)) deallocate(Cl_tensor)
          allocate(Cl_tensor(lmin:CP%Max_l_tensor, CP%InitPower%nn, CT_Temp:CT_Cross))
          Cl_tensor = 0
        end if

        end subroutine Init_Cls

        subroutine output_cl_files(ScalFile,TensFile, TotFile, LensFile, LensTotFile, factor)
        implicit none
        integer in,il
        character(LEN=*) ScalFile, TensFile, TotFile, LensFile, LensTotFile
        real(dl), intent(in), optional :: factor
        real(dl) fact
        integer last_C

        if (present(factor)) then
          fact = factor
        else
          fact =1
        end if

         if (CP%WantScalars .and. ScalFile /= '') then
           last_C=min(C_PhiTemp,C_last)
           open(unit=fileio_unit,file=ScalFile,form='formatted',status='replace')
           do in=1,CP%InitPower%nn
             do il=lmin,min(10000,CP%Max_l)
               write(fileio_unit,trim(numcat('(1I6,',last_C))//'E15.5)')il ,fact*Cl_scalar(il,in,C_Temp:last_C)
             end do
             do il=10100,CP%Max_l, 100
               write(fileio_unit,trim(numcat('(1E15.5,',last_C))//'E15.5)') real(il),&
                       fact*Cl_scalar(il,in,C_Temp:last_C)
             end do
            end do
            close(fileio_unit)
         end if

       if (CP%WantTensors .and. TensFile /= '') then
           open(unit=fileio_unit,file=TensFile,form='formatted',status='replace')
            do in=1,CP%InitPower%nn
             do il=lmin,CP%Max_l_tensor
               write(fileio_unit,'(1I6,4E15.5)')il, fact*Cl_tensor(il, in, CT_Temp:CT_Cross)
             end do
            end do
           close(fileio_unit)
        end if

        if (CP%WantTensors .and. CP%WantScalars .and. TotFile /= '') then
           open(unit=fileio_unit,file=TotFile,form='formatted',status='replace')
           do in=1,CP%InitPower%nn
             do il=lmin,CP%Max_l_tensor

                write(fileio_unit,'(1I6,4E15.5)')il, fact*(Cl_scalar(il, in, C_Temp:C_E)+ Cl_tensor(il,in, C_Temp:C_E)), &
                   fact*Cl_tensor(il,in, CT_B), fact*(Cl_scalar(il, in, C_Cross) + Cl_tensor(il, in, CT_Cross))
             end do
             do il=CP%Max_l_tensor+1,CP%Max_l
                  write(fileio_unit,'(1I6,4E15.5)')il ,fact*Cl_scalar(il,in,C_Temp:C_E), 0._dl, fact*Cl_scalar(il,in,C_Cross)
             end do
           end do
           close(fileio_unit)
        end if

        if (CP%WantScalars .and. CP%DoLensing .and. LensFile /= '') then
           open(unit=fileio_unit,file=LensFile,form='formatted',status='replace')
            do in=1,CP%InitPower%nn
             do il=lmin, lmax_lensed
               write(fileio_unit,'(1I6,4E15.5)')il, fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
             end do
            end do
           close(fileio_unit)
        end if


       if (CP%WantScalars .and. CP%WantTensors .and. CP%DoLensing .and. LensTotFile /= '') then
           open(unit=fileio_unit,file=LensTotFile,form='formatted',status='replace')
           do in=1,CP%InitPower%nn
             do il=lmin,min(CP%Max_l_tensor,lmax_lensed)
                write(fileio_unit,'(1I6,4E15.5)')il, fact*(Cl_lensed(il, in, CT_Temp:CT_Cross)+ Cl_tensor(il,in, CT_Temp:CT_Cross))
             end do
             do il=min(CP%Max_l_tensor,lmax_lensed)+1,lmax_lensed
                write(fileio_unit,'(1I6,4E15.5)')il, fact*Cl_lensed(il, in, CT_Temp:CT_Cross)
             end do
           end do

        end if
        end subroutine output_cl_files

        subroutine output_lens_pot_files(LensPotFile, factor)
      !Write out L TT EE BB TE PP PT PE where P is the lensing potential, all unlensed
      !This input supported by LensPix from 2010
        implicit none
        integer in,il
        real(dl), intent(in), optional :: factor
        real(dl) fact, scale, BB, TT, TE, EE
        character(LEN=*) LensPotFile
         !output file of dimensionless [l(l+1)]^2 C_phi_phi/2pi and [l(l+1)]^(3/2) C_phi_T/2pi
         !This is the format used by Planck_like but original LensPix uses scalar_output_file.

         !(Cl_scalar and scalar_output_file numbers are instead l^4 C_phi and l^3 C_phi
         ! - for historical reasons)

        if (present(factor)) then
          fact = factor
        else
          fact =1
        end if

        if (CP%WantScalars .and. CP%DoLensing .and. LensPotFile/='') then

           open(unit=fileio_unit,file=LensPotFile,form='formatted',status='replace')
           do in=1,CP%InitPower%nn
             do il=lmin,min(10000,CP%Max_l)

               TT = Cl_scalar(il, in, C_Temp)
               EE = Cl_scalar(il, in, C_E)
               TE = Cl_scalar(il, in, C_Cross)
               if (CP%WantTensors .and. il <= CP%Max_l_tensor) then
                TT= TT+Cl_tensor(il,in, CT_Temp)
                EE= EE+Cl_tensor(il,in, CT_E)
                TE= TE+Cl_tensor(il,in, CT_Cross)
                BB= Cl_tensor(il,in, CT_B)
               else
                BB=0
               end if
               scale = (real(il+1)/il)**2/OutputDenominator !Factor to go from old l^4 factor to new

               write(fileio_unit,'(1I6,7E15.5)') il , fact*TT, fact*EE, fact*BB, fact*TE, scale*Cl_scalar(il,in,C_Phi),&
                   (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*Cl_scalar(il,in,C_PhiTemp:C_PhiE)

             end do
             do il=10100,CP%Max_l, 100
               scale = (real(il+1)/il)**2/OutputDenominator
               write(fileio_unit,'(1E15.5,7E15.5)') real(il), fact*Cl_scalar(il,in,C_Temp:C_E),0.,fact*Cl_scalar(il,in,C_Cross), &
                    scale*Cl_scalar(il,in,C_Phi),&
                   (real(il+1)/il)**1.5/OutputDenominator*sqrt(fact)*Cl_scalar(il,in,C_PhiTemp:C_PhiE)
             end do
            end do
            close(fileio_unit)
         end if
        end subroutine output_lens_pot_files


        subroutine output_veccl_files(VecFile, factor)
        implicit none
        integer in,il
        character(LEN=*) VecFile
        real(dl), intent(in), optional :: factor
        real(dl) fact


        if (present(factor)) then
          fact = factor
        else
          fact =1
        end if


       if (CP%WantVectors .and. VecFile /= '') then
           open(unit=fileio_unit,file=VecFile,form='formatted',status='replace')
            do in=1,CP%InitPower%nn
             do il=lmin,CP%Max_l
               write(fileio_unit,'(1I5,4E15.5)')il, fact*Cl_vector(il, in, CT_Temp:CT_Cross)
             end do
            end do

           close(fileio_unit)
        end if

        end subroutine output_veccl_files


        subroutine output_COBElikelihood
          integer in
          do in=1, CP%InitPower%nn
             write(*,*)'COBE Likelihood relative to CP%flat=',COBElikelihoods(in)
          end do
        end  subroutine output_COBElikelihood


      subroutine NormalizeClsAtL(lnorm)
        implicit none
        integer, intent(IN) :: lnorm
        integer in
        real(dl) Norm

         do in=1,CP%InitPower%nn

             if (CP%WantScalars) then
                Norm=1/Cl_scalar(lnorm,in, C_Temp)
                Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_Cross) = Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_Cross) * Norm
             end if

             if (CP%WantTensors) then
                  if (.not.CP%WantScalars) Norm = 1/Cl_tensor(lnorm,in, C_Temp)
                  !Otherwise Norm already set correctly
                  Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) =  &
                    Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) * Norm
             end if
       end do

      end  subroutine NormalizeClsAtL

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine COBEnormalize
        use precision
        use ModelParams


        integer in
        real(dl) xlog10
        real(dl) c10, d1,d2,d3,d4,d5,d6,d7, xlogl, COBE_scale
        real(dl) x1, x2,x3,x4,x5,x6,x7,sy,s,sx,sxy,sxx,delt,d1pr,d1ppr
        real(dl) Ctot(lmin:20)


           if (allocated(COBElikelihoods)) deallocate(COBElikelihoods)
           if (allocated(COBE_scales)) deallocate(COBE_scales)
           allocate(COBElikelihoods(CP%InitPower%nn))
           allocate(COBE_scales(CP%InitPower%nn))



        xlog10=log(10._dl)


! COBE normalization
! fit the spectrum to a quadratic around C_10 with equal weights in logl

        do in=1,CP%InitPower%nn

           if (CP%WantTensors) then
              Ctot =  Cl_tensor(lmin:20, in, C_Temp)
           else
              Ctot = 0
           end if
           if (CP%WantScalars) then
              Ctot=Ctot + Cl_scalar(lmin:20, in, C_Temp)

           end if
           c10=Ctot(10)

           d1=(Ctot(3))/c10-1._dl
           d2=(Ctot(4))/c10-1._dl
           d3=(Ctot(6))/c10-1._dl
           d4=(Ctot(8))/c10-1._dl
           d5=(Ctot(12))/c10-1._dl
           d6=(Ctot(15))/c10-1._dl
           d7=(Ctot(20))/c10-1._dl


           x1=log(3._dl)/xlog10-1._dl
           x2=log(4._dl)/xlog10-1._dl
           x3=log(6._dl)/xlog10-1._dl
           x4=log(8._dl)/xlog10-1._dl
           x5=log(12._dl)/xlog10-1._dl
           x6=log(15._dl)/xlog10-1._dl
           x7=log(20._dl)/xlog10-1._dl
           sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7
           s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7
           sx=x1**3+x2**3+x3**3+x4**3+x5**3+x6**3+x7**3
           sxy=x1**2*d1+x2**2*d2+x3**2*d3+x4**2*d4+ &
              x5**2*d5+x6**2*d6+x7**2*d7
           sxx=x1**4+x2**4+x3**4+x4**4+x5**4+x6**4+x7**4
           delt=s*sxx-sx*sx
           d1pr=(sxx*sy-sx*sxy)/delt
           d1ppr=2._dl*(s*sxy-sx*sy)/delt

! Bunn and White fitting formula
           c10=(0.64575d0+0.02282d0*d1pr+0.01391d0*d1pr*d1pr &
           -0.01819d0*d1ppr-0.00646d0*d1pr*d1ppr &
           +0.00103d0*d1ppr*d1ppr)/c10
! logl
           xlogl=-0.01669d0+1.19895d0*d1pr-0.83527d0*d1pr*d1pr &
                 -0.43541d0*d1ppr-0.03421d0*d1pr*d1ppr &
                 +0.01049d0*d1ppr*d1ppr
          ! write(*,*)'COBE Likelihood relative to CP%flat=',exp(xlogl)
           COBElikelihoods(in) = exp(xlogl)

! density power spectrum normalization;

           COBE_scale=c10/OutputDenominator*1.1d-9
           COBE_scales(in)=COBE_scale

!!$!delta^2 = k^4*(tf)^2*ScalarPower(k,in)*COBE_scale where (tf) is output in the transfer function file
!!$!delta^2 = 4*pi*k^3 P(k)


! C_l normalization; output l(l+1)C_l/twopi
           c10=c10*2.2d-9/fourpi

           if (CP%WantScalars) Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_last) = &
                        Cl_scalar(lmin:CP%Max_l, in, C_Temp:C_last)*c10
           if (CP%WantTensors) Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross) = &
                                    Cl_tensor(lmin:CP%Max_l_tensor, in, CT_Temp:CT_Cross)*c10

          end do !in
         end subroutine COBEnormalize

         subroutine ModelData_Free

             call Free_ClTransfer(CTransScal)
             call Free_ClTransfer(CTransVec)
             call Free_ClTransfer(CTransTens)
             if (allocated(Cl_vector)) deallocate(Cl_vector)
             if (allocated(Cl_tensor)) deallocate(Cl_tensor)
             if (allocated(Cl_scalar)) deallocate(Cl_scalar)
             if (allocated(Cl_lensed)) deallocate(Cl_lensed)
             if (allocated(COBElikelihoods)) deallocate(COBElikelihoods)
             if (allocated(COBE_scales)) deallocate(COBE_scales)

         end subroutine ModelData_Free

        end module ModelData


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    module MassiveNu
      use precision
      use ModelParams
      implicit none
        private

          real(dl), parameter  :: const  = 7._dl/120*pi**4 ! 5.68219698_dl
             !const = int q^3 F(q) dq = 7/120*pi^4
          real(dl), parameter  :: const2 = 5._dl/7/pi**2   !0.072372274_dl
          real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
          real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
          real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl

          integer, parameter  :: nrhopn=2000
          real(dl), parameter :: am_min = 0.01_dl  !0.02_dl
            !smallest a*m_nu to integrate distribution function rather than using series
          real(dl), parameter :: am_max = 600._dl
            !max a*m_nu to integrate

          real(dl),parameter  :: am_minp=am_min*1.1
          real(dl), parameter :: am_maxp=am_max*0.9

          real(dl) dlnam

          real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1,ddr1

          !Sample for massive neutrino momentum
          !These settings appear to be OK for P_k accuate at 1e-3 level
          integer, parameter :: nqmax0=80 !maximum array size of q momentum samples
          real(dl) :: nu_q(nqmax0), nu_int_kernel(nqmax0)

          integer nqmax !actual number of q modes evolves

       public const,Nu_Init,Nu_background, Nu_rho, Nu_drho,  nqmax0, nqmax, &
           nu_int_kernel, nu_q
       contains
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine Nu_init

!  Initialize interpolation tables for massive neutrinos.
!  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.

         integer i
         real(dl) dq,dlfdlq, q, am, rhonu,pnu
         real(dl) spline_data(nrhopn)

!  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
!  Get number density n of neutrinos from
!  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
!  then m = Omega_nu/N_nu rho_crit /n
!  Error due to velocity < 1e-5

        do i=1, CP%Nu_mass_eigenstates
         nu_masses(i)=const/(1.5d0*zeta3)*grhom/grhor*CP%omegan*CP%Nu_mass_fractions(i) &
               /CP%Nu_mass_degeneracies(i)
        end do

        if (allocated(r1)) return
        allocate(r1(nrhopn),p1(nrhopn),dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn))


        nqmax=3
        if (AccuracyBoost >1) nqmax=4
        if (AccuracyBoost >2) nqmax=5
        if (AccuracyBoost >3) nqmax=nint(AccuracyBoost*10)
          !note this may well be worse than the 5 optimized points

        if (nqmax > nqmax0) call MpiStop('Nu_Init: qmax > nqmax0')

        !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
        !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
        !see CAMB notes
        if (nqmax==3) then
          !Accurate at 2e-4 level
          nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
          nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)

        else if (nqmax==4) then
          !This seems to be very accurate (limited by other numerics)
           nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
           nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)

        else if (nqmax==5) then
        !exact for n=-4,-2..3
        !This seems to be very accurate (limited by other numerics)
         nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)
         nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/)

        else
         dq = (12 + nqmax/5)/real(nqmax)
         do i=1,nqmax
            q=(i-0.5d0)*dq
            nu_q(i) = q
            dlfdlq=-q/(1._dl+exp(-q))
            nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)

         end do
        end if
        nu_int_kernel=nu_int_kernel/const

        dlnam=-(log(am_min/am_max))/(nrhopn-1)


        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
        !$OMP & PRIVATE(am, rhonu,pnu)
        do i=1,nrhopn
          am=am_min*exp((i-1)*dlnam)
          call nuRhoPres(am,rhonu,pnu)
          r1(i)=log(rhonu)
          p1(i)=log(pnu)
        end do
        !$OMP END PARALLEL DO


        call splini(spline_data,nrhopn)
        call splder(r1,dr1,nrhopn,spline_data)
        call splder(p1,dp1,nrhopn,spline_data)
        call splder(dr1,ddr1,nrhopn,spline_data)


        end subroutine Nu_init

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine nuRhoPres(am,rhonu,pnu)
!  Compute the density and pressure of one eigenstate of massive neutrinos,
!  in units of the mean density of one flavor of massless neutrinos.

        real(dl),  parameter :: qmax=30._dl
        integer, parameter :: nq=100
        real(dl) dum1(nq+1),dum2(nq+1)
        real(dl), intent(in) :: am
        real(dl), intent(out) ::  rhonu,pnu
        integer i
        real(dl) q,aq,v,aqdn,adq


!  q is the comoving momentum in units of k_B*T_nu0/c.
!  Integrate up to qmax and then use asymptotic expansion for remainder.
        adq=qmax/nq
        dum1(1)=0._dl
        dum2(1)=0._dl
        do  i=1,nq
          q=i*adq
          aq=am/q
          v=1._dl/sqrt(1._dl+aq*aq)
          aqdn=adq*q*q*q/(exp(q)+1._dl)
          dum1(i+1)=aqdn/v
          dum2(i+1)=aqdn*v
        end do
        call splint(dum1,rhonu,nq+1)
        call splint(dum2,pnu,nq+1)
!  Apply asymptotic corrrection for q>qmax and normalize by relativistic
!  energy density.
        rhonu=(rhonu+dum1(nq+1)/adq)/const
        pnu=(pnu+dum2(nq+1)/adq)/const/3._dl

        end subroutine nuRhoPres

!cccccccccccccccccccccccccccccccccccccccccc
       subroutine Nu_background(am,rhonu,pnu)
        use precision
        use ModelParams
        real(dl), intent(in) :: am
        real(dl), intent(out) :: rhonu, pnu

!  Compute massive neutrino density and pressure in units of the mean
!  density of one eigenstate of massless neutrinos.  Use cubic splines to
!  interpolate from a table.

        real(dl) d
        integer i

        if (am <= am_minp) then
          rhonu=1._dl + const2*am**2
          pnu=(2-rhonu)/3._dl
          return
        else if (am >= am_maxp) then
          rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
          pnu = 900._dl/120._dl/const*(zeta5-63._dl/4*Zeta7/am**2)/am
          return
        end if


        d=log(am/am_min)/dlnam+1._dl
        i=int(d)
        d=d-i

!  Cubic spline interpolation.
          rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
               -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
          pnu=p1(i)+d*(dp1(i)+d*(3._dl*(p1(i+1)-p1(i))-2._dl*dp1(i) &
               -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2._dl*(p1(i)-p1(i+1)))))
          rhonu=exp(rhonu)
          pnu=exp(pnu)

        end subroutine Nu_background

!cccccccccccccccccccccccccccccccccccccccccc
       subroutine Nu_rho(am,rhonu)
        use precision
        use ModelParams
        real(dl), intent(in) :: am
        real(dl), intent(out) :: rhonu

!  Compute massive neutrino density in units of the mean
!  density of one eigenstate of massless neutrinos.  Use cubic splines to
!  interpolate from a table.

        real(dl) d
        integer i

        if (am <= am_minp) then
          rhonu=1._dl + const2*am**2
          return
        else if (am >= am_maxp) then
          rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
          return
        end if

        d=log(am/am_min)/dlnam+1._dl
        i=int(d)
        d=d-i

!  Cubic spline interpolation.
        rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
               -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
        rhonu=exp(rhonu)
       end subroutine Nu_rho

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function Nu_drho(am,adotoa,rhonu) result (rhonudot)
        use precision
        use ModelParams

!  Compute the time derivative of the mean density in massive neutrinos
!  and the shear perturbation.
        real(dl) adotoa,rhonu,rhonudot
        real(dl) d
        real(dl), intent(IN) :: am
        integer i

        if (am< am_minp) then

           rhonudot = 2*const2*am**2*adotoa

        else if (am>am_maxp) then

           rhonudot = 3/(2*const)*(zeta3*am - (15*zeta5)/2/am)*adotoa

        else

           d=log(am/am_min)/dlnam+1._dl
           i=int(d)
           d=d-i
           !  Cubic spline interpolation for rhonudot.
           rhonudot=dr1(i)+d*(ddr1(i)+d*(3._dl*(dr1(i+1)-dr1(i)) &
                -2._dl*ddr1(i)-ddr1(i+1)+d*(ddr1(i)+ddr1(i+1) &
                +2._dl*(dr1(i)-dr1(i+1)))))

           rhonudot=rhonu*adotoa*rhonudot/dlnam
        end if

        end function Nu_drho

      end module MassiveNu

! wrapper function to avoid cirular module references
      subroutine init_massive_nu(has_massive_nu)
        use MassiveNu
        use ModelParams
        implicit none
        logical, intent(IN) :: has_massive_nu

        if (has_massive_nu) then
             call Nu_Init
        else
             nu_masses = 0
        end if
      end subroutine init_massive_nu


!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module Transfer
        use ModelData
        use Errors
        implicit none
        public
        integer, parameter :: Transfer_kh =1, Transfer_cdm=2,Transfer_b=3,Transfer_g=4, &
                              Transfer_r=5, Transfer_nu = 6,  & !massless and massive neutrino
                              Transfer_tot=7

        integer, parameter :: Transfer_max = Transfer_tot

        logical :: transfer_interp_matterpower  = .true. !output regular grid in log k
         !set to false to output calculated values for later interpolation

        integer :: transfer_power_var = Transfer_tot
         !What to use to calulcate the output matter power spectrum and sigma_8
         !Transfer_tot uses total matter perturbation

        Type MatterTransferData
         !Computed data
         integer   ::  num_q_trans   !    number of steps in k for transfer calculation
         real(dl), dimension (:), pointer :: q_trans => NULL()
         real(dl), dimension (:,:), pointer ::  sigma_8 => NULL()
         real, dimension(:,:,:), pointer :: TransferData => NULL()
         !TransferData(entry,k_index,z_index) for entry=Tranfer_kh.. Transfer_tot
        end Type MatterTransferData

        Type MatterPowerData
         !everything is a function of k/h
          integer   ::  num_k, num_z
          real(dl), dimension(:), pointer :: log_kh => NULL(), redshifts => NULL()
          !matpower is log(P_k)
          real(dl), dimension(:,:), allocatable :: matpower, ddmat 
          !if NonLinear, nonlin_ratio =  sqrt(P_nonlinear/P_linear)
          !function of k and redshift NonLinearScaling(k_index,z_index)
          real(dl), dimension(:,:), pointer :: nonlin_ratio => NULL()
        end Type MatterPowerData

        Type (MatterTransferData), save :: MT

      contains

        subroutine Transfer_GetMatterPowerData(MTrans, PK_data, in, itf_only)
         !Does *NOT* include non-linear corrections
          !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
          !sepctrum is generated to beyond the CMB k_max
          Type(MatterTransferData), intent(in) :: MTrans
          Type(MatterPowerData) :: PK_data
          integer, intent(in) :: in
          integer, intent(in), optional :: itf_only
          real(dl) h, kh, k, power
          integer ik
          integer nz,itf, itf_start, itf_end

          if (present(itf_only)) then
              itf_start=itf_only
              itf_end = itf_only
              nz = 1
          else
              itf_start=1
              nz= size(MTrans%TransferData,3)
              itf_end = nz
          end if
          PK_data%num_k = MTrans%num_q_trans
          PK_Data%num_z = nz

          allocate(PK_data%matpower(PK_data%num_k,nz))
          allocate(PK_data%ddmat(PK_data%num_k,nz))
          allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
          allocate(PK_data%log_kh(PK_data%num_k))
          allocate(PK_data%redshifts(nz))
          PK_data%redshifts = CP%Transfer%Redshifts(itf_start:itf_end)

          h = CP%H0/100

          do ik=1,MTrans%num_q_trans
                 kh = MTrans%TransferData(Transfer_kh,ik,1)
                 k = kh*h
                 PK_data%log_kh(ik) = log(kh)
                 power = ScalarPower(k,in)
                 if (global_error_flag/=0) then
                     call MatterPowerdata_Free(PK_data)
                     return
                 end if
                 do itf = 1, nz
                   PK_data%matpower(ik,itf) = &
                    log(MTrans%TransferData(transfer_power_var,ik,itf_start+itf-1)**2*k &
                                   *pi*twopi*h**3*power)
                 end do
          end do

          call MatterPowerdata_getsplines(PK_data)

        end subroutine Transfer_GetMatterPowerData

        subroutine MatterPowerData_Load(PK_data,fname)
          !Loads in kh, P_k from file for one redshiftr and one initial power spectrum
          !Not redshift is not stored in file, so not set correctly
          !Also note that output _matterpower file is already interpolated, so re-interpolating is probs not a good idea

          !Get total matter power spectrum in units of (h Mpc^{-1})^3 ready for interpolation.
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          use AmlUtils
          character(LEN=*) :: fname
          Type(MatterPowerData) :: PK_data
          real(dl)kh, Pk
          integer ik
          integer nz


          nz = 1
          call openTxtFile(fname, fileio_unit)

          PK_data%num_k = FileLines(fileio_unit)
          PK_Data%num_z = 1

          allocate(PK_data%matpower(PK_data%num_k,nz))
          allocate(PK_data%ddmat(PK_data%num_k,nz))
          allocate(PK_data%nonlin_ratio(PK_data%num_k,nz))
          allocate(PK_data%log_kh(PK_data%num_k))

          allocate(PK_data%redshifts(nz))
          PK_data%redshifts = 0

          do ik=1,PK_data%num_k
              read (fileio_unit,*) kh, Pk
              PK_data%matpower(ik,1) = log(Pk)
              PK_data%log_kh(ik) = log(kh)
          end do

          call MatterPowerdata_getsplines(PK_data)

        end subroutine MatterPowerData_Load


        subroutine MatterPowerdata_getsplines(PK_data)
          Type(MatterPowerData) :: PK_data
          integer i
          real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl

          do i = 1,PK_Data%num_z

           call spline(PK_data%log_kh,PK_data%matpower(1,i),PK_data%num_k,&
                               cllo,clhi,PK_data%ddmat(1,i))
          end do

        end subroutine MatterPowerdata_getsplines

        subroutine MatterPowerdata_MakeNonlinear(PK_data)
          Type(MatterPowerData) :: PK_data

          call NonLinear_GetRatios(PK_data)
          PK_data%matpower = PK_data%matpower +  2*log(PK_data%nonlin_ratio)
          call MatterPowerdata_getsplines(PK_data)

        end subroutine MatterPowerdata_MakeNonlinear

        subroutine MatterPowerdata_Free(PK_data)
          Type(MatterPowerData) :: PK_data
          integer i

          deallocate(PK_data%log_kh,stat=i)
          deallocate(PK_data%matpower,stat=i)
          deallocate(PK_data%ddmat,stat=i)
          deallocate(PK_data%nonlin_ratio,stat=i)
          deallocate(PK_data%redshifts,stat=i)
          nullify(PK_data%log_kh, & !PK_data%matpower,PK_data%ddmat 
             PK_data%nonlin_ratio,PK_data%redshifts)

        end subroutine MatterPowerdata_Free

        function MatterPowerData_k(PK,  kh, itf) result(outpower)
         !Get matter power spectrum at particular k/h by interpolation
          Type(MatterPowerData) :: PK
          integer, intent(in) :: itf
          real (dl), intent(in) :: kh
          real(dl) :: logk
          integer llo,lhi
          real(dl) outpower, dp
          real(dl) ho,a0,b0
          integer, save :: i_last = 1

           logk = log(kh)
           if (logk < PK%log_kh(1)) then
              dp = (PK%matpower(2,itf) -  PK%matpower(1,itf)) / &
                 ( PK%log_kh(2)-PK%log_kh(1) )
              outpower = PK%matpower(1,itf) + dp*(logk - PK%log_kh(1))
           else if (logk > PK%log_kh(PK%num_k)) then
            !Do dodgy linear extrapolation on assumption accuracy of result won't matter

             dp = (PK%matpower(PK%num_k,itf) -  PK%matpower(PK%num_k-1,itf)) / &
                 ( PK%log_kh(PK%num_k)-PK%log_kh(PK%num_k-1) )
             outpower = PK%matpower(PK%num_k,itf) + dp*(logk - PK%log_kh(PK%num_k))
           else

            llo=min(i_last,PK%num_k)
            do while (PK%log_kh(llo) > logk)
               llo=llo-1
            end do
            do while (PK%log_kh(llo+1)< logk)
               llo=llo+1
            end do
            i_last =llo
            lhi=llo+1
            ho=PK%log_kh(lhi)-PK%log_kh(llo)
            a0=(PK%log_kh(lhi)-logk)/ho
            b0=1-a0

            outpower = a0*PK%matpower(llo,itf)+ b0*PK%matpower(lhi,itf)+&
                  ((a0**3-a0)* PK%ddmat(llo,itf) &
                       +(b0**3-b0)*PK%ddmat(lhi,itf))*ho**2/6

          end if

          outpower = exp(max(-30._dl,outpower))

        end function MatterPowerData_k


        subroutine Transfer_GetMatterPower(MTrans,outpower, itf, in, minkh, dlnkh, npoints)
          !Allows for non-smooth priordial spectra
          !if CP%Nonlinear/ = NonLinear_none includes non-linear evolution
          !Get total matter power spectrum at logarithmically equal intervals dlnkh of k/h starting at minkh
          !in units of (h Mpc^{-1})^3.
          !Here there definition is < Delta^2(x) > = 1/(2 pi)^3 int d^3k P_k(k)
          !We are assuming that Cls are generated so any baryonic wiggles are well sampled and that matter power
          !sepctrum is generated to beyond the CMB k_max
          Type(MatterTransferData), intent(in) :: MTrans
          Type(MatterPowerData) :: PK

          integer, intent(in) :: itf, in, npoints
          real, intent(out) :: outpower(npoints)
          real, intent(in) :: minkh, dlnkh
          real(dl), parameter :: cllo=1.e30_dl,clhi=1.e30_dl
          integer ik, llo,il,lhi,lastix
          real(dl) matpower(MTrans%num_q_trans), kh, kvals(MTrans%num_q_trans), ddmat(MTrans%num_q_trans)
          real(dl) atransfer,xi, a0, b0, ho, logmink,k, h


          if (npoints < 2) stop 'Need at least 2 points in Transfer_GetMatterPower'

!         if (minkh < MTrans%TransferData(Transfer_kh,1,itf)) then
!            stop 'Transfer_GetMatterPower: kh out of computed region'
!          end if
          if (minkh*exp((npoints-1)*dlnkh) > MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf) &
                .and. FeedbackLevel > 0 ) &
                    write(*,*) 'Warning: extrapolating matter power in Transfer_GetMatterPower'


          if (CP%NonLinear/=NonLinear_None) then
           call Transfer_GetMatterPowerData(MTrans, PK, in, itf)
           call NonLinear_GetRatios(PK)
          end if

          h = CP%H0/100
          logmink = log(minkh)
          do ik=1,MTrans%num_q_trans
             kh = MTrans%TransferData(Transfer_kh,ik,itf)
             k = kh*h
             kvals(ik) = log(kh)
             atransfer=MTrans%TransferData(transfer_power_var,ik,itf)
             if (CP%NonLinear/=NonLinear_None) &
                 atransfer = atransfer* PK%nonlin_ratio(ik,1) !only one element, this itf
             matpower(ik) = log(atransfer**2*k*pi*twopi*h**3)
                 !Put in power spectrum later: transfer functions should be smooth, initial power may not be
          end do

          call spline(kvals,matpower,MTrans%num_q_trans,cllo,clhi,ddmat)

            llo=1
            lastix = npoints + 1
            do il=1, npoints
               xi=logmink + dlnkh*(il-1)
               if (xi < kvals(1)) then
                 outpower(il)=-30.
                 cycle
               end if
               do while ((xi > kvals(llo+1)).and.(llo < MTrans%num_q_trans))
                  llo=llo+1
                  if (llo >= MTrans%num_q_trans) exit
               end do
               if (llo == MTrans%num_q_trans) then
                   lastix = il
                   exit
               end if
               lhi=llo+1
               ho=kvals(lhi)-kvals(llo)
               a0=(kvals(lhi)-xi)/ho
               b0=(xi-kvals(llo))/ho

               outpower(il) = a0*matpower(llo)+ b0*matpower(lhi)+((a0**3-a0)* ddmat(llo) &
                       +(b0**3-b0)*ddmat(lhi))*ho**2/6

            end do

            do while (lastix <= npoints)
               !Do linear extrapolation in the log
               !Obviouly inaccurate, non-linear etc, but OK if only using in tails of window functions
               outpower(lastix) = 2*outpower(lastix-1) - outpower(lastix-2)
               lastix = lastix+1
            end do

            outpower = exp(max(-30.,outpower))

            do il = 1, npoints
               k = exp(logmink + dlnkh*(il-1))*h
               outpower(il) = outpower(il) * ScalarPower(k,in)
               if (global_error_flag /= 0) exit
            end do

          if (CP%NonLinear/=NonLinear_None) call MatterPowerdata_Free(PK)

        end subroutine Transfer_GetMatterPower

        subroutine Transfer_Get_sigma8(MTrans, sigr8)
          use MassiveNu
          Type(MatterTransferData) :: MTrans
          integer ik, itf, in
          real(dl) kh, k, h, x, win, delta
          real(dl) lnk, dlnk, lnko
          real(dl) dsig8, dsig8o, sig8, sig8o, powers
          real(dl), intent(IN) :: sigr8

          !Calculate MTrans%sigma_8^2 = int dk/k win**2 T_k**2 P(k), where win is the FT of a spherical top hat
          !of radius sigr8 h^{-1} Mpc

           if (global_error_flag /= 0) return

         H=CP%h0/100._dl
         do in = 1, CP%InitPower%nn
          do itf=1,CP%Transfer%num_redshifts
            lnko=0
            dsig8o=0
            sig8=0
            sig8o=0
          do ik=1, MTrans%num_q_trans
               kh = MTrans%TransferData(Transfer_kh,ik,itf)
               if (kh==0) cycle
               k = kh*H

               delta = k**2*MTrans%TransferData(transfer_power_var,ik,itf)
               !if (CP%NonLinear/=NonLinear_None) delta= delta* MTrans%NonLinearScaling(ik,itf)
               !sigma_8 defined "as though it were linear"

               x= kh *sigr8
               win =3*(sin(x)-x*cos(x))/x**3
               lnk=log(k)
               if (ik==1) then
                  dlnk=0.5_dl
                 !Approx for 2._dl/(CP%InitPower%an(in)+3)  [From int_0^k_1 dk/k k^4 P(k)]
                 !Contribution should be very small in any case
               else
                  dlnk=lnk-lnko
               end if
               powers = ScalarPower(k,in)
               dsig8=(win*delta)**2*powers
               sig8=sig8+(dsig8+dsig8o)*dlnk/2
               dsig8o=dsig8
               lnko=lnk


          end do

          MTrans%sigma_8(itf,in) = sqrt(sig8)
          end do
         end do

        end subroutine Transfer_Get_sigma8

        subroutine Transfer_output_Sig8(MTrans)
           Type(MatterTransferData), intent(in) :: MTrans

           integer in, j

           do in=1, CP%InitPower%nn
            if (CP%InitPower%nn>1)  write(*,*) 'Power spectrum : ', in
            do j = 1, CP%Transfer%num_redshifts
               write(*,*) 'at z = ',real(CP%Transfer%redshifts(j)), ' sigma8 (all matter)=', real(MTrans%sigma_8(j,in))
            end do
           end do

         end subroutine Transfer_output_Sig8


        subroutine Transfer_output_Sig8AndNorm(MTrans)
           Type(MatterTransferData), intent(in) :: MTrans
           integer in, j

           do in=1, CP%InitPower%nn
             write(*,*) 'Power spectrum ',in, ' COBE_scale = ',real(COBE_scales(in))
            do j = 1, CP%Transfer%num_redshifts
               write(*,*) 'at z = ',real(CP%Transfer%redshifts(j)), ' sigma8(all matter) = ', &
                    real(MTrans%sigma_8(j,in)*sqrt(COBE_scales(in)))
            end do
           end do

         end subroutine Transfer_output_Sig8AndNorm


        subroutine Transfer_Allocate(MTrans)
         Type(MatterTransferData) :: MTrans
         integer st

          deallocate(MTrans%q_trans, STAT = st)
          deallocate(MTrans%TransferData, STAT = st)
          deallocate(MTrans%sigma_8, STAT = st)
          allocate(MTrans%q_trans(MTrans%num_q_trans))
          allocate(MTrans%TransferData(Transfer_max,MTrans%num_q_trans,CP%Transfer%num_redshifts))
          allocate(MTrans%sigma_8(CP%Transfer%num_redshifts, CP%InitPower%nn))

        end  subroutine Transfer_Allocate

        subroutine Transfer_Free(MTrans)
          Type(MatterTransferData):: MTrans
          integer st

          deallocate(MTrans%q_trans, STAT = st)
          deallocate(MTrans%TransferData, STAT = st)
          deallocate(MTrans%sigma_8, STAT = st)
          nullify(MTrans%q_trans)
          nullify(MTrans%TransferData)
          nullify(MTrans%sigma_8)

        end subroutine Transfer_Free

       subroutine Transfer_SetForNonlinearLensing(P)
          Type(TransferParams) :: P
          integer i
          real maxRedshift

          P%kmax = 5*AccuracyBoost
          P%k_per_logint  = 0
          maxRedshift = 10
          P%num_redshifts =  nint(10*AccuracyBoost)
          if (HighAccuracyDefault) then
           !only notionally more accuracy, more stable for RS
           maxRedshift =15
           P%num_redshifts = P%num_redshifts *3
          end if
          if (P%num_redshifts > max_transfer_redshifts) &
                stop 'Transfer_SetForNonlinearLensing: Too many redshifts'
          do i=1,P%num_redshifts
           P%redshifts(i) = real(P%num_redshifts-i)/(P%num_redshifts/maxRedshift)

          end do

       end subroutine Transfer_SetForNonlinearLensing



        subroutine Transfer_SaveToFiles(MTrans,FileNames)
          use IniFile
          Type(MatterTransferData), intent(in) :: MTrans
          integer i,ik
          character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)

          do i=1, CP%Transfer%num_redshifts
            if (FileNames(i) /= '') then
            open(unit=fileio_unit,file=FileNames(i),form='formatted',status='replace')
             do ik=1,MTrans%num_q_trans
                if (MTrans%TransferData(Transfer_kh,ik,i)/=0) then
                 write(fileio_unit,'(7E14.6)') MTrans%TransferData(Transfer_kh:Transfer_max,ik,i)
                end if
             end do
            close(fileio_unit)
            end if
          end do


        end subroutine Transfer_SaveToFiles

        subroutine Transfer_SaveMatterPower(MTrans, FileNames)
          use IniFile
          !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
          Type(MatterTransferData), intent(in) :: MTrans
          character(LEN=Ini_max_string_len), intent(IN) :: FileNames(*)
          integer itf,in,i
          integer points
          real, dimension(:,:), allocatable :: outpower
          character(LEN=80) fmt
          real minkh,dlnkh
          Type(MatterPowerData) :: PK_data


          write (fmt,*) CP%InitPower%nn+1
          fmt = '('//trim(adjustl(fmt))//'E15.5)'
          do itf=1, CP%Transfer%num_redshifts
            if (FileNames(itf) /= '') then


             if (.not. transfer_interp_matterpower ) then

             points = MTrans%num_q_trans
             allocate(outpower(points,CP%InitPower%nn))

                 do in = 1, CP%InitPower%nn

                   call Transfer_GetMatterPowerData(MTrans, PK_data, in, itf)

                  if (CP%NonLinear/=NonLinear_None) call MatterPowerdata_MakeNonlinear(PK_Data)

                   outpower(:,in) = exp(PK_data%matpower(:,1))
                   call MatterPowerdata_Free(PK_Data)
                 end do

                 open(unit=fileio_unit,file=FileNames(itf),form='formatted',status='replace')
                 do i=1,points
                  write (fileio_unit, fmt) MTrans%TransferData(Transfer_kh,i,1),outpower(i,1:CP%InitPower%nn)
                 end do
                 close(fileio_unit)

             else


             minkh = 1e-4
             dlnkh = 0.02
             points = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/dlnkh+1
!             dlnkh = log(MTrans%TransferData(Transfer_kh,MTrans%num_q_trans,itf)/minkh)/(points-0.999)
             allocate(outpower(points,CP%InitPower%nn))
             do in = 1, CP%InitPower%nn
              call Transfer_GetMatterPower(MTrans,outpower(1,in), itf, in, minkh,dlnkh, points)
              if (CP%OutputNormalization == outCOBE) then
                 if (allocated(COBE_scales)) then
                  outpower(:,in) = outpower(:,in)*COBE_scales(in)
                 else
                  if (FeedbackLevel>0) write (*,*) 'Cannot COBE normalize - no Cls generated'
                 end if
             end if
             end do

             open(unit=fileio_unit,file=FileNames(itf),form='formatted',status='replace')
             do i=1,points
              write (fileio_unit, fmt) minkh*exp((i-1)*dlnkh),outpower(i,1:CP%InitPower%nn)
             end do
             close(fileio_unit)

             end if

             deallocate(outpower)

            end if
          end do

        end subroutine Transfer_SaveMatterPower

        end module Transfer


!ccccccccccccccccccccccccccccccccccccccccccccccccccc

        module ThermoData
        use ModelData
        implicit none
        private
        integer,parameter :: nthermo=20000

        real(dl) tb(nthermo),cs2(nthermo),xe(nthermo)
        real(dl) dcs2(nthermo)
        real(dl) dotmu(nthermo), ddotmu(nthermo)
        real(dl) sdotmu(nthermo),emmu(nthermo)
        real(dl) demmu(nthermo)
        real(dl) dddotmu(nthermo),ddddotmu(nthermo)
        real(dl) winlens(nthermo),dwinlens(nthermo), scalefactor(nthermo)
        real(dl) tauminn,dlntau,Maxtau
        real(dl), dimension(:), allocatable :: vis,dvis,ddvis,expmmu,dopac, opac, lenswin
        logical, parameter :: dowinlens = .false.

        real(dl) :: tight_tau, actual_opt_depth
         !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
        real(dl) :: matter_verydom_tau
        real(dl) :: r_drag0, z_star, z_drag  !!JH for updated BAO likelihood.

        public thermo,inithermo,vis,opac,expmmu,dvis,dopac,ddvis,lenswin, tight_tau,&
               Thermo_OpacityToTime,matter_verydom_tau, ThermoData_Free,&
               z_star, z_drag  !!JH for updated BAO likelihood.
       contains

        subroutine thermo(tau,cs2b,opacity, dopacity)
        !Compute unperturbed sound speed squared,
        !and ionization fraction by interpolating pre-computed tables.
        !If requested also get time derivative of opacity
        implicit none
        real(dl) tau,cs2b,opacity
        real(dl), intent(out), optional :: dopacity

        integer i
        real(dl) d

        d=log(tau/tauminn)/dlntau+1._dl
        i=int(d)
        d=d-i
        if (i < 1) then
        !Linear interpolation if out of bounds (should not occur).
          cs2b=cs2(1)+(d+i-1)*dcs2(1)
          opacity=dotmu(1)+(d-1)*ddotmu(1)
          stop 'thermo out of bounds'
        else if (i >= nthermo) then
          cs2b=cs2(nthermo)+(d+i-nthermo)*dcs2(nthermo)
          opacity=dotmu(nthermo)+(d-nthermo)*ddotmu(nthermo)
          if (present(dopacity)) then
             dopacity = 0
             stop 'thermo: shouldn''t happen'
           end if
        else
        !Cubic spline interpolation.
          cs2b=cs2(i)+d*(dcs2(i)+d*(3*(cs2(i+1)-cs2(i))  &
              -2*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)  &
              +2*(cs2(i)-cs2(i+1)))))
          opacity=dotmu(i)+d*(ddotmu(i)+d*(3*(dotmu(i+1)-dotmu(i)) &
              -2*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
              +2*(dotmu(i)-dotmu(i+1)))))

         if (present(dopacity)) then

          dopacity=(ddotmu(i)+d*(dddotmu(i)+d*(3*(ddotmu(i+1)  &
              -ddotmu(i))-2*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
              +dddotmu(i+1)+2*(ddotmu(i)-ddotmu(i+1))))))/(tau*dlntau)

         end if
        end if
        end subroutine thermo



       function Thermo_OpacityToTime(opacity)
         real(dl), intent(in) :: opacity
         integer j
         real(dl) Thermo_OpacityToTime
         !Do this the bad slow way for now..
          !The answer is approximate
         j =1
         do while(dotmu(j)> opacity)
            j=j+1
         end do

         Thermo_OpacityToTime = exp((j-1)*dlntau)*tauminn

       end function Thermo_OpacityToTime

     subroutine inithermo(taumin,taumax)
!  Compute and save unperturbed baryon temperature and ionization fraction
!  as a function of time.  With nthermo=10000, xe(tau) has a relative
! accuracy (numerical integration precision) better than 1.e-5.
        use constants
        use precision
        use ModelParams
        use MassiveNu
        real(dl) taumin,taumax


        real(dl) tau01,adot0,a0,a02,x1,x2,barssc,dtau
        real(dl) xe0,tau,a,a2
        real(dl) adot,tg0,ahalf,adothalf,fe,thomc,thomc0,etc,a2t
        real(dl) dtbdla,vfi,cf1,maxvis, vis
        integer ncount,i,j1,j2,iv,ns
        real(dl) spline_data(nthermo)
        real(dl) last_dotmu
        real(dl) dtauda  !diff of tau w.CP%r.t a and integration
        external dtauda
        real(dl) a_verydom
        real(dl) awin_lens1,awin_lens2,dwing_lens, rs, DA
        real(dl) rombint
        external rombint

        call Recombination_Init(CP%Recomb, CP%omegac, CP%omegab,CP%Omegan, & 
        CP%Omegav, CP%h0,CP%tcmb,CP%yhe,CP%Num_Nu_massless + CP%Num_Nu_massive)
          !almost all the time spent here
        if (global_error_flag/=0) return
        Maxtau=taumax
        tight_tau = 0
        actual_opt_depth = 0
        ncount=0
        z_star=0.d0
        z_drag=0.d0
        thomc0= Compton_CT * CP%tcmb**4
        r_drag0 = 3.d0/4.d0*CP%omegab*grhom/grhog
        !thomc0=5.0577d-8*CP%tcmb**4

        tauminn=0.05d0*taumin
        dlntau=log(CP%tau0/tauminn)/(nthermo-1)
        last_dotmu = 0

        matter_verydom_tau = 0
        a_verydom = AccuracyBoost*5*(grhog+grhornomass)/(grhoc+grhob)

!  Initial conditions: assume radiation-dominated universe.
        tau01=tauminn
        adot0=adotrad
        a0=adotrad*tauminn
        a02=a0*a0
!  Assume that any entropy generation occurs before tauminn.
!  This gives wrong temperature before pair annihilation, but
!  the error is harmless.
        tb(1)=CP%tcmb/a0
        xe0=1._dl
        x1=0._dl
        x2=1._dl
        xe(1)=xe0+0.25d0*CP%yhe/(1._dl-CP%yhe)*(x1+2*x2)
        barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(1))
        cs2(1)=4._dl/3._dl*barssc*tb(1)
        dotmu(1)=xe(1)*akthom/a02
        sdotmu(1)=0

          do i=2,nthermo
          tau=tauminn*exp((i-1)*dlntau)
          dtau=tau-tau01
!  Integrate Friedmann equation using inverse trapezoidal rule.

          a=a0+adot0*dtau
          scaleFactor(i)=a
          a2=a*a

          adot=1/dtauda(a)

          if (matter_verydom_tau ==0 .and. a > a_verydom) then
             matter_verydom_tau = tau
          end if

          a=a0+2._dl*dtau/(1._dl/adot0+1._dl/adot)
!  Baryon temperature evolution: adiabatic except for Thomson cooling.
!  Use  quadrature solution.
! This is redundant as also calculated in REFCAST, but agrees well before reionization
          tg0=CP%tcmb/a0
          ahalf=0.5d0*(a0+a)
          adothalf=0.5d0*(adot0+adot)
!  fe=number of free electrons divided by total number of free baryon
!  particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
!  more accuracy is required (unlikely) then this can be iterated with
!  the solution of the ionization equation.
          fe=(1._dl-CP%yhe)*xe(i-1)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i-1))
          thomc=thomc0*fe/adothalf/ahalf**3
          etc=exp(-thomc*(a-a0))
          a2t=a0*a0*(tb(i-1)-tg0)*etc-CP%tcmb/thomc*(1._dl-etc)
          tb(i)=CP%tcmb/a+a2t/(a*a)

! If there is re-ionization, smoothly increase xe to the
! requested value.
          if (CP%Reion%Reionization .and. tau > CP%ReionHist%tau_start) then
             if(ncount == 0) then
                ncount=i-1
             end if

            xe(i) = Reionization_xe(a, tau, xe(ncount))
            !print *,1/a-1,xe(i)
            if (CP%AccurateReionization .and. FeedbackLevel > 0) then
                dotmu(i)=(Recombination_xe(a) - xe(i))*akthom/a2

                if (last_dotmu /=0) then
                 actual_opt_depth = actual_opt_depth - 2._dl*dtau/(1._dl/dotmu(i)+1._dl/last_dotmu)
                end if
                last_dotmu = dotmu(i)
            end if

          else
            xe(i)=Recombination_xe(a)
          end if


!  Baryon sound speed squared (over c**2).
          dtbdla=-2._dl*tb(i)-thomc*adothalf/adot*(a*tb(i)-CP%tcmb)
          barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*xe(i))
          cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3._dl)


! Calculation of the visibility function
          dotmu(i)=xe(i)*akthom/a2

          if (tight_tau==0 .and. 1/(tau*dotmu(i)) > 0.005) tight_tau = tau !0.005
           !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)

          if (tau < 0.001) then
             sdotmu(i)=0
          else
             sdotmu(i)=sdotmu(i-1)+2._dl*dtau/(1._dl/dotmu(i)+1._dl/dotmu(i-1))
          end if

          a0=a
          tau01=tau
          adot0=adot
          end do !i

          if (CP%Reion%Reionization .and. (xe(nthermo) < 0.999d0)) then
             write(*,*)'Warning: xe at redshift zero is < 1'
             write(*,*) 'Check input parameters an Reionization_xe'
             write(*,*) 'function in the Reionization module'
          end if

        do j1=1,nthermo
           if (sdotmu(j1) - sdotmu(nthermo)< -69) then
           emmu(j1)=1.d-30
           else
           emmu(j1)=exp(sdotmu(j1)-sdotmu(nthermo))
           if (.not. CP%AccurateReionization .and. &
               actual_opt_depth==0 .and. xe(j1) < 1e-3) then
              actual_opt_depth = -sdotmu(j1)+sdotmu(nthermo)
           end if
           if (CP%AccurateReionization .and. z_star==0.d0) then
              if (sdotmu(nthermo)-sdotmu(j1) - actual_opt_depth < 1) then
                tau01=1-(sdotmu(nthermo)-sdotmu(j1) - actual_opt_depth)
                tau01=tau01*(1._dl/dotmu(j1)+1._dl/dotmu(j1-1))/2
                z_star = 1/(scaleFactor(j1)- tau01/dtauda(scaleFactor(j1))) -1
              end if
           end if
          end if
        end do

        if (CP%AccurateReionization .and. FeedbackLevel > 0 .and. CP%DerivedParameters) then
         write(*,'("Reion opt depth      = ",f7.4)') actual_opt_depth
        end if


        iv=0
        vfi=0._dl
! Getting the starting and finishing times for decoupling and time of maximum visibility
        if (ncount == 0) then
           cf1=1._dl
           ns=nthermo
           else
              cf1=exp(sdotmu(nthermo)-sdotmu(ncount))
              ns=ncount
           end if
         maxvis = 0
         do j1=1,ns
           vis = emmu(j1)*dotmu(j1)
           tau = tauminn*exp((j1-1)*dlntau)
           vfi=vfi+vis*cf1*dlntau*tau
           if ((iv == 0).and.(vfi > 1.0d-7/AccuracyBoost)) then
              taurst=9._dl/10._dl*tau
              iv=1
           elseif (iv == 1) then
               if (vis > maxvis) then
                maxvis=vis
                tau_maxvis = tau
               end if
               if (vfi > 0.995) then
                taurend=tau
                iv=2
                exit
               end if
           end if
         end do

           if (iv /= 2) then
             call GlobalError('inithermo: failed to find end of recombination',error_reionization)
             return
           end if

          if (dowinlens) then
             vfi=0
             awin_lens1=0
              awin_lens2=0
              winlens=0
              do j1=1,nthermo-1
                   vis = emmu(j1)*dotmu(j1)
                   tau = tauminn*exp((j1-1)*dlntau)
                   vfi=vfi+vis*cf1*dlntau*tau
                    if (vfi < 0.995) then
                          dwing_lens =  vis*cf1*dlntau*tau / 0.995

                          awin_lens1 = awin_lens1 + dwing_lens
                          awin_lens2 = awin_lens2 + dwing_lens/(CP%tau0-tau)
                      end if
                      winlens(j1)= awin_lens1/(CP%tau0-tau) - awin_lens2
                 end do
         end if

! Calculating the timesteps during recombination.

           if (CP%WantTensors) then
              dtaurec=min(dtaurec,taurst/160)/AccuracyBoost
           else
              dtaurec=min(dtaurec,taurst/40)/AccuracyBoost
              if (do_bispectrum .and. hard_bispectrum) dtaurec = dtaurec / 4
           end if

           if (CP%Reion%Reionization) taurend=min(taurend,CP%ReionHist%tau_start)

         if (DebugMsgs) then
           write (*,*) 'taurst, taurend = ', taurst, taurend
         end if

        call splini(spline_data,nthermo)
        call splder(cs2,dcs2,nthermo,spline_data)
        call splder(dotmu,ddotmu,nthermo,spline_data)
        call splder(ddotmu,dddotmu,nthermo,spline_data)
        call splder(dddotmu,ddddotmu,nthermo,spline_data)
        call splder(emmu,demmu,nthermo,spline_data)
        if (dowinlens) call splder(winlens,dwinlens,nthermo,spline_data)

        call SetTimeSteps

        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC)
        do j2=1,TimeSteps%npoints
             call DoThermoSpline(j2,TimeSteps%points(j2))
        end do
         !$OMP END PARALLEL DO


        if ((CP%want_zstar .or. CP%DerivedParameters) .and. z_star==0.d0) call find_z(optdepth,z_star)
        if (CP%want_zdrag .or. CP%DerivedParameters) call find_z(dragoptdepth,z_drag)

        if (CP%DerivedParameters) then
            rs =rombint(dsound_da,1d-8,1/(z_star+1),1d-6)
            DA = AngularDiameterDistance(z_star)/(1/(z_star+1))

            ThermoDerivedParams( derived_zstar ) = z_star
            ThermoDerivedParams( derived_rstar ) = rs
            ThermoDerivedParams( derived_thetastar ) = 100*rs/DA
            ThermoDerivedParams( derived_zdrag ) = z_drag
            rs =rombint(dsound_da,1d-8,1/(z_drag+1),1d-6)
            ThermoDerivedParams( derived_rdrag ) = rs
            ThermoDerivedParams( derived_kD ) =  sqrt(1.d0/(rombint(ddamping_da, 1d-8, 1/(z_star+1), 1d-6)/6))
            ThermoDerivedParams( derived_thetaD ) =  100*pi/ThermoDerivedParams( derived_kD )/DA
            ThermoDerivedParams( derived_zEQ ) = (grhob+grhoc)/(grhog+grhornomass+sum(grhormass(1:CP%Nu_mass_eigenstates))) -1
            ThermoDerivedParams( derived_thetaEQ ) = 100*timeOfz( ThermoDerivedParams( derived_zEQ ))/DA
           ! ThermoDerivedParams( derived_mnu ) = 0
           ! if (CP%Num_Nu_Massive > 0) then
           !         do nu_i=1, CP%Nu_mass_eigenstates
           !             ThermoDerivedParams( derived_mnu ) =  &
           !         ThermoDerivedParams( derived_mnu ) + CP%nu_mass_degeneracies(nu_i)*1.68e-4*nu_masses(nu_i)
           !         end do
           ! end if

            if (FeedbackLevel > 0) then

                        write(*,'("zstar                = ",f8.2)') ThermoDerivedParams( derived_zstar )
                        write(*,'("r_s(zstar)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rstar )
                        write(*,'("100*theta            = ",f9.6)') ThermoDerivedParams( derived_thetastar )

                        write(*,'("zdrag                = ",f8.2)') ThermoDerivedParams( derived_zdrag )
                        write(*,'("r_s(zdrag)/Mpc       = ",f7.2)') ThermoDerivedParams( derived_rdrag )

                        write(*,'("k_D(zstar) Mpc       = ",f7.4)') ThermoDerivedParams( derived_kD )
                        write(*,'("100*theta_D          = ",f9.6)') ThermoDerivedParams( derived_thetaD )

                        write(*,'("z_EQ (if v_nu=1)     = ",f8.2)') ThermoDerivedParams( derived_zEQ )
                        write(*,'("100*theta_EQ         = ",f9.6)') ThermoDerivedParams( derived_thetaEQ )
             end if
        end if

     end subroutine inithermo


        subroutine SetTimeSteps
        real(dl) dtau0
        integer nri0, nstep

         call Ranges_Init(TimeSteps)

         call Ranges_Add_delta(TimeSteps, taurst, taurend, dtaurec)

        ! Calculating the timesteps after recombination
           if (CP%WantTensors) then
              dtau0=max(taurst/40,Maxtau/2000._dl/AccuracyBoost)
           else
              dtau0=Maxtau/500._dl/AccuracyBoost
              if (do_bispectrum) dtau0 = dtau0/3
             !Don't need this since adding in Limber on small scales
              !  if (CP%DoLensing) dtau0=dtau0/2
              !  if (CP%AccurateBB) dtau0=dtau0/3 !Need to get C_Phi accurate on small scales
           end if

         call Ranges_Add_delta(TimeSteps,taurend, CP%tau0, dtau0)

         if (CP%Reion%Reionization) then

              nri0=int(Reionization_timesteps(CP%ReionHist)*AccuracyBoost)
                !Steps while reionization going from zero to maximum
              call Ranges_Add(TimeSteps,CP%ReionHist%tau_start,CP%ReionHist%tau_complete,nri0)

         end if

!Create arrays out of the region information.
        call Ranges_GetArray(TimeSteps)
        nstep = TimeSteps%npoints

        if (allocated(vis)) then
           deallocate(vis,dvis,ddvis,expmmu,dopac, opac)
           if (dowinlens) deallocate(lenswin)
        end if
        allocate(vis(nstep),dvis(nstep),ddvis(nstep),expmmu(nstep),dopac(nstep),opac(nstep))
        if (dowinlens) allocate(lenswin(nstep))

        if (DebugMsgs .and. FeedbackLevel > 0) write(*,*) 'Set ',nstep, ' time steps'

        end subroutine SetTimeSteps


        subroutine ThermoData_Free
         if (allocated(vis)) then
           deallocate(vis,dvis,ddvis,expmmu,dopac, opac)
           if (dowinlens) deallocate(lenswin)
         end if
         call Ranges_Free(TimeSteps)

        end subroutine ThermoData_Free

!cccccccccccccc
        subroutine DoThermoSpline(j2,tau)
        integer j2,i
        real(dl) d,ddopac,tau

!     Cubic-spline interpolation.
           d=log(tau/tauminn)/dlntau+1._dl
           i=int(d)

           d=d-i
           if (i < nthermo) then
          opac(j2)=dotmu(i)+d*(ddotmu(i)+d*(3._dl*(dotmu(i+1)-dotmu(i)) &
              -2._dl*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1) &
              +2._dl*(dotmu(i)-dotmu(i+1)))))
          dopac(j2)=(ddotmu(i)+d*(dddotmu(i)+d*(3._dl*(ddotmu(i+1)  &
              -ddotmu(i))-2._dl*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i) &
              +dddotmu(i+1)+2._dl*(ddotmu(i)-ddotmu(i+1))))))/(tau &
              *dlntau)
          ddopac=(dddotmu(i)+d*(ddddotmu(i)+d*(3._dl*(dddotmu(i+1) &
              -dddotmu(i))-2._dl*ddddotmu(i)-ddddotmu(i+1)  &
              +d*(ddddotmu(i)+ddddotmu(i+1)+2._dl*(dddotmu(i) &
              -dddotmu(i+1)))))-(dlntau**2)*tau*dopac(j2)) &
              /(tau*dlntau)**2
          expmmu(j2)=emmu(i)+d*(demmu(i)+d*(3._dl*(emmu(i+1)-emmu(i)) &
              -2._dl*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1) &
              +2._dl*(emmu(i)-emmu(i+1)))))

          if (dowinlens) then
          lenswin(j2)=winlens(i)+d*(dwinlens(i)+d*(3._dl*(winlens(i+1)-winlens(i)) &
              -2._dl*dwinlens(i)-dwinlens(i+1)+d*(dwinlens(i)+dwinlens(i+1) &
              +2._dl*(winlens(i)-winlens(i+1)))))
          end if
          vis(j2)=opac(j2)*expmmu(j2)
          dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
          ddvis(j2)=expmmu(j2)*(opac(j2)**3+3*opac(j2)*dopac(j2)+ddopac)
          else
          opac(j2)=dotmu(nthermo)
          dopac(j2)=ddotmu(nthermo)
          ddopac=dddotmu(nthermo)
          expmmu(j2)=emmu(nthermo)
          vis(j2)=opac(j2)*expmmu(j2)
          dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
          ddvis(j2)=expmmu(j2)*(opac(j2)**3+3._dl*opac(j2)*dopac(j2)+ddopac)

          end if
        end subroutine DoThermoSpline


        function ddamping_da(a)
          real(dl) :: ddamping_da
          real(dl), intent(in) :: a
          real(dl) :: R
          real(dl) :: dtauda
          external dtauda

          R=r_drag0*a
          !ignoring reionisation, not relevant for distance measures
          ddamping_da = (R**2 + 16*(1+R)/15)/(1+R)**2*dtauda(a)*a**2/(Recombination_xe(a)*akthom)

        end function ddamping_da


!!!!!!!!!!!!!!!!!!!
!JH: functions and subroutines for calculating z_star and z_drag

        function doptdepth_dz(z)
          real(dl) :: doptdepth_dz
          real(dl), intent(in) :: z
          real(dl) :: a
          real(dl) :: dtauda
          external dtauda

          a = 1._dl/(1._dl+z)

          !ignoring reionisation, not relevant for distance measures
          doptdepth_dz = Recombination_xe(a)*akthom*dtauda(a)

        end function doptdepth_dz

        function optdepth(z)
          real(dl) :: rombint2
          external rombint2
          real(dl) optdepth
          real(dl),intent(in) :: z

          optdepth = rombint2(doptdepth_dz, 0.d0, z, 1d-5, 20, 100)

        end function optdepth


        function ddragoptdepth_dz(z)
          real(dl) :: ddragoptdepth_dz
          real(dl), intent(in) :: z
          real(dl) :: a
          real(dl) :: dtauda
          external dtauda

          a = 1._dl/(1._dl+z)
          ddragoptdepth_dz = doptdepth_dz(z)/r_drag0/a

        end function ddragoptdepth_dz


        function dragoptdepth(z)
          real(dl) :: rombint2
          external rombint2
          real(dl) dragoptdepth
          real(dl),intent(in) :: z

          dragoptdepth =  rombint2(ddragoptdepth_dz, 0.d0, z, 1d-5, 20, 100)

        end function dragoptdepth


       subroutine find_z(func,zout)  !find redshift at which (photon/drag) optical depth = 1
          real(dl), external :: func
          real(dl), intent(out) :: zout
          real(dl) :: try1,try2,diff,avg
          integer :: i

          try1 = 0.d0
          try2 = 10000.d0

          i=0
          diff = 10.d0
         do while (diff .gt. 1d-3)
             i=i+1
             if (i .eq. 100) then
               call GlobalError('optical depth redshift finder did not converge',error_reionization)
               zout=0
               return
            end if

             diff = func(try2)-func(try1)
             avg = 0.5d0*(try2+try1)
             if (func(avg) .gt. 1.d0) then
                try2 = avg
             else
                try1 = avg
             end if
          end do

          zout = avg

        end subroutine find_z

!!!!!!!!!!!!!!!!!!! end JH

      end module ThermoData
