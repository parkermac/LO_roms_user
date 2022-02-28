      MODULE biology_mod
!
! **********************************************************************
! NPZDO ecosystem model, adapted from fennel.h. (The architecture of the other
! ROMS NPZD models is a little different.) Not part of the ROMS distribution.
!
! v1.1	updated for ROMS 3.6 by Neil
! v1.0  minor edits by Neil (NPZDO_BANAS_DEBUG, COAG_PHYTOS).
! This is the version used in the PNWTOX papers.
! v0.4  Large detritus pool & Coag added by Sam,
!         Optimal Uptake Kinetics added by Kristen, 
!		  cppdefs for optional light attenuation scheme
!         and for creating a bio BC out of Salish Sea
! v0.3  grazing and Z losses rewritten
!         some renaming (like iSDeN -> iDetr)
! v0.2  chlorophyll, NH4, and large detritus removed (with associated fluxes)
!         phytoplankton growth rewritten without temp dependence
!         diagnostic fluxes now saved divided by timestep
!         AttS included in light attenuation calculation
!         sediment-interaction parameterization removed; now detritus just
!         sinks out of the system
! v0.1  CARBON code removed
! **********************************************************************
!
!
      implicit none
!
      PRIVATE
      PUBLIC  :: biology
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE biology (ng,tile)
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=MyFile
!        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                   GRID(ng) % rmask_full,                         &
# endif

#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
#if defined OXYGEN|| defined CARBON

# ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
# else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
# endif

#endif
#ifdef CARBON
     &                   OCEAN(ng) % pH,                                &
!    &                   OCEAN(ng) % pHfull,                            &
#endif

#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                          &
     &                   DIAGS(ng) % DiaBio3d,                          &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif

      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                         rmask_full,                              &
# endif
#endif
     &                         Hz, z_r, z_w, srflx,                     &
#if defined OXYGEN|| defined CARBON

# ifdef BULK_FLUXES
     &                         Uwind, Vwind,                            &
# else
     &                         sustr, svstr,                            &
# endif
#endif
#ifdef CARBON
     &                         pH,                                      &
!    &                         pH,pHfull,                               &
#endif
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d, DiaBio3d,                      &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:,LBj:)
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# if defined OXYGEN || defined CARBON

#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  endif

# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:,LBj:)
!     real(r8), intent(inout):: pHfull(LBi:,LBj:,UBk:) 

# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else


# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:UBi,LBj:UBj)
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# if defined OXYGEN
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!



# ifdef CARBON
      integer, parameter :: Nsink = 3        ! only detritus sinks.
#else
      integer, parameter :: Nsink = 2        ! only detritus sinks.
                                             ! but the machinery is retained
                                             ! to let more things sink.
#endif
#ifdef OXYGEN
      real(r8) :: u10squ
#endif

      integer :: Iter, i, ibio, isink, itrc, ivar, j, k, ks

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

#ifdef OXYGEN
      real(r8), parameter :: OA0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: OA1 = 3.22014_r8       ! saturation
      real(r8), parameter :: OA2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: OA3 = 4.94457_r8
      real(r8), parameter :: OA4 =-0.256847_r8
      real(r8), parameter :: OA5 = 3.88767_r8
      real(r8), parameter :: OB0 =-0.00624523_r8
      real(r8), parameter :: OB1 =-0.00737614_r8
      real(r8), parameter :: OB2 =-0.0103410_r8
      real(r8), parameter :: OB3 =-0.00817083_r8
      real(r8), parameter :: OC0 =-0.000000488682_r8
      real(r8), parameter :: rOxNO3= 8.625_r8       ! 138/16
      real(r8) :: l2mol = 1000.0_r8/22.9316_r8      ! liter to mol
#endif
#ifdef CARBON
      integer :: iday, month, year

      integer, parameter :: DoNewton = 1            ! pCO2 solver

      real(r8), parameter :: Acoef = 2073.1_r8      ! Schmidt
      real(r8), parameter :: Bcoef = 125.62_r8      ! number
      real(r8), parameter :: Ccoef = 3.6276_r8      ! transfer
      real(r8), parameter :: Dcoef = 0.043219_r8    ! coefficients

      real(r8), parameter :: A1 = -60.2409_r8       ! surface
      real(r8), parameter :: A2 = 93.4517_r8        ! CO2
      real(r8), parameter :: A3 = 23.3585_r8        ! solubility
      real(r8), parameter :: B1 = 0.023517_r8       ! coefficients
      real(r8), parameter :: B2 = -0.023656_r8
      real(r8), parameter :: B3 = 0.0047036_r8

      real(r8) :: pmonth                         ! months since Jan 1951
      real(r8) :: pCO2air_secular
      real(r8) :: yday, hour

      real(r8), parameter :: pi2 = 6.2831853071796_r8

      real(r8), parameter :: D0 = 282.6_r8          ! coefficients
      real(r8), parameter :: D1 = 0.125_r8          ! to calculate
      real(r8), parameter :: D2 =-7.18_r8           ! secular trend in
      real(r8), parameter :: D3 = 0.86_r8           ! atmospheric pCO2
      real(r8), parameter :: D4 =-0.99_r8
      real(r8), parameter :: D5 = 0.28_r8
      real(r8), parameter :: D6 =-0.80_r8
      real(r8), parameter :: D7 = 0.06_r8

     real(r8):: fco2(N(ng)), co2(N(ng)),ppco2(N(ng))
     real(r8):: sili2(N(ng)),po4(N(ng))
     real(r8):: hco3(N(ng)), co3(N(ng)),lat(N(ng))
     real(r8) :: OmegaA(N(ng)), OmegaC(N(ng)), BetaD(N(ng)) 
     real(r8) :: rhoSW(N(ng)), p(N(ng)), tempis(N(ng))
!     REAL(r8) :: p80, sw_temp 
!     REAL(r8) :: sw_ptmp, sw_adtg
!     EXTERNAL  p80, sw_temp
!     EXTERNAL  sw_ptmp,sw_adtg
#endif
#ifdef CARBON
! Zooplankton Carbon:Nitrogen ratio [mole_C/mole_N], {5.0d0}.

        real(r8), parameter ::ZooCN = 6.625_r8
! Phytoplankton Carbon:Nitrogen ratio [mole_C/mole_N] , {6.625d0}.

       real(r8), parameter ::PhyCN = 6.625_r8
! Phytoplankton Carbon:Nitrogen ratio [mole_C/mole_N] , {6.625d0}.
    
        real(r8), parameter ::DetCN = 7.3125_r8
        real(r8), parameter ::CaS = 10280.0_r8 !10.27 mmol/kg*1026 kg/m3
        real(r8), parameter ::kdiss = 0.38_r8 !per day [0.01, 0.38 H&E, 1]


#endif

      real(r8) :: dtdays, mu, alphaE, PAR, Att

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5
      real(r8) :: fac1, fac2, fac3
      real(r8) :: cffL, cffR, cu, dltL, dltR

      real(r8) :: total_N, NO3loss

#ifdef OXYGEN
      real(r8) :: SchmidtN_Ox, O2satu, O2_Flux
      real(r8) :: TS, AA
#endif
#ifdef CARBON
      real(r8) :: CaCO3prod
      real(r8) :: CO2_Flux, CO2_sol, SchmidtN, TempK,Ksp
      real(r8) :: SqrtSalt, invTempk, logTempk,CaCO3diss,omega
#endif

      real(r8) :: N_Flux_Grazing       ! assim + egest + excret
      real(r8) :: N_Flux_PProd
      real(r8) :: N_Flux_Pmort, N_Flux_Zmort, N_Flux_Excret
      real(r8) :: N_Flux_Remin, N_Flux_Remin2
      real(r8) :: N_Flux_CoagD, N_Flux_CoagP
      real(r8) :: totalN, C_Flux_ReminL, C_Flux_ReminS


      real(r8), dimension(Nsink) :: Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur

#ifdef CARBON
      real(r8), dimension(IminS:ImaxS) :: pCO2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: CO3_test
#endif
#ifdef DIAGNOSTICS_BIO
      real(r8) :: fiter
#endif
 


      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

#include "set_bounds.h"
#ifdef DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                            &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO ivar=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,ivar)=0.0_r8
            END DO
          END DO
        END DO
        DO ivar=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,ivar)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)

#ifdef DIAGNOSTICS_BIO
 !
 !  A factor to account for the number of iterations in accumulating
 !  diagnostic rate variables.
 !
	fiter=1.0_r8/REAL(BioIter(ng),r8)
#endif

!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iDetr
      idsink(2)=iLDetr

      idsink(3)=iCaCO3

!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=detWsink(ng)
      Wbio(2)=LdetWsink(ng)


      Wbio(3)=LdetWsink(ng)

!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
       totalN = 0.0_r8
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
#ifdef NPZDO_BANAS_DEBUG
               IF (j.EQ.2.0_r8) THEN
		         IF (i.EQ.2.0_r8) THEN
			       IF (k.EQ.30) THEN
			         print *, ibio,'_beginning', Bio(i,k,ibio)
                                   IF (ibio.LT.7) THEN
					totalN = Bio(i,k,ibio)+totalN
				      	print *, 'TotalN', totalN
                                   END IF
                                END IF
            	         END IF
               END IF
#endif
            END DO
          END DO
        END DO

#ifdef CARBON
        DO k=1,N(ng)
          DO i=Istr,Iend
!            Bio_old(i,k,iTIC_)=MIN(Bio_old(i,k,iTIC_),3000.0_r8)
!            Bio_old(i,k,iTIC_)=MAX(Bio_old(i,k,iTIC_),400.0_r8)
!            Bio(i,k,iTIC_)=Bio_old(i,k,iTIC_)
          END DO
        END DO
#endif
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong s the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!
!  In the implicit algorithm, we have for example (N: nitrate,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            Att=0.0_r8
            IF (PARsur(i).gt.0.0_r8) THEN
              DO k=N(ng),1,-1
!
!  Attenuate the light to the center of the grid cell.
!  AttS is a salinity dependent correction to AttSW; AttSW is defined at
!  a reference salinity of 32.
!
                Att=EXP(-0.5_r8*(AttSW(ng)+                             &
     &                           AttP(ng)*Bio(i,k,iPhyt)+               &
     &                           AttS(ng)*(Bio(i,k,isalt)-32_r8)) *     &
     &                  (z_w(i,j,k)-z_w(i,j,k-1)))
                PAR=PAR*Att
!
!  Phytoplankton growth.
!
!
                alphaE = PAR * phyAlpha(ng)
                mu = phyMu0(ng) * alphaE /                              &
     &                 SQRT(phyMu0(ng)*phyMu0(ng) + alphaE*alphaE)
!
!  Clamping growth in the Salish Sea and Columbia River
#ifdef NO_GROW_SALISH
!  Turning off growth in Salish Sea (3 boxes which include all of Salish and east side of Strait to -123.5degE
		IF (j.gt.340_r8.AND.i.gt.46_r8) THEN
		      mu=0d0
		      Bio(i,k,iPhyt)=0.01d0
		      Bio(i,k,iZoop)=0.01d0
		ELSE IF (j.lt.340_r8.AND.j.gt.309_r8.AND.i.gt.106_r8) THEN
		      mu=0d0
		      Bio(i,k,iPhyt)=0.01d0
		      Bio(i,k,iZoop)=0.01d0
		ELSE IF (j.lt.309_r8.AND.j.gt.199_r8.AND.i.gt.129_r8) THEN
		      mu=0d0
		      Bio(i,k,iPhyt)=0.01d0
		      Bio(i,k,iZoop)=0.01d0
		END IF
!  Turning off growth in Columbia River (1 box)
		IF (i.gt.125_r8.AND.j.gt.96_r8.AND.j.lt.147_r8) THEN
			mu=0d0
			Bio(i,k,iPhyt)=0.01d0
			Bio(i,k,iZoop)=0.01d0
		END IF

#endif
!
!  Implement the Optimal Uptake Model for nutrient uptake (default is M-M kinetics)
#ifdef OPT_UPTAKE
                cff = mu * dtdays * Bio(i,k,iPhyt) /                    &
      &                 (phyKs(ng) + Bio(i,k,iNO3_)                     &
      &                  + 2.0_r8*SQRT(phyKs(ng)*Bio(i,k,iNO3_)))
#else
                cff = mu * dtdays * Bio(i,k,iPhyt) /                    &
      &                 (phyKs(ng) + Bio(i,k,iNO3_))
#endif				   
                Bio(i,k,iNO3_) = Bio(i,k,iNO3_) / (1.0_r8+cff)
                N_Flux_PProd = Bio(i,k,iNO3_) * cff
                Bio(i,k,iPhyt) = Bio(i,k,iPhyt) + N_Flux_PProd
!
#ifdef DIAGNOSTICS_BIO
                DiaBio3d(i,j,k,iPPro)=DiaBio3d(i,j,k,iPPro)+            &
      &                                N_Flux_PProd*fiter 
#endif
#ifdef OXYGEN
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg) + N_Flux_PProd*rOxNO3    
#endif
#ifdef CARBON
!
!  Total inorganic carbon (CO2) uptake during phytoplankton growth.
!
                cff1=PhyCN*(N_Flux_PProd)
                !CaCO3prod=cff1*0.05_r8 !PIC/POC ratio
                CaCO3prod=0.0_r8 !PIC/POC ratio



		Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-cff1
		Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-CaCO3prod 
! is this necessary? or is it already in there....

		Bio(i,k,iCaCO3)=Bio(i,k,iCaCO3)+CaCO3prod

!
!  Account for the uptake of N on total alkalinity.
                Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+N_Flux_PProd
		Bio(i,k,iTAlk)=Bio(i,k,iTAlk)-2.0_r8*CaCO3prod

#endif
!
!  Attenuate the light to the bottom of the grid cell.
!
                PAR=PAR*Att
#ifdef NPZDO_BANAS_DEBUG
		        IF (j.EQ.2.0_r8) THEN
		          IF (i.EQ.2.0_r8) THEN
			        IF (k.EQ.30) THEN
                 print *, 'N_Flux_PProd', N_Flux_PProd!/dtdays
                 print *, 'Oxy_Flux_PProd', N_Flux_PProd*rOxNO3!/dtdays
                 print *, 'DIC_Flux_PProd', N_Flux_PProd*PhyCN!/dtdays
                 print *, 'ALK_Flux_PProd', N_Flux_PProd!/dtdays
                 print *, 'CaCO3_Prod', CaCO3prod!/dtdays
                 print *, 'loc', i, j, k
                    END IF
                  END IF
                END IF
#endif
              END DO
            END IF
          END DO
!
!-----------------------------------------------------------------------
!  grazing and mortality
!-----------------------------------------------------------------------
!
          fac1=dtdays*zooI0(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
! Total grazing (quadratic)
!
              cff = fac1 * Bio(i,k,iZoop) * Bio(i,k,iPhyt) /            &
     &            (zooKs(ng)*zooKs(ng) + Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              Bio(i,k,iPhyt) = Bio(i,k,iPhyt) / (1.0_r8 + cff)
              N_Flux_Grazing = cff * Bio(i,k,iPhyt)
#ifdef DIAGNOSTICS_BIO
              DiaBio3d(i,j,k,iGraz) = DiaBio3d(i,j,k,iGraz) +           &
     &                                 N_Flux_Grazing*fiter 
#endif

!
! Partition total grazing among Z, D, and N
!
              Bio(i,k,iZoop) = Bio(i,k,iZoop) +                         &
     &            zooEps(ng) * N_Flux_Grazing
              Bio(i,k,iDetr) = Bio(i,k,iDetr) +                         &
     &            (1-zooEps(ng)) * zooFegest(ng) * N_Flux_Grazing
              N_Flux_Excret =                                           &
     &            (1-zooEps(ng)) * (1-zooFegest(ng)) * N_Flux_Grazing
              Bio(i,k,iNO3_) = Bio(i,k,iNO3_) + N_FLux_Excret

        
#ifdef OXYGEN
              Bio(i,k,iOxyg) = Bio(i,k,iOxyg) - rOxNO3 * N_Flux_Excret
!	      cff = N_Flux_Excret * dtdays * Bio(i,k,iOxyg) /        &
!     &                 (grzKs(ng) + Bio(i,k,iOxyg))				   
!              Bio(i,k,iOxyg) = Bio(i,k,iOxyg) / (1.0_r8+cff)
#endif
#ifdef CARBON
!              Bio(i,k,iCaCO3)=Bio(i,k,iCaCO3)+                          &
!     &                       0.05_r8*(ZooCN*(N_Flux_Grazing))

              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       ZooCN*(N_Flux_Excret)
!	      Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-                            &
!     &                       (0.05_r8*(ZooCN*(N_Flux_Grazing)))
!       	      Bio(i,k,iTAlk)=Bio(i,k,iTAlk)                             &
!     &                       - 2.0_r8*                                  &
!     &                       (0.05*ZooCN*(N_Flux_Grazing))
 	      Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+                            &
     &                       (N_Flux_Excret)
#endif

!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
              N_Flux_Pmort = dtdays * phyM(ng) *                        &
     &              MAX(Bio(i,k,iPhyt)-phyMin(ng),0.0_r8)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt) - N_Flux_Pmort
              Bio(i,k,iDetr)=Bio(i,k,iDetr) + N_Flux_Pmort
!
!  zooplankton mortality (quadratic, limited by a zooplankton minimum)
!
              N_Flux_Zmort = dtdays * zooZeta(ng) * Bio(i,k,iZoop) *    &
     &              MAX(Bio(i,k,iZoop)-zooMin(ng),0.0_r8)
              Bio(i,k,iZoop) = Bio(i,k,iZoop) - N_Flux_Zmort
              Bio(i,k,iDetr) = Bio(i,k,iDetr) + N_Flux_Zmort
#ifdef CARBON

!	      Bio(i,k,iCaCO3)=Bio(i,k,iCaCO3)+                          &
!     &                       0.05_r8*(DetCN*(N_Flux_Zmort))+            &
!     &                       0.05_r8*(DetCN*(N_Flux_Pmort))

! 	      Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-                            &
!     &                       0.05_r8*(DetCN*(N_Flux_Zmort))-            &
!     &                       0.05_r8*(DetCN*(N_Flux_Pmort))
!       	      Bio(i,k,iTAlk)=Bio(i,k,iTAlk)-                            &
!     &                       2.0_r8*(0.05_r8*DetCN*(N_Flux_Zmort))-     &
!     &                       2.0_r8*(0.05_r8*DetCN*(N_Flux_Pmort))
#endif
#ifdef NPZDO_BANAS_DEBUG
	      IF (j.EQ.2.0_r8) THEN
		     IF (i.EQ.2.0_r8) THEN
			IF (k.EQ.30.0_r8) THEN
        	         print *, 'N_Flux_Excret', N_Flux_Excret
       	                 print *, 'DIC_Flux_Excret', N_Flux_Excret*ZooCN
        	         print *, 'D_Flux_grazing', N_Flux_Grazing      &
     &                   *(1-zooEps(ng)) * (1-zooFegest(ng))
                         print *,'Flux_Graz',N_Flux_Grazing*0.05*ZooCN
			 print *, 'Z_Flx_Graz',N_Flux_Grazing*zooEps(ng)
			 print *, 'Flux_Pmort', 0.05*DetCN*N_Flux_Pmort
                         print *, 'Flux_Zmort', 0.05*DetCN*N_Flux_Zmort
		        END IF
		    END IF
	      END IF
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Coagulation of small Detritus (and maybe Phytoplankton) into Large Detritus
!-----------------------------------------------------------------------
!
          fac1=dtdays*CoagR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
#ifdef COAG_PHYTOS
              cff1=fac1*(Bio(i,k,iDetr)+Bio(i,k,iPhyt))
              cff2=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
              Bio(i,k,iDetr)=Bio(i,k,iDetr)*cff2
              N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
              Bio(i,k,iLDetr)=Bio(i,k,iLDetr)+N_Flux_CoagP
#else
              cff1=fac1*Bio(i,k,iDetr)
              cff2=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iDetr)=Bio(i,k,iDetr)*cff2
              N_Flux_CoagP=0.0_r8		
#endif
              N_Flux_CoagD=Bio(i,k,iDetr)*cff1
              Bio(i,k,iLDetr)=Bio(i,k,iLDetr)+N_Flux_CoagD
              
#ifdef NPZDO_BANAS_DEBUG
              IF (j.EQ.2.0_r8) THEN
		        IF (i.EQ.2.0_r8) THEN
			      IF (k.EQ.30.0_r8) THEN
			        print *, 'N_Flux_CoagP', N_Flux_CoagP
			        print *, 'N_Flux_CoagD', N_Flux_CoagD
                  END IF
		        END IF
  	          END IF
#endif
             END DO
           END DO
!

!
!-----------------------------------------------------------------------
!  Detritus remineralization
!-----------------------------------------------------------------------
!

          cff = dtdays * detRemin(ng)
          fac2 = 1.0_r8   ! used for O2 dependence
          DO k=1,N(ng)
            DO i=Istr,Iend
#ifdef OXYGEN
              ! saturation curve with threshhold for O2 dependence
!             fac1 = MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8)
!              fac2 = MAX(fac1/(3.0_r8+fac1),0.0_r8)
              
              Bio(i,k,iDetr) = Bio(i,k,iDetr) / (1.0_r8 + cff*fac2)
              Bio(i,k,iLDetr) = Bio(i,k,iLDetr) / (1.0_r8 + cff*fac2)


              N_Flux_Remin = Bio(i,k,iDetr) * cff*fac2
              N_Flux_Remin2 = Bio(i,k,iLDetr) * cff*fac2


	      IF(((N_Flux_Remin+N_Flux_Remin2)*rOxNO3).GT.Bio(i,k,iOxyg))THEN

                 Bio(i,k,iNO3_) = Bio(i,k,iNO3_) - N_Flux_Remin      &
     &                    - N_Flux_Remin2
	      ELSE
                Bio(i,k,iNO3_) = Bio(i,k,iNO3_) + N_Flux_Remin       &
     &                   + N_Flux_Remin2

	            
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg) - N_Flux_Remin * rOxNO3 
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg) - N_Flux_Remin2 * rOxNO3
	     ENDIF
#else
             
              Bio(i,k,iDetr) = Bio(i,k,iDetr) / (1.0_r8 + cff*fac2)
              Bio(i,k,iLDetr) = Bio(i,k,iLDetr) / (1.0_r8 + cff*fac2)
              N_Flux_Remin = Bio(i,k,iDetr) * cff*fac2
              N_Flux_Remin2 = Bio(i,k,iLDetr) * cff*fac2
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Remin+N_Flux_Remin2


#endif
        
          				     
#ifdef CARBON
!
!	  DO k=1,N(ng)
!            DO i=Istr,Iend
              C_Flux_ReminS= N_Flux_Remin*PhyCN 
              C_Flux_ReminL= N_Flux_Remin2*PhyCN
              Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                            &
     &                       C_Flux_ReminS+C_Flux_ReminL
	      Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+N_Flux_Remin+N_Flux_Remin2	
#ifdef NPZDO_BANAS_DEBUG
               IF (j.EQ.2.0_r8) THEN
		     IF (i.EQ.2.0_r8) THEN
                         IF (k.eq.30) THEN
                 	  print *, 'N_Flux_Remin', N_Flux_Remin !*rOxNO3
                 	  print *, 'N_Flux_Remin2', N_Flux_Remin2!*rOxNO3
			  print *, 'C_Flux_Remin', C_Flux_ReminS 
                 	  print *, 'C_Flux_Remin2', C_Flux_ReminL
			  print *, 'CaCO3diss', CaCO3diss
                         END IF
                    END IF
               END IF
#endif
            
           END DO
          END DO
#endif
#ifdef DIAGNOSTICS_BIO
            IF(((N_Flux_Remin+N_Flux_Remin2)*rOxNO3).GT.Bio(i,k,iOxyg)) THEN
	      DiaBio3d(i,j,k,iDNIT)=DiaBio3d(i,j,k,iDNIT)+          &
    &             (N_Flux_Remin+N_Flux_Remin2) *fiter
            ELSE            
              DiaBio3d(i,j,k,iRem)=DiaBio3d(i,j,k,iRem)+            &
    &             (N_Flux_Remin + N_Flux_Remin2) *fiter
            ENDIF
#endif
	     
#ifdef CARBON
!
!-----------------------------------------------------------------------
!  CaCO3 Dissolution
!-----------------------------------------------------------------------

!               DO k=1,N(ng)
!                DO i=Istr,Iend			  
!                sili2(k)=0.0_r8
!                po4(k)= 0.0_r8			  
!                 lat(k)= 48.0_r8
		 !pHfull(i,j,k)=8.0d0
!	    CO3_test(i,k)= (Bio(i,k,iTAlk)                             &
!     &                     -0.96_r8* Bio(i,k,iTIC_))/2.0_r8

!	if(Bio(i,k,iCaCO3).gt.0.0_r8) then
!	    if (CO3_test(i,k).lt.100.0d0) then

!       CALL vars(pHfull(i,j,:), ppco2, fco2, co2, hco3, co3,           &
!     &           OmegaA, OmegaC, BetaD, rhoSW, p, tempis,              &
#  ifdef MASKING
!     &                        rmask,                                    &
#  endif
!     &		Bio(i,:,itemp), Bio(i,:,isalt), Bio(i,:,iTAlk),        &
!     &           Bio(i,:,iTIC_), sili2(:), po4(:),                     &
!     &	        -z_w(i,j,:), lat, N(ng),                               &
!     &           "mol/m3","Tpot","m","l10", "m10", "dg")
	      

!!	    omega=(CaS*CO3(i,k))/(Ksp*(10.0_r8**12.0_r8)*               &
!!     &             (1026.0_r8/1000.0_r8))

	    
!	     IF (OmegaA(k).LT.1.0_r8) THEN
!	       CaCO3diss=Bio(i,k,iCaCO3)*kdiss*(1.0_r8-OmegaA(k))*dtdays !!!
!	       Bio(i,k,iCaCO3)=Bio(i,k,iCaCO3)-CaCO3diss
!	       Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+CaCO3diss                            
!	       Bio(i,k,iTAlk)=Bio(i,k,iTAlk)+2.0_r8*CaCO3diss

!            END IF
!	   ELSE
	       CaCO3diss = 0.d0
!           endif
!        endif

#ifdef NPZDO_BANAS_DEBUG
               IF (j.EQ.2.0_r8) THEN
		     IF (i.EQ.2.0_r8) THEN
                         IF (k.eq.30) THEN
			  print *, 'CaCO3diss', CaCO3diss
                         END IF
                    END IF
               END IF
#endif

!                END DO
!             END DO
#endif

#ifdef CARBON
!
!-----------------------------------------------------------------------
!  Surface CO2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute equilibrium partial pressure inorganic carbon (ppmv) at the
!  surface.
!
          k=N(ng)
# ifdef pCO2_RZ
          CALL pCO2_water_RZ (Istr, Iend, LBi, UBi, LBj, UBj,           &
     &                        IminS, ImaxS, j, DoNewton,                &
#  ifdef MASKING
     &                        rmask,                                    &
#  endif
     &                        Bio(IminS:,k,itemp), Bio(IminS:,k,isalt), &
     &                        Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk), &
     &                        pH, pCO2)
# else
          CALL pCO2_water (Istr, Iend, LBi, UBi, LBj, UBj,              &
     &                     IminS, ImaxS, j, DoNewton,                   &
#  ifdef MASKING
     &                     rmask,                                       &
#  endif
     &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
     &                     Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk),    &
     &                     0.0_r8, 0.0_r8, pH, pCO2)
	    
# endif
!      if(pCO2(i).lt.0.0_r8)then
!	print *, 'pco2', pCO2(i), 'TEMP',Bio(i,k,itemp),               &
!     &   Bio(i,k,isalt), 'ALK',Bio(i,k,iTAlk), 'TIC',Bio(i,k,iTIC_)
!      endif			
!
!  Compute surface CO2 gas exchange.
!
          cff1=rho0*550.0_r8
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
          DO i=Istr,Iend
!
!  Compute CO2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)**2+Vwind(i,j)**2
# else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
            SchmidtN=Acoef-                                             &
     &               Bio(i,k,itemp)*(Bcoef-                             &
     &                               Bio(i,k,itemp)*(Ccoef-             &
     &                               Bio(i,k,itemp)*Dcoef))
            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN)
!
!  Calculate CO2 solubility [mol/(kg.atm)] using Weiss (1974) formula.
!
            TempK=0.01_r8*(Bio(i,k,itemp)+273.15_r8)
            CO2_sol=EXP(A1+                                             &
     &                  A2/TempK+                                       &
     &                  A3*LOG(TempK)+                                  &
     &                  Bio(i,k,isalt)*(B1+TempK*(B2+B3*TempK)))
!
!  Add in CO2 gas exchange.
!
!            CALL caldate (r_date, tdays(ng), year, yday, month, iday,   &
!     &                    hour)
!            pmonth=2003.0_r8-1951.0_r8+yday/365.0_r8
!!          pCO2air_secular=D0+D1*pmonth*12.0_r8+                       &
!!   &                         D2*SIN(pi2*pmonth+D3)+                   &
!!   &                         D4*SIN(pi2*pmonth+D5)+                   &
!!   &                         D6*SIN(pi2*pmonth+D7)
!!          CO2_Flux=cff3*CO2_sol*(pCO2air_secular-pCO2(i))

	    if(pCO2(i).gt.0.0_r8)then
                  CO2_Flux=cff3*CO2_sol*(pCO2air(ng)-pCO2(i))
	    else
                  CO2_Flux = 0.0_r8
            endif
            Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+                              &
     &                     CO2_Flux*Hz_inv(i,k)
# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iCOfx)=DiaBio2d(i,j,iCOfx)+                    &
#  ifdef WET_DRY
     &                          rmask_full(i,j)*                          &
#  endif
     &                          CO2_Flux*fiter
            DiaBio2d(i,j,ipCO2)=pCO2(i)
#  ifdef WET_DRY
            DiaBio2d(i,j,ipCO2)=DiaBio2d(i,j,ipCO2)*rmask_full(i,j)
#  endif
# endif
        

#ifdef NPZDO_BANAS_DEBUG
               IF (j.EQ.2.0_r8) THEN
		     IF (i.EQ.2.0_r8) THEN
                        IF (k.eq.30) THEN
                 	  print *, 'pCO2flux', CO2_Flux*Hz_inv(i,k)
                         END IF
                    END IF
               END IF
#endif
     
     END DO
#endif
          

#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!
          cff1=rho0*550.0_r8
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
          k=N(ng)
          DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!

# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)


# endif

# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
            SchmidtN_Ox=1638.0_r8-                                      &
     &                  Bio(i,k,itemp)*(81.83_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (1.483_r8-                      &
     &                                   Bio(i,k,itemp)*0.008004_r8))
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
            SchmidtN_Ox=1953.4_r8-                                      &
     &                  Bio(i,k,itemp)*(128.0_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (3.9918_r8-                     &
     &                                   Bio(i,k,itemp)*0.050091_r8))
# endif

            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
            TS=LOG((298.15_r8-Bio(i,k,itemp))/                          &
     &             (273.15_r8+Bio(i,k,itemp)))
            AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+          &
     &             Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &             OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
            O2satu=l2mol*EXP(AA)
!
!  Add in O2 gas exchange.
!
            O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
     &                     O2_Flux*Hz_inv(i,k)
# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iO2fx)=DiaBio2d(i,j,iO2fx)+                    &
     &                          O2_Flux *fiter
# endif
#ifdef NPZDO_BANAS_DEBUG
               IF (j.EQ.2.0_r8) THEN
		     IF (i.EQ.2.0_r8) THEN
			IF (k.EQ.30.0_r8) THEN
 				   print *, 'O2_flux', O2_Flux*Hz_inv(i,k)
			END IF
                     END IF
                END IF
#endif
          END DO
#endif

!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                cffR=2.0_r8*dltR
                cffL=2.0_r8*dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO

#ifdef BIO_SEDIMENT
! **** insert a new sediment-interface parameterization here ****
! Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.

            cff2=4.0_r8/16.0_r8
# ifdef OXYGEN
            cff3=115.0_r8/16.0_r8
            cff4=108.0_r8/16.0_r8 !including oxid from ammonia to nitrate changes from 106 to 108
# endif
	    IF ((ibio.eq.iDetr).or.                                     &
     &          (ibio.eq.iLDetr)  )THEN
              DO i=Istr,Iend
		 cff1=(FC(i,0)*Hz_inv(i,1))
                 NO3loss=1.2_r8*dtdays*Hz_inv(i,1)
 
#  ifdef OXYGEN
		 IF(FC(i,0).GT.Bio(i,1,iOxyg))THEN
  		    
 	             Bio(i,1,iNO3_)=Bio(i,1,iNO3_)-cff1		     
                 ELSE
		     Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
 	             Bio(i,1,iNO3_)=Bio(i,1,iNO3_)+cff1
                 ENDIF
#  else                 
  		   Bio(i,1,iNO3_)=Bio(i,1,iNO3_)+cff1
#endif
		   IF(cff1.gt.NO3loss) THEN
		      Bio(i,1,iNO3_)=Bio(i,1,iNO3_)-NO3loss
	           END IF
!! the above addition is from unpublished results fuchsman et al (submitted)
!! from Al devol which show a constant flux of nitrate into the sediments at 
!!all depths down to 1000 m or so from oregon. 1.2 mmol NO3/m2 day from figure 4
# ifdef DIAGNOSTICS_BIO
		   IF(FC(i,0).GT.Bio(i,1,iOxyg))THEN
		      DiaBio3d(i,j,1,iDNIT)=DiaBio3d(i,j,1,iDNIT)       &
     &                                       +cff1*fiter  
                   ELSE
                      DiaBio2d(i,j,iDsed)=DiaBio2d(i,j,iDsed)+          &
     &                     cff1*cff4*Hz(i,j,1)*fiter 
                   END IF
# endif
# endif
              END DO
            END IF

# ifdef CARBON
!! STILL need to add denit effect on ALK
#ifdef CARBON

	    IF  (ibio.eq.iCaCO3) THEN 

              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iTAlk)=Bio(i,1,iTAlk)+2*cff1
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+cff1
              END DO
            END IF

             IF ((ibio.eq.iDetr).or.                                     &
     &          (ibio.eq.iLDetr)  )THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iTIC_)=Bio(i,1,iTIC_)+cff1*DetCN
              END DO
            END IF
# endif
#endif


          END DO SINK_LOOP
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of N and
!  C even when advection causes tracer concentration to go negative.
!  (J. Wilkin and H. Arango, Apr 27, 2012)
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
#ifdef NPZDO_BANAS_DEBUG
		      IF (j.EQ.2.0_r8) THEN
		        IF (i.EQ.2.0_r8) THEN
			      IF (k.EQ.30) THEN
			        print *, 'N_end', Bio(i,k,iNO3_)
			        print *, 'P_end', Bio(i,k,iPhyt)
			        print *, 'Z_end', Bio(i,k,iZoop)
                    print *, 'D_end', Bio(i,k,iDetr)
                    print *, 'LD_end', Bio(i,k,iLDetr)
			        print *, 'O_end', Bio(i,k,iOxyg)
				print *, 'DIC_end', Bio(i,k,iTIC_)
				print *, 'ALK_end', Bio(i,k,iTAlk)
				print *, 'CaCO3_end', Bio(i,k,iCaCO3)
                  END IF
                END IF
              END IF
#endif
            END DO
          END DO
        END DO
      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile


# ifdef pCO2_RZ
      SUBROUTINE pCO2_water_RZ (Istr, Iend,                             &
     &                          LBi, UBi, LBj, UBj, IminS, ImaxS,       &
     &                          j, DoNewton,                            &
#  ifdef MASKING
     &                          rmask,                                  &
#  endif
     &                          T, S, TIC, TAlk, pH, pCO2)
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4=0,            !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2= ppmv  (DoNewton=0)                              !
!                pCO2= ppmv  (DoNewton=1)                              !
!                                                                      !
!  This subroutine was adapted by Katja Fennel (Nov 2005) from         !
!  Zeebe and Wolf-Gladrow (2001).                                      !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Zeebe, R.E. and D. Wolf-Gladrow,  2005:  CO2 in Seawater:         !
!      Equilibrium, kinetics, isotopes, Elsevier Oceanographic         !
!      Series, 65, pp 346.                                             !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: TIC(IminS:)
      real(r8), intent(in) :: TAlk(IminS:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
!      real(r8), intent(inout) :: CO3(LBi:,LBj:)

#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
!      real(r8), intent(inout) :: CO3(LBi:UBi,LBj:UBj)
#  endif

      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30

      integer :: Hstep, Ibrack, Inewton, i

      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: scl, sqrtS
      real(r8) :: borate, alk, dic
      real(r8) :: ff, K1, K2, K12, Kb, Kw
      real(r8) :: p5, p4, p3, p2, p1, p0
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid
      real(r8) :: CO2star, Htotal, Htotal2
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif

        Tk=T(i)+273.15_r8
        centiTk=0.01_r8*Tk
        invTk=1.0_r8/Tk
        logTk=LOG(Tk)
        sqrtS=SQRT(S(i))
        scl=S(i)/1.80655_r8

        alk= TAlk(i)*0.000001_r8
        dic = TIC(i)*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
        ff=EXP(-162.8301_r8+                                            &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carboinic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  From Millero (1995; page 664) using Mehrbach et al. (1973) data on
!  seawater scale.
!-----------------------------------------------------------------------
!
        K1=10.0_r8**(62.008_r8-                                         &
     &               invTk*3670.7_r8-                                   &
     &               logTk*9.7944_r8+                                   &
     &               S(i)*(0.0118_r8-                                   &
     &                     S(i)*0.000116_r8))
        K2=10.0_r8**(-4.777_r8-                                         &
     &               invTk*1394.7_r8+                                   &
     &               S(i)*(0.0184_r8-                                   &
     &                     S(i)*0.000118_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
        Kb=EXP(-invTk*(8966.90_r8+                                      &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
        Kw=EXP(148.9652_r8-                                             &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974).
!-----------------------------------------------------------------------
!
        borate=0.000232_r8*scl/10.811_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Solve for h in fifth-order polynomial. First calculate
!  polynomial coefficients.
!
        K12 = K1*K2

        p5 = -1.0_r8;
        p4 = -alk-Kb-K1;
        p3 = dic*K1-alk*(Kb+K1)+Kb*borate+Kw-Kb*K1-K12
        p2 = dic*(Kb*K1+2*K12)-alk*(Kb*K1+K12)+Kb*borate*K1             &
     &       +(Kw*Kb+Kw*K1-Kb*K12)
        p1 = 2.0_r8*dic*Kb*K12-alk*Kb*K12+Kb*borate*K12                 &
     &       +Kw*Kb*K1+Kw*K12
        p0 = Kw*Kb*K12;
!
!  Set first guess and brackets for [H+] solvers.
!
        pH_guess=pH(i,j)         ! Newton-Raphson
        pH_hi=10.0_r8            ! high bracket/bisection
        pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
        X_guess=10.0_r8**(-pH_guess)
        X_lo=10.0_r8**(-pH_hi)
        X_hi=10.0_r8**(-pH_lo)
        X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
        IF (DoNewton.eq.1) THEN
          X=X_guess
!
          DO Inewton=1,InewtonMax
!
!  Evaluate f([H+]) = p5*x^5+...+p1*x+p0
!
            fn=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
!
!  Evaluate derivative, df([H+])/dx:
!
!     df= d(fn)/d(X)
!
            df=(((5*p5*X+4*p4)*X+3*p3)*X+2*p2)*X+p1
!
!  Evaluate increment in [H+].
!
            deltaX=-fn/df
!
!  Update estimate of [H+].
!
            X=X+deltaX
          END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
        ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
          BRACK_IT: DO Ibrack=1,IbrackMax
            DO Hstep=1,3
              IF (Hstep.eq.1) X=X_hi
              IF (Hstep.eq.2) X=X_lo
              IF (Hstep.eq.3) X=X_mid
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
              fni(Hstep)=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
            END DO
!
!  Now, bracket solution within two of three.
!
            IF (fni(3).eq.0) THEN
               EXIT BRACK_IT
            ELSE
               ftest=fni(1)/fni(3)
               IF (ftest.gt.0) THEN
                 X_hi=X_mid
               ELSE
                 X_lo=X_mid
               END IF
               X_mid=0.5_r8*(X_lo+X_hi)
            END IF
          END DO BRACK_IT
!
! Last iteration gives value.
!
          X=X_mid
        END IF
!-----------------------------------------------------------------------
!  Determine [CO3--].
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO3--] (mole/m3) as defined in Zeebe, chapter 1
!
!	CO3(i,j)=dic/(1+(Htotal/K2)+(Htotal2/(K1*K2)))!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
        CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
        pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
        pCO2(i)=CO2star*1000000.0_r8/ff

#  ifdef MASKING
      ELSE
        pH(i,j)=0.0_r8
        pCO2(i)=0.0_r8
      END IF
#  endif

      END DO I_LOOP

      RETURN
      END SUBROUTINE pCO2_water_RZ
# else
     SUBROUTINE pCO2_water (Istr, Iend,                                &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, j, DoNewton,                 &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
     &                       T, S, TIC, TAlk, PO4, SiO3, pH, pCO2)
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     PO4        Inorganic phosphate (millimol/m3).                    !
!     SiO3       Inorganic silicate (millimol/m3).                     !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4=0,            !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2=0.35074945E+03 ppmv  (DoNewton=0)                !
!                pCO2=0.35073560E+03 ppmv  (DoNewton=1)                !
!                                                                      !
!  This subroutine was adapted by Mick Follows (Oct 1999) from OCMIP2  !
!  code CO2CALC. Modified for ROMS by Hernan Arango (Nov 2003).        !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: TIC(IminS:)
      real(r8), intent(in) :: TAlk(IminS:)
      real(r8), intent(inout) :: pH(LBi:,LBj:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: PO4
      real(r8), intent(in) :: SiO3

      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30

      integer :: Hstep, Ibrack, Inewton, i

      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: SO4, scl, sqrtS, sqrtSO4
      real(r8) :: alk, dic, phos, sili
      real(r8) :: borate, sulfate, fluoride
      real(r8) :: ff, K1, K2, K1p, K2p, K3p, Kb, Kf, Ks, Ksi, Kw
      real(r8) :: K12, K12p, K123p, invKb, invKs, invKsi
      real(r8) :: A, A2, B, B2, C, dA, dB
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid,St
      real(r8) :: CO2star, Htotal, Htotal2,is,sqrtis,is2
      real(r8) :: pK1o,ma1,mb1,mc1,pK1,pK2o,ma2,mb2,mc2,pK2
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
      if (TIC(i).gt.100.0_r8.AND.TAlk(i).gt.100.0_r8) then
           
			    
        Tk=T(i)+273.15_r8
        centiTk=0.01_r8*Tk
        invTk=1.0_r8/Tk
        logTk=LOG(Tk)
        sqrtS=SQRT(S(i))
        SO4=19.924_r8*S(i)/(1000.0_r8-1.005_r8*S(i))
        sqrtSO4=SQRT(SO4)
        scl=S(i)/1.80655_r8
	St=0.02824_r8*(S(i)/35.0_r8)
!       Ionic strength:
	is = 19.924d0*S(i)/(1000.0d0 - 1.005d0*S(i))
        is2 = is* is
        sqrtis = SQRT(is)

        alk=TAlk(i)*0.000001_r8
        dic=TIC(i)*0.000001_r8
        phos=PO4*0.000001_r8
        sili=SiO3*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
        ff=EXP(-162.8301_r8+                                            &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carboinic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  From Millero (1995; page 664) using Mehrbach et al. (1973) data on
!  seawater scale.
!  Updated with Millero(2010) for freshwater. on total scale
!
!
!-----------------------------------------------------------------------
!
!        K1=10.0_r8**(62.008_r8-                                         &
!     &               invTk*3670.7_r8-                                   &
!     &               logTk*9.7944_r8+                                   &
!     &               S(i)*(0.0118_r8-                                   &
!     &                     S(i)*0.000116_r8))
!        K2=10.0_r8**(-4.777_r8-                                         &
!     &               invTk*1394.7_r8+                                   &
!     &               S(i)*(0.0184_r8-                                   &
!     &                     S(i)*0.000118_r8))

        pK1o = 6320.813d0*invTk + 19.568224_r8*logTk -126.34048_r8
	ma1 = 13.4051_r8*sqrtS + 0.03185_r8*S(i) - (5.218e-5)*(S(i)*S(i))
	mb1 = -531.095_r8*sqrtS - 5.7789_r8*S(i)
        mc1 = -2.0663_r8*sqrtS
        pK1 = pK1o + ma1 + mb1*invTk + mc1*logTk
        K1 = 10.0_r8**(-pK1) 

        pK2o = 5143.692_r8*invTk + 14.613358_r8*logTk -90.18333_r8
	ma2 = 21.5724_r8*sqrtS + 0.1212_r8*S(i) - (3.714e-4)*(S(i)*S(i))
	mb2 = -798.292_r8*sqrtS - 18.951_r8*S(i)
        mc2 = -3.403_r8*sqrtS
        pK2 = pK2o + ma2 + mb2*invTk + mc2*logTk
        K2 = 10.0_r8**(-pK2)


!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
        Kb=EXP(-invTk*(8966.90_r8+                                      &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute first (K1p), second (K2p), and third (K3p) dissociation
!  constant of phosphoric acid:
!
!           K1p = [H][H2PO4]/[H3PO4]
!           K2p = [H][HPO4]/[H2PO4]
!           K3p = [H][PO4]/[HPO4]
!
!  From DOE (1994) equations 7.2.20, 7.2.23, and 7.2.26, respectively.
!  With footnote using data from Millero (1974).
!-----------------------------------------------------------------------
!
        K1p=EXP(115.525_r8-                                             &
     &          invTk*4576.752_r8-                                      &
     &          logTk*18.453_r8+                                        &
     &          sqrtS*(0.69171_r8-invTk*106.736_r8)-                    &
     &          S(i)*(0.01844_r8+invTk*0.65643_r8))
        K2p=EXP(172.0883_r8-                                            &
     &          invTk*8814.715_r8-                                      &
     &          logTk*27.927_r8+                                        &
     &          sqrtS*(1.3566_r8-invTk*160.340_r8)-                     &
     &          S(i)*(0.05778_r8-invTk*0.37335_r8))
        K3p=EXP(-18.141_r8-                                             &
     &          invTk*3070.75_r8+                                       &
     &          sqrtS*(2.81197_r8+invTk*17.27039_r8)-                   &
     &          S(i)*(0.09984_r8+invTk*44.99486_r8))
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of silica, Ksi=[H][SiO(OH)3]/[Si(OH)4].
!  From Millero (1995; page 671) using data from Yao and Millero (1995).
!-----------------------------------------------------------------------
!
        Ksi=EXP(117.385_r8-                                             &
     &          invTk*8904.2_r8-                                        &
     &          logTk*19.334_r8+                                        &
     &          sqrtSO4*(3.5913_r8-invTk*458.79_r8)-                    &
     &          SO4*(1.5998_r8-invTk*188.74_r8-                         &
     &               SO4*(0.07871_r8-invTk*12.1652_r8))+                &
     &          LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
        Kw=EXP(148.9652_r8-                                             &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!------------------------------------------------------------------------
!  Compute salinity constant of hydrogen sulfate, Ks = [H][SO4]/[HSO4].
!  From Dickson (1990, J. chem. Thermodynamics 22, 113)
!------------------------------------------------------------------------
!
!       Ks = [H][SO4]/[HSO4]
!       (free scale)
!       Dickson (1990, J. chem. Thermodynamics 22, 113)
        Ks = EXP(-4276.1_r8*invTk + 141.328_r8 - 23.093_r8*logTk         &
     &           + (-13856.0_r8*invTk + 324.57_r8 - 47.986_r8*logTk)     &
     &           *sqrtis                                                 &
     &           + (35474.0_r8*invTk - 771.54_r8 + 114.723_r8*logTk) * is&
     &           - 2698.0_r8*invTk*is**1.5_r8                            &
     &           + 1776.0_r8*invTk*is2                                   & 
     &           + LOG(1.0_r8 - 0.001005_r8*S(i)))
!!ORIGINAL, no ionic strength terms
!        Ks=EXP(141.328_r8-                                              &
!     &         invTk*4276.1_r8-                                         &
!     &         logTk*23.093_r8+                                         &
!     &         sqrtSO4*(324.57_r8-invTk*13856.0_r8-logTk*47.986_r8-     &
!     &                  SO4*invTk*2698.0_r8)-                           &
!     &         SO4*(771.54_r8-invTk*35474.0_r8-logTk*114.723_r8-        &
!     &              SO4*invTk*1776.0_r8)+                               &
!     &         LOG(1.0_r8-0.001005_r8*S(i)))
!
!-----------------------------------------------------------------------
!  Compute stability constant of hydrogen fluorid, Kf = [H][F]/[HF].
!  From Dickson and Riley (1979) -- change pH scale to total.
!-----------------------------------------------------------------------
!The following is from the vars subroutine below for 'dg' option
         Kf = EXP(1590.2_r8*invTk - 12.641_r8 + 1.525_r8*sqrtis +        &
      &             LOG(1.0_r8 - 0.001005_r8*S(i)) +                     &
      &             LOG(1.0_r8 + St/Ks))
!! ORIGINAL
!        Kf=EXP(-12.641_r8+                                              &
!     &         invTk*1590.2_r8+                                         &
!     &         sqrtSO4*1.525_r8+                                        &
!     &         LOG(1.0_r8-0.001005_r8*S(i))+                            &
!     &         LOG(1.0_r8+0.1400_r8*scl/(96.062_r8*Ks)))
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974), sulfate (Morris
! and Riley, 1966), and fluoride (Riley, 1965).
!-----------------------------------------------------------------------
!
        borate=0.000232_r8*scl/10.811_r8
        sulfate=0.14_r8*scl/96.062_r8
        fluoride=0.000067_r8*scl/18.9984_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Set first guess and brackets for [H+] solvers.
!
        pH_guess=pH(i,j)         ! Newton-Raphson
        pH_hi=10.0_r8            ! high bracket/bisection
        pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
        X_guess=10.0_r8**(-pH_guess)
        X_lo=10.0_r8**(-pH_hi)
        X_hi=10.0_r8**(-pH_lo)
        X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
        IF (DoNewton.eq.1) THEN

          X=X_guess
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          DO Inewton=1,InewtonMax
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
            X2=X*X
            X3=X2*X
            invX=1.0_r8/X
            invX2=1.0_r8/X2

            A=X*(K12p+X*(K1p+X))
            B=X*(K1+X)+K12
            C=1.0_r8/(1.0_r8+sulfate*invKs)

            A2=A*A
            B2=B*B
            dA=X*(2.0_r8*K1p+3.0_r8*X)+K12p
            dB=2.0_r8*X+K1
!
!  Evaluate f([H+]):
!
!     fn=HCO3+CO3+borate+OH+HPO4+2*PO4+H3PO4+silicate+Hfree+HSO4+HF-TALK
!
            fn=dic*K1*(X+2.0_r8*K2)/B+                                  &
     &         borate/(1.0_r8+X*invKb)+                                 &
     &         Kw*invX+                                                 &
     &         phos*(K12p*X+2.0_r8*K123p-X3)/A+                         &
     &         sili/(1.0_r8+X*invKsi)-                                  &
     &         X*C-                                                     &
     &         sulfate/(1.0_r8+Ks*invX*C)-                              &
     &         fluoride/(1.0_r8+Kf*invX)-                               &
     &         alk
!
!  Evaluate derivative, f(prime)([H+]):
!
!     df= d(fn)/d(X)
!
            df=dic*K1*(B-dB*(X+2.0_r8*K2))/B2-                          &
     &         borate/(invKb*(1.0+X*invKb)**2)-                         &
     &         Kw*invX2+                                                &
     &         phos*(A*(K12p-3.0_r8*X2)-dA*(K12p*X+2.0_r8*K123p-X3))/A2-&
     &         sili/(invKsi*(1.0_r8+X*invKsi)**2)+                      &
     &         C+                                                       &
     &         sulfate*Ks*C*invX2/((1.0_r8+Ks*invX*C)**2)+              &
     &         fluoride*Kf*invX2/((1.0_r8+Kf*invX)**2)
!
!  Evaluate increment in [H+].
!
            deltaX=-fn/df
!
!  Update estimate of [H+].
!
            X=X+deltaX
          END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
        ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
          K12=K1*K2
          K12p=K1p*K2p
          K123p=K12p*K3p
          invKb=1.0_r8/Kb
          invKs=1.0_r8/Ks
          invKsi=1.0_r8/Ksi
!
          BRACK_IT: DO Ibrack=1,IbrackMax
            DO Hstep=1,3
              IF (Hstep.eq.1) X=X_hi
              IF (Hstep.eq.2) X=X_lo
              IF (Hstep.eq.3) X=X_mid
!
!  Set some common combinations of parameters used in the iterative [H+]
!  solver.
!
              X2=X*X
              X3=X2*X
              invX=1.0_r8/X

              A=X*(K12p+X*(K1p+X))+K123p
              B=X*(K1+X)+K12
              C=1.0_r8/(1.0_r8+sulfate*invKs)

              A2=A*A
              B2=B*B
              dA=X*(K1p*2.0_r8+3.0_r8*X2)+K12p
              dB=2.0_r8*X+K1
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
              fni(Hstep)=dic*(K1*X+2.0_r8*K12)/B+                       &
     &                   borate/(1.0_r8+X*invKb)+                       &
     &                   Kw*invX+                                       &
     &                   phos*(K12p*X+2.0_r8*K123p-X3)/A+               &
     &                   sili/(1.0_r8+X*invKsi)-                        &
     &                   X*C-                                           &
     &                   sulfate/(1.0_r8+Ks*invX*C)-                    &
     &                   fluoride/(1.0_r8+Kf*invX)-                     &
     &                   alk
            END DO
!
!  Now, bracket solution within two of three.
!
            IF (fni(3).eq.0.0_r8) THEN
              EXIT BRACK_IT
            ELSE
              ftest=fni(1)/fni(3)
              IF (ftest.gt.0.0) THEN
                X_hi=X_mid
              ELSE
                X_lo=X_mid
              END IF
              X_mid=0.5_r8*(X_lo+X_hi)
            END IF
          END DO BRACK_IT
!
! Last iteration gives value.
!
          X=X_mid
        END IF
!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
        Htotal=X
        Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
        CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
        pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
        pCO2(i)=CO2star*1000000.0_r8/ff

				    

#  ifdef MASKING
      ELSE
        pH(i,j)=0.0_r8
        pCO2(i)=0.0_r8
      END IF
#  endif
	
  ELSE
        pH(i,j)=0.0_r8
        pCO2(i)=0.0_r8
      END IF

      END DO I_LOOP

      RETURN
      END SUBROUTINE pCO2_water
# endif
      END MODULE biology_mod

