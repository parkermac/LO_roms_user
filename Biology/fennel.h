      SUBROUTINE biology (ng,tile)
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
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
        BIONAME(iNLM)=__FILE__
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
     &                   OCEAN(ng) % pHfull,                            &
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
     &                         pH,pHfull,                               &
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
     real(r8), intent(inout):: pHfull(LBi:,LBj:,UBk:) 

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

     SUBROUTINE vars(ph, ppco2, fco2, co2, hco3, co3,      &
     &           OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
     &           tempot, sal, alk, dic, sil, phos,         &
     &           depth, lat, N,                            &
     &           optRHO, optT, optP,                       &
     &           optB, optK1K2, optKf)
				     
!   Purpose:
!     Computes other standard carbonate system variables (pH, CO2*, HCO3- and CO32-, OmegaA, OmegaC, R)
!     as 1D arrays
!     FROM:
!     temperature, salinity, pressure,
!     total alkalinity (ALK), dissolved inorganic carbon (DIC),
!     silica and phosphate concentrations (all 1-D arrays)

!     INPUT variables:
!     ================
!     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
!             = pressure [db] (with optP='db')
!     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
!             = dummy array (unused when optP='db')
!     tempot  = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not temp)
!             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
!     sal     = salinity in [psu]
!     alk     = total alkalinity in [eq/m^3] with optRHO = 'mol/m3'
!             =               [eq/kg]  with optRHO = 'mol/kg'
!     dic     = dissolved inorganic carbon [mol/m^3] with optRHO = 'mol/m3'
!             =                            [mol/kg]  with optRHO = 'mol/kg'
!     sil     = silica    [mol/m^3] with optRHO = 'mol/m3'
!             =           [mol/kg]  with optRHO = 'mol/kg'
!     phos    = phosphate [mol/m^3] with optRHO = 'mol/m3'
!             =           [mol/kg]  with optRHO = 'mol/kg'

!     optRHO: choose input concentration units - mol/kg (data) vs. mol/m^3 (models)
!     -----------
!       -> 'mol/kg' for DIC and ALK given on mokal scale, i.e., in mol/kg  (std DATA units)
!       -> 'mol/m3' for DIC and ALK given in mol/m^3 (std MODEL units)

!     optT: choose in situ vs. potential temperature as input
!     ---------
!     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
!       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
!       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)

!     optP: choose depth (m) vs pressure (db) as input
!     ---------
!       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
!       -> 'db' means "depth" input is already in situ pressure [db], not m (thus p = depth)

!     optB: choose total boron formulation - Uppström (1974) vs. Lee et al. (2010)
!     ---------
!       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
!       -> 'l10' means use newer   formulation of Lee et al. (2010) for total Boron

!     optK1K2:
!     ---------
!       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
!                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
!       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
!                **** Valid for 0 < T < 50°C and 1 < S < 50 psu

!     optKf:
!     ----------
!       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
!               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
!       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)

!     OUTPUT variables:
!     =================
!     ph   = pH on total scale
!     pco2 = CO2 partial pressure (uatm)
!     fco2 = CO2 fugacity (uatm)
!     co2  = aqueous CO2 concentration [mol/m^3]
!     hco3 = bicarbonate (HCO3-) concentration [mol/m^3]
!     co3  = carbonate (CO3--) concentration [mol/m^3] 
!     OmegaA = Omega for aragonite, i.e., the aragonite saturation state
!     OmegaC = Omega for calcite, i.e., the   calcite saturation state
!     BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
!     rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
!     p = pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
!     tempis  = in-situ temperature [degrees C]

   USE mod_kinds
  IMPLICIT NONE

! Input variables
!>     number of records
  INTEGER, INTENT(in) :: N
!> either <b>in situ temperature</b> (when optT='Tinsitu', typical data) 
!! OR <b>potential temperature</b> (when optT='Tpot', typical models) <b>[degree C]</b>
  REAL(R8), INTENT(in),    DIMENSION(N):: tempot
!> salinity <b>[psu]</b>
  REAL(R8), INTENT(in), DIMENSION(N) :: sal
!> total alkalinity in <b>[eq/m^3]</b> (when optRHO = 'mol/m3') OR in <b>[eq/kg]</b>  (when optRHO = 'mol/kg')
  REAL(r8), INTENT(in), DIMENSION(N) :: alk
!> dissolved inorganic carbon in <b>[mol/m^3]</b> (when optRHO = 'mol/m3') OR in <b>[mol/kg]</b> (when optRHO = 'mol/!kg')
  REAL(r8), INTENT(in), DIMENSION(N) :: dic
!> SiO2 concentration in <b>[mol/m^3]</b> (when optRHO = 'mol/m3') OR in <b>[mol/kg]</b> (when optRHO = 'mol/kg')
  REAL(r8), INTENT(in), DIMENSION(N) :: sil
!> phosphate concentration in <b>[mol/m^3]</b> (when optRHO = 'mol/m3') OR in <b>[mol/kg]</b> (when optRHO = 'mol/kg'!)
  REAL(r8), INTENT(in), DIMENSION(N) :: phos
!f2py optional , depend(sal) :: n=len(sal)
!> depth in \b meters (when optP='m') or \b decibars (when optP='db')
  REAL(r8), INTENT(in),    DIMENSION(N) :: depth
!> latitude <b>[degrees north]</b>
  REAL(r8), INTENT(in),    DIMENSION(N) :: lat


  !> choose either \b 'mol/kg' (std DATA units) or \b 'mol/m3' (std MODEL units) to select 
  !! concentration units for input (for alk, dic, sil, phos) & output (co2, hco3, co3)
  CHARACTER(*), INTENT(in) :: optRHO
  !> choose \b 'Tinsitu' for in situ temperature or \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
  CHARACTER(*), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(*), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20), 
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
  CHARACTER(*), INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
  CHARACTER(*), INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
  CHARACTER(*), INTENT(in) :: optK1K2

! Output variables:
  !> pH on the <b>total scale</b>
  REAL(r8), INTENT(out), DIMENSION(N) :: ph
  !> CO2 partial pressure <b>[uatm]</b>
  REAL(r8), INTENT(out), DIMENSION(N) :: ppco2
  !> CO2 fugacity <b>[uatm]</b>
  REAL(r8), INTENT(out), DIMENSION(N) :: fco2
  !> aqueous CO2* concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg</b>] depending on choice for optRHO
  REAL(r8), INTENT(out), DIMENSION(N) :: co2
  !> bicarbonate ion (HCO3-) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optRHO
  REAL(r8), INTENT(out), DIMENSION(N) :: hco3
  !> carbonate ion (CO3--) concentration, either in <b>[mol/m^3]</b> or <b>[mol/kg]</b> depending on choice for optRHO
  REAL(r8), INTENT(out), DIMENSION(N) :: co3
  !> Omega for aragonite, i.e., the aragonite saturation state
  REAL(r8), INTENT(out), DIMENSION(N) :: OmegaA
  !> Omega for calcite, i.e., the calcite saturation state
  REAL(r8), INTENT(out), DIMENSION(N) :: OmegaC
  !> Revelle factor, i.e., dpCO2/pCO2 / dDIC/DIC
  REAL(r8), INTENT(out), DIMENSION(N) :: BetaD
  !> in-situ density of seawater; rhoSW = f(s, t, p) in <b>[kg/m3]</b>
  REAL(r8), INTENT(out), DIMENSION(N) :: rhoSW
  !> pressure <b>[decibars]</b>; p = f(depth, latitude) if computed from depth [m] (when optP='m') OR p = depth [db] (when optP='db')
  REAL(r8), INTENT(out), DIMENSION(N) :: p
  !> in-situ temperature \b <b>[degrees C]</b>
  REAL(r8), INTENT(out), DIMENSION(N) :: tempis

! Local variables
! REAL(kind=r4) :: tempot
  REAL(r8) :: ssal, salk, sdic, ssil, sphos

  REAL(r8) :: tempis68, tempot68
  REAL(r8) :: dfCO2, dpCO2
  REAL(r8) :: drho
  REAL(r8) :: invtk,dlogtk

  REAL(r8) :: Kh, K1, K2, Kb, Kw, Ks, Kf, Kspc
  REAL(r8) :: Kspa, K1p, K2p, K3p, Ksi
  REAL(r8) :: St, Ft, Bt

  REAL(r8), DIMENSION(1) :: aKh, aK1, aK2, aKb, aKw, aKs, aKf, aKspc
  REAL(r8), DIMENSION(1) :: aKspa, aK1p, aK2p, aK3p, aKsi
  REAL(r8), DIMENSION(1) :: aSt, aFt, aBt

  REAL(r8) :: DD,A,B,C,PhiD

  INTEGER :: i, icount

  REAL(r8) :: t, tk, prb
  REAL(r8) :: s
  REAL(r8) :: tc, ta
  REAL(r8) :: sit, pt
  REAL(r8) :: ah1

  REAL(r8) :: HSO4, HF, HSI, HPO4
  REAL(r8) :: ab, aw, ac, ah2, erel

  REAL(r8) :: cu, cb, cc
  INTEGER :: kcomp

!  REAL(r8) :: p80, sw_temp, rho
!  EXTERNAL  p80, sw_temp, rho

  icount = 0
  DO i = 1, N


     icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Model temperature tracer is usually "potential temperature"
!    * Model vertical grid is usually in meters
!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry,
!    convert these 2 model vars (input to this routine)
!    - depth [m] => convert to pressure [db]
!    - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     IF (optP .eq. 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p(i) = p80(depth(i), lat(i))
     ELSEIF (optP .eq. 'db') THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p(i) = depth(i)
     ELSE
        PRINT *,"optP must be either 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
       IF (optT .eq. "Tpot"   &
      &     .AND. tempot(i) .gt. -3.0_r8.AND. tempot(i) .lt. 100.0_r8) THEN
!       This is the case for most models and some data
!       a) Convert the pot. temp on todays "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot(i) - 0.0002_r8) / 0.99975_r8
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p(i), 0.0_r8 )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis(i) = 0.99975_r8*tempis68 + 0.0002_r8
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (optT .eq. 'Tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis(i) = tempot(i)

     ELSE
        PRINT *,"optT must be either 'Tpot' or 'Tinsitu'"
        STOP
     ENDIF

!    ================================================================
!    Carbonate chemistry computations
!    ================================================================
     IF (dic(i) > 0. .AND. dic(i) < 1.0e+4) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (       sal(i) < 0.   &
      &            .OR.  alk(i) < 0.   &
      &       .OR.  dic(i) < 0.   &
      &       .OR.  sil(i) < 0.   &
      &       .OR. phos(i) < 0.   &
      &       .OR.  sal(i) > 1e+3 &
      &       .OR.  alk(i) > 1e+4 &
      &       .OR.  dic(i) > 1e+4 &
      &       .OR.  sil(i) > 1e+3 &
      &       .OR. phos(i) > 1e+3) THEN
!           PRINT *, 'i, icount,' , &
!      &          'tempot(i), sal(i), alk(i), dic(i), sil(i), phos(i) =', &
!      &          i, icount, 
!      &          tempot(i), sal(i), alk(i), dic(i), sil(i), phos(i)
        ENDIF
!       Zero out any negative salinity, phosphate, silica, dic, and alk
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF
        IF (phos(i) < 0.0) THEN
           sphos = 0.0
        ELSE
           sphos = phos(i)
        ENDIF
        IF (sil(i) < 0.0) THEN
           ssil = 0.0
        ELSE
           ssil = sil(i)
        ENDIF
        IF (dic(i) < 0.0) THEN
          sdic = 0.0
        ELSE
          sdic = dic(i)
        ENDIF
        IF (alk(i) < 0.0) THEN
          salk = 0.0
        ELSE
          salk = alk(i)
        ENDIF

!       Absolute temperature (Kelvin) & related variables
        t = DBLE(tempis(i))
        tk = 273.15d0 + t
        invtk=1.0d0/tk
        dlogtk=LOG(tk)

!       Rressure effect (prb is in bars)
        prb = DBLE(p(i)) / 10.0d0

!       Salinity (equivalent array in double precision)
        s = DBLE(ssal)

!       Get all equilibrium constants and total concentrations of SO4, F, B
        CALL constants(aKh, aK1, aK2, aKb, aKw, aKs, aKf, aKspc, aKspa,  &
      &                 aK1p, aK2p, aK3p, aKsi,                           &
      &                 aSt, aFt, aBt,                                    &
      &                 tempot(i), sal(i),                                &
      &                 depth(i), lat(i), 1,                              &
      &                 optT, optP, optB, optK1K2, optKf)

!       Unlike f77, in F90 we cant assign an array (dimen=1) to a scalar in a routine argument
!       Thus, set scalar constants equal to array (dimension=1) values required as arguments
        Kh = aKh(1) ; K1 = aK1(1) ; K2 = aK2(1) ; Kb = aKb(1) ; Kw = aKw(1) 
        Ks = aKs(1) ; Kf = aKs(1) ; Kspc = aKspc(1) ; Kspa = aKspa(1) 
        K1p = aK1p(1) ; K2p = aK2p(1) ; K3p = aK3p(1) ; Ksi = aKsi(1)
        St = aSt(1) ; Ft = aFt(1) ; Bt = aBt(1)


!       Compute in-situ density [kg/m^3]
        rhoSW(i) =  rho_calc(ssal, tempis68, REAL(prb))

!       Either convert units of DIC and ALK (MODEL case) or not (DATA case)
        IF     (optRHO == 'mol/kg') THEN
!          No conversion:
!          print *,'DIC and ALK already given in mol/kg (std DATA units)'
           drho = 1.
        ELSEIF (optRHO .eq. 'mol/m3') THEN
!          Do conversion:
!          print *,"DIC and ALK given in mol/m^3 (std MODEL units)"
           drho = DBLE(rhoSW(i))
        ELSE
           PRINT *,"optRHO must be either 'mol/kg' or 'mol/m3'"
           STOP
        ENDIF

        tc = DBLE(sdic/1000)/drho
        ta = DBLE(salk/1000)/drho
        sit = DBLE(ssil/1000)/drho
        pt  = DBLE(sphos/1000)/drho


!       Computation loop for [H+]
        ah1 = 1.e-8_r8
        kcomp = 0
2       kcomp = kcomp+1
!       print *, "kcomp  =  ", kcomp
        HSO4 = St/(1.0d0 + Ks/(ah1/(1.0d0 + St/Ks)))
        HF = 1.0d0/(1.0d0 + Kf/ah1)
        HSI = 1.0d0/(1.0d0 + ah1/Ksi)
        HPO4 = (K1p*K2p*(ah1 + 2.*K3p) - ah1**3)/   &
     &               (ah1**3 + K1p*ah1**2 + K1p*K2p*ah1 +   &
     &          K1p*K2p*K3p)

        ab = Bt/(1.0d0 + ah1/Kb)
        aw = Kw/ah1 - ah1/(1.0d0 + St/Ks)
        ac = ta + hso4 - sit*hsi - ab - aw + Ft*hf - pt*hpo4
        ah2 = SQRT((tc - ac)**2 + 4.0d0*(ac*K2/K1)*(2.0d0*tc - ac))
        ah2 = 0.5d0*K1/ac*((tc - ac) + ah2)
        erel = (ah2 - ah1)/ah2

        IF (ABS(erel) >= 1.e-7_r8) THEN
           ah1 = ah2
           IF (kcomp < 25) GOTO 2
        END IF

        IF (ah1 > 0.d0) THEN
           ph(i) = REAL(-LOG10(ah1))
        ELSE
           ph(i) = 1.e20_r4
        ENDIF

!       Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
        cu = (2.0d0 * tc - ac) / (2.0d0 + K1 / ah1)
        cb = K1 * cu / ah1
        cc = K2 * cb / ah1

!       If optRHO = 'mol/m3', then:
!       convert output var concentrations from mol/kg to mol/m^3
!       e.g., for case when drho = 1028, multiply by [1.028 kg/L  x  1000 L/m^3])
        co2(i)  = REAL(cu * drho)
        hco3(i) = REAL(cb * drho)
        co3(i)  = REAL(cc * drho)

!       Determine CO2 pressure and fugacity (in microatm)
!       NOTE: equation below for pCO2 requires CO2* in mol/kg
        dfCO2 = cu   * 1.e6_r8/Kh

        B=(-1636.75d0+12.0408d0*tk-0.0327957d0*(tk*tk)+0.0000316528d0* &
      &       (tk*tk*tk))*1e-6

        dpCO2= dfCO2 / EXP(100000.0d0 &
       &      *(B+ 2.0d0*(57.7d0 - 0.118d0*tk)*1e-6_r8) &
       &      / (8.314d0*tk) &
       &      )

        fCO2(i) = REAL(dfCO2)
        ppCO2(i) = REAL(dpCO2)

!       Determine Omega Calcite et Aragonite
     OmegaA(i) = REAL(((0.01028d0*s/35.0d0)*cc)/Kspa)
        OmegaC(i) = REAL(((0.01028d0*s/35.0d0)*cc)/Kspc)



!       Determine Revelle Factor (BetaD w/ code from Seacarb software)
        DD = -((-Kb*Bt)/((ah1 + Kb)*(ah1+Kb))) - (-Kw/(ah1*ah1)) + 1.0d0
        A = (2.0d0*K2*(2.0d0*cc + cb) + ah1*(ah1 + 2.0d0*K2)*DD)             &
     &       / ((ah1 + 2.0d0*K2)*(ah1 + 2.0d0*K2))
        B = ( ( (2.0d0*cc + cb) * ah1)/((ah1 + 2.0d0*K2)*K1) + (ah1/K1)* A )
        C = (-K2*(2.0d0*cc + cb) + K2*(2.0d0*K2 + ah1)*DD)                   &
     &       / ((ah1+2.0d0*K2)*(ah1+2.0d0*K2))
        PhiD = -1.0d0/(ah1*LOG(10.0d0) * ( B+A+C ) )
        BetaD(i) = REAL(-ah1*LOG(10.0d0)*tc/cu*B*PhiD)

     ELSE

        ph(i)     = 1.e20_r8
        ppco2(i)   = 1.e20_r8
        fco2(i)   = 1.e20_r8
        co2(i)    = 1.e20_r8
        hco3(i)   = 1.e20_r8
        co3(i)    = 1.e20_r8
        OmegaA(i) = 1.e20_r8
        OmegaC(i) = 1.e20_r8
        BetaD(i)  = 1.e20_r8
        rhoSW(i)  = 1.e20_r8
        p(i)      = 1.e20_r8
        tempis(i) = 1.e20_r8
 
      ENDIF

      END DO

      RETURN
      END SUBROUTINE vars

!> Compute thermodynamic constants
!! FROM temperature, salinity, and pressure (1D arrays)
!! @item aggr information about the aggregates
!! @item Handle special case
SUBROUTINE constants(Kh, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,  &
     &                K1p, K2p, K3p, Ksi,                      &
     &                St, Ft, Bt,                              &
     &                tempot, sal,                             &
     &                depth, lat, N,                           &
     &                optT, optP, optB, optK1K2, optKf)

  !   Purpose:
  !     Compute thermodynamic constants
  !     FROM: temperature, salinity, and pressure (1D arrays)

  !     INPUT variables:
  !     ================
  !     depth   = depth [m]     (with optP=0, i.e., for a z-coordinate model vertical grid is depth, not pressure)
  !             = pressure [db] (with optP=1)
  !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP=0)
  !             = dummy array (unused when optP=1)
  !     tempot  = potential temperature [degrees C] (with optT=0, i.e., models carry tempot, not temp)
  !             = in situ   temperature [degrees C] (with optT=1, e.g., for data)
  !     sal     = salinity in [psu]

  !     optT: choose in situ vs. potential temperature as input
  !     ---------
  !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
  !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
  !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)

  !     optP: choose depth (m) vs pressure (db) as input
  !     ---------
  !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
  !       -> 'db' means "depth" input is already in situ pressure [db], not m (p = depth)

  !     optB:
  !     ---------
  !       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
  !       -> 'l10' means use   new formulation of Lee et al. (2010) for total Boron

  !     optK1K2:
  !     ---------
  !       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
  !                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
  !       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
  !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu

  !     optKf:
  !     ----------
  !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
  !               **** BUT Valid only for  9 < T < 33°C and 10 < S < 40.
  !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)

  !     OUTPUT variables:
  !     =================
  !     Kh, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
  !     St, Ft, Bt

  USE mod_kinds
  IMPLICIT NONE

! Input variables
  !>     number of records
  INTEGER, INTENT(in) :: N
  !> in <b>situ temperature</b> (when optT='Tinsitu', typical data) 
  !! OR <b>potential temperature</b> (when optT='Tpot', typical models) [degree C]
  REAL(kind=r8), INTENT(in),    DIMENSION(N) :: tempot
  !> depth in <b>meters</b> (when optP='m') or <b>decibars</b> (when optP='db')
  REAL(kind=r8), INTENT(in),    DIMENSION(N) :: depth
  !> latitude <b>[degrees north]</b>
  REAL(kind=r8), INTENT(in),    DIMENSION(N) :: lat
  !> salinity <b>[psu]</b>
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: sal
!f2py optional , depend(sal) :: n=len(sal)

  !> for tempot input, choose \b 'Tinsitu' for in situ Temp or 
  !! \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
  CHARACTER(*), INTENT(in) :: optT
  !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
  CHARACTER(*), INTENT(in) :: optP
  !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
  !! The 'l10' formulation is based on 139 measurements (instead of 20),
  !! uses a more accurate method, and
  !! generally increases total boron in seawater by 4% 
  CHARACTER(*), INTENT(in) :: optB
  !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
  CHARACTER(*), INTENT(in) :: optKf
  !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) 
  CHARACTER(*), INTENT(in) :: optK1K2

! Ouput variables
  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kh
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K2
  !> equilibrium constant for dissociation of boric acid 
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ksi
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Ft
  !> total boron
  !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: Bt

! Local variables
  REAL(kind=r8) :: ssal
  REAL(kind=r8) :: p
  REAL(kind=r8) :: tempis68, tempot68
  REAL(kind=r8) :: tempis
  REAL(kind=r8) :: is, invtk, dlogtk, is2, s2, sqrtis
  REAL(kind=r8) :: Ks_0p, Kf_0p
  REAL(kind=r8) :: total2free, free2SWS, total2SWS, SWS2total
  REAL(kind=r8) :: total2free_0p, free2SWS_0p, total2SWS_0p
! REAL(kind=r8) :: free2SWS, free2SWS_0p

  REAL(kind=r8) :: R

  REAL(kind=r8) :: pK1o, ma1, mb1, mc1, pK1
  REAL(kind=r8) :: pK2o, ma2, mb2, mc2, pK2

  REAL(kind=r8), DIMENSION(12) :: a0, a1, a2, b0, b1, b2
  REAL(kind=r8), DIMENSION(12) :: deltav, deltak, lnkpok0
  REAL(kind=r8) :: tmp, nKhwe74

  INTEGER :: i, icount, ipc

  REAL(kind=r8) :: t, tk, prb
  REAL(kind=r8) :: s, sqrts, s15, scl

!  REAL ::   p80, sw_temp
!  EXTERNAL p80, sw_temp

  ! CONSTANTS
  ! =========
  ! Constants in formulation for Pressure effect on Ks (Millero, 95)
  ! with corrected coefficients for Kb, Kw, Ksi, etc.

  ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
  !            9) K1P, 10) K2P, 11) K3P, 12) Ksi

  DATA a0 /-25.5_r8, -15.82_r8, -29.48_r8, -20.02_r8, &
          -18.03_r8,  -9.78_r8, -48.76_r8, -46._r8, &
          -14.51_r8, -23.12_r8, -26.57_r8, -29.48_r8/
  DATA a1 /0.1271_r8, -0.0219_r8, 0.1622_r8, 0.1119_r8, &
           0.0466_r8, -0.0090_r8, 0.5304_r8, 0.5304_r8, &
           0.1211_r8, 0.1758_r8, 0.2020_r8, 0.1622_r8/
  DATA a2 /     0.0_r8,       0.0_r8, -2.608e-3_r8, -1.409e-3_r8, &
           0.316e-3_r8, -0.942e-3_r8,  0.0_r8,       0.0_r8, &
          -0.321e-3_r8, -2.647e-3_r8, -3.042e-3_r8, -2.6080e-3_r8/
  DATA b0 /-3.08e-3_r8, 1.13e-3_r8,  -2.84e-3_r8,   -5.13e-3_r8, &
           -4.53e-3_r8, -3.91e-3_r8, -11.76e-3_r8, -11.76e-3_r8, &
           -2.67e-3_r8, -5.15e-3_r8,  -4.08e-3_r8,  -2.84e-3_r8/
  DATA b1 /0.0877e-3_r8, -0.1475e-3_r8, 0.0_r8,       0.0794e-3_r8, &
           0.09e-3_r8,    0.054e-3_r8,  0.3692e-3_r8, 0.3692e-3_r8, &
           0.0427e-3_r8,  0.09e-3_r8,   0.0714e-3_r8, 0.0_r8/
  DATA b2 /12*0.0_r8/

  R = 83.14472_r8

  icount = 0
  DO i = 1, N
     icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Modeled temperature tracer is "potential temperature"
!    * Modeled vertical grid is in meters

!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry, if appropriate,
!    convert these 2 model vars (input to this routine)
!     - depth [m] => convert to pressure [db]
!     - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
     IF (optP .eq. 'm' ) THEN
!       Compute pressure [db] from depth [m] and latitude [degrees]
        p = p80(depth(i), lat(i))
     ELSEIF (optP .eq. 'db' ) THEN
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
        p = depth(i)
     ELSE
        PRINT *,"optP must be 'm' or 'db'"
        STOP
     ENDIF

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
     IF (optT .eq. "Tpot" &
     &     .AND. tempot(i) .gt. -3.0_r8 .AND.tempot(i) .lt. 100.0_r8) THEN
!       This is the case for most models and some data
!       a) Convert the pot. temp on todays "ITS 90" scale to older IPTS 68 scale
!       (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
        tempot68 = (tempot(i) - 0.0002_r8) / 0.99975_r8
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
        tempis68 = sw_temp(sal(i), tempot68, p, 0._r8)
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
        tempis = 0.99975_r8*tempis68 + 0.0002_r8
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
     ELSEIF (optT .eq. 'Tinsitu') THEN
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
        tempis = tempot(i)

     ELSE
        PRINT *,"optT must be 'Tpot' or 'Tinsitu'"
        STOP
     ENDIF

!    Compute constants
     IF (tempot(i) >= -2. .AND. tempot(i) < 1.0e+2) THEN
!       Test to indicate if any of input variables are unreasonable
        IF (      sal(i) < 0.  .OR.  sal(i) > 1e+3) THEN
           PRINT *, 'i, icount, tempot(i), sal(i) =', &
      &          i, icount, tempot(i), sal(i)
        ENDIF
!       Zero out negative salinity (prev case for OCMIP2 model w/ slightly negative S in some coastal cells)
        IF (sal(i) < 0.0) THEN
           ssal = 0.0
        ELSE
           ssal = sal(i)
        ENDIF

!       Compute absolute temperature (Kelvin) and related values
        t = DBLE(tempis)
        tk = 273.15d0 + t
        invtk=1.0d0/tk
        dlogtk=LOG(tk)

!       Set pressure effect (prb is in bars)
        prb = DBLE(p) / 10.0d0

!       Salinity and simply related values
        s = DBLE(ssal)
        s2=s*s
        sqrts=SQRT(s)
        s15=s**1.5d0
        scl=s/1.80655d0

!       Ionic strength:
        is = 19.924d0*s/(1000.0d0 - 1.005d0*s)
        is2 = is*is
        sqrtis = SQRT(is)

!       Total concentrations for sulfate, fluoride, and boron

!       Sulfate: Morris & Riley (1966)
        St(i) = 0.14d0 * scl/96.062d0

!       Fluoride:  Riley (1965)
        Ft(i) = 0.000067d0 * scl/18.9984d0

!       Boron:
        IF (optB .eq. 'l10') THEN
!          New formulation from Lee et al (2010)
           Bt(i) = 0.0002414d0 * scl/10.811d0
        ELSEIF (optB .eq. 'u74') THEN
!          Classic formulation from Uppström (1974)
           Bt(i) = 0.000232d0  * scl/10.811d0
        ELSE
           PRINT *,"optB must be 'l10' or 'u74'"
           STOP
        ENDIF

!       Kh (K Henry)
!       CO2(g) <-> CO2(aq.)
!       Kh  = [CO2]/ fCO2
!       Weiss (1974)   [mol/kg/atm]
        tmp = 9345.17d0*invtk - 60.2409d0 + 23.3585d0 * LOG(tk/100.0d0)
        nKhwe74 = tmp + s*(0.023517d0 - 0.00023656d0*tk &
     &                      + 0.0047036e-4_r8                       &
     &           *tk*tk)
        Kh(i) = EXP(nKhwe74)

!       K1 = [H][HCO3]/[H2CO3]
!       K2 = [H][CO3]/[HCO3]
        IF (optK1K2 .eq. 'l') THEN
!         Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total pH scale)
          K1(i) = 10.0d0**(-1.0d0*(3633.86d0*invtk - 61.2172d0       &
     &            + 9.6777d0*dlogtk                                  &
     &             - 0.011555d0*s + 0.0001152d0*s2))

          K2(i) = 10.0d0**(-1*(471.78d0*invtk + 25.9290d0 - 3.16967d0&
     &                                                  *dlogtk      &
     &             - 0.01781d0*s + 0.0001122d0*s2))
        ELSEIF (optK1K2 .eq. 'm10') THEN
!         Millero (2010, Mar. Fresh Wat. Res.) (total pH scale)
          pK1o = 6320.813d0*invtk + 19.568224d0*dlogtk -126.34048d0
          ma1 = 13.4051d0*sqrts + 0.03185d0*s - (5.218e-5)*s2
          mb1 = -531.095d0*sqrts - 5.7789d0*s
          mc1 = -2.0663d0*sqrts
          pK1 = pK1o + ma1 + mb1*invtk + mc1*dlogtk
          K1(i) = 10.0d0**(-pK1) 

          pK2o = 5143.692d0*invtk + 14.613358d0*dlogtk -90.18333d0
          ma2 = 21.5724d0*sqrts + 0.1212d0*s - (3.714e-4)*s2
          mb2 = -798.292d0*sqrts - 18.951d0*s
          mc2 = -3.403d0*sqrts
          pK2 = pK2o + ma2 + mb2*invtk + mc2*dlogtk
          K2(i) = 10.0d0**(-pK2)
        ELSE
           PRINT *, "optK1K2 must be either 'l' or 'm10'"
           STOP
        ENDIF

!       Kb = [H][BO2]/[HBO2]
!       (total scale)
!       Millero p.669 (1995) using data from Dickson (1990)
        Kb(i) = EXP((-8966.90d0 - 2890.53d0*sqrts - 77.942d0*s +  &
    &            1.728d0*s15 - 0.0996d0*s2)*invtk +              &
    &            (148.0248d0 + 137.1942d0*sqrts + 1.62142d0*s) +   &
    &            (-24.4344d0 - 25.085d0*sqrts - 0.2474d0*s) *      &
    &            dlogtk + 0.053105d0*sqrts*tk)

!       K1p = [H][H2PO4]/[H3PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 65
!       Use Millero equations 115.540 constant instead of 115.525 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K1p(i) = EXP(-4576.752d0*invtk + 115.540d0 - 18.453d0*dlogtk +  &
    &             (-106.736d0*invtk + 0.69171d0) * sqrts +             &
    &             (-0.65643d0*invtk - 0.01844d0) * s)

!       K2p = [H][HPO4]/[H2PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!       Millero (1995), p.670, eq. 66
!       Use Millero equations 172.1033 constant instead of 172.0833 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K2p(i) = EXP(-8814.715d0*invtk + 172.1033d0 - 27.927d0*dlogtk +  &
    &             (-160.340d0*invtk + 1.3566d0)*sqrts +                 &
    &             (0.37335d0*invtk - 0.05778d0)*s)

!       K3p = [H][PO4]/[HPO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 67
!       Use Millero equations 18.126 constant instead of 18.141 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        K3p(i) = EXP(-3070.75d0*invtk - 18.126d0 +            &
     &            (17.27039d0*invtk + 2.81197d0) *             &
     &            sqrts + (-44.99486d0*invtk - 0.09984d0) * s)

!       Ksi = [H][SiO(OH)3]/[Si(OH)4]
!       (seawater scale)
!       Millero (1995), p.671, eq. 72
!       Use Millero equations 117.400 constant instead of 117.385 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Ksi(i) = EXP(-8904.2d0*invtk  + 117.400d0 - 19.334d0*dlogtk +  &
     &            (-458.79d0*invtk + 3.5913d0) * sqrtis +             &
     &            (188.74d0*invtk - 1.5998d0) * is +                  &
     &            (-12.1652d0*invtk + 0.07871d0) * is2 +              &
     &            LOG(1.0 - 0.001005d0*s))

!       Kw = [H][OH]
!       (seawater scale)
!       Millero (1995) p.670, eq. 63 from composite data
!       Use Millero equations 148.9802 constant instead of 148.9652 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
        Kw(i) = EXP(-13847.26d0*invtk + 148.9802d0 - 23.6521d0*dlogtk +  &
     &          (118.67d0*invtk - 5.977d0 + 1.0495d0 * dlogtk) *          &
     &          sqrts - 0.01615d0 * s)

!       Ks = [H][SO4]/[HSO4]
!       (free scale)
!       Dickson (1990, J. chem. Thermodynamics 22, 113)
        Ks_0p = EXP(-4276.1d0*invtk + 141.328d0 - 23.093d0*dlogtk      &
     &           + (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk)      &
     &           + (35474.d0*invtk - 771.54 + 114.723d0*dlogtk) * is   &
     &           * sqrtis- 2698.d0*invtk*is**1.5 + 1776.d0*invtk*is2   &
     &           + LOG(1.0d0 - 0.001005d0*s))

!       Kf = [H][F]/[HF]
!       (total scale)
        IF (optKf .eq. 'dg') THEN
!          Dickson and Riley (1979) -- change pH scale to total (following Dickson & Goyet, 1994)
           Kf_0p = EXP(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis +  &
      &             LOG(1.0d0 - 0.001005d0*s) +                     &
      &             LOG(1.0d0 + St(i)/Ks_0p))
        ELSEIF (optKf .eq. 'pf') THEN
!          Perez and Fraga (1987) - Already on Total scale (no need for last line above)
!          Formulation as given in Dickson et al. (2007)
           Kf_0p = EXP(874.d0*invtk - 9.68d0 + 0.111d0*sqrts)
        ELSE
           PRINT *, "optKf must be either 'dg' or 'pf'"
           STOP
        ENDIF

!       Kspc (calcite) - apparent solubility product of calcite
!       (no scale)
!       Kspc = [Ca2+] [CO32-] when soln is in equilibrium w/ calcite
!       Mucci 1983 mol/kg-soln
        Kspc(i) = 10d0**(-171.9065d0 - 0.077993d0*tk + 2839.319d0/tk    &
      &           + 71.595d0*LOG10(tk)                             &
      &           + (-0.77712d0 + 0.0028426d0*tk + 178.34d0/tk)*sqrts  &
      &           -0.07711d0*s + 0.0041249d0*s15 )

!       Kspa (aragonite) - apparent solubility product of aragonite
!       (no scale)
!       Kspa = [Ca2+] [CO32-] when soln is in equilibrium w/ aragonite
!       Mucci 1983 mol/kg-soln
        Kspa(i) = 10.d0**(-171.945d0 - 0.077993d0*tk + 2903.293d0/tk &
      &       +71.595d0*LOG10(tk) &
      &       +(-0.068393d0 + 0.0017276d0*tk + 88.135d0/tk)*sqrts &
      &       -0.10018d0*s + 0.0059415d0*s15 )

!       Pressure effect on Ks (based on Millero, (1995)
!           index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
!                  K1p(9), K2p(10), K3p(11), Ksi(12)
        DO ipc = 1, 12
           deltav(ipc)  =  a0(ipc) + a1(ipc) *t + a2(ipc) *t*t
           deltak(ipc)   = (b0(ipc)  + b1(ipc) *t + b2(ipc) *t*t)
           lnkpok0(ipc)  = (-(deltav(ipc)) &
      &          +(0.5d0*deltak(ipc) * prb) &
      &          )                         * prb/(R*tk)
        END DO

!       Pressure correction on Ks (Free scale)
        Ks(i) = Ks_0p*EXP(lnkpok0(5))
!       Conversion factor total -> free scale
        total2free     = 1.d0/(1.d0 + St(i)/Ks(i))      ! Kfree = Ktotal*total2free
!       Conversion factor total -> free scale at pressure zero
        total2free_0p  = 1.d0/(1.d0 + St(i)/Ks_0p)   ! Kfree = Ktotal*total2free

!       Pressure correction on Kf
!       Kf must be on FREE scale before correction
        Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
        Kf(i) = Kf_0p * EXP(lnkpok0(6)) !Pressure correction (on Free scale)
        Kf(i) = Kf(i)/total2free        !Convert back from Free to Total scale

!       Convert between seawater and total scales (taken from seacarb)
        free2SWS  = 1 + St(i)/Ks(i) + Ft(i)/(Kf(i)*total2free)  ! using Kf on free scale
        total2SWS = total2free * free2SWS                       ! KSWS = Ktotal*total2SWS
        SWS2total = 1 / total2SWS

!       Conversion at pressure zero
        free2SWS_0p  = 1 + St(i)/Ks_0p + Ft(i)/(Kf_0p)  !using Kf on free scale
        total2SWS_0p = total2free_0p * free2SWS_0p ! KSWS = Ktotal*total2SWS

!       Convert from Total to Seawater scale before pressure correction
!       Must change to SEAWATER scale: K1, K2, Kb
        K1(i)  = K1(i)*total2SWS_0p
        K2(i)  = K2(i)*total2SWS_0p
        Kb(i)  = Kb(i)*total2SWS_0p

!       Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw
!       K1p(i) = K1p(i)*total2SWS_0p 
!       K2p(i) = K2p(i)*total2SWS_0p
!       K3p(i) = K3p(i)*total2SWS_0p
!       Ksi(i) = Ksi(i)*total2SWS_0p
!       Kw(i)  = Kw(i)*total2SWS_0p

!       Keep other contants on same scale:
!          - Ks (already on Free scale)
!          - Kspc, Kspa (no scale)

!       Perform actual pressure correction (on seawater scale)
        K1(i)   = K1(i)*EXP(lnkpok0(1))
        K2(i)   = K2(i)*EXP(lnkpok0(2))
        Kb(i)   = Kb(i)*EXP(lnkpok0(3))
        Kw(i)   = Kw(i)*EXP(lnkpok0(4))
        Kspc(i) = Kspc(i)*EXP(lnkpok0(7))
        Kspa(i) = Kspa(i)*EXP(lnkpok0(8))
        K1p(i)  = K1p(i)*EXP(lnkpok0(9))
        K2p(i)  = K2p(i)*EXP(lnkpok0(10))
        K3p(i)  = K3p(i)*EXP(lnkpok0(11))
        Ksi(i)  = Ksi(i)*EXP(lnkpok0(12))

!       Convert back to original total scale:
        K1(i)  = K1(i)*SWS2total
        K2(i)  = K2(i)*SWS2total
        K1p(i) =  K1p(i)*SWS2total
        K2p(i) = K2p(i)*SWS2total
        K3p(i) = K3p(i)*SWS2total
        Kb(i)  = Kb(i)*SWS2total
        Ksi(i) = Ksi(i)*SWS2total
        Kw(i)  = Kw(i)*SWS2total

     ELSE

        Kh(i)   = 1.e20_r8
        K1(i)   = 1.e20_r8
        K2(i)   = 1.e20_r8
        Kb(i)   = 1.e20_r8
        Kw(i)   = 1.e20_r8
        Ks(i)   = 1.e20_r8
        Kf(i)   = 1.e20_r8
        Kspc(i) = 1.e20_r8
        Kspa(i) = 1.e20_r8
        K1p(i)  = 1.e20_r8
        K2p(i)  = 1.e20_r8
        K3p(i)  = 1.e20_r8
        Ksi(i)  = 1.e20_r8
        Bt(i)   = 1.e20_r8
        Ft(i)   = 1.e20_r8
        St(i)   = 1.e20_r8

     ENDIF

  END DO

  RETURN
END SUBROUTINE constants


!>     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db).
!!     This subroutine is needed because rho is a function (cannot accept arrays)
SUBROUTINE rhoinsitu(salt, tempis, pdbar, N, rhois)

  !     Purpose:
  !     Compute in situ density from salinity (psu), in situ temperature (C), & pressure (db)
  !     Needed because rho is a function (cannot accept arrays)

  USE mod_kinds
  IMPLICIT NONE

  INTEGER :: N

! INPUT variables
  ! salt   = salinity [psu]
  ! tempis = in situ temperature [C]
  ! pdbar  = pressure [db]

  !> salinity [psu]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: salt
  !> in situ temperature [C]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: tempis
  !> pressure [db]
  REAL(kind=r4), INTENT(in), DIMENSION(N) :: pdbar
!f2py optional , depend(salt) :: n=len(salt)

! OUTPUT variables:
  ! rhois  = in situ density

  !> in situ density [kg/m3]
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: rhois

! Local variables
  INTEGER :: i

  !REAL(kind=r8) ::  rho
  !EXTERNAL rho

  DO i = 1,N
     rhois(i) = rho_calc(salt(i), tempis(i), pdbar(i)/10._r4)
  END DO

  RETURN
END SUBROUTINE rhoinsitu
!>     Compute pressure [db] from depth [m] & latitude [degrees north].
!!     This subroutine is needed because p80 is a function (cannot accept arrays)
SUBROUTINE depth2press(depth, lat, pdbar, N)

  !     Purpose:
  !     Compute pressure [db] from depth [m] & latitude [degrees north].
  !     Needed because p80 is a function (cannot accept arrays)

  USE mod_kinds
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> depth [m]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: depth
  !> latitude [degrees]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: lat
!f2py optional , depend(depth) :: n=len(depth)

! OUTPUT variables:
  !> pressure [db]
  REAL(kind=r4), INTENT(out), DIMENSION(N) :: pdbar

  !     Local variables
  INTEGER :: i

!  REAL(kind=r8) ::  p80
!  EXTERNAL p80

  DO i = 1,N
     pdbar(i) = p80(depth(i), lat(i))
  END DO

  RETURN
END SUBROUTINE depth2press
!>    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
!!    This subroutine is needed because sw_ptmp is a function (cannot accept arrays)
SUBROUTINE tpot(salt, tempis, press, pressref, N, tempot)
  !    Purpose:
  !    Compute potential temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_ptmp is a function (cannot accept arrays)

  USE mod_kinds
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> salinity [psu]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: salt
  !> in situ temperature [C]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: tempis
  !> pressure [db]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: press
!f2py optional , depend(salt) :: n=len(salt)
  !> pressure reference level [db]
  REAL(kind=r8), INTENT(in) :: pressref

! OUTPUT variables:
  !> potential temperature [C] for pressref
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: tempot

  REAL(kind=r8) :: dsalt, dtempis, dpress, dpressref
  REAL(kind=r8) :: dtempot

  INTEGER :: i

!  REAL(kind=r8) :: sw_ptmp
!  EXTERNAL sw_ptmp

  DO i = 1,N
     dsalt     = DBLE(salt(i))
     dtempis   = DBLE(tempis(i))
     dpress    = DBLE(press(i))
     dpressref = DBLE(pressref)

     dtempot   = sw_ptmp(dsalt, dtempis, dpress, dpressref)

     tempot(i) = REAL(dtempot)
  END DO

  RETURN
END SUBROUTINE tpot

!>    Compute in situ temperature from arrays of potential temp, salinity, and pressure.
!!    This subroutine is needed because sw_temp is a function (cannot accept arrays)
SUBROUTINE tis(salt, tempot, press, pressref, N, tempis)
  !    Purpose:
  !    Compute in situ temperature from arrays of in situ temp, salinity, and pressure.
  !    Needed because sw_temp is a function (cannot accept arrays)

  USE mod_kinds
  IMPLICIT NONE

  !> number of records
  INTEGER, intent(in) :: N

! INPUT variables
  !> salinity [psu]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: salt
  !> potential temperature [C]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: tempot
  !> pressure [db]
  REAL(kind=r8), INTENT(in), DIMENSION(N) :: press
!f2py optional , depend(salt) :: n=len(salt)
  !> pressure reference level [db]
  REAL(kind=r8), INTENT(in) :: pressref

! OUTPUT variables:
  !> in situ temperature [C] 
  REAL(kind=r8), INTENT(out), DIMENSION(N) :: tempis

! REAL(kind=r8) :: dsalt, dtempis, dpress, dpressref
! REAL(kind=r8) :: dtempot

  INTEGER :: i

! REAL(kind=r8) :: sw_temp
!  REAL(kind=r8) :: sw_temp
!  EXTERNAL sw_temp

  DO i = 1,N
    !dsalt     = DBLE(salt(i))
    !dtempot   = DBLE(tempot(i))
    !dpress    = DBLE(press(i))
    !dpressref = DBLE(pressref)
    !dtempis   = sw_temp(dsalt, dtempot, dpress, dpressref)
    !tempis(i) = REAL(dtempis)

     tempis   = sw_temp(salt(i), tempot(i), press(i), pressref)
  END DO

  RETURN
END SUBROUTINE tis

!> Function to compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)
FUNCTION rho_calc (salt, temp, pbar)

  ! Compute in situ density from salinity (psu), in situ temperature (C), & pressure (bar)

  USE mod_kinds
  IMPLICIT NONE

  !> salinity [psu]
  REAL(kind=r8) :: salt
  !> in situ temperature (C)
  REAL(kind=r8) :: temp
  !> pressure (bar) [Note units: this is NOT in decibars]
  REAL(kind=r4) :: pbar

  REAL(kind=r8) :: s, t, p
! REAL(kind=r8) :: t68
  REAL(kind=r8) :: X
  REAL(kind=r8) :: rhow, rho0
  REAL(kind=r8) :: a, b, c
  REAL(kind=r8) :: Ksbmw, Ksbm0, Ksbm
  REAL(kind=r8) :: drho

  REAL(kind=r8) :: rho_calc

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = in situ temperature [degree C (IPTS-68)]
  !     p  = pressure            [bar] !!!!  (not in [db]

  s = DBLE(salt)
  t = DBLE(temp)
  p = DBLE(pbar)

! Convert the temperature on todays "ITS 90" scale to older IPTS 68 scale
! (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
! According to Best-Practices guide, line above should be commented & 2 lines below should be uncommented
! Guides answer of rho (s=35, t=25, p=0) = 1023.343 is for temperature given on ITPS-68 scale
! t68 = (T - 0.0002) / 0.99975
! X = t68
! Finally, dont do the ITS-90 to IPTS-68 conversion (T input var now already on IPTS-68 scale)
  X = T

! Density of pure water
  rhow = 999.842594d0 + 6.793952e-2_r8*X          &
     &       -9.095290e-3_r8*X*X + 1.001685e-4_r8*X**3  &
     &  -1.120083e-6_r8*X**4 + 6.536332e-9_r8*X**5

! Density of seawater at 1 atm, P=0
  A = 8.24493e-1_r8 - 4.0899e-3_r8*X                         &
     &  + 7.6438e-5_r8*X*X - 8.2467e-7_r8*X**3 + 5.3875e-9_r8*X**4
  B = -5.72466e-3_r8 + 1.0227e-4_r8*X - 1.6546e-6_r8*X*X
  C = 4.8314e-4_r8

  rho0 = rhow + A*S + B*S*SQRT(S) + C*S**2.0d0

! Secant bulk modulus of pure water
! The secant bulk modulus is the average change in pressure
! divided by the total change in volume per unit of initial volume.
  Ksbmw = 19652.21d0 + 148.4206d0*X - 2.327105d0*X*X &
     &  + 1.360477e-2_r8*X**3 - 5.155288e-5_r8*X**4

! Secant bulk modulus of seawater at 1 atm
  Ksbm0 = Ksbmw + S*( 54.6746d0 - 0.603459d0*X + 1.09987e-2_r8*X**2 &
     &  - 6.1670e-5_r8*X**3) &
     &  + S*SQRT(S)*( 7.944e-2_r8 + 1.6483e-2_r8*X - 5.3009e-4_r8*X**2)

! Secant bulk modulus of seawater at S,T,P
  Ksbm = Ksbm0 &
     & + P*(3.239908d0 + 1.43713e-3_r8*X + 1.16092e-4_r8*X**2 &
     &  - 5.77905e-7_r8*X**3) &
     & + P*S*(2.2838e-3_r8 - 1.0981e-5_r8*X - 1.6078e-6_r8*X**2) &
     & + P*S*SQRT(S)*1.91075e-4_r8 &
     &  + P*P*(8.50935e-5_r8 - 6.12293e-6_r8*X + 5.2787e-8_r8*X**2) &
     &  + P*P*S*(-9.9348e-7_r8 + 2.0816e-8_r8*X + 9.1697e-10_r8*X**2)

! Density of seawater at S,T,P
  drho = rho0/(1.0d0 - P/Ksbm)
  rho_calc = REAL(drho)

  RETURN
END FUNCTION rho_calc

!>     Function to compute pressure from depth using Saunders (1981) formula with eos80.
FUNCTION p80(dpth,xlat) 

  !     Compute Pressure from depth using Saunders (1981) formula with eos80.

  !     Reference:
  !     Saunders, Peter M. (1981) Practical conversion of pressure
  !     to depth, J. Phys. Ooceanogr., 11, 573-574, (1981)

  !     Coded by:
  !     R. Millard
  !     March 9, 1983
  !     check value: p80=7500.004 dbars at lat=30 deg., depth=7321.45 meters

  !     Modified (slight format changes + added ref. details):
  !     J. Orr, 16 April 2009

  USE mod_kinds
  IMPLICIT NONE

! Input variables:
  !> depth [m]
  REAL(kind=r8), INTENT(in) :: dpth
  !> latitude [degrees]
  REAL(kind=r8), INTENT(in) :: xlat

! Output variable:
  !> pressure [db]
  REAL(kind=r8) :: p80

! Local variables:
  REAL(kind=r8) :: pi
  REAL(kind=r8) :: plat, d, c1

  pi=3.141592654

  plat = ABS(xlat*pi/180.)
  d  = SIN(plat)
  c1 = 5.92e-3+d**2 * 5.25e-3

  p80 = ((1-c1)-SQRT(((1-c1)**2)-(8.84e-6*dpth))) / 4.42e-6

  RETURN
END FUNCTION p80

!> Function to calculate potential temperature [C] from in-situ Temperature
FUNCTION sw_ptmp(s,t,p,pr)

  !     ==================================================================
  !     Calculates potential temperature [C] from in-situ Temperature [C]
  !     From UNESCO 1983 report.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================

  !     Input arguments:
  !     -------------------------------------
  !     s  = salinity            [psu      (PSS-78) ]
  !     t  = temperature         [degree C (IPTS-68)]
  !     p  = pressure            [db]
  !     pr = reference pressure  [db]

  USE mod_kinds
  IMPLICIT NONE

! Input arguments
  !> salinity [psu (PSS-78)]
  REAL(kind=r8) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=r8) :: t
  !> pressure [db]
  REAL(kind=r8) :: p
  !> reference pressure  [db]  
  REAL(kind=r8) :: pr

! local arguments
  REAL(kind=r8) :: del_P ,del_th, th, q
  REAL(kind=r8) :: onehalf, two, three
  PARAMETER (onehalf = 0.5d0, two = 2.d0, three = 3.d0 )

 ! REAL(kind=r8) :: sw_adtg
 ! EXTERNAL sw_adtg

! Output 
  REAL(kind=r8) :: sw_ptmp

  ! theta1
  del_P  = PR - P
  del_th = del_P*sw_adtg(S,T,P)
  th     = T + onehalf*del_th
  q      = del_th

  ! theta2
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 - 1.d0/SQRT(two))*(del_th - q)
  q      = (two-SQRT(two))*del_th + (-two+three/SQRT(two))*q

  ! theta3
  del_th = del_P*sw_adtg(S,th,P+onehalf*del_P)
  th     = th + (1.d0 + 1.d0/SQRT(two))*(del_th - q)
  q      = (two + SQRT(two))*del_th + (-two-three/SQRT(two))*q

  ! theta4
  del_th = del_P*sw_adtg(S,th,P+del_P)
  sw_ptmp     = th + (del_th - two*q)/(two*three)

  RETURN
END FUNCTION sw_ptmp

!> Function to compute in-situ temperature [C] from potential temperature [C]
FUNCTION sw_temp( s, t, p, pr )
  !     =============================================================
  !     SW_TEMP
  !     Computes in-situ temperature [C] from potential temperature [C]
  !     Routine available in seawater.f (used for MIT GCM)
  !     Downloaded seawater.f (on 17 April 2009) from
  !     http://ecco2.jpl.nasa.gov/data1/beaufort/MITgcm/bin/
  !     =============================================================

  !     REFERENCES:
  !     Fofonoff, P. and Millard, R.C. Jr
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !     Eqn.(31) p.39

  !     Bryden, H. 1973.
  !     "New Polynomials for thermal expansion, adiabatic temperature gradient
  !     and potential temperature of sea water."
  !     DEEP-SEA RES., 1973, Vol20,401-408.
  !     =============================================================

  !     Simple modifications: J. C. Orr, 16 April 2009
  !     - combined fortran code from MITgcm site & simplification in
  !       CSIRO code (matlab equivalent) from Phil Morgan

  USE mod_kinds
  IMPLICIT NONE

  !     Input arguments:
  !     -----------------------------------------------
  !     s  = salinity              [psu      (PSS-78) ]
  !     t  = potential temperature [degree C (IPTS-68)]
  !     p  = pressure              [db]
  !     pr = reference pressure    [db]

  !> salinity [psu (PSS-78)]
  REAL(kind=r8) ::   s
  !> potential temperature [degree C (IPTS-68)]
  REAL(kind=r8) ::   t
  !> pressure [db]
  REAL(kind=r8) ::   p
  !> reference pressure [db]
  REAL(kind=r8) ::   pr

  REAL(kind=r8) ::  ds, dt, dp, dpr
  REAL(kind=r8) :: dsw_temp

  REAL(kind=r8) ::   sw_temp
!  EXTERNAL sw_ptmp
!  REAL(kind=r8) ::   sw_ptmp

  ds = DBLE(s)
  dt = DBLE(t)
  dp = DBLE(p)
  dpr = DBLE(pr)

  !    Simple solution
  !    (see https://svn.mpl.ird.fr/us191/oceano/tags/V0/lib/matlab/seawater/sw_temp.m)
  !    Carry out inverse calculation by swapping P_ref (pr) and Pressure (p)
  !    in routine that is normally used to compute potential temp from temp
  dsw_temp = sw_ptmp(ds, dt, dpr, dp)
  sw_temp = REAL(dsw_temp)

  !    The above simplification works extremely well (compared to Table in 1983 report)
  !    whereas the sw_temp routine from MIT GCM site does not seem to work right

  RETURN
END FUNCTION sw_temp

!>  Function to calculate adiabatic temperature gradient as per UNESCO 1983 routines.
FUNCTION sw_adtg(s,t,p)

  !     ==================================================================
  !     Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================
  USE mod_kinds
  IMPLICIT NONE
  !> salinity [psu (PSU-78)]
  REAL(kind=r8) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=r8) :: t
  !> pressure [db]
  REAL(kind=r8) :: p

  REAL(kind=r8) :: a0,a1,a2,a3,b0,b1,c0,c1,c2,c3,d0,d1,e0,e1,e2
  REAL(kind=r8) :: sref

  REAL(kind=r8) :: sw_adtg

  sref = 35.d0
  a0 =  3.5803d-5
  a1 = +8.5258d-6
  a2 = -6.836d-8
  a3 =  6.6228d-10

  b0 = +1.8932d-6
  b1 = -4.2393d-8

  c0 = +1.8741d-8
  c1 = -6.7795d-10
  c2 = +8.733d-12
  c3 = -5.4481d-14

  d0 = -1.1351d-10
  d1 =  2.7759d-12

  e0 = -4.6206d-13
  e1 = +1.8676d-14
  e2 = -2.1687d-16

  sw_adtg =  a0 + (a1 + (a2 + a3*T)*T)*T &
     &       + (b0 + b1*T)*(S-sref) &
     &       + ( (c0 + (c1 + (c2 + c3*T)*T)*T) &
     &       + (d0 + d1*T)*(S-sref) )*P &
     &  + (  e0 + (e1 + e2*T)*T )*P*P
RETURN
END FUNCTION sw_adtg

