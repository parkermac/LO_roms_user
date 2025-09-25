! Next two lines taken from ROMS/Nonlinear/Biology/red_tide.h, although in the
!   case of red_tide.h the NEUMANN option is set to #undef. Stock et al. 2005
!   implemented a gradual reduction in ascent velocity toward zero at the sea
!   surface. The closest analog in the ROMS implementation is to #define
!   NEUMANN (which I interpret as a no-flux condition at the sea surface and
!   seabed). Algal concentrations will thus accumulate in the top/bottom levels.
#include "cppdefs.h"
#define NEUMANN
#define LIMIT_INTERIOR
!
      MODULE mod_oae
#ifdef OAE_BGC
!
! Fortran 95 implementation of OAE module by K.E.Hinson (kyle.hinson@pnnl.gov)
!   with multiple customizations, notably the addition of secondary precipitation
!   of carbonates The format as well as some blocks of code are
!   borrowed from Rutgers' `fennel.h' to maximize compatibility with the Rutgers
!   branch of ROMS. This module is meant to be a `add-on' to an existing BGC
!   codebase and can be activated/deactivated at will with #define/undef OAE_BGC.
! It adds 2 new state variables (CO3/HCO3 if missing from TIC) and utilizes a sinking
!   particle "sp_crystal" that is saved between time steps but otherwise shares
!   the other state variables of the original BGC code
!  (inorganic nutrients, organic matter, etc).
! The module is self-contained inside the present file but a number of
!   modifications were made to ECB, notably the introduction of the 4 new state
!   variables inside ecb_mod.h (inside a #ifdef MARGALEX block).
!
! CAREFUL, you are currently testing individual components of the code
!   by commenting/uncommenting them!
! TO-DO:
!      --implementation of secondary precipitate sinking is preliminary,
!      --cysts at the bottom (not a state variable at the moment;
!          could be an infinite reservoir, or a 2-D field like tacc),
!      --...
!
      USE mod_kinds, only: r8
      implicit none
!
      PRIVATE
      PUBLIC oae
      PUBLIC fix_cff3


! DEFINE SIZE OF ROMS GRID, UNIQUE TO ALL DOMAINS
! This must be changed to match the size of eta_rho xi_rho, and s_rho in the user's
! grid file
      real(r8), save :: sp_crystal( 0:335 , 0:563 , 1:20 ) = 0.
!
!     We assume a grid dimension Lm=334, Mm=562, N=20, which is double-checked
!       at the beginning of precip_tile.
!       The fixed size is required to use the `save' attribute.
!     E.g., sp_crystal( - 2 : Lm+3, - 2 : Mm+3, [# vertical levels] ) = 0.
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE oae(ng, tile, myit)
!***********************************************************************
!
      USE mod_param
# ifdef DIAGNOSTICS_BIO
      USE mod_diags
# endif
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
      integer, intent(in) :: ng, tile, myit
!

!     real( r8 ), intent(  out) ::                                      &
!     &  aragss(  IminS : ImaxS        )

#include "tile.h"

!
      CALL oae_tile( ng, tile, myit,                                    &
     &                    LBi, UBi, LBj, UBj, N(ng), NT(ng),            &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng), nnew(ng),                           &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  ifdef WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
#  endif
# endif
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % z_r,                               &
     &                    GRID(ng) % z_w,                               &
     &                    FORCES(ng) % srflx,                           &
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % bvstr,                           &
     &                    GRID(ng) % h,                                 &
#ifdef CARBON
     &                    OCEAN(ng) % pH,                               &
#endif
     &                    OCEAN(ng) % t,                                &
     &                    sp_crystal(LBi:UBi,LBj:UBj,:) )
!
      RETURN
      END SUBROUTINE oae
!
!-----------------------------------------------------------------------
      SUBROUTINE oae_tile( ng, tile, myit,                              &
     &                          LBi, UBi, LBj, UBj, UBk, UBt,           &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp, nnew,                             &
# ifdef MASKING
     &                          rmask,                                  &
#  if defined WET_DRY
     &                          rmask_wet,                              &
#  endif
# endif
     &                          Hz, z_r, z_w, srflx,                    &
     &                          bustr, bvstr,                           &
     &                          h,                                      &
     &                          pH,                                     &
     &                          t,                                      &
     &                          sp_crystal_tile )
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      USE dateclock_mod, ONLY : caldate
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, myit
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(      LBi :, LBj : )
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(  LBi :, LBj : )
#   endif
#  endif
      real(r8), intent(in) :: Hz(      LBi :, LBj :, :       )
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: h(       LBi :, LBj :          )
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:,LBj:)
!      real(r8), intent(inout) :: TAlk(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: t( LBi : , LBj : , :, :, : )
      real(r8), intent(inout) :: sp_crystal_tile(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(      LBi : UBi, LBj : UBj )
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(  LBi : UBi, LBj : UBj )
#   endif
#  endif
      real(r8), intent(in) :: Hz(      LBi : UBi, LBj : UBj, UBk    )
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(       LBi : UBi, LBj : UBj         )
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: t( LBi : UBi, LBj : UBj, UBk, 3, UBt  )
      real(r8), intent(inout) :: sp_crystal_tile(LBi:UBi,LBj:UBj,UBk)
# endif
!
!  Local variable declarations.
!
      integer, parameter :: DoNewton = 0

      integer :: Iter, i, ibio, iswim, itrc, ivar, j, k, ks, year,      &
     &  hr_day_utc, jcla
      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6

#ifdef LOSTCO3
      real(r8) :: CO3( LBi : UBi, LBj : UBj, UBk )
!     real(r8), dimension(IminS:, ImaxS:, JminS:, JmaxS:, KminS:, KmaxS:) :: CO3
      real(r8) :: HCO3( LBi : UBi, LBj : UBj, UBk )
      real(r8) :: aragss( LBi : UBi, LBj : UBj, UBk )
!#else
!      real(r8), intent(inout) :: CO3(LBi:,LBj:)
!      real(r8), intent(inout) :: HCO3(LBi:,LBj:)
# endif

! Create new level-specific pH variable for aragonite calculations at all depth levels
      real(r8) :: pH_kLev( LBi : UBi, LBj : UBj )
      real(r8), dimension(IminS:ImaxS) :: pCO2

! Create new array to return redissolved (recycled) sp_crystal to CO3 pool 
      real(r8) :: REC_CO3( LBi : UBi, LBj : UBj, UBk )

!     Calcite/Aragonite saturation state
!      real(r8) :: cass, ksp_arag, ksp_calc, zed_factor, temp_factor

      real(r8) :: arag_threshold
      real(r8) :: diss_threshold
      integer :: mcloc
      real(r8), parameter :: spsink = 2.0            ! m d-1

#ifdef CARBON
      integer, parameter :: Nsink = 6
#else
      integer, parameter :: Nsink = 4
#endif

      integer :: isink

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

!      CASE ('idTvar(iSPC_)')
!        idTvar(iSPC_)=varid

!     Stand in cass values
!      real(r8), parameter :: calcite = 41.6  ! Z & W-G (2001)
!     real(r8), parameter :: aragonite = 63  ! Z & W-G (2001)
!     Calcium ion concentration in seawater - from  Z & W-G (2001)
!     real(r8), parameter :: ca_conc = 0.01028  ! mol kg-1
!     Sinking rate of secondary precipitate particle
!     Runaway precipitation surface area parameter
!     real(r8), parameter :: spagg = 1.2             ! dimensionless
      

      real(r8) :: fac1, fac2, fac3
      real(r8) :: cffL, cffR, cu, dltL, dltR


# ifdef DIAGNOSTICS_BIO
!      real(r8) :: fiter
# endif
      real(dp) :: yday
      real(r8) :: dtdays

      real(r8) :: alpha, deltaL, deltaR, dz, wdt

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

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
!
!     Double-check the grid dimensions for the grid of time-accumulated
!       temperature
!
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
# ifdef DIAGNOSTICS_BIO
!
!  A factor to account for the number of iterations in accumulating
!  diagnostic rate variables.
!
!      fiter=1.0_r8/REAL(BioIter(ng),r8)
# endif

!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
      J_LOOP : DO j=Jstr,Jend
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              Bio(i,k,ibio)=Bio_old(i,k,ibio)
            END DO
          END DO
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO

#ifdef OAE_BGC
!
!-----------------------------------------------------------------------
!  Calculate aragonite saturation state based on pH and other factors
!-----------------------------------------------------------------------

        DO k=1,N(ng)
!         First, recompute pH to get the value at each vertical level
          CALL pCO2_water_RZ (Istr, Iend, LBi, UBi, LBj, UBj,           &
     &                        IminS, ImaxS, j, DoNewton,                &
     &                        rmask,                                    &
     &                        Bio(IminS:,k,itemp), Bio(IminS:,k,isalt), &
     &                        Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk), &
     &                        pH_kLev, pCO2)

!         write(*,*) 'min/max pH_kLev,j,k=',minval(pH_kLev),maxval(pH_kLev),j,k
!         IF (Bio(15,k,iTAlk) > 3000.0_r8) THEN
!           write(*,*) 'alk(15,j,k),j,k=',Bio(15,k,iTAlk),j,k
!           write(*,*) 'pH(15,j),j,k=',pH_kLev(15,j),j,k
!         END IF
!         IF (t(15,j,k,nstp,itemp) > 0.1_r8) THEN
!         IF (Bio(15,k,itemp) > 0.1_r8) THEN
!           write(*,*) 'pH_temp(15,j,k)=',Bio(15,k,itemp),j,k
!           write(*,*) 'pH_salt(15,j,k)=',Bio(15,k,isalt),j,k
!           write(*,*) 'pH_TIC(15,j,k)=',Bio(15,k,iTIC_),j,k
!           write(*,*) 'pH_TAlk(15,j,k)=',Bio(15,k,iTAlk),j,k
!           write(*,*) 'pH_kLev(15,j,k)=',pH_kLev(15,j),j,k
!         END IF

          CALL calc_CO3 (                                               &
     &           Istr, Iend, LBi, UBi, LBj, UBj,                        &
     &           IminS, ImaxS, j, rmask,                                &
!     &           h,                                                     &
     &           z_r(IminS:,j,k),                                       &
     &           Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),              &
     &           Bio(IminS:,k,iTIC_),                                   &
     &           pH_kLev, CO3(IminS:,j,k) )
!         write(*,*) 'minCO3(15,j,k)=',minval(CO3(:,j,k)),j,k
!         write(*,*) 'maxCO3(15,j,k)=',maxval(CO3(:,j,k)),j,k

!         Compute aragonite saturation state from pH,temp,salt,TIC,pressure:
!       DO k=1,N(ng)
          call aragonite_saturation_state(                              &
     &           Istr, Iend, LBi, UBi, LBj, UBj,                        &
     &           IminS, ImaxS, j, rmask,                                &
!    &           h,                                                     &
     &           z_r(IminS:,j,k),                                       &
     &           Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),              &
     &           Bio(IminS:,k,iTIC_),                                   &
     &           pH_kLev, aragss(IminS:,j,k) )

!          IF (maxval(Bio(IminS:,k,iTAlk)) > 5000.0_r8) THEN
!             mcloc = maxloc(Bio(IminS:,k,iTAlk),DIM=1)
!             write(*,*) 'pH_temp(:,j,k)=',Bio(:,k,itemp),j,k
!             write(*,*) 'pH_salt(:,j,k)=',Bio(:,k,isalt),j,k
!             write(*,*) 'pH_depth(:,j)=',z_r(:,j,k),j,k
!             write(*,*) 'pH_TIC(:,j,k)=',Bio(:,k,iTIC_),j,k
!             write(*,*) 'pH_TAlk(:,j,k)=',Bio(:,k,iTAlk),j,k
!             write(*,*) 'pH_kLev(:,j,k)=',pH_kLev(:,j),j,k
!             write(*,*) 'pH_CO3(:,j,k)=',CO3(:,j,k),j,k
!             write(*,*) 'aragss(:,j,k),j,k=',aragss(:,j,k),j,k
!!            write(*,*) 'min/max aragss,j,k=',minval(aragss(:,j,k)),maxval(aragss(:,j,k)),j,k
!           END IF
        END DO

        ! pH(IminS:,JminS)
!
! Determine if threshold of aragonite_saturation_state exceeds 7 (chosen
! based on reported values in Moras et al., 2022), and if so then begin
! secondary precipitation of calcium carbonate. This operation is a sink
! for CO3, and transforms 50 percent of the available CO3 in a given cell
! to a sinking particle that can redissolve under favorable conditions.
! These conditions will be met when aragonite_saturation_state reaches a
! value of 2, where the carbonate ion in the precipitate crystal will be
! returned to the dissolved CO3 pool.

        arag_threshold = 7.0_r8
        diss_threshold = 1.0_r8
        do i = Istr, Iend
          do k = 1, N(ng)
            IF (aragss(i,j,k) > arag_threshold) THEN
              write(*,*) 'Here is i,j,k,CO3(i,j,k)=',CO3(i,j,k),i,j,k
              sp_crystal_tile(i,j,k) = CO3(i,j,k) * 0.5_r8
              sp_crystal_tile(i,j,k) = max( sp_crystal_tile(i,j,k), 0.0_r8 )
              IF (CO3(i,j,k) < 0 ) THEN
                write(*,*) 'Something bad happened, negative CO3'
              END IF
              write(*,*) 'Here is i,j,k,aragss(i,j,k)=',aragss(i,j,k),i,j,k
              write(*,*) 'Threshold met,i,j,k,sp_crystal_tile(i,j,k)=',i,j,k,sp_crystal_tile(i,j,k)
              CO3(i,j,k) = CO3(i,j,k) * 0.5_r8
              write(*,*) 'Here is updated i,j,k,CO3(i,j,k)=',CO3(i,j,k),i,j,k
              REC_CO3(i,j,k) = 0.0_r8
            ELSE IF (aragss(i,j,k) <= diss_threshold) THEN
              CO3(i,j,k) = CO3(i,j,k) + sp_crystal_tile(i,j,k)
              REC_CO3(i,j,k) = sp_crystal_tile(i,j,k)
              sp_crystal_tile(i,j,k) = 0.0_r8             
            ELSE
              sp_crystal_tile(i,j,k) = 0.0_r8 ! neither precipitating nor dissolving
              REC_CO3(i,j,k) = 0.0_r8
            END IF
          end do
        end do

! Now we sink our secondary precipitate crystal


!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=sp_crystal_tile(i,j,k)
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
            cff=dtdays*ABS(spsink)
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
                sp_crystal_tile(i,j,k)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO


!       END DO ITER_LOOP
! end if statement for OAE flag
#endif
!      end do J_LOOP


!
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
!       DO itrc=1,NBT
!          ibio=idbio(itrc)
!          DO k=1,N(ng)
!            DO i=Istr,Iend
!              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
!# ifdef MASKING
!              cff=cff*rmask(i,j)
!#  ifdef WET_DRY
!              cff=cff*rmask_wet(i,j)
!#  endif
!# endif
!              t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
!            END DO
!          END DO
!        END DO
        DO i=Istr,Iend
          DO k=1,N(ng)
!            t(i,j,k,nnew,iTIC_) = -9999.0_r8
            t(i,j,k,nnew,iTIC_) = t(i,j,k,nnew,iTIC_) - sp_crystal_tile(i,j,k) * Hz(i,j,k)            & 
           &                        + REC_CO3(i,j,k) * Hz(i,j,k)
!           t(i,j,k,nnew,iTAlk) = -9999.0_r8
!           write(*,*) 'Original t(i,j,k,nnew,iTAlk),i,j,k=',t(i,j,k,nnew,iTAlk),i,j,k
           t(i,j,k,nnew,iTAlk) = t(i,j,k,nnew,iTAlk) - (sp_crystal_tile(i,j,k) * 2.0_r8 * Hz(i,j,k)) & 
            &                        + (REC_CO3(i,j,k) * 2.0_r8 * Hz(i,j,k))
!           write(*,*) 'Updated t(i,j,k,nnew,iTAlk),i,j,k=',t(i,j,k,nnew,iTAlk),i,j,k
          END DO
        END DO


      END DO J_LOOP


      RETURN
      end SUBROUTINE oae_tile

!      SUBROUTINE SEC_PREC(TAlk,pH)
!              real(r8), intent(inout) :: TAlk
!              real(r8), intent(in) :: pH
!
!              IF (pH>=9.75) THEN
!                TAlk=0.8575_r8*TAlk
!              ELSE
!                TAlk=TAlk
!             END IF
!     end SUBROUTINE SEC_PREC

      subroutine calc_CO3(                                              &
     &       Istr, Iend, LBi, UBi, LBj, UBj, IminS, ImaxS,              &
     &       j, rmask,                                                  &
     &       zvec, TempC, Sal, TIC, pH, CO3 )
      implicit none
!     Check value: DIC   = 1474.29 mmol-C/m3 <=> 1474.29 umol-C/kgSW,
!                  pH_sws= 8.62941, salinity= 19.9825, temp= 6.8575dec.C,
!                  press = 5 dbars, Omega_calcite = 5.443654.
!
      integer,    intent(in   ) ::                                      &
     &  Istr, Iend, LBi, UBi, LBj, UBj, IminS, ImaxS, j
!
      real( r8 ), intent(in   ) ::                                      &
     &  rmask( LBi : UBi, LBj : UBj ),                                  &
!     &  hbot(  LBi : UBi, LBj : UBj ),                                  &
     &  zvec(  IminS : ImaxS        ),                                  &
     &  TempC( IminS : ImaxS        ),                                  &
     &  Sal(   IminS : ImaxS        ),                                  &
     &  TIC(   IminS : ImaxS        ),                                  &
     &  pH(    LBi : UBi, LBj : UBj )

      real( r8 ), intent(  out) ::                                      &
     &  CO3(  IminS : ImaxS        )
!
      real( r8 ) ::                                                     &
     &  TempK, logTempK, sqrSal,                                        &
     &  H, Pbar, RT, F1, F2, K1,                                        &
     &  K2, pK1, pK2, fH, TC, deltaV, Kappa, lnK1fac, lnK2fac,          &
     &  K1fac, K2fac

      integer    :: i

      I_LOOP: do i = Istr, Iend
        if ( rmask(i, j) .le. tiny(0._r8) ) then ! Dry cell.
          CO3( i ) = 0._r8
          continue
        end if
        TempK    = TempC( i ) + 273.15_r8 ! Convert Celsius to Kelvin
        logTempK = log( TempK )
        sqrSal   = sqrt( Sal( i ) )
!
!       CO2SYS uses as input the hydrostatic pressure (zero at sea surface).
!         Here we neglect differences between the sea floor depth and the
!         position of the rho point, we assume an undisturbed sea surface, and
!        we assume that 1meter of seawater equates to a pressure of 1decibar.
!       Note that hbot is bounded to zero in case that WET_DRY is activated.
!        Pbar      = max(hbot(i, j), 0._r8) * 1.e-1_r8 ! In bars.
        Pbar      = max(zvec(i), 0._r8) * 1.e-1_r8 ! In bars.
!       Assume an ideal gas constant `R' of 83.1451 mL bar-1 K-1 mol-1.
        RT        = 83.1451_r8 * TempK
!
        F1      = 200.1_r8   / TempK + 0.3220_r8
        pK1     = 3404.71_r8 / TempK + 0.032786_r8 * TempK -14.8435_r8  &
     &          - 0.071692_r8  * F1 * sqrt(Sal( i ))                    &
     &          + 0.0021487_r8 * Sal( i )
        fH      = 1.2948_r8 - 0.002036_r8 * TempK                       &
     &          + (0.0004607_r8 - 0.000001475_r8 * TempK) * Sal( i )**2
        K1      = 10._r8**( - pK1 ) / fH
        F2      = - 129.24_r8 / TempK + 1.4381_r8
        pK2     = 2902.39_r8  / TempK + 0.02379_r8 * TempK - 6.4980_r8  &
     &          - 0.3191_r8 * F2 * sqrt(Sal(i)) + 0.0198_r8 * Sal(i)
        K2      = 10._r8**( - pK2      ) / fH
!
        deltaV  =  - 25.5_r8 + 0.1271_r8 * TempC(i)
        Kappa   = (- 3.08_r8 + 0.0877_r8 * TempC(i)) * 1.e-3_r8
        lnK1fac = (- deltaV + 0.5_r8 * Kappa * Pbar) * Pbar / RT
!
        deltaV  = - 15.82_r8 - 0.0219_r8 * TempC(i)
        Kappa   = (  1.13_r8 - 0.1475_r8 * TempC(i)) * 1.e-3_r8
        lnK2fac = (- deltaV + 0.5_r8 * Kappa * Pbar) * Pbar / RT
        K1fac   = exp( lnK1fac )
        K2fac   = exp( lnK2fac )
        K1      = K1 * K1fac
        K2      = K2 * K2fac

        H       = 10._r8**( - pH(i, j) )
!       Convert DIC from mmol-C/m3 to mol-C/kg-seawater (approximately):
        TC      = TIC( i ) * 1.e-6_r8
        CO3(i)  = (TC * K1 * K2 / (K1 * H + H**2 + K1 * K2)) * 1.e6_r8
        CO3(i)  = max( CO3(i), 0.0_r8 )
        IF ( CO3(i) > 10000._r8 ) THEN
          write(*,*) 'Here is TC=',TC,i,j
          write(*,*) 'Here is pH=',pH(i,j)
          write(*,*) 'Here is TempC=',TempC(i)
          write(*,*) 'Here is Sal=',Sal(i)
          write(*,*) 'Here is CO3=',CO3(i)
!              write(*,*) 'Here is i,j,k,salt(i,j,k)=',t(i,j,k,nstp,isalt),i,j,k
!              write(*,*) 'Here is i,j,k,TIC(i,j,k)=',t(i,j,k,nstp,iTIC_),i,j,k
!              write(*,*) 'Here is i,j,k,pH_kLev(i,j)=',pH_kLev(i,j),i,j
!              write(*,*) 'Here is i,j,k,h(i,j,k)=',h(i,j),i,j
!              write(*,*) 'Here is i,j,k,aragss(i,j,k)=',aragss(i,j,k),i,j,k
        END IF
      end do I_LOOP
!     write(*,*) 'maxCO3=',maxval(CO3)
      end subroutine calc_CO3


      
! Calculate aragonite saturation state
      subroutine aragonite_saturation_state(                            &
     &       Istr, Iend, LBi, UBi, LBj, UBj, IminS, ImaxS,              &
     &       j, rmask,                                                  &
     &       zvec, TempC, Sal, TIC, pH, aragValue)
      implicit none
!     This subroutine is based on CO2SYS.m version 2.0 (20 Dec.2016). The
!       implementation uses the dissociation constants of Cai & Wang 1998 and
!       pH as input (assumed to be on the seawater scale). The subroutine also
!       takes DIC (in millimoles-C/m3) as input and assumes for simplicity that
!       1 mmol-C/m3 <=> 1 umol-C/kg-seawater (i.e., we accept the <3% error).
!     If used with double precision arithmetic (r8=dp), the subroutine produces
!       the same calcite saturation state as the aforementioned CO2SYS code.
!     Check value: DIC   = 1474.29 mmol-C/m3 <=> 1474.29 umol-C/kgSW,
!                  pH_sws= 8.62941, salinity= 19.9825, temp= 6.8575dec.C,
!                  press = 5 dbars, Omega_calcite = 5.443654.
!
      integer,    intent(in   ) ::                                      &
     &  Istr, Iend, LBi, UBi, LBj, UBj, IminS, ImaxS, j
!
      real( r8 ), intent(in   ) ::                                      &
     &  rmask( LBi : UBi, LBj : UBj ),                                  &
!    &  hbot(  LBi : UBi, LBj : UBj ),                                  &
     &  zvec(  IminS : ImaxS        ),                                  &
     &  TempC( IminS : ImaxS        ),                                  &
     &  Sal(   IminS : ImaxS        ),                                  &
     &  TIC(   IminS : ImaxS        ),                                  &
     &  pH(    LBi : UBi, LBj : UBj )
!
      real( r8 ), intent(  out) ::                                      &
     &  aragValue(  IminS : ImaxS        )
!
      real( r8 ) ::                                                     &
     &  Ca, TempK, logTempK, sqrSal, logKAr, KAr, deltaVKCa, KappaKCa,  &
     &  deltaVKAr, KappaKAr, lnKCafac, H, lCO3, Pbar, RT, F1, F2, K1,   &
     &  K2, pK1, pK2, fH, TC, deltaV, Kappa, lnK1fac, lnK2fac,          &
     &  K1fac, K2fac, lnKArfac
      integer    :: i
!
      I_LOOP: do i = Istr, Iend
        if ( rmask(i, j) .le. tiny(0._r8) ) then ! Dry cell.
          aragValue( i ) = 0._r8
          continue
        end if
!
        Ca       = 0.02128_r8 / 40.087_r8  * (Sal( i ) / 1.80655_r8)
!
        TempK    = TempC( i ) + 273.15_r8 ! Convert Celsius to Kelvin
        logTempK = log( TempK )
        logKAr   = - 171.945_r8 - 0.077993_r8 * TempK                   &
     &                           + 2903.293_r8 / TempK
        sqrSal   = sqrt( Sal( i ) )
        logKAr   = logKAr + 71.595_r8 * logTempK / log( 10._r8 )
        logKAr   = logKAr + ( - 0.068393_r8 + 0.0017276_r8 * TempK      &
     &                                     + 88.135_r8    / TempK )     &
     &           * sqrSal
        logKAr   = logKAr - 0.10018_r8   * Sal( i )                     &
     &           + 0.0059415_r8 * sqrSal * Sal( i )
        KAr      = 10._r8**logKAr ! (mol/kg-SW)^2

!       CO2SYS uses as input the hydrostatic pressure (zero at sea surface).
!         Here we neglect differences between the sea floor depth and the
!         position of the rho point, we assume an undisturbed sea surface, and
!        we assume that 1meter of seawater equates to a pressure of 1decibar.
!       Note that hbot is bounded to zero in case that WET_DRY is activated.
!        Pbar      = max(hbot(i, j), 0._r8) * 1.e-1_r8 ! In bars.
        Pbar      = max(zvec(i), 0._r8) * 1.e-1_r8 ! In bars.
!       Assume an ideal gas constant `R' of 83.1451 mL bar-1 K-1 mol-1.
        RT        = 83.1451_r8 * TempK
        deltaVKCa = - 48.76_r8 + 0.5304_r8 * TempC( i ) ! New bug found?
        deltaVKAr = deltaVKCa + 2.8_r8
        KappaKAr  = (- 11.76_r8 + 0.3692_r8 * TempC( i )) * 1.e-3_r8
!       lnKCafac  = (- deltaVKCa + 0.5 * KappaKCa * Pbar) * Pbar / RT
        lnKArfac  = (- deltaVKAr + 0.5_r8 * KappaKAr * Pbar) * Pbar / RT
        KAr       = KAr * exp( lnKArfac )
!
        F1      = 200.1_r8   / TempK + 0.3220_r8
        pK1     = 3404.71_r8 / TempK + 0.032786_r8 * TempK -14.8435_r8  &
     &          - 0.071692_r8  * F1 * sqrt(Sal( i ))                    &
     &          + 0.0021487_r8 * Sal( i )
        fH      = 1.2948_r8 - 0.002036_r8 * TempK                       &
     &          + (0.0004607_r8 - 0.000001475_r8 * TempK) * Sal( i )**2
        K1      = 10._r8**( - pK1 ) / fH
        F2      = - 129.24_r8 / TempK + 1.4381_r8
        pK2     = 2902.39_r8  / TempK + 0.02379_r8 * TempK - 6.4980_r8  &
     &          - 0.3191_r8 * F2 * sqrt(Sal(i)) + 0.0198_r8 * Sal(i)
        K2      = 10._r8**( - pK2      ) / fH
!
        deltaV  =  - 25.5_r8 + 0.1271_r8 * TempC(i)
        Kappa   = (- 3.08_r8 + 0.0877_r8 * TempC(i)) * 1.e-3_r8
        lnK1fac = (- deltaV + 0.5_r8 * Kappa * Pbar) * Pbar / RT
!
        deltaV  = - 15.82_r8 - 0.0219_r8 * TempC(i)
        Kappa   = (  1.13_r8 - 0.1475_r8 * TempC(i)) * 1.e-3_r8
        lnK2fac = (- deltaV + 0.5_r8 * Kappa * Pbar) * Pbar / RT
        K1fac   = exp( lnK1fac )
        K2fac   = exp( lnK2fac )
        K1      = K1 * K1fac
        K2      = K2 * K2fac
!
        H       = 10._r8**( - pH(i, j) )
!       Convert DIC from mmol-C/m3 to mol-C/kg-seawater (approximately):
        TC      = TIC( i ) * 1.e-6_r8
        lCO3     = TC * K1 * K2 / (K1 * H + H**2 + K1 * K2)
!       Compute ARAGonite Saturation State (aragValue):
        aragValue(i) = lCO3 * Ca / KAr
      end do I_LOOP
      end subroutine aragonite_saturation_state



      subroutine pCO2_water_RZ( Istr, Iend,                             &
     &                    LBi, UBi, LBj, UBj, IminS, ImaxS,             &
     &                    j, DoNewton,                                  &
#ifdef MASKING
     &                    rmask,                                        &
#endif
     &                    T, S, TIC, TAlk, pH_kLev, pCO2)
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
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(inout) :: pH_kLev(LBi:UBi,LBj:UBj)
!
      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30
!
      integer :: Hstep, Ibrack, Inewton, i
!
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
      real(r8) :: pkcw, fhcw
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
        IF (rmask(i,j).gt.0.0_r8) THEN
          Tk=T(i)+273.15_r8
          centiTk=0.01_r8*Tk
          invTk=1.0_r8/Tk
          logTk=LOG(Tk)
          sqrtS=SQRT(S(i))
          scl=S(i)/1.80655_r8
!
          alk= TAlk(i)*0.000001_r8
          dic = TIC(i)*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
          ff=EXP(-162.8301_r8+                                          &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carbonic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  Use K1,K2 of Cai and Wang, 1998, Limnology and Oceanography, p. 661.
!  The coefficients are on NBS scale, need to use "fhcw" to convert to SWS.
!  "fhcw" calculation is copied from CO2SYS.m (Takahashi et al, 1982)
!-----------------------------------------------------------------------
!
          fhcw =  1.2948_r8   - 2.036e-3_r8 * Tk                        &
     &         + (4.607e-4_r8 - 1.475e-6_r8 * Tk) * S(i)**2
!
          pkcw = 200.1_r8     * invTk + 0.322_r8
          pkcw = 3404.71_r8   * invTk + 0.032786_r8 * Tk - 14.8435_r8   &
     &         - 7.1692e-2_r8 * pkcw  * sqrtS       + 2.1487e-3_r8 *S(i)
          K1   = 10._r8**(- pkcw) ! This is on NBS scale
          K1   = K1 / fhcw        ! Convert to SWS scale
!
          pkcw = - 129.24_r8 * invTk + 1.4381_r8
          pkcw = 2902.39_r8  * invTk + 0.02379_r8 * Tk - 6.498_r8       &
     &         - 0.3191_r8   * pkcw  * sqrtS      + 0.0198_r8 * S(i)
          K2   = 10._r8**(- pkcw) ! This is on NBS scale
          K2   = K2 / fhcw        ! Convert to SWS scale
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
          Kb=EXP(-invTk*(8966.90_r8+                                    &
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
          Kw=EXP(148.9652_r8-                                           &
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
!
          p5 = -1.0_r8;
          p4 = -alk-Kb-K1;
          p3 = dic*K1-alk*(Kb+K1)+Kb*borate+Kw-Kb*K1-K12
          p2 = dic*(Kb*K1+2*K12)-alk*(Kb*K1+K12)+Kb*borate*K1           &
     &       +(Kw*Kb+Kw*K1-Kb*K12)
          p1 = 2.0_r8*dic*Kb*K12-alk*Kb*K12+Kb*borate*K12               &
     &       +Kw*Kb*K1+Kw*K12
          p0 = Kw*Kb*K12;
!
!  Set first guess and brackets for [H+] solvers.
!
          pH_guess=pH_kLev(i,j)    ! Newton-Raphson
!         pH_hi=10.0_r8            ! high bracket/bisection
          pH_hi=14.0_r8            ! high bracket/bisection
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
!             IF ( abs( fni(3) ) < 1.e-37_r8 ) THEN
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
          pH_kLev(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
          pCO2(i)=CO2star*1000000.0_r8/ff
        ELSE
          pH_kLev(i,j)=0.0_r8
          pCO2(i)=0.0_r8
        END IF
      END DO I_LOOP
      end subroutine pCO2_water_RZ

!***********************************************************************
! Initiate subroutine to modify air-sea gas exchange rates
!***********************************************************************
!
      SUBROUTINE FIX_CFF3(cff2p5,cff3,pH,dtdays)
              real(r8), intent(in) :: cff2p5
              real(r8), intent(inout) :: cff3
              real(r8), intent(in) :: pH
              real(r8), intent(in) :: dtdays

               IF (pH>=9) THEN
!                cff3=MAX(cff2p5,0.00111_r8) ! where to modify CO2 flux KH
                 cff3=MAX(cff2p5,4._r8*24._r8/100._r8*dtdays) ! where to modify CO2 flux KH
               ELSE
                 cff3=cff2p5
              END IF
!             cff3=MAX(cff2p5,0.00111_r8) ! where to modify CO2 flux KH
      end SUBROUTINE FIX_CFF3

#endif
      end MODULE mod_oae

