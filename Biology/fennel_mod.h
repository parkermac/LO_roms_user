!npzdo_Banas_mod.h created by K. Davis 2011-03-31, based on fennel_mod.h
!svn $Id: fennel_mod.h 556 2011-04-26 22:34:37Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
 USE mod_param

 implicit none
 
!  Set biological tracer identification indices.
!
  integer, allocatable :: idbio(:)  ! Biological tracers
  integer :: iNO3_                  ! Nitrate concentration
  integer :: iPhyt                  ! Phytoplankton concentration
  integer :: iZoop                  ! Zooplankton concentration
  integer :: iDetr                  ! Detritus concentration
  integer :: iLDetr                  ! Detritus concentration
  integer :: iOxyg				  ! Oxygen Concentration
  integer :: iTIC_				  ! DIC Concentration
  integer :: iTAlk				  ! ALK Concentration
  integer :: iCaCO3				  ! CaCO3 Concentration

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Biological 2D diagnostic variable IDs.
 
      integer, allocatable :: iDbio2(:)       ! 2D biological terms
 
      integer  :: iCOfx                       ! air-sea CO2 flux
!      integer  :: iDNIT                       ! denitrification flux 
      integer  :: ipCO2                       ! partial pressure of CO2
integer  :: iO2fx                        ! air-sea O2 flux
integer  :: iDsed                        ! sed remin flux
!
!  Biological 3D diagnostic variable IDs.
!
      integer, allocatable :: iDbio3(:)       ! 3D biological terms
 
      integer  :: iPPro = 1                   ! primary productivity
      integer  :: iRem = 2                   ! remineralization
      integer  :: iGraz = 3                   ! total grazing
      integer  :: iDNIT                      ! denitrification

#endif

!
!  Biological parameters.
!
 integer, allocatable :: BioIter(:)
 
 real(r8), allocatable :: PARfrac(:)
 real(r8), allocatable :: AttSW(:)
 real(r8), allocatable :: AttP(:)
 real(r8), allocatable :: AttS(:)
 real(r8), allocatable :: phyAlpha(:)
 real(r8), allocatable :: phyMu0(:)
 real(r8), allocatable :: phyKs(:)
 real(r8), allocatable :: phyMin(:)
 real(r8), allocatable :: phyM(:)
 real(r8), allocatable :: zooI0(:)
 real(r8), allocatable :: zooKs(:)
 real(r8), allocatable :: zooEps(:)
 real(r8), allocatable :: zooFegest(:)
 real(r8), allocatable :: zooMin(:)
 real(r8), allocatable :: zooZeta(:)
 real(r8), allocatable :: detRemin(:)
 real(r8), allocatable :: CoagR(:)
 real(r8), allocatable :: detWsink(:)
 real(r8), allocatable :: LdetWsink(:)
 real(r8), allocatable :: pCO2air(:)            ! ppmv


	CONTAINS

	SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
	integer :: i, ic
!
!-----------------------------------------------------------------------
!  Set number of biological tracers.
!-----------------------------------------------------------------------
!
#ifdef CARBON
# ifdef OXYGEN
      NBT=9
# else
      NBT=8 
!with CaCO3 this is 8
# endif
#else
# ifdef OXYGEN
      NBT=6
# else
      NBT=5
# endif
#endif
!!	NBT=5

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
      NDbio3d=3
      NDbio2d=0
# ifdef OXYGEN
      NDbio2d=NDbio2d+2
      NDbio3d=NDbio3d+1
# endif
#ifdef CARBON
      NDbio2d=NDbio2d+4
      NDbio3d=NDbio3d+1
#endif
!
!  Allocate biological diagnostics vectors
!
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
      END IF
      IF (.not.allocated(iDbio3)) THEN
        allocate ( iDbio3(NDbio3d) )
      END IF
!
!  Initialize biology diagnostic indices.
!
      ic=0
# ifdef OXYGEN
      iO2fx=ic+1
      iDsed=ic+2
      iDNIT=ic+4
# endif
#ifdef CARBON
       iCOfx =ic+5
       ipCO2 =ic+6
#endif
#endif

!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
	IF (.not.allocated(BioIter)) THEN
	 allocate ( BioIter(Ngrids) )
	END IF
	IF (.not.allocated(PARfrac)) THEN
	 allocate ( PARfrac(Ngrids) )
	END IF
	IF (.not.allocated(AttSW)) THEN
	 allocate ( AttSW(Ngrids) )
	END IF
	IF (.not.allocated(AttP)) THEN
	 allocate ( AttP(Ngrids) )
	END IF
	IF (.not.allocated(AttS)) THEN
	 allocate ( AttS(Ngrids) )
	END IF
	IF (.not.allocated(phyAlpha)) THEN
	 allocate ( phyAlpha(Ngrids) )
	END IF
	IF (.not.allocated(phyMu0)) THEN
	 allocate ( phyMu0(Ngrids) )
	END IF
	IF (.not.allocated(phyKs)) THEN
	 allocate ( phyKs(Ngrids) )
	END IF
	IF (.not.allocated(phyMin)) THEN
	 allocate ( phyMin(Ngrids) )
	END IF
	IF (.not.allocated(phyM)) THEN
	 allocate ( phyM(Ngrids) )
	END IF
	IF (.not.allocated(zooI0)) THEN
	 allocate ( zooI0(Ngrids) )
	END IF
	IF (.not.allocated(zooKs)) THEN
	 allocate ( zooKs(Ngrids) )
	END IF
	IF (.not.allocated(zooEps)) THEN
	 allocate ( zooEps(Ngrids) )
	END IF
	IF (.not.allocated(zooFegest)) THEN
	 allocate ( zooFegest(Ngrids) )
	END IF
	IF (.not.allocated(zooMin)) THEN
	 allocate ( zooMin(Ngrids) )
	END IF
	IF (.not.allocated(zooZeta)) THEN
	 allocate ( zooZeta(Ngrids) )
	END IF
	IF (.not.allocated(detRemin)) THEN
	 allocate ( detRemin(Ngrids) )
	END IF
	IF (.not.allocated(CoagR)) THEN
	 allocate ( CoagR(Ngrids) )
	END IF
	IF (.not.allocated(detWsink)) THEN
	 allocate ( detWsink(Ngrids) )
	END IF
	IF (.not.allocated(LdetWsink)) THEN
	 allocate ( LdetWsink(Ngrids) )
	END IF
        IF (.not.allocated(pCO2air)) THEN
         allocate ( pCO2air(Ngrids) )
        END IF
!

!
!  Allocate biological tracer vector.
!
	IF (.not.allocated(idbio)) THEN
	allocate ( idbio(NBT) )
	END IF
!
!  Set identification indices.
!
	ic=NAT+NPT+NCS+NNS
	DO i=1,NBT
	idbio(i)=ic+i
	END DO
	iNO3_=ic+1
	iPhyt=ic+2
	iZoop=ic+3
	iDetr=ic+4
        iLDetr=ic+5
	ic=ic+5
#ifdef OXYGEN
	iOxyg=ic+1
	ic=ic+1
#endif
#ifdef CARBON
	iTIC_=ic+1
	iTAlk=ic+2


	iCaCO3=ic+3
        ic=ic+3

#endif

	RETURN
	END SUBROUTINE initialize_biology

