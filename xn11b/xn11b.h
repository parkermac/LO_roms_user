/*
** svn $Id: upwelling.h 1054 2021-03-06 19:47:12Z arango $
*******************************************************************************
** Copyright (c) 2002-2021 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Upwelling Test.
**
** Application flag:   Capitalized name of this folder
** Input script:       liveocean.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define SALINITY
#define SOLVE3D

/* LiveOcean specific choices */
#define NONLIN_EOS
#define MASKING
#define SPHERICAL
#define SOLAR_SOURCE
#define BULK_FLUXES
#define LONGWAVE_OUT
#define EMINUSP
#define DEFLATE
#define HDF5
#undef AVERAGES
#define WET_DRY
#define OMEGA_IMPLICIT

#undef RADIATION_2D
#undef SSH_TIDES
#undef UV_TIDES
#undef ADD_FSOBC
#undef ADD_M2OBC

#define ANA_BTFLUX
#define ANA_BSFLUX

/* LiveOcean bio choices */
#define BIO_FENNEL
#ifdef BIO_FENNEL
#  define BIO_SEDIMENT
#  define DENITRIFICATION
#  undef RIVER_DON
#  define OXYGEN
#  define CARBON
#  define pCO2_RZ
#  define PCO2AIR_SECULAR
#  define TALK_NONCONSERV
#  define ANA_SPFLUX
#  define ANA_BPFLUX
#  define INCREASE_ATTSW  /*jx*/
#  undef INCREASE_NO3LOSS /*jx*/
#  define BURY           /*jx*/
#endif

#define GLS_MIXING
#if defined GLS_MIXING
# define CANUTO_A
# define N2S2_HORAVG
# define RI_SPLINES
#endif

#undef PERFECT_RESTART
#ifdef PERFECT_RESTART
# undef  AVERAGES
# undef  DIAGNOSTICS_BIO
# undef  DIAGNOSTICS_TS
# undef  DIAGNOSTICS_UV
# undef  OUT_DOUBLE
#endif
