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
#define TS_DIF2
#define  MIX_GEO_TS
#define SALINITY
#define SOLVE3D

#undef AVERAGES
#undef DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV

/* LiveOcean specific choices */
#define NONLIN_EOS
#define MASKING
#define SPHERICAL
#define SOLAR_SOURCE
#define BULK_FLUXES
#define LONGWAVE_OUT
#define DEFLATE
#define HDF5
#define RADIATION_2D
#undef RAMP_TIDES
#define SSH_TIDES
#define UV_TIDES
#define ADD_FSOBC
#define ADD_M2OBC

/* LiveOcean bio choices */
#define BIO_FENNEL
#ifdef BIO_FENNEL
#  define BIO_SEDIMENT
#  define DENITRIFICATION
#  undef RIVER_DON
#  define OXYGEN
#  define CARBON
#  define pCO2_RZ
#  define TALK_NONCONSERV
#  define ANA_SPFLUX
#  define ANA_BPFLUX
#endif

#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX

#define GLS_MIXING
#if defined GLS_MIXING
# define CANUTO_A
# define N2S2_HORAVG
# define RI_SPLINES
#endif

#ifdef PERFECT_RESTART
# undef  AVERAGES
# undef  DIAGNOSTICS_BIO
# undef  DIAGNOSTICS_TS
# undef  DIAGNOSTICS_UV
# define OUT_DOUBLE
#endif
