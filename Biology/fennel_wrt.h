/*
** svn $Id: npzdo_Banas_wrt.h 523 2011-01-05 03:21:38Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2011 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes npzdo_Banas ecosystem model input parameters into **
**  output NetCDF files. It is included in routine "wrt_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Write out NPZDO_Banas biological model parameters.
!

      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',              &
     &                      BioIter(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'PARfrac',              &
     &                      PARfrac(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'AttSW',                &
     &                      AttSW(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'AttP',                 &
     &                      AttP(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'AttS',                 &
     &                      AttS(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'phyAlpha',             &
     &                      phyAlpha(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'phyMu0',               &
     &                      phyMu0(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'phyKs',                &
     &                      phyKs(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'phyMin',               &
     &                      phyMin(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'phyM',                 &
     &                      phyM(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'zooI0',                &
     &                      zooI0(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'zooKs',                &
     &                      zooKs(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'zooEps',               &
     &                      zooEps(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'zooFegest',            &
     &                      zooFegest(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'zooMin',               &
     &                      zooMin(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'zooZeta',              &
     &                      zooZeta(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'detRemin',             &
     &                      detRemin(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'CoagR',                &
     &                      CoagR(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'detWsink',             &
     &                      detWsink(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'LdetWsink',            &
     &                      LdetWsink(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'pCO2air',              &
     &                      pCO2air(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
