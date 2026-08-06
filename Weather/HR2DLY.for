!=================================================================
!  HR2DLY.for, Subroutine, Fabio Oliveira
!  Calculates daily values of temperature, solar radiation,
!  and rain from hourly input file with FlexibleIO
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  10/06/2025  FO Written
!=======================================================================

      SUBROUTINE HR2DLY(YRDOY,                             !Input
     &    SRAD, TMAX, TMIN, RAIN, PAR)                      !Output

!-----------------------------------------------------------------------
      USE flexibleio
      IMPLICIT NONE
      INTEGER H,TINCR, YRDOY
      LOGICAL PARFND

      ! Daily (D) and (H) weather variables
      REAL HSRAD, SRADMJ, SRADJ, HTMAX, HTMIN, HRAIN, HPARV
      REAL DSRAD, DTMAX, DTMIN, DRAIN, DPAR

      ! Weather variables converted
      REAL SRAD, TMAX, TMIN, RAIN, PAR

      PARAMETER (TINCR=24)
!-----------------------------------------------------------------------
!     Initialize
      DSRAD  = 0.0
      DRAIN  = 0.0
      DPAR   = 0.0
      PARFND = .FALSE.

!     Loop to compute daily weather data over hourly data.
      DO H = 1, TINCR
        CALL fio % get('WTH', YRDOY, H, 'SRAD',HSRAD)
        CALL fio % get('WTH', YRDOY, H, 'TMAX',HTMAX)
        CALL fio % get('WTH', YRDOY, H, 'TMIN',HTMIN)
        CALL fio % get('WTH', YRDOY, H, 'RAIN',HRAIN)
        CALL fio % get('WTH', YRDOY, H, 'PAR', HPARV)

        !SRAD in hourly WTH file is in MJ/m2/h; sum directly for daily total
        SRADMJ = HSRAD
        !SRADJ must be in W/m2 (J/m2/s) for HMET/HPAR; convert MJ/m2/h -> W/m2
        SRADJ  = HSRAD * 1.0E6 / 3600.0
        CALL fio % set('WTH', YRDOY, H, 'SRADJ',SRADJ)

        DSRAD = DSRAD + SRADMJ
        DRAIN = DRAIN + HRAIN
        !PAR column is optional in the hourly WTH file (mol/m2/h, same
        !per-hour convention as SRAD in MJ/m2/h). Missing hours return
        !-99 (FlexibleIO sentinel); only sum hours where data was found.
        IF (HPARV .GT. -90.0) THEN
          DPAR   = DPAR + HPARV
          PARFND = .TRUE.
        ENDIF
        IF(H .EQ. 1) THEN
            DTMAX = HTMAX
            DTMIN = HTMIN
        ELSE
            DTMAX = MAX(DTMAX,HTMAX)
            DTMIN = MIN(DTMIN,HTMIN)
        ENDIF
      ENDDO

      !Assign the values converted
      SRAD = DSRAD
      TMAX = DTMAX
      TMIN = DTMIN
      RAIN = DRAIN
      !If no PAR column was present for this day, return -99 (sentinel)
      !so downstream HPAR falls back to its SRAD-based PAR estimate,
      !exactly as before this change.
      IF (PARFND) THEN
        PAR = DPAR
      ELSE
        PAR = -99.0
      ENDIF

      RETURN
      END SUBROUTINE HR2DLY