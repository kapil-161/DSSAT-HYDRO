C=======================================================================
C  HYDROPONICS, Subroutine
C
C  Main driver of the hydroponic solution module. Called from LAND in
C  every DYNAMIC phase, after SOIL and before SPAM, in the same way that
C  SOIL drives the soil water and soil nutrient routines. It calls the
C  routines that describe the nutrient solution:
C     HYDRO_FERT   scheduled nutrient additions
C     HYDRO_WATER  solution volume (water balance)
C     SOLEC        electrical conductivity and EC stress factors
C     SOLPH        solution pH and pH stress factors
C     SOLO2        dissolved oxygen and O2 stress factor
C     OPSOL        daily output file Solution.OUT
C
C  Nutrient uptake from the solution (HYDRO_NUTRIENT, SOLPi, SOLKi) is
C  called by the plant N, P and K uptake routines, in the same way that
C  soil N uptake is computed in NUPTAK from the soil N pools.
C
C  Order of calls within each phase reproduces the previous sequence
C  inside SPAM, so results are unchanged.
C-----------------------------------------------------------------------
C  Revision history
C  11/09/2026 KB  Created; hydroponic calls moved out of SPAM and
C                 CROPGRO into a single driver called from LAND.
C-----------------------------------------------------------------------
C  Called from: LAND
C  Calls:       HYDRO_FERT, HYDRO_WATER, SOLEC, SOLPH, SOLO2, OPSOL
C=======================================================================

      SUBROUTINE HYDROPONICS(CONTROL, ISWITCH,
     &    EP, ES, RWUMX, WEATHER,                         !Input
     &    TRWU, TRWUP)                                    !Output

!-----------------------------------------------------------------------
      USE ModuleDefs
      USE ModuleData
      IMPLICIT NONE
      EXTERNAL HYDRO_FERT, HYDRO_WATER, SOLEC, SOLPH, SOLO2, OPSOL
      SAVE

      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH
      TYPE (WeatherType) WEATHER

!     Input
      REAL EP          ! Actual plant transpiration (mm/d), from SPAM
      REAL RWUMX       ! Max water uptake per unit root length, from PLANT

      REAL ES          ! Solution evaporation (mm/d), from SPAM
!     Output
      REAL TRWU        ! Actual root water uptake (cm/d)
      REAL TRWUP       ! Potential root water uptake (cm/d)

!     Local
      INTEGER DYNAMIC
      LOGICAL WATBAL                           ! Solution water balance on
      REAL NO3_CONC, NH4_CONC, P_CONC, K_CONC  ! Solution conc. (mg/L)
      REAL EC_CALC, EC_TARGET                  ! EC (dS/m)
      REAL PH_CALC, PH_TARGET                  ! pH
      REAL DO2_CALC, DO2_SAT                   ! Dissolved O2 (mg/L)
      REAL UNO3_HYDRO, UNH4_HYDRO              ! N uptake (kg/ha/d)
      REAL UPO4_HYDRO, UK_HYDRO                ! P, K uptake (kg/ha/d)
      REAL ROOT_RESP                           ! Root respiration
      REAL PLTPOP                              ! Plants per m2

      IF (ISWITCH % ISWHYDRO .NE. 'Y') RETURN

      DYNAMIC = CONTROL % DYNAMIC

!     WATER=N: the solution water balance is not simulated. The daily
!     volume update (HYDRO_WATER, INTEGR) is skipped, so the volume
!     stays at its initial value. The potential root water uptake
!     (HYDRO_WATER, RATE) is a root property used by ETPHOT and is
!     still calculated. All other solution routines still run.
      WATBAL = ISWITCH % ISWWAT .NE. 'N'

!***********************************************************************
!     Run initialization
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
        CALL OPSOL(CONTROL, ISWITCH)

!***********************************************************************
!     Seasonal initialization
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
!       Solution volume
        CALL HYDRO_WATER(CONTROL, ISWITCH,
     &    0.0, ES, RWUMX,                                 !Input
     &    TRWUP, TRWU)                                    !Output

!       Initial solution concentrations
        CALL GET('HYDRO','NO3_CONC',NO3_CONC)
        CALL GET('HYDRO','NH4_CONC',NH4_CONC)
        CALL GET('HYDRO','P_CONC',P_CONC)
        CALL GET('HYDRO','K_CONC',K_CONC)
        IF (P_CONC .LT. 0.1) THEN
          P_CONC = 60.0      ! Default P concentration (mg/L)
          CALL PUT('HYDRO','P_CONC',P_CONC)
        ENDIF
        IF (K_CONC .LT. 0.1) THEN
          K_CONC = 240.0     ! Default K concentration (mg/L)
          CALL PUT('HYDRO','K_CONC',K_CONC)
        ENDIF

        CALL GET('PLANT','PLTPOP',PLTPOP)
        IF (PLTPOP .LT. 0.1) PLTPOP = 10.0

        CALL SOLEC(CONTROL, ISWITCH,
     &    NO3_CONC, NH4_CONC, P_CONC, K_CONC,             !Input
     &    EC_CALC, EC_TARGET)                             !Output
        CALL SOLPH(CONTROL, ISWITCH,
     &    0.0, 0.0,                                       !Input
     &    PH_CALC, PH_TARGET)                             !Output
        ROOT_RESP = 0.0
        CALL SOLO2(CONTROL, ISWITCH, WEATHER,
     &    PLTPOP, ROOT_RESP,                              !Input
     &    DO2_CALC, DO2_SAT)                              !Output

        UNO3_HYDRO = 0.0
        UNH4_HYDRO = 0.0
        UPO4_HYDRO = 0.0
        UK_HYDRO   = 0.0

        CALL HYDRO_FERT(CONTROL, ISWITCH)
        CALL OPSOL(CONTROL, ISWITCH)

!***********************************************************************
!     Daily rate calculations
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!       Potential root water uptake from the solution
        CALL HYDRO_WATER(CONTROL, ISWITCH,
     &    0.0, ES, RWUMX,                                 !Input
     &    TRWUP, TRWU)                                    !Output

!       EC and EC stress factors (read by SPAM for transpiration and
!       by the plant routines for uptake and growth)
        CALL GET('HYDRO','NO3_CONC',NO3_CONC)
        CALL GET('HYDRO','NH4_CONC',NH4_CONC)
        CALL GET('HYDRO','P_CONC',P_CONC)
        CALL GET('HYDRO','K_CONC',K_CONC)
        CALL GET('PLANT','PLTPOP',PLTPOP)
        CALL SOLEC(CONTROL, ISWITCH,
     &    NO3_CONC, NH4_CONC, P_CONC, K_CONC,             !Input
     &    EC_CALC, EC_TARGET)                             !Output

!       pH from the previous day's N uptake (stored by NUPTAK)
        UNO3_HYDRO = 0.0
        UNH4_HYDRO = 0.0
        UPO4_HYDRO = 0.0
        UK_HYDRO   = 0.0
        CALL GET('HYDRO','UNO3',UNO3_HYDRO)
        CALL GET('HYDRO','UNH4',UNH4_HYDRO)
        CALL GET('HYDRO','UPO4',UPO4_HYDRO)
        CALL GET('HYDRO','UK',UK_HYDRO)
        CALL SOLPH(CONTROL, ISWITCH,
     &    UNO3_HYDRO, UNH4_HYDRO,                         !Input
     &    PH_CALC, PH_TARGET)                             !Output

!       Dissolved oxygen (root respiration from ROOTS, one day lag)
        CALL GET('HYDRO','ROOT_RESP',ROOT_RESP)
        IF (ROOT_RESP .LT. 0.0) ROOT_RESP = 0.0
        CALL SOLO2(CONTROL, ISWITCH, WEATHER,
     &    PLTPOP, ROOT_RESP,                              !Input
     &    DO2_CALC, DO2_SAT)                              !Output

!***********************************************************************
!     Daily integration
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. INTEGR) THEN
!       Scheduled nutrient additions
        CALL HYDRO_FERT(CONTROL, ISWITCH)

!       Solution volume from today's transpiration (EP from SPAM RATE)
        IF (WATBAL) CALL HYDRO_WATER(CONTROL, ISWITCH,
     &    EP, ES, RWUMX,                                  !Input
     &    TRWUP, TRWU)                                    !Output

!       Solution properties after water loss and additions
        CALL SOLEC(CONTROL, ISWITCH,
     &    NO3_CONC, NH4_CONC, P_CONC, K_CONC,             !Input
     &    EC_CALC, EC_TARGET)                             !Output
        CALL SOLPH(CONTROL, ISWITCH,
     &    UNO3_HYDRO, UNH4_HYDRO,                         !Input
     &    PH_CALC, PH_TARGET)                             !Output
        CALL SOLO2(CONTROL, ISWITCH, WEATHER,
     &    PLTPOP, ROOT_RESP,                              !Input
     &    DO2_CALC, DO2_SAT)                              !Output

!***********************************************************************
!     Daily and seasonal output
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT .OR. DYNAMIC .EQ. SEASEND) THEN
        CALL OPSOL(CONTROL, ISWITCH)

      ENDIF

      RETURN
      END SUBROUTINE HYDROPONICS

!-----------------------------------------------------------------------
! DO2_CALC   Dissolved oxygen concentration (mg/L)
! EC_CALC    Electrical conductivity of the solution (dS/m)
! EP         Actual plant transpiration (mm/d)
! ES         Evaporation from the solution at the plant holes (mm/d)
! PH_CALC    Solution pH
! PLTPOP     Plant population (plants/m2)
! ROOT_RESP  Root respiration, used for O2 consumption
! RWUMX      Max water uptake per unit root length (cm3/cm root/d)
! TRWU       Actual root water uptake (cm/d)
! TRWUP      Potential root water uptake (cm/d)
! UNO3_HYDRO NO3 uptake from the solution (kg N/ha/d)
! UNH4_HYDRO NH4 uptake from the solution (kg N/ha/d)
!-----------------------------------------------------------------------
