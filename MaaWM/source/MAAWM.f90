!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                       !!!
!!!  MaaWM Simulation Program                                             !!!
!!!  ------------------------                                             !!!
!!!                                                                       !!!
!!!  This program simulates the dynamics of a fixed wing aircraft         !!!
!!!  It also simulates a quadcopter (Lisa Schibelius 2016)                !!!
!!!  The code is written in the FORTRAN 95 programming language.          !!!
!!!                                                                       !!!
!!!  The program has been developed by:                                   !!!
!!!                                                                       !!!
!!!  Carlos Montalvo                                                      !!!
!!!  University of South Alabama                                          !!!
!!!  Mobile, Alabama 36608                                                !!!
!!!  Electronic Mail: cmontalvo@southalabama.edu                          !!!
!!!                                                                       !!!
!!!  Code Adapted from TAPAS version 2.0 developed by:                    !!!
!!!                                                                       !!!
!!!  Mark Costello - Georgia Institute of Technology                      !!!
!!!                                                                       !!!
!!!  Version 1.0 Released November 2014                                   !!!
!!!                                                                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!! MODULE MAAWMDATATYPES!!!!!!!!!!!!!!!!!!!!!!!!!!

module MAAWMDATATYPES
 
 integer,parameter :: MAXNALT = 100                ! Units: 'nd', Desc: 'Maximum Number of Atmosphere Altitude Table Points'
 integer,parameter :: MAXXYDIM = 100               ! Units: 'nd', Desc: 'Maximum Number of Dimension for Atmosphere Model
 integer,parameter :: MAXNLSE = 20                 ! Units: 'nd', Desc: 'Maximum Number of Lifting Surface Elements'
 integer,parameter :: MAXPROP = 20                 ! Units: 'nd', Desc: 'Maximum Number of Propellor Table size
 integer,parameter :: MAXNAOA = 100                ! Units: 'nd', Desc: 'Maximum Number of Aerodynamic Angle of Attack Table Points'
 integer,parameter :: MAXX = 20                  ! Units: 'nd', Desc: 'Maximum Number of System States'
 integer,parameter :: NOACTUATORS = 1              ! Units: 'nd', Desc: Number of Actuators
 integer,parameter :: MAXAIRCRAFT = 10             ! Units: 'nd', Desc: Number of Aircraft
 integer,parameter :: MAXCOMMAND = 1000            ! Units: 'nd', Desc: Maximum number of waypoints allowed
 integer,parameter :: NDATAMAX = 100000              ! Units: 'nd', Desc: Number of data points to sample before surrogen is run
 real*8,parameter :: FT2M = 0.3048                ! Conversion from Feet to Meters
 real*8,parameter :: M2FT = 3.28084               ! Conversion from Meters to Feet
 real*8,parameter :: PI = 3.14159265358979323846   ! Units: 'nd', Desc: 'Pi'

!!!!!!!!!!!!!!!!!!!!!!ATMOSPHERE STRUCTURE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 type ATMOSPHERESTRUCTURE
    integer :: markX = 1                             ! X coordinate Marker in interpolation routine 
    integer :: markY = 1                             ! Y coordinate Marker in interpolation routine 
    integer :: markZ = 1                             ! Z coordinate Marker in interpolation routine 
    integer :: markT = 1                             ! T coordinate Marker in interpolation routine 
    integer :: markXT = 1                            ! X coordinate Marker in turbulence interpolation routine 
    integer :: markYT = 1                            ! Y coordinate Marker in turbulence interpolation routine 
    integer :: markZT = 1                            ! Z coordinate Marker in turbulenceinterpolation routine 
    integer :: bounds = 0                            ! Flag indicating wether or not you have gone out of WRF bounds
    integer :: boundflag = 1                         ! Flag indicating wether or not you have gone out of WRF bounds
    integer :: boundsT = 0                           ! Flag indicating wether or not you have gone out of TURB bounds
    integer :: boundflagT = 1                        ! Flag indicating wether or not you have gone out of TURB bounds
    integer :: parameters(5) = 0                     ! Vector of a few parameters 
    integer :: parametersSEVEN(7) = 0                     ! Vector of a few parameters 
    integer :: parametersFIVE(5) = 0                     ! Vector of a few parameters 
    integer :: tcoord(601) = 0                       ! Vector of time
    integer :: dimX = 0                              ! Dimension of Spatial wind data matrices
    integer :: dimY = 0                              ! Dimension of Spatial wind data matrices
    integer :: dimZ = 0                              ! Dimension of Spatial wind data matrices
    integer :: dimT = 500                            ! Dimension of Turbulence wind data matrices
    integer :: tlength = 601                         ! Dimension of time vector
    integer :: OFFON = 0                             ! Units: 'nd', Desc: 'Off/On Switch'
    integer :: TABSIZE = 0                           ! Units: 'nd', Desc: 'Table Size'
    integer :: MODNO = 1                             ! Units: 'nd', Desc: 'Model Number'
    integer :: IP = 1                                ! Units: 'nd', Desc: 'Table Pointer'
    integer :: DQFLAG = 0                            ! Units: 'nd', Desc: 'Data Quality Flag (0=Data Not Loaded Successfully, 1=Data Loaded Successfully)'
    integer :: TIMEVARYING = 0                       ! Units: 'nd', Desc: 'Static Wind field or not
    integer :: INTERPTYPE = 0
    real*8 :: zcoord(40) = 0                         ! Z coordinates of wind data in z (m)
    real*8 :: xcoord(MAXXYDIM) = 0                         ! X coordinates of wind data in x (m)
    real*8 :: ycoord(MAXXYDIM) = 0                         ! Y coordinates of wind data in y (m)
    real*8 :: xcoordT(500) = 0                       ! X coordinates of turbulence data in x (m)
    real*8 :: ycoordT(500) = 0                       ! Y coordinates of turbulence data in y (m)
    !real*8 :: terrain(40,40) = 0                     ! Terrain height in meters
    !!!REVISIT REVISIT REVIST - Always make sure this is back to normal
    ! real*8 :: U0(1,1,1),V0(1,1,1),W0(1,1,1)          ! U,V,W velocity at timestep t (m/s)
    ! real*8 :: Udt(1,1,1),Vdt(1,1,1),Wdt(1,1,1)       ! U,V,W velocity at timestep t+dt (m/s)
    ! real*8 :: UTURB(1,1),VTURB(1,1),WTURB(1,1)       ! U,V,W turbulence at timestep t (m/s)
    real*8 :: U0(MAXXYDIM,MAXXYDIM,MAXXYDIM),V0(MAXXYDIM,MAXXYDIM,MAXXYDIM),W0(MAXXYDIM,MAXXYDIM,MAXXYDIM) ! U,V,W velocity at timestep t (m/s)
    real*8 :: Udt(MAXXYDIM,MAXXYDIM,MAXXYDIM),Vdt(MAXXYDIM,MAXXYDIM,MAXXYDIM),Wdt(MAXXYDIM,MAXXYDIM,MAXXYDIM) ! U,V,W velocity at timestep t+dt (m/s)
    real*8 :: UTURB(500,500),VTURB(500,500),WTURB(500,500) ! U,V,W turbulence at timestep t (m/s)
    real*8 :: dx = 0                                 ! X resolution (m)
    real*8 :: dxT = 1                                ! XY resolution of turbulence (m)
    real*8 :: dyT = 1                                ! XY resolution of turbulence (m)
    real*8 :: dy = 0                                 ! X resolution (m)
    real*8 :: ztop = 0                               ! Highest altitude to sample from (m)
    real*8 :: IWINDSCALE = 1                         ! Scale of WRF winds
    real*8 :: TURBLEVEL = 1                          ! Scale of Turbulence 
    real*8 :: Vtrim = 65.6                             ! Trim velocity of craft flying through winds (ft/s)
    real*8 :: TIMESTEP = 0.0001                      ! Timestep of Simulation
    real*8 :: WINDGUST(3) = 0                        ! Vector holding WINDGUST values at current timestep (ft/s)
    real*8 :: VWAKE(3) = 0                           ! Vector holding wake data values at current timestep (ft/s)
    real*8 :: WRFX = 0                               ! Wrf winds in x (ft/s)
    real*8 :: WRFY = 0                               ! Wrf winds in y (ft/s)
    real*8 :: WRFZ = 0                               ! Wrf winds in z (ft/s)
    real*8 :: ALT = 0.0                              ! Units: 'm', Desc: 'Altitude'
    real*8 :: DEN = 0.002363                          ! Units: 'slugs/ft^3', Desc: 'Density'
    real*8 :: WINDSPEED = 0.0                        ! Units: 'm/s', Desc: 'Wind Speed'
    real*8 :: WINDDIR = 0.0                          ! Units: 'rad', Desc: 'Wind Direction'
    real*8 :: WINDELEV = 0.0                         ! Units: 'rad', Desc: 'Wind Elevation'
    real*8 :: XI = 0.0                               ! Units: 'm', Desc: 'Inertial X Position'
    real*8 :: YI = 0.0                               ! Units: 'm', Desc: 'Inertial Y Position'
    real*8 :: ZI = 0.0                               ! Units: 'm', Desc: 'Inertial Z Position'
    real*8 :: VXWIND = 0.0                           ! Units: 'm/s', Desc: 'Atmospheric Wind Along Inertial I Axis'
    real*8 :: VYWIND = 0.0                           ! Units: 'm/s', Desc: 'Atmospheric Wind Along Inertial J Axis'
    real*8 :: VZWIND = 0.0                           ! Units: 'm/s', Desc: 'Atmospheric Wind Along Inertial K Axis'
    real*8 :: ALTTAB(MAXNALT) = 0.0                  ! Units: 'm', Desc: 'Density Altitude Table'
    real*8 :: DENTAB(MAXNALT) = 0.0                  ! Units: 'kg/m^3', Desc: 'Density Table'
    real*8 :: VXWINDTAB(MAXNALT) = 0.0               ! Units: 'm/s', Desc: 'VXWIND Wind Table'
    real*8 :: VYWINDTAB(MAXNALT) = 0.0               ! Units: 'm/s', Desc: 'VYWIND Wind Table'
    real*8 :: VZWINDTAB(MAXNALT) = 0.0               ! Units: 'm/s', Desc: 'VZWIND Wind Table'
    real*8 :: xshift = 0                             ! Units: 'm' , 'Desc: shift in x coordinate for interpolation
    real*8 :: yshift = 0                             ! Units: 'm' , 'Desc: shift in y coordinate for interpolation
    real*8 :: zshift = 0                             ! Units: 'm' , 'Desc: shift in y coordinate for interpolation
    real*8 :: xshiftT = 0                            ! Units: 'm' , 'Desc: shift in x coordinate for interpolation
    real*8 :: yshiftT = 0                            ! Units: 'm' , 'Desc: shift in y coordinate for interpolation
    real*8 :: zshiftT = 0                            ! Units: 'm' , 'Desc: shift in y coordinate for interpolation
    real*8 :: PSIOFFSET = 0                          ! Units: 'rad' ,'Desc: used to rotate the dryden and WRF model'
    real*8 :: WAVESPEED(2) = 0                       ! Units: 'm/s', Desc: Speed that the static waves propagate to simulate time varying component
    real*8 :: RAMPTIME = 0                           ! Units: 'sec', Desc: Time to ramp in winds
    character*128 PATH
    character*256 U0name,V0name,W0name,UTurbname,Vturbname,Wturbname
    character*256 Udtname,Vdtname,Wdtname
 end type ATMOSPHERESTRUCTURE

!!!!!!!!!!!!!!!!!!!!!!!!!!!! QUADCOPTER STRUCTURE !!!!!!!!!!!!!!!!!!!!!!!!!!

include 'coptermodule'

!!!!!!!!!!!!!!!!!!!!!!!!!!! AIRCRAFT STRUCTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 type AIRCRAFTSTRUCTURE
  integer :: DYNOFFON = 0                          ! Units: 'nd', Desc: 'Dynamics Flag (0=Off, 1=On)'
  integer :: GRAVOFFON = 0                         ! Units: 'nd', Desc: 'Gravity Flag (0=Off, 1=On)'
  integer :: AEROOFFON = 0                         ! Units: 'nd', Desc: 'Aerodynamics Flag (0=Off, 1=On)'
  integer :: CONTOFFON = 0                         ! Units: 'nd', Desc: 'Contact Flag (0=Off, 1=On)'
  integer :: DQFLAG = 0                            ! Units: 'nd', Desc: 'Data Quality Flag (0=Data Not Loaded Successfully, 1=Data Loaded Successfully)'
  integer :: NLSE = 0                              ! Units: 'nd', Desc: 'Number of Lifting Surface Elements'
  integer :: IP(MAXNLSE) = 1                       ! Units: 'nd', Desc: 'Table Look-up Pointer'
  integer :: TABSIZE(MAXNLSE) = 1                  ! Units: 'nd', Desc: 'Table Size'
  integer :: IPROP(2) = 1                          ! Prop table look up pointer
  integer :: NPROPTAB = 0                          ! Propellor table size
  integer :: THRUSTMOD = 0                         ! 0 = linear , 1 = table look up
  integer :: COMMANDFLAG(MAXCOMMAND) = 1           ! 0 or 1 to advance waypoint command
  integer :: WAYPOINT = 1                          ! Waypoint counter
  integer :: BFLAG = 0                             ! Outside the boundary
  integer :: NUMWAYPOINTS = 0                      ! Number of Waypoints
  real*8 :: TURNRADIUS = 0                         ! Turn Radius of aircraft
  real*8 :: MINWAYPOINT = 0                        ! Minimum waypoint radius 
  real*8 :: C_T_TAB(MAXPROP) = 0                   ! Thrust table
  real*8 :: OMEGA_TAB(MAXPROP) = 0                 ! Rotational speed table
  real*8 :: ADVANCERATIO_TAB(MAXPROP) = 0          ! Advance ratio table
  real*8 :: DELTHRUST_TAB(MAXPROP) = 0             ! control input table
  real*8 :: SAREA = 0                              ! Reference area of aircraft (m^2)
  real*8 :: B = 0                                  ! Wingspan of aircraft(m)
  real*8 :: C_BAR = 0                              ! Mean chord of aircraft(m)
  real*8 :: V_T = 0                                ! Trim velocity of aircraft (m/s)
  real*8 :: PD = 0                                 ! propellor diameter (m)
  real*8 :: C_L_0 = 0                              ! zero lift slope
  real*8 :: C_L_ALPHA = 0                          ! lift curve slope
  real*8 :: C_L_U = 0                              ! increase in lift with speed
  real*8 :: C_L_Q = 0                              ! increase in lift with q
  real*8 :: C_L_DE = 0                             ! increase in lift with elevator
  real*8 :: C_Y_BETA = 0                           ! side force w.r.t sideslip
  real*8 :: C_Y_P = 0                              ! side force w.r.t roll rate
  real*8 :: C_Y_R = 0                              ! side force w.r.t. yaw rate
  real*8 :: C_Y_DR = 0                             ! side force w.r.t. rudder
  real*8 :: C_Y_DA = 0                             ! side force w.r.t. aileron
  real*8 :: C_D_0 = 0                              ! zero lift drag
  real*8 :: C_D_ALPHA2 = 0                         ! drag polar
  real*8 :: C_D_U = 0                              ! increase in drag with speed
  real*8 :: C_D_DE = 0                             ! inrease in drag with elevator
  real*8 :: C_L_BETA = 0                           ! roll moment w.r.t beta
  real*8 :: C_L_P = 0                              ! roll moment w.r.t roll rate
  real*8 :: C_L_R = 0                              ! roll moment w.r.t. yaw rate
  real*8 :: C_L_DR = 0                             ! roll moment w.r.t. rudder
  real*8 :: C_L_DA = 0                             ! roll moment w.r.t. aileron
  real*8 :: C_M_0 = 0                              ! zero lift moment
  real*8 :: C_M_ALPHA = 0                          ! pitch moment curve
  real*8 :: C_M_U = 0                              ! increase in pitch moment w.r.t speed
  real*8 :: C_M_Q = 0                              ! pitch moment w.r.t. q
  real*8 :: C_M_DE = 0                             ! pitch moment w.r.t elevator
  real*8 :: C_N_BETA = 0                           ! yaw moment w.r.t. beta
  real*8 :: C_N_P = 0                              ! yaw moment w.r.t. roll rate
  real*8 :: C_N_R = 0                              ! yaw moment w.r.t. yaw rate
  real*8 :: C_N_DR = 0                             ! yaw moment w.r.t. rudder
  real*8 :: C_N_DA = 0                             ! yaw moment w.r.t. aileron
  real*8 :: C_X_DT = 0                             ! Thrust w.r.t. delthrust
  real*8 :: ALC = 0                                ! Quad aero parameter
  real*8 :: ALS = 0                                ! Quad aero parameter
  real*8 :: DXD = 0                                ! Quad aero parameter
  real*8 :: DYD = 0                                ! Quad aero parameter
  real*8 :: RNEW = 0                               ! Quad aero parameter
  real*8 :: C_T = 0                                ! Quad aero parameter
  real*8 :: C_TAU = 0                              ! Quad aero parameter
  real*8 :: LPHI = 0                               ! Quad aero parameter
  real*8 :: LTHETA = 0                             ! Quad aero parameter
  real*8 :: OMEGAMAX = 0                           ! Quad aero parameter
  real*8 :: OMEGAVEC(4,1) = 0                      ! Quad aero parameter
  real*8 :: KT = 0                                 ! Quad Aero Parameter
  real*8 :: OMEGA0 = 0                             ! Quad Aero Parameter
  real*8 :: IRR = 0                                ! Quad Aero Parameter
  real*8 :: MASS = 0.0                             ! Units: 'kg', Desc: 'Mass'
  real*8 :: WEIGHT = 0.0                           ! Units: 'N', Desc: 'Weight'  
  real*8 :: SLCG = 0.0                             ! Units: 'm', Desc: 'Stationline of Mass Center'
  real*8 :: BLCG = 0.0                             ! Units: 'm', Desc: 'Buttline of Mass Center'
  real*8 :: WLCG = 0.0                             ! Units: 'm', Desc: 'Waterline of Mass Center'
  real*8 :: IXX = 0.0                              ! Units: 'kg m^2', Desc: 'Ixx of Inertia Matrix'
  real*8 :: IYY = 0.0                              ! Units: 'kg m^2', Desc: 'Iyy of Inertia Matrix'
  real*8 :: IZZ = 0.0                              ! Units: 'kg m^2', Desc: 'Izz of Inertia Matrix'
  real*8 :: IXY = 0.0                              ! Units: 'kg m^2', Desc: 'Ixy of Inertia Matrix'
  real*8 :: IXZ = 0.0                              ! Units: 'kg m^2', Desc: 'Ixz of Inertia Matrix'
  real*8 :: IYZ = 0.0                              ! Units: 'kg m^2', Desc: 'Iyz of Inertia Matrix'
  real*8 :: IXXI = 0.0                             ! Units: '1/(kg m^2)', Desc: 'Ixx Inverse of Inertia Matrix'
  real*8 :: IYYI = 0.0                             ! Units: '1/(kg m^2)', Desc: 'Iyy Inverse of Inertia Matrix'
  real*8 :: IZZI = 0.0                             ! Units: '1/(kg m^2)', Desc: 'Izz Inverse of Inertia Matrix'
  real*8 :: IXYI = 0.0                             ! Units: '1/(kg m^2)', Desc: 'Ixy Inverse of Inertia Matrix'
  real*8 :: IXZI = 0.0                             ! Units: '1/(kg m^2)', Desc: 'Ixz Inverse of Inertia Matrix'
  real*8 :: IYZI = 0.0                             ! Units: '1/(kg m^2)', Desc: 'Iyz Inverse of Inertia Matrix'
  real*8 :: TIA(3,3) = 0.0                         ! Units: 'nd', Desc: 'aircraft to Inertial Frame Transformation Matrix'
  real*8 :: TAI(3,3) = 0.0                         ! Units: 'nd', Desc: 'Inertial to aircraft Frame Transformation Matrix'
  real*8 :: ELEVATOR = 0                           ! Units: 'rad', Desc: Elevator Command
  real*8 :: RUDDER = 0                             ! Units: 'rad', Desc: RUDDER Command
  real*8 :: AILERON = 0                            ! Units: 'rad', Desc: AILERON Commandn
  real*8 :: DELTHRUST = 0                          ! Units: 'rad', Desc: DELTHRUST Command (Total Thrust = Cxdt*DELTHRUST + T0)
  real*8 :: SLLSE(MAXNLSE) = 0.0                   ! Units: 'm', Desc: 'Stationline of Lifting Surface Element'
  real*8 :: BLLSE(MAXNLSE) = 0.0                   ! Units: 'm', Desc: 'Waterline of Lifting Surface Element'
  real*8 :: WLLSE(MAXNLSE) = 0.0                   ! Units: 'm', Desc: 'Buttline of Lifting Surface Element'
  real*8 :: AREALSE(MAXNLSE) = 0.0                 ! Units: 'm', Desc: 'Lifting Surface Element Reference Area'
  real*8 :: ALPHALSE(MAXNLSE) = 0.0                ! Units: 'rad', Desc: 'Aerodynamic Angle of Attack of Lifting Surface Element'
  real*8 :: PHILSE(MAXNLSE) = 0.0                  ! Units: 'rad', Desc: 'Phi Orientation Angle of Lifting Surface Element'
  real*8 :: GAMLSE(MAXNLSE) = 0.0                  ! Units: 'rad', Desc: 'Gamma Orientation Angle of Lifting Surface Element'
  real*8 :: DELTALSE(MAXNLSE) = 0.0                ! Units: 'rad', Desc: 'Delta Orientation Angle of Lifting Surface Element'
  real*8 :: FXLSE(MAXNLSE) = 0.0                   ! Units: 'nd', Desc: 'Fx Force of Lifting Surface Element in aircraft Reference Frame'
  real*8 :: FYLSE(MAXNLSE) = 0.0                   ! Units: 'nd', Desc: 'Fy Force of Lifting Surface Element in aircraft Reference Frame'
  real*8 :: FZLSE(MAXNLSE) = 0.0                   ! Units: 'nd', Desc: 'Fz Force of Lifting Surface Element in aircraft Reference Frame'
  real*8 :: MXLSE(MAXNLSE) = 0.0                   ! Units: 'nd', Desc: 'Mx Moment of Lifting Surface Element in aircraft Reference Frame'
  real*8 :: MYLSE(MAXNLSE) = 0.0                   ! Units: 'nd', Desc: 'My Moment of Lifting Surface Element in aircraft Reference Frame'
  real*8 :: MZLSE(MAXNLSE) = 0.0                   ! Units: 'nd', Desc: 'Mz Moment of Lifting Surface Element in aircraft Reference Frame'
  real*8 :: ALPHALSETAB(MAXNLSE,MAXNAOA) = 0.0     ! Units: 'rad', Desc: 'Aerodynamic Angle of Attack Table for Lifting Surface Element'
  real*8 :: CLLSETAB(MAXNLSE,MAXNAOA) = 0.0        ! Units: 'rad', Desc: 'Lift Coefficient Table for Lifting Surface Element'
  real*8 :: CDLSETAB(MAXNLSE,MAXNAOA) = 0.0        ! Units: 'rad', Desc: 'Drag Coefficient Table for Lifting Surface Element'
  real*8 :: FXGRAV = 0.0                           ! Units: 'N', Desc: 'X Gravity Forces in Body Frame'
  real*8 :: FYGRAV = 0.0                           ! Units: 'N', Desc: 'Y Gravity Forces in Body Frame'
  real*8 :: FZGRAV = 0.0                           ! Units: 'N', Desc: 'Z Gravity Forces in Body Frame'
  real*8 :: MXGRAV = 0.0                           ! Units: 'N m', Desc: 'X Gravity Moment About CG in Body Frame'
  real*8 :: MYGRAV = 0.0                           ! Units: 'N m', Desc: 'Y Gravity Moment About CG in Body Frame'
  real*8 :: MZGRAV = 0.0                           ! Units: 'N m', Desc: 'Z Gravity Moment About CG in Body Frame'
  real*8 :: FXAERO = 0.0                           ! Units: 'N', Desc: 'X Aerodynamic Force in Body Frame'
  real*8 :: FYAERO = 0.0                           ! Units: 'N', Desc: 'Y Aerodynamic Force in Body Frame'
  real*8 :: FZAERO = 0.0                           ! Units: 'N', Desc: 'Z Aerodynamic Force in Body Frame'
  real*8 :: MXAERO = 0.0                           ! Units: 'N m', Desc: 'X Aerodynamic Moment About CG in Body Frame'
  real*8 :: MYAERO = 0.0                           ! Units: 'N m', Desc: 'Y Aerodynamic Moment About CG in Body Frame'
  real*8 :: MZAERO = 0.0                           ! Units: 'N m', Desc: 'Z Aerodynamic Moment About CG in Body Frame'
  real*8 :: FXCONT = 0.0                           ! Units: 'N', Desc: 'X Contact Force in Body Frame'
  real*8 :: FYCONT = 0.0                           ! Units: 'N', Desc: 'Y Contact Force in Body Frame'
  real*8 :: FZCONT = 0.0                           ! Units: 'N', Desc: 'Z Contact Force in Body Frame'
  real*8 :: MXCONT = 0.0                           ! Units: 'N m', Desc: 'X Contact Moment in Body Frame'
  real*8 :: MYCONT = 0.0                           ! Units: 'N m', Desc: 'Y Contact Moment in Body Frame'
  real*8 :: MZCONT = 0.0                           ! Units: 'N m', Desc: 'Z Contact Moment in Body Frame'
  real*8 :: FXTOTAL = 0.0                          ! Units: 'N', Desc: 'X Total Applied Force in Body Frame'
  real*8 :: FYTOTAL = 0.0                          ! Units: 'N', Desc: 'Y Total Applied Force in Body Frame'
  real*8 :: FZTOTAL = 0.0                          ! Units: 'N', Desc: 'Z Total Applied Force in Body Frame'
  real*8 :: MXTOTAL = 0.0                          ! Units: 'N m', Desc: 'X Total Applied Moment About Mass Center in Body Frame'
  real*8 :: MYTOTAL = 0.0                          ! Units: 'N m', Desc: 'Y Total Applied Moment About Mass Center in Body Frame'
  real*8 :: MZTOTAL = 0.0                          ! Units: 'N m', Desc: 'Z Total Applied Moment About Mass Center in Body Frame'
  real*8 :: INITIALPHI = 0.0                       ! Units: 'rad', Desc: 'Initial Euler Roll Angle of aircraft'
  real*8 :: INITIALTHETA = 0.0                     ! Units: 'rad', Desc: 'Initial Euler Pitch Angle of aircraft'
  real*8 :: INITIALPSI = 0.0                       ! Units: 'rad', Desc: 'Initial Euler Yaw Angle of aircraft'
  real*8 :: INITIALSTATE(12) = 0.0                 ! Units: 'vd', Desc: 'Initial aircraft State Vector'
  real*8 :: PHI = 0                                ! Units: 'rad', Desc: aircraft roll angle
  real*8 :: THETA = 0                              ! Units: 'rad', Desc: aircraft roll angle
  real*8 :: PSI = 0                                ! Units: 'rad', Desc: aircraft roll angle
  real*8 :: STATE(12) = 0.0                        ! Units: 'vd', Desc: 'aircraft State Vector'
  real*8 :: STATEPREV(12) = 0.0                    ! Units: 'vd', Desc: 'aircraft State Vector before camera snapshot'
  real*8 :: STATEDOT(12) = 0.0                     ! Units: 'vd', Desc: 'aircraft State Vector Derivative'
  real*8 :: TARGETPOSITIONFROMCAMERA(3) = -999.0   ! Units: 'nd', sensed target position in aircraft frame
  real*8 :: TARGETROTATIONFROMCAMERA(3) = -999.0   ! Units: 'nd', sensed target rotation in aircraft frame
  real*8 :: TARGETBEARINGFROMCAMERA(2) = -999.0    ! Units: 'nd', sensed bearing angle from long range sensor
  real*8 :: PSICOMMAND = 0                         ! Units: 'nd', slope of heading path to track
  real*8 :: UCOMMAND = 20.0                        ! Units: 'm/s', velocity command
  real*8 :: XCOMMAND(MAXCOMMAND) = 500                         ! Units: 'm', X Waypoint command
  real*8 :: YCOMMAND(MAXCOMMAND) = 500                         ! Units: 'm', Y Waypoint command
  real*8 :: ZCOMMAND(MAXCOMMAND) = -200                        ! Units: 'm', Altitude command
  real*8 :: NOMINALSTATE(MAXX) = 0.0               ! Units: 'md', Nominal state vector of aircraft
  real*8 :: KRKBODY(MAXX,4) = 0.0                  ! Units: 'md', RK4 vector
 end type AIRCRAFTSTRUCTURE

!!!!!!!!!!!!!!!!!!!!!!!!!! SIMULATION STRUCTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 type SIMULATIONSTRUCTURE
  character(128) :: SIMINPUTFILE = ' '             ! Units: 'nd', Desc: 'Simulation Input File'
  character(128) :: RESTARTFILE = ' '         ! Units: 'nd', Desc: 'Simulation Input File'
  character(128) :: STATEDES(MAXX) = ' '           ! Units: 'nd', Desc: 'State Description'
  integer :: NOSTATES = 1                          ! Units: 'nd', Desc: 'Number of States'
  integer :: NOQUADSTATES = 1                          ! Units: 'nd', Desc: 'Number of States'
  integer :: IOUTSKIP = 1                          ! Units: 'nd', Desc: 'Output Skip Parameter'
  integer :: ICS = 0                               ! Units: 'nd', Desc: 'Control System Off On Flag'
  integer :: IDEBUG = 0                            ! Units: 'nd', Desc: 'Debug Flag'
  integer :: IDX = 0                               ! Units: 'nd', Desc: 'Integration Index'
  integer :: IDXOUT = 0                            ! Units: 'nd', Desc: ''
  integer :: IECHOOFFON = 1                        ! Units: 'nd', Desc: 'Input Data Echo Flag'
  integer :: DQFLAG = 0                            ! Units: 'nd', Desc: 'Data Quality Flag (0=Data Not Loaded Successfully, 1=Data Loaded Successfully)'
  integer :: CREATERESTART = 0                     ! Units: 'nd', Desc: 'Create Restart Point'
  integer :: RESTART = 0                           ! Units: 'nd', Desc: 'Use a restart file'
  integer :: IMPACT = 0                            ! Units: 'nd', Desc: 'model increase of mass or impact'
  integer :: ACTUATORONOFF = 0                     ! Units: 'nd', Desc: Turn on actuator dynamics or not
  integer :: NOAIRCRAFT = 1.0                      ! Units: 'nd', Desc: 'Number of aircraft simulated'
  integer :: NOAIRPLANES = 1.0                     ! Units: 'nd', Desc: 'Number of airplanes simulated'
  integer :: NOQUADCOPTERS = 1.0                   ! Units: 'nd', Desc: 'Number of quadcopters simulated'
  integer :: CONTROLTYPE = 0                       ! Units: 'nd', Desc: 'Wander Aimlessly or follow waypoints'
  integer :: NDATA = 0                             ! Units: 'nd', Desc: Number of wind data points
  integer :: RBFDATA = 4000                        ! Units: 'nd', Desc: Number of wind data points to sample before RBF routine
  integer :: IKIND = 1                             ! Units: 'nd', Desc: Integer between 1-7 that changes the order of the RBF
  integer :: K35 = 1                               ! Units: 'nd', Desc: Integer between 1-4 that changes the polynomial fits from linear to full quartic
  integer :: DYNSKIP = 0                           ! Units: 'nd', Desc: 'Flag set to skip dynamic model or not
  integer :: NCENS = 0                             ! Units: 'nd', Desc: 'Number of centers computed'
  real*8 :: BOUNDARY = 0                           ! Units: 'm', Desc: Boundary of world
  real*8 :: INCREASETIME = 0                       ! Units: 'nd', Desc: 'time to increase mass or model impulse'
  real*8 :: INCREASEMASS = 0                       ! Units: 'nd', Desc: 'amount of mass to increase'
  real*8 :: IMPACTVELOCITY = 0                     ! Units: 'nd', Desc: 'impact velocity'
  real*8 :: TIME = 0.0                             ! Units: 's', Desc: 'Time'
  real*8 :: RESTARTTIME = 0                        ! Units: 's', Desc: 'Time to create restart time
  real*8 :: DELTATIME = 0.0                        ! Units: 's', Desc: 'Delta Time'
  real*8 :: INITIALTIME = 0.0                      ! Units: 's', Desc: 'Initial Time'
  real*8 :: CPUTIMEUSER = 0.0                      ! Units: 's', Desc: 'Elapsed CPU Time by User Process'
  real*8 :: CPUTIMESYSTEM = 0.0                    ! Units: 's', Desc: 'Elapsed CPU Time by System Process'
  real*8 :: CPUTIMETOTAL = 0.0                     ! Units: 's', Desc: 'Elapsed CPU Time by User and System Processes'
  real*8 :: FINALTIME = 0.0                        ! Units: 's', Desc: 'Final Time'
  !real*8 :: STATE(MAXX) = 0.0                      ! Units: 'vd', Desc: 'State'
  real*8 :: STATEPREV(MAXX) = 0.0                  ! Units: 'vd', Desc: 'State before camera snapshot'
  real*8 :: INITIALSTATE(MAXX) = 0.0               ! Units: 'vd', Desc: 'Initial State'
  !real*8 :: STATEDOT(MAXX) = 0.0                   ! Units: 'vd', Desc: 'State Dot'
  real*8 :: INITIALSTATEDOT(MAXX) = 0.0            ! Units: 'vd', Desc: 'Initial State Dot'
  real*8 :: ACTUATOR(NOACTUATORS) = 0.0            ! Units: 'vd', Desc: 'Vector of Actuators
  real*8 :: ACTUATORDOT(NOACTUATORS) = 0.0         ! Units: 'vd', Desc: 'Vector of Actuators Derivatives
  
 end type SIMULATIONSTRUCTURE

!!!!!!!!!!!!!!!!!!!!!!!!!! CONTROL SYSTEM STRUCTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 type CONTROLSYSTEMSTRUCTURE
  integer :: DQFLAG = 0                            ! Units: 'nd', Desc: 'Data Quality Flag (0=Data Not Loaded Successfully, 1=Data Loaded Successfully)'
  integer :: TETHERCONTROLOFFON = 0                ! Units: 'nd', Desc: 'Tether Line Length Control Flag (0=Off, 1=On)'
  integer :: PARAFOILCONTROLOFFON = 0              ! Units: 'nd', Desc: 'Parafoil Brake Control Flag (0=Off, 1=On)'
  integer :: AIRCRAFTCONTROLOFFON = 0              ! Units: 'nd', Desc: 'aircraft controller off on (0=Off, 1=On)'
  real*8 :: SENSORERRORS = 0                      ! Units: 'nd', Desc: 'Sensor Errors off on (0=Off, 1=On)'
  real*8 :: KU = 0                                 ! Units: 'nd', Desc: airplane proportional gain on speed
  real*8 :: KPP = 0                                ! Units: 'nd', Desc: airplane proportional gain on roll angle
  real*8 :: KDP = 0                                ! Units: 'nd', Desc: airplane derivative gain on roll angle
  real*8 :: KPT = 0                                ! Units: 'nd', Desc: airplane proportional gain on pitch angle
  real*8 :: KPZ = 0                                ! Units: 'nd', Desc: airplane proportional gain on altitude
  real*8 :: KDZ = 0                                ! Units: 'nd', Desc: airplane derivative gain on altitude
  real*8 :: KPY = 0                                ! Units: 'nd', Desc: airplane proportional gain on crossrange
  real*8 :: KDY = 0                                ! Units: 'nd', Desc: airplane derivative gain on crossrange
  real*8 :: KV = 0                                 ! Units: 'nd', Desc: airplane proportional gain on sideslip
  real*8 :: KPSI = 0                               ! Units: 'nd', Desc: airplane proportional gain on heading
  real*8 :: KPPHI = 0                              ! Units: 'nd', Desc: quadcopter proportional gain on roll 
  real*8 :: KDPHI = 0                              ! Units: 'nd', Desc: quadcopter derivative gain on roll
  real*8 :: KPTHETA = 0                            ! Units: 'nd', Desc: quadcopter proportional gain on pitch
  real*8 :: KDTHETA = 0                            ! Units: 'nd', Desc: quadcopter derivative gain on pitch
  real*8 :: KPPSI = 0                              ! Units: 'nd', Desc: quadcopter proportional gain on heading
  real*8 :: KDPSI = 0                              ! Units: 'nd', Desc: quadcopter derivative gain on heading
  real*8 :: KPXQUAD = 0                            ! Units: 'nd', Desc: quadcopter proportional gain on x-position
  real*8 :: KDXQUAD = 0                            ! Units: 'nd', Desc: quadcopter derivative gain on x-position
  real*8 :: KPYQUAD = 0                            ! Units: 'nd', Desc: quadcopter proportional gain on y-position
  real*8 :: KDYQUAD = 0                            ! Units: 'nd', Desc: quadcopter derivative gain on y-position
  real*8 :: KPZQUAD = 0                            ! Units: 'nd', Desc: quadcopter proportional gain on z-position
  real*8 :: KDZQUAD = 0                            ! Units: 'nd', Desc: quadcopter derivative gain on z-position
  real*8 :: KTETHER = 0                            ! Units: 'nd', Desc: Tether proportional gain on length
  real*8 :: KDTETHER = 0                           ! Units: 'nd', Desc: Tether derivative gain on length
  real*8 :: KITETHER = 0                            ! Units: 'nd', Desc: Tether integral gain on length
  real*8 :: KPARAFOIL = 0                          ! Units: 'nd', Desc: Parafoil proportional gain on 
  real*8 :: KDPARAFOIL = 0                         ! Units: 'nd', Desc: Parafoil derivative gain on 
  real*8 :: UPDATERATE = 0.0                       ! Units: 's', Desc: 'Rate for Control System Updates'
  real*8 :: UPDATETIME = 0.0                       ! Units: 's', Desc: 'Time for Next Control System Update'
  real*8 :: LEN = 0.0                              ! Units: 'm', Desc: Feedback of length of tether
  real*8 :: THETA = 0.0                            ! Units: 'rad', Desc: Feedback of pitch angle of tether
  real*8 :: PSI = 0.0                              ! Units: 'rad', Desc: Feedback of yaw angle of tether
  real*8 :: PSI_PURE = 0.0                         ! Units: 'rad', Desc: Unpolluted Feedback of yaw angle of tether
  real*8 :: DELPREV = 0.0                          ! Units: 'nd', Desc: Previous Dely value used for derivative filter
  real*8 :: DELDOTPREV = 0.0                       ! Units: 'nd', Desc: Previous Delydot value used for derivative filter
  real*8 :: BIRDROLL = 0.0                         ! Units: 'rad', Desc: Bird Roll Angle from feedback sensors
  real*8 :: BIRDYAW = 0.0                          ! Units: 'rad', Desc: Bird Yaw Angle from feedback sensors
  real*8 :: RISERROLL = 0.0                        ! Units: 'rad', Desc: Angle of Risers from Bird to Canopy
  real*8 :: TENSIONBIRD = 0.0                      ! Units: 'N', Desc: Tension At Bird from Feedback Sensor
  real*8 :: TENSIONWINCH = 0.0                     ! Units: 'N', Desc: Tension At Winch from Feedback Sensor
  real*8 :: BIRDROLLERROR(3) = 0.0                 ! Units: 'md', Desc: Bias, Scale Factor, and Noise
  real*8 :: RISERROLLERROR(3) = 0.0                ! Units: 'md', Desc: Bias, Scale Factor, and Noise
  real*8 :: TENSIONERROR(3) = 0.0                  ! Units: 'md', Desc: Bias, Scale Factor, and Noise
 end type CONTROLSYSTEMSTRUCTURE

!!!!!!!!!!!!!!!!!!!!!!!!!! MAAWM STRUCTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 type MAAWMSTRUCTURE
  character(128) :: FILEINPUTFILE = ' '            ! Units: 'nd', Desc: 'File of Files Input File'
  character(128) :: ATMOSPHEREINPUTFILE = ' '      ! Units: 'nd', Desc: 'Atmosphere Input File'
  character(128) :: AIRPLANEINPUTFILE = ' '        ! Units: 'nd', Desc: 'AIRPLANE Input File'
  character(128) :: CSINPUTFILE = ' '              ! Units: 'nd', Desc: 'Control System Input File'
  character(128) :: SIMINPUTFILE = ' '             ! Units: 'nd', Desc: 'Simulation Input File'
  character(128) :: ICSINPUTFILE = ' '             ! Units: 'nd', Desc: 'Initial Condition Input File'
  character(128) :: WAYINPUTFILE = ' '              ! Units: 'nd', Desc:  'Waypoint Input File'
  character(128) :: RUNLOGFILE = ' '               ! Units: 'nd', Desc: 'Run Log File'
  character(128) :: STATEOUTPUTFILE = ' '          ! Units: 'nd', Desc: 'State Output File'
  character(128) :: MISCOUTPUTFILE = ' '           ! Units: 'nd', Desc: 'Miscellaneous Output File'
  character(128) :: CONTROLOUTPUTFILE = ' '        ! Units: 'nd', Desc: 'Control Output File'
  character(128) :: FORCEOUTPUTFILE = ' '          ! Units: 'nd', Desc: 'Force Output File'
  character(128) :: ERROROUTPUTFILE = ' '          ! Units: 'nd', Desc: 'Error Output File'
  character(128) :: WINDOUTPUTFILE = ' '           ! Units: 'nd', Desc: 'Wind Output File
  character(24) :: DATEANDTIMEOFRUN = ' '          ! Units: 'nd', Desc: 'Date and Time of Run'
  real*8 :: GRAVITY = 32.2                         ! Units: 'ft/s^2', Desc: 'Gravity'
  type(ATMOSPHERESTRUCTURE) :: ATM
  type(AIRCRAFTSTRUCTURE) :: AC(MAXAIRCRAFT)
  type(COPTERSTRUCTURE) :: COPTER(MAXAIRCRAFT)
  type(CONTROLSYSTEMSTRUCTURE) :: CS
  type(SIMULATIONSTRUCTURE) :: SIM
 end type MAAWMSTRUCTURE

end module MAAWMDATATYPES

!!!!!!!!!!!!!!!!!!!!!!!! PROGRAM MAAWM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM MAAWM
 use MAAWMDATATYPES
 implicit none
 integer openflag,readflag,LENGTH,k,ierr,i
 real*8 x(NDATAMAX,4),u(NDATAMAX,3),ac
 character(128) inputfilename
 character(12) inputfiletype
 type(MAAWMSTRUCTURE) T

 readflag = 0
 openflag = 0
 
 ! Input File of Files

 call getarg(1,T%FILEINPUTFILE)
 open(unit=93,file=T%FILEINPUTFILE,status='old',iostat=openflag)
 if (len(trim(T%FILEINPUTFILE)) .lt. 1) then
    write(*,*) 'No input ifiles file given as input argument'; STOP
 else if (openflag .ne. 0) then
  write(*,*) 'Error Opening MAAWM Input File: ',T%FILEINPUTFILE;  STOP
 end if
 rewind(93)
 
 do while (readflag .eq. 0)
  inputfiletype = ' '; inputfilename = ' '
  read(unit=93,fmt=*,iostat=readflag) inputfiletype,inputfilename
  write(*,*) 'Input file: ', inputfilename
  if ((inputfiletype.eq.'ATM') .or. (inputfiletype.eq.'atm') .or. (inputfiletype.eq.'Atm')) then
   T%ATMOSPHEREINPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'AC') .or. (inputfiletype.eq.'ac') .or. (inputfiletype.eq.'Ac')) then
   T%AIRPLANEINPUTFILE = inputfilename;  
  end if
  if ((inputfiletype .eq. 'QUAD') .or. (inputfiletype .eq. 'quad') .or. (inputfiletype .eq. 'Quad')) then
    T%COPTER(1)%INPUTFILE = inputfilename;
  end if
  if ((inputfiletype.eq.'CS') .or. (inputfiletype.eq.'cs') .or. (inputfiletype.eq.'Cs')) then
   T%CSINPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'SIM') .or. (inputfiletype.eq.'sim') .or. (inputfiletype.eq.'Sim')) then
   T%SIMINPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'ICS') .or. (inputfiletype.eq.'ics') .or. (inputfiletype.eq.'Ics')) then
   T%ICSINPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'WAY') .or. (inputfiletype.eq.'way') .or. (inputfiletype.eq.'Way')) then
   T%WAYINPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'LOG') .or. (inputfiletype.eq.'log') .or. (inputfiletype.eq.'log')) then
   T%RUNLOGFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'SOUT') .or. (inputfiletype.eq.'sout') .or. (inputfiletype.eq.'Sout')) then
   T%STATEOUTPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'MOUT') .or. (inputfiletype.eq.'mout') .or. (inputfiletype.eq.'Mout')) then
   T%MISCOUTPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'COUT') .or. (inputfiletype.eq.'cout') .or. (inputfiletype.eq.'Cout')) then
     T%CONTROLOUTPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'FOUT') .or. (inputfiletype.eq.'fout') .or. (inputfiletype.eq.'Fout')) then
     T%FORCEOUTPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'EOUT') .or. (inputfiletype.eq.'eout') .or. (inputfiletype.eq.'Eout')) then
     T%ERROROUTPUTFILE = inputfilename;  
  end if
  if ((inputfiletype.eq.'WOUT') .or. (inputfiletype.eq.'wout') .or. (inputfiletype.eq.'Wout')) then
     T%WINDOUTPUTFILE = inputfilename;  
  end if
 end do

 close(93)

!!!!!!!!!!!!!!!!!!! Load Data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call ATMOSPHERE(T,1)
 call AIRCRAFT(T,1,1)
 call COPTER(T%COPTER(1),1) !load copter input file
 call CONTROL(T,1)
 call SIMULATION(T,1)

!!!!!!!!!!!!!!!!!!!!!!!!! Echo Data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 openflag = 0
 open(unit=25,file=T%RUNLOGFILE,iostat=openflag)
 if (openflag .ne. 0) then
  write(*,*) 'Error Opening Run Log File: ',T%RUNLOGFILE;  STOP
 end if

 write(25,*) ' '
 write(25,*) '**********************************************************************************************************'
 write(25,*) '**********************************************************************************************************'
 write(25,*) '                                                                                                          '
 write(25,*) 'Multi-Agent Atmospheric Wind Mapper                                                                                                       '
 write(25,*) '                                                                                                          '
 write(25,*) 'Version 1.0, November 2014                                                                                 '
 write(25,*) '                                                                                                          '
 write(25,*) '**********************************************************************************************************'
 write(25,*) '**********************************************************************************************************'
 write(25,*) ' '
 write(25,*) 'MAAWM Input File Names File: ',trim(T%FILEINPUTFILE)
 write(25,*) 'Atmosphere Input File: ',trim(T%ATMOSPHEREINPUTFILE)
 write(25,*) 'Airplane Input File: ',trim(T%AIRPLANEINPUTFILE)
 write(25,*) 'Quadcopter Input File: ',trim(T%COPTER(1)%INPUTFILE)
 write(25,*) 'Control System Input File: ',trim(T%CSINPUTFILE)
 write(25,*) 'Simulation Input File: ',trim(T%SIMINPUTFILE)
 write(25,*) 'Initial Condition Input File: ',trim(T%ICSINPUTFILE)
 write(25,*) 'Waypoint Input File: ',trim(T%WAYINPUTFILE)
 write(25,*) 'Run Log Output File: ',trim(T%RUNLOGFILE)
 write(25,*) 'State Output File: ',trim(T%STATEOUTPUTFILE)
 write(25,*) 'Miscellaneous Output File: ',trim(T%MISCOUTPUTFILE)
 write(25,*) 'Control Output File: ',trim(T%CONTROLOUTPUTFILE)
 write(25,*) 'Force Output File: ',trim(T%FORCEOUTPUTFILE)
 write(25,*) 'Error Output File: ',trim(T%ERROROUTPUTFILE)
 write(25,*) 'Wind Output File: ',trim(T%WINDOUTPUTFILE)
 write(25,*) ' '

 !Echo All other Data
 write(*,*) 'Echoing atmosphere...'
 call ATMOSPHERE(T,2)
 write(*,*) 'Echoing aircraft...'
 call AIRCRAFT(T,2,1)
 write(*,*) 'Quad handshake init...'
 do i = 1,T%SIM%NOQUADCOPTERS
    write(*,*) 'Initial Handshake...',i,' out of ',T%SIM%NOQUADCOPTERS
    call QUAD_HANDSHAKE_INIT(T,i) !Pass some specific things to quad
 end do
 write(*,*) 'Echoing Copter...'
 call COPTER(T%COPTER(1),2)
 write(*,*) 'Echoing Control...'
 call CONTROL(T,2)
 write(*,*) 'Echoing Simulation...'
 call SIMULATION(T,2)

!!!!!!!!!!!!!!!!!!! Compute Simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (T%SIM%DYNSKIP .eq. 0) then
    call SIMULATION(T,3)
    write(*,*) 'Running Surrogen with RBFDATA = ',T%SIM%RBFDATA
    if (T%SIM%RBFDATA .ne. 0) then
       close(23)
       if (T%SIM%RBFDATA .lt. NDATAMAX) then
          write(*,*) 'Simulation Complete, Computing RBF Fit'
          write(*,*) 'Data Points Sampled = ',T%SIM%NDATA
          T%SIM%RBFDATA = T%SIM%NDATA
          call COMPUTENCENS(T)
          call SURROGEN(T%SIM%IKIND,T%SIM%K35,T%SIM%TIME,T%SIM%NCENS)
       else
          write(*,*) 'Maximum Number of data points allowed for RBF routine exceeded'
          write(*,*) 'Time Simulation = ',T%SIM%TIME
          write(*,*) 'Number of Aircraft = ',T%SIM%NOAIRCRAFT
          write(*,*) 'Number of Quadcopters = ',T%SIM%NOQUADCOPTERS
          STOP
       end if
    else
       write(*,*) 'RBFDATA is zero so most likely it is turned off'
    end if
 else
    write(*,*) 'Dynamic Model Skipped thus using Wind.OUT File instead of simulating data sampling'
    if (T%SIM%RBFDATA .lt. NDATAMAX) then
       !write(*,*) 'Dynamic Model Skipped thus using Wind.OUT File instead of simulating data sampling'
       !Unfortunately we will need to open data file and determine the number of data points
       open(unit=47,file='Output_Files/Winds.OUT',status='old')
       k = 1
       do while (1 .eq. 1) !!!Change this for loop to a while loop
          read(unit=47,fmt=*,iostat=ierr,end=555) x(k,1),x(k,2),x(k,3),x(k,4),u(k,1),u(k,2),u(k,3),ac
          k = k + 1 !!!!ADD THIS k = k + 1 down here
        enddo  
        !!!!!!!!ADD THESE LINES OF CODE
 555    if (ierr .ne. 0) then
           T%SIM%RBFDATA = k-1
           ierr = 0
           close(47)
        endif
        call COMPUTENCENS(T)
       call SURROGEN(T%SIM%IKIND,T%SIM%K35,T%SIM%TIME,T%SIM%NCENS)
    else
       write(*,*) 'Maximum Number of data points allowed for RBF routine exceeded'
    end if
 end if

!!!Or Compute Linearization

 !call LINEARIZEAC(T)

!!!!!!!!!!!!!!!!! Close Output Files !!!!!!!!!!!!!!!!!!!!!!!

 close(25)
 
 STOP

END PROGRAM MAAWM

!!!!!!!!!!!!!!!!!!!! SUBROUTINE SIMULATION !!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SIMULATION(T,iflag)
 use MAAWMDATATYPES
 implicit none
 integer iflag,openflag,readflag
 integer i,j,k,npts,stateindex,a,numcom
 real*8 readreal,xslope,yslope,zslope,idx
 real*8 sum,nominaltime,nominalstate(MAXX),rkalfa(4),krkbody(MAXX,4)
 real*4 elapsed(2)
 real tictotal,ticuser,ticsystem,toctotal,tocuser,tocsystem
 real etime
 real*8 krkactuator(NOACTUATORS,4),nomactuator(NOACTUATORS)
 type(MAAWMSTRUCTURE) T

 !!!!!!!!!!!!!!! COMPUTE iflag = 3 !!!!!!!!!!!!!!!!!!!!!!1

 if (iflag .eq. 3) then


    !Debug Wind Model
    ! T%ATM%XI = -1731.7916154971
    ! T%ATM%YI = 146.1588682476
    ! T%ATM%ZI = -655.8709426696
    ! call ATMOSPHERE(T,3)
    ! write(*,*) T%ATM%VXWIND,T%ATM%VYWIND,T%ATM%VZWIND
    ! STOP

 ! Define Constants

  write(*,*) 'Running Simulation...'

  rkalfa(1) = 1.0; rkalfa(2) = 2.0; rkalfa(3) = 2.0; rkalfa(4) = 1.0

  ! Initial State Vector

  call SYSTEMDERIVATIVES(T,3)
  write(*,*) 'System Derivatives Initialized....'
  call SYSTEMDERIVATIVES(T,2) !Echo data
  write(*,*) 'System Derivatives Echoed'

  write(*,*) 'Initializing Control....'
  if (T%SIM%ICS .eq. 1) then
     call CONTROL(T,3) !Initialize control
  end if

  ! Open Time Simulation Output Files

  open(unit=92,file=T%STATEOUTPUTFILE)
  write(25,*) ' '
  write(25,*) 'STATE VECTOR OUTPUT FILE CREATED: ',trim(T%STATEOUTPUTFILE)
  write(25,*) ' '
  rewind(92)
  open(unit=96,file=T%MISCOUTPUTFILE)
  write(25,*) ' '
  write(25,*) 'MISCELLANEOUS VECTOR OUTPUT FILE CREATED: ',trim(T%MISCOUTPUTFILE)
  write(25,*) ' '
  rewind(96)
  open(unit=91,file=T%CONTROLOUTPUTFILE)
  write(25,*) ' '
  write(25,*) 'CONTROL VECTOR OUTPUT FILE CREATED: ',trim(T%CONTROLOUTPUTFILE)
  write(25,*) ' '
  rewind(91)
  open(unit=83,file=T%FORCEOUTPUTFILE)
  write(25,*) ' '
  write(25,*) 'FORCE VECTOR OUTPUT FILE CREATED: ',trim(T%FORCEOUTPUTFILE)
  write(25,*) ' '
  rewind(83)
  open(unit=42,file=T%ERROROUTPUTFILE)
  write(25,*) ' '
  write(25,*) 'ERROR VECTOR OUTPUT FILE CREATED: ',trim(T%ERROROUTPUTFILE)
  write(25,*) ' '
  rewind(42)
  open(unit=23,file=T%WINDOUTPUTFILE)
  write(25,*) ' '
  write(25,*) 'WIND VECTOR OUTPUT FILE CREATED: ',trim(T%WINDOUTPUTFILE)
  write(25,*) ' '
  rewind(23)
  !write(23,fmt=*) T%SIM%RBFDATA,T%SIM%RBFDATA,T%SIM%RBFDATA,T%SIM%RBFDATA,T%SIM%RBFDATA,T%SIM%RBFDATA,T%SIM%RBFDATA

  write(*,*) 'Output Vectors Created'
  
  ! Integrate Equations of Motion

  tictotal = ETIME(elapsed)
  ticuser = elapsed(1)
  ticsystem = elapsed(2)
  npts = nint((T%SIM%FINALTIME-T%SIM%INITIALTIME)/T%SIM%DELTATIME)

  write(*,*) 'Simulation Start'
  
  do i=1,npts

   T%SIM%IDX = i
  
   ! Output Data to File

   !write(*,*) 'Outputting files....'

   !!OUTPUTS 

   if (mod(i-1,T%SIM%IOUTSKIP) .eq. 0) then 
      !State.OUT File, Misc.OUT File, and Actuator File,Force File
      write(92,fmt='(1000F30.10)',advance='no') T%SIM%TIME
      do a = 1,T%SIM%NOAIRPLANES
         write(92,fmt='(1000F30.10)',advance='no') T%AC(a)%STATE(1:T%SIM%NOSTATES)
      end do
      do a = 1,T%SIM%NOQUADCOPTERS-1
         write(92,fmt='(1000F30.10)',advance='no') T%COPTER(a)%STATE(1:T%SIM%NOSTATES)
      end do
      write(92,fmt='(1000F30.10)') T%COPTER(T%SIM%NOQUADCOPTERS)%STATE(1:T%SIM%NOSTATES)
      write(96,fmt='(1000F30.10)') T%SIM%TIME !Miscellaneous File
      !Control File - outputting aircraft controls of 1 because if NOAIRPLANE is non-zero gauranteed the 1st module is an airplane
      !Omegavec is outputting NOAIRCRAFT because if there are any quads gauranteed the last aircraft is a quad
      write(91,fmt='(1000F30.10)') T%SIM%TIME,T%AC(1)%DELTHRUST,T%AC(1)%RUDDER,T%AC(1)%ELEVATOR,T%AC(1)%AILERON,T%AC(T%SIM%NOAIRCRAFT)%OMEGAVEC
      write(83,fmt='(1000F30.10)') T%SIM%TIME !Force File 
      write(42,fmt='(1000F30.10)') T%SIM%TIME !Error File
   end if 
   
   if (mod(i,npts/10) .eq. 1) then 
      idx = i
      write(*,*) 'Time Simulation: ',int(100*idx/npts)+10, ' % Complete'
   end if

   ! Store Nominal State Values

   nominaltime = T%SIM%TIME
   do a = 1,T%SIM%NOAIRPLANES
      T%AC(a)%NOMINALSTATE(1:T%SIM%NOSTATES) = T%AC(a)%STATE(1:T%SIM%NOSTATES)
   end do

   do a = 1,T%SIM%NOQUADCOPTERS
      T%COPTER(a)%NOMINALSTATE(1:T%SIM%NOQUADSTATES) = T%COPTER(a)%STATE(1:T%SIM%NOQUADSTATES)
   end do

   !Store Nominal Actuators
   if (T%SIM%ACTUATORONOFF .eq. 1) then
      !nomactuator(1) = T%SIM%ACTUATOR(1) 
   end if
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Numerical Integration of Equations of Motion

   !write(*,*) 'T = ',T%SIM%TIME

   do j=1,4

    ! State Values to Evaluate Derivatives

    if (j .ne. 1) then
       T%SIM%TIME = nominaltime + T%SIM%DELTATIME/rkalfa(j)
       do k=1,T%SIM%NOSTATES
          T%AC(1:T%SIM%NOAIRPLANES)%STATE(k) = T%AC(1:T%SIM%NOAIRPLANES)%NOMINALSTATE(k) + T%AC(1:T%SIM%NOAIRPLANES)%KRKBODY(k,j-1)/rkalfa(j)
       end do
       do k = 1,T%SIM%NOQUADSTATES
          T%COPTER(1:T%SIM%NOQUADCOPTERS)%STATE(k) = T%COPTER(1:T%SIM%NOQUADCOPTERS)%NOMINALSTATE(k) + T%COPTER(1:T%SIM%NOQUADCOPTERS)%KRKBODY(k,j-1)/rkalfa(j)
       end do
       !Actuator Dynamics
       if (T%SIM%ACTUATORONOFF .eq. 1) then
          do k = 1,NOACTUATORS
             T%SIM%ACTUATOR(k) = nomactuator(k) + krkactuator(k,j-1)/rkalfa(j)
          end do
       end if
    end if

    ! Compute Derivatives

    !write(*,*) 'Calling Sys Derivs...'

    call SYSTEMDERIVATIVES(T,3)

    ! Runge-Kutta Constants

    do k=1,T%SIM%NOSTATES
     T%AC(1:T%SIM%NOAIRPLANES)%KRKBODY(k,j) = T%SIM%DELTATIME*T%AC(1:T%SIM%NOAIRPLANES)%STATEDOT(k)
    end do
    do k = 1,T%SIM%NOQUADSTATES
       T%COPTER(1:T%SIM%NOQUADCOPTERS)%KRKBODY(k,j) = T%SIM%DELTATIME*T%COPTER(1:T%SIM%NOQUADCOPTERS)%STATEDOT(k)
    end do

    
    if (T%SIM%ACTUATORONOFF .eq. 1) then
       do k=1,NOACTUATORS
          krkactuator(k,j) = T%SIM%DELTATIME*T%SIM%ACTUATORDOT(k)
       end do
    end if

 end do

  !write(*,*) 'RK4 Complete'

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!END RK4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Step Time

   T%SIM%TIME = nominaltime + T%SIM%DELTATIME

   ! Step States
   
   do a=1,T%SIM%NOAIRPLANES
      do j=1,T%SIM%NOSTATES
         sum = 0.0
         do k=1,4
            sum = sum + rkalfa(k)*T%AC(a)%KRKBODY(j,k)
         end do
         T%AC(a)%STATE(j) = T%AC(a)%NOMINALSTATE(j) + sum/6.0
      end do
   end do

   do a=1,T%SIM%NOQUADCOPTERS
      do j=1,T%SIM%NOQUADSTATES
         sum = 0.0
         do k=1,4
            sum = sum + rkalfa(k)*T%COPTER(a)%KRKBODY(j,k)
         end do
         T%COPTER(a)%STATE(j) = T%COPTER(a)%NOMINALSTATE(j) + sum/6.0
      end do
   end do

   !write(*,*) 'Stepped states....'

   ! Step Actuators
   if (T%SIM%ACTUATORONOFF .eq. 1) then
      do j=1,NOACTUATORS
         sum = 0.0
         do k=1,4
            sum = sum + rkalfa(k)*krkactuator(j,k)
         end do
         T%SIM%ACTUATOR(j) = nomactuator(j) + sum/6.0
      end do
   end if

   ! Step State Derivatives

   !write(*,*) 'Clean up Derivs Call'

   call SYSTEMDERIVATIVES(T,3)

   ! Step Control System

   !write(*,*) 'Clean up Control call'

   if (T%SIM%ICS .eq. 1) then
      call CONTROL(T,3)
   end if

  ! State Limits

   ! call STATELIMITS(T)

  end do !End of routine

  toctotal = ETIME(elapsed)
  tocuser = elapsed(1)
  tocsystem = elapsed(2)

  T%SIM%CPUTIMEUSER = tocuser - ticuser
  T%SIM%CPUTIMESYSTEM = tocsystem - ticsystem
  T%SIM%CPUTIMETOTAL = toctotal - tictotal

  write(25,*) ' '
  write(25,*) 'SIMULATION CPU TIME USER (sec): ',T%SIM%CPUTIMEUSER,tocuser,ticuser
  write(25,*) 'SIMULATION CPU TIME SYSTEM (sec): ',T%SIM%CPUTIMESYSTEM,tocsystem,ticsystem
  write(25,*) 'SIMULATION CPU TIME TOTAL (sec): ',T%SIM%CPUTIMETOTAL,toctotal,tictotal
  write(25,*) ' '
 
  write(*,*) ' '
  write(*,*) 'TIME SIMULATION COMPLETE'
  write(*,*) ' '
  write(*,*) 'DATA POINTS = ',T%SIM%NDATA
  write(25,*) ' '
  write(25,*) 'TIME SIMULATION COMPLETE'
  write(25,*) ' '

  ! Close Output Files

  close(92)
  close(96)
  close(91)
  close(83)
  close(42)
  close(23)

  RETURN
  
 end if
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ECHO DATA iflag = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 if (iflag .eq. 2) then
 
  write(25,*) ' '
  write(25,*) 'SSSS IIII M   M U  U LL   AAAA TTTT IIII OOOO N  N'
  write(25,*) 'SS    II  MM MM U  U LL   A  A  TT   II  O  O NN N'
  write(25,*) 'SSSS  II  M M M U  U LL   AAAA  TT   II  O  O N NN'
  write(25,*) '  SS  II  M   M U  U LL   A  A  TT   II  O  O N  N'
  write(25,*) 'SSSS IIII M   M UUUU LLLL A  A  TT  IIII OOOO N  N'
  write(25,*) ' '
  write(25,*) 'Simulation Input File: ',trim(T%SIMINPUTFILE)
  write(25,*) ' '
  write(25,*) 'Skip Dynamic Model On or Off: ',T%SIM%DYNSKIP
  write(25,*) 'Initial Time (s): ',T%SIM%INITIALTIME
  write(25,*) 'Final Time (s): ',T%SIM%FINALTIME
  write(25,*) 'Delta Time (s): ',T%SIM%DELTATIME
  write(25,*) 'Output Skip Parameter (nd): ',T%SIM%IOUTSKIP
  write(25,*) 'Control System Flag (0=Off, 1=On): ',T%SIM%ICS
  write(25,*) 'Waypoint System Flag (0=Off, 1=On): ',T%SIM%CONTROLTYPE
  if (T%SIM%CONTROLTYPE .eq. 0) then
     write(25,*) 'Boundary of World (m): ',T%SIM%BOUNDARY
  end if
  write(25,*) 'Data Points to Sample Before RBF Routine: ',T%SIM%RBFDATA
  write(25,*) 'RBF Order: ',T%SIM%IKIND
  write(25,*) 'Polynomial Order: ',T%SIM%K35
  write(25,*) 'Debug Flag (0=Off, 1=On): ',T%SIM%IDEBUG
  write(25,*) 'Debug Index (nd): ',T%SIM%IDXOUT
  write(25,*) ' '  
  write(25,*) 'Initial Condition Input File: ',trim(T%ICSINPUTFILE)
  write(25,*) ' '  
  write(25,*) 'Initial State Vector'
  write(25,*) '--------------------'
  write(25,*) 'Aircraft States'
  do j = 1,T%SIM%NOAIRPLANES
     write(25,*) 'Aircraft ',j
     do i=1,T%SIM%NOSTATES
        write(25,fmt='(a1,i4,e18.8)') ' ',i,T%AC(j)%INITIALSTATE(i)
     end do
  end do
  do j = 1,T%SIM%NOQUADCOPTERS
     write(25,*) 'Quad ',j
     do i=1,T%SIM%NOQUADSTATES
        write(25,fmt='(a1,i4,e18.8)') ' ',i,T%COPTER(j)%INITIALSTATE(i)
     end do
  end do
  if (T%SIM%CONTROLTYPE .eq. 1) then 
     write(25,*) 'Waypoint Input File: ',trim(T%WAYINPUTFILE)
     write(*,*) 'Waypoint Input File: ',trim(T%WAYINPUTFILE)
     write(25,*) ' '  
     write(25,*) 'Waypoint Vector'
     write(25,*) '--------------------'
     do j = 1,T%SIM%NOAIRPLANES
        write(25,*) 'Aircraft ',j
        numcom = 1
        do while (T%AC(j)%COMMANDFLAG(numcom) .eq. 1)
           write(25,fmt='(1000F30.10)') T%AC(j)%XCOMMAND(numcom),T%AC(j)%YCOMMAND(numcom),T%AC(j)%ZCOMMAND(numcom)
           numcom = numcom + 1
        end do
        T%AC(j)%NUMWAYPOINTS = numcom
        write(*,*) 'Number of Waypoints = ',T%AC(j)%NUMWAYPOINTS
     end do
     do j = 1,T%SIM%NOQUADCOPTERS
        write(*,*) 'Quad ',j
        write(25,*) 'Quad ',j
        numcom = 1
        do while (T%COPTER(j)%COMMANDFLAG(numcom) .eq. 1)
           write(25,fmt='(1000F30.10)') T%COPTER(j)%XCOM(numcom),T%COPTER(j)%YCOM(numcom),T%COPTER(j)%ZCOM(numcom)
           numcom = numcom + 1
        end do
        T%COPTER(j)%NUMWAYPOINTS = numcom
        write(*,*) 'Number of Waypoints = ',T%COPTER(j)%NUMWAYPOINTS
     end do
  end if
  write(25,*) 'Data Quality Flag (nd, 0=Data Not Loaded Successfully, 1=Data Loaded Successfully): ',T%SIM%DQFLAG

  write(*,*) 'Simulation Log Complete'
  
  RETURN
  
 end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOAD DATA iflag = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 if (iflag .eq. 1) then

  !!!12 = aircraft
    T%SIM%NOSTATES = 12
    T%SIM%NOQUADSTATES = 20
  
  open(unit=90,file=T%SIMINPUTFILE,status='old',iostat=openflag)
  if (openflag .ne. 0) then
   write(*,*) 'Error Opening Simulation Input File: ',T%SIMINPUTFILE;  STOP
  end if
  rewind(90)
  
  read(unit=90,fmt=*,iostat=readflag) T%SIM%DYNSKIP
  read(unit=90,fmt=*,iostat=readflag) T%SIM%INITIALTIME
  read(unit=90,fmt=*,iostat=readflag) T%SIM%FINALTIME
  read(unit=90,fmt=*,iostat=readflag) T%SIM%DELTATIME
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%IOUTSKIP = int(readreal)
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%ICS = int(readreal)
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%CONTROLTYPE = int(readreal)
  read(unit=90,fmt=*,iostat=readflag) T%SIM%BOUNDARY
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%RBFDATA = int(readreal)
  if (T%SIM%RBFDATA .gt. NDATAMAX) then
     write(*,*) 'Maximum Number of data points allowed for RBF routine exceeded - check the SIM inputfile'
     STOP
  end if
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%IKIND = int(readreal)
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%K35 = int(readreal)

  !Compute K35 
  if (T%SIM%K35 .eq. 4) then
     T%SIM%K35 = 35
  else if (T%SIM%K35 .eq. 3) then
     T%SIM%K35 = 20
  else if (T%SIM%K35 .eq. 2) then
     T%SIM%K35 = 10
  else if (T%SIM%K35 .eq. 1) then
     T%SIM%K35 = 4
  else
     T%SIM%K35 = 0
  end if

  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%ACTUATORONOFF = int(readreal)
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%IDEBUG = int(readreal)
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%NCENS = int(readreal)

  !Initialize Time Vector
  T%SIM%TIME = T%SIM%INITIALTIME

  close(90)

  !!!Read Initial Aircraft States
  open(unit=90,file=T%ICSINPUTFILE,status='old',iostat=openflag)
  if (openflag .ne. 0) then
   write(*,*) 'Error Opening Initial Condition Input File: ',T%ICSINPUTFILE;  STOP
  end if
  rewind(90)

  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%NOAIRCRAFT = int(readreal)

  if (T%SIM%NOAIRCRAFT .gt. MAXAIRCRAFT) then
     write(*,*) 'Number of aircraft requested = ',T%SIM%NOAIRCRAFT
     write(*,*) 'Maximum number of aircraft allowed = ',MAXAIRCRAFT
     STOP;
  end if

  !Ok here we are going to add the number of aircraft and quads
  read(unit=90,fmt=*,iostat=readflag) readreal; T%SIM%NOAIRPLANES = int(readreal)

  if (T%SIM%NOAIRPLANES .gt. T%SIM%NOAIRCRAFT) then
     write(*,*) 'Truncating Airplanes to NOAIRCRAFT and Quadcopters to zero'
     write(*,*) 'Re-check NOAIRCRAFT and NOAIRPLANES'
     T%SIM%NOAIRPLANES = T%SIM%NOAIRCRAFT
     T%SIM%NOQUADCOPTERS = 0
  else
     T%SIM%NOQUADCOPTERS = T%SIM%NOAIRCRAFT - T%SIM%NOAIRPLANES
  end if

  write(*,*) 'Number of Airplanes = ',T%SIM%NOAIRPLANES
  write(*,*) 'Number of Quadcopters = ',T%SIM%NOQUADCOPTERS

  !Copy Airplane Structures into other structure for now (Assume all mass and geometry properties are identical
  if (T%SIM%NOAIRPLANES .gt. 1) then
     do i = 2,T%SIM%NOAIRPLANES
        T%AC(i) = T%AC(1)
     end do
  end if

  !Do the same for the COPTER structure
  if (T%SIM%NOQUADCOPTERS .gt. 1) then
     do i = 2,T%SIM%NOQUADCOPTERS
        T%COPTER(i) = T%COPTER(1)
     end do
  end if
  
  !!!!Read Aircraft Initial States!!!
  do i = 1,T%SIM%NOAIRPLANES
     read(unit=90,fmt=*,iostat=readflag) readreal; !Dummy variable
     do j=1,12
        read(unit=90,fmt=*,iostat=readflag) T%AC(i)%INITIALSTATE(j)
     end do
     T%AC(i)%STATE = T%AC(i)%INITIALSTATE
  end do
  !So we need to all all aircraft first and then all quad states
  do i = 1,T%SIM%NOQUADCOPTERS
     read(unit=90,fmt=*,iostat=readflag) readreal; !Dummy variable
     do j=1,12
        read(unit=90,fmt=*,iostat=readflag) T%COPTER(i)%INITIALSTATE(j)
     end do
     T%COPTER(i)%STATE = T%COPTER(i)%INITIALSTATE
  end do
  close(90)

  !!!Read Waypoint File
  if (T%SIM%CONTROLTYPE .eq. 1) then
     open(unit=90,file=T%WAYINPUTFILE,status='old',iostat=openflag)
     if (openflag .ne. 0) then
        write(*,*) 'Error Opening Waypoint Input File: ',T%WAYINPUTFILE;  STOP
     end if
     rewind(90)

     read(unit=90,fmt=*,iostat=readflag) readreal; !Dummy Row

     do a = 1,T%SIM%NOAIRPLANES
        read(unit=90,fmt=*,iostat=readflag) readreal; i = int(readreal); !Aircraft Number
        j = 1 !Initialize row counter
        read(unit=90,fmt=*,iostat=readflag) T%AC(i)%XCOMMAND(j),T%AC(i)%YCOMMAND(j),T%AC(i)%ZCOMMAND(j),T%AC(i)%COMMANDFLAG(j)
        do while (T%AC(i)%COMMANDFLAG(j) .eq. 1)
           j = j+1
           if (j .gt. MAXCOMMAND) then
              write(*,*) 'Maximum number of waypoints exceeded. Maximum allowed is ',MAXCOMMAND
              write(*,*) 'Check your Input File: ',T%WAYINPUTFILE; STOP;
           end if
           read(unit=90,fmt=*,iostat=readflag) T%AC(i)%XCOMMAND(j),T%AC(i)%YCOMMAND(j),T%AC(i)%ZCOMMAND(j),T%AC(i)%COMMANDFLAG(j)
        end do
     end do

     do a = 1,T%SIM%NOQUADCOPTERS
        read(unit=90,fmt=*,iostat=readflag) readreal; i = int(readreal); !Quad Number
        i = i - T%SIM%NOAIRPLANES
        j = 1 !Initialize row counter
        read(unit=90,fmt=*,iostat=readflag) T%COPTER(i)%XCOM(j),T%COPTER(i)%YCOM(j),T%COPTER(i)%ZCOM(j),T%COPTER(i)%COMMANDFLAG(j)
        do while (T%COPTER(i)%COMMANDFLAG(j) .eq. 1)
           j = j+1
           if (j .gt. MAXCOMMAND) then
              write(*,*) 'Maximum number of waypoints exceeded. Maximum allowed is ',MAXCOMMAND
              write(*,*) 'Check your Input File: ',T%WAYINPUTFILE; STOP;
           end if
           read(unit=90,fmt=*,iostat=readflag) T%COPTER(i)%XCOM(j),T%COPTER(i)%YCOM(j),T%COPTER(i)%ZCOM(j),T%COPTER(i)%COMMANDFLAG(j)
        end do
     end do
     
     close(90)
  else
     !This is the wander aimlessly algorithm
     do a = 1,T%SIM%NOAIRPLANES
        T%AC(a)%XCOMMAND(1) = T%AC(a)%INITIALSTATE(1) + 4*T%SIM%BOUNDARY*cos(T%AC(a)%INITIALSTATE(6))
        T%AC(a)%YCOMMAND(1) = T%AC(a)%INITIALSTATE(2) + 4*T%SIM%BOUNDARY*sin(T%AC(a)%INITIALSTATE(6))
     end do
     do a = 1,T%SIM%NOQUADCOPTERS
        T%COPTER(a)%XCOM(1) = T%COPTER(a)%INITIALSTATE(1)
        T%COPTER(a)%YCOM(1) = T%COPTER(a)%INITIALSTATE(2)
     end do
  end if
  
  T%SIM%DQFLAG = 1

  ! State Limits
  
  ! call STATELIMITS(T)

  write(*,*) 'SIMULATION Load Complete'

  RETURN
 
 end if
 
 RETURN
END SUBROUTINE SIMULATION

!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINE SYSTEMDERIVATIVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SYSTEMDERIVATIVES(T,iflag)
 use MAAWMDATATYPES
 implicit none
 integer stateindex,iflag,i,j
 type(MAAWMSTRUCTURE) T

!!!!!!!!!!!!!!!!!!!!!!!! COMPUTE iflag = 3!!!!!!!!!!!!!!!!!

if (iflag .eq. 3) then
   
    stateindex = 0 
    
    ! Aircraft State Derivatives
    
    do i = 1,T%SIM%NOAIRPLANES
       T%AC(i)%STATEDOT = 0.0
       if (T%AC(i)%DYNOFFON .eq. 1) then
          call AIRCRAFT(T,3,i)
       end if
    end do

    ! Quadcopter State Derivatives

    do i = 1,T%SIM%NOQUADCOPTERS
       T%COPTER(i)%STATEDOT = 0.0
       if (T%COPTER(i)%DYNOFFON .eq. 1) then
          call QUAD_HANDSHAKE(T,i)
          call COPTER(T%COPTER(i),3)
       end if
    end do

    RETURN
    
 end if

!!!!!!!!!!!!!!!!!!!!ECHO DATA iflag = 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (iflag .eq. 2) then
   
   write(25,*) ' '
   write(25,*) 'STATE AT T = ',T%SIM%TIME
   write(25,*) ' '
   do j = 1,T%SIM%NOAIRPLANES
      write(25,*) 'Airplane ',j
      do i=1,T%SIM%NOSTATES
         write(25,fmt='(a1,i4,e18.8)') ' ',i,T%AC(j)%STATE(i)
      end do
   end do
   do j = 1,T%SIM%NOQUADCOPTERS
      write(25,*) 'Quadcopter ',j
      do i=1,T%SIM%NOQUADSTATES
         write(25,fmt='(a1,i4,e18.8)') ' ',i,T%COPTER(j)%STATE(i)
      end do
   end do
   write(25,*) ' '
   write(25,*) 'STATE DERIVATIVES AT T = ',T%SIM%TIME
   write(25,*) ' '
   do j = 1,T%SIM%NOAIRPLANES
      write(25,*) 'Airplane ',j
      do i=1,T%SIM%NOSTATES
         write(25,fmt='(a1,i4,e18.8)') ' ',i,T%AC(j)%STATEDOT(i)
      end do
   end do
   do j = 1,T%SIM%NOQUADCOPTERS
      write(25,*) 'Quadcopter ',j
      do i=1,T%SIM%NOQUADSTATES
         write(25,fmt='(a1,i4,e18.8)') ' ',i,T%COPTER(j)%STATEDOT(i)
      end do
   end do
   write(25,*) ' '

   RETURN
            
end if
 
 RETURN
END SUBROUTINE SYSTEMDERIVATIVES

!!!!!!!!!!!!!!!!!!!1!!! SUBROUTINE STATELIMITS!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SUBROUTINE STATELIMITS(T)
!  use MAAWMDATATYPES
!  implicit none
!  integer i,stateindex,counter
!  real*8 q0,q1,q2,q3,norm
!  type(MAAWMSTRUCTURE) T

 ! Normalize Aircraft Quaternions

 ! if (T%AC%DYNOFFON .eq. 1) then
 !  T%AC%STATE(1:12) = T%SIM%STATE(1:13)
 !  q0 = T%AC%STATE(4)
 !  q1 = T%AC%STATE(5)
 !  q2 = T%AC%STATE(6)
 !  q3 = T%AC%STATE(7)
 !  norm = sqrt(q0**2+q1**2+q2**2+q3**2)
 !  T%AC%STATE(4) = q0/norm
 !  T%AC%STATE(5) = q1/norm
 !  T%AC%STATE(6) = q2/norm
 !  T%AC%STATE(7) = q3/norm  
 !  T%SIM%STATE(1:13) = T%AC%STATE(1:13)
 ! end if

!  RETURN
! END SUBROUTINE STATELIMITS

!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINE CONTROL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CONTROL(T,iflag)
 use MAAWMDATATYPES
 implicit none
 integer iflag,openflag,readflag,i,stateindex,now(3),ac,j
 real*8 readreal,phicommand,deltheta,delpsi,zcommand,rAS_I(3,1),vAS_I(3,1),thetacommand
 real*8 rAS_A(3,1),vAS_A(3,1),lencommand,rGS_I(3,1),vGS_I(3,1),psicp,rGS_P(3,1),vGS_P(3,1),tension,len
 real*8 angles(2,1),delphi,delphidot,ldotnom,ldot,n1,x,y,z,q0,q1,q2,q3,phi,theta,psi,u,v,p,xdot,ydot,zdot
 real*8 normAC,rcx,rcy,lbarx,lbary,flightdirection,planedirection,xs,ys
 real*8 spsi,cpsi,spsic,cpsic,delx,dely,delydot,xcommand,ycommand,delz,dist,dx,dy
 real*8 xother,yother,zother,center(2,1),thetacirc,ddotcirc(2,1),radii,ways
 real*8 omegaNot,omegaRight,omegaOpp,omegaBack,omegaFront,omegaDiag,domegaLeft,domegaFront
 real*8 domegaDiag,addyaw,addroll,addpitch,q,r
 type(MAAWMSTRUCTURE) T
 
!!!!!!!!!!!!!!!!!!!!!!! COMPUTE iflag = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!
  
 if (iflag .eq. 3) then

    !write(*,*) 'Compute controls....'

    !!!Feedback from Sensors
    call FEEDBACK(T) 

    !write(*,*) 'Feedback routine complete'

    !!!AIRCRAFT CONTROL SYSTEM
    if (T%CS%AIRCRAFTCONTROLOFFON .eq. 1) then
       do i = 1,T%SIM%NOAIRPLANES

          !!Extract State vector
          x = T%AC(i)%STATE(1)
          y = T%AC(i)%STATE(2)
          z = T%AC(i)%STATE(3)
          phi = T%AC(i)%STATE(4)
          theta = T%AC(i)%STATE(5)
          psi = T%AC(i)%STATE(6)
          u = T%AC(i)%STATE(7)
          v = T%AC(i)%STATE(8)
          p = T%AC(i)%STATE(10)
          q = T%AC(i)%STATE(11)
          r = T%AC(i)%STATE(12)

          !!Extract Statedot Vector
          xdot = T%AC(i)%STATEDOT(1)
          ydot = T%AC(i)%STATEDOT(2)
          zdot = T%AC(i)%STATEDOT(3)
          
          !!COMMANDS
          
          !!!The Controltype parameter takes two values. 0 = Wander Aimlessly
          !!!1 = Follow Waypoints
          if (T%SIM%CONTROLTYPE .eq. 1) then

             !!!!WAYPOINTS 
             xcommand = T%AC(i)%XCOMMAND(T%AC(i)%WAYPOINT)
             ycommand = T%AC(i)%YCOMMAND(T%AC(i)%WAYPOINT)
             zcommand = T%AC(i)%ZCOMMAND(T%AC(i)%WAYPOINT)
             
             delx = xcommand - T%AC(i)%STATE(1)
             dely = ycommand - T%AC(i)%STATE(2)
             delz = zcommand - T%AC(i)%STATE(3)
             
             !!!Check for Advancing Waypoint 
             if (sqrt(delx**2+dely**2+delz**2) .lt. T%AC(i)%MINWAYPOINT) then
                ways = 0
                do ac = 1,T%SIM%NOAIRCRAFT
                   ways = ways + T%AC(ac)%COMMANDFLAG(T%AC(ac)%WAYPOINT)
                end do
                if ((ways .eq. 0) .and. (T%SIM%RBFDATA .eq. -1)) then
                   write(*,*) 'Last Waypoint Reached'
                   T%SIM%RBFDATA = T%SIM%NDATA+1
                end if
                !write(*,*) 'Waypoint ',T%AC(i)%WAYPOINT,' out of ',T%AC(i)%NUMWAYPOINTS,' Reached'
                T%AC(i)%WAYPOINT = T%AC(i)%WAYPOINT + T%AC(i)%COMMANDFLAG(T%AC(i)%WAYPOINT)
             end if
          else
             !!Wander Aimlessly
             if (T%AC(i)%BFLAG .eq. 1) then
                !I am currently out of bounds check to see if I finally am in bounds again
                if ((abs(T%AC(i)%STATE(1)) .lt. T%SIM%BOUNDARY-T%AC(i)%TURNRADIUS) .and. (abs(T%AC(i)%STATE(2)) .lt. T%SIM%BOUNDARY-T%AC(i)%TURNRADIUS)) then
                   T%AC(i)%BFLAG = 0
                end if
                !Double check in the event of something bad
                if ((abs(T%AC(i)%STATE(2)) .gt. T%SIM%BOUNDARY) .or. (abs(T%AC(i)%STATE(1)) .gt. T%SIM%BOUNDARY)) then
                   T%AC(i)%XCOMMAND(1) = 0
                   T%AC(i)%YCOMMAND(1) = 0
                end if
                if ((abs(T%AC(i)%STATE(2)) .gt. T%SIM%BOUNDARY-T%AC(i)%TURNRADIUS) .and. (abs(T%AC(i)%STATE(1)) .gt. T%SIM%BOUNDARY-T%AC(i)%TURNRADIUS)) then
                   T%AC(i)%XCOMMAND(1) = 0
                   T%AC(i)%YCOMMAND(1) = 0
                end if
             else
                !Check if something bad happened
                if ((T%AC(i)%XCOMMAND(1) .eq. 0) .and. (T%AC(i)%YCOMMAND(1) .eq. 0)) then
                   !Recompute waypoint
                   T%AC(i)%XCOMMAND(1) = T%AC(i)%STATE(1) + 4*T%SIM%BOUNDARY*cos(T%AC(i)%STATE(6))
                   T%AC(i)%YCOMMAND(1) = T%AC(i)%STATE(2) + 4*T%SIM%BOUNDARY*sin(T%AC(i)%STATE(6))
                end if
                !I am currently in bounds, check and see if I travel out of bounds
                if (abs(T%AC(i)%STATE(1)) .gt. T%SIM%BOUNDARY-T%AC(i)%TURNRADIUS) then
                   T%AC(i)%BFLAG = 1
                   !Set a new waypoint
                   T%AC(i)%XCOMMAND(1) = T%AC(i)%STATE(1) - 4*T%SIM%BOUNDARY*cos(T%AC(i)%STATE(6))
                   T%AC(i)%YCOMMAND(1) = T%AC(i)%STATE(2) + 4*T%SIM%BOUNDARY*sin(T%AC(i)%STATE(6))
                end if
                if (abs(T%AC(i)%STATE(2)) .gt. T%SIM%BOUNDARY-T%AC(i)%TURNRADIUS) then
                   T%AC(i)%BFLAG = 1
                   !Set a new waypoint
                   T%AC(i)%XCOMMAND(1) = T%AC(i)%STATE(1) + 4*T%SIM%BOUNDARY*cos(T%AC(i)%STATE(6))
                   T%AC(i)%YCOMMAND(1) = T%AC(i)%STATE(2) - 4*T%SIM%BOUNDARY*sin(T%AC(i)%STATE(6))
                end if
             end if

             xcommand = T%AC(i)%XCOMMAND(1)
             ycommand = T%AC(i)%YCOMMAND(1)
             zcommand = T%AC(i)%INITIALSTATE(3)

             delx = xcommand - T%AC(i)%STATE(1)
             dely = ycommand - T%AC(i)%STATE(2)
             delz = zcommand - T%AC(i)%STATE(3)

          end if

          !!Here is the avoidance algorithm
          !if (T%SIM%CONTROLTYPE .gt. -1) then
          if (T%SIM%CONTROLTYPE .eq. -99) then !Turn off avoidance for now
             !!Check to make sure you aren't to close to another aircraft
             do ac = 1,T%SIM%NOAIRCRAFT
                if (ac .ne. i) then
                   xother = T%AC(ac)%STATE(1)
                   yother = T%AC(ac)%STATE(2)
                   zother = T%AC(ac)%STATE(3)
                   dx = xother - x 
                   dy = yother - y
                   dist = sqrt((dx)**2 + (dy)**2)
                   if ((dist .lt. T%AC(i)%TURNRADIUS) .and. (abs(z-zother) .lt. 10)) then
                      center(1,1) = (x+xother)/2
                      center(2,1) = (y+yother)/2
                      thetacirc = atan2(y-center(2,1),x-center(1,1))
                      ! T%AC(i)%XCOMMAND(1) = x + 4*T%SIM%BOUNDARY*cos(thetacirc)
                      ! T%AC(i)%YCOMMAND(1) = y + 4*T%SIM%BOUNDARY*sin(thetacirc)
                      ! xcommand = T%AC(i)%XCOMMAND(1)
                      ! ycommand = T%AC(i)%YCOMMAND(1)
                      xcommand = x + 4*T%SIM%BOUNDARY*cos(thetacirc)
                      ycommand = y + 4*T%SIM%BOUNDARY*sin(thetacirc)                      
                      zcommand = T%AC(i)%INITIALSTATE(3)

                      delx = xcommand - T%AC(i)%STATE(1)
                      dely = ycommand - T%AC(i)%STATE(2)
                      delz = zcommand - T%AC(i)%STATE(3)
                      
                      !Recompute waypoint so that when avoidance is over you can keep going straight.
                      !I'm pretty sure this is only used for the wander aimlessly algorithm
                      !Ok actually this is used during any algorithm and causes some wierd things if aircraft come in contact.
                      !However when this was written it was not indended to be used during waypoint following
                      T%AC(i)%XCOMMAND(1) = T%AC(i)%STATE(1) + 4*T%SIM%BOUNDARY*cos(T%AC(i)%STATE(6))
                      T%AC(i)%YCOMMAND(1) = T%AC(i)%STATE(2) + 4*T%SIM%BOUNDARY*sin(T%AC(i)%STATE(6))
                      ! if (i .eq. 1) then
                      !    write(*,*) 'T,Xcommand,Ycommand = ',T%SIM%TIME,T%AC(i)%XCOMMAND(1), T%AC(i)%YCOMMAND(1)
                      ! end if
                   end if
                end if
             end do
          end if

          !Airplane control
          !!Velocity Command - This is only for the Aircraft but it doesn't really
          !matter because the quadcopter won't use it
          T%AC(i)%UCOMMAND = T%AC(i)%V_T
          T%AC(i)%PSICOMMAND = atan(dely,delx)
          !!Thrust Controller
          T%AC(i)%DELTHRUST = T%CS%KU*(T%AC(i)%UCOMMAND-u)
          !!!Sideslip Controller - 
          T%AC(i)%RUDDER = T%CS%KV*v
          !!!Altitude Hold
          deltheta = T%CS%KPZ*(zcommand-z) + T%CS%KDZ*zdot
          thetacommand = theta+deltheta
          if (abs(thetacommand) .gt. 15*PI/180) then
             thetacommand = sign(15*PI/180,thetacommand)
          end if
          T%AC(i)%ELEVATOR = T%CS%KPT*(thetacommand-theta)
          !!!Heading Controller
          spsi = sin(T%AC(i)%STATE(6))
          cpsi = cos(T%AC(i)%STATE(6))
          spsic = sin(T%AC(i)%PSICOMMAND)
          cpsic = cos(T%AC(i)%PSICOMMAND)
          delpsi = atan2(spsi*cpsic-cpsi*spsic,cpsi*cpsic+spsi*spsic);
          if (abs(delpsi) .gt. 45*pi/180) then
             delpsi = sign(45*pi/180,delpsi)
          end if
          phicommand = T%CS%KPSI*delpsi
          if (abs(phicommand) .gt. 45*PI/180) then
             phicommand = sign(45*PI/180,phicommand)
          end if
          T%AC(i)%AILERON = T%CS%KPP*(phi-phicommand)+T%CS%KDP*p
          !Check for saturation of control surfaces
          if (abs(T%AC(i)%RUDDER) .gt. 15*PI/180) then
             T%AC(i)%RUDDER = sign(15*PI/180,T%AC(i)%RUDDER);
          end if
          if (abs(T%AC(i)%ELEVATOR) .gt. 20*PI/180) then
             T%AC(i)%ELEVATOR = sign(20*PI/180,T%AC(i)%ELEVATOR);
          end if
          if (abs(T%AC(i)%AILERON) .gt. 20*PI/180) then
             T%AC(i)%AILERON = sign(20*PI/180,T%AC(i)%AILERON)
          end if
       end do !End do NOAIRPLANES

       !Now call the control system for the quad. For now we will
       !just have them follow waypoint
       ways = 0
       do i = 1,T%SIM%NOQUADCOPTERS
          !write(*,*) 'Current Waypoint = ',T%COPTER(i)%XCOM(T%COPTER(i)%WAYPOINT),T%COPTER(i)%YCOM(T%COPTER(i)%WAYPOINT),T%COPTER(i)%ZCOM(T%COPTER(i)%WAYPOINT)
          
          delx = T%COPTER(i)%XCOM(T%COPTER(i)%WAYPOINT) - T%COPTER(i)%STATE(1)
          dely = T%COPTER(i)%YCOM(T%COPTER(i)%WAYPOINT) - T%COPTER(i)%STATE(2)
          delz = T%COPTER(i)%ZCOM(T%COPTER(i)%WAYPOINT) - T%COPTER(i)%STATE(3)
             
          !!!Check for Advancing Waypoint 
          if (sqrt(delx**2+dely**2+delz**2) .lt. T%COPTER(i)%MINWAYPOINT) then
             T%COPTER(i)%XINTEGRAL = 0.0
             T%COPTER(i)%YINTEGRAL = 0.0
             T%COPTER(i)%ZINTEGRAL = 0.0
             !write(*,*) 'Waypoint ',T%AC(i)%WAYPOINT,' out of ',T%AC(i)%NUMWAYPOINTS,' Reached'
             T%COPTER(i)%WAYPOINT = T%COPTER(i)%WAYPOINT + T%COPTER(i)%COMMANDFLAG(T%COPTER(i)%WAYPOINT)
          end if
          ways = ways + T%COPTER(i)%COMMANDFLAG(T%COPTER(i)%WAYPOINT)
          T%COPTER(i)%XCOMMAND = T%COPTER(i)%XCOM(T%COPTER(i)%WAYPOINT)
          T%COPTER(i)%YCOMMAND = T%COPTER(i)%YCOM(T%COPTER(i)%WAYPOINT)
          T%COPTER(i)%ZCOMMAND = T%COPTER(i)%ZCOM(T%COPTER(i)%WAYPOINT)
          call COPTER_CONTROL(T%COPTER(i))
       end do
       if (T%SIM%NOQUADCOPTERS .gt. 0) then
          if ((ways .eq. 0) .and. (T%SIM%RBFDATA .eq. -1)) then
             write(*,*) 'Last Quadcopter Waypoint Reached'
             T%SIM%RBFDATA = T%SIM%NDATA+1
          end if
       end if
    end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RETURN
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ECHO DATA iflag = 2 !!!!!!!!!!!!!!!!!!!!!!!!

 if (iflag .eq. 2) then 
  write(25,*) ' '
  write(25,*) 'CCCC OOOO N  N TTTT RRR  OOOO LL   SSSS'
  write(25,*) 'CC   O  O NN N  TT  RR R O  O LL   SS  '
  write(25,*) 'CC   O  O N NN  TT  RRR  O  O LL   SSSS'
  write(25,*) 'CC   O  O N  N  TT  RRR  O  O LL     SS'
  write(25,*) 'CCCC OOOO N  N  TT  RR R OOOO LLLL SSSS'
  write(25,*) ' '
  write(25,*) 'Control System Input File: '
  write(25,*) trim(T%CSINPUTFILE)
  write(25,*) ' '
  write(25,*) 'Airplane Control (0=Off, 1=On): ',T%CS%AIRCRAFTCONTROLOFFON
  write(25,*) 'Airplane proportional gain on speed: ',T%CS%KU 
  write(25,*) 'Airplane proportional gain on roll angle: ',T%CS%KPP
  write(25,*) 'Airplane derivative gain on roll angle: ',T%CS%KDP
  write(25,*) 'Airplane proportional gain on pitch angle: ',T%CS%KPT
  write(25,*) 'Airplane proportional gain on altitude: ',T%CS%KPZ
  write(25,*) 'Airplane derivative gain on altitude: ',T%CS%KDZ
  write(25,*) 'Airplane proportional gain on crossrange: ',T%CS%KPY
  write(25,*) 'Airplane derivative gain on crossrange: ',T%CS%KDY
  write(25,*) 'Airplane proportional gain on sideslip: ',T%CS%KV 
  write(25,*) 'Airplane proportional gain on heading: ',T%CS%KPSI
  write(25,*) 'Feedback Update Rate: ',T%CS%UPDATERATE
  write(25,*) 'Sensor Errors On (1) off (0): ', T%CS%SENSORERRORS
  write(25,*) 'Data Quality Flag (nd, 0=Data Not Loaded Successfully, 1=Data Loaded Successfully): ',T%CS%DQFLAG

  RETURN
  
 end if
  
!!!!!!!!!!!!!!!!!!!!!!!!LOAD DATA iflag = 1!!!!!!!!!!!!!!!!!!!!!!!!!
 
 if (iflag .eq. 1) then
 
  open(unit=94,file=T%CSINPUTFILE,status='old',iostat=openflag)
  if (openflag .ne. 0) then
   write(*,*) 'Error Opening Control System Input File: ',T%CSINPUTFILE;  STOP
  end if
  rewind(94)

  read(unit=94,fmt=*,iostat=readflag) readreal; T%CS%AIRCRAFTCONTROLOFFON = readreal
  read(unit=94,fmt=*,iostat=readflag) T%CS%KPY
  read(unit=94,fmt=*,iostat=readflag) T%CS%KDY
  read(unit=94,fmt=*,iostat=readflag) T%CS%KPSI
  read(unit=94,fmt=*,iostat=readflag) T%CS%KPP
  read(unit=94,fmt=*,iostat=readflag) T%CS%KDP
  read(unit=94,fmt=*,iostat=readflag) T%CS%KPZ
  read(unit=94,fmt=*,iostat=readflag) T%CS%KDZ
  read(unit=94,fmt=*,iostat=readflag) T%CS%KPT
  read(unit=94,fmt=*,iostat=readflag) T%CS%KV 
  read(unit=94,fmt=*,iostat=readflag) T%CS%KU
  read(unit=94,fmt=*,iostat=readflag) T%CS%UPDATERATE
  !!Sensor Errors
  read(unit=94,fmt=*,iostat=readflag) T%CS%SENSORERRORS

  call itime(now)

  do i = 1,now(1)*3600+now(2)*60+now(3)
     call RandUniform(n1)  
  end do

  ! T%CS%TENSIONERROR(1) = T%CS%TENSIONERROR(1)*n1
  ! call RandUniform(n1)  
  ! T%CS%TENSIONERROR(2) = T%CS%TENSIONERROR(2)*n1
  ! call RandUniform(n1)  
  ! T%CS%BIRDROLLERROR(1) = T%CS%BIRDROLLERROR(1)*n1
  ! call RandUniform(n1)  
  ! T%CS%BIRDROLLERROR(2) = T%CS%BIRDROLLERROR(2)*n1
  ! call RandUniform(n1)  
  ! T%CS%RISERROLLERROR(1) = T%CS%RISERROLLERROR(1)*n1
  ! call RandUniform(n1)  
  ! T%CS%RISERROLLERROR(2) = T%CS%RISERROLLERROR(2)*n1
  
  close(94) 

  write(*,*) 'CONTROL SYSTEM Load Complete'
  
  T%CS%DQFLAG = 1

  RETURN
 
 end if

 RETURN
END SUBROUTINE CONTROL

!!!!!!!!!!!!!!!!!!!1!!! SUBROUTINE ATMOSPHERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ATMOSPHERE(T,iflag)
 use MAAWMDATATYPES
 implicit none
 integer i,iflag,ifind,openflag,readflag,ii,jj,ierr
 real*8 m,readreal,wind,dirnominal,windnominal,winddir
 real*8 dim,rx,test(3),ramp
 character*11 HeightFile
 character*14 ParametersFile
 character*1 letter
 character*10 ZcoordFile,number
 character*15 TimesFile
 character*256 inParameters,inZcoord,inTimes,inHeight,temp
 type(MAAWMSTRUCTURE) T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTE iflag = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (iflag .eq. 3) then

    !write(*,*) 'Atmosphere at (x,y,z) = ',T%ATM%XI,T%ATM%YI,T%ATM%ZI

    !!Reset all winds to zero
    T%ATM%VXWIND = 0
    T%ATM%VYWIND = 0
    T%ATM%VZWIND = 0
    T%ATM%WRFX = 0
    T%ATM%WRFY = 0
    T%ATM%WRFZ = 0
    T%ATM%WINDGUST(1) = 0
    T%ATM%WINDGUST(2) = 0
    T%ATM%WINDGUST(3) = 0

    T%ATM%ALT = -T%ATM%ZI

  ! Constant Model

  !!!!!REVISIT REVISIT REVISIT!!!!!!
    dirnominal = T%ATM%WINDDIR
    windnominal = 0.0
    if ((T%ATM%WINDDIR .ne. dirnominal) .or. (T%ATM%WINDSPEED .ne. windnominal)) then
       if (T%SIM%TIME .gt. 1) then
          winddir = T%ATM%WINDDIR
          wind = T%ATM%WINDSPEED
       else
          wind = windnominal + (T%SIM%TIME/1)*(T%ATM%WINDSPEED-windnominal)
          winddir = dirnominal + (T%SIM%TIME/1)*(T%ATM%WINDDIR-dirnominal)     
       endif
    else
       wind = windnominal
       winddir = dirnominal
    end if
    
    if (T%ATM%MODNO .eq. 1) then
       T%ATM%VXWIND = wind*cos(winddir)*cos(T%ATM%WINDELEV)
       T%ATM%VYWIND = wind*sin(winddir)*cos(T%ATM%WINDELEV)
       T%ATM%VZWIND = -wind*sin(T%ATM%WINDELEV)
       T%ATM%DEN = 0.002363;
    end if

  ! Equation Density and Constant Wind Model

  if (T%ATM%MODNO .eq. 2) then
     T%ATM%DEN = 0.002363*(1.00000000-0.0000225696709*T%ATM%ALT)**4.258
     T%ATM%VXWIND = T%ATM%WINDSPEED*cos(T%ATM%WINDDIR)
     T%ATM%VYWIND = T%ATM%WINDSPEED*sin(T%ATM%WINDDIR)
     T%ATM%VZWIND = 0.0
  end if

  ! Table Look-Up Model 

  if (T%ATM%MODNO .eq. 3) then

   ifind = 0

   ! Position Pointer  

   if (T%ATM%ALT .le. T%ATM%ALTTAB(T%ATM%IP)) then 
    ifind = -1 
    do while ((ifind.ne.0) .and. (T%ATM%IP.gt.1))
     T%ATM%IP = T%ATM%IP - 1 
     if (T%ATM%ALTTAB(T%ATM%IP)   .le. T%ATM%ALT) then 
     if (T%ATM%ALTTAB(T%ATM%IP+1) .gt. T%ATM%ALT) then 
      ifind = 0
     end if
     end if 
    end do 
   end if
   if (T%ATM%ALT .gt. T%ATM%ALTTAB(T%ATM%IP+1)) then 
    ifind = 1
    do while ((ifind.ne.0) .and. (T%ATM%IP.lt.T%ATM%TABSIZE-1))
     T%ATM%IP = T%ATM%IP + 1 
     if (T%ATM%ALTTAB(T%ATM%IP)   .le. T%ATM%ALT) then 
     if (T%ATM%ALTTAB(T%ATM%IP+1) .gt. T%ATM%ALT) then 
      ifind = 0 
     end if 
     end if 
    end do 
   end if
   if (ifind .eq. 0) then
    m = (T%ATM%ALT-T%ATM%ALTTAB(T%ATM%IP))/(T%ATM%ALTTAB(T%ATM%IP+1)-T%ATM%ALTTAB(T%ATM%IP))
   else if (ifind .eq. -1) then
    m = 0.0
   else if (ifind .eq. 1) then
    m = 1.0
   end if

   ! Interpolate

   T%ATM%DEN    = T%ATM%DENTAB(T%ATM%IP)    + m*(T%ATM%DENTAB(T%ATM%IP+1)-T%ATM%DENTAB(T%ATM%IP))
   T%ATM%VXWIND = T%ATM%VXWINDTAB(T%ATM%IP) + m*(T%ATM%VXWINDTAB(T%ATM%IP+1)-T%ATM%VXWINDTAB(T%ATM%IP))
   T%ATM%VYWIND = T%ATM%VYWINDTAB(T%ATM%IP) + m*(T%ATM%VYWINDTAB(T%ATM%IP+1)-T%ATM%VYWINDTAB(T%ATM%IP))
   T%ATM%VZWIND = T%ATM%VZWINDTAB(T%ATM%IP) + m*(T%ATM%VZWINDTAB(T%ATM%IP+1)-T%ATM%VZWINDTAB(T%ATM%IP))

  end if

  if (T%ATM%MODNO .eq. 4) then
     !write(*,*) 'Calling WRF and Turbulence Model'
     T%ATM%DEN = 0.002363; !Keep density constant at least for now
     call WRFMODEL(T) !Call WRF model. All velocities needed are set inside this subroutine
     !//Add in Full Field Dryden Gust model
     call TURBULENCE(T)
  end if

  !write(*,*) T%ATM%WRFX

  !!Add all winds
  T%ATM%VXWIND = T%ATM%VXWIND + T%ATM%WINDGUST(1) + T%ATM%WRFX
  T%ATM%VYWIND = T%ATM%VYWIND + T%ATM%WINDGUST(2) + T%ATM%WRFY
  T%ATM%VZWIND = T%ATM%VZWIND + T%ATM%WINDGUST(3) + T%ATM%WRFZ

  !Ramp Winds in
  if (T%ATM%MODNO .eq. 4) then
     ramp = 0
     if (T%SIM%TIME .lt. T%ATM%RAMPTIME) then
        ramp = T%SIM%TIME/T%ATM%RAMPTIME
     else
        ramp = 1
     end if
     T%ATM%VXWIND = ramp*T%ATM%VXWIND
     T%ATM%VYWIND = ramp*T%ATM%VYWIND
     T%ATM%VZWIND = ramp*T%ATM%VZWIND
  end if

  !write(*,*) T%ATM%VXWIND

  RETURN

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ECHO DATA iflag = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (iflag .eq. 2) then 
  write(25,*) ' '
  write(25,*) 'AAAA TTTT M   M OOOO SSSS PPPP H  H EEEE RRR  EEEE'
  write(25,*) 'A  A  TT  MM MM O  O SS   PP P H  H EE   RR R EE  '
  write(25,*) 'AAAA  TT  M M M O  O SSSS PPPP HHHH EEEE RRRR EEEE'
  write(25,*) 'A  A  TT  M   M O  O   SS PP   H  H EE   RRR  EE  '
  write(25,*) 'A  A  TT  M   M OOOO SSSS PP   H  H EEEE RR R EEEE'
  write(25,*) ' '
  write(25,*) 'Atmosphere Input File: '
  write(25,*) trim(T%ATMOSPHEREINPUTFILE)
  write(25,*) ' '
  write(25,*) 'Model Number (0=Constant, 1=Equation, 2=Table): ',T%ATM%MODNO
  write(25,*) ' '
  if (T%ATM%MODNO .eq. 1) then
   write(25,*) 'Density (kg/m^3): ',T%ATM%DEN
   write(25,*) 'Wind Speed (m/s): ',T%ATM%WINDSPEED
   write(25,*) 'Wind Direction (deg): ',57.3*T%ATM%WINDDIR
   write(25,*) 'Wind Elevation (deg): ',57.3*T%ATM%WINDELEV
   write(25,*) ' '
  end if
  if (T%ATM%MODNO .eq. 2) then
   write(25,*) 'Wind Speed (m/s): ',T%ATM%WINDSPEED
   write(25,*) 'Wind Direction (deg): ',57.3*T%ATM%WINDDIR
   write(25,*) 'Wind Elevation (deg): ',57.3*T%ATM%WINDELEV
   write(25,*) ' '
  end if
  if (T%ATM%MODNO .eq. 3) then
   write(25,*) 'Altitude (m), Density (kg/m^3), VX Wind Speed (m/s), VY Wind Speed (m/s), VZ Wind Speed (m/s)'
   write(25,*) '---------------------------------------------------------------------------------------------'
   do i=1,T%ATM%TABSIZE  
    write(25,fmt='(5e18.8)') T%ATM%ALTTAB(i),T%ATM%DENTAB(i),T%ATM%VXWINDTAB(i),T%ATM%VYWINDTAB(i),T%ATM%VZWINDTAB(i)
   end do
   write(25,*) ' '
  end if
  if (T%ATM%MODNO .eq. 4) then
     write(25,*) 'Using Wind File = ',T%ATM%PATH 
     write(25,*) 'Grid Size of WRF model = ',T%ATM%dx,' m'
     write(25,*) 'Maximum Height of WRF model = ',T%ATM%ztop,' m'
     write(25,*) 'Grid Size of Turbulence(m) = ', T%ATM%dxT
     write(25,*) 'Wind scale = ',T%ATM%IWINDSCALE
     write(25,*) 'Turbulence Scale = ',T%ATM%TURBLEVEL
     write(25,*) 'Heading Offset in WRF model(rad) = ',T%ATM%PSIOFFSET
     write(25,*) 'Wavespeed in X and Y (m/s) = ',T%ATM%WAVESPEED(1),T%ATM%WAVESPEED(2)
     write(25,*) 'Ramp Time (sec) = ',T%ATM%RAMPTIME
  end if
  write(25,*) 'Data Quality Flag (nd, 0=Data Not Loaded Successfully, 1=Data Loaded Successfully): ',T%ATM%DQFLAG

  RETURN
  
end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!! LOAD DATA iflag = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 if (iflag .eq. 1) then

    write(*,*) 'Loading Atmosphere File'
    
  open(unit=94,file=T%ATMOSPHEREINPUTFILE,status='old',iostat=openflag)
  if (openflag .ne. 0) then
   write(*,*) 'Error Opening Atmosphere Input File: ',T%ATMOSPHEREINPUTFILE;  STOP
  else
     write(*,*) 'Atmosphere File opened'
  end if
  rewind(94)
  
  read(unit=94,fmt=*,iostat=readflag) readreal; T%ATM%MODNO = readreal
  if (T%ATM%MODNO .eq. 1) then
   read(unit=94,fmt=*,iostat=readflag) T%ATM%DEN
   read(unit=94,fmt=*,iostat=readflag) T%ATM%WINDSPEED
   read(unit=94,fmt=*,iostat=readflag) T%ATM%WINDDIR
   read(unit=94,fmt=*,iostat=readflag) T%ATM%WINDELEV
  end if
  if (T%ATM%MODNO .eq. 2) then
   read(unit=94,fmt=*,iostat=readflag) T%ATM%WINDSPEED
   read(unit=94,fmt=*,iostat=readflag) T%ATM%WINDDIR
   read(unit=94,fmt=*,iostat=readflag) T%ATM%WINDELEV
  end if
  if (T%ATM%MODNO .eq. 3) then
   read(unit=94,fmt=*,iostat=readflag) readreal; T%ATM%TABSIZE = readreal
   do i=1,T%ATM%TABSIZE  
    read(unit=94,fmt=*,iostat=readflag) T%ATM%ALTTAB(i),T%ATM%DENTAB(i),T%ATM%VXWINDTAB(i),T%ATM%VYWINDTAB(i),T%ATM%VZWINDTAB(i)
   end do
  end if
  if (T%ATM%MODNO .eq. 4) then
     !Atmospheric Density
     read(unit=94,fmt=*,iostat=readflag) T%ATM%DEN
     read(unit=94,fmt=*,iostat=readflag) T%ATM%IWINDSCALE
     read(unit=94,fmt=*,iostat=readflag) T%ATM%TURBLEVEL
     read(unit=94,fmt=*,iostat=readflag) readreal; T%ATM%INTERPTYPE = int(readreal)
     read(unit=94,fmt=*,iostat=readflag) T%ATM%PSIOFFSET
     read(unit=94,fmt=*,iostat=readflag) T%ATM%RAMPTIME
     read(unit=94,fmt=*,iostat=readflag) T%ATM%WAVESPEED(1)
     read(unit=94,fmt=*,iostat=readflag) T%ATM%WAVESPEED(2)
     !Read PATH
     read(unit=94,fmt=*,iostat=readflag) T%ATM%PATH

     !%%%%%%%%Import Extra Parameters File%%%%%%%%% */

     ParametersFile = 'Parameters.txt';
     inParameters = trim(T%ATM%PATH)//ParametersFile
    open(unit=78,file=inParameters,status='old',iostat=ierr)
     if (ierr .ne. 0) then
        write(*,*) 'Parameters.txt File defined incorrectly. Check to make sure your path is set correctly'
        write(*,*) 'Looking in Path: ',T%ATM%PATH
        STOP;
     endif
     !Zero out parameters
     T%ATM%parametersSEVEN(:) = 0
     T%ATM%parametersFIVE(:) = 0
     !Ok so we can't just read the entire file because some files only have 5 inputs 
     !It's possible that the length of parameters is only 5. 
     !However since we zero out parameters before we read 
     !in the parameters file we can simply check this.
     read(unit=78,fmt=*,iostat=readflag) T%ATM%parametersSEVEN
     close(78)
     if (readflag .eq. -1) then
        write(*,*) 'Error in Parameter File'
        write(*,*) 'Defaulting to 5 Parameter File'
        open(unit=78,file=inParameters,status='old',iostat=ierr)
        read(unit=78,fmt=*) T%ATM%parametersFIVE
        close(78)
        T%ATM%parameters = T%ATM%parametersFIVE
        T%ATM%dimZ = T%ATM%parameters(5);
        T%ATM%dimX = T%ATM%dimZ
        T%ATM%dimY = T%ATM%dimZ
     else
        write(*,*) 'Seven Parameter File Found'
        T%ATM%parameters = T%ATM%parametersSEVEN(1:5)
        T%ATM%dimX = T%ATM%parametersSEVEN(6) 
        T%ATM%dimY = T%ATM%parametersSEVEN(7)
     end if
     T%ATM%dx = T%ATM%parameters(1);
     T%ATM%dy = T%ATM%parameters(2);
     T%ATM%ztop = T%ATM%parameters(3);
     T%ATM%dimZ = T%ATM%parameters(5);
  
     ZcoordFile = 'Zcoord.txt';
     inZcoord = trim(T%ATM%PATH)//ZcoordFile
     open(unit=78,file=inZcoord,status='old',iostat=ierr)
     if (ierr .ne. 0) then
        write(*,*) 'Zcoord.txt File defined incorrectly. Check to make sure your path is set correctly'
        write(*,*) 'Looking in Path: ',T%ATM%PATH
        STOP;
     endif
     read(78,*) T%ATM%zcoord
     close(78)
  
     TimesFile = 'SampleTimes.txt';
     inTimes = trim(T%ATM%PATH)//TimesFile
     open(unit=78,file=inTimes,status='old',iostat=ierr)
     if (ierr .ne. 0) then
        write(*,*) 'SampleTimes.txt File defined incorrectly. Check to make sure your path is set correctly'
        write(*,*) 'Looking in Path: ',T%ATM%PATH
        STOP;
     endif
     read(78,*) T%ATM%tcoord
     close(78)

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

     !%%%%%%%%%%%Create xcoord and ycoord%%%%%%% */
     !I think it's safe to assume that dimX=dimY but just in case
     dim = T%ATM%dimX
     do ii = 1,dim
        T%ATM%xcoord(ii) = -T%ATM%dx*(dim-1)/2 + (ii-1)*T%ATM%dx
     enddo
     dim = T%ATM%dimY
     do ii = 1,dim
        T%ATM%ycoord(ii) = -T%ATM%dy*(dim-1)/2 + (ii-1)*T%ATM%dy
     enddo


     !%%%%%%%Import Initial UVW matrices%%%%%%%%% */
  
     write(number, '(i1)' )  T%ATM%tcoord(1)
     letter = trim('U')
     T%ATM%U0name = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
     letter = trim('V')
     T%ATM%V0name = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
     letter = trim('W')
     T%ATM%W0name = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
     write(number, '(i1)' )  T%ATM%tcoord(2)
     T%ATM%Wdtname = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
     letter = trim('U')
     T%ATM%Udtname = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
     letter = trim('V')
     T%ATM%Vdtname = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'

     call IMPORTWIND(T,T%ATM%U0,T%ATM%U0name);
     call IMPORTWIND(T,T%ATM%Udt,T%ATM%Udtname);
     call IMPORTWIND(T,T%ATM%V0,T%ATM%V0name);
     call IMPORTWIND(T,T%ATM%Vdt,T%ATM%Vdtname);
     call IMPORTWIND(T,T%ATM%W0,T%ATM%W0name);
     call IMPORTWIND(T,T%ATM%Wdt,T%ATM%Wdtname);
  
     !write(*,*) 'UVW Initialization Complete'

     if (T%ATM%TURBLEVEL .gt. 0) then

        write(*,*) 'Dryden Initialization '

        !%%%%%%%%%%%Initialize Variables%%%%%%%%%%%%% */

        !Allocate memory for WRF interpolation
        T%ATM%boundsT = 0;
        T%ATM%boundflagT = 1;

        !%%%%%%%%%%Unwrap Parameters%%%%%%%%%%%%%%% 

        !%%%%%%%%%%%Create xcoord and ycoord%%%%%%%
        do ii = 1,T%ATM%dimT
           T%ATM%xcoordT(ii) = -T%ATM%dxT * (T%ATM%dimT - 1) / 2 + ii * T%ATM%dxT;
           T%ATM%ycoordT(ii) = -T%ATM%dyT * (T%ATM%dimT - 1) / 2 + ii * T%ATM%dyT;
        end do
        !%%%%%%%Import Initial UVW matrices%%%%%%%%%

        T%ATM%Uturbname = trim(T%ATM%PATH)//trim('Uturb.txt')
        T%ATM%Vturbname = trim(T%ATM%PATH)//trim('Vturb.txt')
        T%ATM%Wturbname = trim(T%ATM%PATH)//trim('Wturb.txt')
        
        call IMPORTTURB(T%ATM%UTURB, T%ATM%UTurbname);
        call IMPORTTURB(T%ATM%VTURB, T%ATM%VTurbname);
        call IMPORTTURB(T%ATM%WTURB, T%ATM%WTurbname);
        
        write(*,*) 'Turbulence Initialization Complete'
     end if

     close(94) 
     write(*,*) 'ATMOSPHERE Load Complete'
  
  end if

  T%ATM%DQFLAG = 1
  
  RETURN
end if
 
RETURN
END SUBROUTINE ATMOSPHERE

!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINE AIRCRAFT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AIRCRAFT(T,iflag,aflag)
 use MAAWMDATATYPES
 implicit none
 integer i,jdx,iflag,ifind,openflag,readflag,aflag
 real*8 m,readreal,deti,c1,c2,c3
 real*8 xcg,ycg,zcg,q0,q1,q2,q3,ub,vb,wb,pb,qb,rb
 real*8 xcgdot,ycgdot,zcgdot,q0dot,q1dot,q2dot,q3dot,ubdot,vbdot,wbdot,pbdot,qbdot,rbdot
 real*8 sphi,cphi,stheta,ctheta,spsi,cpsi,tphi(3,3),ttheta,tpsi(3,3),tib(3,3),vATM_A(3,1),vATM_I(3,1)
 real*8 rxlse,rylse,rzlse,xpt,ypt,zpt,upt,vpt,wpt,uptaero,vptaero,wptaero,C_roll,q_inf,salfa,calfa
 real*8 sphilse,cphilse,sgamlse,cgamlse,ulse,vlse,wlse,cllse,cdlse,liftlse,draglse,sangle,cangle
 real*8 fxloc,fyloc,fzloc,ct,J,T_0,T_A,Q_A,V_A,omega,rps,alfa,beta,C_D,C_L,C_Y,C_n,C_m,uaero,vaero,waero
 real*8 phidot,thetadot,psidot,phi,theta,psi,sumomega,thrust,forcevec(4,1),bquad,tauaerovec(3,4),omegar,Gammavec(3,1)
 type(MAAWMSTRUCTURE) T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTE iflag = 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 if (iflag .eq. 3) then  
 
  ! Unwrap State Vector 
  
  xcg = T%AC(aflag)%STATE(1)
  ycg = T%AC(aflag)%STATE(2)
  zcg = T%AC(aflag)%STATE(3)
  phi = T%AC(aflag)%STATE(4) 
  theta = T%AC(aflag)%STATE(5)
  psi = T%AC(aflag)%STATE(6)
  ub = T%AC(aflag)%STATE(7)
  vb = T%AC(aflag)%STATE(8)
  wb = T%AC(aflag)%STATE(9)
  pb = T%AC(aflag)%STATE(10)
  qb = T%AC(aflag)%STATE(11) 
  rb = T%AC(aflag)%STATE(12)

  ! Aircraft to Inertial Transformation Matrix
  ctheta = cos(theta);
  stheta = sin(theta);
  ttheta = stheta / ctheta;
  cphi = cos(phi);
  sphi = sin(phi);
  spsi = sin(psi);
  cpsi = cos(psi);
  T%AC(aflag)%TIA(1,1) = ctheta * cpsi;
  T%AC(aflag)%TIA(2,1) = ctheta * spsi;
  T%AC(aflag)%TIA(3,1) = -stheta;
  T%AC(aflag)%TIA(1,2) = sphi * stheta * cpsi - cphi * spsi;
  T%AC(aflag)%TIA(2,2) = sphi * stheta * spsi + cphi * cpsi;
  T%AC(aflag)%TIA(3,2) = sphi * ctheta;
  T%AC(aflag)%TIA(1,3) = cphi * stheta * cpsi + sphi * spsi;
  T%AC(aflag)%TIA(2,3) = cphi * stheta * spsi - sphi * cpsi;
  T%AC(aflag)%TIA(3,3) = cphi * ctheta;

  ! Inertial to Aircraft Transformation Matrix 
  
  T%AC(aflag)%TAI = transpose(T%AC(aflag)%TIA)
  
  ! Gravity Forces and Moments
  
  T%AC(aflag)%FXGRAV = 0.0; T%AC(aflag)%FYGRAV = 0.0; T%AC(aflag)%FZGRAV = 0.0;
  T%AC(aflag)%MXGRAV = 0.0; T%AC(aflag)%MYGRAV = 0.0; T%AC(aflag)%MZGRAV = 0.0;
  if (T%AC(aflag)%GRAVOFFON .eq. 1) then
   T%AC(aflag)%FXGRAV = T%AC(aflag)%TAI(1,3)*T%AC(aflag)%WEIGHT
   T%AC(aflag)%FYGRAV = T%AC(aflag)%TAI(2,3)*T%AC(aflag)%WEIGHT
   T%AC(aflag)%FZGRAV = T%AC(aflag)%TAI(3,3)*T%AC(aflag)%WEIGHT
  end if
  
  ! Aerodynamic Forces and Moments
  
  T%AC(aflag)%FXAERO = 0.0; T%AC(aflag)%FYAERO = 0.0; T%AC(aflag)%FZAERO = 0.0;
  T%AC(aflag)%MXAERO = 0.0; T%AC(aflag)%MYAERO = 0.0; T%AC(aflag)%MZAERO = 0.0;

  if (T%AC(aflag)%AEROOFFON .eq. 1) then

    !Compute Atmopsheric density and winds

    T%ATM%XI = xcg
    T%ATM%YI = ycg
    T%ATM%ZI = zcg
    call ATMOSPHERE(T,3) !T%ATM%DEN

    vATM_I(1,1) = T%ATM%VXWIND
    vATM_I(2,1) = T%ATM%VYWIND
    vATM_I(3,1) = T%ATM%VZWIND
    vATM_A = matmul(T%AC(aflag)%TAI,vATM_I)

    !Add in atmospheric winds

    uaero = ub - vATM_A(1,1)
    vaero = vb - vATM_A(2,1)
    waero = wb - vATM_A(3,1)

    !Compute total velocity

    V_A = sqrt(uaero**2 + vaero**2 + waero**2)

    !!Dynamic pressure
    q_inf = 0.5*T%ATM%DEN*(V_A**2)*T%AC(aflag)%SAREA
    
    !!Angle of attack and sideslip
    
    if (uaero .gt. 0.0D0) then
       alfa = atan2(waero,uaero)
    else
       alfa = 0.0D0
    end if
    if (V_A .gt. 0.0D0) then
       beta = asin(vaero/V_A)
    else
       beta = 0.0D0
    end if
    calfa = cos(alfa)
    salfa = sin(alfa)

    !!THRUST MODEL
    
    if (T%AC(aflag)%THRUSTMOD .eq. 1) then
       ! Use DELTHRUST to obtain omega
       ifind = 0
       if (T%AC(aflag)%DELTHRUST .le. T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1))) then 
          ifind = -1 
          do while ((ifind.ne.0) .and. (T%AC(aflag)%IPROP(1).gt.1))
             T%AC(aflag)%IPROP(1) = T%AC(aflag)%IPROP(1) - 1 
             if (T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1))   .le. T%AC(aflag)%DELTHRUST) then 
                if (T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1)+1) .gt. T%AC(aflag)%DELTHRUST) then 
                   ifind = 0
                end if
             end if
          end do
       end if
       if (T%AC(aflag)%DELTHRUST .ge. T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1)+1)) then 
          ifind = 1
          do while ((ifind.ne.0) .and. (T%AC(aflag)%IPROP(1).lt.T%AC(aflag)%NPROPTAB-1))
             T%AC(aflag)%IPROP(1) = T%AC(aflag)%IPROP(1) + 1 
             if (T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1))   .le. T%AC(aflag)%DELTHRUST) then 
                if (T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1)+1) .gt. T%AC(aflag)%DELTHRUST) then 
                   ifind = 0 
                end if
             end if
          end do
       end if
       if (ifind .eq. 0) then
          m = (T%AC(aflag)%OMEGA_TAB(T%AC(aflag)%IPROP(1)+1)-T%AC(aflag)%OMEGA_TAB(T%AC(aflag)%IPROP(1)))/(T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1)+1)-T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1)))
       else if (ifind .eq. -1) then
          m = 0.0
       else if (ifind .eq. 1) then
          m = 1.0
       end if

       ! Interpolate

       omega = T%AC(aflag)%OMEGA_TAB(T%AC(aflag)%IPROP(1)) + m*(T%AC(aflag)%DELTHRUST-T%AC(aflag)%DELTHRUST_TAB(T%AC(aflag)%IPROP(1))) !rad/s

       !Compute Advance Ratio

       J = PI*T%AC(aflag)%V_T/(omega*T%AC(aflag)%PD/2)

       !Use Advance Ratio to obtain ct

       ifind = 0
       if (J .le. T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(2))) then 
          ifind = -1 
          do while ((ifind.ne.0) .and. (T%AC(aflag)%IPROP(2).gt.1))
             T%AC(aflag)%IPROP(2) = T%AC(aflag)%IPROP(2) - 1 
             if (T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(2))   .le. J) then 
                if (T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(2)+1) .gt. J) then 
                   ifind = 0
                end if
             end if
          end do
       end if
       if (J .ge. T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(2)+1)) then 
          ifind = 1
          do while ((ifind.ne.0) .and. (T%AC(aflag)%IPROP(2).lt.T%AC(aflag)%NPROPTAB-1))
             T%AC(aflag)%IPROP(2) = T%AC(aflag)%IPROP(2) + 1 
             if (T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(1))   .le. J) then 
                if (T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(1)+1) .gt. J) then 
                   ifind = 0 
                end if
             end if
          end do
       end if
       if (ifind .eq. 0) then
          m = (T%AC(aflag)%C_T_TAB(T%AC(aflag)%IPROP(2)+1)-T%AC(aflag)%C_T_TAB(T%AC(aflag)%IPROP(2)))/(T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(2)+1)-T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(2)))
       else if (ifind .eq. -1) then
          m = 0.0
       else if (ifind .eq. 1) then
          m = 1.0
       end if

       !Now interpolate 
       
       ct = T%AC(aflag)%C_T_TAB(T%AC(aflag)%IPROP(2)) + m*(J-T%AC(aflag)%ADVANCERATIO_TAB(T%AC(aflag)%IPROP(2))) !rad/s
       rps = omega/(2*PI) !!Convert rad/s to rev/s

       !Thrust

       T_A = T%ATM%DEN*(rps**2)*(T%AC(aflag)%PD**4)*ct
       Q_A = T_A*V_A/omega
    else if (T%AC(aflag)%THRUSTMOD .eq. 0) then
       T_0 = 5.0882; !Obtained from Thrust equals drag - Only works for UAV Decathalon
       T_A = T_0 + T%AC(aflag)%C_X_DT*T%AC(aflag)%DELTHRUST
       Q_A = 0
    end if

    !!Lift Drag and Side force

    if (V_A .gt. 0.0D0) then 
       C_L = T%AC(aflag)%C_L_0 + T%AC(aflag)%C_L_ALPHA*alfa + T%AC(aflag)%C_L_U*uaero/V_A + T%AC(aflag)%C_L_Q*qb*T%AC(aflag)%C_BAR/(2*V_A) + T%AC(aflag)%C_L_DE*T%AC(aflag)%ELEVATOR
       C_Y = T%AC(aflag)%C_Y_BETA*beta + T%AC(aflag)%C_Y_P*pb*T%AC(aflag)%B/(2*V_A) + T%AC(aflag)%C_Y_R*rb*T%AC(aflag)%B/(2*V_A) + T%AC(aflag)%C_Y_DR*T%AC(aflag)%RUDDER + T%AC(aflag)%C_Y_DA*T%AC(aflag)%AILERON
       C_D = T%AC(aflag)%C_D_0 + T%AC(aflag)%C_D_ALPHA2*alfa**2 + T%AC(aflag)%C_D_U*uaero/V_A + T%AC(aflag)%C_D_DE*T%AC(aflag)%ELEVATOR

       !!Roll,pitch and yaw coefficients

       C_roll = T%AC(aflag)%C_L_BETA*beta + T%AC(aflag)%C_L_P*pb*T%AC(aflag)%B/(2*V_A) + T%AC(aflag)%C_L_R*rb*T%AC(aflag)%B/(2*V_A) + T%AC(aflag)%C_L_DR*T%AC(aflag)%RUDDER + T%AC(aflag)%C_L_DA*T%AC(aflag)%AILERON
       C_m = T%AC(aflag)%C_M_0 + T%AC(aflag)%C_M_ALPHA*alfa + T%AC(aflag)%C_M_U*uaero/V_A + T%AC(aflag)%C_M_Q*qb*T%AC(aflag)%C_BAR/(2*V_A) + T%AC(aflag)%C_M_DE*T%AC(aflag)%ELEVATOR
       C_n = T%AC(aflag)%C_N_BETA*beta + T%AC(aflag)%C_N_P*pb*T%AC(aflag)%B/(2*V_A) + T%AC(aflag)%C_N_R*rb*T%AC(aflag)%B/(2*V_A) + T%AC(aflag)%C_N_DR*T%AC(aflag)%RUDDER + T%AC(aflag)%C_N_DA*T%AC(aflag)%AILERON

       T%AC(aflag)%FXAERO = -q_inf*(calfa*(C_D) - salfa*C_L) + T_A
       T%AC(aflag)%FYAERO = q_inf*C_Y
       T%AC(aflag)%FZAERO = -q_inf*(salfa*(C_D) + calfa*C_L)
       T%AC(aflag)%MXAERO = q_inf*T%AC(aflag)%B*C_roll
       T%AC(aflag)%MYAERO = q_inf*T%AC(aflag)%C_BAR*C_m
       T%AC(aflag)%MZAERO = q_inf*T%AC(aflag)%B*C_n
    end if
  end if !AEROFORCES

  !At this point we should have F(XYZ)AERO and M(XYZ)AERO populated

  ! Contact Forces and Moments

  T%AC(aflag)%FXCONT = 0.0; T%AC(aflag)%FYCONT = 0.0; T%AC(aflag)%FZCONT = 0.0;
  T%AC(aflag)%MXCONT = 0.0; T%AC(aflag)%MYCONT = 0.0; T%AC(aflag)%MZCONT = 0.0;
  if (T%AC(aflag)%CONTOFFON .eq. 1) then
   T%AC(aflag)%FXCONT = 0.0
   T%AC(aflag)%FYCONT = 0.0
   T%AC(aflag)%FZCONT = 0.0
   T%AC(aflag)%MXCONT = 0.0
   T%AC(aflag)%MYCONT = 0.0
   T%AC(aflag)%MZCONT = 0.0
  end if

  ! Total Forces and Moments
  
  T%AC(aflag)%FXTOTAL = T%AC(aflag)%FXGRAV + T%AC(aflag)%FXAERO + T%AC(aflag)%FXCONT
  T%AC(aflag)%FYTOTAL = T%AC(aflag)%FYGRAV + T%AC(aflag)%FYAERO + T%AC(aflag)%FYCONT 
  T%AC(aflag)%FZTOTAL = T%AC(aflag)%FZGRAV + T%AC(aflag)%FZAERO + T%AC(aflag)%FZCONT
  T%AC(aflag)%MXTOTAL = T%AC(aflag)%MXGRAV + T%AC(aflag)%MXAERO + T%AC(aflag)%MXCONT
  T%AC(aflag)%MYTOTAL = T%AC(aflag)%MYGRAV + T%AC(aflag)%MYAERO + T%AC(aflag)%MYCONT
  T%AC(aflag)%MZTOTAL = T%AC(aflag)%MZGRAV + T%AC(aflag)%MZAERO + T%AC(aflag)%MZCONT
  
  ! State Derivatives
  
  xcgdot = T%AC(aflag)%TIA(1,1)*ub + T%AC(aflag)%TIA(1,2)*vb + T%AC(aflag)%TIA(1,3)*wb
  ycgdot = T%AC(aflag)%TIA(2,1)*ub + T%AC(aflag)%TIA(2,2)*vb + T%AC(aflag)%TIA(2,3)*wb
  zcgdot = T%AC(aflag)%TIA(3,1)*ub + T%AC(aflag)%TIA(3,2)*vb + T%AC(aflag)%TIA(3,3)*wb  

  phidot = pb + sphi * ttheta * qb + cphi * ttheta * rb;
  thetadot = cphi * qb - sphi * rb;
  psidot = (sphi / ctheta) * qb + (cphi / ctheta) * rb;
  ubdot = T%AC(aflag)%FXTOTAL/T%AC(aflag)%MASS + rb*vb - qb*wb
  vbdot = T%AC(aflag)%FYTOTAL/T%AC(aflag)%MASS + pb*wb - rb*ub 
  wbdot = T%AC(aflag)%FZTOTAL/T%AC(aflag)%MASS + qb*ub - pb*vb 
  c1 = T%AC(aflag)%MXTOTAL - pb*(qb*T%AC(aflag)%IXZ-rb*T%AC(aflag)%IXY) - qb*(qb*T%AC(aflag)%IYZ-rb*T%AC(aflag)%IYY) - rb*(qb*T%AC(aflag)%IZZ-rb*T%AC(aflag)%IYZ)
  c2 = T%AC(aflag)%MYTOTAL - pb*(rb*T%AC(aflag)%IXX-pb*T%AC(aflag)%IXZ) - qb*(rb*T%AC(aflag)%IXY-pb*T%AC(aflag)%IYZ) - rb*(rb*T%AC(aflag)%IXZ-pb*T%AC(aflag)%IZZ)
  c3 = T%AC(aflag)%MZTOTAL - pb*(pb*T%AC(aflag)%IXY-qb*T%AC(aflag)%IXX) - qb*(pb*T%AC(aflag)%IYY-qb*T%AC(aflag)%IXY) - rb*(pb*T%AC(aflag)%IYZ-qb*T%AC(aflag)%IXZ)
  pbdot = T%AC(aflag)%IXXI*c1 + T%AC(aflag)%IXYI*c2 + T%AC(aflag)%IXZI*c3
  qbdot = T%AC(aflag)%IXYI*c1 + T%AC(aflag)%IYYI*c2 + T%AC(aflag)%IYZI*c3
  rbdot = T%AC(aflag)%IXZI*c1 + T%AC(aflag)%IYZI*c2 + T%AC(aflag)%IZZI*c3
  
  ! Wrap State Derivatives
  
  T%AC(aflag)%STATEDOT(1) = xcgdot
  T%AC(aflag)%STATEDOT(2) = ycgdot
  T%AC(aflag)%STATEDOT(2) = ycgdot
  T%AC(aflag)%STATEDOT(3) = zcgdot
  T%AC(aflag)%STATEDOT(4) = phidot
  T%AC(aflag)%STATEDOT(5) = thetadot
  T%AC(aflag)%STATEDOT(6) = psidot
  T%AC(aflag)%STATEDOT(7) = ubdot
  T%AC(aflag)%STATEDOT(8) = vbdot
  T%AC(aflag)%STATEDOT(9) = wbdot
  T%AC(aflag)%STATEDOT(10) = pbdot 
  T%AC(aflag)%STATEDOT(11) = qbdot 
  T%AC(aflag)%STATEDOT(12) = rbdot
   
  RETURN
  
 end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ECHO DATA iflag = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 if (iflag .eq. 2) then 
  write(25,*) ' '
  write(25,*) 'AAAA IIII RRR  CCCC RRR  AAAA FFFF TTTT'
  write(25,*) 'A  A  II  RR R CC   RR R A  A FF    TT '
  write(25,*) 'AAAA  II  RRRR CC   RRRR AAAA FFFF  TT '
  write(25,*) 'A  A  II  RRR  CC   RRR  A  A FF    TT '
  write(25,*) 'A  A IIII RR R CCCC RR R A  A FF    TT '
  write(25,*) ' '
  do i = 1,T%SIM%NOAIRPLANES
     write(25,*) 'Airplane Input File: '
     write(25,*) trim(T%AIRPLANEINPUTFILE)
     write(25,*) ' '  
     write(25,*) 'Dynamics Flag (0=Off, 1=On): ',T%AC(i)%DYNOFFON
     write(25,*) 'Gravity Flag (0=Off, 1=On): ',T%AC(i)%GRAVOFFON
     write(25,*) 'Aerodynamics Flag (0=Off, 1=On): ',T%AC(i)%AEROOFFON
     write(25,*) 'Contact Flag (0=Off, 1=On): ',T%AC(i)%CONTOFFON
     write(25,*) 'Mass (kg): ',T%AC(i)%MASS
     write(25,*) 'Weight (N): ',T%AC(i)%WEIGHT
     write(25,*) 'Stationline of Mass Center (m): ',T%AC(i)%SLCG
     write(25,*) 'Buttline of Mass Center (m): ',T%AC(i)%BLCG
     write(25,*) 'Waterline of Mass Center (m): ',T%AC(i)%WLCG
     write(25,*) 'Ixx (kg m^2): ',T%AC(i)%IXX
     write(25,*) 'Iyy (kg m^2): ',T%AC(i)%IYY
     write(25,*) 'Izz (kg m^2): ',T%AC(i)%IZZ
     write(25,*) 'Ixy (kg m^2): ',T%AC(i)%IXY
     write(25,*) 'Ixz (kg m^2): ',T%AC(i)%IXZ
     write(25,*) 'Iyz (kg m^2): ',T%AC(i)%IYZ
     write(25,*) 'Ixx Inverse (1/(kg m^2)): ',T%AC(i)%IXXI
     write(25,*) 'Iyy Inverse (1/(kg m^2)): ',T%AC(i)%IYYI
     write(25,*) 'Izz Inverse (1/(kg m^2)): ',T%AC(i)%IZZI
     write(25,*) 'Ixy Inverse (1/(kg m^2)): ',T%AC(i)%IXYI
     write(25,*) 'Ixz Inverse (1/(kg m^2)): ',T%AC(i)%IXZI
     write(25,*) 'Iyz Inverse (1/(kg m^2)): ',T%AC(i)%IYZI
     write(25,*) 'Turn Radius of AC(m): ',  T%AC(i)%TURNRADIUS
     write(25,*) 'Reference area of aircraft (m^2): ',  T%AC(i)%SAREA
     write(25,*) 'Wingspan of aircraft(m): ',  T%AC(i)%B 
     write(25,*) 'Mean chord of aircraft(m): ',  T%AC(i)%C_BAR
     write(25,*) 'Trim Velocity of Aircraft(m/s): ',  T%AC(i)%V_T
     write(25,*) 'propellor diameter (m): ',  T%AC(i)%PD
     write(25,*) 'zero lift slope: ',  T%AC(i)%C_L_0
     write(25,*) 'lift curve slope: ',  T%AC(i)%C_L_ALPHA
     write(25,*) 'increase in lift with speed: ',  T%AC(i)%C_L_U 
     write(25,*) 'increase in lift with q: ',  T%AC(i)%C_L_Q 
     write(25,*) 'increase in lift with elevator: ',  T%AC(i)%C_L_DE
     write(25,*) 'side force w.r.t sideslip: ',  T%AC(i)%C_Y_BETA
     write(25,*) 'side force w.r.t roll rate: ',  T%AC(i)%C_Y_P 
     write(25,*) 'side force w.r.t. yaw rate: ',  T%AC(i)%C_Y_R 
     write(25,*) 'side force w.r.t. rudder: ',  T%AC(i)%C_Y_DR
     write(25,*) 'side force w.r.t. aileron: ',  T%AC(i)%C_Y_DA
     write(25,*) 'zero lift drag: ',  T%AC(i)%C_D_0 
     write(25,*) 'drag polar: ',  T%AC(i)%C_D_ALPHA2
     write(25,*) 'increase in drag with speed: ',  T%AC(i)%C_D_U 
     write(25,*) 'inrease in drag with elevator: ',  T%AC(i)%C_D_DE
     write(25,*) 'roll moment w.r.t beta: ',  T%AC(i)%C_L_BETA
     write(25,*) 'roll moment w.r.t roll rate: ',  T%AC(i)%C_L_P 
     write(25,*) 'roll moment w.r.t. yaw rate: ',  T%AC(i)%C_L_R 
     write(25,*) 'roll moment w.r.t. rudder: ',  T%AC(i)%C_L_DR
     write(25,*) 'roll moment w.r.t. aileron: ',  T%AC(i)%C_L_DA
     write(25,*) 'zero lift moment: ',  T%AC(i)%C_M_0 
     write(25,*) 'pitch moment curve: ',  T%AC(i)%C_M_ALPHA
     write(25,*) 'increase in pitch moment w.r.t speed: ',  T%AC(i)%C_M_U
     write(25,*) 'pitch moment w.r.t. q: ',  T%AC(i)%C_M_Q
     write(25,*) 'pitch moment w.r.t elevator: ',  T%AC(i)%C_M_DE
     write(25,*) 'yaw moment w.r.t. beta: ',  T%AC(i)%C_N_BETA
     write(25,*) 'yaw moment w.r.t. roll rate: ',  T%AC(i)%C_N_P 
     write(25,*) 'yaw moment w.r.t. yaw rate: ',  T%AC(i)%C_N_R 
     write(25,*) 'yaw moment w.r.t. rudder: ',  T%AC(i)%C_N_DR
     write(25,*) 'yaw moment w.r.t. aileron: ',  T%AC(i)%C_N_DA
     write(25,*) 'Thrust model 0 = linear , 1 = table look up: ', T%AC(i)%THRUSTMOD
     if (T%AC(i)%THRUSTMOD .eq. 0) then
        write(25,*) 'Thrust w.r.t. delthrust', T%AC(i)%C_X_DT
     else if (T%AC(i)%THRUSTMOD .eq. 1 ) then
        write(25,*) 'Propellor table size: ',T%AC(i)%NPROPTAB
        write(25,*) 'Delthrust         Omega(rad/s)       Advance Ratio      C_T'
        write(25,*) '-----------------------------------------------------------'
        do jdx=1,T%AC(i)%NPROPTAB
           write(25,fmt='(4e18.8)') T%AC(i)%DELTHRUST_TAB(jdx),T%AC(i)%OMEGA_TAB(jdx),T%AC(i)%ADVANCERATIO_TAB(jdx),T%AC(i)%C_T_TAB(jdx)
        end do
     end if
     write(25,*) ' '
  end do
  write(25,*) 'Data Quality Flag (nd, 0=Data Not Loaded Successfully, 1=Data Loaded Successfully): ',T%AC(1)%DQFLAG

  RETURN

 end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!! LOAD DATA iflag = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 if (iflag .eq. 1) then

  !!!!!LOAD THE AIRPLANE PARAMETERS
 
  open(unit=94,file=T%AIRPLANEINPUTFILE,status='old',iostat=openflag)
  if (openflag .ne. 0) then
   write(*,*) 'Error Opening Aircraft Input File: ',T%AIRPLANEINPUTFILE;  STOP
  end if
  rewind(94)
  
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%DYNOFFON
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%GRAVOFFON
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%AEROOFFON
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%CONTOFFON
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%MASS
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%SLCG
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%BLCG
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%WLCG
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%IXX
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%IYY
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%IZZ
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%IXY
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%IXZ
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%IYZ
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%TURNRADIUS
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%SAREA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%B 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_BAR
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%V_T
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%PD
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_0
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_ALPHA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_U 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_Q 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_DE
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_Y_BETA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_Y_P 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_Y_R 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_Y_DR
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_Y_DA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_D_0 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_D_ALPHA2
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_D_U 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_D_DE
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_BETA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_P 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_R 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_DR
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_L_DA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_M_0 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_M_ALPHA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_M_U
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_M_Q
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_M_DE
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_N_BETA
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_N_P 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_N_R 
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_N_DR
  read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_N_DA
  read(unit=94,fmt=*,iostat=readflag) readreal; T%AC(1)%THRUSTMOD = readreal
  if (T%AC(1)%THRUSTMOD .eq. 0) then !Linear
     read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_X_DT
  else if (T%AC(1)%THRUSTMOD .eq. 1) then !Table look up
     read(unit=94,fmt=*,iostat=readflag) readreal; T%AC(1)%NPROPTAB = readreal
     do i = 1,T%AC(1)%NPROPTAB
        read(unit=94,fmt=*,iostat=readflag) T%AC(1)%C_T_TAB(i)
     end do
     do i = 1,T%AC(1)%NPROPTAB
        read(unit=94,fmt=*,iostat=readflag) T%AC(1)%ADVANCERATIO_TAB(i)
     end do
     do i = 1,T%AC(1)%NPROPTAB
        read(unit=94,fmt=*,iostat=readflag) T%AC(1)%OMEGA_TAB(i)
     end do
     do i = 1,T%AC(1)%NPROPTAB
        read(unit=94,fmt=*,iostat=readflag) T%AC(1)%DELTHRUST_TAB(i)
     end do
  end if

  close(94) 
  write(*,*) 'AIRPLANE Load Complete'

  !!!DO SOME CALCULATIONS ON AIRPLANE
  T%AC(1)%WEIGHT = T%GRAVITY*T%AC(1)%MASS
  deti = + T%AC(1)%IXX*(T%AC(1)%IYY*T%AC(1)%IZZ-T%AC(1)%IYZ*T%AC(1)%IYZ) - T%AC(1)%IXY*(T%AC(1)%IXY*T%AC(1)%IZZ-T%AC(1)%IYZ*T%AC(1)%IXZ) + T%AC(1)%IXZ*(T%AC(1)%IXY*T%AC(1)%IYZ-T%AC(1)%IYY*T%AC(1)%IXZ)
  T%AC(1)%IXXI = (T%AC(1)%IYY*T%AC(1)%IZZ-T%AC(1)%IYZ*T%AC(1)%IYZ)/deti
  T%AC(1)%IXYI = (T%AC(1)%IYZ*T%AC(1)%IXZ-T%AC(1)%IXY*T%AC(1)%IZZ)/deti
  T%AC(1)%IXZI = (T%AC(1)%IXY*T%AC(1)%IYZ-T%AC(1)%IYY*T%AC(1)%IXZ)/deti
  T%AC(1)%IYYI = (T%AC(1)%IXX*T%AC(1)%IZZ-T%AC(1)%IXZ*T%AC(1)%IXZ)/deti
  T%AC(1)%IYZI = (T%AC(1)%IXY*T%AC(1)%IXZ-T%AC(1)%IXX*T%AC(1)%IYZ)/deti
  T%AC(1)%IZZI = (T%AC(1)%IXX*T%AC(1)%IYY-T%AC(1)%IXY*T%AC(1)%IXY)/deti

  !!Minimum Waypoint radius for aircraft - Crap. These are fixed and would probably change depending on the Aircraft. ugh.
  T%AC(1)%MINWAYPOINT = 50.0D0 
  !T%AC(2)%MINWAYPOINT = 10.0D0 !10 feet seems reasonable for a quad
    
  T%AC(1)%DQFLAG = 1

  RETURN
 
 end if
 
 RETURN
END SUBROUTINE AIRCRAFT

SUBROUTINE LUDCMP(T,a,n,indx,d)
 use MAAWMDATATYPES
 integer i,idxmax,j,k,n,indx(n)
 real*8 aamax,dum,sum,small,d,a(n,n),vv(n)
 type(MAAWMSTRUCTURE) T
 
 small = 1.0e-20
 
 d = 1.0
 do i=1,n
  aamax = 0.0
  do j=1,n
   if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
  end do
  if (aamax .eq. 0.0) then
     !close(90)
   !pause
   write(*,*) "singular matrix in ludcmp"
   !Turn off all simulations
   T%AC(1:T%SIM%NOAIRCRAFT)%DYNOFFON = 0
   STOP;
   !stop
  end if
  vv(i) = 1.0/aamax
 end do
 
 do j=1,n
  do i=1,j-1
   sum = a(i,j)
   do k=1,i-1
    sum = sum - a(i,k)*a(k,j)
   end do
   a(i,j) = sum
  end do
  aamax = 0.0
  do i=j,n
   sum = a(i,j)
   do k=1,j-1
    sum = sum - a(i,k)*a(k,j)
   end do
   a(i,j) = sum
   dum = vv(i)*abs(sum)
   if (dum .ge. aamax) then
    idxmax = i
    aamax = dum
   end if
  end do
  if (j .ne. idxmax)then
   do k=1,n
    dum = a(idxmax,k)
    a(idxmax,k) = a(j,k)
    a(j,k) = dum
   end do
   d = - d
   vv(idxmax) = vv(j)
  end if
  indx(j) = idxmax
  if(a(j,j) .eq. 0.0) then
   a(j,j) = small
  end if
  if (j .ne. n) then
   dum = 1.0/a(j,j)
   do i=j+1,n
    a(i,j) = a(i,j)*dum
   end do
  end if
 end do
      
 RETURN
END SUBROUTINE LUDCMP

SUBROUTINE LUBKSB(T,a,n,indx,b)
 use MAAWMDATATYPES
 integer i,ii,j,ll,n,indx(n)
 real*8 sum,a(n,n),b(n)
 type(MAAWMSTRUCTURE) T
 
 ii = 0
 do i=1,n
  ll = indx(i)
  sum = b(ll)
  b(ll) = b(i)
  if (ii.ne.0) then
   do j=ii,i-1
    sum = sum - a(i,j)*b(j)
   end do
  else if (sum .ne. 0.0) then
   ii = i
  end if
  b(i) = sum
 end do
 do i=n,1,-1
  sum = b(i)
  do j=i+1,n
   sum = sum - a(i,j)*b(j)
  end do
  if (a(i,i) .ne. 0) then
      b(i) = sum/a(i,i)
  else
      b(i) = 0
  end if
 end do
     
 RETURN
END SUBROUTINE LUBKSB

SUBROUTINE LINEARIZEAC(T)
  use MAAWMDATATYPES
  implicit none
  integer stateindex,i,ierr
  real*8 C_Ftether_P(3,1),C_Ftether_I(3,1),A(12,12),statedot0(12,1),del
  type(MAAWMSTRUCTURE) T
  
  !Run Test equilibrium to get theta steady state I guess I could code it here 
  !but we're going for speed at the moment

  !Make sure Atm is off
  T%ATM%MODNO = 1
  T%ATM%WINDSPEED = 0

  !Set contact values of aircraft to zero
  T%AC(1:T%SIM%NOAIRCRAFT)%FXCONT = 0
  T%AC(1:T%SIM%NOAIRCRAFT)%FYCONT = 0
  T%AC(1:T%SIM%NOAIRCRAFT)%FZCONT = 0

  !Echo states
  !T%AC%STATE(1:12) = T%SIM%STATE(1:12)
  write(*,*) 'Aircraft States'
  do i=1,12
     write(*,fmt='(a1,i4,e18.8)') ' ',i,T%AC(1)%STATE(i)
  end do
  
  !Aircraft Derivatives

  T%AC(1)%STATEDOT = 0.0
  call AIRCRAFT(T,3,1)

  !Echo derivatives states
  write(*,*) 'Aircraft Derivatives'
  do i=1,12
     write(*,fmt='(a1,i4,e18.8)') ' ',i,T%AC(1)%STATEDOT(i)
  end do

  !At this point the derivatives should be zero
  !now we can do our finite difference calls and make our A matrix
  statedot0(1:12,1) = T%AC(1)%STATEDOT(1:12)
  del = 1e-8
  do i = 1,12
     !perturb state
     T%AC(1)%STATE(i) = T%AC(1)%STATE(i) + del
     call AIRCRAFT(T,3,1)
     A(1:12,i) = (T%AC(1)%STATEDOT(1:12)-statedot0(1:12,1))/del
     !return state
     T%AC(1)%STATE(i) = T%AC(1)%STATE(i) - del
  end do

  !Print A matrix
  open(unit=99,file='Output_Files/A.txt',iostat=ierr)
  if (ierr .ne. 0) then
     open(unit=99,file='A.txt',iostat=ierr)
  end if
  !write(99,*) 'A matrix'
  do i = 1,12
     write(99,*) A(i,1:12)
  end do

END SUBROUTINE LINEARIZEAC

SUBROUTINE RandGaussian(xg)
  implicit none
  real*8 u1,u2,r,fac
  real*8 xu1,xu2,xg
  
  !c.... Convert Uniform Random Numbers to Gaussian Random Numbers

  r = 2.00000000
  do while (r .ge. 1.00000000)

     !c..... Obtain Uniform Random Numbers Between Zero and One

     call RandUniform(xu1)
     call RandUniform(xu2)

     !c..... Generate Uniform Random Numbers Between Minus One and Plus One
     
     u1 = 2.00000000*xu1 - 1.00000000
     u2 = 2.00000000*xu2 - 1.00000000

     !.. Check if Numbers are in Unit Circle

     r = u1**2 + u2**2
  enddo

  !. Convert Uniform Random Numbers to Gaussian Random Numbers

  fac = sqrt(-2.00000000*log(r)/r)
  xg = u1*fac

  RETURN
END SUBROUTINE

!. SUBROUTINE RandUniform

SUBROUTINE RandUniform(xu)
  implicit none
  integer i,j,iff,m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3,ix1,ix2,ix3
  real*8 rm1,rm2,xu
  real*8 r(97)
  save

  !. Basic Algorithm Data

  parameter(m1=259200,ia1=7141,ic1=54773,rm1=1./m1)
  parameter(m2=134456,ia2=8121,ic2=28411,rm2=1./m2)
  parameter(m3=243000,ia3=4561,ic3=51349)
  data iff /0/

  !. Initialization

  if (iff .eq. 0) then
     iff = 1
     ix1 = mod(ic1-1,m1)
     ix1 = mod(ia1*ix1+ic1,m1)
     ix2 = mod(ix1,m2)
     ix1 = mod(ia1*ix1+ic1,m1)
     ix3 = mod(ix1,m3)
     do j=1,97
        ix1 = mod(ia1*ix1+ic1,m1)
        ix2 = mod(ia2*ix2+ic2,m2)
        r(j) = (float(ix1)+float(ix2)*rm2)*rm1
     enddo
  endif

  !. Generate Uniform Random Number

  do i=1,2
     ix1 = mod(ia1*ix1+ic1,m1)
     ix2 = mod(ia2*ix2+ic2,m2)
     ix3 = mod(ia3*ix3+ic3,m3)
     j = 1 + (97*ix3)/m3
     if ((j.gt.97) .or. (j.lt.1)) then
        j = 43
        write(*,*) ' '
        write(*,*) 'ERROR IN RANDUNIFORM'
        write(*,*) ' '
     endif
     xu = r(j)
     r(j) = (float(ix1)+float(ix2)*rm2)*rm1
  enddo

  RETURN
END SUBROUTINE RandUniform

!!Import routine to load data files placed in text files
SUBROUTINE IMPORTWIND(T,mat,filename)
  use MAAWMDATATYPES
  implicit none
  integer ii,jj,kk,nii,njj,ierr;
  real*8 mat(MAXXYDIM,MAXXYDIM,40),tempmat40(40),tempmat100(100)
  real*8 tempmat(MAXXYDIM)
  character*256 filename
  type(MAAWMSTRUCTURE) T

  write(*,*) 'Importing: ',filename
  
  open(unit=78,file=filename,status='old',iostat=ierr)
  if (ierr .ne. 0) then 
     write(*,*) 'Wind File defined incorrectly'
     write(*,*) filename
     STOP;
  endif
  do jj=1,T%ATM%dimZ
     do ii=1,T%ATM%dimY
        if (T%ATM%dimX .eq. 40) then
           read(78,*) tempmat40
           tempmat(1:40) = tempmat40
        end if
        if (T%ATM%dimX .eq. 100) then
           read(78,*) tempmat100
           tempmat = tempmat100
        end if
        do kk=1,T%ATM%dimX
           mat(ii,kk,jj) = tempmat(kk)
        enddo
     enddo
  enddo
  close(78)

END SUBROUTINE IMPORTWIND

SUBROUTINE IMPORTTURB(outmat,filename) 
  integer ii, jj, kk, nii, njj;
  character*256 filename
  real*8 outmat(500,500),temp(500)

  open(unit=78,file=filename,status='old',iostat=ierr)
  if (ierr .ne. 0) then
     write(*,*) 'Turbulence Data File defined incorrectly'
     write(*,*) filename
     STOP;
  end if
  do ii = 1,500
     read(78,*) temp
     do jj = 1,500
        outmat(ii,jj) = temp(jj)
     end do
  end do
  close(78)

END SUBROUTINE IMPORTTURB

!!!FIND FUNCTIONS FOR ATMOSPHERE MODEL

SUBROUTINE FINDGE(invec,row,val,counter)
  implicit none
  integer idx,counter,row;
  real*8 invec(400),val
  idx = 1;
  counter = row;
  do while (idx .lt. row)
     if ((invec(idx) .ge. val)) then
        counter = idx-1;
        idx = row + 1;
     endif
     idx = idx + 1
  enddo
END SUBROUTINE FINDGE

SUBROUTINE FIND(invec,row,val,counter)
  implicit none
  integer idx,counter,row;
  real*8 invec(400),val
  idx = 1;
  counter = row;
  do while (idx .lt. row)
     if ((invec(idx) .gt. val)) then
        counter = idx-1;
        idx = row + 1;
     endif
     idx = idx + 1
  enddo
END SUBROUTINE FIND

SUBROUTINE FIND2(invec,row,val,counter)
  implicit none
  integer idx,counter,row,invec(601)
  real*8 val
  idx = 1;
  counter = row;
  do while (idx .lt. row)
     if ((invec(idx) .gt. val)) then
        counter = idx-1;
        idx = row + 1;
     endif
     idx = idx + 1
  enddo
END SUBROUTINE FIND2

SUBROUTINE WRFMODEL(T)
  use MAAWMDATATYPES
  implicit none
  integer stepX,stepY,stepZ,stepT,extrapX,extrapY,extrapZ,extrapT,cord2(2)
  integer markX,markY,markZ,markT,gust,body
  integer tinterp,x1,x2,y1,y2,z1,z2,tt,ii,coord1(4),coord2(4),cord1(2)
  real*8 uvw(3,2),xpts2(2),ypts2(2),zpts2(2),zpts1,xpts1,ypts1,rx;
  real*8 u8(8),v8(8),w8(8),u4(4),v4(4),w4(4),uslope,vslope,wslope;
  real*8 u2(2),v2(2),w2(2),u,v,w,tpts(2),Lu,Lv,Lw,sigw,sigu,sigv,tstar;
  real*8 ugo,vgo,wgo,vatm(3),xstar,ystar,zstar,xtemp,vtemp,counter,xwidth,ywidth
  real*8 xi,yi,zi
  character*1 letter
  character*10 number
  type(MAAWMSTRUCTURE) T

  !   /*   %%This function will take in x,y,z(m),t(sec) and location and  */
  !   /*   %return u,v,w(m/s). This uses a fast quad-linear interpolation */
  !   /*   %so many globals must be defined. location is a string that */
  !   /*   %contains the location of the data to be interpolated. */

  xi = T%ATM%XI
  yi = T%ATM%YI
  zi = T%ATM%ZI

  xtemp = xi
  xi = xtemp*cos(T%ATM%PSIOFFSET) + yi*sin(T%ATM%PSIOFFSET);
  yi = -xtemp*sin(T%ATM%PSIOFFSET) + yi*cos(T%ATM%PSIOFFSET);

  !!Extra term is from phase offset to add time-varying component
  xi = xi - T%ATM%WAVESPEED(1)*(T%SIM%TIME-T%SIM%INITIALTIME)
  yi = yi - T%ATM%WAVESPEED(2)*(T%SIM%TIME-T%SIM%INITIALTIME)

  !Now shift the WRF Grid so that XI,YI and ZI fall in side the cube
  T%ATM%xshift = 0
  T%ATM%yshift = 0

  xi = xi*FT2M
  yi = yi*FT2M
  zi = zi*FT2M

  xwidth = T%ATM%xcoord(T%ATM%dimX)-T%ATM%xcoord(1)
  ywidth = T%ATM%ycoord(T%ATM%dimY)-T%ATM%ycoord(1)

  !%Find markX
  if (T%ATM%INTERPTYPE .eq. 0) then
     if (xi .gt. T%ATM%xcoord(T%ATM%dimX)) then
        T%ATM%xshift = -floor(abs(xi-T%ATM%xcoord(1))/xwidth)
     else if (xi .lt. T%ATM%xcoord(1)) then 
        T%ATM%xshift = floor(abs(xi-T%ATM%xcoord(T%ATM%dimX))/xwidth)
     end if

     !%Find markY
     if (yi .gt. T%ATM%ycoord(T%ATM%dimY)) then
        T%ATM%yshift = -floor(abs(yi-T%ATM%ycoord(1))/ywidth)
     else if (yi .lt. T%ATM%ycoord(1)) then 
        T%ATM%yshift = floor(abs(yi-T%ATM%ycoord(T%ATM%dimY))/ywidth)
     end if
  end if

  xstar = xi + T%ATM%xshift*xwidth
  ystar = yi + T%ATM%yshift*ywidth
  zstar = -zi

  ! write(*,*) xi,yi,zi

  ! if (abs(xstar) .gt. T%ATM%xcoord(T%ATM%dimX)) then
  !    write(*,*) 'x = ',xi,xstar
  !    PAUSE
  ! end if
  ! if (abs(ystar) .gt. T%ATM%ycoord(T%ATM%dimY)) then
  !    write(*,*) 'y = ',yi,ystar
  !    PAUSE
  ! end if

  !write(*,*) xi,xstar

  if (T%ATM%IWINDSCALE .gt. 0) then
     if (T%ATM%TIMEVARYING .eq. 0) then
        !write(*,*) 'Constant Time'
        tstar = 0
     else
        tstar = T%SIM%TIME
     end if
     stepX = 1;stepY = 1;stepZ = 1;stepT = 1;
     extrapX = 0;extrapY = 0;extrapZ = 0;extrapT = 0;
     if (zstar .lt. 0) then
        zstar = -zstar
     endif
     tinterp = 2;

     uvw(1,1)=0;uvw(2,1)=0;uvw(3,1)=0;
     uvw(1,2)=0;uvw(2,2)=0;uvw(3,2)=0;

     markX = T%ATM%markX
     markY = T%ATM%markY
     markZ = T%ATM%markZ
     markT = T%ATM%markT

     ! write(*,*) 'Current Markers',markX,markY,markZ
     ! write(*,*) 'xcoord = ',T%ATM%xcoord(1),T%ATM%xcoord(T%ATM%dimX),xstar
     !%%Check X
     if (markX .eq. T%ATM%dimX) then
        markX = markX - 1;
     end if
     if ((xstar .ge. T%ATM%xcoord(markX)) .and. (xstar .le. T%ATM%xcoord(markX+1))) then
        !%%Your in between the markers so keep going
     else
        if (xstar .gt. T%ATM%xcoord(T%ATM%dimX)) then
           !use endpt
           ! write(*,*) 'Using Endpt'
           markX = T%ATM%dimX
           stepX = -1
           extrapX = 1
        else if (xstar .lt. T%ATM%xcoord(1)) then
           ! write(*,*) 'Out of bounds'
           markX = 1
           stepX = 1
           extrapX = 1
        else
           ! write(*,*) 'Searching'
           call FIND(T%ATM%xcoord,T%ATM%dimX,xstar,markX)
           if (markX .eq. T%ATM%dimX) then
              markX = markX - 1;
           else if (markX .le. 0) then
              markX = 1;
           end if
        end if
     end if
     !%%Check Y
     if (markY .eq. T%ATM%dimY) then
        markY = markY - 1;
     end if
     if ((ystar .ge. T%ATM%ycoord(markY)) .and. (ystar .le. T%ATM%ycoord(markY+1))) then
        !%%Your in between the markers so keep going
     else
        if (ystar .gt. T%ATM%ycoord(T%ATM%dimY)) then
           !use endput
           markY = T%ATM%dimY
           stepY = -1
           extrapY = 1
        else if (ystar .lt. T%ATM%ycoord(1)) then
           !use start pt
           markY = 1
           stepY = 1
           extrapY = 1
        else
           call FIND(T%ATM%ycoord,T%ATM%dimY,ystar,markY)
           if (markY .eq. T%ATM%dimY) then
              markY = markY - 1;
           else if (markY .le. 0) then
              markY = 1;
           end if
        end if
     end if
     !%%Check Z
     if (markZ .eq. T%ATM%dimZ) then
        markZ = markZ - 1;
     end if

     if ((zstar .ge. T%ATM%zcoord(markZ)) .and. (zstar .le. T%ATM%zcoord(markZ+1))) then
        !%%Your in between the markers so keep going
     else
        !%Find markZ
        if (zstar .gt. T%ATM%zcoord(T%ATM%dimZ)) then
           !%use endpt
           markZ = T%ATM%dimZ;
           stepZ = -1;
           extrapZ = 1;
        else if (zstar .lt. T%ATM%zcoord(1)) then
           markZ = 1;
           stepZ = 1;
           extrapZ = 1;
        else
           call FIND(T%ATM%zcoord,T%ATM%dimZ,zstar,markZ)
           if (markZ .eq. T%ATM%dimZ) then
              markZ = markZ - 1;
           else if (markZ .eq. 0) then
              markZ = 1;
           end if
        end if
     end if
     !%%Check T
     if (markT .eq. T%ATM%tlength) then
        markT = markT - 1;
     end if
     if ((tstar .ge. T%ATM%tcoord(markT)) .and. (tstar .le. T%ATM%tcoord(markT+1))) then
        !%%Your in between the markers so keep going
     else
        if (T%ATM%TIMEVARYING .eq. 1) then
           !%Find markT
           if (tstar .gt. T%ATM%tcoord(T%ATM%tlength)) then
              !%use endpt
              markT = T%ATM%tlength;
              extrapT = 1;
           else if (tstar .lt. T%ATM%tcoord(1)) then
              !%use start pt
              markT = 1;
              extrapT = 1;
           else
              call FIND2(T%ATM%tcoord,T%ATM%tlength,tstar,markT)
              if (markT .eq. T%ATM%tlength) then
                 markT = markT - 1;
              else if (markT .eq. 0) then
                 markT = 1;
              end if
           end if

           !%%Import U,V,W maps since markT changed
           if (T%ATM%tcoord(markT) .lt. 10) then
              write(number, '(i1)' )  T%ATM%tcoord(markT)
           else
              if (T%ATM%tcoord(markT) .le. 99) then
                 write(number, '(i2)' )  T%ATM%tcoord(markT)
              else
                 write(number, '(i3)' )  T%ATM%tcoord(markT)
              endif
           endif
           letter = trim('U')
           T%ATM%U0name = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
           letter = trim('V')
           T%ATM%V0name = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
           letter = trim('W')
           T%ATM%W0name = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
           !%%only import at markT
           call IMPORTWIND(T,T%ATM%U0,T%ATM%U0name);
           call IMPORTWIND(T,T%ATM%V0,T%ATM%V0name);
           call IMPORTWIND(T,T%ATM%W0,T%ATM%W0name);
           !%U0 = U0(end:-1:1,end:-1:1,:);
           !%V0 = V0(end:-1:1,end:-1:1,:);
           !%W0 = W0(end:-1:1,end:-1:1,:);
           if (extrapT .eq. 1) then
              tinterp = 1;
           else
              !%%import markT + 1
              if (T%ATM%tcoord(markT+1) .lt. 10) then
                 write(number, '(i1)' )  T%ATM%tcoord(markT+1)
              else if (T%ATM%tcoord(markT+1) .le. 99) then
                 write(number, '(i2)' )  T%ATM%tcoord(markT+1)
              else
                 write(number, '(i3)' )  T%ATM%tcoord(markT+1)
              endif
              letter = trim('U')
              T%ATM%Udtname = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
              letter = trim('V')
              T%ATM%Vdtname = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
              letter = trim('W')
              T%ATM%Wdtname = trim(T%ATM%PATH)//trim(letter)//trim(number)//'.txt'
              call IMPORTWIND(T,T%ATM%Udt,T%ATM%Udtname);
              call IMPORTWIND(T,T%ATM%Vdt,T%ATM%Vdtname);
              call IMPORTWIND(T,T%ATM%Wdt,T%ATM%Wdtname);
           end if !(extrapT eq. 1 ) 
        else !Not time varyine
           !write(*,*) 'Not Time Varying'
           tinterp = 1
        end if
     end if

     !%%Interpolation Scheme
     !write(*,*) 'tinterp = ',tinterp
     ! write(*,*) 'MarkX = ',markX,stepX,extrapX
     do tt = 1,tinterp 
        !%Interpolate Spatially

        !%%To start we have 8 discrete point (8 corners of a cube)
        xpts2(1) = T%ATM%xcoord(markX)
        xpts2(2) = T%ATM%xcoord(markX+stepX);
        ypts2(1) = T%ATM%ycoord(markY)
        ypts2(2) = T%ATM%ycoord(markY+stepY);
        zpts2(1) = T%ATM%zcoord(markZ)
        zpts2(2) = T%ATM%zcoord(markZ+stepZ);
        ! write(*,*) xpts2,ypts2,zpts2
        x1 = markX;x2 = markX+stepX;
        y1 = markY;y2 = (markY+stepY);
        z1 = markZ;z2 = markZ+stepZ;
        if (tt .eq. 1) then
           !%%Use U0,V0,W0
           u8(1) = T%ATM%U0(y1,x1,z1);
           u8(2) = T%ATM%U0(y1,x2,z1);
           u8(3) = T%ATM%U0(y2,x2,z1);
           u8(4) = T%ATM%U0(y2,x1,z1);
           u8(5) = T%ATM%U0(y1,x1,z2);
           u8(6) = T%ATM%U0(y1,x2,z2);
           u8(7) = T%ATM%U0(y2,x2,z2);
           u8(8) = T%ATM%U0(y2,x1,z2);
           v8(1) = T%ATM%V0(y1,x1,z1);
           v8(2) = T%ATM%V0(y1,x2,z1);
           v8(3) = T%ATM%V0(y2,x2,z1);
           v8(4) = T%ATM%V0(y2,x1,z1);
           v8(5) = T%ATM%V0(y1,x1,z2);
           v8(6) = T%ATM%V0(y1,x2,z2);
           v8(7) = T%ATM%V0(y2,x2,z2);
           v8(8) = T%ATM%V0(y2,x1,z2);
           w8(1) = T%ATM%W0(y1,x1,z1);
           w8(2) = T%ATM%W0(y1,x2,z1);
           w8(3) = T%ATM%W0(y2,x2,z1);
           w8(4) = T%ATM%W0(y2,x1,z1);
           w8(5) = T%ATM%W0(y1,x1,z2);
           w8(6) = T%ATM%W0(y1,x2,z2);
           w8(7) = T%ATM%W0(y2,x2,z2);
           w8(8) = T%ATM%W0(y2,x1,z2);
        else
           !%%Use Udt,Vdt,Wdt
           u8(1) = T%ATM%Udt(y1,x1,z1);
           u8(2) = T%ATM%Udt(y1,x2,z1);
           u8(3) = T%ATM%Udt(y2,x2,z1);
           u8(4) = T%ATM%Udt(y2,x1,z1);
           u8(5) = T%ATM%Udt(y1,x1,z2);
           u8(6) = T%ATM%Udt(y1,x2,z2);
           u8(7) = T%ATM%Udt(y2,x2,z2);
           u8(8) = T%ATM%Udt(y2,x1,z2);
           v8(1) = T%ATM%Vdt(y1,x1,z1);
           v8(2) = T%ATM%Vdt(y1,x2,z1);
           v8(3) = T%ATM%Vdt(y2,x2,z1);
           v8(4) = T%ATM%Vdt(y2,x1,z1);
           v8(5) = T%ATM%Vdt(y1,x1,z2);
           v8(6) = T%ATM%Vdt(y1,x2,z2);
           v8(7) = T%ATM%Vdt(y2,x2,z2);
           v8(8) = T%ATM%Vdt(y2,x1,z2);
           w8(1) = T%ATM%Wdt(y1,x1,z1);
           w8(2) = T%ATM%Wdt(y1,x2,z1);
           w8(3) = T%ATM%Wdt(y2,x2,z1);
           w8(4) = T%ATM%Wdt(y2,x1,z1);
           w8(5) = T%ATM%Wdt(y1,x1,z2);
           w8(6) = T%ATM%Wdt(y1,x2,z2);
           w8(7) = T%ATM%Wdt(y2,x2,z2);
           w8(8) = T%ATM%Wdt(y2,x1,z2);
        end if

        ! write(*,*) u8
        ! write(*,*) v8
        ! write(*,*) w8


        !%%%%%interpZ%%%%%%%%%%%%

        if (extrapZ .eq. 1) then
           !%%You don't need to interpolate on z and you can just use
           !%%the values at markZ or z1
           zpts1 = zpts2(1);
           u4(1) = u8(1);
           u4(2) = u8(2);
           u4(3) = u8(3);
           u4(4) = u8(4);
           v4(1) = v8(1);
           v4(2) = v8(2);
           v4(3) = v8(3);
           v4(4) = v8(4);
           w4(1) = w8(1);
           w4(2) = w8(2);
           w4(3) = w8(3);
           w4(4) = w8(4);
           T%ATM%bounds = 1;
        else
           !%%Interpolate Between Z points(interpolate pts 1-4 and 5-8)
           !%Pts 1,5 : 2,6 : 3,7 : 4,8
           coord1 = (/1,2,3,4/);
           coord2 = (/5,6,7,8/);
           do ii = 1,4
              uslope = (u8(coord2(ii))-u8(coord1(ii)))/(zpts2(2)-zpts2(1));
              vslope = (v8(coord2(ii))-v8(coord1(ii)))/(zpts2(2)-zpts2(1));
              wslope = (w8(coord2(ii))-w8(coord1(ii)))/(zpts2(2)-zpts2(1));
              u4(ii) = uslope*(zstar-zpts2(1))+u8(coord1(ii));
              v4(ii) = vslope*(zstar-zpts2(1))+v8(coord1(ii));
              w4(ii) = wslope*(zstar-zpts2(1))+w8(coord1(ii));
           end do
           zpts1 = zstar;
        end if

        !%%%%%interpY%%%%%%%%%%%

        if (extrapY .eq. 1) then
           !%%You don't need to interpolate on y
           ypts1 = ypts2(1);
           u2(1) = u4(1);
           u2(2) = u4(2);
           v2(1) = v4(1);
           v2(2) = v4(2);
           w2(1) = w4(1);
           w2(2) = w4(2);
           T%ATM%bounds = 1;
        else
           !%%Interpolate between Y points(interpolate pts 1-2 and 3-4)
           !%%Pts 1,4 : 2,3
           cord1 = (/1,2/);
           cord2 = (/4,3/);
           do ii = 1,2
              uslope = (u4(cord2(ii))-u4(cord1(ii)))/(ypts2(2)-ypts2(1));
              vslope = (v4(cord2(ii))-v4(cord1(ii)))/(ypts2(2)-ypts2(1));
              wslope = (w4(cord2(ii))-w4(cord1(ii)))/(ypts2(2)-ypts2(1));
              u2(ii) = uslope*(ystar-ypts2(1))+u4(cord1(ii));
              v2(ii) = vslope*(ystar-ypts2(1))+v4(cord1(ii));
              w2(ii) = wslope*(ystar-ypts2(1))+w4(cord1(ii));
           end do
           ypts1 = ystar;
        end if

        !%%%%interpX%%%%%%%%%%%%
        if (extrapX .eq. 1) then
           !%%You don't need to interpolate on x
           xpts1 = xpts2(1);
           u = u2(1);
           v = v2(1);
           w = w2(1);
           T%ATM%bounds = 1;
        else
           !%%Interpolate between X points
           uslope = (u2(2)-u2(1))/(xpts2(2)-xpts2(1));
           vslope = (v2(2)-v2(1))/(xpts2(2)-xpts2(1));
           wslope = (w2(2)-w2(1))/(xpts2(2)-xpts2(1));
           u = uslope*(xstar-xpts2(1))+u2(1);
           v = vslope*(xstar-xpts2(1))+v2(1);
           w = wslope*(xstar-xpts2(1))+w2(1);
           xpts1 = xstar;
        end if

        !%%%%Save wind values%%%%%

        uvw(1,tt) = u;
        uvw(2,tt) = v;
        uvw(3,tt) = w;

     end do

     ! write(*,*) uvw

     if (T%ATM%TIMEVARYING .eq. 1) then
        if (extrapT .eq. 1) then
           !//Answer is just first entry of uvw
           vatm(1) = uvw(1,1);
           vatm(2) = uvw(2,1);
           vatm(3) = uvw(3,1);
        else
           !%%Interpolate on T
           tpts(1) = T%ATM%tcoord(markT)
           tpts(2) = T%ATM%tcoord(markT+1);
           u2(1) = uvw(1,1);
           u2(2) = uvw(1,2);
           v2(1) = uvw(2,1);
           v2(2) = uvw(2,2);
           w2(1) = uvw(3,1);
           w2(2) = uvw(3,2);
           uslope = (u2(2)-u2(1))/(tpts(2)-tpts(1));
           vslope = (v2(2)-v2(1))/(tpts(2)-tpts(1));
           wslope = (w2(2)-w2(1))/(tpts(2)-tpts(1));  
           u = uslope*(tstar-tpts(1))+u2(1);
           v = vslope*(tstar-tpts(1))+v2(1);
           w = wslope*(tstar-tpts(1))+w2(1);
           vatm(1) = u;
           vatm(2) = v;
           vatm(3) = w;
        end if
     else
        !//Answer is just first entry of uvw
        vatm(1) = uvw(1,1);
        vatm(2) = uvw(2,1);
        vatm(3) = uvw(3,1);
     end if

     ! write(*,*) vatm
     
     !Rotate by PSIOFFSET
     vtemp = vatm(1)
     vatm(1) = vtemp*cos(T%ATM%PSIOFFSET) + vatm(2)*sin(T%ATM%PSIOFFSET);
     vatm(2) = -vtemp*sin(T%ATM%PSIOFFSET) + vatm(2)*cos(T%ATM%PSIOFFSET);

     !//Multiply by scale
     vatm(1) = T%ATM%IWINDSCALE*vatm(1)
     vatm(2) = T%ATM%IWINDSCALE*vatm(2)
     vatm(3) = T%ATM%IWINDSCALE*vatm(3)

     if ((T%ATM%bounds .eq. 1) .and. (T%ATM%boundflag .eq. 1)) then
        !write(*,*) 'You went out of bounds at T = ',tstar,' XYZ = ',xstar,ystar,zstar
        T%ATM%boundflag = 0;
     end if
  else
     vatm(1) = 0
     vatm(2) = 0
     vatm(3) = 0
  endif

  ! if (vatm(3) .gt. 4) then
  !    write(*,*) xi,xstar,T%ATM%xshift
  !    write(*,*) yi,ystar,T%ATM%yshift
  !    write(*,*) zi,zstar,T%ATM%zshift
  !    write(*,*) vatm(1),vatm(2),vatm(3)
  !    write(*,*) 'Ypts2 = ',ypts2
  !    write(*,*) u2
  !    write(*,*) v2
  !    write(*,*) w2
  !    STOP
  ! end if

  ! write(*,*) vatm
  
  T%ATM%WRFX = vatm(1)*M2FT
  T%ATM%WRFY = vatm(2)*M2FT
  T%ATM%WRFZ = vatm(3)*M2FT

  ! write(*,*) T%ATM%WRFX,T%ATM%WRFY,T%ATM%WRFZ

  T%ATM%markX = markX
  T%ATM%markY = markY
  T%ATM%markZ = markZ
  T%ATM%markT = markT

END SUBROUTINE WRFMODEL

SUBROUTINE TURBULENCE(T)
  use MAAWMDATATYPES
  implicit none
  integer stepX,stepY,stepZ,stepT,extrapX,extrapY,extrapZ,extrapT,cord2(2)
  integer markX,markY,markZ,markT,gust,body
  integer tinterp,x1,x2,y1,y2,z1,z2,tt,ii,coord1(4),coord2(4),cord1(2)
  real*8 uvw(3,2),xpts2(2),ypts2(2),zpts2(2),zpts1,xpts1,ypts1,rx;
  real*8 u8(8),v8(8),w8(8),u4(4),v4(4),w4(4),uslope,vslope,wslope;
  real*8 u2(2),v2(2),w2(2),u,v,w,tpts(2),Lu,Lv,Lw,sigw,sigu,sigv,tstar;
  real*8 ugo,vgo,wgo,vatm(3),xstar,ystar,zstar,xtemp,vtemp,counter
  real*8 xi,yi,zi
  character*1 letter
  character*10 number
  type(MAAWMSTRUCTURE) T

  !   /*   %%This function will take in x,y,z(m),t(sec) and location and  */
  !   /*   %return u,v,w(m/s). This uses a fast quad-linear interpolation */
  !   /*   %so many globals must be defined. location is a string that */
  !   /*   %contains the location of the data to be interpolated. */

  xtemp = T%ATM%XI
  xi = xtemp*cos(T%ATM%PSIOFFSET) + T%ATM%YI*sin(T%ATM%PSIOFFSET);
  yi = -xtemp*sin(T%ATM%PSIOFFSET) + T%ATM%YI*cos(T%ATM%PSIOFFSET);

  xi = xi*FT2M !Convert from feet to meters
  yi = yi*FT2M 
  zi = zi*FT2M 

  xstar = xi + T%ATM%xshiftT
  ystar = yi + T%ATM%yshiftT
  zstar = zi + T%ATM%zshiftT

  if (T%ATM%TURBLEVEL .gt. 0) then
     stepX = 1;stepY = 1;stepZ = 1;stepT = 1;
     extrapX = 0;extrapY = 0;extrapZ = 0;extrapT = 0;
     if (zstar .lt. 0) then
        zstar = -zstar
     endif
     tinterp = 2;

     uvw(1,1)=0;uvw(2,1)=0;uvw(3,1)=0;
     uvw(1,2)=0;uvw(2,2)=0;uvw(3,2)=0;

     markX = T%ATM%markXT
     markY = T%ATM%markYT

     !%%Check X
     if (markX .eq. T%ATM%dimT) then
        markX = markX - 1;
     end if
     if ((xstar .ge. T%ATM%xcoordT(markX)) .and. (xstar .le. T%ATM%xcoordT(markX+1))) then
        !%%Your in between the markers so keep going
     else
        !%Find markX
        if (xstar .gt. T%ATM%xcoordT(T%ATM%dimT)) then
           !%use endpt
           markX = T%ATM%dimT
           stepX = -1
           extrapX = 1
           counter = 0
           do while ((T%ATM%xshiftT + xi .gt. T%ATM%xcoordT(T%ATM%dimT)) .and. (counter .lt. 20))
              T%ATM%xshiftT = -1.9 * T%ATM%xcoordT(T%ATM%dimT) + T%ATM%xshiftT
              counter = counter + 1
           end do
        else if (xstar .lt. T%ATM%xcoordT(1)) then 
           !%use starpt
           markX = 1;
           stepX = 1;
           extrapX = 1;    
           counter = 0
           do while ((T%ATM%xshiftT + xi .lt. T%ATM%xcoordT(1)) .and. (counter .lt. 20))
              T%ATM%xshiftT = 1.9 * T%ATM%xcoordT(T%ATM%dimT) + T%ATM%xshiftT
              counter = counter + 1
           end do
        else
           call FIND(T%ATM%xcoordT,T%ATM%dimT,xstar,markX)
           if (markX .eq. T%ATM%dimT) then
              markX = markX - 1;
           else if (markX .le. 0) then
              markX = 1;
           end if
        end if
     end if
     !%%Check Y
     if (markY .eq. T%ATM%dimT) then
        markY = markY - 1;
     end if
     if ((ystar .ge. T%ATM%ycoordT(markY)) .and. (ystar .le. T%ATM%ycoordT(markY+1))) then
        !%%Your in between the markers so keep going
     else
        !%Find markY
        if (ystar .gt. T%ATM%ycoordT(T%ATM%dimT)) then
           !%use endpt
           markY = T%ATM%dimT;
           stepY = -1;
           extrapY = 1;
           counter = 0
           do while ((T%ATM%yshiftT + yi .gt. T%ATM%ycoordT(T%ATM%dimT)) .and. (counter .lt. 20))
              T%ATM%yshiftT = -1.9 * T%ATM%ycoordT(T%ATM%dimT) + T%ATM%yshiftT
              counter = counter + 1
           end do
        else if (ystar .lt. T%ATM%ycoordT(1)) then
           markY = 1;
           stepY = 1;
           extrapY = 1;
           counter = 0
           do while ((T%ATM%yshiftT + yi .lt. T%ATM%ycoordT(1)) .and. (counter .lt. 20))
              T%ATM%yshiftT = 1.9 * T%ATM%ycoordT(T%ATM%dimT) + T%ATM%yshiftT
              counter = counter + 1;
           end do
        else
           call FIND(T%ATM%ycoordT,T%ATM%dimT,ystar,markY)
           if (markY .eq. T%ATM%dimT) then
              markY = markY - 1;
           else if (markY .le. 0) then
              markY = 1;
           end if
        end if
     end if

     !We start with 8 points
     !%%To start we have 4 discrete point (4 corners of a square)
     xpts2(1) = T%ATM%xcoordT(markX)
     xpts2(2) = T%ATM%xcoordT(markX+stepX);
     ypts2(1) = T%ATM%ycoordT(markY)
     ypts2(2) = T%ATM%ycoordT(markY+stepY);
     x1 = markX;x2 = markX+stepX;
     y1 = markY;y2 = (markY+stepY);
     !%%Use U0,V0,W0
     u4(1) = T%ATM%UTURB(y1,x1)
     u4(2) = T%ATM%UTURB(y1,x2)
     u4(3) = T%ATM%UTURB(y2,x2)
     u4(4) = T%ATM%UTURB(y2,x1)
     v4(1) = T%ATM%VTURB(y1,x1)
     v4(2) = T%ATM%VTURB(y1,x2)
     v4(3) = T%ATM%VTURB(y2,x2)
     v4(4) = T%ATM%VTURB(y2,x1)
     w4(1) = T%ATM%WTURB(y1,x1)
     w4(2) = T%ATM%WTURB(y1,x2)
     w4(3) = T%ATM%WTURB(y2,x2)
     w4(4) = T%ATM%WTURB(y2,x1)

     !%%%%%interpY%%%%%%%%%%%
     if (extrapY .eq. 1) then
        !%%You don't need to interpolate on y
        ypts1 = ypts2(1);
        u2(1) = u4(1);
        u2(2) = u4(2);
        v2(1) = v4(1);
        v2(2) = v4(2);
        w2(1) = w4(1);
        w2(2) = w4(2);
        T%ATM%boundsT = 1;
     else
        !%%Interpolate between Y points(interpolate pts 1-2 and 3-4)
        !%%Pts 1,4 : 2,3
        cord1 = (/1,2/);
        cord2 = (/4,3/);
        do ii = 1,2
           uslope = (u4(cord2(ii))-u4(cord1(ii)))/(ypts2(2)-ypts2(1));
           vslope = (v4(cord2(ii))-v4(cord1(ii)))/(ypts2(2)-ypts2(1));
           wslope = (w4(cord2(ii))-w4(cord1(ii)))/(ypts2(2)-ypts2(1));
           u2(ii) = uslope*(ystar-ypts2(1))+u4(cord1(ii));
           v2(ii) = vslope*(ystar-ypts2(1))+v4(cord1(ii));
           w2(ii) = wslope*(ystar-ypts2(1))+w4(cord1(ii));
        end do
        ypts1 = ystar;
     end if

     !%%%%interpX%%%%%%%%%%%%
     if (extrapX .eq. 1) then
        !%%You don't need to interpolate on x
        xpts1 = xpts2(1);
        u = u2(1);
        v = v2(1);
        w = w2(1);
        T%ATM%boundsT = 1;
     else
        !%%Interpolate between X points
        uslope = (u2(2)-u2(1))/(xpts2(2)-xpts2(1));
        vslope = (v2(2)-v2(1))/(xpts2(2)-xpts2(1));
        wslope = (w2(2)-w2(1))/(xpts2(2)-xpts2(1));
        u = uslope*(xstar-xpts2(1))+u2(1);
        v = vslope*(xstar-xpts2(1))+v2(1);
        w = wslope*(xstar-xpts2(1))+w2(1);
        xpts1 = xstar;
     end if

     !%%%%Save wind values%%%%%

     vatm(1) = u
     vatm(2) = v
     vatm(3) = w

     !Rotate by PSIOFFSET
     vtemp = vatm(1)
     vatm(1) = vtemp*cos(T%ATM%PSIOFFSET) + vatm(2)*sin(T%ATM%PSIOFFSET);
     vatm(2) = -vtemp*sin(T%ATM%PSIOFFSET) + vatm(2)*cos(T%ATM%PSIOFFSET);

     !//Multiply by scale
     vatm(1) = T%ATM%TURBLEVEL*vatm(1)
     vatm(2) = T%ATM%TURBLEVEL*vatm(2)
     vatm(3) = T%ATM%TURBLEVEL*vatm(3)

     if ((T%ATM%boundsT .eq. 1) .and. (T%ATM%boundflagT .eq. 1)) then
        !write(*,*) 'You went out of turbulence bounds at T = ',tstar,' XYZ = ',xstar,ystar,zstar
        T%ATM%boundflagT = 0;
     end if
  else
     vatm(1) = 0
     vatm(2) = 0
     vatm(3) = 0
  endif

  T%ATM%WINDGUST(1) = vatm(1)*M2FT  !Convert from m/s back to ft/s
  T%ATM%WINDGUST(2) = vatm(2)*M2FT
  T%ATM%WINDGUST(3) = vatm(3)*M2FT

  T%ATM%markXT = markX
  T%ATM%markYT = markY
  T%ATM%markZT = markZ

END SUBROUTINE TURBULENCE

SUBROUTINE TUSTINDERIVATIVEFILTER(T,delydot,dely)
  use MAAWMDATATYPES
  implicit none
  real*8 angles(2,1),rGS_I(3,1),theta,psi,norm,delz,dely,dt,tau,wc
  real*8 alfa,y(2),delydot
  type(MAAWMSTRUCTURE) T
  
  dt = T%SIM%DELTATIME
  wc = 10 !rad/s
  tau = 1/wc;
  alfa = 2*tau/dt;
  delydot = (2*(dely-T%CS%DELPREV)+T%CS%DELDOTPREV*(2*tau-dt))/(2*tau+dt);

END SUBROUTINE TUSTINDERIVATIVEFILTER

SUBROUTINE COMPUTENCENS(T)
  use MAAWMDATATYPES
  implicit none
  type(MAAWMSTRUCTURE) T

  !T%SIM%NCENS = sqrt(T%SIM%RBFDATA/1.4D0)
  !T%SIM%NCENS = 16 !256
  !T%SIM%NCENS = 16 !225
  write(*,*) 'Number of Data Points Sampled = ',T%SIM%NDATA
  write(*,*) 'Number of Centers Per Axis = ',T%SIM%NCENS

END SUBROUTINE COMPUTENCENS

SUBROUTINE FEEDBACK(T)
 use MAAWMDATATYPES
 implicit none
 integer i,stateindex,ncenin
 real*8 angles(2,1),rGS_I(3,1),theta,psi,norm,delz,dely,n,rGS_S(3,1)
 type(MAAWMSTRUCTURE) T

 if (T%SIM%TIME .gt. T%CS%UPDATETIME) then

    !Step UPDATETIME
    T%CS%UPDATETIME = T%CS%UPDATETIME + 1/T%CS%UPDATERATE

    !Have airplanes measure the atmosphere and output the data to a file
    do i = 1,T%SIM%NOAIRPLANES
       T%ATM%XI = T%AC(i)%STATE(1)
       T%ATM%YI = T%AC(i)%STATE(2)
       T%ATM%ZI = T%AC(i)%STATE(3)
       !write(*,*) 'Windspeed measurements for airplane = ',i
       call SURROGENCLEANUP(T,i)
    end do

    !Have Quadcopters Measure the atmosphere and output the data to a file
    do i = 1,T%SIM%NOQUADCOPTERS
       T%ATM%XI = T%COPTER(i)%STATE(1)
       T%ATM%YI = T%COPTER(i)%STATE(2)
       T%ATM%ZI = T%COPTER(i)%STATE(3)
       !write(*,*) 'Windspeed measurements for quadcopter = ',i
       call SURROGENCLEANUP(T,i+T%SIM%NOAIRPLANES)
    end do

 endif

END SUBROUTINE FEEDBACK

SUBROUTINE SURROGENCLEANUP(T,i)
  use MAAWMDATATYPES
  implicit none
  integer i
  real*8 n1,n2,n3
  type(MAAWMSTRUCTURE) T

  !write(*,*) T%ATM%XI,T%ATM%YI,T%ATM%ZI

  call ATMOSPHERE(T,3)

  !write(*,*) T%ATM%XI,T%ATM%YI,T%ATM%ZI
  
  !write(*,*) 'Windspeed computed'
  
  T%SIM%NDATA = T%SIM%NDATA + 1
  
  if (T%SIM%NDATA .eq. NDATAMAX) then
     write(*,*) 'Maximum Number of data points allowed for RBF routine exceeded in SURROGENCLEANUP'
     write(*,*) 'Time Simulation = ',T%SIM%TIME
     write(*,*) 'Number of Aircraft = ',T%SIM%NOAIRCRAFT
     write(*,*) 'Number of Quadcopters = ',T%SIM%NOQUADCOPTERS
     STOP
  end if

  !Add Sensor Errors here
  if (T%CS%SENSORERRORS .gt. 0) then
     call RandUniform(n1) !Number b/t 0 and 1
     call RandUniform(n2) !Number b/t 0 and 1
     call RandUniform(n3) !Number b/t 0 and 1
  else
     n1 = 0.5
     n2 = 0.5
     n3 = 0.5
  end if
  
  write(23,fmt='(1000F30.10)') T%ATM%XI,T%ATM%YI,T%ATM%ZI,T%SIM%TIME,(1+T%CS%SENSORERRORS*(n1-0.5))*T%ATM%VXWIND, (1+T%CS%SENSORERRORS*(n2-0.5))*T%ATM%VYWIND, (1+T%CS%SENSORERRORS*(n3-0.5))*T%ATM%VZWIND,real(i)
  !write(*,*) T%ATM%XI,T%ATM%YI,T%ATM%ZI,T%SIM%TIME,T%ATM%VXWIND, T%ATM%VYWIND, T%ATM%VZWIND,real(i)
  
  if (T%SIM%RBFDATA .gt. 0) then
     if (T%SIM%NDATA .eq. T%SIM%RBFDATA) then
          close(23)
          write(*,*) 'Number of Data Points Equals Number of RBF Data Points'
          !Compute number of centers along each axis
          !ncentotal = T%SIM%NDATA/1.4D0
          call COMPUTENCENS(T)
          call SURROGEN(T%SIM%IKIND,T%SIM%K35,T%SIM%TIME,T%SIM%NCENS)
          STOP
       end if
   end if
   
END SUBROUTINE SURROGENCLEANUP

SUBROUTINE QUAD_HANDSHAKE_INIT(T,i)
  use MAAWMDATATYPES
  implicit none
  integer i
  type(MAAWMSTRUCTURE) T
  write(*,*) 'Handshaking...',i
  !Pass some tether dynamics stuff to QUAD
  T%COPTER(i)%THR_DYNOFFON = 0
  T%COPTER(i)%THR_ELASOFFON = 0
  !Basically pass all things that don't change
  T%COPTER(i)%DELTATIME = T%SIM%DELTATIME
  T%COPTER(i)%DYNOFFON = T%COPTER(i)%OFFON
  T%COPTER(i)%MINWAYPOINT = 10.0D0
  call QUAD_HANDSHAKE(T,i)
end SUBROUTINE QUAD_HANDSHAKE_INIT

SUBROUTINE QUAD_HANDSHAKE(T,i)
  use MAAWMDATATYPES
  implicit none
  integer i
  type(MAAWMSTRUCTURE) T
  !Pass Time to COPTER Module
  !T%COPTER%TIME = T%SIM%TIME - T%COPTER%TIMEON
  T%COPTER(i)%TIME = T%SIM%TIME
  !In order to compute the aero model properly you need
  !to pass x,y,z to the ATMOSPHERE model
  if (T%COPTER(i)%AEROOFFON .eq. 1) then
     T%ATM%XI = T%COPTER(i)%STATE(1)
     T%ATM%YI = T%COPTER(i)%STATE(2)
     T%ATM%ZI = T%COPTER(i)%STATE(3)
     !Compute Atmopsheric density and winds - Same for Quad
     call ATMOSPHERE(T,3) !T%ATM%DEN
     T%COPTER(i)%VXWIND = T%ATM%VXWIND
     T%COPTER(i)%VYWIND = T%ATM%VYWIND
     T%COPTER(i)%VZWIND = T%ATM%VZWIND
     T%COPTER(i)%DEN = T%ATM%DEN
  end if
  !Pass some tether dynamics stuff to QUAD
  T%COPTER(i)%FTETHERX = 0
  T%COPTER(i)%FTETHERY = 0
  T%COPTER(i)%FTETHERZ = 0
end SUBROUTINE QUAD_HANDSHAKE
