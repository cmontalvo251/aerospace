!!!!SIXDOF CODE USING FORTRAN 90 STANDARD

!!!!!!!!!!!!!!!! COPTER STRUCTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BODYDATATYPES
  IMPLICIT NONE
  integer,parameter :: MAXNBODY = 10             ! Units: 'nd', Desc: 'Maximum Number of BODY Motion Table Points'
  integer,parameter :: MAXWP = 1000             ! Units: 'nd', Desc: 'Maximum Number of BODY Motion Table Points'
  integer,parameter :: IMAX = 1                    ! Units: 'nd', Desc: 'Maximum Airwake index in X
  integer,parameter :: JMAX = 1                    ! Units: 'nd', Desc: 'Maximum Airwake index in Y
  integer,parameter :: KMAX = 1                    ! Units: 'nd', Desc: 'Maximum Airwake index in Z
  integer,parameter :: NTIMES = 1                  ! Units: 'nd', Desc: 'Maximum Airwake index for t
  integer,parameter :: MAXNBEADS = 100              ! Units: 'nd', Desc: 'Maximum Number of Tether Beads'
  integer,parameter :: MAXNALT = 200                ! Units: 'nd', Desc: 'Maximum Number of Atmosphere Altitude Table Points'
  integer,parameter :: MAXNLSE = 20                 ! Units: 'nd', Desc: 'Maximum Number of Lifting Surface Elements'
  integer,parameter :: MAXPROP = 20                 ! Units: 'nd', Desc: 'Maximum Number of Propellor Table size
  integer,parameter :: MAXNAOA = 100                ! Units: 'nd', Desc: 'Maximum Number of Aerodynamic Angle of Attack Table Points'
  integer,parameter :: MAXX = 1000                  ! Units: 'nd', Desc: 'Maximum Number of System States'
  integer,parameter :: NOACTUATORS = 9              ! Units: 'nd', Desc: Number of Actuators
  integer,parameter :: SWEEPTABLE = 101             ! Units: 'nd', Desc: Sweep Table length
  integer,parameter :: FT2M = 0.3048                ! Conversion from Feet to Meters
  integer,parameter :: M2FT = 3.28084               ! Conversion from Meters to Feet
  integer,parameter  :: NUMSTATES = 33
  real*8,parameter  :: PI = 3.14159265358979323846  ! Units: 'nd', Desc: 'Pi'
  real*8,parameter  :: qPI = 3.14159265358979323846  ! Units: 'nd', Desc: 'Pi' - using qPI to not conflict with any other usage of PI
  type BODYSTRUCTURE
     integer :: OFFON = 0                             ! Units: 'nd', Desc: 'Off/On Switch'
     integer :: DYNOFFON = 0                             ! Units: 'nd', Desc: 'Off/On Switch'
     integer :: GRAVOFFON = 0                         ! Units: 'nd', Desc: 'Gravity Flag (0=Off, 1=On)'
     integer :: AEROOFFON = 0                         ! Units: 'nd', Desc: 'Aerodynamics Flag (0=Off, 1=On)'
     integer :: CONTOFFON = 0                         ! Units: 'nd', Desc: 'Contact Flag (0=Off, 1=On)'
     integer :: PHASE = 1                               ! Phasing for autopilot
     integer :: MODNO = 1                             ! Units: 'nd', Desc: 'Model Number Flag'
     integer :: IP = 1                                ! Units: 'nd', Desc: 'Time Table Look-Up Pointer'
     integer :: TABSIZE = 1                           ! Units: 'nd', Desc: 'Time Table Size'
     integer :: DQFLAG = 0                            ! Units: 'nd', Desc: 'Data Quality Flag (0=Data Not Loaded Successfully, 1=Data Loaded Successfully)'
     integer :: AIRWAKE = 0                           ! Units: 'nd', Desc: 'Air wake on or off'
     integer :: markX = 1                             ! Units: 'nd', X index placeholder
     integer :: markY = 1                             ! Units: 'nd', X index placeholder
     integer :: markZ = 1                             ! Units: 'nd', X index placeholder
     integer :: markT = 1                             ! Units: 'nd', X index placeholder
     integer :: DOWNWASHONOFF = 0                     ! Units: 'nd', X index placeholder
     integer :: THR_DYNOFFON = 0
     integer :: THR_ELASOFFON = 0
     integer :: CONTROLOFFON = 0
     integer :: WAYPOINT = 1                          ! Defaults to 1 but apparently can change in input file
     integer :: NUMWAYPOINTS = 0                         ! Defaults to 1 but apparently can change in input file
     integer :: COMMANDFLAG(MAXWP) = 1           ! 0 or 1 to advance waypoint command     real*8 :: COMMANDFLAG(MAXWP) =
     real*8 :: ALC = 0                                ! Body aero parameter
     real*8 :: ALS = 0                                ! Body aero parameter
     real*8 :: DXD = 0                                ! Body aero parameter
     real*8 :: DYD = 0                                ! Body aero parameter
     real*8 :: RNEW = 0                               ! Body aero parameter
     real*8 :: C_T = 0                                ! Body aero parameter
     real*8 :: C_Q = 0                              ! Body aero parameter
     real*8 :: LPHI12 = 0                             ! Body aero parameter
     real*8 :: LPHI34 = 0                             ! Body aero parameter
     real*8 :: LTHETA12 = 0                           ! Body aero parameter
     real*8 :: LTHETA34 = 0                           ! Body aero parameter
     real*8 :: OMEGAMAX = 0                           ! Body aero parameter
     real*8 :: OMEGAVEC(9,1) = 0                      ! Body aero parameter
     real*8 :: THRUSTVEC(9,1) = 0                     ! Body aero parameter
     real*8 :: MUVEC(9,1) = 0                         ! Body aero parameter
     real*8 :: KT = 0                                 ! Body Aero Parameter
     real*8 :: OMEGA0 = 0                             ! Body Aero Parameter
     real*8 :: IRR = 0                                ! Body Aero Parameter
     real*8 :: PSILOC = 0
     real*8 :: MINWAYPOINT = 1.0D0
     real*8 :: MS_MIN = 0.0                           ! Maximum PWM signal (microseconds)
     real*8 :: MS_MAX = 0.0                           ! Minimum PWM signal (microseconds)
     real*8 :: DOWNWASH = 0.0                         ! Units: 'nd', Magnitude of Downwash
     real*8 :: DIAMETER = 0.0                         ! Units: 'nd', Diameter of Downwash (ft)
     real*8 :: SLCG = 0.0                             ! Units: 'ft', Desc: 'Stationline of Mass Center of BODY'
     real*8 :: BLCG = 0.0                             ! Units: 'ft', Desc: 'Buttline of Mass Center of BODY'
     real*8 :: WLCG = 0.0                             ! Units: 'ft', Desc: 'Waterline of Mass Center of BODY'
     real*8 :: SLREEL = 0.0                           ! Units: 'ft', Desc: 'Stationline of Tether Reel Point on BODY'
     real*8 :: BLREEL = 0.0                           ! Units: 'ft', Desc: 'Buttline of Tether Reel Point on BODY'
     real*8 :: WLREEL = 0.0                           ! Units: 'ft', Desc: 'Waterline of Tether Reel Point on BODY'
     real*8 :: SLAIRWAKE = 0.0                        ! Units: 'ft', Desc: 'Stationline of Airwake grid start on BODY
     real*8 :: BLAIRWAKE = 0.0                        ! Units: 'ft', Desc: 'Buttline of Airwake grid start on BODY
     real*8 :: WLAIRWAKE = 0.0                        ! Units: 'ft', Desc: 'Waterline of Airwake grid start on BODY
     real*8 :: SPEED = 0.0                            ! Units: 'ft/s', Desc: 'Speed of BODY'
     real*8 :: FINALSPEED = 0.0                       ! Units: 'ft/s', Desc: 'Initial Speed of BODY'
     real*8 :: RESTARTSPEED = -999                    ! Units: 'ft/s', Desc: 'Restart speed of BODYBODY' 
     real*8 :: DIRECTION = 0.0                        ! Units: 'rad', Desc: 'Direction of Motion of BODY'
     real*8 :: XCGINITIAL = 0.0                       ! Units: 'ft', Desc: 'Initial XCG Inertial Position of BODY'
     real*8 :: YCGINITIAL = 0.0                       ! Units: 'ft', Desc: 'Initial YCG Inertial Position of BODY'
     real*8 :: ZCGINITIAL = 0.0                       ! Units: 'ft', Desc: 'Initial ZCG Inertial Position of BODY'
     real*8 :: TIME = 0.0                             ! Units: 's', Desc: 'Time'
     real*8 :: TIMEON = 0                             ! Units: 's', Desc: Time for BODY to come on
     real*8 :: XCG = 0.0                              ! Units: 'ft', Desc: 'XCG Inertial Position of BODY'
     real*8 :: YCG = 0.0                              ! Units: 'ft', Desc: 'YCG Inertial Position of BODY'
     real*8 :: ZCG = 0.0                              ! Units: 'ft', Desc: 'ZCG Inertial Position of BODY'
     real*8 :: PHI = 0.0                              ! Units: 'rad', Desc: 'BODY Euler Roll Angle'
     real*8 :: THETA = 0.0                            ! Units: 'rad', Desc: 'BODY Euler Pitch Angle'
     real*8 :: PSI = 0.0                              ! Units: 'rad', Desc: 'BODY Euler Yaw Angle'
     real*8 :: PSIPREV = 0.0                          ! Units: 'rad', Desc: 'BODY Euler Yaw Angle'
     real*8 :: UB = 0.0                               ! Units: 'ft/s', Desc: 'UB Body Frame Mass Center Velocity of BODY'
     real*8 :: VB = 0.0                               ! Units: 'ft/s', Desc: 'VB Body Frame Mass Center Velocity of BODY'
     real*8 :: WB = 0.0                               ! Units: 'ft/s', Desc: 'WB Body Frame Mass Center Velocity of BODY'
     real*8 :: PB = 0.0                               ! Units: 'rad/s', Desc: 'PB Roll Rate of BODY'
     real*8 :: QB = 0.0                               ! Units: 'rad/s', Desc: 'QB Pitch Rate of BODY'
     real*8 :: RB = 0.0                               ! Units: 'rad/s', Desc: 'RB Yaw Rate of BODY'
     real*8 :: XDOT = 0                               ! Units: 'ft/s', Desc: Inertial Frame X velocity of BODY
     real*8 :: YDOT = 0                               ! Units: 'ft/s', Desc: Inertial Frame X velocity of BODY
     real*8 :: ZDOT = 0                               ! Units: 'ft/s', Desc: Inertial Frame X velocity of BODY
     real*8 :: XREEL = 0.0                            ! Units: 'ft', Desc: 'XREEL Inertial Position of Reel'
     real*8 :: YREEL = 0.0                            ! Units: 'ft', Desc: 'YREEL Inertial Position of Reel'
     real*8 :: ZREEL = 0.0                            ! Units: 'ft', Desc: 'ZREEL Inertial Position of Reel'
     real*8 :: XREELDOT = 0.0                         ! Units: 'ft/s', Desc: 'XREEL DOT Inertial Velocity of Reel'
     real*8 :: YREELDOT = 0.0                         ! Units: 'ft/s', Desc: 'YREEL DOT Inertial Velocity of Reel'
     real*8 :: ZREELDOT = 0.0                         ! Units: 'ft/s', Desc: 'ZREEL DOT Inertial Velocity of Reel'
     real*8 :: UREEL = 0.0                            ! Units: 'ft', Desc: 'UREEL Inertial Velcity of Reel in BODY Reference Frame'
     real*8 :: VREEL = 0.0                            ! Units: 'ft', Desc: 'VREEL Inertial Velcity of Reel in BODY Reference Frame'
     real*8 :: WREEL = 0.0                            ! Units: 'ft', Desc: 'WREEL Inertial Velcity of Reel in BODY Reference Frame'
     !REVISIT REVISIT REVISIT
     ! real*8 :: UBODY(1,1,1),VBODY(1,1,1),WBODY(1,1,1) ! U,V,W BODY airwake velocity at timestep t (ft/s)
     ! real*8 :: UBODYDT(1,1,1),VBODYDT(1,1,1),WBODYDT(1,1,1) ! U,V,W airwake velocity at timestep t+dt (ft/s)
     ! real*8 :: XBODY(1,1,1),YBODY(1,1,1),ZBODY(1,1,1) ! U,V,W airwake velocity at timestep t+dt (ft/s)
     real*8 :: UBODY(IMAX,JMAX,KMAX),VBODY(IMAX,JMAX,KMAX),WBODY(IMAX,JMAX,KMAX) ! U,V,W BODY airwake velocity at timestep t (ft/s)
     real*8 :: UBODYDT(IMAX,JMAX,KMAX),VBODYDT(IMAX,JMAX,KMAX),WBODYDT(IMAX,JMAX,KMAX) ! U,V,W airwake velocity at timestep t+dt (ft/s)
     real*8 :: XBODY(IMAX,JMAX,KMAX),YBODY(IMAX,JMAX,KMAX),ZBODY(IMAX,JMAX,KMAX) ! U,V,W airwake velocity at timestep t+dt (ft/s)
     real*8 :: TCOORD(NTIMES)                            ! Units: 's' Desc: Airwake time intervals
     real*8 :: XCOORD(IMAX)                            ! Units: 's' Desc: Airwake time intervals
     real*8 :: YCOORD(JMAX)                            ! Units: 's' Desc: Airwake time intervals
     real*8 :: ZCOORD(KMAX)                            ! Units: 's' Desc: Airwake time intervals
     real*8 :: TIMETAB(MAXNBODY) = 0.0           ! Units: 's', Desc: 'Time Table'
     real*8 :: XCGTAB(MAXNBODY) = 0.0            ! Units: 'ft', Desc: 'XCG Inertial Position of BODY Table'
     real*8 :: YCGTAB(MAXNBODY) = 0.0            ! Units: 'ft', Desc: 'YCG Inertial Position of BODY Table'
     real*8 :: ZCGTAB(MAXNBODY) = 0.0            ! Units: 'ft', Desc: 'ZCG Inertial Position of BODY Table'
     real*8 :: PHITAB(MAXNBODY) = 0.0            ! Units: 'rad', Desc: 'BODY Euler Roll Angle Table'
     real*8 :: THETATAB(MAXNBODY) = 0.0          ! Units: 'rad', Desc: 'BODY Euler Pitch Angle Table'
     real*8 :: PSITAB(MAXNBODY) = 0.0            ! Units: 'rad', Desc: 'BODY Euler Yaw Angle Table'
     real*8 :: UBTAB(MAXNBODY) = 0.0             ! Units: 'ft/s', Desc: 'UB Body Frame Mass Center Velocity of BODY Table'
     real*8 :: VBTAB(MAXNBODY) = 0.0             ! Units: 'ft/s', Desc: 'VB Body Frame Mass Center Velocity of BODY Table'
     real*8 :: WBTAB(MAXNBODY) = 0.0             ! Units: 'ft/s', Desc: 'WB Body Frame Mass Center Velocity of BODY Table'
     real*8 :: PBTAB(MAXNBODY) = 0.0             ! Units: 'rad/s', Desc: 'PB Roll Rate of BODY Table'
     real*8 :: QBTAB(MAXNBODY) = 0.0             ! Units: 'rad/s', Desc: 'QB Pitch Rate of BODY Table'
     real*8 :: RBTAB(MAXNBODY) = 0.0             ! Units: 'rad/s', Desc: 'RB Yaw Rate of BODY Table'
     real*8 :: INITIALSTATE(NUMSTATES) = 0                    ! Units: 'varies', Desc: 'Initial State Vector'
     real*8 :: NOMINALSTATE(NUMSTATES) = 0                    ! Units: 'varies', Desc: 'Initial State Vector'
     real*8 :: STATE(NUMSTATES) = 0                           ! Units: 'varies', Desc: 'State Vector'
     real*8 :: STATEDOT(NUMSTATES) = 0                        ! Units: 'varies', Desc: 'State Dot Vector'
     real*8 :: KRKBODY(NUMSTATES,4) = 0.0                  ! Units: 'md', RK4 vector
     real*8 :: TIS(3,3) = 0.0                          ! Units: 'nd', Desc: 'Transformation from BODY to Inertial Reference Frame'
     real*8 :: TIC(3,3) = 0.0                          ! Units: 'nd', Desc: 'Transformation from BODY to Inertial Reference Frame'
     real*8 :: TCI(3,3) = 0.0                          ! Units: 'nd', Desc: 'Transformation from BODY to Inertial Reference Frame'
     real*8 :: TBI(2,2) = 0.0
     real*8 :: XDDOTNOISE = 0.0                        ! Units: 'ft/s^2', Desc: 'Acceleration Noise of BODYBODY'
     real*8 :: YDDOTSCALE = 0.0                        ! Units: 'ft/s^2', Desc: 'Acceleration Noise of BODYBODY'
     real*8 :: YDDOTPERIOD = 0.0                        ! Units: 'ft/s^2', Desc: 'Acceleration Noise of BODYBODY'
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
     real*8 :: MASS = 0.0                             ! Units: 'kg', Desc: 'Mass'
     real*8 :: WEIGHT = 0.0                           ! Units: 'lbf', Desc: 'Weight'
     real*8 :: TURNRADIUS = 0.0                       ! Units: 'ft', Desc: Minimum turn radius of bodyBODY
     real*8 :: FXTOTAL = 0.0                          ! Units: 'lbf', Desc: 'X Total Applied Force in Body Frame'
     real*8 :: FYTOTAL = 0.0                          ! Units: 'lbf', Desc: 'Y Total Applied Force in Body Frame'
     real*8 :: FZTOTAL = 0.0                          ! Units: 'lbf', Desc: 'Z Total Applied Force in Body Frame'
     real*8 :: MXTOTAL = 0.0                          ! Units: 'N m', Desc: 'X Total Applied Moment About Mass Center in Body Frame'
     real*8 :: MYTOTAL = 0.0                          ! Units: 'N m', Desc: 'Y Total Applied Moment About Mass Center in Body Frame'
     real*8 :: MZTOTAL = 0.0                          ! Units: 'N m', Desc: 'Z Total Applied Moment About Mass Center in Body Frame'
     real*8 :: FXGRAV = 0.0                           ! Units: 'lbf', Desc: 'X Gravity Forces in Body Frame'
     real*8 :: FYGRAV = 0.0                           ! Units: 'lbf', Desc: 'Y Gravity Forces in Body Frame'
     real*8 :: FZGRAV = 0.0                           ! Units: 'lbf', Desc: 'Z Gravity Forces in Body Frame'
     real*8 :: MXGRAV = 0.0                           ! Units: 'N m', Desc: 'X Gravity Moment About CG in Body Frame'
     real*8 :: MYGRAV = 0.0                           ! Units: 'N m', Desc: 'Y Gravity Moment About CG in Body Frame'
     real*8 :: MZGRAV = 0.0                           ! Units: 'N m', Desc: 'Z Gravity Moment About CG in Body Frame'
     real*8 :: FXAERO = 0.0                           ! Units: 'lbf', Desc: 'X Aerodynamic Force in Body Frame'
     real*8 :: FYAERO = 0.0                           ! Units: 'lbf', Desc: 'Y Aerodynamic Force in Body Frame'
     real*8 :: FZAERO = 0.0                           ! Units: 'lbf', Desc: 'Z Aerodynamic Force in Body Frame'
     real*8 :: MXAERO = 0.0                           ! Units: 'N m', Desc: 'X Aerodynamic Moment About CG in Body Frame'
     real*8 :: MYAERO = 0.0                           ! Units: 'N m', Desc: 'Y Aerodynamic Moment About CG in Body Frame'
     real*8 :: MZAERO = 0.0                           ! Units: 'N m', Desc: 'Z Aerodynamic Moment About CG in Body Frame'
     real*8 :: FXCONT = 0.0                           ! Units: 'lbf', Desc: 'X Contact Force in Body Frame'
     real*8 :: FYCONT = 0.0                           ! Units: 'lbf', Desc: 'Y Contact Force in Body Frame'
     real*8 :: FZCONT = 0.0                           ! Units: 'lbf', Desc: 'Z Contact Force in Body Frame'
     real*8 :: MXCONT = 0.0                           ! Units: 'N m', Desc: 'X Contact Moment in Body Frame'
     real*8 :: MYCONT = 0.0                           ! Units: 'N m', Desc: 'Y Contact Moment in Body Frame'
     real*8 :: MZCONT = 0.0                           ! Units: 'N m', Desc: 'Z Contact Moment in Body Frame'
     real*8 :: VXWIND = 0.0                           ! Units: 'ft/s', Desc: 'Atmospheric Wind Along Inertial I Axis'
     real*8 :: VYWIND = 0.0                           ! Units: 'ft/s', Desc: 'Atmospheric Wind Along Inertial J Axis'
     real*8 :: VZWIND = 0.0                           ! Units: 'ft/s', Desc: 'Atmospheric Wind Along Inertial K Axis'
     real*8 :: DEN = 0.002363                         ! Units: 'slug/ft^3', Desc: 'Density'
     real*8 :: FTETHERX = 0.0
     real*8 :: FTETHERY = 0.0
     real*8 :: FTETHERZ = 0.0
     real*8 :: VWAKE(3) = 0.0
     real*8 :: AIRWAKEPOSITION(3) = 0.0
     real*8 :: GRAVITY = 32.2                       ! Units: 'ft/s this is the default but you can change it in the input file
     real*8 :: XCOMMAND = 0.0
     real*8 :: YCOMMAND = 0.0
     real*8 :: ZCOMMAND = 0.0
     real*8 :: PHICOMMAND = 0.0
     real*8 :: THETACOMMAND = 0.0
     real*8 :: PSICOMMAND = 0.0
     real*8 :: KPXBODY = 0                            ! Units: 'nd', Desc: bodyBODY proportional gain on x-position
     real*8 :: KIXBODY = 0                            ! Units: 'nd', Desc: bodyBODY integral gain on x-position
     real*8 :: KDXBODY = 0                            ! Units: 'nd', Desc: bodyBODY derivative gain on x-position
     real*8 :: KPYBODY = 0                            ! Units: 'nd', Desc: bodyBODY proportional gain on y-position
     real*8 :: KIYBODY = 0                            ! Units: 'nd', Desc: bodyBODY integral gain on y-position
     real*8 :: KDYBODY = 0                            ! Units: 'nd', Desc: bodyBODY derivative gain on y-position
     real*8 :: KPZBODY = 0                            ! Units: 'nd', Desc: bodyBODY proportional gain on z-position
     real*8 :: KIZBODY = 0                            ! Units: 'nd', Desc: bodyBODY integral gain on z-position 
     real*8 :: KDZBODY = 0                            ! Units: 'nd', Desc: bodyBODY derivative gain on z-position
     real*8 :: KPPSI = 0                              ! Units: 'nd', Desc: bodyBODY proportional gain on heading
     real*8 :: KIPSI = 0                              ! Units: 'nd', Desc: bodyBODY integral gain on heading
     real*8 :: KDPSI = 0                              ! Units: 'nd', Desc: bodyBODY derivative gain on heading
     real*8 :: KPPHI = 0                              ! Units: 'nd', Desc: bodyBODY proportional gain on roll 
     real*8 :: KIPHI = 0                              ! Units: 'nd', Desc: bodyBODY integral gain on roll 
     real*8 :: KDPHI = 0                              ! Units: 'nd', Desc: bodyBODY derivative gain on roll
     real*8 :: KPTHETA = 0                            ! Units: 'nd', Desc: bodyBODY proportional gain on pitch
     real*8 :: KITHETA = 0                            ! Units: 'nd', Desc: bodyBODY integral gain on pitch
     real*8 :: KDTHETA = 0                            ! Units: 'nd', Desc: bodyBODY derivative gain on pitch
     real*8 :: XINTEGRAL = 0
     real*8 :: YINTEGRAL = 0
     real*8 :: ZINTEGRAL = 0
     real*8 :: PHIINTEGRAL = 0.0                      ! Units: 'ft', Desc: 'PHI integral term in PID Controller'
     real*8 :: THETAINTEGRAL = 0.0                    ! Units: 'ft', Desc: 'THETA integral term in PID Controller'
     real*8 :: PSIINTEGRAL = 0.0                      ! Units: 'ft', Desc: 'PSI integral term in PID Controller'
     real*8 :: MS_0 = 0.0                             ! Units: 'us', Desc: Microsecond impulse to keep body hovering // REVISIT can take this out later
     real*8 :: MS_ROLL = 0.0                          ! Units: 'us', Desc: Microsecond impulse to roll  // REVISIT can take this out later
     real*8 :: MS_PITCH = 0.0                         ! Units: 'us', Desc: Microsecond impulse to pitch // REVISIT can take this out later
     real*8 :: MS_YAW = 0.0                           ! Units: 'us', Desc: Microsecond impulse to yaw   // REVISIT can take this out later
     real*8 :: PWM2F(9,1) = 0.0                       ! Units: 'lbf',Desc: Force converted as a function of input microsecond pulse
     real*8 :: DELTATIME = 0.0                        ! Simulation timestep
     real*8 :: XCOM(MAXWP) = 500                         ! Units: 'm', X Waypoint command
     real*8 :: YCOM(MAXWP) = 500                         ! Units: 'm', Y Waypoint command
     real*8 :: ZCOM(MAXWP) = -200                        ! Units: 'm', Altitude command
     character*256 U0name,V0name,W0name
     character*256 Udtname,Vdtname,Wdtname,AIRWAKEPATH
     character(128) :: INPUTFILE = ' '                ! Units: 'nd', Desc: 'Body Input File'     chara
  end type BODYSTRUCTURE
end module BODYDATATYPES

PROGRAM SIXDOF
 use BODYDATATYPES
 !USE Math_Module, ONLY: Math_Init - Archaic from AREA-I days
 implicit none
 integer iflag,openflag,readflag,LENGTH
 integer i,j,k,npts,stateindex,SWEEPS,IOUTSKIP
 real*8 readreal,xslope,yslope,zslope,idx
 real*8 sum,nominaltime,nominalstate(MAXX),rkalfa(4),krkbody(MAXX,4)
 real*4 tictotal,ticuser,ticsystem,toctotal,tocuser,tocsystem,elapsed(2)
 real*4 etime,aoa,unominal,wnominal,ALFAMAX
 real*8 zint,INITIALTIME,DELTATIME,FINALTIME
 real*8 CPUTIMESYSTEM,CPUTIMETOTAL,CPUTIMEUSER,rk4time,T
 character(128) FILEINPUTFILE
 type(BODYSTRUCTURE) BODY
 !CHARACTER(len=256)::AircraftDataFile !Aircraft-data file name (Moved to main Data structure - CM 8/16/2015
 CHARACTER(len=128)::rec,default
 !INTEGER:: numNodesNLL = 100 (Moved to TOWED Data structure - CM 8/16/2015)
 !REAL::Rde,CL0,Sw,bw,cw,Vflt,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,hx,hy,hz

 readflag = 0
 openflag = 0

 !call Math_Init !Initialize module - Functionality removed as of 1/2/2017
 
 ! Input File of Files

 call getarg(1,FILEINPUTFILE)
 if (len(trim(FILEINPUTFILE)) .lt. 1) then
    default = 'Input_Files/Parameters.Body'
    write(*,*) 'Using Default Parameter File = ',trim(default)
    FILEINPUTFILE = default
 endif
 
 BODY%INPUTFILE = trim(FILEINPUTFILE)
 call Load_Data(BODY)
 call Echo_Data(BODY)

 !Set the initial conditions of the quadcopter
 BODY%STATE(:) = 0
 BODY%STATE(3) = - 10. !Initial altitude

 !Run the RK4 Routine

 ! Define Constants
 rkalfa(1) = 1.0; rkalfa(2) = 2.0; rkalfa(3) = 2.0; rkalfa(4) = 1.0

 ! Initial State Vector
 call DERIVATIVES(T,BODY)
 
 ! Open Time Simulation Output Files

  open(unit=92,file='Output_Files/State.OUTPUT')
  rewind(92)

  !!HARD CODE RK4 STUFF
  FINALTIME = 100.0
  INITIALTIME = 0
  DELTATIME = 0.01D0
  IOUTSKIP = 1
  T = INITIALTIME 
 
  ! Integrate Equations of Motion
  
  tictotal = ETIME(elapsed)
  ticuser = elapsed(1)
  ticsystem = elapsed(2)
  npts = nint((FINALTIME-INITIALTIME)/DELTATIME)

  write(*,*) 'Simulation Start'

  ! THIS IS THE MAIN SIMULATION LOOP
  do i=1,npts

     ! Output Data to File
     !!OUTPUTS OUTPUTFILES STANDARD OUTS FILE OUTPUTS OUTPUT DATA
     !!OUTPUTS OUTPUTFILES STANDARD OUTS FILE OUTPUTS
     if (mod(i-1,IOUTSKIP) .eq. 0) then 
        !State.OUT File
        ! REVISIT subtracted 8 to get rid of thrust states
        write(92,fmt='(1000F30.10)',advance='no') T
        do j = 1,NUMSTATES
           write(92,fmt='(A1)',advance='no') ','
           write(92,fmt='(1000F30.10)',advance='no') BODY%STATE(j)
        end do
        ! do j = 1,NUMSTATES
        !    write(92,fmt='(A1)',advance='no') ','
        !    write(92,fmt='(1000F30.10)',advance='no') BODY%STATEDOT(j)
        ! end do
        do j = 1,9
           write(92,fmt='(A1)',advance='no') ','
           write(92,fmt='(1000F30.10)',advance='no') BODY%MUVEC(j,1)
        end do
        write(92,fmt='(A1)') ' '
     end if
     

     if (npts .gt. 100) then
        if (mod(i,npts/100) .eq. 1) then 
           idx = i
           write(*,*) 'Time Simulation: ',int(100*idx/npts)+1, ' Percent Complete'
        end if
     else
        write(*,*) 'Time Simulation: ',i, ' out of ',npts
     end if

     ! Store Nominal State Values

     nominaltime = T
     nominalstate(1:NUMSTATES) = BODY%STATE(1:NUMSTATES)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Numerical Integration of Equations of Motion

     do j=1,4

        ! State Values to Evaluate Derivatives

        if (j .ne. 1) then
           rk4time = nominaltime + DELTATIME/rkalfa(j)
           do k=1,NUMSTATES
              BODY%STATE(k) = nominalstate(k) + krkbody(k,j-1)/rkalfa(j)
           end do
        end if

        if (j .ne. 1) then
           rk4time = nominaltime + DELTATIME/rkalfa(j)
           do k=1,NUMSTATES
              BODY%STATE(k) = nominalstate(k) + krkbody(k,j-1)/rkalfa(j)
           end do
        end if

        ! Compute Derivatives
        
        call Derivatives(rk4time,BODY)

        ! Runge-Kutta Constants
        
        do k=1,NUMSTATES
           krkbody(k,j) = DELTATIME*BODY%STATEDOT(k)
        end do
    
     end do

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Step Time

     T = nominaltime + DELTATIME

     ! Step States
     
     do j=1,NUMSTATES
        sum = 0.0
        do k=1,4
           sum = sum + rkalfa(k)*krkbody(j,k)
        end do
        BODY%STATE(j) = nominalstate(j) + sum/6.0
     end do

     call Derivatives(T,BODY)

  end do

  toctotal = ETIME(elapsed)
  tocuser = elapsed(1)
  tocsystem = elapsed(2)

  CPUTIMEUSER = tocuser - ticuser
  CPUTIMESYSTEM = tocsystem - ticsystem
  CPUTIMETOTAL = toctotal - tictotal

  write(*,*) ' '
  write(*,*) 'SIMULATION CPU TIME USER (sec): ',CPUTIMEUSER,tocuser,ticuser
  write(*,*) 'SIMULATION CPU TIME SYSTEM (sec): ',CPUTIMESYSTEM,tocsystem,ticsystem
  write(*,*) 'SIMULATION CPU TIME TOTAL (sec): ',CPUTIMETOTAL,toctotal,tictotal
  write(*,*) ' '
 
  write(*,*) ' '
  write(*,*) 'TIME SIMULATION COMPLETE'
  write(*,*) ' '

  ! Close Output Files

  close(92)
  
end PROGRAM SIXDOF

SUBROUTINE Derivatives(T,BODY)
  use BODYDATATYPES
  implicit none
  integer i,iflag,ifind,openflag,readflag
  real*8 m,readreal,ctheta,stheta,sphi,cphi,spsi,cpsi,ttheta,dspeed,psiprev,omegar,T
  real*8 tend,slope,accel,noise,freq,vATM_I(3,1),vATM_A(3,1),uaero,vaero,waero,V_A,sumomega
  real*8 xcg,ycg,zcg,phi,theta,psi,ub,vb,wb,pb,qb,rb,rCF_B(3,1),forcevec(9,1),thrust
  real*8 Gammavec(3,1),kq,C_Ftether_I(3,1),C_Ftether_B(3,1),S_rCF_B(3,3),C_Mtether_B(3,1)
  real*8 xcgdot,ycgdot,zcgdot,phidot,thetadot,psidot,ubdot,vbdot,wbdot,c1,c2,c3,pbdot,qbdot,rbdot
  real*8 rReel_I(3,1),rCG_I(3,1),v_CG_I(3,1),S_wt_B(3,3),v_Reel_I(3,1),deti,yc,ycdot,kpy,kdy,ey,eydot
  real*8 xc,xcdot,kpx,kdx,exdot,ex,thetac,psicdot,kds,kps,kiz,kix,kiy
  real*8 thetacdot,epitch,epitchdot,angtheta,angpsi,eyaw,eyawdot,psic
  real*8 exbody, eybody, exbodydot, eybodydot, thrust_pusher,mupusher
  real*8 TVEC(9),TDOTVEC(9),TDBLDOTVEC(9),zc,zcdot,ez,ezdot,altitude,kp,kd,angphi,munominal
  real*8 sigmaF,omegaF,zetaF,C1F(9),C2F(9),C3F(9),idx,W2Tpwm(9,1),W0,j, psuedox, psuedoy 
  real*8 kpa,kda,phic,phicdot,eroll,erolldot,gammaz,gammazdot,gammax,gammaxdot,gammay,gammaydot
  character*256 xgridname,ygridname,zgridname
  character*1 letter
  character*10 number
  type(BODYSTRUCTURE) BODY

  !Extract everything from state vector
  xcg = BODY%STATE(1)
  ycg = BODY%STATE(2)
  zcg = BODY%STATE(3)
  phi = BODY%STATE(4) 
  theta = BODY%STATE(5)
  psi = BODY%STATE(6)
  ub = BODY%STATE(7)
  vb = BODY%STATE(8)
  wb = BODY%STATE(9)
  pb = BODY%STATE(10)
  qb = BODY%STATE(11) 
  rb = BODY%STATE(12)
  TVEC(1)    = BODY%STATE(13)
  TDOTVEC(1) = BODY%STATE(14)
  TVEC(2)    = BODY%STATE(15)
  TDOTVEC(2) = BODY%STATE(16)
  TVEC(3)    = BODY%STATE(17)
  TDOTVEC(3) = BODY%STATE(18)
  TVEC(4)    = BODY%STATE(19)
  TDOTVEC(4) = BODY%STATE(20)
  TVEC(5)    = BODY%STATE(21)!!!!!!!!!!!ADDED TO THIS TO TRY 8 ROTORS
  TDOTVEC(5) = BODY%STATE(22)
  TVEC(6)    = BODY%STATE(23)
  TDOTVEC(6) = BODY%STATE(24)
  TVEC(7)    = BODY%STATE(25)
  TDOTVEC(7) = BODY%STATE(26)
  TVEC(8)    = BODY%STATE(27)
  TDOTVEC(8) = BODY%STATE(28)
  TVEC(9) = BODY%STATE(29)
  TDOTVEC(9) = BODY%STATE(30)
  gammaz = BODY%STATE(31)
  gammay = BODY%STATE(32)
  gammax = BODY%STATE(33)
  ! Aircraft to Inertial Transformation Matrix - Same for BODYBODY
  ctheta = cos(theta);
  stheta = sin(theta);
  ttheta = stheta / ctheta;
  cphi = cos(phi);
  sphi = sin(phi);
  spsi = sin(psi);
  cpsi = cos(psi);
  BODY%TIC(1,1) = ctheta * cpsi;
  BODY%TIC(2,1) = ctheta * spsi;
  BODY%TIC(3,1) = -stheta;
  BODY%TIC(1,2) = sphi * stheta * cpsi - cphi * spsi;
  BODY%TIC(2,2) = sphi * stheta * spsi + cphi * cpsi;
  BODY%TIC(3,2) = sphi * ctheta;
  BODY%TIC(1,3) = cphi * stheta * cpsi + sphi * spsi;
  BODY%TIC(2,3) = cphi * stheta * spsi - sphi * cpsi;
  BODY%TIC(3,3) = cphi * ctheta;

  BODY%TBI(1,1) = cpsi;
  BODY%TBI(2,1) = spsi;
  BODY%TBI(1,2) = -spsi;
  BODY%TBI(2,2) = cpsi;
  ! Inertial to Aircraft Transformation Matrix - Same for BODYCOPTER

  BODY%TCI = transpose(BODY%TIC)

  ! State Derivatives - Kinematics 
  xcgdot = BODY%TIC(1,1)*ub + BODY%TIC(1,2)*vb + BODY%TIC(1,3)*wb
  ycgdot = BODY%TIC(2,1)*ub + BODY%TIC(2,2)*vb + BODY%TIC(2,3)*wb
  zcgdot = BODY%TIC(3,1)*ub + BODY%TIC(3,2)*vb + BODY%TIC(3,3)*wb  
  phidot = pb + sphi * ttheta * qb + cphi * ttheta * rb;
  thetadot = cphi * qb - sphi * rb;
  psidot = (sphi / ctheta) * qb + (cphi / ctheta) * rb;

  ! Gravity Forces and Moments - Same for Quadcopter

  BODY%FXGRAV = 0.0; BODY%FYGRAV = 0.0; BODY%FZGRAV = 0.0;
  BODY%MXGRAV = 0.0; BODY%MYGRAV = 0.0; BODY%MZGRAV = 0.0;
  if (BODY%GRAVOFFON .eq. 1) then
     BODY%FXGRAV = BODY%TIC(3,1)*BODY%WEIGHT
     BODY%FYGRAV = BODY%TIC(3,2)*BODY%WEIGHT
     BODY%FZGRAV = BODY%TIC(3,3)*BODY%WEIGHT
  end if

  ! Aerodynamic Forces and Moments - Still Same for Quadcopter

  BODY%FXAERO = 0.0; BODY%FYAERO = 0.0; BODY%FZAERO = 0.0;
  BODY%MXAERO = 0.0; BODY%MYAERO = 0.0; BODY%MZAERO = 0.0;

  if (BODY%AEROOFFON .eq. 1) then
     vATM_I(1,1) = BODY%VXWIND
     vATM_I(2,1) = BODY%VYWIND
     vATM_I(3,1) = BODY%VZWIND
     vATM_A = matmul(BODY%TCI,vATM_I)

     !Add in atmospheric winds

     uaero = ub - vATM_A(1,1)
     vaero = vb - vATM_A(2,1)
     waero = wb - vATM_A(3,1)

     !Compute total velocity

     V_A = sqrt(uaero**2 + vaero**2 + waero**2)

     !Quadcopter Aerodynamic Model written by Lisa Schibelius - 12/2016

     !Recompute KT
     BODY%KT = BODY%C_T*((BODY%DEN*qPI*(BODY%RNEW**4)/4))

     !Compute Thrust

     !sigmaF = 0.000437554764978899 !0.000437554764978899
     sigmaF = 0.00133 !This is based off of reading Nghia Huynh's draft journal
     omegaF = 45.42   !18.65
     zetaF  = 0.942     !0.8533

     !!!!!!!!!!!OUTER LOOP CONTROL!!!!!!!!!!!!!
     munominal = (1630.0 - 1100.0)/2+1100.0

     !because the body frame moves and we always want to fly foward 
     ! we need psuedo directions to calculate where the aircraft actually needs to fly
     !!!!!!!!!!!!!!!!roll control!!!!!!!!!!!!
     
     yc = -25.0
     ycdot = 0.0
     kpy = -0.04
     kdy = -0.05
     kiy = 0.001 
     ey = ycg - yc
     eydot = ycgdot - ycdot
    !!!!!!!!!!!!!!!!!pitch control!!!!!!!!!!!!!!!
     xc = 30.0
     xcdot = 0.0
     kpx = 0.035
     kdx = 0.05
     kix = 0.001
     ex = xcg - xc
     exdot = xcgdot - xcdot

     exbody = cpsi * ex + spsi *ey 
     eybody = -spsi * ex + cpsi *ey
     exbodydot = cpsi *exdot +spsi* eydot + psidot *(-spsi *ex+cpsi*ey) 
     eybodydot = -spsi *exdot + cpsi *eydot + psidot * (-cpsi * ex- spsi *ey)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!yaw control!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!yaw!!!!!!!!!!!!!!
    kps = -60.0
    kds = -120.0
    psicdot = 0
    eyawdot = psidot - psicdot
    !!!!!!!!!!!!set phi and theta to 0 and use psi to control Direction

    !!!Altitude Stuff
    zc = -50.0
    kp = 60.0
    kd = 200.0
    kiz = 2.0     
    zcdot = 0
    ez = zcg - zc
    gammazdot = ez
    ezdot = zcgdot - zcdot
    altitude = kp*ez  + kd*ezdot + kiz*gammaz

    mupusher = 1100.0

    !!PHASING
     if (BODY%PHASE == 1) then
        phic = 0
        thetac = 0
        psic = 0
        eyaw = psi - psic
        if (abs(ez) < 1) then
            BODY%PHASE = 2
            write(*,*) 'Switching to Phase 2'
        end if   
        gammaydot = 0
        gammaxdot = 0
     else if (BODY%PHASE == 2) then
        phic = 0
        thetac = 0        
        psic = atan(ey/ex)!this is in radians 
        eyaw = psi - psic   
        if (abs(eyaw) < 2.*PI/180.) then
            BODY%PHASE = 3
            write(*,*) 'Switching to Phase 3 = ',psi,eyaw,psic
            BODY%PSILOC = psi
            !PAUSE
        end if
        gammaydot = 0
        gammaxdot = 0
     else if (BODY%PHASE == 3) then
        phic = kpy*eybody + kdy*eybodydot + kiy*gammay        
        !phic = 0.0
        !thetac = kpx*exbody + kdx*exbodydot + kix*gammax
        thetac = 0.0
        !write(*,*) 'exbody = ',exbody
        mupusher = (-70./40.)*exbody + 1100.0
        !psic = BODY%PSILOC 
        if (abs(ex) .lt. 0.1) then
            psic = SIGN(PI/2,ex)
        else
            psic = atan(ey/ex)!this is in radians 
        end if
        !write(*,*) psic
        eyaw = psi - psic   
        gammaydot = eybody      
        gammaxdot = exbody 
        if ((sqrt(exbody**2 + eybody**2) < 1) .and. (sqrt(exbodydot**2 + eybodydot**2) < 0.1)) then
            BODY%PHASE = 4
            write(*,*) 'Switching to Phase 4',T
            BODY%PSILOC = psi
        end if
    else
        mupusher = 1100.0
        phic = 0
        thetac = 0
        psic = BODY%PSILOC
        eyaw = psi-psic
        gammaydot = 0
        gammaxdot = 0        
        zc = -1
        ez = zcg - zc
        ezdot = zcgdot
        altitude = 10*ezdot + 2*ez
        !write(*,*) zcgdot,ezdot,kd,altitude
        !PAUSE
     end if

    !!//Saturation Block
    if (abs(phic) .gt. 30*PI/180.0) then
        phic = SIGN(30*PI/180.0,phic)
    end if

    if (abs(thetac) .gt. 30*PI/180.0) then
        thetac = SIGN(30*PI/180.0,thetac)
     end if

     !!!!!!!!!!!!!!!altitude control!!!!!!!!!!

     if ((munominal + altitude) .gt. 1900.) then
        altitude = 1800.-munominal
        gammazdot = 0.
     endif  

     !write(*,*) theta*PI/180.0
     !!!!!!!INNER LOOP CONTROL!!!!!!!!!!!!!!!!
    !!!!!!!!!!!roll!!!!!!!!!!!!
     kpa = -60.0
     kda = -30.0
     !phic = 20.0*PI/180.0
     phicdot = 0
     !phic = 0 
     eroll = phi - phic
     erolldot = phidot - phicdot
     angphi = kpa*eroll + kda*erolldot

     !!!!!!!!!pitch!!!!!!!!!!!!!!
     thetacdot = 0
     !thetac = 0
     epitch = theta - thetac
     epitchdot = thetadot - thetacdot
     angtheta = kpa*epitch + kda*epitchdot   

     !!!Yaw     
     angpsi = kps*eyaw + kds*eyawdot
    
     !!! According to dynamic equations, a positive roll will have rotors 1,4,5,8 > 2,3,6,7. This was previously 2,3>1,4
     BODY%MUVEC(1,1) = munominal + altitude + angtheta + angphi + angpsi !cw
     BODY%MUVEC(5,1) = munominal + altitude + angtheta + angphi - angpsi !ccw
     BODY%MUVEC(2,1) = munominal + altitude + angtheta - angphi - angpsi !ccw
     BODY%MUVEC(6,1) = munominal + altitude + angtheta - angphi + angpsi !cw
     BODY%MUVEC(3,1) = munominal + altitude - angtheta - angphi + angpsi !cw
     BODY%MUVEC(7,1) = munominal + altitude - angtheta - angphi - angpsi !ccw
     BODY%MUVEC(4,1) = munominal + altitude - angtheta + angphi - angpsi !ccw
     BODY%MUVEC(8,1) = munominal + altitude - angtheta + angphi + angpsi !cw
     BODY%MUVEC(9,1) = mupusher

     do idx = 1,9
        if (BODY%MUVEC(idx,1) .lt. 1100.0) then
           BODY%MUVEC(idx,1) = 1100.0
        end if
        if (BODY%MUVEC(idx,1) .gt. 1900.0) then
           BODY%MUVEC(idx,1) = 1900.0
        end if
     end do

     !write(*,*) BODY%MUVEC

     !!! Second order filter
     do idx = 1,9
        C1F(idx) = -2*zetaF*TDOTVEC(idx)
        C2F(idx) = -(omegaF)*TVEC(idx)
        C3F(idx) = (BODY%MUVEC(idx,1)-1100.0)*sigmaF*omegaF   ! replaced sigma with force
        TDBLDOTVEC(idx) = omegaF*(C1F(idx) + C2F(idx) + C3F(idx))
        BODY%THRUSTVEC(idx,1) = TVEC(idx)
        if (BODY%THRUSTVEC(idx,1) .gt. 0) then
            BODY%OMEGAVEC(idx,1) = sqrt(BODY%THRUSTVEC(idx,1)/BODY%KT)
        else
            BODY%OMEGAVEC(idx,1) = 0.0
        end if
     end do


     sumomega = sum(BODY%OMEGAVEC)-BODY%OMEGAVEC(9,1)

     !!! Make sure angular velocities of rotor does not go beyond the limit
     if (sumomega .ge. BODY%OMEGAMAX*8) then
        do j = 1,9
           if (BODY%OMEGAVEC(j,1) .gt. BODY%OMEGAMAX) then
              BODY%OMEGAVEC(j,1) = BODY%OMEGAMAX
           end if
           if (BODY%OMEGAVEC(j,1) .lt. 0.00D0) then
              BODY%OMEGAVEC(j,1) = 0.00D0
           end if
        end do
        sumomega = sum(BODY%OMEGAVEC)-BODY%OMEGAVEC(9,1)
        BODY%THRUSTVEC = BODY%KT*BODY%OMEGAVEC**2
        do j = 1,9
           TVEC(idx) = BODY%THRUSTVEC(idx,1)
        end do
     end if
     forcevec = BODY%THRUSTVEC
     thrust = sum(BODY%THRUSTVEC)-BODY%THRUSTVEC(9,1)
     thrust_pusher = BODY%THRUSTVEC(9,1)

     !write(*,*) 'Rotor Thrust = ',forcevec,thrust

     !Aerodynamic Forces
     if (sumomega .gt. 1e-2) then
        BODY%FXAERO = -thrust*(((BODY%ALC/(sumomega*BODY%RNEW))+BODY%DXD)*uaero - ((BODY%ALS)/(sumomega*BODY%RNEW))*vaero)
        BODY%FYAERO = -thrust*(((BODY%ALS)/(sumomega*BODY%RNEW))*uaero + (((BODY%ALC)/(sumomega*BODY%RNEW))+BODY%DYD)*vaero)
        BODY%FZAERO = -thrust
     end if

     !!Add Pusher Thrust
     BODY%FXAERO = BODY%FXAERO + thrust_pusher

     !write(*,*) BODY%FZAERO,thrust,sumomega,BODY%OMEGAVEC(5,1),BODY%MUVEC(5,1),mupusher

     omegar = BODY%OMEGAVEC(1,1) - BODY%OMEGAVEC(5,1) - BODY%OMEGAVEC(2,1) + BODY%OMEGAVEC(6,1) + BODY%OMEGAVEC(3,1) - BODY%OMEGAVEC(7,1) + BODY%OMEGAVEC(8,1) - BODY%OMEGAVEC(4,1)
     Gammavec(1,1) = (BODY%IRR * omegar * qb)
     Gammavec(2,1) = (-BODY%IRR * omegar * pb)
     Gammavec(3,1) = 0
     !gotodynamics

     !!!!!!!!! Aerodynamics
     !bquad = BODY%C_TAU*((BODY%DEN*qPI*(BODY%RNEW**5)/4))
     kq = BODY%DEN*PI*(BODY%RNEW**5)*BODY%C_Q !this is the lumped balde power parameter 
     !!! According to dynamic equations, a positive roll will have rotors 1,4 > 2,3. This was previously 2,3>1,4
     ! Since T3 = T1*Ltheta_front/Ltheta_back we're just going to do this for simplicity
     BODY%MXAERO = Gammavec(1,1) + (BODY%LPHI12*(TVEC(1) +TVEC(5) - TVEC(2) - TVEC(6)) + BODY%LPHI34*(TVEC(4)+TVEC(8) - TVEC(3)-TVEC(7)))
     BODY%MYAERO = Gammavec(2,1) + BODY%LTHETA12*(TVEC(1) +TVEC(5) +TVEC(2) + TVEC(6) - TVEC(3) -TVEC(4) - TVEC(7) - TVEC(8))
     BODY%MZAERO = Gammavec(3,1) + kq*(BODY%OMEGAVEC(1,1)**2 - BODY%OMEGAVEC(5,1)**2 - BODY%OMEGAVEC(2,1)**2 + BODY%OMEGAVEC(6,1)**2+ BODY%OMEGAVEC(3,1)**2 - BODY%OMEGAVEC(7,1)**2 - BODY%OMEGAVEC(4,1)**2 + BODY%OMEGAVEC(8,1)**2)
  else
     BODY%FXAERO = 0
     BODY%FYAERO = 0
     BODY%FZAERO = 0
     BODY%MXAERO = 0
     BODY%MYAERO = 0
     BODY%MZAERO = 0
  end if
     

  ! Total Forces and Moments

  BODY%FXTOTAL = BODY%FXGRAV + BODY%FXAERO
  BODY%FYTOTAL = BODY%FYGRAV + BODY%FYAERO
  BODY%FZTOTAL = BODY%FZGRAV + BODY%FZAERO
  BODY%MXTOTAL = BODY%MXGRAV + BODY%MXAERO
  BODY%MYTOTAL = BODY%MYGRAV + BODY%MYAERO
  BODY%MZTOTAL = BODY%MZGRAV + BODY%MZAERO


  ubdot = BODY%FXTOTAL/BODY%MASS + rb*vb - qb*wb
  vbdot = BODY%FYTOTAL/BODY%MASS + pb*wb - rb*ub
  !write(*,*) vbdot,BODY%FYTOTAL,BODY%FYGRAV,BODY%FYAERO
  wbdot = BODY%FZTOTAL/BODY%MASS + qb*ub - pb*vb

  c1 = BODY%MXTOTAL - pb*(qb*BODY%IXZ-rb*BODY%IXY) - qb*(qb*BODY%IYZ-rb*BODY%IYY) - rb*(qb*BODY%IZZ-rb*BODY%IYZ)
  c2 = BODY%MYTOTAL - pb*(rb*BODY%IXX-pb*BODY%IXZ) - qb*(rb*BODY%IXY-pb*BODY%IYZ) - rb*(rb*BODY%IXZ-pb*BODY%IZZ)
  c3 = BODY%MZTOTAL - pb*(pb*BODY%IXY-qb*BODY%IXX) - qb*(pb*BODY%IYY-qb*BODY%IXY) - rb*(pb*BODY%IYZ-qb*BODY%IXZ)
  pbdot = BODY%IXXI*c1 + BODY%IXYI*c2 + BODY%IXZI*c3
  qbdot = BODY%IXYI*c1 + BODY%IYYI*c2 + BODY%IYZI*c3
  rbdot = BODY%IXZI*c1 + BODY%IYZI*c2 + BODY%IZZI*c3

  ! Wrap State Derivatives

  BODY%STATEDOT(1) = xcgdot
  BODY%STATEDOT(2) = ycgdot
  BODY%STATEDOT(2) = ycgdot
  BODY%STATEDOT(3) = zcgdot
  BODY%STATEDOT(4) = phidot
  BODY%STATEDOT(5) = thetadot
  BODY%STATEDOT(6) = psidot
  BODY%STATEDOT(7) = ubdot
  BODY%STATEDOT(8) = vbdot
  BODY%STATEDOT(9) = wbdot
  BODY%STATEDOT(10) = pbdot 
  BODY%STATEDOT(11) = qbdot 
  BODY%STATEDOT(12) = rbdot
  BODY%STATEDOT(13) = TDOTVEC(1)
  BODY%STATEDOT(14) = TDBLDOTVEC(1)
  BODY%STATEDOT(15) = TDOTVEC(2)
  BODY%STATEDOT(16) = TDBLDOTVEC(2)
  BODY%STATEDOT(17) = TDOTVEC(3)
  BODY%STATEDOT(18) = TDBLDOTVEC(3)
  BODY%STATEDOT(19) = TDOTVEC(4)
  BODY%STATEDOT(20) = TDBLDOTVEC(4)
  BODY%STATEDOT(21) = TDOTVEC(5)
  BODY%STATEDOT(22) = TDBLDOTVEC(5)
  BODY%STATEDOT(23) = TDOTVEC(6)
  BODY%STATEDOT(24) = TDBLDOTVEC(6)
  BODY%STATEDOT(25) = TDOTVEC(7)
  BODY%STATEDOT(26) = TDBLDOTVEC(7)
  BODY%STATEDOT(27) = TDOTVEC(8)
  BODY%STATEDOT(28) = TDBLDOTVEC(8)
  BODY%STATEDOT(29) = TDOTVEC(9)
  BODY%STATEDOT(30) = TDBLDOTVEC(9)
  BODY%STATEDOT(31) = gammazdot 
  BODY%STATEDOT(32) = gammaydot 
  BODY%STATEDOT(33) = gammaxdot

  ! write(*,fmt='(4e18.2)') BODY%STATE
  ! write(*,fmt='(4e18.2)') BODY%STATEDOT
  ! PAUSE

end SUBROUTINE Derivatives

SUBROUTINE Echo_Data(BODY)
  use BODYDATATYPES
  implicit none
  integer i,iflag,ifind,openflag,readflag
  real*8 m,readreal,ctheta,stheta,sphi,cphi,spsi,cpsi,ttheta,dspeed,psiprev,omegar
  real*8 tend,slope,accel,noise,freq,vATM_I(3,1),vATM_A(3,1),uaero,vaero,waero,V_A,sumomega
  real*8 xcg,ycg,zcg,phi,theta,psi,ub,vb,wb,pb,qb,rb,rCF_B(3,1),forcevec(9,1),thrust
  real*8 Gammavec(3,1),kq,C_Ftether_I(3,1),C_Ftether_B(3,1),S_rCF_B(3,3),C_Mtether_B(3,1)
  real*8 xcgdot,ycgdot,zcgdot,phidot,thetadot,psidot,ubdot,vbdot,wbdot,c1,c2,c3,pbdot,qbdot,rbdot
  real*8 rReel_I(3,1),rCG_I(3,1),v_CG_I(3,1),S_wt_B(3,3),v_Reel_I(3,1),deti
  real*8 TVEC(4),TDOTVEC(4),TDBLDOTVEC(4)
  real*8 sigmaF,omegaF,zetaF,C1F(4),C2F(4),C3F(4),idx,W2Tpwm(4,1),W0,j
  character*256 xgridname,ygridname,zgridname
  character*1 letter
  character*10 number
  type(BODYSTRUCTURE) BODY

  write(*,*) 'Echoing Data...'

  open(unit=25,file='Output_Files/Parameters.Check',iostat=openflag)

  write(25,*) 'THIS FILE IS WRONG'
  write(25,*) 'Body Input File: '
  write(25,*) trim(BODY%INPUTFILE)
  write(25,*) ' '
  write(25,*) 'Module off or on: ',BODY%OFFON
  write(25,*) 'Model (0=Integration,1=Constant, 2=Table,3=Quad 6DOF Model): ',BODY%MODNO
  if (BODY%MODNO .eq. 3) then
     !Output quadcopter stuff
     write(25,*) 'Gravity Flag (0=Off, 1=On): ',BODY%GRAVOFFON
     write(25,*) 'Aerodynamics Flag (0=Off, 1=On): ',BODY%AEROOFFON
     write(25,*) 'Contact Flag (0=Off, 1=On): ',BODY%CONTOFFON
     write(25,*) 'Mass (kg): ',BODY%MASS
     write(25,*) 'Weight (N): ',BODY%WEIGHT
     write(25,*) 'Stationline of Mass Center (m): ',BODY%SLCG
     write(25,*) 'Buttline of Mass Center (m): ',BODY%BLCG
     write(25,*) 'Waterline of Mass Center (m): ',BODY%WLCG
     write(25,*) 'Ixx (kg m^2): ',BODY%IXX
     write(25,*) 'Iyy (kg m^2): ',BODY%IYY
     write(25,*) 'Izz (kg m^2): ',BODY%IZZ
     write(25,*) 'Ixy (kg m^2): ',BODY%IXY
     write(25,*) 'Ixz (kg m^2): ',BODY%IXZ
     write(25,*) 'Iyz (kg m^2): ',BODY%IYZ
     write(25,*) 'Ixx Inverse (1/(kg m^2)): ',BODY%IXXI
     write(25,*) 'Iyy Inverse (1/(kg m^2)): ',BODY%IYYI
     write(25,*) 'Izz Inverse (1/(kg m^2)): ',BODY%IZZI
     write(25,*) 'Ixy Inverse (1/(kg m^2)): ',BODY%IXYI
     write(25,*) 'Ixz Inverse (1/(kg m^2)): ',BODY%IXZI
     write(25,*) 'Iyz Inverse (1/(kg m^2)): ',BODY%IYZI
     write(25,*) 'Turn Radius of AC(m): ',  BODY%TURNRADIUS
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%ALC
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%ALS
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%DXD
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%DYD
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%RNEW
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%C_T
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%C_Q
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%LPHI12
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%LPHI34
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%LTHETA12
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%LTHETA34
     write(25,*) 'Quadcopter Aero Parameter: ', BODY%OMEGAMAX
  else
     write(25,*) ' '
     write(25,*) 'Stationline of Tether Reel Point on Copter: ', BODY%SLREEL
     write(25,*) 'Buttline of Tether Reel Point on Copter: ', BODY%BLREEL
     write(25,*) 'Waterline of Tether Reel Point on Copter: ', BODY%WLREEL
     write(25,*) 'Stationline of Airwake grid start on Copter: ', BODY%SLAIRWAKE
     write(25,*) 'Buttline of Airwake grid start on Copter: ', BODY%BLAIRWAKE
     write(25,*) 'Waterline of Airwake grid start on Copter: ', BODY%WLAIRWAKE
     if ((BODY%MODNO .eq. 1) .or. (BODY%MODNO .eq. 0 )) then
        write(25,*) 'Copter Speed from .COPTER File(ft/s): ',BODY%FINALSPEED
        write(25,*) 'Restart Speed from RESTART File (ft/s): ',BODY%RESTARTSPEED
        write(25,*) 'Copter Azimuthal Direction (deg): ',57.3*BODY%PSI
        write(25,*) 'Copter Initial X (ft): ',BODY%XCGINITIAL
        write(25,*) 'Copter Initial Y (ft): ',BODY%YCGINITIAL
        write(25,*) 'Copter Initial Z (ft): ',BODY%ZCGINITIAL
        write(25,*) 'Copter Noise X (ft/s^2): ',BODY%XDDOTNOISE
        write(25,*) 'Copter Noise Y (ft/s^2): ',BODY%YDDOTSCALE
        write(25,*) 'Copter Noise Z (ft/s^2): ',BODY%YDDOTPERIOD
        write(25,*) ' '
     end if
     if (BODY%MODNO .eq. 2) then
        write(25,*) 'Time (s),      Xcg (ft),      Ycg (ft),      Zcg (ft)'
        write(25,*) '----------------------------------------------------'
        do i=1,BODY%TABSIZE  
           write(25,fmt='(4e18.8)') BODY%TIMETAB(i),BODY%XCGTAB(i),BODY%YCGTAB(i),BODY%ZCGTAB(i)
        end do
        write(25,*) ' '
        write(25,*) 'Time (s),    Phi (deg),  Theta (deg),    Psi (deg)'
        write(25,*) '----------------------------------------------------'
        do i=1,BODY%TABSIZE  
           write(25,fmt='(4e18.8)') BODY%TIMETAB(i),57.3*BODY%PHITAB(i),57.3*BODY%THETATAB(i),57.3*BODY%PSITAB(i)
        end do
        write(25,*) ' '
        write(25,*) 'Time (s),     Ub (ft/s),     Vb (ft/s),     Wb (ft/s)'
        write(25,*) '----------------------------------------------------'
        do i=1,BODY%TABSIZE   
           write(25,fmt='(4e18.8)') BODY%TIMETAB(i),BODY%UBTAB(i),BODY%VBTAB(i),BODY%WBTAB(i)
        end do
        write(25,*) ' '
        write(25,*) 'Time (s),     Pb (r/s),      Qb (r/s),    Rb (r/s)'
        write(25,*) '----------------------------------------------------'
        do i=1,BODY%TABSIZE  
           write(25,fmt='(4e18.8)') BODY%TIMETAB(i),BODY%PBTAB(i),BODY%QBTAB(i),BODY%RBTAB(i)
        end do
        write(25,*) ' '
     end if
  end if
  write(25,*) 'Data Quality Flag (nd, 0=Data Not Loaded Successfully, 1=Data Loaded Successfully): ',BODY%DQFLAG

  write(*,*) 'Echoing Complete'
  
  close(25)
  
end SUBROUTINE Echo_Data

SUBROUTINE Load_Data(BODY)
  use BODYDATATYPES
  implicit none
  integer i,iflag,ifind,openflag,readflag
  real*8 m,readreal,ctheta,stheta,sphi,cphi,spsi,cpsi,ttheta,dspeed,psiprev,omegar
  real*8 tend,slope,accel,noise,freq,vATM_I(3,1),vATM_A(3,1),uaero,vaero,waero,V_A,sumomega
  real*8 xcg,ycg,zcg,phi,theta,psi,ub,vb,wb,pb,qb,rb,rCF_B(3,1),forcevec(9,1),thrust
  real*8 Gammavec(3,1),kq,C_Ftether_I(3,1),C_Ftether_B(3,1),S_rCF_B(3,3),C_Mtether_B(3,1)
  real*8 xcgdot,ycgdot,zcgdot,phidot,thetadot,psidot,ubdot,vbdot,wbdot,c1,c2,c3,pbdot,qbdot,rbdot
  real*8 rReel_I(3,1),rCG_I(3,1),v_CG_I(3,1),S_wt_B(3,3),v_Reel_I(3,1),deti
  real*8 TVEC(4),TDOTVEC(4),TDBLDOTVEC(4)
  real*8 sigmaF,omegaF,zetaF,C1F(4),C2F(4),C3F(4),idx,W2Tpwm(4,1),W0,j
  character*256 xgridname,ygridname,zgridname
  character*1 letter
  character*10 number
  type(BODYSTRUCTURE) BODY

  open(unit=94,file=trim(BODY%INPUTFILE),status='old',iostat=openflag)
  if (openflag .ne. 0) then
     write(*,*) BODY%INPUTFILE
     write(*,*) 'Error Opening Body Input File => ',trim(BODY%INPUTFILE),' <= ';STOP
  else
     write(*,*) 'Successfully Opened Body Input File =>',trim(BODY%INPUTFILE)
  end if
   
  rewind(94)
  read(unit=94,fmt=*,iostat=readflag) readreal; BODY%OFFON = readreal
  read(unit=94,fmt=*,iostat=readflag) readreal; BODY%MODNO = readreal

  read(unit=94,fmt=*,iostat=readflag) BODY%GRAVOFFON
  read(unit=94,fmt=*,iostat=readflag) BODY%AEROOFFON
  read(unit=94,fmt=*,iostat=readflag) BODY%CONTOFFON
  read(unit=94,fmt=*,iostat=readflag) BODY%WEIGHT
  read(unit=94,fmt=*,iostat=readflag) BODY%GRAVITY
  read(unit=94,fmt=*,iostat=readflag) BODY%SLCG
  read(unit=94,fmt=*,iostat=readflag) BODY%BLCG
  read(unit=94,fmt=*,iostat=readflag) BODY%WLCG
  read(unit=94,fmt=*,iostat=readflag) BODY%SLREEL
  read(unit=94,fmt=*,iostat=readflag) BODY%BLREEL
  read(unit=94,fmt=*,iostat=readflag) BODY%WLREEL
  read(unit=94,fmt=*,iostat=readflag) BODY%IXX
  read(unit=94,fmt=*,iostat=readflag) BODY%IYY
  read(unit=94,fmt=*,iostat=readflag) BODY%IZZ
  read(unit=94,fmt=*,iostat=readflag) BODY%IXY
  read(unit=94,fmt=*,iostat=readflag) BODY%IXZ
  read(unit=94,fmt=*,iostat=readflag) BODY%IYZ
  read(unit=94,fmt=*,iostat=readflag) BODY%TURNRADIUS
  read(unit=94,fmt=*,iostat=readflag) BODY%ALC
  read(unit=94,fmt=*,iostat=readflag) BODY%ALS
  read(unit=94,fmt=*,iostat=readflag) BODY%DXD
  read(unit=94,fmt=*,iostat=readflag) BODY%DYD
  read(unit=94,fmt=*,iostat=readflag) BODY%RNEW
  read(unit=94,fmt=*,iostat=readflag) BODY%C_T
  read(unit=94,fmt=*,iostat=readflag) BODY%C_Q
  read(unit=94,fmt=*,iostat=readflag) BODY%LPHI12
  read(unit=94,fmt=*,iostat=readflag) BODY%LPHI34
  read(unit=94,fmt=*,iostat=readflag) BODY%LTHETA12
  read(unit=94,fmt=*,iostat=readflag) BODY%LTHETA34
  read(unit=94,fmt=*,iostat=readflag) BODY%OMEGAMAX
  read(unit=94,fmt=*,iostat=readflag) BODY%IRR
  read(unit=94,fmt=*,iostat=readflag) BODY%MS_MIN
  read(unit=94,fmt=*,iostat=readflag) BODY%MS_MAX
  read(unit=94,fmt=*,iostat=readflag) readreal; BODY%CONTROLOFFON = int(readreal)
  read(unit=94,fmt=*,iostat=readflag) BODY%KPXBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KIXBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KDXBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KPYBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KIYBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KDYBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KPZBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KIZBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KDZBODY
  read(unit=94,fmt=*,iostat=readflag) BODY%KPPHI
  read(unit=94,fmt=*,iostat=readflag) BODY%KIPHI
  read(unit=94,fmt=*,iostat=readflag) BODY%KDPHI
  read(unit=94,fmt=*,iostat=readflag) BODY%KPTHETA
  read(unit=94,fmt=*,iostat=readflag) BODY%KITHETA
  read(unit=94,fmt=*,iostat=readflag) BODY%KDTHETA
  read(unit=94,fmt=*,iostat=readflag) BODY%KPPSI
  read(unit=94,fmt=*,iostat=readflag) BODY%KIPSI
  read(unit=94,fmt=*,iostat=readflag) BODY%KDPSI
  read(unit=94,fmt=*,iostat=readflag) BODY%XINTEGRAL 
  read(unit=94,fmt=*,iostat=readflag) BODY%YINTEGRAL 
  read(unit=94,fmt=*,iostat=readflag) BODY%ZINTEGRAL
  read(unit=94,fmt=*,iostat=readflag) BODY%PHIINTEGRAL 
  read(unit=94,fmt=*,iostat=readflag) BODY%THETAINTEGRAL
  read(unit=94,fmt=*,iostat=readflag) BODY%PSIINTEGRAL
  read(unit=94,fmt=*,iostat=readflag) BODY%XCOMMAND
  read(unit=94,fmt=*,iostat=readflag) BODY%YCOMMAND
  read(unit=94,fmt=*,iostat=readflag) BODY%ZCOMMAND
  read(unit=94,fmt=*,iostat=readflag) BODY%WAYPOINT

!!!DO SOME CALCULATIONS ON Quadcopter
  BODY%MASS = BODY%WEIGHT/BODY%GRAVITY 
  deti = + BODY%IXX*(BODY%IYY*BODY%IZZ-BODY%IYZ*BODY%IYZ) - BODY%IXY*(BODY%IXY*BODY%IZZ-BODY%IYZ*BODY%IXZ) + BODY%IXZ*(BODY%IXY*BODY%IYZ-BODY%IYY*BODY%IXZ)
  BODY%IXXI = (BODY%IYY*BODY%IZZ-BODY%IYZ*BODY%IYZ)/deti
  BODY%IXYI = (BODY%IYZ*BODY%IXZ-BODY%IXY*BODY%IZZ)/deti
  BODY%IXZI = (BODY%IXY*BODY%IYZ-BODY%IYY*BODY%IXZ)/deti
  BODY%IYYI = (BODY%IXX*BODY%IZZ-BODY%IXZ*BODY%IXZ)/deti
  BODY%IYZI = (BODY%IXY*BODY%IXZ-BODY%IXX*BODY%IYZ)/deti
  BODY%IZZI = (BODY%IXX*BODY%IYY-BODY%IXY*BODY%IXY)/deti

  close(94) 
  write(*,*) 'Body Load Complete'

  BODY%DQFLAG = 1

  RETURN
   
END SUBROUTINE LOAD_DATA

