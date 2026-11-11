import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# --- 1. System Parameters ---
mass = 4.5       # Mass in kg
length = 0.3     # Dimension along X in meters (30 cm)
width = 0.1      # Dimension along Y in meters (10 cm)

# Moment of Inertia for a uniform rectangular plate about its center
I = (1.0 / 12.0) * mass * (length**2 + width**2)
print(f"Calculated Moment of Inertia (I): {I:.5f} kg*m^2")

# --- 2. Initial Conditions & Simulation Time ---
theta_0 = 0.0                   # Initial roll angle (radians)
omega_0 = np.radians(10.0)      # Initial angular velocity (10 deg/s converted to rad/s)
initial_state = [theta_0, omega_0]

# --- 2. Simulation & Fixed-Step Settings ---
t_start = 0.0
t_end = 10.0
dt = 0.01                        # Fixed time step size (10 ms / 100 Hz loop)
t = np.arange(t_start, t_end + dt, dt)
n_steps = len(t)

# --- 3. Controller Parameters (PID) ---
Kp = 0.5                        # Proportional gain
Kd = 0.1                        # Derivative gain
theta_target = 0*45*np.pi/180.0             # Desired target roll angle (radians)
RCS = False
RCSTIME = -1e20

# --- 4. Differential Equation Model ---
def cubesat_dynamics(theta,omega, t, I, Kp, Kd, theta_target):
    global current,angular_acceleration_RW, RCS, RCSTIME,Thrust
    #theta, omega = state
    
    # PID Control Torque calculation inside the model
    error = theta - theta_target
    error_rate = omega - 0.0  # target velocity is 0
    tau_desired = -(Kp * error + Kd * error_rate)
    
    """
    ###Let's Assume MagTs
    #Maximum dipole Moment = 1000 A-m^2
    #Actual torque = n*current*A * magnetic field
    beta = 30000e-9 #T
    A = length*width
    max_current = 1.0 #Amps
    n = 1000.0/(max_current*A)
    max_torque = n*max_current*A * beta
    if abs(tau_desired) > max_torque:
        tau_control = np.sign(tau_desired)*max_torque
    else:
        tau_control = tau_desired
    current = tau_control/(beta*n*A)
    """
    
    
    #Let's do this for Reaction wheels
    max_torque_RW = 0.004
    if abs(tau_desired) > max_torque_RW:
        tau_control = np.sign(tau_desired)*max_torque_RW
    else:
        tau_control = tau_desired
    mRW = 0.13
    RRW = 0.042
    IRW = mRW*RRW**2 /2.0
    angular_acceleration_RW = tau_control / IRW

    """    
    ##Alright let's do Thrusters
    if RCS == False and (t-RCSTIME) > 2/1000.0:
        #Turn it on
        Thrust = 25/1000.0*np.sign(tau_desired)
        RCSTIME = t
    if RCS == True and (t-RCSTIME) > 2/1000.0:
        RCS = False
        Thrust = 0.0
        RCSTIME = t
    tau_control = Thrust*length
    """
    
    # Equations of Motion
    dtheta_dt = omega
    domega_dt = tau_control / I
    
    return dtheta_dt, domega_dt

# --- 5. Run Simulation ---
# Storage arrays for fixed-step simulation
theta_sim = np.zeros(n_steps)
omega_sim = np.zeros(n_steps)
tau_control_sim = np.zeros(n_steps)
omega_sim[0] = omega_0
# --- 5. Fixed-Step RK4 Integration Loop ---
for i in range(1, n_steps):
    th = theta_sim[i-1]
    om = omega_sim[i-1]
    
    # Runge-Kutta 4th Order (RK4) steps
    # k1
    dth1, dom1 = cubesat_dynamics(th, om,t[i-1],I, Kp, Kd, theta_target)
    
    # k2
    th2 = th + 0.5 * dt * dth1
    om2 = om + 0.5 * dt * dom1
    dth2, dom2 = cubesat_dynamics(th2, om2,t[i-1],I, Kp, Kd, theta_target)
    
    # k3
    th3 = th + 0.5 * dt * dth2
    om3 = om + 0.5 * dt * dom2
    dth3, dom3 = cubesat_dynamics(th3, om3,t[i-1],I, Kp, Kd, theta_target)
    
    # k4
    th4 = th + dt * dth3
    om4 = om + dt * dom3
    dth4, dom4 = cubesat_dynamics(th4, om4,t[i-1],I, Kp, Kd, theta_target)
    
    # Update states using weighted average of increments
    theta_sim[i] = th + (dt / 6.0) * (dth1 + 2*dth2 + 2*dth3 + dth4)
    omega_sim[i] = om + (dt / 6.0) * (dom1 + 2*dom2 + 2*dom3 + dom4)

# Post-calculate tau_control history for plotting
current_sim = 0*theta_sim
alpha_RW = 0*theta_sim
omega_RW = 0*theta_sim
thrust_sim = 0*theta_sim
RCS = False
RCSTIME = -1e20
for i in range(0,len(theta_sim)):
    cubesat_dynamics(theta_sim[i],omega_sim[i], t[i], I, Kp, Kd, theta_target)

    ##MagTs
    #current_sim[i] = current

    ##RWs
    alpha_RW[i] = angular_acceleration_RW
    #Let's also use Euler's method to simulate the RW angular velocity
    omega_RW[i] = omega_RW[i-1] + dt * angular_acceleration_RW

    ##RCS
    #thrust_sim[i] = Thrust

# --- 6. Plot Results (3 Subplots) ---
plt.figure(figsize=(10, 10))

# Subplot 1: Roll Angle
plt.subplot(3, 1, 1)
plt.plot(t, np.degrees(theta_sim), label='Roll Angle (deg)', color='blue', linewidth=2)
plt.axhline(0, color='black', linestyle='--', alpha=0.5, label='Target Angle')
plt.ylabel('Roll Angle (deg)')
plt.title('2D CubeSat Attitude Stabilization & Reaction Wheel Torque')
plt.grid(True)
plt.legend(loc='upper right')

# Subplot 2: Angular Velocity
plt.subplot(3, 1, 2)
plt.plot(t, np.degrees(omega_sim), label='Angular Velocity (deg/s)', color='red', linewidth=2)
plt.axhline(0, color='black', linestyle='--', alpha=0.5)
plt.ylabel('Angular Velocity (deg/s)')
plt.grid(True)
plt.legend(loc='upper right')

# Subplot 3: Control Value
plt.subplot(3, 1, 3)
plt.plot(t, alpha_RW, label='Angular Acceleration(rad/s^2)', color='green', linewidth=2)
plt.plot(t, omega_RW, label='Angular Velocity RW(rad/s)', color='orange', linewidth=2)
plt.axhline(0, color='black', linestyle='--', alpha=0.5)
plt.xlabel('Time (s)')
#plt.ylabel('Current (Amps)')
plt.ylabel('Reaction Wheels')
#plt.ylabel('Thrust RCS (N)')
plt.grid(True)
plt.legend(loc='upper right')

plt.tight_layout()
plt.show()
