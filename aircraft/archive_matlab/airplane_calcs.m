purge


%%%Center of mass calcs

m_structure = 320; %%grams
x_structure = 13.5; %%Inches

%%%%Servos
m_servo = 14; %%grams
x_servo1 = 8; %%%aileron
x_servo2 = 20; %%%elev
x_servo3 = 0.0; %%%rudder -- removed

%%%ESC
m_esc = 26; %%grams
x_esc = 3.0;

%%%Motor 
m_motor = 14+80; %%%66.6 grams to kilos to ounces
x_motor = -0.5;

%%%Reveiver
m_rec = 12; %%%grams
x_rec = 10.5;

%%%Arduino MEGA Flight Controller
m_ard = 80; %%grams
x_ard = 11.5; %%inches

%%%%2-Cell MEGA
m_2cell = 25;
x_2cell = 15.5;

%%%%Battery - 3 cell - 11.1 V
m_battery = 55; 
x_battery = 4.5;

%%%%Total weight
m_total = m_ard + m_2cell + m_structure + 2*m_servo + m_esc + m_motor + m_battery + m_rec

%%%%Center of mass all up
x_cg = (m_rec*x_rec + m_ard*x_ard + m_2cell*x_2cell + m_structure*x_structure + m_servo*(x_servo1 + x_servo2 + x_servo3) + m_esc*x_esc + m_motor*x_motor + m_battery*x_battery)/m_total

%%%Center of pressure
chord = 12.5;
cop_main_wing = 5.75 + chord/4;
chord_tail = 10;
cop_tail = 27 + chord_tail/4;

%%%Area of wing and tail
wingspan = 21.5*2; %%Inches
tip_chord = 6; %%%Inches
S_wing = 0.5*(chord+tip_chord)*wingspan; %%inches

wingspan_tail = 9*2; %%%inches
tip_chord_tail = 5; %%inches
S_tail = 0.5*(chord+tip_chord_tail)*wingspan_tail;

%%%TOtal Area
S_ft = (S_wing + S_tail)/144

%%%Center of pressure of entire airplane
x_cop = (cop_main_wing*S_wing + cop_tail*S_tail)/(S_wing+S_tail)

%%%Wing Loading
WL = 0.035274*m_total/ (S_wing+S_tail)*144

%%%OUNCES TO GRAMS - 1 gram = 0.035274 ounces

%%% THESE ARE IN UNITS OF OUNCES TO SQ FT
%Wing loading ballpark figures
%Slow flyers and thermal gliders – under 4 - 
%Trainers, park flyers, 3D – 5 to 7
%General sport and scale aerobatics – 7 to 10
%Sport and scale models – 10 to 13
%Warbirds and racers – 13 and over

%%%Thrust to weight
TW = 1.2;

Thrust_grams = TW*m_total

Thrust_lbs = Thrust_grams/1000*2.2

%%%Weight of Vehicle
Weight = m_total/1000*2.2

%%%Trim Velocity
Vtrim = 45; %%From pitot probe

%%Density
rho = 0.00238; %%%slugs/ft^3

CLtrim = 2*Weight/(rho*Vtrim^2*S_ft)

alfa_max = 12*pi/180
CLmax = 2*pi*alfa_max*0.8
Vstall = sqrt(2*Weight/(rho*S_ft*CLmax))

Vstall_ms = Vstall / 3.28