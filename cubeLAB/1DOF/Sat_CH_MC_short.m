clear
clc
close all

A = tf([1.3],[1,1.3])
G = tf([1],[.00303,0,0])
H = tf([1000],[1,1000])

GOL = A*G*H;
bode(GOL)
margin(GOL)
figure()
rlocus(GOL)

 C1 = tf([1,1],[1,650])
 C2 = tf([1,5],[1,300])
 C = C1*C2;
 CGOL = C*GOL

bode(CGOL)
margin(CGOL)

figure()
rlocus(CGOL)

K = 5;
GCL = K*CGOL/(1+K*CGOL);

figure()
step(GCL)
figure()
pzmap(GCL)


