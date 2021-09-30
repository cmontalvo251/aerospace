clear
clc
close all

A = tf([1.316],[1,1.316])
G = tf([1],[.00303,0,0])
H = tf([1000],[1,1000])

bode(G)
margin(G)
GOL = A*G*H;
bode(GOL)
margin(GOL)
figure()
rlocus(GOL)

C1 = tf([1,1],[1,500])
C2 = tf([1,10],[1,200])
C = C1*C2;
CGOL = C*GOL
pzmap(CGOL)
bode(CGOL)
margin(CGOL)

figure()
rlocus(CGOL)
