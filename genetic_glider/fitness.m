function J = fitness(x)

%x represents the state vector of the panels
%x = [(l1,w1,d1,g1,l2,w2,d2,g2,length)];
[m,I,xbody] = convertstate(x)
%xbody = [panel1,panel2,panel3,panel4,W,I] %distances in body frame
%panel1 = [x1,y1,z1,phi1,theta1,psi1]
%note panel1 = panel3 just on the opposite side

J = m;