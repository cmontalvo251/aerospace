function breeze=RBFhack1(P,xknot,alfcough);
NNN=size(alfcough);
Nknot=NNN(2);
breeze=0;
for Jay=1:Nknot
Pwork=(P(1)-xknot(Jay,1))^2+(P(2)-xknot(Jay,2))^2+(P(3)-xknot(Jay,3))^2;
Pwork=sqrt(Pwork);
breeze=breeze+alfcough(Jay)*Pwork;
end;
