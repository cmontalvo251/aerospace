function breeze=RBFhack2(P,xknot,alfcough,ikind);
NNN=size(alfcough);
Nknot=NNN(2);
breeze=0;
for Jay=1:Nknot
Pwork=(P(1)-xknot(Jay,1))^2+(P(2)-xknot(Jay,2))^2+(P(3)-xknot(Jay,3))^2;
Pwork=sqrt(Pwork);
if(ikind==1)
breeze=breeze+alfcough(Jay)*Pwork;
end
if(ikind==3)
breeze=breeze+alfcough(Jay)*Pwork^3;
end
if(ikind==5)
breeze=breeze+alfcough(Jay)*Pwork^5;
end
if(ikind==7)
breeze=breeze+alfcough(Jay)*Pwork^7;
end
%! Big if for behavior of Pwork near zero.
if Pwork > 0.01
if(ikind==2)
breeze=breeze+alfcough(Jay)*Pwork^2*log(Pwork);
end
if(ikind==4)
breeze=breeze+alfcough(Jay)*Pwork^4*log(Pwork);
end
if(ikind==6)
breeze=breeze+alfcough(Jay)*Pwork^6*log(Pwork);
end
end
if Pwork < 0.01
if(ikind==2)
breeze=breeze+alfcough(Jay)*Pwork*log(Pwork^Pwork);
end
if(ikind==4)
breeze=breeze+alfcough(Jay)*Pwork^3*log(Pwork^Pwork);
end
if(ikind==6)
breeze=breeze+alfcough(Jay)*Pwork^5*log(Pwork^Pwork);
end
end
end;
