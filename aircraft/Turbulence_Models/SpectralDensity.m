purge

%Regenerate S
N = [100,100];
wupper = [1,1].*0.5;
wlower = [-1,-1].*0.5;
dw = (wupper-wlower)./N;

S11 = zeros(N(1),N(2));
S12 = S11;
S22 = S11;

w1 = linspace(wlower(1),wupper(1),N(1));
w2 = linspace(wlower(2),wupper(2),N(2));
[ww1,ww2] = meshgrid(w1,w2);

sigu = 1;
Lu = 5000/3.28;
alfa = 4*(sigu^2)*(Lu^5)/(pi^2);
w3 = 0.5;
clip = 10;
for ii = 1:N(1)
  w1k1 = wlower(1) + (ii-1/2)*dw(1);
  for jj = 1:N(2)
    w2k2 = wlower(2) + (jj-1/2)*dw(2);
    w = [w1k1;w2k2;w3];
    W = sqrt(w(1)^2+w(2)^2+w(3)^2);
    wp = w;
    e = 1/(1+(Lu*W)^2)^3;
    f = alfa*e;
    %S
    W2 = W^2;
    S11(ii,jj) = f*(W2-wp(1)^2);
    if (abs(S11(ii,jj))>clip) S11(ii,jj) = clip*sign(S11(ii,jj)); end
    S12(ii,jj) = f*-wp(1)*wp(2);
    if (abs(S12(ii,jj))>clip) S12(ii,jj) = clip*sign(S12(ii,jj)); end
    S22(ii,jj) = f*(W2-wp(2)^2);
    if (abs(S22(ii,jj))>clip) S22(ii,jj) = clip*sign(S22(ii,jj)); end
    S31(ii,jj) = f*-wp(1)*wp(3);
    if (abs(S31(ii,jj))>clip) S31(ii,jj) = clip*sign(S31(ii,jj)); end
    S32(ii,jj) = f*-wp(2)*wp(3);
    if (abs(S32(ii,jj))>clip) S32(ii,jj) = clip*sign(S32(ii,jj)); end
    S33(ii,jj) = f*(W2-wp(3)^2);
    if (abs(S33(ii,jj))>clip) S33(ii,jj) = clip*sign(S33(ii,jj)); end
  end
end

plottool(1,'Name',12,'w1','w2','S11');
mesh(ww1,ww2,S11)

plottool(1,'Name',12,'w1','w2','S12');
mesh(ww1,ww2,S12)

plottool(1,'Name',12,'w1','w2','S22');
mesh(ww1,ww2,S22)

plottool(1,'Name',12,'w1','w2','S31');
mesh(ww1,ww2,S31)

plottool(1,'Name',12,'w1','w2','S32');
mesh(ww1,ww2,S32)

plottool(1,'Name',12,'w1','w2','S33');
mesh(ww1,ww2,S33)



drawnow

