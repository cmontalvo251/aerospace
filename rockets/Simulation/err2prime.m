function out = err2prime(Tk,time,V)


N = length(V);

dCd = 0.01;

[veuler2, teuler2] = eulerbeta(time,V,Tk+2*dCd);

vinterp2 = interp1(teuler2,veuler2,time);

err2 = (1/N)*sum((V-vinterp2));

[veuler1, teuler1] = eulerbeta(time,V,Tk+dCd);

vinterp1 = interp1(teuler1,veuler1,time);

err1 = (1/N)*sum((V-vinterp1));

[veuler0, teuler0] = eulerbeta(time,V,Tk);

vinterp0 = interp1(teuler0,veuler0,time);

err0 = (1/N)*sum((V-vinterp0));

out = (err2-2*err1+err0)/dCd^2;
