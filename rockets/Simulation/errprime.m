function out = errprime(Tk,V,time)


N = length(V);

dCd = 0.01;

[veuler2, teuler2] = eulerbeta(time,V,Tk+dCd);

vinterp2 = interp1(teuler2,veuler2,time);

err2 = (1/N)*sum((V-vinterp2));

[veuler1, teuler1] = eulerbeta(time,V,Tk);

vinterp1 = interp1(teuler1,veuler1,time);

err1 = (1/N)*sum((V-vinterp1));

out = (err2-err1)/dCd;


