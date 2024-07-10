% Representative example for building velocity from two sources.
function veecomp = Vbuild(P,veeze,k35)
veecomp=0;

VEK=vecpoly(P);
for k=1:k35
    veecomp=veecomp+VEK(k)*veeze(k,1);  % 2 for the y-component, etc.
end;
