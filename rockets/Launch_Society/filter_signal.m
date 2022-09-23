function out = filter_signal(in,sig)

%sig = 0.05;
out = 0*in;
out(1) = in(1);
for idx = 2:length(in)
    out(idx) = out(idx-1)*(1-sig) + in(idx)*sig;
end