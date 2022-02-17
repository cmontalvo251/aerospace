function WINDout = WINDsim(WINDin,UVW,NOISETYPE)

[r,c] = size(WINDin);


%%If this code is run it assumes that UVW already
%%contains errors in it.
%%if this isn't true then you should run the code below

%%%RUN THIS IF UVW DOESN'T ALREADY HAVE ERRORS IN IT


%UVW = UVWsim(UVW,NOISETYPE);


%%%TO BACK OUT wind disturbance we need to first add in
%%pitot tube errors

WINDout = WINDin;

if NOISETYPE == 1 %%white noise
  windNoise = 2; %%m/s (std deviation)
  WINDout(1:3,:) = WINDin(1:3,:) - windNoise + (2*windNoise).*rand(3,c);
end

%%Now that pitot-static errors have been added,
%%we can back out the actual wind disturbance

WINDout = WINDout - UVW;
  