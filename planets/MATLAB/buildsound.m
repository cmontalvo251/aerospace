input = outvec;

% input = [0 64;
%          .04 40;
%          .04 44;
%          .04 47;
%          .1 59;
%          .20 64;
%          .24 40;
%          .24 44;
%          .24 47;
%          .3 59]
fileout = 'run_2';
time = input(:,1)*tscale;

notes = input(:,2);

y_sound = zeros(100,2);
M=size(input,1);
cclim = 0;

[M,c] = size(outvec);

for ii=1:M
    if ii/M*100 >=cclim
        fprintf('%6.2f%% Complete\n',cclim);
        cclim = cclim+1;
    end

    t = time(ii);

    strinG=['data/' num2str(notes(ii)) '.mat'];
    load(strinG);

    tip = round(t*Fs);

    pad = size([zeros(tip-1,2); y],1)-size(y_sound,1);
%     size(padarray(y_sound,[pad,0],'post'));
%     size([zeros(tip-1,2); y]);
    y_sound = padarray(y_sound,[pad,0],'post')+[zeros(tip-1,2); y];
%     pause
end

save('soundfile.mat','y_sound','Fs','time')

t = [1:size(y_sound,1)]/Fs;
figure;
plot(t(1:length(y_sound)),y_sound);
