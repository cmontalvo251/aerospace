function hgfhgf = tunes2pianosound(m)

% m - is a vector with tunes from pianonotes 1 to 88 played at the same time.
% the lenght of the tunes can't be controlled

M=length(m);
for i=1:M
    strinG=['data/' num2str(m(i)) '.mat'];
    load(strinG)
    try
        y_sound=y_sound+y;
    catch
        y_sound=y;
    end
end
soundsc(y_sound,Fs)
