function savefigures2pdf(outputfilename)

if isunix
    system('rm *.pdf');
elseif ispc
    system('del *.pdf');
end

ii = 50;
while ~isempty(get(0,'CurrentFigure') )
    f = getfilename(ii,5);
    saveas(gcf,['Frame_',f,'.pdf'])
    close(gcf)
    ii = ii - 1;
end
%%%On linux you can combine all pdfs into 1
if isunix
    disp('Combining pdfs')
    system('combinepdfs');
    system(['mv combined.pdf PDFs/',outputfilename]);
end

function filename = getfilename(i,L)
%%%filename = getfilename(i,L)
%%%%L = max number of digits
%%%i = current frame number
%%%use this in conjunction with save_frame(getfilename(framenumber,maxdigits))
numdx = [];
for ii = 1:L
  numdx = [numdx,'0'];
end
numi = num2str(i);
filename = [numdx(1:end-length(numi)),numi];
% Copyright - Carlos Montalvo 2015
% You may freely distribute this file but please keep my name in here
% as the original owner
