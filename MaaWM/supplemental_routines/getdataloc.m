function dataloc = getdataloc(ifiles)

fid = fopen(ifiles);

if (fid==-1) 
    disp(['sorry ',ifiles,' not found'])
    dataloc = [];
end

filetype = 'START';

while strcmp(filetype,'ATM') == 0
    line = fgetl(fid);
    line = line(2:end); %%%Delete first tick mark
    loc = find(line=='''',1)-1;
    filetype = line(1:loc);
end
line = line(loc+2:end); %%%Delete filetype
line = line(1:find(line=='!')-1); %%%Delete comment in the file
line(line == '''') = []; %%%Delete tick marks
line(line == ' ') = []; %%%Delete space
%%%Get rid of tabs
temp = double(line); %%%Convert to ASCII equivalent
temp(temp == 9) = []; %%%9 is ASCII for tab
line = char(temp); %%Convert back to characters

ATMFILE = line;

%%%Open ATM file
fid = fopen(ATMFILE);
if (fid == -1)
    disp(['sorry ',ATMFILE,' not found']);
end

%%%Line we care about is on line 9
for idx = 1:10
    line = fgetl(fid);
end

line = line(2:find(line=='!')-1); %%%Get rid of first tick mark and comment
line(line=='''')=[]; %%%Get rid of tick makrs
line(line==' ') = []; %%%Get rid of spaces

dataloc = line;