function [DYNSKIP,BOUNDARY,IKIND,K35] = read_boundary(filename)

fid = fopen(filename);

dynskipstr = fgetl(fid);

L = find(dynskipstr=='!')-1;
DYNSKIP = str2num(dynskipstr(1:L));

for ii = 1:7
    BSTR = fgetl(fid);
end

L = find(BSTR == '!')-1;

BOUNDARY = str2num(BSTR(1:L));

datastr = fgetl(fid);
ikindstr = fgetl(fid);

L = find(ikindstr=='!')-1;

IKIND = str2num(ikindstr(1:L));

K35STR = fgetl(fid);

L = find(K35STR =='!')-1;

K35 = str2num(K35STR(1:L));

switch K35
    case 4
        K35 = 35;
    case 3
        K35 = 20;
    case 2
        K35 = 10;
    case 1
        K35 = 4;
    otherwise
        K35 = 0;
end

fclose(fid);
