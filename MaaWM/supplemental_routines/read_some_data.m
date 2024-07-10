function [dataflag,xknot,alfu_uno,alfv_dos,alfw_tres,uveeze,vveeze,wveeze] = read_some_data(root_directory,K35);

dataflag = 0;

Xcenfile = [root_directory,'Xcen.dat'];
Ucenfile = [root_directory,'Ualf.dat'];
Vcenfile = [root_directory,'Valf.dat'];
Wcenfile = [root_directory,'Walf.dat'];
Uveezefile = [root_directory,'Uveeze.dat'];
Vveezefile = [root_directory,'Vveeze.dat'];
Wveezefile = [root_directory,'Wveeze.dat'];

fidcen=fopen(Xcenfile);
fidalf1=fopen(Ucenfile);
fidalf2=fopen(Vcenfile);
fidalf3=fopen(Wcenfile);
if fidcen > 0
    NRANK2014=fscanf(fidcen, '%i',1);
else
    disp([Xcenfile,' not found'])
    return
end
NRANkdum1=fscanf(fidalf1, '%i',1);
NRANkdum2=fscanf(fidalf2, '%i',1);
NRANkdum3=fscanf(fidalf3, '%i',1);

for kayfabe=1:NRANK2014
xknot(kayfabe,1)=fscanf(fidcen, '%f',1);
xknot(kayfabe,2)=fscanf(fidcen, '%f',1);
xknot(kayfabe,3)=fscanf(fidcen, '%f',1);
end;

for kayfabe=1:NRANK2014
alfu_uno(kayfabe)=fscanf(fidalf1, '%f',1);
end;

for kayfabe=1:NRANK2014
alfv_dos(kayfabe)=fscanf(fidalf2, '%f',1);
end;

for kayfabe=1:NRANK2014
alfw_tres(kayfabe)=fscanf(fidalf3, '%f',1);
end;

dataflag = 1;

if K35
    uveeze = dlmread(Uveezefile);
    vveeze = dlmread(Vveezefile);
    wveeze = dlmread(Wveezefile);
else
    uveeze = 0;
    vveeze = 0;
    wveeze = 0;
end
