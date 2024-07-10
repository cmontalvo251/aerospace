%%%%
purge

%%%Generate a windfield
%dataloc = '/home/carlos/Documents/GitLab_Repos/Research/Data25dx_HF_0/';
dataloc = '/home/carlos/Files/GitLab_Repos/Research/MaaWM/Auto_Generated_WindField/Wave1/';

%%%To generate a windfield use the function
a0 = (1-2*[rand rand rand])*5;
%%%For starters let's make a wind field with a wavelength of 200 meters
tau = 500;
w0 = pi/tau;
N = 20; %%%10 terms in this fourier series
N_vec = 1:N;
wavelengths = 2*pi./(N_vec*w0);

% 1000, 500, 333.33, 250, 200, 166.67, 142.86, 125, 111.11, 100, 90.909, 83.333

% 76.923, 71.429, 66.667, 62.5, 58.824, 55.556, 52.632, 50

%%%Sine Terms
bnu = zeros(1,length(wavelengths));
bnu(2) = 3;  %%%5th term is 200 wavelength. %%2nd term is 500
bnv = bnu*rand*2;
bnw = bnv*rand*2;
phi = (1-2*rand(2,3))*pi;

%%%Cosine Terms
anu = 0*bnu;
anv = anu;
anw = anv;
theta = (1-2*rand(2,3))*pi;

an = [anu;anv;anw];
bn = [bnu;bnv;bnw];
%%%The last number is the dimension of the matrix. Nominal or WRF is 40.
%%%The first number is the spatial wavelength
Auto_Generate_WindField(dataloc,100,an,bn,a0,phi,theta,100); 

%break

%%%%Test wind field on both axes
[wavelengthx,Euvw_avgx] = wind_freq(dataloc,1);
[wavelengthy,Euvw_avgy] = wind_freqY(dataloc,1);

plottool(1,'Error ALL',12,'Wavelength (m)','Normalized Error')
plot(wavelengthx,Euvw_avgx,'b-','LineWidth',2)
plot(wavelengthy,Euvw_avgy,'r-','LineWidth',2)
EFuvw_avg = 0.5*Euvw_avgx + 0.5*Euvw_avgy;
plot(wavelengthx,EFuvw_avg,'g-','LineWidth',2)
legend('E_N','F_N','EF_N')

