function Mdipole = dipole(Mrad,B,B500)

dipole_moment = Mrad/(B500*1e-9);
Mdipole = B.*1e-9.*dipole_moment;