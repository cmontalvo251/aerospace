IGRF -- Version 11 ----- March 14, 2010 (updated April 30, 2010)

The International Geomagnetic Reference Field (IGRF) model is the empirical 
representation of the Earth's magnetic field recommended for scientific use 
by a special Working Group of the International Association of Geomagnetism 
and Aeronomy (IAGA). The IGRF model represents the main (core) field without 
external sources. The model employs the usual spherical harmonics expansion 
of the scalar potential in geocentric coordinates. The IGRF model coefficients 
are based on all available data sources including geomagnetic measurements 
from observatories, ships, aircrafts and satellites.

The IGRF model consists of sets of coefficients for a global representation 
of the Earth magnetic field for the years 1945, 1950, 1955, etc. There are 
definitive coefficient sets (DGRF####.DAT; #### = year) for which no further 
revisions are anticipated and IGRF####.DAT and IGRF####S.DAT for which future 
updates are expected. IGRF####S.DAT provides the first time derivatives of 
the coefficients for extrapolation into the future. The 11th generation of the 
IGRF model (IGRF-11) consists of definitive coefficients sets for 1900 thru 
2005 (DGRF1945 thru DGRF2005) and preliminary sets for 1900 to 1940 and for 
2010 (IGRF2010) and for extrapolating from 2005 to 2010 (IGRF2010s.DAT).

In combination with the IGRF coefficient sets different subroutines have been 
used to determine the components of the magnetic field vector and the L-value 
at a given location. This sofware package uses the subroutines FELDG (magnetic 
field vector) and SHELLG (L shell) developed by G. Kluge at the European Space 
Operations Center (ESOC). His use of inverse cartesian co- ordinates simplifies 
the computation. The IGRF subroutines were developed by A. Zunde of the U.S. 
Geological Survey (USGS). FELDG and SHELLG are included togehter with other 
funbctions and subroutines in the file SHELLIG.FOR. The program BILCAL.FOR 
produces tables of the geomagnetic field strength, vector components (B-abs., 
B-north, B-east, B-down, declination, inclination), equatorial/minimum field 
strength (B0), dipole moment, and L-value in latitude, longitude (geodetic), 
altitude, or year (decimal). In addition to the main driver program BILCAL.FOR
this distribution also includes a subroutine for computation of IGRF parameters
in the file IGRF_SUB.FOR. 

The IGRF homepage is at http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html.
