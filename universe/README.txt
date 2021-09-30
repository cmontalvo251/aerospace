=====================================================================
  These data are to be used as described in the related document
  titled "Keplerian Elements for Approximate Positions of the
  Major Planets" by E.M. Standish (JPL/Caltech) available from
  the JPL Solar System Dynamics web site (http://ssd.jpl.nasa.gov/).
=====================================================================

Hey if the government isn't shutdown you can access the JPL network

telnet horizons.jpl.nasa.gov 6775

You can then search by planet to see alot of the data. Pretty neat

Table 1.

Keplerian elements and their rates, with respect to the mean ecliptic
and equinox of J2000, valid for the time-interval 1800 AD - 2050 AD.

Right Ascension is the "yaw" angle measured from the vernal equinox in the celestial equator I'm assuming the sun?
Declination is the angle from the plane of the celestial equator check page 56 in Fundamentals of Astrodynamics
I'm assuming that the parameters below are these classic orbital elements which are defined on page 58 of fundamentals of astro dynamics

Ok actually these 6 orbital elements are slightly different. I recommend reading the pdf (aprx_pos_planets.pdf)

If you'd like I would recommend taking a look at page 59 of the orbital mechanics book there is a great figure there

EM Bary means the Earth Moon Barycenter which is an approximation of the center of mass of the Earth Moon System.

Julian Day - This is the integer day from the start of the gregorian calendar.
On January 1st 2000 this number was - 2 451 545
J2000 is then the reference frame used at Julian Day 2000 where apparently the tilt of the earth was about 23.4 degrees
Apparently the J2000 frame is an international celestial reference frame - ICRF and in 1998 they created the ICRF1 frame and in 2009 it was updated to the ICRF2 frame.
J2000 is also apparently an ECI or Earth centered inertial frame

AU is astronomical unit - 149597870700 meters (wikipedia) - 149597870.700 km (JPL) - These are the same. lol.
/Cy is actually per century

T is the number of centuries past J2000  --- T = (Teph - 2 451 545 ) / 36525

where 2451545 is the Julian day Jan 1 2000 and 36525 is number of days in a century.

So if you want a at a particular date in time you need to first know the ephemeris day Teph and then compute the actual value of a like this

a = a0 + adot * T 

makes sense? You can also just count the number of days past Jan 1st 2000 and divide that number by 36525 days which is the number of days in a century (365.25 * 100 = 36525)

According to the PDF titled aprx_pos_planets.pdf the orbital elements below are defined as follows

a = semi major axis AU
e = eccentricity
i = inclination
L = mean longitude 
long.peri. = longitude of the perihelion (wbar)
long.node. = longitude of the ascending node (OMEGA)

Again take a look at page 59 in fundamentals of astrodynamics for this to make more sense.

However, this pdf from JPL, aprx_pos_planets.pdf has a pretty good derivation of how to get the cartesian coordinates of every planet using the 6 orbital elements below. 
Remember you have 6 degrees of freedom so you need 6 coordinates. Converting these quantities to 3 cartesian coordinates is a reduction in order but that's also because you can get velocity from these orbital elements as well but that's not so important.

These are only valid until 2050 AD btw.

               a              e               I                L            long.peri. (wbar)     long.node. (OMEGA)
           AU, AU/Cy     rad, rad/Cy     deg, deg/Cy      deg, deg/Cy      deg, deg/Cy            deg, deg/Cy
-----------------------------------------------------------------------------------------------------------
Mercury   0.38709927      0.20563593      7.00497902      252.25032350     77.45779628     48.33076593
          0.00000037      0.00001906     -0.00594749   149472.67411175      0.16047689     -0.12534081
Venus     0.72333566      0.00677672      3.39467605      181.97909950    131.60246718     76.67984255
          0.00000390     -0.00004107     -0.00078890    58517.81538729      0.00268329     -0.27769418
EM Bary   1.00000261      0.01671123     -0.00001531      100.46457166    102.93768193      0.0
          0.00000562     -0.00004392     -0.01294668    35999.37244981      0.32327364      0.0
Mars      1.52371034      0.09339410      1.84969142       -4.55343205    -23.94362959     49.55953891
          0.00001847      0.00007882     -0.00813131    19140.30268499      0.44441088     -0.29257343
Jupiter   5.20288700      0.04838624      1.30439695       34.39644051     14.72847983    100.47390909
         -0.00011607     -0.00013253     -0.00183714     3034.74612775      0.21252668      0.20469106
Saturn    9.53667594      0.05386179      2.48599187       49.95424423     92.59887831    113.66242448
         -0.00125060     -0.00050991      0.00193609     1222.49362201     -0.41897216     -0.28867794
Uranus   19.18916464      0.04725744      0.77263783      313.23810451    170.95427630     74.01692503
         -0.00196176     -0.00004397     -0.00242939      428.48202785      0.40805281      0.04240589
Neptune  30.06992276      0.00859048      1.77004347      -55.12002969     44.96476227    131.78422574
          0.00026291      0.00005105      0.00035372      218.45945325     -0.32241464     -0.00508664
Pluto    39.48211675      0.24882730     17.14001206      238.92903833    224.06891629    110.30393684
         -0.00031596      0.00005170      0.00004818      145.20780515     -0.04062942     -0.01183482

Correction terms for Jupiter through Pluto are below.

Planet  b   c   s   f
Jupiter -0.00012452 0.06064060 -0.35635438 38.35125
Saturn 0.00025899 -0.13434469 0.87320147 38.35125
Uranus 0.00058331 -0.97731848 0.17689245 7.67025
Neptune -0.0041348 0.68346318 -0.10162547 7.67025
Pluto -0.01262724 0.0 0.0 0.0


Alright so using the PDF aprx_pos_planets.pdf you can use the elements above....

0.) Compute T for the current day
T = (Teph - 2 451 545 ) / 36525

1.) First compute each of the planets 6 elements for the particular day of the year
a = a0 + adot*T

2.) Compute the argument of the perihelion, w, and the mean anomaly, M

w = wbar - OMEGA
M = L - wbar + bT**2 + c*cos(fT) + s*sin(f*T)

where the last 3 terms must be added for Jupiter, Saturn, Uranus, Neptune and Pluto up to 3000 AD. I'll be long dead by then anyway.

3.) Modulus M such that M is in the set of +-180 degrees and then solve the kepler equation below to solve for the Eccentric Anomaly

M = E - estar*sin(E)

where estar is e in degrees such that

estar = 180/pi*e

4.) Compute planet's heliocentric coordinates in it's orbital plane

xprime = a*(cos(E)-e)
yprime = a*sqrt(1-e**2)*sin(E)
zprime = 0

5.) Compute the coordinates in the J2000 ecliptic plane

x = (cos(w)*cos(OMEGA) - sin(w)*sin(OMEGA)*cos(I))*xprime + (-sin(w)*cos(OMEGA)-cos(w)*sin(OMEGA)*cos(I))*yprime

y = (cos(w)*sin(OMEGA)+sin(w)*cos(OMEGA)*cos(I))*xprime + (-sin(w)*sin(OMEGA)+cos(w)*cos(OMEGA)*cos(I))*yprime

z = (sin(w)*sin(I))*xprime + (cos(w)*sin(I))*yprime

When you are testing this code there are two images in this folder titled orbit_plot_inner and orbit_plot_outer. These were taken on Jan 1 2018. Note that the Julian day for Jan 1st 2018 is 2 458 120  

Since I will probably end up using this code in 2019 here is the Julian Day for Jan 1 2019 - 2 458 485

Oh and I added a Julian_Day.png which shows the Julian Day per year. The next leap year is 2020. 
