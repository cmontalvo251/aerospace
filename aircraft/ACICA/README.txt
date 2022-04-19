Hey Carlos,

Ok, here are all the files.

I used airfoildata.m to help me decide what airfoil data to use. It just takes a bunch of files, which I've included for Reynolds numbers of 100,000 and 200,000, and plots things like the lift curve and lift to drag ratio. I got the data from here if you want to choose your own -- there are lots I didn't look at.
Here are the airfoils tested: http://www.ae.illinois.edu/m-selig/uiuc_lsat_airfoilsTested.html
Here is the data for them: http://www.ae.illinois.edu/m-selig/pd.html
These are all low Reynolds number tests so if you're aircraft is going to be significantly bigger/faster you'll have to get data from elsewhere.

StabilityDerivs.m calculates all of the aero coefficients/stability derivatives and control derivatives. Just beware that there are some variables that I calculated in other codes and copied and pasted in, like 'xcgstruct'.  I also have some of the masses of the innards of the aircraft that I pulled from the GenMAV (the old airplane that we were using before from the folks at Eglin AFB). Some of those you might want to keep and other not, like the shocks are the articulated wings so you won't need those and if you're aircraft is gas powered then you won't need the battery, etc. Note that it also reads in an airfoil data file.

TractorDesign.m is more of a design code. There are spots when you can comment or uncomment to make plots of cruise velocity and wing loading for different span, root chord, and wing taper ratio. Again, it requires airfoil data files for the wing and horizontal and vertical tails. It will also plot required thrust vs. velocity so you can find your most efficient cruise speed. It also calculates your required angle of attack and elevator angle to trim so you can make sure trim is feasible and also calculates the rudder angle necessary for max crosswind landing to make sure that's feasible.

stingdesign.m and Inertias.m do the cg and inertia calculations. Sting design also has other useless stuff for figuring out how the sting can be from the aircraft to the load cell but you don't need that.


Sorry it's really messy and kind of allover the place.


XFOIL AND NACA STUFF AND DISTRIBUTED LOAD STUFF
https://web.iit.edu/sites/web/files/departments/academic-affairs/Academic%20Resource%20Center/pdfs/Distributed_Loading.pdf
http://web.mit.edu/16.unified/www/SPRING/fluids/xfoil_info.txt

LOAD airfoil_file  (skip this if airfoil_file argument is given to xfoilexec)

NACA 0010          (can create a NACA airfoil if airfoil_file is not used)

OPER               (enter OPER menu)
?                  (show full OPER menu; commands below are from this menu)

A    5             (specify alpha; can also type full command ALFA)
ASEQ 0 8 2         (specify alpha sequence  0 2 4 6 8 )
CPV                (make spiffy Cp vector plot of last point computed, 8 deg)
CPX                (make usual Cp(x) plot again)
HARD               (write current plot to hardcopy PostScript file plot.ps; optional)
CPWR filename      (write Cp(x) distribution to text file; optional)

