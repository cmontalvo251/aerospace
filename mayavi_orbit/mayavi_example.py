import numpy as np
from mayavi import mlab
import scipy.spatial.transform as strans

def make_scene(time=0):
    dphi, dtheta = np.pi/250.0, np.pi/250.0
    [phi,theta] = np.mgrid[0:np.pi+dphi*1.5:dphi,0:2*np.pi+dtheta*1.5:dtheta]

    phi0 = np.arange(0, np.pi+dphi*1.5, dphi) # azimuth
    theta0 = np.arange(0, 2*np.pi+dtheta*1.5, dtheta) #

    orbit_radius = 5
    decl = 20 * np.pi / 180

    periodmoon = 27.3*24
    periodearth = 24


    moonrot = strans.Rotation.from_euler('zyx',
        [np.pi/2 + time / periodmoon * np.pi * 2,
        -decl, 0], degrees=False)

    def apply_surf(x, y, z, rot):
        shape = np.shape(x)
        X = np.vstack([x.flat, y.flat, z.flat])
        Xnew = rot.apply(X.T).T

        x = Xnew[0, :].reshape(shape)
        y = Xnew[1, :].reshape(shape)
        z = Xnew[2, :].reshape(shape)
        return x, y, z


    # Moon orbit


    xorb = orbit_radius * np.cos(theta0)
    yorb = orbit_radius * np.sin(theta0)
    zorb = 0 * np.cos(theta0)

    X = moonrot.apply(np.vstack([xorb, yorb, zorb]).T).T

    #s = mlab.plot3d(X[0, :], X[1, :], X[2, :])
    s = mlab.plot3d(xorb,yorb,zorb)

    # Moon

    r = 0.4
    x = r*np.cos(phi)*np.sin(theta) + orbit_radius
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)

    x,y,z = apply_surf(x, y, z, moonrot)

    #s = mlab.mesh(x, y, z, color=(0.7, 0.7, 0.7))

    # Earth
    r = 1
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)

    #s = mlab.mesh(x, y, z, opacity=1, color=(0., 0.9, 0.0))

    # Equator
    x = r*np.cos(theta0)
    y = r*np.sin(theta0)
    z = r*np.cos(theta0) * 0

    #s = mlab.plot3d(x, y, z, )

    # Axis
    x = np.zeros(10)
    y = np.zeros(10)
    z = np.linspace(-r*1.3, r*1.3, 10)
    #s = mlab.plot3d(x, y, z, color=(0, 0, 0))

    # Victoria
    ph = time * np.pi * 2 / periodearth
    x = r * np.cos(ph) * np.cos(49*np.pi/180)
    y = r * np.sin(ph) * np.cos(49*np.pi/180)
    z = r* np.sin(49*np.pi/180)

    #s = mlab.points3d(x, y, z, mode='sphere', scale_factor=.4)

    # parallel at Victoria

    x = r * np.cos(theta0) * np.cos(49*np.pi/180)
    y = r * np.sin(theta0) * np.cos(49*np.pi/180)
    z = r* np.sin(49*np.pi/180) + 0 * theta0

    #s = mlab.plot3d(x, y, z)

    # Tidal Bulge

    x = 1.7 * r * np.cos(phi)*np.sin(theta)
    y = 1.2 * r * np.sin(phi)*np.sin(theta)
    z = 1.2 * r * np.cos(theta)

    x, y, z = apply_surf(x, y, z, moonrot)

    #s = mlab.mesh(x, y, z, opacity=0.5, color=(0., 0.5, 1.0))

    X,Y = np.mgrid[-10:10:0.2,-10:10:0.2]

    # Bulge Equator

    x = 1.7 * r * np.cos(theta0)
    y = 1.2*r*np.sin(theta0)
    z = 0 * r*np.cos(theta0)

for time in np.arange(0, 49, 0.5) + 0*24:
    make_scene(time=time)
    mlab.view(azimuth=-0, elevation=0, distance=20)
    mlab.orientation_axes()
    mlab.savefig(f'moonearth/Above{int(time*10):04d}.png')
    mlab.clf()

mlab.show()
