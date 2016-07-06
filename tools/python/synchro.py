import read as read
import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.interpolate import griddata
from congrid import congrid
import vtk as v
#import vtktonumpy as ah
#from tvtk import array_handler as ah
import numpy_support as ah
import img_scale
from scipy import spatial
from scipy.ndimage.filters import gaussian_filter
from matplotlib import pyplot
import matplotlib as mpl
import pickle, os

coord='cartesian'

def get_pointdata(offset,filenameout='data',type='pvtu',attribute_mode='cell',mirror=None):

    if type == 'vtu':
        filename=''.join([filenameout,repr(offset).zfill(4),'.vtu'])
        datareader = v.vtkXMLUnstructuredGridReader()
    elif type == 'pvtu':
        filename=''.join([filenameout,repr(offset).zfill(4),'.pvtu'])
        datareader = v.vtkXMLPUnstructuredGridReader()
        
    print '=== Reading ',filename,' ==='
        
    datareader.SetFileName(filename)
    datareader.Update()
    data = datareader.GetOutput()

    if attribute_mode == 'cell':
        c2p = v.vtkCellDataToPointData()
        c2p.SetInput(data)
        pointdata=c2p.GetOutput()
        pointdata.Update()
    else:
        pointdata = data 

    if mirror != None:
        vr=v.vtkReflectionFilter()
        vr.SetInput(pointdata)
        vr.SetPlane(mirror)
        pointdata=vr.GetOutput()
        vr.Update()
        pointdata.Update()

    vtk_points=pointdata.GetPoints().GetData()
    points=ah.vtk2array(vtk_points).astype(np.float64)
 
    print '=== Done with reading data! ==='
    return {'pointdata': pointdata, 'points': points}


def rot2D(offset,nphi,filenameout='data',type='pvtu'):
    """
    creates a mock vtu datastructure with rotated values.
    No connectivity is set, vtu is just used as a container for pointdata     since
    the following routines expect to get a vtu unstructured datatype.
    """

    data=get_pointdata(offset,filenameout=filenameout,type=type)
    points2D=data.get('points')
    data3D = v.vtkUnstructuredGrid()

    nelem2D = len(points2D)
    nelem3D = nelem2D*nphi

    print 'preparing ',nelem3D, 'points'

    points3D = np.empty((nelem3D,3))

    phis=np.linspace(0.,2.*np.pi*float(nphi)/float(nphi+1),nphi)

    i3D = 0
    for i2D in range(nelem2D):
# pick next 2D point:
        point = points2D[i2D]
        r = point[0]
        z = point[1]
        for phi in phis:
            sinphi = np.sin(phi)
            cosphi = np.cos(phi)
            points3D[i3D,0] = r*cosphi
            points3D[i3D,1] = r*sinphi
            points3D[i3D,2] = z
            i3D=i3D+1

    vtkpoints=ah.array2vtkPoints(points3D)
    data3D.SetPoints(vtkpoints)

# rotate variables:
    ivar = 0
    while ivar < data.get('pointdata').GetPointData().GetNumberOfArrays():
        var = data.get('pointdata').GetPointData().GetArrayName(ivar)
        print 'rotating variable:',ivar,var
        i3D=0
# Treat vectors:
        if var == 'u1' or var == 'b1' or var == 'v1':
            array2D_vec = np.empty((nelem2D,3))
            array3D_vec = np.empty((nelem3D,3))
            bx = np.empty(nelem3D)
            by = np.empty(nelem3D)
            bz = np.empty(nelem3D)
            array2D_vec[:,0] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar),attribute_mode='point')
            array2D_vec[:,1] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar+1),attribute_mode='point')
            array2D_vec[:,2] = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar+2),attribute_mode='point')
            for i2D in range(nelem2D):
                for phi in phis:
# Should have three components (makes no sense otherwise):
# Input in cylindrical coordinates, convention:
# b1==br; b2==bz; b3==bphi  and similar for velocities.
                    array3D_vec[i3D,0] = array2D_vec[i2D,0]
                    array3D_vec[i3D,1] = array2D_vec[i2D,1]
                    array3D_vec[i3D,2] = array2D_vec[i2D,2]
# Now create the rotated vectors:
                    sinphi = np.sin(phi)
                    cosphi = np.cos(phi)
                    bz[i3D] = array3D_vec[i3D,1]  
                    bx[i3D] = array3D_vec[i3D,0] * cosphi - array3D_vec[i3D,2] * sinphi
                    by[i3D] = array3D_vec[i3D,2] * cosphi + array3D_vec[i3D,0] * sinphi
                    i3D=i3D+1

            vtkarray3D = ah.array2vtk(bx)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar))
            data3D.GetPointData().AddArray(vtkarray3D)

            vtkarray3D = ah.array2vtk(by)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar+1))
            data3D.GetPointData().AddArray(vtkarray3D)
 
            vtkarray3D = ah.array2vtk(bz)
            vtkarray3D.SetName(data.get('pointdata').GetPointData().GetArrayName(ivar+2))
            data3D.GetPointData().AddArray(vtkarray3D)
  
            ivar = ivar+3
        else:
# Treat scalars:
            array2D = np.empty(nelem2D)
            array3D = np.empty(nelem3D)
            array2D = read.extract(data.get('pointdata'),data.get('pointdata').GetPointData().GetArrayName(ivar),attribute_mode='point')
            for i2D in range(nelem2D):
                for phi in phis:
                    array3D[i3D] = array2D[i2D]
                    i3D=i3D+1
            ivar = ivar+1

            vtkarray3D = ah.array2vtk(array3D)
            vtkarray3D.SetName(var)
            data3D.GetPointData().AddArray(vtkarray3D)
    
    data3D.Update()
    return {'pointdata': data3D, 'points': points3D}


def simpleEmiss(data,phi,theta,var='p',delta=0):

    return {'emiss': read.extract(data.get('pointdata'),var,attribute_mode='point'),'points':data.get('points'),'phi':phi,'theta':theta,'delta':delta,'nu':1e12,'alpha':0.6,'recipe':1}

def emiss(data,phi,theta,nu,alpha,recipe=1,polarized=1,delta=0,noswing=None,epsb=0.):

    print '=== Getting emissivities... ==='

# peel data:
    pointdata=data.get('pointdata')

# Synchrotron constants:
    
    c  = 2.99792458e10 # cm s^-1
    c1 = 6.27e18 # cm^-7/2 g^-5/2 s^4
    c2 = 2.368e-3 # g^-2 cm^-1 s^3
    me = 9.1093897e-28 # g

# Euclidian norm:
    norm = lambda x: np.sqrt(
    np.square(x[:,0]) + np.square(x[:,1]) + np.square(x[:,2])
    )
# Scalar product:
    scalar = lambda x,y: x[:,0]*y[:,0] + x[:,1]*y[:,1] + x[:,2]*y[:,2]

# Gamma from alpha:
    Gamma = 2.*alpha + 1.

# points:
    ux=read.extract(pointdata,'u1',attribute_mode='point')
    uy=read.extract(pointdata,'u2',attribute_mode='point')
    uz=read.extract(pointdata,'u3',attribute_mode='point')
    lfac=read.extract(pointdata,'lfac',attribute_mode='point')

    nelem = len(ux)


# get line of sight vector:
# phi and theta should be given in deg:
    n = np.empty((nelem,3))
    sintheta = np.sin(np.pi/180.*theta)
    costheta = np.cos(np.pi/180.*theta)
    sinphi   = np.sin(np.pi/180.*phi)
    cosphi   = np.cos(np.pi/180.*phi)
    n[:,0] = sintheta*cosphi
    n[:,1] = sintheta*sinphi
    n[:,2] = costheta


# get doppler factor:
    beta = np.empty((nelem,3))
    beta[:,0] = ux/lfac/c
    beta[:,1] = uy/lfac/c
    beta[:,2] = uz/lfac/c

    D = np.empty(nelem)
    D = 1./(
    lfac*(1.-(scalar(n,beta)))
    )

    print 'min and max Doppler factors: ', D.min(), D.max()

# get ndash:
    ndash = np.empty((nelem,3))
    ndash[:,0] = D*n[:,0] - (D+1.) * lfac/(lfac+1) * beta[:,0]
    ndash[:,1] = D*n[:,1] - (D+1.) * lfac/(lfac+1) * beta[:,1]
    ndash[:,2] = D*n[:,2] - (D+1.) * lfac/(lfac+1) * beta[:,2]

    print 'got ndash min and max:', ndash[:,0].min(), ndash[:,0].max(), ndash[:,1].min(), ndash[:,1].max(), ndash[:,2].min(), ndash[:,2].max()

# get bdash:
    bdash = get_bdash(pointdata)

    if epsb > 0. :
# keep direction of magnetic field, but assume magnitude from equipartition 
# with thermal energy.  
# Equipartition fraction is epsb, such that E_b = epsb * 3 * p
        p = read.extract(pointdata,'p',attribute_mode='point')
        normbdash = norm(bdash)
        e = 3.*p + normbdash**2/(8.*np.pi)
        underCutoff = normbdash**2/(8.*np.pi) < epsb * e

#        bdash[underCutoff,0] = bdash[underCutoff,0] * powerorzero(normbdash[underCutoff],-1) * np.sqrt(8.*np.pi*e[underCutoff]/(1./epsb+1.))
#        bdash[underCutoff,1] = bdash[underCutoff,1] * powerorzero(normbdash[underCutoff],-1) * np.sqrt(8.*np.pi*e[underCutoff]/(1./epsb+1.))
#        bdash[underCutoff,2] = bdash[underCutoff,2] * powerorzero(normbdash[underCutoff],-1) * np.sqrt(8.*np.pi*e[underCutoff]/(1./epsb+1.))
        bdash[underCutoff,0] = bdash[underCutoff,0] * powerorzero(normbdash[underCutoff],-1) * np.sqrt(8.*np.pi*epsb*e[underCutoff])
        bdash[underCutoff,1] = bdash[underCutoff,1] * powerorzero(normbdash[underCutoff],-1) * np.sqrt(8.*np.pi*epsb*e[underCutoff])
        bdash[underCutoff,2] = bdash[underCutoff,2] * powerorzero(normbdash[underCutoff],-1) * np.sqrt(8.*np.pi*epsb*e[underCutoff])



    print 'got bdash min and max:', bdash[:,0].min(), bdash[:,0].max(), bdash[:,1].min(), bdash[:,1].max(), bdash[:,2].min(), bdash[:,2].max()

# bdash perp:
    bdashp = np.empty(nelem)
    tmp    = np.empty((nelem,3))
    tmp    = np.cross(ndash,bdash)
    bdashp = norm(tmp)

    print 'got bdashp min and max: ', bdashp.min(), bdashp.max()

# electron densities and cutoff:
    eps_inf=read.extract(pointdata,'eps_inf',attribute_mode='point')

# emissivity:
    if recipe ==1:
        rhoe=read.extract(pointdata,'rhoe',attribute_mode='point')
        rhoe_0=read.extract(pointdata,'rhoe_0',attribute_mode='point')
    elif recipe ==2:
        rhoe=read.extract(pointdata,'n',attribute_mode='point')
        rhoe_0=read.extract(pointdata,'n0',attribute_mode='point')

    if recipe == 1 or recipe ==2:
        emiss = (c2*rhoe_0/me*c2*np.power(c1,(Gamma-3.)/2.)/(8.*np.pi) *
                 np.power(D,2.+(Gamma-1)/2.) *
                 powerorzero(bdashp,(Gamma+1.)/2.)*
                 np.power(rhoe_0/rhoe,-((Gamma+2.)/3.)) *
                 np.power(nu,(1.-Gamma)/2.) *
                 powerorzero((1. - np.sqrt(nu)/(np.sqrt(c1*bdashp)*eps_inf)),Gamma-2.)
             )

    if recipe == 3:
        emiss = (norm(bdash)**2 *
                 np.power(D,2.+(Gamma-1)/2.) *
                 powerorzero(bdashp,(Gamma+1.)/2.)*
                 np.power(nu,(1.-Gamma)/2.)
             )


    if recipe == 4:
        p = read.extract(pointdata,'p',attribute_mode='point')
        emiss = (p *
                 np.power(D,2.+(Gamma-1)/2.) *
                 powerorzero(bdashp,(Gamma+1.)/2.)*
                 np.power(nu,(1.-Gamma)/2.)
             )

# mask out the wind zone:
    m = lfac >9.
    emiss = np.ma.filled(np.ma.masked_array(emiss,m),0.)

    if polarized !=1:
        print '=== Done with emissivities! ==='
        return {'emiss': emiss, 'points': data.get('points'),
                'phi': phi, 'theta': theta, 'delta': delta, 'nu': nu, 'alpha': alpha, 'recipe': recipe, 'epsb': epsb}
# The following only if polarized emission is requested:

# direction of emission:
    q       = np.empty((nelem,3))
    absb    = np.empty(nelem)
    b       = np.empty((nelem,3))
    ehat    = np.empty((nelem,3))
    absehat = np.empty((nelem))
    lhat    = np.empty((nelem,3))
    coschie = np.empty(nelem)
    sinchie = np.empty(nelem)
    cos2chie= np.empty(nelem)
    sin2chie= np.empty(nelem)
     
    b[:,0] =read.extract(pointdata,'b1',attribute_mode='point')
    b[:,1] =read.extract(pointdata,'b2',attribute_mode='point')
    b[:,2] =read.extract(pointdata,'b3',attribute_mode='point')

    absb = norm(b)

# to test the effect of the swing you can use the noswing parameter:
    if noswing==None:
        q = b + np.cross(n,np.cross(beta,b))
    else:
        q = b

    q[:,0] = q[:,0]/absb
    q[:,1] = q[:,1]/absb
    q[:,2] = q[:,2]/absb

    ehat      = np.cross(n,q) 
    absehat   = norm(ehat)
    ehat[:,0] = ehat[:,0]/absehat
    ehat[:,1] = ehat[:,1]/absehat
    ehat[:,2] = ehat[:,2]/absehat

    lhat1 = -costheta*cosphi
    lhat2 = -costheta*sinphi
    lhat3 = sintheta

    if delta!=0:
        cosa = np.cos(delta*np.pi/180.)
        sina = np.sin(delta*np.pi/180.)

        n1 = sintheta*cosphi
        n2 = sintheta*sinphi
        n3 = costheta

        omca  = 1.-cosa
        lhatd1 = (n1**2*omca+cosa)*lhat1 + (n1*n2*omca-n3*sina)*lhat2 + (n1*n3*omca+n2*sina)*lhat3
        lhatd2 = (n2*n1*omca+n3*sina)*lhat1 + (n2**2*omca+cosa)*lhat2 + (n2*n3*omca-n1*sina)*lhat3
        lhatd3 = (n3*n1*omca-n2*sina)*lhat1 + (n3*n2*omca+n1*sina)*lhat2 + (n3**2*omca+cosa)*lhat3

        lhat1=lhatd1
        lhat2=lhatd2
        lhat3=lhatd3
        
    lhat[:,0] = lhat1
    lhat[:,1] = lhat2
    lhat[:,2] = lhat3


    coschie = scalar(lhat,ehat)
    sinchie = scalar(n,np.cross(lhat,ehat))

# complex phase with Doppelwinkelsatz:
    sin2chie = 2.*sinchie*coschie
    cos2chie = np.square(coschie)-np.square(sinchie)

# assemble complex polarized emissivity:
    piem = (Gamma+1)/(Gamma+7./3.)
    emisspol = np.ma.filled(np.ma.masked_array(piem * emiss * (cos2chie + 1j* sin2chie),absb==0),0.)

    print '=== Done with emissivities! ==='
    return {'emiss': emiss, 'emisspol': emisspol, 'points': data.get('points'),
            'phi': phi, 'theta': theta, 'delta': delta, 'nu': nu, 'alpha': alpha, 'recipe': recipe, 'epsb': epsb}



def powerorzero(x,y):
    m = np.logical_or(x<=0., np.isnan(x))
    masked = np.ma.masked_array(np.power(x,y),m)
    return np.ma.filled(masked,0.)



def get_bdash(data):

    c  = 2.99792458e10 # cm s^-1

    bx=read.extract(data,'b1',attribute_mode='point')
    by=read.extract(data,'b2',attribute_mode='point')
    bz=read.extract(data,'b3',attribute_mode='point') 
    lfac=read.extract(data,'lfac',attribute_mode='point')
    ux=read.extract(data,'u1',attribute_mode='point')/lfac
    uy=read.extract(data,'u2',attribute_mode='point')/lfac
    uz=read.extract(data,'u3',attribute_mode='point')/lfac
    
    bdash = np.empty((len(bx),3))
    bdash[:,0] = bx/lfac + lfac/(lfac+1.)*ux*(bx*ux+by*uy+bz*uz)/c/c
    bdash[:,1] = by/lfac + lfac/(lfac+1.)*uy*(bx*ux+by*uy+bz*uz)/c/c
    bdash[:,2] = bz/lfac + lfac/(lfac+1.)*uz*(bx*ux+by*uy+bz*uz)/c/c

    return bdash


def show_emiss(inputdata,zi,nregrid,xrange=[-1.8e18,1.8e18],yrange=[-1.8e18,1.8e18]):
   
    fontsize=18
    nlevels=256

    points = inputdata.get('points')
    emiss = inputdata.get('emiss')

    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

   # this is needed because griddata freaks out otherwise:
    smalldouble=1.e-8

    tmp=(xrange[1]-xrange[0])*smalldouble
    xi = np.linspace(xrange[0]-tmp,xrange[1]+tmp,nregrid[0])
    tmp=(yrange[1]-yrange[0])*smalldouble
    yi = np.linspace(yrange[0]-tmp,yrange[1]+tmp,nregrid[1])

    nelem = nregrid[0]*nregrid[1]
    xij=np.empty(nelem)
    yij=np.empty(nelem)
    zij=np.empty(nelem)
    zij[:]=zi
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            xij[i+nregrid[0]*j] = xi[i]
            yij[i+nregrid[0]*j] = yi[j]

    data = griddata((x,y,z), emiss, (xij,yij,zij), method='nearest')

    data_2D = np.empty((nregrid[0],nregrid[1]))
    for i in range(nregrid[0]):
        for j in range(nregrid[1]):
            data_2D[j,i] = data[i+nregrid[0]*j]

    plotrange=np.array([np.log10(data_2D.max())-4.,np.log10(data_2D.max())])

    plt.figure()
    plt.clf()
    cs = plt.contourf(xi,yi,np.log10(data_2D),levels=np.arange(nlevels)/float(nlevels-1)*(plotrange[1]-plotrange[0])+plotrange[0])
#    cs = plt.contourf(xi,yi,data_2D,nlevels)
    cb = plt.colorbar()

    ax = plt.gca()
    ax.set_aspect('equal')
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

    for tick in cb.ax.get_yticklabels():
        tick.set_fontsize(fontsize)


    plt.ylabel('y [cm]',size=fontsize)
    plt.xlabel('x [cm]',size=fontsize)
 
    plt.show()

def get_ray(data,nregrid,beg,end):


    points = data.get('points')
    emiss = data.get('emiss')

    x=points[:,0]
    y=points[:,1]
    z=points[:,2]

    xi = np.linspace(beg[0],end[0],nregrid)
    yi = np.linspace(beg[1],end[1],nregrid)
    zi = np.linspace(beg[2],end[2],nregrid)

    data = griddata((x,y,z), emiss, (xi,yi,zi), method='nearest')    

    return data



def emissFrom2D(data,theta=0,npixel=[200,200,200],L=[100,100,100],x0=[0,0,0],nu=1,alpha=0.6,delta=0,recipe=2,polarized=0,dog=None,mirror=None):

# Synchrotron constants:
    c  = 2.99792458e10 # cm s^-1
#    c  = 1.

    c1 = 6.27e18 # cm^-7/2 g^-5/2 s^4
    c2 = 2.368e-3 # g^-2 cm^-1 s^3
    me = 9.1093897e-28 # g

# Euclidian norm:
    norm = lambda x: np.sqrt(
    np.square(x[:,0]) + np.square(x[:,1]) + np.square(x[:,2])
    )

    # Scalar product:
    scalar = lambda x,y: x[:,0]*y[:,0] + x[:,1]*y[:,1] + x[:,2]*y[:,2]

    # Gamma from alpha:
    Gamma = 2.*alpha + 1.

    # Since this is from an axisymmetric run, we put phi=0
    phi=0
    grid=make_grid(phi,theta,npixel,L,delta=delta,x0=x0)

    # get back the points:
    points=data.get('points')

# mask out the wind zone:
    lfac = read.extract(data.get('pointdata'),'lfac',attribute_mode='point')
    m = lfac > 9.

    # get velocities for Doppler factor:
    ur = ontoGrid(extract_masked(data,'u1',m)
                  ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)
    uz = ontoGrid(extract_masked(data,'u2',m)
                  ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)
    uphi = ontoGrid(extract_masked(data,'u3',m)
                  ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)

    # get magnetic field:
    br = ontoGrid(extract_masked(data,'b1',m)
                  ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)
    bz = ontoGrid(extract_masked(data,'b2',m)
                  ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)
    bphi = ontoGrid(extract_masked(data,'b3',m)
                  ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)

# mirror the variables if needed:
    if (mirror == 1):
        uz[grid[2]<0]   = - uz[grid[2]<0]
        br[grid[2]<0]   = - br[grid[2]<0]
        bphi[grid[2]<0] = - bphi[grid[2]<0]        



    # Cartesian velocity components:
    rcyl = np.sqrt(grid[0]**2+grid[1]**2)
    cosphi = grid[0] / rcyl
    sinphi = grid[1] / rcyl

    ux = (ur*cosphi - uphi*sinphi)
    uy = (ur*sinphi + uphi*cosphi)

    # Lorentz factor:
    lfac = np.sqrt(1.+(ux**2+uy**2+uz**2)/c**2)
    print '=== Minimum and maximum Lorentz factor: ', lfac.min(), lfac.max(), ' ==='

    nelem = len(ux)

# get line of sight vector:
# phi and theta should be given in deg:
    n = np.empty((nelem,3))
    sintheta = np.sin(np.pi/180.*theta)
    costheta = np.cos(np.pi/180.*theta)
    sinphi_n   = 0.
    cosphi_n   = 1.
    n[:,0] = sintheta*cosphi_n
    n[:,1] = sintheta*sinphi_n
    n[:,2] = costheta


    beta = np.empty((nelem,3))
    beta[:,0] = ux/lfac/c
    beta[:,1] = uy/lfac/c
    beta[:,2] = uz/lfac/c

    # Doppler factor:
    D = np.empty(nelem)
    D = 1./(
    lfac*(1.-(scalar(n,beta)))
    )

    print '=== Minimum and maximum Doppler factor: ', D.min(), D.max(), ' ==='
# get ndash:
    ndash = np.empty((nelem,3))
    ndash[:,0] = D*n[:,0] - (D+1.) * lfac/(lfac+1.) * beta[:,0]
    ndash[:,1] = D*n[:,1] - (D+1.) * lfac/(lfac+1.) * beta[:,1]
    ndash[:,2] = D*n[:,2] - (D+1.) * lfac/(lfac+1.) * beta[:,2]


    print 'got ndash min and max:', ndash[:,0].min(), ndash[:,0].max(), ndash[:,1].min(), ndash[:,1].max(), ndash[:,2].min(), ndash[:,2].max()


    bx = (br*cosphi - bphi*sinphi)
    by = (br*sinphi + bphi*cosphi)

#    bdash = np.empty((nelem,3))
#    bdash[:,0] = bx/lfac + 1./(lfac+1.)*ux*(bx*ux+by*uy+bz*uz)/c/c/lfac
#    bdash[:,1] = by/lfac + 1./(lfac+1.)*uy*(bx*ux+by*uy+bz*uz)/c/c/lfac
#    bdash[:,2] = bz/lfac + 1./(lfac+1.)*uz*(bx*ux+by*uy+bz*uz)/c/c/lfac


# bdash perp:
#    bdashp = np.empty(nelem)
#    tmp    = np.empty((nelem,3))
#    tmp    = np.cross(ndash,bdash)
#    bdashp = norm(tmp)

    b       = np.empty((nelem,3))
    b[:,0] = bx
    b[:,1] = by
    b[:,2] = bz
    bdashp = 1./lfac* np.sqrt((bx**2+by**2+bz**2)-D**2*scalar(b,n)**2+2.*lfac*D*scalar(b,n)*scalar(beta,b))

    print 'got bdashp min and max: ', bdashp.min(), bdashp.max()


    if (recipe == 1 or recipe == 2):
# electron densities and cutoff:
        eps_inf = read.extract(data.get('pointdata'),'eps_inf',attribute_mode='point')
        eps_inf = np.ma.filled(np.ma.masked_array(eps_inf,m),0.)
        eps_inf = ontoGrid(
            eps_inf,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)

# emissivity:
    if recipe ==1:
        rhoe = ontoGrid(
            np.ma.filled(
                np.ma.masked_array(
                    read.extract(data.get('pointdata'),'rhoe',attribute_mode='point'),m),
                0.)
            ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)
        rhoe_0=ontoGrid(
            np.ma.filled(
                np.ma.masked_array(
            read.extract(data.get('pointdata'),'rhoe_0',attribute_mode='point'),m),
                0.)
            ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)

    elif recipe ==2:
        rhoe = ontoGrid(
            np.ma.filled(
                np.ma.masked_array(
                    read.extract(data.get('pointdata'),'n',attribute_mode='point'),m),
                0.)
            ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)
        rhoe_0=ontoGrid(
            np.ma.filled(
                np.ma.masked_array(
                    read.extract(data.get('pointdata'),'n0',attribute_mode='point'),m),
                0.)
            ,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)

    if recipe == 1 or recipe ==2:
        emiss = (c2*rhoe_0/me*c2*powerorzero(c1,(Gamma-3.)/2.)/(8.*np.pi) *
                 powerorzero(D,2.+(Gamma-1)/2.) *
                 powerorzero(bdashp,(Gamma+1.)/2.)*
                 powerorzero(rhoe_0/rhoe,-((Gamma+2.)/3.)) *
                 powerorzero(nu,(1.-Gamma)/2.) *
                 powerorzero((1. - np.sqrt(nu)/(np.sqrt(c1*bdashp)*eps_inf)),Gamma-2.)
             )
    if recipe == 3:
        emiss = powerorzero(D,2.+(Gamma-1)/2.)*powerorzero(bdashp,(Gamma+1.)/2.)

    if recipe == 4:
        rho = read.extract(data.get('pointdata'),'rho',attribute_mode='point')
        rho = np.ma.filled(np.ma.masked_array(rho,m),0.)
        rho = ontoGrid(
            rho,theta,npixel,L,points,delta=delta,x0=x0,dog=dog)
        emiss = rho * powerorzero(D,2.+(Gamma-1)/2.) * powerorzero(bdashp,(Gamma+1.)/2.)


    if polarized !=1:
        return {'emiss':emiss,'points':points,
                'phi': phi, 'theta': theta, 'delta': delta, 'nu': nu, 'alpha': alpha, 'recipe': recipe,
                'L':L,'npixel':npixel,'x0':x0}

# direction of emission:
    q       = np.empty((nelem,3))
    absb    = np.empty(nelem)
    ehat    = np.empty((nelem,3))
    absehat = np.empty((nelem))
    lhat    = np.empty((nelem,3))
    coschie = np.empty(nelem)
    sinchie = np.empty(nelem)
    cos2chie= np.empty(nelem)
    sin2chie= np.empty(nelem)
     

    q = b + np.cross(n,np.cross(beta,b))

    absb = norm(b)

    ehat      = np.cross(n,q) 
    absehat   = norm(ehat)
    ehat[:,0] = ehat[:,0]/absehat
    ehat[:,1] = ehat[:,1]/absehat
    ehat[:,2] = ehat[:,2]/absehat

    lhat1 = -costheta*cosphi_n
    lhat2 = -costheta*sinphi_n
    lhat3 = sintheta

    if delta!=0:
        cosa = np.cos(delta*np.pi/180.)
        sina = np.sin(delta*np.pi/180.)

        n1 = sintheta*cosphi_n
        n2 = sintheta*sinphi_n
        n3 = costheta

        omca  = 1.-cosa
        lhatd1 = (n1**2*omca+cosa)*lhat1 + (n1*n2*omca-n3*sina)*lhat2 + (n1*n3*omca+n2*sina)*lhat3
        lhatd2 = (n2*n1*omca+n3*sina)*lhat1 + (n2**2*omca+cosa)*lhat2 + (n2*n3*omca-n1*sina)*lhat3
        lhatd3 = (n3*n1*omca-n2*sina)*lhat1 + (n3*n2*omca+n1*sina)*lhat2 + (n3**2*omca+cosa)*lhat3

        lhat1=lhatd1
        lhat2=lhatd2
        lhat3=lhatd3
        
    lhat[:,0] = lhat1
    lhat[:,1] = lhat2
    lhat[:,2] = lhat3


    coschie = scalar(lhat,ehat)
    sinchie = scalar(n,np.cross(lhat,ehat))

# complex phase with Doppelwinkelsatz:
    sin2chie = 2.*sinchie*coschie
    cos2chie = np.square(coschie)-np.square(sinchie)

# assemble complex polarized emissivity:
    piem = (Gamma+1)/(Gamma+7./3.)
    emisspol = np.ma.filled(np.ma.masked_array(piem * emiss * (cos2chie + 1j* sin2chie),absb==0),0.)

    print '=== Done with emissivities! ==='
    return {'emiss': emiss, 'emisspol': emisspol, 'points': data.get('points'),
            'phi': phi, 'theta': theta, 'delta': delta, 'nu': nu, 'alpha': alpha, 'recipe': recipe,
                'L':L,'npixel':npixel,'x0':x0}



def extract_masked(data,var,mask,maskval=0.,attribute_mode='point'):

    mydata=np.ma.filled(
        np.ma.masked_array(
            read.extract(
                data.get('pointdata'),var,attribute_mode=attribute_mode)
            ,mask),
        maskval)

    return mydata


def ontoGrid(var,theta,npixel,L,points,delta=0,x0=[0,0,0],dog=None):

    print '=== Obtaining variable on raygrid ==='

    phi=0
    x=points[:,0]
    y=points[:,1]
    grid=make_grid(phi,theta,npixel,L,delta=delta,x0=x0)

    # Resolution of raygrid:
    drg = np.min((L[0]/np.double(npixel[0]), L[1]/np.double(npixel[1]), L[2]/np.double(npixel[2])))

    # Resolution of original grid:
    Lx = x.max() - x.min()
    Ly = y.max() - y.min()
    if dog==None:
        nelem = len(var)
        Lx = x.max() - x.min()
        Ly = y.max() - y.min()
        ny = np.int(np.sqrt(nelem * Ly/Lx))
        nx = np.int(nelem/ny)
        dog = np.min((Ly/ny,Lx/nx))
    else:
        nx = Lx/dog
        ny = Ly/dog

    xi,yi = np.mgrid[0:x.max():nx*1j, 
                     y.min():y.max():ny*1j]

#    varGrid = griddata((x,y),var,(xi,yi),method='nearest')
#    varGridSmooth = gaussian_filter(varGrid,drg/dog,mode='nearest')

    # Regrid the smoothed version onto grid with resolution drg (might not be necessary):
#    nx = np.int(Lx / drg)
#    ny = np.int(Ly / drg)
#    xiC,yiC = np.mgrid[0:x.max():nx*1j, 
#                     y.min():y.max():ny*1j]
#    varCoarseGrid = griddata((xi.ravel(),yi.ravel()),varGridSmooth.ravel(),(xiC,yiC),method='nearest')
    
#    data = griddata((xi.ravel(),yi.ravel()),varGridSmooth.ravel(),(np.sqrt(grid[0]**2+grid[1]**2),grid[2]),method='nearest')
    data = griddata((x,y),var,(np.sqrt(grid[0]**2+grid[1]**2),grid[2]),method='nearest')

    return data

def raygun(inputdata):

    print '=== Calculating radiation map... ==='
    data     = inputdata.get('emiss')
    npixel   = inputdata.get('npixel')
    L        = inputdata.get('L')

    phi      = inputdata.get('phi')
    theta    = inputdata.get('theta')
    delta    = inputdata.get('delta')

    map    = np.zeros((npixel[0],npixel[1]))
    if 'emisspol' in inputdata:
        mapp_r = np.zeros((npixel[0],npixel[1]))
        mapp_i = np.zeros((npixel[0],npixel[1]))
        datap_r = inputdata.get('emisspol').real
        datap_i = inputdata.get('emisspol').imag
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                map[i,j] = map[i,j] + data[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
                mapp_r[i,j] = mapp_r[i,j] + datap_r[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
                mapp_i[i,j] = mapp_i[i,j] + datap_i[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
    else:
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                map[i,j] = map[i,j] + data[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()

    if 'emisspol' in inputdata:
        return {'Inu': map*L[2]/npixel[2],
                'Inupol': (mapp_r + 1j * mapp_i) * L[2]/npixel[2],
                'phi': phi, 'theta': theta, 'delta': delta, 
                'nu': inputdata.get('nu'), 'alpha': inputdata.get('alpha'), 'recipe': inputdata.get('recipe'), 
                'npixel': npixel, 'L': L}
    else:
        return {'Inu': map*L[2]/npixel[2],
                'phi': phi, 'theta': theta, 'delta': delta, 
                'nu': inputdata.get('nu'), 'alpha': inputdata.get('alpha'), 'recipe': inputdata.get('recipe'), 
                'npixel': npixel, 'L': L}

def cutEmiss(data,grid,npixel,points,coord):
    
    if coord == 'cartesian':
        xrange=[points[:,0].min(),points[:,0].max()]
        yrange=[points[:,1].min(),points[:,1].max()]
        zrange=[points[:,2].min(),points[:,2].max()]
        data[grid[0]<xrange[0]] = 0.
        data[grid[0]>xrange[1]] = 0.
        data[grid[1]<yrange[0]] = 0.
        data[grid[1]>yrange[1]] = 0.
        data[grid[2]<zrange[0]] = 0.
        data[grid[2]>zrange[1]] = 0.
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                for k in range(npixel[2]):
                    ii=j*npixel[2]+i*npixel[2]*npixel[1]+k
                    if data[ii] != 0. and data[ii+1] ==0.:
                       data[ii+1:(j+1)*npixel[2]+i*npixel[2]*npixel[1]-1]=0.
    elif coord == 'spherical':
        points_r     = np.sqrt(points[:,0]**2+points[:,1]**2+points[:,2]**2)
        points_theta = np.arccos(points[:,2]/points_r)
        points_phi   = np.arctan2(points[:,1],points[:,0])
        points_phi[points_phi<0.] = 2.*np.pi + points_phi[points_phi<0.]
        
        r_range=[points_r.min(),points_r.max()]
        theta_range=[points_theta.min(),points_theta.max()]
        phi_range=[points_phi.min(),points_phi.max()]

        r     = np.sqrt(grid[0]**2+grid[1]**2+grid[2]**2)
        theta = np.arccos(grid[2]/r)
        phi   = np.arctan2(grid[1],grid[0])
        phi[phi<0.] = 2.*np.pi + phi[phi<0.]

        data[r<r_range[0]]          = 0.
        data[r>r_range[1]]          = 0.
        data[theta<theta_range[0]]  = 0.
        data[theta>theta_range[1]]  = 0.
        data[phi<phi_range[0]]      = 0.
        data[phi>phi_range[1]]      = 0.
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                for k in range(npixel[2]):
                    ii=j*npixel[2]+i*npixel[2]*npixel[1]+k
                    if data[ii] != 0. and data[ii+1] ==0.:
                       data[ii+1:(j+1)*npixel[2]+i*npixel[2]*npixel[1]-1]=0.
    else:
        print 'COORDINATES', coord, 'not implemented yet'

    return data

def shotgun(inputdata,npixel,L,x0=None,coordinate='cartesian',nocut=True):

    print '=== Calculating radiation map... ==='
    points   = inputdata.get('points')
    emiss    = inputdata.get('emiss')
    if 'emisspol' in inputdata:
        emisspol = inputdata.get('emisspol')
    phi      = inputdata.get('phi')
    theta    = inputdata.get('theta')
    delta    = inputdata.get('delta')
    x=points[:,0]
    y=points[:,1]
    z=points[:,2]
    if x0 == None:
        x0=[(x.min()+x.max())/2.,(y.min()+y.max())/2.,(z.min()+z.max())/2.] 

    grid = make_grid(phi,theta,npixel,L,delta=delta,x0=x0)

    data    = griddata((x,y,z),emiss, (grid[0],grid[1],grid[2]),method='nearest')
    
# remove emissivity outside of original data:
    if nocut==False:
        data = cutEmiss(data,grid,npixel,points,coordinate)

    if 'emisspol' in inputdata:
        datap_r = griddata((x,y,z),emisspol.real, (grid[0],grid[1],grid[2]),method='nearest')
        datap_i = griddata((x,y,z),emisspol.imag, (grid[0],grid[1],grid[2]),method='nearest')
    print '=== Done with regridding, now integrating los ==='
    map    = np.zeros((npixel[0],npixel[1]))
    
    if 'emisspol' in inputdata:
        mapp_r = np.zeros((npixel[0],npixel[1]))
        mapp_i = np.zeros((npixel[0],npixel[1]))
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                map[i,j] = map[i,j] + data[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
                mapp_r[i,j] = mapp_r[i,j] + datap_r[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
                mapp_i[i,j] = mapp_i[i,j] + datap_i[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()
    else:
        for i in range(npixel[0]):
            for j in range(npixel[1]):
                map[i,j] = map[i,j] + data[np.arange(npixel[2])+j*npixel[2]+i*npixel[2]*npixel[1]].sum()

    print '=== Done with map! ==='
    if 'emisspol' in inputdata:
        return {'Inu': map*L[2]/npixel[2],'Inupol': (mapp_r + 1j * mapp_i) * L[2]/npixel[2],
                'phi': phi, 'theta': theta, 'delta': delta, 
            'nu': inputdata.get('nu'), 'alpha': inputdata.get('alpha'), 'recipe': inputdata.get('recipe'), 
            'npixel': npixel, 'L': L}
    else:
        return {'Inu': map*L[2]/npixel[2],
                'phi': phi, 'theta': theta, 'delta': delta, 
            'nu': inputdata.get('nu'), 'alpha': inputdata.get('alpha'), 'recipe': inputdata.get('recipe'), 
            'npixel': npixel, 'L': L}

def shoot(inputdata,phi,theta,npixel,xgrid,ygrid):

# Euclidian norm:
    norm = lambda x: np.sqrt(np.square(x).sum())

    points = inputdata.get('points')
    emiss  = inputdata.get('emiss')


    dy = (xgrid[1]-xgrid[0])/npixel[0]
    dz = dy*np.sqrt(2.)

    beg = np.array([-1.e18,-1.e18,-2.e18])
    end = np.array([1.e18,-1.e18,0.])

    map = np.empty((npixel[0],npixel[1]))

    for i in range(npixel[0]):
        for j in range(npixel[1]):
            map[i,j] = get_ray(emiss,points,512,beg,end).sum()*norm(end-beg)/512.
            print 'got i,j:', i,j, i*npixel[0]+j,'map:',map[i,j]
            beg[1] = beg[1]+dy
            end[1] = end[1]+dy
        beg[1] = beg[1]-npixel[1]*dy
        end[1] = end[1]-npixel[1]*dy
        beg[2] = beg[2]+dz
        end[2] = end[2]+dz
        

    return map

def phase(map):
    '''map should be a 2D array of complex numbers'''
    phase = np.empty((map.shape[0],map.shape[1]))

    x = map.real
    y = map.imag

    for i in range(map[0].shape[0]):
        for j in range(map[0].shape[1]):
            if x[i,j] > 0.:
                phase[i,j] = np.arctan(y[i,j]/x[i,j])
            elif x[i,j] < 0. and y[i,j] >= 0.:
                phase[i,j] = np.arctan(y[i,j]/x[i,j]) + np.pi
            elif x[i,j] < 0. and y[i,j] < 0.:
                phase[i,j] = np.arctan(y[i,j]/x[i,j]) - np.pi
            elif x[i,j] == 0. and y[i,j] > 0.:
                 phase[i,j] = np.pi/2.
            elif x[i,j] == 0. and y[i,j] < 0.:
                phase[i,j] = -np.pi/2.
            else:
                phase[i,j] = float('nan') 

    return phase

def beam(input,sigma):
    
    nsig = input.get('npixel')[0]*sigma/(input.get('L')[0])
    map = mapshow = np.copy(input.get('Inu'))
    mapSmooth = gaussian_filter(map,nsig,mode='nearest')

    output = np.copy(input)
    output['Inu'] = mapSmooth

    return output

def make_grid(phi,theta,npixel,L,delta=0,x0=[0,0,0]):

    xi,yi,zi=np.mgrid[0:npixel[0],0:npixel[1],0:npixel[2]]

    xi=xi/(npixel[0]-1.)
    yi=yi/(npixel[1]-1.)
    zi=zi/(npixel[2]-1.)

    xi = xi*L[0] - L[0]/2.
    yi = yi*L[1] - L[1]/2.
    zi = zi*L[2] - L[2]/2.

    xi = xi.flatten()
    yi = yi.flatten()
    zi = zi.flatten()

    nelem = len(xi)

# Since the standard phi counts from x (and not -y as we do by rotating around x first)
# add 90 to the internal phi.  
    phii = phi + 90.
# Now rotate the points around theta and then phi:

    cosphi = np.cos(phii*np.pi/180.)
    sinphi = np.sin(phii*np.pi/180.)
    costheta = np.cos(theta*np.pi/180.)
    sintheta = np.sin(theta*np.pi/180.)

    xid = xi
    yid = np.empty(nelem)
    zid = np.empty(nelem)
# Rotate by theta around the x-axis:
    yid = costheta*yi - sintheta*zi
    zid = sintheta*yi + costheta*zi

    xidd = np.empty(nelem)
    yidd = np.empty(nelem)
    zidd = zid
# Rotate by phi around the z-axis:
    xidd =  cosphi*xid - sinphi*yid
    yidd =  sinphi*xid + cosphi*yid
# Optional: rotate around delta to orient view:
    if delta==0:
        return [xidd+x0[0],yidd+x0[1],zidd+x0[2]]
    else:
        cosa = np.cos(delta*np.pi/180.)
        sina = np.sin(delta*np.pi/180.)
        # make sure to use the original phi angle: 
        n1 = sintheta*np.cos(phi*np.pi/180.)
        n2 = sintheta*np.sin(phi*np.pi/180.)
        n3 = costheta
        
        xiddd = np.empty(nelem)
        yiddd = np.empty(nelem)
        ziddd = np.empty(nelem)
        
        omca  = 1.-cosa
        xiddd = (n1**2*omca+cosa)*xidd + (n1*n2*omca-n3*sina)*yidd + (n1*n3*omca+n2*sina)*zidd
        yiddd = (n2*n1*omca+n3*sina)*xidd + (n2**2*omca+cosa)*yidd + (n2*n3*omca-n1*sina)*zidd
        ziddd = (n3*n1*omca-n2*sina)*xidd + (n3*n2*omca+n1*sina)*yidd + (n3**2*omca+cosa)*zidd
        return [xiddd+x0[0],yiddd+x0[1],ziddd+x0[2]]

def show_map(inputdata,npol='none',filename='none',function=None,background=None,vmin=None,vmax=None,colpol='w',cmap=cm.gist_heat_r,ncbarticks=5,orientation='horizontal'):

    dpi=300
    fontsize=12
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'


    L   = inputdata.get('L')
# square root filter:
    map = inputdata.get('Inu')
    if function == None:
        mapshow = np.copy(inputdata.get('Inu'))
    else:
        mapshow = function(np.copy(inputdata.get('Inu')))


    if 'Inupol' in inputdata:
        mappol = inputdata.get('Inupol')

    fig1 = plt.figure(figsize=(4,4))
    ax  = fig1.add_subplot(1,1,1)
    ax.set_xlim(-L[0]/2.,L[0]/2.)
    ax.set_ylim(-L[1]/2.,L[1]/2.)

    extent=(-L[0]/2.,L[0]/2.,-L[1]/2.,L[1]/2.)

    palette = cmap
    if background != None:
        palette.set_over(background, 1.0)
        palette.set_under(background, 1.0)
        palette.set_bad(background, 1.0)

    if vmin==None:
        vmin=mapshow.min()
    if vmax==None:
        vmax=mapshow.max()
    if background == None:
        map_clip = np.clip(mapshow,vmin,vmax)
    else:
        map_clip = np.ma.masked_where(mapshow<vmin, mapshow)
        map_clip = np.ma.masked_where(mapshow>vmax, map_clip)

    normColors=mpl.colors.Normalize(vmin=vmin,vmax=vmax,clip = False)

    cs = ax.imshow(map_clip.transpose(),origin='lower',cmap=palette,extent=extent,norm=normColors)
#    ax.patch.set_facecolor(background)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

#    for tick in cb.ax.get_yticklabels():
#        tick.set_fontsize(fontsize)

    plt.plot(0,0,'w+',markersize=10)

    if npol != 'none' and 'Inupol' in inputdata:
        xbeg=0.5*L[0]*float(npol[0]-1)/float(npol[0])
        ybeg=0.5*L[1]*float(npol[1]-1)/float(npol[1])
        x=np.linspace(-xbeg,xbeg,npol[0])
        y=np.linspace(-ybeg,ybeg,npol[1])
        x2D=np.empty((npol[0],npol[1]))
        y2D=np.empty((npol[0],npol[1]))
        for i in range(npol[1]):
            x2D[:,i] = x
        for i in range(npol[0]):
            y2D[i,:] = y
#        Q=congrid(mappol.real,npol,method='linear',centre=True)
#        U=congrid(mappol.imag,npol,method='linear',centre=True)
#        evec = get_polvector(Q,U)

        ex=-congrid(np.sqrt(mappol).imag,npol,method='linear',centre=True)
        ey=congrid(np.sqrt(mappol).real,npol,method='linear',centre=True) 
        pi=congrid(abs(mappol)/map,npol,method='linear',centre=True)
        
        norm=np.sqrt(ex**2+ey**2)
        ex=ex/norm * pi
        ey=ey/norm * pi
# evectors:
#        plt.quiver(x2D,y2D,ex,ey,width=0.002,headlength=0,headwidth=0,headaxislength=0,pivot='middle',color=colpol, edgecolors=(colpol))
# bvectors:
        plt.quiver(x2D,y2D,-ey,ex,width=0.002,headlength=0,headwidth=0,headaxislength=0,pivot='middle',color=colpol,scale=npol[0]/4.,scale_units='inches', edgecolors=(colpol))
        
# make a little inset showing the scaling of the polarization
        alpha = inputdata.get('alpha')
        pimax = (alpha+1.)/(alpha+5./3.)
#        if colpol == 'w':
#            plt.fill((-7.5e16,-6e16,-6e16,-7.5e16),(-7.5e16,-7.5e16,-6e16,-6e16),'k')
#        else:
#            plt.fill((-7.5e16,-6e16,-6e16,-7.5e16),(-7.5e16,-7.5e16,-6e16,-6e16),'w')            
#        plt.quiver(-6.75e16,-6.75e16,pimax,0,width=0.002,headlength=0,headwidth=0,headaxislength=0,pivot='middle',color=colpol,zorder=2,scale_units='inches',scale=npol[0]/4)


# Now for the colorbar (plotted separately):
#    figcb = pyplot.figure(figsize=(4,0.25))
#    axcb = figcb.add_subplot(1,1,1)
#    cbarticks = np.linspace(vmin, vmax, ncbarticks, endpoint=True)
#    if vmax == mapshow.max():
#        cb = mpl.colorbar.ColorbarBase(axcb, cmap=palette,
#                                       norm=normColors,
#                                       ticks=cbarticks,
#                                    orientation='horizontal')
#    else:
#        cb = mpl.colorbar.ColorbarBase(axcb, cmap=palette,
#                                       norm=normColors,extend='max',
#                                       ticks=cbarticks,
#                                    orientation='horizontal')
#        
#    for tick in cb.ax.get_yticklabels():
#        tick.set_fontsize(fontsize)
#    cb.ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
#    cb.ax.yaxis.get_offset_text().set_fontsize(fontsize-2)
#    cb.solids.set_rasterized(True) 


    if filename == 'none':
        plt.show()
    else: 
        fig1.savefig(''.join([filename,'.pdf']), format="pdf", transparent=True,dpi=dpi,bbox_inches='tight')
        figcb.savefig(''.join([filename,'_colorbar','.pdf']), format="pdf",
                                   transparent=True,dpi=dpi,bbox_inches='tight')
        plt.close()


def show_polmap(inputdata,npol='none',filename='none',background='w',vmin=None,vmax=None,colpol='w',cmap=cm.gist_heat_r):

    dpi=300
    fontsize=12
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'


    L   = inputdata.get('L')
    mappol = inputdata.get('Inupol')
    map    = inputdata.get('Inu')

    fig1 = plt.figure(figsize=(4,4))
    ax  = fig1.add_subplot(1,1,1)
    ax.set_xlim(-L[0]/2.,L[0]/2.)
    ax.set_ylim(-L[1]/2.,L[1]/2.)

    extent=(-L[0]/2.,L[0]/2.,-L[1]/2.,L[1]/2.)

    if vmin==None:
        vmin=map.min()
    if vmax==None:
        vmax=map.max()
    map_clip = np.clip(np.abs(mappol),vmin,vmax)
    cs = plt.imshow(map_clip.transpose(),origin='lower',cmap=cmap,extent=extent,vmin=vmin,vmax=vmax)
#    cb = plt.colorbar()
#    for tick in cb.ax.get_yticklabels():
#        tick.set_fontsize(fontsize)
#    cb.ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
#    cb.ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

#    ax.patch.set_facecolor(background)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

#    for tick in cb.ax.get_yticklabels():
#        tick.set_fontsize(fontsize)

    plt.plot(0,0,'w+',markersize=10)

#    plt.ylabel('b [cm]',size=fontsize)
#    plt.xlabel('a [cm]',size=fontsize)
#    plt.title(''.join((['phi=',str(inputdata.get('phi'))[0:4],' ; theta=',str(inputdata.get('theta'))[0:4]
#                        ,' ; delta=',str(inputdata.get('delta'))[0:4]])))

    if npol != 'none' and 'Inupol' in inputdata:
        xbeg=0.5*L[0]*float(npol[0]-1)/float(npol[0])
        ybeg=0.5*L[1]*float(npol[1]-1)/float(npol[1])
        x=np.linspace(-xbeg,xbeg,npol[0])
        y=np.linspace(-ybeg,ybeg,npol[1])
        x2D=np.empty((npol[0],npol[1]))
        y2D=np.empty((npol[0],npol[1]))
        for i in range(npol[1]):
            x2D[:,i] = x
        for i in range(npol[0]):
            y2D[i,:] = y
        ex=-congrid(np.sqrt(mappol).imag,npol,method='linear',centre=True)
        ey=congrid(np.sqrt(mappol).real,npol,method='linear',centre=True) 
        pi=congrid(abs(mappol)/map,npol,method='linear',centre=True)
        
        norm=np.sqrt(ex**2+ey**2)
        ex=ex/norm * pi
        ey=ey/norm * pi
# evectors:
#        plt.quiver(x2D,y2D,ex,ey,width=0.002,headlength=0,headwidth=0,headaxislength=0,pivot='middle',color=colpol)
# bvectors:
        plt.quiver(x2D,y2D,-ey,ex,width=0.002,headlength=0,headwidth=0,headaxislength=0,pivot='middle',color=colpol,scale=npol[0]/4.,scale_units='inches')
        
# make a little inset showing the scaling of the polarization
        alpha = inputdata.get('alpha')
        pimax = (alpha+1.)/(alpha+5./3.)
        plt.fill((-7.5e17,-6e17,-6e17,-7.5e17),(-7.5e17,-7.5e17,-6e17,-6e17),'w')
        plt.quiver(-6.75e17,-6.75e17,pimax,0,width=0.002,headlength=0,headwidth=0,headaxislength=0,pivot='middle',color=colpol,zorder=2,scale_units='inches',scale=npol[0]/4)


    if filename == 'none':
        plt.show()
    else: 
        plt.savefig(''.join([filename,'.pdf']), format="pdf", transparent=False,dpi=dpi,bbox_inches='tight')
        plt.close()



def show_pol(inputdata,npol='none',filename='none',background='w',vmin=None,vmax=None,cmap=cm.gist_heat,ncbarticks=5,orientation='horizontal'):

    dpi=300
    fontsize=12
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'


    L   = inputdata.get('L')
    pi = np.abs(inputdata.get('Inupol'))/inputdata.get('Inu')

    fig1 = plt.figure(figsize=(4,4))
    ax  = fig1.add_subplot(1,1,1)
    ax.set_xlim(-L[0]/2.,L[0]/2.)
    ax.set_ylim(-L[1]/2.,L[1]/2.)

    extent=(-L[0]/2.,L[0]/2.,-L[1]/2.,L[1]/2.)

    if vmin==None:
        vmin=pi.min()
    if vmax==None:
        vmax=pi.max()
    pi_clip = np.clip(pi,vmin,vmax)
    cs = plt.imshow(pi_clip.transpose(),origin='lower',cmap=cmap,extent=extent,vmin=vmin,vmax=vmax)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize-2)

    plt.plot(0,0,'k+',markersize=10)

# Now for the colorbar (plotted separately):
    figcb = pyplot.figure(figsize=(4,0.25))
#    axcb = figcb.add_axes([0.05, 0.80, 0.9, 0.15])
    axcb = figcb.add_subplot(1,1,1)

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cbarticks = np.linspace(vmin, vmax, ncbarticks, endpoint=True)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap,
                                   norm=norm,
                                   orientation=orientation,
                                    ticks=cbarticks)
# adjust the size of font:
    for tick in cb.ax.get_yticklabels():
        tick.set_fontsize(fontsize)
    cb.ax.xaxis.get_offset_text().set_fontsize(fontsize-2)
    cb.ax.yaxis.get_offset_text().set_fontsize(fontsize-2)


    if filename == 'none':
        plt.show()
    else: 
        fig1.savefig(''.join([filename,'.pdf']), format="pdf", transparent=False,dpi=dpi,bbox_inches='tight')
        figcb.savefig(''.join([filename,'_colorbar','.pdf']), format="pdf",
                                   transparent=False,dpi=dpi,bbox_inches='tight')
        plt.close()



def get_polvector(Q,U):
    
    evec = np.empty((Q.shape[0],Q.shape[1],2))

    chie = 0.5* np.arctan(U/Q)

    evec[Q>0,0] = - np.sin(chie[Q>0])
    evec[Q>0,1] = + np.cos(chie[Q>0])

    evec[Q<=0,0] = + np.sin(chie[Q<=0])
    evec[Q<=0,1] = + np.cos(chie[Q<=0])

    return evec

def sweep_theta(offset,phi,theta_beg,theta_end,ntheta,nu,alpha,recipe=2,filenameout='data',type='pvtu'):

    npixel=[400,400,400]
    L=[1e19,1e19,1e19]

    theta=np.linspace(theta_beg,theta_end,ntheta)

#    data = rot2D(offset,128,type=type)
    data = get_pointdata(offset,filenameout=filenameout,type=type)
    for theta_i in theta:
#        em = emiss(data,phi,theta_i,nu,alpha,polarized=0,recipe=recipe)
        em = emissFrom2D(data,theta=theta_i,L=L,npixel=npixel,nu=nu)
#        map = shotgun(em,npixel,L)
        map = raygun(em)
        dump(map,dumpfile=''.join(['map',str(offset).zfill(4),'_phi',str(phi)[0:4],'theta',str(theta_i).zfill(5),'.dat']))
        del em, map


def sweep_phi(offset,phi_beg,phi_end,theta,nphi,nu,alpha,filenameout='data',type='pvtu'):

    npixel=[512,512,512]
    L=[1e18,1e18,2e18]

    phi=np.linspace(phi_beg,phi_end,nphi)

    data = get_pointdata(offset,filenameout=filenameout,type=type)
    i = 0
    for phi_i in phi:
        em = emiss(data,phi_i,theta,nu,alpha)
        map = shotgun(em,npixel,L)
        show_map(map,filename=''.join(['movie/t',str(offset).zfill(4),'frame',str(i).zfill(4),'phi',str(phi_i)[0:4],'theta',str(theta)[0:4]]))
        show_map(map,filename=''.join(['movie/t',str(offset).zfill(4),'frame',str(i).zfill(4),'_pol_phi',str(phi_i)[0:4],'theta',str(theta)[0:4]]),npol=[32,32])
        i = i+1

def sweep_time(offset_beg,offset_end,phi,theta,nu,alpha,recipe=1,delta=0,cmap=cm.hot,filenameout='data',type='pvtu',attribute_mode='cell',vmin=None,vmax=None,
path='./',
npixel=[200,200,200],
L=[2e18,2e18,2e18],
foldername='xray'):
    os.system(''.join(['mkdir ',path,foldername,str(recipe)]))
    

    offset=np.arange(offset_end-offset_beg+1)+offset_beg

    for offset_i in offset:
        data = get_pointdata(offset_i,filenameout=filenameout,type=type,attribute_mode=attribute_mode)
        em = emiss(data,phi,theta,nu,alpha,delta=delta,recipe=recipe)
        map = shotgun(em,npixel,L)
#        show_map(map,cmap=cmap,vmin=vmin,vmax=vmax,filename=''.join([path,foldername,str(recipe),'/t',str(offset_i).zfill(4),'_phi',str(phi)[0:4],'theta',str(theta)[0:4]]))
#        show_map(map,cmap=cmap,vmin=vmin,vmax=vmax,filename=''.join([path,foldername,str(recipe),'/t',str(offset_i).zfill(4),'_pol_phi',str(phi)[0:4],'theta',str(theta)[0:4]]),npol=[32,32])
        dump(map,dumpfile=''.join([path,foldername,str(recipe),'/t',str(offset_i).zfill(4),'_phi',str(phi)[0:4],'theta',str(theta)[0:4],'.dat']))
        del data, em, map


def sweep_composite(offset_beg,offset_end,phi,theta,nu_r,nu_g,nu_b,alpha,rotate=0.,filenameout='data',type='pvtu'):

    npixel=[256,256,512]
    L=[3e18,3e18,3e18]

    offset=np.arange(offset_end-offset_beg+1)+offset_beg

    i = 0
    for offset_i in offset:
        data = get_pointdata(offset_i,filenameout=filenameout,type=type)
        em_r = emiss(data,phi,theta,nu_r,alpha,polarized=0)
        em_g = emiss(data,phi,theta,nu_g,alpha,polarized=0)
        em_b = emiss(data,phi,theta,nu_b,alpha,polarized=0)

        map_r = shotgun(em_r,npixel,L)
        map_g = shotgun(em_g,npixel,L)
        map_b = shotgun(em_b,npixel,L)

        composite(map_r,map_g,map_b,rotate=rotate,filename=''.join(['movie_composite/t',str(offset_i).zfill(4),'frame',str(i).zfill(4),'phi',str(phi)[0:4],'theta',str(theta)[0:4]]))

        i = i+1

def composite(map_r,map_g,map_b,rotate=0.0,filename='none'):
    
    img=np.zeros((map_r.get('Inu').shape[0],map_r.get('Inu').shape[1],3))

    img[:,:,0] = img_scale.sqrt(map_r.get('Inu'))
    img[:,:,1] = img_scale.sqrt(map_g.get('Inu'))
    img[:,:,2] = img_scale.sqrt(map_b.get('Inu'))

    img = scipy.ndimage.rotate(img, rotate, reshape=False,order=1) 

    fig=plt.figure(frameon=False)
    plt.imshow(img, aspect='equal',origin='lower')
    ax=fig.add_subplot(1,1,1)
    ax.set_axis_off()
    ax.set_aspect('equal')
    
    if filename == 'none':
        plt.show()
    else: 
        fig.set_size_inches(4,4)
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        plt.savefig(filename,bbox_inches=extent)
        plt.close()

def getFlux(offset,phi,theta,nu,alpha,recipe=1,nphi=None,filenameout='data',type='pvtu'):

    npixel=[128,128,128]
    L=[5e17,5e17,5e17]

    if nphi==None:
        data = get_pointdata(offset,filenameout=filenameout,type=type)
    else:
        data = rot2D(offset,nphi,filenameout=filenameout,type=type)
    em   = emiss(data,phi,theta,nu,alpha,recipe=recipe,polarized=0)
    map  = shotgun(em,npixel,L)

    flux = map.get('Inu').sum()

    return flux


def dump(data,dumpfile='dumpfile.dat'):

    file=open(dumpfile,'wb') 
    pickle.dump(data,file, pickle.HIGHEST_PROTOCOL)
    file.close()


def restore(dumpfile='dumpfile.dat'):

    file=open(dumpfile,'rb')
    data = pickle.load(file)
    file.close()
    return data
