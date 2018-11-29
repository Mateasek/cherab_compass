from cherab.compass.equilibrium.equilibrium import COMPASSEquilibrium
from raysect.primitive.mesh.mesh import Mesh
from raysect.primitive.mesh.stl import export_stl
from raysect.primitive import Cylinder, Sphere, export_vtk, Box
from raysect.optical.material import AbsorbingSurface
from raysect.optical.material.emitter import UnityVolumeEmitter, UnitySurfaceEmitter, UniformVolumeEmitter
from raysect.optical.spectralfunction import ConstantSF
from raysect.optical.spectrum import Spectrum

from raysect.core.math.transform import translate 
from raysect.primitive.csg import Subtract
from raysect.optical.observer import MeshCamera, MeshPixel, PowerPipeline0D, PowerPipeline1D, MonoAdaptiveSampler1D
from raysect.optical import World
from raysect.core import Point3D
import numpy as np
import matplotlib.pyplot as plt
from raysect.core import SerialEngine

def makemesh(polygon, n_segments, torangle, closed=False):
    polygon_len = polygon.shape[0]
    vertices = []
    triangles = []

    for i in range(n_segments):
        for j in range(polygon_len):
            vertices.append([polygon[j,0] * np.cos(i * torangle / n_segments),
                            polygon[j,0] * np.sin(i * torangle / n_segments),
                            polygon[j,1]])

    n_vertices = len(vertices)
    close_surf = True

    if close_surf:
        segrange = n_segments
    else:
        segrange = n_segments-1


    for i in range(segrange):
        for j in range(polygon_len):
            s0 = i * polygon_len
            s1 = (s0+polygon_len)%n_vertices
            triangles.append([s0 + j, (s1 + j)%n_vertices, s1+(j+1)%polygon_len])
            triangles.append([s1 + (j + 1)%polygon_len, s0 + (j + 1)%polygon_len, s0 + j])
            
    return vertices, triangles
  
 
 #get equilibrium
path = "data/efit_17636.h5"
shot_time = 1.125

# get equilibrium
equilibrium = COMPASSEquilibrium(path=path)
equilibrium_slice = equilibrium.time(shot_time)

world=World()

n_segments = 180
torangle = np.pi/2

limiter = equilibrium_slice.limiter_polygon
vertices, triangles = makemesh(limiter, n_segments, torangle, closed = True)

pfcs = Mesh(vertices=vertices, triangles=triangles, parent = world, material=AbsorbingSurface(), closed=False)
#export_stl(pfcs,"pfc.stl")

a = 0.25/(np.sqrt(2*np.pi*0.56))
cylinder1 = Cylinder(radius= 0.56 + a/2, height = a,
                     transform =translate(0, 0, -a/2))

cylinder2 = Cylinder(radius= 0.56 - a/2, height = 1, transform =translate(0, 0, -0.5))

#box1 = Box(lower=Point3D(-1.5, -1.5, -1.5), upper=Point3D(1.5,0,1.5))
#box2 = Box(lower=Point3D(-1.5, -1.5, -1.5), upper=Point3D(0, 1.5, 1.5))
box1 = Box(lower=Point3D(-1.5,0,-1), upper=Point3D(1.5, 0, 1))
box2 = Box(lower=Point3D(0,-1,-1), upper=Point3D(0, 1, 1))


#box = Box(lower=Point3D(0, 0, 0), upper=Point3D(1,1,1))
radiation = UniformVolumeEmitter(emission_spectrum=ConstantSF(1.0), scale=1.5e5)
radiator1 = Subtract(cylinder1, cylinder2)
radiator2 = Subtract(radiator1, box1)
radiatior3 = Subtract(radiator2, box2, parent = world, material=UnityVolumeEmitter())

min_wl = 200
max_wl = 1000
samples = 100000

power = PowerPipeline1D()
sampler = MonoAdaptiveSampler1D(power, fraction=0.2, ratio=25.0, min_samples=1000, cutoff=0.1)
camera = MeshCamera(
    pfcs,
    surface_offset=1e-6,  # launch rays 1mm off surface to avoid intersection with absorbing mesh
    pipelines=[power],
    frame_sampler=sampler,
    parent=world,
    spectral_bins=1,
    min_wavelength=650,
    max_wavelength=651,
    pixel_samples=250
)
camera.render_engine = SerialEngine()


camera.observe()
