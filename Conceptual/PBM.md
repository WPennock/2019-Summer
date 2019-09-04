# Model-Based PBM

## Continuous, different collision times
Take time step, for each bin, calculate collision time and divide time step by
collision time. Multiply this fraction by current population in bin and add half
this to next bin.
## Imports
```python
from aguaclara.play import *
import aguaclara.research.floc_model as floc
import pdb
import pickle as pkl
import matplotlib.animation as anim
import types
import html
from numpngw import AnimatedPNGWriter

plt.rcParams["animation.convert_path"] = "C:\Program Files\ImageMagick-7.0.6-Q16\magick.exe"
plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]
params = {'text.usetex' : True,
          'font.size' : 10,
          'font.family' : 'serif',
          'text.latex.unicode': True,
          'axes.facecolor': 'white',
          'axes.labelcolor': 'black',
          'savefig.facecolor': 'white',
          'axes.edgecolor': 'black',
          'savefig.edgecolor': 'black',
          'xtick.color': 'black',
          'ytick.color': 'black',
          'grid.linestyle': 'None',
          #'lines.markersize': 1
          }
plt.rcParams.update(params)
```

### Load Data
```python
PSD = pkl.load(open('PSD.pkl','rb'))*90*u.NTU
Bins = pkl.load(open('Bins.pkl','rb'))
```

```python
# Figuring out bins
Dose = 2.65*u.mg/u.L
d_Max = 100*u.um
d_0 = 7*u.um
floc.DIM_FRACTAL
NumP = (d_Max/d_0)**(floc.DIM_FRACTAL)
NBin = np.ceil(np.log2(NumP))
NBin
Coll = np.arange(0,NBin+1)
Coll
BinD = floc.diam_fractal(floc.DIM_FRACTAL,d_0,Coll)
BinD
C_0 = 90*u.NTU
BinC = np.zeros(len(BinD))*u.NTU
BinC[0] = C_0
BinC_0 = BinC

def frac_vol_floc_initial(ConcAluminum, ConcClay, coag, material):
    """Return the volume fraction of flocs initially present, accounting for both suspended particles and coagulant precipitates.
    :param ConcAluminum: Concentration of aluminum in solution
    :type ConcAluminum: float
    :param ConcClay: Concentration of particle in suspension
    :type ConcClay: float
    :param coag: Type of coagulant in solution
    :type coag: float
    :param material: Type of particles in suspension, e.g. floc_model.Clay
    :type material: floc_model.Material
    :return: Volume fraction of particles initially present
    :rtype: float
    """
    return ((floc.conc_precipitate(ConcAluminum, coag)/coag.PrecipDensity) + (ConcClay / material.Density)).to(u.dimensionless)

phi_0 = frac_vol_floc_initial(Dose,C_0,floc.PACl,floc.Clay)
def gamma_coag(ConcClay, ConcAluminum, coag, material,
               DiamTube, RatioHeightDiameter):
    """Return the coverage of clay with nanoglobs.
    This function accounts for loss to the tube flocculator walls
    and a poisson distribution on the clay given random hits by the
    nanoglobs. The poisson distribution results in the coverage only
    gradually approaching full coverage as coagulant dose increases.
    :param ConcClay: Concentration of clay in suspension
    :type ConcClay: float
    :param ConcAluminum: Concentration of aluminum in solution
    :type ConcAluminum: float
    :param coag: Type of coagulant in solution, e.g. floc_model.PACl
    :type coag: floc_model.Material
    :param material: Type of clay in suspension, e.g. floc_model.Clay
    :type material: floc_model.Material
    :param DiamTube: Diameter of flocculator tube (assumes tube flocculator for calculation of reactor surface area)
    :type DiamTube: float
    :param RatioHeightDiameter: Dimensionless ratio of clay height to clay diameter
    :type RatioHeightDiameter: float
    :return: Fraction of the clay surface area that is coated with coagulant precipitates
    :rtype: float
    """
    return (1 - np.exp((
                       (-frac_vol_floc_initial(ConcAluminum, 0*u.kg/u.m**3, coag, material)
                         * material.Diameter)
                        / (frac_vol_floc_initial(0*u.kg/u.m**3, ConcClay, coag, material)
                           * coag.Diameter))
                       * (1 / np.pi)
                       * (floc.ratio_area_clay_total(ConcClay, material,
                                                DiamTube, RatioHeightDiameter)
                          / floc.ratio_clay_sphere(RatioHeightDiameter))
))
Gamma = gamma_coag(C_0,Dose,floc.PACl,floc.Clay,1.25*u.inch,floc.RATIO_HEIGHT_DIAM)
Gamma

def time_col_laminar(EnergyDis, Temp, ConcAl, ConcClay, coag, material,
                     DiamTarget, DiamTube, DIM_FRACTAL, RatioHeightDiameter):
    """Calculate single collision time for laminar flow mediated collisions.
    Calculated as a function of floc size.
    """
    return (((1/6) * ((6/np.pi)**(1/3))
             * frac_vol_floc_initial(ConcAl, ConcClay, coag, material) ** (-2/3)
             * (pc.viscosity_kinematic(Temp) / EnergyDis) ** (1 / 2)
             * (DiamTarget / material.Diameter) ** (2*DIM_FRACTAL/3 - 2)
             )  # End of the numerator
            / (gamma_coag(ConcClay, ConcAl, coag, material, DiamTube,
                          RatioHeightDiameter)
               )  # End of the denominator
            ).to(u.s)

time_col_laminar(25*u.mW/u.kg, 22*u.degC, Dose, C_0, floc.PACl, floc.Clay, 10*u.um, 1.25*u.inch, floc.DIM_FRACTAL, 0.1)

alpha = 2*Gamma-Gamma**2*frac_vol_floc_initial(Dose,C_0,floc.PACl,floc.Clay)
def t_c(Dose,conc,coag,material,EDR,temp):
    return ((np.pi*6/frac_vol_floc_initial(Dose,C_0,coag,material))**(2/3)/np.pi*(pc.viscosity_kinematic(temp)/EDR)**(1/2)).to(u.s)
t_c(Dose,C_0,floc.PACl,floc.Clay,25*u.mW/u.kg,22*u.degC)

t = 0
T = 387
dt = 0.1
(T/dt)

# Set up GIF
ims = []
ts = []

while t<=T:
    for i in range(0,len(BinD)-1):
        tc = t_c(Dose,BinC[i],floc.PACl,floc.Clay,24*u.mW/u.kg,22*u.degC).magnitude
        dC = alpha*dt/tc*BinC[i]
        BinC[i] = BinC[i] - dC
        BinC[i+1] = BinC[i+1] + dC
    if round(t,1)%10 == 0:
        ims.append(BinC.copy())    
        ts.append(t)
    t = t + dt

pkl.dump(ims,open("PSD_265.pkl","wb"))
pkl.dump(ts,open("t_265.pkl","wb"))

# Create Animation function
plt.clf(),plt.close('all')
fig, ax = plt.subplots()

# Plot data graphs
ind = 5
width1 = np.zeros(len(Bins[ind]))
width2 = np.zeros(len(Bins[ind]))
len(Bins[0])
for j in range(0,len(Bins[ind])-1):
    width1[j] = Bins[ind][j+1]-Bins[ind][j]
    width1[6] = 20
    width2[j] = Bins[ind+6][j+1]-Bins[ind+6][j]
    width2[6] = 20

ax.bar(Bins[ind],PSD[ind],align='edge',width=width1,color='b',alpha=0.5,label=r"1st 2.65 mg/L")
ax.bar(Bins[ind+6],PSD[ind+6],align='edge',width=width2,color='r',alpha=0.5,label=r"2nd 2.65 mg/L")
ax.xlabel(r'')
ax.xlabel(r'')

# ax.axis([0,100,0,90])

# Plot simulation graphs
width = np.zeros(len(BinD))
len(BinD)
for i in range(0,len(BinD)-1):
    width[i+1] = -(BinD[i+1] - BinD[i]).magnitude
    width[0] = -BinD[0].magnitude
width    
Graph = ax.bar(BinD,BinC_0,color='g',alpha = 0.5,width=width,align='edge',label='Simulation')
label = ax.text(0.42,0.95,'t = 0 s',transform=ax.transAxes)
plt.bar(BinD,BinC_0,color='g',alpha = 0.5,width=width,align='edge')

fig.legend()

def init():
    for rect in Graph:
        rect.set_height(0)
    return Graph
def anibar(i):
    for rect, y in zip(Graph,ims[i]):
        rect.set_height(y.magnitude)
    label.set_text('t = ' + str(round(ts[i],1)) + ' s')
    return Graph
ani = anim.FuncAnimation(fig,anibar,init_func=init,frames=len(ims),blit=True,interval=500,repeat_delay=500)
# ani.to_html5_video()

# html.HTML(ani.to_html5_video())
# ani.save('PSD_265.html',writer='imagemagick')
ani.save('PSD_265.gif',writer='imagemagick' ,extra_args="convert")

# ani.save('PSD_265_v.png',dpi=50,writer=AnimatedPNGWriter(fps=1))

np.savetxt('PSD_265.csv', (BinD.magnitude,BinC.magnitude), delimiter=',')

```
