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
from scipy.interpolate import interp1d

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

## Load Data
```python
PSD = pkl.load(open('PSD.pkl','rb'))*90*u.NTU

CDF = np.zeros(PSD.shape)

for i in range(0,len(PSD)):
    for j in range(0,len(PSD[i])):
        CDF[i][j] = np.sum(PSD[i][0:j+1].magnitude)

CDF = CDF*u.NTU        


Bins = pkl.load(open('Bins.pkl','rb'))

Temps = pkl.load(open('Temps.pkl','rb'))*u.degC
Tempsp = []
for i in range(0,len(Temps)):
    Tempsp.append(np.append(Temps[i]*u.degC,Temps[i][-1]*u.degC))

Doses = pkl.load(open('Doses.pkl','rb'))*u.mg/u.L
Dosesp = []
for i in range(0,len(Doses)):
    Dosesp.append(np.append(Doses[i]*u.mg/u.L,Doses[i][-1]*u.mg/u.L))

BinsDens = np.zeros([12,7])

for i in range(0,len(BinsDens)):
    for j in range(0,len(BinsDens[i])):
        BinsDens[i][j] = dens_floc(Dosesp[i][j]*u.mg/u.L,C_0,Dim,Bins[i][j]*u.um,floc.PACl,floc.Clay,Tempsp[i][j]*u.degC).magnitude

nPP = nP(CDF,BinsDens*u.kg/u.m**3,Bins*u.um)
nPC = np.zeros(nPP.shape)

for i in range(0,len(nPC)):
    for j in range(0,len(nPC[i])):
        nPC[i][j] = np.sum(nPP[i][0:j+1].magnitude)

nPC[0]/nP_0    
```

## Define Functions
```python
def frac_vol_floc_initial(ConcAluminum, ConcClay, coag, material):
    return ((floc.conc_precipitate(ConcAluminum, coag)/coag.PrecipDensity) + (ConcClay / material.Density)).to(u.dimensionless)


def conc_floc(ConcAluminum, concClay, coag):
    return floc.conc_precipitate(ConcAluminum, coag) + concClay

def dens_floc_init(ConcAluminum, ConcClay, coag, material):
    return (conc_floc(ConcAluminum, ConcClay, coag)
            / frac_vol_floc_initial(ConcAluminum, ConcClay, coag, material)
            )

def dens_floc(ConcAl, ConcClay, DIM_FRACTAL, DiamTarget, coag, material, Temp):
    WaterDensity = pc.density_water(Temp)
    return ((dens_floc_init(ConcAl, ConcClay, coag, material)
             - WaterDensity
             )
            * (material.Diameter / DiamTarget)**(3 - DIM_FRACTAL)
            + WaterDensity
            ).to(u.kg/u.m**3)

def frac_vol_clay(ConcClay, Density):
    return ((ConcClay / Density)).to(u.dimensionless)

def gamma_coag(ConcClay, ConcAluminum, coag, material,
               DiamTube, RatioHeightDiameter):
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

# def time_col_laminar(EnergyDis, Temp, ConcAl, ConcClay, coag, material,
                     # DiamTarget, DiamTube, DIM_FRACTAL, RatioHeightDiameter):
    # return (((1/6) * ((6/np.pi)**(1/3))
             # * frac_vol_floc_initial(ConcAl, ConcClay, coag, material) ** (-2/3)
             # * (pc.viscosity_kinematic(Temp) / EnergyDis) ** (1 / 2)
             # * (DiamTarget / material.Diameter) ** (2*DIM_FRACTAL/3 - 2)
             # )  # End of the numerator
            # / (gamma_coag(ConcClay, ConcAl, coag, material, DiamTube,
                          # RatioHeightDiameter)
               # )  # End of the denominator
            # ).to(u.s)


def dnP(k,nP,Diam,EDR,temp,alpha,dt):
    return (k*np.pi*alpha*nP**(5/3)*Diam**2*(EDR/pc.viscosity_kinematic(temp))**(1/2)*dt).to(1/u.m**3)

def nP(conc,Density,Diameter):
    return(6*conc.to(u.kg/u.m**3)/(np.pi*Diameter.to(u.m)**3*Density.to(u.kg/u.m**3))).to(u.m**-3)

# def t_c(conc,Density,EDR,temp):
    # return ((np.pi/(6*frac_vol_clay(conc,Density)))**(2/3)/np.pi*(pc.viscosity_kinematic(temp)/EDR)**(1/2)).to(u.s)
#
# while t<=T:
    # dC = np.zeros(len(BinD)-1)
    # for i in range(0,len(BinD)-1):
        # tc = t_c(BinC[i],DensD[i],EDR,Temp).magnitude
        # dC[i] = kf*alpha*dt/tc*(BinC[i].to(u.NTU)).magnitude
    # for i in range(0,len(BinD)-1):
        # BinC[i] = BinC[i] - dC[i]*u.NTU
        # BinC[i+1] = BinC[i+1] + dC[i]*u.NTU
    # if round(t,1)%10 == 0:
    # ims.append(BinC.magnitude.copy())    
    # ts.append(t)
    # t = t + dt
def IntegratednP(t,material,k,EDR,temp,alpha,nP0):
    return (((2/3*((EDR/pc.viscosity_kinematic(temp))**(1/2)).to(1/u.s)*t*k*alpha*np.pi*material.Diameter**2)+nP0**(-2/3))**(-3/2)).to(1/u.m**3)
IntegratednP(T,floc.Clay,0.1,EDR,Temp,alpha,)
```

## Simulation
```python
# Set up simulation

# Conditions
ind = 4
Dose = 2.04*u.mg/u.L
EDR = 22.8*u.mW/u.kg
Temp = 23.2*u.degC
C_0 = 90*u.NTU
nP_0 = floc.particle_number_concentration(C_0,floc.Clay)

Gamma = gamma_coag(C_0,Dose,floc.PACl,floc.Clay,1.25*u.inch,floc.RATIO_HEIGHT_DIAM)
Gamma
alpha = 2*Gamma-Gamma**2
alpha

# Set up simulation bins
d_Max = 110*u.um
d_0 = 7*u.um
Dim = floc.DIM_FRACTAL
NumP = (d_Max/d_0)**Dim
NBin = np.ceil(np.log2(NumP))
NBin
Coll = np.arange(0,NBin+1)
Coll
BinD = floc.diam_fractal(Dim,d_0,Coll)
BinD
BinnP = np.zeros(len(BinD))/u.m**3
BinnP[0] = nP_0
BinnP_0 = BinnP
BinD
DensD = dens_floc(Dose,C_0,Dim,BinD,floc.PACl,floc.Clay,Temp)
DensD

# Set up time steps
t = 0
T = 397
dt = 0.1
(T/dt)
# kf = 1
# kfit = pkl.load(open("k_041.pkl",'rb'))
# kf = 5*np.exp(-BinD[:-1]/BinD[0])
# kf = 2**(-Coll[:-1]*(1-1/Dim)) # Theory
kf = 10*(BinD[:-1]/BinD[0])**(Dim-1)
kf

# Clear previous results
ims = []
ts = []

while t<=T:
    d_nP = dnP(kf,BinnP[:-1],BinD[:-1],EDR,Temp,alpha,dt*u.s)
    for i in range(0,len(BinD)-1):
        BinnP[i] = BinnP[i] - d_nP[i]
        BinnP[i+1] = BinnP[i+1] + d_nP[i]
    # if round(t,1)%10 == 0:
    ims.append(BinnP.magnitude.copy())    
    ts.append(t)
    t = t + dt

# Quick Quality Check
len(ims)
ims[0]
ims[-1]
```

## Examining behavior of primary particles
```python
prim = np.zeros(len(ims))
for i in range(0, len(prim)):
    prim[i]=ims[i][0]

plt.loglog(ts,prim,label='Simulation')
plt.loglog(np.array(ts[1:]),np.array(ts[1:])**(-3/2),label=r'$t^{-3/2}$')
# plt.loglog(ts,Integrated(ts*u.s,floc.Clay,kf,EDR,Temp,alpha),label='Model')
# plt.plot(np.array(ts[1:]),.9*np.array(ts[1:])**(-3/2),label=r'$9t^{-3/2}$')
plt.xlabel(r'$t$')
plt.ylabel(r'$C_P$')
plt.legend()
# plt.savefig('log1.png')    
plt.show()

```

## Save Results
```python
kf

pkl.dump(ims,open("PSD_204.pkl","wb"))
pkl.dump(ts,open("t_204.pkl","wb"))
```

## Trim Data
```python
# # Load results (optional)
# ims = pkl.load(open("PSD_041.pkl","rb"))
# ts = pkl.load(open("t_041.pkl","rb"))
#
# for i in range(0,len(ims)):
#     ims[i] = ims[i].magnitude
# ims = ims*u.NTU

ims10 = []
t10 = []

for i in range(0,len(ims)):
    if round(ts[i],2)%10 == 0.0:
        ims10.append(ims[i].copy())
        t10.append(ts[i])
    elif i == len(ims)-1:    
        ims10.append(ims[i].copy())
        t10.append(ts[i])
ims10 = ims10*u.m**(-3)

SIM = np.zeros(ims10.shape)

for i in range(0,len(ims10)):
    for j in range(0,len(ims10[i])):
        SIM[i][j] = np.sum(ims10[i][0:j+1].magnitude)
SIM = SIM*u.m**(-3)
```

## Animation
```python
plt.clf(),plt.close('all')
fig, ax = plt.subplots()

# Plot data graphs
nPC1 = nPC[ind]/(nP_0)
nPC2 = nPC[ind+6]/(nP_0)
len(Bins[ind])
len(nPC1)
ax.plot(Bins[ind][:-1],nPC1[:-1],'b+',label=r"1st 2.04 mg/L")
ax.plot(Bins[ind+6][:-1],nPC2[:-1],'r+',label=r"2nd 2.04 mg/L")
ax.set_xlabel(r'Particle Diameter ($\mathrm{\mu m}$)')
ax.set_ylabel(r'$\sum\frac{n_\mathrm{P}}{n_\mathrm{P_0}}$')

# Plot simulation graphs

nPCS = []
for i in range(0,len(SIM)):
    nPCS.append(SIM[i].copy()/(nP_0))
nPCS
splines = []
for i in range(0,len(nPCS)):
    splines.append(interp1d(BinD[:-1], nPCS[i][:-1], kind='cubic'))    

nPC_0 = np.ones(len(BinD))*nP_0
nPCS_0 = nPC_0/(nP_0)

spline_0 = interp1d(BinD,nPCS_0,kind='cubic')

Diams = np.arange(7,BinD[-2].magnitude,0.1)*u.um

# Graph, = ax.plot(Diams,spline_0(Diams),'g-',label='Simulation, k = '+str(round(kf,3)))
Graph, = ax.plot(Diams,spline_0(Diams),'g-',label=r'Simulation, k = 1')
 # \left(\frac{d}{d_0}\right)^{1-D}$

ax.legend(loc=4)
label = ax.text(0.42,0.96,'t = 0 s',transform=ax.transAxes)

ax.set_ylim(1e-4,1e-3)
ax.set_yscale('log')

def init():
    Graph.set_data([],[])
    return Graph,
def aniline(i):
    y = splines[i](Diams)
    Graph.set_data(Diams,y)
    label.set_text('t = ' + str(round(t10[i],0)) + ' s')
    return Graph,
ani = anim.FuncAnimation(fig,aniline,init_func=init,frames=len(CDFS),blit=True,interval=750,repeat_delay=3000,repeat=True)

ani.save('Lambda_204_fractal.gif',writer='imagemagick' ,extra_args="convert") # run 3 times
```

## Miscellaneous Calculations
```python
from scipy.optimize import curve_fit

Dar = np.array([0.41,0.87,1.41,2.04,2.65])
kar = np.array([0.1,0.067,0.049,0.038,0.028])
def expfit(x, a, b):
    return a * np.exp(-b * x)

expopt, expcov = curve_fit(expfit, Dar, kar)

def invfit(x, a):
    return a/x

invopt, invcov = curve_fit(invfit, Dar, kar)

plt.plot(Dar,kar,'+')
plt.plot(Dar,expfit(Dar,expopt[0],expopt[1]),label='exp')
plt.plot(Dar,invfit(Dar,invopt[0]),label='inv')
plt.legend()

Gar = gamma_coag(C_0,Dar*u.mg/u.L,floc.PACl,floc.Clay,1.25*u.inch,floc.RATIO_HEIGHT_DIAM)
plt.plot(Dar,Gar)
aar = 2*Gar - Gar**2
plt.plot(Dar,aar)

```
