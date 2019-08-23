# Attempt at learning generalized doubling rules

## Discrete, same time for collision
```python
from aguaclara.play import *
import aguaclara.research.floc_model as floc
import pdb
pop = np.ones(1000000)
frac = 0.5
collisions = 11
for i in range(1,collisions+1):
    for j in range(0,i):
        # pdb.set_trace()
        sizej = int(2**j)
        popj = np.where(pop==sizej)
        coll = int(np.floor(frac*0.5*np.size(popj)/2)*2)
        if coll > 1:
            pop[popj[0][0:coll]] = 2*sizej
            pop[popj[0][-coll:]] = 0
pop = pop[pop!=0]
bins = np.arange(0,np.log2(np.max(pop)))

plt.hist(pop,bins=bins,align='mid')
bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
plt.xticks(np.arange(min(bins)+bin_w/2, max(bins), bin_w), bins)
plt.xlim(bins[0], bins[-1])
plt.show()


```
## Continuous, different collision times
Take time step, for each bin, calculate collision time and divide time step by
collision time. Multiply this fraction by current population in bin and add half
this to next bin.

```python

# Figuring out bins
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


alpha = 2*Gamma-Gamma**2
frac_vol_floc_initial(Dose,C_0,floc.PACl,floc.Clay)
def t_c(Dose,conc,coag,material,EDR,temp):
    return ((np.pi*6/frac_vol_floc_initial(Dose,C_0,coag,material))**(2/3)/np.pi*(pc.viscosity_kinematic(temp)/EDR)**(1/2)).to(u.s)
t_c(Dose,C_0,floc.PACl,floc.Clay,25*u.mW/u.kg,22*u.degC)

C_0 = 90*u.NTU
BinC = np.zeros(len(BinD))*u.NTU
BinC[0] = C_0
Dose = 2.65*u.mg/u.L
t = 0*u.s
T = 387*u.s
dt = 0.1*u.s
(T/dt).to(u.dimensionless)
while t<=T:
    for i in range(0,len(BinD)-1):
        tc = t_c(Dose,BinC[i],floc.PACl,floc.Clay,24*u.mW/u.kg,22*u.degC)
        dC = alpha*dt/tc*BinC[i]
        BinC[i] = BinC[i] - dC
        BinC[i+1] = BinC[i+1] + dC
    t = t + dt
BinC
plt.bar(BinD,BinC)
plt.savefig('PSD_265.png')
plt.bar(BinD,BinC/BinC[0])
np.savetxt('PSD_265.csv', (BinD.magnitude,BinC.magnitude), delimiter=',')
# B = np.array([95,100,105,80,300])
# np.sqrt(np.mean(B))
# np.mean(np.sqrt(B))


# (np.sqrt(8*0.06/.025*u.g_0*(.05)*1*u.mm)).to(u.cm/u.s)
vs = 0.5*u.ft/u.min
vs.to(u.mm/u.s)
```
