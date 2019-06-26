#Design Example for Vertically-Baffled Hydraulic Flocculator
```python
from aguaclara.play import *
from aguaclara.
```
## Set Initial Conditions
```python
Q = 10*u.L/u.s
h_L = 40*u.cm
H = 2*u.m
T_Des = 0*u.degC
Gtheta = 3.7e4
```
## Find General Design Parameters
```python
nu = pc.viscosity_kinematic(T_Des)
G = (u.g_0*h_L/(nu*Gtheta)).to(1/u.s)
G
theta = Gtheta/G
theta
V_Floc = (Q*theta).to(u.m**3)
V_Floc
```
## Flocculator Dimensions
```python
```
