#Design Example for 100 L/s Vertically-Baffled Hydraulic Flocculator
```python
from aguaclara.play import *
```

## Set Initial Conditions
```python
Q = 100*u.L/u.s
h_L = 40*u.cm
H_Min = 2*u.m
T_Des = 0*u.degC
T_Des = T_Des.to(u.degK)
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

### Channel Length
```python
W_Human = 45*u.cm
L_V_Max = (V_Floc/(2*W_Human*H_Min)).to(u.m)
L_V_Max
L_Sed_Max = 7*u.m
L = np.min(np.array([L_V_Max.magnitude,L_Sed_Max.magnitude]))*u.m
L

# After computing W

# L_V_Max_New = (V_Floc/(2*W_Min_0*H_Min)).to(u.m)
# L_V_Max_New
```

### Channel Width
```python
# K = 2.56
# def W_MIN_HYD(Q,H,K,nu,G):
#     return ((Q*3/H)*(K/(2*H*nu*G**2))**(1/3)).to(u.m)
# W_Min_Hyd = W_MIN_HYD(Q,H_Min,K,nu,G)
# W_Min_Hyd
# W_Min = np.max([(W_Human.to(u.m)).magnitude,W_Min_Hyd.magnitude])*u.m
# W_Min
#
# W_Floc = V_Floc/(H_Min*L)
# W_Floc
#
# n_Channel = np.floor(W_Floc/W_Min)
# n_Channel
# n_Channel = 6 # Need even number
# W = W_Floc/n_Channel
# W
```
### Channel Depth
```python
h_Free = 10*u.cm
D_Min = H_Min + h_L + h_Free
D_Min
```

## Baffle Module Dimensions
```python
# def S_DES(n_Obs,n_Channel,L,K,Q,h_L,W):
#     return (np.cbrt((((n_Obs+1)*n_Channel*L*K*Q**2)/(2*u.g_0*h_L*W**2)).to(u.m**3)))*u.m
#
# S_Des = S_DES(0,n_Channel,L,K,Q,h_L,W)
# S_Des
# v_Scour = 15*u.cm/u.s
# def S_MAX(H,Q,v_Scour,W):
#     return np.min([((H/3).to(u.m)).magnitude,((Q/(v_Scour*W)).to(u.m)).magnitude])*u.m
#
# S_Max = S_MAX(H_Min,Q,v_Scour,W_0)
# S_Max
# S_Des<S_Max
#
# Pi_Des = H_Min/S_Des
# Pi_Des
# 3<Pi_Des<7
#
# # Need zero obstacles.
# S = S_Des
# H_e = Pi_Des*S
# H_e
# # Check minor loss
# n_Obs = 0
# h_L_Act = (n_Obs+1)*n_Channel*(L/S)*K*Q**2/(2*u.g_0*W**2*S**2)
# h_L_Act.to(u.cm)# Check theta
# theta_Act = (n_Channel*L*W*(H_Min+h_L/2))/Q
# theta_Act.to(u.s)
# # Check major loss
# n_space = np.ceil(L/S_Des)
# n_space
# Length = n_space*H_Min
# rough = (S_Des*0.3*u.mm+W*0.002*u.mm)/(S_Des+W)
# rough
# bonus = pc.headloss_fric_rect(Q,W,S_Des,Length,nu,rough,False)
# bonus.to(u.m)
```
# Original Design Method
```python
# Channel Dimensions
W_Floc = V_Floc/(H_Min*L)
W_Floc
K = 2.56

def WMINHYD_0(Q,H,K,nu,G):
    return ((3*Q/H)*(K/(2*H*nu*G**2))**(1/3)).to(u.m)

W_Min_Hyd_0 = WMINHYD_0(Q,H_Min,K,nu,G)
W_Min_Hyd_0
W_Min_0 = np.max(np.array([(W_Human.to(u.m)).magnitude,(W_Min_Hyd_0.to(u.m)).magnitude]))*u.m
W_Min_0
n_Channel_0 = np.floor(W_Floc/(2*W_Min_0))*2
n_Channel_0
W_0 = W_Floc/n_Channel_0
W_0
# Baffle Dimensions
## He
def HEMAX(K,nu,G,Q,W):
  return ((K/(2*nu*G**2)*(Q*6/W)**3)**(1/4)).to(u.m)

HeMax = HEMAX(K,nu,G,Q,W_0)  
HeMax
nExp = np.ceil(H_Min/HeMax)
nExp
He = H_Min/nExp
He
nObs = H_Min/He-1
nObs
## S
def S(K,He,G,nu,Q,W):
  return ((K/(2*He*G**2*nu))**(1/3)*(Q/W)).to(u.m)
S_0 = S(K,He,G,nu,Q,W_0)
S_0
S_0<S_Max

Pi_0 = He/S_0
Pi_0
3<Pi_0<6
(Q/(W_0*S_0)).to(u.m/u.s)
S_0
# Baffle Heights
H_Bot = H_Min-S_0
H_Bot
H_Free = 0.1*u.m
H_Top = H_Min-S_0+h_L+0.5*H_Free
H_Top
# Check minor loss
h_L_Act_0 = (nObs+1)*n_Channel_0*(L/S_0-1)*K*Q**2/(2*u.g_0*W_0**2*S_0**2)
h_L_Act_0.to(u.cm)
```
