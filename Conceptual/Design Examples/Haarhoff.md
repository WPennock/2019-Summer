```python
from aguaclara.play import *

Q = 5*u.L/u.s
Temp = 0*u.degC
nu = pc.viscosity_kinematic_water(Temp.to(u.K))
h_L = 40*u.cm
Gtheta = 3.7e4
G = (u.g_0*h_L/(nu*Gtheta)).to(1/u.s)
G
t = Gtheta/G
# K = 3.2
K = 2.56
p = 1
w = 0.1*u.m
r = 10

# t_turn = np.arange(1,40,1)*u.s
# N = t/t_turn
# N
N = np.arange(2,300,1)

B = (((N-1)*Q**2*K/(2*nu*G**2*t*r**2))**(1/4)).to(u.m)
q = (Q*t-(N-1)*r*B**2*p*w)/(N*r*B**3)-2*p
q
q[172] #119
Index = 172
N[Index]
B[Index]
Width = (q[Index]+2*p)*B[Index]
Width
Depth = r*B[Index]
Depth
Length = N[Index]*B[Index]
Length
Length_Actual = Length + (N[Index]-1)*w
Length_Actual  
Area = Length*Width
Area
Area_Actual = Width*Length_Actual
Area_Actual
Volume = Length*Width*Depth
Volume
(1800*u.ft**3).to(u.m**3)

P = (100*u.L/u.s*40*u.cm*u.g_0*1000*u.kg/u.m**3).to(u.kW)
(P/0.6).to(u.hp)
Energy = (P*1*u.year).to(u.kWh)
Energy/0.6
Energy/0.6*0.10*u.USD

2*0.45*u.m*1.7*u.m + 0.1*u.m*1.7*u.m
# t_turn[4]
G_Exam = np.array([40, 70, 35, 20, 50, 35, 25, 50, 50, 50])
t_Exam = np.array([1800, 500, 500, 500, 420, 420, 420, 1000, 1000, 1000])
Gt_Exam = G_Exam*t_Exam
Gt_Exam
np.mean(Gt_Exam)
np.min(Gt_Exam)
np.max(Gt_Exam)
```
