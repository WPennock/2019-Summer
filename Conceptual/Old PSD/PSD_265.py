from aguaclara.play import *
import aguaclara.research.floc_model as floc
import pdb
import pickle as pkl
import matplotlib.animation as anim
1*u.NTU
PSD_265 = pkl.load(open("PSD_265.pkl","rb"))

# Create Animation function
plt.clf(),plt.close('all')
fig, ax = plt.subplots()
