'''

create_reference_list.py

'''

import os 
import numpy as np
import matplotlib.pyplot as plt
from hst_fun import *
from linear6d import *

######################

dir = os.path.dirname(os.path.realpath(__file__))

datstr = "_HOROLOGIUM-I_F814W_oc"

flags, ra, dec, xr, yr, flux, c_star, magr, id, xc, yc, xo, yo = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12

brightcut = 1000
pixtol = 5

filelist = ["jdan21l8q", "jdan21laq", "jdan21lhq", "jdan21llq"]

dither0 = np.loadtxt(dir+"/Data/"+filelist[0]+datstr+".dat")
dither1 = np.loadtxt(dir+"/Data/"+filelist[1]+datstr+".dat")
dither2 = np.loadtxt(dir+"/Data/"+filelist[2]+datstr+".dat")
dither3 = np.loadtxt(dir+"/Data/"+filelist[3]+datstr+".dat")


dither0 = dither0[dither0[:,magr].argsort()]
dither1 = dither1[dither1[:,magr].argsort()]
dither2 = dither2[dither2[:,magr].argsort()]
dither3 = dither3[dither3[:,magr].argsort()]

print(dither0[:10, magr])

bright_d0 = dither0[:brightcut]
bright_d1 = dither1[:brightcut]
bright_d2 = dither2[:brightcut]
bright_d3 = dither3[:brightcut]


np.savetxt(dir+"/Data/"+filelist[0]+datstr+"_bright.reg", bright_d0[:,[xr, yr]], fmt="%1.6f")
np.savetxt(dir+"/Data/"+filelist[1]+datstr+"_bright.reg", bright_d1[:,[xr, yr]], fmt="%1.6f")
np.savetxt(dir+"/Data/"+filelist[2]+datstr+"_bright.reg", bright_d2[:,[xr, yr]], fmt="%1.6f")
np.savetxt(dir+"/Data/"+filelist[3]+datstr+"_bright.reg", bright_d3[:,[xr, yr]], fmt="%1.6f")


plt.clf()
plt.scatter(bright_d1[:,xo], bright_d1[:,yo], s=5, color="tomato")
plt.scatter(bright_d2[:,xo], bright_d2[:,yo], s=4, color="forestgreen")
plt.scatter(bright_d3[:,xo], bright_d3[:,yo], s=3, color="royalblue")
plt.scatter(bright_d0[:,xo], bright_d0[:,yo], s=2, color="black")

plt.savefig(dir+"/Plots/offset_bright.png", dpi=300)


####
   # Start first loop to establish masterlist
####

masterid, b1_matchid = matchlistid(bright_d0, bright_d1, pixtol, xo, yo, xo, yo)

master = bright_d0[masterid]

print(len(master))

masterid, b2_matchid = matchlistid(master, bright_d2, pixtol, xo, yo, xo, yo)

master = master[masterid]

print(len(master))

masterid, b3_matchid = matchlistid(master, bright_d3, pixtol, xo, yo, xo, yo)

master = master[masterid]

print(len(master))

####
   # Now match the final master list to each of the other three dithers and 
   # trim down the dither lists
####

b3_match = bright_d3[b3_matchid]

tempid, b1_matchid = matchlistid(master, bright_d1, pixtol, xo, yo, xo, yo)

b1_match = bright_d1[b1_matchid]

tempid, b2_matchid = matchlistid(master, bright_d2, pixtol, xo, yo, xo, yo)

b2_match = bright_d2[b2_matchid]


####
   # Create weights and perform the 6D linear transformations
####

weights = np.zeros((len(master)))
weights.fill(1.0)


new_b1, new_d1 = test_linear(b1_match[:,xo], b1_match[:,yo], master[:,xo], master[:,yo], weights, weights, dither1[:,xo], dither1[:,yo])

new_b2, new_d2 = test_linear(b2_match[:,xo], b2_match[:,yo], master[:,xo], master[:,yo], weights, weights, dither2[:,xo], dither2[:,yo])

new_b3, new_d3 = test_linear(b3_match[:,xo], b3_match[:,yo], master[:,xo], master[:,yo], weights, weights, dither3[:,xo], dither3[:,yo])


####
   # Check on the differences between the original lists and the transformed positions
####

rn = 6 #rounding number

# Dither 1
print("")
print("Dither 1 stats: ")
print("Average offset (x, y): " + str(np.round(np.average(b1_match[:,xo] - master[:,xo]), rn))+", "+str(np.round(np.average(b1_match[:,yo]-master[:,yo]), rn)))

print("Sigma of offset (x, y): " + str(np.round(np.std(b1_match[:,xo] - master[:,xo]), rn)) + ", " + str(np.round(np.std(b1_match[:,yo]-master[:,yo]), rn)))

print("Average offset (x, y): " + str(np.round(np.average(new_b1[:,0] - master[:,xo]), rn)) + ", " + str(np.round(np.average(new_b1[:,1] - master[:,yo]), rn)))

print("Sigma of offset (x, y): " + str(np.round(np.std(new_b1[:,0] - master[:,xo]), rn)) + ", " + str(np.round(np.std(new_b1[:,1] - master[:,yo]), rn)))


# Dither 2

print("")
print("Dither 2 stats: ")
print("Average offset (x, y): " + str(np.round(np.average(b2_match[:,xo] - master[:,xo]), rn))+", "+str(np.round(np.average(b2_match[:,yo]-master[:,yo]), rn)))

print("Sigma of offset (x, y): " + str(np.round(np.std(b2_match[:,xo] - master[:,xo]), rn)) + ", " + str(np.round(np.std(b2_match[:,yo]-master[:,yo]), rn)))

print("Average offset (x, y): " + str(np.round(np.average(new_b2[:,0] - master[:,xo]), rn)) + ", " + str(np.round(np.average(new_b2[:,1] - master[:,yo]), rn)))

print("Sigma of offset (x, y): " + str(np.round(np.std(new_b2[:,0] - master[:,xo]), rn)) + ", " + str(np.round(np.std(new_b2[:,1] - master[:,yo]), rn)))

# Dither 3

print("")
print("Dither 3 stats: ")
print("Average offset (x, y): " + str(np.round(np.average(b3_match[:,xo] - master[:,xo]), rn))+", "+str(np.round(np.average(b3_match[:,yo]-master[:,yo]), rn)))

print("Sigma of offset (x, y): " + str(np.round(np.std(b3_match[:,xo] - master[:,xo]), rn)) + ", " + str(np.round(np.std(b3_match[:,yo]-master[:,yo]), rn)))

print("Average offset (x, y): " + str(np.round(np.average(new_b3[:,0] - master[:,xo]), rn)) + ", " + str(np.round(np.average(new_b3[:,1] - master[:,yo]), rn)))

print("Sigma of offset (x, y): " + str(np.round(np.std(new_b3[:,0] - master[:,xo]), rn)) + ", " + str(np.round(np.std(new_b3[:,1] - master[:,yo]), rn)))


####
   # Create new plot to look at new stellar positions
####

plt.clf()
plt.scatter(new_b1[:,0], new_b1[:,1], s=5, color="tomato")
plt.scatter(new_b2[:,0], new_b2[:,1], s=4, color="forestgreen")
plt.scatter(new_b3[:,0], new_b3[:,1], s=3, color="royalblue")
plt.scatter(master[:,xo], master[:,yo], s=2, color="black")

plt.savefig(dir+"/Plots/offset_bright_transformed.png", dpi=300)


## Check the all position output

plt.clf()
plt.scatter(dither0[:,xo], dither0[:,yo], s=5, color="black")
plt.scatter(new_d1[:,0], new_d1[:,1], s=2, color="tomato")
plt.savefig(dir+"/Plots/dither0_dither1_transformed.png", dpi=300)


## Compare to the original positions

plt.clf()
plt.scatter(dither0[:,xo], dither0[:,yo], s=5, color="black")
plt.scatter(dither1[:,xo], dither1[:,yo], s=2, color="tomato")
plt.savefig(dir+"/Plots/dither0_dither1.png", dpi=300)


####
   # Now test to see how many matches the original dither 0 and 1 versus now
   # WARNING: This will take ~10x as long as the rest of the script
####

'''

all0_matchid, all1_matchid = matchlistid(dither0, dither1, 2.0, xo, yo, xo, yo)

new0_matchid, new1_matchid = matchlistid(dither0, new_d1, 2.0, xo, yo, 0, 1)

# Print results

print("")
print("Number of matches when matching all dither 0 to all dither 1: ")
print(len(all0_matchid))
print("Number of matches when mattching all dither 0 to transformed dither 1: ")
print(len(new0_matchid))

'''
