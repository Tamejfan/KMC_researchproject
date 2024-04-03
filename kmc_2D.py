import numpy as np
import matplotlib.pyplot as plt
from random import uniform, expovariate

# D = k_tot * d^2/number of k's,
#defining transitions for sytstem
transitions1 = [(3,10),(2,10),(4,10),(7,10)]
transitions2 = [(1,10),(3,10),(5,10),(8,10)]
transitions3 = [(2,10),(1,10),(6,10),(9,10)]
transitions4 = [(6,10),(5,10),(7,10),(1,10)]
transitions5 = [(4,10),(6,10),(8,10),(2,10)]
transitions6 = [(5,10),(4,10),(9,10),(3,10)]
transitions7 = [(9,10),(8,10),(1,10),(4,10)]
transitions8 = [(7,10),(9,10),(2,10),(5,10)]
transitions9 = [(8,10),(7,10),(3,10),(6,10)]
#k_total for each site need to be created (should be improved) 
k_totals = []
jumppossibilies = len(transitions1) #only works for sytems with same amount of transtions on each site
sum_k = 0
for k in transitions1:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions2:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions3:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions4:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions5:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions6:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions7:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions8:
    sum_k += k[1]
k_totals.append(sum_k)
sum_k = 0
for k in transitions9:
    sum_k += k[1]
k_totals.append(sum_k)
#change box size for the system
boxlength = 3
#define postions for sites
site1 = (np.array([0.0,0.0,0.0]), transitions1, k_totals[0])
site2 = (np.array([1.0,0.0,0.0]), transitions2, k_totals[1])
site3 = (np.array([2.0,0.0,0.0]), transitions3, k_totals[2])
site4 = (np.array([0.0,1.0,0.0]), transitions4, k_totals[3])
site5 = (np.array([1.0,1.0,0.0]), transitions5, k_totals[4])
site6 = (np.array([2.0,1.0,0.0]), transitions6, k_totals[5])
site7 = (np.array([0.0,2.0,0.0]), transitions7, k_totals[6])
site8 = (np.array([1.0,2.0,0.0]), transitions8, k_totals[7])
site9 = (np.array([2.0,2.0,0.0]), transitions9, k_totals[8])
sites = [site1, site2, site3, site4, site5, site6, site7, site8, site9]
nsteps = 100
nsim = 1000
curr_site = 1
def run_sim(trajectories: list[list[int]],timess: list[list[int]], sumtimes: list[list[int]], distances: list[list[int]]):
    curr_site = 1
    #creation of trajectory and times
    trajectory = [curr_site] 
    times = [0]
    t = 0
    sumtime =[0]
    for i in range(nsteps):
        k_tot = sites[curr_site - 1][2]
        x = uniform(0, k_tot)
        t_step = expovariate(k_tot)
        k_trans = 0
        for transitions in sites[curr_site - 1][1]:
            k_trans += transitions[1]
            if x <= k_trans:
                curr_site = transitions[0]
                break
        trajectory.append (curr_site)
        t += t_step
        times.append(t_step)
        sumtime.append(t)
    trajectories.append(trajectory)
    timess.append(times)
    sumtimes.append(sumtime)
    movement = 0
    distance=[0]
    #calculation of movement from the trajectory
    for steps in range(nsteps):
        curr_site = sites[trajectory[steps] -1][0]
        next_site = sites[trajectory[steps + 1] -1][0]
        length = np.linalg.norm(next_site - curr_site)
        if length > 1:#!!!
            length = -1 * (boxlength - np.linalg.norm(next_site - curr_site))
        elif length < -1:#!!!
            length = 1 * (boxlength - np.linalg.norm(next_site - curr_site))
        if trajectory[steps] < trajectory[steps + 1]:
            movement += length
        elif trajectory[steps] > trajectory[steps + 1]:
            movement -= length
        distance.append(movement)
    distances.append(distance)
trajectories = []
timess =[]
sumtimes = []
distances = []
for i in range (nsim):
     run_sim( trajectories, timess, sumtimes, distances)


MSDs = []
#claclulation of MSD
MSD = np.zeros(len(distances[0]))
for i in range(nsim):
    MSDs.append(np.array(distances[i])**2)
    MSD+=MSDs[i]
MSD/=nsim

#averaging of times
avgtimes = np.zeros(len(sumtimes[0]))
for i in range(nsim):
    avgtimes+=sumtimes[i]
avgtimes/=nsim
A = np.vstack([avgtimes, np.ones(len(avgtimes))]).T
k, m = np.linalg.lstsq(A, MSD,rcond=None)[0]

#print(trajectories)
#print(k_totals)
#print(t)
#print(sumtimes)
#print(distances)
#print(timess)
#print(MSD)
#print(k)
#print(avgtimes)
_ = plt.plot(avgtimes, MSD, 'o', label='Data', markersize=10)
_ = plt.plot(avgtimes, k*np.array(avgtimes) + m, 'r', label='linear Fit')
_ = plt.legend()
plt.show()
D = k/jumppossibilies
print(D)