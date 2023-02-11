import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 8})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams["font.sans-serif"] = "Arial"
import csv

#script for generating confinement quench of a G1-sizer model
#code for Fig. 5F
folder = 'confinement_quench'

#set quench colony numbers
initial_area = np.array([5000,12000,47000])/10
radius_labels = ['.5','.75','1.5']

run_num = len(initial_area)

def reaction_rate(area):

    #division rate as a function of area

    if area > 1:
        r = division_s
    else:
        r = 0

    return r

#define timesteps and dt
num_steps = 800
dt = .005
unit_time = int(1/dt)

#choose number of replicates (varies with colony size)
replicates = np.array([20,10,5])

#we set the division rate to obtain a desired average rate total cell cycle length
#apprxoimately consistent with experiment
division_s = 3

#initializing storage variables for parameter sweep
pop_num_all = []
avg_area_all = []
var_area_all = []
Sfrac_all = []
confinement_all = []

#perform parameter sweep
for run in np.arange(run_num):

    #we keep track of population number,
    #total population biomass size, average area
    #variance in area, confinement, and fraction of cells in S/G2/M

    pop_num = np.zeros((replicates[run],num_steps))
    area_tot = np.zeros((replicates[run],num_steps))
    avg_area = np.zeros((replicates[run],num_steps))
    var_area = np.zeros((replicates[run],num_steps))
    confinement = np.zeros((replicates[run],num_steps))
    Sfrac = np.zeros((replicates[run],num_steps))

    for j in np.arange(replicates[run]):

        record_counter = 0

        print(j)

        #initialize colony size to 4 cells
        #with total area summing to 8, and an estimate of the correct S fraction
        N_0 = 4
        a_0 = 2.
        A_0 = N_0*a_0
        a = np.random.uniform(.5,1.,size=N_0)
        a = a/np.sum(a)*a_0*N_0

        #the phase variable represents cell cycle phase:
        #phase < 0 indicates G1, phase > 0 indicates S/G2, division occurs when phase = 0
        phase = np.random.randint(1,unit_time,size=N_0)*np.random.binomial(1,.3*1/(1+1/division_s),size=N_0) - 1
        since_div = phase

        #initialize unconfined growth till the required quench colony size is reached
        while A_0 < initial_area[run]:

            #each cell grows with random addition of area
            #cell phase counter decreases with 1 each time step
            cell_num = len(a)
            a += np.random.uniform(.9,1.1,size=cell_num)*dt
            phase += -1
            since_div += -1

            for cell in np.arange(cell_num):

                #we generate a random number to see if cells enter S
                r = reaction_rate(a[cell])*dt
                if r > np.random.uniform():
                    
                    #a cell is allowed to enter S if it is in G1 phase (i.e. phase < 0)
                    #if a cell enters S its phase is set to a counter to enforce deterministic S length
                    #this amount of time is 1 unit of simulation time, which sets the timescale
                    #for comparison to experiment
                    if phase[cell] < 0:
                        phase[cell] = unit_time


                #when a cell's phase reaches 0, it divides
                if phase[cell] == 0:

                    #cell area is halved
                    a[cell] = a[cell]/2.

                    #a new cell is created
                    a = np.append(a,a[cell])
                    phase = np.append(phase,phase[cell])
                    since_div[cell] = 0
                    since_div = np.append(since_div,since_div[cell])
            A_0 = np.sum(a)

        print(A_0)
        print(A_0/len(a))
        R = np.sqrt(A_0*1/np.pi)

        #once quench colony size is reached, begin simulation under confinement
        for i in np.arange(num_steps):

            area_tot[j,i] = np.sum(a)
            avg_area[j,i] = np.sum(a)/len(a)
            var_area[j,i] = np.std(a)
            Sfrac[j,i] = len(np.where(phase > 0)[0])/len(phase)

            cell_num = len(a)
            R_prev = R

            #colony radius expands outwards at fixed rate
            #this provides a limited amount of area for the growing cells in the colony
            R += 0.*dt
            A  = np.pi*R**2
            A_prev = np.pi*R_prev**2
            confinement[j,i] = 1 - (A-A_prev)/dt/A
            mean_growth = (A-A_prev)/cell_num

            #new area is allocated uniformly across all existing cells with random fluctuation
            growth_disorder = np.random.uniform(-1,1,size=len(a))*mean_growth
            growth_disorder -= np.mean(growth_disorder) 
            a += mean_growth + 1.*growth_disorder

            phase += -1
            since_div += -1


            for cell in np.arange(cell_num):

                #we generate a random number to see if cells enter S
                r = reaction_rate(a[cell])*dt
                if r > np.random.uniform():
                    if phase[cell] < 0:

                        #a cell is allowed to enter S if it is in G1 phase (i.e. phase < 0)
                        #if a cell enters S its phase is set to a counter to enforce deterministic S length
                        #this amount of time is 1 unit of simulation time, which sets a timescale
                        #for comparison to experiment
                        phase[cell] = unit_time


                #when a cell's phase variable reaches 0, it divides
                if phase[cell] == 0:

                    #cell area is halved
                    a[cell] = a[cell]/2.

                    #a new cell is created
                    a = np.append(a,a[cell])
                    phase = np.append(phase,phase[cell])
                    since_div[cell] = 0
                    since_div = np.append(since_div,since_div[cell])

            pop_num[j,i] = len(a)

    pop_num_all.append(pop_num)
    avg_area_all.append(avg_area)
    var_area_all.append(var_area)
    Sfrac_all.append(Sfrac)
    confinement_all.append(confinement)

radius_labels = ['.5 mm','.75 mm','1.5 mm']
colors = ['black','pink','goldenrod']

#when we plot, we scale simulations to approximate experimental values
#1 unit sim time = 9.2 hrs
#1 unit sim volume = 1200 um3

#plot population number as function of time
plt.figure()
for run in np.arange(run_num):
    plt.plot(9.2*dt*np.arange(num_steps),pop_num_all[run].T,linewidth=3,alpha=.1,color=colors[run])
    plt.plot(9.2*dt*np.arange(num_steps),np.mean(pop_num_all[run].T,axis=1),
             linewidth=3,color=colors[run],label=radius_labels[run])
plt.xlabel('time (hrs)',fontsize=18)
plt.xticks([0,10,20,30],[0,10,20,30],fontsize=18)
plt.ylabel('total cells',fontsize=18)
plt.yticks([0,2000,4000,6000],fontsize=18)
plt.legend(fontsize=14)
plt.savefig(folder+'_popnum.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_popnum.pdf',format = 'pdf',bbox_inches="tight")
plt.show()

#plot volume as function of time
plt.figure()
for run in np.arange(run_num):
    plt.plot(9.2*dt*np.arange(num_steps),150*8*avg_area_all[run].T,linewidth=3,alpha=.1,color=colors[run])
    plt.plot(9.2*dt*np.arange(num_steps),150*8*np.mean(avg_area_all[run].T,axis=1),
            linewidth=3,color=colors[run],label=radius_labels[run])
plt.xlabel('time (hrs)',fontsize=18)
plt.xticks([0,10,20,30],[0,10,20,30],fontsize=18)
plt.ylabel('mean volume (in um3)',fontsize=18)
plt.yticks([0,1200,2400],[0,1200,2400],fontsize=18)
plt.legend(fontsize=14)
plt.savefig(folder+'_avgarea.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_avgarea.pdf',format = 'pdf',bbox_inches="tight")
plt.show()

#plot average division fraction as function of time
plt.figure()
trim = 1
for run in np.arange(run_num):
    plt.plot(9.2*dt*(trim+np.arange(num_steps-2*trim)),
             Sfrac_all[run].T[trim:-trim],linewidth=3,alpha=.1,color=colors[run])
    plt.plot(9.2*dt*(trim+np.arange(num_steps-2*trim)),
             np.mean(Sfrac_all[run].T,axis=1)[trim:-trim],
             linewidth=3,color=colors[run],label=radius_labels[run])
plt.xlabel('time (hrs)',fontsize=18)
plt.xticks([0,10,20,30],[0,10,20,30],fontsize=18)
plt.ylabel('S/G2/M fraction',fontsize=18)
plt.yticks([0,.25,.5,.75],[0,.25,.5,.75],fontsize=18)
plt.legend(fontsize=14)
plt.savefig(folder+'_avgdivfrac.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_avgdivfrac.pdf',format = 'pdf',bbox_inches="tight")
plt.show()

#plot volume deviation as a function of time
plt.figure()
for run in np.arange(run_num):
    plt.plot(9.2*dt*np.arange(num_steps),150*8*var_area_all[run].T,linewidth=3,alpha=.1,color=colors[run])
    plt.plot(9.2*dt*np.arange(num_steps),150*8*np.mean(var_area_all[run].T,axis=1),
            linewidth=3,color=colors[run],label=radius_labels[run])
plt.xlabel('time (hrs)',fontsize=18)
plt.xticks([0,10,20,30],[0,10,20,30],fontsize=18)
plt.ylabel('volume deviation (in um3)',fontsize=18)
plt.yticks([0,200,400,600],[0,200,400,600],fontsize=18)
plt.legend(fontsize=14)
plt.savefig(folder+'_stdarea.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_stdarea.pdf',format = 'pdf',bbox_inches="tight")
plt.show()