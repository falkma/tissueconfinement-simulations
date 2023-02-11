import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import linregress
matplotlib.rcParams.update({'font.size': 8})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams["font.sans-serif"] = "Arial"

#script for generating single-cell statistics of a G1-sizer model
#at varying growth rates
#code for SI Fig 7

folder = 'single_cell_comparison_to_theory'

#define timesteps and dt
num_steps = 8000
dt = .005
unit_time = int(1/dt)

#choose number of replicates
replicates = 100

#set initial growth of cells to provide baseline
initial_growth = 1.

#set growth rates for parameter sweep
growth_rates = .5*2.**np.array([-5, -4.5, -4, -3.5, -3,
                          -1, -.5,  0,  .5,   1, 
                          3,  3.5,  4,  4.5,  5])
growth_num = len(growth_rates)

def reaction_rate(area,k):

    #division rate as a function of area

    if area > 1:
        r = k
    else:
        r = 0

    return r

def evaluate_behavior(area_rate,division_rate,initial_area_rate):

    #accumulate statistics for a given division rate and growth rate
    #we keep track of size/area of the cell
    #cell cycle phase
    #mother size and daughter size

    area_traces = []
    phase_traces = []
    mother_size = []
    daughter_size = []

    for _ in np.arange(replicates):

        #cells are initialized with random size
        #and a phase = -1, which means they have just divided

        a = [initial_area_rate*(1 + np.random.uniform(0,1))]
        
        #the phase variable represents cell cycle phase:
        #phase < 0 indicates G1, phase > 0 indicates S/G2, division occurs when phase = 0
        phase = [-1]

        for j in np.arange(num_steps):

            #we now let cells grow and divie for num_steps

            if j > 10*unit_time:

                #only after a certain time period do we allow cells
                #to drop to our desired growth rate
                a.append(a[-1] + area_rate*dt)

            else:

                #initially we equilibrate cells via a high initial growth rate
                a.append(a[-1] + initial_area_rate*dt)

            #every time step, the phase decreases by 1
            phase.append(phase[-1] - 1)

            #we generate a random number to see if cells enter S
            r = reaction_rate(a[-1],division_rate)*dt
            if r > np.random.uniform():
                if phase[-1] < 0:

                    #a cell is allowed to enter S if it is in G1 phase (i.e. phase < 0)
                    #if a cell enters S its phase is set to a counter to enforce deterministic S length
                    #this amount of time is 1 unit of simulation time, which sets the timescale
                    #for comparison to experiment
                    phase[-1] = unit_time

            if phase[-1] == 0:

                #when a cell's phase variable reaches 0, it divides

                if j > 20*unit_time:

                    #after equilibrating from the growth quench
                    #we begin to record the cell size pre- and post-division
                    mother_size.append(150*8*a[-1])
                    a[-1] = a[-1]/2.
                    daughter_size.append(150*8*a[-1])
                else:
                    a[-1] = a[-1]/2.

        area_traces.append(a)
        phase_traces.append(phase)

    return area_traces, phase_traces, mother_size, daughter_size


#initializing storage variables for the growth rate sweep
area_sweep = []
mother_sweep = []
daughter_sweep = []
confinement = np.zeros(growth_num)
divtime = np.zeros(growth_num)

#we set the division rate to obtain an average rate of S entrance for unconfined growth
division_s = 3

#perform the parameter sweep over growth rates
for i in np.arange(growth_num):

    print(i)
    area_s = growth_rates[i]
    a_trace, p_trace, m_sweep, d_sweep = evaluate_behavior(area_s,division_s,initial_growth)

    #we use the phase data after quench in order to count the number of division events
    p_trace = np.array(p_trace)[:,10*unit_time+1000:]
    count = len(np.where(p_trace == 0)[0])

    #we estimate confinement based on a guess of the unconfined rate compared to the measured
    #rate of division
    confinement[i] = 1-count/replicates/(p_trace.shape[1]/unit_time/(1+1/division_s))
    divtime[i] = p_trace.shape[1]/unit_time/(count/replicates)

    a_trace = np.array(a_trace)[:,10*unit_time+1000:]

    area_sweep.append(a_trace)
    mother_sweep.append(m_sweep)
    daughter_sweep.append(d_sweep)

#all plots in non-dimensional simulation units

#plot confinement as a function of growth
#compare to analytic forms for different growth regimes
fig = plt.figure(figsize=(7,2))
plt.plot(np.log10(growth_rates),confinement,color='black',label='numerical')
plt.plot(np.log10(growth_rates[6:]),np.zeros(len(growth_rates[6:])),linestyle='--',color='fuchsia',label='theory limits')
plt.plot(np.log10(growth_rates[0:10]),1-(1+1/division_s)/(1/2+1/2/division_s+1/2/growth_rates[0:10]),linestyle='--',color='fuchsia')
plt.ylabel('confinement',fontsize=13)
plt.xlabel('log(growth rate)',fontsize=13)
plt.savefig(folder+'_confinement.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_confinement.pdf',format = 'pdf',bbox_inches="tight")
plt.show()

#plot average division time as a function of growth
#compare to analytic forms for different growth regimes
fig = plt.figure(figsize=(7,2))
plt.plot(np.log10(growth_rates),np.log10(divtime),color='black',label='numerical')
plt.plot(np.log10(growth_rates[5:]),np.log10((1+1/division_s)*np.ones(len(growth_rates[5:]))),linestyle='--',color='fuchsia',label='theory limits')
plt.plot(np.log10(growth_rates[0:11]),np.log10(1/2+1/2/division_s+1/2/growth_rates[0:11]),linestyle='--',color='fuchsia')
plt.ylabel(r'log($\tau$)',fontsize=13)
plt.xlabel('log(growth rate)',fontsize=13)
plt.savefig(folder+'_avgdivtime.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_avgdivtime.pdf',format = 'pdf',bbox_inches="tight")
plt.show()

#plot average cell size as a function of growth
#compare to analytic forms for different growth regimes
fig = plt.figure(figsize=(7,2))
plt.plot(np.log10(growth_rates),np.log10(np.mean(np.mean(area_sweep,axis=2),axis=1)),color='black',label='numerical')
plt.plot(np.log10(growth_rates[5:]),np.log10(1.5*(1+1/division_s)*growth_rates[5:]),linestyle='--',color='fuchsia',label='theory limits')
plt.plot(np.log10(growth_rates[0:11]),np.log10(.75*(1+(1+1/division_s)*growth_rates[0:11])),linestyle='--',color='fuchsia')
plt.ylabel('log(mean area/a)',fontsize=13)
plt.xlabel('log(growth rate)',fontsize=13)
plt.legend()
plt.savefig(folder+'_avgarea.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_avgarea.pdf',format = 'pdf',bbox_inches="tight")
plt.show()