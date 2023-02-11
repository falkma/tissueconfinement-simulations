import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import linregress
matplotlib.rcParams.update({'font.size': 8})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams["font.sans-serif"] = "Arial"

#script for generating single-cell birth size-division size correlation
#statistics of a G1-sizer model at varying (linear) growth rates
#code for Fig 5E
folder = 'birth_div_correlation_linear_growth'

#define timesteps and dt
num_steps = 8000
dt = .005
unit_time = int(1/dt)

#choose number of replicates
replicates = 400

#set initial growth of cells to provide baseline
initial_growth = 1.

#set growth rates for parameter sweep
growth_rates = np.array([1.,.9*3./4.,.9*1./4.,1./20.])
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
    #we keep track of size of the cell,
    #cell cycle phase, mother size, and daughter size

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
                    #convert approx to experimental units - 1 unit sim volume = 1200 um3
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
num_divisions = np.zeros(growth_num)

#we set the division rate to obtain a desired average rate total cell cycle length
#consistent with experiment
division_s = 3

#perform the parameter sweep over growth rates
for i in np.arange(growth_num):

    print(i)
    area_s = growth_rates[i]
    a_trace, p_trace, m_sweep, d_sweep = evaluate_behavior(area_s,division_s,initial_growth)

    #we use the phase data after quench in order to count the number of division events
    p_trace = np.array(p_trace)[:,20*unit_time:]
    count = len(np.where(p_trace == 0)[0])

    #we compare divisions to the number we would expect from unconfined growth
    num_divisions[i] = count/replicates/(p_trace.shape[1]/unit_time/(1+1/division_s))

    area_sweep.append(a_trace)
    mother_sweep.append(m_sweep)
    daughter_sweep.append(d_sweep)

colors = ['black','purple','fuchsia','orange']

#estimate confinement to first decimal
confinement = np.round(10*(1-num_divisions))/10
fig = plt.figure(figsize=(4.5,6))
for i in np.arange(len(growth_rates)):

    #linear fit to scatter of birth size vs division size
    coeffs = linregress(daughter_sweep[i][:-1],mother_sweep[i][1:])[0:2]

    #plot birth size vs division size scatter
    plt.scatter(daughter_sweep[i][:-1],mother_sweep[i][1:],alpha=.1,
             color=colors[i],s=4)

    size_range = np.linspace(0.,np.max(daughter_sweep[0][:-1]))

    #output data to csv
    a = np.array([daughter_sweep[i][:-1],mother_sweep[i][1:]]).T
    np.savetxt(folder+'_confinement_%s.csv' % confinement[i], a)

    #plot linear fit on top of scatter
    plt.plot(size_range,coeffs[0]*size_range+coeffs[1],
             label='confinement = %s \nslope = %s' % (confinement[i],np.round(10*coeffs[0])/10),
             color=colors[i],
             linewidth=2,linestyle='--')

#plot all growth rates in the same figure
plt.ylabel(r'division volume $(um^3)$',fontsize=13)
plt.xlabel(r'birth volume $(um^3)$',fontsize=13)
plt.xlim([0.,3000])
plt.ylim([0.,4000])
plt.xticks([0,1000,2000,3000],[0,1000,2000,3000],fontsize=13)
plt.yticks([0,1000,2000,3000,4000],[0,1000,2000,3000,4000],fontsize=13)
plt.legend(fontsize=10,loc='lower right')
plt.savefig(folder+'_birthdivsize.png',format = 'png',bbox_inches="tight")
plt.savefig(folder+'_birthdivsize.pdf',format = 'pdf',bbox_inches="tight")
plt.show()