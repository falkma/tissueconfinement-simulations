## Microscopic manipulation with capillary machines (numerics)

The code uploaded here is for replicating the simulation results of Devaney et al., "Tissue confinement regulates cell growth and size in epithelia". All code is written in Python (v3.7).

The scripts in this repository implement agent-based simulations for the G1-Sizer model discussed in the paper. A brief description is provided below; please find more detail in the Methods section of the paper.

### G1-Sizer model overview

Briefly, our phenomenological model of cell size control is a “G1 Sizer”, which posits that cells' exit from G1 into S/G2/M is controlled by a stochastic, size-dependent rate function. Following G1-exit, division occurs after a deterministic period of time. Growth can occur across the entire cell cycle. These two elements – size-dependent S entrance, and growth of cells – constitute the core of how we advance single-cell trajectories through time. 

We simulate this model with an agent-based approach - each cell carries an index i, as well as two quantities a_i and p_i, representing the size and cell cycle phase respectively. Time and size are in non-dimensional units during the simulation and are converted to dimensional units for analysis after the simulation. We do this by associating 1 unit of simulation time to the length of a typical S phase (~10 hours), and 1 unit of simulation size to be the volume at which it is seen experimentally that cells transition from a size-dependent to a size-independent division rate (~1200 um3).

We implement two categories of simulations – ensemble simulations, and single-cell simulations. In single-cell simulations, we track the trajectory of only one cell, following only one daughter cell after division. In ensemble simulations, we track a whole population of cells, and the growth rate can therefore depend on quantities like N, the total number of cells in the population.

### Single-cell simulations

- single_cell_comparison_to_theory.py: compares statistics of single-cell trajectories to predictions from analytical expressions.
- birth_div_correlation_linear_growth.py: computes correlation between birth size and division size for single-cell trajectories which grow with a linear growth rate.
- birth_div_correlation_exp_growth.py: computes correlation between birth size and division size for single-cell trajectories which grow with a exponential growth rate.
- birth_div_correlation_smooth_rates.py: computes correlation between birth size and division size for single-cell trajectories which exit G1 with a rate function derived from a fit to experimental data.

### Ensemble simulations

- confinement_quench.py: computes statistics for colonies of cells which experience a quench in confinement at varying initial colony sizes. The colony first grows in unconfined conditions, and after being quenched, grows with colony radius expanding outwards at a constant radial velocity.