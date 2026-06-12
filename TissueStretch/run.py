# ###############################################################
#
#  Performs simulations of horizontal tissue stretch
#
#  Author: Rastko Sknepnek, (c) 2023
#  Modified by: Yuan He, (c) 2025
#  Date: 09-Dec-2023 (Original), <Modification Date 22-Oct-2025>
#
# ###############################################################


# ###############################################################
#
# Load standard Python modules
#
# ###############################################################
import argparse 

# ###############################################################
#
# Load AJM modules
#
# ###############################################################
from VMToolkit.VM import *
#from tissue_plot import make_plotter

# ###############################################################
#
# Read command line arguments 
#
# ###############################################################
parser = argparse.ArgumentParser()  
parser.add_argument('--input', dest = 'input', type = str, default = 'honeycomb.json' , help = 'input fule')
parser.add_argument('--dt', dest = 'dt', type = float, default = 0.1, help = 'timestep')
parser.add_argument('--seed', dest = 'seed', type = int, default = None, help = 'random number generator seed')
parser.add_argument('--nrun', dest = 'nrun', type = int, default = 4000, help = 'number of run steps')
parser.add_argument('--loading-increment', dest = 'loading_increment', type = float, default = 0.03, help = 'boundary displacement applied after each relaxed state is saved')
args = parser.parse_args()


# ###############################################################
#
# Set parameters
#
# ###############################################################
freq = int(round(10.0/args.dt))   # This makes sure that we output data once per unit of time

# Cell mechanics parameters
kappa = 1.0     # area stiffness
gamma = 0.16     # perimeter stiffness
lam = 0.56 #0.592 # lambda parameter
boundary_tension = 0
# ################################################################
#
# Set up simulation objects
#
# ################################################################

if args.seed != None:
    seed = args.seed
else:
    seed = -1        # Use current time as the seed

tissue  = Tissue()                                               # initialise mesh
sim_sys = System(tissue)                                         # base object for the system
forces = Force(sim_sys)                                          # handles all types of forces
integrators = Integrate(sim_sys, forces, seed)                   # handles all integrators
topology = Topology(sim_sys, forces)                             # handles all topology changes (T1, division, ingression)
dumps = Dump(sim_sys, forces)                                    # handles all data output 
simulation = Simulation(sim_sys, integrators, forces, topology)  # simulation object
#IntegratorBrownian(sim_sys,forces,seed)

# #################################################################
#
# Create the initial configuration and read it
#
# #################################################################

sim_sys.read_input(args.input)           # read input configuration


# #################################################################
#
# Add forces to the system
#
# #################################################################


forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)


# Set parameters for each cell type
forces.set_params('area', 'passive' , {'kappa' : kappa})
forces.set_params('perimeter', 'passive', {
    'gamma': gamma,
    'lambda': lam,
    'boundary_tension': boundary_tension,
})
   


# #################################################################
#
# Set conditions for the T1 transition
#
# #################################################################

topology.set_params({'min_edge_len': 0.1, 'new_edge_len': 0.12}) 


# #################################################################
#
# Add Brownian integrator that will handle mechanical part
#
# #################################################################

integrators.add('brownian')    
integrators.set_params('brownian', {'loading_increment': args.loading_increment})

# #################################################################
#
# Simulation starts here
#
# #################################################################

integrators.set_dt(args.dt) # set time step

step = 0       # Step counter in terms of time units

print('Pull for {:d} steps'.format(args.nrun))


for i in range(args.nrun):
    if i < 6:
        n = 4000#50000
    else:
        n = 3000

    if simulation.run(n):
        #print(f'Simulation  at step {i+1} .')
        #dumps.dump_junctions(f'junctions_fpull_{args.fpull:.4f}_step_{i:08d}.vtp')
        dumps.dump_cells(f'cells_step_{i:08d}.vtp')
    # if i == 200:
    #     dumps.dump_json('200 steps.json')
    # if i == 400:
    #     dumps.dump_json('400 steps.json')
    if i == 600:
        dumps.dump_json('600 steps.json')
    if i == 1200:
        dumps.dump_json('1200 steps.json')
    if i == 1800:
        dumps.dump_json('1800 steps.json')
    if i == 2400:
        dumps.dump_json('2400 steps.json')
    if i == 3000:
        dumps.dump_json('3000 steps.json')
   
dumps.dump_json(f'final steps.json')
