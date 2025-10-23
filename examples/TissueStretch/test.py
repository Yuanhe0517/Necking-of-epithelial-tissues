from mpi4py import MPI
import sys
sys.path.append('/home/ubuntu/Desktop/VMTutorial')
from VMToolkit.VM import *

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
print("[Python] Hello from machine " + name + ", MPI rank " + str( rank ) + " out of " + str( size ))


seed = -1
tissue  = Tissue()                                               # initialise mesh
sim_sys = System(tissue)                                         # base object for the system
forces = Force(sim_sys)                                          # handles all types of forces
integrators = Integrate(sim_sys, forces, seed)                   # handles all integrators

#print(comm)
integrators.say_hi()