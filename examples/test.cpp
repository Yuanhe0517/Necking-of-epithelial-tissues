#include <mpi.h>
#include <VMToolkit.VM>

int main(int argc, char **argv)
{   

    // # ###############################################################
    // #
    // # Read command line arguments 
    // #
    // # ###############################################################
    // parser = argparse.ArgumentParser();
    // parser.add_argument('--input', dest = 'input', type = str, default = 'honeycomb.json' , help = 'input fule');
    // parser.add_argument('--fpull', dest = 'fpull', type = float, default = 0.1, help = 'pulling force on left and right sides');
    // parser.add_argument('--dt', dest = 'dt', type = float, default = 0.05, help = 'timestep');
    // parser.add_argument('--seed', dest = 'seed', type = int, default = None, help = 'random number generator seed');
    // parser.add_argument('--dumpfreq', dest = 'dumpfreq', type = int, default = 100, help = 'how often to produce output');
    // parser.add_argument('--nrun', dest = 'nrun', type = int, default =1400, help = 'number of run steps');
    // parser.add_argument('--n',dest = 'n', tyype = int, default = 100, help = 'number of vertices');
    // args = parser.parse_args();
    double dt = 0.1;
    double nrun = 1400;
    // # ##########################
    // # MPI setting
    // # ##########################
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    //
    //  set up the data partitioning across processors
    //
    // int vertices_per_proc = (n + n_proc - 1) / n_proc;

    // // # ###############################################################
    // // #
    // // # Set parameters
    // // #
    // // # ###############################################################
    // int freq = int(round(10.0/dt))   # This makes sure that we output data once per unit of time

    // # Cell mechanics parameters
    // double kappa = 1.0;     # area stiffness
    // double gamma = 0.16;     # perimeter stiffness
    // double lam = 0.592;    # lambda parameter
    
    // // # ################################################################
    // // #
    // // # Set up simulation objects
    // // #
    // // # ################################################################

    // int seed = -1;        //# Use current time as the seed

    // Mymesh& tissue  = Tissue();                                               //# initialise mesh
    // System& sim_sys = System(tissue);                                         //# base object for the system
    // ForceCompute& forces = Force(sim_sys);                                          //# handles all types of forces
    // Integrate& integrators = Integrate(sim_sys, forces, seed, n_proc,rank);                   //# handles all integrators
    // Topology& topology = Topology(sim_sys, forces);                             //# handles all topology changes (T1, division, ingression)
    // Dump& dumps = Dump(sim_sys, forces);                                    //# handles all data output 
    // Simulation& simulation = Simulation(sim_sys, integrators, forces, topology);  //# simulation object

    // sim_sys.read_input('honeycomb.json');           //# read input configuration

    // forces.add('area');         //# add area force form term E = 0.5*kappa*(A-A0)^2
    // forces.add('perimeter');    //# add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)

    // //# Set parameters for each cell type
    // forces.set_params('area', 'passive' , {'kappa' : kappa});
    // forces.set_params('perimeter', 'passive' ,  {'gamma': gamma, 'lambda': lam};)

    // topology.set_params({'min_edge_len': 0.1, 'new_edge_len': 0.15}); 

    // integrators.add('brownian');  

    // integrators.set_dt(dt); //# set time step

    // int step = 0;       //# Step counter in terms of time units

    // for (int i in range(nrun))
    //     if (simulation.run(20000))
    //         //dumps.dump_junctions(f'junctions_fpull_{args.fpull:.4f}_step_{i:08d}.vtp')
    //         dumps.dump_cells(f"cells_fpull_{args.fpull:.4f}_step_{i:08d}.vtp");
    //     // if i == 250:
    //     //     dumps.dump_json('250 steps.json')
    
    // MPI.Finalize()
}