
/*!
 * \file simulation.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Simulation class 
 */ 

#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#define XSTR(s) STR(s)
#define STR(s) #s


#include <iostream>
#include <string>
#include <vector>

#include "system.hpp"
#include "force_compute.hpp"
#include "integrate.hpp"
#include "topology.hpp"
#include "dump.hpp"
#include "version.hpp"

using std::cout;
using std::string;
using std::endl;
using std::to_string;
using std::vector;

namespace VMTutorial
{

  struct Simulation
  {
    Simulation(System& sys, Integrate& integ, ForceCompute& f, Topology& t) : _sys{sys}, 
                                                                              _integ{integ}, 
                                                                              _force_compute{f},
                                                                              _topology{t},
                                                                              print_freq{100},
                                                                              num_zeros{8},
                                                                              bar_width{40},
                                                                              sim_step{0},
                                                                              t{0},
                                                                              t1_last_increment{0}
                                                                              //print_eq{false}                                                                              
    { 

    }
    bool run(int, bool = true, bool = true);
    const string print_version() { return  "branch : "+static_cast<string>(XSTR(GIT_BRANCH))+" commit : "+static_cast<string>(XSTR(GIT_COMMIT_HASH)); }
    void progress_bar(double, const string&);
    void store_rng_state();
    void appendValueToCSV(int k);

    // variables and parameters
    System& _sys;
    Integrate& _integ;
    ForceCompute& _force_compute;
    Topology& _topology;
    //Dump& _dump = Dump(_sys,_force_compute);

    //int rescale_freq;
    int print_freq;
    int num_zeros;
    int bar_width;
    int sim_step;
    int t;                           // cumulative T1 count
    int t1_last_increment;           // T1 increment at the most recent step
    vector<int> t1_increment_history;
    vector<int> t1_total_history;
    //bool print_eq;
    
    //void set_print_eq(bool eq) {print_eq = eq;}   

    //void set_d(Dump& _d){_dump = _d;} 
  };

  void export_Simulation(py::module&);


}

#endif
