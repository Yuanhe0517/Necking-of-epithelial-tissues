/*!
 * \file integrate.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \modified by Yuan He, yuanhe@westlake.edu.cn
 * \date 30-Nov-2023 <Modification Date 22-Oct-2025>
 * \brief Integrate class 
 */
 
#ifndef __INTEGRATE_COMPUTE_HPP__
#define __INTEGRATE_COMPUTE_HPP__

#include <exception>
#include <algorithm>
#include <vector>
#include <string>
#include <mpi.h>
#include <omp.h>

#include "system.hpp"
#include "force_compute.hpp"
#include "integrator.hpp"
#include "integrator_brownian.hpp"

using std::runtime_error;
using std::transform;
using std::vector;
using std::string;
using std::map;


namespace VMTutorial
{

  class Integrate : public ClassFactory<Integrator>
  {
    public:

      Integrate(System& sys, ForceCompute& fc, int seed) : _sys{sys},
                                                           _force_compute{fc},
                                                           _dt{0.01}, 
                                                           //_alpha{0.1},
                                                           //_dt_max{1},
                                                           //_power{0},
                                                           _seed{seed},
                                                           _meq{false},
                                                           mpi_size(4),
                                                           mpi_rank(4)
                                                           { 
          
                                                           } 
      ~Integrate() = default; 

      Integrate(const ForceCompute&) = delete;

      //Integrate() : comm_global(MPI_COMM_WORLD){}

      void apply()
      {
        // int rank, size;
        // MPI_Comm_rank(comm_global, &rank);
        // MPI_Comm_size(comm_global, &size);
        //cout<<rank<<endl;
        for (auto i : this->integ_order)
          if (this->factory_map[i]->is_enabled())
          {
            //cout<<"step"<<endl;
            this->factory_map[i]->step();
            //cout<< "finish step" <<endl;
            if (this->factory_map[i]->meq)
            {
              _meq = true;
              this->factory_map[i]->set_meq(false);
            }
          }
        _sys.time_step()++;
        _sys.simulation_time() += _dt;
      }

      void set_params(const string& iname, const params_type& params)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_params(params);
        else
          throw runtime_error("set_params: Integrator type " + iname + " is not used in this simulation.");
      }

      void set_type_params(const string& iname, const string& type, const params_type& params)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_type_params(type, params);
        else
          throw runtime_error("set_type_params: Integrator type " + iname + " is not used in this simulation.");
      }

      void set_string_params(const string& iname, const string_params_type& params)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_string_params(params);
        else
          throw runtime_error("set_string_params: Integrator type " + iname + " is not used in this simulation.");
      }

      void set_external_force(const string& iname, const string& vtype, const Vec& f)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_external_force(vtype, f);
        else
          throw runtime_error("set_external_force: Integrator type " + iname + " is not used in this simulation.");
      }


      void set_flag(const string& iname, const string& flag)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->set_flag(flag);
        else
          throw runtime_error("set_flag: Integrator type " + iname + " is not used in this simulation.");
      }

      void enable(const string& iname)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->enable();
        else
          throw runtime_error("enable: Integrator type " + iname + " is not used in this simulation.");
      }

      void disable(const string& iname)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->disable();
        else
          throw runtime_error("disable: Integrator type " + iname + " is not used in this simulation.");
      }

      void set_dt(double dt) 
      { 
        if (this->factory_map.size() < 1)
          throw runtime_error("There are no integrators defined. Time step cannot be changed.");
        _dt = dt; 
        for (auto& integ : this->factory_map)
          integ.second->set_dt(dt);
      }

      void enable_constraint(const string& iname, bool enable)
      {
        if (this->factory_map.find(iname) != this->factory_map.end())
          this->factory_map[iname]->enable_constraint(enable);
        else
          throw runtime_error("enable_constraint: Integrator type " + iname + " is not used in this simulation.");
      }

      
      void add_integrator(const string& iname)
      {
        string name = iname; 
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        if (name == "brownian")
          this->add<IntegratorBrownian, System&, ForceCompute&, int, int, int>(name, _sys, _force_compute, _seed, mpi_size,mpi_rank);
        else 
          throw runtime_error("Unknown integrator type : " + name + ".");
        integ_order.push_back(name);
      }

      void set_eq(bool eq){_meq=eq;}

      void set_rank(int rank){mpi_rank = rank;}

      void set_size(int size){mpi_size = size;}

      void say_hi() {
        int argc = 0;
        char** argv = nullptr;
        MPI_Init(&argc, &argv);
        //cout<<"hello2"<<endl;
        //cout<<comm_global<<endl;
        int world_size;
        MPI_Comm_size(comm_global, &world_size);
        int world_rank;
        MPI_Comm_rank(comm_global, &world_rank);
        char processor_name[MPI_MAX_PROCESSOR_NAME] = "localhost";
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);
        printf("[C++]    Hello from machine %s, MPI rank %d out of %d\n",
        processor_name,
        world_rank,
        world_size);
      }

      bool _meq;
      double _dt;

    private: 

      System& _sys;
      ForceCompute& _force_compute;
      
      int _seed;
      MPI_Comm comm_global;

      vector<string> integ_order;  // Keeps track in which order integrators were added

      int mpi_rank;
      int mpi_size;

  };

  void export_Integrate(py::module&);

}
#endif


  
