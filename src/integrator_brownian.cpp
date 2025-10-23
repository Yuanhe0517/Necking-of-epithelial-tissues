/*!
 * \file integrator.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \modified by Yuan He, yuanhe@westlake.edu.cn
 * \date 30-Nov-2023 <Modification Date 22-Oct-2025>
 * \brief IntegratorBrownian class 
*/
#include <thread>
#include "integrator_brownian.hpp"

namespace py = pybind11;
using pymod = pybind11::module;

namespace VMTutorial
{
  void IntegratorBrownian::step()
  {
    double mu = 1.0 / _gamma;    // mobility 
    double B = sqrt(2.0*mu*_T);
    double sqrt_dt = sqrt(_dt);
    bool at_tensile_timestep = false;
    //double dt_max = 10 * _dt;  
    //double alpha = 0.1;  

    // Compute force on each vertex
    //cout<< "compute force" <<endl;

    for (auto& v : _sys.mesh().vertices()){
      if (!v.erased)
      {
        //cout<< "calculate"<<endl;
        //cout<< v.id<<endl;
        _force_compute.compute(v);
      }

    }
    
    // Compte energy of all faces in root process
    if (mpi_rank == 0) {
      //cout<< "energy"<<endl;
      double energy = _force_compute.total_energy();
      cout<<energy<<endl;

      // Initial Equilibrium
      if (mstep==0 && abs(energy - mEnergy) < 0.0000001)//0.000000001)
      {
        mstep = mstep + 1;
        at_tensile_timestep = true;
        std::cout<<"Equilibrium"<<mstep<<endl;
        meq = true;
        // Broadcast meq, at_tensile_timestep to all processes
        //MPI_Bcast(&meq, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //MPI_Bcast(&at_tensile_timestep, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }

      if (!mstep==0 && abs(energy - mEnergy) < 0.0000001)
      {
        mstep = mstep + 1;
        at_tensile_timestep = true;
        std::cout<<"Equilibrium"<<mstep<<endl;
        meq = true;
        // Broadcast meq, at_tensile_timestep to all processes
        //MPI_Bcast(&meq, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //MPI_Bcast(&at_tensile_timestep, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }

      mEnergy = energy;

      }  

      
      //MPI_Barrier(MPI_COMM_WORLD);
      // cout<<"apply force"<<endl;
      for (auto& v : _sys.mesh().vertices())
      {
        if (!v.erased)
        {
          // add external force for tensile?
          Vec f = v.data().force + 0*_constant_force[v.data().vert_type];
        
        
          // apply constraint after initial equilibrium, aka, mstep > 0
          if (_constraint_enabled && !mstep==0)
          {
            if (v.data().vert_type == 2 || v.data().vert_type == 0 )//&& !v.data().constraint == "fixed") 
            {
              v.data().constraint = "tensile";
              //if (!(v.data().vert_type == 0 && f.x > 0))
              f.x=0;
              //f.y=0;
            }
            else f = _constrainer->apply_vertex(v, f);
          }


          Vec rold = v.r;
          if (at_tensile_timestep)
          {
            // Apply tensile stretch
            if (v.data().vert_type == 0) v.r += Vec(0.0102,0.0);

            //if (v.data().vert_type == 2) v.r += Vec(-0.03,0.0);
            //cout<<v.r<<endl;
          }
          else
          {
            //v.r += _dt*f;   // deterministic part of the integrator step
            // Vervelt
            //Vec velocity = f;

            //velocity = (1-_alpha)*v.r+_alpha*f.unit()*v.r.len();
            v.r += _dt*f;
            //_power += velocity.dot(f);
          }


          if (_T > 1.0)
          {
            Vec ffr(B*_rng.gauss_rng(), B*_rng.gauss_rng());  // random noise contribution to force
            Vec fr = ffr;
            if (_constraint_enabled) 
              fr = _constrainer->apply_vector(v, ffr);
            v.r += sqrt_dt*fr*0.01;  // update vertex position due to noise
          }


          // v.data().vel = (1.0 / _dt) * (v.r - rold);  
        } 
      }


  }


}