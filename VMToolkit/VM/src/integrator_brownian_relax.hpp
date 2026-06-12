/*!
 * \file integrator_brownian_relax.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief IntegratorBrownianRelax class
*/

#ifndef __INTEGRATOR_BROWNIAN_RELAX_HPP__
#define __INTEGRATOR_BROWNIAN_RELAX_HPP__

#include "constrainer.hpp"
#include "constraint_none.hpp"
#include "constraint_fixed.hpp"


#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace std::chrono;
using std::map;
using std::make_unique;
using std::string;

namespace VMTutorial
{

  using ConstrainerType = unique_ptr<Constrainer>;

  class IntegratorBrownianRelax : public Integrator
  {

    public:

      IntegratorBrownianRelax(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc, seed},
                                                                         _T{0.0},
                                                                         _gamma{1.0},
                                                                         _omega0{0.0},
                                                                         _J{0.0},
                                                                         _Tphi{0.0},
                                                                         _A0_amp{0.0},
                                                                         _Dr{0.0},
                                                                         _update_n{false}

      { 
        map<string,int>& vert_types = _sys.vert_types();
        map<string,int>& cell_types = _sys.cell_types();
        for (int i = 0; i < vert_types.size(); i++)  _constant_force.push_back(Vec(0.0,0.0));

        _constrainer = make_unique<Constrainer>();
        _constrainer->add<ConstraintNone>("none");
        _constrainer->add<ConstraintFixed>("fixed");
      }

      void step() override;
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "T")
            _T = p.second;
          if (p.first == "gamma")
            _gamma = p.second;
          if (p.first == "omega0")
            _omega0 = p.second;
          if (p.first == "J")
            _J = p.second;
          if (p.first == "Tphi")
            _Tphi = p.second;
          if (p.first == "A0_amp")
            _A0_amp = p.second;
          if (p.first == "Dr")
            _Dr = p.second;
        }
      };
      void set_type_params(const string& type, const params_type& params) override { }
      void set_string_params(const string_params_type& params) override
      {
        for (auto& p : params)
        {
          if (p.first != "phi_init")
            throw runtime_error("Brownian relax integrator: Unknown string parameter : " + p.first + ".");

          if (p.second != "zero" && p.second != "random")
            throw runtime_error("Brownian relax integrator: phi_init must be either 'zero' or 'random'.");

          _phi_init_mode = p.second;
          _phi_initialized = false;
          _cell_phi_initialized = false;
        }
      }
      void set_external_force(const string& vtype, const Vec& f) override 
      { 
        _constant_force[_sys.vert_types()[vtype]] = f;
      }
      
      void set_flag(const string& flag) override 
      {  
        if (flag == "update_n")
          _update_n = true;
        else
          throw runtime_error("Brownian relax integrator: Unknown flag : " + flag + ".");
      }

      
    private:

      ConstrainerType _constrainer; // Apply various constraints
      double _T;                 // temperature
      double _gamma;             // friction 
      double _omega0;            // phase drift frequency
      double _J;                 // Kuramoto phase coupling
      double _Tphi;              // phase noise strength
      double _A0_amp;            // amplitude of cell target-area oscillation
      vector<Vec> _constant_force;
      double _Dr;                // rotational diffusion constant
      bool _update_n;            // If true, update direction of the cell self-propulsion direction
      string _phi_init_mode = "zero";   // initialize edge phase to zero or random
      bool _phi_initialized = false;    // edge phases are initialized lazily at the first step
      bool _cell_phi_initialized = false; // cell phases are initialized lazily at the first step
      
    };

}

#endif
