/*!
 * \file integrator_brownian_relax.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief IntegratorBrownianRelax class
*/

#include <cmath>
#include <unordered_set>
#include <vector>

#include "integrator_brownian_relax.hpp"

namespace VMTutorial
{
  void IntegratorBrownianRelax::step()
  {
    constexpr double TWO_PI = 2.0 * 3.14159265358979323846;
    double mu = 1.0 / _gamma;    // mobility 
    double B = sqrt(2.0*mu*_T);
    double sqrt_dt = sqrt(_dt);
    double phase_noise = sqrt(2.0*_Tphi*_dt);

    if (!_phi_initialized)
    {
      if (_sys.edge_phi_loaded())
      {
        _phi_initialized = true;
        _sys.set_edge_phi_loaded(false);
      }
      else
      {
        for (auto& e : _sys.mesh().edges())
        {
          if (!e.erased)
          {
            if (_phi_init_mode == "random")
              e.data().phi = TWO_PI * _rng.drnd();
            else
              e.data().phi = 0.0;
          }
        }
        _phi_initialized = true;
      }
    }

    std::vector<double> phi_old(_sys.mesh().edges().size(), 0.0);
    for (auto& e : _sys.mesh().edges())
      if (!e.erased)
        phi_old[e.idx()] = e.data().phi;

    for (auto& e : _sys.mesh().edges())
    {
      if (!e.erased)
      {
        std::unordered_set<int> neighbour_ids;
        double coupling_sum = 0.0;

        for (auto he_n : e.he()->from()->circulator())
        {
          int nidx = he_n.edge()->idx();
          if (nidx != e.idx() && !he_n.edge()->erased)
            neighbour_ids.insert(nidx);
        }

        for (auto he_n : e.he()->to()->circulator())
        {
          int nidx = he_n.edge()->idx();
          if (nidx != e.idx() && !he_n.edge()->erased)
            neighbour_ids.insert(nidx);
        }

        for (int nidx : neighbour_ids)
          coupling_sum += std::sin(phi_old[nidx] - phi_old[e.idx()]);

        e.data().phi = phi_old[e.idx()] + _omega0 * _dt + _J * _dt * coupling_sum + phase_noise * _rng.gauss_rng();
        e.data().phi = fmod(e.data().phi, TWO_PI);
        if (e.data().phi < 0.0)
          e.data().phi += TWO_PI;
      }
    }

    if (!_cell_phi_initialized)
    {
      if (_sys.cell_phi_loaded())
      {
        _cell_phi_initialized = true;
        _sys.set_cell_phi_loaded(false);
      }
      else
      {
        for (auto& f : _sys.mesh().faces())
        {
          if (!(f.erased || f.outer))
          {
            if (_phi_init_mode == "random")
              f.data().phi = TWO_PI * _rng.drnd();
            else
              f.data().phi = 0.0;
          }
        }
        _cell_phi_initialized = true;
      }
    }

    std::vector<double> cell_phi_old(_sys.mesh().faces().size(), 0.0);
    for (auto& f : _sys.mesh().faces())
      if (!(f.erased || f.outer))
        cell_phi_old[f.id] = f.data().phi;

    for (auto& f : _sys.mesh().faces())
    {
      if (!(f.erased || f.outer))
      {
        double coupling_sum = 0.0;

        for (auto he : f.circulator())
        {
          const Face<Property>& fj = *(he.pair()->face());
          if (!(fj.erased || fj.outer))
            coupling_sum += std::sin(cell_phi_old[fj.id] - cell_phi_old[f.id]);
        }

        f.data().phi = cell_phi_old[f.id] + _omega0 * _dt + _J * _dt * coupling_sum + phase_noise * _rng.gauss_rng();
        f.data().phi = fmod(f.data().phi, TWO_PI);
        if (f.data().phi < 0.0)
          f.data().phi += TWO_PI;

        f.data().A0 = 1.0 + _A0_amp * std::sin(f.data().phi);
      }
    }

    // Compute force on each vertex
    for (auto& v : _sys.mesh().vertices())
      if (!v.erased)
        _force_compute.compute(v);
    
    // This is actual integrator 
    for (auto& v : _sys.mesh().vertices())
    {
      if (!v.erased)
      {
        // add external force 
        Vec f = v.data().force + _constant_force[v.data().vert_type];
        
        // apply constraint
        if (_constraint_enabled)
          f = _constrainer->apply_vertex(v, f);

        Vec rold = v.r;
        v.r += _dt*mu*f;  // deterministic part of the integrator step
        if (_T > 0.0)
        {
          Vec ffr(B*_rng.gauss_rng(), B*_rng.gauss_rng());  // random noise contribution to force
          Vec fr = ffr;
          if (_constraint_enabled) 
            fr = _constrainer->apply_vector(v, ffr);
          v.r += sqrt_dt*fr;  // update vertex position due to noise
        }
        v.data().vel = (1.0 / _dt) * (v.r - rold);  
      }
    }

    // Update direction of cell director using simple Brownian dynamics
    if (_update_n)
    {
      double stoch_coeff = sqrt(_Dr*_dt);
      for (auto& f : _sys.mesh().faces())
      {
        if (!(f.erased || f.outer))
        {
          double dtheta = stoch_coeff*_rng.gauss_rng();
          Vec n = f.data().n;
          // Sine and cosine of the rotation angle
          double c = cos(dtheta), s = sin(dtheta);
          // Rotation matrix around z axis
          double Rxx = c, Rxy = -s;
          double Ryx = s, Ryy = c;
          // Apply rotation matrix
          double nx = Rxx*n.x + Rxy*n.y;
          double ny = Ryx*n.x + Ryy*n.y;
          double len = sqrt(nx*nx + ny*ny);
          f.data().n = Vec(nx/len, ny/len);  // Rotation does not change length of a vector, but numerical errors can accumulate, so we normalize it
        }
      }
    }

  }

}
