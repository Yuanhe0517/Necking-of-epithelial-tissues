/*!
 * \file force_active_edge_tension.hpp
 * \author OpenAI
 * \date 01-Apr-2026
 * \brief Active edge tension force with edge phase dynamics handled by the integrator
 */

#ifndef __FORCE_ACTIVE_EDGE_TENSION_HPP__
#define __FORCE_ACTIVE_EDGE_TENSION_HPP__

#include <cmath>

#include "force.hpp"

namespace VMTutorial
{

  class ForceActiveEdgeTension : public Force
  {
    public:

      ForceActiveEdgeTension(System& sys) : Force{sys},
                                            _gamma0{0.0},
                                            _alpha{0.0}
      {
      }

      virtual ~ForceActiveEdgeTension() { }

      Vec compute(const Vertex<Property>&, const HalfEdge<Property>&) override;

      double tension(const HalfEdge<Property>&) override;

      double energy(const Face<Property>&) override
      {
        return 0.0;
      }

      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "gamma0" && p.first != "alpha")
            throw runtime_error("Unknown parameter " + p.first + ".");

        if (params.find("gamma0") == params.end())
          throw runtime_error("ActiveEdgeTension force requires parameter gamma0.");
        if (params.find("alpha") == params.end())
          throw runtime_error("ActiveEdgeTension force requires parameter alpha.");

        try
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("For active edge tension: Cell type " + cell_type + " is not defined.");
          _gamma0 = params.at("gamma0");
          _alpha = params.at("alpha");
        }
        catch (const exception& e)
        {
          cerr << "Problem with setting active edge tension parameters. Exception: " << e.what() << '\n';
          throw;
        }
      }

      void set_vec_params(const string&, const vec_params_type&) override { }

      void set_flag(const string&) override { }

    private:

      double edge_gamma(const Edge<Property>& e) const
      {
        //return _gamma0 * (1.0 + _alpha * std::sin(e.data().phi));
        return  _alpha * std::sin(e.data().phi);
      }

      double _gamma0;
      double _alpha;
  };

}

#endif
