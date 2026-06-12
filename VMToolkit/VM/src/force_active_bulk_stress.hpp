/*!
 * \file force_active_bulk_stress.hpp
 * \author OpenAI
 * \date 01-Apr-2026
 * \brief Active bulk stress force using sigma^(act) = -beta Q
 */

#ifndef __FORCE_ACTIVE_BULK_STRESS_HPP__
#define __FORCE_ACTIVE_BULK_STRESS_HPP__

#include "force.hpp"

namespace VMTutorial
{

  class ForceActiveBulkStress : public Force
  {
    public:

      ForceActiveBulkStress(System& sys) : Force{sys}
      {
        _beta.resize(_sys.cell_types().size(), 0.0);
      }

      virtual ~ForceActiveBulkStress() { }

      Vec compute(const Vertex<Property>&, const HalfEdge<Property>&) override;

      double tension(const HalfEdge<Property>&) override
      {
        return 0.0;
      }

      double energy(const Face<Property>&) override
      {
        return 0.0;
      }

      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "beta")
            throw runtime_error("Unknown parameter " + p.first + ".");

        if (params.find("beta") == params.end())
          throw runtime_error("ActiveBulkStress force requires parameter beta.");

        try
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("For active bulk stress: Cell type " + cell_type + " is not defined.");
          if (_beta.size() < _sys.get_num_cell_types())
            _beta.resize(_sys.get_num_cell_types(), 0.0);
          int ct = _sys.cell_types()[cell_type];
          _beta[ct] = params.at("beta");
        }
        catch (const exception& e)
        {
          cerr << "Problem with setting active bulk stress parameters. Exception: " << e.what() << '\n';
          throw;
        }
      }

      void set_vec_params(const string&, const vec_params_type&) override { }

      void set_flag(const string&) override { }

    private:

      void shape_tensor(const Face<Property>& f, double& Qxx, double& Qxy, double& Qyy) const;

      vector<double> _beta;
  };

}

#endif
