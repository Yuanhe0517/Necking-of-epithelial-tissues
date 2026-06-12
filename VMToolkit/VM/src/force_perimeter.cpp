/*!
 * \file force_perimeter.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief ForcePerimeter class 
*/ 

#include "force_perimeter.hpp"

namespace VMTutorial
{
  namespace
  {
    template <typename P>
    bool boundary_tension_applies(const HalfEdge<P>& he)
    {
      const auto& v_from = *he.from();
      const auto& v_to = *he.to();
      return he.edge()->boundary &&
             v_from.boundary && v_to.boundary &&
             (v_from.data().type_name == "regular" || v_to.data().type_name == "regular");
    }
  }

  // Vec ForcePerimeter::compute(const Vertex<Property>& v, const HalfEdge<Property>& he)
  // {
  //   Vec l = he.to()->r - v.r;                    // vector along the junction pointing away from the vertex
  //   const Face<Property>& f   = *(he.face());         // cell to the right of the half edge
  //   const Face<Property>& fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
  //   double P1 = _sys.mesh().perim(f);
  //   double P2 = _sys.mesh().perim(fp);
  //   double lambda_1, lambda_2;
    
  //   double  gamma_1 = (f.outer)    ? 0.0 : _gamma[f.data().face_type];
  //   double  gamma_2 = (fp.outer)   ? 0.0 : _gamma[fp.data().face_type];
  //   if (!_lambda_P0)
  //   {
  //     lambda_1 = (f.outer)   ? 0.0 : _lambda[f.data().face_type];
  //     lambda_2 = (fp.outer)  ? 0.0 : _lambda[fp.data().face_type];
  //   }
  //   else 
  //   {
  //     lambda_1 = gamma_1*f.data().P0;
  //     lambda_2 = gamma_2*fp.data().P0;
  //   }

  //   double lambda = lambda_1 + lambda_2;
  //   double fedges = gamma_1*P1 + gamma_2*P2 - lambda;

  //   return fedges*l.unit();
  // }

  Vec ForcePerimeter::compute(const Vertex<Property>& v, const HalfEdge<Property>& he)
  {
    Vec l = he.to()->r - v.r;    
    Vec l1 = he.prev()->from()->r - v.r ;                // vector along the junction pointing away from the vertex
    const Face<Property>& f   = *(he.face());         // cell to the right of the half edge

    double P1 = _sys.mesh().perim(f);

    double lambda_1;//, lambda_2;
    
    double  gamma_1 = (f.outer)    ? 0.0 : _gamma[f.data().face_type];
    double boundary_tension = (f.outer || !_boundary_tension_enabled) ? 0.0 : _boundary_tension;

    if (!_lambda_P0)
    {
      lambda_1 = (f.outer)   ? 0.0 : _lambda[f.data().face_type];
    }
    else 
    {
      lambda_1 = gamma_1*f.data().P0;
    }

    double fedges = gamma_1*P1  - lambda_1;

    Vec fvec = fedges*(l.unit()+l1.unit());

    if (boundary_tension_applies(he))
    {
      fvec += boundary_tension * l.unit();
    }

    if (boundary_tension_applies(*he.prev()))
    {
      fvec += boundary_tension * l1.unit();
    }

    return fvec;
  }

  double ForcePerimeter::tension(const HalfEdge<Property>& he)
  {
    const Face<Property>& f   = *(he.face());         // cell to the right of the half edge
    const Face<Property>& fp  = *(he.pair()->face()); // pair cell (opposite side of the same junction)
    
    double lambda_1, lambda_2;
    
    double gamma_1 = (f.outer)    ? 0.0 : _gamma[f.data().face_type];
    double gamma_2 = (fp.outer)   ? 0.0 : _gamma[fp.data().face_type];
    if (!_lambda_P0)
    {
      lambda_1 = (f.outer)   ? 0.0 : _lambda[f.data().face_type];
      lambda_2 = (fp.outer) ? 0.0 : _lambda[fp.data().face_type];
    }
    else
    {
      lambda_1 = gamma_1*f.data().P0;
      lambda_2 = gamma_2*fp.data().P0;
    }

    double lambda = lambda_1 + lambda_2;
  
    return gamma_1*_sys.mesh().perim(f) + gamma_2*_sys.mesh().perim(fp) - lambda;
  }

  double ForcePerimeter::energy(const Face<Property>& f)
  {
    if (f.outer || f.erased)
      return 0.0;

    double P = _sys.mesh().perim(f);
    
    double lambda = 0.0;
    double gamma = _gamma[f.data().face_type];
    if (!_lambda_P0)
        lambda = _lambda[f.data().face_type];
    
    double P0;
    if (_lambda_P0)
      P0 = f.data().P0;
    else
      P0 = lambda/gamma;
    
    double dP = P - P0;

    double E = 0.5*gamma*dP*dP;

        
    if (_boundary_tension_enabled)
    {
      for (auto he : f.circulator())
      {
        if (boundary_tension_applies(he))
          E += _boundary_tension * (he.to()->r - he.from()->r).len();
      }
    }

    return E;
    
  }

}
