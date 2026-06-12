/*!
 * \file force_active_edge_tension.cpp
 * \author OpenAI
 * \date 01-Apr-2026
 * \brief Active edge tension force with per-edge phase
 */

#include "force_active_edge_tension.hpp"

namespace VMTutorial
{

  Vec ForceActiveEdgeTension::compute(const Vertex<Property>& v, const HalfEdge<Property>& he)
  {
    const Edge<Property>& e = *(he.edge());
    if (e.erased)
      return Vec(0.0, 0.0);

    Vec l = he.to()->r - v.r;
    double len = l.len();
    if (len < SMALL_NUMBER)
      return Vec(0.0, 0.0);

    return edge_gamma(e) * l.unit();
  }

  double ForceActiveEdgeTension::tension(const HalfEdge<Property>& he)
  {
    const Edge<Property>& e = *(he.edge());
    if (e.erased)
      return 0.0;
    return edge_gamma(e);
  }

}
