/*!
 * \file force_active_bulk_stress.cpp
 * \author OpenAI
 * \date 01-Apr-2026
 * \brief Active bulk stress force using sigma^(act) = -beta Q
 */

#include "force_active_bulk_stress.hpp"
#include <cmath>

namespace VMTutorial
{

  void ForceActiveBulkStress::shape_tensor(const Face<Property>& f, double& Qxx, double& Qxy, double& Qyy) const
  {
    Qxx = 0.0;
    Qxy = 0.0;
    Qyy = 0.0;

    if (f.outer || f.erased)
      return;

    double P = _sys.mesh().perim(f);
    if (P < SMALL_NUMBER)
      return;

    for (auto he : f.circulator())
    {
      Vec e = he.to()->r - he.from()->r;
      double lk = e.len();
      if (lk < SMALL_NUMBER)
        continue;

      Vec t = e.unit();
      Qxx += lk * (t.x * t.x - 0.5);
      Qxy += lk * (t.x * t.y);
      Qyy += lk * (t.y * t.y - 0.5);
    }

    Qxx /= P;
    Qxy /= P;
    Qyy /= P;
  }

  Vec ForceActiveBulkStress::compute(const Vertex<Property>& v, const HalfEdge<Property>& he)
  {
    const Face<Property>& f = *(he.face());
    if (f.outer || f.erased)
      return Vec(0.0, 0.0);

    double beta = _beta[f.data().face_type];
    if (std::fabs(beta) < SMALL_NUMBER)
      return Vec(0.0, 0.0);

    double Qxx, Qxy, Qyy;
    shape_tensor(f, Qxx, Qxy, Qyy);

    // Tlili et al. scheme:
    // R_i^(T) = 1/2 (r_{i+1} - r_{i-1}) x zhat
    // Here ez_cross_v() gives zhat x v, hence the minus sign.
    Vec d = he.to()->r - he.prev()->from()->r;
    Vec R = -0.5 * d.ez_cross_v();

    // sigma^(act) = -beta Q, and F_i = -R_i dot sigma.
    // Since sigma is symmetric, this is equivalent to -sigma * R.
    double sxx = -beta * Qxx;
    double sxy = -beta * Qxy;
    double syy = -beta * Qyy;

    return Vec(-(sxx * R.x + sxy * R.y),
               -(sxy * R.x + syy * R.y),
               v.r.box);
  }

}
