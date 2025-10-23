/*!
 * \file force_area.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceArea class
 */

#include "force_area.hpp"

namespace VMTutorial
{
	// Vec ForceArea::compute(const Vertex<Property> &v, const HalfEdge<Property> &he)
	// {
	// 	Vec l = he.to()->r - v.r;
	// 	Vec l1 = v.r - he.prev()->from()->r;
	// 	//cout<<v.r<<v.id<<endl;
	// 	//cout<<l<<l1<<endl;					// vector along the junction pointing away from the vertex
	// 	const Face<Property> &f = *(he.face());			// cell to the right of the half edge
	// 	const Face<Property> &fp = *(he.pair()->face()); // pair cell (opposite side of the same junction)
	// 	double A1 = _sys.mesh().area(f);
	// 	double A2 = _sys.mesh().area(fp);
	// 	double A0_1 = f.data().A0;
	// 	double A0_2 = fp.data().A0;

	// 	double kappa_1 = (f.outer) ? 0.0 : _kappa[f.data().face_type];
	// 	double kappa_2 = (fp.outer) ? 0.0 : _kappa[f.data().face_type];
	// 	//cout<<kappa_1<<";"<<kappa_2<<endl;

	// 	Vec farea_vec = 0.5 * (kappa_1 * (A1 - A0_1) - kappa_2 * (A2 - A0_2)) * l.ez_cross_v();
	// 	cout<<v.id<<endl;

	// 	cout<<farea_vec<<endl;

	// 	return farea_vec;
	// }

	Vec ForceArea::compute(const Vertex<Property> &v, const HalfEdge<Property> &he)
	{
		Vec l = he.to()->r - v.r;	
		Vec l1 = v.r - he.prev()->from()->r;	
		Vec l2 = -l1-l;			// vector along the junction pointing away from the vertex
		const Face<Property> &f = *(he.face());			// cell to the right of the half edge
		double A1 = _sys.mesh().area(f);
		double A0_1 = f.data().A0;
		// if (!f.outer)
		// {
		// 	return 0;
		// }
		// else
		// {

		// }
		double kappa_1 = (f.outer) ? 0.0 : _kappa[f.data().face_type];
		//double kappa_2 = (fp.outer) ? 0.0 : _kappa[f.data().face_type];

		Vec farea_vec = -0.5 * (kappa_1 * (A1 - A0_1)) * l2.ez_cross_v();

		// cout<<v.id<<endl;

		// cout<<l<<l1<<l2.ez_cross_v()<<endl;

		// cout<<farea_vec<<endl;

		return farea_vec;
	}


	double ForceArea::energy(const Face<Property> &f)
	{
		double A = _sys.mesh().area(f);
		double A0 = f.data().A0;
		
		if (f.outer || f.erased)
			return 0.0;

		double kappa = _kappa[f.data().face_type];

		double dA = A - A0;
		return 0.5 * kappa * dA * dA;
	}
}
