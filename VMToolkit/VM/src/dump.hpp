/*!
 * \file dump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Dump class
 */

#ifndef __DUMP_HPP__
#define __DUMP_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <exception>
#include <iomanip>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sys/stat.h>
#include <regex>

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

#include "json.hpp"

#include "system.hpp"
#include "force_compute.hpp"

using std::back_inserter;
using std::cerr;
using std::copy;
using std::ifstream;
using std::istream_iterator;
using std::istringstream;
using std::map;
using std::move;
using std::ofstream;
using std::runtime_error;
using std::setprecision;
using std::setw;
using std::stod;
using std::stoi;
using std::string;
using std::stringstream;
using std::to_string;
using std::tolower;
using std::transform;
using std::unique_ptr;
using std::vector;

namespace VMTutorial
{

	class Dump
	{
	public:
		Dump(System &sys, ForceCompute &fc) : _sys{sys}, _force_compute{fc}, _sfc(0.95) {}

		void dump_cells(const string &, bool = false, bool = false);
		void dump_junctions(const string &, bool = false);
		void dump_mesh(const string &, bool = false);
		void dump_json(const string &);

		void set_sfc(double sfc) { _sfc = sfc; }

		double get_face_stress_xx(const Face<Property> &f)
		{
			double stress_xx = 0.0; 
			double area = _sys.mesh().area(f);
			for (auto he : f.circulator())
			{
				VertexHandle<Property> v_handle = he.from();
        		Vertex<Property> &v = *v_handle; // 解引用 VertexHandle 获取 Vertex 对象

        		HalfEdge<Property> &he_obj = he; // 解引用 HEHandle 获取 HalfEdge 对象

				Vec node_location_ralate_to_centroid = v_handle->r - _sys.mesh().get_face_centre(f);

				// he_obj? 
				//v.circulator();f.circulator()
				Vec Force = _force_compute.compute(v, he_obj);//he);//v->data().force;
				stress_xx += -1.0/area*node_location_ralate_to_centroid.x*Force.x;
			}

			return stress_xx;
		}

		// double get_nominal_face_stress_xx(const Face<Property> &f)
		// {
		// 	double nominal_stress_xx = 0.0; 
		// 	double area = _sys.mesh().area(f);
		// 	for (auto he : f.circulator())
		// 	{
		// 		VertexHandle<Property> v_handle = he.from();
        // 		Vertex<Property> &v = *v_handle; // 解引用 VertexHandle 获取 Vertex 对象

        // 		HalfEdge<Property> &he_obj = he; // 解引用 HEHandle 获取 HalfEdge 对象

		// 		Vec node_location_ralate_to_centroid = v_handle->r - _sys.mesh().get_face_centre(f);

		// 		// he_obj? 
		// 		//v.circulator();f.circulator()
		// 		Vec Force = _force_compute.compute(v, he_obj);//he);//v->data().force;
		// 		stress_xx += -1.0/area*abs(node_location_ralate_to_centroid.x*Force.x);
		// 	}

		// 	return stress_xx;
		// }
		double get_face_stress_yy(const Face<Property> &f)
		{
			double stress_yy = 0.0; 
			double area = _sys.mesh().area(f);
			for (auto he : f.circulator())
			{
				VertexHandle<Property> v_handle = he.from();
        		Vertex<Property> &v = *v_handle; // 解引用 VertexHandle 获取 Vertex 对象

        		HalfEdge<Property> &he_obj = he; // 解引用 HEHandle 获取 HalfEdge 对象

				Vec node_location_ralate_to_centroid = v_handle->r - _sys.mesh().get_face_centre(f);

				// he_obj? 
				//v.circulator();f.circulator()
				Vec Force = _force_compute.compute(v, he_obj);//he);//v->data().force;
				stress_yy += -1.0/area*node_location_ralate_to_centroid.y*Force.y;
			}

			return stress_yy;
		}


		double get_face_stress_xy(const Face<Property> &f)
		{
			double stress_xy = 0.0; 
			double area = _sys.mesh().area(f);
			// for (auto he : f.circulator())
			// {
			// VertexHandle<Property> v = he.from();
			// Vec node_location_ralate_to_centroid = v->r - this->get_face_centre(f);
			// Vec Force = v->data().force;
			// stress_xy += -1.0/area*node_location_ralate_to_centroid.x*Force.y;
			// }

			for (auto he : f.circulator())
			{
				VertexHandle<Property> v_handle = he.from();
        		Vertex<Property> &v = *v_handle; // 解引用 VertexHandle 获取 Vertex 对象

        		HalfEdge<Property> &he_obj = he; // 解引用 HEHandle 获取 HalfEdge 对象

				Vec node_location_ralate_to_centroid = v_handle->r - _sys.mesh().get_face_centre(f);

				// he_obj? 
				//v.circulator();f.circulator()
				Vec Force = _force_compute.compute(v, he_obj);//he);//v->data().force;
				stress_xy += -1.0/area*node_location_ralate_to_centroid.x*Force.y;
			}


			return stress_xy;
		}

		void get_face_shape_tensor(const Face<Property> &f, double &Qxx, double &Qxy, double &Qyy)
		{
			Qxx = 0.0;
			Qxy = 0.0;
			Qyy = 0.0;

			if (f.outer || f.erased)
				return;

			double P = _sys.mesh().perim(f);
			if (P < 1e-12)
				return;

			for (auto he : f.circulator())
			{
				Vec e = he.to()->r - he.from()->r;
				double lk = e.len();
				if (lk < 1e-12)
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

		double get_face_shape_qxx(const Face<Property> &f)
		{
			double Qxx, Qxy, Qyy;
			get_face_shape_tensor(f, Qxx, Qxy, Qyy);
			return Qxx;
		}

		double get_face_shape_qxy(const Face<Property> &f)
		{
			double Qxx, Qxy, Qyy;
			get_face_shape_tensor(f, Qxx, Qxy, Qyy);
			return Qxy;
		}

		double get_face_shape_anisotropy(const Face<Property> &f)
		{
			double Qxx, Qxy, Qyy;
			get_face_shape_tensor(f, Qxx, Qxy, Qyy);
			return std::sqrt(2.0 * (Qxx * Qxx + 2.0 * Qxy * Qxy + Qyy * Qyy));
		}

		double get_face_shape_theta(const Face<Property> &f)
		{
			double Qxx, Qxy, Qyy;
			get_face_shape_tensor(f, Qxx, Qxy, Qyy);
			if (std::fabs(Qxx) < 1e-12 && std::fabs(Qxy) < 1e-12 && std::fabs(Qyy) < 1e-12)
				return 0.0;
			return 0.5 * std::atan2(2.0 * Qxy, Qxx - Qyy);
		}

		void get_face_gyration_long_axis(const Face<Property> &f, double &g1, double &g2, double &aspect, double &theta, double &vx, double &vy)
		{
			g1 = 0.0;
			g2 = 0.0;
			aspect = 0.0;
			theta = 0.0;
			vx = 1.0;
			vy = 0.0;

			if (f.outer || f.erased)
				return;

			Vec rc = _sys.mesh().get_face_centroid(f);
			double Gxx = 0.0, Gxy = 0.0, Gyy = 0.0;
			int n = 0;

			for (auto he : f.circulator())
			{
				Vec dr = he.from()->r - rc;
				Gxx += dr.x * dr.x;
				Gxy += dr.x * dr.y;
				Gyy += dr.y * dr.y;
				n++;
			}

			if (n == 0)
				return;

			Gxx /= static_cast<double>(n);
			Gxy /= static_cast<double>(n);
			Gyy /= static_cast<double>(n);

			double trace = Gxx + Gyy;
			double diff = Gxx - Gyy;
			double discr = std::sqrt(diff * diff + 4.0 * Gxy * Gxy);
			g1 = 0.5 * (trace + discr);
			g2 = 0.5 * (trace - discr);
			aspect = (g2 > 1e-12) ? g1 / g2 : 0.0;

			if (discr > 1e-12)
				theta = 0.5 * std::atan2(2.0 * Gxy, diff);

			vx = std::cos(theta);
			vy = std::sin(theta);
		}

		void get_face_hexatic_psi6(const Face<Property> &f, double &psi6_re, double &psi6_im)
		{
			psi6_re = 0.0;
			psi6_im = 0.0;

			if (f.outer || f.erased)
				return;

			Vec rj = _sys.mesh().get_face_centroid(f);
			int Nj = 0;
			for (auto he : f.circulator())
			{
				const Face<Property>& fk = *(he.pair()->face());
				if (fk.outer || fk.erased)
					continue;

				Vec dr = _sys.mesh().get_face_centroid(fk) - rj;
				double theta = std::atan2(dr.y, dr.x);
				psi6_re += std::cos(6.0 * theta);
				psi6_im += std::sin(6.0 * theta);
				Nj += 1;
			}

			if (Nj > 0)
			{
				psi6_re /= static_cast<double>(Nj);
				psi6_im /= static_cast<double>(Nj);
			}
		}

		double get_face_hexatic_psi6_re(const Face<Property> &f)
		{
			double psi6_re, psi6_im;
			get_face_hexatic_psi6(f, psi6_re, psi6_im);
			return psi6_re;
		}

		double get_face_hexatic_psi6_im(const Face<Property> &f)
		{
			double psi6_re, psi6_im;
			get_face_hexatic_psi6(f, psi6_re, psi6_im);
			return psi6_im;
		}

		double get_face_hexatic_psi6_abs(const Face<Property> &f)
		{
			double psi6_re, psi6_im;
			get_face_hexatic_psi6(f, psi6_re, psi6_im);
			return std::sqrt(psi6_re * psi6_re + psi6_im * psi6_im);
		}

		//double get_face_centre


	private:
		System &_sys;
		ForceCompute &_force_compute;
		double _sfc; // Scaling factor for junction output
	};

	void to_json(json &, const HalfEdge<Property> &);
	void to_json(json &, const Edge<Property> &);
	void to_json(json &, const Vertex<Property> &);
	void to_json(json &, const Face<Property> &);
	void to_json(json &, const Box &);

	void from_json(const json &, HalfEdge<Property> &);
	void from_json(const json &, Edge<Property> &);
	void from_json(const json &, Vertex<Property> &);
	void from_json(const json &, Face<Property> &);

	vector<string> split(const std::string &, char);

	void export_Dump(py::module &);


}

#endif
