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