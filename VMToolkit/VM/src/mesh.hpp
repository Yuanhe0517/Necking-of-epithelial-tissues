
/*!
 * \file mesh.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Jun-2017
 * \brief Mesh class
 */

#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <algorithm>
#include <vector>
#include <iterator>
#include <iostream>
#include <exception>
#include <memory>
#include <list>
#include <math.h>

#include <pybind11/stl.h>

#include "types.hpp"
#include "box.hpp"

using std::find_if;
using std::list;
using std::next;
using std::prev;
using std::runtime_error;
using std::shared_ptr;
using std::vector;

namespace py = pybind11;

namespace VMTutorial
{

	template <typename Property>
	class Mesh
	{

	public:
		Mesh() : _n_faces(0), _box(nullptr) {}

		//! Mesh setup functions
		void add_vertex(const Vertex<Property> &v)
		{
			_vertices.push_back(v);
			if (v.erased)
				_erased_vertices.push_back(v.id);
		}
		void add_edge(const Edge<Property> &e) { _edges.push_back(e); }
		void add_halfedge(const HalfEdge<Property> &he) { _halfedges.push_back(he); }
		void set_box(const shared_ptr<Box> &box) { _box = box; }
		void add_face(const vector<int> &, bool = false);


		void wipe()
		{
			_halfedges.clear();
			_vertices.clear();
			_edges.clear();
			_faces.clear();
			_erased_edges.clear();
			_erased_halfedges.clear();
			_erased_vertices.clear();
			_n_faces = 0;
			_box = nullptr;
		}

		// accessor functions
		HEHandle<Property> get_mesh_he(int);
		VertexHandle<Property> get_mesh_vertex(int);
		EdgeHandle<Property> get_mesh_edge(int);
		FaceHandle<Property> get_mesh_face(int);

		Vertex<Property> &get_vertex(int);
		HalfEdge<Property> &get_halfedge(int);
		Edge<Property> &get_edge(int);
		Face<Property> &get_face(int);

		vector<HalfEdge<Property>> &halfedges() { return _halfedges; }
		vector<Vertex<Property>> &vertices() { return _vertices; }
		vector<Edge<Property>> &edges() { return _edges; }
		vector<Face<Property>> &faces() { return _faces; }

		int num_vert() { return _vertices.size(); }
		int num_faces() { return _faces.size(); }

		// mesh manipulation functions
		
		bool T1(Edge<Property> &, double);

		// mesh info functions

		double area(const Face<Property> &) const;
		double perim(const Face<Property> &) const;
		double len(const Edge<Property> &);
		int coordination(const Vertex<Property> &);
		int face_sides(const Face<Property> &);
		bool is_boundary_face(const Face<Property> &);

		void tidyup();	  // get the mesh in order (set boundary edges, outer faces, etc. )
		Vec get_centre(); // compute geometric centre of the mesh
		Vec get_face_centre(const Face<Property> &);
		Vec get_face_centroid(const Face<Property> &);
		Vec get_face_direction(const Face<Property> &);
		double get_face_stress_xx(const Face<Property> &);
		double get_face_stress_yy(const Face<Property> &);
		double get_face_stress_xy(const Face<Property> &);

		const shared_ptr<Box> &box() const { return _box; }

	private:
		vector<HalfEdge<Property>> _halfedges;
		vector<Vertex<Property>> _vertices;
		vector<Edge<Property>> _edges;
		vector<Face<Property>> _faces;

		list<int> _erased_vertices;
		list<int> _erased_edges;
		list<int> _erased_halfedges;
		list<int> _erased_faces;

		int _n_faces;		  // number of faces;
		shared_ptr<Box> _box; // simulation box
	};

	template <typename Property>
	Vertex<Property> &Mesh<Property>::get_vertex(int i)
	{
		if ((i < 0) || (i > _vertices.size()))
			throw runtime_error("Vertex index out of bounds.");
		else
			return _vertices[i];
	}

	template <typename Property>
	HalfEdge<Property> &Mesh<Property>::get_halfedge(int i)
	{
		if ((i < 0) || (i > _halfedges.size()))
			throw runtime_error("Junction index out of bounds.");
		else
			return _halfedges[i];
	}

	template <typename Property>
	Edge<Property> &Mesh<Property>::get_edge(int i)
	{
		if ((i < 0) || (i > _edges.size()))
			throw runtime_error("Junction index out of bounds.");
		else
			return _edges[i];
	}

	template <typename Property>
	Face<Property> &Mesh<Property>::get_face(int i)
	{
		if ((i < 0) || (i > _faces.size()))
			throw runtime_error("Face index out of bounds.");
		else
			return *(find_if(_faces.begin(), _faces.end(), [i](const Face<Property> &f) -> bool
							 { return (f.id == i); }));
	}

	// accessor functions
	template <typename Property>
	HEHandle<Property> Mesh<Property>::get_mesh_he(int i)
	{
		if ((i >= 0) && (i < _halfedges.size()))
		{
			return std::next(_halfedges.begin(), i);
		}
		else
		{
			throw runtime_error("HalfEdge index out of bounds.");
		}
	}

	template <typename Property>
	VertexHandle<Property> Mesh<Property>::get_mesh_vertex(int i)
	{
		if ((i >= 0) && (i < _vertices.size()))
		{
			return std::next(_vertices.begin(), i);
		}
		else
		{
			throw runtime_error("Vertex index out of bounds.");
		}
	}

	template <typename Property>
	EdgeHandle<Property> Mesh<Property>::get_mesh_edge(int i)
	{
		if ((i >= 0) && (i < _edges.size()))
		{
			return std::next(_edges.begin(), i);
		}
		else
		{
			throw runtime_error("Edge index out of bounds.");
		}
	}

	template <typename Property>
	FaceHandle<Property> Mesh<Property>::get_mesh_face(int i)
	{
		if ((i >= 0) && (i < _faces.size()))
		{
			return std::next(_faces.begin(), i);
		}
		else
		{
			throw runtime_error("Face index out of bounds.");
		}
	}
	// end of accessor functions

	// mesh setup functions
	template <typename Property>
	void Mesh<Property>::add_face(const vector<int> &vert_ids, bool erased)
	{
		_faces.push_back(Face<Property>(_n_faces++, erased, *this));
		if (!erased)
		{
			FaceHandle<Property> fh = prev(_faces.end());
			HEHandle<Property> he;
			int prev_he;
			int first_he;
			for (int i = 0; i < vert_ids.size(); i++)
			{
				int v1_id = vert_ids[i], v2_id = vert_ids[(i == (vert_ids.size() - 1)) ? 0 : i + 1];
				VertexHandle<Property> vh_from = this->get_mesh_vertex(v1_id);
				VertexHandle<Property> vh_to = this->get_mesh_vertex(v2_id);
				_halfedges.push_back(HalfEdge<Property>(*this));
				_halfedges.back().set_idx(_halfedges.size() - 1);
				he = this->get_mesh_he(_halfedges.size() - 1);
				if (i == 0)
					first_he = he->idx();
				he->_from = vh_from->id;
				he->_to = vh_to->id;
				vh_from->_he = he->idx();
				EdgeHandle<Property> eh = find_if(_edges.begin(), _edges.end(), [v1_id, v2_id](const Edge<Property> &e) -> bool
												  { return (v1_id == e.i && v2_id == e.j) || (v1_id == e.j && v2_id == e.i); });

				if (eh == _edges.end())
				{
					_edges.push_back(Edge<Property>(v1_id, v2_id, *this));
					_edges.back().set_idx(_edges.size() - 1);
					eh = prev(_edges.end());
					eh->_he = he->idx();
				}
				else
				{
					eh->he()->_pair = he->idx();
					he->_pair = eh->he()->idx();
					eh->boundary = false;
				}
				he->_edge = eh->idx();
				he->_face = fh->id;
				if (i > 0)
				{
					he->_prev = _halfedges[prev_he].idx();
					_halfedges[prev_he]._next = he->idx();
				}
				prev_he = he->idx();
			}
			_halfedges[first_he]._prev = he->idx();
			he->_next = _halfedges[first_he].idx();
			fh->_he = _halfedges[first_he].idx(); // Make sure that the first half edge is the face half-edge. This makes life easier when reading in "per edge" data.
			fh->nsides = this->face_sides(*fh);
		}
		else
			_erased_faces.push_back(_n_faces - 1);
	}

	// Mesh manipulation functions
	//! Implements actual T1
	template <typename Property>
	bool Mesh<Property>::T1(Edge<Property> &e, double edge_len)
	{
		if (e.boundary)
			return false;
		//cout<<e._he<<"_he"<<endl;
	
		HEHandle<Property> he = e.he();		 // Half-edge belonging to the end
		HEHandle<Property> hep = he->pair(); // Half-edge pair to he


		VertexHandle<Property> v1 = he->from(); // We define v1 and the vertex he points from
		VertexHandle<Property> v2 = he->to();	// We define v2 and the vertex he points to
		
		//cout<<"1 v1 he"<<v1->he()->from()->r<<v1->he()->to()->r<<endl;
		//cout<<"1 v2 hep"<<v2->he()->from()->r<<v2->he()->to()->r<<endl;

		Vec AB = v2->r - v1->r;
		Vec l = 0.5 * AB; // Vector l points from v1 towards the geometric centre of the v1-v2 line

		Vec rc = v1->r + l;

		// Vec rot_l = Vec(-l.y*1.5*0.1/AB.len(),l.x*1.5*0.1/AB.len());

		// v1->r = rc - 0.5 * rot_l;
		// v2->r = rc + 0.5 * rot_l;

		Vec rot_l = Vec(-l.y, l.x).unit();

		v1->r = rc - 0.5 * edge_len * rot_l;
		v2->r = rc + 0.5 * edge_len * rot_l;

		// cout<<"2 v1 he"<<v1->he()->from()->r<<v1->he()->to()->r<<endl;
		// cout<<"2 v2 hep"<<v2->he()->from()->r<<v2->he()->to()->r<<endl;

		v1->set_he(he->_idx);
		v2->set_he(hep->_idx);

		// cout<<"he"<<he->from()->r<<he->to()->r<<endl;
		// cout<<"hep"<<hep->from()->r<<hep->to()->r<<endl;

		// cout<<"3 v1 he"<<v1->he()->from()->r<<v1->he()->to()->r<<endl;
		// cout<<"3 v2 hep"<<v2->he()->from()->r<<v2->he()->to()->r<<endl;

		//cout<<he.idx()<<endl;

		HEHandle<Property> he1 = he->prev();
		HEHandle<Property> he2 = he->next();
		HEHandle<Property> he3 = hep->prev();
		HEHandle<Property> he4 = hep->next();

		// cout<<"he1"<<he1->from()->r<<he1->to()->r<<endl;
		// cout<<"he2"<<he2->from()->r<<he2->to()->r<<endl;
		// cout<<"he3"<<he3->from()->r<<he3->to()->r<<endl;
		// cout<<"he4"<<he4->from()->r<<he4->to()->r<<endl;

		// Change he relation
		he1->set_next(he2->_idx);
		he2->set_prev(he1->_idx);
		he3->set_next(he4->_idx);
		he4->set_prev(he3->_idx);
		// cout<<"he1_next"<<he1->next()->from()->r<<he1->next()->to()->r<<endl;

		// he1->next() = he2;
		// he2->prev() = he1;
		// he3->next() = he4;
		// he4->prev() = he3;

		he->set_next(he1->pair()->_idx);
		he->set_prev(he4->pair()->_idx);
		hep->set_next(he3->pair()->_idx);
		hep->set_prev(he2->pair()->_idx);


		// he->next() = he1->pair();
		// he->prev() = he4->pair();
		// hep->next() = he3->pair();
		// hep->prev() = he2->pair();

		he1->pair()->set_prev(he->_idx);
		he2->pair()->set_next(hep->_idx);
		he3->pair()->set_prev(hep->_idx);
		he4->pair()->set_next(he->_idx);

		// he1->pair()->prev() = he;
		// he2->pair()->next() = hep;
		// he3->pair()->prev() = hep;
		// he4->pair()->next() = he;

		// he1->to() = v2;
		// he1->pair()->from() = v2;
		//cout<<"he1_to"<<he1->_to<<endl;
		//cout<<"he1_from"<<he1->_from<<endl;

		he1->set_to(v2->id);
		he1->pair()->set_from(v2->id);

		// cout<<"new_he1"<<he1->from()->r<<he1->to()->r<<endl;
		// cout<<"new_he1_pair"<<he1->pair()->from()->r<<he1->pair()->to()->r<<endl;

		// he3->to() = v1;
		// he3->pair()->from() = v1;

		he3->set_to(v1->id);
		he3->pair()->set_from(v1->id);

		// cout<<"new_he3"<<he3->from()->r<<he3->to()->r<<endl;
		// cout<<"new_he3_pair"<<he3->pair()->from()->r<<he3->pair()->to()->r<<endl;
		// cout<<"he"<<he->from()->r<<he->to()->r<<endl;
		// he2->to() = v1;
		// he1->pair()->from() = v1;

		// he3->to() = v2;
		// he3->pair()->from() = v2;

		he->face()->set_he(he2->_idx);
		hep->face()->set_he(he4->_idx);

		he->set_face(he1->pair()->face()->id);
		hep->set_face(he2->pair()->face()->id);

		he->face()->set_nsides(this->face_sides(*(he->face())));
		hep->face()->set_nsides(this->face_sides(*(hep->face())));
		he2->face()->set_nsides(this->face_sides(*(he2->face())));
		he4->face()->set_nsides(this->face_sides(*(he4->face())));

		//he->face()->set_nsides(this->face_sides(*(he->face())));
		//hep->face()->set_nsides(this->face_sides(*(hep->face())));
		// cout<<"he_face"<<he->face()->nsides<<endl;
		// cout<<"hep_face"<<hep->face()->nsides<<endl;
		// cout<<"he2_face"<<he2->face()->nsides<<endl;
		// cout<<"he4_face"<<he4->face()->nsides<<endl;

		// for (auto i : he->face()->circulator())
		// {
		// 	cout<<i.from()->r<<i.to()->r<<endl;
		// }

		// for (auto i : he2->face()->circulator())
		// {
		// 	cout<<i.from()->r<<i.to()->r<<endl;
		// }
		

		//he->face()->he() = he2;
		//hep->face()->he() = he4;

		// he->face() = he1->pair()->face();
		// hep->face() = he2->pair()->face();

		//he->face()->nsides = this->face_sides(*(he->face()));
		//hep->face()->nsides = this->face_sides(*(hep->face())); 

		return true;
	}

	template <typename Property>
	void Mesh<Property>::tidyup()
	{
		//cout<<"begin"<<endl;
		for (auto& v : _vertices)
		{	//cout<<"process vertices"<<endl;
			if (!v.erased)
				v.coordination = this->coordination(v);
		}

		for (auto& e : _edges)
		{
			//cout<<"process edge"<<endl;
			e.boundary = false;
			if (e.he()->from()->boundary && e.he()->to()->boundary)
				e.boundary = true;
		}
	}

	// Mesh info functions
	template <typename Property>
	double Mesh<Property>::area(const Face<Property> &f) const
	{
		if (f.outer)
			return 0.0;

		Vec r0 = f.he()->from()->r;
		double A = 0.0;
		for (auto he : f.circulator())
		{
			Vec r1 = he.from()->r - r0; // this takes care of the boundary conditions
			Vec r2 = he.to()->r - r0;
			A += r1.x * r2.y - r2.x * r1.y;
		}
		return 0.5 * fabs(A);
	}

	template <typename Property>
	double Mesh<Property>::perim(const Face<Property> &f) const
	{
		if (f.outer)
			return 0.0;

		double P = 0.0;
		for (auto he : f.circulator())
			P += (he.to()->r - he.from()->r).len();

		return P;
	}

	template <typename Property>
	double Mesh<Property>::len(const Edge<Property> &e)
	{
		return (e.he()->to()->r - e.he()->from()->r).len();
	}

	template <typename Property>
	int Mesh<Property>::coordination(const Vertex<Property> &v)
	{
		//cout<<"coordination"<<endl;
		if (v.erased)
			return -1;
		int i = 0;
		for (auto he : v.circulator())
			i++;
		//cout<<i<<endl;
		return i;
	}

	template <typename Property>
	int Mesh<Property>::face_sides(const Face<Property> &f)
	{
		int i = 0;
		for (auto he : f.circulator())
			i++;

		return i;
	}

	template <typename Property>
	bool Mesh<Property>::is_boundary_face(const Face<Property> &f)
	{
		for (auto he : f.circulator())
		{
			if (he.from()->boundary)
				return true;
		}
		return false;
	}

	// Compute geometric centre of the mesh by tracing positions of boundary vertices
	template <typename Property>
	Vec Mesh<Property>::get_centre()
	{
		Vec cm(0.0, 0.0);
		for (auto v : _vertices)
			cm += v.r;
		return static_cast<double>(1.0 / _vertices.size()) * cm;
	}

	// Compute centre of a face
	template <typename Property>
	Vec Mesh<Property>::get_face_centre(const Face<Property> &f)
	{
		if (f.outer)
			return Vec(0.0, 0.0);
		Vec r0 = f.he()->from()->r;
		Vec rc(0.0, 0.0);
		for (auto he : f.circulator())
		{
			Vec dr = he.from()->r - r0;
			rc += Vec(dr.x, dr.y);
		}
		Vec Rc = (1.0 / f.nsides) * rc + r0;
		return Vec(Rc.x, Rc.y, this->_box);
	}

	// Compute position of the face centroid
	template <typename Property>
	Vec Mesh<Property>::get_face_centroid(const Face<Property> &f)
	{
		if (f.outer)
			return Vec(0.0, 0.0);

		Vec r0 = f.he()->from()->r;
		Vec rc(0.0, 0.0);
		for (auto he : f.circulator())
		{
			Vec ri = he.from()->r - r0;
			Vec rj = he.to()->r - r0;
			double fact = ri.x * rj.y - ri.y * rj.x;
			rc.x += (ri.x + rj.x) * fact;
			rc.y += (ri.y + rj.y) * fact;
		}
		Vec Rc = (1.0 / (6 * this->area(f))) * rc + r0;
		return Vec(Rc.x, Rc.y, this->_box);
	}

	template <typename Property>
	double Mesh<Property>::get_face_stress_xx(const Face<Property> &f)
	{
		double stress_xx = 0.0; 
		double area = this->area(f);
		for (auto he : f.circulator())
		{
			VertexHandle<Property> v = he.from();
			Vec node_location_ralate_to_centroid = v->r - this->get_face_centre(f);

			Vec Force = v->data().force;
			stress_xx += -1.0/area*node_location_ralate_to_centroid.x*Force.x;
		}

		return stress_xx;
	}

	template <typename Property>
	double Mesh<Property>::get_face_stress_yy(const Face<Property> &f)
	{
		double stress_yy = 0.0; 
		double area = this->area(f);
		for (auto he : f.circulator())
		{
			VertexHandle<Property> v = he.from();
			Vec node_location_ralate_to_centroid = v->r- this->get_face_centre(f);
			Vec Force = v->data().force;
			stress_yy += -1.0/area*node_location_ralate_to_centroid.y*Force.y;
		}

		return stress_yy;
	}

	template <typename Property>
	double Mesh<Property>::get_face_stress_xy(const Face<Property> &f)
	{
		double stress_xy = 0.0; 
		double area = this->area(f);
		for (auto he : f.circulator())
		{
			VertexHandle<Property> v = he.from();
			Vec node_location_ralate_to_centroid = v->r - this->get_face_centre(f);
			Vec Force = v->data().force;
			stress_xy += -1.0/area*node_location_ralate_to_centroid.x*Force.y;
		}

		return stress_xy;
	}

};

#endif
