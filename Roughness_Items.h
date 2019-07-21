#ifndef Roughness_ITEMS_H
#define Roughness_ITEMS_H

/*!
 * \file Roughness_Items.h
 * \brief Enrichement of the Facets, Halfedges and vertices of the polyhedra
 * \author \author Delmotte Arnaud, NAIST, based on code from Guillaume Lavoue,  INSA-Lyon, LIRIS UMR 5205, M2DISCO.
    \date	2018
 */
 
#include <polyhedron_enrichment_define.h>
#ifdef BUILD_component_Roughness

#include "Polyhedron/polyhedron_shared_items.h"

/**
 \class	Roughness_Facet

 \brief	Enriches the Facets of a Polyhedra 

 */
template <class Refs, class T, class P, class Norm, class Plane>
class Roughness_Facet : virtual public MEPP_Common_Facet<Refs, T, Norm>
{
	public:

        Roughness_Facet() {}
};


/*!
 * \class Roughness_Halfedge
 * \brief Enriches the Halfedges of a Polyhedra
 */
template <class Refs, class Tprev, class Tvertex, class Tface, class Norm>
class Roughness_Halfedge : virtual public MEPP_Common_Halfedge<Refs,Tprev,Tvertex,Tface,Norm>
{
	public:

        Roughness_Halfedge() {}
};


/*!
 * \class Roughness_Vertex
 * \brief Enriches the Vertices of a Polyhedra
 */
template <class Refs, class T, class P, class Norm>
class Roughness_Vertex : virtual public MEPP_Common_Vertex<Refs, T, P, Norm>
{
	public:

        Roughness_Vertex() {}


        void Roughness(double val)
        {
            _Roughness = val;
        }

        double Roughness()
        {
            return _Roughness;
        }

        void VKmin(Norm val)
        {
            _VKmin = val;
        }

        void VKmax(Norm val)
        {
            _VKmax = val;
        }

        Norm VKmin()
        {
            return _VKmin;
        }

        Norm VKmax()
        {
            return _VKmax;
        }

        void Kmin(double val)
        {
            _Kmin = val;
        }

        void Kmax(double val)
        {
            _Kmax = val;
        }

        double Kmin()
        {
            return _Kmin;
        }

        double Kmax()
        {
            return _Kmax;
        }

        double CourbureMoyenne;
        double CourbureVariance;
        double _Roughness;
        Norm _VKmin;
        Norm _VKmax;
        double _Kmin;
        double _Kmax;

};

/*!
 * \class Roughness_Polyhedron
 * \brief Enriches the Polyhedron
 */
template <class kernel, class items>
class Roughness_Polyhedron : virtual public MEPP_Common_Polyhedron<kernel,items>
{
	public:

        Roughness_Polyhedron() {  }

        void MinNrmMinCurvature(double val)
        {
            _MinNrmMinCurvature = val;
        }

        void MinNrmMaxCurvature(double val)
        {
            _MinNrmMaxCurvature = val;
        }

        void MaxNrmMinCurvature(double val)
        {
            _MaxNrmMinCurvature = val;
        }

        void MaxNrmMaxCurvature(double val)
        {
            _MaxNrmMaxCurvature = val;
        }

        double MinNrmMinCurvature()
        {
            return _MinNrmMinCurvature;
        }

        double MinNrmMaxCurvature()
        {
            return _MinNrmMaxCurvature;
        }

        double MaxNrmMinCurvature()
        {
            return _MaxNrmMinCurvature;
        }

        double MaxNrmMaxCurvature()
        {
            return _MaxNrmMaxCurvature;
        }

        void MinNrmRoughCurvature(double val)
        {
            _MinNrmRoughCurvature = val;
        }

        void MaxNrmRoughCurvature(double val)
        {
            _MaxNrmRoughCurvature = val;
        }

        double MinNrmRoughCurvature()
        {
            return _MinNrmRoughCurvature;
        }

        double MaxNrmRoughCurvature()
        {
            return _MaxNrmRoughCurvature;
        }

        double _MinNrmRoughCurvature, _MaxNrmRoughCurvature;
        double _MinNrmMinCurvature, _MinNrmMaxCurvature, _MaxNrmMinCurvature, _MaxNrmMaxCurvature;
};

#endif

#endif // Curvature_ITEMS_H
