

#ifndef ROUGHNESS_H
#define ROUGHNESS_H


//#include <CGAL/circulator.h>
#include <CGAL/HalfedgeDS_decorator.h>
//#include "Tools_Polyhedron_decorator.h"
#include <CGAL/HalfedgeDS_face_base.h>
#include<set>
#include<map>
#include<list>
#include<vector>
#include<stack>

#include <polyhedron_enrichment_define.h>
#include "Polyhedron/polyhedron.h"


#include "normal_cycle.h"





template <class _Poly>
class CRoughness
{
	typedef _Poly                                        Polyhedron;

	typedef typename Polyhedron::Point_3                 Point_3;
	typedef typename Polyhedron::Vector                Vector_3;

  typedef typename Polyhedron::Halfedge_data_structure HDS;
  typedef typename Polyhedron::Vertex                  Vertex;
  typedef typename Polyhedron::Halfedge                Halfedge;
  typedef typename Polyhedron::Facet                   Facet;

  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Vertex_iterator         Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron::Edge_iterator           Edge_iterator;
  typedef typename Polyhedron::Facet_iterator          Facet_iterator;

  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator
                                            Halfedge_around_vertex_circulator;
  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Polyhedron::FT                           FT;
  typedef typename Polyhedron::Vector                       Vector;
  typedef typename Polyhedron::HalfedgeDS             HalfedgeDS;



public:

Polyhedron *m_Polyhedron;

CRoughness(Polyhedron * _Poly1)
{
    m_Polyhedron=_Poly1;
}

void Taubin_smooth(Polyhedron * m_Poly)
	{
		int MatSize=m_Poly->size_of_vertices();

		double* X=new double[MatSize];
		double* Y=new double[MatSize];
		double* Z=new double[MatSize];
		////TAubin smoothing

		int Ind=0;

		for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
				pVertex	!= m_Poly->vertices_end();
				pVertex++)
		{

				int valence=0;
				double SommeWi=0;
				X[Ind]=0;
				Y[Ind]=0;
				Z[Ind]=0;

				Halfedge_around_vertex_circulator	pHalfEdge	=	pVertex->vertex_begin();
				Halfedge_around_vertex_circulator	d	=	pHalfEdge;
				CGAL_For_all(pHalfEdge,d)
				{
					valence++;
					SommeWi+=1;

					X[Ind]+=pHalfEdge->opposite()->vertex()->point().x();
					Y[Ind]+=pHalfEdge->opposite()->vertex()->point().y();
					Z[Ind]+=pHalfEdge->opposite()->vertex()->point().z();

				}

				X[Ind]=X[Ind]/SommeWi-pVertex->point().x();
				Y[Ind]=Y[Ind]/SommeWi-pVertex->point().y();
				Z[Ind]=Z[Ind]/SommeWi-pVertex->point().z();

				X[Ind]=pVertex->point().x()+0.6307*X[Ind];
				Y[Ind]=pVertex->point().y()+0.6307*Y[Ind];
				Z[Ind]=pVertex->point().z()+0.6307*Z[Ind];

				Ind++;

		}
		Ind=0;
		for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
				pVertex	!= m_Poly->vertices_end();
				pVertex++)
		{


			Point p=Point(X[Ind],Y[Ind],Z[Ind]);
			pVertex->point()=p;
			Ind=Ind+1;

		}

		Ind=0;
		for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
				pVertex	!= m_Poly->vertices_end();
				pVertex++)
		{

				int valence=0;
				double SommeWi=0;
				X[Ind]=0;
				Y[Ind]=0;
				Z[Ind]=0;

				Halfedge_around_vertex_circulator	pHalfEdge	=	pVertex->vertex_begin();
				Halfedge_around_vertex_circulator	d	=	pHalfEdge;
				CGAL_For_all(pHalfEdge,d)
				{
					valence++;
					SommeWi+=1;

					X[Ind]+=pHalfEdge->opposite()->vertex()->point().x();
					Y[Ind]+=pHalfEdge->opposite()->vertex()->point().y();
					Z[Ind]+=pHalfEdge->opposite()->vertex()->point().z();

				}

				X[Ind]=X[Ind]/SommeWi-pVertex->point().x();
				Y[Ind]=Y[Ind]/SommeWi-pVertex->point().y();
				Z[Ind]=Z[Ind]/SommeWi-pVertex->point().z();

				X[Ind]=pVertex->point().x()-0.6352*X[Ind];
				Y[Ind]=pVertex->point().y()-0.6352*Y[Ind];
				Z[Ind]=pVertex->point().z()-0.6352*Z[Ind];

				Ind++;

		}



		Ind=0;
		for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
				pVertex	!= m_Poly->vertices_end();
				pVertex++)
		{


			Point p=Point(X[Ind],Y[Ind],Z[Ind]);
			pVertex->point()=p;
			Ind=Ind+1;

		}




	}

#define DEBUG_MOD 100000
void debug_iteration(int i, int n, int depth) {
    if ((i % DEBUG_MOD) == 0) {
        printf("\e[1;1m\e[38;5;087m");
        for (int j = 0; j < depth; j++) {
             printf("│");
        }
        printf("\e[0m %d/%d\n", i, n);
    }
}

void Taubin_smooth_multi_scale(Polyhedron * m_Poly, double radius)
{
	int MatSize=m_Poly->size_of_vertices();

	double* X=new double[MatSize];
	double* Y=new double[MatSize];
	double* Z=new double[MatSize];

	double* XBUF=new double[MatSize];
	double* YBUF=new double[MatSize];
	double* ZBUF=new double[MatSize];

	int Ind=0;
	for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
			pVertex	!= m_Poly->vertices_end();
			pVertex++)
	{


		XBUF[Ind]=pVertex->point().x();
		YBUF[Ind]=pVertex->point().y();
		ZBUF[Ind]=pVertex->point().z();
        Ind++;

	}


	 Ind=0;
     printf("\e[1;1m\e[38;5;087m││┌\e[0m First step\n");
	for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
			pVertex	!= m_Poly->vertices_end();
			pVertex++)
	{
		X[Ind]=0;
		Y[Ind]=0;
		Z[Ind]=0;
        debug_iteration(Ind + 1, MatSize, 3);
		Taubin_smooth_multi_scale_per_vertex(m_Poly,(&(*pVertex)),X[Ind],Y[Ind],Z[Ind],radius,0.6307,XBUF[Ind],YBUF[Ind],ZBUF[Ind]);
		Ind++;
	}
     printf("\e[1;1m\e[38;5;087m││└\e[0m\n");

	Ind=0;
	for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
			pVertex	!= m_Poly->vertices_end();
			pVertex++)
	{


		XBUF[Ind]=X[Ind];
		YBUF[Ind]=Y[Ind];
		ZBUF[Ind]=Z[Ind];
        Ind++;

	}

	Ind=0;

     printf("\e[1;1m\e[38;5;087m││┌\e[0m Second step\n");
	for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
			pVertex	!= m_Poly->vertices_end();
			pVertex++)
	{
        debug_iteration(Ind + 1, MatSize, 3);
		Taubin_smooth_multi_scale_per_vertex(m_Poly,(&(*pVertex)),X[Ind],Y[Ind],Z[Ind],radius,-0.6732,XBUF[Ind],YBUF[Ind],ZBUF[Ind]);//-0.6352
		Ind++;
	}
     printf("\e[1;1m\e[38;5;087m││└\e[0m\n");

	Ind=0;
	for(Vertex_iterator	pVertex	=	m_Poly->vertices_begin();
			pVertex	!= m_Poly->vertices_end();
			pVertex++)
	{


        if (isnan(X[Ind]) || isnan(Y[Ind]) || isnan(Z[Ind])) {
            fprintf(stderr, "\e[38;5;226m┌\e[1;1m Warning:\e[0m Smoothing vertex %d gave the smoothed vertex::\n", Ind + 1);
            fprintf(stderr, "\e[38;5;226m└\e[1;1m (%e, %e, %e) which contains NaNs, ignoring the smoothing of this vertex.\n", X[Ind], Y[Ind], Z[Ind]);
        } else {
            Point p=Point(X[Ind],Y[Ind],Z[Ind]);
            pVertex->point()=p;
        }
        Ind++;

	}


	delete []X;
	delete []Y;
	delete []Z;
	delete []XBUF;
	delete []YBUF;
	delete []ZBUF;



}

bool sphere_clip_vector(Point &O, double r,const Point &P, Vector &V)
{

	Vector W = P - O ;
	double a = (V*V);
	double b = 2.0 * V * W ;
	double c = (W*W) - r*r ;
	double delta = b*b - 4*a*c ;



	if( a==0)
		return true ;

	if(delta < 0) {
		// Should not happen, but happens sometimes (numerical precision)

		return true ;
	}
	double t = (- b + std::sqrt(delta)) / (2.0 * a) ;
	if(t < 0.0) {

		// Should not happen, but happens sometimes (numerical precision)
		return true ;
	}
	if(t >= 1.0) {
		// Inside the sphere
		return false ;
	}

	if(t==0)
	{

		t=0.01;
	}

	V=V*t;

	return true ;
}


void Taubin_smooth_multi_scale_per_vertex(Polyhedron *m_Poly,Vertex* pVertex, double &x, double &y, double &z,double radius, double coef,double &xb, double &yb, double &zb)
{

	std::set<Vertex*> vertices ;
	Point O = pVertex->point() ;
	std::stack<Vertex*> S ;
	S.push(pVertex) ;
	vertices.insert(pVertex) ;
	int iter=0;

	x=y=z=0.;
    double SommeWi; double wi;
	SommeWi=wi=0.;

	while(!S.empty())
	{
		Vertex* v = S.top() ;
		S.pop() ;
		Point P = v->point() ;
		Halfedge_around_vertex_circulator h = v->vertex_begin();
		Halfedge_around_vertex_circulator pHalfedgeStart = h;
		CGAL_For_all(h,pHalfedgeStart)
		{
			Point p1 = h->vertex()->point();
			Point p2 = h->opposite()->vertex()->point();
			Vector V = (p2-p1);
			if(v==pVertex || V * (P - O) > 0.0)
			{
				double len_old = std::sqrt(V*V);


				bool isect = sphere_clip_vector(O, radius, P, V) ;



				if(!isect)
				{

					Vertex_iterator w=h->opposite()->vertex();
					if(vertices.find(&(*w)) == vertices.end())
					{
						vertices.insert(&(*w)) ;
						S.push(&(*w)) ;
					}
				}
				else
				{
					p2=Point(V.x()+p1.x(),V.y()+p1.y(),V.z()+p1.z());
					double Wi=1;///(sqrt(V*V));///Wi=1
                    SommeWi+=Wi;

					x+=Wi*p2.x();
					y+=Wi*p2.y();
					z+=Wi*p2.z();
				}
			}

		}
		iter++;
	}


	if(SommeWi==0)
	{
		x=xb;
		y=yb;
		z=zb;

	}
	else
	{


		x=x/SommeWi-xb;
		y=y/SommeWi-yb;
		z=z/SommeWi-zb;


		x=xb+coef*x;
		y=yb+coef*y;
		z=zb+coef*z;
	}




}

void Processroughness_per_vertex_curve(Polyhedron * PolyUsed,Vertex* pVertex,double radius,std::vector<double> &TabDistance,std::vector<Point> &TabPoint,double & moyenneRet,bool IsGauss = false)
	{
		std::set<Vertex*> vertices ;
        Point O = pVertex->point() ;
        std::stack<Vertex*> S ;
        S.push(pVertex) ;
        vertices.insert(pVertex) ;


		int NbSommetInSphere=0;
		double SommeDistance=0;


        while(!S.empty())
		{
			Vertex* v = S.top() ;
            S.pop() ;
            Point P = v->point() ;
            Halfedge_around_vertex_circulator h = v->vertex_begin();
			Halfedge_around_vertex_circulator pHalfedgeStart = h;
			CGAL_For_all(h,pHalfedgeStart)
			{
                Point p1 = h->vertex()->point();
				Point p2 = h->opposite()->vertex()->point();
				Vector V = (p2-p1);
                if(v==pVertex || V * (P - O) > 0.0)
				{
					double len_old = std::sqrt(V*V);
					bool isect = sphere_clip_vector(O, radius, P, V) ;
					double len_edge = std::sqrt(V*V);

					NbSommetInSphere++;
					///ici on prend en compte la distance map des sommets

					double DistancePondere;
					Point PPondere;
					if(len_old!=0)
					{
						DistancePondere=(1-len_edge/len_old)*h->vertex()->Kmax()+len_edge/len_old*h->opposite()->vertex()->Kmax();
						PPondere=p1+V;
					}
					else
					{
						DistancePondere=h->opposite()->vertex()->Kmax();
						PPondere=p2;
					}

					TabDistance.push_back(DistancePondere);
					TabPoint.push_back(PPondere);

					SommeDistance+=DistancePondere;

					if(!isect)
					{

						Vertex_iterator w=h->opposite()->vertex();
                        if(vertices.find(&(*w)) == vertices.end())
						{
                            vertices.insert(&(*w)) ;
                            S.push(&(*w)) ;
                        }
                    }

				}

			}

		}

		double moyenne=0;
		double variance=0;

		if(IsGauss==false)
		{
			moyenne=SommeDistance/(double)NbSommetInSphere;




			variance=0;
			if(NbSommetInSphere!=0)
			{
				for(int i=0;i<NbSommetInSphere;i++)
					variance+=pow(TabDistance[i]-moyenne,2);

				variance=variance/(double)NbSommetInSphere;
				variance=sqrt(variance);
			}
		}
		else
		{//variance = 0.008
			SommeDistance=0;
			double SommeWi=0;
			for(int i=0;i<NbSommetInSphere;i++)
			{
				Vector DistancePt=TabPoint[i]-pVertex->point();
				double distPt=sqrt(DistancePt*DistancePt);
				double wi=1/0.008/sqrt(2*3.141592)*exp(-(distPt*distPt)/2/0.008/0.008);
				SommeWi+=wi;
				SommeDistance+=TabDistance[i]*wi;
			}
			moyenne=SommeDistance/(double)SommeWi;


            for(int i=0;i<NbSommetInSphere;i++)
			{
				Vector DistancePt=TabPoint[i]-pVertex->point();
				double distPt=sqrt(DistancePt*DistancePt);
				double wi=1/0.008/sqrt(2*3.141592)*exp(-(distPt*distPt)/2/0.008/0.008);

				variance+=wi*pow(TabDistance[i]-moyenne,2);
			}

				variance=variance/SommeWi;
				variance=sqrt(variance);

		}



			pVertex->CourbureMoyenne=moyenne;
			pVertex->CourbureVariance=variance;


	}


	void Processroughness_curve(Polyhedron *PolyUsed,double radius)
	{

		double somme_roughness=0;
		int NbVert=0;

		for(Vertex_iterator	pVertex	=	PolyUsed->vertices_begin();
					pVertex	!= PolyUsed->vertices_end();
					pVertex++)
		{
				std::vector<double> TabDistance;
				std::vector<Point> TabPoint;
				double moyenne;

                debug_iteration(NbVert + 1, PolyUsed->size_of_vertices(), 1);
				Processroughness_per_vertex_curve(PolyUsed,(&(*pVertex)),radius,TabDistance,TabPoint,moyenne,true);
				NbVert++;

		}

	}

void compute_Roughness(double radius, double SmoothRadius, double CurvatureRadius, double scaling)
{

	Polyhedron SmoothPoly=*m_Polyhedron;

	m_Polyhedron->MinNrmRoughCurvature(10000);
	m_Polyhedron->MaxNrmRoughCurvature(-10000);



    printf("\e[1;1m\e[38;5;087m┌ Info:\e[0m Adaptive smoothing with 2-step Taubin filter with radius %lf\n", fabs(SmoothRadius));
    if (SmoothRadius < 0) {
    	bool ok = SmoothPoly.load_mesh("smooth.off");
        if (!ok) {
            fprintf(stderr, "Failed to read mesh from smooth.off");
            exit(EXIT_FAILURE);
        }
        printf("\e[1;1m\e[38;5;087m└\e[0m read from smooth.off\n");
    } else {
        int n = 5;
    	for(int i=0;i<n;i++) {
            printf("\e[1;1m\e[38;5;087m│┌\e[0m %d/%d\n", i + 1, n);
    		Taubin_smooth_multi_scale(&SmoothPoly,SmoothRadius);
        }
    	SmoothPoly.write_off("smooth.off", false, false);
        printf("\e[1;1m\e[38;5;087m└\e[0m written to smooth.off\n");
    }

/*
    printf("Laplacian smoothing with 2-step Taubin filter\n");
	for(int i=0;i<10;i++)
	{
        printf("%d\n", i);
		Taubin_smooth(&SmoothPoly);
		Taubin_smooth(&SmoothPoly);
		Taubin_smooth(&SmoothPoly);
		Taubin_smooth(&SmoothPoly);
		Taubin_smooth(&SmoothPoly);
		Taubin_smooth(&SmoothPoly);
		Taubin_smooth(&SmoothPoly);
	}*/


	SmoothPoly.compute_normals();

	Normal_cycle<Polyhedron> estimateur;
    printf("\e[1;1m\e[38;5;087m┌ Info:\e[0m Compute maximum curvature kmax of each vertex of the smooth mesh\n");
	estimateur.principal_curvature(SmoothPoly,true,CurvatureRadius,scaling);
    FILE *smooth_kmax = fopen("smooth_kmax.txt", "w");
    for (Vertex_iterator pVertex = SmoothPoly.vertices_begin(); pVertex != SmoothPoly.vertices_end(); pVertex++)
    {
        fprintf(smooth_kmax, "%lf\n", pVertex->Kmax());
    }
    fclose(smooth_kmax);
    printf("\e[1;1m\e[38;5;087m└\e[0m written to smooth_kmax.txt\n");

    printf("\e[1;1m\e[38;5;087m┌ Info:\e[0m Compute maximum curvature kmax of each vertex of the original mesh\n");
	estimateur.principal_curvature(*m_Polyhedron,true,CurvatureRadius,scaling);
    FILE *original_kmax = fopen("original_kmax.txt", "w");
    for (Vertex_iterator pVertex = m_Polyhedron->vertices_begin(); pVertex != m_Polyhedron->vertices_end(); pVertex++)
    {
        fprintf(original_kmax, "%lf\n", pVertex->Kmax());
    }
    fclose(original_kmax);
    printf("\e[1;1m\e[38;5;087m└\e[0m written to original_kmax.txt\n");

    printf("\e[1;1m\e[38;5;087m┌ Info:\e[0m Compute average curvature kav of each vertex of the smooth mesh\n");
	Processroughness_curve(&SmoothPoly,radius);
    FILE *smooth_kav = fopen("smooth_kav.txt", "w");
    for (Vertex_iterator pVertex = SmoothPoly.vertices_begin(); pVertex != SmoothPoly.vertices_end(); pVertex++)
    {
        fprintf(smooth_kav, "%lf\n", pVertex->CourbureMoyenne);
    }
    fclose(smooth_kav);
    printf("\e[1;1m\e[38;5;087m└\e[0m written to smooth_kav.txt\n");
    printf("\e[1;1m\e[38;5;087m┌ Info:\e[0m Compute average curvature kav of each vertex of the original mesh\n");
	Processroughness_curve(m_Polyhedron,radius);
    FILE *original_kav = fopen("original_kav.txt", "w");
    for (Vertex_iterator pVertex = m_Polyhedron->vertices_begin(); pVertex != m_Polyhedron->vertices_end(); pVertex++)
    {
        fprintf(original_kav, "%lf\n", pVertex->CourbureMoyenne);
    }
    fclose(original_kav);

    printf("\e[1;1m\e[38;5;087m┌ Info:\e[0m Compute roughness for each vertex\n");
	Vertex_iterator	pVertex2=SmoothPoly.vertices_begin();
	for(Vertex_iterator	pVertex	=	m_Polyhedron->vertices_begin();
				pVertex	!= m_Polyhedron->vertices_end();
				pVertex++)
	{
		double offset=0;
		if(pVertex->CourbureMoyenne>pVertex2->CourbureMoyenne)
			offset=fabs(pVertex->CourbureMoyenne-pVertex2->CourbureMoyenne);

		pVertex->Roughness(offset);
        if(offset<m_Polyhedron->MinNrmRoughCurvature())
					m_Polyhedron->MinNrmRoughCurvature(offset);

		if(offset>m_Polyhedron->MaxNrmRoughCurvature())
			m_Polyhedron->MaxNrmRoughCurvature(offset);

		pVertex2++;
	}
    FILE *roughness = fopen("roughness.txt", "w");
    for (Vertex_iterator pVertex = m_Polyhedron->vertices_begin(); pVertex != m_Polyhedron->vertices_end(); pVertex++)
    {
        fprintf(roughness, "%lf\n", pVertex->Roughness());
    }
    fclose(roughness);
    printf("\e[1;1m\e[38;5;087m└\e[0m written to roughness.txt\n");
}



};
//SURFLAB_END_NAMESPACE
#endif
