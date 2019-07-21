

#ifndef NORMAL_CYCLE_H
#define NORMAL_CYCLE_H


//#include <CGAL/circulator.h>
#include <CGAL/HalfedgeDS_decorator.h>
//#include "Tools_Polyhedron_decorator.h"
#include <CGAL/HalfedgeDS_face_base.h> 
#include<set>
#include<map>
#include<list>
#include<vector>
#include<stack>
#include <CGAL/Origin.h>
#include <CGAL/Kernel/global_functions_3.h>

#include"extract_Vpropres.h"


//SURFLAB_BEGIN_NAMESPACE
const static  double m_pi = 3.14159265359;


template <class _Poly>
class Normal_cycle
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
	
	
  
public:

	
void principal_curvature(Polyhedron &pMesh,bool IsGeod,double radius)
{
	pMesh.MinNrmMinCurvature (100000);
	pMesh.MaxNrmMinCurvature (-100000);

	pMesh.MinNrmMaxCurvature (100000);
	pMesh.MaxNrmMaxCurvature (-100000);


  Vertex_iterator pVertex = NULL;
  for(pVertex = pMesh.vertices_begin();
      pVertex != pMesh.vertices_end();
      pVertex++)
  {
    double ppMatrix_sum[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	double eigenvalues[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	double eigenvectors[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	
    if(IsGeod==true)//voisinage geodesique
		geodes_principal_curvature_per_vert((&(*pVertex)),ppMatrix_sum,radius);
	else//voisinage topologique 1-ring
		principal_curvature_per_vert(*pVertex,ppMatrix_sum);

	//valeurs propres
			double **CovMat=(double**)malloc((4)*sizeof(double*));
			double **VectPro=(double**)malloc((4)*sizeof(double*));
			double **Valpro=(double**)malloc((4)*sizeof(double*));
			for(int i=0;i<(4);i++)
			{
				CovMat[i]=(double*)malloc((4)*sizeof(double));
				VectPro[i]=(double*)malloc((4)*sizeof(double));
				Valpro[i]=(double*)malloc((4)*sizeof(double));
			}

			CovMat[1][1]=ppMatrix_sum[0][0];
			CovMat[1][2]=ppMatrix_sum[0][1];
			CovMat[1][3]=ppMatrix_sum[0][2];
			CovMat[2][1]=ppMatrix_sum[1][0];
			CovMat[2][2]=ppMatrix_sum[1][1];
			CovMat[2][3]=ppMatrix_sum[1][2];
			CovMat[3][1]=ppMatrix_sum[2][0];
			CovMat[3][2]=ppMatrix_sum[2][1];
			CovMat[3][3]=ppMatrix_sum[2][2];

			//la matrice n'est elle pas déja diagonale?
			if(ppMatrix_sum[0][1]==0 && ppMatrix_sum[0][2]==0 &
				ppMatrix_sum[1][0]==0 && ppMatrix_sum[1][2]==0 &
				ppMatrix_sum[2][1]==0 && ppMatrix_sum[2][0]==0)
			{
				for(int i=1;i<4;i++)
					for(int j=1;j<4;j++)
					{
						Valpro[i][j]=CovMat[i][j];
						if(i==j)		
						VectPro[i][j]=1;
						else
						VectPro[i][j]=0;


					}



			}
			else
			{
				//recherche des vecteurs et valeurs propres
				if(ValPro(3,CovMat,1e-15,10000.,VectPro,Valpro)==-1)
				{
					pVertex->VKmax(CGAL::NULL_VECTOR);
					pVertex->VKmin(CGAL::NULL_VECTOR);
					return;
				}
			}
				//  Call the Jacovi subroutine 
			for(int u=0;u<4;u++)
				for(int v=0;v<4;v++)
				{
					Valpro[u][v]=fabs(Valpro[u][v]);

				}
			EigSrt(Valpro,VectPro,3);
			Vector_3 VKmax(VectPro[1][2],VectPro[2][2],VectPro[3][2]); 
			Vector_3 VKmin(VectPro[1][1],VectPro[2][1],VectPro[3][1]); 


			eigenvalues[0][0]=Valpro[1][1];
			eigenvalues[0][1]=Valpro[1][2];
			eigenvalues[0][2]=Valpro[1][3];
			eigenvalues[1][0]=Valpro[2][1];
			eigenvalues[1][1]=Valpro[2][2];
			eigenvalues[1][2]=Valpro[2][3];
			eigenvalues[2][0]=Valpro[3][1];
			eigenvalues[2][1]=Valpro[3][2];
			eigenvalues[2][2]=Valpro[3][3];

			eigenvectors[0][0]=VectPro[1][1];
			eigenvectors[0][1]=VectPro[1][2];
			eigenvectors[0][2]=VectPro[1][3];
			eigenvectors[1][0]=VectPro[2][1];
			eigenvectors[1][1]=VectPro[2][2];
			eigenvectors[1][2]=VectPro[2][3];
			eigenvectors[2][0]=VectPro[3][1];
			eigenvectors[2][1]=VectPro[3][2];
			eigenvectors[2][2]=VectPro[3][3];


			
				
			pVertex->VKmax(VKmax);
			pVertex->VKmin(VKmin);

			pVertex->Kmax(Valpro[1][1]);
			pVertex->Kmin(Valpro[2][2]);

			
			for(int i=0;i<(3);i++)
			{
				free(CovMat[i]);
				free(VectPro[i]);
				free(Valpro[i]);
			}
			free(CovMat);
			free(VectPro);
			free(Valpro);


			
			pMesh.MinNrmMinCurvature(std::min(pMesh.MinNrmMinCurvature(),pVertex->Kmin()));
			pMesh.MaxNrmMinCurvature(std::max(pMesh.MaxNrmMinCurvature(),pVertex->Kmin()));

			pMesh.MinNrmMaxCurvature(std::min(pMesh.MinNrmMaxCurvature(),pVertex->Kmax()));
			pMesh.MaxNrmMaxCurvature(std::max(pMesh.MaxNrmMaxCurvature(),pVertex->Kmax()));
   
  }

  

 
}

//**********************************************
// principal curvature for a vertex
//**********************************************
void principal_curvature_per_vert(Vertex pVertex,
                                       double ppMatrix_sum[3][3])
{

	double area=0;
	
  // iterate over all edges
  Halfedge_around_vertex_circulator pHalfedge = pVertex.vertex_begin();
  Halfedge_around_vertex_circulator pHalfedgeStart = pHalfedge;
  CGAL_For_all(pHalfedge,pHalfedgeStart)
  {
	  
    // build edge vector and comput its norm
	  Point p1 = pHalfedge->vertex()->point();
	  Point p2 = pHalfedge->opposite()->vertex()->point();
	  Vector edge = (p1-p2);
		double len_edge = std::sqrt(edge*edge);
		if(len_edge == 0) // avoid divide by zero
		continue;

    // compute (signed) angle between two incident faces, if exists
    Facet_handle pFacet1 = pHalfedge->facet();
    Facet_handle pFacet2 = pHalfedge->opposite()->facet();
    CGAL_assertion(pFacet1 != pFacet2);
    if(pFacet1 == NULL || pFacet2 == NULL)
      continue; // border edge

	area+=AreaFacetTriangle(pHalfedge->facet());

	Vector normal1 = pFacet1->normal();
	Vector normal2 = pFacet2->normal();

    double sine = (CGAL::cross_product(normal1,normal2)*edge)/len_edge;
    double beta = fix_sine(sine);
  
    // compute edge * edge^t * coeff, and add it to current matrix
    double pVector_edge[3] = {edge.x(),edge.y(),edge.z()};
    double ppMatrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    vector_times_transpose_mult(pVector_edge,ppMatrix,beta/len_edge);
    add(ppMatrix,ppMatrix_sum);
  }
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			ppMatrix_sum[i][j]/=area;
}

void geodes_principal_curvature_per_vert(Vertex *pVertex,
                                       double ppMatrix_sum[3][3], double radius)
    {
		
       std::set<Vertex*> vertices ;
        Point O = pVertex->point() ;
        std::stack<Vertex*> S ;
        S.push(pVertex) ;
        vertices.insert(pVertex) ;
		int iter=0;
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
                    
                    if(!h->is_border_edge()) 
					{
						double len_edge = std::sqrt(V*V);
                         // compute (signed) angle between two incident faces, if exists
						Facet_handle pFacet1 = h->facet();
						Facet_handle pFacet2 = h->opposite()->facet();
						CGAL_assertion(pFacet1 != pFacet2);
						if(pFacet1 == NULL || pFacet2 == NULL)
							continue; // border edge
						Vector normal1 = pFacet1->normal();
						Vector normal2 = pFacet2->normal();

						double sine = (CGAL::cross_product(normal1,normal2)*V)/len_edge;
						double beta = fix_sine(sine);

                         // compute edge * edge^t * coeff, and add it to current matrix
						

						double pVector_edge[3] = {V.x(),V.y(),V.z()};
						double ppMatrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
						vector_times_transpose_mult(pVector_edge,ppMatrix,beta/len_edge);
						add(ppMatrix,ppMatrix_sum);
                    }


                    if(!isect) {
						
						Vertex_iterator w=h->opposite()->vertex();
                        if(vertices.find(&(*w)) == vertices.end()) {
                            vertices.insert(&(*w)) ;
                            S.push(&(*w)) ;
                        }
                    }
				}
                
			}
			iter++;
		}

		double area=m_pi*radius*radius;
		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			ppMatrix_sum[i][j]/=area;
    }


//**********************************************
// compute v x v^t
//**********************************************
void vector_times_transpose_mult(double pVector[3],
                                               double ppMatrix[3][3],
                                               double coeff)
{
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      ppMatrix[i][j] = coeff * pVector[i] * pVector[j];
}

//**********************************************
// add two matrices
//**********************************************
void add(double pMatrix[3][3],
                       double pMatrixSum[3][3])
{
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      pMatrixSum[i][j] += pMatrix[i][j];
}

//**********************************************
// fix sine
//**********************************************
double fix_sine(double sine)
{
  if(sine >= 1)
    return m_pi/2;
  else
    if(sine <= -1)
      return -m_pi/2;
    else
      return std::asin(sine);
}

 


 double AreaFacetTriangle(const Facet_handle &f)
	{
		Halfedge_around_facet_circulator pHalfedge = f->facet_begin();
		Point P = pHalfedge->vertex()->point();
		Point Q = pHalfedge->next()->vertex()->point();
		Point R = pHalfedge->next()->next()->vertex()->point();

		Vector PQ=Q-P;
		Vector PR=R-P;
		Vector QR=R-Q;

	
		Vector normal	=	CGAL::cross_product(PQ,QR);
		double area=0.5*sqrt(normal*normal);

		return area;

	}


	bool sphere_clip_vector(Point &O, double r,const Point &P, Vector &V)
    {

        Vector W = P - O ;
        double a = (V*V);
        double b = 2.0 * V * W ;
        double c = (W*W) - r*r ;
        double delta = b*b - 4*a*c ;
        if(delta < 0) {
            // Should not happen, but happens sometimes (numerical precision)
            return true ;
        }
        double t = (- b + ::sqrt(delta)) / (2.0 * a) ;
        if(t < 0.0) {
            // Should not happen, but happens sometimes (numerical precision)
            return true ;
        }
        if(t >= 1.0) {
            // Inside the sphere
            return false ;
        }

        V=V*t;

        return true ;
    }


};
//SURFLAB_END_NAMESPACE
#endif
