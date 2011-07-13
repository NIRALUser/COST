#ifndef SphereIkosohedron_txx
#define SphereIkosohedron_txx

#include <math.h>
#include <fstream>

#include "SphereIkosahedron.h"

#define X .525731112119133606
#define Z .850650808352039932

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

//Vertices, triangles, edges of a single icosahedron
static double vert[12][3] = {
	{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
	{0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
	{Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}
};
static short edge[30][2] = {
	{0,1}, {0,4}, {0,6}, {0,9}, {0,11}, {1,4}, {1,6}, {1,8}, {1,10}, {2,3},
	{2,5}, {2,7}, {2,9}, {2,11}, {3,5}, {3,7}, {3,8}, {3,10}, {4,5}, {4,8},
	{4,9}, {5,8}, {5,9}, {6,7}, {6,10}, {6,11}, {7,10}, {7,11}, {8,10}, {9,11}
};
static short triang[20][3] = {
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
	{8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
};

using namespace itk;
using namespace std;

//************************************************************************************************

template < typename T >
unsigned short SphereIkosahedron< T >::GetNumberOfVertices()
{
	return m_NumberOfVertices;
}

//************************************************************************************************

template < typename T >
void SphereIkosahedron< T >::ComputeSubdivision()
{
	int i;
	
	m_NumberOfVertices = 12 + 30*m_SubdivisionLevel + 20*m_SubdivisionLevel*(m_SubdivisionLevel-1)/2;//12 is the number of vertices of the simple iko, 30: number of edge, 20: number of triangles
	m_NumberOfTriangles = 20*(m_SubdivisionLevel+1)*(m_SubdivisionLevel+1);

	unsigned short subdiv = m_SubdivisionLevel + 1;
	int m,n;
	double x1, y1, z1, x2, y2, z2, x3, y3, z3, dx12, dy12, dz12, dx23, dy23, dz23, length;
	for(i=0; i<30; i++) 
	{
		//For each edge, we divide it into 2 equal ones
		x1 = vert[edge[i][0] ][0];
		y1 = vert[edge[i][0] ][1];
		z1 = vert[edge[i][0] ][2];
		x2 = vert[edge[i][1] ][0];
		y2 = vert[edge[i][1] ][1];
		z2 = vert[edge[i][1] ][2];
		dx12 = (x2 - x1)/subdiv;
		dy12 = (y2 - y1)/subdiv;
		dz12 = (z2 - z1)/subdiv;
		for(n=1; n<subdiv; n++) 
		{
			//For each level of subdivision, we subdivide each edge into 2 smaller ones
			VectorType Point, PhiTheta;
			Point.push_back(x1 + n*dx12);
			Point.push_back(y1 + n*dy12);
			Point.push_back(z1 + n*dz12);
			
			length = sqrt((double) Point[0]*Point[0]+
					Point[1]*Point[1]+
					Point[2]*Point[2]);
			//Normalisation (to make sure each new point is at a distance of 1 from the center so that we approximate a sphere)
			Point[0] /= length;
			Point[1] /= length;
			Point[2] /= length;
			m_CoordinateTable.push_back(Point);
			//From cartesian to polar coordinates
			double temp_phi = atan2(Point[1] , Point[0]);
			while(temp_phi >= M_PI)
			{
				temp_phi -= 2*M_PI;
			}
			while(temp_phi < -M_PI)
			{
				temp_phi += 2*M_PI;
			}
			PhiTheta.push_back(temp_phi);//phi
			PhiTheta.push_back(acos(Point[2]));//theta
		
			m_PhiThetaTable.push_back(PhiTheta);
		}
	}

	if(subdiv > 2) 
	{
		//If the subdivision level is > 2, we have to create points inside each basic triangle
		for(i=0; i<20; i++) 
		{
			//For each triangle, we subdivide each edge.
			x1 = vert[triang[i][0] ][0];
			y1 = vert[triang[i][0] ][1];
			z1 = vert[triang[i][0] ][2];
			x2 = vert[triang[i][1] ][0];
			y2 = vert[triang[i][1] ][1];
			z2 = vert[triang[i][1] ][2];
			x3 = vert[triang[i][2] ][0];
			y3 = vert[triang[i][2] ][1];
			z3 = vert[triang[i][2] ][2];
			dx12 = (x2 - x1)/subdiv;
			dy12 = (y2 - y1)/subdiv;
			dz12 = (z2 - z1)/subdiv;
			dx23 = (x3 - x2)/subdiv;
			dy23 = (y3 - y2)/subdiv;
			dz23 = (z3 - z2)/subdiv;

			n = 1;
			do 
			{
				for(m=1; m<=n; m++) 
				{
					//Then, we subdivide each distance between 2 opposite points from 2 edges of the same triangle
					VectorType Point, PhiTheta;
					Point.push_back(x1 + (n+1)*dx12 + m*dx23);
					Point.push_back(y1 + (n+1)*dy12 + m*dy23);
					Point.push_back(z1 + (n+1)*dz12 + m*dz23);
					//Normalisation
					length = sqrt((double) Point[0]*Point[0]+
							Point[1]*Point[1]+
							Point[2]*Point[2]);
					Point[0] /= length;
					Point[1] /= length;
					Point[2] /= length;
					m_CoordinateTable.push_back(Point);
					//From cartesian to polar coordinates
					double temp_phi = atan2(Point[1] , Point[0]);
					while(temp_phi >= M_PI)
					{
						temp_phi -= 2*M_PI;
					}
					while(temp_phi < -M_PI)
					{
						temp_phi += 2*M_PI;
					}
					PhiTheta.push_back(temp_phi);//phi
					PhiTheta.push_back(acos(Point[2]));//theta
		
					m_PhiThetaTable.push_back(PhiTheta);
				}
				n++;
			}while( n<=(subdiv-2) );
		}
	}
	//Create the triangle from the new made vertices
	if (subdiv > 1) 
	{
		for(i=0; i<20; i++) 
		{
			x1 = vert[triang[i][0] ][0];
			y1 = vert[triang[i][0] ][1];
			z1 = vert[triang[i][0] ][2];
			x2 = vert[triang[i][1] ][0];
			y2 = vert[triang[i][1] ][1];
			z2 = vert[triang[i][1] ][2];
			x3 = vert[triang[i][2] ][0];
			y3 = vert[triang[i][2] ][1];
			z3 = vert[triang[i][2] ][2];
			dx12 = (x2 - x1)/subdiv;
			dy12 = (y2 - y1)/subdiv;
			dz12 = (z2 - z1)/subdiv;
			dx23 = (x3 - x2)/subdiv;
			dy23 = (y3 - y2)/subdiv;
			dz23 = (z3 - z2)/subdiv;

			n = 1;
			do 
			{
				for(m=1; m<=n; m++) 
				{
					// Draw lower triangle
					VectorType triangle1, triangle2, triangle3;
					std::vector< VectorType > zero_triangs;
					
					triangle1.push_back(x1 + n*dx12 + m*dx23);
					triangle1.push_back(y1 + n*dy12 + m*dy23);
					triangle1.push_back(z1 + n*dz12 + m*dz23);
					length = sqrt( triangle1[0]*triangle1[0] + triangle1[1]*triangle1[1] + triangle1[2]*triangle1[2]);
					triangle1[0] /= length;
					triangle1[1] /= length;
					triangle1[2] /= length;
					zero_triangs.push_back(triangle1);
					
					triangle2.push_back(x1 + (n-1)*dx12 + (m-1)*dx23);
					triangle2.push_back(y1 + (n-1)*dy12 + (m-1)*dy23);
					triangle2.push_back(z1 + (n-1)*dz12 + (m-1)*dz23);
					length = sqrt( triangle2[0]*triangle2[0] + triangle2[1]*triangle2[1] + triangle2[2]*triangle2[2]);
					triangle2[0] /= length;
					triangle2[1] /= length;
					triangle2[2] /= length;
					zero_triangs.push_back(triangle2);
					
					triangle3.push_back(x1 + n*dx12 + (m-1)*dx23);
					triangle3.push_back(y1 + n*dy12 + (m-1)*dy23);
					triangle3.push_back(z1 + n*dz12 + (m-1)*dz23);
					length = sqrt( triangle3[0]*triangle3[0] + triangle3[1]*triangle3[1] + triangle3[2]*triangle3[2]);
					triangle3[0] /= length;
					triangle3[1] /= length;
					triangle3[2] /= length;
					zero_triangs.push_back(triangle3);
					
					m_all_triangs.push_back(zero_triangs);
					
					if ( m != n ) 
					{
						// Draw lower left triangle
						VectorType triangle4, triangle5, triangle6;
						std::vector< VectorType > one_triangs;
						
						triangle4.push_back(x1 + n*dx12 + m*dx23);
						triangle4.push_back(y1 + n*dy12 + m*dy23);
						triangle4.push_back(z1 + n*dz12 + m*dz23);
						length = sqrt( triangle4[0]*triangle4[0] + triangle4[1]*triangle4[1] + triangle4[2]*triangle4[2]);
						triangle4[0] /= length;
						triangle4[1] /= length;
						triangle4[2] /= length;
						one_triangs.push_back(triangle4);
						
						triangle5.push_back(x1 + (n-1)*dx12 + m*dx23);
						triangle5.push_back(y1 + (n-1)*dy12 + m*dy23);
						triangle5.push_back(z1 + (n-1)*dz12 + m*dz23);
						length = sqrt( triangle5[0]*triangle5[0] + triangle5[1]*triangle5[1] + triangle5[2]*triangle5[2]);
						triangle5[0] /= length;
						triangle5[1] /= length;
						triangle5[2] /= length;
						one_triangs.push_back(triangle5);
						
						triangle6.push_back(x1 + (n-1)*dx12 + (m-1)*dx23);
						triangle6.push_back(y1 + (n-1)*dy12 + (m-1)*dy23);
						triangle6.push_back(z1 + (n-1)*dz12 + (m-1)*dz23);
						length = sqrt( triangle6[0]*triangle6[0] + triangle6[1]*triangle6[1] + triangle6[2]*triangle6[2]);
						triangle6[0] /= length;
						triangle6[1] /= length;
						triangle6[2] /= length;
						one_triangs.push_back(triangle6);
						
						m_all_triangs.push_back(one_triangs);
					}
				}
				n++;
			} while( n<=subdiv );
		}
	}
	double epsilon = 0.00001;
	
	for (i = 0 ; i < m_NumberOfTriangles ; i++)
	{
		std::vector<short> triangs;
		for( int l = 0; l < 3 ; l++ )
		{
			triangs.push_back(-1);
		}
		m_OrdinatedTriangles.push_back(triangs);
	}
	// find indexes
	for(i = 0; i < m_NumberOfVertices; i++) 
	{
		for (int j = 0; j < m_NumberOfTriangles ; j++) 
		{
			for(int k = 0 ; k < 3 ; k++)
			{
				if (m_OrdinatedTriangles[j][k] < 0)
				{
					if ( (fabs(m_CoordinateTable[i][0] - m_all_triangs[j][k][0]) < epsilon) && 
										(fabs(m_CoordinateTable[i][1] - m_all_triangs[j][k][1]) < epsilon) && 
										(fabs(m_CoordinateTable[i][2] - m_all_triangs[j][k][2]) < epsilon ) ) 
					{
						m_OrdinatedTriangles[j][k] = i;
					}
				}
			}
		}
	}
}

//************************************************************************************************

template < typename T >
void SphereIkosahedron< T >::Initialisation_tables()
{
	int i = 0;
	
	for( i = 0 ; i < 12 ; i++ )
	{
		VectorType PhiTheta, Coordinate;
		
		Coordinate.push_back(vert[i][0]);
		Coordinate.push_back(vert[i][1]);
		Coordinate.push_back(vert[i][2]);
		m_CoordinateTable.push_back(Coordinate);
		
		double temp_phi = atan2(Coordinate[1] , Coordinate[0]);
		while(temp_phi >= M_PI)
		{
			temp_phi -= 2*M_PI;
		}
		while(temp_phi < -M_PI)
		{
			temp_phi += 2*M_PI;
		}
		PhiTheta.push_back(temp_phi);//phi
		PhiTheta.push_back(acos(Coordinate[2]));//theta
		m_PhiThetaTable.push_back(PhiTheta);
	}
	if(m_SubdivisionLevel > 0)
	{
		ComputeSubdivision();
		CreateVTKFile();
	}
	else
	{
		m_NumberOfVertices = 12;
//CreateVTKFile () ;
	}
}

//************************************************************************************************


template< typename T >
int SphereIkosahedron< T >::PhiThetaToIndex(double Phi, double Theta)
{
	int indice = -1;
	int i,j;
	
	for(i = 0 ; i < m_NumberOfVertices ; i++)
	{
		if( Phi == m_PhiThetaTable[i][0] && Theta == m_PhiThetaTable[i][1])
		{
			indice = i;
		}
	}
	if(indice == -1)
	{
		cout<<"Values of Phi and/or Theta not correct, index not found."<<endl;
	}
	return indice;
}

//************************************************************************************************

template < typename T >
std::vector<double> SphereIkosahedron< T >::GetPhiThetaTableatIndex(short index)
{
	return m_PhiThetaTable[index];
}

//************************************************************************************************

template < typename T >
void SphereIkosahedron< T >::CreateVTKFile()
{
	ofstream myfile;
	myfile.open ("../VTKFiles/myVTKFile.vtk");
	myfile << "# vtk DataFile Version 5.2\n";
	//myfile << "Ikosahedron subdivision. Subdivision level: "<<m_SubdivisionLevel<<"\n";
        myfile << "vtk output" << std::endl ;
	myfile << "ASCII\n";
	myfile << "DATASET POLYDATA\n";
	myfile << "POINTS "<<m_NumberOfVertices<<" float\n";
	for(int i = 0 ; i < m_NumberOfVertices ; i++)
	{
		myfile<<m_CoordinateTable[i][0]<<" "<<m_CoordinateTable[i][1]<<" "<<m_CoordinateTable[i][2]<<"\n";
	}
	myfile << "POLYGONS "<<m_NumberOfTriangles<<" "<<m_NumberOfTriangles*4<<"\n";
	for(int j = 0 ; j < m_NumberOfTriangles ; j++)
	{
		myfile<<"3 "<<m_OrdinatedTriangles[j][0]<<" "<<m_OrdinatedTriangles[j][1]<<" "<<m_OrdinatedTriangles[j][2]<<"\n";
	}
	myfile << "\nCELL_DATA "<<m_NumberOfTriangles<<"\n";
	myfile << "POINT_DATA "<<m_NumberOfVertices<<"\n\n";
	myfile.close();
}

//************************************************************************************************

template < typename T >
std::vector<double> SphereIkosahedron< T >::GetCoordinateTableatIndex(short index)
{
	return m_CoordinateTable[index];
}

//************************************************************************************************

#endif
