#ifndef SphereIkosohedron_h
#define SphereIkosohedron_h

#include <iostream>
#include <itkDataObject.h>
#include <itkObjectFactory.h>
#include <vector>

namespace itk {
template < typename T >
class ITK_EXPORT SphereIkosahedron : public DataObject
{
	public:
		typedef std::vector< double >		VectorType;
		typedef SphereIkosahedron 		Self;
		typedef DataObject 			Superclass;
		typedef SmartPointer<Self> 		Pointer;
		typedef SmartPointer<const Self> 	ConstPointer;
		
		itkNewMacro(Self);
		itkTypeMacro(SphereIkosahedron, DataObject);
		itkSetMacro(SubdivisionLevel, unsigned short);
		itkGetMacro(SubdivisionLevel, unsigned short);
		
		void 					Initialisation_tables();
		void 					CreateVTKFile();
		unsigned short 				GetNumberOfVertices();
		VectorType 				GetCoordinateTableatIndex(short);
		VectorType 				GetPhiThetaTableatIndex(short);
		int 					PhiThetaToIndex(double, double);
		
	protected:
		void 					ComputeSubdivision();
		SphereIkosahedron(){}
		~SphereIkosahedron(){}
		
		std::vector< VectorType > 		m_CoordinateTable;
		std::vector< VectorType > 		m_PhiThetaTable;//Phi, first coordinate, theta second one
		std::vector< std::vector< short > >	m_OrdinatedTriangles;
		std::vector<std::vector< VectorType > >	m_all_triangs;
		unsigned short 				m_SubdivisionLevel;
		unsigned short 				m_NumberOfVertices;
		short 					m_NumberOfTriangles;
		T *					m_FeatureTable;
		
	private:
		SphereIkosahedron(const Self&); //purposely not implemented
		void 	operator=(const Self&); //purposely not implemented
};
}//end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "SphereIkosahedron.txx"
#endif

#endif
