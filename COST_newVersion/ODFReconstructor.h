#ifndef DEF_ODF_RECONSTRUCTOR
#define DEF_ODF_RECONSTRUCTOR

#include <iostream>
#include <vector>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkRGBPixel.h>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"

using namespace itk;

class ODFReconstructor
{
   public:
      typedef double 				ODFType;
      typedef itk::VectorImage <ODFType, 3>     ODFImageType;
      typedef itk::Image <ODFType, 3>           ImageType ;
      typedef vnl_vector_fixed<int,2>		LmVector;
      typedef vnl_vector_fixed<double,3>	CoordinateType;
      typedef itk::RGBPixel<unsigned char>      RGBPixelType ;
      typedef itk::Image <RGBPixelType, 3>      RGBImageType ;

      ODFReconstructor ( ODFImageType::Pointer, long unsigned, long unsigned, std::string ) ;
      ~ODFReconstructor() {}
		
      ODFImageType::Pointer 	        	ReconstructODFImage();
      ImageType::Pointer                        GetFAImage () ;
      ODFReconstructor::RGBImageType::Pointer   GetColorFAImage () ;
      std::vector < CoordinateType >            GetCoordinateTable () 
	{
	  return m_CoordinateTable ;
	} 

   protected:
      double 					Y(long int, double, double);//ODF conversion
      double 					K(long int, long int);//ODF conversion
      double 					LegendreP( int, int, double);//ODF conversion
      LmVector	         			GetLM(long unsigned int);//ODF conversion
		
      unsigned long 				m_XDimension, m_YDimension, m_ZDimension;
      unsigned int                              m_NumberOfVertices ;
      unsigned int 				m_NumberOfSphericalHarmonics;
      ODFImageType::Pointer 		        m_ODFImage, m_CoefsImage;
      ODFType 		 	        	*m_ODFArray, *m_CoefsArray ;
      vnl_matrix<double> 			m_RSHBasisMatrix;
      std::vector < CoordinateType >            m_PhiThetaTable, m_CoordinateTable ;

      double m_GlobalMin, m_GlobalMax, m_GlobalRange ;

      void InitializeFromFile ( std::string ) ;
      void MakeSymmetric () ;
      void WritePhiTheta () ;
      void CheckOrder () ;
      ODFReconstructor::RGBPixelType DirToRGB ( unsigned long, double ) ;
};

#endif
