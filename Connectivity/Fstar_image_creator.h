#ifndef DEF_FSTAR_IMAGE_CREATOR
#define DEF_FSTAR_IMAGE_CREATOR

#include <iostream>
#include <vector>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkRGBPixel.h>
#include <itkImageFileWriter.h>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"

#include "SphereIkosahedronImage.h"
#include "Fstar_prob3D.h"

using namespace itk;

class Fstar_image_creator
{
	public:
		typedef double 					FeatureType;
		typedef double 					ImageElementType;
		typedef SphereIkosahedronImage<FeatureType> 	VectorImageType;
		typedef itk::Image < ImageElementType , 3 > 	ImageType;
		typedef vnl_vector_fixed<int,2>			LmVector;
		typedef vnl_vector_fixed<double,3>		VectorType;
		
		typedef itk::RGBPixel<unsigned char> RGBPixelType;
		typedef itk::Image<RGBPixelType,3> RGBImageType;
		typedef itk::ImageFileWriter<RGBImageType> RGBImageWriterType;


		Fstar_image_creator(ImageType::Pointer, long unsigned, double);//with synthetic labeled image
		Fstar_image_creator(VectorImageType::Pointer, long unsigned, double);//with ODF
		~Fstar_image_creator() {}
		
		void 					Image_creator();
		VectorImageType::Pointer 		GetOutputImage();
		VectorImageType::Pointer 		CreateVectorImage();
		
	protected:
		double 					GetDirectionSimilarity( VectorType, unsigned long);
		double 					GetAngleCoefficient(double, double, double);
		double 					EvaluateBasis(long int, double, double);
		double 					Y(long int, double, double);//ODF conversion
		double 					K(long int, long int);//ODF conversion
		double 					LegendreP( int, int, double);//ODF conversion
		void 					ComputeODFSumImage();//ODF conversion
		LmVector				GetLM(long unsigned int);//ODF conversion
		
		FeatureType *				m_OneCoefficientImageArray;
		RGBPixelType *				m_ColorFAImageArray;
		unsigned long 				m_XDimension, m_YDimension, m_ZDimension;
		unsigned long				m_x_it, m_y_it, m_z_it;
		unsigned long 				m_LINE, m_SLICE;
		unsigned long 				m_Direction, m_NumberOfVertices;
		long 					m_SubdivisionLevel;
		int 					m_Order;
		int 					m_NumberOfSphericalHarmonics;
		int					m_CurrentIndex, m_OutputIndex;
		double					m_ThresholdCoefficient;
		SphereIkosahedron<double>::Pointer	m_SphereIkosahedronObject;
		ImageElementType			m_CurrentFeature;
		ImageElementType * 			m_LabeledImageArray;
		ImageType::Pointer			m_SumImage;
		ImageType::Pointer			m_OneCoefficientImage;
		RGBImageType::Pointer			m_ColorFAImage;
		VectorImageType::Pointer 		m_IkoOutputImage;
		FeatureType *		 		m_ODFImageArray;
		FeatureType *				m_SumImageArray;
		FeatureType * 				m_IkoOutputArray;
		vnl_matrix<double> 			m_RSHBasisMatrix;
};

#endif
