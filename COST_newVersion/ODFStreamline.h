#ifndef DEF_ODF_STREAMLINE
#define DEF_ODF_STREAMLINE

#include <iostream>
#include <vector>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"

class ODFStreamline
{
   public:
      typedef double 				 ODFType;
      typedef itk::VectorImage <ODFType, 3>      ODFImageType;
      typedef itk::Image <ODFType, 3>            DoubleImageType ;
      typedef itk::Image <unsigned long, 3 >     LabelImageType ;
      typedef vnl_vector_fixed<double, 3>        CoordinateType ;

      ODFStreamline ( ODFImageType::Pointer, unsigned long ) ;
      ~ODFStreamline() {} ;
		
      LabelImageType::Pointer Streamline ( LabelImageType::Pointer, unsigned long, bool ) ;

      void SetFAImage ( DoubleImageType::Pointer fa )
      {
	m_FAImage = fa ;
	m_FAArray = fa->GetPixelContainer()->GetBufferPointer() ;
      }
      void SetCoordinateTable ( std::vector < CoordinateType > vec ) 
      {
          m_CoordinateTable = vec ;
      }

   protected:
      unsigned long 	         m_XDim, m_YDim, m_ZDim ;
      unsigned long              m_NumberOfDirs ;
      ODFImageType::Pointer      m_ODFImage ;
      LabelImageType::Pointer    m_StreamlineImage ;
      DoubleImageType::Pointer   m_FAImage ;
      unsigned long              m_NumberOfSourceVoxels ;
      ODFType                    *m_ODFArray ;
      unsigned long              *m_StreamlineArray, *m_SourceArray ;
      double                     *m_FAArray ;
      std::vector < CoordinateType >  m_SourceVoxels ;
      std::vector < CoordinateType > m_CoordinateTable ;
      double m_StepSize ;

      void StreamlineFromVoxel ( unsigned long, unsigned long, unsigned long, bool ) ;
      void FindSourceVoxels () ;
      void WriteSourceVoxels () ;
      bool FacingAway ( int i, int j ) ;
      int GetOpposite ( unsigned int i ) ;
      long unsigned GetPeakDir ( unsigned long ) ;
};

#endif


