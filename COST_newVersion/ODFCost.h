#ifndef DEF_ODF_COST
#define DEF_ODF_COST

#include <iostream>
#include <vector>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"
#include "math.h"

#define INF 99999
#define PI 3.1415926

//Double image definition
typedef itk::Image < double, 3 >                DoubleImageType ;
typedef itk::ImageFileWriter < DoubleImageType > DoubleWriterType ;

class ODFCost
{
   public:
      typedef double 				 ODFType;
      typedef itk::VectorImage <ODFType, 3>      ODFImageType;
      typedef itk::Image <ODFType, 3>            DoubleImageType ;
      typedef itk::Image <unsigned long, 3 >     LabelImageType ;
      typedef vnl_vector_fixed<double, 3>        CoordinateType ;

      //ODFCost ( ODFImageType::Pointer, unsigned long,double ) ;
      //modified by Wenyu
      ODFCost ( ODFImageType::Pointer,DoubleImageType::Pointer, unsigned long,double ) ;

      ~ODFCost() {} ;
		
      DoubleImageType::Pointer Cost ( LabelImageType::Pointer ) ;

      void SetFAImage ( DoubleImageType::Pointer fa )
      {
	m_FAImage = fa ;
	m_FAArray = fa->GetPixelContainer()->GetBufferPointer() ;
      }
      void SetCoordinateTable ( std::vector < CoordinateType > vec ) 
      {
          m_CoordinateTable = vec ;
	  this->PrecomputeAnglePenalty () ;
      }

      DoubleImageType::Pointer GetOrigImage ()
	{
	  return m_FinalOrigImage ;
	}

      DoubleImageType::Pointer GetLengthImage ()
	{
	  return m_FinalLengthImage ;
	}

   protected:
      unsigned long 	         m_XDim, m_YDim, m_ZDim, m_Slice ;
      unsigned long              m_NumberOfDirs ;
      ODFImageType::Pointer      m_ODFImage, m_ODFQCImage ;
      ODFImageType::Pointer      m_CostImage, m_LengthImage, m_AvgCostImage, m_OrigImage ;
      DoubleImageType::Pointer   m_FinalCostImage, m_FinalOrigImage, m_FinalLengthImage ;
      DoubleImageType::Pointer   m_FAImage ;

      //type is changed by wenyu from unsigned long to unsiged long long
     
      unsigned long long          m_NumberOfSourceVoxels ;
      ODFType                    *m_ODFArray, *m_ODFQCArray ;
      unsigned long              *m_SourceArray ;
      double                     *m_FAArray ;
      ODFType                    *m_CostArray, *m_LengthArray, *m_AvgCostArray, *m_OrigArray ;
      double                     *m_FinalCostArray, *m_FinalOrigArray, *m_FinalLengthArray ; 
      double                      m_alpha ;
 
  
      double                       m_startx,m_starty,m_startz;
      ODFType                     *m_Cost;
      DoubleImageType::Pointer     m_OrigFAImage ;
      double                      *m_OrigFAArray;
      ODFType                     *m_OrigInDir;
      ODFType                     *m_OrigOutDir;
      DoubleImageType::Pointer      m_OrigInDirImage, m_OrigOutDirImage ;

      std::vector < CoordinateType >  m_SourceVoxels ;
      std::vector < CoordinateType > m_CoordinateTable ;
      std::vector < double > m_DNTable ;
      std::vector < std::vector < double > > m_AnglePenaltyTable ;
      std::vector < std::vector < double > > m_NeighborAnglePenaltyTable ;
      std::vector < std::vector < double > > m_NeighborPenaltyTable ;
      
      std::vector <double> m_FinslerTable ;
      std::vector < std::vector < int > > m_DirectionTable;

      std::vector < std::vector < unsigned long > > m_NeighborIndexTable ;
      std::vector < std::vector < bool > > m_NeighborValidTable ;

      bool m_Converged ;
      unsigned long m_NChanged ;

      void FindSourceVoxels () ;

      void MakeOnePass () ;
      void MakeOnePassFstar () ;
      void VisitVoxel ( unsigned long, unsigned long, unsigned long ) ;
      void InitializeSourceVoxels () ;
      void CreateCostMaps () ;
      void ComputeFinalCosts () ;

      unsigned long GetNeighbor ( unsigned long, unsigned long, unsigned long, unsigned short, bool &) ;
     
      void UpdateCost ( unsigned long, unsigned long, double, double, double, unsigned int, unsigned int, unsigned int, unsigned int) ;

      double f_ODF ( unsigned long ) ;
      double AnglePenalty ( unsigned, unsigned ) ;
      double NeighborAnglePenalty ( unsigned, unsigned ) ;
      double stepLength ( unsigned ) ;
      void PrecomputeStepLength () ;
      void PrecomputeAnglePenalty () ;
      void PrecomputeNeighbors () ;
      void NormalizeODF () ;

      void CreateFinslerTable();
      double DirectionPenalty(unsigned long , unsigned long, unsigned long, unsigned);
      void PrecomputeDirection();
      
      double NeighborPenalty( unsigned, unsigned );

      void FstarTraversalUpdates () ;
      long int m_fstar_updates[8][14] ;
     
      unsigned short m_fstar_direction ;

};


#endif


