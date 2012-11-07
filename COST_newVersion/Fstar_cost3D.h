#ifndef DEF_FSTAR_COST3D
#define DEF_FSTAR_COST3D

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>

#include "Fstar_image_creator.h"

class Fstar_cost3D
{
	public:
		typedef double 						CostType;
		typedef double                                          LengthType ;
		//typedef SphereIkosahedronImage< CostType > 		VectorImageType;
		typedef std::vector<double> 				VectorType;
		typedef double 						ImageElementType;
		typedef itk::Image < ImageElementType , 3 > 		ImageType;
		typedef itk::Image < LengthType, 3 >                    LengthImageType ;
		typedef itk::ImageRegionIterator < Fstar_cost3D::VectorImageType >    IteratorType;
		typedef itk::ImageFileWriter < Fstar_cost3D::ImageType > WriterType;
		typedef itk::ImageFileWriter < LengthImageType > 	LengthWriterType;

		Fstar_cost3D(VectorImageType::Pointer,ImageType::Pointer, long unsigned int, double);//overloaded constructor
		~Fstar_cost3D();

		//Return the FINAL RESULT connectivity map
		void 				                        GetCostMap(std::string mincost, std::string lengthmap, std::string originmap, std::string avgcost, std::string mincostperminlength);
		
	protected:
		void 							Initialisation();//Launch the precomputation of SphereIkosahedron stuff
		void PrecomputeCosts () ;
		void FstarInitialization () ;
		void 							ComputeCostMap();//Main algorithm, also count the number of changes
		void 							ComputeFinalMinImage( std::string mincost_filename, std::string length_filename, std::string orig_filename, std::string avgcostfilename, std::string mincostperminlengthfilename);
		void 							Fstar_Cost_map_algorithm_3D();//F-star style algorithm
		void 							Voxel_algorithm_3D();
		void 							ChangeFlagOptimization();
		void							AlgorithmProcessing();
		double                                                  f_ODF ( double odf ) ;
		double GetToNeighbourDirectionSimilarity ( double x, double y, double z, unsigned short direction_index ) ;
		double InOutDirectionSimilarity ( unsigned short in, unsigned short out ) ;
		double GetAngleCoefficient ( double AngleSimilarityTemp, double mu, double sigma) ;
		void CheckSanity () ;

		SphereIkosahedron<CostType>::Pointer        		m_SphereIkosahedronObject;
		VectorImageType::Pointer 				m_CurrentImage;
		VectorImageType::Pointer 				m_CurrentLengthImage;
		VectorImageType::Pointer 				m_CurrentAvgCostImage;
		VectorImageType::Pointer                                m_CurrentOrigImage ;
		CostType *						m_CostMapArray;//the image array storing the data at time "n-1"
		CostType *						m_LengthMapArray;//the image array storing the data at time "n-1"
		CostType *						m_AvgCostMapArray;//the image array storing the data at time "n-1"
		CostType *                                              m_OrigMapArray ;
		CostType *						m_ODFArray;//the input/data image array
		CostType 						m_CurrentCost;//feature of the current voxel
		ImageElementType *					m_LabelImageArray;
		std::vector<short>					m_ChangeFlag;
		std::vector<short>					m_PropagationFlag;

		unsigned long int 						m_NumberOfLoops;
		unsigned long int 						m_CurrentIndex;//Index of the current voxel
		unsigned long int 						m_NeighbourIndex;//index of current neighbour of current voxel
		unsigned long int 						m_CostMapChange;//how many changes are done per iteration
		unsigned long 						m_DirectionIterator, m_DirectionNeighbourIterator;//Direction iterators
		unsigned long 						m_DirectionToNeighbourIterator;
		unsigned long 						m_SLICE, m_LINE ;
		unsigned long                                           m_ZFactor, m_YFactor ;
		unsigned long 						m_XDimension, m_YDimension, m_ZDimension;//dimensions of the image
		unsigned long 						m_NumberOfDirections;//total number of directions available
		 long 							m_x_it, m_y_it, m_z_it;//iterators
		long 							m_x_temp_it, m_y_temp_it, m_z_temp_it;//temporary iterators
		double 							m_WeightCoefficient;//neighbourhood participation + angle coeff
		double 							m_Epsilon;//error allowed to know if we stop the algorithm
		double							m_ThresholdCoefficient;
		bool 							m_ChangingFlag;
		VectorImageType::SpacingType 				m_Spacing;
		std::vector < std::vector < std::vector < std::vector < double > > > > m_NeighborDirCosts ;
		std::vector < std::vector < double > > m_DirDirCosts ;
		CostType m_NeighborDist[3][3][3] ;
		long int m_fstar_updates[8][14] ;
		long int m_neighbors[27][3] ;
		long int m_fstar_direction ;
		unsigned long int m_MinCost ;
};


#endif
