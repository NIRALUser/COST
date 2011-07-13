#ifndef DEF_FSTAR_PROB3D
#define DEF_FSTAR_PROB3D

#include <itkVector.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>

#include "SphereIkosahedronImage.h"
#include "Fstar_image_creator.h"

class Fstar_prob3D
{
	public:
		typedef double 						ProbType;
		typedef double                                          LengthType ;
		typedef SphereIkosahedronImage< ProbType > 		VectorImageType;
		typedef std::vector<double> 				VectorType;
		typedef double 						ImageElementType;
		typedef itk::Image < ImageElementType , 3 > 		ImageType;
		typedef itk::Image < LengthType, 3 >                    LengthImageType ;
		typedef itk::ImageRegionIterator < Fstar_prob3D::VectorImageType >    IteratorType;
		typedef itk::ImageFileWriter < Fstar_prob3D::ImageType > WriterType;
		typedef itk::ImageFileWriter < LengthImageType > 	LengthWriterType;

		Fstar_prob3D(VectorImageType::Pointer,ImageType::Pointer, long unsigned int, double);//overloaded constructor
		~Fstar_prob3D();

		//Return the FINAL RESULT connectivity map
		void 				                        GetProbMap(std::string mincost ) ;
		
	protected:
		void 							Initialisation();//Launch the precomputation of SphereIkosahedron stuff
		void PrecomputeProbs () ;
		void FstarInitialization () ;
		void 							ComputeProbMap();//Main algorithm, also count the number of changes
		void 							ComputeFinalProbImage( std::string mincost_filename ) ;

		void 							Fstar_Prob_map_algorithm_3D();//F-star style algorithm
		void 							Voxel_algorithm_3D();
		void Normalize () ;
		void							AlgorithmProcessing();
		double GetToNeighbourDirectionSimilarity ( short x, short y, short z, unsigned short direction_index ) ;
		double InOutDirectionSimilarity ( unsigned short in, unsigned short out ) ;
		double GetAngleCoefficient ( double AngleSimilarityTemp, double mu, double sigma) ;

		SphereIkosahedron<ProbType>::Pointer        		m_SphereIkosahedronObject;
		VectorImageType::Pointer 				m_CurrentImage;

		ProbType *						m_ProbMapArray;//the image array storing the data at time "n-1"

		ProbType *						m_ODFArray;//the input/data image array
		ProbType 						m_CurrentProb;//feature of the current voxel
		ImageElementType *					m_LabelImageArray;
		std::vector<short>					m_PropagationFlag;

		unsigned long int 						m_NumberOfLoops;
		unsigned long int 						m_CurrentIndex;//Index of the current voxel
		unsigned long int 						m_NeighbourIndex;//index of current neighbour of current voxel
		unsigned long int 						m_ProbMapChange;//how many changes are done per iteration
		unsigned long 						m_out, m_DirectionNeighbourIterator;//Direction iterators
		unsigned long 						m_DirectionToNeighbourIterator;
		unsigned long 						m_SLICE, m_LINE ;
		unsigned long                                           m_ZFactor, m_YFactor ;
		unsigned long 						m_XDimension, m_YDimension, m_ZDimension;//dimensions of the image
		unsigned long 						m_NumberOfDirections;//total number of directions available
		 long 							m_x_it, m_y_it, m_z_it;//iterators
		double 							m_WeightCoefficient;//neighbourhood participation + angle coeff
		double 							m_Epsilon;//error allowed to know if we stop the algorithm
		double							m_ThresholdCoefficient;
		bool 							m_ChangingFlag;
		VectorImageType::SpacingType 				m_Spacing;
		std::vector < std::vector < std::vector < std::vector < double > > > > m_NeighborDirProbs ;
		std::vector < std::vector < double > > m_DirDirProbs ;
		std::vector < std::vector < double > > m_ParticipationToNeighbor ;
		ProbType m_NeighborDist[3][3][3] ;
		long int m_neighbors[27][3] ;
		int m_reverse[27] ;
		std::vector < std::vector < unsigned long int > > m_neighborIndex ;
		std::vector < std::vector < bool > > m_exists ;
		std::vector < ProbType > m_prevValues ;
};


#endif
