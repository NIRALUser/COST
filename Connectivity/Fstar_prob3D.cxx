#include "Fstar_prob3D.h"
#include <itkImageDuplicator.h>
#include <math.h>
#include <time.h>

//on input image
#define SEED 1

#define INF 99999999

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

//***************************************************************************************

Fstar_prob3D::Fstar_prob3D(VectorImageType::Pointer Input_image , ImageType::Pointer InputLabelImage, unsigned long int subdivision_level, double thresholdCoefficient)
{
  m_ThresholdCoefficient = thresholdCoefficient ;
  m_XDimension = InputLabelImage->GetLargestPossibleRegion().GetSize()[0] ;
  m_YDimension = InputLabelImage->GetLargestPossibleRegion().GetSize()[1] ;
  m_ZDimension = InputLabelImage->GetLargestPossibleRegion().GetSize()[2] ;
  m_SLICE = m_XDimension*m_YDimension ;
  m_LINE = m_XDimension ;
  cout<<"xdim: "<<m_XDimension<<"   // ydim: "<<m_YDimension<<"   // zdim: "<<m_ZDimension<<endl ;
	
  //Icosahedron subdivision
  m_SphereIkosahedronObject = SphereIkosahedron<ProbType>::New() ;
  m_SphereIkosahedronObject->SetSubdivisionLevel(subdivision_level) ;
  m_SphereIkosahedronObject->Initialisation_tables();
  m_NumberOfDirections = m_SphereIkosahedronObject->GetNumberOfVertices();
  cout<<"Number of vertices/directions: "<<m_NumberOfDirections<<endl;
  m_ZFactor = m_SLICE * m_NumberOfDirections ;
  m_YFactor = m_LINE * m_NumberOfDirections ;
  m_Spacing = Input_image->GetSpacing();
	
  //The image that will be modified/worked on
  m_CurrentImage = VectorImageType::New();
  VectorImageType::IndexType start;
  itk::VariableLengthVector< ProbType > f( m_NumberOfDirections );
  itk::VariableLengthVector< LengthType > l( m_NumberOfDirections );
  VectorImageType::SizeType  size;
  for( unsigned int i=0; i<m_NumberOfDirections; i++ ) { f[i] = 0; l[i] = -1 ; }
  start[0] = start[1] = start[2] = 0;
  size[0] = m_XDimension;
  size[1] = m_YDimension;
  size[2] = m_ZDimension;
  VectorImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  m_CurrentImage->SetVectorLength( m_NumberOfDirections );
  m_CurrentImage->SetRegions( region );
  m_CurrentImage->Allocate();
  m_CurrentImage->FillBuffer( f );

  m_LabelImageArray = InputLabelImage->GetPixelContainer()->GetBufferPointer();
  m_ODFArray = Input_image->GetPixelContainer()->GetBufferPointer();
  m_ProbMapArray = m_CurrentImage->GetPixelContainer()->GetBufferPointer();

  m_Epsilon = 0;
  m_ChangingFlag = true;

  m_prevValues.resize ( m_NumberOfDirections ) ;

  PrecomputeProbs () ;
  FstarInitialization () ;
}

//***************************************************************************************

void Fstar_prob3D::Initialisation()
{
  unsigned long int sourceCounter = 0 ;
  m_ParticipationToNeighbor.resize ( m_ZDimension * m_YDimension * m_XDimension ) ;
  for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      unsigned long int counter = 0;
	      m_CurrentIndex = m_x_it + m_y_it*m_LINE + m_z_it*m_SLICE;
	      if ( m_LabelImageArray[m_CurrentIndex] == SEED )
		{
		  sourceCounter++ ;
		}

	      m_ParticipationToNeighbor[m_CurrentIndex].resize ( 27 ) ;
	      for ( int neighbor = 0 ; neighbor < 27 ; neighbor++ )
		m_ParticipationToNeighbor[m_CurrentIndex][neighbor] = 0 ;

	      for( m_out = 0 ; m_out < m_NumberOfDirections ; m_out++ )
		{
		  unsigned long int CurrentVoxelIndex = m_out + m_x_it*m_NumberOfDirections + m_y_it*m_YFactor + m_z_it*m_ZFactor;
					
		  if(m_ODFArray[CurrentVoxelIndex] > 0)
		    {
			//To check wether there is a propagation along that voxel, if not, no computation along it
		      counter++;
		    }
		  if(m_LabelImageArray[m_CurrentIndex] == SEED )
		    {
		    m_ProbMapArray[CurrentVoxelIndex] = 1.0/m_NumberOfDirections;//If voxel is labeled as source, set it to 0 along each direction
		    }
		  else
		    {
		    m_ProbMapArray[CurrentVoxelIndex] = 0;
		    }
		}
	      if(counter > 0)
		m_PropagationFlag.push_back(1);//propagation allowed
	      else
		m_PropagationFlag.push_back(0);//propagation denied 
	    }
	}
    }
  std::cout << "The seed was separated in " << sourceCounter << " for tracking purposes." << std::endl ;
}

void Fstar_prob3D::FstarInitialization ()
{
  short int x, y, z, temp, reverse ;
  for ( int i = 0 ; i < 27 ; i++ )
    {
      x = i % 3 ;
      m_neighbors[i][0] = x ;

      temp = ( i - x ) / 3 ;
      y = temp % 3 ;
      m_neighbors[i][1] = y ;

      z = ( temp - y ) / 3 ;
      m_neighbors[i][2] = z ;

      x = 2 - x ;
      y = 2 - y ;
      z = 2 - z ;
      reverse = z*9+y*3+x ;
      m_reverse[i] = reverse ;
    }

  short int nz, ny, nx ;
  long int z_dest, y_dest, x_dest, z_cur, y_cur, x_cur ;
  unsigned long int index = 0 ;
  m_exists.resize ( m_XDimension * m_YDimension * m_ZDimension ) ;
  m_neighborIndex.resize ( m_XDimension * m_YDimension * m_ZDimension ) ;
  for ( z_cur = 0 ; z_cur < m_ZDimension ; z_cur++ )
    {
      for ( y_cur = 0 ; y_cur < m_YDimension ; y_cur++ )
	{
	  for ( x_cur = 0 ; x_cur < m_XDimension ; x_cur++ )
	    {
	      index = x_cur + y_cur * m_LINE + z_cur * m_SLICE ;
	      m_exists[index].resize ( 27 ) ;
	      m_neighborIndex[index].resize ( 27 ) ;
	      for ( unsigned long int neighbor = 0 ; neighbor < 27 ; neighbor++ )
		{
		  m_exists[index][neighbor] = true ;
		  nz = m_neighbors[neighbor][2] ;
		  z_dest = z_cur + nz - 1 ;
		  if ( ( z_dest > m_ZDimension ) || ( z_dest <= 0 ) ) m_exists[index][neighbor] = false ;

		  ny = m_neighbors[neighbor][1] ;
		  y_dest = y_cur + ny - 1 ;
		  if ( ( y_dest > m_YDimension ) || ( y_dest <= 0 ) )  m_exists[index][neighbor] = false ;

		  nx = m_neighbors[neighbor][0] ;
		  x_dest = x_cur + nx - 1 ;
		  if ( ( x_dest > m_XDimension ) || ( x_dest <= 0 ) )  m_exists[index][neighbor] = false ;

		  if ( ( nx == 1 ) && ( ny == 1 ) && ( nz == 1 ) )  m_exists[index][neighbor] = false ;
		  if ( m_exists[index][neighbor] ) 
		    {
		      m_neighborIndex[index][neighbor] = z_dest * m_SLICE + y_dest * m_LINE + x_dest ;
		    }
		    
		}
	    }
	}
    }
  std::cout << "Initialization complete" << std::endl ;
}

void Fstar_prob3D::Voxel_algorithm_3D ()
{
  const unsigned long CurrentVoxelIndex = m_out + m_CurrentIndex * m_NumberOfDirections ;
  m_prevValues[m_out] = m_ProbMapArray[CurrentVoxelIndex] ;
  const double fODFCurrent = m_ODFArray[CurrentVoxelIndex] ;

  unsigned long int baseIndex ;
  ProbType neighborProb, stepProb ;
  ProbType currentProb = 0 ;
  short int nx, ny, nz ;
  // for each neighbor 
  for ( unsigned short int neighbor = 0 ; neighbor < 27 ; neighbor++ )
    {
      if ( ! m_exists[m_CurrentIndex][neighbor] ) continue ;
      nx = m_neighbors[neighbor][0] ;
      ny = m_neighbors[neighbor][1] ;
      nz = m_neighbors[neighbor][2] ;

      baseIndex = ( m_neighborIndex[m_CurrentIndex][neighbor] ) * m_NumberOfDirections ;
      ProbType stepLength = m_NeighborDist[nz][ny][nx] ;

      // for each direction in
      for( unsigned long int in = 0 ; in < m_NumberOfDirections ; in++ )
	{
	  m_NeighbourIndex = in + baseIndex ;
	  neighborProb = m_ProbMapArray[m_NeighbourIndex] ;
	  if ( neighborProb == 0 ) continue ;

	  stepProb =  stepLength * m_DirDirProbs[m_out][in] * m_NeighborDirProbs[nz][ny][nx][in] ;
	  currentProb += neighborProb * stepProb ;
	}
    }

  currentProb *= fODFCurrent ;
  if ( currentProb > 1 ) currentProb = 1 ;
  m_ProbMapArray[CurrentVoxelIndex] = currentProb ;
  return ;
}

void Fstar_prob3D::Normalize () 
{
  double IncomingTotal = 0, OutgoingTotal = 0 ;
  const unsigned long CurrentVoxelIndex = m_CurrentIndex * m_NumberOfDirections ;

  for ( int neighbor = 0 ; neighbor < 27 ; neighbor++ )
  {
    double tempOut = 0 ;
    int nx, ny, nz ;
    nx = m_neighbors[neighbor][0] ;
    ny = m_neighbors[neighbor][1] ;
    nz = m_neighbors[neighbor][2] ;
    // if this is a valid neighbor
    if ( !m_exists[m_CurrentIndex][neighbor] ) continue ;
    for ( m_out = 0 ; m_out < m_NumberOfDirections ; m_out++ )
      {
	tempOut += m_ProbMapArray[CurrentVoxelIndex+m_out] * m_NeighborDirProbs[nz][ny][nx][m_out] ;
      }

    m_ParticipationToNeighbor[m_CurrentIndex][neighbor] = tempOut ;
    OutgoingTotal += tempOut ;
  }

  for ( int neighbor = 0 ; neighbor < 27 ; neighbor++ )
    {
      int nx, ny, nz ;
      nx = m_neighbors[neighbor][0] - 1;
      ny = m_neighbors[neighbor][1] - 1;
      nz = m_neighbors[neighbor][2] - 1;

      // if this is a valid neighbor
      if ( !m_exists[m_CurrentIndex][neighbor] ) continue ;
      int neighbourIndex = m_neighborIndex[m_CurrentIndex][neighbor] ;

      int reverseNeighbor = m_reverse[neighbor] ;
      IncomingTotal += m_ParticipationToNeighbor[neighbourIndex][reverseNeighbor];
   }

  double normalizationFactor = IncomingTotal / OutgoingTotal ;

  if ( std::isinf ( normalizationFactor ) || std::isnan ( normalizationFactor ) || ( normalizationFactor >= 1 ) )
      return ;

  for ( m_out = 0 ; m_out < m_NumberOfDirections ; m_out++ )
    {
      unsigned long currentDirIndex = CurrentVoxelIndex+m_out ;
      m_ProbMapArray[currentDirIndex] *= normalizationFactor ;
    }
}

//***************************************************************************************

void Fstar_prob3D::AlgorithmProcessing()
{
  if ( m_PropagationFlag[m_CurrentIndex] == 1 )
    {
      m_CurrentProb = m_LabelImageArray[m_CurrentIndex] ;
      if ( m_CurrentProb != SEED )
	{
	  for( m_out = 0 ; m_out < m_NumberOfDirections ; m_out++ )
	    {
		Voxel_algorithm_3D () ;
	    }
	  Normalize () ;
	

      unsigned long int current ;
      for ( m_out = 0, current = m_CurrentIndex * m_NumberOfDirections ; m_out < m_NumberOfDirections ; m_out++, current++ )
	{
	  if ( ( m_prevValues[m_out] && ( fabs ( m_ProbMapArray[current] - m_prevValues[m_out] ) / m_prevValues[m_out] > 0.01 ) ) || ( !m_prevValues[m_out] && m_ProbMapArray[current] ) )
	    m_ProbMapChange++ ;
	}
	}
    }
}

//***************************************************************************************

void Fstar_prob3D::Fstar_Prob_map_algorithm_3D ()
{
  cout<<"Full algorithm starting"<<endl ;
  for( m_z_it = 0, m_CurrentIndex = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing() ;
	      m_CurrentIndex++ ;
	    }
	  for( m_x_it = m_XDimension - 1 ; m_x_it >= 0 ; m_x_it-- )
	    {
	      m_CurrentIndex-- ;
	      AlgorithmProcessing() ;
	    }
	  m_CurrentIndex += m_LINE ;
	}
      for( m_y_it = m_YDimension - 1 ; m_y_it >= 0 ; m_y_it-- )
	{
	  m_CurrentIndex -= m_LINE ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing() ;
	      m_CurrentIndex++ ;
	    }
	  for( m_x_it = m_XDimension - 1 ; m_x_it >= 0 ; m_x_it-- )
	    {
	      m_CurrentIndex-- ;
	      AlgorithmProcessing();
	    }
	}
      m_CurrentIndex += m_SLICE ;
    }
  cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_ProbMapChange<<endl;

  m_ProbMapChange = 0;
  m_NumberOfLoops++;
  for( m_z_it = m_ZDimension - 1 ; m_z_it >= 0 ; m_z_it-- )
    {
      m_CurrentIndex -= m_SLICE ;
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing();
	      m_CurrentIndex++ ;
	    }
	  for( m_x_it = m_XDimension - 1 ; m_x_it >= 0 ; m_x_it-- )
	    {
	      m_CurrentIndex-- ;
	      AlgorithmProcessing();
	    }
	  m_CurrentIndex += m_LINE ;
	}
      for( m_y_it = m_YDimension - 1 ; m_y_it >= 0 ; m_y_it-- )
	{
	  m_CurrentIndex -= m_LINE ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing();
	      m_CurrentIndex++ ;
	    }
	  for( m_x_it = m_XDimension - 1 ; m_x_it >= 0 ; m_x_it-- )
	    {
	      m_CurrentIndex -- ;
	      AlgorithmProcessing();
	    }
	}
    }
}

//***************************************************************************************
void Fstar_prob3D::ComputeFinalProbImage( std::string mincost_filename )
{
  ImageType::Pointer SumImage = ImageType::New() ;

  ImageType::IndexType startsum ;
  ImageType::SizeType  sizesum ;
  startsum[0] = startsum[1] = startsum[2] = 0 ;
  sizesum[0] = m_XDimension ;
  sizesum[1] = m_YDimension;
  sizesum[2] = m_ZDimension;
  ImageType::RegionType regionsum;
  regionsum.SetSize( sizesum );
  regionsum.SetIndex( startsum );

  SumImage->SetRegions( regionsum );
  SumImage->SetSpacing(m_Spacing);
  SumImage->Allocate();

  ProbType *SumImageArray = SumImage->GetPixelContainer()->GetBufferPointer();
  ProbType temp ;
  unsigned long int index ;

  for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      m_CurrentIndex = m_x_it + m_y_it*m_LINE + m_z_it*m_SLICE;
	      SumImageArray[m_CurrentIndex] = 0 ;
	      for ( m_out = 0 ; m_out < m_NumberOfDirections ; m_out++ )
		{
		  index = m_CurrentIndex*m_NumberOfDirections + m_out ;
		  temp = m_ProbMapArray[index] ;
		  SumImageArray[m_CurrentIndex] += temp ;
		}
	    }
	}
    }

  WriterType::Pointer Writer = WriterType::New();
  cout<<"Writing the final output image of the sum prob at each voxel: " << mincost_filename << endl;
  Writer->SetFileName ( mincost_filename.c_str() );
  Writer->SetInput( SumImage );
	
  try
    {
      Writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      cerr<<"Problem writing the output file"<<endl;
      cerr<<excp<<endl;
    }

}

//***************************************************************************************

void Fstar_prob3D::ComputeProbMap()
{
  time_t seconds_1 = time (NULL), seconds_2, diff;
  m_NumberOfLoops = 1;
  cout<<"Starting Initialisation."<<endl;

  Initialisation();

  cout<<"Starting the Algorithm."<<endl;
  m_ProbMapChange = 0;
  Fstar_Prob_map_algorithm_3D ();
  cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_ProbMapChange<<endl;

  m_NumberOfLoops++;
  do
    {
      do
	{
	  m_ProbMapChange = 0;
	  Fstar_Prob_map_algorithm_3D ();
	  cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_ProbMapChange<<endl;

	  m_NumberOfLoops++;
	}
      while( m_ProbMapChange > 2 );
		
      cout<<"Extra loop to check if there is really no change."<<endl;

      unsigned long int currentIndexBaseZ = 0, currentIndexBaseY = 0 ;
      for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
	{
	  currentIndexBaseZ = m_z_it * m_SLICE ;
	  for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	    {
	      currentIndexBaseY = currentIndexBaseZ + m_y_it*m_LINE ;
	      for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
		{
		  m_CurrentIndex = m_x_it + currentIndexBaseY ;
		}
	    }
	}
      m_ProbMapChange = 0;
      Fstar_Prob_map_algorithm_3D ();
      cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_ProbMapChange<<endl;

      m_NumberOfLoops++;
    }
  while( m_ProbMapChange > 2 );
	
  seconds_2 = time (NULL);
  diff = seconds_2 - seconds_1;
  cout<<"Number of seconds for execution: "<<diff<<endl;
  cout<<"Number of minutes for execution: "<<diff/60<<endl;
}

//***************************************************************************************

void Fstar_prob3D::GetProbMap(std::string mincost)
{
  ComputeProbMap();
  ComputeFinalProbImage( mincost ) ;
}

//***************************************************************************************

Fstar_prob3D::~Fstar_prob3D()
{
}

double Fstar_prob3D::GetAngleCoefficient(double AnglesSimilarityTemp, double mu, double sigma)
{
	//We apply an angle function that gives a coefficient between 0 (direction non similar) and 1 (similar)
	//This function is centered on mu and has a width of sigma around mu
	//The use of the error function makes it non linear
	double anglecoefficient;
	
	//tempresult = ( mu - AnglesSimilarityTemp )/sigma;
	anglecoefficient = 0.5*( 1 - erf( AnglesSimilarityTemp ));
	
	return anglecoefficient;
}

double Fstar_prob3D::InOutDirectionSimilarity ( unsigned short in, unsigned short out ) 
{
	//dot product is v1.v2 = |v1|*|v2|*cos(angle)
	VectorType v1, v2;
	v1 = m_SphereIkosahedronObject->GetCoordinateTableatIndex(in);
	v2 = m_SphereIkosahedronObject->GetCoordinateTableatIndex(out);

	double anglesimilarity;
	double x1, y1, z1, x2, y2, z2;
	x1 = v1[0];
	y1 = v1[1];
	z1 = v1[2];
	x2 = v2[0];
	y2 = v2[1];
	z2 = v2[2];
	
	//Norm of the vectors
	double v1_magnitude = sqrt(x1*x1 + y1*y1 + z1*z1);
	double v2_magnitude = sqrt(x2*x2 + y2*y2 + z2*z2);
	
	//Normalisation of the vectors
	x1 /= v1_magnitude;
	x2 /= v2_magnitude;
	y1 /= v1_magnitude;
	y2 /= v2_magnitude;
	z1 /= v1_magnitude;
	z2 /= v2_magnitude;
	
	double sum = x1*x2 + y1*y2 + z1*z2;
	if(sum <= -1)
	{
		sum = -1;
	}
	else if(sum >= 1)
	{
		sum = 1;
	}
	//anglesimilarity = acos(sum);
	anglesimilarity = - sum ;
	return GetAngleCoefficient ( anglesimilarity, cos ( M_PI/8 ), cos ( M_PI/16 ) ) ;
}

double Fstar_prob3D::GetToNeighbourDirectionSimilarity(short x, short y, short z, unsigned short direction_index)
{
	//dot product is v1.v2 = |v1|*|v2|*cos(angle)
	VectorType v1, v2;
	
	v1.push_back(x);
	v1.push_back(y);
	v1.push_back(z);
	v2 = m_SphereIkosahedronObject->GetCoordinateTableatIndex(direction_index);
	
	double anglesimilarity;
	double x1, y1, z1, x2, y2, z2;
	x1 = v1[0];
	y1 = v1[1];
	z1 = v1[2];
	x2 = v2[0];
	y2 = v2[1];
	z2 = v2[2];
	
	//Norm of the vectors
	double v1_magnitude = sqrt(x1*x1 + y1*y1 + z1*z1);
	double v2_magnitude = sqrt(x2*x2 + y2*y2 + z2*z2);
	
	//Normalisation of the vectors
	x1 /= v1_magnitude;
	x2 /= v2_magnitude;
	y1 /= v1_magnitude;
	y2 /= v2_magnitude;
	z1 /= v1_magnitude;
	z2 /= v2_magnitude;
	
	double sum = x1*x2 + y1*y2 + z1*z2;
	if(sum <= -1)
	{
		sum = -1;
	}
	else if(sum >= 1)
	{
		sum = 1;
	}
	//anglesimilarity = acos(sum);
	anglesimilarity = -sum ;
	return GetAngleCoefficient ( anglesimilarity, cos ( M_PI/8 ) , cos ( M_PI/16 ) ) ;
}

void Fstar_prob3D::PrecomputeProbs () 
{
  double length ;
  m_NeighborDirProbs.resize ( 3 ) ;
  long int xdir, ydir, zdir ;
  for ( short z = 0 ; z < 3 ; z++ )
    {
      zdir = z - 1 ;
      m_NeighborDirProbs[z].resize ( 3 ) ;
      for ( short y = 0 ; y < 3 ; y++ )
	{
	  ydir = y - 1 ;
	  m_NeighborDirProbs[z][y].resize ( 3 ) ;
	  for ( short x = 0 ; x < 3 ; x++ )
	    {
	      xdir = x - 1 ;
	      m_NeighborDirProbs[z][y][x].resize ( m_NumberOfDirections ) ;

	      length = abs ( xdir ) + abs ( ydir ) + abs ( zdir ) ;

	      if ( length == 1 ) length = 1 ;
	      else if ( length == 2 ) length = 1. / sqrt ( 2 ) ;
	      else if ( length == 3 ) length = 1. / sqrt ( 3 ) ;
	      else if ( length != 0 ) std::cout << "This is never supposed to happen!!!!" << std::endl ;

	      m_NeighborDist[z][y][x] = length ;
	      for ( unsigned short dir = 0 ; dir < m_NumberOfDirections ; dir++ )
		{
		  m_NeighborDirProbs[z][y][x][dir] = GetToNeighbourDirectionSimilarity ( xdir, ydir, zdir, dir ) ;
		}
	    }
	}
    }

  std::cout << InOutDirectionSimilarity (0, 0 ) << " should be one" << std::endl ;
  m_DirDirProbs.resize ( m_NumberOfDirections ) ;
  for ( unsigned long i = 0 ; i < m_NumberOfDirections ; i++ )
    {
      m_DirDirProbs[i].resize ( m_NumberOfDirections ) ;
      for ( unsigned long j = 0 ; j < m_NumberOfDirections ; j++ )
	{
	  m_DirDirProbs[i][j] = InOutDirectionSimilarity ( i, j ) ;
	  if ( m_DirDirProbs[i][j] < 0 ) { std::cout << "Dirdir" << std::endl ; exit ( 0 ) ;}
	  if ( m_DirDirProbs[i][j] > 1 ) { std::cout << "Dirdir+" << std::endl ; exit ( 0 ) ;}
	}
    }
}
