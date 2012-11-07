#include "Fstar_cost3D.h"
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

Fstar_cost3D::Fstar_cost3D(VectorImageType::Pointer Input_image , ImageType::Pointer InputLabelImage, unsigned long int subdivision_level, double thresholdCoefficient)
{
  m_ThresholdCoefficient = thresholdCoefficient ;
  m_XDimension = InputLabelImage->GetLargestPossibleRegion().GetSize()[0] ;
  m_YDimension = InputLabelImage->GetLargestPossibleRegion().GetSize()[1] ;
  m_ZDimension = InputLabelImage->GetLargestPossibleRegion().GetSize()[2] ;
  m_SLICE = m_XDimension*m_YDimension ;
  m_LINE = m_XDimension ;
  cout<<"xdim: "<<m_XDimension<<"   // ydim: "<<m_YDimension<<"   // zdim: "<<m_ZDimension<<endl ;
	
  //Icosahedron subdivision
  m_SphereIkosahedronObject = SphereIkosahedron<CostType>::New() ;
  m_SphereIkosahedronObject->SetSubdivisionLevel(subdivision_level) ;
  m_SphereIkosahedronObject->Initialisation_tables(true);
  m_NumberOfDirections = m_SphereIkosahedronObject->GetNumberOfVertices();
  cout<<"Number of vertices/directions: "<<m_NumberOfDirections<<endl;
  m_ZFactor = m_SLICE * m_NumberOfDirections ;
  m_YFactor = m_LINE * m_NumberOfDirections ;
  m_Spacing = Input_image->GetSpacing();
	
  //The image that will be modified/worked on
  m_CurrentImage = VectorImageType::New();
  m_CurrentLengthImage = VectorImageType::New();
  m_CurrentOrigImage = VectorImageType::New () ;
  m_CurrentAvgCostImage = VectorImageType::New () ;
  VectorImageType::IndexType start;
  itk::VariableLengthVector< CostType > f( m_NumberOfDirections );
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

  m_CurrentLengthImage->SetVectorLength ( m_NumberOfDirections ) ;
  m_CurrentLengthImage->SetRegions ( region ) ;
  m_CurrentLengthImage->Allocate () ;
  m_CurrentLengthImage->FillBuffer ( l ) ;

  m_CurrentOrigImage->SetVectorLength ( m_NumberOfDirections ) ;
  m_CurrentOrigImage->SetRegions ( region ) ;
  m_CurrentOrigImage->Allocate () ;
  m_CurrentOrigImage->FillBuffer ( l ) ;
			
  m_CurrentAvgCostImage->SetVectorLength ( m_NumberOfDirections ) ;
  m_CurrentAvgCostImage->SetRegions ( region ) ;
  m_CurrentAvgCostImage->Allocate () ;
  m_CurrentAvgCostImage->FillBuffer ( f ) ;

  m_LabelImageArray = InputLabelImage->GetPixelContainer()->GetBufferPointer();
  m_ODFArray = Input_image->GetPixelContainer()->GetBufferPointer();
  m_CostMapArray = m_CurrentImage->GetPixelContainer()->GetBufferPointer();
  m_LengthMapArray = m_CurrentLengthImage->GetPixelContainer()->GetBufferPointer();
  m_OrigMapArray = m_CurrentOrigImage->GetPixelContainer()->GetBufferPointer() ;
  m_AvgCostMapArray = m_CurrentAvgCostImage->GetPixelContainer()->GetBufferPointer() ;
  m_Epsilon = 0;
  m_ChangingFlag = true;

  PrecomputeCosts () ;
  FstarInitialization () ;
}

//***************************************************************************************

void Fstar_cost3D::Initialisation()
{
  unsigned long int sourceCounter = 0 ;
  for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      unsigned long int counter = 0;
	      m_ChangeFlag.push_back(-1);
	      m_CurrentIndex = m_x_it + m_y_it*m_LINE + m_z_it*m_SLICE;
	      if ( m_LabelImageArray[m_CurrentIndex] == SEED )
		{
		  sourceCounter++ ;
		}

	      for( m_DirectionIterator = 0 ; m_DirectionIterator < m_NumberOfDirections ; m_DirectionIterator++ )
		{
		  unsigned long int CurrentVoxelIndex = m_DirectionIterator + m_x_it*m_NumberOfDirections + m_y_it*m_YFactor + m_z_it*m_ZFactor;
					
		  if(m_ODFArray[CurrentVoxelIndex] > 0)
		    {
			//To check wether there is a propagation along that voxel, if not, no computation along it
		      counter++;
		    }
		  if(m_LabelImageArray[m_CurrentIndex] == SEED )
		    {
		    m_CostMapArray[CurrentVoxelIndex] = 0;//If voxel is labeled as source, set it to 0 along each direction
		    m_AvgCostMapArray[CurrentVoxelIndex] = 0;//If voxel is labeled as source, set it to 0 along each direction
		    m_LengthMapArray[CurrentVoxelIndex] = 0 ;
		    m_OrigMapArray[CurrentVoxelIndex] = sourceCounter ;
		    }
		  else
		    {
		    m_CostMapArray[CurrentVoxelIndex] = INF;
		    m_AvgCostMapArray[CurrentVoxelIndex] = INF;
		    m_LengthMapArray[CurrentVoxelIndex] = INF ;
		    m_OrigMapArray[CurrentVoxelIndex] = -1 ;
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

void Fstar_cost3D::FstarInitialization ()
{
  m_fstar_updates[0][0] = 13 ;
  m_fstar_updates[0][1] = 0 ;
  m_fstar_updates[0][2] = 1 ;
  m_fstar_updates[0][3] = 2 ;
  m_fstar_updates[0][4] = 3 ;
  m_fstar_updates[0][5] = 4 ;
  m_fstar_updates[0][6] = 5 ;
  m_fstar_updates[0][7] = 6 ;
  m_fstar_updates[0][8] = 7 ;
  m_fstar_updates[0][9] = 8 ;
  m_fstar_updates[0][10] = 9 ;
  m_fstar_updates[0][11] = 10 ;
  m_fstar_updates[0][12] = 11 ;
  m_fstar_updates[0][13] = 12 ;

  m_fstar_updates[1][0] = 1 ;
  m_fstar_updates[1][1] = 14 ;

  m_fstar_updates[2][0] = 4 ;
  m_fstar_updates[2][1] = 12 ;
  m_fstar_updates[2][2] = 15 ;
  m_fstar_updates[2][3] = 16 ;
  m_fstar_updates[2][4] = 17 ;

  m_fstar_updates[3][0] = 1 ;
  m_fstar_updates[3][1] = 14 ;

  m_fstar_updates[4][0] = 13 ;
  m_fstar_updates[4][1] = 14 ;
  m_fstar_updates[4][2] = 15 ;
  m_fstar_updates[4][3] = 16 ;
  m_fstar_updates[4][4] = 17 ;
  m_fstar_updates[4][5] = 18 ;
  m_fstar_updates[4][6] = 19 ;
  m_fstar_updates[4][7] = 20 ;
  m_fstar_updates[4][8] = 21 ;
  m_fstar_updates[4][9] = 22 ;
  m_fstar_updates[4][10] = 23 ;
  m_fstar_updates[4][11] = 24 ;
  m_fstar_updates[4][12] = 25 ;
  m_fstar_updates[4][13] = 26 ;

  m_fstar_updates[5][0] = 1 ;
  m_fstar_updates[5][1] = 14 ;

  m_fstar_updates[6][0] = 4 ;
  m_fstar_updates[6][1] = 12 ;
  m_fstar_updates[6][2] = 15 ;
  m_fstar_updates[6][3] = 16 ;
  m_fstar_updates[6][4] = 17 ;

  m_fstar_updates[7][0] = 1 ;
  m_fstar_updates[7][1] = 14 ;

  short int x, y, z, temp ;
  for ( int i = 0 ; i < 27 ; i++ )
    {
      x = i % 3 ;
      m_neighbors[i][2] = x ;

      temp = ( i - x ) / 3 ;
      y = temp % 3 ;
      m_neighbors[i][1] = y ;

      z = ( temp - y ) / 3 ;
      m_neighbors[i][0] = z ;
    }
}

void Fstar_cost3D::Voxel_algorithm_3D ()
{
  const unsigned long CurrentVoxelIndex = m_DirectionIterator + m_CurrentIndex * m_NumberOfDirections ;
  const double fODFCurrent = f_ODF ( m_ODFArray[CurrentVoxelIndex] ) ;

  long int neighbor, x_neighbour, y_neighbour, z_neighbour ;
  unsigned long int baseIndex ;
  CostType minAvgCost = INF, minRealAvgCost = INF, neighborCost, stepCost, minCost = INF;
  CostType currentCost, currentAvgCost, neighborAvgCost ;
  CostType minLength = -1, neighborLength = -1, stepLength, currentLength ;
  CostType origin = -1 ;

  // for each neighbor j
  unsigned long int nNeighbors = m_fstar_updates[m_fstar_direction][0] ;
  for ( unsigned long int j = 0 ; j < nNeighbors ; j++ )
    {
      neighbor = m_fstar_updates[m_fstar_direction][j+1] ;
      m_z_temp_it = m_neighbors[neighbor][0] ; 
      z_neighbour = m_z_it + m_z_temp_it - 1 ;
      if ( ( z_neighbour >= m_ZDimension ) || ( z_neighbour < 0 ) ) continue ;

      m_y_temp_it = m_neighbors[neighbor][1] ;
      y_neighbour = m_y_it + m_y_temp_it - 1 ;
      if ( ( y_neighbour >= m_YDimension ) || ( y_neighbour < 0 ) ) continue ;

      m_x_temp_it = m_neighbors[neighbor][2] ;
      x_neighbour = m_x_it + m_x_temp_it - 1 ;
      if ( ( x_neighbour >= m_XDimension ) || ( x_neighbour < 0 ) ) continue ;

      baseIndex = ( x_neighbour + y_neighbour*m_LINE + z_neighbour*m_SLICE ) * m_NumberOfDirections ;
      stepLength = m_NeighborDist[m_z_temp_it][m_y_temp_it][m_x_temp_it] ;

      // for each direction d
      for( unsigned long int d = 0 ; d < m_NumberOfDirections ; d++ )
	{
	  m_NeighbourIndex = d + baseIndex ;
	  neighborAvgCost = m_AvgCostMapArray[m_NeighbourIndex] ;
	  //if ( m_fstar_direction > 3 ) std::cout << "C" << std::endl ;
	  if ( neighborAvgCost == INF ) continue ;
	  neighborLength = m_LengthMapArray[m_NeighbourIndex] ; 
	  neighborCost = m_CostMapArray[m_NeighbourIndex] ;

	  // we consider whether it's cheaper to go to the current voxel from this neighbour
	  stepCost =  f_ODF ( m_ODFArray[m_NeighbourIndex] ) + fODFCurrent + m_DirDirCosts[m_DirectionIterator][d] + m_NeighborDirCosts[m_z_temp_it][m_y_temp_it][m_x_temp_it][m_DirectionIterator] ;
	  currentCost = neighborCost + stepCost * stepLength ;
	  currentLength = neighborLength + stepLength ;
	  currentAvgCost = currentCost / currentLength ;

	  if ( currentAvgCost < minRealAvgCost ) 
	    {
	      minRealAvgCost = currentAvgCost ;
	      if ( currentAvgCost < neighborAvgCost ) 
		{
		  currentAvgCost = neighborAvgCost ;
		  currentCost = currentAvgCost * currentLength ;
		}
	      origin = m_OrigMapArray[m_NeighbourIndex] ;
	      minLength = currentLength ;
	      minAvgCost = currentAvgCost ;
	      minCost = currentCost ;
	    }
	}
    }

  //check if there was a change or not
  if ( m_AvgCostMapArray[CurrentVoxelIndex] > minAvgCost ) 
    {
      m_CostMapArray[CurrentVoxelIndex] = minCost ;
      m_AvgCostMapArray[CurrentVoxelIndex] = minAvgCost ;
      m_LengthMapArray[CurrentVoxelIndex] = minLength  ;
      m_OrigMapArray[CurrentVoxelIndex] = origin ;
      m_CostMapChange++ ;
    }
  return ;
}

//***************************************************************************************

void Fstar_cost3D::ChangeFlagOptimization()
{
  //For faster convergence, if no change we upgrade the flag to a less likely computation stage (for 10, the voxel will be computed only every 10 iterations)
  if( m_NumberOfLoops > 2)
    {
      switch(m_ChangeFlag[m_CurrentIndex])
	{
	case -1:
	  m_ChangeFlag[m_CurrentIndex] = 1;
	  break;
	case 1:
	  m_ChangeFlag[m_CurrentIndex] = 2;
	  break;
	case 2:
	  m_ChangeFlag[m_CurrentIndex] = 5;
	  break;
	case 5:
	  m_ChangeFlag[m_CurrentIndex] = 10;
	  break;
	default:
	  break;
	}
    }
}
//***************************************************************************************

void Fstar_cost3D::AlgorithmProcessing()
{
  if( m_PropagationFlag[m_CurrentIndex] == 1 )
    {
      m_CurrentCost = m_LabelImageArray[m_CurrentIndex] ;

      if( m_CurrentCost != SEED )
	{
	  for( m_DirectionIterator = 0 ; m_DirectionIterator < m_NumberOfDirections ; m_DirectionIterator++ )
	    {
		Voxel_algorithm_3D () ;
	    }
	}
    }
}

//***************************************************************************************

void Fstar_cost3D::Fstar_Cost_map_algorithm_3D ()
{
  cout<<"Full algorithm starting"<<endl ;
  for( m_z_it = 0, m_CurrentIndex = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  m_fstar_direction = 0 ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing() ;
	      m_CurrentIndex++ ;
	    }
	  m_fstar_direction = 1 ;
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
	  m_fstar_direction = 2 ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing() ;
	      m_CurrentIndex++ ;
	    }
	  m_fstar_direction = 3 ;
	  for( m_x_it = m_XDimension - 1 ; m_x_it >= 0 ; m_x_it-- )
	    {
	      m_CurrentIndex-- ;
	      AlgorithmProcessing();
	    }
	}
      m_CurrentIndex += m_SLICE ;
    }
  cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_CostMapChange<<endl;

  m_CostMapChange = 0;
  m_NumberOfLoops++;
  for( m_z_it = m_ZDimension - 1 ; m_z_it >= 0 ; m_z_it-- )
    {
      m_CurrentIndex -= m_SLICE ;
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  m_fstar_direction = 4 ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing();
	      m_CurrentIndex++ ;
	    }
	  m_fstar_direction = 5 ;
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
	  m_fstar_direction = 6 ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      AlgorithmProcessing();
	      m_CurrentIndex++ ;
	    }
	  m_fstar_direction = 7 ;
	  for( m_x_it = m_XDimension - 1 ; m_x_it >= 0 ; m_x_it-- )
	    {
	      m_CurrentIndex -- ;
	      AlgorithmProcessing();
	    }
	}
    }
}

//***************************************************************************************
void Fstar_cost3D::ComputeFinalMinImage( std::string mincost_filename, std::string length_filename, std::string orig_filename, std::string avgcostfilename, std::string mincostperminlengthfilename)
{
  ImageType::Pointer SumImage = ImageType::New() ;
  ImageType::Pointer AvgCostImage = ImageType::New() ;
  LengthImageType::Pointer LengthImage = LengthImageType::New() ;
  LengthImageType::Pointer OrigImage = LengthImageType::New() ;

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
  LengthImage->SetRegions ( regionsum ) ;
  LengthImage->SetSpacing ( m_Spacing ) ;
  LengthImage->Allocate () ;
  OrigImage->SetRegions ( regionsum ) ;
  OrigImage->SetSpacing ( m_Spacing ) ;
  OrigImage->Allocate () ;
  AvgCostImage->SetRegions ( regionsum ) ;
  AvgCostImage->SetSpacing ( m_Spacing ) ;
  AvgCostImage->Allocate () ;

  CostType *SumImageArray = SumImage->GetPixelContainer()->GetBufferPointer();
  CostType temp ;
  LengthType *LengthImageArray = LengthImage->GetPixelContainer()->GetBufferPointer();
  unsigned long int index ;
  LengthType *OrigImageArray = OrigImage->GetPixelContainer()->GetBufferPointer() ;
  CostType *AvgCostImageArray = AvgCostImage->GetPixelContainer()->GetBufferPointer () ;

  for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      m_CurrentIndex = m_x_it + m_y_it*m_LINE + m_z_it*m_SLICE;
	      SumImageArray[m_CurrentIndex] = INF ;
	      LengthImageArray[m_CurrentIndex] = INF ;
	      OrigImageArray[m_CurrentIndex] = INF ;
	      AvgCostImageArray[m_CurrentIndex] = INF ;
	      for ( m_DirectionIterator = 0 ; m_DirectionIterator < m_NumberOfDirections ; m_DirectionIterator++ )
		{
		  index = m_CurrentIndex*m_NumberOfDirections + m_DirectionIterator ;
		  temp = m_AvgCostMapArray[index] ;
		  if ( AvgCostImageArray[m_CurrentIndex] > temp )
		    {
		      SumImageArray[m_CurrentIndex] = m_CostMapArray[index] ;
		      LengthImageArray[m_CurrentIndex] = (LengthType) m_LengthMapArray[index] ;
		      OrigImageArray[m_CurrentIndex] = (LengthType) m_OrigMapArray[index] ;
		      AvgCostImageArray[m_CurrentIndex] = temp ;
		    }
		}
	      if ( AvgCostImageArray[m_CurrentIndex] == INF )
		{
		  SumImageArray[m_CurrentIndex] = -1 ;
		  LengthImageArray[m_CurrentIndex] = -1 ;
		  OrigImageArray[m_CurrentIndex] = -1 ;
		  AvgCostImageArray[m_CurrentIndex] = -1 ;
		}

	      // make length + 1 to be able to visualize it as labels
	      LengthImageArray[m_CurrentIndex]++ ;
	      OrigImageArray[m_CurrentIndex]++ ;
	    }
	}
    }

  WriterType::Pointer Writer = WriterType::New();
  cout<<"Writing the final output image of the min cost at each voxel: " << mincost_filename << endl;
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

  LengthWriterType::Pointer LengthWriter = LengthWriterType::New();
  cout<<"Writing the final output image of the length to each voxel: " << length_filename << endl;
  LengthWriter->SetFileName ( length_filename.c_str() );
  LengthWriter->SetInput( LengthImage );
	
  try
    {
      LengthWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      cerr<<"Problem writing the output file"<<endl;
      cerr<<excp<<endl;
    }

  LengthWriterType::Pointer OrigWriter = LengthWriterType::New();
  cout<<"Writing the final output image of the origin that led to each voxel: " << orig_filename << endl;
  OrigWriter->SetFileName ( orig_filename.c_str() );
  OrigWriter->SetInput( OrigImage );
	
  try
    {
      OrigWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      cerr<<"Problem writing the output file"<<endl;
      cerr<<excp<<endl;
    }

  WriterType::Pointer AvgCostWriter = WriterType::New();
  cout<<"Writing the final output image of the avg cost to each voxel: " << avgcostfilename << endl;
  AvgCostWriter->SetFileName ( avgcostfilename.c_str() );
  AvgCostWriter->SetInput( AvgCostImage );
	
  try
    {
      AvgCostWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      cerr<<"Problem writing the output file"<<endl;
      cerr<<excp<<endl;
    }

}

//***************************************************************************************

void Fstar_cost3D::ComputeCostMap()
{
  time_t seconds_1 = time (NULL), seconds_2, diff;
  m_NumberOfLoops = 1;
  cout<<"Starting Initialisation."<<endl;

  Initialisation();

  cout<<"Starting the Algorithm."<<endl;
  m_CostMapChange = 0;
  Fstar_Cost_map_algorithm_3D ();
  cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_CostMapChange<<endl;

  m_NumberOfLoops++;
  do
    {
      do
	{
	  m_CostMapChange = 0;
	  Fstar_Cost_map_algorithm_3D ();
	  cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_CostMapChange<<endl;

	  m_NumberOfLoops++;
	}
      while( m_CostMapChange > 2 );
		
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
      m_CostMapChange = 0;
      Fstar_Cost_map_algorithm_3D ();
      cout<<"loop done: #"<<m_NumberOfLoops<<" - changes: "<<m_CostMapChange<<endl;

      m_NumberOfLoops++;
    }
  while( m_CostMapChange > 2 );
	
  seconds_2 = time (NULL);
  diff = seconds_2 - seconds_1;
  cout<<"Number of seconds for execution: "<<diff<<endl;
  cout<<"Number of minutes for execution: "<<diff/60<<endl;
}

//***************************************************************************************

void Fstar_cost3D::GetCostMap(std::string mincost, std::string lengthmap, std::string originmap, std::string avgcost, std::string mincostperminlength)
{
  ComputeCostMap();
  ComputeFinalMinImage( mincost, lengthmap, originmap, avgcost, mincostperminlength ) ;
}

//***************************************************************************************

Fstar_cost3D::~Fstar_cost3D()
{
}

double Fstar_cost3D::f_ODF ( double odf ) 
{
  //if ( odf > 0.01 )  std::cout << odf << std::endl ;
  //return m_alpha * ( 1. - odf ) ;
  return 1. - odf ;
}

double Fstar_cost3D::GetAngleCoefficient(double AnglesSimilarityTemp, double mu, double sigma)
{
	//We apply an angle function that gives a coefficient between 0 (direction non similar) and 1 (similar)
	//This function is centered on mu and has a width of sigma around mu
	//The use of the error function makes it non linear
	double anglecoefficient;
	
	//tempresult = ( mu - AnglesSimilarityTemp )/sigma;
	anglecoefficient = 0.5*( 1 + erf( AnglesSimilarityTemp ));
	
	return anglecoefficient;
}

double Fstar_cost3D::InOutDirectionSimilarity ( unsigned short in, unsigned short out ) 
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

void Fstar_cost3D::CheckSanity () 
{
  double dirdir, xyzdir ;
  double x, y, z, vmag ;
  for ( unsigned short i = 0 ; i < m_NumberOfDirections ; i++ )
    {
      for ( unsigned short j = 0 ; j < m_NumberOfDirections ; j++ )
	{
	  dirdir = m_DirDirCosts[i][j] ;
	  VectorType v = m_SphereIkosahedronObject->GetCoordinateTableatIndex(i) ;
	  x = v[0] ; 
	  y = v[1] ;
	  z = v[2] ;
	  xyzdir = GetToNeighbourDirectionSimilarity ( x, y, z, j ) ;

	  if ( dirdir != xyzdir ) std::cout << dirdir << " " << xyzdir << std::endl ;
	}
    }
  exit ( 0 ) ;
}

double Fstar_cost3D::GetToNeighbourDirectionSimilarity(double x1, double y1, double z1, unsigned short direction_index)
{
	//dot product is v1.v2 = |v1|*|v2|*cos(angle)
	VectorType  v2;
	
	v2 = m_SphereIkosahedronObject->GetCoordinateTableatIndex(direction_index);
	
	double anglesimilarity;
	double  x2, y2, z2;
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

void Fstar_cost3D::PrecomputeCosts () 
{
  double length ;
  m_NeighborDirCosts.resize ( 3 ) ;
  long int xdir, ydir, zdir ;
  for ( short z = 0 ; z < 3 ; z++ )
    {
      zdir = z - 1 ;
      m_NeighborDirCosts[z].resize ( 3 ) ;
      for ( short y = 0 ; y < 3 ; y++ )
	{
	  ydir = y - 1 ;
	  m_NeighborDirCosts[z][y].resize ( 3 ) ;
	  for ( short x = 0 ; x < 3 ; x++ )
	    {
	      xdir = x - 1 ;
	      m_NeighborDirCosts[z][y][x].resize ( m_NumberOfDirections ) ;

	      length = abs ( xdir ) + abs ( ydir ) + abs ( zdir ) ;

	      if ( length == 1 ) length = 1 ;
	      else if ( length == 2 ) length = sqrt ( 2 ) ;
	      else if ( length == 3 ) length = sqrt ( 3 ) ;
	      else if ( length != 0 ) std::cout << "This is never supposed to happen!!!!" << std::endl ;

	      m_NeighborDist[z][y][x] = length ;
	      for ( unsigned short dir = 0 ; dir < m_NumberOfDirections ; dir++ )
		{
		  m_NeighborDirCosts[z][y][x][dir] = GetToNeighbourDirectionSimilarity ( xdir, ydir, zdir, dir ) ;
		}
	    }
	}
    }

  std::cout << InOutDirectionSimilarity (0, 0 ) << " should be zero" << std::endl ;
  m_DirDirCosts.resize ( m_NumberOfDirections ) ;
  for ( unsigned long i = 0 ; i < m_NumberOfDirections ; i++ )
    {
      m_DirDirCosts[i].resize ( m_NumberOfDirections ) ;
      for ( unsigned long j = 0 ; j < m_NumberOfDirections ; j++ )
	{
	  m_DirDirCosts[i][j] = InOutDirectionSimilarity ( i, j ) ;
	  if ( m_DirDirCosts[i][j] < 0 ) { std::cout << "Dirdir" << std::endl ; exit ( 0 ) ;}
	  if ( m_DirDirCosts[i][j] > 1 ) { std::cout << "Dirdir+" << std::endl ; exit ( 0 ) ;}
	}
    }

  //CheckSanity () ;
}
