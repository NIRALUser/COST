#include "ODFCost.h"

static unsigned int iterationNum = 0; 

// initialize internal variables
ODFCost::ODFCost ( ODFCost::ODFImageType::Pointer odf, DoubleImageType::Pointer origfa, unsigned long verts ,double alpha)
{ 
   this->m_ODFImage = odf ;

   //read FA file
   this->m_OrigFAImage = origfa;

   this->m_alpha = alpha;
   this->m_NumberOfDirs = verts ;
   this->m_NumberOfSourceVoxels = 0 ;
   this->m_ODFArray = odf->GetPixelContainer()->GetBufferPointer();

   
   this->m_OrigFAArray = origfa->GetPixelContainer()->GetBufferPointer();

   this->PrecomputeStepLength () ;
   this->FstarTraversalUpdates () ;
   std::cout << "initialization done" << std::endl ;
}      

// main function - public
ODFCost::DoubleImageType::Pointer ODFCost::Cost ( ODFCost::LabelImageType::Pointer source ) 
{
  this->m_Converged = false ;
  this->m_NChanged = 0 ;
  this->m_SourceArray = source->GetPixelContainer()->GetBufferPointer();

  this->CreateCostMaps () ;
  
  this->FindSourceVoxels () ;
  
  this->InitializeSourceVoxels () ;
 
  this->NormalizeODF () ;
  
  this->CreateFinslerTable();

  this->PrecomputeDirection();

  std::cout << "Entering main loop" << std::endl ;
  unsigned i = 0 ;

  while ( !this->m_Converged ) 
    {
      iterationNum ++ ;
      this->MakeOnePassFstar () ;
      i++ ;
      std::cout << "Pass " << i << " complete. " << this->m_NChanged << " voxels updated in this pass." << std::endl ;
      // i get impatient when debugging - but this needs to be taken out so that the algorithm runs to convergence
      
      if(this->m_NChanged == 0)
	{
	  this->m_Converged = true;
	}
      else
	{
	  this->m_Converged = false;
	}    
      //---------------------------------------
      this->m_NChanged = 0 ;
    }

  std::cout << "Compute final costs" << std::endl ;
  this->ComputeFinalCosts () ;

  return this->m_FinalCostImage ;
}

// visit every voxel in the image
void ODFCost::MakeOnePass ()
{
  
  unsigned long int x, y, z ;

  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	}
   }
}


// visit every voxel according to F* traversal algorithm
void ODFCost::MakeOnePassFstar ()
{
  short x, y, z ;

  for( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for( y = 0 ; y < this->m_YDim ; y++ )
	{
	  this->m_fstar_direction = 0 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	  this->m_fstar_direction = 1 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	}
      for( y = this->m_YDim - 1 ; y >= 0 ; y-- )
	{
	  this->m_fstar_direction = 2 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	  this->m_fstar_direction = 3 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	}
    }

  std::cout << "Direction flip. " << this->m_NChanged << " voxels updated so far in this pass." << std::endl ;
  for( z = this->m_ZDim - 1 ; z >= 0 ; z-- )
    {
      for( y = 0 ; y < this->m_YDim ; y++ )
	{
	  this->m_fstar_direction = 4 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	  this->m_fstar_direction = 5 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	}
      for( y = this->m_YDim - 1 ; y >= 0 ; y-- )
	{
	  this->m_fstar_direction = 6 ;
	  for( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	  this->m_fstar_direction = 7 ;
	  for( x = this->m_XDim - 1 ; x >= 0 ; x-- )
	    {
	      this->VisitVoxel ( x, y, z ) ;
	    }
	}
    }
}

// locate all voxels that are in the source region, store their coordinates in a table
void ODFCost::FindSourceVoxels ()
{
  unsigned long index ;
  ODFCost::CoordinateType currentVoxel ;
  this->m_NumberOfSourceVoxels = 0 ;
  this->m_SourceVoxels.clear () ;

  for ( unsigned long z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( unsigned long y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( unsigned long x = 0 ; x < this->m_XDim ; x++ )
	    {
	       
	      index = x + y * this->m_XDim + z * this->m_Slice ;
	      if ( this->m_SourceArray[index] )
		{
		  currentVoxel[0] = x ;
		  currentVoxel[1] = y ;
		  currentVoxel[2] = z ;
		  this->m_SourceVoxels.push_back ( currentVoxel ) ;
		  this->m_NumberOfSourceVoxels ++ ;
		}
	    }
	}
    }
}

// allocate and initialize cost images
void ODFCost::CreateCostMaps () 
{
   //create empty output image
   this->m_CostImage = ODFImageType::New () ;
   this->m_AvgCostImage = ODFImageType::New () ;
   this->m_FinalCostImage = DoubleImageType::New () ;
   this->m_LengthImage = ODFImageType::New () ;
   this->m_OrigImage = ODFImageType::New () ;
   this->m_FinalOrigImage = DoubleImageType::New () ;
   this->m_FinalLengthImage = DoubleImageType::New () ;

   //added by Wenyu

   this->m_OrigDirImage = DoubleImageType::New () ;
   this->m_OrigNeighborImage = DoubleImageType::New () ;
   this->m_OrigInDirImage = DoubleImageType::New () ;
   this->m_OrigOutDirImage = DoubleImageType::New () ;

//    this->m_FinslerImageX = DoubleImageType::New();
//    this->m_FinslerImageY = DoubleImageType::New();
//    this->m_FinslerImageZ = DoubleImageType::New();

   itk::VariableLengthVector< ODFType > f ( this->m_NumberOfDirs ) ;
   for ( unsigned long i = 0 ; i < this->m_NumberOfDirs ; i++ ) 
     { 
       f[i] = INF; 
     }

   ODFImageType::IndexType start;
   ODFImageType::SizeType  size;
   start[0] = 0;   start[1] = 0;   start[2] = 0;
   size = this->m_ODFImage->GetLargestPossibleRegion().GetSize() ;

   this->m_XDim = size[0] ;
   this->m_YDim = size[1] ;
   this->m_ZDim = size[2] ;
   this->m_Slice = this->m_YDim * this->m_XDim ;

   ODFImageType::SpacingType spacing;
   spacing = this->m_ODFImage->GetSpacing () ;
   this->m_CostImage->SetSpacing ( spacing ) ;
   this->m_FinalCostImage->SetSpacing ( spacing ) ;
   this->m_LengthImage->SetSpacing ( spacing ) ;
   this->m_AvgCostImage->SetSpacing ( spacing ) ;
   this->m_OrigImage->SetSpacing ( spacing ) ;
   this->m_FinalOrigImage->SetSpacing ( spacing ) ;
   this->m_FinalLengthImage->SetSpacing ( spacing ) ;

   //added by Wenyu
   this->m_OrigDirImage->SetSpacing ( spacing ) ;
   this->m_OrigNeighborImage->SetSpacing ( spacing ) ;
   this->m_OrigInDirImage->SetSpacing ( spacing ) ;
   this->m_OrigOutDirImage->SetSpacing ( spacing ) ;

//    this->m_FinslerImageX->SetSpacing( spacing );
//    this->m_FinslerImageY->SetSpacing( spacing );
//    this->m_FinslerImageZ->SetSpacing( spacing );

   ODFImageType::RegionType region;
   region.SetSize ( size ) ;
   region.SetIndex ( start ) ;
   this->m_CostImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_CostImage->SetRegions ( region ) ;
   this->m_CostImage->Allocate () ;
   this->m_CostImage->FillBuffer ( f ) ;

   this->m_AvgCostImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_AvgCostImage->SetRegions ( region ) ;
   this->m_AvgCostImage->Allocate () ;
   this->m_AvgCostImage->FillBuffer ( f ) ;

   this->m_LengthImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_LengthImage->SetRegions ( region ) ;
   this->m_LengthImage->Allocate () ;
   this->m_LengthImage->FillBuffer ( f ) ;

   this->m_OrigImage->SetVectorLength ( this->m_NumberOfDirs ) ;
   this->m_OrigImage->SetRegions ( region ) ;
   this->m_OrigImage->Allocate () ;
   this->m_OrigImage->FillBuffer ( f ) ;

   this->m_FinalCostImage->SetRegions ( region ) ;
   this->m_FinalCostImage->Allocate () ;
   this->m_FinalCostImage->FillBuffer ( 0 ) ;

   this->m_FinalOrigImage->SetRegions ( region ) ;
   this->m_FinalOrigImage->Allocate () ;
   this->m_FinalOrigImage->FillBuffer ( 0 ) ;

   this->m_FinalLengthImage->SetRegions ( region ) ;
   this->m_FinalLengthImage->Allocate () ;
   this->m_FinalLengthImage->FillBuffer ( 0 ) ;

   //added by Wenyu

   this->m_OrigDirImage->SetRegions ( region ) ;
   this->m_OrigDirImage->Allocate () ;
   this->m_OrigDirImage->FillBuffer ( 0 ) ;

   this->m_OrigNeighborImage->SetRegions ( region ) ;
   this->m_OrigNeighborImage->Allocate () ;
   this->m_OrigNeighborImage->FillBuffer ( 0 ) ;


   this->m_OrigInDirImage->SetRegions ( region ) ;
   this->m_OrigInDirImage->Allocate () ;
   this->m_OrigInDirImage->FillBuffer ( 0 ) ;


   this->m_OrigOutDirImage->SetRegions ( region ) ;
   this->m_OrigOutDirImage->Allocate () ;
   this->m_OrigOutDirImage->FillBuffer ( 0 ) ;

//    this->m_FinslerImageX->SetRegions ( region ) ;
//    this->m_FinslerImageX->Allocate () ;
//    this->m_FinslerImageX->FillBuffer ( 0 ) ;
//    this->m_FinslerImageY->SetRegions ( region ) ;
//    this->m_FinslerImageY->Allocate () ;
//    this->m_FinslerImageY->FillBuffer ( 0 ) ;
//    this->m_FinslerImageZ->SetRegions ( region ) ;
//    this->m_FinslerImageZ->Allocate () ;
//    this->m_FinslerImageZ->FillBuffer ( 0 ) ;

   this->m_CostArray = this->m_CostImage->GetPixelContainer()->GetBufferPointer();
   this->m_AvgCostArray = this->m_AvgCostImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalCostArray = this->m_FinalCostImage->GetPixelContainer()->GetBufferPointer() ;
   this->m_LengthArray = this->m_LengthImage->GetPixelContainer()->GetBufferPointer();
   this->m_OrigArray = this->m_OrigImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalOrigArray = this->m_FinalOrigImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalLengthArray = this->m_FinalLengthImage->GetPixelContainer()->GetBufferPointer();

 //added by Wenyu
   this->m_OrigDir = this->m_OrigDirImage->GetPixelContainer()->GetBufferPointer();
   this->m_OrigNeighbor = this->m_OrigNeighborImage->GetPixelContainer()->GetBufferPointer();
   this->m_OrigInDir = this->m_OrigInDirImage->GetPixelContainer()->GetBufferPointer();
   this->m_OrigOutDir = this->m_OrigOutDirImage->GetPixelContainer()->GetBufferPointer();

//    this->m_FinslerArrayX = this->m_FinslerImageX->GetPixelContainer()->GetBufferPointer();
//    this->m_FinslerArrayY = this->m_FinslerImageY->GetPixelContainer()->GetBufferPointer();
//    this->m_FinslerArrayZ = this->m_FinslerImageZ->GetPixelContainer()->GetBufferPointer();
   this->m_Cost = this->m_CostImage->GetPixelContainer()->GetBufferPointer();
}

// source voxels have 0 cost
void ODFCost::InitializeSourceVoxels () 
{
  ODFCost::CoordinateType currentVoxel ;
  unsigned long index, indexTimesNDirs ;
  unsigned int counter = 1 ;

  
  unsigned long int x, y, z;
  index = 0, indexTimesNDirs = 0 ;
  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      index =  x + y * this->m_XDim + z * this->m_Slice ;

	      
	      
	      this->m_OrigDir[index] = this->m_NumberOfDirs;
	      this->m_OrigNeighbor[index] = 13;
	      this->m_OrigInDir[index] = this->m_NumberOfDirs;
	      this->m_OrigOutDir[index] = this->m_NumberOfDirs;
	      
	      indexTimesNDirs = index * this->m_NumberOfDirs ;
	      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
		{
		  this->m_CostArray[indexTimesNDirs + dir] = INF;
		  this->m_AvgCostArray[indexTimesNDirs + dir] = INF;
		  this->m_LengthArray[indexTimesNDirs+dir] = 0   ;
		  this->m_OrigArray[indexTimesNDirs+dir] = this->m_NumberOfDirs ; 		
		}
	      
	      this->m_FinalCostArray[index] = INF ;
	    }
	}
    } 
  
  for ( unsigned long i = 0 ; i < this->m_NumberOfSourceVoxels ; i++ )
    {
      currentVoxel = this->m_SourceVoxels[i] ;
      
      index = currentVoxel[0] + (currentVoxel[1]) * this->m_XDim + (currentVoxel[2]) * this->m_Slice ;
      indexTimesNDirs = index * this->m_NumberOfDirs ;
      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir ++ )
	{
	  this->m_CostArray[indexTimesNDirs+dir] = 0 ;
	  this->m_AvgCostArray[indexTimesNDirs+dir] = 0 ;
	  this->m_LengthArray[indexTimesNDirs+dir] = 0 ;
	  this->m_OrigArray[indexTimesNDirs+dir] = counter ;

	}
      this->m_FinalCostArray[index] = 0 ;
      counter++ ;
    }
  
  
}

// visit voxel, updating cost
void ODFCost::VisitVoxel ( unsigned long x, unsigned long y, unsigned long z )
{
  
  unsigned long voxelIndex = x + y * this->m_XDim + z * this->m_Slice ;
  unsigned long voxelIndexTimesNDirs = voxelIndex * this->m_NumberOfDirs ;
  //if ( this->m_FAArray[voxelIndex] < 0.4 || this->m_SourceArray[voxelIndex] )
  if ( this->m_OrigFAArray[voxelIndex] < 0.05 ||  this->m_SourceArray[voxelIndex] )
    {
      // not worth the computation time TODO: add as command line argument
      return ;
    }

  unsigned long neighborIndex, neighborIndexTimesNDirs ;
  unsigned short neighborIterator ;
  bool validNeighbor ;
  double newStepCost ;
  ODFCost::CoordinateType idir, jdir ,proj, ndir, projn;
  double dn;
  double ndirCost;
  unsigned int origInDir = this->m_OrigInDir[voxelIndex];
  unsigned int origOutDir = this->m_OrigOutDir[voxelIndex];
  unsigned int origNeighbor = this->m_OrigNeighbor[voxelIndex];
  
    // instead of all 26 neighbors, only visit the neighbors as indicated by fstar for our current traversal direction
  for ( unsigned short n = 0 ; n < this->m_fstar_updates[this->m_fstar_direction][0] ; n++ )
    {
     
      neighborIterator = this->m_fstar_updates[this->m_fstar_direction][n+1] ;
      neighborIndex = this->GetNeighbor ( x, y, z, neighborIterator, validNeighbor ) ;
      if ( ! validNeighbor ) continue ;

      neighborIndexTimesNDirs = neighborIndex * this->m_NumberOfDirs ;
      double dn = this->m_DNTable[neighborIterator] ;
      
      // for each incoming dir
      for ( unsigned int inDir = 0 ; inDir < this->m_NumberOfDirs ; inDir++ )
     
	{
	  // for each outgoing dir
	   for ( unsigned int outDir = 0 ; outDir < this->m_NumberOfDirs ; outDir++ )
	    {
	      //I use this constraints for mouse data, it's unneccessary for synthetic data
	      //if(this->m_DirectionTable[inDir][origInDir] && this->m_DirectionTable[outDir][origOutDir])
	     	{
		  newStepCost =  dn*this->stepCost ( inDir, outDir, neighborIterator, voxelIndexTimesNDirs, neighborIndexTimesNDirs, neighborIndex );
		  double f_odf1 = this->m_FinslerTable[ voxelIndexTimesNDirs + inDir];
		  double f_odf2 = this->m_FinslerTable[ neighborIndexTimesNDirs + outDir ];
		  
		  this->UpdateCost ( voxelIndexTimesNDirs+inDir, voxelIndex, newStepCost + this->m_CostArray[neighborIndexTimesNDirs+outDir], this->m_AvgCostArray[neighborIndexTimesNDirs+outDir], dn+this->m_LengthArray[neighborIndexTimesNDirs+outDir], this->m_OrigArray[neighborIndexTimesNDirs+outDir], inDir, outDir, neighborIterator ) ;
		   
		}
	    }
	}
    }
}

//for 5x5x5 neighbors
/*
void ODFCost::PrecomputeDirection()
{
  std::cout<<"Compute direction table"<<std::endl;
  this->m_DirectionTable.resize ( this->m_NumberOfDirs + 1) ;
  for(unsigned int i = 0; i < this->m_NumberOfDirs ; i ++)
    {
      this->m_DirectionTable[i].resize ( this->m_NumberOfDirs + 1) ;
      for(unsigned int j=0; j < this->m_NumberOfDirs; j++)
	{
	  ODFCost::CoordinateType jdir = this->m_CoordinateTable[j];
	  ODFCost::CoordinateType idir = this->m_CoordinateTable[i];  
	  double magi = sqrt(idir[0]*idir[0]+idir[1]*idir[1]+idir[2]*idir[2]);
	  double magj = sqrt(jdir[0]*jdir[0]+jdir[1]*jdir[1]+jdir[2]*jdir[2]);
	  double sum = idir[0]*jdir[0]+idir[1]*jdir[1]+idir[2]*jdir[2];
	  double angle = acos(sum/(magi*magj));
	  if( angle < 0.5*PI )
	    this->m_DirectionTable[i][j] = 1;
	  else
	    this->m_DirectionTable[i][j] = 0; 
	}
      this->m_DirectionTable[i][this->m_NumberOfDirs] = 0;
    }
  this->m_DirectionTable[this->m_NumberOfDirs].resize ( this->m_NumberOfDirs + 1) ;
  for(unsigned int k = 0; k <= this->m_NumberOfDirs; k ++ )
    {
      this->m_DirectionTable[this->m_NumberOfDirs][k] = 0;
    }
  
  this->m_NeighborDirectionTable.resize ( this->m_NumberOfDirs + 1 ) ;
  for ( unsigned i = 0 ; i < this->m_NumberOfDirs ; i++ )
    {
      
      this->m_NeighborDirectionTable[i].resize ( 124 ) ;
      //std::cout<<"phase 2"<<std::endl;
      for ( unsigned j = 0 ; j < 124 ; j++ )
	{
	  ODFCost::CoordinateType idir, jdir ;
	  idir = this->m_CoordinateTable[i] ;
	  
	  short nx, ny, nz ;
	  nx = j % 5  ;
	  short neighborTemp = ( j - nx ) / 5;
	  nx = nx - 2 ;
	  ny = neighborTemp % 5 ;
	  nz = ( neighborTemp - ny ) / 5 - 2 ;
	  ny = ny - 2 ;
  
	  jdir[0] = nx / this->m_DNTable[j] ;
	  jdir[1] = ny / this->m_DNTable[j] ;
	  jdir[2] = nz / this->m_DNTable[j] ;
  
	  double magi = sqrt(idir[0]*idir[0]+idir[1]*idir[1]+idir[2]*idir[2]);
	  double magj = sqrt(jdir[0]*jdir[0]+jdir[1]*jdir[1]+jdir[2]*jdir[2]);
	  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
	  double angle = acos(sum/(magi*magj));
	  //std::cout<<"phase 3 "<<std::endl;
	  if(angle < 0.5*PI)
	    this->m_NeighborDirectionTable[i][j] = 1 ;
	  else
	    this->m_NeighborDirectionTable[i][j] = 0 ;
	}
    }
  this->m_NeighborDirectionTable[this->m_NumberOfDirs].resize ( 124 ) ;
  for( unsigned int p = 0; p < 124; p++)
    {
      this->m_NeighborDirectionTable[this->m_NumberOfDirs][p] = 0;
    }
}
*/

void ODFCost::PrecomputeDirection()
{
  std::cout<<"Compute direction table"<<std::endl;
  this->m_DirectionTable.resize ( this->m_NumberOfDirs + 1) ;
  for(unsigned int i = 0; i < this->m_NumberOfDirs ; i ++)
    {
      this->m_DirectionTable[i].resize ( this->m_NumberOfDirs + 1) ;
      for(unsigned int j=0; j < this->m_NumberOfDirs; j++)
	{
	  ODFCost::CoordinateType jdir = this->m_CoordinateTable[j];
	  ODFCost::CoordinateType idir = this->m_CoordinateTable[i];  
	  double magi = sqrt(idir[0]*idir[0]+idir[1]*idir[1]+idir[2]*idir[2]);
	  double magj = sqrt(jdir[0]*jdir[0]+jdir[1]*jdir[1]+jdir[2]*jdir[2]);
	  double sum = idir[0]*jdir[0]+idir[1]*jdir[1]+idir[2]*jdir[2];
	  double angle = acos(sum/(magi*magj));
	  if( angle < 0.5*PI )
	    this->m_DirectionTable[i][j] = 1;
	  else
	    this->m_DirectionTable[i][j] = 0; 
	}
      this->m_DirectionTable[i][this->m_NumberOfDirs] = 1;
    }
  this->m_DirectionTable[this->m_NumberOfDirs].resize ( this->m_NumberOfDirs + 1) ;
  for(unsigned int k = 0; k <= this->m_NumberOfDirs; k ++ )
    {
      this->m_DirectionTable[this->m_NumberOfDirs][k] = 1;
    }
  
  this->m_NeighborDirectionTable.resize ( this->m_NumberOfDirs + 1 ) ;
  for ( unsigned i = 0 ; i < this->m_NumberOfDirs ; i++ )
    {
      
      this->m_NeighborDirectionTable[i].resize ( 27 ) ;
      
      for ( unsigned j = 0 ; j < 27 ; j++ )
	{
	  ODFCost::CoordinateType idir, jdir ;
	  idir = this->m_CoordinateTable[i] ;

	  short nx, ny, nz ;
	  nx = j % 3  ;
	  short neighborTemp = ( j - nx ) / 3;
	  nx = nx - 1 ;
	  ny = neighborTemp % 3 ;
	  nz = ( neighborTemp - ny ) / 5 - 1 ;
	  ny = ny - 1 ;  

	  jdir[0] = nx / this->m_DNTable[j] ;
	  jdir[1] = ny / this->m_DNTable[j] ;
	  jdir[2] = nz / this->m_DNTable[j] ;
  
	  double magi = sqrt(idir[0]*idir[0]+idir[1]*idir[1]+idir[2]*idir[2]);
	  double magj = sqrt(jdir[0]*jdir[0]+jdir[1]*jdir[1]+jdir[2]*jdir[2]);
	  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
	  double angle = acos(sum/(magi*magj));
	  
	  if(angle < 0.5*PI)
	    this->m_NeighborDirectionTable[i][j] = 1 ;
	  else
	    this->m_NeighborDirectionTable[i][j] = 0 ;
	}
    }
  this->m_NeighborDirectionTable[this->m_NumberOfDirs].resize ( 27 ) ;
  for( unsigned int p = 0; p < 27; p++)
    {
      this->m_NeighborDirectionTable[this->m_NumberOfDirs][p] = 1;
    }

  
  this->m_NeighborNeighborTable.resize ( 28 ) ;
  for ( unsigned i = 0 ; i < 28 ; i++ )
    {
      
      this->m_NeighborNeighborTable[i].resize ( 28 ) ;
      
      for ( unsigned j = 0 ; j < 28 ; j++ )
	{
	  
	  if(i==27 || j == 27 || i == 13 || j == 13)
	    {
	      this->m_NeighborNeighborTable[i][j] = 1;
	    }
	  else
	    {
	      ODFCost::CoordinateType idir, jdir ;
      
	      short nx1, ny1, nz1 ;
	      nx1 = i % 3  ;
	      short neighborTemp1 = ( i - nx1 ) / 3;
	      nx1-- ;
	      ny1 = neighborTemp1 % 3 ;
	      nz1 = ( neighborTemp1 - ny1 ) / 3 - 1 ;
	      ny1-- ;
      
	      idir[0] = nx1 / this->m_DNTable[i] ;
	      idir[1] = ny1 / this->m_DNTable[i] ;
	      idir[2] = nz1 / this->m_DNTable[i] ;
	      
	      short nx2, ny2, nz2 ;
	      nx2 = j % 3 ;
	      short neighborTemp2 = ( j - nx2 ) / 3;
	      nx2-- ;
	      ny2 = neighborTemp2 % 3 ;
	      nz2 = ( neighborTemp2 - ny2 ) / 3 - 1 ;
	      ny2-- ;
	      
	      jdir[0] = nx2 / this->m_DNTable[j] ;
	      jdir[1] = ny2 / this->m_DNTable[j] ;
	      jdir[2] = nz2 / this->m_DNTable[j] ;

	      double magi = sqrt(idir[0]*idir[0]+idir[1]*idir[1]+idir[2]*idir[2]);
	      double magj = sqrt(jdir[0]*jdir[0]+jdir[1]*jdir[1]+jdir[2]*jdir[2]);
	      double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
	      double angle = acos(sum/(magi*magj));	      
	      if(angle < 0.4*PI)
		this->m_NeighborNeighborTable[i][j] = 1 ;
	      else
		this->m_NeighborNeighborTable[i][j] = 0 ;
	    }
	}
    }
  
}

//added by Wenyu
void ODFCost::CreateFinslerTable()
{
  std::cout<<"Compute Finsler Table"<<std::endl;
  unsigned long int x, y, z, index = 0, indexTimesNDirs = 0 ;
  ODFCost::CoordinateType idir, jdir ;
  unsigned long curdir; 
//  long finslerDirX = -1, finslerDirY = -1, finslerDirZ = -1 ;
  double f_ODF;

//   for ( curdir = 0 ; curdir < this->m_NumberOfDirs ; curdir++ )
//     {
//       ODFCost::CoordinateType dir = this->m_CoordinateTable[curdir];
//       if ( dir[0] == 1 ) finslerDirX = curdir ;
//       if ( dir[1] == 1 ) finslerDirY = curdir ;
//       if ( dir[2] == 1 ) finslerDirZ = curdir ;
//     }
//   if ( ( ( finslerDirX == -1 ) || ( finslerDirY == -1 ) || ( finslerDirZ == -1 ) ) && ( this->m_NumberOfDirs == 6 ) )
//     {
//       std::cout << "couldnt find dirs for writing out finsler images" << std::endl ;
//       exit ( 0 ) ;
//     }
//   else
//     {
//       std::cout << finslerDirX << " " << finslerDirY << " " << finslerDirZ << std::endl ;
//    }

  double min_fodf = 1, max_fodf = 0 ;
  index = -1 ;

  indexTimesNDirs = -this->m_NumberOfDirs ;
  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      index++ ;
	      
	      indexTimesNDirs += this->m_NumberOfDirs ;
	      for ( curdir = 0 ; curdir < this->m_NumberOfDirs ; curdir++ )
		{
		  idir = this->m_CoordinateTable[curdir] ;		  	  

		  double curODF = this->m_ODFArray [ indexTimesNDirs+curdir ] ; 
		  if ( curODF != 0 ) 
		    {
		      double sumODF_proj = 0 ;
		      curODF *= 0.25*this->m_NumberOfDirs;
		      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
			{
			  jdir = this->m_CoordinateTable [dir] ;
			  double dirODF = this->m_ODFArray [indexTimesNDirs+dir];
			  double cosAlpha = idir[0]*jdir[0]+idir[1]*jdir[1]+idir[2]*jdir[2] ;
			  double sinAlpha2 = 1 - cosAlpha * cosAlpha ;
			  double sinAlpha = sinAlpha2 < 0 ? 0 : sqrt ( sinAlpha2 ) ;
			  sumODF_proj += dirODF * sinAlpha ;
			}

		      double finsler = curODF/sumODF_proj ;
		      if ( finsler < 0 ) 
			{
			  std::cout << "Something is very wrong. This should never be negative. Check your ODF file. " << std::endl ;
			  exit ( 0 ) ;
			}
		      else if ( finsler > 1 ) 
			{
			  std::cout << "Something is very wrong. This should never be higher than 1. Check your ODF file and/or the normalization. " << std::endl ;
			  exit ( 0 ) ;
			}
		      f_ODF = 1 - finsler ;

		    }
		  else // curodf== 0
		    {
		      f_ODF = 1 ;
		    }

		  if ( min_fodf > f_ODF )
		    min_fodf = f_ODF ; 
		  if ( max_fodf < f_ODF )
		    max_fodf = f_ODF ;
		  //f_ODF = f_ODF*f_ODF*f_ODF*f_ODF*f_ODF*f_ODF;
		  f_ODF = f_ODF*f_ODF*f_ODF*f_ODF*f_ODF*f_ODF;
		  this->m_FinslerTable.push_back ( f_ODF ) ;
		 
// 		  if( curdir == finslerDirX ) this->m_FinslerArrayX[index] = f_ODF ;
// 		  if( curdir == finslerDirY ) this->m_FinslerArrayY[index] = f_ODF ;
// 		  if( curdir == finslerDirZ ) this->m_FinslerArrayZ[index] = f_ODF ;
		}
	     
	    }
	}
    }
  
  std::cout<<"min_fodf is "<<min_fodf<<","<<"max_fodf is "<<max_fodf<<std::endl;

//   if ( ( finslerDirX == -1 ) || ( finslerDirY == -1 ) || ( finslerDirZ == -1 ) )
//     return ;

//   //print finsler into a finsler file
//   std::cout<<"Start writing finsler files"<<std::endl;
//   std::string fileNameX = "finslerX.nrrd";
//   DoubleWriterType::Pointer finWriterX = DoubleWriterType::New();
//   finWriterX->SetFileName(fileNameX);
//   finWriterX->SetInput( this->m_FinslerImageX );
//   try
//     {
//       finWriterX->Update();
//     }
//   catch( itk::ExceptionObject & excp )
//     {
//       std::cerr << "Problem writing the finsler file " << fileNameX << std::endl ;
//       std::cerr << excp << std::endl;
//     }

//   std::string fileNameY = "finslerY.nrrd";
//   DoubleWriterType::Pointer finWriterY = DoubleWriterType::New();
//   finWriterY->SetFileName(fileNameY);
//   finWriterY->SetInput( this->m_FinslerImageY );
//   try
//     {
//       finWriterY->Update();
//     }
//   catch( itk::ExceptionObject & excp )
//     {
//       std::cerr << "Problem writing the finsler file " << fileNameY << std::endl ;
//       std::cerr << excp << std::endl;
//     }

//   std::string fileNameZ = "finslerZ.nrrd";
//   DoubleWriterType::Pointer finWriterZ = DoubleWriterType::New();
//   finWriterZ->SetFileName(fileNameZ);
//   finWriterZ->SetInput( this->m_FinslerImageZ );
//   try
//     {
//       finWriterZ->Update();
//     }
//   catch( itk::ExceptionObject & excp )
//     {
//       std::cerr << "Problem writing the finsler file " << fileNameZ << std::endl ;
//       std::cerr << excp << std::endl;
//     }

//   exit ( 0 ) ;
}

// for each voxel, the final cost is the minimum cost among the costs associated with all the directions
void ODFCost::ComputeFinalCosts ()
{
  unsigned long int x, y, z, index = 0, indexTimesNDirs = 0 ;
  double currentCost, minCost, minLength ;
  unsigned int minOrig;
  for ( z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( x = 0 ; x < this->m_XDim ; x++ )
	    {
	      minCost = INF ;
	      minOrig = -1 ;
	      minLength = -1 ;
	      unsigned long mindir = 0;
	      double minODF = 2;
	      double odfcost;
	      unsigned int odfdir;
	      double sum = 0;
	      for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
		{
		 
		  currentCost = this->m_AvgCostArray[indexTimesNDirs+dir];
		  sum += currentCost;		    

		  if ( ( currentCost < minCost ) && ( currentCost >= 0 ) )
		    {
		      mindir = dir;
		      minCost = currentCost ;
		      minOrig = this->m_OrigArray[indexTimesNDirs+dir] ;
		      minLength = this->m_LengthArray[indexTimesNDirs+dir] ;
		    }

		  if(minODF > this->m_FinslerTable[indexTimesNDirs+ dir])
		    {
		      minODF =  this->m_FinslerTable[indexTimesNDirs+ dir];
		      odfdir = dir;
		      odfcost = currentCost;
		    }
		  
		} 

	      if ( minCost == INF ) minCost = -1 ;
	     
	      this->m_FinalCostArray[index] = minCost ;
	      
	      this->m_FinalOrigArray[index] = minOrig ;
	      this->m_FinalLengthArray[index] = minLength ;
	      index++ ;
	      indexTimesNDirs += this->m_NumberOfDirs ;
	    }
	}
    }

}


 //for 5x5x5 neighbors
/*
unsigned long ODFCost::GetNeighbor ( unsigned long x, unsigned long y, unsigned long z, unsigned short neighbor, bool &valid) 
{
  short nx, ny, nz ;
  nx = neighbor % 5;
  neighbor = ( neighbor - nx ) / 5;
  nx = nx - 2 ;
  ny = neighbor % 5 ;
  nz = ( neighbor - ny ) / 5 - 2 ;
  ny = ny - 2 ;

  short neighborX = nx + x ;
  short neighborY = ny + y ;
  short neighborZ = nz + z ;

  valid = true ;
  if ( neighborX < 0 || neighborX >= this->m_XDim || neighborY < 0 || neighborY >= this->m_YDim || neighborZ < 0 || neighborZ >= this->m_ZDim )
    valid = false ;

  return neighborX + neighborY * this->m_XDim + neighborZ * this->m_Slice ;
}

*/

// get the neighbor index (x+y*xDim+z*xDim*yDim) given current voxel and 0 <= neighbor < 27
// if the neighbor is out of the image volume, "valid" becomes false
unsigned long ODFCost::GetNeighbor ( unsigned long x, unsigned long y, unsigned long z, unsigned short neighbor, bool &valid) 
{
  short nx, ny, nz ;
  nx = neighbor % 3;
  neighbor = ( neighbor - nx) / 3;
  nx-- ;
  ny = neighbor % 3 ;
  nz = ( neighbor - ny ) / 3 - 1 ;
  ny-- ;
  
  short neighborX = nx + x ;
  short neighborY = ny + y ;
  short neighborZ = nz + z ;
  
  valid = true ;
  if ( neighborX < 0 || neighborX >= this->m_XDim || neighborY < 0 || neighborY >= this->m_YDim || neighborZ < 0 || neighborZ >= this->m_ZDim )
    valid = false ;

  return neighborX + neighborY * this->m_XDim + neighborZ * this->m_Slice ;
}


// if the current path is cheaper than the current best path, update
//modified by Wenyu
//void ODFCost::UpdateCost ( unsigned long index, double newCost, double newLength, unsigned int newOrig, unsigned int dir)
void ODFCost::UpdateCost ( unsigned long index, unsigned long voxelIndex, double newCost,double oldAvgCost, double newLength, unsigned int newOrig, unsigned int indir,unsigned int outdir, unsigned int neighbor)
{

  //static int update_counter = 0;
  double currentAvgCost = this->m_AvgCostArray[index] ;
  double newAvgCost = newCost / newLength ;
  
  if ( newAvgCost < oldAvgCost)
    {
      newAvgCost = oldAvgCost;
    
    } 
  
  if ( newAvgCost*1.05 < currentAvgCost ) 
    {
      this->m_CostArray[index] = newCost ;
      this->m_Converged = false ;
      this->m_NChanged++ ;
      this->m_LengthArray[index] = newLength ;
      this->m_AvgCostArray[index] = newAvgCost ;
      this->m_OrigArray[index] = newOrig;
      //added by Wenyu
      //this->m_OrigDir[voxelIndex] = indir;
      this->m_OrigInDir[voxelIndex] = indir;
      this->m_OrigOutDir[voxelIndex] = outdir;
      //std::cout<<"Orig dir is "<<dir<<std::endl;
      this->m_OrigNeighbor[voxelIndex] = neighbor;
      // std::cout<<"Orig Neighbor is "<<neighbor<<std::endl;
      //std::cout<<"avg cost is "<<newAvgCost<<std::endl;
    }
}

// the cost of one step, given current voxel, neighbor, in and out directions
double ODFCost::stepCost ( unsigned inIndex, unsigned outIndex, unsigned neighborIndex, unsigned long voxelIndexTimesNDirs, unsigned long neighborVoxelIndexTimesNDirs ,unsigned long neighborVoxel ) 
{
 
  unsigned int origDir =  this->m_OrigDir[ neighborVoxel ];
  unsigned int origNeighbor = this->m_OrigNeighbor[neighborVoxel ];
  unsigned int origInDir = this->m_OrigInDir[ neighborVoxel ];
  unsigned int origOutDir = this->m_OrigOutDir[ neighborVoxel];
 
  double preAngleCost =  this->m_AnglePenaltyTable[inIndex][origInDir]+this->m_AnglePenaltyTable[outIndex][origOutDir];

   double f_odf1 =  this->m_FinslerTable[ voxelIndexTimesNDirs + inIndex];

  double f_odf2 = this->m_FinslerTable[neighborVoxelIndexTimesNDirs + outIndex ];

  double cost = (this->m_NeighborAnglePenaltyTable[inIndex][neighborIndex] + this->m_NeighborAnglePenaltyTable[outIndex][neighborIndex]+this->m_AnglePenaltyTable[inIndex][outIndex])+preAngleCost;

  //For synthetic data, we do not need the preAngleCost part
  //double cost = (this->m_NeighborAnglePenaltyTable[inIndex][neighborIndex] + this->m_NeighborAnglePenaltyTable[outIndex][neighborIndex]+this->m_AnglePenaltyTable[inIndex][outIndex]);

  double f_odf = f_odf1 + f_odf2;

  return this->m_alpha*f_odf + cost;

}

// the penalty associated with the angle discrepancy between two directions (arguments are indices from 0 to nDirs)
double ODFCost::AnglePenalty ( unsigned u, unsigned v )
{
  ODFCost::CoordinateType idir, jdir ;
  idir = this->m_CoordinateTable[u] ;
  jdir = this->m_CoordinateTable[v] ;
  
  
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2];
  double angle = acos(sum);
  double e = 0.5 * ( 1. - sum ) ;
 
  return e;
  
}

//5x5x5 neighbors
// penalty associated with the angle discrepancy between two directions ( 0 <= u < nDirs, 0 <= neighbor < 124 )
/*
double ODFCost::NeighborAnglePenalty ( unsigned u, unsigned neighbor )
{
  ODFCost::CoordinateType idir, jdir ;
  idir = this->m_CoordinateTable[u] ;

  short nx, ny, nz ;
  nx = neighbor % 5  ;
  short neighborTemp = ( neighbor - nx ) / 5;
  nx = nx - 2 ;
  ny = neighborTemp % 5 ;
  nz = ( neighborTemp - ny ) / 5 - 2 ;
  ny = ny - 2 ;
  
  jdir[0] = nx / this->m_DNTable[neighbor] ;
  jdir[1] = ny / this->m_DNTable[neighbor] ;
  jdir[2] = nz / this->m_DNTable[neighbor] ;
  
  
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
  double e = 0.5 * ( 1. - sum ) ;
  double angle = acos(sum);
 
  return e ;
  
}
*/

 //3x3x3 neighbors
// penalty associated with the angle discrepancy between two directions ( 0 <= u < nDirs, 0 <= neighbor < 27 )
double ODFCost::NeighborAnglePenalty ( unsigned u, unsigned neighbor )
{
  ODFCost::CoordinateType idir, jdir, kdir ;
  idir = this->m_CoordinateTable[u] ;

  //std::cout<<"coornate of dir "<<u<<" : "<<idir[0]<<" "<<idir[1]<<" "<<idir[2]<<std::endl;

  short nx, ny, nz ;
  nx = neighbor % 3  ;
  short neighborTemp = ( neighbor - nx ) / 3;
  nx-- ;
  ny = neighborTemp % 3 ;
  nz = ( neighborTemp - ny ) / 3 - 1 ;
  ny-- ;
  
  jdir[0] = nx / this->m_DNTable[neighbor] ;
  jdir[1] = ny / this->m_DNTable[neighbor] ;
  jdir[2] = nz / this->m_DNTable[neighbor] ;
  
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
  
  double e = 0.5 * ( 1. - sum ) ;
  return e ;
  
}

double ODFCost::NeighborPenalty ( unsigned neighbor1, unsigned neighbor2 )
{
  ODFCost::CoordinateType idir, jdir ;

  short nx1, ny1, nz1 ;
  nx1 = neighbor1 % 3  ;
  short neighborTemp1 = ( neighbor1 - nx1 ) / 3;
  nx1-- ;
  ny1 = neighborTemp1 % 3 ;
  nz1 = ( neighborTemp1 - ny1 ) / 3 - 1 ;
  ny1-- ;
  
  idir[0] = nx1 / this->m_DNTable[neighbor1] ;
  idir[1] = ny1 / this->m_DNTable[neighbor1] ;
  idir[2] = nz1 / this->m_DNTable[neighbor1] ;

  short nx2, ny2, nz2 ;
  nx2 = neighbor2 % 3 ;
  short neighborTemp2 = ( neighbor2 - nx2 ) / 3;
  nx2-- ;
  ny2 = neighborTemp2 % 3 ;
  nz2 = ( neighborTemp2 - ny2 ) / 3 - 1 ;
  ny2-- ;
  
  jdir[0] = nx2 / this->m_DNTable[neighbor2] ;
  jdir[1] = ny2 / this->m_DNTable[neighbor2] ;
  jdir[2] = nz2 / this->m_DNTable[neighbor2] ;
 
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
  
  double e = 0.5 * ( 1. - sum ) ;
  return e ;
 
}

//steplegnth for 5x5x5 neighbors
/*
double ODFCost::stepLength ( unsigned neighbor ) 
{
  short nx, ny, nz ;
  nx = neighbor % 5  ;
  neighbor = ( neighbor - nx ) / 5; 
  nx = nx - 2 ;
  ny = neighbor % 5 ;
  nz = ( neighbor - ny ) / 5 - 2 ;
  ny = ny - 2 ;

  return sqrt(nx*nx+ny*ny+nz*nz);
}
*/

// a face-neighbor has a step length of 1, edge-neighbor \sqrt(2), vertex-neighbor \sqrt(3)
// 0 <= neighbor < 27
double ODFCost::stepLength ( unsigned neighbor ) 
{
  short nx, ny, nz ;
  nx = neighbor % 3  ;
  neighbor = ( neighbor - nx ) / 3;
  nx-- ;
  ny = neighbor % 3 ;
  nz = ( neighbor - ny ) / 3 - 1 ;
  ny-- ;

  short dn = abs ( nx ) + abs ( ny ) + abs ( nz ) ;
  if ( dn == 0 ) return 0 ;
  if ( dn == 1 ) return 1 ;
  if ( dn == 2 ) return sqrt ( 2 ) ;
  return sqrt ( 3 ) ;
}

//5x5x5 neighbors
/*
void ODFCost::PrecomputeStepLength () 
{
  this->m_DNTable.resize ( 124 ) ;
  for ( unsigned i = 0 ; i < 124 ; i++ )
    {    
      this->m_DNTable[i] = this->stepLength ( i ) ;
    }
}
*/

// compute step lengths and store for faster access
void ODFCost::PrecomputeStepLength () 
{
  this->m_DNTable.resize ( 27 ) ;
  for ( unsigned i = 0 ; i < 27 ; i++ )
    {    
      this->m_DNTable[i] = this->stepLength ( i ) ;
    }
}


// compute angle penalties and store for faster access
void ODFCost::PrecomputeAnglePenalty () 
{
  std::cout<<"precompute angle penalty"<<std::endl;
  this->m_AnglePenaltyTable.resize ( this->m_NumberOfDirs+1 ) ;
  for ( unsigned i = 0 ; i <= this->m_NumberOfDirs ; i++ )
    {
      this->m_AnglePenaltyTable[i].resize ( this->m_NumberOfDirs + 1 ) ;
      for ( unsigned j = 0 ; j <= this->m_NumberOfDirs ; j++ )
	{
	  if(i ==  this->m_NumberOfDirs || j ==  this->m_NumberOfDirs)
	    {
	      this->m_AnglePenaltyTable[i][j] = 0;
	    }
	  else
	    {
	      this->m_AnglePenaltyTable[i][j] = this->AnglePenalty ( i, j ) ;
	    }
	 
	}
    }

  this->m_NeighborAnglePenaltyTable.resize ( this->m_NumberOfDirs+1 ) ;
  
  for ( unsigned i = 0 ; i <= this->m_NumberOfDirs ; i++ )
    {
      //5x5x5 neighbors
      /*
      this->m_NeighborAnglePenaltyTable[i].resize ( 124 ) ;
      for ( unsigned j = 0 ; j < 124 ; j++ )
	{
	  this->m_NeighborAnglePenaltyTable[i][j] = this->NeighborAnglePenalty ( i, j ) ;
	}
      */
      
      this->m_NeighborAnglePenaltyTable[i].resize ( 27 ) ;
      for ( unsigned j = 0 ; j < 27 ; j++ )
	{
	  if(i ==  this->m_NumberOfDirs)
	    {
	      this->m_NeighborAnglePenaltyTable[i][j] = 0 ;
	    }
	  else
	    {
	      this->m_NeighborAnglePenaltyTable[i][j] = this->NeighborAnglePenalty ( i, j ) ;
	      if(j == 13)
		{
		  this->m_NeighborAnglePenaltyTable[i][j] = 0;
		}
	     }
	}
      
    }
  this->m_NeighborPenaltyTable.resize ( 27 ) ;
  
  for ( unsigned i = 0 ; i < 27 ; i++ )
    {
      
      this->m_NeighborPenaltyTable[i].resize ( 27 ) ;
      for ( unsigned j = 0 ; j < 27 ; j++ )
	{
		 
	  this->m_NeighborPenaltyTable[i][j] = this->NeighborPenalty ( i, j ) ;
	  if(i==13 || j == 13)
	    {
	      this->m_NeighborPenaltyTable[i][j] = 0;
	    }
	  
	}
      
    }
}

// we want the odf values to be 0..1
void ODFCost::NormalizeODF () 
{
   double min = INF, max = -INF, current ;
   for ( unsigned long index = 0 ; index < this->m_ZDim * this->m_YDim * this->m_XDim * this->m_NumberOfDirs ; index++ )
     {
       current = this->m_ODFArray[index] ;
       if ( current < min ) min = current ;
       if ( current > max ) max = current ;
     }
   std::cout << "ODF min-max: " << min << " " << max << std::endl ;
   double range = max - min ;

  double sum ;
  unsigned long index, indexTimesNDir = 0 ;
  for ( index = 0 ; index < this->m_ZDim * this->m_YDim * this->m_XDim ; index++ )
    {
      sum = 0 ;
      for ( unsigned dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
	{
	  sum+= this->m_ODFArray[indexTimesNDir+dir] - min ;
	}
      for ( unsigned dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
	{
	  current = this->m_ODFArray[indexTimesNDir+dir] - min ;
	  current /= sum ;
	  this->m_ODFArray[indexTimesNDir+dir] = current ;
	}
      indexTimesNDir += this->m_NumberOfDirs ;
    }
}
/*
 //added by Wenyu
 //a new 5x5x5 f star traversal
 //62 is the voxel itself
void ODFCost::FstarTraversalUpdates ()
{
  int i,j;
  this->m_fstar_updates[0][0] = 49;
  this->m_fstar_updates[0][1] = 1;
  this->m_fstar_updates[0][2] = 3;
  this->m_fstar_updates[0][3] = 5;
  this->m_fstar_updates[0][4] = 6;
  this->m_fstar_updates[0][5] = 7;
  this->m_fstar_updates[0][6] = 8;
  this->m_fstar_updates[0][7] = 9;
  this->m_fstar_updates[0][8] = 11;
  this->m_fstar_updates[0][9] = 13;
  this->m_fstar_updates[0][10] = 15;
  this->m_fstar_updates[0][11] = 16;
  this->m_fstar_updates[0][12] = 17;
  this->m_fstar_updates[0][13] = 18;
  this->m_fstar_updates[0][14] = 19;
  this->m_fstar_updates[0][15] = 21;
  this->m_fstar_updates[0][16] = 23;

  for(i = 17,j = 25;i <= 41; i ++, j++)
    {
      this->m_fstar_updates[0][i] = j;
    }

  this->m_fstar_updates[0][42] = 51;
  this->m_fstar_updates[0][43] = 53;
  this->m_fstar_updates[0][44] = 55;
  this->m_fstar_updates[0][45] = 56;
  this->m_fstar_updates[0][46] = 57;
  this->m_fstar_updates[0][47] = 58;
  this->m_fstar_updates[0][48] = 59;
  this->m_fstar_updates[0][49] = 61;  

  this->m_fstar_updates[1][0] = 1 ;
  this->m_fstar_updates[1][1] = 63;

  this->m_fstar_updates[2][0] = 8 ;
  this->m_fstar_updates[2][1] = 61 ;
  this->m_fstar_updates[2][2] = 65 ;
  this->m_fstar_updates[2][3] = 66 ;
  this->m_fstar_updates[2][4] = 67 ;
  this->m_fstar_updates[2][5] = 68 ;
  this->m_fstar_updates[2][6] = 69 ;
  this->m_fstar_updates[2][7] = 71 ;
  this->m_fstar_updates[2][8] = 73 ;
  this->m_fstar_updates[2][1] = 71 ;

  this->m_fstar_updates[3][0] = 1 ;
  this->m_fstar_updates[3][1] = 63 ;

  this->m_fstar_updates[4][0] = 49 ;
  this->m_fstar_updates[4][1] = 63 ;
  this->m_fstar_updates[4][2] = 65 ;
  this->m_fstar_updates[4][3] = 66 ;
  this->m_fstar_updates[4][4] = 67 ;
  this->m_fstar_updates[4][5] = 68 ;
  this->m_fstar_updates[4][6] = 69 ;
  this->m_fstar_updates[4][7] = 71 ;
  this->m_fstar_updates[4][8] = 73 ;

  for(i = 9,j = 75; i <= 33; i++, j++)
    {
      this->m_fstar_updates[4][i] = j ;
    }
  this->m_fstar_updates[4][34] = 101 ;
  this->m_fstar_updates[4][35] = 103 ;
  this->m_fstar_updates[4][36] = 105 ;
  this->m_fstar_updates[4][37] = 106;
  this->m_fstar_updates[4][38] = 107;
  this->m_fstar_updates[4][39] = 108;
  this->m_fstar_updates[4][40] = 109;
  this->m_fstar_updates[4][41] = 111;
  this->m_fstar_updates[4][42] = 113;
  this->m_fstar_updates[4][43] = 115;
  this->m_fstar_updates[4][44] = 116;
  this->m_fstar_updates[4][45] = 117;
  this->m_fstar_updates[4][46] = 118;
  this->m_fstar_updates[4][47] = 119;
  this->m_fstar_updates[4][48] = 121;
  this->m_fstar_updates[4][49] = 123;
   
  this->m_fstar_updates[5][0] = 1 ;
  this->m_fstar_updates[5][1] = 63 ;

  this->m_fstar_updates[6][0] = 8 ;
  this->m_fstar_updates[6][1] = 61 ;
  this->m_fstar_updates[6][2] = 65 ;
  this->m_fstar_updates[6][3] = 66 ;
  this->m_fstar_updates[6][4] = 67 ;
  this->m_fstar_updates[6][5] = 68 ;
  this->m_fstar_updates[6][6] = 69 ;
  this->m_fstar_updates[6][7] = 71 ;
  this->m_fstar_updates[6][8] = 73 ;  

  this->m_fstar_updates[7][0] = 1 ;
  this->m_fstar_updates[7][1] = 63;
}
*/

// initialize the neighbors that need visiting along each fstar traversal direction
// this table is structured as follows:
// m_fstar_updates[i][0] = number (n) of neighbors that need visiting if we are on the i-direction (0<=i<8)
// m_fstar_updates[i][1..n] = the neighbors that need visiting
// the neighbors are indexed 0..27 (so 13 would be the voxel itself)
void ODFCost::FstarTraversalUpdates ()
{
  this->m_fstar_updates[0][0] = 13 ;
  this->m_fstar_updates[0][1] = 0 ;
  this->m_fstar_updates[0][2] = 1 ;
  this->m_fstar_updates[0][3] = 2 ;
  this->m_fstar_updates[0][4] = 3 ;
  this->m_fstar_updates[0][5] = 4 ;
  this->m_fstar_updates[0][6] = 5 ;
  this->m_fstar_updates[0][7] = 6 ;
  this->m_fstar_updates[0][8] = 7 ;
  this->m_fstar_updates[0][9] = 8 ;
  this->m_fstar_updates[0][10] = 9 ;
  this->m_fstar_updates[0][11] = 10 ;
  this->m_fstar_updates[0][12] = 11 ;
  this->m_fstar_updates[0][13] = 12 ;

  this->m_fstar_updates[1][0] = 1 ;
  this->m_fstar_updates[1][1] = 14 ;

  this->m_fstar_updates[2][0] = 4 ;
  this->m_fstar_updates[2][1] = 12 ;
  this->m_fstar_updates[2][2] = 15 ;
  this->m_fstar_updates[2][3] = 16 ;
  this->m_fstar_updates[2][4] = 17 ;

  this->m_fstar_updates[3][0] = 1 ;
  this->m_fstar_updates[3][1] = 14 ;

  this->m_fstar_updates[4][0] = 13 ;
  this->m_fstar_updates[4][1] = 14 ;
  this->m_fstar_updates[4][2] = 15 ;
  this->m_fstar_updates[4][3] = 16 ;
  this->m_fstar_updates[4][4] = 17 ;
  this->m_fstar_updates[4][5] = 18 ;
  this->m_fstar_updates[4][6] = 19 ;
  this->m_fstar_updates[4][7] = 20 ;
  this->m_fstar_updates[4][8] = 21 ;
  this->m_fstar_updates[4][9] = 22 ;
  this->m_fstar_updates[4][10] = 23 ;
  this->m_fstar_updates[4][11] = 24 ;
  this->m_fstar_updates[4][12] = 25 ;
  this->m_fstar_updates[4][13] = 26 ;

  this->m_fstar_updates[5][0] = 1 ;
  this->m_fstar_updates[5][1] = 14 ;

  this->m_fstar_updates[6][0] = 4 ;
  this->m_fstar_updates[6][1] = 12 ;
  this->m_fstar_updates[6][2] = 15 ;
  this->m_fstar_updates[6][3] = 16 ;
  this->m_fstar_updates[6][4] = 17 ;

  this->m_fstar_updates[7][0] = 1 ;
  this->m_fstar_updates[7][1] = 14 ;
}


