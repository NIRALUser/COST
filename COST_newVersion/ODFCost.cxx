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

  //this->PrecomputeDirection();

  this->PrecomputeNeighbors () ;

  std::cout << "Entering main loop" << std::endl ;
  unsigned i = 0 ;

  while ( !this->m_Converged ) 
    {
      iterationNum ++ ;
      this->MakeOnePassFstar () ;
      i++ ;
      std::cout << "Pass " << i << " complete. " << this->m_NChanged << " voxels updated in this pass." << std::endl ;
      
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

   //this->m_OrigInDirImage = DoubleImageType::New () ;
   //this->m_OrigOutDirImage = DoubleImageType::New () ;

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
   //this->m_OrigInDirImage->SetSpacing ( spacing ) ;
   //this->m_OrigOutDirImage->SetSpacing ( spacing ) ;

   ODFImageType::PointType origin ;
   origin = this->m_ODFImage->GetOrigin () ;
   this->m_CostImage->SetOrigin ( origin ) ;
   this->m_FinalCostImage->SetOrigin ( origin ) ;
   this->m_LengthImage->SetOrigin ( origin ) ;
   this->m_AvgCostImage->SetOrigin ( origin ) ;
   this->m_OrigImage->SetOrigin ( origin ) ;
   this->m_FinalOrigImage->SetOrigin ( origin ) ;
   this->m_FinalLengthImage->SetOrigin ( origin ) ;
   //this->m_OrigInDirImage->SetOrigin ( origin ) ;
   //this->m_OrigOutDirImage->SetOrigin ( origin ) ;

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

   //this->m_OrigInDirImage->SetRegions ( region ) ;
   //this->m_OrigInDirImage->Allocate () ;
   //this->m_OrigInDirImage->FillBuffer ( 0 ) ;


   //this->m_OrigOutDirImage->SetRegions ( region ) ;
   //this->m_OrigOutDirImage->Allocate () ;
   //this->m_OrigOutDirImage->FillBuffer ( 0 ) ;

   this->m_CostArray = this->m_CostImage->GetPixelContainer()->GetBufferPointer();
   this->m_AvgCostArray = this->m_AvgCostImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalCostArray = this->m_FinalCostImage->GetPixelContainer()->GetBufferPointer() ;
   this->m_LengthArray = this->m_LengthImage->GetPixelContainer()->GetBufferPointer();
   this->m_OrigArray = this->m_OrigImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalOrigArray = this->m_FinalOrigImage->GetPixelContainer()->GetBufferPointer();
   this->m_FinalLengthArray = this->m_FinalLengthImage->GetPixelContainer()->GetBufferPointer();

 //added by Wenyu
   //this->m_OrigInDir = this->m_OrigInDirImage->GetPixelContainer()->GetBufferPointer();
   //this->m_OrigOutDir = this->m_OrigOutDirImage->GetPixelContainer()->GetBufferPointer();

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
	      
	      //this->m_OrigInDir[index] = this->m_NumberOfDirs;
	      //this->m_OrigOutDir[index] = this->m_NumberOfDirs;
	      
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
  if ( this->m_OrigFAArray[voxelIndex] < 0.05 ||  this->m_SourceArray[voxelIndex] )
    {
      // not worth the computation time TODO: add as command line argument
      return ;
    }

  unsigned long neighborIndex, neighborIndexTimesNDirs ;
  unsigned short neighborIterator ;
  double newStepCost ;
  //unsigned int origInDir ;
  //unsigned int origOutDir ;
  
  // instead of all 26 neighbors, only visit the neighbors as indicated by fstar for our current traversal direction
  for ( unsigned short n = 0 ; n < this->m_fstar_updates[this->m_fstar_direction][0] ; n++ )
    {
      neighborIterator = this->m_fstar_updates[this->m_fstar_direction][n+1] ;

      if ( ! this->m_NeighborValidTable[voxelIndex][neighborIterator] ) continue ;
      neighborIndex = this->m_NeighborIndexTable[voxelIndex][neighborIterator] ;

      neighborIndexTimesNDirs = neighborIndex * this->m_NumberOfDirs ;
      double dn = this->m_DNTable[neighborIterator] ;

      //origInDir = this->m_OrigInDir[ neighborIndex ];
      //origOutDir = this->m_OrigOutDir[ neighborIndex];
      
      // for each incoming dir
      for ( unsigned int inDir = 0 ; inDir < this->m_NumberOfDirs ; inDir++ )  
	{
	  unsigned int voxelAndIn = voxelIndexTimesNDirs + inDir ;

	  // for each outgoing dir
	   for ( unsigned int outDir = 0 ; outDir < this->m_NumberOfDirs ; outDir++ )
	    {
	      unsigned int neighborAndOut = neighborIndexTimesNDirs+outDir ;

	      //I use this constraints for mouse data, it's unneccessary for synthetic data
	      //if(this->m_DirectionTable[inDir][origInDir] && this->m_DirectionTable[outDir][origOutDir])
	     	{
 		  //double preAngleCost =  this->m_AnglePenaltyTable[inDir][origInDir] + this->m_AnglePenaltyTable[outDir][origOutDir] ;

		  double f_odf1 =  this->m_FinslerTable[voxelAndIn];  
		  double f_odf2 = this->m_FinslerTable[neighborAndOut];
		  double f_odf = f_odf1 + f_odf2;

		  //For synthetic data, we do not need the preAngleCost part
		  //double cost = this->m_NeighborAnglePenaltyTable[inDir][neighborIterator] + this->m_NeighborAnglePenaltyTable[outDir][neighborIterator]+this->m_AnglePenaltyTable[inDir][outDir]+preAngleCost;
		  double cost = this->m_NeighborAnglePenaltyTable[inDir][neighborIterator] + this->m_NeighborAnglePenaltyTable[outDir][neighborIterator]+this->m_AnglePenaltyTable[inDir][outDir];

		  newStepCost = (  this->m_alpha*f_odf + cost ) * dn ;

		  this->UpdateCost ( voxelAndIn, voxelIndex, newStepCost + this->m_CostArray[neighborAndOut], this->m_AvgCostArray[neighborAndOut], dn+this->m_LengthArray[neighborAndOut], this->m_OrigArray[neighborAndOut], inDir, outDir, neighborIterator ) ;	   
		}
	    }
	}
    }
}

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
}

//added by Wenyu
void ODFCost::CreateFinslerTable()
{
  std::cout<<"Compute Finsler Table"<<std::endl;
  unsigned long int x, y, z, index = 0, indexTimesNDirs = 0 ;
  ODFCost::CoordinateType idir, jdir ;
  unsigned long curdir; 
  double f_ODF;

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
		      curODF *= 0.25*this->m_NumberOfDirs; // NOT SURE WHAT THE APPROPRIATE COEFFICIENT HERE SHOULD REALLY BE
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
			  std::cout << curODF << " " << sumODF_proj << std::endl ;
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
		  f_ODF = f_ODF*f_ODF*f_ODF*f_ODF*f_ODF*f_ODF;
		  this->m_FinslerTable.push_back ( f_ODF ) ;
		}
	    }
	}
    }
  
  std::cout<<"min_fodf is "<<min_fodf<<","<<"max_fodf is "<<max_fodf<<std::endl;
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


void ODFCost::PrecomputeNeighbors ()
{
  std::cout << "Precompute neighbors" << std::endl ;
  unsigned long nVoxels = this->m_XDim * this->m_YDim * this->m_ZDim ;   
  this->m_NeighborIndexTable.resize ( nVoxels ) ;
  this->m_NeighborValidTable.resize ( nVoxels ) ;

  unsigned long neighborIndex ;
  bool neighborValid ;

  for ( unsigned int z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( unsigned int y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( unsigned int x = 0 ; x < this->m_XDim ; x++ )
	    {
	      unsigned long currentVoxel = x + y * this->m_XDim + z * this->m_Slice ;
	      this->m_NeighborIndexTable[currentVoxel].resize ( 27 ) ;
	      this->m_NeighborValidTable[currentVoxel].resize ( 27 ) ;

	      for ( unsigned int neighbor = 0 ; neighbor < 27 ; neighbor++ )
		{
		  neighborIndex = this->GetNeighbor ( x, y, z, neighbor, neighborValid ) ;
		  this->m_NeighborIndexTable[currentVoxel][neighbor] = neighborIndex ;
		  this->m_NeighborValidTable[currentVoxel][neighbor] = neighborValid ;
		}
	    }
	}
   }
}

// if the current path is cheaper than the current best path, update
//modified by Wenyu
//void ODFCost::UpdateCost ( unsigned long index, double newCost, double newLength, unsigned int newOrig, unsigned int dir)
void ODFCost::UpdateCost ( unsigned long index, unsigned long voxelIndex, double newCost,double oldAvgCost, double newLength, unsigned int newOrig, unsigned int indir,unsigned int outdir, unsigned int neighbor)
{
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
      //this->m_OrigInDir[voxelIndex] = indir;
      //this->m_OrigOutDir[voxelIndex] = outdir;
    }
}

// the penalty associated with the angle discrepancy between two directions (arguments are indices from 0 to nDirs)
double ODFCost::AnglePenalty ( unsigned u, unsigned v )
{
  ODFCost::CoordinateType idir, jdir ;
  idir = this->m_CoordinateTable[u] ;
  jdir = this->m_CoordinateTable[v] ;
  
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2];
  //double angle = acos(sum);
  double e = 0.5 * ( 1. - sum ) ;
 
  return e;
}

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

