#include "ODFStreamline.h"

// initialize internal structures
ODFStreamline::ODFStreamline ( ODFStreamline::ODFImageType::Pointer odf, unsigned long verts )
{ 
   this->m_ODFImage = odf ;
   this->m_NumberOfDirs = verts ;
   this->m_NumberOfSourceVoxels = 0 ;

   //create empty output image
   this->m_StreamlineImage = LabelImageType::New () ;

   LabelImageType::IndexType start;
   LabelImageType::SizeType  size;
   start[0] = 0;   start[1] = 0;   start[2] = 0;
   size = odf->GetLargestPossibleRegion().GetSize() ;
   this->m_XDim = size[0] ;
   this->m_YDim = size[1] ;
   this->m_ZDim = size[2] ;
   std::cout << "Source image dimensions: " << size[0] << " " << size[1] << " " << size[2] << std::endl ;

   LabelImageType::SpacingType spacing;
   spacing = odf->GetSpacing () ;
   this->m_StreamlineImage->SetSpacing ( spacing ) ;
   LabelImageType::RegionType region;
   region.SetSize ( size ) ;
   region.SetIndex ( start ) ;
   this->m_StreamlineImage->SetRegions ( region ) ;
   this->m_StreamlineImage->Allocate () ;

   this->m_ODFArray = this->m_ODFImage->GetPixelContainer()->GetBufferPointer();
   this->m_StreamlineArray = this->m_StreamlineImage->GetPixelContainer()->GetBufferPointer();
   this->m_StepSize = 1 ;
}

// main functionality - public
ODFStreamline::LabelImageType::Pointer ODFStreamline::Streamline ( ODFStreamline::LabelImageType::Pointer source, unsigned long maxLength, bool singleTract ) 
{
  std::cout << "Starting streamline algorithm" << std::endl ;
  m_SourceArray = source->GetPixelContainer()->GetBufferPointer();

  this->FindSourceVoxels () ;
  std::cout << "Number of source voxels: " << m_NumberOfSourceVoxels << std::endl ;
  if ( singleTract ) 
    this->m_NumberOfSourceVoxels = 1 ;

  for ( unsigned long int i = 0 ; i < m_NumberOfSourceVoxels ; i++ )
    {
      // go both 'forward' and 'backward' from each source voxel
      this->StreamlineFromVoxel ( i, maxLength, 1, true ) ;
      this->StreamlineFromVoxel ( i, maxLength, 2, false ) ;
    }

  return this->m_StreamlineImage ;
}

// trace the streamline from a single voxel 
void ODFStreamline::StreamlineFromVoxel ( unsigned long sourceVoxelId, unsigned long maxLength, unsigned long label, bool forward )
{
  ODFStreamline::CoordinateType currentVoxel ;
  unsigned long slice = this->m_XDim * this->m_YDim ;
  unsigned long index ;
  unsigned long maxODFindex ;
  LabelImageType::SpacingType spacing = this->m_StreamlineImage->GetSpacing () ;
  ODFStreamline::CoordinateType oneStep ;

  oneStep[0] = spacing[0] * this->m_StepSize ;
  oneStep[1] = spacing[1] * this->m_StepSize ;
  oneStep[2] = spacing[2] * this->m_StepSize ;

  int lastStep = -1 ;
  if ( ! forward ) lastStep = -2 ;

  currentVoxel = this->m_SourceVoxels[sourceVoxelId] ;

  for ( unsigned long step = 0 ; step < maxLength ; step++ )
    {
      index = round ( currentVoxel[0] ) + round ( currentVoxel[1] ) * this->m_XDim + round ( currentVoxel[2] ) * slice ;
      //std::cout << "Visiting " << currentVoxel[0] << " " << currentVoxel[1] << " " << currentVoxel[2] << " " << index << " " ;

      this->m_StreamlineArray[index] = label ;

      if ( this->m_FAArray[index] < 0.4 ) return ;

      // find the peak diffusion direction at current voxel
      maxODFindex = this->GetPeakDir ( index ) ;
      //std::cout << "Walking along: " << maxODFindex << " "  ;

      // take a step
      if ( this->FacingAway ( maxODFindex, lastStep ) )
	{
	  maxODFindex = this->GetOpposite ( maxODFindex ) ;
	  //std::cout << "Turn around: " << maxODFindex << " " ;
	}
      lastStep = maxODFindex ;

      ODFStreamline::CoordinateType maxDir = this->m_CoordinateTable[maxODFindex] ;
      currentVoxel[0] += oneStep[0] * maxDir[0] ;  
      currentVoxel[1] += oneStep[1] * maxDir[1] ;  
      currentVoxel[2] += oneStep[2] * maxDir[2] ;  
    }
}

// scan the source file to identify all the seed voxels and store coordinates in a table for faster access
void ODFStreamline::FindSourceVoxels ()
{
  unsigned long index ;
  unsigned long slice = this->m_XDim * this->m_YDim ;
  ODFStreamline::CoordinateType currentVoxel ;
  this->m_NumberOfSourceVoxels = 0 ;
  this->m_SourceVoxels.clear () ;

  for ( unsigned long z = 0 ; z < this->m_ZDim ; z++ )
    {
      for ( unsigned long y = 0 ; y < this->m_YDim ; y++ )
	{
	  for ( unsigned long x = 0 ; x < this->m_XDim ; x++ )
	    {
	      index = x + y * this->m_XDim + z * slice ;
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
  //this->WriteSourceVoxels () ;
}

// debug
void ODFStreamline::WriteSourceVoxels () 
{
  std::cout << "Number of source voxels: " << this->m_NumberOfSourceVoxels << std::endl ;
  for ( unsigned long i = 0 ; i < this->m_NumberOfSourceVoxels ; i++ )
    {
      std::cout << this->m_SourceVoxels[i][0] << " " << this->m_SourceVoxels[i][1] << " " << this->m_SourceVoxels[i][2] << std::endl ;
    }
}

// get the index of the direction diametrically opposite to i
int ODFStreamline::GetOpposite ( unsigned int i )
{
  if ( i <= this->m_NumberOfDirs / 2 )
    return i + this->m_NumberOfDirs/2 ;

  return i - this->m_NumberOfDirs/2 ;
}

// are two directions facing away from each other (<180?)
bool ODFStreamline::FacingAway ( int i, int j )
{
  if ( i == -1 || j == -1  ) return false ;
  if ( i == -2 || j == -2 ) return true ;

  ODFStreamline::CoordinateType idir, jdir ;
  idir = this->m_CoordinateTable[i] ;
  jdir = this->m_CoordinateTable[j] ;
  double sum = idir[0] * jdir[0] + idir[1] * jdir[1] + idir[2] * jdir[2] ;
  if ( sum > 0 )
    return false ;
  else 
    return true ;
}

// get direction index that has the maximum diffusion at the current voxel
unsigned long ODFStreamline::GetPeakDir ( unsigned long index )
{
  unsigned long maxODFindex = 0 ;
  double maxODF = -1 , currentODF ;
  unsigned long indexTimesNDirs = index * this->m_NumberOfDirs ;

  for ( unsigned long dir = 0 ; dir < this->m_NumberOfDirs ; dir++ )
    {
      currentODF = this->m_ODFArray[indexTimesNDirs+dir] ;
      if ( currentODF > maxODF )
	{
	  maxODF = currentODF ;
	  maxODFindex = dir ;
	}
    }
  return maxODFindex ;
}
