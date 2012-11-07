#include <itkImageDuplicator.h>
#include <itkImageRegionIterator.h>
#include "vtkMath.h"
#include <itkImageFileWriter.h>

#include "ODFReconstructor.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

typedef itk::ImageDuplicator< ODFReconstructor::ODFImageType >    DuplicatorType;
typedef itk::ImageRegionIterator<ODFReconstructor::ODFImageType>    IteratorType;

//************************************************************************************************

ODFReconstructor::ODFReconstructor(ODFImageType::Pointer Input_image, unsigned long numberOfSamplesOnSphere, unsigned long numberOfSpharm, std::string filename )
{
   this->m_NumberOfVertices = numberOfSamplesOnSphere ;
   this->m_NumberOfSphericalHarmonics = numberOfSpharm ;
   this->m_CoefsImage = Input_image ;

   m_XDimension = Input_image->GetLargestPossibleRegion().GetSize()[0];
   m_YDimension = Input_image->GetLargestPossibleRegion().GetSize()[1];
   m_ZDimension = Input_image->GetLargestPossibleRegion().GetSize()[2];
   std::cout << "Image Dimensions: " <<  m_XDimension << " " << m_YDimension << " " << m_ZDimension << std::endl ;
   
   this->CheckOrder () ;
   this->InitializeFromFile ( filename ) ;
   
   //m_RSHBasisMatrix = new vnl_matrix< double >(numberOfSamplesOnSphere*2, m_NumberOfSphericalHarmonics) ;
   m_RSHBasisMatrix.set_size ( numberOfSamplesOnSphere * 2, m_NumberOfSphericalHarmonics ) ;
   
   //Creation of the output image 
   this->m_ODFImage = ODFImageType::New();
   
   ODFImageType::IndexType start;
   itk::VariableLengthVector< ODFType > f ( numberOfSamplesOnSphere * 2) ;

   ODFImageType::SizeType  size;
   for( unsigned int i=0; i<numberOfSamplesOnSphere * 2; i++ ) 
     { 
       f[i] = 0; 
     }
   start[0] = 0;   start[1] = 0;   start[2] = 0;
   size[0] = m_XDimension;
   size[1] = m_YDimension;
   size[2] = m_ZDimension;
   ODFImageType::SpacingType spacing;
   spacing = Input_image->GetSpacing () ;
   m_ODFImage->SetSpacing ( spacing ) ;
   ODFImageType::RegionType region;
   region.SetSize ( size ) ;
   region.SetIndex ( start ) ;
   m_ODFImage->SetVectorLength ( numberOfSamplesOnSphere * 2) ;
   m_ODFImage->SetRegions ( region ) ;
   m_ODFImage->Allocate () ;
   m_ODFImage->FillBuffer ( f ) ;

   //To manipulate easily the data
   m_CoefsArray = m_CoefsImage->GetPixelContainer()->GetBufferPointer();
   m_ODFArray = m_ODFImage->GetPixelContainer()->GetBufferPointer();
}

//************************************************************************************************

double ODFReconstructor::LegendreP(  int n, int m, double x )
{
   //Handle the case of negative m
   if (m < 0)
   {
      //Determine (-1)^(m)
      int sign = 1;
      if ( vcl_abs(m) % 2 == 1) sign = -1;

      // Compute factorial(n+m) / factorial(n-m) 
      double f = 1.0;
      // m < 0 so (n-m) > (n+m)
      // (n+m)!/(n-m)! = Prod
      for (long int i = (n-m); i>(n+m) ; i--)
      {
         f /= i;
      }
      
      return sign * f * LegendreP(n,-m,x);
   }
  
   if (m == 0 && n == 0)
   {
      return 1.0;
   }
  
   //Compute P_{n-2)_m
   //initialize with P_m^m = (-1)^m (2m-1)!! (1-x^2)^(m/2)
   double pm2 = 1;
   if ( m % 2 == 1 )
   {
      pm2 = -1;
   }

   //(-1)^m (2m-1)!!
   for ( long int i=1 ; i<2*m ; i=i+2 )
   {
      pm2 *= i;
   }
  
   pm2 *= vcl_pow(vcl_sqrt( (1+x)*(1-x) ), m);

   if (m==n) return pm2;

   //Compute P_(m+1)^m(x)   = x (2m+1)P_m^m(x).
   double pm1 = x * (2 * m + 1) * pm2;

   if (n==m+1) return pm1;
  
   // Iterate (n-m) P_n^m(x) = x(2n-1) P_(n-1)^m(x) - (n+m-1) P_(n-2)^m(x).
   double pn = 0;
   for (long int nn = m+2; nn<=n; nn++)
   {
      pn = (x * (2 * nn -1) * pm1 - (nn+m-1) * pm2) / (nn - m);
      pm2 = pm1;
      pm1 = pn;
   }
   return pn;

}

//************************************************************************************************

double ODFReconstructor::Y(long int c, double theta, double phi)
{
   LmVector vec = GetLM(c);
   const int l = vec[0];
   const int m = vec[1];

   if( m == 0 ) /// Y_l^0
      return K(l,0) * LegendreP(l,m,vcl_cos(theta));
   else if( m < 0 ) /// sqrt2 re(y_l^m)
      return vnl_math::sqrt2 * K(l,m) * vcl_cos(m*phi) * LegendreP(l,m,vcl_cos(theta));
   else ///(m > 0) sqrt2 im(y_l^m)
      return vnl_math::sqrt2* K(l,m) * vcl_sin(m*phi) * LegendreP(l,m,vcl_cos(theta));
}

//************************************************************************************************

double ODFReconstructor::K(  long int l, long int m )
{
   double f = 1; //if m=0
   if (m > 0)
   {
      for(long int i=l-m+1; i<l+m+1; i++)
      {
         f /= i; 
      }
   }
   else
   {
      for( long int i=l+m+1; i<l-m+1; i++)
      {
         f *= i; 
      }
   }
   return vcl_sqrt( ( (2*l+1) / ( 4*( vnl_math::pi ) ) * f ) );
}

//************************************************************************************************

ODFReconstructor::LmVector ODFReconstructor::GetLM(unsigned long int j)
{
   const int l = 2 * (int) ( ((1 + vcl_sqrt(8 * j - 7)) / 2) / 2);
   const int m = j - 1 - l * (l+1) / 2;
   LmVector retVal;
   
   retVal[0] = l;
   retVal[1] = m;
   
   return retVal;
}

//************************************************************************************************

ODFReconstructor::ODFImageType::Pointer ODFReconstructor::ReconstructODFImage()
{
  std::cout << "Number of vertices: " << m_NumberOfVertices << std::endl ;

  vnl_matrix < double > R (m_NumberOfVertices, m_NumberOfSphericalHarmonics ) ;

  for ( unsigned int v = 0 ; v < m_NumberOfVertices ; v++)
  {
    //Get the phi and the theta associated with the direction 
    ODFReconstructor::CoordinateType phitheta = m_PhiThetaTable[v] ;
    for ( unsigned int c = 0; c < m_NumberOfSphericalHarmonics; c++)
    {
      //Get the Matrix of the ODF
      //m_RSHBasisMatrix[v][c]  = Y ( c+1, phitheta[1], phitheta[0] ) ;
      R(v,c)  = Y ( c+1, phitheta[1], phitheta[0] ) ;
      //std::cout << m_RSHBasisMatrix[v][c] << " " ;
//      std::cout << R[v][c] << " " ;
    }
 //   std::cout << std::endl ;
  }

  std::cout << "Finished initializing RSH Basis Matrix" << std::endl ;
//  exit ( 0 ) ;

  unsigned int x, y, z, d, dir ;
  unsigned long int index, indexTimesNSPHARM, indexTimesNVerts ;
  unsigned long outputIndex ;
  unsigned long int slice = m_XDimension * m_YDimension ;
  //  bool flag = false ;
  for( z = 0 ; z < m_ZDimension ; z++ )
    {
      for( y = 0 ; y < m_YDimension ; y++ )
      {
	  for( x = 0 ; x < m_XDimension ; x++ )
	    {
	      vnl_vector< double > Coefficients ( m_NumberOfSphericalHarmonics ) ;
	      index = x + y * m_XDimension + z * slice ;
	      indexTimesNSPHARM = index * m_NumberOfSphericalHarmonics ;
	      indexTimesNVerts = index * m_NumberOfVertices ;
	      //flag = false ;
	      for( d = 0 ; d < m_NumberOfSphericalHarmonics ; d++ )
		{
		  outputIndex = d + indexTimesNSPHARM ;
		  //get the coefficients of the ODF at each voxel
		  Coefficients[d] = m_CoefsArray[outputIndex];
		  /*		  if ( fabs ( Coefficients[d] ) > 0.2 && fabs ( Coefficients[d] - 1 ) > 0.1 ) 
		    {
		      flag = true ;
		      std::cout << Coefficients[d] << " " ;
		      }*/
		}
	      //if ( flag ) std::cout << std::endl ;

	      //Multiplication of Matrix by Coefficients to get the sampled values of the ODF
	      //vnl_vector < double > Values = m_RSHBasisMatrix * Coefficients;
	      vnl_vector < double > Values = R * Coefficients;
	      double sum = 0 ;
	      for ( dir = 0 ; dir < this->m_NumberOfVertices ; dir++ )
		{
		  sum += Values[dir] ;
		}				
	      sum /= 2 * this->m_NumberOfVertices ;

  	     for( dir = 0 ; dir < this->m_NumberOfVertices ; dir++)
		{
		  outputIndex = dir + indexTimesNVerts ;
		  if ( sum != 0 ) 
		      m_ODFArray[outputIndex] = Values[dir] / sum ;
		  else
		      m_ODFArray[outputIndex] = 0 ;
		}
	    }
	}
    }
  std::cout << "Done reconstructing ODF image" << std::endl ;
  return m_ODFImage;
}

//************************************************************************************************

void ODFReconstructor::InitializeFromFile ( std::string filename ) 
{
  std::ifstream myfile;
  myfile.open ( filename.c_str() ) ;
  ODFReconstructor::CoordinateType coordinates, PhiTheta ;
  m_GlobalMin = 9999999 ;
  m_GlobalMax = m_GlobalRange = -1 ;

  unsigned long int skip = 27 ;
  for ( unsigned int i = 6 ; i < m_NumberOfVertices ; i++ )
  {
     skip += i+1 ;
  }
  char name[256] ;
  for ( unsigned long int s = 0 ; s < skip ; s++ )
  {
     myfile.getline ( name, 256 ) ;
  }

  for ( unsigned int i = 0 ; i < m_NumberOfVertices ; i++)
  {
    myfile >> coordinates[0] >> coordinates[1] >> coordinates[2];
    m_CoordinateTable.push_back ( coordinates ) ;
    // compute phi-theta
    double temp_phi = atan2 ( coordinates[1] , coordinates[0] ) ;
    while ( temp_phi >= M_PI )
    {
      temp_phi -= 2*M_PI;
    }
    while(temp_phi < -M_PI)
    {
       temp_phi += 2*M_PI;
    }
    PhiTheta[0] = temp_phi ;//phi
    PhiTheta[1] = acos ( coordinates[2] ) ; //theta		
    m_PhiThetaTable.push_back ( PhiTheta ) ;
  }
  myfile.close () ;

  this->MakeSymmetric () ;
  //this->WritePhiTheta () ;
}

void ODFReconstructor::MakeSymmetric()
{
  std::vector < ODFReconstructor::CoordinateType > symmetricPhiTheta, symmetricCoordinates ;

   // pick half		
   for ( unsigned int i = 0 ; i < m_NumberOfVertices ; i++ )
   {
      symmetricPhiTheta.push_back ( m_PhiThetaTable[i] ) ;
      symmetricCoordinates.push_back ( m_CoordinateTable[i] ) ;
   } 
 
   // create opposites
   ODFReconstructor::CoordinateType PhiTheta, Coordinates ;
   for ( unsigned int i = 0 ; i < m_NumberOfVertices ; i++ )
   {
      PhiTheta = symmetricPhiTheta[i] ;      
      PhiTheta[0] = M_PI + PhiTheta[0] ;
      PhiTheta[1] = M_PI - PhiTheta[1] ;
      symmetricPhiTheta.push_back ( PhiTheta ) ;
      Coordinates[0] = sin ( PhiTheta[1] ) * cos ( PhiTheta[0] ) ;
      Coordinates[1] = sin ( PhiTheta[1] ) * sin ( PhiTheta[0] ) ;
      Coordinates[2] = cos ( PhiTheta[1] ) ;
      symmetricCoordinates.push_back ( Coordinates ) ;
   }

   m_NumberOfVertices *= 2 ;

   // update tables
   m_PhiThetaTable.clear () ;
   m_CoordinateTable.clear () ;
   for ( unsigned int i = 0 ; i < m_NumberOfVertices ; i++ )
   {
      m_PhiThetaTable.push_back ( symmetricPhiTheta[i] ) ;
      m_CoordinateTable.push_back ( symmetricCoordinates[i] ) ;
   }
}

void ODFReconstructor::WritePhiTheta () 
{
  std::cout << "Phi-Theta Table Dump" << std::endl ;
  for ( unsigned int i = 0 ; i < m_NumberOfVertices ; i++ )
    {
      std::cout << i << " " << m_PhiThetaTable[i][0] << " " << m_PhiThetaTable[i][1] << std::endl ;
    }
}

void ODFReconstructor::CheckOrder()
{
  int order = 2 * (int) (((1 + vcl_sqrt(8 * this->m_NumberOfSphericalHarmonics - 7)) / 2) / 2);

  switch ( order )
  {
    case 0:    case 2:    case 4:    case 6:    case 8:
    case 10:    case 12:    case 14:    case 16:    case 18:
    case 20:
    {
      std::cout << "RSH order check passed" << std::endl ;
       return ;
       break;
    }
   default:
   {
      std::cout << "Unsupported RSH ORDER : " << m_NumberOfSphericalHarmonics << std::endl;
      exit ( -1 ) ;
      //EXCEPTION!!!
   }
  }
}

ODFReconstructor::ImageType::Pointer ODFReconstructor::GetFAImage () 
{
   ODFReconstructor::ImageType::Pointer FAImage = ODFReconstructor::ImageType::New();
   
   ODFReconstructor::ImageType::IndexType start;
   ODFReconstructor::ImageType::SizeType  size;
   start[0] = 0;   start[1] = 0;   start[2] = 0;
   size[0] = m_XDimension;
   size[1] = m_YDimension;
   size[2] = m_ZDimension;
   ODFReconstructor::ImageType::SpacingType spacing;
   spacing = this->m_CoefsImage->GetSpacing () ;
   FAImage->SetSpacing ( spacing ) ;
   ODFReconstructor::ImageType::RegionType region;
   region.SetSize ( size ) ;
   region.SetIndex ( start ) ;
   FAImage->SetRegions ( region ) ;
   FAImage->Allocate () ;
      
   //To manipulate easily the data
   double *FAArray = FAImage->GetPixelContainer()->GetBufferPointer();

   unsigned long slice = m_XDimension * m_YDimension ;
   unsigned long index, indexTimesNVerts ;

   for ( unsigned long z = 0 ; z < m_ZDimension ; z++ )
   {
     for ( unsigned long y = 0 ; y < m_YDimension ; y++ )
       {
	 for ( unsigned long x = 0 ; x < m_XDimension ; x++ )
	   {
	     index = x + y * m_XDimension + z * slice ;
	     indexTimesNVerts = index * m_NumberOfVertices ;
	     double max = -1 ;
	     for ( unsigned long v = 0 ; v < m_NumberOfVertices ; v++ )
	       {
		 double current = this->m_ODFArray[indexTimesNVerts+v] ;
		 if ( current > max )
		   max = current ;
		 if ( current > m_GlobalMax ) 
		   m_GlobalMax = current ;
		 if ( current < m_GlobalMin ) 
		   m_GlobalMin = current ;
	       }
	     FAArray[index] = max ;
	   }
       }
   }

   // normalize
   index = 0 ;
   m_GlobalRange = m_GlobalMax - m_GlobalMin ;

   for ( unsigned int z = 0 ; z < m_ZDimension ; z++ )
   {
     for ( unsigned int y = 0 ; y < m_YDimension ; y++ )
       {
	 for ( unsigned int x = 0 ; x < m_XDimension ; x++ )
	   {
	     double current = FAArray[index] ;
	     FAArray[index] = ( current - m_GlobalMin ) / m_GlobalRange ;
	     index++ ;
	   }
       }
   }
  std::cout << m_GlobalMax << " " << m_GlobalMin << " " << m_GlobalRange << std::endl ;
  return FAImage ;
}

ODFReconstructor::RGBImageType::Pointer ODFReconstructor::GetColorFAImage () 
{
  if ( m_GlobalMax == -1 )
    this->GetFAImage () ;

   ODFReconstructor::RGBImageType::Pointer FAImage = ODFReconstructor::RGBImageType::New();
   
   ODFReconstructor::RGBImageType::IndexType start;
   ODFReconstructor::RGBImageType::SizeType  size;
   start[0] = 0;   start[1] = 0;   start[2] = 0;
   size[0] = m_XDimension;
   size[1] = m_YDimension;
   size[2] = m_ZDimension;
   ODFReconstructor::RGBImageType::SpacingType spacing;
   spacing = this->m_CoefsImage->GetSpacing () ;
   FAImage->SetSpacing ( spacing ) ;
   ODFReconstructor::RGBImageType::RegionType region;
   region.SetSize ( size ) ;
   region.SetIndex ( start ) ;
   FAImage->SetRegions ( region ) ;
   FAImage->Allocate () ;
      
   //To manipulate easily the data
   ODFReconstructor::RGBPixelType *FAArray = FAImage->GetPixelContainer()->GetBufferPointer() ;

   unsigned long slice = m_XDimension * m_YDimension ;
   unsigned long index, indexTimesNVerts ;

   for ( unsigned int z = 0 ; z < m_ZDimension ; z++ )
   {
     for ( unsigned int y = 0 ; y < m_YDimension ; y++ )
       {
	 for ( unsigned int x = 0 ; x < m_XDimension ; x++ )
	   {
	     index = x + y * m_XDimension + z * slice ;
	     indexTimesNVerts = index * m_NumberOfVertices ;
	     double max = -1 ;
	     int maxIndex = -1 ;
	     for ( unsigned int v = 0 ; v < m_NumberOfVertices ; v++ )
	       {
		 double current = this->m_ODFArray[indexTimesNVerts+v] ;
		 if ( current > max )
		   {
		     max = current ;
		     maxIndex = v ;
		   }
	       }
	     FAArray[index] = this->DirToRGB ( maxIndex, max ) ;
	   }
       }
   }
  return FAImage ;
}

ODFReconstructor::RGBPixelType ODFReconstructor::DirToRGB ( unsigned long dirIndex, double fa )
{
  RGBPixelType rgb ;
  CoordinateType direction = this->m_CoordinateTable[dirIndex] ;
  fa = ( fa - m_GlobalMin ) * m_GlobalRange * 192 ;
  rgb[0] = fabs(direction[0]) * fa ;
  rgb[1] = fabs(direction[1]) * fa ;
  rgb[2] = fabs(direction[2]) * fa ;

  return rgb ;
}

