#include <itkImageDuplicator.h>
#include <itkImageRegionIterator.h>
#include "vtkMath.h"
#include <itkImageFileWriter.h>

#include "Fstar_image_creator.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

using namespace std;

typedef itk::ImageDuplicator< Fstar_image_creator::VectorImageType > 	DuplicatorType;
typedef itk::ImageRegionIterator<Fstar_image_creator::VectorImageType> 	IteratorType;
typedef itk::ImageFileWriter < Fstar_prob3D::ImageType > 		WriterType;

//************************************************************************************************
// For synthetic image creation
// if label = 1, we set probability to one to directions pointing to the superior
// if label = 2, we set probability to one to directions pointing to the posterior
// if label = 4, we set probability to one to directions pointing to the left
//************************************************************************************************

Fstar_image_creator::Fstar_image_creator(ImageType::Pointer Input_image, unsigned long m_SubdivisionLevel, double thresholdCoefficient)//with synthetic labeled image
{
	m_XDimension = Input_image->GetLargestPossibleRegion().GetSize()[0];
	m_YDimension = Input_image->GetLargestPossibleRegion().GetSize()[1];
	m_ZDimension = Input_image->GetLargestPossibleRegion().GetSize()[2];
	m_SLICE = m_XDimension*m_YDimension;
	m_LINE = m_XDimension;

	m_ThresholdCoefficient = thresholdCoefficient;
	
	m_SphereIkosahedronObject = SphereIkosahedron<double>::New();
	
	m_SphereIkosahedronObject->SetSubdivisionLevel(m_SubdivisionLevel);
	m_SphereIkosahedronObject->Initialisation_tables();
	m_NumberOfVertices = m_SphereIkosahedronObject->GetNumberOfVertices();
	
	m_IkoOutputImage = VectorImageType::New();
	
	VectorImageType::IndexType start;
	itk::VariableLengthVector< ImageElementType > f( m_NumberOfVertices );
	VectorImageType::SizeType  size;
	for( unsigned long int i=0; i<m_NumberOfVertices; i++ ) { f[i] = 0; }
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	size[0] = m_XDimension;
	size[1] = m_YDimension;
	size[2] = m_ZDimension;
	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	m_IkoOutputImage->SetVectorLength( m_NumberOfVertices );
	m_IkoOutputImage->SetRegions( region );
	m_IkoOutputImage->Allocate();
	m_IkoOutputImage->FillBuffer( f );
	
	m_LabeledImageArray = Input_image->GetPixelContainer()->GetBufferPointer();
	m_IkoOutputArray = m_IkoOutputImage->GetPixelContainer()->GetBufferPointer();
}

//************************************************************************************************

Fstar_image_creator::Fstar_image_creator(VectorImageType::Pointer Input_image, unsigned long m_SubdivisionLevel, double thresholdCoefficient)//with ODF
{
	m_XDimension = Input_image->GetLargestPossibleRegion().GetSize()[0];
	m_YDimension = Input_image->GetLargestPossibleRegion().GetSize()[1];
	m_ZDimension = Input_image->GetLargestPossibleRegion().GetSize()[2];
	std::cout << m_XDimension << " " << m_YDimension << " " << m_ZDimension << std::endl ;
	m_SLICE = m_XDimension*m_YDimension;
	m_LINE = m_XDimension;
	
	m_ThresholdCoefficient = thresholdCoefficient;
	
	m_NumberOfSphericalHarmonics = 15;
	m_Order = 2 * (int) (((1 + vcl_sqrt(8 * m_NumberOfSphericalHarmonics - 7)) / 2) / 2);
	std::cout << "Order: " << m_Order << std::endl ;

	m_SphereIkosahedronObject = SphereIkosahedron<double>::New();
	
	m_SphereIkosahedronObject->SetSubdivisionLevel(m_SubdivisionLevel);
	m_SphereIkosahedronObject->Initialisation_tables();//Icosahedron subdivision
	m_NumberOfVertices = m_SphereIkosahedronObject->GetNumberOfVertices();
	
	m_RSHBasisMatrix.set_size(m_NumberOfVertices, m_NumberOfSphericalHarmonics);
	
	//Creation of the output image 
	m_IkoOutputImage = VectorImageType::New();
	
	VectorImageType::IndexType start;
	itk::VariableLengthVector< ImageElementType > f( m_NumberOfVertices );
	VectorImageType::SizeType  size;
	for( unsigned int i=0; i<m_NumberOfVertices; i++ ) { f[i] = 0; }
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	size[0] = m_XDimension;
	size[1] = m_YDimension;
	size[2] = m_ZDimension;
	VectorImageType::SpacingType spacing;
	spacing = Input_image->GetSpacing();
	m_IkoOutputImage->SetSpacing(spacing);
	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	m_IkoOutputImage->SetVectorLength( m_NumberOfVertices );
	m_IkoOutputImage->SetRegions( region );
	m_IkoOutputImage->Allocate();
	m_IkoOutputImage->FillBuffer( f );
	
	//To manipulate easily the data
	m_ODFImageArray = Input_image->GetPixelContainer()->GetBufferPointer();
	m_IkoOutputArray = m_IkoOutputImage->GetPixelContainer()->GetBufferPointer();
	
	//Creation of the sum image (visualisation)
	m_SumImage = ImageType::New();
	ImageType::IndexType startsum;
	ImageType::SizeType  sizesum;
	startsum[0] = 0;
	startsum[1] = 0;
	startsum[2] = 0;
	sizesum[0] = m_XDimension;
	sizesum[1] = m_YDimension;
	sizesum[2] = m_ZDimension;
	ImageType::RegionType regionsum;
	regionsum.SetSize( sizesum );
	regionsum.SetIndex( startsum );
	m_SumImage->SetRegions( regionsum );
	m_SumImage->Allocate();
	
	//Creation of the unique coefficient image (visualisation)
	m_OneCoefficientImage = ImageType::New();
	ImageType::IndexType startOneCoefficient;
	ImageType::SizeType  sizeOneCoefficient;
	startOneCoefficient[0] = 0;
	startOneCoefficient[1] = 0;
	startOneCoefficient[2] = 0;
	sizeOneCoefficient[0] = m_XDimension;
	sizeOneCoefficient[1] = m_YDimension;
	sizeOneCoefficient[2] = m_ZDimension;
	ImageType::RegionType regionOneCoefficient;
	regionOneCoefficient.SetSize( sizeOneCoefficient );
	regionOneCoefficient.SetIndex( startOneCoefficient );
	m_OneCoefficientImage->SetRegions( regionOneCoefficient );
	m_OneCoefficientImage->Allocate();
	
	m_ColorFAImage = RGBImageType::New () ;
	RGBImageType::IndexType startColorFA ;
	RGBImageType::SizeType sizeColorFA ;
	startColorFA[0] = startColorFA[1] = startColorFA[2] = 0 ;
	sizeColorFA[0] = m_XDimension ;
	sizeColorFA[1] = m_YDimension ;
	sizeColorFA[2] = m_ZDimension ;
	RGBImageType::RegionType regionColorFA ;
	regionColorFA.SetSize ( sizeColorFA ) ;
	regionColorFA.SetIndex ( startColorFA ) ;
	m_ColorFAImage->SetRegions ( regionColorFA ) ;
	m_ColorFAImage->Allocate () ;

	m_SumImageArray = m_SumImage->GetPixelContainer()->GetBufferPointer();
	m_OneCoefficientImageArray = m_OneCoefficientImage->GetPixelContainer()->GetBufferPointer();
	m_ColorFAImageArray = m_ColorFAImage->GetPixelContainer()->GetBufferPointer() ;
}

//************************************************************************************************

double Fstar_image_creator::LegendreP(  int n, int m, double x )
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

double Fstar_image_creator::Y(long int c, double theta, double phi)
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

double Fstar_image_creator::K(  long int l, long int m )
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

Fstar_image_creator::LmVector Fstar_image_creator::GetLM(unsigned long int j)
{
	const int l = 2 * (int) ( ((1 + vcl_sqrt(8 * j - 7)) / 2) / 2);
	const int m = j - 1 - l * (l+1) / 2;
	LmVector retVal;
	
	retVal[0] = l;
	retVal[1] = m;
	
	return retVal;
}

//************************************************************************************************

double Fstar_image_creator::EvaluateBasis( long int c, double theta, double phi)
{
	switch ( m_Order )
	{
		case 0:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 2:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 4:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 6:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 8:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 10:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 12:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 14:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 16:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 18:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		case 20:
		{
			return Y(c+1,static_cast<double>(theta),static_cast<double>(phi));
			break;
		}
		default:
		{
			std::cout << "Unsupported RSH ORDER : " << m_Order << std::endl;
			return -1.0;
      			//EXCEPTION!!!
		}
	}
}

//************************************************************************************************

Fstar_image_creator::VectorImageType::Pointer Fstar_image_creator::CreateVectorImage()
{
  std::cout << "Number of vertices: " << m_NumberOfVertices << std::endl ;
  for(m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
    {
      //Get the phi and the theta associated with the direction 
      std::vector<double> phitheta = m_SphereIkosahedronObject->GetPhiThetaTableatIndex(m_Direction) ;
      for ( long int c = 0; c < m_NumberOfSphericalHarmonics; c++)
	{
	  //Get the Matrix of the ODF
	  m_RSHBasisMatrix[m_Direction][c]  = EvaluateBasis(c, phitheta[1], phitheta[0]);
	}
    }

  unsigned long int zbase, ybase, xbase ;
  unsigned long int zbase2, ybase2, xbase2 ;
  unsigned long int zindex, yindex, xindex ;

  const unsigned long int SliceTimesNSphericalHarmonics = m_SLICE * m_NumberOfSphericalHarmonics ;
  const unsigned long int SliceTimesNVertices = m_SLICE * m_NumberOfVertices ;
  const unsigned long int LineTimesNSphericalHarmonics = m_LINE * m_NumberOfSphericalHarmonics ;
  const unsigned long int LineTimesNVertices = m_LINE * m_NumberOfVertices ;

  zbase = zbase2 = zindex = 0 ;
  for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      yindex = zindex ;
      ybase = zbase ;
      ybase2 = zbase2 ;
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  xindex = yindex ;
	  xbase = ybase ;
	  xbase2 = ybase2 ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      vnl_vector< double > Coefficients(m_NumberOfSphericalHarmonics) ;
	      for(  long int d = 0 ; d < m_NumberOfSphericalHarmonics ; d++)
		{
		  m_OutputIndex = d + xbase ;
		  //get the coefficients of the ODF at each voxel
		  Coefficients[d] = m_ODFImageArray[m_OutputIndex];
		}
	      
	      //Multiplication of Matrix by Coefficients to get the sampled values of the ODF
	      vnl_vector< double > Values = m_RSHBasisMatrix * Coefficients;
				
	      //Visualisation
	      double max = -1 ;	
	      unsigned long maxIndex = -1 ;
	     for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
		{
		  m_OutputIndex = m_Direction + xbase2 ;
		  //Fill the output image with the sampled ODF values
		  m_IkoOutputArray[m_OutputIndex] = Values[m_Direction];
		  if ( Values[m_Direction] > max ) 
		    {
		      max = Values[m_Direction] ;
		      maxIndex = m_Direction ;
		    }
		  }

	      m_OneCoefficientImageArray[xindex] = max;
	      RGBPixelType rgb ;	
	      vector < double > maxDirection = m_SphereIkosahedronObject->GetCoordinateTableatIndex(maxIndex) ;

	      rgb[0] = fabs(maxDirection[0]) * 255 * max ;
	      rgb[1] = fabs(maxDirection[1]) * 255 * max ;
	      rgb[2] = fabs(maxDirection[2]) * 255 * max ;
	      	      
	      m_ColorFAImageArray[xindex] = rgb ;
	      xindex++ ;
	      xbase += m_NumberOfSphericalHarmonics ;
	      xbase2 += m_NumberOfVertices ;
	    }

	  yindex += m_LINE ;
	  ybase += LineTimesNSphericalHarmonics ;
	  ybase2 += LineTimesNVertices ;
	}
      zindex += m_SLICE ;
      zbase += SliceTimesNSphericalHarmonics ;
      zbase2 += SliceTimesNVertices ;
    }

	WriterType::Pointer CoeffWriter = WriterType::New();
	cout<<"Writing the image of the coefficient 0 of ODF at each voxel."<<endl;
	CoeffWriter->SetFileName ( "../images/PseudoFA.nrrd" );
	CoeffWriter->SetInput( m_OneCoefficientImage );

        RGBImageWriterType::Pointer RGBimageWriter = RGBImageWriterType::New();	
	RGBimageWriter->SetFileName("../images/PseudoColorFA.nrrd");
	RGBimageWriter->SetInput( m_ColorFAImage );
	RGBimageWriter->Update();

	try
	{
		CoeffWriter->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		cerr<<"Problem writing the output file"<<endl;
		cerr<<excp<<endl;
	}
    
  //First Normalization
  //Uniforme probability along each direction
  const double ProbUniform = 2.0 / double(m_NumberOfVertices);
  cout<<"ProbUniform = "<<ProbUniform<<endl;
  //To information display only
  double max = 0;

  double temp ;
  zindex = zbase = 0 ;
  for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      yindex = zindex ;
      ybase = zbase ;
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  xindex = yindex ;
	  xbase = ybase ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      double NormFactor = 0;
	      //Sum the values for each direction at each voxel
	      for(m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
		{
		  m_OutputIndex = m_Direction + xbase;
		  NormFactor += m_IkoOutputArray[m_OutputIndex];
		  if(m_IkoOutputArray[m_OutputIndex] > max)
		    {
		      max = m_IkoOutputArray[m_OutputIndex];
		    }
		}
	      NormFactor /= 2;
	      if(NormFactor > 0)
		{
		  for(m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
		    {
		      m_OutputIndex = m_Direction + xbase ;
		      //Normalisation of each direction (so that the sum = 2, or = 1 on each half sphere)
		      m_IkoOutputArray[m_OutputIndex] /= NormFactor;
		      //Change values to be centred on 1
		      temp = 1 + (m_IkoOutputArray[m_OutputIndex] - ProbUniform)/ProbUniform;
		      //Increase the gap between < 1 and > 1 values to sharpen the ODF
		      m_IkoOutputArray[m_OutputIndex] = temp*temp*temp*temp*temp*temp;
						
		      //If value is < to the threshold, we don't deal with it (set to 0)
		      if( m_IkoOutputArray[m_OutputIndex] < m_ThresholdCoefficient )
			{
			  m_IkoOutputArray[m_OutputIndex] = 0.0;
			}
		    }
		}
	      xindex++ ;
	      xbase += m_NumberOfVertices ;
	    }
	  yindex += m_LINE ;
	  ybase += LineTimesNVertices ;
	}
      zindex += m_SLICE ;
      zbase += SliceTimesNVertices ;
    }
  cout<<"Maximum is "<<max<<endl;

  //Second normalization
  zbase = 0 ;
  for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
    {
      ybase = zbase ;
      for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
	{
	  xbase = ybase ;
	  for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
	    {
	      double NormFactor = 0;

	      for(m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
		{
		  m_OutputIndex = m_Direction + xbase ;
					
		  NormFactor += m_IkoOutputArray[m_OutputIndex];
		}
	      if(NormFactor > 0)
		{
		  NormFactor /= 2;
		  for(m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
		    {
		      m_OutputIndex = m_Direction + xbase ;
		      m_IkoOutputArray[m_OutputIndex] /= NormFactor;
		    }
		}
	      xbase += m_NumberOfVertices ;
	    }
	  ybase += LineTimesNVertices ;
	}
      zbase += SliceTimesNVertices ;
    }
	
  ComputeODFSumImage();
	
  return m_IkoOutputImage;
}

//************************************************************************************************

double Fstar_image_creator::GetDirectionSimilarity(VectorType StaticVector , unsigned long int m_Direction)
{
	//dot product is v1.v2 = |v1|*|v2|*cos(angle)
	//We want the angle
	std::vector<double> v1 = m_SphereIkosahedronObject->GetCoordinateTableatIndex(m_Direction);//Cartesian direction vector
	
	double anglesimilarity;
	double x1, y1, z1, x2, y2, z2;
	x1 = v1[0];
	y1 = v1[1];
	z1 = v1[2];
	x2 = StaticVector[0];
	y2 = StaticVector[1];
	z2 = StaticVector[2];
	
	//Norm of the vectors
	double v1_magnitude = sqrt(x1*x1 + y1*y1 + z1*z1);
	double v2_magnitude = sqrt(x2*x2 + y2*y2 + z2*z2);
	
	//Normalisation
	x1 /= v1_magnitude;
	x2 /= v2_magnitude;
	y1 /= v1_magnitude;
	y2 /= v2_magnitude;
	z1 /= v1_magnitude;
	z2 /= v2_magnitude;
	
	//Get the angle between the two vectors (directions)
	anglesimilarity = acos(x1*x2 + y1*y2 + z1*z2);
	
	return anglesimilarity;
}

//************************************************************************************************

double Fstar_image_creator::GetAngleCoefficient(double AnglesSimilarityTemp, double mu, double sigma)
{
	double anglecoefficient, tempresult;
	
	tempresult = ( mu - AnglesSimilarityTemp )/sigma;
	anglecoefficient = 0.25*( 1 + erf( tempresult ));
	
	return anglecoefficient;
}

//************************************************************************************************

void Fstar_image_creator::Image_creator()
{
	//Creation of reference vectors
	VectorType Vertical_up, Vertical_down;
	Vertical_up[0] = 0; Vertical_up[1] = 0; Vertical_up[2] = 1;
	Vertical_down[0] = 0; Vertical_down[1] = 0; Vertical_down[2] = -1;
	
	VectorType Horiz_right, Horiz_left;
	Horiz_right[0] = -1; Horiz_right[1] = 0; Horiz_right[2] = 0;
	Horiz_left[0] = 1; Horiz_left[1] = 0; Horiz_left[2] = 0;
	
	VectorType Horiz_posterior, Horiz_anterior;
	Horiz_posterior[0] = 0; Horiz_posterior[1] = 1; Horiz_posterior[2] = 0;
	Horiz_anterior[0] = 0; Horiz_anterior[1] = -1; Horiz_anterior[2] = 0;
	
	for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
	{
		for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
		{	
			for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
			{
				m_CurrentIndex = m_x_it + m_y_it*m_LINE + m_z_it*m_SLICE;
				m_CurrentFeature = m_LabeledImageArray[m_CurrentIndex];//Label of the voxel
				
				//Details of the labels are given at the beginning of the file
				if( m_CurrentFeature == 1 )
				{
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						
						double AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_up , m_Direction);//angle between direction and ref vector
						double anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);//Apply the Angle function
						
						if(anglecoeff > 0)
						{
							//Direction are similar
							m_IkoOutputArray[m_OutputIndex] = anglecoeff;
						}
						else
						{
							//Direction are not similar, we test with opposite direction
							AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_down , m_Direction);
							m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						}
					}
				}
				else if( m_CurrentFeature == 2 )
				{
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						
						double AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_posterior , m_Direction);
						double anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
						if(anglecoeff > 0)
						{
							m_IkoOutputArray[m_OutputIndex] = anglecoeff;
						}
						else
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_anterior , m_Direction);
							m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						}
					}
				}
				else if( m_CurrentFeature == 3 )//Both label 1 and 2
				{
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						
						double AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_posterior , m_Direction);
						double anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
						if(anglecoeff > 0)
						{
							m_IkoOutputArray[m_OutputIndex] = anglecoeff;
						}
						else
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_anterior , m_Direction);
							m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						}
						if(m_IkoOutputArray[m_OutputIndex] <= 0.00001)
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_up , m_Direction);
							anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
							if(anglecoeff > 0)
							{
								m_IkoOutputArray[m_OutputIndex] = anglecoeff;
							}
							else
							{
								AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_down , m_Direction);
								m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
							}
						}
					}
				}
				else if( m_CurrentFeature == 4 )
				{
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						
						double AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_left , m_Direction);
						double anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
						if(anglecoeff > 0)
						{
							m_IkoOutputArray[m_OutputIndex] = anglecoeff;
						}
						else
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_right , m_Direction);
							m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						}
					}
				}
				else if( m_CurrentFeature == 5 )//Both label 1 and 4
				{
					//cout<<"m_z_it = "<<m_z_it<<" // m_y_it = "<<m_y_it<<" // m_x_it = "<<m_x_it<<endl;
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						
						double AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_left , m_Direction);
						double anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
						if(anglecoeff > 0)
						{
							m_IkoOutputArray[m_OutputIndex] = anglecoeff;
						}
						else
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_right , m_Direction);
							m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						}
						if(m_IkoOutputArray[m_OutputIndex] <= 0.00001)
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_up , m_Direction);
							anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
							
							if(anglecoeff > 0)
							{
								m_IkoOutputArray[m_OutputIndex] = anglecoeff;
							}
							else
							{
								AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_down , m_Direction);
								m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
							}
						}
					}
				}
				else if( m_CurrentFeature == 6 )//Both label 2 and 4
				{
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						
						double AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_posterior , m_Direction);
						double anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
						if(anglecoeff > 0)
						{
							m_IkoOutputArray[m_OutputIndex] = anglecoeff;
						}
						else
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_anterior , m_Direction);
							m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						}
						if(m_IkoOutputArray[m_OutputIndex] <= 0.00001)
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_left , m_Direction);
							anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
							if(anglecoeff > 0)
							{
								m_IkoOutputArray[m_OutputIndex] = anglecoeff;
							}
							else
							{
								AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_right , m_Direction);
								m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
							}
						}
					}
				}
				else if( m_CurrentFeature == 7 )//Labels 1, 2 and 4
				{
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						
						double AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_posterior , m_Direction);
						double anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
						if(anglecoeff > 0)
						{
							m_IkoOutputArray[m_OutputIndex] = anglecoeff;
						}
						else
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_anterior , m_Direction);
							m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						}
						if(m_IkoOutputArray[m_OutputIndex] <= 0.00001)
						{
							AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_left , m_Direction);
							anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
							if(anglecoeff > 0)
							{
								m_IkoOutputArray[m_OutputIndex] = anglecoeff;
							}
							else
							{
								AnglesSimilarityTemp = GetDirectionSimilarity( Horiz_right , m_Direction);
								m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
							}
							if(m_IkoOutputArray[m_OutputIndex] <= 0.00001)
							{
								AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_up , m_Direction);
								anglecoeff = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
						
								if(anglecoeff > 0)
								{
									m_IkoOutputArray[m_OutputIndex] = anglecoeff;
								}
								else
								{
									AnglesSimilarityTemp = GetDirectionSimilarity( Vertical_down , m_Direction);
									m_IkoOutputArray[m_OutputIndex] = GetAngleCoefficient(AnglesSimilarityTemp, M_PI/8, M_PI/16);
								}
							}
						}
					}
				}
				else//If no lqbel, set each direction to 0
				{
					for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
					{
						m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
						m_IkoOutputArray[m_OutputIndex] = 0;
					}
				}
			}
		}
	}
}

//************************************************************************************************

void Fstar_image_creator::ComputeODFSumImage()
{
	const char * sum_file;
	
	for( m_z_it = 0 ; m_z_it < m_ZDimension ; m_z_it++ )
	{
		for( m_y_it = 0 ; m_y_it < m_YDimension ; m_y_it++ )
		{
			for( m_x_it = 0 ; m_x_it < m_XDimension ; m_x_it++ )
			{
				m_CurrentIndex = m_x_it + m_y_it*m_LINE + m_z_it*m_SLICE;
				double sumatvoxel = 0;
				
				for( m_Direction = 0 ; m_Direction < m_NumberOfVertices ; m_Direction++)
				{
					m_OutputIndex = m_Direction + m_x_it*m_NumberOfVertices + m_y_it*m_LINE*m_NumberOfVertices + m_z_it*m_SLICE*m_NumberOfVertices;
					sumatvoxel += m_IkoOutputArray[m_OutputIndex];
				}
				
				m_SumImageArray[m_CurrentIndex] = sumatvoxel;//Visualisation
			}
		}
	}
	WriterType::Pointer Writer = WriterType::New();
	sum_file = "../images/SumODF.nrrd";
	cout<<"Writing the image of the sum of ODF at each voxel."<<endl;
	Writer->SetFileName ( sum_file );
	Writer->SetInput( m_SumImage );
	
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

//************************************************************************************************

Fstar_image_creator::VectorImageType::Pointer Fstar_image_creator::GetOutputImage()
{
	return m_IkoOutputImage;
}

//************************************************************************************************
