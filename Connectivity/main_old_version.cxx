#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>

#include <iostream>

#include "SphereIkosahedronImage.h"
#include "Fstar_prob3D.h"
#include "Fstar_image_creator.h"
#include "Fstar_cost3D.h"

//Vector image definition
typedef double 						CostType;
typedef SphereIkosahedronImage<CostType> 		VectorImageType;
typedef itk::ImageFileReader < VectorImageType > 	VectorReaderType;
typedef itk::ImageFileWriter < VectorImageType > 	VectorWriterType;

//Normal image definition
typedef double 					ImageElementType;
typedef itk::Image < ImageElementType , 3 > 	ImageType;
typedef itk::ImageFileReader < ImageType > 	ReaderType;

int main(int argc, char * argv[])
{
  std::cout << "Usage: " << argv[0] << " basis_image output_image source_label_image icosahedron_subdivision threshold_coefficient vector_or_odf prob_or_cost alpha" << std::endl ;
  std::cout << "vector_or_odf: 0 for vector (synthethic), 1 for odf" << std::endl ;
  std::cout << "prob_or_cost: 0 for probability, 1 for cost" << std::endl ;
  std::cout << "alpha: weight of odf versus angle similarity" << std::endl ;
  if ( argc != 8 ) return 0 ;
  //Vector image pointers
  VectorImageType::Pointer VectorBasisImage 	= VectorImageType::New();
  VectorImageType::Pointer VectorInputImage 	= VectorImageType::New();
  VectorReaderType::Pointer VectorReader       	= VectorReaderType::New();
  VectorWriterType::Pointer VectorWriter 	= VectorWriterType::New();
  VectorImageType::Pointer VectorOutputImage 	= VectorImageType::New();
	
  //Normal image pointers
  ImageType::Pointer InputLabelImage 		= ImageType::New();
  ImageType::Pointer BasisImage        		= ImageType::New();
  ReaderType::Pointer Reader_1 			= ReaderType::New();
  ReaderType::Pointer Reader_2 			= ReaderType::New();
	
  const char * file_basis_image = argv [1];//Input file
  const char * file_final = argv[2];//Output file
  const char * file_input_label = argv[3];//Label file
	
  short icosahedronSubdivision = atoi ( argv[4] ) ;//Icosahedron subdivision level
  double thresholdCoefficient = atof ( argv[5] ) ;
  short vector_or_odf = atoi ( argv[6] ) ;
  short prob_or_cost = atoi ( argv[7] ) ;
  //double alpha = atof ( argv[8] ) ;
  //double step_cost = atof ( argv[9] ) ;

  Reader_1->SetFileName(file_input_label);
  try
    {
      Reader_1->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Problem reading the input 'basis' file" << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
  InputLabelImage = Reader_1->GetOutput();
	
  if( vector_or_odf )
    {
      std::cout << "odf" << std::endl ;
      VectorReader->SetFileName ( file_basis_image );
      try
	{
	  VectorReader->Update();
	}
      catch( itk::ExceptionObject & excp )
	{
	  std::cerr << "Problem reading the input 'basis' file" << std::endl;
	  std::cerr << excp << std::endl;
	  return EXIT_FAILURE;
	}
      VectorInputImage = VectorReader->GetOutput();
      Fstar_image_creator Fstar_image_creator_object(VectorInputImage, icosahedronSubdivision, thresholdCoefficient);//make the ODF usable
      VectorBasisImage = Fstar_image_creator_object.CreateVectorImage();
    }
  else
    {	
      std::cout << "synthetic" << std::endl ;
      Reader_2->SetFileName ( file_basis_image );
      try
	{
	  Reader_2->Update();
	}
      catch( itk::ExceptionObject & excp )
	{
	  std::cerr << "Problem reading the input 'basis' file" << std::endl;
	  std::cerr << excp << std::endl;
	  return EXIT_FAILURE;
	}
      BasisImage = Reader_2->GetOutput();
      Fstar_image_creator Fstar_image_creator_object(BasisImage, icosahedronSubdivision, thresholdCoefficient);//compute a synthesis vector image
      Fstar_image_creator_object.Image_creator();
      VectorBasisImage = Fstar_image_creator_object.GetOutputImage();
      //print out min/max here
    }
  exit(0) ;
  // choose and run main algorithm
  if ( prob_or_cost == 0) 
    {
      std::cout << "Probability" << std::endl ;
      Fstar_prob3D Fstar_prob3D_object(VectorBasisImage, InputLabelImage, icosahedronSubdivision, thresholdCoefficient) ;
      Fstar_prob3D_object.GetProbMap ( file_final ) ;
    }
  else 
    {
      std::cout << "Cost" << std::endl ;
      Fstar_cost3D Fstar_cost3D_object ( VectorBasisImage, InputLabelImage, icosahedronSubdivision, thresholdCoefficient ) ;

      std::string mincost = file_final ;
      mincost.insert ( mincost.rfind ( "." ) , "_mincost" ) ;
      std::string lengthmap = file_final ;
      lengthmap.insert ( lengthmap.rfind ( "." ), "_length" ) ;
      std::string origmap = file_final ;
      origmap.insert ( origmap.rfind ( "." ), "_orig" ) ;
      std::string avgcostmap = file_final ;
      avgcostmap.insert ( avgcostmap.rfind ( "." ), "_avgcost" ) ;
      std::string mincostperminlengthmap = file_final ;
      mincostperminlengthmap.insert ( mincostperminlengthmap.rfind ( "." ), "_mincost_minlength" ) ;

      Fstar_cost3D_object.GetCostMap ( mincost, lengthmap, origmap, avgcostmap, mincostperminlengthmap ) ;
    }
	
  return 0;
}
