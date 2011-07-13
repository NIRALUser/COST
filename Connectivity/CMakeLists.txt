cmake_minimum_required(VERSION 2.6)

project(Prob_Project)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.") 
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
   MESSAGE(FATAL_ERROR "VTK not found. Please set VTK_DIR.") 
ENDIF(VTK_FOUND)

IF (USE_VTK_FILE)
   INCLUDE(${USE_VTK_FILE})
ENDIF (USE_VTK_FILE)

SET(TXX_FILES 
  SphereIkosahedron.txx
  SphereIkosahedronImage.txx
 )

SET(CXX_FILES 
  main.cxx
  Fstar_prob3D.cxx
  Fstar_cost3D.cxx
  Fstar_image_creator.cxx
)

SET(H_FILES
  SphereIkosahedron.h
  Fstar_prob3D.h
  SphereIkosahedronImage.h
  Fstar_image_creator.h
)

add_executable( 
	Prob_Project
	${TXX_FILES}
	${CXX_FILES}
	${H_FILES}
)
TARGET_LINK_LIBRARIES(Prob_Project ${ITK_LIBRARIES} vtkHybrid)

find_package(GenerateCLP REQUIRED)
include(${GenerateCLP_USE_FILE})

set(main_SOURCE main.cxx )
GENERATECLP(main_SOURCE main.xml)
link_directories(${vtkITK_LIB_DIR})
add_executable(main ${main_SOURCE} 
	Prob_Project
	${TXX_FILES}
	${CXX_FILES}
	${H_FILES}
)
target_link_libraries(main vtkImaging vtkGraphics vtkIO ${ITK_LIBRARIES} vtkHybrid)
include_directories(${vtkITK_SOURCE_DIR} ${vtkITK_BINARY_DIR})
