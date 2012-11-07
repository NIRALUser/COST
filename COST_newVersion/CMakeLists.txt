cmake_minimum_required(VERSION 2.6)
project(COST)

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

SET(CXX_FILES 
  main.cxx
  ODFReconstructor.cxx
  ODFStreamline.cxx
  ODFCost.cxx
)

SET(H_FILES
  ODFReconstructor.h
  ODFStreamline.h
  ODFCost.h
  ${CMAKE_CURRENT_BINARY_DIR}/define.h # Added by Adrien Kaiser
)

find_package(GenerateCLP REQUIRED)
include(${GenerateCLP_USE_FILE})

GENERATECLP(CXX_FILES Cost.xml)

# Added by Adrien Kaiser
SET (GK_Table_dir ${CMAKE_CURRENT_BINARY_DIR} CACHE FILEPATH "Table containing tensors computed by Gordon Kindlman" )
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/define.h.in ${CMAKE_CURRENT_BINARY_DIR}/define.h ) # configure and copy
include_directories(  ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}  )

add_executable( 
       COST
	${CXX_FILES}
	${H_FILES}
)

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/tensorGK.dat ${CMAKE_CURRENT_BINARY_DIR}/tensorGK.dat ) # Copy
TARGET_LINK_LIBRARIES(COST ${ITK_LIBRARIES} vtkHybrid)
