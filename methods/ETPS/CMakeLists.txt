PROJECT(spixel)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(OpenCV REQUIRED)
FIND_PACKAGE(OpenMP) 
FIND_PACKAGE(PNG REQUIRED)

FIND_PATH(png++_INCLUDE_DIR
	NAMES  png++/color.hpp  png++/config.hpp png++/consumer.hpp png++/convert_color_space.hpp png++/end_info.hpp  png++/error.hpp  png++/ga_pixel.hpp png++/generator.hpp png++/gray_pixel.hpp png++/image.hpp png++/image_info.hpp png++/index_pixel.hpp png++/info.hpp png++/info_base.hpp  png++/io_base.hpp  png++/packed_pixel.hpp png++/palette.hpp png++/pixel_buffer.hpp  png++/pixel_traits.hpp png++/png.hpp png++/reader.hpp png++/require_color_space.hpp png++/rgb_pixel.hpp png++/rgba_pixel.hpp png++/streaming_base.hpp  png++/tRNS.hpp png++/types.hpp png++/writer.hpp)
  
IF(OPENMP_FOUND)
    MESSAGE("OpenMP Enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(OPENMP_FOUND)

IF(CMAKE_COMPILER_IS_GNUCXX)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -msse4.2") 
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

SET(HEADER_FILES src/functions.h src/segengine.h src/stdafx.h src/structures.h src/targetver.h src/tsdeque.h src/utils.h src/contrib/SGMStereo.h)
SET(SOURCE_FILES src/stdafx.cpp src/utils.cpp src/functions.cpp src/structures.cpp src/segengine.cpp src/spixel.cpp src/contrib/SGMStereo.cpp)

INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS} ${PNG_INCLUDE_DIR} ${png++_INCLUDE_DIR})
ADD_EXECUTABLE(spixel ${SOURCE_FILES} ${HEADER_FILES})
TARGET_LINK_LIBRARIES(spixel ${OpenCV_LIBS})

INSTALL(TARGETS spixel DESTINATION bin)
