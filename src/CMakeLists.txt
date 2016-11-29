cmake_minimum_required(VERSION 2.8)
project(Sniffles)

include_directories (../lib/bamtools-2.3.0/src)
include_directories(../lib/tclap-1.2.1/include)

configure_file( Version.h.in ${CMAKE_SOURCE_DIR}/src/Version.h )

add_executable(sniffles 	
                    tree/Breakpoint_Tree.cpp
                    Genotyper/Genotyper.cpp
					Alignment.cpp
					BamParser.cpp
					Sniffles.cpp
					Ignore_Regions.cpp
					tree/Intervall_bed.cpp
					sub/Detect_Breakpoints.cpp
					sub/Breakpoint.cpp
					tree/IntervallTree.cpp
					tree/IntervallList.cpp
					realign/SWCPU.cpp
					realign/Realign.cpp
					print/VCFPrinter.cpp
					print/BedpePrinter.cpp
					print/IPrinter.cpp
					tree/BinTree.cpp
					print/NGMPrinter.cpp
					plane-sweep/PlaneSweep_slim.cpp
					cluster/Cluster_SVs.cpp
					)										
					
#target_link_libraries(ngm-core pthread)	
TARGET_LINK_LIBRARIES(sniffles BamTools-static)
TARGET_LINK_LIBRARIES(sniffles zlibstatic)
	
add_executable(sniffles-debug 	
                    tree/Breakpoint_Tree.cpp
                    Genotyper/Genotyper.cpp
					Alignment.cpp
					BamParser.cpp
					Sniffles.cpp
					Ignore_Regions.cpp
					tree/Intervall_bed.cpp
					tree/IntervallList.cpp
					sub/Detect_Breakpoints.cpp
					sub/Breakpoint.cpp
					tree/IntervallTree.cpp
					realign/SWCPU.cpp
					realign/Realign.cpp
					print/VCFPrinter.cpp
					print/BedpePrinter.cpp
					print/IPrinter.cpp
					tree/BinTree.cpp
					print/NGMPrinter.cpp
					plane-sweep/PlaneSweep_slim.cpp
					cluster/Cluster_SVs.cpp
					)	

SET_TARGET_PROPERTIES(sniffles-debug PROPERTIES COMPILE_FLAGS "-g3 -O0")

#target_link_libraries(sniffles-debug pthread)
TARGET_LINK_LIBRARIES(sniffles-debug BamTools-static)
TARGET_LINK_LIBRARIES(sniffles-debug zlibstatic)

