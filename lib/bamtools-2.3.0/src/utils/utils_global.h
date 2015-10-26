// ***************************************************************************
// utils_global.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides macros for exporting & importing BamTools-utils library symbols
// ***************************************************************************

#ifndef UTILS_GLOBAL_H
#define UTILS_GLOBAL_H

#include "shared/bamtools_global.h"

#ifdef BAMTOOLS_UTILS_LIBRARY
#  define UTILS_EXPORT BAMTOOLS_LIBRARY_EXPORT
#else
#  define UTILS_EXPORT BAMTOOLS_LIBRARY_IMPORT
#endif

#endif // UTILS_GLOBAL_H
