/*
    This file is part of LibImportance library that provides a technique for guiding
    transport paths towards the important places in the scene. This is a direct implementation
    of the method described in the paper "On-line Learning of Parametric Mixture 
    Models for Light Transport Simulation", ACM Trans. Graph. (SIGGRAPH 2014) 33, 4 (2014).
   
    Copyright (c) 2014 by Jiri Vorba, Ondrej Karlik, Martin Sik.

    LibImportance library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    LibImportance library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once


//#define SINGLE_PRECISION

// CONFIGURATION ////////////////////////////////////////////////////////////////////////////////////////////////////

//#define FAST_FLOAT
//#define CUSTOM_RANDOM

//#define IMPORTANCE_DEBUG
#ifndef DOUBLE_PRECISION
    #define IMPORTANCE_SINGLE_PRECISION
#endif

namespace Importance {
    //const int MAX_KNN_PARTICLES = 20000;
    const int ENVIRO_DIR_SEARCH_COUNT = 100;
    const int MAX_KNN_PARTICLES = 1000;
    const int PARTICLE_TREE_LEAF_SIZE = 20;
    const int CACHE_TREE_LEAF_SIZE = 5;

    const int MAX_FITTED_LOBES = 64;
    //const int MAX_FITTED_LOBES = 2000;
    const int MAX_CACHE_KNN = 3;
    //const int MAX_BRDF_LOBES = 2;

    const int ENVIRO_DIR_MAX_LOBES = 100;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <limits>

#define IMPORTANCE_INLINE __forceinline

//#define IMPORTANCE_ASSERT(x) _ASSERTE(x)
//#define IMPORTANCE_ASSERT(cond) do { \
//    if ( !(cond) ) { \
//        if (1 /*IsDebuggerPresent()*/) { \
//            __debugbreak(); \
//        } \
//        else { \
//            throw std::runtime_error(#cond); \
//        } \
//    } \
//} while(0)

//#define LIBIMP_STATS
//#define LIBIMP_GATHER_PARTICLES
//#define LIBIMP_LOG

#if defined(LIBIMP_DEBUG) || defined( _DEBUG )
    #define IMPORTANCE_ASSERT(cond) if ( !(cond) ) { __debugbreak(); }
    #undef LIBIMP_DEBUG
    #define LIBIMP_DEBUG
    //#pragma message( "Importance library debugging is ON." )
#else
    //#pragma message( "Importance library debugging is OFF." )
    #define IMPORTANCE_ASSERT(cond) /*if ( !(cond) ) { __debugbreak(); }*///{ throw std::runtime_error(#cond); }
#endif

namespace Importance {

    class Vector3d;
    class Vector3f;    

    typedef _int64 StatsInt64;

#ifdef IMPORTANCE_SINGLE_PRECISION
    typedef float Float;
    typedef Vector3f Vector3;
#else
    typedef double Float;
    typedef Vector3d Vector3;    
#endif

    const float FLOAT_NAN = std::numeric_limits<float>::signaling_NaN();
    const float FLOAT_INFINITY = std::numeric_limits<float>::infinity();
    const double DOUBLE_NAN = std::numeric_limits<double>::signaling_NaN();
    const double DOUBLE_INFINITY = std::numeric_limits<double>::infinity();
    const Float NAN = std::numeric_limits<Float>::signaling_NaN();
    const Float INFINITY = std::numeric_limits<Float>::infinity();
    const Float ONE = 1.f;
    const Float ZERO = 0.f;
    const Float NORMAL_TO_NORMAL_LIMIT      = 1e-2f;
    const Float INCID_DIR_TO_NORMAL_LIMIT   = 1e-3f;

    const Float PI = Float(3.141592653589793238462643383279);    

    struct NullObject { };
}

#ifdef LIBIMP_STATS
#define IMPORTANCE_INC_STATS(type) ((type)++)
#define IMPORTANCE_INC_STATS_MORE(type, value) ((type) += (value))
#else
#define IMPORTANCE_INC_STATS(x)
#define IMPORTANCE_INC_STATS_MORE(type, value)
#endif
