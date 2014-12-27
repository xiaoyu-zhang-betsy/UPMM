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


#if !defined(__IMPORTANCE_TIMER_H)
#define __IMPORTANCE_TIMER_H

#include <string>
#ifndef WIN32
#include <sys/time.h>
#endif

namespace Importance {
    class Timer {
    public:
	    /// Create a new timer and reset it
        Timer();

	    /// Reset the timer
	    void reset();

	    /// Return the milliseconds which have passed since the last reset
        unsigned int getMilliseconds() const;

	    /// Return the microseconds which have passed since the last reset
        unsigned int getMicroseconds() const;

	    /// Return a string representation
        std::string toString() const;

	    /// Virtual destructor
        virtual ~Timer() {}
    private:
    #ifdef WIN32
	    __int64 m_start;
	    __int64 m_frequency;
    #else
	    struct timeval m_start;
    #endif
    };
}

#endif /* __IMPORTANCE_TIMER_H */