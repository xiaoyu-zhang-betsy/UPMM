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


#include <windows.h>
#include "Timer.h"
#include <iostream>
#include <sstream>

namespace Importance { 

    /// Create a new timer and reset it
    Timer::Timer() { reset(); }

    /// Reset the timer
    void Timer::reset() {
    #ifdef WIN32
        QueryPerformanceFrequency((PLARGE_INTEGER) &m_frequency);
        QueryPerformanceCounter((PLARGE_INTEGER) &m_start);
    #else
        gettimeofday(&m_start, NULL);
    #endif
    }

    /// Return the milliseconds which have passed since the last reset
    unsigned int Timer::getMilliseconds() const {
    #ifdef WIN32
        LARGE_INTEGER current;
        QueryPerformanceCounter(&current);

        return (int) ((current.QuadPart - m_start) * 1000 / m_frequency);
    #else
        struct timeval current;

        gettimeofday(&current, NULL);

        return (current.tv_sec - m_start.tv_sec) * 1000 +
            (current.tv_usec - m_start.tv_usec) / 1000;
    #endif
    }

    /// Return the microseconds which have passed since the last reset
    unsigned int Timer::getMicroseconds() const {
    #ifdef WIN32
        LARGE_INTEGER current;
        QueryPerformanceCounter(&current);

        return (int) ((current.QuadPart - m_start) * 1000000 / m_frequency);
    #else
        struct timeval current;

        gettimeofday(&current, NULL);

        return (current.tv_sec - m_start.tv_sec) * 1000000 +
            (current.tv_usec - m_start.tv_usec);
    #endif
    }

    /// Return a string representation
    std::string Timer::toString() const {
        std::ostringstream oss;
        oss << "Timer[ms=" << getMilliseconds() << "]";
        return oss.str();
    }
}