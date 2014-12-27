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


#ifndef ___READER_WRITER_LOCK___
#define ___READER_WRITER_LOCK___

#include <windows.h>

#ifdef min
    #undef min
#endif
#ifdef max
    #undef max
#endif

#include "../shared/Config.h"

namespace Importance {
    /// \brief Abstraction of a reader-writer lock. It is used for synchronization of data 
    ///        structures which can be read by multiple readers at once, but only one thread 
    ///        can modify it at any given time.
    class ReaderWriterLock {
    protected:

        /// \brief win API data structure for the lock
        SRWLOCK lock;

    public:

        ReaderWriterLock() {
            memset(&lock, 0, sizeof(lock));
            InitializeSRWLock(&lock);
        }

        /// \brief Locks for reading. Will not block other threads from reading.
        IMPORTANCE_INLINE void lockRead() {
            AcquireSRWLockShared(&lock);
        }

        /// \brief Unlocks reading for current thread
        IMPORTANCE_INLINE void unlockRead() {
            ReleaseSRWLockShared(&lock);
        }

        /// \brief Locks writing and reading
        IMPORTANCE_INLINE void lockWrite() {
            AcquireSRWLockExclusive(&lock);
        }

        /// \brief Unlocks writing and reading
        IMPORTANCE_INLINE void unlockWrite() {
            ReleaseSRWLockExclusive(&lock);
        }
    };
}

#endif