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

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdarg>
#include "Config.h"

#pragma warning(push)
#pragma warning(disable:4127)

namespace Importance {

    enum ELogLevel {                
        ETrace,
        EDebug,
        EInfo,
        EWarn,
        EError
    };

    class SimpleLogger {
    public:
        static SimpleLogger * getInstance() {
            if ( s_logger == NULL ) {
                s_logger = new SimpleLogger();
            }

            return s_logger;
        }

        static void initStatic() {
            s_impLoggerClosed = false;
        }

        static void close() {           
            if ( s_logger != NULL ) {
                if ( SimpleLogger::s_impLoggerClosed ) {
                    s_logger->log( EWarn, "Importance logger was once closed during run of your program! It is possible that previously written data were lost. Ensure that \"Importance::close()\" method is called only once in your host code.", __LINE__, __FILE__ );
                }
                delete s_logger;
                s_logger = NULL;
            }
            SimpleLogger::s_impLoggerClosed = true;
        }

        void log( ELogLevel level, std::string message, int line, const char * file ) {
            if ( level < m_logLevel ) {
                return;
            }

            log( level, message.c_str(), line, file );
        }

        void log( ELogLevel level, const char * message, int line, const char * file, ... ) {
            if ( level < m_logLevel ) {
                return;
            }                        

            char tmp[512], *msg = tmp;
            va_list iterator;

            va_start(iterator, file);
            size_t size = _vscprintf(message, iterator) + 1;

            if (size >= sizeof(tmp))
                msg = new char[size];

            vsnprintf_s(msg, size, size-1, message, iterator);
            va_end(iterator);

            std::ostringstream ostr;
            ostr << file << " (" << line << "): " << msg << std::endl;

            std::cout << ostr.str();

            m_file.open( "importancelib.log", std::ofstream::app );
            if ( !m_file.fail() ) {
                m_file    << ostr.str();
                //m_file.flush();
			    m_file.close();
            }

            if (msg != tmp)
                delete[] msg;

            if ( level == EError ) {
                throw std::runtime_error( ostr.str().c_str() );
            }
        }

        void setLogLevel( ELogLevel level ) {
            m_logLevel = level;
        }

       ~SimpleLogger() {
           m_file.close();
       }
    private:
        SimpleLogger() {
            m_logLevel = EInfo;
            m_file.open( "importancelib.log" );
            m_file.close();
        }
    private:
        static SimpleLogger * s_logger;
        static bool s_impLoggerClosed;
        ELogLevel m_logLevel;
        std::ofstream m_file;
    };    
}

#ifdef LIBIMP_LOG
    #define ILog( level, msg, ... ) Importance::SimpleLogger::getInstance()->log( level, msg, __LINE__, __FILE__, ## __VA_ARGS__ );
#else
    #define ILog( level, msg, ... )
#endif


#pragma warning(pop)