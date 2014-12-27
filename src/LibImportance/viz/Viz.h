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
#include "../LibImportance.h"
#pragma warning(disable:4996)
#define __WXMSW__
#define WXUSINGDLL
#include <wx/wx.h>
#include "MainWindow.h"

#include <windows.h>
#include <stdio.h>
#include <fcntl.h>
#include <io.h>
#include <iostream>
#include <fstream>

class UseConsole {
    public:
        UseConsole() {
            m_good = !!AttachConsole(ATTACH_PARENT_PROCESS);
            if (m_good)
                RedirectIOToConsole();
        }

        ~UseConsole() {
            if (m_good)
                FreeConsole();
        }
    private:
        // The following function is taken nearly verbatim from
        // http://www.halcyon.com/~ast/dload/guicon.htm
        void RedirectIOToConsole() {
            int hConHandle;
            long lStdHandle;
            FILE *fp;

            // redirect unbuffered STDOUT to the console
            lStdHandle = (long)GetStdHandle(STD_OUTPUT_HANDLE);
            hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
            fp = _fdopen( hConHandle, "w" );
            *stdout = *fp;
            setvbuf( stdout, NULL, _IONBF, 0 );

            // redirect unbuffered STDIN to the console
            lStdHandle = (long)GetStdHandle(STD_INPUT_HANDLE);
            hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
            fp = _fdopen( hConHandle, "r" );
            *stdin = *fp;
            setvbuf( stdin, NULL, _IONBF, 0 );

            // redirect unbuffered STDERR to the console
            lStdHandle = (long)GetStdHandle(STD_ERROR_HANDLE);
            hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
            fp = _fdopen( hConHandle, "w" );
            *stderr = *fp;
            setvbuf( stderr, NULL, _IONBF, 0 );

            // make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog
            // point to console as well
            std::ios::sync_with_stdio();
        }

    private:
        bool    m_good;
};

namespace Importance {
    class VizApp : public wxApp {
    private:
        UseConsole useConsole;

        virtual bool OnInit() {
            return true;
        }
    };


    template<class TDistributionModel>
    void showVizWindow(
        const std::string & caption, 
        SimpleSampler<TDistributionModel>* sampler, 
        IEnviroSampler * enviroSampler,
        const std::vector<Hit> * markers,
        const bool initWxWidgets = true,
        TriangleIterator* geometry = NULL, 
        bool* closedFlag = NULL, 
        const InitCamera* camera = NULL, 
        const bool runEventLoop = false, 
        const Config * impCfg = NULL ) {

            wxApp::SetInstance(new VizApp());
			wxInitializer * wxInit = new wxInitializer();
            MainWindow<TDistributionModel>* window = new MainWindow<TDistributionModel>( caption );
            wxTheApp->SetTopWindow(window);
            bool closed;
            if(impCfg) {
                window->setImportanceConfig( *impCfg );
            }
            window->show(sampler, enviroSampler, geometry, (closedFlag == NULL || runEventLoop) ? &closed : closedFlag, camera, markers );

            if(runEventLoop) {
                MSG msg;
                while(!closed) {
                    int ret = GetMessage(&msg, NULL, 0, 0); ret = ret;
                    TranslateMessage(&msg);
                    DispatchMessage(&msg);
                }
            }
			delete wxInit;
    }

    inline void showVizWindow(
        const std::string & caption, 
        Sampler * sampler, 
        IEnviroSampler * enviroSampler,
        const Config & cfg,
        const std::vector<Hit> * markers,
        const bool initWxWidgets = true,
        TriangleIterator* geometry = NULL, 
        bool* closedFlag = NULL, 
        const InitCamera* camera = NULL, 
        const bool runEventLoop = false ) {

        if ( sampler == NULL ) {
            /* This is the case when we need to visualize only environment sampler */
            showVizWindow<GaussianMixtureSSEHemisphere>( caption, NULL, enviroSampler, markers, initWxWidgets, geometry, closedFlag, camera, runEventLoop, &cfg );
            return;
        }

            // our "brilliant" software engineering stuff
        switch ( cfg.distType ) {
            case EPharr:
                {
                    PharrSimpleSampler * tmp = dynamic_cast<PharrSimpleSampler*>( sampler );
                    if ( tmp != NULL ) {
                        showVizWindow( caption, tmp, enviroSampler,  markers, initWxWidgets, geometry, closedFlag, camera, runEventLoop, &cfg );
                    }
                }
                break;
            case EHey:
                {
                    HeySimpleSampler * tmp = dynamic_cast<HeySimpleSampler*>( sampler );
                    if ( tmp != NULL ) {
                        showVizWindow( caption, tmp, enviroSampler,  markers, initWxWidgets, geometry, closedFlag, camera, runEventLoop, &cfg );
                    }                
                }
                break;
            case EGaussianMixture: 
                {
                    GaussianMixtureSimpleSampler * tmp = dynamic_cast<GaussianMixtureSimpleSampler*>( sampler );
                    if ( tmp != NULL ) {
                        showVizWindow( caption, tmp, enviroSampler,  markers, initWxWidgets, geometry, closedFlag, camera, runEventLoop, &cfg );
                    }                
                }
                break;
            case EJensen:
                {
                    JensenSimpleSampler * tmp = dynamic_cast<JensenSimpleSampler*>( sampler );
                    if ( tmp != NULL ) {
                        showVizWindow( caption, tmp, enviroSampler,  markers, initWxWidgets, geometry, closedFlag, camera, runEventLoop, &cfg );
                    }                
                }
                break;
        }
    }


}
