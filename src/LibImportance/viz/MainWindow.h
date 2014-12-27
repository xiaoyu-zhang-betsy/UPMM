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

#pragma warning(disable:4996)
#define __WXMSW__
#define WXUSINGDLL
#include <wx/wx.h>
#include <wx/spinctrl.h>

#include "GlWindow.h"
#include "State.h"
#include "RecordViz.h"

namespace Importance {

    enum Ids {
        ID_WHARRGARBL = wxID_HIGHEST + 1,
        ID_WHARRGARBL2,
        ID_DISPLAY_PHOTONS,
        ID_DISPLAY_RECORDS,
        ID_NTH_PHOTON,
        ID_NTH_RECORD,
        ID_DISPLAY_GEOMETRY,
        ID_GEOMETRY_BGCOLOR,
        ID_RECORDS_SAMECOLOR,
        ID_INTERPOLATE_ON_CLICK,
        ID_PHOTONS_IN_RECORD,
        ID_GEOMETRY_IN_RECORD,
        ID_RESET_VIEW,
        ID_MAPPING_TYPE,
        ID_SHIRLEY_MAPPING,
        ID_SHOW_RADII,
        ID_PHOTON_LENGTH,
        ID_RELATIVE_COLOR_MAPPING,
        ID_PHOTONS_FALSE_COLOR,
        ID_DISPLAY_ACTIVE_PHOTONS,
        ID_SHOW_CAMERA_RATIO,
        ID_SHOW_INTERPOLATION_METRIC,
        ID_METRIC_NONE,
        ID_METRIC_KL,
        ID_METRIC_IRRADIANCE,
        ID_METRIC_D_IRRADIANCE,
        ID_METRIC_TEST,
		ID_METRIC_GRADIENT,
        ID_SHOW_MARKS,
        ID_SCREENSHOT,
        ID_SHOWENVIRONMENTONLY,
		ID_GAMMA,
    };

    template<class TDistributionModel>
    class MainWindow : public wxFrame {
#ifdef LIBIMP_GATHER_PARTICLES
        Jensen * jensenDist;        
#endif

        CacheViz<TDistributionModel>* glCanvas;
        RecordViz<TDistributionModel>* recordCanvas;
        wxStaticText* recordInfo;
        wxBoxSizer* sideSizer;
        SamplersFacade<TDistributionModel> * m_samplers;
        wxMenu * cacheMenu;

        void OnPaint(wxPaintEvent&) {
        }

        void OnClose(wxCloseEvent&) {
            Show(false);

            if ( m_samplers != NULL ) {
                delete m_samplers;
            }
        }
        void OnSpin(wxSpinEvent& evt) {
            switch(evt.GetId()) {
            case ID_NTH_RECORD:
                state.everyNthRecord = evt.GetValue();
                break;
            case ID_NTH_PHOTON:
                state.everyNthPhoton = evt.GetValue();
                break;
            case ID_PHOTON_LENGTH:
                state.photonLength = evt.GetValue()/100.f;
                break;
			}
            Refresh();
        }
		void OnSpinDouble(wxSpinDoubleEvent& evt) {
			switch (evt.GetId()) {
			case ID_GAMMA:
				state.recordViz.gamma = evt.GetValue();
				break;
			}
			Refresh();
		}
        void OnCommand(wxCommandEvent& evt) {
            const bool checked = evt.IsChecked();
            switch(evt.GetId()) {
            case ID_DISPLAY_PHOTONS:
                state.showPhotons = checked;
                break;
            case ID_DISPLAY_ACTIVE_PHOTONS:
                state.showActivePhotons = checked;
                break;
            case ID_DISPLAY_RECORDS:
                state.showRecords = checked;
                break;
            case ID_DISPLAY_GEOMETRY:
                state.showGeometry = checked;
                break;
            case ID_GEOMETRY_BGCOLOR:
                state.isGeometryBgColor = checked;
                break;
            case ID_RECORDS_SAMECOLOR:
                state.isRecordsSameColor = checked;
                break;
            case ID_INTERPOLATE_ON_CLICK:
                state.interpolateOnClick = checked;
                break;
            case ID_GEOMETRY_IN_RECORD: 
                if(!checked || state.triangles.size() < 3000 || wxMessageDialog(this, "Really show geometry? You have a shitton of polygons in the scene, so you should go make a coffee if you click 'yes'.", "Confirmation", wxOK|wxCANCEL|wxCANCEL_DEFAULT|wxCENTRE|wxICON_QUESTION|wxSTAY_ON_TOP).ShowModal() == wxID_OK) {
                    state.recordViz.showGeometry = checked;
                }
                break;
            case ID_PHOTONS_IN_RECORD:
                state.recordViz.showPhotons = checked;
                break;
            case ID_RESET_VIEW:
                this->glCanvas->resetView();
                break;
            case ID_MAPPING_TYPE:
                state.recordViz.mapping = (evt.GetSelection() == 0) ? MAPPING_HEMISPHERE : MAPPING_SHIRLEY;
                break;
            case ID_SHIRLEY_MAPPING:
                state.recordViz.mapping = checked ? MAPPING_SHIRLEY : MAPPING_HEMISPHERE;
                break;
            case ID_SHOW_RADII:
                state.showRadii = checked;
                break;
            case ID_RELATIVE_COLOR_MAPPING:
                state.recordViz.relativeColorMapping = checked;
                break;
            case ID_PHOTONS_FALSE_COLOR:
                state.isPhotonsFalseColor = checked;
                break;
            case ID_SHOW_CAMERA_RATIO:
                state.showCameraRatio = checked;
                break;
            case ID_SHOW_MARKS:
                state.isShowMarks = checked;
                break;
            case ID_METRIC_NONE:
                state.displayedMetric = METRIC_NONE;
                break;
            case ID_METRIC_TEST:
                state.displayedMetric = METRIC_TEST;
                break;
            case ID_METRIC_KL:
                state.displayedMetric = METRIC_KL_DIVERGENCE;
                break;
            case ID_METRIC_IRRADIANCE:
                state.displayedMetric = METRIC_IRRADIANCE;
                break;
            case ID_METRIC_D_IRRADIANCE:
                state.displayedMetric = METRIC_D_IRRADIANCE;
                break;
			case ID_METRIC_GRADIENT:
                state.displayedMetric = METRIC_GRADIENT;
                break;
            case ID_SCREENSHOT: {
                bool res = recordCanvas->getBitmap().SaveFile( "recordViz.bmp", wxBITMAP_TYPE_BMP );
                                }
                break;
            case ID_SHOWENVIRONMENTONLY:
                state.showEnvironmentOnly = !state.showEnvironmentOnly;
                break;
            }
            Refresh();
        }

        bool* closedFlag;
    public:

        State state;

        MainWindow( const std::string & caption ) : wxFrame(NULL, wxID_ANY, caption, wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE|wxFULL_REPAINT_ON_RESIZE ) {

#ifdef LIBIMP_GATHER_PARTICLES
            jensenDist = NULL;
#endif

            wxMenuBar* menu = new wxMenuBar;
            cacheMenu = new wxMenu;
            wxMenu* recordMenu = new wxMenu;
            cacheMenu->Append(ID_RESET_VIEW, "Reset view");
            cacheMenu->AppendSeparator();
            cacheMenu->AppendCheckItem(ID_DISPLAY_PHOTONS, wxT("Display &Photons"));
            cacheMenu->Check(ID_DISPLAY_PHOTONS, state.showPhotons);
            cacheMenu->AppendCheckItem(ID_DISPLAY_RECORDS, wxT("Display &Records"));
            cacheMenu->Check(ID_DISPLAY_RECORDS, state.showRecords);
            cacheMenu->AppendCheckItem(ID_DISPLAY_GEOMETRY, wxT("Display &Geometry"));
            cacheMenu->Check(ID_DISPLAY_GEOMETRY, state.showGeometry);
            cacheMenu->AppendCheckItem(ID_GEOMETRY_BGCOLOR, "Geometry in Bg color");
            cacheMenu->Check(ID_GEOMETRY_BGCOLOR, state.isGeometryBgColor);
            cacheMenu->AppendCheckItem(ID_RECORDS_SAMECOLOR, "Records in same color");
            cacheMenu->Check(ID_RECORDS_SAMECOLOR, state.isRecordsSameColor);
            cacheMenu->AppendCheckItem(ID_DISPLAY_ACTIVE_PHOTONS, "Display KNN Photons");
            cacheMenu->Check(ID_DISPLAY_ACTIVE_PHOTONS, state.showActivePhotons);
            cacheMenu->AppendCheckItem(ID_SHOW_RADII, "Display all radii");
            cacheMenu->Check(ID_SHOW_RADII, state.showRadii);
            cacheMenu->AppendCheckItem(ID_PHOTONS_FALSE_COLOR, "Paint photons in false colors" );
            cacheMenu->Check(ID_PHOTONS_FALSE_COLOR, state.isPhotonsFalseColor);
            cacheMenu->AppendCheckItem(ID_SHOW_CAMERA_RATIO, "Show camera pixel->world" );
            cacheMenu->Check(ID_SHOW_CAMERA_RATIO, state.showCameraRatio);
            cacheMenu->AppendCheckItem(ID_SHOW_MARKS, "Show custom marks");
            cacheMenu->Check(ID_SHOW_MARKS, state.isShowMarks);
            cacheMenu->AppendCheckItem(ID_SHOWENVIRONMENTONLY, "Show environment only");
            cacheMenu->Check(ID_SHOWENVIRONMENTONLY, state.showEnvironmentOnly);            


   //         wxMenu* submenu = new wxMenu;
   //         submenu->AppendRadioItem(ID_METRIC_NONE, "None");
   //         submenu->AppendRadioItem(ID_METRIC_KL, "KL divergence");
   //         submenu->AppendRadioItem(ID_METRIC_IRRADIANCE, "Irradiance");
   //         submenu->AppendRadioItem(ID_METRIC_D_IRRADIANCE, "dIrradiance");
   //         submenu->AppendRadioItem(ID_METRIC_TEST, "Test");
			//submenu->AppendRadioItem(ID_METRIC_GRADIENT, "Gradient");
   //         cacheMenu->AppendSubMenu(submenu, "Displayed metric");            

            recordMenu->AppendCheckItem(ID_INTERPOLATE_ON_CLICK, wxT("&Interpolate on click"));
            recordMenu->Check(ID_INTERPOLATE_ON_CLICK, state.interpolateOnClick);
            recordMenu->AppendCheckItem(ID_SHIRLEY_MAPPING, "Use Shirley mapping");
            recordMenu->Check(ID_SHIRLEY_MAPPING, state.recordViz.mapping == MAPPING_SHIRLEY);
            recordMenu->AppendCheckItem(ID_RELATIVE_COLOR_MAPPING, "Relative false color mapping");
            recordMenu->Check(ID_RELATIVE_COLOR_MAPPING, state.recordViz.relativeColorMapping);
            recordMenu->AppendSeparator();
            recordMenu->AppendCheckItem(ID_PHOTONS_IN_RECORD, "Show photons");
            recordMenu->Check(ID_PHOTONS_IN_RECORD, state.recordViz.showPhotons);
            //recordMenu->AppendCheckItem(ID_GEOMETRY_IN_RECORD, "Show geometry");
            //recordMenu->Check(ID_GEOMETRY_IN_RECORD, state.recordViz.showGeometry);

            
            menu->Append(cacheMenu, wxT("&Cache view"));
            menu->Append(recordMenu, wxT("&Record view"));
            SetMenuBar(menu);

            wxSizer* control = new wxGridSizer(2);
            control->Add(new wxStaticText(this, wxID_ANY, "Every Nth photon"));
            control->Add(new wxSpinCtrl(this, ID_NTH_PHOTON, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxALIGN_RIGHT, 1, 999, state.everyNthPhoton, "Every nth photon"));

            control->Add(new wxStaticText(this, wxID_ANY, "Every Nth Record"));
            control->Add(new wxSpinCtrl(this, ID_NTH_RECORD, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxALIGN_RIGHT, 1, 999, state.everyNthRecord, "Every nth record"));

            control->Add(new wxStaticText(this, wxID_ANY, "Photon length mult."));
            control->Add(new wxSpinCtrl(this, ID_PHOTON_LENGTH, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS|wxALIGN_RIGHT, 0, 100, 100*state.photonLength, "PhotonLengthMult"));

			control->Add(new wxStaticText(this, wxID_ANY, "PDF display gamma"));
			control->Add(new wxSpinCtrlDouble(this, ID_GAMMA, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS | wxALIGN_RIGHT, 0.01f, 99.f, state.recordViz.gamma, 0.1f, "PDF display gamma"));


            control->Add( new wxButton( this, ID_SCREENSHOT, "Screeshot" ) );



            this->glCanvas = new CacheViz<TDistributionModel>(this);
            this->glCanvas->createGlContext();
            this->recordCanvas = new RecordViz<TDistributionModel>(this, &state);
            sideSizer = new wxBoxSizer(wxVERTICAL);
            sideSizer->Add(control);
            this->recordInfo = new wxStaticText(this, wxID_ANY, "\n\n\n\n\n\n\n\n\n\n\n\n\n\n", wxDefaultPosition, wxSize(300, wxDefaultCoord));
            wxFont smallFont(8, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
            this->recordInfo->SetFont(smallFont);
            sideSizer->Add(recordInfo);

            sideSizer->Add(recordCanvas);
            SetSizerAndFit(sideSizer);

            Connect(wxEVT_CLOSE_WINDOW,             wxCloseEventHandler(MainWindow::OnClose));
            Connect(wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(MainWindow::OnCommand));
            Connect(wxEVT_COMMAND_SPINCTRL_UPDATED, wxSpinEventHandler(MainWindow::OnSpin));
			Connect(wxEVT_COMMAND_SPINCTRLDOUBLE_UPDATED, wxSpinDoubleEventHandler(MainWindow::OnSpinDouble));
            Connect(wxEVT_PAINT,                    wxPaintEventHandler(MainWindow::OnPaint));
            Connect(wxEVT_CLOSE_WINDOW,             wxCommandEventHandler(MainWindow::OnQuit));
            Connect(wxEVT_COMMAND_MENU_SELECTED,    wxCommandEventHandler(MainWindow::OnCommand));
            Connect(wxEVT_SIZE,                     wxSizeEventHandler(MainWindow::OnResize));
            Connect(wxEVT_COMMAND_COMBOBOX_SELECTED, wxCommandEventHandler(MainWindow::OnCommand));
            Connect(wxEVT_BUTTON,                   wxCommandEventHandler(MainWindow::OnCommand));
            this->SetSize(1024, 768);
            recursiveSetup(this);
        }

        void setImportanceConfig( const Config & cfg ) {
            state.config = cfg;
        }

        void OnResize(wxSizeEvent& evt) {
            wxSize size = evt.GetSize();
            const int glSize = size.GetX()*2/3;
            const int recordSize = size.GetX() - glSize;
            this->glCanvas->SetSize(glSize, size.GetY());
            this->glCanvas->SetPosition(wxPoint(recordSize, 0));
            this->recordCanvas->SetSize(recordSize, recordSize);
            Refresh(true);
        }

        void recursiveSetup( wxWindow* window) {
            IMPORTANCE_ASSERT(window);

            typedef CacheViz<TDistributionModel> CacheType;
            // connect key events
            window->Connect(wxID_ANY, wxEVT_KEY_UP, wxKeyEventHandler(CacheType::OnKey), (wxObject*) NULL, this->glCanvas);
            window->Connect(wxEVT_MOUSEWHEEL, wxMouseEventHandler(CacheType::OnScroll), NULL, glCanvas);

            wxWindowListNode* child = window->GetChildren().GetFirst();
            while(child) {
                this->recursiveSetup(child->GetData());
                child = child->GetNext();
            }
        }

        void OnQuit(wxCommandEvent&) {
            if(this->closedFlag) {
                *this->closedFlag = true;
            }
        }

        ~MainWindow() {
            if(this->closedFlag) {
                *this->closedFlag = true;
            }
        }

        struct Pos {
            Pos( int _x, int _y ) : x(_x), y(_y) {}
            int x, y;
        };

        void plotJensens( const std::vector<Pos> & points, SimpleSampler<TDistributionModel> * sampler ) {
            if ( sampler->m_gatheredParticles.size() == points.size() ) {
                for ( size_t i = 0; i < points.size(); ++i ) {
                    int res[2] = {32, 64};
                    for ( size_t j = 0; j < 2; ++j ) {
                        if ( jensenDist != NULL ) {
                            jensenDist->release();
                            delete jensenDist;
                        }
                        jensenDist = new Importance::Jensen();
                        jensenDist->init(  res[j], res[j], sampler->m_gatherPoints[i].normal  );
                        for ( auto p = sampler->m_gatheredParticles[i].begin(); p < sampler->m_gatheredParticles[i].end(); ++p ) {
                            jensenDist->add( *p, Vector3(), -1.f );
                        }


                        // Write it to a bitmap and to the file
                        bool relativeColorMapping = false;                
                        wxSize size;
                        size.x = 619; size.y = 619;
                        wxBitmap drawbuffer(size, 24);
                        wxNativePixelData data(drawbuffer);
                        wxNativePixelData::Iterator p(data);

                        ImBitmap<float> pdfs(size.x, size.y);
                        float maxPdf = 1e-6f;
                        for(int x = 0; x < size.x; ++x) {
                            for(int y = 0; y < size.y; ++y) {
                                const Vector3 dir = state.coordsToDir(x, y, size.x, size.y, sampler->m_gatherPoints[ i ].normal );
                                pdfs(x, y) = jensenDist->pdf(dir);
                                maxPdf = std::max(maxPdf, pdfs(x, y));
                            }
                        }
                        const float maxMappedIntensity = recordCanvas->mapIntensity(relativeColorMapping ? maxPdf : 255.f);

                        for(int x = 0; x < size.x; ++x) {
                            for(int y = 0; y < size.y; ++y) {
                                p.MoveTo(data, x, y);
                                unsigned char r,g,b;
                                Float w = std::min(1.f, recordCanvas->mapIntensity(pdfs(x, y))/maxMappedIntensity);
                                Vector3 color;
                                VizTools::getColor2(w, &color);
                                color *= 255.f;
                                r = color.x; g = color.y; b = color.z;                        
                                p.Red() = r;
                                p.Green() = g;
                                p.Blue() = b;
                            }
                        }

                        std::ostringstream ostr;
                        ostr << "dist_" << points[ i ].x << "_" << points[ i ].y << "_ref_res" << res[ j ] << ".png";
                        drawbuffer.SaveFile(ostr.str().c_str(), wxBITMAP_TYPE_PNG);
                    }
                }
            }        
        }


        void show(SimpleSampler<TDistributionModel>* sampler, 
                  IEnviroSampler * enviroSampler, TriangleIterator* triangles, 
                  bool* closedFlag, const InitCamera* camera, const std::vector<Hit> * markers ) {
            this->state.marks = markers;
            this->closedFlag = closedFlag;
            if(closedFlag) {
                *closedFlag = false;
            }            
            this->m_samplers = new SamplersFacade<TDistributionModel>( enviroSampler, sampler );
            this->recordCanvas->setSamplers( m_samplers );            
            this->glCanvas->setSamplers( m_samplers );
            this->glCanvas->setup(&state, triangles, camera);
            this->glCanvas->computeSceneStatistics();
            this->cacheMenu->Enable(ID_SHOWENVIRONMENTONLY, m_samplers->getEnviroSampler() != NULL);
            Iconize(false);
			this->Maximize(true);
            Show(true);

#ifdef LIBIMP_GATHER_PARTICLES  // pro vytvareni figure distribuci do paperu
			this->state.recordViz.showPhotons = false;
			this->state.recordViz.relativeColorMapping = false;
			wxMouseEvent evt;

            // generating teaser image
            std::vector<Pos> points;
            //table
            //points.push_back( Pos( 886, 2088 ) );
            points.push_back( Pos( 1071, 878 ) ); //table right down
            points.push_back( Pos( 1015, 929 ) ); //table inner corner
            points.push_back( Pos( 1062, 958 ) );
            points.push_back( Pos( 985, 916 ) ); //under the table
            points.push_back( Pos( 974, 891 ) );
            points.push_back( Pos( 981, 892 ) );



            //books
            points.push_back( Pos( 470, 790 ) ); //shelf front
            points.push_back( Pos( 463, 821 ) ); //shelf front 2
            points.push_back( Pos( 371, 813 ) ); //book
            points.push_back( Pos( 299, 765 ) ); //shelf inside (noisy place in insets)
            points.push_back( Pos( 374, 746 ) );

            plotJensens( points, sampler );

            std::vector<Hit> gatherHitPoints;
            for ( auto p = points.cbegin(); p < points.cend(); ++p ) {
                evt.SetX(p->x);
                evt.SetY(p->y);
                this->glCanvas->leftClickEvent(evt);
                if ( state.lastClick.position.isReal() ) {
                    gatherHitPoints.push_back( state.lastClick );
                }
                this->recordCanvas->Update();
                std::ostringstream ostr;
                ostr << "dist_" << p->x << "_" << p->y << ".png";
                this->recordCanvas->getBitmap().SaveFile(ostr.str().c_str(), wxBITMAP_TYPE_PNG);
            }
            
            sampler->saveGatherPoints( gatherHitPoints );

            //Ondra's for generating cbox figures
			//evt.SetX(612);
			//evt.SetY(331);
			//this->glCanvas->leftClickEvent(evt);
			//this->recordCanvas->Update();
			//this->recordCanvas->getBitmap().SaveFile("dist1.png", wxBITMAP_TYPE_PNG);

			//evt.SetX(374);
			//evt.SetY(238);
			//this->glCanvas->leftClickEvent(evt);
			//this->recordCanvas->Update();
			//this->recordCanvas->getBitmap().SaveFile("dist2.png", wxBITMAP_TYPE_PNG);
#endif
        }
    };
}
