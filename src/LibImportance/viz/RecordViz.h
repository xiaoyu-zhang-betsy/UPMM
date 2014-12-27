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

#pragma once

#define __WXMSW__
#define WXUSINGDLL
#include <wx/wx.h>
#include <wx/rawbmp.h>
#include <wx/glcanvas.h>
#include <gl/GLU.h>

#include "State.h"

namespace Importance {


    template<class TDistributionModel>
    class RecordViz : public wxWindow {
    public:

        IMPORTANCE_INLINE void setSamplers( SamplersFacade<TDistributionModel> * samplers ) {
            m_samplers = samplers;
        }

        RecordViz(wxWindow* parent, State* state ) : wxWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize/*,  wxSize(320, 320)*/) {
            this->state = state;            
            SetBackgroundStyle(wxBG_STYLE_CUSTOM);
            Connect(wxEVT_PAINT, wxPaintEventHandler(RecordViz::OnPaint));
            Connect(wxEVT_LEFT_UP, wxMouseEventHandler(RecordViz::OnClick));
            Connect(wxEVT_MOTION, wxMouseEventHandler(RecordViz::OnMouseMotion));
        }

        wxBitmap getBitmap() {            
            return drawbuffer;
        }

        IMPORTANCE_INLINE float mapIntensity(const float intensity) const {
            //return log(log(intensity+1)+1)*0.5f;
            return pow(log(intensity + 1)*0.5f, 1.f/state->recordViz.gamma);
        }

        Particle getParticle( int index ) const {            
            return m_samplers->getParticle( index, this->state );
        } 

        void click(const int x, const int y) {
            state->recordViz.clickedX = x;
            state->recordViz.clickedY = y;
            const Vector3 dir = state->coordsToDir(state->recordViz.clickedX, state->recordViz.clickedY, GetSize().x, GetSize().y);
            state->recordViz.clickedPdf = state->pdf(dir);
            std::cout << "----Clicked PDF: " << state->recordViz.clickedPdf << std::endl;
            GetParent()->Refresh();
        }

    protected:
        State* state;
        SamplersFacade<TDistributionModel> * m_samplers;

        wxBitmap drawbuffer;
        int lastRecord;        

        IMPORTANCE_INLINE SimpleSampler<TDistributionModel> * getSampler() {
            return m_samplers->getSampler();
        }

        IMPORTANCE_INLINE const SimpleSampler<TDistributionModel> * getSampler() const { return m_samplers->getSampler(); }

        IMPORTANCE_INLINE IEnviroSampler * getEnviroSampler() { return m_samplers->getEnviroSampler(); }

        IMPORTANCE_INLINE const IEnviroSampler * getEnviroSampler() const { return m_samplers->getEnviroSampler(); }

        void OnClick(wxMouseEvent& evt) {
            click(evt.GetX(), evt.GetY());
            
        }
        void OnMouseMotion(wxMouseEvent& evt) {
            if(evt.LeftIsDown()) {
                click(evt.GetX(), evt.GetY());
            }
        }

        void OnPaint(wxPaintEvent&) {
            wxPaintDC dc(this);
            wxSize size = dc.GetSize();
            if(drawbuffer.GetSize() != size) {
                drawbuffer = wxBitmap(size, 24);
            }
            Hit recordGeometry = state->selectedRecordHit();

            wxNativePixelData data(drawbuffer);
            wxNativePixelData::Iterator p(data);
            if(state->coordsToDir(0, 0, size.x, size.y).isReal()) {
				
                ImBitmap<float> pdfs(size.x, size.y);
                float maxPdf = 1e-6f;
                if(!state->recordViz.showGeometry) {
                    for(int x = 0; x < size.x; ++x) {
                        for(int y = 0; y < size.y; ++y) {
                            const Vector3 dir = state->coordsToDir(x, y, size.x, size.y);
                            pdfs(x, y) = state->pdf(dir);
                            maxPdf = std::max(maxPdf, pdfs(x, y));
                        }
                    }
                }

				//static int outputCount = 0;
				//std::ofstream dumpStream(std::string(wxString::Format("pdfviz-%d.txt", outputCount++)));
				//dumpStream << pdfs.getWidth() << " " << pdfs.getHeight() << std::endl;
				//for (int i = 0; i < pdfs.pixelCount(); ++i) {
				//	dumpStream << pdfs.getNth(i) << std::endl;
				//}
				//dumpStream.flush();
				//dumpStream.close();

				const float maxMappedIntensity = mapIntensity(state->recordViz.relativeColorMapping ? maxPdf : 255.f);


                for(int x = 0; x < size.x; ++x) {
                    for(int y = 0; y < size.y; ++y) {
                        p.MoveTo(data, x, y);
                        unsigned char r,g,b;
                        if(state->recordViz.showGeometry) {
                            Vector3 color;
                            Hit hit;
                            const Vector3 dir = state->coordsToDir(x, y, size.x, size.y);
                            if(dot(dir, recordGeometry.normal) > 0.f) {
                                int primType;
                                hit = state->rayTrace(recordGeometry.position, dir, color, primType);
                                if(hit.position.isReal()) {
                                    color *= 255.f * abs(dot(hit.normal, dir));
                                    r = color.x;
                                    g = color.y;
                                    b = color.z;
                                } else {
                                    r = g = b = 50;
                                }
                            } else {
                                r = g = b = 30;
                            }
                        } else {
                            Float w = std::min(1.f, mapIntensity(pdfs(x, y))/maxMappedIntensity);
                            Vector3 color;
                            VizTools::getColor2(w, &color);
                            color *= 255.f;
                            r = color.x; g = color.y; b = color.z;
                        }
                        
                        p.Red() = r;
                        p.Green() = g;
                        p.Blue() = b;
                        if(state->recordViz.clickedX == x && state->recordViz.clickedY == y) {
                            p.Red() = p.Blue() = p.Green() = 255;
                        }
                    }
                }
                if( state->recordViz.showPhotons ) {
                    const Hit hit = state->selectedRecordHit();
                    ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> particles;                   
                    const int found = m_samplers->nnquery( hit, particles, this->state );

                    if ( !state->isDisplayUsedParticles ) {
                        for(int i = 0; i < found; ++i) {
                            const Vector3 & incidentDir = m_samplers->getParticleRecVizDirection( particles[i].index, this->state );
                            plotDirection(incidentDir, size, p, data);
                        }
                    } else {
                        IStack<Particle> usedParticles;                            
                        m_samplers->getUsedParticles( particles, found, hit, state, usedParticles );
                        for ( auto it = usedParticles.begin(); it < usedParticles.end(); ++it ) {
                            plotDirection( it->incidentDir, size, p, data );
                        }
                    }
                }
            }

            dc.DrawBitmap(this->drawbuffer, wxPoint(0, 0));
        }    

        void plotDirection( const Vector3 & incidentDir, const wxSize & size, wxNativePixelData::Iterator & p, wxNativePixelData & data ) {
            int x, y;
            state->dirToCoords(-incidentDir, x, y, size.x, size.y);
            if(x >= 0 && y >= 0 && x < size.x && y < size.y) {
                p.MoveTo(data, x, y);
                p.Red() = p.Green() = p.Blue() = 255;
            }
        }
    };
}