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

#define __WXMSW__
#define WXUSINGDLL
#include <wx/wx.h>
#include <wx/glcanvas.h>
#include <gl/GLU.h>
#include "State.h"
#include "../shared/Triangle.h"
#include "sphere.h"
#include "../shared/tools.h"
#include "vizapi.h"
#include "samplersfacade.h"
#include "../shared/Utils.h"

    const float FOV = 55.f;

namespace Importance {

    //enum {
    //    ID_REFRESH,
    //};

    template<class TDistributionModel>
    class CacheViz : public wxGLCanvas {
    public:

        CacheViz(wxWindow* parent) : wxGLCanvas(parent, wxID_ANY, NULL, wxDefaultPosition, wxDefaultSize /*wxSize(800, 600)*/) {
            m_samplerVizAPI         = NULL;
            m_samplers              = NULL;
            m_enviroSamplerVizAPI   = NULL;            
            myGlContext             = NULL;   
            sunFrom = sunTo = Vector3( 0.f );
            sunMode = false;
            
            //wxTimer* timer = new wxTimer(this, ID_REFRESH);
            //timer->Start(100);
            Connect(wxEVT_PAINT, wxPaintEventHandler(CacheViz::OnPaint));
            Connect(wxEVT_MOTION, wxMouseEventHandler(CacheViz::OnMotion));
            Connect(wxEVT_KEY_UP, wxKeyEventHandler(CacheViz::OnKey));
            Connect(wxEVT_MOUSEWHEEL, wxMouseEventHandler(CacheViz::OnScroll));
            Connect(wxEVT_SCROLL_LINEDOWN, wxMouseEventHandler(CacheViz::OnScroll));
            Connect(wxEVT_SCROLL_LINEUP, wxMouseEventHandler(CacheViz::OnScroll));
            Connect(wxEVT_LEFT_DOWN, wxMouseEventHandler(CacheViz::OnClick));
            Connect(wxEVT_LEFT_UP, wxMouseEventHandler(CacheViz::OnClick));
            Connect(wxEVT_RIGHT_DOWN, wxMouseEventHandler(CacheViz::OnClick));
            Connect(wxEVT_RIGHT_UP, wxMouseEventHandler(CacheViz::OnClick));
            Connect(wxEVT_MIDDLE_DOWN, wxMouseEventHandler(CacheViz::OnClick));
            Connect(wxEVT_MIDDLE_UP, wxMouseEventHandler(CacheViz::OnClick));
            Connect(wxEVT_MOUSE_CAPTURE_LOST, wxMouseEventHandler(CacheViz::OnCaptureLost));
        }
        ~CacheViz() {
            if ( m_enviroSamplerVizAPI != NULL ) {
                delete m_enviroSamplerVizAPI;
            }

            if ( m_samplerVizAPI != NULL ) {
                delete m_samplerVizAPI;
            }

            delete myGlContext;
        }

        bool createGlContext() {            
            this->myGlContext = new wxGLContext(this);
            return myGlContext != NULL;
        }

        void setup(State* state, TriangleIterator* geometry, const InitCamera* camera) {
            this->state = state;

            if ( getSampler() != NULL ) {
                m_samplerVizAPI = getSampler()->getVizApi();
            } else {
                m_samplerVizAPI = new SimpleSamplerVizAPI<TDistributionModel>( NULL );
            }

            if ( getEnviroSampler() != NULL ) {
                m_enviroSamplerVizAPI = getEnviroSampler()->getVizApi();
            } else {
                m_enviroSamplerVizAPI = new EnviroSamplerVizAPI<TDistributionModel>( NULL );
            }

            if(geometry) {
                Vector3f v0, v1, v2, n0, n1, n2, color;

                while(geometry->next(v0, v1, v2, n0, n1, n2, color)) {
                    state->triangles.push(VizTriangle());
                    state->triangles.back().v0 = v0;
                    state->triangles.back().v1 = v1;
                    state->triangles.back().v2 = v2;
                    //const Vector3 colorNormal = (n0+n1+n2).getNormalized()/2 + Vector3f(0.5f);
                    const Vector3 colorNormal = ((n0+n1+n2).getNormalized() + Vector3f( 1.f )) * 0.25f;
                    state->triangles.back().color =  colorNormal;/**0.5f + Vector3f(color.l1Norm()/3.f*0.5f);*/
                    state->triangles.back().n0 = n0;
                    state->triangles.back().n1 = n1;
                    state->triangles.back().n2 = n2;

                    this->state->bbox += v0;
                    this->state->bbox += v1;
                    this->state->bbox += v2;
                }
            }            

            if ( getSampler() != NULL ) {
                //state->lock();
                for(int i = 0; i < getSampler()->particles.size(); i += state->everyNthPhoton) {
                    const Particle& particle = getSampler()->particles[i];
                    this->state->bbox += particle.position;
                }
            }

            //state->unlock();
            if(camera) {
                this->camera.origin = camera->origin;
                this->camera.dir = (camera->target-camera->origin).getNormalized();
                this->camera.roll = camera->roll.getNormalized();
                this->camera.xScale = camera->xScale;
                this->initCameraDistance = (camera->target - camera->origin).size();                
            } else {
                this->camera.origin = this->state->bbox.getCenter() + 1.5f*this->state->bbox.size();
                this->camera.dir = (this->state->bbox.getCenter()-this->camera.origin).getNormalized();
                this->camera.roll = Vector3(0, 0, 1);
                this->camera.xScale = 1.f;
                this->initCameraDistance = this->state->bbox.size().size();                
            }
            this->initialCamera = this->camera;

            if ( getEnviroSampler() != NULL && getEnviroSampler()->isInitialized() ) {
                //Float r = (this->state->bbox.point2 - this->state->bbox.point1).size();
                Float r = getEnviroSampler()->getDiscRadius();
                this->state->enviroDummySph = Sphere( getEnviroSampler()->getSceneCenter(), r, true );
            }

            SetFocus();
        }

        void OnKey(wxKeyEvent& evt) {
            Vector3 movement(0.f);
            switch(evt.GetKeyCode()) {
            case 'I':
            case WXK_UP:
                movement = camera.dir;
                break;
            case 'K':
            case WXK_DOWN:
                movement = -camera.dir;
                break;
            case 'J':
            case WXK_LEFT:
                movement = cross(camera.roll, camera.dir).getNormalized();
                break;
            case 'L':
            case WXK_RIGHT:
                movement = -cross(camera.roll, camera.dir).getNormalized();
                break;
            }
            movement *= initCameraDistance/40.f;
            this->camera.origin += movement;

            switch( evt.GetKeyCode() ) {
            case 'F':
                this->camera = this->initialCamera;
                break;
            case 'B':
                this->camera = this->initialCamera;
                this->camera.dir *= -1.f;
                break;
            case 'T': {
                this->camera = this->initialCamera;
                Vector3 axis = Vector3( 1.f, 0.f, 0.f );
                this->camera.dir = Transform::rotation( axis, IMP_PI * 0.5f ) * this->camera.dir; }
                break;
                      
            case 'D': {
                this->camera = this->initialCamera;
                Vector3 axis = Vector3( 1.f, 0.f, 0.f );
                this->camera.dir = Transform::rotation( axis, - IMP_PI * 0.5f ) * this->camera.dir; }
                break;

            case 'S': 
                this->sunMode = !this->sunMode;
                break;

            case 'N':
                this->state->isDisplayNormals = !this->state->isDisplayNormals;
                break;

            /* asterix */
            case 387:
                this->state->isDisplayUsedParticles = !this->state->isDisplayUsedParticles;
                GetParent()->Refresh();
                break;
            }
             
            Refresh();
        }

        void resetView() {
            this->camera = this->initialCamera;
            Refresh();
        }

        void OnScroll(wxMouseEvent& evt) {
            float SCALE = initCameraDistance / 4000.f;
            if(evt.ShiftDown()) {
                SCALE *= 8;
            }
            this->camera.origin += this->camera.dir * evt.GetWheelRotation()*SCALE;
            Refresh();
        }

        IMPORTANCE_INLINE void setSamplers( SamplersFacade<TDistributionModel> * samplers ) {
            m_samplers = samplers;
        }

        IMPORTANCE_INLINE SimpleSampler<TDistributionModel> * getSampler() { return m_samplers->getSampler(); }

        IMPORTANCE_INLINE const SimpleSampler<TDistributionModel> * getSampler() const { return m_samplers->getSampler(); }

        IMPORTANCE_INLINE IEnviroSampler * getEnviroSampler() { return m_samplers->getEnviroSampler(); }

        IMPORTANCE_INLINE const IEnviroSampler * getEnviroSampler() const { return m_samplers->getEnviroSampler(); }

        Particle getParticle( int index ) const {            
            return m_samplers->getParticle( index, this->state );
        }      

        void computeSceneStatistics() {
            double & maxw = state->sceneStats.particleMaxWeight;
            double & minw = state->sceneStats.particleMinWeight;
            double & avgw = state->sceneStats.particleAvgWeight;

            maxw = -std::numeric_limits<double>::max();
            minw = std::numeric_limits<double>::max();
            avgw = 0.0;   

            int count = 0;
            if ( getSampler() != NULL ) {
                for( int i = 0; i < getSampler()->particles.size(); ++i ) {
                    const Particle& particle = getSampler()->particles[i];
                    maxw = std::max( maxw, (double) particle.weight );
                    minw = std::min( minw, (double) particle.weight );
                    avgw += particle.weight;
                }
                count += (int) getSampler()->particles.size();
            }     

            if ( getEnviroSampler() != NULL ) {
                for( int i = 0; i < getEnviroSampler()->getParticles().size(); ++i ) {
                    const Particle& particle = getEnviroSampler()->getParticles()[i];
                    maxw = std::max( maxw, (double) particle.weight );
                    minw = std::min( minw, (double) particle.weight );
                    avgw += particle.weight;
                }
                count += (int) getEnviroSampler()->getParticles().size();
            }
            avgw /= count;
        }

    public:
        void leftClickEvent(wxMouseEvent& evt) {
			std::cout << "CLICK: " << evt.GetX() << " " << evt.GetY() << std::endl;
            std::string label = state->getRecordDesc();
            //recordInfo->SetLabel(label.c_str());

            const float unitWidth = 2*::tan(FOV*PI/180.f/2.f) / GetSize().GetY();
            const Vector3 xVect = cross(camera.dir, camera.roll).getNormalized() * unitWidth;
            const Vector3 yVect = cross(camera.dir, xVect).getNormalized()*unitWidth;
            const Vector3 rayDir = camera.dir + camera.xScale * xVect * (evt.GetX()-GetSize().GetX()/2) + yVect * (evt.GetY()-GetSize().GetY()/2);            

            int primType;
            state->lastClick = state->rayTrace(camera.origin, rayDir, Vector3(), primType);


            if(state->lastClick.position.isReal()) {
                if( this->sunMode ) {
                    /// we are in setting guideline mode
                    if ( evt.ControlDown() ) {
                        this->sunTo = state->lastClick.position;
                    } else {                    
                        this->sunFrom = state->lastClick.position;
                    }
                }

                // resolve whether we hit environment form inside or regular surface
                if( primType != State::PT_SPHERE ) {
                    state->recordViz.vRec = m_samplerVizAPI->getVizHitRecord( state->lastClick, state->recordViz.wharrgarbl, *state );
                    state->recordViz.isEnviro = false;

                    if ( getSampler() != NULL ) {
                        Hit hit = state->recordViz.vRec.toHit();
                        ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> qres;
                        int found = getSampler()->nnquery( hit, qres );
                        dumpNNParticles( hit, qres, getSampler()->particles, found );
                    }
                } else {
                    state->recordViz.vRec = m_enviroSamplerVizAPI->getVizHitRecord( state->lastClick, state->recordViz.wharrgarbl, *state );
                    state->recordViz.isEnviro = true;
                }

                // sample statistics
                ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> particles; 
                Hit hit = state->selectedRecordHit();                
                const int found = nnquery( hit, particles );
                IStack<Particle> nearestParticles( found );
                for ( int i = 0; i < found; ++i ) {
                    nearestParticles[ i ] = getParticle( particles[ i ].index );
                }                
                computeSampleStatistics( nearestParticles, state->smplStats );
                
                
                //////////////////////////////////////////////////////////////////////////
                // Debug print 

                state->print( std::cout, label );
                std::stringstream temp;
                printInfo( temp );
                std::cout << temp.str();
                std::cout.flush();
                OutputDebugStringA(temp.str().c_str());

                std::ofstream file( "vizconsole.txt" );
                printInfo( file );            
                file.close();
            }

            GetParent()->Refresh();
        }


        void OnClick(wxMouseEvent& evt) {
            if(evt.ButtonDown()) {
                CaptureMouse();
                if(evt.GetButton() == wxMOUSE_BTN_LEFT) {
                    leftClickEvent(evt);
                }
            } else {
                IMPORTANCE_ASSERT(evt.ButtonUp());
                ReleaseCapture();
            }
        }
        void OnCaptureLost(wxMouseEvent&) {
            ReleaseCapture();
        }
        void OnMotion(wxMouseEvent& evt) {
            static int lastX = -1;
            static int lastY = -1;
            const int newX = evt.GetX();
            const int newY = evt.GetY();
            if(evt.LeftIsDown()) {
                leftClickEvent(evt);
            } else if(evt.ControlDown() && evt.MiddleIsDown()) {
                this->camera.origin += camera.dir * (lastY-newY)*initCameraDistance/1000.f;
            } else if((evt.AltDown() && evt.MiddleIsDown()) || evt.RightIsDown()) {
                Frame rotator(camera.roll);
                const Vector3 localDir = rotator.toLocal(camera.dir);
                float rotation, elevation;
                uniformSphereInverse(localDir, rotation, elevation);
                rotation += (newX-lastX)/5000.f;
                elevation += (newY-lastY)/1000.f;
                elevation = std::max(0.001f, std::min(0.999f, elevation));
                const Vector3 target = camera.origin+camera.dir*initCameraDistance;
                camera.dir = rotator.toWorld(uniformSphere(rotation, elevation));
                IMPORTANCE_ASSERT(camera.dir.isReal());
                camera.origin = target-camera.dir*initCameraDistance;
            } else if (evt.MiddleIsDown()) {
                const Vector3 xDir = cross(camera.dir, camera.roll).getNormalized();
                const Vector3 yDir = cross(camera.dir, -xDir).getNormalized();
                float SCALE = initCameraDistance / 1000.f;
                if ( evt.ShiftDown() ) {
                    SCALE *= 8;
                }
                camera.origin += -xDir*this->camera.xScale*(newX-lastX)*SCALE + yDir*(newY-lastY)*SCALE;
            }
            lastX = newX;
            lastY = newY;
            Refresh();
        }

        void paintMetric( const int width, const int height, wxPaintDC& paintDc, const CachedSampler<TDistributionModel> * sampler ) {
            if ( sampler == NULL ) {
                return;
            }

            const float unitWidth = 2*::tan(FOV*PI/180.f/2.f) / GetSize().GetY();
            const Vector3 xVect = cross(camera.dir, camera.roll).getNormalized()*unitWidth;
            const Vector3 yVect = cross(camera.dir, xVect).getNormalized()*unitWidth;
            
            wxBitmap drawbuffer(width, height, 24);
            wxNativePixelData data(drawbuffer);

            //#pragma omp parallel
            for(int x = 0; x < width; ++x) {
                if(x % 10) {
                    std::cout << ".";
                }
                wxNativePixelData::Iterator p(data);
                for(int y = 0; y < height; ++y) {
                    p.MoveTo(data, x, y);
                    const Vector3 rayDir = (camera.dir + xVect*(x-width/2) + yVect*(y-height/2)).getNormalized();
                    int primType;
                    Hit hit = state->rayTrace(camera.origin, rayDir, Vector3(), primType);
                    float val;
                    if(hit.normal.isReal() && (val = sampler->getInterpolationMetric(hit, state->displayedMetric)) > 0) {
                        val = std::min(val, 1.f);
                        Vector3 color;
                        VizTools::getColor2(val, &color);
                        color *= 255.f;
                        p.Red() = color.x;
                        p.Green() = color.y;
                        p.Blue() = color.z;
                    } else {
                        p.Red() = 0;
                        p.Green() = 0;
                        p.Blue() = 0;
                    }
                }
            }
            paintDc.DrawBitmap(drawbuffer, wxPoint(0, 0));
        }

        void paintAllPhotons( const SimpleSampler<TDistributionModel> * sampler ) {
            if ( sampler == NULL ) {
                return;
            }

            const float PHOTON_COLOR = 0.8f;
            Vector3 color;
            glColor3f(PHOTON_COLOR, PHOTON_COLOR, PHOTON_COLOR);
            glPointSize(1.0f);
            glEnable(GL_POINT_SMOOTH);

            glBegin(GL_POINTS);
            for(int i = 0; i < sampler->particles.size(); i += state->everyNthPhoton) {
                const Particle& particle = sampler->particles[i];
                if ( state->isPhotonsFalseColor ) {
                    VizTools::getColor2( particle.weight / state->sceneStats.particleMaxWeight, &color );
                    glColor3f( color.x, color.y, color.z );
                }
                glVertex3f(particle.position.x, particle.position.y, particle.position.z);
            }
            glEnd();
        }

        void paintMarks() {
            if ( state->marks == NULL ) {
                return;
            }

            const float MARK_COLOR = 0.8f;
            Vector3 color;
            glColor3f(MARK_COLOR, MARK_COLOR, 0.f);
            glPointSize(2.0f);
            glEnable(GL_POINT_SMOOTH);

            glBegin(GL_POINTS);
            for( auto mark = state->marks->begin(); mark < state->marks->end(); ++mark ) {
                const Vector3 & pos = mark->position;
                glVertex3f( pos.x, pos.y, pos.z );
            }
            glEnd();
        }

        void paintAllEnviroImportons( const IEnviroSampler * sampler ) {
            if ( sampler == NULL ) {
                return;
            }

            const float IMPORTON_COLOR = 0.8f;
            Vector3 color;
            glColor3f(IMPORTON_COLOR, IMPORTON_COLOR, IMPORTON_COLOR);
            glPointSize(1.0f);
            glEnable(GL_POINT_SMOOTH);

            glBegin(GL_POINTS);
            float r = state->enviroDummySph.getRadius();
            for(int i = 0; i < sampler->getParticles().size(); i += state->everyNthPhoton) {
                const Particle& particle = sampler->getParticles()[ i ];
                if ( state->isPhotonsFalseColor ) {
                    VizTools::getColor2( particle.weight / state->sceneStats.particleMaxWeight, &color );
                    glColor3f( color.x, color.y, color.z );
                }
                //Vector3 pos = sampler->getSceneCenter() - r * particle.position;
                Vector3 pos = sampler->getIncidentPosition( particle );
                glVertex3f(pos.x, pos.y, pos.z);
            }
            glEnd();
        }

        void paintAABB() {
            GLfloat left, right, front, back, top, bottom;
            left    = this->state->bbox.point1.x;
            right   = this->state->bbox.point2.x;
            top     = this->state->bbox.point1.z;
            bottom  = this->state->bbox.point2.z;
            back    = this->state->bbox.point1.y;
            front   = this->state->bbox.point2.y;

            glBegin(GL_LINES);
                glColor3f( 0.5f, 0.5f, 0.5f );
                glVertex3f( left, front, top );
                glVertex3f( left, back, top );                
                glVertex3f( left, back, top );
                glVertex3f( left, back, bottom );
                glVertex3f( left, back, bottom );
                glVertex3f( left, front, bottom );
                glVertex3f( left, front, bottom );
                glVertex3f( left, front, top );

                glVertex3f( right, front, top );
                glVertex3f( right, back, top );                
                glVertex3f( right, back, top );
                glVertex3f( right, back, bottom );
                glVertex3f( right, back, bottom );
                glVertex3f( right, front, bottom );
                glVertex3f( right, front, bottom );
                glVertex3f( right, front, top );

                glVertex3f( right, front, top );
                glVertex3f( right, back, top );
                glVertex3f( right, back, top );
                glVertex3f( left,  back, top );
                glVertex3f( left,  back, top );
                glVertex3f( left,  front, top );
                glVertex3f( left,  front, top );
                glVertex3f( right, front, top );

                glVertex3f( right, front, bottom );
                glVertex3f( right, back, bottom );
                glVertex3f( right, back, bottom );
                glVertex3f( left,  back, bottom );
                glVertex3f( left,  back, bottom );
                glVertex3f( left,  front, bottom );
                glVertex3f( left,  front, bottom );
                glVertex3f( right, front, bottom );


                glVertex3f( right, front, bottom );
                glVertex3f( right, front, top );
                glVertex3f( right, front, top );
                glVertex3f( left,  front, top );
                glVertex3f( left,  front, top );
                glVertex3f( left,  front, bottom );
                glVertex3f( left,  front, bottom );
                glVertex3f( right, front, bottom );

                glVertex3f( right, back, bottom );
                glVertex3f( right, back, top );
                glVertex3f( right, back, top );
                glVertex3f( left,  back, top );
                glVertex3f( left,  back, top );
                glVertex3f( left,  back, bottom );
                glVertex3f( left,  back, bottom );
                glVertex3f( right, back, bottom );
            glEnd();
        }

        int nnquery( const Hit & hit, ImStaticArray<KdQueryResult,MAX_KNN_PARTICLES> & particles ) {
            return m_samplers->nnquery( hit, particles, this->state );
        }

        void paintUsedParticles() {
            ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> particles; 
            Hit hit = state->selectedRecordHit();
            if(hit.position.isReal() && state->showActivePhotons) {
                const int found = nnquery( hit, particles );
                IStack<Particle> usedParticles;
                m_samplers->getUsedParticles( particles, found, hit, state, usedParticles );

                glBegin(GL_POINTS);
                glColor3f(1, 1, 1);
                for ( auto it = usedParticles.begin(); it < usedParticles.end(); ++it ) {
                    const Particle & p = *it;
                    glVertex3f(p.position.x, p.position.y, p.position.z);
                }
                glEnd();

                if ( state->photonLength > 0.f ) {
                    glLineWidth(0.5f);
                    glBegin(GL_LINES);
                    glColor3f(1, 1, 1);
                    Vector3 color;
                    Vector3 avgDir(0.f);
                    Float avgDist = 0.f;
                    for(auto it = usedParticles.begin(); it < usedParticles.end(); ++it ) {                        
                        const Particle & p = *it;
                        avgDir += p.incidentDir;
                        avgDist += p.distance;
                        if ( state->isPhotonsFalseColor ) {
                            VizTools::getColor2( p.weight / state->sceneStats.particleMaxWeight, &color );
                            glColor3f( color.x, color.y, color.z );
                        }
                        glVertex3f(p.position.x, p.position.y, p.position.z);
                        const Vector3 pEnd = p.position-p.incidentDir*state->photonLength*p.distance;
                        glVertex3f(pEnd.x, pEnd.y, pEnd.z);
                    }
                    avgDir = avgDir / found;
                    avgDist /= found;
                    const Vector3 pEnd = hit.position-avgDir*4*state->photonLength*avgDist;
                    glLineWidth( 1.5f );
                    glColor3f(1.f, 0.f, 0.f);
                    glVertex3f(hit.position.x, hit.position.y, hit.position.z);
                    glVertex3f(pEnd.x, pEnd.y, pEnd.z);
                    glEnd();
                }
            }
        }

        void paintActivePhotons() {
            ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> particles; 
            Hit hit = state->selectedRecordHit();
            if(hit.position.isReal() && state->showActivePhotons) {
                const int found = nnquery( hit, particles );                

                glBegin(GL_POINTS);
                glColor3f(1, 1, 1);
                IStack<Particle> debugParticles;
                for(int i = 0; i < found; ++i) {
                    const Particle photon = getParticle( particles[i].index );
                    debugParticles.push( photon );
                    glVertex3f(photon.position.x, photon.position.y, photon.position.z);
                }
                glEnd();
                if(state->photonLength > 0.f) {
                    glLineWidth(0.5f);
                    glBegin(GL_LINES);
                    glColor3f(1, 1, 1);
                    Vector3 color;
                    Vector3 avgDir(0.f);
                    Float avgDist = 0.f;
                    for(int i = 0; i < found; ++i) {
                        const Particle particle = getParticle( particles[i].index );
                        avgDir += particle.incidentDir;
                        avgDist += particle.distance;
                        if ( state->isPhotonsFalseColor ) {
                            VizTools::getColor2( particle.weight / state->sceneStats.particleMaxWeight, &color );
                            glColor3f( color.x, color.y, color.z );
                        }
                        glVertex3f(particle.position.x, particle.position.y, particle.position.z);
                        const Vector3 pEnd = particle.position-particle.incidentDir*state->photonLength*particle.distance;
                        glVertex3f(pEnd.x, pEnd.y, pEnd.z);
                    }
                    //avgDir = avgDir / found;
                    //avgDist /= found;
                    //const Vector3 pEnd = hit.position-avgDir*4*state->photonLength*avgDist;
                    //glLineWidth( 1.5f );
                    //glColor3f(1.f, 0.f, 0.f);
                    //glVertex3f(hit.position.x, hit.position.y, hit.position.z);
                    //glVertex3f(pEnd.x, pEnd.y, pEnd.z);
                    glEnd();
                }
            }
        }

        void OnPaint(wxPaintEvent&) {
            
            wxPaintDC paintDc(this);                            

            wglMakeCurrent(paintDc.GetHDC(), myGlContext->GetGLRC());
            const Vector3 & c = this->state->bgColor;
            glClearColor(c.x, c.y, c.z, 0.0f);
            glEnable(GL_DEPTH_TEST);
            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

            glLoadIdentity();

            glViewport(0, 0, (GLint)GetSize().x, (GLint)GetSize().y);
            gluPerspective(FOV, GetSize().GetX()/float(GetSize().GetY()), this->state->bbox.size().size()*0.01f, this->state->bbox.size().size()*5.f);
            const Vector3 target = camera.origin + camera.dir;
            glScalef( camera.xScale, 1.f, 1.f );
            gluLookAt(camera.origin.x, camera.origin.y, camera.origin.z, target.x, target.y, target.z, camera.roll.x, camera.roll.y, camera.roll.z);

            if(state->showGeometry && !state->showEnvironmentOnly) {                
                glBegin(GL_TRIANGLES);
                for(int i = 0; i < state->triangles.size(); ++i) {
                    const VizTriangle& act = state->triangles[i];                    
                    if ( !state->isGeometryBgColor ) {
                        glColor3f(act.color.x, act.color.y, act.color.z);
                    } else {
                        const Vector3 & c = this->state->bgColor;
                        glColor3f(c.x, c.y, c.z);
                    }
                    glVertex3f(act.v0.x, act.v0.y, act.v0.z);
                    glVertex3f(act.v1.x, act.v1.y, act.v1.z);
                    glVertex3f(act.v2.x, act.v2.y, act.v2.z);
                }
                glEnd();
            }

            if(state->showCameraRatio && state->selectedRecordHit().position.isReal()) {
                Hit hit = state->selectedRecordHit();
                Float r = state->config.cache.maxRadius * state->config.pixel2WorldRatio(hit, getSampler()->camera);
                VizAPI::drawRadius(hit.position, hit.normal, r, Vector3(0,1,0));
            }

            if(state->showPhotons && !state->showEnvironmentOnly) {
                paintAllPhotons( getSampler() );
            }
            
            if( state->interpolateOnClick && state->recordViz.vRec.isValid ) {
                if ( !state->isDisplayUsedParticles) {
                    paintActivePhotons();
                } else {
                    paintUsedParticles();
                }
            }

            if ( !state->showEnvironmentOnly ) {
                m_samplerVizAPI->draw( *state );
            }
            m_enviroSamplerVizAPI->draw( *state );

            if ( state->showPhotons ) {
                paintAllEnviroImportons( getEnviroSampler() );
            }

            if ( state->isShowMarks ) {
                paintMarks();
            }

            paintAABB();

            if(!isNaN(state->lastClick.position.x)) {
                glPointSize(4.f);
                glColor3f(0, 0, 1);
                glBegin(GL_POINTS);
                glVertex3f(state->lastClick.position.x, state->lastClick.position.y, state->lastClick.position.z);
                glEnd();
                glBegin(GL_LINES);
                glVertex3f(state->lastClick.position.x, state->lastClick.position.y, state->lastClick.position.z);
                glVertex3f(state->lastClick.position.x+state->lastClick.normal.x, 
                    state->lastClick.position.y+state->lastClick.normal.y, 
                    state->lastClick.position.z+state->lastClick.normal.z);
                glEnd();
            }

            if ( this->sunMode ) {
                glColor3f(1, 0, 0);
                glPointSize(4.f);
                glBegin(GL_LINES);
                glVertex3f(this->sunFrom.x, this->sunFrom.y, this->sunFrom.z);
                glVertex3f(this->sunTo.x, this->sunTo.y, this->sunTo.z);
                glEnd();
            }

            glFlush();      
            ::SwapBuffers(paintDc.GetHDC());        
         }      

    private:
        void printInfo( std::ostream & out ) const {
            int nEnvParticles = getEnviroSampler() != NULL ? getEnviroSampler()->getParticles().size() : 0;

            out << "-------------------------------------------------------------------------------------------" << std::endl << std::endl;

            std::string tmp = state->interpolateOnClick ? "Click" : "Record";
            out << tmp << " position: " << state->recordViz.vRec.position.toString() << std::endl;
            out << tmp << " normal: " << state->recordViz.vRec.normal.toString() << std::endl;
            out << "Number of all particles (including environment): " << ( getSampler()->particles.size() + nEnvParticles ) << std::endl;
            out << "Number of environment particles (out of all): " << nEnvParticles << std::endl;
            out << "Particle max weight (all particles): " << state->sceneStats.particleMaxWeight << std::endl;
            out << "Particle min weight (all particles): " << state->sceneStats.particleMinWeight << std::endl;
            out << "Particle avg weight (all particles): " << state->sceneStats.particleAvgWeight << std::endl;
            out << std::endl;

            out << "Statistics of selected particles:" << std::endl << state->smplStats.toString() << std::endl;        

            if( !state->interpolateOnClick && state->recordViz.vRec.isValid ) {                
                const TDistributionModel * dist = NULL;
                ImmediateDistribution<TDistributionModel> * tmp = dynamic_cast<ImmediateDistribution<TDistributionModel>*>( state->recordViz.vRec.pDistr );
                if ( tmp != NULL ) {
                    dist = &tmp->model;
                }

                if ( dist != NULL ) {
                    out << dist->toString() << std::endl;
                    const EmInfo * info = dist->getInfo();				    
                    if ( info != NULL ) {
                        out << "Iterations: " << info->nIteration << std::endl;
                        out << "Photons found: " << info->nPhotonsFound << std::endl;
                        out << "Photons used: " << info->nPhotonsUsed << std::endl;
                        if ( info->lkhFunction.size() > 0 ) {
                            out << "Likelihood: " << info->lkhFunction[ info->lkhFunction.size() - 1 ];
                        }
                        out << "(Particles used for fitting)" << std::endl;
                        out << "Clamped weights: " << info->clampedWeights << std::endl;
                        out << "Min weight: " << info->minWeight << std::endl;
                        out << "Max weight: " << info->maxWeight << std::endl;
                        out << "Avg weight: " << info->avgWeight << std::endl;
                    }
                    const GaussianMixtureSSEHemisphere * gmmDist = dynamic_cast<const GaussianMixtureSSEHemisphere*>( dist );
                    if ( gmmDist != NULL ) {
                        out << "Progressive k: " << gmmDist->k << std::endl;
                        out << "Photons used to fit: " << gmmDist->m_cacheStats.photonCount << std::endl;
                        out << "Photons search radius: " << gmmDist->m_cacheStats.maxDistance << std::endl;
                    }
                }
            }

            out << "Visualization guideline: " << (this->sunTo - this->sunFrom).getNormalized().toString() << std::endl;
        }

    protected:
        wxGLContext* myGlContext;
        State* state;       

        struct Camera {
            Vector3 origin;
            Vector3 dir;
            Vector3 roll;
            float xScale;
        };
        Camera camera, initialCamera;
        
        float initCameraDistance;
        bool sunMode;
        Vector3 sunFrom, sunTo;

    private:
        SamplersFacade<TDistributionModel> * m_samplers;
        VizAPI * m_enviroSamplerVizAPI;
        VizAPI * m_samplerVizAPI;
    };
}