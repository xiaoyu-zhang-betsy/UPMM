/*
    This file is part of a demo implementation of an importance sampling technique
    described in the "On-line Learning of Parametric Mixture Models for Light Transport Simulation"
    (SIGGRAPH 2014) paper.
    The implementation is based on Mitsuba, a physically based rendering system.

    Copyright (c) 2014 by Jiri Vorba, Ondrej Karlik, Martin Sik.
    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once 

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/guiding_config.h>
#include <mitsuba/render/guided_particletracing.h>
#include "LibImportance/LibImportance.h"

#define MTS_VECTOR(v) Vector((v).x, (v).y, (v).z)
#define IMP_VECTOR3(v) Importance::Vector3((v).x, (v).y, (v).z)

MTS_NAMESPACE_BEGIN

/* Translate Mitsuba Intersection to importance lib Hit */
inline Importance::Hit itsToHit( const Intersection & its ) {
    Importance::Hit hit;
    hit.position.x = its.p.x;
    hit.position.y = its.p.y;
    hit.position.z = its.p.z;
    hit.normal.x = its.shFrame.n.x;
    hit.normal.y = its.shFrame.n.y;
    hit.normal.z = its.shFrame.n.z;
    return hit;
}

/************************************************************************/
/*  Importance library stuff                                            */
/************************************************************************/

/// Conversion of photons from Mitsuba into the depth of Importance sampling library
class PhotonsIterator : public Importance::InputIterator {
public:
    PhotonsIterator( const ref<PhotonMap> photonMap ) 
        : m_photonMap( photonMap ), m_pos( 0 )  {};

    virtual Importance::Vector3 position() const {
        const Point & pos = (*m_photonMap)[ m_pos ].getPosition();
        return Importance::Vector3( pos.x, pos.y, pos.z );
    }

    virtual Importance::Vector3 incidentDirection() const {
        const Vector & dir = (*m_photonMap)[ m_pos ].getDirection();
        return Importance::Vector3( dir.x, dir.y, dir.z );
    }

    virtual Float weight() const {
        //return (*m_photonMap)[ m_pos ].getPower().average();
        IMPORTANCE_ASSERT((*m_photonMap)[ m_pos ].getPower().max() > 0.0f );
        return (*m_photonMap)[ m_pos ].getPower().max();        
    }

    virtual Importance::Vector3 normal() const {
        const Normal & n = (*m_photonMap)[ m_pos ].getNormal();
        return Importance::Vector3( n.x, n.y, n.z );
    }

    virtual Float distance() const {
        return (*m_photonMap)[ m_pos ].getDistance();
    }

    virtual void next() {
        SAssert( m_pos < m_photonMap->size() );
        m_pos++;
    }

    virtual bool isValid() const {
        return m_pos < m_photonMap->size();
    }
private:
    ref<PhotonMap> m_photonMap;
    size_t m_pos;
};

class EnviroMapInterface : public Emitter, public Importance::EnviroMap {
public:

    /************************************************************************/
    /*     Importance::EnviroMap                                            */
    /************************************************************************/

	/** Returns env. map bounding sphere center */
	virtual Importance::Vector3 getSceneCenter() const {
		return IMP_VECTOR3( this->getEnviroSphere().center );
	}

	/** Returns env. map bounding sphere radius */
	virtual float getDiscRadius() const {
		return this->getEnviroSphere().radius;
	}

	/** Returns a pdf of direction from environment */
	virtual float _dirPdf( const Importance::Vector3 dirFromEnviro ) const {
		Vector dir = MTS_VECTOR( dirFromEnviro );
		PositionSamplingRecord pRec;
		pRec.time = 0.f;
		DirectionSamplingRecord dRec( dir );
		return this->pdfDirection( dRec, pRec );
	}    

    /************************************************************************/
    /*       Emitter                                                        */
    /************************************************************************/

    /** Returns env. map bounding sphere */
    virtual const BSphere & getEnviroSphere() const = 0;

    /************************************************************************/
    /*   class EnviroMapInterface                                           */
    /************************************************************************/

	/** Destructor */
	virtual ~EnviroMapInterface() {}

	/** Sets env. map environment sampler */
	void setEnviroSampler( Importance::IEnviroSampler * enviroSampler) const {
		m_enviroSampler = enviroSampler;
	}

protected:
	EnviroMapInterface(const Properties &props) : Emitter(props), m_enviroSampler(0) {}
	EnviroMapInterface(Stream *stream, InstanceManager *manager) : Emitter(stream,manager), m_enviroSampler(0) {}

	/** Importance library sampler of environment map */
	mutable Importance::IEnviroSampler * m_enviroSampler;
};

/// Information for Importance caching inside the Importance lib
class MitsubaImportanceCamera : public Importance::Camera {
private:
    ref<const mitsuba::PerspectiveCamera> m_cam;
    /// Average of pixel edge sizes in the image plane
    Float m_avgPixelEdgeImgPlane;
    Point m_position;
public:
    MitsubaImportanceCamera( ref<const mitsuba::Sensor> cam ) {
        m_cam = dynamic_cast<const mitsuba::PerspectiveCamera*>( cam.get() );
        IMPORTANCE_ASSERT( m_cam != NULL );   
        const mitsuba::Film * film = m_cam->getFilm();
        Float invAspect = (Float) film->getSize().y / (Float) film->getSize().x;
        m_avgPixelEdgeImgPlane = ( std::tan( degToRad( m_cam->getXFov() ) * 0.5f ) * (1.f + invAspect ) ) /
            film->getSize().x;
    }

    // This is not used in default setting so it does not need to be implemented
    virtual Float pixelToWorldRatio( const Importance::Hit& hit ) const {   
        Float time = 0.f;
        Point center = m_cam->getWorldTransform()->eval( time ).transformAffine( Point( 0.0f ) );
        Importance::Vector3 dir = IMP_VECTOR3( center ) - hit.position;
        Float dist = dir.length();
        dir /= dist;        
        Float cosWorldNormal = std::abs( Importance::dot( hit.normal, dir ) );
        cosWorldNormal = std::max( 0.2f, cosWorldNormal );
        /// it is just the approx. - dist si not used correctly
        Float res = (m_avgPixelEdgeImgPlane * dist ) / std::sqrt( cosWorldNormal );
        IMPORTANCE_ASSERT( res > 0.f );        
        return res;
    }
};

class EnviroImportonsIterator : public Importance::InputIterator {
public:
    EnviroImportonsIterator( const BgParticlesVec * particlesVec ) 
        : m_particlesVec( particlesVec ), m_pos( 0 ) {}

    virtual Importance::Vector3 position() const {
        const Point & pos = m_particlesVec->m_data[ m_pos ].getPosition();
        return Importance::Vector3( pos.x, pos.y, pos.z );
    }

    virtual Importance::Vector3 incidentDirection() const {
        const Vector & dir = m_particlesVec->m_data[ m_pos ].getDirection();
        return Importance::Vector3( dir.x, dir.y, dir.z );
    }

    virtual Float weight() const {
        //return (*m_photonMap)[ m_pos ].getPower().average();
        IMPORTANCE_ASSERT(m_particlesVec->m_data[ m_pos ].getPower().max() > 0.0f );
        return m_particlesVec->m_data[ m_pos ].getPower().max();        
    }

    virtual Importance::Vector3 normal() const {
		const Vector & pos = m_particlesVec->m_data[ m_pos ].getNormal();
		return Importance::Vector3( pos.x, pos.y, pos.z );
    }

    virtual Float distance() const {
        return m_particlesVec->m_data[ m_pos ].getDistance();
    }

    virtual void next() {
        SAssert( m_pos < m_particlesVec->m_data.size() );
        m_pos++;
    }

    virtual bool isValid() const {
        return m_pos < m_particlesVec->m_data.size();
    }
private:
    ref<const BgParticlesVec> m_particlesVec;
    size_t m_pos;
};


MTS_NAMESPACE_END