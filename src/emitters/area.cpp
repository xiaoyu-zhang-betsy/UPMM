/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include <mitsuba/render/emitter.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{area}{Area light}
 * \icon{emitter_area}
 * \order{2}
 * \parameters{
 *     \parameter{radiance}{\Spectrum}{
 *         Specifies the emitted radiance in units of
 *         power per unit area per unit steradian.
 *     }
 *     \parameter{samplingWeight}{\Float}{
 *         Specifies the relative amount of samples
 *         allocated to this emitter. \default{1}
 *     }
 * }
 *
 * This plugin implements an area light, i.e. a light source that emits
 * diffuse illumination from the exterior of an arbitrary shape.
 * Since the emission profile of an area light is completely diffuse, it
 * has the same apparent brightness regardless of the observer's viewing
 * direction. Furthermore, since it occupies a nonzero amount of space, an
 * area light generally causes scene objects to cast soft shadows.
 *
 * When modeling scenes involving area lights, it is preferable
 * to use spheres as the emitter shapes, since they provide a
 * particularly good direct illumination sampling strategy (see
 * the \pluginref{sphere} plugin for an example).
 *
 * To create an area light source, simply instantiate the desired
 * emitter shape and specify an \code{area} instance as its child:
 *
 * \vspace{4mm}
 * \begin{xml}
 * <!-- Create a spherical light source at the origin -->
 * <shape type="sphere">
 *     <emitter type="area">
 *         <spectrum name="radiance" value="1"/>
 *     </emitter>
 * </shape>
 * \end{xml}
 */

class AreaLight : public Emitter {
public:
	AreaLight(const Properties &props) : Emitter(props) {
		m_type |= EOnSurface;

		if (props.hasProperty("toWorld"))
			Log(EError, "Found a 'toWorld' transformation -- this is not "
				"allowed -- the area light inherits this transformation from "
				"its parent shape");

		m_radiance = props.getSpectrum("radiance", Spectrum::getD65());
		m_power = Spectrum(0.0f); /// Don't know the power yet
	}

	AreaLight(Stream *stream, InstanceManager *manager)
		: Emitter(stream, manager) {
		m_radiance = Spectrum(stream);
		m_power = Spectrum(stream);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Emitter::serialize(stream, manager);
		m_radiance.serialize(stream);
		m_power.serialize(stream);
	}

	Spectrum samplePosition(PositionSamplingRecord &pRec,
			const Point2 &sample, const Point2 *extra) const {
		m_shape->samplePosition(pRec, sample);
		return m_power;
	}

	Spectrum evalPosition(const PositionSamplingRecord &pRec) const {
		return m_radiance * M_PI;
	}

	Spectrum eval(const Intersection &its, const Vector &d) const {
		if (dot(its.shFrame.n, d) <= 0)
			return Spectrum(0.0f);
		else
			return m_radiance;
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {
		return m_shape->pdfPosition(pRec);
	}

	Spectrum sampleDirection(DirectionSamplingRecord &dRec,
			PositionSamplingRecord &pRec,
			const Point2 &sample, const Point2 *extra) const {
		Vector local = Warp::squareToCosineHemisphere(sample);
		dRec.d = Frame(pRec.n).toWorld(local);
		dRec.pdf = Warp::squareToCosineHemispherePdf(local);
		dRec.measure = ESolidAngle;
		return Spectrum(1.0f);
	}

	Spectrum evalDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {
		Float dp = dot(dRec.d, pRec.n);

		if (dRec.measure != ESolidAngle || dp < 0)
			dp = 0.0f;

		return Spectrum(INV_PI * dp);
	}

	Float pdfDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {
		Float dp = dot(dRec.d, pRec.n);

		if (dRec.measure != ESolidAngle || dp < 0)
			dp = 0.0f;

		return INV_PI * dp;
	}

	Spectrum sampleRay(Ray &ray,
			const Point2 &spatialSample,
			const Point2 &directionalSample,
			Float time) const {
		PositionSamplingRecord pRec(time);
		m_shape->samplePosition(pRec, spatialSample);
		Vector local = Warp::squareToCosineHemisphere(directionalSample);
		ray.setTime(time);
		ray.setOrigin(pRec.p);
		ray.setDirection(Frame(pRec.n).toWorld(local));
		return m_power;
	}

	Spectrum sampleDirect(DirectSamplingRecord &dRec,
			const Point2 &sample) const {
		m_shape->sampleDirect(dRec, sample);

		/* Check that the emitter and reference position are oriented correctly
		   with respect to each other. Note that the >= 0 check
		   for 'refN' is intentional -- those sampling requests that specify
		   a reference point within a medium or on a transmissive surface
		   will set dRec.refN = 0, hence they should always be accepted. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0 && dRec.pdf != 0) {
			return m_radiance / dRec.pdf;
		} else {
			dRec.pdf = 0.0f;
			return Spectrum(0.0f);
		}
	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
		/* Check that the emitter and receiver are oriented correctly
		   with respect to each other. */
		if (dot(dRec.d, dRec.refN) >= 0 && dot(dRec.d, dRec.n) < 0) {
			return m_shape->pdfDirect(dRec);
		} else {
			return 0.0f;
		}
	}

	void setParent(ConfigurableObject *parent) {
		Emitter::setParent(parent);

		if (parent->getClass()->derivesFrom(MTS_CLASS(Shape))) {
			Shape *shape = static_cast<Shape *>(parent);
			if (m_shape == shape || shape->isCompound())
				return;

			if (m_shape != NULL)
				Log(EError, "An area light cannot be parent of multiple shapes");

			m_shape = shape;
			m_shape->configure();
			m_power = m_radiance * M_PI * m_shape->getSurfaceArea();
		} else {
			Log(EError, "An area light must be child of a shape instance");
		}
	}

	AABB getAABB() const {
		return m_shape->getAABB();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "AreaLight[" << endl
			<< "  radiance = " << m_radiance.toString() << "," << endl
			<< "  samplingWeight = " << m_samplingWeight << "," << endl
			<< "  surfaceArea = ";
		if (m_shape)
			oss << m_shape->getSurfaceArea();
		else
			oss << "<no shape attached!>";
		oss << "," << endl
		    << "  medium = " << indent(m_medium.toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	Float gatherAreaPdf(PositionSamplingRecord pRec, Point gatherPosition, Float gatherRadius, Vector4 &bbox, Vector4 *bboxd) const{
// 		Vector local = Warp::squareToCosineHemisphere(sample);
// 		dRec.d = Frame(pRec.n).toWorld(local);
// 		dRec.pdf = Warp::squareToCosineHemispherePdf(local);
// 		dRec.measure = ESolidAngle;
// 		return Spectrum(1.0f);

		bbox = Vector4(0.f, 0.5f * M_PI, 0.f, 2.f * M_PI);

		Vector wo = Frame(pRec.n).toLocal(gatherPosition - pRec.p);		
		Vector dir = wo;
		Float dis = dir.length();
		if (dis < gatherRadius) return 1.f;
		dir /= dis;
		Float dTheta = acos(sqrt(dis * dis - gatherRadius * gatherRadius) / dis);
		Float theta = acos(dir.z);
		Float theta0 = theta - dTheta;
		Float theta1 = std::min(0.5f * M_PI, theta + dTheta);
		Float prob = 0.f;
		if (theta0 < 0.f){
			// sample the full sphere cap of polar
			Float cos0 = cos(2.f * theta0);
			Float cos1 = cos(2.f * theta1);
			prob = 0.5f * (cos0 - cos1);
			bbox.x = theta0; bbox.y = theta1;
		}
		else{
			// sample the bbox of the cone
			Vector dirProj = Vector(wo.x, wo.y, 0.f);
			Float disProj = dirProj.length();
			Float dPhi = acos(sqrt(disProj * disProj - gatherRadius * gatherRadius) / disProj);
			Float phi = atan2(dirProj.y, dirProj.x);
			Float cos0 = cos(2.f * theta0);
			Float cos1 = cos(2.f * theta1);
			prob = 0.5f * dPhi * (cos0 - cos1) / M_PI;
			bbox.x = theta0; bbox.y = theta1;
			bbox.z = phi - dPhi; bbox.w = phi + dPhi;
		}
		return prob;
	}
	Vector sampleGatherArea(DirectionSamplingRecord &dRec, PositionSamplingRecord pRec, Point gatherPosition, Float gatherRadius, Point2 sample, Vector4 bbox, Vector4 bboxd) const{
		Vector dir;
		if (bbox.x == 0.f && bbox.y == 0.5f * M_PI && bbox.z == 0.f && bbox.w == 2.f * M_PI){
			// uniform sampling
			dir = Warp::squareToCosineHemisphere(sample);
		}
		else if (bbox.z == 0.f && bbox.w == 2.f * M_PI){
			// Sampling the whole sphere cap
			Point2 smp = sample;
			Float r0 = sin(std::max((Float)0.0, bbox.x));
			Float r1 = sin(bbox.y);
			smp.x = smp.x * 2.f - 1.f;
			smp.y = smp.y * 2.f - 1.f;
			Float baser = r0;
			if (smp.x * smp.x > smp.y * smp.y){
				Float r2r1 = smp.y / smp.x;
				if (smp.x < 0.f) baser = -baser;
				smp.x = smp.x * (r1 - r0) + baser;
				smp.y = r2r1 * smp.x;
			}
			else{
				Float r1r2 = smp.x / smp.y;
				if (smp.y < 0.f) baser = -baser;
				smp.y = smp.y * (r1 - r0) + baser;
				smp.x = r1r2 * smp.y;
			}
			smp.x = (smp.x + 1.f) * 0.5f;
			smp.y = (smp.y + 1.f) * 0.5f;
			dir = Warp::squareToCosineHemisphere(smp);
		}
		else{
			// sampling a bbox in theta-phi space
			Float sin0 = sin(std::max((Float)0.0, bbox.x));
			Float sin1 = sin(bbox.y);
			Float phi0 = bbox.z;
			Float phi1 = bbox.w;
			Point2 smp = sample;
			Float isin = smp.x * (sin1 - sin0) + sin0;
			Float iphi = smp.y * (phi1 - phi0) + phi0;
			Float r = isin;
			Float z = sqrt(1.f - isin * isin);
			Float cosPhi, sinPhi;
			math::sincos(iphi, &sinPhi, &cosPhi);
			if (EXPECT_NOT_TAKEN(z == 0))
				z = 1e-10f;
			dir = Vector(r * cosPhi, r* sinPhi, z);
		}
		dir = Frame(pRec.n).toWorld(dir);
		return dir;
	}

	Float gatherAreaPdf(Vector wi, Vector wo, Float gatherRadius, Vector4 &bbox, Vector4 *bboxd) const{
		
	}

	Vector sampleGatherArea(Vector wi, Vector wo, Float gatherRadius, Point2 sample, Vector4 bbox, Vector4 bboxd) const{
		
	}

	MTS_DECLARE_CLASS()
protected:
	Spectrum m_radiance, m_power;
};

// ================ Hardware shader implementation ================

class AreaLightShader : public Shader {
public:
	AreaLightShader(Renderer *renderer, const Spectrum &radiance)
		: Shader(renderer, EEmitterShader), m_radiance(radiance) {
	}

	void resolve(const GPUProgram *program, const std::string &evalName,
			std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_radiance", false));
	}

	void generateCode(std::ostringstream &oss, const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_radiance;" << endl
			<< endl
			<< "vec3 " << evalName << "_area(vec2 uv) {" << endl
			<< "    return " << evalName << "_radiance * pi;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_dir(vec3 wo) {" << endl
			<< "    if (cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return vec3(inv_pi);" << endl
			<< "}" << endl;
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_radiance);
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum m_radiance;
};

Shader *AreaLight::createShader(Renderer *renderer) const {
	return new AreaLightShader(renderer, m_radiance);
}

MTS_IMPLEMENT_CLASS(AreaLightShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(AreaLight, false, Emitter)
MTS_EXPORT_PLUGIN(AreaLight, "Area light");
MTS_NAMESPACE_END
