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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

static StatsCounter uniformShootRatio("Unbiased photon mapping", "Percentage of uniform sampled shoot (diffuse)", EPercentage);
static StatsCounter thetaShootRatio("Unbiased photon mapping", "Percentage of theta bounded shoot (diffuse)", EPercentage);
static StatsCounter phiShootRatio("Unbiased photon mapping", "Percentage of theta-phi bounded shoot (diffuse)", EPercentage);

/*!\plugin{diffuse}{Smooth diffuse material}
 * \order{1}
 * \icon{bsdf_diffuse}
 * \parameters{
 *     \parameter{reflectance}{\Spectrum\Or\Texture}{
 *       Specifies the diffuse albedo of the
 *       material \default{0.5}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Homogeneous reflectance, see \lstref{diffuse-uniform}}
 *         {bsdf_diffuse_plain}
 *     \rendering{Textured reflectance, see \lstref{diffuse-textured}}
 *         {bsdf_diffuse_textured}
 * }
 *
 * The smooth diffuse material (also referred to as ``Lambertian'')
 * represents an ideally diffuse material with a user-specified amount of
 * reflectance. Any received illumination is scattered so that the surface
 * looks the same independently of the direction of observation.
 *
 * Apart from a  homogeneous reflectance value, the plugin can also accept
 * a nested or referenced texture map to be used as the source of reflectance
 * information, which is then mapped onto the shape based on its UV
 * parameterization. When no parameters are specified, the model uses the default
 * of 50% reflectance.
 *
 * Note that this material is one-sided---that is, observed from the
 * back side, it will be completely black. If this is undesirable,
 * consider using the \pluginref{twosided} BRDF adapter plugin.
 * \vspace{4mm}
 *
 * \begin{xml}[caption={A diffuse material, whose reflectance is specified
 *     as an sRGB color}, label=lst:diffuse-uniform]
 * <bsdf type="diffuse">
 *     <srgb name="reflectance" value="#6d7185"/>
 * </bsdf>
 * \end{xml}
 *
 * \begin{xml}[caption=A diffuse material with a texture map,
 *     label=lst:diffuse-textured]
 * <bsdf type="diffuse">
 *     <texture type="bitmap" name="reflectance">
 *         <string name="filename" value="wood.jpg"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */
class SmoothDiffuse : public BSDF {
public:
	SmoothDiffuse(const Properties &props)
		: BSDF(props) {
		/* For better compatibility with other models, support both
		   'reflectance' and 'diffuseReflectance' as parameter names */
		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance"
				: "diffuseReflectance", Spectrum(.5f)));
	}

	SmoothDiffuse(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	void configure() {
		/* Verify the input parameter and fix them if necessary */
		m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

		m_components.clear();
		if (m_reflectance->getMaximum().max() > 0)
			m_components.push_back(EDiffuseReflection | EFrontSide
				| (m_reflectance->isConstant() ? 0 : ESpatiallyVarying));
			m_usesRayDifferentials = m_reflectance->usesRayDifferentials();

		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_reflectance->eval(its);
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		return m_reflectance->eval(bRec.its)
			* (INV_PI * Frame::cosTheta(bRec.wo));
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		return Warp::squareToCosineHemispherePdf(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.wo = Warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		return m_reflectance->eval(bRec.its);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.wo = Warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = Warp::squareToCosineHemispherePdf(bRec.wo);
		return m_reflectance->eval(bRec.its);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
				&& (name == "reflectance" || name == "diffuseReflectance")) {
			m_reflectance = static_cast<Texture *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_reflectance.get());
	}

	Float getRoughness(const Intersection &its, int component) const {
		return std::numeric_limits<Float>::infinity();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothDiffuse[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  reflectance = " << indent(m_reflectance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	Float gatherAreaPdf(Vector wi, Vector wo, Float gatherRadius, 
		std::vector<Vector2> &componentCDFs, std::vector<Vector4> &componentBounds) const{
		if (Frame::cosTheta(wi) <= 0) return 0.f;
		Vector4 bbox = Vector4(0.f, 0.5f * M_PI, 0.f, 2.f * M_PI);		
		Vector dir = wo;
		Float dis = dir.length();
		if (dis < gatherRadius){
			int numNode = 1;
			int ptrBound = -componentBounds.size();
			componentCDFs.push_back(Vector2(1.f, *(float*)&numNode));		// level root node
			componentCDFs.push_back(Vector2(1.f, *(float*)&ptrBound));		// single CDF node
			componentBounds.push_back(bbox);
			return 1.f;
		}
		dir /= dis;
		Float dTheta = acos(sqrt(dis * dis - gatherRadius * gatherRadius) / dis);
		Float theta = acos(dir.z);
		Float theta0 = theta - dTheta;
		Float theta1 = std::min(0.5f * M_PI, theta + dTheta);		
		Float prob = 0.f;
		if (theta0 < 0.f){
			// sample the full sphere cap of polar
			Float cos1 = cos(2.f * theta1);
			prob = 0.5f * (1.f - cos1);
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
		int numNode = 1;
		int ptrBound = -componentBounds.size();
		componentCDFs.push_back(Vector2(prob, *(float*)&numNode));		// level root node
		componentCDFs.push_back(Vector2(prob, *(float*)&ptrBound));			// single CDF node
		componentBounds.push_back(bbox);
		return prob;		
	}
	
	Vector sampleGatherArea(Vector wi, Vector wo, Float gatherRadius, Point2 sample, 
		int ptrTree, std::vector<Vector2> componentCDFs, std::vector<Vector4> componentBounds) const{
		if (Frame::cosTheta(wi) <= 0) return Vector(0.f);		
		uniformShootRatio.incrementBase();
		thetaShootRatio.incrementBase();
		phiShootRatio.incrementBase();

		// sample CDF tree
		Vector2 rootnode = componentCDFs[ptrTree];
		Float invTotalPdf = 1.f / rootnode.x;
		int numNode = *(int*)&rootnode.y;
		Float cdfi = 0.f;
		int chosenLobe = -1;
		Vector4 bbox;
		for (int i = 0; i < numNode; i++){
			Vector2 nodei = componentCDFs[ptrTree + 1 + i];
			Float pdfi = nodei.x;
			if (sample.x <= (cdfi + pdfi) * invTotalPdf){
				// choose this component
				chosenLobe = i;
				int ptrBound = -*(int*)&nodei.y;
				bbox = componentBounds[ptrBound];
				sample.x = (sample.x - cdfi * invTotalPdf) / (pdfi * invTotalPdf);
				break;
			}
			cdfi += pdfi;
		}

		Vector dir;
		if (bbox.x == 0.f && bbox.y == 0.5f * M_PI && bbox.z == 0.f && bbox.w == 2.f * M_PI){
			// uniform sampling
			dir = Warp::squareToCosineHemisphere(sample);
			++uniformShootRatio;
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
			++thetaShootRatio;
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
			++phiShootRatio;
		}
		return dir;
	}	

	Float getBandwidth() const{
		return 0.f;
	}


	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_reflectance;
};

// ================ Hardware shader implementation ================

class SmoothDiffuseShader : public Shader {
public:
	SmoothDiffuseShader(Renderer *renderer, const Texture *reflectance)
		: Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
		m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_reflectance.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_reflectanceShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_reflectance;
	ref<Shader> m_reflectanceShader;
};

Shader *SmoothDiffuse::createShader(Renderer *renderer) const {
	return new SmoothDiffuseShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(SmoothDiffuseShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothDiffuse, "Smooth diffuse BRDF")
MTS_NAMESPACE_END
