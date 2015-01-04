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
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

static StatsCounter uniformShootRatio("Unbiased photon mapping", "Percentage of uniform sampled shoot (phong)", EPercentage);
static StatsCounter thetaShootRatio("Unbiased photon mapping", "Percentage of theta bounded shoot (phong)", EPercentage);
static StatsCounter phiShootRatio("Unbiased photon mapping", "Percentage of theta-phi bounded shoot (phong)", EPercentage);

/*!\plugin{phong}{Modified Phong BRDF}
 * \order{14}
 * \parameters{
 *     \parameter{exponent}{\Float\Or\Texture}{
 *         Specifies the Phong exponent \default{30}.
 *     }
 *     \parameter{specular\showbreak Reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the weight of the specular reflectance component.\default{0.2}}
 *     \parameter{diffuse\showbreak Reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the weight of the diffuse reflectance component\default{0.5}}
 * }
 * \renderings{
 *     \rendering{Exponent$\,=60$}{bsdf_phong_60}
 *     \rendering{Exponent$\,=300$}{bsdf_phong_300}
 * }

 * This plugin implements the modified Phong reflectance model as described in
 * \cite{Phong1975Illumination} and \cite{Lafortune1994Using}. This heuristic
 * model is mainly included for historical reasons---its use in new scenes is
 * discouraged, since significantly more realistic models have been developed
 * since 1975.
 *
 * If possible, it is recommended to switch to a BRDF that is based on
 * microfacet theory and includes knowledge about the material's index of
 * refraction. In Mitsuba, two good alternatives to \pluginref{phong} are
 * the plugins \pluginref{roughconductor} and \pluginref{roughplastic}
 * (depending on the material type).
 *
 * When using this plugin, note that the diffuse and specular reflectance
 * components should add up to a value less than or equal to one (for each
 * color channel). Otherwise, they will automatically be scaled appropriately
 * to ensure energy conservation.
 */
class Phong : public BSDF {
public:
	Phong(const Properties &props)
		: BSDF(props) {
		m_diffuseReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(0.2f)));
		m_exponent = new ConstantFloatTexture(
			props.getFloat("exponent", 30.0f));
		m_specularSamplingWeight = 0.0f;
	}

	Phong(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_exponent = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide |
			((!m_specularReflectance->isConstant()
			  || !m_exponent->isConstant()) ? ESpatiallyVarying : 0));
		m_components.push_back(EDiffuseReflection | EFrontSide
			| (m_diffuseReflectance->isConstant() ? 0 : ESpatiallyVarying));

		/* Verify the input parameters and fix them if necessary */
		std::pair<Texture *, Texture *> result = ensureEnergyConservation(
			m_specularReflectance, m_diffuseReflectance,
			"specularReflectance", "diffuseReflectance", 1.0f);
		m_specularReflectance = result.first;
		m_diffuseReflectance = result.second;

		/* Compute weights that steer samples towards
		   the specular or diffuse components */
		Float dAvg = m_diffuseReflectance->getAverage().getLuminance(),
			  sAvg = m_specularReflectance->getAverage().getLuminance();
		m_specularSamplingWeight = sAvg / (dAvg + sAvg);

		//hack to force single component
// 		if (m_specularSamplingWeight < 0.5f){
// 			m_specularSamplingWeight = 0.f;
// 			m_specularReflectance = new ConstantSpectrumTexture(Spectrum(0.f));
// 		}
// 		else{
// 			m_specularSamplingWeight = 1.f;
// 			m_diffuseReflectance = new ConstantSpectrumTexture(Spectrum(0.f));
// 		}

		m_usesRayDifferentials =
			m_diffuseReflectance->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials() ||
			m_exponent->usesRayDifferentials();

		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_diffuseReflectance->eval(its);
	}

	Spectrum getSpecularReflectance(const Intersection &its) const {
		return m_specularReflectance->eval(its);
	}

	/// Reflection in local coordinates
	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
			return Spectrum(0.0f);

		bool hasSpecular = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		Spectrum result(0.0f);
		if (hasSpecular) {
			Float alpha    = dot(bRec.wo, reflect(bRec.wi)),
				  exponent = m_exponent->eval(bRec.its).average();

			if (alpha > 0.0f) {
				result += m_specularReflectance->eval(bRec.its) *
					((exponent + 2) * INV_TWOPI * std::pow(alpha, exponent));
			}
		}

		if (hasDiffuse)
			result += m_diffuseReflectance->eval(bRec.its) * INV_PI;

		return result * Frame::cosTheta(bRec.wo);
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
			return 0.0f;

		bool hasSpecular = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		Float diffuseProb = 0.0f, specProb = 0.0f;

		if (hasDiffuse)
			diffuseProb = Warp::squareToCosineHemispherePdf(bRec.wo);

		if (hasSpecular) {
			Float alpha    = dot(bRec.wo, reflect(bRec.wi)),
				  exponent = m_exponent->eval(bRec.its).average();
			if (alpha > 0)
				specProb = std::pow(alpha, exponent) *
					(exponent + 1.0f) / (2.0f * M_PI);
		}

		if (hasDiffuse && hasSpecular)
			return m_specularSamplingWeight * specProb +
				   (1-m_specularSamplingWeight) * diffuseProb;
		else if (hasDiffuse)
			return diffuseProb;
		else if (hasSpecular)
			return specProb;
		else
			return 0.0f;
	}

	inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasSpecular = (bRec.typeMask & EGlossyReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
				&& (bRec.component == -1 || bRec.component == 1);

		if (!hasSpecular && !hasDiffuse)
			return Spectrum(0.0f);

		bool choseSpecular = hasSpecular;

		if (hasDiffuse && hasSpecular) {
			if (sample.x <= m_specularSamplingWeight) {
				sample.x /= m_specularSamplingWeight;
			} else {
				sample.x = (sample.x - m_specularSamplingWeight)
					/ (1-m_specularSamplingWeight);
				choseSpecular = false;
			}
		}

		if (choseSpecular) {
			Vector R = reflect(bRec.wi);
			Float exponent = m_exponent->eval(bRec.its).average();

			/* Sample from a Phong lobe centered around (0, 0, 1) */
			Float sinAlpha = std::sqrt(1-std::pow(sample.y, 2/(exponent + 1)));
			Float cosAlpha = std::pow(sample.y, 1/(exponent + 1));
			Float phi = (2.0f * M_PI) * sample.x;
			Vector localDir = Vector(
				sinAlpha * std::cos(phi),
				sinAlpha * std::sin(phi),
				cosAlpha
			);

			/* Rotate into the correct coordinate system */
			bRec.wo = Frame(R).toWorld(localDir);
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyReflection;

			if (Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);
		} else {
			bRec.wo = Warp::squareToCosineHemisphere(sample);
			bRec.sampledComponent = 0;
			bRec.sampledType = EDiffuseReflection;
		}
		bRec.eta = 1.0f;

		_pdf = pdf(bRec, ESolidAngle);

		if (_pdf == 0)
			return Spectrum(0.0f);
		else
			return eval(bRec, ESolidAngle) / _pdf;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		Float pdf;
		return Phong::sample(bRec, pdf, sample);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "exponent")
				m_exponent = static_cast<Texture *>(child);
			else if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else if (name == "diffuseReflectance")
				m_diffuseReflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_diffuseReflectance.get());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_exponent.get());
	}

	Float getRoughness(const Intersection &its, int component) const {
		Assert(component == 0 || component == 1);
		/* Find the Beckmann-equivalent roughness */
		if (component == 0)
			return std::sqrt(2 / (2+m_exponent->eval(its).average()));
		else
			return std::numeric_limits<Float>::infinity();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "Phong[" << endl
   			<< "  id = \"" << getID() << "\"," << endl
			<< "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
			<< "  diffuseSamplingWeight = " << (1-m_specularSamplingWeight) << "," << endl
			<< "  exponent = " << indent(m_exponent->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	Float gatherAreaPdf(Vector wi, Vector wo, Float gatherRadius, 
		std::vector<Vector2> &componentCDFs, std::vector<Vector4> &componentBounds) const{
		if (Frame::cosTheta(wi) <= 0)
			return 0.f;		

		// initial spec component
		Vector4 bbox = Vector4(1.f, 0.f, 0.f, 1.f);		
		Vector4 bboxd = Vector4(0.f, 0.5f * M_PI, 0.f, 2.f * M_PI);		
		
		Float dis = wo.length();
		if (dis < gatherRadius){
			int numNode = 2;			
			componentCDFs.push_back(Vector2(1.f, *(float*)&numNode));		// level root node
			int ptrBound = -componentBounds.size();										// specular node
			componentCDFs.push_back(Vector2(m_specularSamplingWeight, *(float*)&ptrBound));
			componentBounds.push_back(bbox);
			ptrBound = -componentBounds.size();											// diffuse node
			componentCDFs.push_back(Vector2(1.f - m_specularSamplingWeight, *(float*)&ptrBound));
			componentBounds.push_back(bboxd);
			return 1.f;
		}
		
		Float dTheta = acos(sqrt(dis * dis - gatherRadius * gatherRadius) / dis);		

		// specular component
		Vector wir = reflect(wi);
		Vector dir = Frame(wir).toLocal(wo / dis);		
		Float exponent = m_exponent->getAverage().average();		
		Float theta = acos(dir.z);
		Float theta0 = theta - dTheta;
		Float theta1 = theta + dTheta;		
		theta0 = std::max(theta0, (Float)0.0);
		theta1 = std::min(theta1, (Float)(0.5 * M_PI));
		Float cos0 = std::max((Float)0.0, std::min((Float)1.0, pow(cos(theta0), exponent + 1.f)));
		Float cos1 = std::max((Float)0.0, std::min((Float)1.0, pow(cos(theta1), exponent + 1.f)));
		Vector dirProj = Vector(dir.x, dir.y, 0.f);
		Float disProj = dirProj.length() * dis;
		Float sqrDisTangent = disProj * disProj - gatherRadius * gatherRadius;
		Float probSpec = 1.f;	
		bbox.x = cos0; bbox.y = cos1;
		if (sqrDisTangent >= 0.f){
			// sample the bbox of the cone				
			Float cosdPhi = sqrt(sqrDisTangent) / disProj;
			Float dPhi = acos(cosdPhi);
			Float phi = atan2(dirProj.y, dirProj.x);
			Float inv2Pi = 0.5f / M_PI;
			bbox.z = (phi - dPhi) * inv2Pi;
			bbox.w = (phi + dPhi) * inv2Pi;
			probSpec = dPhi / M_PI * (cos0 - cos1);
			if (probSpec < 0.f){
				float fuck = 1.f;
			}
		}
		else{
			// sample the full sphere cap of polar			
			probSpec = 1.f - cos1;
		}

		// diffuse component
		dir = wo / dis;
		theta = acos(dir.z);
		theta0 = theta - dTheta;
		theta1 = theta + dTheta;
		theta0 = std::max(theta0, (Float)0.0);
		theta1 = std::min(theta1, (Float)(0.5 * M_PI));
		dirProj = Vector(dir.x, dir.y, 0.f);
		disProj = dirProj.length() * dis;
		sqrDisTangent = disProj * disProj - gatherRadius * gatherRadius;
		Float probDiff = 1.f;		
		bboxd.x = theta0; bboxd.y = theta1;
		if (sqrDisTangent >= 0.f){
			// sample the bbox of the cone			
			Float cosdPhi = sqrt(sqrDisTangent) / disProj;
			Float dPhi = acos(cosdPhi);
			Float phi = atan2(dirProj.y, dirProj.x);
			cos0 = cos(2.f * theta0);
			cos1 = cos(2.f * theta1);
			probDiff = 0.5f * dPhi * (cos0 - cos1) / M_PI;
			bboxd.z = phi - dPhi; bboxd.w = phi + dPhi;
		} else {
			// sample the full sphere cap of polar
			cos1 = cos(2.f * theta1);
			probDiff = 0.5f * (1.f - cos1);
		}
		Float prob = probSpec * m_specularSamplingWeight + probDiff * (1.f - m_specularSamplingWeight);

		int numNode = 2;
		componentCDFs.push_back(Vector2(prob, *(float*)&numNode));		// level root node
		int ptrBound = -componentBounds.size();										// specular node
		componentCDFs.push_back(Vector2(probSpec * m_specularSamplingWeight, *(float*)&ptrBound));
		componentBounds.push_back(bbox);
		ptrBound = -componentBounds.size();											// diffuse node
		componentCDFs.push_back(Vector2(probDiff * (1.f - m_specularSamplingWeight), *(float*)&ptrBound));
		componentBounds.push_back(bboxd);
		
		return prob;
	}
	
	Vector sampleGatherArea(Vector wi, Vector wo, Float gatherRadius, Point2 sample, 
		int ptrTree, std::vector<Vector2> componentCDFs, std::vector<Vector4> componentBounds) const{
		if (Frame::cosTheta(wi) <= 0)
			return Vector(0.f);

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
		if (chosenLobe == 0){
			/* Update statistics */			
			if (bbox.x == 1.f && bbox.y == 0.f && bbox.z == 0.f && bbox.w == 1.f)
				++uniformShootRatio;
			else if (bbox.z == 0.f && bbox.w == 1.f)
				++thetaShootRatio;
			else
				++phiShootRatio;
			/* Sample from a Phong lobe centered around (0, 0, 1) */
			Float exponent = m_exponent->getAverage().average();
			sample.y = sample.y * (bbox.x - bbox.y) + bbox.y;
			sample.x = sample.x * (bbox.w - bbox.z) + bbox.z;
			Vector R = reflect(wi);
			Float sinAlpha = std::sqrt(1 - std::pow(sample.y, 2 / (exponent + 1)));
			Float cosAlpha = std::pow(sample.y, 1 / (exponent + 1));
			Float phi = (2.0f * M_PI) * sample.x;
			Vector localDir = Vector(
				sinAlpha * std::cos(phi),
				sinAlpha * std::sin(phi),
				cosAlpha
				);
			/* Rotate into the correct coordinate system */
			dir = Frame(R).toWorld(localDir);
			if (Frame::cosTheta(dir) <= 0)
				return Vector(0.0f);
		}
		else{			
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
		}
		return dir;
	}

	Float getBandwidth() const{
		if (m_specularSamplingWeight == 0.f)
			return 0.f;
		else{
			return m_exponent->getAverage().average();;
		}
	}	

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_diffuseReflectance;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_exponent;
	Float m_specularSamplingWeight;
};

// ================ Hardware shader implementation ================

/**
 * The GLSL implementation clamps the exponent to 30 so that a
 * VPL renderer will able to handle the material reasonably well.
 */
class PhongShader : public Shader {
public:
	PhongShader(Renderer *renderer, const Texture *exponent,
			const Texture *diffuseColor, const Texture *specularColor)
		  : Shader(renderer, EBSDFShader),
			m_exponent(exponent),
			m_diffuseReflectance(diffuseColor),
			m_specularReflectance(specularColor) {
		m_exponentShader = renderer->registerShaderForResource(m_exponent.get());
		m_diffuseReflectanceShader = renderer->registerShaderForResource(m_diffuseReflectance.get());
		m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
	}

	bool isComplete() const {
		return m_exponentShader.get() != NULL &&
			   m_diffuseReflectanceShader.get() != NULL &&
			   m_specularReflectanceShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_exponentShader.get());
		deps.push_back(m_diffuseReflectanceShader.get());
		deps.push_back(m_specularReflectanceShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_exponent.get());
		renderer->unregisterShaderForResource(m_diffuseReflectance.get());
		renderer->unregisterShaderForResource(m_specularReflectance.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) <= 0.0 || cosTheta(wo) <= 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    vec3 R = vec3(-wi.x, -wi.y, wi.z);" << endl
			<< "    float specRef = 0.0, alpha = dot(R, wo);" << endl
			<< "    float exponent = min(30.0, " << depNames[0] << "(uv)[0]);" << endl
			<< "    if (alpha > 0.0)" << endl
			<< "    	specRef = pow(alpha, exponent) * " << endl
			<< "      (exponent + 2) * 0.15915;" << endl
			<< "    return (" << depNames[1] << "(uv) * inv_pi" << endl
			<< "           + " << depNames[2] << "(uv) * specRef) * cosTheta(wo);" << endl
			<< "}" << endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (wi.z <= 0.0 || wo.z <= 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[1] << "(uv) * (inv_pi * cosTheta(wo));" << endl
			<< "}" << endl;
	}


	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_exponent;
	ref<const Texture> m_diffuseReflectance;
	ref<const Texture> m_specularReflectance;
	ref<Shader> m_exponentShader;
	ref<Shader> m_diffuseReflectanceShader;
	ref<Shader> m_specularReflectanceShader;
};

Shader *Phong::createShader(Renderer *renderer) const {
	return new PhongShader(renderer, m_exponent.get(),
		m_diffuseReflectance.get(), m_specularReflectance.get());
}

MTS_IMPLEMENT_CLASS(PhongShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Phong, false, BSDF)
MTS_EXPORT_PLUGIN(Phong, "Modified Phong BRDF");
MTS_NAMESPACE_END
