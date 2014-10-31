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

#pragma once
#if !defined(__MITSUBA_BIDIR_PATHSAMPLER_H_)
#define __MITSUBA_BIDIR_PATHSAMPLER_H_

#include <mitsuba/bidir/path.h>
#include <boost/function.hpp>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/kdtree.h>

MTS_NAMESPACE_BEGIN

/*
*	Misc for VCM
*/
enum EMisTech {
	EVCM = 0,
	EVC = 1,
	EVM = 2,
	EMisTechs = 3
};
struct MTS_EXPORT_BIDIR MisState{
	float state[EMisTechs];
	MisState(){
		state[EVCM] = 0.f;
		state[EVC] = 0.f;
		state[EVM] = 0.f;
	}
	float& operator[](EMisTech tech){
		return state[tech];
	}
};
struct LightVertex{
	Spectrum importanceWeight;
	Vector wo;
	MisState emitterState;
	inline  LightVertex(const PathVertex *vs, const PathVertex *vsPred,
		const MisState &state, Spectrum _wgt) :
		emitterState(state), importanceWeight(_wgt){
		wo = normalize(vsPred->getPosition() - vs->getPosition());
	}
};
struct LightVertexExt{
	int depth;
	const Shape *shape;
	Point3 position;
	Frame shFrame;
	Frame geoFrame;
	uint8_t measure;
	uint8_t type;

	bool hasVsPred;
	Point3 posPred;
	Float pdfPred;
	uint8_t typePred;
	uint8_t measPred;

	inline LightVertexExt(const PathVertex *vs, const PathVertex *vsPred, int _depth){
		const Intersection &its = vs->getIntersection();
		depth = _depth;
		shape = its.shape;
		position = its.p;
		shFrame = its.shFrame;
		geoFrame = its.geoFrame;
		measure = vs->measure;
		type = vs->type;
		if (vsPred != NULL){
			hasVsPred = true;
			posPred = vsPred->getPosition();
			pdfPred = vsPred->pdf[EImportance];
			typePred = vsPred->type;
			measPred = vsPred->measure;
		}
		else{
			hasVsPred = false;
		}
	}
};
struct LightPathNodeData{
	int depth;
	size_t vertexIndex;
};
struct LightPathNode : public SimpleKDNode < Point, LightPathNodeData > {

	inline LightPathNode(){}

	inline  LightPathNode(const Point3 p, size_t vertexIndex, int depth){
		position = p;
		data.vertexIndex = vertexIndex;
		data.depth = depth;
	}
	inline LightPathNode(Stream *stream) {
		// TODO
	}
	void serialize(Stream *stream) const {
		// TODO
	}
};
typedef PointKDTree<LightPathNode>		LightPathTree;
typedef LightPathTree::IndexType     IndexType;
typedef LightPathTree::SearchResult SearchResult;


/*
*	Misc for CMLT
*/
enum EConnectionFlags {
	EConnectVisibility = 1,
	EConnectGeometry = 2,
	EConnectBRDF = 4,
	EConnectMis = 8,
	EConnectImportance = 16,
	EConnectRadiance = 32,
	EConnectAll = 64
};
struct MTS_EXPORT_BIDIR SplatListImp {
	/// Represents a screen-space splat produced by a path sampling technique
	typedef std::pair<Point2, Spectrum> Splat;

	/// A series of splats associated with the current sample
	std::vector<Splat> splats;
	std::vector<Spectrum> importances;
	/// Combined importance of all splats in this sample
	Float importance;
	/// Total number of samples in the splat list
	int nSamples;

	inline SplatListImp() : importance(1.0f), nSamples(0) { }

	/// for arbitrary Metropolis importance functions	
	/// Appends a splat entry to the list
	inline void append(const Point2 &samplePos, const Spectrum &value, const Spectrum &imp) {
		splats.push_back(std::make_pair(samplePos, value));
		importances.push_back(imp);
		importance += imp.getLuminance();
		++nSamples;
	}

	/// Increases the contribution of an existing splat
	inline void accum(size_t i, const Spectrum &value, const Spectrum &imp) {
		splats[i].second += value;
		importances[i] += imp;
		importance += imp.getLuminance();
		++nSamples;
	}

	/// Returns the number of contributions
	inline size_t size() const {
		return splats.size();
	}

	/// Clear the splat list
	inline void clear() {
		importance = 1;
		nSamples = 0;
		splats.clear();
		importances.clear();
	}

	/// Return the position associated with a splat in the list
	inline const Point2 &getPosition(size_t i) const { return splats[i].first; }

	/// Return the spectral contribution associated with a splat in the list
	inline const Spectrum &getValue(size_t i) const { return splats[i].second; }

	inline const Spectrum &getImportance(size_t i) const { return importances[i]; }

	/**
	* \brief Normalize the splat list
	*
	* This function divides all splats so that they have unit
	* luminance (though it leaves the \c luminance field untouched).
	* When given an optional importance map in 2-stage MLT approaches,
	* it divides the splat values by the associated importance
	* map values
	*/
	void normalize(const Bitmap *importanceMap = NULL){
		if (importanceMap) {
			BDAssert(false); // not implemented yet
		}
		if (importance > 0) {
			/* Normalize the contributions */
			Float invImportance = 1.0f / importance;
			for (size_t i = 0; i < splats.size(); ++i){
				splats[i].second *= invImportance;
				importances[i] *= invImportance;
			}
		}
	}

	/// Return a string representation
	std::string toString() const;	
};

/* ==================================================================== */
/*                         Work result for UPM                        */
/* ==================================================================== */
#define UPM_DEBUG 1
class UPMWorkResult : public WorkResult {
public:
	UPMWorkResult(const int width, const int height, const int maxDepth, const ReconstructionFilter *rfilter){
		/* Stores the 'camera image' -- this can be blocked when
		spreading out work to multiple workers */
		Vector2i blockSize = Vector2i(width, height);

		m_block = new ImageBlock(Bitmap::ESpectrum, blockSize, rfilter);
		m_block->setOffset(Point2i(0, 0));
		m_block->setSize(blockSize);

		/* When debug mode is active, we additionally create
		full-resolution bitmaps storing the contributions of
		each individual sampling strategy */
#if UPM_DEBUG == 1
		m_debugBlocks.resize(
			maxDepth*(5 + maxDepth) / 2);

		for (size_t i = 0; i<m_debugBlocks.size(); ++i) {
			m_debugBlocks[i] = new ImageBlock(
				Bitmap::ESpectrum, blockSize, rfilter);
			m_debugBlocks[i]->setOffset(Point2i(0, 0));
			m_debugBlocks[i]->setSize(blockSize);
		}
#endif
		sampleCount = 0;
	}

	// Clear the contents of the work result
	void clear(){
#if UPM_DEBUG == 1
		for (size_t i = 0; i < m_debugBlocks.size(); ++i)
			m_debugBlocks[i]->clear();
#endif
		m_block->clear();
		sampleCount = 0;
	}

	/// Fill the work result with content acquired from a binary data stream
	virtual void load(Stream *stream){
#if UPM_DEBUG == 1
		for (size_t i = 0; i < m_debugBlocks.size(); ++i)
			m_debugBlocks[i]->load(stream);
#endif
		m_block->load(stream);
	}

	/// Serialize a work result to a binary data stream
	virtual void save(Stream *stream) const{
#if UPM_DEBUG == 1
		for (size_t i = 0; i < m_debugBlocks.size(); ++i)
			m_debugBlocks[i]->save(stream);
#endif
		m_block->save(stream);
	}

	/// Aaccumulate another work result into this one
	void put(const UPMWorkResult *workResult){
#if UPM_DEBUG == 1
		for (size_t i = 0; i < m_debugBlocks.size(); ++i)
			m_debugBlocks[i]->put(workResult->m_debugBlocks[i].get());
#endif
		m_block->put(workResult->m_block.get());
		sampleCount += workResult->getSampleCount();
	}

#if UPM_DEBUG == 1
	/* In debug mode, this function allows to dump the contributions of
	the individual sampling strategies to a series of images */
	void dump(const int width, const int height, const int maxDepth,
		const fs::path &prefix, const fs::path &stem) const {
		Float weight = 1.f / (Float)sampleCount;
		for (int k = 1; k <= maxDepth; ++k) {
			for (int t = 0; t <= k + 1; ++t) {
				size_t s = k + 1 - t;
				Bitmap *bitmap = const_cast<Bitmap *>(m_debugBlocks[strategyIndex(s, t)]->getBitmap());
				ref<Bitmap> ldrBitmap = bitmap->convert(Bitmap::ERGB, Bitmap::EFloat, -1, weight);
				fs::path filename =
					prefix / fs::path(formatString("%s_upm_k%02i_s%02i_t%02i.pfm", stem.filename().string().c_str(), k, s, t));
				ref<FileStream> targetFile = new FileStream(filename,
					FileStream::ETruncReadWrite);
				ldrBitmap->write(Bitmap::EPFM, targetFile, 1);
			}
		}
	}

	inline void putDebugSample(int s, int t, const Point2 &sample,
		const Spectrum &spec) {
		m_debugBlocks[strategyIndex(s, t)]->put(sample, (const Float *)&spec);
	}
#endif

	inline void putSample(const Point2 &sample, const Float *value) {
		m_block->put(sample, value);
	}

	// 	inline void putLightSample(const Point2 &sample, const Spectrum &spec) {
	// 		m_lightImage->put(sample, spec, 1.0f);
	// 	}

	inline const ImageBlock *getImageBlock() const {
		return m_block.get();
	}
	inline ImageBlock *getImageBlock() {
		return m_block.get();
	}

	// 	inline const ImageBlock *getLightImage() const {
	// 		return m_lightImage.get();
	// 	}

	inline void setSize(const Vector2i &size) {
		m_block->setSize(size);
	}

	inline void setOffset(const Point2i &offset) {
		m_block->setOffset(offset);
	}

	void accumSampleCount(size_t count){
		sampleCount += count;
	}
	size_t getSampleCount() const{
		return sampleCount;
	}

	/// Return a string representation
	std::string toString() const {
		return m_block->toString();
	}

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~UPMWorkResult(){}

	inline int strategyIndex(int s, int t) const {
		int above = s + t - 2;
		return s + above*(5 + above) / 2;
	}
protected:
#if UPM_DEBUG == 1
	ref_vector<ImageBlock> m_debugBlocks;
#endif
	size_t sampleCount;
	ref<ImageBlock> m_block; // , m_lightImage;
};

/**
 * \brief Implements a sampling strategy that is able to produce paths using
 * bidirectional path tracing or unidirectional volumetric path tracing.
 *
 * This versatile class does the heavy lifting under the hood of Mitsuba's
 * PSSMLT implementation. It is also used to provide the Veach-MLT
 * implementation with a luminance estimate and seed paths.
 *
 * \author Wenzel Jakob
 * \ingroup libbidir
 */
class MTS_EXPORT_BIDIR PathSampler : public Object {
public:
	/**
	 * \brief Callback type for use with \ref samplePaths()
	 *
	 * The arguments are (s, t, weight, path) where \c s denotes the number
	 * of steps from the emitter, \c is the number of steps from the sensor,
	 * and \c weight contains the importance weight associated with the sample.
	 */
	typedef boost::function<void (int, int, Float, Path &)> PathCallback;

	/// Specifies the sampling algorithm that is internally used
	enum ETechnique {
		/// Bidirectional path tracing
		EBidirectional,

		/// Unidirectional path tracing (via the 'volpath' plugin)
		EUnidirectional
	};

	/**
	 * Construct a new path sampler
	 *
	 * \param technique
	 *     What path generation technique should be used (unidirectional
	 *     or bidirectional path tracing?)
	 *
	 * \param scene
	 *     \ref A pointer to the underlying scene
	 *
	 * \param emitterSampler
	 *     Sample generator that should be used for the random walk
	 *     from the emitter direction
	 *
	 * \param sensorSampler
	 *     Sample generator that should be used for the random walk
	 *     from the sensor direction
	 *
	 * \param directSampler
	 *     Sample generator that should be used for direct sampling
	 *     strategies (or \c NULL, when \c sampleDirect=\c false)
	 *
	 * \param maxDepth
	 *     Maximum path depth to be visualized (-1==infinite)
	 *
	 * \param rrDepth
	 *     Depth to begin using russian roulette
	 *
	 * \param excludeDirectIllum
	 *     If set to true, the direct illumination component will
	 *     be ignored. Note that this parameter is unrelated
	 *     to the next one (\a sampleDirect) although they are
	 *     named similarly.
	 *
	 * \param sampleDirect
	 *     When this parameter is set to true, specialized direct
	 *     sampling strategies are used for s=1 and t=1 paths.
	 *
	 * \param lightImage
	 *    Denotes whether or not rendering strategies that require a 'light image'
	 *    (specifically, those with <tt>t==0</tt> or <tt>t==1</tt>) are included
	 *    in the rendering process.
	 */
	PathSampler(ETechnique technique, const Scene *scene, Sampler *emitterSampler,
		Sampler *sensorSampler, Sampler *directSampler, int maxDepth, int rrDepth,
		bool excludeDirectIllum, bool sampleDirect, bool lightImage = true,
		Float gatherRadius = 0.0f, Sampler *lightPathSampler = NULL);

	/**
	 * \brief Generate a sample using the configured sampling strategy
	 *
	 * The result is stored as a series of screen-space splats (pixel position
	 * and spectral value pairs) within the parameter \c list. These can be
	 * used to implement algorithms like Bidirectional Path Tracing or Primary
	 * Sample Space MLT.
	 *
	 * \param offset
	 *    Specifies the desired integer pixel position of the sample. The special
	 *    value <tt>Point2i(-1)</tt> results in uniform sampling in screen space.
	 *
	 * \param list
	 *    Output parameter that will receive a list of splats
	 */
	void sampleSplats(const Point2i &offset, SplatList &list);

	/**
	 * \brief Sample a series of paths and invoke the specified callback
	 * function for each one.
	 *
	 * This function is similar to \ref sampleSplats(), but instead of
	 * returning only the contribution of the samples paths in the form of
	 * screen-space "splats", it returns the actual paths by invoking a
	 * specified callback function multiple times.
	 *
	 * This function is currently only implemented for the bidirectional
	 * sampling strategy -- i.e. it cannot be used with the unidirectional
	 * path tracer.
	 *
	 * \param offset
	 *    Specifies the desired integer pixel position of the sample. The special
	 *    value <tt>Point2i(-1)</tt> results in uniform sampling in screen space.
	 *
	 * \param pathCallback
	 *    A callback function that will be invoked once for each
	 *    path generated by the BDPT sampling strategy. The first argument
	 *    specifies the path importance weight.
	 */
	void samplePaths(const Point2i &offset, PathCallback &callback);

	/**
	 * \brief Generates a sequence of seeds that are suitable for
	 * starting a MLT Markov Chain
	 *
	 * This function additionally computes the average luminance
	 * over the image plane.
	 *
	 * \param sampleCount
	 *     The number of luminance samples that will be taken
	 * \param seedCount
	 *     The desired number of MLT seeds (must be > \c sampleCount)
	 * \param fineGrained
	 *     This parameter only matters when the technique is set to
	 *     \ref EBidirectional. It specifies whether to generate \ref PathSeed
	 *     records at the granularity of entire sensor/emitter subpaths or at
	 *     the granularity of their constituent sampling strategies.
	 * \param seeds
	 *     A vector of resulting MLT seeds
	 * \return The average luminance over the image plane
	 */
	Float generateSeeds(size_t sampleCount, size_t seedCount,
			bool fineGrained, const Bitmap *importanceMap,
			std::vector<PathSeed> &seeds);

	/**
	 * \brief Compute the average luminance over the image plane
	 * \param sampleCount
	 *     The number of luminance samples that will be taken
	 */
	Float computeAverageLuminance(size_t sampleCount);

	/**
	 * \brief Reconstruct a path from a \ref PathSeed record
	 *
	 * Given a \ref PathSeed data structure, this function rewinds
	 * the random number stream of the underlying \ref ReplayableSampler
	 * to the indicated position and recreates the associated path.
	 */
	void reconstructPath(const PathSeed &seed, 
		const Bitmap *importanceMap, Path &result);

	/// Return the underlying memory pool
	inline MemoryPool &getMemoryPool() { return m_pool; }

	/// for Connection MLT
	void sampleSplatsConnection(const Point2i &offset, SplatListImp &list, const int connectionFlag);
	Float generateSeedsConnection(size_t sampleCount, size_t seedCount,
		bool fineGrained, const Bitmap *importanceMap,
		std::vector<PathSeed> &seeds, const int connectionFlag);
	int getConnectionFlag(bool connectionImportance, bool connectionRadiance, bool connectionVisibility,
		bool connectionMIS, bool connectionBSDFs, bool connectionGeometry, bool connectionFull);

	/// for VCM
	void gatherLightPaths(const bool useVC, const bool useVM, const float gatherRadius, const int nsample, ImageBlock* lightImage = NULL);
	void sampleSplatsVCM(const bool useVC, const bool useVM, const float gatherRadius, const Point2i &offset, const size_t cameraPathIndex, SplatList &list);
	/// for UPM
	void sampleSplatsUPM(UPMWorkResult *wr, const float gatherRadius, const Point2i &offset, const size_t cameraPathIndex, SplatList &list);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~PathSampler();
protected:
	ETechnique m_technique;
	ref<const Scene> m_scene;
	ref<SamplingIntegrator> m_integrator;
	ref<Sampler> m_emitterSampler;
	ref<Sampler> m_sensorSampler;
	ref<Sampler> m_directSampler;
	int m_maxDepth;
	int m_rrDepth;
	bool m_excludeDirectIllum;
	bool m_sampleDirect;
	bool m_lightImage;
	int m_emitterDepth, m_sensorDepth;
	Path m_emitterSubpath, m_sensorSubpath;
	Path m_connectionSubpath, m_fullPath;
	MemoryPool m_pool;

	// VCM
	size_t m_lightPathNum;
	ref<Sampler> m_lightPathSampler;	
	LightPathTree m_lightPathTree;
	std::vector<LightVertex> m_lightVertices;
	std::vector<LightVertexExt> m_lightVerticesExt;
	std::vector<size_t> m_lightPathEnds;
};

/**
 * \brief Stores information required to re-create a seed path (e.g. for MLT)
 *
 * This class makes it possible to transmit a path over the network or store
 * it locally, while requiring very little storage to do so. This is done by
 * describing a path using an index into a random number stream, which allows
 * to generate it cheaply when needed.
 */
struct PathSeed {
	size_t sampleIndex; ///< Index into a rewindable random number stream
	Float luminance;    ///< Luminance value of the path (for sanity checks)
	int s;              ///< Number of steps from the luminaire
	int t;              ///< Number of steps from the eye

	inline PathSeed() { }

	inline PathSeed(size_t sampleIndex, Float luminance, int s = 0, int t = 0)
		: sampleIndex(sampleIndex), luminance(luminance), s(s), t(t) { }

	inline PathSeed(Stream *stream) {
		sampleIndex = stream->readSize();
		luminance = stream->readFloat();
		s = stream->readInt();
		t = stream->readInt();
	}

	void serialize(Stream *stream) const {
		stream->writeSize(sampleIndex);
		stream->writeFloat(luminance);
		stream->writeInt(s);
		stream->writeInt(t);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "PathSeed[" << endl
			<< "  sampleIndex = " << sampleIndex << "," << endl
			<< "  luminance = " << luminance << "," << endl
			<< "  s = " << s << "," << endl
			<< "  t = " << t << endl
			<< "]";
		return oss.str();
	}
};

/**
 * MLT work unit -- wraps a \ref PathSeed into a
 * \ref WorkUnit instance.
 */
class SeedWorkUnit : public WorkUnit {
public:
	inline void set(const WorkUnit *wu) {
		m_id = static_cast<const SeedWorkUnit *>(wu)->m_id;
		m_seed = static_cast<const SeedWorkUnit *>(wu)->m_seed;
		m_timeout = static_cast<const SeedWorkUnit *>(wu)->m_timeout;
	}	

	inline const PathSeed &getSeed() const {
		return m_seed;
	}

	inline void setSeed(const PathSeed &seed) {
		m_seed = seed;
	}

	inline int getTimeout() const {
		return m_timeout;
	}

	inline void setTimeout(int timeout) {
		m_timeout = timeout;
	}

	inline void load(Stream *stream) {
		m_seed = PathSeed(stream);
		m_timeout = stream->readInt();
	}

	inline void save(Stream *stream) const {
		m_seed.serialize(stream);
		stream->writeInt(m_timeout);
	}

	inline std::string toString() const {
		return "SeedWorkUnit[]";
	}

	inline const int getID() const{
		return m_id;
	}
	inline const void setID(int id){
		m_id = id;
	}

	MTS_DECLARE_CLASS()
private:
	int m_id;
	PathSeed m_seed;
	int m_timeout;
};

/**
 * \brief List storage for the image-space contributions ("splats") of a
 * path sample generated using \ref PathSampler.
 */
struct MTS_EXPORT_BIDIR SplatList {
	/// Represents a screen-space splat produced by a path sampling technique
	typedef std::pair<Point2, Spectrum> Splat;

	/// A series of splats associated with the current sample
	std::vector<Splat> splats;
	/// Combined luminance of all splats in this sample
	Float luminance;
	/// Total number of samples in the splat list
	int nSamples;

	inline SplatList() : luminance(0.0f), nSamples(0) { }

	/// Appends a splat entry to the list
	inline void append(const Point2 &samplePos, const Spectrum &value) {
		splats.push_back(std::make_pair(samplePos, value));
		luminance += value.getLuminance();
		++nSamples;
	}

	/// Increases the contribution of an existing splat
	inline void accum(size_t i, const Spectrum &value) {
		splats[i].second += value;
		luminance += value.getLuminance();
		++nSamples;
	}

	/// Returns the number of contributions
	inline size_t size() const {
		return splats.size();
	}

	/// Clear the splat list
	inline void clear() {
		luminance = 0;
		nSamples = 0;
		splats.clear();
	}

	/// Return the position associated with a splat in the list
	inline const Point2 &getPosition(size_t i) const { return splats[i].first; }

	/// Return the spectral contribution associated with a splat in the list
	inline const Spectrum &getValue(size_t i) const { return splats[i].second; }

	/**
	 * \brief Normalize the splat list
	 *
	 * This function divides all splats so that they have unit
	 * luminance (though it leaves the \c luminance field untouched).
	 * When given an optional importance map in 2-stage MLT approaches,
	 * it divides the splat values by the associated importance
	 * map values
	 */
	void normalize(const Bitmap *importanceMap = NULL);

	/// Return a string representation
	std::string toString() const;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_BIDIR_PATHSAMPLER_H_ */
