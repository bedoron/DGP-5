//=============================================================================
//
//  CLASS HarmonicMapViewer
//
//=============================================================================


#ifndef HARMONIC_MAP_VIEWER_HH
#define HARMONIC_MAP_VIEWER_HH


//== INCLUDES =================================================================

#include <gmm.h>
#include "MeshViewer.hh"



//== CLASS DEFINITION =========================================================

	      

class HarmonicMapViewer : public MeshViewer
{
public:
  enum ParameterizationMode {NoParameterization, Uniform, Harmonic, Mean}; 

  typedef gmm::dense_matrix<double>  gmmMatrix;
	typedef std::vector<double>        gmmVector;

   
  /// default constructor
  HarmonicMapViewer(const char* _title, int _width, int _height, int iRepeat);

  // destructor
  ~HarmonicMapViewer();

  /// open mesh
  virtual bool open_mesh(const char* _meshfilename);


private:
  void calc_uniform_parameterization();
  void calc_harmonic_parameterization();
  void calc_mean_parameterization();

  void ComputeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iReapeats, ParameterizationMode imode);

  void calc_distortion(ParameterizationMode imode);

  // solve linear system _A * _x = _b
	void solve_linear_system( gmmMatrix& _A, 
														gmmVector& _b, 
														gmmVector& _x );

  virtual void init();
  virtual void draw(const std::string& _draw_mode);

  virtual void keyboard(int key, int x, int y);

private:
  int _Repeat;

  OpenMesh::VPropHandleT<OpenMesh::Vec2f>   vparam_u_, vparam_h_, texcoord_u_, texcoord_h_  ;
	OpenMesh::VPropHandleT<OpenMesh::Vec2f> vparam_m_, texcoord_m_;
  
  OpenMesh::VPropHandleT<Mesh::Point>       vpos_;

  GLuint  textureID_;
  float *_TextureCoordinates_u, *_TextureCoordinates_h, *_TextureCoordinates_m; // Not sure 

  Mesh::Point _bbMin3D, _bbMax3D, _bbMin2D, _bbMax2D;

  ParameterizationMode _ParameterizationMode_;

	float *TextureCoordinates() {
		switch(_ParameterizationMode_) {
		case Uniform: return _TextureCoordinates_u;
		case Harmonic: return _TextureCoordinates_h;
		case Mean: return _TextureCoordinates_m;
		default: throw std::string("Moo");
		}
	}

	OpenMesh::VPropHandleT<OpenMesh::Vec2f> &texcoord() {
		switch(_ParameterizationMode_) {
		case Uniform: return texcoord_u_;
		case Harmonic: return texcoord_h_;
		case Mean: return texcoord_m_;
		default: throw std::string("Moo");
		}
	}
	OpenMesh::VPropHandleT<OpenMesh::Vec2f> &vparam() {
		switch(_ParameterizationMode_) {
		case Uniform: return vparam_u_;
		case Harmonic: return vparam_h_;
		case Mean: return vparam_m_;
		default: throw std::string("Moo");
		}
	}

	inline Vec2f &texcoord(const Mesh::VertexHandle& vh) {
		return mesh_.property(texcoord(), vh);
	}

	inline Vec2f &vparam(const Mesh::VertexHandle& vh) {
		return mesh_.property(vparam(), vh);
	}

	
};


//=============================================================================
#endif // HARMONIC_MAP_VIEWER_HH defined
//=============================================================================

