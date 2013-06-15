//=============================================================================
//
//  CLASS HarmonicMapViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "HarmonicMapViewer.hh"
#include <windows.h>
#include <vector>

using std::vector;

//== IMPLEMENTATION ========================================================== 


HarmonicMapViewer::
	HarmonicMapViewer(const char* _title, int _width, int _height, int iRepeat)
	: MeshViewer(_title, _width, _height)
{ 
	_Repeat = iRepeat;

	mesh_.request_vertex_colors();

	mesh_.add_property(vpos_);
	mesh_.add_property(vparam_u_);
	mesh_.add_property(vparam_h_);

	mesh_.add_property(texcoord_u_);
	mesh_.add_property(texcoord_h_);

	mesh_.add_property(texcoord_m_);
	mesh_.add_property(vparam_m_);

	add_draw_mode("UV Domain");
	add_draw_mode("Textured mesh");

	_TextureCoordinates_u=NULL;
	_TextureCoordinates_h=NULL;
	_TextureCoordinates_m=NULL;
	_ParameterizationMode_ = NoParameterization;

	init();
}

HarmonicMapViewer::
	~HarmonicMapViewer()
{ 
	if (glIsTexture(textureID_))  
		glDeleteTextures( 1, &textureID_);

	if (_TextureCoordinates_u)
		delete [] _TextureCoordinates_u;

	if (_TextureCoordinates_h)
		delete [] _TextureCoordinates_h;

	if (_TextureCoordinates_m)
		delete [] _TextureCoordinates_m;
}

void
	HarmonicMapViewer::
	init()
{
	// base class first
	MeshViewer::init();

	// generate checkerboard-like image
	GLubyte tex[256][256][3];
	int index=0;
	for (int x=0; x<256; ++x)
	{
		for (int y=0; y<256; ++y)
		{
			if ((x<128&&y<128) ||(x>128&&y>128))
			{
				tex[x][y][0] = 0;
				tex[x][y][1] = 255;
				tex[x][y][2] = 0;
			}
			else
			{
				tex[x][y][0] = 255;
				tex[x][y][1] = 255;
				tex[x][y][2] = 255;
			}
		}
	}
	// generate texture
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 
	if (!glIsTexture(textureID_))
		glGenTextures(1, &textureID_);
	glBindTexture(GL_TEXTURE_2D, textureID_);

	// copy texture to GL
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, 256, 256,
		0, GL_RGB, GL_UNSIGNED_BYTE, tex);
}

bool
	HarmonicMapViewer::
	open_mesh(const char* _meshfilename)
{
	// load mesh
	if (MeshViewer::open_mesh(_meshfilename))
	{
		// store vertex initial positions and 3D mesh bounding box
		Mesh::VertexIter v_it=mesh_.vertices_begin(), v_end(mesh_.vertices_end());
		_bbMin3D = _bbMax3D = mesh_.point(v_it);
		for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
		{
			mesh_.property(vpos_,v_it) = mesh_.point(v_it);
			_bbMin3D.minimize(mesh_.point(v_it));
			_bbMax3D.maximize(mesh_.point(v_it));
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------


void
	HarmonicMapViewer::solve_linear_system( gmmMatrix& _M, gmmVector& _b, gmmVector& _x)
{
	unsigned int N = _b.size();
	_x.resize(N);
	std::vector< size_t >  ipvt(N);
	gmm::lu_factor( _M, ipvt );
	gmm::lu_solve( _M, ipvt, _x, _b );
} 


//-----------------------------------------------------------------------------

void HarmonicMapViewer::
	calc_uniform_parameterization()
{
	// ------------- IMPLEMENT HERE ---------
	// TASK 5.1 Uniform map computation:
	// Search and order boundary vertices
	// Compute boundary parameters
	// Solve the linear system for internal vertex parameters using solve_linear_system()
	// Store the parameters in the vparam_u_ vertex property
	// ------------- IMPLEMENT HERE ---------
	int nSize = mesh_.n_vertices();
	gmmMatrix W(nSize,nSize);
	gmm::clear(W);

	for(Mesh::HalfedgeIter heIter = mesh_.halfedges_begin(); heIter != mesh_.halfedges_end(); ++heIter) {
		Mesh::VertexHandle from = mesh_.from_vertex_handle(heIter);
		Mesh::VertexHandle to = mesh_.to_vertex_handle(heIter);

		int fromIdx = from.idx();
		int toIdx = to.idx();

		if(mesh_.is_boundary(from)) continue;

		if(!mesh_.is_boundary(from)) {
			W(fromIdx, toIdx) = 1;
		}
	}

	for(int diagonal=0; diagonal < nSize; ++diagonal) {
		Mesh::VertexHandle vertex = mesh_.vertex_handle(diagonal);
		if(mesh_.is_boundary(vertex)) {
			W(diagonal, diagonal)  = 1;
		} else {
			int sumLine = 0;
			for(int col=0; col < nSize; ++col) {
				if(col == diagonal) continue; 
				sumLine += W(diagonal, col);
			}
			W(diagonal,diagonal) = -sumLine;
		} 
	}

	// Find first Boundary
	Mesh::HalfedgeHandle firstBoundary, currentBoundary;
	for(Mesh::HalfedgeIter heIter = mesh_.halfedges_begin(); heIter != mesh_.halfedges_end(); ++heIter) {
		if(mesh_.is_boundary(heIter)) {
			firstBoundary = heIter;
			currentBoundary = mesh_.next_halfedge_handle(heIter.handle());
			break;
		}
	}

	gmmVector Bx(nSize);
	gmmVector By(nSize);

	gmmVector x(nSize);
	gmmVector y(nSize);

	gmm::clear(Bx); gmm::clear(By);
	gmm::clear(x); gmm::clear(y);

	float arcLength = mesh_.calc_edge_length(firstBoundary);
	int firstIdx = mesh_.from_vertex_handle(firstBoundary).idx();
	Bx[firstIdx] = 0;

	while(currentBoundary!=firstBoundary) {
		double currentLength = mesh_.calc_edge_length(currentBoundary);
		int currentIdx = mesh_.from_vertex_handle(currentBoundary).idx();
		Bx[currentIdx] = arcLength;

		arcLength += currentLength;
		currentBoundary = mesh_.next_halfedge_handle(currentBoundary);
	}

	for(int i=0; i < nSize; ++i) {
		Mesh::VertexHandle vertex = mesh_.vertex_handle(i);
		if(mesh_.is_boundary(vertex)) {
			double theta = (Bx[i]/arcLength)*2*M_PI;
			Bx[i] = cos(theta);
			By[i] = sin(theta);
		}
	}

	gmmMatrix W1(W); // I hate you gmm
	solve_linear_system(W, Bx, x);
	solve_linear_system(W1, By, y);

	for(Mesh::VertexIter viter(mesh_.vertices_begin()); viter != mesh_.vertices_end(); viter++) {
		int idx = viter.handle().idx();
		Vec2f paramCoord(x[idx], y[idx]);
		mesh_.property(vparam_u_, viter) = paramCoord;
	}

}

//-----------------------------------------------------------------------------

void HarmonicMapViewer::
	calc_harmonic_parameterization()
{
	// ------------- IMPLEMENT HERE ---------
	// TASK 5.2 harmonic map computation:
	// Search and order boundary vertices
	// Compute boundary parameters
	// Solve the linear system for internal vertex parameters using solve_linear_system()
	// Store the parameters in the vparam_h_ vertex property
	// ------------- IMPLEMENT HERE ---------
	int nSize = mesh_.n_vertices();
	gmmMatrix W(nSize,nSize);
	gmm::clear(W);

	for(Mesh::HalfedgeIter heIter = mesh_.halfedges_begin(); heIter != mesh_.halfedges_end(); ++heIter) {
		Mesh::VertexHandle from = mesh_.from_vertex_handle(heIter);
		Mesh::VertexHandle to = mesh_.to_vertex_handle(heIter);

		int fromIdx = from.idx();
		int toIdx = to.idx();

		if(mesh_.is_boundary(from)) 
			continue;

		float weight;

		Mesh::HalfedgeHandle opposite = mesh_.opposite_halfedge_handle(heIter);

		Mesh::VertexHandle va = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heIter));
		Mesh::VertexHandle vb = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(opposite));

		Mesh::Point p1 = mesh_.point(from);
		Mesh::Point p2 = mesh_.point(to);
		Mesh::Point p3 = mesh_.point(va);
		Mesh::Point p4 = mesh_.point(vb);

		Vec3f va1 = (p1-p3).normalize();
		Vec3f vb1 = (p1-p4).normalize();
		Vec3f va2 = (p2-p3).normalize();
		Vec3f vb2 = (p2-p4).normalize();

		float cosa = max(-0.99, min(0.99, va1 |va2));
		float cosb = max(-0.99, min(0.99, vb1 |vb2));

		float cota = (cosa) / sqrtf(1 - cosa*cosa);
		float cotb = (cosb) / sqrtf(1 - cosb*cosb);

		weight = (cota + cotb) / 2.0;
		if(weight<0)
			weight = 0;
		
		W(fromIdx, toIdx) = weight;
		
	}

	for(int diagonal=0; diagonal < nSize; ++diagonal) {
		Mesh::VertexHandle vertex = mesh_.vertex_handle(diagonal);
		if(mesh_.is_boundary(vertex)) {
			W(diagonal, diagonal)  = 1;
		} else {
			float sumLine = 0;
			for(int col=0; col < nSize; ++col) {
				if(col == diagonal) continue; 
				sumLine += W(diagonal, col);
			}
			W(diagonal,diagonal) = -sumLine;
		} 
	}

	// Find first Boundary
	Mesh::HalfedgeHandle firstBoundary, currentBoundary;
	for(Mesh::HalfedgeIter heIter = mesh_.halfedges_begin(); heIter != mesh_.halfedges_end(); ++heIter) {
		if(mesh_.is_boundary(heIter)) {
			firstBoundary = heIter;
			currentBoundary = mesh_.next_halfedge_handle(heIter.handle());
			break;
		}
	}

	gmmVector Bx(nSize);
	gmmVector By(nSize);

	gmmVector x(nSize);
	gmmVector y(nSize);

	gmm::clear(Bx); gmm::clear(By);
	gmm::clear(x); gmm::clear(y);

	float arcLength = mesh_.calc_edge_length(firstBoundary);
	int firstIdx = mesh_.from_vertex_handle(firstBoundary).idx();
	Bx[firstIdx] = 0;

	while(currentBoundary!=firstBoundary) {
		double currentLength = mesh_.calc_edge_length(currentBoundary);
		int currentIdx = mesh_.from_vertex_handle(currentBoundary).idx();
		Bx[currentIdx] = arcLength;

		arcLength += currentLength;
		currentBoundary = mesh_.next_halfedge_handle(currentBoundary);
	}

	for(int i=0; i < nSize; ++i) {
		Mesh::VertexHandle vertex = mesh_.vertex_handle(i);
		if(mesh_.is_boundary(vertex)) {
			double theta = (Bx[i]/arcLength)*2*M_PI;
			Bx[i] = cos(theta);
			By[i] = sin(theta);
		}
	}

	gmmMatrix W1(W); // I hate you gmm
	solve_linear_system(W, Bx, x);
	solve_linear_system(W1, By, y);

	for(Mesh::VertexIter viter(mesh_.vertices_begin()); viter != mesh_.vertices_end(); viter++) {
		int idx = viter.handle().idx();
		Vec2f paramCoord(x[idx], y[idx]);
		mesh_.property(vparam_h_, viter) = paramCoord;
	}
}

//-----------------------------------------------------------------------------

void HarmonicMapViewer::
	calc_mean_parameterization()
{
	int nSize = mesh_.n_vertices();
	gmmMatrix W(nSize,nSize);
	gmm::clear(W);

	for(Mesh::HalfedgeIter heIter = mesh_.halfedges_begin(); heIter != mesh_.halfedges_end(); ++heIter) {
		Mesh::VertexHandle from = mesh_.from_vertex_handle(heIter);
		Mesh::VertexHandle to = mesh_.to_vertex_handle(heIter);

		int fromIdx = from.idx();
		int toIdx = to.idx();

		if(mesh_.is_boundary(from)) 
			continue;

		float weight;

		Mesh::HalfedgeHandle opposite = mesh_.opposite_halfedge_handle(heIter);

		Mesh::VertexHandle va = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heIter));
		Mesh::VertexHandle vb = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(opposite));

		Mesh::Point p1 = mesh_.point(from);
		Mesh::Point p2 = mesh_.point(to);
		Mesh::Point p3 = mesh_.point(va);
		Mesh::Point p4 = mesh_.point(vb);

		/*Vec3f va1 = (p1-p3).normalize();
		Vec3f vb1 = (p1-p4).normalize();
		Vec3f va2 = (p2-p3).normalize();
		Vec3f vb2 = (p2-p4).normalize();*/
		Vec3f m = (p2-p1).normalize();
		Vec3f a = (p4-p1).normalize();
		Vec3f b = (p3-p1).normalize();

		float cosa = max(-0.99, min(0.99, m |a));
		float cosb = max(-0.99, min(0.99, m |b));

		float gamma = acos(cosa);
		float delta = acos(cosb);

		float htgg = (tan(gamma/2.));
		float htgd = (tan(delta/2.));

		weight = (htgg + htgd)/(2.0*(p2-p1).length());

		W(fromIdx, toIdx) = weight;

	}

	for(int diagonal=0; diagonal < nSize; ++diagonal) {
		Mesh::VertexHandle vertex = mesh_.vertex_handle(diagonal);
		if(mesh_.is_boundary(vertex)) {
			W(diagonal, diagonal)  = 1;
		} else {
			float sumLine = 0;
			for(int col=0; col < nSize; ++col) {
				if(col == diagonal) continue; 
				sumLine += W(diagonal, col);
			}
			W(diagonal,diagonal) = -sumLine;
		} 
	}

	// Find first Boundary
	Mesh::HalfedgeHandle firstBoundary, currentBoundary;
	for(Mesh::HalfedgeIter heIter = mesh_.halfedges_begin(); heIter != mesh_.halfedges_end(); ++heIter) {
		if(mesh_.is_boundary(heIter)) {
			firstBoundary = heIter;
			currentBoundary = mesh_.next_halfedge_handle(heIter.handle());
			break;
		}
	}

	gmmVector Bx(nSize);
	gmmVector By(nSize);

	gmmVector x(nSize);
	gmmVector y(nSize);

	gmm::clear(Bx); gmm::clear(By);
	gmm::clear(x); gmm::clear(y);

	float arcLength = mesh_.calc_edge_length(firstBoundary);
	int firstIdx = mesh_.from_vertex_handle(firstBoundary).idx();
	Bx[firstIdx] = 0;

	while(currentBoundary!=firstBoundary) {
		double currentLength = mesh_.calc_edge_length(currentBoundary);
		int currentIdx = mesh_.from_vertex_handle(currentBoundary).idx();
		Bx[currentIdx] = arcLength;

		arcLength += currentLength;
		currentBoundary = mesh_.next_halfedge_handle(currentBoundary);
	}

	for(int i=0; i < nSize; ++i) {
		Mesh::VertexHandle vertex = mesh_.vertex_handle(i);
		if(mesh_.is_boundary(vertex)) {
			double theta = (Bx[i]/arcLength)*2*M_PI;
			Bx[i] = cos(theta);
			By[i] = sin(theta);
		}
	}

	gmmMatrix W1(W); // I hate you gmm
	solve_linear_system(W, Bx, x);
	solve_linear_system(W1, By, y);

	for(Mesh::VertexIter viter(mesh_.vertices_begin()); viter != mesh_.vertices_end(); viter++) {
		int idx = viter.handle().idx();
		Vec2f paramCoord(x[idx], y[idx]);
		mesh_.property(vparam_m_, viter) = paramCoord;
	}
}

//-----------------------------------------------------------------------------


void 
	HarmonicMapViewer::ComputeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iRepeats, ParameterizationMode imode)
{
	// ------------- IMPLEMENT HERE ---------
	// TASK 5.3 Compute texture coordinates for textured mesh
	// rendering using a texture image of dimension iTextureWidth*iTextureHeights and iRepeats repeats
	// and store them in a mesh property.
	// If imode is equals to Uniform, compute the texture coordinates using the
	// parameters stored in vparam_u_ and store the result in texcoord_u_.
	// If imode is equals to Harmonic, compute the texture coordinates using the
	// parameters stored in vparam_h_ and store the result in texcoord_h_.
	// ------------- IMPLEMENT HERE ---------

	OpenMesh::VPropHandleT<OpenMesh::Vec2f> *vvparam, *ttexcoord;
	vvparam = &vparam();				//(imode == Harmonic)?&vparam_h_:&vparam_u_;
	ttexcoord = &texcoord();//(imode==Harmonic)?&texcoord_h_:&texcoord_u_;

	for(Mesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		Vec2f coord = ((vparam(viter) + Vec2f(1,1))/2.0) * (iRepeats);
		texcoord(viter) = coord;
	}

}

//-----------------------------------------------------------------------------

void HarmonicMapViewer::
	calc_distortion(ParameterizationMode imode)
{
	float angle_distortion=0., area_distortion=0.;

	// ------------- IMPLEMENT HERE ---------
	// TASK 5.4 Compute distortion of triangle areas and angles
	// and print it in the output window.
	// If imode is equals to Uniform, uniform map distortion has to be 
	// computed.
	// If imode is equals to Harmonic, harmonic map distortion has to be
	// computed.
	// ------------- IMPLEMENT HERE ---------

	float total_origin_area = 0;
	float total_param_area = 0;

	OpenMesh::VPropHandleT<OpenMesh::Vec2f> *vvparam;
	// Collect areas per face, first=3D Area, second=2D Area
	std::vector<std::pair<float, float>> face_to_areas;//(mesh_.n_faces());

	vvparam = &vparam(); // (imode == Harmonic)?&vparam_h_:&vparam_u_;
	// Collect angles
	for(Mesh::FaceIter fiter = mesh_.faces_begin(); fiter != mesh_.faces_end(); ++fiter) {
		Mesh::FVIter fvit(mesh_.fv_iter(fiter));
		Mesh::Point p0(mesh_.point(fvit));
		Mesh::VertexHandle h0 = fvit;;
		++fvit;
		Mesh::Point p1(mesh_.point(fvit));
		Mesh::VertexHandle h1 = fvit;
		++fvit;
		Mesh::Point p2(mesh_.point(fvit));
		Mesh::VertexHandle h2 = fvit;
		++fvit;
		
		Vec2f j0 = mesh_.property(*vvparam, h0);
		Vec2f j1 = mesh_.property(*vvparam, h1);
		Vec2f j2 = mesh_.property(*vvparam, h2);

		Vec3f p0p1 = p0 - p1;
		Vec3f p1p2 = p1 - p2;
		Vec3f p2p0 = p2 - p0;
		
		Vec2f j0j1 = j0 - j1;
		Vec2f j1j2 = j1 - j2;
		Vec2f j2j0 = j2 - j0;

		float originArea = (p0p1%(-p2p0)).length()/2.0;
		float paramArea = (
			Vec3f(0, j0j1[0], j0j1[1]) %
			(-Vec3f(0, j2j0[0], j2j0[1]))).length()/2.0;

		total_origin_area += originArea;
		total_param_area += paramArea;
		face_to_areas.push_back(std::pair<float, float>(originArea, paramArea));

		float oa = acos((p0p1.normalize()|(-p2p0).normalize()));
		float pa = acos(j0j1.normalize()|(-j2j0).normalize());

		float ob = acos((-p0p1).normalize()|p1p2.normalize());
		float pb = acos((-j0j1).normalize()|j1j2.normalize());

		float oc = acos(p2p0.normalize()|(-p1p2).normalize());
		float pc = acos(j2j0.normalize()|(-j1j2).normalize());

		angle_distortion += (oa-pa)*(oa-pa) + (ob-pb)*(ob-pb) + (oc-pc)*(oc-pc);
	}
	
	for(int i=0; i < face_to_areas.size(); ++i) {
		float originArea = face_to_areas[i].first;
		float paramArea = face_to_areas[i].second;
		float dist = (originArea/total_origin_area) - (paramArea/total_param_area);
		area_distortion += dist*dist;
	}

	cout << "Parameterization distortion: \n * ";
	switch(imode) {
	case Uniform: cout << "Uniform"; break;
	case Harmonic:cout << "Harmonic"; break;
	case Mean:cout << "Mean"; break;
	}
	cout << "\n";
	cout << " ---- Angle distortion: " << angle_distortion << " Area distortion: " << area_distortion << "\n";
}

//-----------------------------------------------------------------------------


static bool _ParameterizationComputed_u = false, _ParameterizationComputed_h = false, _ParameterizationComputed_m = false; 
static bool _BoundingBox2DComputed = false;

void 
	HarmonicMapViewer::
	draw(const std::string& _draw_mode)
{
	if (indices_.empty())
	{
		MeshViewer::draw(_draw_mode);
		return;
	}
	if (_draw_mode == "UV Domain")
	{
		if (_ParameterizationMode_!=NoParameterization)
		{
			OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = vparam();//_ParameterizationMode_ == Uniform? vparam_u_: vparam_h_);
			Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
			for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
			{
				OpenMesh::Vec2f UVCoord = mesh_.property(TexCoordHandle,v_it);
				mesh_.set_point(v_it, Mesh::Point(-UVCoord[0], UVCoord[1], 0.));
			}
			mesh_.update_normals();

			if (!_BoundingBox2DComputed)
			{
				Mesh::ConstVertexIter  v_it(mesh_.vertices_begin()), 
					v_end(mesh_.vertices_end());
				_bbMin2D = _bbMax2D = mesh_.point(v_it);
				for (; v_it!=v_end; ++v_it)
				{
					_bbMin2D.minimize(mesh_.point(v_it));
					_bbMax2D.maximize(mesh_.point(v_it));
				}
			}

			set_scene( (Vec3f)(_bbMin2D + _bbMax2D)*0.5, 0.5*(_bbMin2D - _bbMax2D).norm());

			MeshViewer::draw("Wireframe");
		}
	}
	else 
	{
		Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
		for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
			mesh_.set_point(v_it, mesh_.property(vpos_,v_it));
		mesh_.update_normals();

		set_scene( (Vec3f)(_bbMin3D + _bbMax3D)*0.5, 0.5*(_bbMin3D - _bbMax3D).norm());

		if (_draw_mode == "Textured mesh" && _ParameterizationMode_!=NoParameterization)
		{
			float *TextureCoord = TextureCoordinates();//(_ParameterizationMode_ == Uniform? _TextureCoordinates_u: _TextureCoordinates_h);
			OpenMesh::VPropHandleT<OpenMesh::Vec2f> & TexCoordHandle = texcoord(); //(_ParameterizationMode_ == Uniform? texcoord_u_: texcoord_h_);
//			if (!TextureCoord)
			{
				int nvertices = mesh_.n_vertices();
				
				if(TextureCoord)
					delete [] TextureCoord;

				TextureCoord = new float[2*nvertices];

				Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
				int index=0;
				for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
				{
					OpenMesh::Vec2f UVParam = mesh_.property(TexCoordHandle,v_it);
					TextureCoord[index++] = UVParam[0];
					TextureCoord[index++] = UVParam[1];
				}

			}
			glEnable( GL_TEXTURE_2D ); 

			glEnable(GL_LIGHTING);
			glShadeModel(GL_SMOOTH);

			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glEnableClientState(GL_TEXTURE_COORD_ARRAY);
			GL::glVertexPointer(mesh_.points());
			GL::glNormalPointer(mesh_.vertex_normals());
			GL::glTexCoordPointer(2, GL_FLOAT, 0, TextureCoord);

			glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			glDisableClientState(GL_TEXTURE_COORD_ARRAY);

			glDisable( GL_TEXTURE_2D );
		}
		else MeshViewer::draw(_draw_mode);
	}
}


//-----------------------------------------------------------------------------
static std::string currentDrawMode;


void
	HarmonicMapViewer::
	keyboard(int key, int x, int y)
{
	switch (toupper(key))
	{ 
	case 'O':
		{
			OPENFILENAME ofn={0};
			char szFileName[MAX_PATH]={0};
			ofn.lStructSize=sizeof(OPENFILENAME);
			ofn.Flags=OFN_ALLOWMULTISELECT|OFN_EXPLORER;
			ofn.lpstrFilter="All Files (*.*)\0*.*\0";
			ofn.lpstrFile=szFileName;
			ofn.nMaxFile=MAX_PATH;
			if(GetOpenFileName(&ofn))
				open_mesh(szFileName);
			_ParameterizationComputed_u = false;
			_ParameterizationComputed_h = false;
			_ParameterizationComputed_m = false;
		}
		break;
	case 'U':

		_ParameterizationMode_ = Uniform;
		if (!_ParameterizationComputed_u)
		{
			calc_uniform_parameterization();
			ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_);
			calc_distortion(_ParameterizationMode_);

			_ParameterizationComputed_u = true;      
		}
		glutPostRedisplay();
		break;

	case 'H':
		_ParameterizationMode_ = Harmonic;
		if (!_ParameterizationComputed_h)
		{
			calc_harmonic_parameterization();
			ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_);
			calc_distortion(_ParameterizationMode_);

			_ParameterizationComputed_h = true;       
		}
		glutPostRedisplay();
		break;

	case 'M':
		_ParameterizationMode_ = Mean;
		if (!_ParameterizationComputed_m)
		{
			calc_mean_parameterization();
			ComputeTextureCoordinates(256, 256, _Repeat, _ParameterizationMode_);
			calc_distortion(_ParameterizationMode_);

			_ParameterizationComputed_m = true;       
		}
		glutPostRedisplay();
		break;

		break;
	case 'R': 
		_Repeat++;
		ComputeTextureCoordinates(256, 256, _Repeat, Uniform);
		ComputeTextureCoordinates(256, 256, _Repeat, Harmonic);
		printf("Number of repeats: %d\n",_Repeat);
		
		glutPostRedisplay();
		break;

	case 'E': 
		if (_Repeat>1) _Repeat--;
		ComputeTextureCoordinates(256, 256, _Repeat, Uniform);
		ComputeTextureCoordinates(256, 256, _Repeat, Harmonic);
		printf("Number of repeats: %d\n",_Repeat);
		glutPostRedisplay();
		break;

	default:
		MeshViewer::keyboard(key, x, y);
		break;

	}
}
//-----------------------------------------------------------------------------

//=============================================================================
