/*===========================================================================*\
 *                                                                           *
 *                                IsoEx                                      *
 *        Copyright (C) 2002 by Computer Graphics Group, RWTH Aachen         *
 *                         www.rwth-graphics.de                              *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>

#include <IsoEx/Implicits/ImplicitSphere.hh>
#include <IsoEx/Implicits/CSG.hh>

#include <IsoEx/Grids/ImplicitGrid.hh>

#include <IsoEx/Extractors/MarchingCubesT.hh>
#include <IsoEx/Extractors/ExtendedMarchingCubesT.hh>

#include <unistd.h>


//-----------------------------------------------------------------------------


using namespace IsoEx;
using namespace OpenMesh;
using namespace OpenMesh::IO;


//-----------------------------------------------------------------------------


// Define the mesh to be used: need vertex normal and status for EMC
struct MyTraits : public DefaultTraits
{
  VertexAttributes   (Attributes::Normal | Attributes::Status);
  HalfedgeAttributes (Attributes::PrevHalfedge);
};
typedef TriMesh_ArrayKernelT<MyTraits>  MyMesh;


//-----------------------------------------------------------------------------



void usage(const char* _argv0)
{
  std::cerr << "\n\nUsage: \n"
	    << _argv0 << "  <-e | -m> <-a angle> <-r resolution> <-o filename> \n\n";
  
  std::cerr << "  -e   Use Extended Marching Cubes (default)\n"
	    << "  -m   Use standard Marching Cubes\n"
	    << "  -a   Feature detection threshold\n"
	    << "  -r   Grid resolution (default is 50)\n"
	    << "  -o   Write result to filename (should be *.{off,obj,stl}), "
	    << "defaults to output.off\n"
	    << "\n";

  exit(1);
}


//-----------------------------------------------------------------------------


int main(int argc, char** argv)
{
  // parameters
  const char*       filename = "output.off";
  unsigned int      res      = 50;
  enum { MC, EMC }  mode     = EMC;
  float             angle    = 30.0;



  // parse command line
  int         c;
  extern char *optarg;
  //extern int  optind;

  while ((c = getopt(argc, argv, "a:ehmo:r:")) != -1)
  {
    switch (c)
    {
      case 'a':
      {
	angle = atof(optarg);
	break;
      }

      case 'e':
      {
	mode = EMC;
	break;
      }

      case 'm':
      {
	mode = MC;
	break;
      }

      case 'o':
      {
	filename = optarg;
	break;
      }

      case 'r':
      {
	res = atoi(optarg);
	break;
      }

      case 'h':
      default:
      {
	usage(argv[0]);
	break;
      }
    }
  }



  // output parameters
  switch (mode)
  {
    case MC:  
      std::cout << "Standard Marching Cubes\n"; 
      break;

    case EMC: 
      std::cout << "Extended Marching Cubes\n"
		<< "Feature detection angle: " << angle 
		<< std::endl;
      break;
  }
  std::cout << "Grid: " << res << "x" << res << "x" << res << std::endl;
  std::cout << "Output: " << filename << std::endl;

  


  // construct union and diff of 3 spheres
  ImplicitSphere     s1(Vec3f(-0.5,  0.0,  0.0), 1.0);
  ImplicitSphere     s2(Vec3f( 0.5,  0.5,  0.3), 0.7);
  ImplicitSphere     s3(Vec3f( 0.1,  0.0,  1.0), 0.5);
  CSG::Union         i1(s1, s2);
  CSG::Difference    i2(i1, s3);



  // define the grid
  ImplicitGrid grid(i2,               // implicit
		    Vec3f(-2,-2,-2),  // origin
		    Vec3f(4,0,0),     // x-axis
		    Vec3f(0,4,0),     // y-axis
		    Vec3f(0,0,4),     // z-axis
		    res, res, res);   // resolution




  // extract 0-level isosurface
  MyMesh  mesh;
  switch (mode)
  {
    case MC:  
      grid.build_scalar_distance_cache();
      marching_cubes(grid, mesh); 
      break;

    case EMC:
      grid.build_is_inside_cache();
      extended_marching_cubes(grid, mesh, angle);
      break;
  }



  // write result
  write_mesh(mesh, filename);


  return 0;
}


//-----------------------------------------------------------------------------
