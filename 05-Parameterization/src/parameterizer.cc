//=============================================================================
#include "HarmonicMapViewer.hh"



int main(int argc, char **argv)
{
  glutInit(&argc, argv);

  HarmonicMapViewer window("Harmonic Map", 512, 512, 5);

  glutMainLoop();
}
