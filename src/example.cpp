#include "dynamicvoronoi.h"

#include <iostream>
#include <fstream>
#include <cstring>

void loadPGM( std::istream &is, int *sizeX, int *sizeY, bool ***map ) {
  std::string tag;
  is >> tag;
  if (tag!="P5") {
    std::cerr << "Awaiting 'P5' in pgm header, found " << tag << std::endl;
    exit(-1);
  }
  
  while (is.peek()==' ' || is.peek()=='\n') is.ignore();
  while (is.peek()=='#') is.ignore(255, '\n');
  is >> *sizeX;
  while (is.peek()=='#') is.ignore(255, '\n');
  is >> *sizeY;
  while (is.peek()=='#') is.ignore(255, '\n');
  is >> tag;
  if (tag!="255") {
    std::cerr << "Awaiting '255' in pgm header, found " << tag << std::endl;
    exit(-1);
  }
  is.ignore(255, '\n');  
  
  *map = new bool*[*sizeX];

  for (int x=0; x<*sizeX; x++) {
    (*map)[x] = new bool[*sizeY];
  }
  for (int y=*sizeY-1; y>=0; y--) {
    for (int x=0; x<*sizeX; x++) {
      int c = is.get();
      if ((double)c<255-255*0.2) (*map)[x][y] = true; // cell is occupied
      else (*map)[x][y] = false; // cell is free
      if (!is.good()) {
        std::cerr << "Error reading pgm map.\n";
        exit(-1);
      }
    }
  }
}


int main( int argc, char *argv[] ) {

  if(argc<2 || argc>3 || (argc==3 && !(strcmp(argv[2],"prune")==0 || strcmp(argv[2],"pruneAlternative")==0))) {
    std::cerr<<"usage: "<<argv[0]<<" <pgm map> [prune|pruneAlternative]\n";
    exit(-1);
  }
  
  bool doPrune = false;
  bool doPruneAlternative = false;
  if (argc==3) doPrune = true;

  if(doPrune && strcmp(argv[2],"pruneAlternative")==0){
    doPrune = false;
    doPruneAlternative = true;
  }


  // LOAD PGM MAP AND INITIALIZE THE VORONOI
  std::ifstream is(argv[1]);
  if (!is) {
    std::cerr << "Could not open_queue_ map file for reading.\n";
    exit(-1);
  }
  
  bool **map=nullptr;
  int sizeX, sizeY;
  loadPGM( is, &sizeX, &sizeY, &map );
  is.close();
  fprintf(stderr, "Map loaded (%dx%d).\n", sizeX, sizeY);

  // create the voronoi object and initialize it with the map
  DynamicVoronoi voronoi;
    voronoi.InitializeMap(sizeX, sizeY, map);
    voronoi.Update(); // update distance map and Voronoi diagram
  if (doPrune) voronoi.Prune();  // prune the Voronoi
  if (doPruneAlternative) voronoi.UpdateAlternativePrunedDiagram();  // prune the Voronoi

    voronoi.Visualize("initial.ppm");
  std::cerr << "Generated initial frame.\n";
  return 0;
}
