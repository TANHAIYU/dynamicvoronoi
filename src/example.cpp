#include "dynamicvoronoi.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

void loadPGM(std::istream &is, int *sizeX, int *sizeY, bool ***map) {
    std::string tag;
    is >> tag;
    if (tag != "P5") {
        std::cerr << "Awaiting 'P5' in pgm header, found " << tag << std::endl;
        exit(-1);
    }

    while (is.peek() == ' ' || is.peek() == '\n') is.ignore();
    while (is.peek() == '#') is.ignore(255, '\n');
    is >> *sizeX;
    while (is.peek() == '#') is.ignore(255, '\n');
    is >> *sizeY;
    while (is.peek() == '#') is.ignore(255, '\n');
    is >> tag;
    if (tag != "255") {
        std::cerr << "Awaiting '255' in pgm header, found " << tag << std::endl;
        exit(-1);
    }
    is.ignore(255, '\n');

    *map = new bool *[*sizeX];

    for (int x = 0; x < *sizeX; x++) {
        (*map)[x] = new bool[*sizeY];
    }
    for (int y = *sizeY - 1; y >= 0; y--) {
        for (int x = 0; x < *sizeX; x++) {
            int c = is.get();
            if ((double) c < 255 - 255 * 0.2) (*map)[x][y] = true; // cell is occupied
            else (*map)[x][y] = false; // cell is free
            if (!is.good()) {
                std::cerr << "Error reading pgm map.\n";
                exit(-1);
            }
        }
    }
}


int main(int argc, char *argv[]) {

    if (argc < 2 || argc > 3 ||
        (argc == 3 && !(strcmp(argv[2], "prune") == 0 || strcmp(argv[2], "pruneAlternative") == 0))) {
        std::cerr << "usage: " << argv[0] << " <pgm map> [prune|pruneAlternative]\n";
        exit(-1);
    }

    bool doPrune = false;
    bool doPruneAlternative = false;
    if (argc == 3) doPrune = true;

    if (doPrune && strcmp(argv[2], "pruneAlternative") == 0) {
        doPrune = false;
        doPruneAlternative = true;
    }


    // LOAD PGM MAP AND INITIALIZE THE VORONOI
    std::ifstream is(argv[1]);
    if (!is) {
        std::cerr << "Could not open_queue_ map file for reading.\n";
        exit(-1);
    }

    bool **map = nullptr;
    int sizeX, sizeY;
    loadPGM(is, &sizeX, &sizeY, &map);
    is.close();
    fprintf(stderr, "Map loaded (%dx%d).\n", sizeX, sizeY);

    // create the voronoi_manager object and initialize it with the map
    DynamicVoronoi voronoi_manager;
    voronoi_manager.InitializeMap(sizeX, sizeY, map);

    auto start = std::chrono::system_clock::now();
    voronoi_manager.Update(); // update distance map and voronoi diagram

    if (doPrune)
        voronoi_manager.Prune();  // prune the voronoi
    if (doPruneAlternative)
        voronoi_manager.UpdateAlternativePrunedDiagram();  // prune the voronoi

    auto init_end   = std::chrono::system_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(init_end - start);
    std::cout <<  "init_stage cost: "
              << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000
              << "ms" << std::endl;


    voronoi_manager.Visualize("initial.ppm");

    auto visu_end   = std::chrono::system_clock::now();
    duration = duration_cast<std::chrono::microseconds>(visu_end - init_end);
    std::cout <<  "init_visu cost: "
              << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000
              << "ms" << std::endl;
    std::cerr << "Generated initial frame.\n";
    // now perform some updates with random obstacles
    int num_pts = 10 + sizeX * sizeY * 0.005;
    for (int frame = 1; frame <= 1; frame++) {
        std::vector<Eigen::Vector2i> new_obstacles;
        for (int i = 0; i < num_pts; i++) {
            double x = 2 + rand() / (double) RAND_MAX * (sizeX - 4);
            double y = 2 + rand() / (double) RAND_MAX * (sizeY - 4);
            new_obstacles.push_back(Eigen::Vector2i(x, y));
        }
        voronoi_manager.ExchangeObstacles(new_obstacles); // register the new obstacles (old ones will be removed)
        voronoi_manager.Update();
        if (doPrune)
            voronoi_manager.Prune();
        std::string filename = "frame" + std::to_string(frame) + ".pgm";
        voronoi_manager.Visualize(filename.data());
        std::cerr << "Performed update with random obstacles.\n";
    }

    auto update_end   = std::chrono::system_clock::now();
    duration = duration_cast<std::chrono::microseconds>(update_end - visu_end);
    std::cout <<  "update cost: "
              << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000
              << "ms" << std::endl;

    // now remove all random obstacles again.
    // final.pgm should be very similar to initial.pgm, except for ambiguous spots
    std::vector<Eigen::Vector2i> empty;
    voronoi_manager.ExchangeObstacles(empty);
    voronoi_manager.Update();
    if (doPrune)
        voronoi_manager.Prune();
    voronoi_manager.Visualize("final.ppm");
    auto end   = std::chrono::system_clock::now();
    duration = duration_cast<std::chrono::microseconds>(end - start);
    std::cout <<  "all cost: "
         << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den * 1000
         << "ms" << std::endl;

    std::cerr << "Done with final update (all random obstacles removed).\n";
    return 0;
}
