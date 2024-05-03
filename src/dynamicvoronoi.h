#ifndef _DYNAMICVORONOI_H_
#define _DYNAMICVORONOI_H_


#include <cstdlib>
#include <cstdio>
#include <climits>
#include <queue>

#include "bucketedqueue.h"

//! A DynamicVoronoi object computes and updates a distance map and Voronoi diagram.
class DynamicVoronoi {

public:

  DynamicVoronoi();
  ~DynamicVoronoi();

  //! Initialization with an empty map
  void InitializeEmpty(int _sizeX, int _sizeY, bool initGridMap=true);
  //! Initialization with a given binary map (false==free, true==occupied)
  void InitializeMap(int _sizeX, int _sizeY, bool** _gridMap);

  //! add an obstacle at the specified cell coordinate
  void OccupyCell(int x, int y);
  //! remove an obstacle at the specified cell coordinate
  void ClearCell(int x, int y);
  //! remove old dynamic obstacles and add the new ones
  void ExchangeObstacles(std::vector<Eigen::Vector2i>& new_obstacles);

  //! update distance map and Voronoi diagram to reflect the changes
  void Update(bool update_real_dist=true);
  //! prune the Voronoi diagram
  void Prune();
  //! prune the Voronoi diagram by globally revisiting all Voronoi nodes. Takes more time but gives a more sparsely pruned Voronoi graph. You need to call this after every call to udpate()
  void UpdateAlternativePrunedDiagram();
  //! retrieve the alternatively pruned diagram. see UpdateAlternativePrunedDiagram()
  int** AlternativePrunedDiagram(){
    return alternative_diagram_;
  };
  //! retrieve the number of neighbors that are Voronoi nodes (4-connected)
  int GetNumVoronoiNeighborsAlternative(int x, int y);
  //! returns whether the specified cell is part of the alternatively pruned diagram. See UpdateAlternativePrunedDiagram.
  bool IsVoronoiAlternative(int x, int y );

  //! returns the obstacle distance at the specified location
  float GetDistance(int x, int y );
  //! returns whether the specified cell is part of the (pruned) Voronoi graph
  bool IsVoronoi(int x, int y );
  //! checks whether the specficied location is occupied
  bool IsOccupied(int x, int y);
  //! write the current distance map and voronoi diagram as ppm file
  void Visualize(const char* filename="result.ppm");

  //! returns the horizontal size of the workspace/map
  unsigned int GetSizeX() const {return sizeX;}
  //! returns the vertical size of the workspace/map
  unsigned int GetSizeY() const {return size_y;}

private:
  struct Cell {
    float dist;
    char voronoi;
    char queueing;
    int obstX;
    int obstY;
    bool needsRaise;
    int sqdist;
  };

  typedef enum {voronoiKeep=-4, freeQueued = -3, voronoiRetry=-2, voronoiPrune=-1, free=0, occupied=1} State;
  typedef enum {fwNotQueued=1, fwQueued=2, fwProcessed=3, bwQueued=4, bwProcessed=1} QueueingState;
  typedef enum {invalidObstData = SHRT_MAX/2} ObstDataState;
  typedef enum {pruned, keep, retry} MarkerMatchResult;



  // methods
  void SetObstacle(int x, int y);
  void RemoveObstacle(int x, int y);
  inline void CheckVoro(int x, int y, int nx, int ny, Cell& c, Cell& nc);
  void ReCheckVoro();
  void CommitAndColorize(bool updateRealDist=true);
  inline void ReviveVoroNeighbors(int &x, int &y);

  static inline bool IsOccupied(int &x, int &y, Cell &c);
  inline MarkerMatchResult MarkerMatch(int x, int y);
  inline bool MarkerMatchAlternative(int x, int y);
  inline int GetVoronoiPruneValence(int x, int y);

  // queues

  BucketPrioQueue<Eigen::Vector2i> open_queue_;
  std::queue<Eigen::Vector2i> prune_queue_;
  BucketPrioQueue<Eigen::Vector2i> sorted_prune_queue_;

  std::vector<Eigen::Vector2i> remove_list_;
  std::vector<Eigen::Vector2i> add_list_;
  std::vector<Eigen::Vector2i> last_obstacles_;

  // maps
  int size_y;
  int sizeX;
  Cell** data_;
  bool** grid_map_;
  bool allocated_grid_map_;

  // parameters
  int padding;
  double doubleThreshold;
  double sqrt2;
  int** alternative_diagram_;
};


#endif

