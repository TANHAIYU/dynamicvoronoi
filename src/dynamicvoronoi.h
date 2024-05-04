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
    void InitializeEmpty(int _sizeX, int _sizeY, bool init_grid_map=true);
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
    //! prune the Voronoi diagram by globally revisiting all Voronoi nodes. Takes more time but gives a more sparsely PRUNED Voronoi graph. You need to call this after every call to udpate()
    void UpdateAlternativePrunedDiagram();
    //! retrieve the alternatively PRUNED diagram. see UpdateAlternativePrunedDiagram()
    int** AlternativePrunedDiagram(){
      return alternative_diagram_;
    };
    //! retrieve the number of neighbors that are Voronoi nodes (4-connected)
    int GetNumVoronoiNeighborsAlternative(int x, int y);
    //! returns whether the specified cell is part of the alternatively PRUNED diagram. See UpdateAlternativePrunedDiagram.
    bool IsVoronoiAlternative(int x, int y );

    //! returns the obstacle distance at the specified location
    float GetDistance(int x, int y );
    //! returns whether the specified cell is part of the (PRUNED) Voronoi graph
    bool IsVoronoi(int x, int y );
    //! checks whether the specficied location is occupied
    bool IsOccupied(int x, int y);
    //! write the current distance map and voronoi diagram as ppm file
    void Visualize(const char* filename="result.ppm");

    //! returns the horizontal size of the workspace/map
    unsigned int GetSizeX() const {return size_x;}
    //! returns the vertical size of the workspace/map
    unsigned int GetSizeY() const {return size_y;}

private:
    struct Cell {
        float dist;
        char voronoi_status;
        char queueing_status;
        int obst_x;
        int obst_y;
        bool needs_raise;
        int square_dist;
    };

    typedef enum {
        VORONOI_KEEP = -4,
        FREE_QUEUED = -3,
        VORONOI_RETRY=-2,
        VORONOI_PRUNE=-1,
        FREE = 0,
        OCCUPIED=1
    } State;

    typedef enum {
        FORWARD_INITIALIZED = 1,
        FORWARD_QUEUED = 2,
        FORWARD_PROCESSED = 3, // process done
        BACKWARD_QUEUED = 4,
        BACKWARD_PROCESSED = 1
    } QueueingState;

    typedef enum {INIT = SHRT_MAX / 2} ObstDataState;

    typedef enum {
        PRUNED = 0,
        KEEP = 1,
        RETRY = 3
    } MarkerMatchResult;



    // methods
    void ProcessRaise(bool update_real_dist, int x, int y, Cell &curr_cell);
    void ProcessLower(bool update_real_dist, int x, int y, Cell &curr_cell);

    bool IsSurrounded(int x, int y) const;
    bool NeedOverwrite(const Cell &curr_cell, int nx, int ny, Cell &neighbor_cell, int &new_square_dis) const;

    void SetObstacle(int x, int y);
    void RemoveObstacle(int x, int y);
    inline void CheckVoro(int x, int y, int nx, int ny, Cell& cell, Cell& neighbor_cell);
    void ReCheckVoro();
    void CommitAndColorize(bool update_real_dist=true);
    inline void ReviveVoroNeighbors(int &x, int &y);

    static inline bool IsOccupied(int &x, int &y, Cell &c);
    inline MarkerMatchResult MarkerMatch(int x, int y);
    inline bool MarkerMatchAlternative(int x, int y);
    inline int GetVoronoiPruneValence(int x, int y);

    // queues
    BucketPrioQueue<Eigen::Vector2i> open_queue_;
    std::queue<Eigen::Vector2i> prune_queue_;
    BucketPrioQueue<Eigen::Vector2i> sorted_prune_queue_;

    // newly become obstacle-free grid
    std::vector<Eigen::Vector2i> remove_list_;
    // newly add obstacle occupied grid
    std::vector<Eigen::Vector2i> add_list_;
    std::vector<Eigen::Vector2i> last_obstacles_;

    // maps
    int size_y;
    int size_x;
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

