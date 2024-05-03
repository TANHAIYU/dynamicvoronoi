#include "dynamicvoronoi.h"

#include <cmath>
#include <iostream>

DynamicVoronoi::DynamicVoronoi() : size_y(), sizeX(), padding(), doubleThreshold() {
  sqrt2 = sqrt(2.0);
    data_ = nullptr;
    grid_map_ = nullptr;
    alternative_diagram_ = nullptr;
    allocated_grid_map_ = false;
}

DynamicVoronoi::~DynamicVoronoi() {
  if (data_) {
    for (int x=0; x<sizeX; x++) delete[] data_[x];
    delete[] data_;
  }
  if (allocated_grid_map_ && grid_map_) {
    for (int x=0; x<sizeX; x++) delete[] grid_map_[x];
    delete[] grid_map_;
  }
}

void DynamicVoronoi::InitializeEmpty(int _sizeX, int _sizeY, bool initGridMap) {
  if (data_) {
    for (int x=0; x<sizeX; x++) delete[] data_[x];
    delete[] data_;
      data_ = nullptr;
  }
  if(alternative_diagram_){
    for (int x=0; x<sizeX; x++) delete[] alternative_diagram_[x];
    delete[] alternative_diagram_;
      alternative_diagram_ = nullptr;
  }
  if (initGridMap) {
    if (allocated_grid_map_ && grid_map_) {
      for (int x=0; x<sizeX; x++) delete[] grid_map_[x];
      delete[] grid_map_;
        grid_map_ = nullptr;
        allocated_grid_map_ = false;
    }
  }


  sizeX = _sizeX;
    size_y = _sizeY;
    data_ = new Cell*[sizeX];
  for (int x=0; x<sizeX; x++) data_[x] = new Cell[size_y];

  if (initGridMap) {
      grid_map_ = new bool*[sizeX];
    for (int x=0; x<sizeX; x++) grid_map_[x] = new bool[size_y];
      allocated_grid_map_ = true;
  }
  
  Cell c = {};
  c.dist = INFINITY;
  c.sqdist = INT_MAX;
  c.obstX = invalidObstData;
  c.obstY = invalidObstData;
  c.voronoi = free;
  c.queueing = fwNotQueued;
  c.needsRaise = false;

  for (int x=0; x<sizeX; x++)
    for (int y=0; y < size_y; y++) data_[x][y] = c;

  if (initGridMap) {
    for (int x=0; x<sizeX; x++) 
      for (int y=0; y < size_y; y++) grid_map_[x][y] = false;
  }
}

void DynamicVoronoi::InitializeMap(int _sizeX, int _sizeY, bool** _gridMap) {
    grid_map_ = _gridMap;
    InitializeEmpty(_sizeX, _sizeY, false);

  for (int x=0; x<sizeX; x++) {
    for (int y=0; y < size_y; y++) {
      if (grid_map_[x][y]) {
        Cell c = data_[x][y];
        if (!IsOccupied(x, y, c)) {
          
          bool isSurrounded = true;
          for (int dx=-1; dx<=1; dx++) {
            int nx = x+dx;
            if (nx<=0 || nx>=sizeX-1) continue;
            for (int dy=-1; dy<=1; dy++) {
              if (dx==0 && dy==0) continue;
              int ny = y+dy;
              if (ny<=0 || ny >= size_y - 1) continue;

              if (!grid_map_[nx][ny]) {
                isSurrounded = false;
                break;
              }
            }
          }
          if (isSurrounded) {
            c.obstX = x;
            c.obstY = y;
            c.sqdist = 0;
            c.dist=0;
            c.voronoi=occupied;
            c.queueing = fwProcessed;
              data_[x][y] = c;
          } else SetObstacle(x, y);
        }
      }
    }
  }
}

void DynamicVoronoi::OccupyCell(int x, int y) {
    grid_map_[x][y] = true;
    SetObstacle(x, y);
}
void DynamicVoronoi::ClearCell(int x, int y) {
    grid_map_[x][y] = false;
    RemoveObstacle(x, y);
}

void DynamicVoronoi::SetObstacle(int x, int y) {
  Cell c = data_[x][y];
  if(IsOccupied(x, y, c)) return;
  
  add_list_.push_back(Eigen::Vector2i(x, y));
  c.obstX = x;
  c.obstY = y;
    data_[x][y] = c;
}

void DynamicVoronoi::RemoveObstacle(int x, int y) {
  Cell c = data_[x][y];
  if(!IsOccupied(x, y, c)) return;

  remove_list_.push_back(Eigen::Vector2i(x, y));
  c.obstX = invalidObstData;
  c.obstY  = invalidObstData;    
  c.queueing = bwQueued;
    data_[x][y] = c;
}

void DynamicVoronoi::ExchangeObstacles(std::vector<Eigen::Vector2i>& new_obstacles) {

  for (auto& lastObstacle : last_obstacles_) {
    int x = lastObstacle.x();
    int y = lastObstacle.y();

    bool v = grid_map_[x][y];
    if (v) continue;
      RemoveObstacle(x, y);
  }  

  last_obstacles_.clear();

  for (auto& point : new_obstacles) {
    int x = point.x();
    int y = point.y();
    bool v = grid_map_[x][y];
    if (v) continue;
      SetObstacle(x, y);
    last_obstacles_.push_back(point);
  }  
}

void DynamicVoronoi::Update(bool update_real_dist) {

    CommitAndColorize(update_real_dist);

  while (!open_queue_.empty()) {
    Eigen::Vector2i p = open_queue_.pop();
    int x = p.x();
    int y = p.y();
    Cell c = data_[x][y];

    if(c.queueing==fwProcessed) continue; 

    if (c.needsRaise) {
      // RAISE
      for (int dx=-1; dx<=1; dx++) {
        int nx = x+dx;
        if (nx<=0 || nx>=sizeX-1) continue;
        for (int dy=-1; dy<=1; dy++) {
          if (dx==0 && dy==0) continue;
          int ny = y+dy;
          if (ny<=0 || ny >= size_y - 1) continue;
          Cell nc = data_[nx][ny];
          if (nc.obstX!=invalidObstData && !nc.needsRaise) {
            if(!IsOccupied(nc.obstX, nc.obstY, data_[nc.obstX][nc.obstY])) {
              open_queue_.push(nc.sqdist, Eigen::Vector2i(nx, ny));
              nc.queueing = fwQueued;
              nc.needsRaise = true;
              nc.obstX = invalidObstData;
              nc.obstY = invalidObstData;
              if (update_real_dist) nc.dist = INFINITY;
              nc.sqdist = INT_MAX;
                data_[nx][ny] = nc;
            } else {
              if(nc.queueing != fwQueued){
                open_queue_.push(nc.sqdist, Eigen::Vector2i(nx, ny));
                nc.queueing = fwQueued;
                  data_[nx][ny] = nc;
              }
            }      
          }
        }
      }
      c.needsRaise = false;
      c.queueing = bwProcessed;
        data_[x][y] = c;
    }
    else if (c.obstX != invalidObstData && IsOccupied(c.obstX, c.obstY, data_[c.obstX][c.obstY])) {

      // LOWER
      c.queueing = fwProcessed;
      c.voronoi = occupied;

      for (int dx=-1; dx<=1; dx++) {
        int nx = x+dx;
        if (nx<=0 || nx>=sizeX-1) continue;
        for (int dy=-1; dy<=1; dy++) {
          if (dx==0 && dy==0) continue;
          int ny = y+dy;
          if (ny<=0 || ny >= size_y - 1) continue;
          Cell nc = data_[nx][ny];
          if(!nc.needsRaise) {
            int distx = nx-c.obstX;
            int disty = ny-c.obstY;
            int newSqDistance = distx*distx + disty*disty;		
            bool overwrite =  (newSqDistance < nc.sqdist);
            if(!overwrite && newSqDistance==nc.sqdist) { 
              if (nc.obstX == invalidObstData || !IsOccupied(nc.obstX, nc.obstY, data_[nc.obstX][nc.obstY])) overwrite = true;
            }
            if (overwrite) {
              open_queue_.push(newSqDistance, Eigen::Vector2i(nx, ny));
              nc.queueing = fwQueued;
              if (update_real_dist) {
                nc.dist = sqrt((double) newSqDistance);
              }
              nc.sqdist = newSqDistance;
              nc.obstX = c.obstX;
              nc.obstY = c.obstY;
            } else {
                CheckVoro(x, y, nx, ny, c, nc);
            }
              data_[nx][ny] = nc;
          }
        }
      }
    }
      data_[x][y] = c;
  }
}

float DynamicVoronoi::GetDistance(int x, int y ) {
  if( (x>0) && (x<sizeX) && (y>0) && (y < size_y)) return data_[x][y].dist;
  else return -INFINITY;
}

bool DynamicVoronoi::IsVoronoi(int x, int y ) {
  Cell c = data_[x][y];
  return (c.voronoi==free || c.voronoi==voronoiKeep);
}

bool DynamicVoronoi::IsVoronoiAlternative(int x, int y) {
  int v = alternative_diagram_[x][y];
  return (v == free || v == voronoiKeep);
}

void DynamicVoronoi::CommitAndColorize(bool updateRealDist) {
  // ADD NEW OBSTACLES
  for (auto p : add_list_) {
    int x = p.x();
    int y = p.y();
    Cell c = data_[x][y];

    if(c.queueing != fwQueued){
      if (updateRealDist) c.dist = 0;
      c.sqdist = 0;
      c.obstX = x;
      c.obstY = y;
      c.queueing = fwQueued;
      c.voronoi = occupied;
        data_[x][y] = c;
      open_queue_.push(0, Eigen::Vector2i(x, y));
    }
  }

  // REMOVE OLD OBSTACLES
  for (auto p : remove_list_) {
    int x = p.x();
    int y = p.y();
    Cell c = data_[x][y];

    if (IsOccupied(x, y, c)) continue; // obstacle was removed and reinserted
    open_queue_.push(0, Eigen::Vector2i(x, y));
    if (updateRealDist) c.dist  = INFINITY;
    c.sqdist = INT_MAX;
    c.needsRaise = true;
      data_[x][y] = c;
  }
  remove_list_.clear();
  add_list_.clear();
}


void DynamicVoronoi::CheckVoro(int x, int y, int nx, int ny, Cell& c, Cell& nc) {

  if ((c.sqdist>1 || nc.sqdist>1) && nc.obstX!=invalidObstData) { 
    if (abs(c.obstX-nc.obstX) > 1 || abs(c.obstY-nc.obstY) > 1) {
      //compute dist from x,y to obstacle of nx,ny	 
      int dxy_x = x-nc.obstX;
      int dxy_y = y-nc.obstY;
      int sqdxy = dxy_x*dxy_x + dxy_y*dxy_y;
      int stability_xy = sqdxy - c.sqdist;
      if (sqdxy - c.sqdist<0) return;

      //compute dist from nx,ny to obstacle of x,y
      int dnxy_x = nx - c.obstX;
      int dnxy_y = ny - c.obstY;
      int sqdnxy = dnxy_x*dnxy_x + dnxy_y*dnxy_y;
      int stability_nxy = sqdnxy - nc.sqdist;
      if (sqdnxy - nc.sqdist <0) return;

      //which cell is added to the Voronoi diagram?
      if(stability_xy <= stability_nxy && c.sqdist>2) {
        if (c.voronoi != free) {
          c.voronoi = free;
            ReviveVoroNeighbors(x, y);
          prune_queue_.push(Eigen::Vector2i(x, y));
        }
      }
      if(stability_nxy <= stability_xy && nc.sqdist>2) {
        if (nc.voronoi != free) {
          nc.voronoi = free;
            ReviveVoroNeighbors(nx, ny);
          prune_queue_.push(Eigen::Vector2i(nx, ny));
        }
      }
    }
  }
}


void DynamicVoronoi::ReviveVoroNeighbors(int &x, int &y) {
  for (int dx=-1; dx<=1; dx++) {
    int nx = x+dx;
    if (nx<=0 || nx>=sizeX-1) continue;
    for (int dy=-1; dy<=1; dy++) {
      if (dx==0 && dy==0) continue;
      int ny = y+dy;
      if (ny<=0 || ny >= size_y - 1) continue;
      Cell nc = data_[nx][ny];
      if (nc.sqdist != INT_MAX && !nc.needsRaise && (nc.voronoi == voronoiKeep || nc.voronoi == voronoiPrune)) {
        nc.voronoi = free;
          data_[nx][ny] = nc;
        prune_queue_.push(Eigen::Vector2i(nx, ny));
      }
    }
  }
}


bool DynamicVoronoi::IsOccupied(int x, int y) {
  Cell c = data_[x][y];
  return (c.obstX==x && c.obstY==y);
}

bool DynamicVoronoi::IsOccupied(int &x, int &y, Cell &c) {
  return (c.obstX==x && c.obstY==y);
}

void DynamicVoronoi::Visualize(const char *filename) {
  // write ppm files

  FILE* F = fopen(filename, "w");
  if (!F) {
    std::cerr << "could not open_queue_ 'result.pgm' for writing!\n";
    return;
  }
  fprintf(F, "P6\n#\n");
  fprintf(F, "%d %d\n255\n", sizeX, size_y);

  for(int y = size_y - 1; y >= 0; y--){
    for(int x = 0; x<sizeX; x++){	
      unsigned char c = 0;
      if (alternative_diagram_ != nullptr && (alternative_diagram_[x][y] == free || alternative_diagram_[x][y] == voronoiKeep)) {
        fputc( 255, F );
        fputc( 0, F );
        fputc( 0, F );
      } else if(IsVoronoi(x, y)){
        fputc( 0, F );
        fputc( 0, F );
        fputc( 255, F );
      } else if (data_[x][y].sqdist == 0) {
        fputc( 0, F );
        fputc( 0, F );
        fputc( 0, F );
      } else {
        float f = 80+(sqrt(data_[x][y].sqdist) * 10);
        if (f>255) f=255;
        if (f<0) f=0;
        c = (unsigned char)f;
        fputc( c, F );
        fputc( c, F );
        fputc( c, F );
      }
    }
  }
  fclose(F);
}


void DynamicVoronoi::Prune() {
  // filler
  while(!prune_queue_.empty()) {
    Eigen::Vector2i p = prune_queue_.front();
    prune_queue_.pop();
    int x = p.x();
    int y = p.y();

    if (data_[x][y].voronoi == occupied) continue;
    if (data_[x][y].voronoi == freeQueued) continue;

      data_[x][y].voronoi = freeQueued;
    sorted_prune_queue_.push(data_[x][y].sqdist, p);

    /* tl t tr
       l c r
       bl b br */

    Cell tr,tl,br,bl;
    tr = data_[x + 1][y + 1];
    tl = data_[x - 1][y + 1];
    br = data_[x + 1][y - 1];
    bl = data_[x - 1][y - 1];

    Cell r,b,t,l;
    r = data_[x + 1][y];
    l = data_[x - 1][y];
    t = data_[x][y + 1];
    b = data_[x][y - 1];

    if (x+2<sizeX && r.voronoi==occupied) { 
      // fill to the right
      if (tr.voronoi!=occupied && br.voronoi!=occupied && data_[x + 2][y].voronoi != occupied) {
        r.voronoi = freeQueued;
        sorted_prune_queue_.push(r.sqdist, Eigen::Vector2i(x + 1, y));
          data_[x + 1][y] = r;
      }
    } 
    if (x-2>=0 && l.voronoi==occupied) { 
      // fill to the left
      if (tl.voronoi!=occupied && bl.voronoi!=occupied && data_[x - 2][y].voronoi != occupied) {
        l.voronoi = freeQueued;
        sorted_prune_queue_.push(l.sqdist, Eigen::Vector2i(x - 1, y));
          data_[x - 1][y] = l;
      }
    } 
    if (y+2 < size_y && t.voronoi == occupied) {
      // fill to the top
      if (tr.voronoi!=occupied && tl.voronoi!=occupied && data_[x][y + 2].voronoi != occupied) {
        t.voronoi = freeQueued;
        sorted_prune_queue_.push(t.sqdist, Eigen::Vector2i(x, y + 1));
          data_[x][y + 1] = t;
      }
    } 
    if (y-2>=0 && b.voronoi==occupied) { 
      // fill to the bottom
      if (br.voronoi!=occupied && bl.voronoi!=occupied && data_[x][y - 2].voronoi != occupied) {
        b.voronoi = freeQueued;
        sorted_prune_queue_.push(b.sqdist, Eigen::Vector2i(x, y - 1));
          data_[x][y - 1] = b;
      }
    } 
  }


  while(!sorted_prune_queue_.empty()) {
    Eigen::Vector2i p = sorted_prune_queue_.pop();
    Cell c = data_[p.x()][p.y()];
    int v = c.voronoi;
    if (v!=freeQueued && v!=voronoiRetry) { // || v>free || v==voronoiPrune || v==voronoiKeep) {
      //      assert(v!=retry);
      continue;
    }

    MarkerMatchResult r = MarkerMatch(p.x(), p.y());
    if (r==pruned) c.voronoi = voronoiPrune;
    else if (r==keep) c.voronoi = voronoiKeep;
    else { // r==retry
      c.voronoi = voronoiRetry;
      //      printf("RETRY %d %d\n", x, sizeY-1-y);
      prune_queue_.push(p);
    }
      data_[p.x()][p.y()] = c;

    if (sorted_prune_queue_.empty()) {
      while (!prune_queue_.empty()) {
        Eigen::Vector2i p = prune_queue_.front();
        prune_queue_.pop();
        sorted_prune_queue_.push(data_[p.x()][p.y()].sqdist, p);
      }
    }
  }
  //  printf("match: %d\nnomat: %d\n", matchCount, noMatchCount);
}

void DynamicVoronoi::UpdateAlternativePrunedDiagram() {

  if(alternative_diagram_ == nullptr){
      alternative_diagram_ = new int*[sizeX];
    for(int x=0; x<sizeX; x++){
        alternative_diagram_[x] = new int[size_y];
    }
  }


  std::queue<Eigen::Vector2i> end_cells;
  BucketPrioQueue<Eigen::Vector2i> sortedPruneQueue;
  for(int x=1; x<sizeX-1; x++){
    for(int y=1; y < size_y - 1; y++){
      Cell& c = data_[x][y];
        alternative_diagram_[x][y] = c.voronoi;
	if(c.voronoi <=free){
	  sortedPruneQueue.push(c.sqdist, Eigen::Vector2i(x,y));
	  end_cells.push(Eigen::Vector2i(x, y));
	}
    }
  }

  for(int x=1; x<sizeX-1; x++){
    for(int y=1; y < size_y - 1; y++){
      if(GetNumVoronoiNeighborsAlternative(x, y) >= 3){
          alternative_diagram_[x][y] = voronoiKeep;
	sortedPruneQueue.push(data_[x][y].sqdist, Eigen::Vector2i(x, y));
	end_cells.push(Eigen::Vector2i(x, y));
      }
    }
  }

  for(int x=1; x<sizeX-1; x++){
    for(int y=1; y < size_y - 1; y++){
      if(GetNumVoronoiNeighborsAlternative(x, y) >= 3){
          alternative_diagram_[x][y] = voronoiKeep;
	sortedPruneQueue.push(data_[x][y].sqdist, Eigen::Vector2i(x, y));
	end_cells.push(Eigen::Vector2i(x, y));
      }
    }
  }


  while (!sortedPruneQueue.empty()) {
    Eigen::Vector2i p = sortedPruneQueue.pop();

    if (MarkerMatchAlternative(p.x(), p.y())) {
        alternative_diagram_[p.x()][p.y()]=voronoiPrune;
    } else {
        alternative_diagram_[p.x()][p.y()]=voronoiKeep;
    }
  }

  // //delete worms
  while (!end_cells.empty()) {
    Eigen::Vector2i p = end_cells.front();
    end_cells.pop();

    if (IsVoronoiAlternative(p.x(), p.y()) && GetNumVoronoiNeighborsAlternative(p.x(), p.y()) == 1) {
        alternative_diagram_[p.x()][p.y()] = voronoiPrune;

      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          if (!(dx || dy) || (dx && dy)) {
            continue;
          }
          int nx = p.x() + dx;
          int ny = p.y() + dy;
          if (nx < 0 || nx >= sizeX || ny < 0 || ny >= size_y) {
            continue;
          }
          if (IsVoronoiAlternative(nx, ny)) {
            if (GetNumVoronoiNeighborsAlternative(nx, ny) == 1) {
              end_cells.push(Eigen::Vector2i(nx, ny));
            }
          }
        }
      }
    }
  }
}

bool DynamicVoronoi::MarkerMatchAlternative(int x, int y) {
// prune if this returns true

  bool f[8];

  int nx, ny;
  int dx, dy;

  int i = 0;
//  int obstacleCount=0;
  int voroCount = 0;
  for (dy = 1; dy >= -1; dy--) {
    ny = y + dy;
    for (dx = -1; dx <= 1; dx++) {
      if (dx || dy) {
        nx = x + dx;
        int v = alternative_diagram_[nx][ny];
        bool b = (v <= free && v != voronoiPrune);
        //	if (v==occupied) obstacleCount++;
        f[i] = b;
        if (v <= free && !(dx && dy))
          voroCount++;
        i++;
      }
    }
  }

  /*
   * 5 6 7
   * 3   4
   * 0 1 2
   */

  {
    //connected horizontal or vertically to only one cell
    if (voroCount == 1 && (f[1] || f[3] || f[4] || f[6])) {
      return false;
    }

    // 4-connected
    if ((!f[0] && f[1] && f[3]) || (!f[2] && f[1] && f[4]) || (!f[5] && f[3] && f[6]) || (!f[7] && f[6] && f[4]))
      return false;

    if ((f[3] && f[4] && !f[1] && !f[6]) || (f[1] && f[6] && !f[3] && !f[4]))
      return false;

  }
  return true;
}

int DynamicVoronoi::GetNumVoronoiNeighborsAlternative(int x, int y) {
  int count = 0;
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      if ((dx == 0 && dy == 0) || (dx != 0 && dy != 0)) {
        continue;
      }

      int nx = x + dx;
      int ny = y + dy;
      if (nx < 0 || nx >= sizeX || ny < 0 || ny >= size_y) {
        continue;
      }
      if (alternative_diagram_[nx][ny] == free || alternative_diagram_[nx][ny] == voronoiKeep) {
        count++;
      }
    }
  }
  return count;
}



DynamicVoronoi::MarkerMatchResult DynamicVoronoi::MarkerMatch(int x, int y) {
  // implementation of connectivity patterns
  bool f[8];

  int nx, ny;
  int dx, dy;

  int i=0;
  int count=0;
  //  int obstacleCount=0;
  int voroCount=0;
  int voroCountFour=0;

  for (dy=1; dy>=-1; dy--) {
    ny = y+dy;
    for (dx=-1; dx<=1; dx++) {
      if (dx || dy) {
        nx = x+dx;
        Cell nc = data_[nx][ny];
        int v = nc.voronoi;
        bool b = (v<=free && v!=voronoiPrune); 
        //	if (v==occupied) obstacleCount++;
        f[i] = b;
        if (b) {
          voroCount++;
          if (!(dx && dy)) voroCountFour++;
        }
        if (b && !(dx && dy) ) count++;
        //	if (v<=free && !(dx && dy)) voroCount++;
        i++;
      }
    }
  }
  if (voroCount<3 && voroCountFour==1 && (f[1] || f[3] || f[4] || f[6])) {
    //    assert(voroCount<2);
    //    if (voroCount>=2) printf("voro>2 %d %d\n", x, y);
    return keep;
  }

  // 4-connected
  if ((!f[0] && f[1] && f[3]) || (!f[2] && f[1] && f[4]) || (!f[5] && f[3] && f[6]) || (!f[7] && f[6] && f[4])) return keep;
  if ((f[3] && f[4] && !f[1] && !f[6]) || (f[1] && f[6] && !f[3] && !f[4])) return keep;
  


  // keep voro cells inside of blocks and retry later
  if (voroCount>=5 && voroCountFour>=3 && data_[x][y].voronoi != voronoiRetry) {
    return retry;
  }

  return pruned;
}
