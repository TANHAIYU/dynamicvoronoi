#include "dynamicvoronoi.h"

#include <cmath>
#include <iostream>
#include <climits>

DynamicVoronoi::DynamicVoronoi() : size_y(), size_x(), padding(), doubleThreshold() {
    sqrt2 = sqrt(2.0);
    data_ = nullptr;
    grid_map_ = nullptr;
    alternative_diagram_ = nullptr;
    allocated_grid_map_ = false;
}

DynamicVoronoi::~DynamicVoronoi() {
    if (data_) {
        for (int x = 0; x < size_x; x++) delete[] data_[x];
        delete[] data_;
    }
    if (allocated_grid_map_ && grid_map_) {
        for (int x = 0; x < size_x; x++) delete[] grid_map_[x];
        delete[] grid_map_;
    }
}

void DynamicVoronoi::InitializeEmpty(int _sizeX, int _sizeY, bool init_grid_map) {
    if (data_) {
        for (int x = 0; x < size_x; x++) {
            delete[] data_[x];
        }
        delete[] data_;
        data_ = nullptr;
    }
    if (alternative_diagram_) {
        for (int x = 0; x < size_x; x++) {
            delete[] alternative_diagram_[x];
        }
        delete[] alternative_diagram_;
        alternative_diagram_ = nullptr;
    }
    if (init_grid_map) {
        if (allocated_grid_map_ && grid_map_) {
            for (int x = 0; x < size_x; x++) {
                delete[] grid_map_[x];
            }
            delete[] grid_map_;
            grid_map_ = nullptr;
            allocated_grid_map_ = false;
        }
    }


    size_x = _sizeX;
    size_y = _sizeY;
    data_ = new Cell *[size_x];
    for (int x = 0; x < size_x; x++) {
        data_[x] = new Cell[size_y];
    }

    if (init_grid_map) {
        grid_map_ = new bool *[size_x];
        for (int x = 0; x < size_x; x++) {
            grid_map_[x] = new bool[size_y];
        }
        allocated_grid_map_ = true;
    }

    Cell c = {};
    c.dist = INFINITY;
    c.sqrt_dist = INT_MAX;
    c.obst_x = invalidObstData;
    c.obst_y = invalidObstData;
    c.voronoi_status = FREE;
    c.queueing_status = FORWARD_INITIALIZED;
    c.needs_raise = false;

    for (int x = 0; x < size_x; x++){
        for (int y = 0; y < size_y; y++) {
            data_[x][y] = c;
        }
    }

    if (init_grid_map) {
        for (int x = 0; x < size_x; x++) {
            for (int y = 0; y < size_y; y++) {
                grid_map_[x][y] = false;
            }
        }
    }
}

void DynamicVoronoi::InitializeMap(int _sizeX, int _sizeY, bool **_gridMap) {
    grid_map_ = _gridMap;
    InitializeEmpty(_sizeX, _sizeY, false);

    for (int x = 0; x < size_x; x++) {
        for (int y = 0; y < size_y; y++) {
            if (grid_map_[x][y]) {
                Cell curr_cell = data_[x][y];
                if (!IsOccupied(x, y, curr_cell)) {
                    if (IsSurrounded(x, y)) {
                        curr_cell.obst_x = x;
                        curr_cell.obst_y = y;
                        curr_cell.sqrt_dist = 0;
                        curr_cell.dist = 0;
                        curr_cell.voronoi_status = OCCUPIED;
                        curr_cell.queueing_status = FORWARD_PROCESSED;
                        data_[x][y] = curr_cell;
                    } else {
                        SetObstacle(x, y);
                    }
                }
            }
        }
    }
}

bool DynamicVoronoi::IsSurrounded(int x, int y) const {
    bool is_surrounded = true;
    for (int dx = -1; dx <= 1; dx++) {
        int nx = x + dx;
        if (nx <= 0 || nx >= size_x - 1)
            continue;
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0)
                continue;

            int ny = y + dy;
            if (ny <= 0 || ny >= size_y - 1)
                continue;

            if (!grid_map_[nx][ny]) {
                is_surrounded = false;
                break;
            }
        }
    }
    return is_surrounded;
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
    if (IsOccupied(x, y, c))
        return;

    add_list_.push_back(Eigen::Vector2i(x, y));
    c.obst_x = x;
    c.obst_y = y;
    data_[x][y] = c;
}

void DynamicVoronoi::RemoveObstacle(int x, int y) {
    Cell c = data_[x][y];
    if (!IsOccupied(x, y, c))
        return;

    remove_list_.push_back(Eigen::Vector2i(x, y));
    c.obst_x = invalidObstData;
    c.obst_y = invalidObstData;
    c.queueing_status = BACKWARD_QUEUED;
    data_[x][y] = c;
}

void DynamicVoronoi::ExchangeObstacles(std::vector<Eigen::Vector2i> &new_obstacles) {

    for (auto &lastObstacle: last_obstacles_) {
        int x = lastObstacle.x();
        int y = lastObstacle.y();

        if (grid_map_[x][y])
            continue;
        RemoveObstacle(x, y);
    }

    last_obstacles_.clear();

    for (auto &point: new_obstacles) {
        int x = point.x();
        int y = point.y();
        if (grid_map_[x][y])
            continue;
        SetObstacle(x, y);
        last_obstacles_.push_back(point);
    }
}

void DynamicVoronoi::Update(bool update_real_dist) {

    CommitAndColorize(update_real_dist);

    while (!open_queue_.empty()) {
        Eigen::Vector2i p = open_queue_.pop();
        Cell curr_cell = data_[p.x()][p.y()];

        if (curr_cell.queueing_status == FORWARD_PROCESSED)
            continue;

        if (curr_cell.needs_raise) {
            // RAISE
            ProcessRaise(update_real_dist, p.x(), p.y(), curr_cell);
        } else if (curr_cell.obst_x != invalidObstData &&
                   IsOccupied(curr_cell.obst_x, curr_cell.obst_y, data_[curr_cell.obst_x][curr_cell.obst_y])) {
            // LOWER
            ProcessLower(update_real_dist, p.x(), p.y(), curr_cell);
        }
        data_[p.x()][p.y()] = curr_cell;
    }
}

void DynamicVoronoi::ProcessLower(bool update_real_dist, int x, int y, DynamicVoronoi::Cell &curr_cell) {
    curr_cell.queueing_status = FORWARD_PROCESSED;
    curr_cell.voronoi_status = OCCUPIED;

    for (int dx = -1; dx <= 1; dx++) {
        int nx = x + dx;
        if (nx <= 0 || nx >= size_x - 1)
            continue;
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0)
                continue;
            int ny = y + dy;
            if (ny <= 0 || ny >= size_y - 1)
                continue;
            Cell neighbor_cell = data_[nx][ny];
            if (!neighbor_cell.needs_raise) {
                int distx = nx - curr_cell.obst_x;
                int disty = ny - curr_cell.obst_y;
                int newSqDistance = distx * distx + disty * disty;
                bool overwrite = (newSqDistance < neighbor_cell.sqrt_dist);
                if (!overwrite && newSqDistance == neighbor_cell.sqrt_dist) {
                    if (neighbor_cell.obst_x == invalidObstData ||
                        !IsOccupied(neighbor_cell.obst_x, neighbor_cell.obst_y, data_[neighbor_cell.obst_x][neighbor_cell.obst_y]))
                        overwrite = true;
                }
                if (overwrite) {
                    open_queue_.push(newSqDistance, Eigen::Vector2i(nx, ny));
                    neighbor_cell.queueing_status = FORWARD_QUEUED;
                    if (update_real_dist) {
                        neighbor_cell.dist = std::sqrt((double) newSqDistance);
                    }
                    neighbor_cell.sqrt_dist = newSqDistance;
                    neighbor_cell.obst_x = curr_cell.obst_x;
                    neighbor_cell.obst_y = curr_cell.obst_y;
                } else {
                    CheckVoro(x, y, nx, ny, curr_cell, neighbor_cell);
                }
                data_[nx][ny] = neighbor_cell;
            }
        }
    }
}

void DynamicVoronoi::ProcessRaise(bool update_real_dist, int x, int y, DynamicVoronoi::Cell &curr_cell) {
    for (int dx = -1; dx <= 1; dx++) {
        int nx = x + dx;
        if (nx <= 0 || nx >= size_x - 1)
            continue;
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0)
                continue;
            int ny = y + dy;
            if (ny <= 0 || ny >= size_y - 1)
                continue;
            Cell neighbor_cell = data_[nx][ny];
            if (neighbor_cell.obst_x != invalidObstData && !neighbor_cell.needs_raise) {
                if (!IsOccupied(neighbor_cell.obst_x, neighbor_cell.obst_y, data_[neighbor_cell.obst_x][neighbor_cell.obst_y])) {
                    open_queue_.push(neighbor_cell.sqrt_dist, Eigen::Vector2i(nx, ny));
                    neighbor_cell.queueing_status = FORWARD_QUEUED;
                    neighbor_cell.needs_raise = true;
                    neighbor_cell.obst_x = invalidObstData;
                    neighbor_cell.obst_y = invalidObstData;
                    if (update_real_dist)
                        neighbor_cell.dist = INFINITY;
                    neighbor_cell.sqrt_dist = INT_MAX;
                    data_[nx][ny] = neighbor_cell;
                } else {
                    if (neighbor_cell.queueing_status != FORWARD_QUEUED) {
                        open_queue_.push(neighbor_cell.sqrt_dist, Eigen::Vector2i(nx, ny));
                        neighbor_cell.queueing_status = FORWARD_QUEUED;
                        data_[nx][ny] = neighbor_cell;
                    }
                }
            }
        }
    }
    curr_cell.needs_raise = false;
    curr_cell.queueing_status = BACKWARD_PROCESSED;
    data_[x][y] = curr_cell;
}

float DynamicVoronoi::GetDistance(int x, int y) {
    if ((x > 0) && (x < size_x) && (y > 0) && (y < size_y)) return data_[x][y].dist;
    else return -INFINITY;
}

bool DynamicVoronoi::IsVoronoi(int x, int y) {
    Cell cell = data_[x][y];
    return (cell.voronoi_status == FREE || cell.voronoi_status == VORONOI_KEEP);
}

bool DynamicVoronoi::IsVoronoiAlternative(int x, int y) {
    int v = alternative_diagram_[x][y];
    return (v == FREE || v == VORONOI_KEEP);
}

void DynamicVoronoi::CommitAndColorize(bool update_real_dist) {
    // ADD NEW OBSTACLES
    for (auto p: add_list_) {
        int x = p.x();
        int y = p.y();
        Cell cell = data_[x][y];

        if (cell.queueing_status != FORWARD_QUEUED) {
            if (update_real_dist)
                cell.dist = 0;
            cell.sqrt_dist = 0;
            cell.obst_x = x;
            cell.obst_y = y;
            cell.queueing_status = FORWARD_QUEUED;
            cell.voronoi_status = OCCUPIED;
            data_[x][y] = cell;
            open_queue_.push(0, Eigen::Vector2i(x, y));
        }
    }

    // REMOVE OLD OBSTACLES
    for (auto p: remove_list_) {
        int x = p.x();
        int y = p.y();
        Cell cell = data_[x][y];

        if (IsOccupied(x, y, cell))
            continue; // obstacle was removed and reinserted
        open_queue_.push(0, Eigen::Vector2i(x, y));
        if (update_real_dist)
            cell.dist = INFINITY;
        cell.sqrt_dist = INT_MAX;
        cell.needs_raise = true;
        data_[x][y] = cell;
    }
    remove_list_.clear();
    add_list_.clear();
}


void DynamicVoronoi::CheckVoro(int x, int y, int nx, int ny, Cell &cell, Cell &neighbor_cell) {

    if ((cell.sqrt_dist > 1 || neighbor_cell.sqrt_dist > 1) && neighbor_cell.obst_x != invalidObstData) {
        if (abs(cell.obst_x - neighbor_cell.obst_x) > 1 || abs(cell.obst_y - neighbor_cell.obst_y) > 1) {
            //compute dist from x,y to obstacle of nx,ny
            int dxy_x = x - neighbor_cell.obst_x;
            int dxy_y = y - neighbor_cell.obst_y;
            int sqdxy = dxy_x * dxy_x + dxy_y * dxy_y;
            int stability_xy = sqdxy - cell.sqrt_dist;
            if (sqdxy - cell.sqrt_dist < 0)
                return;

            //compute dist from nx,ny to obstacle of x,y
            int dnxy_x = nx - cell.obst_x;
            int dnxy_y = ny - cell.obst_y;
            int sqdnxy = dnxy_x * dnxy_x + dnxy_y * dnxy_y;
            int stability_nxy = sqdnxy - neighbor_cell.sqrt_dist;
            if (sqdnxy - neighbor_cell.sqrt_dist < 0)
                return;

            //which cell is added to the Voronoi diagram?
            if (stability_xy <= stability_nxy && cell.sqrt_dist > 2) {
                if (cell.voronoi_status != FREE) {
                    cell.voronoi_status = FREE;
                    ReviveVoroNeighbors(x, y);
                    prune_queue_.push(Eigen::Vector2i(x, y));
                }
            }
            if (stability_nxy <= stability_xy && neighbor_cell.sqrt_dist > 2) {
                if (neighbor_cell.voronoi_status != FREE) {
                    neighbor_cell.voronoi_status = FREE;
                    ReviveVoroNeighbors(nx, ny);
                    prune_queue_.push(Eigen::Vector2i(nx, ny));
                }
            }
        }
    }
}


void DynamicVoronoi::ReviveVoroNeighbors(int &x, int &y) {
    for (int dx = -1; dx <= 1; dx++) {
        int nx = x + dx;
        if (nx <= 0 || nx >= size_x - 1)
            continue;
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0)
                continue;
            int ny = y + dy;
            if (ny <= 0 || ny >= size_y - 1)
                continue;
            Cell neighbor_cell = data_[nx][ny];
            if (neighbor_cell.sqrt_dist != INT_MAX &&
                !neighbor_cell.needs_raise &&
                (neighbor_cell.voronoi_status == VORONOI_KEEP || neighbor_cell.voronoi_status == VORONOI_PRUNE)) {
                neighbor_cell.voronoi_status = FREE;
                data_[nx][ny] = neighbor_cell;
                prune_queue_.push(Eigen::Vector2i(nx, ny));
            }
        }
    }
}


bool DynamicVoronoi::IsOccupied(int x, int y) {
    Cell c = data_[x][y];
    return (c.obst_x == x && c.obst_y == y);
}

bool DynamicVoronoi::IsOccupied(int &x, int &y, Cell &c) {
    return (c.obst_x == x && c.obst_y == y);
}

void DynamicVoronoi::Visualize(const char *filename) {
    // write ppm files

    FILE *F = fopen(filename, "w");
    if (!F) {
        std::cerr << "could not open_queue_ 'result.pgm' for writing!\n";
        return;
    }
    fprintf(F, "P6\n#\n");
    fprintf(F, "%d %d\n255\n", size_x, size_y);

    for (int y = size_y - 1; y >= 0; y--) {
        for (int x = 0; x < size_x; x++) {
            unsigned char c = 0;
            if (alternative_diagram_ != nullptr &&
                (alternative_diagram_[x][y] == FREE || alternative_diagram_[x][y] == VORONOI_KEEP)) {
                fputc(255, F);
                fputc(0, F);
                fputc(0, F);
            } else if (IsVoronoi(x, y)) {
                fputc(0, F);
                fputc(0, F);
                fputc(255, F);
            } else if (data_[x][y].sqrt_dist == 0) {
                fputc(0, F);
                fputc(0, F);
                fputc(0, F);
            } else {
                float f = 80 + (sqrt(data_[x][y].sqrt_dist) * 10);
                if (f > 255) f = 255;
                if (f < 0) f = 0;
                c = (unsigned char) f;
                fputc(c, F);
                fputc(c, F);
                fputc(c, F);
            }
        }
    }
    fclose(F);
}


void DynamicVoronoi::Prune() {
    // filler
    while (!prune_queue_.empty()) {
        Eigen::Vector2i p = prune_queue_.front();
        prune_queue_.pop();
        int x = p.x();
        int y = p.y();

        if (data_[x][y].voronoi_status == OCCUPIED)
            continue;
        if (data_[x][y].voronoi_status == FREE_QUEUED)
            continue;

        data_[x][y].voronoi_status = FREE_QUEUED;
        sorted_prune_queue_.push(data_[x][y].sqrt_dist, p);

        /* tl t tr
           l c r
           bl b br */

        Cell tr, tl, br, bl;
        tr = data_[x + 1][y + 1];
        tl = data_[x - 1][y + 1];
        br = data_[x + 1][y - 1];
        bl = data_[x - 1][y - 1];

        Cell r, b, t, l;
        r = data_[x + 1][y];
        l = data_[x - 1][y];
        t = data_[x][y + 1];
        b = data_[x][y - 1];

        if (x + 2 < size_x && r.voronoi_status == OCCUPIED) {
            // fill to the right
            if (tr.voronoi_status != OCCUPIED && br.voronoi_status != OCCUPIED && data_[x + 2][y].voronoi_status != OCCUPIED) {
                r.voronoi_status = FREE_QUEUED;
                sorted_prune_queue_.push(r.sqrt_dist, Eigen::Vector2i(x + 1, y));
                data_[x + 1][y] = r;
            }
        }
        if (x - 2 >= 0 && l.voronoi_status == OCCUPIED) {
            // fill to the left
            if (tl.voronoi_status != OCCUPIED && bl.voronoi_status != OCCUPIED && data_[x - 2][y].voronoi_status != OCCUPIED) {
                l.voronoi_status = FREE_QUEUED;
                sorted_prune_queue_.push(l.sqrt_dist, Eigen::Vector2i(x - 1, y));
                data_[x - 1][y] = l;
            }
        }
        if (y + 2 < size_y && t.voronoi_status == OCCUPIED) {
            // fill to the top
            if (tr.voronoi_status != OCCUPIED && tl.voronoi_status != OCCUPIED && data_[x][y + 2].voronoi_status != OCCUPIED) {
                t.voronoi_status = FREE_QUEUED;
                sorted_prune_queue_.push(t.sqrt_dist, Eigen::Vector2i(x, y + 1));
                data_[x][y + 1] = t;
            }
        }
        if (y - 2 >= 0 && b.voronoi_status == OCCUPIED) {
            // fill to the bottom
            if (br.voronoi_status != OCCUPIED && bl.voronoi_status != OCCUPIED && data_[x][y - 2].voronoi_status != OCCUPIED) {
                b.voronoi_status = FREE_QUEUED;
                sorted_prune_queue_.push(b.sqrt_dist, Eigen::Vector2i(x, y - 1));
                data_[x][y - 1] = b;
            }
        }
    }


    while (!sorted_prune_queue_.empty()) {
        Eigen::Vector2i p = sorted_prune_queue_.pop();
        Cell cell = data_[p.x()][p.y()];
        auto v = cell.voronoi_status;
        if (v != FREE_QUEUED && v != VORONOI_RETRY) { // || v>free || v==VORONOI_PRUNE || v==VORONOI_KEEP) {
            //      assert(v!=RETRY);
            continue;
        }

        MarkerMatchResult r = MarkerMatch(p.x(), p.y());
        if (r == PRUNED) cell.voronoi_status = VORONOI_PRUNE;
        else if (r == KEEP) cell.voronoi_status = VORONOI_KEEP;
        else { // r==RETRY
            cell.voronoi_status = VORONOI_RETRY;
            //      printf("RETRY %d %d\n", x, sizeY-1-y);
            prune_queue_.push(p);
        }
        data_[p.x()][p.y()] = cell;

        if (sorted_prune_queue_.empty()) {
            while (!prune_queue_.empty()) {
                Eigen::Vector2i p = prune_queue_.front();
                prune_queue_.pop();
                sorted_prune_queue_.push(data_[p.x()][p.y()].sqrt_dist, p);
            }
        }
    }
}

void DynamicVoronoi::UpdateAlternativePrunedDiagram() {

    if (alternative_diagram_ == nullptr) {
        alternative_diagram_ = new int *[size_x];
        for (int x = 0; x < size_x; x++) {
            alternative_diagram_[x] = new int[size_y];
        }
    }


    std::queue<Eigen::Vector2i> end_cells;
    BucketPrioQueue<Eigen::Vector2i> sortedPruneQueue;
    for (int x = 1; x < size_x - 1; x++) {
        for (int y = 1; y < size_y - 1; y++) {
            Cell &c = data_[x][y];
            alternative_diagram_[x][y] = c.voronoi_status;
            if (c.voronoi_status <= FREE) {
                sortedPruneQueue.push(c.sqrt_dist, Eigen::Vector2i(x, y));
                end_cells.push(Eigen::Vector2i(x, y));
            }
        }
    }

    for (int x = 1; x < size_x - 1; x++) {
        for (int y = 1; y < size_y - 1; y++) {
            if (GetNumVoronoiNeighborsAlternative(x, y) >= 3) {
                alternative_diagram_[x][y] = VORONOI_KEEP;
                sortedPruneQueue.push(data_[x][y].sqrt_dist, Eigen::Vector2i(x, y));
                end_cells.push(Eigen::Vector2i(x, y));
            }
        }
    }

    for (int x = 1; x < size_x - 1; x++) {
        for (int y = 1; y < size_y - 1; y++) {
            if (GetNumVoronoiNeighborsAlternative(x, y) >= 3) {
                alternative_diagram_[x][y] = VORONOI_KEEP;
                sortedPruneQueue.push(data_[x][y].sqrt_dist, Eigen::Vector2i(x, y));
                end_cells.push(Eigen::Vector2i(x, y));
            }
        }
    }


    while (!sortedPruneQueue.empty()) {
        Eigen::Vector2i p = sortedPruneQueue.pop();

        if (MarkerMatchAlternative(p.x(), p.y())) {
            alternative_diagram_[p.x()][p.y()] = VORONOI_PRUNE;
        } else {
            alternative_diagram_[p.x()][p.y()] = VORONOI_KEEP;
        }
    }

    // //delete worms
    while (!end_cells.empty()) {
        Eigen::Vector2i p = end_cells.front();
        end_cells.pop();

        if (IsVoronoiAlternative(p.x(), p.y()) && GetNumVoronoiNeighborsAlternative(p.x(), p.y()) == 1) {
            alternative_diagram_[p.x()][p.y()] = VORONOI_PRUNE;

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    if (!(dx || dy) || (dx && dy)) {
                        continue;
                    }
                    int nx = p.x() + dx;
                    int ny = p.y() + dy;
                    if (nx < 0 || nx >= size_x || ny < 0 || ny >= size_y) {
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
    bool f[8];

    int nx, ny;
    int dx, dy;

    int i = 0;
    int voro_count = 0;
    for (dy = 1; dy >= -1; dy--) {
        ny = y + dy;
        for (dx = -1; dx <= 1; dx++) {
            if (dx || dy) {
                nx = x + dx;
                int v = alternative_diagram_[nx][ny];
                bool b = (v <= FREE && v != VORONOI_PRUNE);
                f[i] = b;
                if (v <= FREE && !(dx && dy))
                    voro_count++;
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
        if (voro_count == 1 && (f[1] || f[3] || f[4] || f[6])) {
            return false;
        }

        // 4-connected
        if ((!f[0] && f[1] && f[3]) ||
            (!f[2] && f[1] && f[4]) ||
            (!f[5] && f[3] && f[6]) ||
            (!f[7] && f[6] && f[4]))
            return false;

        if ((f[3] && f[4] && !f[1] && !f[6]) ||
            (f[1] && f[6] && !f[3] && !f[4]))
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
            if (nx < 0 || nx >= size_x || ny < 0 || ny >= size_y) {
                continue;
            }
            if (alternative_diagram_[nx][ny] == FREE ||
                alternative_diagram_[nx][ny] == VORONOI_KEEP) {
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

    int i = 0;
    int count = 0;
    int voro_count = 0;
    int voro_count_four = 0;

    for (dy = 1; dy >= -1; dy--) {
        ny = y + dy;
        for (dx = -1; dx <= 1; dx++) {
            if (dx || dy) {
                nx = x + dx;
                Cell neighbor_cell = data_[nx][ny];
                auto v = neighbor_cell.voronoi_status;
                bool b = (v <= FREE && v != VORONOI_PRUNE);
                //	if (v==occupied) obstacleCount++;
                f[i] = b;
                if (b) {
                    voro_count++;
                    if (!(dx && dy)) {
                        voro_count_four++;
                    }
                }
                if (b && !(dx && dy)) count++;
                //	if (v<=free && !(dx && dy)) voro_count++;
                i++;
            }
        }
    }
    if (voro_count < 3 && voro_count_four == 1 && (f[1] || f[3] || f[4] || f[6])) {
        //    assert(voro_count<2);
        //    if (voro_count>=2) printf("voro>2 %d %d\n", x, y);
        return KEEP;
    }

    // 4-connected
    if ((!f[0] && f[1] && f[3]) ||
        (!f[2] && f[1] && f[4]) ||
        (!f[5] && f[3] && f[6]) ||
        (!f[7] && f[6] && f[4]))
        return KEEP;
    if ((f[3] && f[4] && !f[1] && !f[6]) ||
        (f[1] && f[6] && !f[3] && !f[4]))
        return KEEP;



    // KEEP voro cells inside of blocks and RETRY later
    if (voro_count >= 5 && voro_count_four >= 3 &&
        data_[x][y].voronoi_status != VORONOI_RETRY) {
        return RETRY;
    }

    return PRUNED;
}
