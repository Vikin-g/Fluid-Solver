
//---------------------------------ChatGPT----------------------------------------
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <rapidobj/rapidobj.hpp>


using ObjectMask = std::vector<std::vector<bool>>;

struct Point2D
{
    double x, y;
};

struct Segment
{
    Point2D a, b;
};

inline ObjectMask FindCrossSection(
    const rapidobj::Result& obj,
    const Bounds3D& bounds,
    double z,
    int nx,
    int ny,
    double cell)
{

    if (nx <= 0 || ny <= 0 || cell <= 0.0)
        throw std::invalid_argument("nx, ny and cell must be positive.");

    ObjectMask mask(nx, std::vector<bool>(ny, false));

    const auto& P = obj.attributes.positions;
    if (P.empty())
        return mask;

    auto vertex = [&](int index)
        {
            const size_t i = 3 * static_cast<size_t>(index);
            return std::array<double, 3>{ P[i], P[i + 1], P[i + 2] };
        };

    // An edge crosses the plane when its endpoints are on different sides.
    // A vertex exactly on the plane counts as "above" (z >= plane), so there is
    // no special case for touching vertices or edges lying in the plane, and
    // shared edges always give matching points in neighbouring triangles.
    auto crossing = [&](const std::array<double, 3>& p,
        const std::array<double, 3>& q,
        Point2D& out)
        {
            const bool pAbove = p[2] >= z;
            const bool qAbove = q[2] >= z;
            if (pAbove == qAbove)
                return false;

            const double t = (z - p[2]) / (q[2] - p[2]);
            out = { p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]) };
            return true;
        };

    // 1. Collect the line segments where triangles meet the plane.
    std::vector<Segment> segments;

    for (const auto& shape : obj.shapes)
    {
        const auto& mesh = shape.mesh;
        size_t offset = 0;

        for (unsigned char faceSize : mesh.num_face_vertices)
        {
            if (faceSize != 3)
            {
                offset += faceSize;
                continue;
            }

            const auto a = vertex(mesh.indices[offset].position_index);
            const auto b = vertex(mesh.indices[offset + 1].position_index);
            const auto c = vertex(mesh.indices[offset + 2].position_index);
            offset += 3;

            Point2D pts[3];
            int n = 0;
            if (crossing(a, b, pts[n])) ++n;
            if (crossing(b, c, pts[n])) ++n;
            if (crossing(c, a, pts[n])) ++n;

            if (n == 2)
                segments.push_back({ pts[0], pts[1] });
        }
    }

    if (segments.empty())
        return mask;

    // 2. Scanline fill (even-odd): for each row of cell centres, find where the
    //    outline crosses that row, sort, and fill between pairs of crossings.
    std::vector<double> xs;

    for (int iy = 0; iy < ny; ++iy)
    {
        const double y = bounds.miny + (iy + 0.5) * cell;
        xs.clear();

        for (const auto& s : segments)
        {
            if ((s.a.y > y) != (s.b.y > y))
                xs.push_back(s.a.x + (y - s.a.y) * (s.b.x - s.a.x) / (s.b.y - s.a.y));
        }

        std::sort(xs.begin(), xs.end());

        for (size_t k = 0; k + 1 < xs.size(); k += 2)
        {
            // first / one-past-last cell whose centre lies in [xs[k], xs[k+1])
            const int i0 = std::max(0, (int)std::ceil((xs[k] - bounds.minx) / cell - 0.5));
            const int i1 = std::min(nx, (int)std::ceil((xs[k + 1] - bounds.minx) / cell - 0.5));

            for (int ix = i0; ix < i1; ++ix)
                mask[ix][iy] = true;
        }
    }

    return mask;
}

ObjectMask ResampleMask(const ObjectMask& src, int nx, int ny, float threshold = 0.5f)
{
    ObjectMask dst(nx, std::vector<bool>(ny, false));
    if (src.empty() || src[0].empty() || nx <= 0 || ny <= 0) return dst;

    const int sw = (int)src.size();
    const int sh = (int)src[0].size();

    // Each destination cell covers this region of the source, in source cells
    const float sx = float(sw) / nx;
    const float sy = float(sh) / ny;

    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            int x0 = (int)std::floor(x * sx);
            int x1 = std::max(x0 + 1, (int)std::ceil((x + 1) * sx));   // at least one source cell
            int y0 = (int)std::floor(y * sy);
            int y1 = std::max(y0 + 1, (int)std::ceil((y + 1) * sy));

            x1 = std::min(x1, sw);
            y1 = std::min(y1, sh);

            int filled = 0, total = 0;
            for (int i = x0; i < x1; i++)
                for (int j = y0; j < y1; j++)
                {
                    filled += src[i][j];
                    total++;
                }

            dst[x][y] = total > 0 && float(filled) / total >= threshold;
        }
    }
    return dst;
}

//----------------------------Claude-------------------------------------------


// ---------------------------------------------------------------------------
// 1. Bounding box of the whole mesh (model units). Compute ONCE, not per slice.
// ---------------------------------------------------------------------------
struct Bounds3D
{
    float minx = FLT_MAX, maxx = -FLT_MAX;
    float miny = FLT_MAX, maxy = -FLT_MAX;
    float minz = FLT_MAX, maxz = -FLT_MAX;

    float width()  const { return maxx - minx; }
    float height() const { return maxy - miny; }
    float depth()  const { return maxz - minz; }
};

inline Bounds3D ComputeBounds(const rapidobj::Result& result)
{
    const auto& P = result.attributes.positions;      // x,y,z,x,y,z,...
    if (P.size() < 3)
        throw std::runtime_error("ComputeBounds: mesh has no vertices");

    Bounds3D b;
    for (size_t i = 0; i + 2 < P.size(); i += 3)
    {
        b.minx = std::min(b.minx, P[i]);     b.maxx = std::max(b.maxx, P[i]);
        b.miny = std::min(b.miny, P[i + 1]); b.maxy = std::max(b.maxy, P[i + 1]);
        b.minz = std::min(b.minz, P[i + 2]); b.maxz = std::max(b.maxz, P[i + 2]);
    }
    return b;
}

// ---------------------------------------------------------------------------
// 2. Metres per model unit, from a known real size along one axis.
//    real_width_m is the real size along model X (e.g. 0.177).
// ---------------------------------------------------------------------------
inline float ComputeScale(const Bounds3D& b, float real_width_m)
{
    if (b.width() <= 0.0f)
        throw std::runtime_error("ComputeScale: zero width along X");
    return real_width_m / b.width();
}

// Sanity check: the other axis should match its known real size.
// Example: ScaleMatches(b, to_m, 0.039f) -> true if Y really is ~39 mm.
inline bool ScaleMatches(const Bounds3D& b, float to_m, float expected_height_m, float tol = 0.05f)
{
    float h = b.height() * to_m;
    return std::fabs(h - expected_height_m) <= tol * expected_height_m;
}

// ---------------------------------------------------------------------------
// 3. Grid size of the object, fixed for every slice.
// ---------------------------------------------------------------------------
struct GridDims { int nx; int ny; };

inline GridDims ComputeObjectGrid(const Bounds3D& b, float to_m, float cell_size_m)
{
    GridDims g;
    g.nx = std::max(1, (int)std::lround(b.width() * to_m / cell_size_m));
    g.ny = std::max(1, (int)std::lround(b.height() * to_m / cell_size_m));
    return g;
}

// Number of cells for a real-world length (use for the whole domain).
inline int CellsForLength(float length_m, float cell_size_m)
{
    return std::max(1, (int)std::lround(length_m / cell_size_m));
}

// ---------------------------------------------------------------------------
// 4. Mapping between grid cells and model coordinates (same for every slice).
//    Gives the model-space position of the CENTRE of cell (i, j).
// ---------------------------------------------------------------------------
inline void CellToModel(const Bounds3D& b, float to_m, float cell_size_m,
    int i, int j, float& model_x, float& model_y)
{
    float cell_model = cell_size_m / to_m;            // cell size in model units
    model_x = b.minx + (i + 0.5f) * cell_model;
    model_y = b.miny + (j + 0.5f) * cell_model;
}

// ---------------------------------------------------------------------------
// 5. Slicing along Z.
// ---------------------------------------------------------------------------
inline int ComputeSliceCount(const Bounds3D& b, float to_m, float slice_spacing_m)
{
    return std::max(1, (int)std::ceil(b.depth() * to_m / slice_spacing_m));
}

// Model-space z of the middle of slice k.
inline float SliceZ(const Bounds3D& b, float to_m, float slice_spacing_m, int k)
{
    return b.minz + (k + 0.5f) * slice_spacing_m / to_m;
}

// ---------------------------------------------------------------------------
// 6. Placement: object x position (cells from inlet) and vertical centring.
// ---------------------------------------------------------------------------
inline void ComputePlacement(int resy, const GridDims& obj, float inlet_distance_m,
    float cell_size_m, int& posx, int& posy)
{
    posx = (int)std::lround(inlet_distance_m / cell_size_m);
    posy = (resy - obj.ny) / 2;
}

// Check the object fits inside the domain at that placement.
inline bool ObjectFits(int resx, int resy, const GridDims& obj, int posx, int posy)
{
    return posx >= 0 && posy >= 0 &&
        posx + obj.nx <= resx - 1 &&               // keep clear of the outlet column
        posy + obj.ny <= resy - 1;
}

// ---------------------------------------------------------------------------
// 7. Mask helpers (only needed if FindCrossSection returns a bitmap
//    that you then have to resize).
// ---------------------------------------------------------------------------
inline bool AnyTrue(const ObjectMask& m)
{
    for (const auto& col : m)
        for (bool v : col)
            if (v) return true;
    return false;
}

// Number of filled cells in column x -> useful to skip slices that are too thin.
inline int MaxThickness(const ObjectMask& m)
{
    int best = 0;
    for (const auto& col : m)
    {
        int c = 0;
        for (bool v : col) c += v;
        best = std::max(best, c);
    }
    return best;
}

// Resample a mask onto nx x ny cells (grow or shrink).
// Pass the WHOLE uncropped bounding-box bitmap, same size for every slice.
inline ObjectMask ResampleMask(const ObjectMask& src, int nx, int ny, float threshold = 0.5f)
{
    ObjectMask dst(nx, std::vector<bool>(ny, false));
    if (src.empty() || src[0].empty() || nx <= 0 || ny <= 0) return dst;

    const int sw = (int)src.size();
    const int sh = (int)src[0].size();
    const float sx = float(sw) / nx;
    const float sy = float(sh) / ny;

    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
        {
            int x0 = (int)std::floor(x * sx);
            int x1 = std::min(sw, std::max(x0 + 1, (int)std::ceil((x + 1) * sx)));
            int y0 = (int)std::floor(y * sy);
            int y1 = std::min(sh, std::max(y0 + 1, (int)std::ceil((y + 1) * sy)));

            int filled = 0, total = 0;
            for (int i = x0; i < x1; i++)
                for (int j = y0; j < y1; j++)
                {
                    filled += src[i][j];
                    total++;
                }
            dst[x][y] = total > 0 && float(filled) / total >= threshold;
        }
    return dst;
}
