#pragma once

//---------------------------------ChatGPT----------------------------------------

#include "rapidobj/rapidobj.hpp"

#include <stdexcept>

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
    int height,
    int resolution)
{
    if (resolution <= 0.0)
        throw std::invalid_argument("Resolution must be positive.");

    const auto& positions = obj.attributes.positions;

    if (positions.empty())
        return {};

    // Find the XY bounds of the complete model.
    double minX = positions[0], maxX = positions[0];
    double minY = positions[1], maxY = positions[1];

    for (size_t i = 0; i + 2 < positions.size(); i += 3)
    {
        minX = std::min(minX, (double)positions[i]);
        maxX = std::max(maxX, (double)positions[i]);
        minY = std::min(minY, (double)positions[i + 1]);
        maxY = std::max(maxY, (double)positions[i + 1]);
    }

    const size_t nx =
        static_cast<size_t>(std::ceil((maxX - minX) / resolution)) + 1;
    const size_t ny =
        static_cast<size_t>(std::ceil((maxY - minY) / resolution)) + 1;

    ObjectMask mask(nx, std::vector<bool>(ny, false));
    std::vector<Segment> segments;

    auto getVertex = [&](int index)
        {
            const size_t i = 3 * static_cast<size_t>(index);
            return std::array<double, 3>{
                positions[i], positions[i + 1], positions[i + 2]
            };
        };

    constexpr double eps = 1e-7;

    // Intersect one triangle edge with z = height.
    auto intersectEdge = [&](const auto& a, const auto& b,
        std::vector<Point2D>& points)
        {
            const double da = a[2] - height;
            const double db = b[2] - height;

            if (std::abs(da) <= eps)
                points.push_back({ a[0], a[1] });

            if ((da < -eps && db > eps) ||
                (da > eps && db < -eps))
            {
                const double t = da / (da - db);

                points.push_back({
                    a[0] + t * (b[0] - a[0]),
                    a[1] + t * (b[1] - a[1])
                    });
            }
        };

    // Extract the intersection segments from every triangle.
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

            const auto a = getVertex(
                mesh.indices[offset].position_index);
            const auto b = getVertex(
                mesh.indices[offset + 1].position_index);
            const auto c = getVertex(
                mesh.indices[offset + 2].position_index);

            offset += faceSize;

            const double da = a[2] - height;
            const double db = b[2] - height;
            const double dc = c[2] - height;

            // Skip triangles coplanar with the cutting plane.
            if (std::abs(da) <= eps &&
                std::abs(db) <= eps &&
                std::abs(dc) <= eps)
                continue;

            std::vector<Point2D> points;
            points.reserve(3);

            intersectEdge(a, b, points);
            intersectEdge(b, c, points);
            intersectEdge(c, a, points);

            // Remove duplicate intersection points.
            std::vector<Point2D> unique;
            for (const auto& p : points)
            {
                bool duplicate = false;

                for (const auto& q : unique)
                {
                    if (std::hypot(p.x - q.x, p.y - q.y) <= eps)
                    {
                        duplicate = true;
                        break;
                    }
                }

                if (!duplicate)
                    unique.push_back(p);
            }

            if (unique.size() == 2 &&
                std::hypot(unique[0].x - unique[1].x,
                    unique[0].y - unique[1].y) > eps)
            {
                segments.push_back({ unique[0], unique[1] });
            }
        }
    }

    // Classify a point using the even-odd ray-crossing rule.
    auto inside = [&](double x, double y)
        {
            bool result = false;

            for (const auto& s : segments)
            {
                const auto& a = s.a;
                const auto& b = s.b;

                if ((a.y > y) != (b.y > y))
                {
                    const double xIntersect =
                        a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y);

                    if (x < xIntersect)
                        result = !result;
                }
            }

            return result;
        };

    // Sample the cross-section at every grid point.
    for (size_t ix = 0; ix < nx; ++ix)
    {
        for (size_t iy = 0; iy < ny; ++iy)
        {
            const double x = minX + ix * resolution;
            const double y = minY + iy * resolution;

            mask[ix][iy] = inside(x, y);
        }
    }

    return mask;
}
