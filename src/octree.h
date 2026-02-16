#pragma once

#include <vector>
#include <cmath>
#include <algorithm>

class BarnesHutTree {
public:
    void build(const float* px, const float* py, const float* pz,
               const float* masses, std::size_t count,
               float theta = 0.7f, float softening = 0.05f)
    {
        px_ = px; py_ = py; pz_ = pz;
        masses_ = masses;
        n_ = count;
        theta2_ = theta * theta;
        eps2_ = softening * softening;
        nodeCount_ = 0;

        if (n_ == 0) return;

        std::size_t estimated = n_ * 3 + 256;
        if (nodes_.size() < estimated) nodes_.resize(estimated);

        float minX = px[0], maxX = px[0];
        float minY = py[0], maxY = py[0];
        float minZ = pz[0], maxZ = pz[0];
        for (std::size_t i = 1; i < n_; i++) {
            if (px[i] < minX) minX = px[i]; if (px[i] > maxX) maxX = px[i];
            if (py[i] < minY) minY = py[i]; if (py[i] > maxY) maxY = py[i];
            if (pz[i] < minZ) minZ = pz[i]; if (pz[i] > maxZ) maxZ = pz[i];
        }

        float cx = (minX + maxX) * 0.5f;
        float cy = (minY + maxY) * 0.5f;
        float cz = (minZ + maxZ) * 0.5f;
        float hs = std::max({maxX - minX, maxY - minY, maxZ - minZ}) * 0.5f + 0.1f;

        allocNode(cx, cy, cz, hs);

        for (std::size_t i = 0; i < n_; i++) {
            insert(0, static_cast<int>(i), 0);
        }

        computeCOM(0);
    }

    void computeAccel(std::size_t i, float& ax, float& ay, float& az) const {
        ax = ay = az = 0;
        if (nodeCount_ == 0) return;

        const float bx = px_[i], by = py_[i], bz = pz_[i];
        const int skipBody = static_cast<int>(i);

        int stack[256];
        int top = 0;
        stack[top++] = 0;

        while (top > 0) {
            const int ni = stack[--top];
            const Node& nd = nodes_[ni];

            if (!nd.isInternal) {
                if (nd.body >= 0 && nd.body != skipBody) {
                    float dx = nd.comX - bx;
                    float dy = nd.comY - by;
                    float dz = nd.comZ - bz;
                    float d2 = dx * dx + dy * dy + dz * dz + eps2_;
                    float inv = 1.0f / std::sqrt(d2);
                    float f = nd.mass * inv * inv * inv;
                    ax += dx * f;
                    ay += dy * f;
                    az += dz * f;
                }
                continue;
            }

            float dx = nd.comX - bx;
            float dy = nd.comY - by;
            float dz = nd.comZ - bz;
            float dist2 = dx * dx + dy * dy + dz * dz;
            float size = nd.halfSize * 2.0f;

            if (size * size < theta2_ * dist2) {
                float d2 = dist2 + eps2_;
                float inv = 1.0f / std::sqrt(d2);
                float f = nd.mass * inv * inv * inv;
                ax += dx * f;
                ay += dy * f;
                az += dz * f;
            } else {
                for (int c = 0; c < 8; c++) {
                    if (nd.children[c] >= 0) {
                        stack[top++] = nd.children[c];
                    }
                }
            }
        }
    }

    std::size_t nodeCount() const { return nodeCount_; }

private:
    static constexpr int MAX_DEPTH = 30;

    struct Node {
        float comX = 0, comY = 0, comZ = 0;
        float mass = 0;
        float cx = 0, cy = 0, cz = 0;
        float halfSize = 0;
        int children[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
        int body = -1;
        bool isInternal = false;
    };

    std::vector<Node> nodes_;
    std::size_t nodeCount_ = 0;

    const float* px_ = nullptr;
    const float* py_ = nullptr;
    const float* pz_ = nullptr;
    const float* masses_ = nullptr;
    std::size_t n_ = 0;
    float theta2_ = 0.49f;
    float eps2_ = 0.0025f;

    int allocNode(float cx, float cy, float cz, float hs) {
        int idx = static_cast<int>(nodeCount_++);
        if (nodeCount_ > nodes_.size()) {
            nodes_.resize(nodes_.size() * 2 + 256);
        }
        Node& nd = nodes_[idx];
        nd.cx = cx; nd.cy = cy; nd.cz = cz;
        nd.halfSize = hs;
        for (int i = 0; i < 8; i++) nd.children[i] = -1;
        nd.body = -1;
        nd.isInternal = false;
        nd.mass = 0;
        nd.comX = nd.comY = nd.comZ = 0;
        return idx;
    }

    static int octant(float cx, float cy, float cz,
                      float px, float py, float pz) {
        return ((px > cx) ? 1 : 0) | ((py > cy) ? 2 : 0) | ((pz > cz) ? 4 : 0);
    }

    int getOrMakeChild(int pi, int oct) {
        if (nodes_[pi].children[oct] == -1) {
            float hs = nodes_[pi].halfSize * 0.5f;
            float ncx = nodes_[pi].cx + ((oct & 1) ? hs : -hs);
            float ncy = nodes_[pi].cy + ((oct & 2) ? hs : -hs);
            float ncz = nodes_[pi].cz + ((oct & 4) ? hs : -hs);
            nodes_[pi].children[oct] = allocNode(ncx, ncy, ncz, hs);
        }
        return nodes_[pi].children[oct];
    }

    void insert(int ni, int bi, int depth) {
        if (!nodes_[ni].isInternal && nodes_[ni].body == -1) {
            nodes_[ni].body = bi;
            return;
        }

        if (depth >= MAX_DEPTH) return;

        if (!nodes_[ni].isInternal) {
            int existing = nodes_[ni].body;
            nodes_[ni].body = -1;
            nodes_[ni].isInternal = true;

            int o = octant(nodes_[ni].cx, nodes_[ni].cy, nodes_[ni].cz,
                           px_[existing], py_[existing], pz_[existing]);
            int child = getOrMakeChild(ni, o);
            insert(child, existing, depth + 1);
        }

        int o = octant(nodes_[ni].cx, nodes_[ni].cy, nodes_[ni].cz,
                       px_[bi], py_[bi], pz_[bi]);
        int child = getOrMakeChild(ni, o);
        insert(child, bi, depth + 1);
    }

    void computeCOM(int ni) {
        Node& nd = nodes_[ni];

        if (!nd.isInternal) {
            if (nd.body >= 0) {
                nd.comX = px_[nd.body];
                nd.comY = py_[nd.body];
                nd.comZ = pz_[nd.body];
                nd.mass = masses_[nd.body];
            }
            return;
        }

        double mx = 0, my = 0, mz = 0, totalMass = 0;
        for (int i = 0; i < 8; i++) {
            int c = nd.children[i];
            if (c >= 0) {
                computeCOM(c);
                double cm = nodes_[c].mass;
                totalMass += cm;
                mx += nodes_[c].comX * cm;
                my += nodes_[c].comY * cm;
                mz += nodes_[c].comZ * cm;
            }
        }

        nd.mass = static_cast<float>(totalMass);
        if (totalMass > 0) {
            double inv = 1.0 / totalMass;
            nd.comX = static_cast<float>(mx * inv);
            nd.comY = static_cast<float>(my * inv);
            nd.comZ = static_cast<float>(mz * inv);
        }
    }
};
