#pragma once

#include "particle_system.h"
#include "basic_thread_pool.h"
#include "octree.h"

template<typename NumericT>
struct PhysicalParticleSystem {
    mrc::ParticleSystem<NumericT> ps;
    std::vector<sc::utils::Vec<NumericT, 3>> speed;
    std::vector<NumericT> weight;
};

template<typename NumericT>
void iterate(PhysicalParticleSystem<NumericT> &particles, NumericT step, NumericT epsilon = NumericT(0.05)) {
    const std::size_t N = particles.speed.size();
    if (N == 0) return;

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i + 1; j < N; ++j) {
            auto diff = particles.ps.positions[j] - particles.ps.positions[i];

            NumericT dist2 = sc::utils::dot(diff, diff) + epsilon * epsilon;
            NumericT invDist = NumericT(1) / std::sqrt(dist2);
            NumericT invDist3 = invDist * invDist * invDist;

            auto forceDir = diff * invDist3;
            particles.speed[i] += forceDir * (particles.weight[j] * step);
            particles.speed[j] += forceDir * (-particles.weight[i] * step);
        }
    }
    for (std::size_t i = 0; i < N; ++i) {
        particles.ps.positions[i] += particles.speed[i] * step;
    }
}

template<typename NumericT>
void iterateMT(PhysicalParticleSystem<NumericT> &particles, NumericT step,
               NumericT epsilon = NumericT(0.05), NumericT theta = NumericT(0.7)) {
    static ThreadPool pool;
    static BarnesHutTree tree;
    static std::vector<float> px, py, pz;

    const std::size_t N = particles.speed.size();
    if (N == 0) return;

    px.resize(N); py.resize(N); pz.resize(N);
    for (std::size_t i = 0; i < N; i++) {
        px[i] = static_cast<float>(particles.ps.positions[i][0]);
        py[i] = static_cast<float>(particles.ps.positions[i][1]);
        pz[i] = static_cast<float>(particles.ps.positions[i][2]);
    }

    tree.build(px.data(), py.data(), pz.data(),
               particles.weight.data(), N,
               static_cast<float>(theta), static_cast<float>(epsilon));

    pool.runConcurrentTask([&](std::size_t tid, std::size_t tc) {
        const std::size_t chunk = (N + tc - 1) / tc;
        const std::size_t start = std::min(chunk * tid, N);
        const std::size_t end = std::min(start + chunk, N);

        for (std::size_t i = start; i < end; i++) {
            float ax, ay, az;
            tree.computeAccel(i, ax, ay, az);

            particles.speed[i][0] += static_cast<NumericT>(ax) * step;
            particles.speed[i][1] += static_cast<NumericT>(ay) * step;
            particles.speed[i][2] += static_cast<NumericT>(az) * step;

            particles.ps.positions[i][0] += particles.speed[i][0] * step;
            particles.ps.positions[i][1] += particles.speed[i][1] * step;
            particles.ps.positions[i][2] += particles.speed[i][2] * step;
        }
    });
}
