#include <random>

#include "main_pipeline.h"
#include "physical_particle.h"
#include "black_body_color.h"

inline sc::utils::Vec<float, 3> circularOrbitVelocity(
    const sc::utils::Vec<float, 3>& position,
    float centralMass,
    float G = 1.f
) {
    float R = std::sqrt(position[0] * position[0] +
                        position[1] * position[1]);

    if (R < 1e-4f)
        return sc::utils::Vec<float, 3>{0.f, 0.f, 0.f};

    float speed = std::sqrt(G * centralMass / R);

    return sc::utils::Vec<float, 3>{
        -position[1] / R * speed,
         position[0] / R * speed,
         0.f
    };
}

inline void generateParticles(
    PhysicalParticleSystem<float>& particles,
    std::size_t count,
    float baseTemperature = 6000.f,
    const sc::utils::Vec<float, 3> offset = {},
    const sc::utils::Vec<float, 3> rotation = {},
    const sc::utils::Vec<float, 3> speedOffset = {}
) {
    std::size_t oldSize = particles.speed.size();
    std::size_t newSize = count + oldSize;

    particles.ps.positions.resize(newSize);
    particles.ps.colors.resize(newSize);
    particles.ps.sizes.resize(newSize);
    particles.ps.enableRender.resize(newSize);
    particles.speed.resize(newSize);
    particles.weight.resize(newSize);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> uniform01(0.f, 1.f);
    std::normal_distribution<float> normal(0.f, 1.f);

    constexpr float galaxyRadius = 2.5f;
    constexpr float diskHeight = 0.25f;
    constexpr float centralMass = 50.f;

    constexpr float armTightness = 3.5f;
    constexpr int armCount = 4;
    constexpr float armWidth = 0.25f;

    constexpr float eccentricityNoise = 0.15f;
    constexpr double temperatureSigma = 3000.0;

    for (std::size_t i = oldSize; i < newSize; ++i) {
        float r = galaxyRadius * std::sqrt(uniform01(gen));

        float armIndex = std::floor(uniform01(gen) * armCount);
        float armPhase = armIndex * (2.f * static_cast<float>(M_PI) / armCount);

        float theta = armTightness * std::log(r + 0.1f) + armPhase;

        theta += normal(gen) * armWidth;

        float x = r * std::cos(theta);
        float y = r * std::sin(theta);
        float z = normal(gen) * diskHeight;

        sc::utils::Vec<float, 3> basePos{x, y, z};

        auto circVel = -circularOrbitVelocity(basePos, centralMass);

        float R = std::sqrt(x*x + y*y);
        sc::utils::Vec<float, 3> radial{
            x / (R + 1e-6f),
            y / (R + 1e-6f),
            0.f
        };

        circVel += radial * (normal(gen) * eccentricityNoise);
        circVel[2] += normal(gen) * 0.05f;

        particles.ps.positions[i] =
            sc::utils::rotateEuler(basePos, rotation) + offset;

        particles.speed[i] =
            sc::utils::rotateEuler(circVel, rotation) + speedOffset;

        particles.weight[i] = 0.00001f;

        particles.ps.sizes[i] =
            0.005f + 0.005f * uniform01(gen);

        std::normal_distribution<double> tempNoise(0.0, temperatureSigma);
        double kelvin = baseTemperature + tempNoise(gen);
        kelvin = std::clamp(kelvin, 2500.0, 20000.0);

        particles.ps.colors[i] = bbc(kelvin);

        particles.ps.enableRender[i] = true;
    }

    particles.weight[oldSize] = centralMass;
    particles.ps.positions[oldSize] = offset;
    particles.ps.sizes[oldSize] = 0.12f;
    particles.speed[oldSize] = speedOffset;
    particles.ps.colors[oldSize] = sc::utils::Vec<float, 3>{1.f, .5f, .5f};
}


int main() {
    sc::Camera<float, sc::VecArray> camera;
    camera.pos()[2] = 2.0f;
    camera.setLen(0.3);
    camera.setRes(sc::utils::Vec<float, 2>{1000, 800});

    PhysicalParticleSystem<float> particles;
    mrc::makeCircleBillboard(particles.ps, 8);

    generateParticles(particles, 10000, 9000.f
        , sc::utils::Vec<float, 3>{-5.f, 5.f, 0.f}
        , sc::utils::Vec<float, 3>{-1.f, 5.f, 0.f}
        , sc::utils::Vec<float, 3>{.1f, -.1f, .5f}
        );
    generateParticles(particles, 10000, 4000.f
        , sc::utils::Vec<float, 3>{5.f, 0.f, -5.f}
        , sc::utils::Vec<float, 3>{1.f, 0.f, -2.f}
        , sc::utils::Vec<float, 3>{-.1f, .1f, -.5f}
        );

    float step = 0.001f;
    bool pause = false;

    auto efmu = [&particles, &step, &pause](std::size_t, std::size_t) {
        if (!pause) iterateMT(particles, step);
    };

    std::vector<std::pair<std::vector<int>, std::function<void()>>> customKeyHandlers = {
        {{GLFW_KEY_Q}, [&step](){ step *= 1.1f; },},
        {{GLFW_KEY_E}, [&step](){ step *= .9f; },},
        {{GLFW_KEY_LEFT_ALT, GLFW_KEY_P}, [&pause](){ pause = !pause; },},
    };

    mrc::initMrcRender(camera, { }, { },
        efmu, { }, customKeyHandlers,
        sc::utils::Vec<int, 2>{-1,-1}, 60, {}, &particles.ps);

    return 0;
}
