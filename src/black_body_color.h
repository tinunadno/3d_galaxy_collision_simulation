#pragma once

#include "utils/vec.h"

inline sc::utils::Vec<float, 3> bbc(double kelvin)
{
    kelvin = std::clamp(kelvin, 1000.0, 40000.0);
    kelvin /= 100.0;
    double red, green, blue;
    if (kelvin <= 66.0)
    {
        red = 255.0;
    }
    else
    {
        red = kelvin - 60.0;
        red = 329.698727446 * std::pow(red, -0.1332047592);
        red = std::clamp(red, 0.0, 255.0);
    }
    if (kelvin <= 66.0)
    {
        green = kelvin;
        green = 99.4708025861 * std::log(green) - 161.1195681661;
    }
    else
    {
        green = kelvin - 60.0;
        green = 288.1221695283 * std::pow(green, -0.0755148492);
    }
    green = std::clamp(green, 0.0, 255.0);

    if (kelvin >= 66.0)
    {
        blue = 255.0;
    }
    else if (kelvin <= 19.0)
    {
        blue = 0.0;
    }
    else
    {
        blue = kelvin - 10.0;
        blue = 138.5177312231 * std::log(blue) - 305.0447927307;
        blue = std::clamp(blue, 0.0, 255.0);
    }

    return sc::utils::Vec<float, 3>{
        static_cast<float>(red / 255.0),
        static_cast<float>(green / 255.0),
        static_cast<float>(blue / 255.0)
    };
}