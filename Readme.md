 # 3d galaxy collision simulation
 
**exmaple**

~40 fps, 20k particles

![Screen Recording 2026-02-17 at 02.30.06.gif](contents/Screen%20Recording%202026-02-17%20at%2002.30.06.gif)

# run with your own galaxy settings

```bash

git clone --recurse-submodules https://github.com/tinunadno/3d_galaxy_collision_simulation.git
cd 3d_galaxy_collision_simulator
# probably change galaxy setup in src/main.cpp
mkdir build && cd build
cmake .. && cmake --build .
./galaxy_sim
```