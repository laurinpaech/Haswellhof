# Haswellhof
Advanced Systems Lab Project 2020

Team 42: Carla Jancik, Sebastian Winberg, Valentin Wolf, Laurin Paech

## Getting started

To run this code do the following steps:

- Open a terminal and navigate to the project.

- Create a new directory `build` in the root project directory:

```
mkdir build
```

- Navigate to this newly created directory `build`:

```
cd build
```

- Run [CMake](https://cmake.org/download/) to create project structure:

```
cmake .. -DCMAKE_BUILD_TYPE=Release
```

Note that with CMake is compatible with many IDEs and can e.g. create an Xcode project by adding the flag `-G Xcode`. For more information have a look at the CMake Documentation.

- Next, you should call the Makefile and compile the program:

```
make
```

- If all went smoothly and everything compiled and linked correctly, you can now execute the program `surf` the following way:

```
./src/surf
```

- Note that everytime you change something in the project structure, you need to call `cmake ..` again. If something isn't working correctly it is sometimes also helpful to delete the current `build` folder and redo all the above steps again.

## Running Benchmarking

- Instead of executing the program, one can also run benchmarking in directory `build`:

```
./src/benchmark
```

- This will run benchmarking with the configuration defined in `src/benchmark/main_benchmark.cpp` and create csv timings of every function in `benchmarking_files`.

## Running Validation

- Instead of executing the program, one can also validate functions:

```
./src/validation IMAGE_PATH desc1
```

- This will run validation with the configuration defined in `src/validation/main_validation.cpp` and will compare the outputs to the baseline.

## Plotting Keypoints and Matchings

For visualising the found keypoints of a single image the `visualize_keypoints.py` Python script is provided in the `/scripts/` folder. First change the `img_fp` the `dest` variables at the top of the script as desired. Then run:

```
cd scripts/
python visualize_keypoints.py
```

To show feature matches between two images use the `visualize_matches.py` script that is also included in the `/scripts/` folder. Fist, change `img1_fp`,  `img2_fp` and 
`dest` variables as desired and run:

```
python visualize_matches.py
```

For matching the OpenCV2 function `cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)` is used.

Note, that by default results from both our base implementation and the openCV SIFT implementation will be created and shown in both scripts. To use a different version of our implementation (i.e. optimized) an other executable can be specified in the scripts. Results from OpenSURF can also be visualized but keypoint files have to be first generated using OpenSURF and then specified in the script. 

