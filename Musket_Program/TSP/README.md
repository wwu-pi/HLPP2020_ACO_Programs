# Run Program

To locally the programm the build-and-run.sh script can be used. It uses the CMake Files included. However, this requires to adjust the paths of the src/aco_start.cu file in line 591-602. Please consider that the build-and-run.sh script includes also very big cities which you might want to exclude when they are run locally. In the current configuration the script runs more than 16 hours.

# Structure
The .musket file is the Musket code. **Include**, **src**, and **CMake Files** are generated. Please keep in mind that reading the files is not included in the Musket code but done manually. 
