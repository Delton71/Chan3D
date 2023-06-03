# Chan3D
Implementation of the 3-d Divide-and-Conquer Convex Hull

## How to use:
* Install [GoogleTest](https://google.github.io/googletest/):
```
sudo apt-get install libgtest-dev libgmock-dev libtbb-dev cmake
cd /usr/src/googletest/
sudo mkdir build
cd build
sudo cmake ..
sudo make
sudo cp lib/* /usr/lib
```
Or use locally using github this way (don't forget to change include paths in .cpp files):
```
git clone https://github.com/google/googletest.git
```
* Clone this repository;
* Create tests (or use my example):
```
mkdir build && cd build
cmake ..
make
./Chan3D
```
