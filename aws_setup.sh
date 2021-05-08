sudo apt install make
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc
sudo apt install gfortran
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install libfftw3-dev
mkdir vtkdata
sudo apt install libblacs-mpi-dev

export OMP_NUM_THREADS=1
export OMPI_MCA_btl=^openib
