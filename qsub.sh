#!/bin/bash

#*** PBS setting when needed
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -j oe
#PBS -m ea
#PBS -V

PBS_O_WORKDIR=/work/a06/yingying/camada/HydroDA/src
cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate data
python run.py

