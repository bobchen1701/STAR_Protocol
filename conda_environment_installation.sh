#!/bin/sh
conda create -n qc_pipe python=3.8
conda activate qc_pipe
conda install numpy pandas matplotlib scipy scikit-learn networkx
conda install -c conda-forge python-igraph leidenalg louvain numba pytables scikit-image pot
pip install --user scanpy
pip install --user ipykernel
python -m ipykernel install --user --name=qc_pipe
