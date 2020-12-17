# STAR_Protocol

### Install Dependencies
The environment for Step 1 of this protocol require Singularity. The details of this installation can be found in the [Step 1 notebook](https://github.com/KenLauLab/STAR_Protocol/blob/master/Step_1_Bioinformatics_Pipeline.ipynb).


The environment for Steps 2 and 3 of this protocol can be installed through the included .yml with the following, which will also initialize a jupyter kernel: 
```
conda env create -f qc_pipe_env.yml
python -m ipykernel install --user --name=qc_pipe
```

Otherwise, the requirements.txt file can be manually installed In a working Python virtual environment using `pip`:
```
pip install -r requirements.txt
```
