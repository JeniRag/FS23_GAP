# FS23_GAP
My workflow of my Master's Thesis "Machine Learning the Equation of State of Dense Hydrogen using Gaussian Approximation Potential".
Still in Progress.

* notebooks contain jupyter notebook folders
	1_PrepareData.ipynb: example to extract frames and store in.xyz file. (example data used in notebook not included)
	2_ReducedGradientModel.ipynb: example to train and test GAP model.
* intermediate_data contains a demo file used in 2_ReducedGradientModel.ipynb


# How to install librascal on Science Cluster
Create Conda environment on Cluster
```
module load anaconda3
conda create --name gap python=3.10
source activate gap
```
(Here I created an environment called gap)
```
git clone https://github.com/lab-cosmo/librascal.git
cd librascal
pip install -r requirements.txt
```
Make sure to build the packages, as described in the [librascal instructions](https://lab-cosmo.github.io/librascal/installation.html).
 
To run simulations using your GAP potential, you need to install i-Pi:
```
pip install git+https://github.com/lab-cosmo/i-pi.git
```

After setting up your environment create kernel, which can be used in you jupyter notebooks.

```
ipython kernel install --user --name gap
```
Deactivate environment if you don't need it anymore.

```
conda deactivate gap
```


