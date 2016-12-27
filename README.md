# OptoMC-heterogeneous
Mesh-based (heterogeneous) Monte Carlo Simulation with Direct Photon Flux Recording Strategy for Optogenetics.

To run the program below you need an additional iso2mesh MATLAB toolbox for meshing and visualization (http://iso2mesh.sourceforge.net).

This version includes a rodent brain model and a four-layer tissue model, and the heterogeneous models are made of mesh.
The four-layer  tissue model is intended to compare our model to the results of MCML.

To run the code for the four-layered tissue model, please excute:  'optomc_fourLayer.m'

To run the code for the rodent brain model, please excute:  'optomc_rat.m'

To visualize precomputed results for the rodent brain model, please excute:  'visualization.m'


See the article below for more details.
Shin, Younghoon, and Hyuk-Sang Kwon. "Mesh-based Monte Carlo method for fibre-optic optogenetic neural stimulation with direct photon flux recording strategy." Physics in medicine and biology 61.6 (2016): 2265.
DOI:  10.1088/0031-9155/61/6/2265

See also the following journal articles.
Shin, Younghoon, et al. "Characterization of fiber-optic light delivery and light-induced temperature changes in a rodent brain for precise optogenetic neuromodulation." Biomedical Optics Express 7.11 (2016): 4450-4471.
DOI:  10.1364/BOE.7.004450
