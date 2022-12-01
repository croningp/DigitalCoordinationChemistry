# Algorithmic Exploration and Novelty Quantification

Contains scripts for the machine-learning led exploration of a search space of metal carboxylate mixed-ligand complexes. Covers the full data science workflow from digital data interpretation, calculation of a mapping function to sampling of new conditions based on random sampling from a uiform distribution, Latin Hypercube Sampling, or Gaussian Proess Bayesian Optimisation. Raw mass spectrometry data from the exploration is also included to permit replication. 

This work accompanies the paper Kowalski, MacGregor, Long, Bell, and Cronin, *Automated Library Generation and Serendipity Quantification enables Diverse Discovery in Coordination Chemistry*, Journal of the American Chemical Society, 2022.

No further development is planned at this time.

## Dependencies

Required packages and versions. May require their own dependencies.
* ```python``` -> 3.7.1
* ```ipython``` -> 7.2.0
* ```matplotlib``` -> 3.4.1
* ```nmrglue``` -> 0.8
* ```numpy``` -> 1.21.5
* ```pandas``` -> 1.2.4
* ```scikit-learn``` -> 0.24.2
* ```scikit-optimize``` -> 0.9.0
* ```summit``` -> 0.8.4

## Included Packages

The source files are included herein and the package can be installed by: (i) downloading the source directories from Git, (ii) on the user's computer, navigating to the directory containing the ```setup.py``` file in Command Prompt, (iii) installing with pip by running the command ```pip install -e .``` in this directory. Packages are required to run much of the exploration code.

### Carboxylate Exploration Package, ```carboxy```

Supporting code for the algorithmic exploration of chemical seach spaces defined by a series of multinuclear, mixed-ligand metal acetate complexes.

### Spectra Visualiser Package, ```spectra_visualiser```

This is a bespoke package used to import analytical data to python and handle mundane visualisation needs. 

## Support

Please contact the authors at lee.cronin@glasgow.ac.uk for support or further information
