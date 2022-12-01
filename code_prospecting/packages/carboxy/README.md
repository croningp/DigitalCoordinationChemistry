# Carboxylate Exploration Package, ```carboxy```

Supporting code for the algorithmic exploration of chemical seach spaces defined by a series of multinuclear, mixed-ligand metal acetate complexes.

This work accompanies the paper "Hybrid Human-Robot Laboratories enable Library Generation and Discovery in Coordination Chemistry", by Kowalski, MacGregor, Long, Bell, and Cronin.

No further development is planned at this time.

## Dependencies
Required packages and versions. May require their own dependencies.
* ```python``` -> 3.7.1
* ```matplotlib``` -> 3.4.1
* ```numpy``` -> 1.21.5
* ```pandas``` -> 1.2.4
* ```scikit-learn``` -> 0.24.2
* ```scikit-optimize``` -> 0.9.0
* ```summit``` -> 0.8.4

## Package Structure

```python
carboxy         # DIGITAL DATA INTERPRETATION / CONDITION GENERATION
│   amf.py          # calculation of analytical mapping function from MS data
│   explore.py      # generation of samples for subsequent iterations
│   msrep.py        # reproducibility calculation though comparison of MS data
│   __init__.py 
│
└───matrix          # VISUALISATION
       explore_vis.py   # matrix plots for input conditions
                        # plots of AMF score vs experiment number
```

## Support

Please contact the authors at lee.cronin@glasgow.ac.uk for support or further information