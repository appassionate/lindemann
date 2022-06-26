# lindemann
lindemann index based on MDAnalysis
## There will be some features.
1. run in parallel: using multiprocessing, the LindemannAna can run successfully in parallel.
2. less memory using: running variance and running average is used to avoid value overflow.


## example
```python
from lindemann.parallel import parallel_exec
from lindemann.lindemann import LindemannAna
from MDAnalysis import Universe

#initialize Universe in mda and load trajectory
u = Universe("./traj.xyz") 
#init lindemann analysis class
lindemann = LindemannAna(u.atoms,# point out the  MDA Universe 
                         verbose=True,# show the analysis progress bar
                          cell=None, #set the trajectory cell
                          )

#parallel execute lindemann analysis
parallel_exec(lindemann.run, )

```