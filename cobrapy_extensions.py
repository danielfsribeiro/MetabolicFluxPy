### Extension to COBRApy 
### Daniel Ribeiro 2020
from cobra.core.gene import Gene

# Add __lt__ method to cobra.core.gene
# This method is required in numpy.ndarray.sort()...
# Called in libsmop.py #def setdiff(a, b) 
#Without it there is an error
"""
TypeError: '<' not supported between instances of 'Gene' and 'Gene'
"""

def __lt__(self, other): # Have to add self since this will become a method
    return self.name < other.name
#Register attibute
setattr(Gene, '__lt__', __lt__)
