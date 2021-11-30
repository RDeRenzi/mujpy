from numpy import sqrt, pi, exp
import numpy as np
from scipy.special import erf

class muprompt(object):
    def _init_(self,x,y,e=1):
        '''
        generic __init__ for all mufunction classes
        x, y, e are numpy arrays
        e is always defined, 
        either provided by the caller 
        or default to np.ones(x.shape[0])
        '''
        if x.shape[0]!=y.shape[0]:
            raise ValueError('x, y have different lengths')
        else:
            self.x = x
            self.y = y
        try: 
            if e==1:
                self.e = np.ones(x.shape[0])
        except:
            if e.shape[0]!=x.shape[0]:
                raise ValueError('x, e have different lengths')           
            else:
                self.e = e
    # ---- end generic __init__
    def f(self,x,a,x0,dx,ak1,ak2): 
        '''
        fit function for PSI t0 fit on the prompt peak, 
        sluggish python version
        a gaussian peak coincident with the edge betwee two plateaus (a constant + an erf)
        parameters: 
        peak height, 
        peak center, 
        peak std width, 
        first plateau, 
        second plateau.
        Needs x as first ture argument because it is also invoked for plot.
        '''
        return a/(sqrt(2.*pi)*dx) * exp(-.5*((x-x0)/dx)**2)+ak2/2.+ak1+ak2/2.*erf((x-x0)/sqrt(2.)/dx)
        
    def __call__(self,a, x0, dx, ak1, ak2):
        '''
        chisquare, normalized if self.e was provided, unnormalized otherwise 
        '''
        return sum(((self.f(self.x,a,x0,dx,ak1,ak2)-self.y)/self.e)**2)

