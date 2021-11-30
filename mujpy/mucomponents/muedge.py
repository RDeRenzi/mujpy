from numpy import sqrt, pi, exp
import numpy as np

class muedge(object):
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

    def f(self,t, t00, N, D):
        '''
        fit function for ISIS t0 fit on edge 
        see mulab function f=mathedge(t,N,Delta,TauPi,TauMu,t0)
        obtained from MuonBeamSimple.nb in Mathematica
        here time in musec (vs mulab in ns)
        '''
        tp = 0.026003 # mus    
        t0 = 0.82*tp-t00
        tm = 2.197 # mus 
        f = (6*N/(D**3*(tm-tp)))*\
             exp(-(((t+t0+D/2)*tm+(t+t0+D/2)*tp)/(tm*\
                                               tp)))*\
             ((-2*exp(((t+t0+D/2)*tm+D*tp)/(tm*tp))*\
             tm**2*(-D/2+tm)+\
             2*exp((D*tm+(t+t0+D/2)*tp)/(tm*tp))*\
             tp**2*(-(D/2)+tp)+\
             exp(((t+t0+D/2)*(tm+tp))/(tm*tp))*(tm-tp)*\
             (t**2+t0**2+2*t0*(t-tm-tp)-2*t*(tm+tp)+\
             2*(-(D**2/8)+tm**2+tm*tp+tp**2)))*\
             ((t+t0-D/2)>0)+\
             (2*exp((t+t0+D/2)/tp)*tm**2*(D/2+tm)-\
              2*exp((t+t0+D/2)/tm)*tp**2*(D/2+tp)-\
              exp(((t+t0+D/2)*(tm+tp))/(tm*tp))*(tm-tp)*\
              (t**2+t0**2+2*t0*(t-tm-tp)-2*t*(tm+tp)+\
              2*(-(D**2/8)+tm**2+tm*tp+tp**2)))*((t+t0+D/2)>0))
        return f
        
    def __call__(self,t00, N, D):
        '''
        chisquare, normalized if self.e was provided, unnormalized otherwise 
        '''
        return sum(((self.f(self.x, t00, N, D)-self.y)/self.e)**2)

