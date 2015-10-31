from rsf.proj import *

septi2d = '../src/septi2d.x'

# ========================================= #
class sepTIop2d:
    '''separable qP-wave propagation'''
    def __init__(self, v, ss, rr, par, custom=''):
        self.v = v
        self.custom = custom
        self.par = 'verb=%(verb)s nb=%(nb)d snap=%(snap)s '%par+custom
        
        self.dep = [self.v]
        self.ss = ''
        if (ss != ''):
            self.ss = ' sou=' + ss + '.rsf '
            self.dep.append(ss)
        
        self.rr = ''
        if (rr != ''):
            self.rr = ' rec=' + rr + '.rsf '
            self.dep.append(rr)

    # ------------------------------------- #
    def FORW(self, m, d):
        Flow(d, [m]+self.dep,
            septi2d 
            + ''' adj=n model=${SOURCES[1]} ''' 
            + self.ss + self.rr + self.par)

    # ------------------------------------- #
    def ADJT(self, m, d):
        Flow(m, [d]+self.dep,
            septi2d 
            + ''' adj=y model=${SOURCES[1]} ''' 
            + self.ss + self.rr + self.par)

