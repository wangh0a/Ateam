from rsf.proj import *

sepvti2d = '../src/sepvti2d.x'
septti2d = '../src/septti2d.x'

# ========================================= #
class sepTIop2d:
    '''separable qP-wave propagation'''
    def __init__(self, v, ss, rr, par, custom=''):
        self.v = v
        self.custom = custom
        self.par = 'verb=%(verb)s nb=%(nb)d snap=%(snap)s '%par+custom
        self.septi2d = sepvti2d
        if par['atype'] == 'tti':
            self.septi2d = septti2d
        
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
            self.septi2d 
            + ''' adj=n model=${SOURCES[1]} ''' 
            + self.ss + self.rr + self.par)

    # ------------------------------------- #
    def ADJT(self, m, d):
        Flow(m, [d]+self.dep,
            self.septi2d 
            + ''' adj=y model=${SOURCES[1]} ''' 
            + self.ss + self.rr + self.par)

