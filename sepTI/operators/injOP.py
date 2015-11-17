#
# Inject operator
#  
try:    from rsf.cluster import *
except: from rsf.proj    import *

INJ2d = '../src/INJOP2D.x'
#INJ2d = 'injop2d'
#INJ3d = 'injop3d'
# ------------------------------------------------------------
class injop2d:
  def __init__(self,c,par,custom=''):
    self.c=c
    self.par=par.copy() # copy the par file instead of referencing it
    self.custom=custom

  def FORW(self,m,d):
    Flow(d,[m,self.c],
         INJ2d + '''
          adj=n
         coo=${SOURCES[1]}
         nz=%(nz)d oz=%(oz)g dz=%(dz)g
         nx=%(nx)d ox=%(ox)g dx=%(dx)g
         '''%self.par+' '+self.custom)

  def ADJT(self,m,d):
    Flow(m,[d,self.c],
         INJ2d + '''
          injop2d adj=y
         coo=${SOURCES[1]} 
         '''+' '+self.custom)

# ------------------------------------------------------------
