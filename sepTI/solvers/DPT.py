#
# Dot product 
#    
try:    from rsf.cluster import *
except: from rsf.proj    import *

class DPT:
    def __init__(self,G,m,d):
        self.G=G
        self.m=m
        self.d=d

    # ------------------------------------------------------------
    def DOT(self,d,i,o):
        Flow(d,[i,o],
             '''
             add mode=p ${SOURCES[1]} |
             stack axis=0 norm=n |
             put o1=0 d1=1 o2=0 d2=1 o3=0 d3=1
             ''')

    # ------------------------------------------------------------
    def TEST(self,name=''):
        
        prfx = '_DPT'+self.m+self.d 

        mi = prfx + self.m + 'i' # input
        mo = prfx + self.m + 'o' # output
        md = prfx + self.m + 'd' # dot product
        
        di = prfx + self.d + 'i' # input
        do = prfx + self.d + 'o' # output
        dd = prfx + self.d + 'd' # dot product

        Flow(mi,self.m,'noise rep=y seed=10101') # random m
        Flow(di,self.d,'noise rep=y seed=01010') # random d

        self.G.FORW(mi,do) # do = G [mi]
        self.G.ADJT(mo,di) # mo = G*[di]

        self.DOT(md,mi,mo) # output name is md
        self.DOT(dd,di,do) # output name is dd

        # report
        Flow( prfx ,[md,dd],
            '''
            math m=${SOURCES[0]} d=${SOURCES[1]} output="abs(d-m)"|
            cat axis=1 space=n ${SOURCES[1]} ${SOURCES[0]} |
            scale axis=123 |
            reverse which=1 |
            disfil format="DPTEST %9.6f" col=1
            ''',stdin=0,stdout=-1)


