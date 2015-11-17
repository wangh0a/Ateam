/* inject/extract in/from 2D wavefield */
#include <rsf.h>
#include "fdutil.h"

int main(int argc, char* argv[])
{
    bool verb,adj,expl;
    int method;
    
    /* I/O files */
    sf_file Ftrc=NULL; /* traces       */
    sf_file Fcoo=NULL; /* coordinates  */
    sf_file Fwfl=NULL; /* wavefield    */
    
    /* cube axes */
    sf_axis at,az,ax,aa,ac;
    
    /* I/O arrays */
    float  *wco  =NULL;  /* traces   */
    pt2d   *coo  =NULL;  /* coordinates   */
    /* scoef2 *csinc=NULL;  [> weights/indices <] */
    scoef2d csinc=NULL;  /* weights/indices */
    lint2d  cow;         /* weights/indices */
    float **wfl  =NULL;  /* wavefield   */
    float **exwfl=NULL;  /* extended wavefield   */
    float  *sum  =NULL;  /* stacked trace for explosive sources adjoint*/
    
    fdm2d fdm=NULL;
    int   iz,ix,it,i;
    int   nz,nx;
    float dz,dx;
    float oz,ox;
    float ssum=0;
    
    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    
    /*------------------------------------------------------------*/
    if(! sf_getint("method",&method)) method=1 ;
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("adj", &adj))   adj=false; /* adjoint flag */
    if(! sf_getbool("expl",&expl)) expl=false; /* Multiple sources, one wvlt */
    
    /*------------------------------------------------------------*/
    /* setup I/O */
    Fcoo = sf_input("coo"); /* coordinates */
    ac = sf_iaxa(Fcoo,2); sf_setlabel(ac,"c"); sf_setunit(ac,"");
    coo = (pt2d*) sf_alloc(sf_n(ac),sizeof(*coo));
    pt2dread1(Fcoo,coo,sf_n(ac),2); /* read (x,z) coordinates */
    
    if(adj) {
        Fwfl = sf_input ("in");  /* wavefield */
        Ftrc = sf_output("out"); /* traces   */
        
        az = sf_iaxa(Fwfl,1); sf_setlabel(az,"z");
        ax = sf_iaxa(Fwfl,2); sf_setlabel(ax,"x");
        at = sf_iaxa(Fwfl,3); sf_setlabel(at,"t");
        
        aa = sf_maxa(1,0,1);
        if(expl) sf_oaxa(Ftrc,aa,1);
        else sf_oaxa(Ftrc,ac,1);
        sf_oaxa(Ftrc,at,2);
        sf_oaxa(Ftrc,aa,3);
        
    } else {
        Ftrc = sf_input ("in" ); /* traces   */
        Fwfl = sf_output("out"); /* wavefield */
        
        at = sf_iaxa(Ftrc,2); sf_setlabel(at,"t");
        
        if(!sf_getint  ("nz",&nz)) nz=1;
        if(!sf_getfloat("oz",&oz)) oz=0.0;
        if(!sf_getfloat("dz",&dz)) dz=1.0;
        az = sf_maxa(nz,oz,dz);
        sf_setlabel(az,"z");
        
        if(!sf_getint  ("nx",&nx)) nx=1;
        if(!sf_getfloat("ox",&ox)) ox=0.0;
        if(!sf_getfloat("dx",&dx)) dx=1.0;
        ax = sf_maxa(nx,ox,dx);
        sf_setlabel(ax,"x");
        
        sf_oaxa(Fwfl,az,1);
        sf_oaxa(Fwfl,ax,2);
        sf_oaxa(Fwfl,at,3);
    }
    
    if(verb) {
        sf_raxa(az);
        sf_raxa(ax);
        sf_raxa(at);
        sf_raxa(ac);
    }
    
    /* allocate wavefield arrays */
    if(expl & !adj) wco = sf_floatalloc(       1);
    else            wco = sf_floatalloc(sf_n(ac));
    wfl = sf_floatalloc2(sf_n(az),sf_n(ax));
    
    if(expl & adj) sum = sf_floatalloc( 1);
    
    /* interpolation coefficients */
    fdm = fdutil_init(verb,'n',az,ax,4,1);
    
    exwfl =  sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    if(verb) sf_warning("nxpad=%d nzpad=%d",fdm->nzpad,fdm->nxpad);
    
    if(method!=1){
    csinc = sinc2d_make(sf_n(ac),coo,fdm);
        
    if(verb) sf_warning("1st trace position %d,%d",csinc[0].ix,csinc[0].iz);
        
    } else cow = lint2d_make(sf_n(ac),coo,fdm);
    
    if(verb) sf_warning("\nmethod = %d ",method);


    /*------------------------------------------------------------*/
    /*
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<sf_n(at); it++) {
        if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);
        
        // Reading and zeroing arrays
        if(adj){
            sf_floatread(wfl[0],sf_n(az)*sf_n(ax),Fwfl);
            
            if(expl) {wco[0]=0.0f;}
            else { for(i=0; i<sf_n(ac); i++) wco[i]=0.0f;}
            
        } else {
            for    (ix=0; ix<fdm->nx; ix++)
                for(iz=0; iz<fdm->nz; iz++)
                    wfl[ix][iz]=0.0f;
            
            if(expl) sf_floatread(wco,1,Ftrc);
            else sf_floatread(wco,sf_n(ac),Ftrc);
        }
        
        expand(wfl,exwfl,fdm);
        
        // Interpolation
        if(method==1){ // linear injection
            if(adj)
                lint2d_extract(exwfl,wco,cow);

            else {
                if(expl)lint2d_inject1(exwfl,wco[0],cow);
                else    lint2d_inject (exwfl,wco   ,cow);
                
            }
        } else {  // kaiser-sinc injection (Hicks, 2002)
            if(adj)
                if(expl) {sum[0]=0.0f; sinc2d_extract1(exwfl,sum,csinc);}
                else     sinc2d_extract(exwfl,wco,csinc);
                
            else {
                if(expl) sinc2d_inject1(exwfl,wco[0],csinc);
                else     sinc2d_inject (exwfl,wco   ,csinc);
            }
        }
        
        //Writing
        if(adj){
            if(expl){
                if(method==1){
                sum[0]=0.0f;
                for (i=0; i<sf_n(ac); i++)
                    sum[0] += wco[i];}
                sf_floatwrite(sum,1,Ftrc);
            }
            else sf_floatwrite(wco,sf_n(ac),Ftrc);
        } else{
            cut2d(exwfl,wfl,fdm,az,ax);
            sf_floatwrite(wfl[0],sf_n(az)*sf_n(ax),Fwfl);
        }

    }	/* end time loop */
    if(verb) fprintf(stderr,"\n");
    free(csinc);
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*exwfl); free(exwfl);
    free(*wfl); free(wfl);
    free(wco);
    free(sum);
    
    /*------------------------------------------------------------*/ 
    /* close files */
    if (Ftrc!=NULL) sf_fileclose(Ftrc); 
    if (Fwfl!=NULL) sf_fileclose(Fwfl);
    if (Fcoo!=NULL) sf_fileclose(Fcoo);
    
    exit (0);
}
