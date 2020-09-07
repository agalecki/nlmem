/* MINMOD: GLUCOSE DISAPPEARANCE MODEL                  --*/
/* SAS/IML code for NLINMIX: A. Galecki, June 2002      --*/

/*--------------------------------------------------------*/
/*--                                                    --*/
/*--    Glucose disappearance model                     --*/
/*--                                                    --*/
/*--                                                    --*/
/*-- glcd      - integrand function for MINMOD model    --*/
/*               and sensitivity eqs                    --*/
/*-- glc_solv  - solves glucose disappearance model     --*/
/*-- ws_odef   - invokes glc_solv                       --*/
/*-- jacn      - returns jacobian (numerical)           --*/
/*-- wghtc     - weight for observations                --*/
/*-- interpolate - interpolation. Auxiliary function    --*/

/*--------------------------------------------------------*/


start glcd(t,y)  global(dd,I_b,G_b,I_t,time);
/* Integrand for glucose disappearance model */
/* time is a vector extracted from data */
/* t is a scalar */
/* function <interpolate> interpolates insulin level at time t */

    d=dd;

    p1=d[1]; /* S_G  Glucose sensitivity*/
    p2=d[2]; /* p_2 */
    p3=d[3]; /* p_3 */
    G0=d[4]; /* G_0  Glucose at time 0*/

    g_t =y[1];
    x_t =y[2];

    
    idiff= interpolate(t,time,I_t) - I_b;

    /* Prevents idiff to be less than zero */
     idiff= choose(idiff>0,idiff,0);    
     *print 'glcd' t p1 p2 p3 G0 I_b G_b idiff;
     dg_dt= -(p1+x_t)*g_t + p1*G_b;
     dx_dt= -p2*x_t       + p3*idiff;
     *print 'glcd' t g_t x_t dg_dt dx_dt;

/* Integrand function for sensitivity equations */
    g_1=y[3];
    g_2=y[4];     x_2=y[7];
    g_3=y[5];     x_3=y[8];
    g_4=y[6];

    
    c=p1+x_t;  
    dg_1=-c*g_1 - (g_t - G_b);
    dg_2=-c*g_2 -  g_t*x_2;
    dg_3=-c*g_3 -  g_t*x_3;
    dg_4=-c*g_4;
    
    dx_2=-p2*x_2-x_t;
    dx_3=-p2*x_3+idiff;

    return(dg_dt//dx_dt //
        dg_1//dg_2//dg_3//dg_4//
        dx_2//dx_3
);
finish glcd;


start glc_solv(d) global(dd,I_t,g_t,time,id);
/* Function returns solution vector for G_t and X_t  and sensitivity eqs.*/
/* Invokes GLCD function */

   dd=d;

  ***h= {1E-12  1. 1E-5};            /* hmin, hmax, hinit step size*/
   h1=10**(-12);
   h2=1;
   h3=10**(-5);
   h=h1||h2||h3;

  c= d[4] // {0, 0,0,0,1,0,0};       /* Initial condition y(0)=c */

  /* Initial condition y(0)=c */
  id1=id[1];
  *print 'Before ODE'  id1 c d time;
  eps=10**(-4);
  call ode(res,"glcd",c,time,h) eps=eps;  
  glc_solv=c`//res`;
  return(glc_solv);
finish glc_solv;


start interpolate(x0,x,y);

/* This is an auxiliary function to perform linear interpolation  */
/* Vectors x,y represent coordinates of several points            */
/* on the plane.                                                  */
/* Elements of vector x are in ascending order                    */
/* Scalar x0 indicates at which point on x-axis to interpolate    */
/* Function returns y0 with interpolated value at x0       */

i2=(x>=x0)[<:>];
i1=max(i2-1,1);
num=x0-x[i1];
den=x[i2]-x[i1];
r=0;
if den>10**(-9) then r=num/den;
y0=y[i1, ]+ r*(y[i2, ]-y[i1, ]);
*print 'int' x0 i1 i2 num den r y0;

return(y0);
finish interpolate;


