/*
INDOMETHACIN DATA.
TWO COMPARTMENTS MODEL.
ODE SPECIFICATION. NUMERICAL SOLUTION.
See "Population PK/PD models using NLINMIX" manuscript

RANDOM    U1 U2
NEXT      =11        
SKIPNLIN  no
DEBUG     blank
APPEND    blank
OPTIONS   ???
CPU time:  secs
*/

options nocenter ls=78;

%let path=../..;
filename indomet  "&path/datasets/indomet_dt.sas";
filename nlinmix  "&path/nlmm800.sas";
filename nlmem    "&path/nlmem.sas";
%include nlinmix / nosource;

%include indomet /nosource;
%inc nlinmix / nosource;
  

%macro useriml;
/* ---- Within subject model ----- */ 

start comp2_integrand(t,q) global(a);
/* Integrand function for two compartments model */
Dq=a*q;
*print t q dq;
return(Dq);
finish comp2_integrand;

start f(d) global(time, a);
 q0  = d[1];  
 k01 = d[2];
 k21 = d[3];
 k12 = d[4];

 a      = j(2,2,.);               /* Compartmental matrix */ 
 a[1,1] = -(k01+k21);
 a[1,2] = k12;
 a[2,1] = k21;
 a[2,2] = -k12;

                       /* Solve ODE */
 c     = q0 // 0;              /* Initial values: q1(0),q2(0) */
 h     = {1e-12 1. 1e-5};      /* Minimum allowable step size,
                                  Maximum allowable step size,
                                  Initial step size to start integration */

 t12  = 0 || t(time);          /* Time intervals for integration */
 call ode(q,"comp2_integrand", c, t12, h) eps=1E-12;
 q1 = q[1,];                   /* Extract q_1 */
 return(t(q1));
finish f;
  
                           /* Between subject model */ 
start d(b,u) global(du);
   du={1 0,
       0 1,
       0 0, 
       0 0};
    d = b + du*u;
 return(d);
finish d;

                           /* Derivatives for between subject model */
start dd(b,u) global(du);
 dd= i(4) || du;
 return(dd);
finish dd;  
%mend useriml;
   

     
Title "ODE: Numerical solution";
   
                           /* NLINMIX invocation */
%nlinmix(data=indomet, 
   next  =11,
   globvar= dpredv,
  
   model=%str(

      b=t(_b);
      u=t(_u);
      d       =    d(b,u);   /* */
      dd      =   dd(b,u);
      par     =   11 // {1,.};
      call nlpfdd(q1, dq1,h, "f", d) par=par;
      predv = q1;  
      dpredv     = dq1*dd;                /* Chain rule */
   ),
      
  derivs=%str(
      _db = dpredv[,1:4];
      _du = dpredv[,5:6]; 
      /* print _db; */
      ),
   parms=%str(q0=3 k01=2 k21=2 k12=2),
   stmts=%str(
      class subject;
      model pseudo_conc = d_q0 d_k01 d_k21 d_k12/ noint notest solution cl;
      random d_u0 d_u01 / subject=subject solution cl;
   ),
   expand   =  eblup,
   options  =  /* debug     skipnlin   notes*/,
   imlfun   =  useriml,
   converge =  2e-5,
   nlmem    =  nlmem
);

