/*----------------------------------------------------------------*/
/*--                                                            --*/
/*-- Title:    FSIVGTT  data.                                   --*/
/*-- Filename: ivgtt_example.sas                                --*/
/*--           Invokes minmod.inc and fsivgtt_dt.sas            --*/
/*-- Date:     August 1, 2002                                   --*/
/*--                                                            --*/
/*-- Ref:                                                       --*/
/*--                                                            --*/
/*--                                                            --*/
/*----------------------------------------------------------------*/


libname lib '.';

%let path=../..;   
filename ivgtt_dt "&path/datasets/fsivgtt_dt.sas";
filename nlinmix  "&path/nlmm800.sas";
filename nlmem    "&path/nlmem.sas";
filename minmod   "minmod.inc";
%include nlinmix / nosource;
%include ivgtt_dt;


%macro useriml;
/* User defined SAS/IML functions */
%include minmod;          /* Minimal model definition */

start bsm(b,u) global(der_bu,age);
/* Between-subject model */
/* Derivs for: p1,p2,p3,G0 
          wrt b and u */

p1      = b[1] + u[1];
logp2   = b[2] + u[2];
logsi   = b[3] + u[3] + (age-70)*b[5]/100; /* Linear model for log(si) */
G0      = b[4];

db      = i(4) || j(4,1,0);
db[3,5] = (age-70)/100;

du      = {1 0 0,
           0 1 0,
           0 0 1,
           0 0 0};


logp3   = logp2 + logsi;
dlogp3b = { 1  1}* db[2:3,];
dlogp3u = { 1  1}* du[2:3,];

p2      = 10**(logp2);
db[2,2] = p2* db[2,2]*log(10);
du[2,2] = p2* du[2,2]*log(10);

p3      = 10**(logp3);
db[3,]  = p3* dlogp3b*log(10);
du[3,]  = p3* dlogp3u*log(10);

der_bu  = db || du;
return(p1//p2//p3//G0);
finish bsm;

start dbsm(b,u) global(der_bu);;
return(der_bu);
finish;

start wghtc(dummy) global(g_t,time, predv,weight);

/* Vector of weights */
g_weight= choose(time>=8,1,0);  /* Zero-weighted measures for 
                                         time<8 mins */
g_weight[6]=50;                 /* Overweighted measure for time=8 */
w= g_weight/(g_t#g_t);
return(w);
finish;


%mend useriml;

%nlinmix(data = fsivgtt,
         next = id,
         xvar = age,
         stmts   =  %str( 
         class id;
         model pseudo_g_t = d_beta1 d_beta2 d_beta3  d_beta4 d_beta5
                               / noint notest solution cl;
         random d_u1 d_u2 d_u3 / subject=id solution cl type=un; 
         weight _weight_;
         ods output Iterhistory=iter;
                      ),
         globvar = predv dpredv dd  g_b i_b  x_solv ins_action,
         model=%str(
           g_b=g_t[1];  /* Baseline glucose  */
           i_b=i_t[1];  /* Baseline insulin */

           /* Between subject model */
           _d=bsm(_b,_u);
           der_bu=dbsm(_b,_u);

           /* Within-subject model */
           *scale = {0.01 0.01 0.00001 100}; 
           scale = {0.01 1 1 100}; 
           d=_d # t(scale);
           g_solv = glc_solv(d);   /* Solving ODE       */
           predv  = g_solv[,1];    /* Predicted glucose */
           x_solv = g_solv[,2];    /* Predicted glucose in a remote compartment */
           jacs   = g_solv[,3:6];  /* Solution for sensitivity eqs */
           jacs   = jacs # scale; 
           dpredv = jacs * der_bu;
           ins_action = 100*x_solv/(d[1]+x_solv);       /* Insulin action calculated */
),
       derivs  =  %str(
                  _db=dpredv[,1:5]; /* Derivatives wrt _b extracted */ 
                  _du=dpredv[,6:8]; /* Derivatives wrt _u extracted */ 
                  ),
       savevar =  x_solv ins_action,
       weight  =  %str( _weight_= wghtc(.);),
       parms   =  %str(beta1=2  beta2=-2  beta3=-3  beta4=2.5 beta5=0),
       expand  =  zero,
       imlfun  =  useriml,         /* Points to a SAS macro with a user defined SAS/IML functions */
       converge=  1e-3,
       options = skipnlin,         /* Skips initial call to a nonlinear regression */
       nlmem   =  nlmem /nosource  /* Points to NLMEM macro source code */
);
