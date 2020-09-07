/*
ORANGE TREE DATA. 
LOGISTIC GROWTH MODEL. 
ODE SPECIFICATION

parms      = b1=190 b2=10 b3=-.003

RANDOM     U1
NEXT       =7
SKIPNLIN   no
DEBUG      yes
USER_DEBUG yes
DERIVS     nlpfdd
APPEND     no
CONVERGE   2e-5

CPU time:  25 secs
*/

/*---Examples for NLINMIX macro---*/

/*---Data on orange trees, Draper and Smith, 1981, p. 524, and
     Lindstrom and Bates, Biometrics, 1990, 673-687---*/

%let path=../..;
filename tree_dt "&path/datasets/tree_dt.sas";
filename nlinmix "&path/nlmm800.sas";
filename nlmem   "&path/nlmem.sas";
%include tree_dt;
%inc  nlinmix / nosource;

%macro user_debug;
if _isubjct=1 then do; 
   print _b _u; 
   print predv _dbu; 
end;
%mend  user_debug;

%macro useriml;

start bsm(bu);
/* Between subject model */
/* Extract _b and _u */
   _b=bu[1:3];
   _u=bu[4];
   bsm = i(3)*_b + {1,0,0}*_u;
 return(bsm);
finish bsm;

start dbsm(bu);
/* Derivatives for between subject model */
 dbsm= i(3) || {1,0,0}; 
 return(dbsm);
finish dbsm;

start der_ode(x,y) global(dd);
/*-- Logistic model. ODE specification  --*/
 d=dd;
 dy_x=-d[3]*y*(1-y/d[1]);
 return(dy_x);
finish;

start _logistic_growth(d) global(dd,x,tree);
 /*-- Stage 1  of the hierarchical model --*/
 /*-- Within subjects variation          --*/
 dd=d;
 h={1E-12 1 1E-5};
 eps=1e-12;
 
 t=0//x;                 /* vector of time points for integration limits*/
 c=d[1]/(1+d[2]);        /* value at t=0*/
 tree1=tree[1];
 call ode(predv,"der_ode",c,t,h) eps=eps;
 return(t(predv));
finish _logistic_growth;

%mend useriml;  
   
/*---logistic model---*/
%nlinmix(data=tree,
   next =7,
   globvar=dbu,
   model=%str(
      bu      = _b||_u;
      bsm     = bsm(bu);
      dbsm    = dbsm(bu);

      par = 7 // {1,.};
      call nlpfdd(predv, dpredvd,h, "_logistic_growth", bsm) par=par;
      dbu     = dpredvd*dbsm;            /* Chain rule */
   ),

  derivs=%str(
      _db=dbu[,1:3];
      _du=dbu[,4];
      /* print _db; */
      ),
   parms=%str(b1=150 b2=10 b3=-.003),
   stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
      random d_u1 / subject=tree solution cl;
   ),
   expand   =  eblup, 
   debug=user_debug,
   options  =  debug    /* skipnlin   notes*/,
   imlfun   =  useriml,
   converge =  2e-5,
   nlmem    =  nlmem / nosource
);

