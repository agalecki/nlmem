/*
ORANGE TREE DATA. 
LOGISTIC GROWTH MODEL. 
ODE SPECIFICATION

PARMS      = b1=150 b2=10 b3=-.001
RANDOM     U1
NEXT       =7
SKIPNLIN   no
IMLNLIN    usernlin , constraints 
DEBUG      yes
USER_DEBUG yes
DERIVS     nlpfdd
APPEND     no
CONVERGE   2e-5

CPU time:  40secs
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
   print _b _u; 
   print predv _dbu; 
%mend  user_debug;

%macro usernlin;
 *print _nobs;
  optn= _nobs ||{0};
 *CON defines constraints;
  con = {100  5 -0.5,
         500 50 -0.00001};
  call nlplm(rc,_bnlin,"_objf",_b,optn) jac="_dobjf" blc=con ;
  print _bnlin;
%mend usernlin;


%macro useriml;

start bsm(bu);
/* Between subject model */
/* Extract _b and _u */
   _b=bu[1:3];
   _u=bu[4];
   bsm = i(3)*_b + {1,0,0}*_u;
 return(bsm);
finish bsm;

start der_ode(x,y) global(dd);
/*-- Logistic model. ODE specification  --*/
 d=dd;
 dy_x=-d[3]*y*(1-y/d[1]);
 return(dy_x);
finish;

start _logist_growth(d) global(dd,x,tree);
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
finish _logist_growth;

%mend useriml;  
   
/*---logistic model---*/
%nlinmix(data=tree,
   next=7, 
   globvar=dbu,
   model=%str(
      bu      = _b||_u;
      bsm     = bsm(bu);
      predv   = _logist_growth(bsm);
    ),

   parms=%str(b1=150 b2=10 b3=-.001),
   stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
      random d_u1 / subject=tree solution cl;
   ),
   expand   =  eblup, 
   options  =  debug    /* skipnlin   notes*/,
   imlfun   =  useriml,
   converge =  2e-5,
   imlnlin  =  usernlin,
   debug    =  user_debug,
   nlmem    =  nlmem / nosource
);


  
