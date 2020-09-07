/*
ORANGE TREE DATA.
LOGISTIC GROWTH MODEL.
ANALYTICAL SPECIFICATION

RANDOM    U1
NEXT      =1         Default
SKIPNLIN  no
DEBUG     blank
DERIVS    blank
APPEND    blank
OPTIONS   DEBUG
CPU time: 6.6 secs
*/


/*---Data on orange trees, Draper and Smith, 1981, p. 524, and
     Lindstrom and Bates, Biometrics, 1990, 673-687---
*/
%let path=../..;   
filename tree_dt "&path/datasets/tree_dt.sas";
filename nlinmix "&path/nlmm800.sas";
filename nlmem   "&path/nlmem.sas";
%include tree_dt;   
%include nlinmix / nosource;

/*---logistic model---*/

%nlinmix(data=tree,
   model=%str(
      num = _b[1]+_u[1];
      e = exp(_b[3]*x);
      den = 1 + _b[2]*e;
      predv = num/den;
   ),

   parms=%str(b1=150 b2=10 b3=-.001),
   stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
      random d_u1 / subject=tree solution cl;
      ),
   expand   = eblup,
   nlmem    = nlmem / nosource
);

