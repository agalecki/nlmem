/*
ORANGE TREE DATA. 
LOGISTIC GROWTH MODEL.
ANALYTICAL SPECIFICATION

RANDOM    U1
NEXT      =7
SKIPNLIN  no
DEBUG     no
DERIVS    nlpfdd
APPEND    no
CPU time: 4.41 secs
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

%macro useriml;
start predv(bu) global(x);
      *print bu;
      _b=bu[1:3];
      _u=bu[4];
      num = _b[1]+_u[1];
      e = exp(_b[3]*x);
      den = 1 + _b[2]*e;
      predv = num/den;
return(predv);
finish;
%mend useriml;  
   
/*---logistic model---*/
%nlinmix(data=tree,
   next =7,
   globvar = dbu,
   model=%str(
      bu  = _b||_u;
      par = {7,1,.};  /* 7 can be replaced with _ni */
      call nlpfdd(predv, dbu,h, "predv", bu) par=par;
      /* print predv*/ /* dbu segmentation violation */;
   ),

  derivs=%str(
      _db=dbu[,1:3];
      _du=dbu[,4];
      /* print _db; */
      ),
  parms=%str(b1=150 b2=10 b3=-.001),
  stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
      random d_u1 / subject=tree solution cl;
      ),
  expand   =  eblup, 
  imlfun   =  useriml,
  nlmem    =  nlmem / nosource
);

