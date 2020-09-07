/* No random effects */

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
%let next=7;  
%let savevar= e den;

start predv(_b) global(x,e,den);
      *print bu;
      num = _b[1];
      e = exp(_b[3]*x);
      den = 1 + _b[2]*e;
      predv = num/den;
return(predv);
finish;
%mend useriml;  
   
/*---logistic model---*/
%nlinmix(data=tree,
   globvar = dbu,
   model=%str(
      par = _ni // {1,.};
      call nlpfdd(predv, dbu,h, "predv", _b) par=par;
      /* print predv*/ /* dbu segmentation violation */;
   ),

  derivs=%str(
      _db=dbu[,1:3];
      /* print _db; */
      ),
   parms=%str(b1=150 b2=10 b3=-.001),
   stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
   ),
   expand     =  eblup, 
   options    =  debug   /* skipnlin  notes*/,
   imlfun     =  useriml,
   nlmem      =  nlmem / nosource
);

