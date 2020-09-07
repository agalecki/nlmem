/*
tree1a.sas illustrates how does:
  switch= argument 
AND
 hybrid estimation work
*/

/*---Examples for NLINMIX macro---*/

/*---Orange trees data, Draper and Smith, 1981, p. 524, and
     Lindstrom and Bates, Biometrics, 1990, 673-687---*/

%let path=../..;
filename tree_dt "&path/datasets/tree_dt.sas";
filename nlinmix "&path/nlmm800.sas";
%include tree_dt;
%inc  nlinmix / nosource;


/*---logistic model---*/
%nlinmix(data=tree,
   model=%str(
      num = b1+u1;
      e = exp(b3*x);
      den = 1 + (b2+u2)*e;
      predv = num/den;
   ),
   derivs=%str(
   d_b1=1/den;
   d_b2=-num/den/den*e;
   d_b3=-num/den/den*b2*x*e;
   d_u1=d_b1;
   d_u2=d_b2;
   ),
   parms=%str(b1=190 b2=10 b3=-.001),
   stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
      random d_u1 d_u2/ subject=tree solution cl type=un;
   ),
   expand= 1 0,
   switch=3, 
   options=debug 
)
run;

proc print data=alldebug;
run;
