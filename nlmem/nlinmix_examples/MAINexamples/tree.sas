/*---Examples for NLINMIX macro---*/

/*---Data on orange trees, Draper and Smith, 1981, p. 524, and
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
      den = 1 + b2*e;
      predv = num/den;
   ),
   parms=%str(b1=150 b2=10 b3=-.001),
   stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
      random d_u1 / subject=tree solution cl;
   ),
   expand=eblup
);


