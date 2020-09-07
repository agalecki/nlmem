libname lib '.';

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
start fun(_bu) global(x,num,e,den);
      _b=_bu[1:3];
      _u=_bu[4:5];
      num = _b[1]+_u[1];
      e = exp(_b[3]*x);
      den = 1 + (_b[2]+_u[2])*e;
      predv = num/den;
return(predv);
finish;
%mend useriml;


/*---logistic model---*/

%nlinmix(data=tree,
   globvar= dbu,
   next =7,
   model=%str(
      _bu= _b||_u;
      par = {7};
      call nlpfdd(predv,dbu,h, "fun",_bu) par=par;
      *print predv /*  dbu segmentation violation */;
      ),

   derivs=%str(
      _db=dbu[,1:3];
      _du=j(_ni,_nu,.); /* derivs _du will be ignored due to dtmplt.
                           Numerical derivatives will be calulated instead */
      ),
   dtmplt = 1 1 1 0 0, 
   parms=%str(b1=150 b2=10 b3=-.001),
   stmts=%str(
      class tree;
      model pseudo_y = d_b1 d_b2 d_b3 / noint notest solution cl;
      random d_u1 d_u2 / subject=tree type=un solution cl;
      ods output IterHistory = iter;
   ),
   expand   =  1 0,
   imlfun   =  useriml,
   converge =  1e-6,
   nlmem    =  nlmem / nosource
);
run;

