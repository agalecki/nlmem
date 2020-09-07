%include "&sysparm";

%macro symcheck(name);
%if %nrquote(&&&name)=%nrstr(&)&name %then %let rc=0;
%else %let rc=1;
%put rc=&rc;
&rc
%mend symcheck;

%macro symvalue(name);
 %if %symcheck(&name)=0 %then %let rvalue=;
 %else %let rvalue=&&&name;
 %put rvalue=&rvalue;
 &rvalue
%mend symvalue;

%macro default(input_value,default_value); 
  %let output_value=&input_value;
  %if %length(&input_value)=0 %then %let output_value=&default_value; 
  &output_value
%mend  default;


%let data    =%symvalue(data);
%let modinit =%symvalue(modinit);
%let model   =%symvalue(model);
%let weight  =%symvalue(weight);
%let derivs  =%symvalue(derivs);
%let tol     =%symvalue(tol);
%let parms   =%symvalue(parms);
%let stmts   =%symvalue(stmts);
%let expand  =%symvalue(expand);     
%let converge=%symvalue(converge);
%let maxit   =%symvalue(maxit);
%let switch  =%symvalue(switch);
%let append  =%symvalue(append);
%let procopt =%symvalue(procopt);
%let nlinopt =%symvalue(nlinopt);
%let nlinstmt=%symvalue(nlinstmt);
%let outdata =%symvalue(outdata);
%let options =%symvalue(options);
%let nlmem   =%symvalue(nlmem);
%let imlfun  =%symvalue(imlfun); 
%let globvar =%symvalue(globvar);
%let retain  =%symvalue(retain);
%let next    =%symvalue(next);
%let xvar    =%symvalue(xvar);
%let tvar    =%symvalue(tvar);
%let savevar =%symvalue(savevar);
%let dtmplt  =%symvalue(dtmplt);
%let debug   =%symvalue(debug);
%let imlnlin =%symvalue(imlnlin);

/* Default values assigned if blank */
%let data     = %default(&data,    _last_);
%let tol      = %default(&tol,     1e-5);
%let expand   = %default(&expand,  zero);
%let converge = %default(&converge,1e-8);
%let maxit    = %default(&maxit,   30);
%let switch   = %default(&switch,  0);
%let outdata  = %default(&outdata, _outdata);
%let next     = %default(&next,    1);
%let imlnlin  = %default(&imlnlin, default_nlin default_objf 
    default_initobjf);


%let path=../..;
filename  sasdata "&path/datasets/&data._dt.sas";
filename  nlinmix "&path/nlmm800.sas";
%include  sasdata;
%include  nlinmix / nosource;


Title "NLINMIX: Data=&data";
%nlinmix(data=&data,   /* _LAST_  */
   modinit=&modinit,
   model=&model,
   weight=&weight,
   derivs=&derivs,
   tol=&tol,           /* default 1e-5  */
   parms=&parms,
   stmts=&stmts,
   expand=&expand,     /* zero */
   converge=&converge, /* 1e-8 */
   maxit=&maxit,       /* 30   */
   switch=&switch,     /* 0    */
   append=&append,
   procopt=&procopt,
   nlinopt=&nlinopt,
   nlinstmt=&nlinstmt,
   outdata=&outdata,   /*_nlinmix*/
   options=&options,
   nlmem=&nlmem,
   imlfun=&imlfun,  
   globvar=&globvar,
   retain=&retain,
   next=&next,        /* 1 */
   xvar=&xvar,
   tvar=&tvar,  
   savevar=&savevar,
   dtmplt= &dtmplt,
   debug = &debug,
   imlnlin=&imlnlin   /* default_nlin default_objf default_initobjf */
);

