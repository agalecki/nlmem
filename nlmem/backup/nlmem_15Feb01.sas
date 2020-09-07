%let nlmemlocal= read_all _resdr iml_catalog  imlnlinw imlnlin_default
      imlstore_catalog imlload_catalog dbtmplt dutmplt;

%macro default_initobjf;
free _objf _dobjf;
%mend  default_initobjf;

%macro default_objf;
/* Cumulates weighted differences into vector _objf 
   and corresponding derivatives into _dobjf matrix
*/
      %if not (%bquote(&weight)= ) %then %do;
/* predv_y0 = predv - _observedy  with some zero elements */
      predv_y0   = sqrt(_weight_) # predv_y0;  
       _db       = sqrt(_weight_) # _db;
      %end;

    _objf       = _objf     // predv_y0; 
    _dobjf      = _dobjf    // _db;    
%mend  default_objf;

%macro default_nlin;
/*
 Default IML code to solve nonlinear least-squares problem
  using NLPLM.

 The following vars and vectors are available:
  _nobs  number of obs in the input dataset
  _b  
*/
 *print _nobs;
  optn= _nobs ||{5};
  call nlplm(rc,_bnlin,"_objf",_b,optn) jac="_dobjf";
  print _bnlin;
%mend  default_nlin;

%macro default_debug;
/*
 Default IML code for  debugging
*/
 *print _nobs;
  optn= _nobs ||{0};
  call nlpfdd(_objfv,g,h,"_objf",_b,optn) jac="_dobjf";
  print _objvf _dobjf;

%mend  default_debug;

%macro nlmem_default2;

/* Redefine IMLNLIN -> IMLNLINW */
%let imlnlinw=;
%let imlnlin_default= default_nlin  default_objf default_initobjf;
%do i=1 %to 3;
%let temp=%scan(&imlnlin,&i);
%if %length(&temp)=0 %then %let temp=%scan(&imlnlin_default,&i);
%let imlnlinw=&imlnlinw &temp; 
%end;

%let _resdr=sub_objf;
%let read_all  =_reads;
%let subject= &next;

data _null_;
length readc $ 40;
readc ="&next";

c1=substr(left(readc),1,1);
if indexc(c1,'0123456789') then do;
      call symput('read_all','_readn');
      call symput('subject','');
   end;
run;
%mend nlmem_default2;

%macro aux(list,oper);
%scan(&list,1)
%let i=1;
%do %while(%length(%scan(&list,&i+1)) );
%let i=%eval(&i+1);
&oper %scan(&list,&i)
%end;
%mend aux;

/* SUB_OBJF and RES_DR2 are _res_dr macros.

   Following vectors are available in _res_dr macros

       _isubjct  is subject #
       _iobs     vector of indices pointing to _isubjct
       _iobs1    first element of _iobs

       _b        vector of fixed effects extracted from data
       _u        vector of random effects extracted from data
       _ni       number of records per subject.

       &dtvar and &xvar from the input data

*/

%macro sub_objf;
/* Called from _objf  function */
   _u= _subject(_iobs);
   predv  = _modelb(_b);
   predv_y0    = choose(_observedy=., 0,  predv - _observedy); 

   %if not (%bquote(&weight) = ) %then %do;
      _weight_    = _weightb(_b);
   %end;

     _db         = _derivsb(_b);
   *print predv _db;

   %let mnlin     = %scan(&imlnlinw,2);
   %&mnlin;

   %if not (%bquote(&retain)= ) %then %do;
     _retain   = _retainb(_b);
   %end;
*print 'Macro sub_objf ends ....';
%mend sub_objf;

%macro res_dr2;
  _u= _subject(_iobs);
  *print 'res_dr2' _u;
  _bu=_b||_u;
  predv  = _model(_bu);
  &pseudoresp = _observedy - predv;
   _pseudoy = &pseudoresp;

   %if not (%bquote(&weight)= ) %then %do;
      _weight_    = _weight_(_bu);
   %end;
 
   _dbu= _derivs(_bu);

      %if (&nb) %then %do; 
         do i=1 to _nb;
            &pseudoresp  = &pseudoresp + _dbu[,i]*_b[i];
         end;   
      %end;

      %if (&nu) & (&expandu=EBLUP) %then %do;            
          do i=1 to _nu; 
            if (_expandn[i]) then  
            &pseudoresp = &pseudoresp + _dbu[,(_nb+i)]*_u[i];
          end;
      %end;
   %if %length(&debug) %then %&debug;

   _appnd=predv || _pseudoy || &pseudoresp || _dbu;


   %if %length(&savevar) %then %do;
    _savevar=_savevar(.);
    _appnd =_appnd || _savevar;
   %end;

   %if not (%bquote(&weight)= ) %then %do;
    _appnd =_appnd || _weight_;
   %end;

   append from _appnd;
   %if not (%bquote(&retain)= ) %then %do;
     _retain   = _retain(_bu);
   %end;

%mend res_dr2;

/*  READ_ALL macros   */

%macro _reads;
/* Reads in subject with various number of records at once */
read var {&subject} into _idk;
_ni=.;
_iobs=1;
_isubjct=0;

do _j=2 to _nobs;
 read var {%unquote(&subject)} into _idj point _j;
  if _idk = _idj then _iobs=_iobs//_j;
   else do;
    /* Next subject */
     _isubjct=_isubjct+1;
     _ni=nrow(_iobs);
     %&_resdr

    /* Initiate new subject */
     _idk=_idj; _iobs=_j;
   end;
end;

/* Last subject */
   _isubjct=_isubjct+1;
   _ni=nrow(_iobs);
   %&_resdr
%mend _reads;

%macro _readn;
_ni=&next;
do _i=_ni to _nobs by _ni;
_isubjct=_i/_ni;
_iobs = (_i-_ni+1):_i;
   %&_resdr
end;
%mend _readn;

%macro _busplit;
/* _bu -> _b _u */
 _b=.;
 _u=.;
 if _nb then _b=_bu[,(1:_nb)];
 if _nu then _u=_bu[,(_nb+1):(_nu+_nb)];
%mend _busplit;

%macro _pseudoder;
start _pseudoder(_b) global(_nb,_nu,_call, _nobs,_observedy,_ni
                           );
  /* _weight_=.;*/
  %if %length(&expandn) %then _expandn ={&expandn};;
  _objf=.;
  create _nlinmx2 var {predv  _pseudoy
                       &pseudoresp  &d_fixed &d_random &savevar
                       %if not (%bquote(&weight) = ) %then _weight_;
  };
  use &outdata var {&dtallvar &random};
  setin &outdata NOBS _nobs;
  %&read_all;
  close _nlinmx2;
return(_objf);
finish;
%mend _pseudoder;

%macro _imlinit;

start _subject(_iobs) global(_iobs1,_u,_newsub,_observedy,
                              %aux(&xvar &tvar, %str(,)) );
 setin &outdata NOBS _nobs;
 _iobs1=_iobs[1];
 _newsub= _iobs=_iobs1;
 %if &next=1 %then %do;
   read var {&xvar &tvar &random} point _iobs1;
 %end;
 %else %do;
   read var {&xvar &random} point _iobs1;
   read var {&tvar} point _iobs;
 %end;

 _u=.;
 %if (&nu) %then _u=%aux(&random,||);;
 _observedy=&response;
return(_u);
finish _subject;

  _nb=0;
  _fixed=.;

  %if (&nb) %then %do;
      _nb=&nb;
      _fixed =j(1,_nb,rowcat({[40] ' '}));
      %do i=1 %to &nb;
      %let fixedi=%scan(&fixed,&i);
      _fixed[&i]="&fixedi";
      %end;
   %end;

   _u=.;
   _nu=0;
   _random=.;
   %if (&nu) %then %do;
      _nu=&nu;
      _random =j(1,_nu,rowcat({[40] ' '}));
      %do i=1 %to &nu;
      %let randomi=%scan(&random,&i);
      _random[&i]="&randomi";
      %end;
   %end; 

/* _b */
   use &outdata var {&fixed} NOBS _nobs;
   %if (&nb) %then %do;
       read next 1;
      _b=%aux(&fixed,||);
   %end;
   close &outdata;
   free &fixed;
   
/* Global vars available to use in imlinit ??? statement below:
   _b0  , _nb, _fixed,  
   _u0=., _nu, _random
  
   _nobs
*/


  %if %length(&imlload_catalog) %then %do;
    %let i=0;
    %do %while (%length(%scan(&imlload_catalog,&i+1,' ')));
    %let i=%eval(&i+1);
    %let item=%scan(&imlload_catalog,&i,' ');
    reset storage=&item;
    load;
    %end;
  %end;
%mend _imlinit;

%macro _nlinobjf;
/* Objective function for nonlinear regression */ 
start _objf(_b) global( _call, _nobs, _observedy, _dobjf,_ni
                      );
 *print '_objf' _b;
 use &outdata var {&dtallvar &random};
  _call=_call+1;
  
/* Examples of nlin_initobjf macros:
  Ex1. free _objf _dobjf;   Default
  Ex2. _objf=0; _dobjf=0;
*/
  %let mnlin=%scan(&imlnlinw,3); /* Default: nlin_initobjf */
  %&mnlin;
  %&read_all;
 close &outdata;
return(_objf);
finish _objf;

start _dobjf(_b) global (_dobjf,_ni);
 return(_dobjf);
finish _dobjf;
%mend _nlinobjf;

%macro _numd;
    _b0=_b;
    _tol=&tol*(1+abs(_b0));                  

  %do i=1 %to &nb;
   %if %scan(&dbtmplt,&i) = 0 %then %do;
    _mcall=_mcall+1;
    _b =_b0;
    _b[&i]=_b0[&i]-_tol[&i];
    &model
    _predvl=predv;
   
    _mcall=_mcall+1;
    _b    =_b0;
    _b[&i]=_b0[&i]+_tol[&i];
    &model
    _predvu = predv;
    _db[,&i]= (_predvu-_predvl) /2/_tol[&i];
    %if not (%bquote(&retain) = ) %then %do;
     predt=predt||_predvl||_predvu;
    %end;
   %end;
  %end;
  _b      =_b0;
%mend _numd;

%macro nlin;
proc iml;
 
start _modelb(_b) global(_u, _ni,_nb,_nu,
                         %aux(&dtallvar &globvar, %str(,))
                         );
    &modinit
    _mcall=1;
    &model
    %if not (%bquote(&retain) = ) %then %do;
        predt=predv;
    %end;
  return(predv);
finish;

%if not (%bquote(&weight) = ) %then %do;

start _weightb(_b) global(_u, _ni,_nb,_nu,
                         %aux(&dtallvar &globvar, %str(,))
                         );
&weight
return(_weight_);
finish;
%end;

start _derivsb(_b) global(_u,_ni,_nb,_nu,
                          %aux(&dtallvar &globvar, %str(,))
                         );
      &derivs
      *print _db;
      
      %if %index(&dbtmplt,0) %then %do;
      %if %bquote(&derivs)= %then _db=j(_ni,_nb,.);;
      _mcall=1;            /* Model call number */
      %_numd
      %end;
return(_db);
finish _derivsb;

%if not (%bquote(&retain) = ) %then %do;
start _retainb(_b) global(_u,  _ni,
                         %aux(&dtallvar &globvar, %str(,))
                         );
  _retain=.;
  &retain
  return(_retain);
finish _retainb;
%end;

 %let _resdr=sub_objf;
 %_nlinobjf;                              
 %_imlinit;

 _call=0;

 %let mnlin=%scan(&imlnlinw,1);
 %&mnlin;   /* By default: default_nls macro */

create _bnlin var {&fixed};
  append from _bnlin;
close _bnlin;
quit;

/* PROC IML: Creates _nlinmix dataset*/
%pseudoder_nlin;

%mend nlin;


 /*----------------------------------------------------------*/
 /*                                                          */
 /*   %pseudoder                                             */
 /*   construct residuals and derivatives   SAS/IML          */
 /*                                                          */
 /*----------------------------------------------------------*/

%macro pseudoder_common;

start _model(_bu) global(_nb,_nu,_ni,
                         %aux(&dtallvar &globvar, %str(,))
                        );
  %_busplit;
  &modinit
  _mcall=1;
  &model
  %if not (%bquote(&retain) = ) %then %do;
    predt=predv;
  %end;
  return(predv);
finish;

%if not (%bquote(&weight) = ) %then %do;

start _weight_(_bu) global(_nu,_nb, _ni,
                           %aux(&dtallvar &globvar, %str(,))
);
%_busplit;
&weight
return(_weight_);
finish;
%end;

start _derivs(_bu) global(_nu,_nb, _ni,
                          %aux(&dtallvar &globvar, %str(,))
                          );
 free _dbu;
 %_busplit;
 &derivs

 _mcall=1;            /* Model call number */
 %if (&nb) & %index(&dbtmplt,0) %then %do;
     %if %bquote(&derivs)= %then %do;
      _db=j(_ni,_nb,.);
     %end;
 
      %_numd
 %end;
 if _nb then  _dbu=_dbu || _db; 

 %if (&nu) & %index(&dutmplt,0) %then %do;
     %if %bquote(&derivs)= %then %do;
      _du=j(_ni,_nu,.);
     %end;

      _u0=_u;
      _tol=&tol*(1+abs(_u0));                  

     %do i=1 %to &nu;
       %if %scan(&dutmplt,&i) = 0 %then %do;
         _mcall=_mcall+1;
         _u=_u0;
         _u[&i]=_u0[&i]-_tol[&i];
         &model
         _predvl=predv;
   
         _mcall=_mcall+1;
         _u=_u0;
         _u[&i]=_u0[&i]+_tol[&i];
         &model
         _predvu = predv;
         _du[,&i]=  (_predvu-_predvl) /2/_tol[&i];
           %if not (%bquote(&retain) = )  %then %do; 
              predt=predt||_predvl||_predvu;
           %end;
       %end;
      %end;
   _u =_u0;
 %end;
 if _nu then _dbu=_dbu||_du;
 %if %index(&options,DEBUG) %then &model;
return(_dbu);
finish;

%if not (%bquote(&retain) = ) %then %do;
start _retain(_bu) global(_nu, _nb,  _ni,
                         %aux(&dtallvar &globvar, %str(,))
                         );
  %_busplit;
  _retain=.;
  &retain
  return(_retain);
finish _retain;
%end;

%if %length(&savevar) %then %do;
start _savevar(_dummy_) global( %aux(&savevar,%str(,)) );
  _savevar=%aux(&savevar,||);
  return(_savevar);
finish;
%end;

 %let _resdr=res_dr2;
 %_pseudoder
 %_imlinit;
%mend   pseudoder_common;

%macro pseudoder_nlin;
 proc iml;
 /* Create _nlinmx2 dataset */
   %pseudoder_common
   use _bnlin var {&fixed} ;
       read next 1;
      _bnlin=%aux(&fixed,||);
   close _bnlin;
   free &fixed;
   _objf = _pseudoder(_bnlin);
quit; /* quit IML */

data &outdata;
 merge &outdata _nlinmx2;
 _tempone=1;
 set _bnlin point=_tempone;
run; 

%if %index(&options,DEBUG) %then %debug(&outdata, -8);
%mend pseudoder_nlin;

%macro pseudoder;
proc iml;
 %pseudoder_common
  _objf=_pseudoder(_b); 
quit; /* quit IML */
   
data &outdata;
 merge &outdata _nlinmx2;
run; 
   
%if %index(&options,DEBUG) %then %do;
   proc means data=&outdata uss;
      var _pseudoy;
   run;
%end;
%mend pseudoder;

%macro nlmem_init;
%let iml_catalog=work._nlmem;

/***** iml_catalog -> imlstore + imlload macro vars ***********/
%let indxc=%index(&iml_catalog,|);
%if  &indxc=1 %then %let iml_catalog=.&iml_catalog;
%let imlstore_catalog=%scan(&iml_catalog,1,'|');
%let imlload_catalog =%scan(&iml_catalog,2,'|');
%if %length(&imlload_catalog)=0 %then %let imlload_catalog=&imlstore_catalog;
%if &imlstore_catalog=. %then %let imlstore_catalog=;

   %if %length(&imlfun) %then %do;
   proc iml;  
   %&imlfun
   _dummy=.;
   %if %length(&imlstore_catalog) %then %do;
         reset storage=&imlstore_catalog;
         store;
   %end;
   quit;
   %end;
      
   %nlmem_default2
   %init;
  

%let dtmpltw=&dtmplt;
%if (%bquote(&derivs)= ) & %length(&dtmplt)=0 %then %do;
%do i=1 %to (&nb+&nu);
%let dtmpltw=&dtmpltw 0;
%end;
%end;

%if  %length(&derivs &dtmplt)=0 %then %do;
%do i=1 %to (&nb+&nu);
%let dtmpltw=&dtmpltw 1;
%end;
%end;

%let dbtmplt=;
    %do i=1 %to &nb;
      %let dbtmplt= &dbtmplt %scan(&dtmpltw,&i);
    %end;

%let dutmplt=;
    %do i=1 %to &nu;  
      %let dutmplt= &dutmplt %scan(&dtmpltw,(&nb+&i));
    %end;
%mend nlmem_init;
