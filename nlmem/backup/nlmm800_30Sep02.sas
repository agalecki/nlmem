/*-----------------------------------------------------------
 TITLE                                                        
 -----                                                        
 NLINMIX: a SAS macro for fitting nonlinear mixed models      
 using PROC NLIN and PROC MIXED. Requires SAS/STAT Version 8.    
                                                              
 SUPPORT                                                      
 -------                                                      
 Russ Wolfinger (russ.wolfinger@sas.com)                          
 Please send email with any suggestions or corrections.       
                                                              
 HISTORY                                                      
 -------                                                      
 initial coding                                     01Jun92 rdw  
 revision                                           21Oct92 rdw  
 removed PROC NLIN from iterations and revised      08Mar94 rdw  
 changed syntax to include random effects                        
    explicitly, added numerical derivatives and                  
    optional Gauss-Newton steps                     19Dec94 rdw  
 added rtype=                                       09Feb95 rdw  
 changed tech= to expand=                           13Feb95 rdw  
 suggestions from Ken Goldberg, Wyeth-Ayerst        29May96 rdw  
 added random2                                      07Jan98 rdw
 more suggestions from Ken Goldberg: covparms   
    only used in iter 0 unless noprev, added   
    covpopt=,nlinopt=,nlinstmt=,outdata=,   
    nestrsub                                        16Nov98 rdw
 fixed problem from Ken Goldberg with random2,   
    eblup, and unequal subject2 block sizes         16Apr99 rdw
 added stmts= and deleted options involving    
    random, repeated, response, and Gauss-Newton    28Feb00 rdw
 added append= from Andrzej Galecki, U Mich.        07Mar00 atg
 additional improvements from Andrzej Galecki       04Apr00 atg
 added nlmem=, enhanced syntax for expand=          15Feb01 atg
 added linear inequality constraints for fixed
                 parameters (neq_constraints=) 
       not adjustable parameters (paramfixed=)
       user defined debugging  (user_debug=)
       nonlinear regression using calls to 
                 PROC MIXED (nlinmethod=)           30Sep02 atg  

 DESCRIPTION                                                  
 -----------                                                  
 Available estimation methods are as follows, all of which    
 are implemented via iterative calls to PROC MIXED using      
 appropriately constructed pseudo data:                       

    1. Expanding the nonlinear function about random effects  
    parameters set equal to zero, which is similar to but not
    the same as Sheiner and Beal's first-order method (NONMEM User's 
    Guide, University of California, San Francisco).  This method is 
    the default for random effects specifications; that is,
    those using RANDOM statements in PROC MIXED.

    2. Expanding the nonlinear function about random effects  
    parameters set equal to their current empirical best      
    linear unbiased predictor (EBLUP), which is Lindstrom and 
    Bates' approximate second-order method (Lindstrom and     
    Bates, 1990, Biometrics 46, 673-687).                     
       The method is implemented using an intitial call to    
    PROC NLIN and then iterative calls to PROC MIXED.  This   
    differs from Lindstrom and Bates' implementation, in      
    which pseudo data are constructed and iterative calls are  
    made alternately to PROC MIXED and PROC NLIN.  This       
    macro's implementation requires much less time and space  
    for larger problems, and it works because solving the     
    mixed-model equations in PROC MIXED is equivalent to      
    taking the first Gauss-Newton step in PROC NLIN           
    (Wolfinger, 1993, Biometrika 80, 791-795).  Use           
    EXPAND=EBLUP to implement this method, which applies for       
    RANDOM statement specifications in PROC MIXED.                       

    3. Iteratively fitting a covariance structure outside of  
    the nonlinear function.  This solves a type of second-order 
    generalized estimating equations (Prentice and    
    Zhao, 1991, Biometrics 47, 825-839). This method uses the 
    REPEATED statement in PROC MIXED.                                  


 SYNTAX                                                       
 ------                                                       
 %nlinmix(                                                    
    data=,                                                    
    model=,                                                   
    modinit=,
    derivs=,                                                  
    tol=,                                                     
    parms=,                                                   
    adjust=,
    stmts=,                                                  
    weight=,
    expand=,                                                  
    converge=,                                                
    maxit=,                                                   
    procopt=,                                                 
    nlinopt=,
    nlinstmt=,
    switch=,
    append=
    outdata=,
    options=,                                                  
    neq-constraints=,
    paramfixed=,
    user_debug=
 )                                                            

 where the arguments are as follows:

  data     specifies the input data set. The default is
           the last created data set.

  model    specifies the nonlinear mixed model in terms of 
           SAS programming statements using parameters of 
           your choice.  You may use multiple SAS  
           statements and auxiliary variables from the
           input data set.  The names of all fixed- and random- 
           effects parameters should be unique and not the
           same as any of the variables in the input data set.       
           You must assign the final value of the model 
           to a variable named PREDV; this is the predicted      
           value for an observation. You should enclose the 
           entire block of code with the %str() macro (see example 
           syntax below). In addition, the fixed-effects       
           parameters must be listed in the PARMS= argument       
           and the random-effects in the RANDOM= argument.    
           In your specification of PREDV, you should scale all 
           parameters so that they are all are around the same order
           of magnitude to avoid instabilities in the
           algorithm. 
                                              

  modinit  specifies modeling statements to be called only    
           once at the beginning of the modeling step; this   
           option is usually used in conjunction with         
           initializing recursively defined models which      
           are to be differentiated numerically.  This code    
           should have no references to the fixed- or random- 
           effects parameters and should be enclosed with
           the %str() macro.                                

  derivs   specifies derivatives of PREDV in the preceding model 
           specification with respect to the fixed- and random-effects 
           parameters.  The derivatives are specified using the same 
           parameter names but with "d_" appended beforehand.  You 
           may use multiple statements, auxiliary variables,  
           and variables defined in the model specification,  
           and you should enclose the entire specification    
           with the %str() macro.  If you do not specify this 
           argument, derivatives are computed numerically using 
           central differences with width tol*(1+abs(parm)), where 
           tol is from the TOL= argument.  You may also specify only
           selected derivatives and the remaining ones will be
           computed numerically.

  tol      specifies the tolerance for numerical derivatives. 
           The default is 1e-5, and the width is computed as         
           tol*(1+abs(parm)), where parm is the parameter     
           being differentiated.                               

  parms    specifies the starting values of the fixed-effects 
           parameters in the form %str(b1=value b2=value ...  
           bN=value).  When not using the SKIPNLIN option,     
           the form can be %str(b1=values b2=values ...       
           bN=values), where the form of the values is        
           defined in the documentation of PROC NLIN's        
           PARAMETER statement.  You may also specify starting
           values for the parameters in the input data set,
           and these take precedence over those in the PARMS=
           specification if they are present.

  stmts    specifies PROC MIXED statements to be evaluated
           at each iteration.  You must include a MODEL
           statement with dependent variable equal to
           the response variable in your input data set, but with
           the prefix "pseudo_" (e.g. "pseudo_y").  The macro 
           creates a new version of this pseudo response variable 
           during each iteration.  The indepedent variables in
           the MODEL statement must consist of the list
           derivatives of PREDV with respect to the fixed-effect 
           parameters (all with prefix "d_"). You must also specify 
           the NOINT and SOLUTION options in your MODEL statement.
           Include the CL option if you want approximate confidence
           limits for the fixed-effects parameters.

           You may also include one or more RANDOM statements
           with effects being one or more of the derivatives of 
           of PREDV with respect to the random-effect parameters (again, 
           all with prefix "d_"). Each derivative variable must be 
           included in one and only one RANDOM statement, and their 
           names must be unique.  You must also specify the SUBJECT=
           option in every RANDOM statement.  Multiple RANDOM 
           statements with different SUBJECT= effects enable you
           to fit nested or crossed hierarchies in your model.  Finally, 
           if you are using EXPAND=EBLUP, then you must specify the 
           SOLUTION option in at least one of the RANDOM statements.  

           Other PROC MIXED statements that you can use include 
           CLASS (good for variables comprising the SUBJECT=
           effects in your RANDOM statements), PARMS (to specify 
           starting values for the covariance parameters for the PROC 
           MIXED calls), and ESTIMATE (to compute estimates of linear 
           combinations of the fixed-effects parameters).  PROC MIXED 
           syntax rules apply throughout.

         
  weight   specifies weights using SAS statements, which may use
           auxiliary variables defined in the model and derivs specifications.
           You assign the final value of the weight to a variable specified
           in WEIGHT statement in STMNTS. You should enclose the
           entire block of code with the %str(),  e.g.
         
           weight=%str(statements;
                  _weight_= expression;
                      ),

  expand   specifies the Taylor series expansion point. The default is
           ZERO (Method 1 above), and EBLUP specifies the current
           EBLUPs (Method 2 above). In addition to EXPAND=ZERO
           or EBLUP, syntax of the form EXPAND = 0 1 0 is allowed.
           Assuming there are 3 random effects the latter means
           that ZERO expansion for the first and third random effect and
           EBLUP expansion for the second random effect is used.
           EXPAND=ZERO is equivalent to EXPAND= 0 0 0. EXPAND=EBLUP is
           equivalent to EXPAND = 1 1 1. Similar feature called
           hybrid expansion is implemented in NONMEM.
           If you have no RANDOM statements in your STMTS= specification, 
           then EXPAND= has no effect.

  converge specifies the convergence criterion for the macro _and_ . 
           for nonlinear regression (using NLINMETHOD=1). The default is CONVERGE=1e-8 
           1e-1.                                       

  maxit    specifies the maximum number of iterations for the 
           macro to converge _and_ maximum number for initial nonlinear regression
           (using NLINMETHOD=1).  The default is MAXIT=30, which means that maximum
           # of iterations is 30 for both NLIN and NLINMIX iterations.            


  procopt  specifies options for the PROC MIXED statement.     

  nlinmethod
           indicates method to perform initial NLIN step.
           By default is = 1 and iteratively calls PROC MIXED.
           If nlinmethod=2 and nlmem is blank (non-IML execution mode) 
             then PROC NLIN is used.
           If nlinmethod=2 and nlmem is not blank(IML execution mode)
              then NLPFDD is used instead of PROC NLIN.
           For other values of nlinmethod argument initial NLIN step is omitted
              (same effect as SKIPNLIN option)

  nlinopt  specifies options for the PROC NLIN statement.

  nlinstmt specifies statements to be included in the PROC NLIN call.

  switch   specifies an iteration number to switch from EXPAND=ZERO
           to EXPAND=EBLUP.


  nlmem    points to a file containing a set of SAS macros, called NLMEM.
           NLMEM retains all the benefits of NLINMIX while allowing the
           systematic part of the model to be specified using IML syntax.   
           In particular it allows us to address advanced population 
           pharmacokinetics and pharmacodynamics models specified 
           by ordinary differential equations.
           NLMEM and examples can be downloaded from URL:
           http://www-personal.umich.edu/~agalecki/
             

  append   provides a means to keep data sets from every iteration.
           The syntax is one or more strings like the following:

              ds -> dsall 

           where ds is the name of a data set created during an
           iteration and dsall is the name of the new data set
           created by appending all of the different instances of
           ds throughout the iteration history.  The name ds must
           be _fit, _soln, _cov, _solnr, or a name that you specify
           for other tables using an ODS statement in your STMTS=
           code.

  outdata  specifies the name of the output working data set
           produced by NLINMIX.  The default name is _NLINMIX.


  neq_constraints
           points to a macro containing linear inequality constraints
           for parameters listed in PARMS argument. Constraints are 
           checked after every execution of PROC MIXED. If any of the constraints
           is not satisfied the step size in NLINMIX iterations is appropriately 
           reduced.
           Example of constraints specification for 3 parameters: b1,b2,b3.
           b1>=0; b2-b1>=0; 1< b3<=10
            
  paramfixed
           Specifies which of model fixed parameters are not adjusted
           during NLINMIX iterations.

  user_debug
           points to a dataset containing a list of user defined macros.
           By default it points to <default_debug> dataset containing
           information on default debug macros.


  options  specifies nlinmix macro options: 

     noprint          suppresses all printing   

     notes            requests printing of SAS notes, date, and page
                      numbers during macro execution.  By default, 
                      the notes, date, and numbers are turned off
                      during macro execution and turned back on 
                      after completion.

     noprev           prevents use previous covariance parameter
                      estimates as starting values for the next
                      iteration                               

     printall         prints all PROC NLIN and MIXED          
                      steps; the final PROC MIXED step        
                      is printed by default                   

     printfirst       prints the first PROC NLIN and          
                      MIXED runs                              

     skipnlin         skips the initial nonlinear regression.
                        (PROC NLIN or iterative calls to PROC MIXED)
                      

 EXAMPLE SYNTAX                                               
 --------------  

 see the nlmm801e.sas file                                             

 OUTPUT
 ------

 NLINMIX prints an iteration history to the log (which you should
 check to verify convergence) and then outputs the final call to PROC
 MIXED. The fixed- and random-effect parameter estimates are prefixed
 by "d_" in the output.  Refer to the PROC MIXED documentation for
 descriptions of each table.  In addition to the OUTDATA= data set,
 there are also a number of other data sets created by the macro 
 (check the notes at the end of the log for their names).
 

 DISCLAIMER                                                   
 -----------                                                  

 This macro is provided by SAS Institute Inc as a service to 
 its users.  It is provided "as is".  There are no 
 warranties, expressed or implied, as to merchantability or   
 fitness for a particular purpose regarding the accuracy of   
 the materials or code contained herein.                      

 -------------------------------------------------------------*/

 /*--------------------------------------------------------------*/
 /*                                                              */
 /*    %expand                                                   */
 /*    Expands argument.                                         */
 /*                                                              */
 /*--------------------------------------------------------------*/

%macro expand;
%let arg=ZERO;
%if (&argx)  %then %let arg=EBLUP;

%if %index(&expandu,&arg) & (&nu) %then %do;
 %let expandn=;
 %do i=1 %to &nu;
  %let expandn=&expandn &argx;
 %end;
%end;
%mend expand;

%macro aux(list,oper);
%scan(&list,1)
%let i=1;
%do %while(%length(%scan(&list,&i+1)) );
%let i=%eval(&i+1);
&oper  %scan(&list,&i)
%end;
%mend aux;

%macro mtstmts;
/*
Macro mtstmts translates stmts -> tstmts
Based on indadjset omits terms in the model statement 
Input macro vars: stmts nb indadjset

*/

%if %length(&stmt2)=0 %then %do; /* No | in STMTS argument */
 %let iv = 1;
 %let tstmts=;
 %let tstmts2=;
 %do %while (%length(%scan(&stmts,&iv,;)));
   %let stmt = %qscan(&stmts,&iv,%str(;));
   %let first = %qscan(&stmt,1);
   %if %index(%qupcase(&first),MODEL) %then %do;
    %let temp1=%qscan(&stmt,1,=);
    %let temp2=%qscan(&stmt,2,=/);
    %let temp3=%qscan(&stmt,2,/);
    %let temp2x=;
    %do i=1 %to &nb;
       %let cx= %scan(&temp2,&i,' ');
     %if %scan(&indadjset,&i,' ') %then
        %let temp2x=&temp2x &cx;
     %end;
     %let tstmts=&tstmts &temp1=&temp2x;

     %if %length(&temp3) %then %let tstmts=&tstmts / &temp3;;
     %let tstmts=&tstmts%str(;);
     %let tstmts2=&tstmts%str(;);
   %end;
    %else  %do;
     %let tstmts=&tstmts&stmt%str(;);
     %if not %index(%qupcase(&first),RANDOM) %then %let tstmts2=&tstmts2&stmt%str(;);
    %end;
   %let iv=%eval(&iv+1);
 %end;

 %if (&nlin) %then %let temp=&tstmts2;
  %else %let temp=&tstmts;
%end;
%else %do;  /* | in STMTS= argument */
 %if (&nlin) %then %let temp=&stmt1&stmt3;
  %else %let temp=&stmt1&stmt2;
%end;
&temp
%mend mtstmts;


 /*--------------------------------------------------------------*/
 /*                                                              */
 /*    %init                                                     */
 /*    Sets the initial values for the iterations.               */
 /*                                                              */
 /*--------------------------------------------------------------*/

%macro init;

%user_debugall;

 /*---determine number of parameters---*/
%let another = 1;
%let nb = 0;
%let prevpos = 1;
%let pos = 1;
%let _error = 0;
%let fixed = ;
%let d_fixed = ;
%let d_fixedn = ;
%do %while(&another = 1);
   %let quoteparms=%quote(&parms);
   %let effectx = %scan(&quoteparms,&pos,=);
   %let ipos=0;
   %do %while(%length(%scan(%quote(&effectx),&ipos+1,' '))); 
     %let ipos = %eval(&ipos + 1);
     %let effect = %scan(%quote(&effectx),&ipos,' ');
   %end;
   %if %length(&effect) & %length(%scan(&quoteparms,&pos+1,=)) %then %do;
      %let fixed = &fixed &effect;
      %let d_fixed = &d_fixed d_&effect;  
      %let d_fixedn = &d_fixedn der.&effect;  
      %let nb = %eval(&nb + 1);
      %let pos = %eval(&pos + 1);
   %end;
   %else %let another = 0;
%end;


/*---response and random effects---*/
%let another = 1;
%let nu = 0;
%let random = ;
%let d_random = ;
%let iv = 1;
%do %while (%length(%scan(&stmts,&iv,;)));
   %let stmt = %qscan(&stmts,&iv,%str(;));
   %let first = %qscan(&stmt,1);
   %if %index(%qupcase(&first),MODEL) %then %do;
      %let pseudoresp = %scan(&stmt,2,' ' =);
      %let i = %index(&pseudoresp,_);
      %let i = %eval(&i + 1);
      %let response = %substr(&pseudoresp,&i);
   %end;
   %else %if %index(%qupcase(&first),RANDOM) %then %do;
      %let i = %index(&stmt,/);
      %let ir = %index(%qupcase(&stmt),RANDOM);
      %let rndeff = %substr(&stmt,%eval(&ir+7),(%eval(&i-&ir-7)));
      %let ive = 1;
      %do %while (%length(%scan(&rndeff,&ive,' ')));
         %let d_effect = %scan(&rndeff,&ive,' ');
         %let effect = %substr(&d_effect,3);
         %let random = &random &effect;
         %let d_random = &d_random &d_effect;
         %let nu = %eval(&nu + 1);
         %let ive = %eval(&ive + 1);
      %end;
   %end;
   %let iv = %eval(&iv + 1);
%end;

%let expandu = %qupcase(&expand);

%if (&nu) %then %do;
  %let expandn=&expand;
  %let argx=0; %expand
  %let argx=1; %expand
  %let argx=0;
  %do i=1 %to &nu;
    %if (%scan(&expandu,&i) =1  ) %then %let argx=1;
  %end;
  %if (&argx) | %index(&expandu,EBLUP) %then %let expandu=EBLUP;
  %else %let expandu=ZERO;
%end;

%if not %index (&expandu,EBLUP) | (&switch > 0) %then   
   %let expandu = ZERO;

/* Adjustable parameters. May 2002. */

  %let notadjset=%qupcase(&paramfixed);
  %let numadjset=;
  %let indadjset=;
  %let adjset=;
  %if (&nb) %then %do;
   %do i=1 %to &nb;
    %let adjfixedi=%scan(&fixed,&i);
    %let j=0;
    %do %while(%length(%scan(&notadjset,&j+1,' ')));
     %let j=%eval(&j+1);
     %let notadjx=%scan(&notadjset,&j,%str( ));
     %if &notadjx=%upcase(&adjfixedi) %then %do;
       %let indadjset=&indadjset 0;
       %goto adjexit;
       %end;
    %end;
     %let indadjset=&indadjset 1;
     %let numadjset=&numadjset &i;
     %let adjset =&adjset &adjfixedi;
     %adjexit:
    %end;
  %end;
  %let adjloop=%aux(&numadjset,%str(,));


 /*---starting values for b and u---*/

%let adjparms=;

data &outdata;
   set %unquote(&data);
   %if (&nb) %then %do;
      array _b{&nb} &fixed;
      %let quoteparms=%quote(&parms);
      %do i = 1 %to &nb;
        %let effectx = %scan(&quoteparms,&i,=);
        %let ipos=0;
         %do %while(%length(%scan(%quote(&effectx),&ipos+1,' '))); 
          %let ipos = %eval(&ipos + 1);
          %let effect = %scan(%quote(&effectx),&ipos,' ');
         %end;

         %let effectx = %scan(&quoteparms,&i+1, =);
         %let ipos=0;
          %let parmvalues =;
          %do %while(%length(%scan(%quote(&effectx),&ipos+2,' '))); 
               %let ipos = %eval(&ipos + 1);
               %let parmvalues = &parmvalues %scan(%quote(&effectx),&ipos,' ');
             %end;
           %if %scan(&indadjset,&i) %then %do;
               %if (&i=&nb) %then  %let parmvalues=%scan(&quoteparms,&i+1,=);
               %let adjparms = &adjparms &effect = &parmvalues;
           %end;
           if _b{&i} =. then &effect= %scan(%quote(&parmvalues),1,' ' %str(,));
      %end;
   %end;
   %if (&nu) %then %do;
      array _u{&nu} &random;
      %do i = 1 %to &nu;
         if _u{&i} =. then %scan(&random,&i,' ') = 0;
      %end;
   %end;
   _one = 1;
run;

 /*---grab all specified variables---*/

/*     xvar specifies list of time independent variables  */
/*     tvar specifies list of time dependent variables    */
/*     dtvar identifies all vars except xvar              */  
 
%let xvar=_one &subject &xvar;

proc contents data=&outdata (drop=&fixed &random &xvar)
         noprint  out=_contentsoutdata(keep=name);
run;
  
data _null_;
   set _contentsoutdata end=last;
   len = length(name)+1;
   lenc + len;
   if last then call symput('lenc',lenc+1);
run;

data _null_;
   set _contentsoutdata end=last;
   length dtvar $ &lenc;   
   retain dtvar '';
   dtvar = trim(dtvar) ||' '|| left(name);
   if last then call symput('dtvar',dtvar);
run;

%if %length(&tvar)=0 %then %let tvar=&dtvar;
%let dtallvar=&xvar &tvar;

%if %index(&options,DEBUG) %then %do;
   %put nb =&nb;
   %put nu =&nu;
   %put model = &model;
   %put derivatives = &derivs;
   %put fixed = &fixed;
   %put d_fixed = &d_fixed;
   %put random = &random;
   %put d_random = &d_random;
   %put options =&options;
   %put expand  =&expand;
   %put expandu =&expandu;
   %put expandn =&expandn;
   %put adjset  =&adjset;
   %put numadjset=&numadjset;
   %put indadjset=&indadjset;
   %put adjloop=&adjloop;
   %put notadjset=&notadjset;
   %put adjparms =&adjparms;
   %put cumloc_no=&cumloc_no;
   %put cummacro_debug=&cummacro_debug;
%end;


%if %index(&options,DEBUG) & %length(&nlmem)>0  %then %do;
   %put subject  = &subject;
   %put next     = &next;
   %put imlnlin  = &imlnlin;
   %put imlnlinw = &imlnlinw;
   %put dtmplt   = &dtmplt;
   %put dbtmplt  = &dbtmplt;
   %put dutmplt  = &dutmplt;
   %put next    = &next;
   %put read    = &read_all;
   %put xvar    = &xvar;
   %put tvar    = &tvar;
   %put dtvar = &dtvar;
%end;
%mend init;


 /*------------------------------------------------------*/
 /*                                                      */
 /*   %numder                                            */
 /*   numerical derivatives                              */
 /*                                                      */
 /*------------------------------------------------------*/

%macro numder;
   %if (&nb) %then %do;
      do _i = 1 to &nb;
         if _db{_i} = . then do;
           _b0 = _b{_i};
           _tol = &tol*(1+abs(_b0));
           _b{_i} = _b0 - _tol;
           &model
           _predl = predv;
           _b{_i} = _b0 + _tol;
           &model
           _predu = predv;
           _db{_i} = (_predu - _predl)/2/_tol; 
           _b{_i} = _b0;
         end;
      end;
     drop _b0;
   %end;
   %if (&nu) & not (&nlin) %then %do;
      do _i = 1 to &nu;
         if _du{_i} = . then do;
           _u0 = _u{_i};
           _tol = &tol*(1+abs(_u0));
           _u{_i} = _u0 - _tol;
           &model
           _predl = predv;
           _u{_i} = _u0 + _tol;
           &model
           _predu = predv;
           _du{_i} = (_predu - _predl)/2/_tol; 
           _u{_i} = _u0;
         end;
      end;
     drop _u0;
   %end;

   %if %index(&options,DEBUG) %then &model;

%mend numder;


 /*------------------------------------------------------*/
 /*                                                      */
 /*   %nlin                                              */
 /*   PROC NLIN call                                     */
 /*                                                      */
 /*------------------------------------------------------*/
%macro nlin;

%if (not %index(&options,NOPRINT)) & (&nu) %then %do;
   %put Calling PROC NLIN to initialize.;
%end;

proc nlin data=&outdata(keep=&dtallvar &notadjset &random) &_printn_ 
   %unquote(&nlinopt); 
   array _b{&nb} &fixed;
   array _db{&nb} &d_fixed;
   array _derb{&nb} &d_fixedn;
   %if (&nu) %then %do;
      array _u{&nu} &random;
      array _du{&nu} &d_random;
   %end;

   /*---set starting values---*/
   parms %unquote(&adjparms);

   /*---compute the nonlinear function and its derivatives---*/
   &modinit
   &model
   model %unquote(&response) = predv;
   &weight
   &derivs
   %numder
   do i = 1 to &nb;
      _derb{i} = _db{i};
   end;

   output out=&outdata parms=&adjset;
   %unquote(&nlinstmt);
   drop i;
run;

%if %index(&options,DEBUG) %then %do;
   %user_debug(888);
   %put Ignore warning message(s) listed above of the type:; 
   %put Variable X was not found on DATA file.;      
%end;

%mend nlin;


 /*----------------------------------------------------------*/
 /*                                                          */
 /*   %pseudoder                                             */
 /*   construct pseudo data and derivatives                  */
 /*                                                          */
 /*----------------------------------------------------------*/
%macro pseudoder;

data &outdata;
   set &outdata(keep=&dtallvar &fixed &random) end=_last;
   %if (&nb) %then %do;
      array _b{&nb} &fixed;
      array _db{&nb} &d_fixed;
   %end;
   %if (&nu) & not (&nlin) %then %do;
      array _u{&nu} &random;
      array _du{&nu} &d_random;
   %end;
   &modinit
   &model
   &pseudoresp = %unquote(&response) - predv;
   %if %index(&options,DEBUG) %then %do;
      _pseudoy = &pseudoresp;
   %end;
   &weight
   &derivs
   %numder
   if (&pseudoresp ne .) then do;
      %if (&nb) %then %do;
         do _i = &adjloop;    /* It was:  1 to nb */
            &pseudoresp = &pseudoresp + _db{_i}*_b{_i};
         end;
      %end;
      %if (&nu) & (&expandu=EBLUP) & not (&nlin) %then %do;
         %do i = 1 %to &nu;
            %if (%scan(&expandn,&i)=1) %then
            &pseudoresp = &pseudoresp + _du{&i}*_u{&i};;
         %end;
      %end;
   end;
   if (_error_ = 1) then do;
      call symput('_error',left(_error_));
      stop;
   end;
   drop _i;
run;
%let _trace=  %str(Pseudoder macro: outdata -> outdata) ;  
%if %length(&postpd)>0 %then %do;
%let _trace=&_trace : POSTPD=&postpd executed;
%&postpd;
%end;
%user_debug(1);

%mend pseudoder;

 /*----------------------------------------------------------*/
 /*                                                          */
 /*   %mixed                                                 */
 /*   PROC MIXED step                                        */
 /*                                                          */
 /*----------------------------------------------------------*/
%macro mixed;

%let ncall = %eval(&ncall + 1);
%if not %index(&options,NOPRINT) %then
   %put %str(   )PROC MIXED call &ncall;

%let _trace=MIXED: %quote(%mtstmts);
%user_debug(1);

%let _trace=MIXED;
proc mixed data=&outdata %unquote(&procopt) ;
   %mtstmts;
   %if (&nu) & not (&nlin) & (&ncall>0) & not %index(&options,NOPREV) & 
      not %index(&stmts,PARMS) & (%sysfunc(exist(_covsave))) %then %do;
         %let _trace = &_trace : parms / pdata=_covsave ;
         %str(parms / pdata=_covsave;) ;
      %end;
   ods output fitstatistics=_fit;
   %if (&nb) %then %do;
      ods output solutionf=_soln;
         %let _trace=  &_trace : _soln created ;  
   %end;
   %if (&nu) & not (&nlin) %then %do;
      ods output covparms=_cov;
      %let _trace=  &_trace : _cov created;  
      %if (&expandu=EBLUP) %then %do;
         ods output solutionr=_solnr;
         %let _trace=  &_trace : _solnr created;  
      %end;
   %end;
run;

data _soln;
   set _soln;
   length param_name $ 38;
   param_name =substr(upcase(effect),3,length(effect)-2);
run;
%user_debug(1);
%user_debug(921);

 /*---check for convergence---*/
%let there = 0;
data _null_;
   set _fit;
   call symput('there',1); 
run;
%if (&there = 0) %then %do;
   %if not %index(&options,NOPRINT) %then   
      %put PROC MIXED did not converge.;
   %let _error = 1;
%end;
%mend mixed;


 /*-----------------------------------------------------------*/
 /*                                                           */
 /*   %merge                                                  */
 /*   merge in new estimates of b and u                       */
 /*                                                           */
 /*-----------------------------------------------------------*/
%macro merge;

 /*---merge in b---*/
%if (&nb) %then %do;
   %user_debug(261);

   data _saveb;
    set _saveb;
    length param_name $38;
    param_name=substr(effect,3,length(effect)-2);
   run;

   proc transpose data=_saveb out=_trbeta;
      var estimate;
      id param_name;
   run;

   data _beta;
      merge _beta _trbeta;
      _one = 1;
      keep _one &fixed;
   run;
   
   data &outdata;
      merge &outdata(drop=&adjset) _beta;
      by _one;
   run;
   %let _trace = MERGE _beta: &adjset;
   %user_debug(1);
   %user_debug(263);
%end;

 /*---merge in u---*/
%if (&nu) & not (&nlin) & (&expandu=EBLUP) & (&step_notadjusted) %then %do;
   %let _trace = MERGE in u;
   /*---loop through RANDOM statements---*/
   %let nui=0;
   %let iv = 1;
   %do %while (%length(%scan(&stmts,&iv,;)));

      %let stmt = %qscan(&stmts,&iv,%str(;));
      %let first = %qscan(&stmt,1);

      %if %index(%qupcase(&first),RANDOM) %then %do;

         /*---peel off pieces---*/
         %let i = %index(&stmt,/);
         %let ir = %index(%qupcase(&stmt),RANDOM);
         %let rndeff = %substr(&stmt,%eval(&ir+7),(%eval(&i-&ir-7)));
         %let rndopt = %substr(&stmt,(%eval(&i+1))); 
         %let i = %index(%qupcase(&rndopt),SUBJECT);
         %let rndsub = %substr(&rndopt,&i);
         %let i = %index(&rndsub,=);
         %let rndsub = %qscan(%substr(&rndsub,(%eval(&i+1))),1);
         %let subvar=;
         %let i=1;
         %do %while(%length( %scan(&rndsub,&i,%str( *())) ));
            %let subvar=&subvar %scan(&rndsub,&i,%str( *()));
            %let i=%eval(&i+1);
         %end;

         /*---loop through effects in this statement---*/
         %let ive = 1;
         %do %while (%length(%scan(&rndeff,&ive,' ')));

            %let d_effect = %scan(&rndeff,&ive,' ');
            %let effect = %substr(&d_effect,3);
            %let nui=%eval(&nui+1);
          %if %scan(&expandn,&nui)=1 %then %do;
            data _eblup;
               set _saveu;
               deff = "&d_effect";
               if (upcase(effect) = upcase(deff));
               /*---convert to numeric---*/
               &effect = estimate;
               keep &subvar &effect; 
            run;

            proc sort data=_eblup;
               by &subvar;
            run;

            proc sort data=&outdata;
               by &subvar;
            run;

            data &outdata;
               merge &outdata(drop=&effect) _eblup;
               by &subvar;
            run;
            %let _trace=&_trace / &d_effect : &subvar;
            %user_debug(265);
          %end;             
          %let ive = %eval(&ive + 1);
         %end;
      %end;

      %let iv = %eval(&iv + 1);
   %end;
   %user_debug(1);
%end;

%mend merge;


 /*-------------------------------------------------------*/
 /*                                                       */
 /*    %iterate                                           */
 /*    Iteration process                                  */
 /*                                                       */
 /*-------------------------------------------------------*/

%macro step_adjust;
/* August, 2002 
   Adjusts  step size for fixed parameters so linear inequality constraints are satisfied. 
*/         

         proc transpose data=_soln(keep=param_name estimate)
             out=_solnt(drop=_name_);
         id param_name;
         quit;
          
         proc transpose data=_solnold(keep=param_name estold) out=_solnoldt(drop=_name_);
         id param_name;   
         quit;
         
         data _solnconstraints;
           set _solnt;
           %do i= 1 %to &n_neqconstraints;
           _solnconstraint=%scan(&cum_neqexpression,&i,%str(|));
           _solnconstraint_satisfied=%scan(&cum_neqconstraints,&i,%str(|));
           output;
           %end;
         run;
            
         data _solnoldconstraints;
           set _solnoldt;
           %do i= 1 %to &n_neqconstraints;
           _solnoldconstraint=%scan(&cum_neqexpression,&i,%str(|));
           output;
           %end;
         run;
   
         data _solnconstraints_combined;
         merge _solnconstraints _solnoldconstraints  _neqconstraints;
         run;

        data _soln_combined;
          merge _soln _solnold;
         run;
         %user_debug(231);
  
         %let maxi=1;
          data _solnconstraints_combined;
           set _solnconstraints_combined;
           retain maxratio 0;
         
         if not _solnconstraint_satisfied  then do;
            delta= _solnconstraint-_solnoldconstraint;
            if delta>0 then xdiff=max_value -_solnoldconstraint;
            if delta<0 then xdiff=min_value -_solnoldconstraint;
            xratio=0;
            if xdiff ne 0 then xratio = delta/xdiff;
            if xratio>=1 then maxratio=max(maxratio,xratio);
         end;
            if  maxratio > 1 then  call symput('maxi', maxratio);
         output;   
         run;
          %let _trace=  STEP_ADJUST: maxi=&maxi;  
         
          %if %eval(&maxi>1) %then %do;
           %let step_reduced=%eval(&step_reduced+1);
           %let step_notadjusted=0;
          data _soln_corrected;
           set _soln_combined;
           delta=estimate-estold;
           original_estimate=estimate;
         
           if &maxi>1 then estimate=estold+ (&divisor)*delta / &maxi;
          run;

             proc datasets lib=work nolist;   
             delete _soln 
              %if %sysfunc(exist(_solnr)) %then _solnr;
              %if %sysfunc(exist(_cov))   %then _cov;
              %if %sysfunc(exist(_fit))   %then _fit;
              ;
              quit;
             
          data _soln(keep=param_name effect estimate StdErr DF tValue Probt Alpha Lower Upper);
           set _soln_corrected;
           StdErr=.;  DF=.; tValue=.; Probt=.; Alpha=.; Lower=.; Upper=.;
          run;
          %let _trace=  &_trace: _soln _solnold -> _soln adjusted;  
          %end;
          %else %do;
           %let step_reduced=0;
           %let _trace =&_trace: _soln NOT adjusted;
          %end;
          %user_debug(1);
          %user_debug(232);
%mend step_adjust;

%macro _neqconstraints_macro;
%let _trace=&_trace : _neqconstraints_macro executed;
/* Initialize linear inequality constraints */
data _neqconstraints;
   set _beta;
   length neqconstraint $ 200;
   %let iv=1;
   %let _neqconstraints=%quote(%&neq_constraints);
   %do %while (%length(%scan(&_neqconstraints,&iv,;)));
   %let stmt = %qscan(&_neqconstraints,&iv,%str(;));
   neqconstraint="&stmt";
   if not (&stmt) then  do;
      put "FATAL: Parameters &fixed=" &fixed " do not satisfy constraint:" neqconstraint;
      call symput('_error','1');
      end;
      fixed="&fixed";
   output;
   %let iv=%eval(&iv+1);
   %end;
   drop &fixed;
run;


data _neqconstraints;
  length fixedi $ 40;
  length word $ 80;
  length expression $ 120;
  length neqconstraint2 ctemp $ 200;
  length cmin_value cmax_value $ 40;
  set _neqconstraints end=last;
      maxc='1E20';
      wordn=0;
      neqs1=index(neqconstraint,'>');
      neqs2=index(neqconstraint,'<');
      do while (scan(neqconstraint,wordn+1,'<>') ne ' ');
       wordn=wordn+1;
       word=scan(upcase(neqconstraint),wordn,'<>');
         do i=1 to &nb;
          fixedi=scan(upcase(fixed),i);
          idx=index(word,fixedi); 

            if (idx) then do; 
              iexpression=wordn;  
              *put '->wordn,word,fixedi,iexpression,idx' wordn word fixedi iexpression idx;
            end;
          end;
      end;

        expression_type=wordn;
        if neqs1*neqs2>0 then do;
         put "FATAL: Inequality signs conflict in " neqconstraint;
         call symput('_error','1');
        end;
        neqconstraint2=neqconstraint;
        if wordn=2 then do;
          if iexpression=1 then do; /* expression ? # */
             word = scan(upcase(neqconstraint),2,'<>');
             expression = scan(upcase(neqconstraint),1,'<>');
             if neqs1 then neqconstraint2=maxc ||' >'|| left(trim(expression)) || '>'||left(trim(word));   /* expression >?= # */ 
             if neqs2 then neqconstraint2='-'||maxc||'<'|| left(trim(expression)) || '<'||left(trim(word));   /* expression <?= # */ 
           end;
          if iexpression=2 then do; /* # ? expression */
             word = scan(upcase(neqconstraint),1,'<>');
             expression = scan(upcase(neqconstraint),2,'<>');
             if neqs1 then neqconstraint2=left(trim(word))|| '>'||left(trim(expression)) || '>-' || maxc;   /* # >?= expression */ 
             if neqs2 then neqconstraint2=left(trim(word))|| '<'||left(trim(expression)) || '<' || maxc;   /* # <?= expression */ 
          end;
        end;

      if neqs1 then do;
      ctemp = neqconstraint2;
      word= scan(ctemp,3,'>');
      neqconstraint2=left(trim(compress(word,'=')))||'<';
      if index(word,'=') then neqconstraint2=left(trim(neqconstraint2))||'=';
      expression=scan(ctemp,2,'>');
      neqconstraint2= left(trim(neqconstraint2))||left(trim(compress(expression,'=')))||'<';
      if index(expression,'=') then neqconstraint2=left(trim(neqconstraint2))||'=';
      word= scan(ctemp,1,'>');
      neqconstraint2=left(trim(neqconstraint2))||left(trim(compress(word,'=')));
      end;
      expression=compress(scan(neqconstraint2,2,'<'),'=');
      cmin_value=compress(scan(neqconstraint2,1,'<'),'=');
      cmax_value=compress(scan(neqconstraint2,3,'<'),'=');
       *length cum_constraints cum_minvalue cum_maxvalue cum_expression $ 1000;
      len_neqconstraints+(length(neqconstraint2)+2);
      len_neqexpression+(length(expression)+2);
      len_neqcmin_value+(length(cmin_value)+2);
      len_neqcmax_value+(length(cmax_value)+2);

      if last then do;
         call symput('n_neqconstraints',left(_n_));
         call symput('len_neqconstraints', left(len_neqconstraints));
         call symput('len_neqexpression',  left(len_neqexpression));
         call symput('len_neqcmin_value',  left(len_neqcmin_value));
         call symput('len_neqcmax_value',  left(len_neqcmax_value));

      end;

      drop fixedi i word idx ctemp;
run;

data _neqconstraints;
   set  _neqconstraints end=last;
   length cum_neqconstraints $ &len_neqconstraints;
   length cum_neqexpression $ &len_neqexpression;
   length cum_neqcmin_value $ &len_neqcmin_value;
   length cum_neqcmax_value $ &len_neqcmax_value;
   retain cum_neqconstraints cum_neqexpression cum_neqcmin_value cum_neqcmax_value;
   cum_neqconstraints=left(trim(cum_neqconstraints)) || left(trim(neqconstraint2)) || '|';
   cum_neqexpression =left(trim(cum_neqexpression))  || left(trim(expression))     || '|';
   cum_neqcmin_value =left(trim(cum_neqcmin_value))  || left(trim(cmin_value))     || '|';
   cum_neqcmax_value =left(trim(cum_neqcmax_value))  || left(trim(cmax_value))     || '|';
  
   if last then do;
     call symput('cum_neqconstraints',cum_neqconstraints);
     call symput('cum_neqexpression', cum_neqexpression);
     call symput('cum_neqcmin_value', cum_neqcmin_value);
     call symput('cum_neqcmax_value', cum_neqcmax_value);
   end;
   drop cum_neqconstraints cum_neqexpression cum_neqcmin_value cum_neqcmax_value;
run;


%if %index(&options,DEBUG) %then %do;
%put n_neqconstraints=&n_neqconstraints;
%put len_neqconstraints=&len_neqconstraints;
%put len_neqexpression=&len_neqexpression;
%put len_neqcmin_value=&len_neqcmin_value;
%put len_neqcmax_value=&len_neqcmax_value;
%put cum_neqconstraints=&cum_neqconstraints;
%put cum_neqexpression=&cum_neqexpression;
%put cum_neqcmin_value=&cum_neqcmin_value;
%put cum_neqcmax_value=&cum_neqcmax_value;
%end;

data _limit_values;
           %do i= 1 %to &n_neqconstraints;
           min_value = %scan(&cum_neqcmin_value,&i,%str(|));    
           max_value = %scan(&cum_neqcmax_value,&i,%str(|));    
           output;
           %end;
run;

data _neqconstraints;
   merge _neqconstraints _limit_values;
run;

%mend _neqconstraints_macro;



%macro iterate;

 /*---initial data set for iterinfo---*/
%let _trace=ITERATE;
%if (&saveb = 0 or &nlinmethod ne 1) %then %do;
  data _beta;
   set &outdata(obs=1);
   %if (&nb) %then %do;
      keep &fixed;
   %end;
  run;
  %let _trace=&_trace : _beta extracted from _nlinmix;

  %if ((&nb) & %length(&neq_constraints)) %then %do;
   %_neqconstraints_macro;
   %if (&_error = 1) %then %do;
      %let _trace=&_trace : GOTO finish_iter;
      %user_debug(1);
      %goto finish_iter;
   %end;

  %user_debug(101); /* DATASETS:  &outdata _beta  */ 
  %end;

  /* Initialize _solnold for step_adjust */

  proc transpose data=_beta(keep=&adjset) out=_solnold (rename=(col1=estold));
  quit;

 data _solnold;
  set _solnold;
  length param_name $38;
  param_name=upcase(_name_);
  drop _name_;
 run;
 %let _trace=&_trace : transpose(_beta) -> _solnold;

%user_debug(105);  
%end;
%user_debug(1);

 /*---initial macro variables---*/
%let crit = .;
%let conv = 0;
%let ni = 0;
%let ncall = -1;
%let step_notadjusted=0;

%if (&nlin) %then %put NLIN step using calls to PROC MIXED;
 
 /*---iterate until convergence---*/
%do %while(&ni <= &maxiter);

   %let _trace=  ITERATION ni=&ni. nlin=&nlin;  
   %if (&ni = 0) and not %index(&options,NOPRINT) %then
      %put Iteratively calling PROC MIXED.;
   
   %if (&switch > 0) and (&ni = &switch) and not (&nlin) %then %do;
      %let expandu = EBLUP;
      %let ni = 0;
      %let switch = 0;
      %put Switching from EXPAND=ZERO to EXPAND=EBLUP;
      %let _trace= &_trace : Switching from EXPAND=0 -> EBLUP;  
   %end;

   /*---save estimates---*/
   %if (&ni ne 0) and (&nu) and (%sysfunc(exist(_covsave))) and not (&nlin) %then %do;
      data _covold;
         set _covsave;
         estold = estimate;
         keep estold;
      run;
      %let _trace=&_trace : _covsave -> _covold ;  
   %end;

   %if (&ni ne 0 or (&ni=0 & &nlinmethod=1 & &saveb)) and (&nb) %then %do;
      data _solnold;
         set _saveb;
         estold = estimate;
         keep estold param_name;
      run;
      %let _trace= &_trace : _saveb -> _solnold ;  
   %end;
   %user_debug(1);
   %user_debug(201);  

   /*---set up pseudo data and compute first step---*/
   %pseudoder
   %user_debug(215);

   %if %index(&options,DEBUG) %then %do;
    %user_debug(880);
   %end;

   %if (&_error = 1) %then %goto finish;
   %mixed
   %user_debug(220);
   %if (&_error = 1) %then %goto finish;

   %let step_notadjusted=1;
   %if ((&nb) & %length(&neq_constraints)) %then %step_adjust;
   %if (&step_notadjusted=0) %then  
       %put %str(     Stepsize for beta parameters reduced (r=&maxi)); 

   %if %length(&postmixed)>0 %then %do;
     %let _trace=POSTMIXED:  &postmixed executed;
     %&postmixed;
   %end;
 
   %if (&nu) & not (&nlin) & (&step_notadjusted) %then %do;
      %let _trace= _&_trace : %str(_cov -> _covsave) ;  
      data _covsave;
         set _cov;
      run;
   %end;
      %user_debug(1);
      %user_debug(233);

   %if (&_error = 1) %then %goto finish;
   /*---check for convergence---*/

   %let _trace=ni ne 0?: ni=&ni:;
   %if (&ni ne 0) %then %do;
      %let crit = 0;

      %if (&nb) %then %do;
         data _null_;
            merge _soln _solnold end=last; 
            retain cr &crit;
            crit = abs(estimate-estold)/
               max(abs(estold),max(abs(estimate),1));
            %if %index(&options,WORST) %then %do;
               if (crit > cr) then do;
                  call symput('worst',left(estold));
               end;
            %end;
            cr = max(cr,crit);
            if last then do;
               call symput('crit',left(cr));
            end;
         run;
         %let _trace=&_trace : _soln and _solnold merged;
      %end;


      %let cov_crit=0;
      %if (&nu) %then %do;
        %if (%sysfunc(exist(_covold))) & (%sysfunc(exist(_cov))) %then %do;
         %let _trace=&_trace _covold exists:;
         %let cov_crit=1;
          data _null_;
            merge _cov _covold end=last;
            retain cr &crit;
            crit = abs(estimate-estold)/
               max(abs(estold),max(abs(estimate),1));
            %if %index(&options,WORST) %then %do;
               if (crit > cr) then do;
                  call symput('worst',left(estold));
               end;
            %end;
            cr = max(cr,crit);
            if last then do;
               call symput('crit',left(cr));
            end;
         run;
         %let _trace=&_trace _cov and _covold merged:;
         %user_debug(235);  
       %end;
      %end;

      %if ((&nlin) or (&cov_crit)) %then %do;
        %let temp= cr < &convrg;
        data _crit_;
         cr = &crit;
         if &temp then call symput('conv',left(1));
        run;
    
         %user_debug(240);  
         %let _trace=&_trace : _crit_; 
      %end;   
      %else %do;
         %let _trace=&_trace NOT; 
         %let crit= &crit  (COVP not used);
      %end;
         %let _trace=&_trace created;
     
   %end;
   %else %if %index(&options,PRINTFIRST) & not (&nlin) %then %do;
      %let _trace=&_trace : ni eq 0;
      %let ods_status=exclude all;
      ods &ods_status;
   %end;
   %user_debug(1);  


   %if not %index(&options,NOPRINT) %then %do;
      %if (&step_notadjusted or &ni=0) %then %do;
        %put %str(   );
        %if (&nlin) %then
              %put NLIN1 iteration = &ni;
        %else %put iteration = &ni;

        %put convergence criterion = &crit;;
        %getbc;

        %put &bstr &cstr;
        %if %index(&options,WORST) %then %do;
          %put worst = &worst;
        %end;
      %end;
       
  %end;

   %let _trace=CONV CHECK;
   %if (&conv = 1) & (&step_notadjusted) %then %do;
      %let maxiter = -1;
      %let _trace = &_trace : maxiter -> -1;
   %end;
   %else %do;
    %if %length(&append) %then %apploop(&appnd,append);
    %if (&step_notadjusted or &ni=0) %then %do;
      %let ni = %eval(&ni + 1);
      %let _trace = &_trace : ni <- ni+1 =&ni;
    %end;
   %end;
   %user_debug(1);
   %user_debug(250);

   /*---save step---*/
   %let _trace = SAVE STEP;
   %if (&nb)  %then %do;
      data _saveb;
         set _soln;
      run;
      %if (&ni=1) %then %let saveb=1;
   %let _trace=&_trace : _soln -> _saveb;
   %end;
   %if (&nu) & (&expandu=EBLUP) & (&step_notadjusted) & not (&nlin) %then %do;
      data _saveu;
         set _solnr;
      run;
   %let _trace=&_trace : _solnr -> _saveu;
   %end;
   %user_debug(1);
   %user_debug(260);
   %merge;

   /*---get rid of fitting data set---*/
   %if (&there = 1) %then %do;
      proc datasets lib=work nolist;
         delete _fit;
      quit;
   %end;
   %user_debug(280);
%end;

 /*---turn on printing and options---*/
  %finish:

  
  %if (&nlin) %then %do;
    %let _trace=GOTO finish_iter;
    %user_debug(1);
    %goto finish_iter;
  %end;

  %let _trace=Turn on printing and options;
  %let ods_status=select all;
  ods &ods_status;
  %let _printn_ = ;
  %let niter = &ni;
  %if not %index(&options,NOTES) %then %do;
   %let notes=on;
   options notes date number;
  %end;
  %user_debug(1);
  %user_debug(300);

%if not %index(&options,NOPRINT) %then %do;
   %if (&conv = 1 & not (&nlin) ) %then %do;
      %put NLINMIX convergence criteria met.;
      /*---compute final results---*/
      %user_debug(310);
      %pseudoder
      %user_debug(315);
      %if %index(&options,DEBUG) %then %do;
        %user_debug(881); 
      %end;
      %if (&expandu = ZERO) %then %let expandu = EBLUP;
      %mixed
      %user_debug(320);
      %if %length(&append) %then %apploop(&appnd,append);
   %end;
   %else %put NLINMIX did not converge.;

   %user_debug(350);
%end;

%finish_iter:
%let _trace=FINISH ITERATE: NLIN1=&nlin;
%user_debug(1);
%mend iterate;


 /*------------------------------------------------------------*/
 /*                                                            */
 /*    %baseinfo                                               */
 /*    Print basic information about the macro                 */
 /*                                                            */
 /*------------------------------------------------------------*/
%macro baseinfo;

%put;
%put %str(                          The NLINMIX Macro);
%put;
%put %str(           Data Set                     : &data);
%put %str(           Response                     : &response);
%if (&nb) %then 
%put %str(           Fixed-Effect Parameters      : &adjset);
%if %length(&notadjset)>0 %then
%put %str(           Not adjustable parameters:   : &notadjset) ;
%if (&nu) %then
%put %str(           Random-Effect Parameters     : &random);
%if (&nu) %then
%put %str(           Expansion Point              : &expand);
%if not %index(&options,SKIPNLIN) %then %do;
  %if &nlinmethod=1 %then 
%put %str(           NLIN step(default method)    : Calls to PROC MIXED);
  %if &nlinmethod=2 %then %do; 
  %if %length(&nlmem)=0 %then
%put %str(           NLIN step(method=2)          : PROC NLIN);
  %if %length(&nlmem)>0 %then
%put %str(           NLIN step(method=2)          : &imlnlin macros);
  %end;
%end;

%if (&nb & %length(&neq_constraints)) %then 
%put %str(           Macro with LinNeq constraints: &neq_constraints);

%put;

%mend baseinfo;

 /*----------------------------------------------------------*/
 /*                                                          */
 /*    %getbc                                                */
 /*    load current estimate of b into macro variables       */
 /*                                                          */
 /*----------------------------------------------------------*/
%macro getbc;

%let bstr = ;
%let cstr = ;

%let _trace=GETBC;
%if (&nb) %then %do;
   data _beta;
      set _beta;
      %do i = 1 %to &nb;
         call symput("bb&i",left(%scan(&fixed,&i,' ')));
      %end;
   run;

   %do i = 1 %to &nb;
      %let bstr = &bstr %scan(&fixed,&i,' ')=&&bb&i;
   %end;
   %let _trace=&_trace beta:  &fixed extracted from _beta;
%end;

%if (&nu) & not (&nlin)  & (&step_notadjusted) %then %do;
   data _null_;
      set _cov nobs=count;
      call symput('ncov',left(put(count,8.)));
   run;

   data _null_;
      set _cov;
      %do i = 1 %to &ncov;
         if (_n_ = &i) then do;
            call symput("cc&i",left(estimate));
         end;
      %end;
   run;

   %do i = 1 %to &ncov;
      %let cstr = &cstr COVP&i=&&cc&i;
   %end;
   %let _trace=&_trace : _cov extracted from _cov;
%end;
%user_debug(1);
%mend getbc;


/*---appending macros---*/
%macro appnd;
   %let tmp=&append;
   %let lentmp=%length(&tmp);
   %let i=%index (&tmp,->);
   %let line=%substr(&tmp,1,%eval(&i-1));
   %let appnd=%str(&line \);
   
   %next:
   %let tmp=%substr(&tmp,%eval(&i+2),%eval(&lentmp-&i-1));
   %let lentmp=%length(&tmp);
   %let i=%index(&tmp,->);
   
   %if not &i %then %goto fin;
   %let line=%substr(&tmp,1,%eval(&i-1));
   %let lenline=%length(&line);
   %let out=%scan(&line,1,%str( ));
   %let appnd=&appnd &out |;
   
   %if &i %then %do;
   %let lenout=%length(&out);
   %let in =%substr(&line,%eval(&lenout+1),%eval(&lenline-&lenout));
   %let appnd=&appnd &in \;
   %goto next;
   %end;
   %fin:
   %let appnd=&appnd &tmp;
%mend appnd;

%macro apploop(appnd,doappnd);
   %let _trace =Apploop: &appnd : &doappnd;
   %user_debug(1);
   %let i=0;
   %do %while (%length(%scan(&appnd,%eval(&i+1),|)));
      %let i=%eval(&i+1);
      %let tmp=%scan(&appnd,&i,|);
      %let in =%scan(&tmp,1,\); %* dsn with options;
      %let out=%scan(&tmp,2,\); %* dsn without options;

      %let in1 =%scan(&in,1,'('); %*dsn without options;
      %let out1=%scan(&out,1,'(');
      
      %let lbin=work;
      %let dtin=&in1;
      %if  %index(&in1,.) %then %do;
         %let lbin=%scan(&in1,1,.);
         %let dtin=%scan(&in1,2,.);
      %end;

      %let lbout=work;
      %let dtout=&out1;
      %if  %index(&out1,.) %then %do;
         %let lbout=%scan(&out1,1,.);
         %let dtout=%scan(&out1,2,.);
      %end;
      %&doappnd;
   %end;
%mend apploop;

%macro appdel;
 %if %sysfunc(exist(&lbout..&dtout)) %then %do;
   proc datasets lib=&lbout nolist;
      delete &dtout;
   quit;
 %end;
%mend appdel;

%macro append;
 %if %sysfunc(exist(&in)) %then %do;
   data _tempout;
      set  &in;
      _switch=&switch;
      _ni = &ni;
      _iter=_ni;
      if _switch=0 then _iter=_ni+&switchadd;
      _call = &ncall+1;
      _step_adjust=&maxi;
      _nlin1 =&nlin; 
      _reduced=&step_reduced;
   run;
   
   proc append base=&out force data=_tempout;
   quit;
 %end;
%mend append;



 /*------------------------------------------------------------*/
 /*                                                            */
 /*    %nlinmix()                                              */
 /*    the main macro                                          */
 /*                                                            */
 /*------------------------------------------------------------*/

%macro nlinmix_init;

   /*---default data set---*/
   %if %bquote(&data)= %then %let data=&syslast;

   %let options = %qupcase(&options);

   %let switchadd=&switch;
   %let maxi=1;

   /*---initialize---*/
   %let notes=on;
   %if not %index(&options,NOTES) %then %do;
      %let notes=off;
      options nonotes nodate nonumber;
   %end;
   proc datasets nolist;
      delete alldebug _crit_ _soln _solnr;
   quit;

%mend nlinmix_init;


%macro nlinmix(data=,modinit=,model=,weight=,derivs=, tol=1e-5,
   parms=,stmts=,expand=zero,converge=1e-8 1e-1 ,maxit=30,switch=0,
   append=,procopt=,nlinopt=,nlinstmt=,outdata=_nlinmix,options=,
   paramfixed =,divisor=0.9, postpd=, postmixed=, nlinmethod=1,
   user_debug=, neq_constraints=,

   /*---NLMEM macro specific arguments---*/
   nlmem=, imlfun=,  globvar=, retain=,
   next=1, xvar=, tvar=, savevar=, dtmplt=, debug=,
   imlnlin= default_nlin default_objf default_initobjf
);

 /*---check for mandatories---*/
%if %bquote(&model)= %then %let missing = MODEL=;
%else %if %bquote(&stmts)= %then %let missing = STMTS=;
%else %let missing =;
%if %length(&missing) %then %do;
   %put ERROR: The NLINMIX &missing argument is not present.;
%end;
%else %do;

   /*---global variables---*/
   %global _printn_ _error;

   /*---local variables---*/

   %local nb nu no ni ncall response pseudoresp fixed d_fixed d_fixedn
      random d_random crit bstr cstr there appnd dtvar nlmemlocal
      dtallvar subject xvar tvar expandn expandu switchadd maxi
      numadjset indadjset adjloop adjset notadjset adjparms
      macro_debug notes step_notadjusted ods_status
      len_neqconstraints n_neqconstraints len_neqexpression len_neqcmin_value len_neqcmax_value
      cum_neqconstraints cum_neqexpression cum_neqcmin_value cum_neqcmax_value
      nmacro_debug cummacro_debug cumloc_no
      maxiter nlin _trace convrg stmt1 stmt2 stmt3 saveb step_reduced
;

%let stmts0=;
%do i=1 %to 3;
 %let stmt&i=%scan(&stmts,&i,%str(|));
 %let stmts0=&stmts0&&stmt&i;
%end;
%let stmts=&stmts0;
%let step_reduced=0;


%if %length(&nlmem)=0 %then %do;
   
   %let tvar=; 
   %let xvar=;
   %let subject=;
   %nlinmix_init
   %init
%end;
%else %do;
  %local nlmem_local;

  /* Include nlmem file with a set of NLMEM macros */
      %include &nlmem;
      %local &nlmemlocal;
      %nlinmix_init
      %nlmem_init
  %end;

   %let nlin=0;
   %if %index(&options,DEBUG) %then %do;
   %pseudoder;

   %user_debug(889); 
   %end;

   %if %length(&append) %then %do;
      %appnd;
      %apploop(&appnd,appdel);
   %end;

   %if %index(&options,PRINTALL) or %index(&options,PRINTFIRST) 
      %then %do;
      %let ods_status=select all;
      %let _printn_ = ;
   %end;
   %else %do;
      %let ods_status=exclude all;
      %let _printn_ = noprint;
   %end;
      ods &ods_status; 
   
   /*---print basic macro information---*/
   %if not %index(&options,NOPRINT) %then %baseinfo;

   /*---run PROC NLIN to get initial estimates---*/
   %let saveb=0;
   %if (&nb) & not %index(&options,SKIPNLIN) %then %do;
      %if not (&nu) %then %let _printn_ = ;
      
      %if &nlinmethod=2 %then %nlin;

      %if &nlinmethod=1 %then %do;
      %let nlin=1;
      %let convrg=%scan(&converge,2,' ');
      %if %length(&convrg)=0 %then %let convrg=%scan(&converge,1,' ');
      %let maxiter=%scan(&maxit,2,' ');
      %if %length(&maxiter)=0 %then %let maxiter=%scan(&maxit,1,' ');
      %iterate;

      proc datasets lib=work;
      delete _fit_
       %if %sysfunc(exist(_cov)) % then _cov;
       %if %sysfunc(exist(_covsave)) % then _covsave;
       %if %sysfunc(exist(_covold)) % then _covold;
            quit;
       %end;
   %end;

   /*---iterate PROC MIXED calls until convergence---*/
   %let maxiter=&maxit;
   %let nlin=0;
   %let convrg=%scan(&converge,1,' ');
   %let maxiter=%scan(&maxit,1,' ');
   %iterate;
%end;

%mend nlinmix;


/****** Debugging macros ******/

   options nonotes nodate nonumber;
   data default_debug;
   length macro_debug $ 80;
   input loc_no macro_debug;
   cards;
   880 debug(&outdata,&ni)
   881 debug(&outdata,&ni+0.1)
   888 debug(&outdata,-8)
   889 debug(&outdata,-9)
   ;
   run;


%macro user_debug(loc_no);

%if &notes=on %then options nonotes nodate nonumber;;
%let ods_status_save=&ods_status;
ods select all;
%let _printn_status=&_printn_;
%let _printn_ = ;
%let macro_debug=;
%if (&nmacro_debug) %then %do;
%do idebug=1 %to &nmacro_debug;
%let iloc_no= %scan(&cumloc_no,&idebug,' ');
%if (&iloc_no=&loc_no) %then %let macro_debug= %scan(&cummacro_debug,&idebug,!);
%end;
%end;
%let temp_debug=&macro_debug;

%if %length(&macro_debug) %then %do;
  %put  %str( -> DEBUG: &temp_debug executed at: &loc_no);
  %&temp_debug;;
%end;
%if &notes=on %then options notes date number;;
%let _printn_ = &_printn_status;
ods &ods_status; 
%mend user_debug;

%macro user_debugall;

data user_debugall;
 length macro_debug $ 80;
 length loc_no 8;
 if _n_<0;
run;


%if %length(&user_debug) %then %do;
  proc append base=user_debugall
              data=&user_debug force;
  quit;
  %end;

%if %index(&options,DEBUG) %then %do;
  proc append base=user_debugall
              data=default_debug;
  quit;
  %end;


%let lenloc_no=1;
%let lenmacro_debug=1;
%let nmacro_debug=0;
data _null_;
 set user_debugall end=last;
 lentemp = length(left(trim(loc_no)))+1;
 lenloc_no + lentemp;
 lentemp = length(macro_debug)+1;
 lenmacro_debug + lentemp; 
 
 if last then do;
    call symput('nmacro_debug',left(trim(_n_)));
    call symput('lenloc_no',left(lenloc_no+1));
    call symput('lenmacro_debug',left(lenmacro_debug+3));
 end;
run;

%let cumloc_no=;
%let cummacro_debug=;
data _null_;
 set user_debugall end=last;
 length cumloc_no  $ &lenloc_no;
 length cummacro_debug $ &lenmacro_debug;
 length ctemp $ 200;
 retain cumloc_no  '';
 retain cummacro_debug "'";
 ctemp=left(loc_no);
 cumloc_no      = trim(cumloc_no)||' '||ctemp;
 ctemp=left(macro_debug);
 cummacro_debug = trim(cummacro_debug)||'!'||trim(ctemp);
 if last then do;
   call symput('cumloc_no',left(cumloc_no));
   call symput('cummacro_debug', trim(cummacro_debug)||"'");
 end;
run;

proc print data=user_debugall;
var loc_no macro_debug;
run;

%let cummacro_debug=%qscan(&cummacro_debug,1,"'");
%mend user_debugall;


%macro debug(data,ni);
   data debug;
      set &data;
      _switch=&switch;
      _ni = &ni;
      _iter=_ni;
      _nlin1=&nlin;
      _reduced=&step_reduced;
      if _switch=0 then _iter=_ni+&switchadd;
   run;

   *proc print data=debug(obs=30);
   run;

   proc append base=alldebug data=debug force;
   run;
%mend debug;

