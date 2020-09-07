/*----------------------------------------------------------------*/
/*--                                                            --*/
/*-- Title:    Examples with phenobarbital data                 --*/
/*-- Filename: pheno.1.sas                                      --*/
/*-- Date:     April 4, 2000                                    --*/
/*--                                                            --*/
/*-- Ref:                                                       --*/
/*--  1.Grasela, T.H., Jr. and Donn, SM                         --*/
/*--    Neonatal population pharmacokinetics for                --*/
/*--    phenobarbital derived from routine clinical data        --*/
/*--    Dev. Pharmacol. Ther. 8(1985), 374-383                  --*/
/*--                                                            --*/
/*--  2.Littel, R.C., Miliken, G.A., Stroup, W.A., Wolfinger R.D.-*/
/*--    SAS System for Mixed Models                             --*/
/*--    SAS Institute, 1996                                     --*/
/*--                                                            --*/
/*----------------------------------------------------------------*/

libname lib '.';
         
%let path=../..;
filename pheno_dt "&path/datasets/pheno_dt.sas";
filename nlinmix  "&path/nlmm800.sas";
filename nlmem    "&path/nlmem.sas";
%include pheno_dt;
%include  nlinmix / nosource;


/* Example 1. Numerical derivatives   */

Title " Pheno data. Example 1. Numerical  derivatives.";

%macro usernlin;
 *print _nobs;
  optn= _nobs ||{5};
  tc =j(1,13,.);
  tc[1]=2;              /* Termination criterion: # iterations */
  call nlplm(rc,_bnlin,"_objf",_b,optn) jac="_dobjf" tc=tc;
  print _bnlin;
%mend usernlin;  

%macro useriml;
  %let next    =  indiv;
  %let imlnlin =  usernlin;
  %let xvar    =  weight apgarlow;
%mend  useriml;

%nlinmix(data=pheno,
         stmts   =  %str(class indiv;
         model pseudo_conc = d_beta1 d_beta2 d_beta3
                               / noint notest solution cl;
         random d_u1 d_u2 / subject=indiv solution cl;
         weight _weight_;
         ods output Iterhistory=iter;
                      ),
         globvar = predv,

         model   = %str(
              clear=_b[1]*weight*exp(_u[1]);
              vol  =_b[2]*weight*(1+_b[3]*apgarlow)*exp(_u[2]);
              predv=j(_ni,1,0);
              predv[1] = dose[1]/vol;
              do j =2 to _ni;
                eterm=exp(-(time[j]-time[j-1])*clear/vol);
                predv[j] =dose[j]/vol + predv[j-1]*eterm;
              end;
              ),

        weight  =  %str(_weight_=1/(predv#predv);),
        parms   =  %str(beta1=.01 beta2=1 beta3=.1),
        expand  =  zero,
        imlfun  =  useriml,
        nlmem   =  nlmem / nosource
);

