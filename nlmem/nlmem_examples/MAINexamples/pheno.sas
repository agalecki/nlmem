/*----------------------------------------------------------------*/
/*--                                                            --*/
/*-- Title:    Examples with phenobarbital data                 --*/
/*-- Filename: pheno.sas                                        --*/
/*--           Invokes pheno_dt.sas                             --*/
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
/*--   Execution time: ~ 3 mins                                 --*/
/*----------------------------------------------------------------*/


libname lib '.';

/* nlinmix macro. Choose appropriate version */

%let path=../..;
filename pheno_dt "&path/datasets/pheno_dt.sas";
filename nlinmix  "&path/nlmm800.sas";
filename nlmem    "&path/nlmem.sas";

%include pheno_dt;
%include  nlinmix / nosource;


/* Example 1. Numerical derivatives   */

Title " Pheno data. Example 1. Numerical  derivatives.";

%nlinmix(data=pheno(drop=tlag),
         stmts   =  %str(class indiv;
         model pseudo_conc = d_beta1 d_beta2 d_beta3
                               / noint notest solution cl;
         random d_u1 d_u2 / subject=indiv solution cl;
         weight _weight_;
         ods output Iterhistory=iter;
                        ),
         globvar = predv  _mcall predt_lag  t_lag predt,
         modinit = %str (
              if (newsub=1) then do;
                      predt_lag=j(1,12,0);
                      t_lag =0;
                  end;
                 ),

         model   = %str(
              predv_lag  =  predt_lag[,_mcall];
              clear      =  _b[1]*weight*exp(_u[1]);
              vol        =  _b[2]*weight*(1+_b[3]*apgarlow)*exp(_u[2]);
              eterm      =  exp(-(time-t_lag)*clear/vol);
              predv      =  dose/vol + predv_lag*eterm;
              ),

       weight  =  %str( _weight_=1/predv**2;),
       retain  =  %str(predt_lag  = predt;
                        t_lag     = time;),

       parms   =  %str(beta1=.01 beta2=1 beta3=.1),
       expand  =  zero,
       nlmem   =  nlmem
);
run;


