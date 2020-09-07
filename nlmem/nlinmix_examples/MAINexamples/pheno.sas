
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
/*----------------------------------------------------------------*/

libname lib '.';

*options date notes number;


%let path=../..;
filename pheno_dt "&path/datasets/pheno_dt.sas";
filename nlinmix "&path/nlmm800.sas";
%include pheno_dt;
%include  nlinmix / nosource;

/* Example 1. Numerical derivatives   */
/* tmp array has been used instead of _temporary_ */

Title " Pheno data. Example 1. Numerical  derivatives.";

%nlinmix(data=pheno,
         stmts   =  %str(class indiv;
             model pseudo_conc = d_beta1 d_beta2 d_beta3
                               / noint notest solution cl;
         random d_u1 d_u2 / subject=indiv solution cl;
         weight _weight_;
         ods output Iterhistory=iter;
                      )
,
       modinit = %str(
                      retain _tmp1 - _tmp12;
                      array  _tmp{12} _tmp1 - _tmp12;

                      if (newsub=1) then
                        do call=1 to 12;
                        _tmp{call}= 0;
                        end;

                      call=0;
                     ),
       model   =  %str(
                       call =call+1;
                       prd0= _tmp{call};

                       clear=beta1*weight*exp(u1);
                       vol = beta2*weight*(1+beta3*apgarlow)*exp(u2);
                       eterm=exp(-(time-tlag)*clear/vol);
                       predv =dose/vol +prd0*eterm;

                       _tmp{call}=predv;
                       ),

       weight  =  %str( _weight_=1/predv**2;),
       parms   =  %str(beta1=.01 beta2=1 beta3=.1),
       expand  =  zero,
       options =  debug
);
run;

