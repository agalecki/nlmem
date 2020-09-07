/*---Examples for NLINMIX macro---*/

%let path=../..;
filename pigs_dt "&path/datasets/pigs_dt.sas";
filename nlinmix "&path/nlmm800.sas";
filename nlmem   "&path/nlmem.sas";

%include pigs_dt / nosource;
%include nlinmix / nosource;


 /* NOTE: The variance parameters on the original scale
    differ by 5 orders of magnitude, so we rescale x to make 
    the problem more stable.  The final estimates on the 
    original scale are close to those given by Lindstrom and Bates, 
    but are not exactly equal. */

 /* The derivatives are commented out so that NLINMIX will
    compute them numerically. */


 /*---log-transformed Michaelis-Menten model---*/

%nlinmix(data=pigs,
   model=%str(
      v =_b[1];
      k =_b[2];
      d =_b[3];
      kz=_u[1];
      dz=_u[2];
      a = k + kz + x;
      ef = v*x/a + (d+dz)*x/100;
      predv = log(ef);
   ),
 /*
   derivs=%str(
      d_v = x/a/ef;
      d_k = -v*x/a/a/ef;
      d_d = x/100/ef;
      d_kz = d_k;
      d_dz = d_d;
   ),
 */
   parms=%str(v=.3 k=5 d=1),
   stmts=%str(
      class pig;
      model pseudo_y = d_v d_k d_d / noint notest solution cl;
      random d_kz d_dz / type=un subject=pig solution;
   ), 
   expand  = eblup,
   nlmem   = nlmem
);
