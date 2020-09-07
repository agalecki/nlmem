./pheno

pheno.1.sas 

     Data for one subject are read in at once (next=indiv).
     Time-independent variables (weight, apgarlow) are declared (xvar=).
     Elements of vector predv are defined recursively in the loop.
     Note: _weight_ = 1/(predv#predv)
     CPU: 25 secs

pheno.1ode.sas

     ODE specification.
     CPU  2 min 6 sec
     

pheno.sas 

     User defined call to NLIN with modified termination criteria
            (Max # iterations=2)
     CPU reduced from  2min 56 sec to 53.44 secs .
     See pheno.sas in ../nlmem_MAINexamples

     Note: savevar illustrated
