./nlmem_tree

1. tree.sas

        Model specification: Analytical 
        Random effects:      1
        Derivatives: numerical
        CPU: 6.58 secs
        Reads 1 record at a time.

2. tree.7.sas 

        Model specification: Analytical
        Random effects:      1
        Derivatives: numerical
        Reads 7 records at a time.
        CPU: 4:65 secs

3. tree1.sas 

        Model specification: Analytical
        Random effects:      1
        Derivatives: numerical specified using NLPFDD
        Reads 7 records at a time.
        CPU: 4:43 secs
        Note: SAS IML code in the user defined macro.
              Global variable dbu used.

4. tree_ode.sas

        Model specification: ODE for within-subject variation
                             Analytic for between-subject variation
        Derivatives:         Numerical for within-subject model using NLPFDD
                             Analytic for between-subject varaiation

        CPU: 25:62 secs
        Note: Convergence criterion 2e-5

5. tree.an.sas

        Model specification: Analytical 
        Random effects:      1
        Derivatives:         Analytical
        Reads 1 record at a time.
        CPU:                 5.99 secs


6. tree0.sas

        Model specification: Analytical 
        Random effects:      None
        Derivatives:         Analytical
        Reads                7 records at a time.
        CPU:                 2.39 secs
        Note:                savevar illustrated
                             savevar and next declared in the body 
                             of the useriml macro.


7. tree5g.sas      

       dtmplt = 1 1 1 0 0
       expand = 1 0
       CPU:                 4.58 sec

8. tree_ode2.sas

        Model specification: ODE for within-subject variation
                             Analytic for between-subject variation
        Derivatives:         Numerical for within-subject model using NLPFDD
                             Analytic for between-subject varaiation
        CPU: 38:59 secs
        Note: Convergence criterion 2e-5
              User specified NLIN step with constraints for parameters

9. tree_ode3.sas

        Model specification: ODE 
        Derivatives:         Numerical
        CPU:                 39:75 secs
        Note: Convergence criterion 2e-5
              User specified NLIN step with constraints for parameters





