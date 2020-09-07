#!/bin/sh
echo 'pheno'
grep -i 'error' pheno.log
grep -i 'warning' pheno.log

echo 'pheno.1'
grep -i 'error' pheno.1.log
grep -i 'warning' pheno.1.log

echo 'pheno.1ode'
grep -i 'error' pheno.1ode.log
grep -i 'warning' pheno.1ode.log

