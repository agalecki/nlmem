#!/bin/sh
echo 'pheno'
grep -i 'error' pheno.log
grep -i 'warning' pheno.log

echo 'pigs'
grep -i 'error' pigs.log
grep -i 'warning' pigs.log

echo 'tree'
grep -i 'error' tree.log
grep -i 'warning' tree.log

