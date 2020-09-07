#!/bin/sh
echo 'pheno1'
grep -i 'error:' pheno1.log
grep -i 'warning' pheno1.log

echo 'tree1'
grep -i 'error:' tree1.log
grep -i 'warning' tree1.log

echo 'tree1a'
grep -i 'error:' tree1a.log
grep -i 'warning' tree1a.log



