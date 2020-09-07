#!/bin/sh

echo 'tree.7'
grep -i 'error' tree.7.log
grep -i 'warning' tree.7.log

echo 'tree.an.sas'
grep -i 'error' tree.an.log
grep -i 'warning' tree.an.log

echo 'tree'
grep -i 'error' tree.log
grep -i 'warning' tree.log

echo 'tree0'
grep -i 'error' tree0.log
grep -i 'warning' tree0.log

echo 'tree1'
grep -i 'error' tree1.log
grep -i 'warning' tree1.log

echo 'tree5g'
grep -i 'error' tree5g.log
grep -i 'warning' tree5g.log

echo 'tree_ode'
grep -i 'error' tree_ode.log
grep -i 'warning' tree_ode.log

echo 'tree_ode2'
grep -i 'error' tree_ode2.log
grep -i 'warning' tree_ode2.log

echo 'tree_ode3'
grep -i 'error' tree_ode3.log
grep -i 'warning' tree_ode3.log


