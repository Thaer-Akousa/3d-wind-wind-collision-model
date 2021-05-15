#!/bin/bash

v_orbs=(470)
phs=(0.11 0.89)
for i in {0..1}
do
	ph="${phs[$i]}"
	v_orb="${v_orbs[0]}"
	
	echo "phase: $ph v_orb: $v_orb"
	./bs $ph $v_orb
done

