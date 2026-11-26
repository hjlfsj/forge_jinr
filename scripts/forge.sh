#!/bin/bash
# crush

./build/bin/crush/crush_trigger  $1
./build/bin/crush/crush_beam $1
./build/bin/crush/crush_t0d1  $1
./build/bin/crush/crush_t0d2 $1
./build/bin/crush/crush_gagg $1

./build/bin/crush_vme 16
./build/bin/sift -x 99 -v 16 

./build/bin/coke -x 99 -v 16
 
./build/bin/smelt -r $1 -t main t0d1 t0d2 gagg beam 
./build/bin/smelt -r $1 -t t1 t0d1 t0d2 gagg beam 

./build/bin/smelt -x 99 -v 16  -t t1 tafd tafcsi t1du t1csiu
