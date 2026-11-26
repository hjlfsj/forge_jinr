#!/bin/bash
# crush

./build/bin/crush/crush_trigger  $1
./build/bin/crush/crush_beam $1
./build/bin/crush/crush_t0d1  $1
./build/bin/crush/crush_t0d2 $1
./build/bin/crush/crush_gagg $1

 
