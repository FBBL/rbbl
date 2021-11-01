# RBBL - RAM-Based BKW for LWE 
[![Master Branch](https://img.shields.io/badge/-master:-gray.svg)](https://github.com/FBBL/fbbl/tree/master) [![Build Status](https://app.travis-ci.com/FBBL/rbbl.svg?branch=master)](https://app.travis-ci.com/FBBL/rbbl) [![Coverage Status](https://coveralls.io/repos/github/FBBL/rbbl/badge.svg?branch=master)](https://coveralls.io/github/FBBL/rbbl?branch=master)

This library implements the BKW algorithm for solving LWE instances presented [here](https://www.mdpi.com/2410-387X/5/4/31). Differently from [fbbl](https://github.com/FBBL/fbbl), rbbl uses only RAM memory. It is therefore much faster than fbbl, but it is also limited by the RAM memory available. It supports smooth-Lazy Modulus Switch and the FWHT-based guessing method (with ot without bruteforce).

As an example, the hardest LWE instance we solved with this library had parameters `n=40`, `q=1601` and standard deviation `σ = 0.005q`. Using 15 cores and 32Gb of RAM, it took <20 seconds. 


##### Build
```
mkdir target
cd target
cmake ..
make
make test
```

Author: Alessandro Budroni
