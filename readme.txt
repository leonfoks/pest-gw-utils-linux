
to install pest groundwater utilites on a linux machine just type "bash compile_unix_gw.sh"

I hope it will work for you as it works for me

It has been tested on a ubuntu server 10.4 64bit machine, using ifort (v 12.0.2).


IMPORTANT NOTE:
please, please, PLEASE: read the header of compile_unix_gw.sh before using it, it contains all what you need to know.

knonw issues:
- this does not work with gfortran. I would really appreciate a pull request if anyone would fix this ;)


A. Borghi, University of Neuchatel, 2011


Leon Foks, Apogee Engineering contracted to USGS, added CMAKE capabilities 10/6/2020

Create a new folder "build" in the repo root folder.  From there run "cmake ../src"

