# Vol-BnC-FCMC
First the paths for libraries need to be changed according to the local directories in makefiles:
[coinUtils](src/coinUtils/makefile)  
[osi](src/osi/makefile)  
[bcp](src/bcp/makefile)  
[mcnd](mcnd/makefile)  

To install:
```
mkdir lib
cd src/coinUtils/
mkdir build
make
cd ../osi/
mkdir build
make
cd ../bcp/
mkdir build
make
cd ../../mcnd/
mkdir build
```

To compile:
``make``
