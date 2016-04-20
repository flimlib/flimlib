swig -java -package loci.slim -outdir target/generated-sources/swig/loci/slim src/main/c/cLibrary.i 
#object files generated in curr directory (so in slim-curve/). To get around this need to cd to desired dir and then cd back here?
cd target/generated-sources/swig/loci/slim
cc -c ../../../../../src/main/c/EcfGlobal.c ../../../../../src/main/c/cLibrary_wrap.c ../../../../../src/main/c/EcfUtil.c ../../../../../src/main/c/EcfSingle.c ../../../../../src/main/c/EcfSPA.c ../../../../../src/main/c/GCI_Phasor.c -I"/Library/Java/JavaVirtualMachines/jdk1.8.0_40.jdk/Contents/Home/include" -I"/Library/Java/JavaVirtualMachines/jdk1.8.0_40.jdk/Contents/Home/include/darwin"
cc -framework JavaVM -bundle EcfGlobal.o cLibrary_wrap.o EcfUtil.o EcfSingle.o EcfSPA.o GCI_Phasor.o -o libEcfGlobalWrapper.jnilib
