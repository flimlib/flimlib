
#object files generated in curr directory (so in slim-curve/). To get around this need to cd to desired dir and then cd back here?
TARGETDIR="target/generated-sources/swig/loci/slim"
SOURCEDIR="../../../../../src/main/c"
JNIDIRMAC="/Library/Java/JavaVirtualMachines/jdk"*"/Contents/Home/include"
JNIDIRLINUX="/usr/lib/jvm/jdk*/include"
HOME=`pwd`
platform="unknown"
unamestr=`uname`
#determine OS
if [[ "$unamestr" == 'Linux' ]]; then
   platform='linux'
elif [[ "$unamestr" == 'Darwin' ]]; then
   platform='mac'
fi
if [[ $platform == 'linux' ]]; then
	cd $JNIDIRLINUX #TODO: Hack to get the expanded directory
	JNIDIRLINUX=`pwd`
	cd $HOME
	#Swig commands to generate swig files for linux
	swig -java -package loci.slim -outdir $TARGETDIR src/main/c/cLibrary.i 
	cd $TARGETDIR
	cc -c -fpic $SOURCEDIR/*.c -I"$JNIDIRLINUX" -I"$JNIDIRLINUX/linux"
	cc -shared *.o -o libEcfGlobalWrapper.so
	rm *.o
elif [[ $platform == 'mac' ]]; then
	cd $JNIDIRMAC #TODO: Hack to get the expanded directory. 
	JNIDIRMAC=`pwd`
	cd $HOME
	#Swig commands to generate swig files for linux
	swig -java -package loci.slim -outdir $TARGETDIR src/main/c/cLibrary.i 
	cd $TARGETDIR
	cc -c $SOURCEDIR/*.c -I"$JNIDIRMAC" -I"$JNIDIRMAC/darwin"
	cc -framework JavaVM -bundle *.o -o libEcfGlobalWrapper.jnilib
	rm *.o
else
	echo "Unsupported OS: $unamestr"
fi

