TARGETDIR="target/generated-sources/swig/loci/slim"
SOURCEDIR="../../../../../src/main/c"
JNIDIRMAC="/Library/Java/JavaVirtualMachines/jdk"*"/Contents/Home/include"
JNIDIRLINUX="/usr/lib/jvm/jdk*/include"
HOME=`pwd`
platform="unknown"
unamestr=`uname`
haspath=""
SWIGPATH="$(which swig)"

#determine OS
if [[ "$unamestr" == 'Linux' ]]; then
   platform='linux'
elif [[ "$unamestr" == 'Darwin' ]]; then
   platform='mac'
fi

#If swig is not on the Eclipse class path, add it.
SWIGPATH="${SWIGPATH%/swig}"
echo $SWIGPATH>~/swigpath.txt
haspath="$(echo "$PATH"|grep -q $SWIGPATH && echo "good")"
echo $haspath>~/haspath.txt
if [[ "$haspath" != 'good' ]]; then
	export PATH=/usr/local/bin:$PATH
fi

#Generate the swig files
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

