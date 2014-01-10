function compileSlimCurve(Cpath)
% compileSlimCurve(Cpath) a function to compile the mxSlimCurve to use with
%   your installation of Matlab.
%
%       Cpath   Directory containing the SlimCurve source code.
%               [Default = '../main/c/']
%
%   The copileSlimCurve function EcfSingle.c, EcfUtil.c, Ecf.h, and 
%   EcfInternal.h from Cpath into its folder and deletes old compiled
%   binaries. Finally, mxSlimCurve is compiled to run in Matlab.
%
%   compileSlimCurve assumes that a functioning mex compiler in your Matlab
%   installation. If the compiler is installed, functional and compatible
%   with mxSlimCurve code, this function should run without any problem. If
%   if fails to compile, make sure a mex compiler is installed and run the
%   following command on the Matlab command line to select the compiler:
%   
%       mex -setup
%
%
%   The following compilers have been tested with mxSlimCurve:
%   Linux 64-bit:       gcc 4.7.2
%   Windows XP 32-bit:  Windows SDK 7.1
%                       gcc 4.7.1 (part of MinGW)
%
%   The following compilers are not compatible with mxSlimCurve:
%   Windows XP 32-bit:  Lcc-win32 2.4.1
%
%   Note for Linux users only: The function copies the mexopts.sh file 
%   from matlabroot directory and removes the instances of "-ansi" from it.
%   Otherwise, mxSlimCurve will fail to compile.
%
% GNU GPL license 3.0
% copyright 2013-2014 Jakub Nedbal



% Assume that C-files are in '../c/' if no folder has been provided
if ~exist('Cpath', 'var')
    Cpath = ['..' filesep 'main' filesep 'c' filesep];
end

% Check if source directory exists
if ~exist(Cpath, 'dir')
    error('Source directory %s for C-files does not exist', Cpath);
end

% Check if SlimCurve source files exist
files = {'EcfSingle.c', 'EcfUtil.c', 'Ecf.h', 'EcfInternal.h'};
for file = files
    if ~exist([Cpath file{1}], 'file')
        error('Cannot find %s in %s.', file{1}, Cpath);
    end
end


% Check is mc.SlimCurve.c exists in the current folder.
if ~exist('mxSlimCurve.c', 'file')
    error('Cannot find mxSlimCurve.c.');
end

% On Linux machines:
% Copy mexopts.sh file from matlabroot and delete all instances of "-ansi"
% from it. Otherwise mxSlimCurve would not compile
if any(strcmpi(computer, {'GLNXA64', 'GLNX86'}))
    copyfile([matlabroot filesep 'bin' filesep 'mexopts.sh'], ...
             'mexopts.sh', 'f');
    perl('replaceANSI.pl', 'mexopts.sh');
end

% Copy all necessary slimcurve source files into the current directory
files = {'EcfSingle.c', 'EcfUtil.c', 'Ecf.h', 'EcfInternal.h'};
for file = files
    copyfile([Cpath file{1}], '.', 'f');
end

% Delete old mxSlimCurve compiled binaries
if exist(['mxSlimCurve.' mexext], 'file')
    delete(['mxSlimCurve.' mexext]);
end

% Compile mxSlimCurve
switch computer
    case {'GLNXA64', 'GLNX86'}
        fprintf('Compiling mxSlimCurve...\n');
        mex -f ./mexopts.sh mxSlimCurve.c EcfUtil.c EcfSingle.c
    case {'PCWIN', 'PCWIN64'}
        fprintf('Compiling mxSlimCurve...\n');
        mex mxSlimCurve.c EcfUtil.c EcfSingle.c
    otherwise
        fprintf('Not sure how to compile mxSlimCurve on your computer ');
        fprintf(['architecture: ' computer '. Attempting...\n']);
        mex mxSlimCurve.c EcfUtil.c EcfSingle.c
end

% Delete the temporary slimcurve source files
for file = files
    delete(file{1});
end

fprintf('Finished compiling mxSlimCurve.\n');