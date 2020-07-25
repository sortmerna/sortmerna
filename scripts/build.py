'''
FILE: build.py
Created: Jul 31, 2019 Wed
'''
import os
import sys
import subprocess
import platform
from optparse import OptionParser
import urllib.request
import tarfile
import re
import fileinput
import yaml
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import zipfile
import time
from shutil import which as shu_which
#from distutils.dir_util import copy_tree

# globals
SMR     = 'sortmerna'
ZLIB    = 'zlib'
ROCKS   = 'rocksdb'
RAPID   = 'rapidjson'
DIRENT  = 'dirent'
CMAKE   = 'cmake'
CONDA   = 'conda'
ALL     = 'all'
CCQUEUE = 'concurrentqueue'

URL_ZLIB   = None
URL_ROCKS  = None
URL_DIRENT = None
URL_RAPID  = None
URL_SMR    = None
URL_CONCURRENTQUEUE = None

MY_OS = None
ENV = None # WIN | WSL | LNX_AWS | LNX_TRAVIS

# define platform
pf = platform.platform()
IS_WIN = 'Windows' in pf
IS_WSL = 'Linux' in pf and 'Microsoft' in pf # Windows Subsystem for Linux (WSL)
IS_LNX = 'Linux' in pf and not 'Microsoft' in pf
if   IS_WIN: MY_OS = 'WIN'
elif IS_WSL: MY_OS = 'WSL'
elif IS_LNX: MY_OS = 'LNX'
else:
    print('Unable to define the platform: {}'.format(pf))
    sys.exit(1)

UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']
CMAKE_GEN   = None
LIB_DIR     = None
DIRENT_DIST = None

SMR_SRC   = None # source root dir
SMR_BUILD = None # build root dir
SMR_DIST  = None # dist root dir

ZLIB_SRC   = None
ZLIB_BUILD = None
ZLIB_DIST  = None

ROCKS_SRC   = None
ROCKS_BUILD = None
ROCKS_DIST  = None

# no binaries, so always build Release only
RAPID_SRC   = None
RAPID_BUILD = None
RAPID_DIST  = None

# Concurrentqueue
CCQUEUE_SRC = None

def test():
    '''
    '''
    print(os.path.dirname(os.path.realpath(__file__)))
    print(os.getcwd())
#END test

def git_clone(url, pdir, force=False):
    '''
    @param url  url to clone
    @param pdir Parent directory where to clone
    '''
    STAMP = '[git_clone]'

    # check clone already exists
    repo = os.path.basename(url).split('.')[0] # e.g. sortmerna
    gitdir = '{}/{}/.git'.format(pdir, repo)
    if os.path.exists(gitdir):
        print('{} Git clone already exists: {}'.format(STAMP, gitdir))
        cmd = ['git', 'status']
        proc_run(cmd, '{}/{}'.format(pdir, repo))
    
    if force or not os.path.exists(gitdir):
        cmd = [ 'git', 'clone', url ]
        proc_run(cmd, pdir)
#END git_clone  

def conda_install(cfg, dir=None, force=False, clean=False):
    '''
    @param MY_OS WIN | WSL | LNX | OSX
    @param cfg   Config dictionary
    @dir         installation root dir. Default: User Home

    pip install -U pyyaml
    pip install -U Jinja2
    '''
    STAMP = '[conda_install]'
    dir = UHOME if not dir else dir
    fsh = 'miniconda.sh'
    #os.chdir(dir)
    # check already installed
    bin_conda = os.path.join(dir, 'miniconda3', 'bin')
    if os.path.exists(bin_conda):
        cmd = ['python', '--version']
        ret = proc_run(cmd, bin_conda, True)
        if not ret['retcode']:
            print('{} Conda already installed: {} : {}'.format(STAMP, bin_conda, ret['stdout']))
            if force:
                print('{} TODO: Force specified - removing the existing installation: {}'.format(STAMP, bin_conda))
            return

    # download the installer if not already present
    if not os.path.exists(os.path.join(dir, fsh)):
        url_conda = cfg[CONDA]['url'][MY_OS]
        if sys.version_info[1] < 7:
            # 'Mozilla/5.0 (Windows NT 6.1; Win64; x64)'
            req = urllib.request.Request(url_conda, headers={'User-Agent': 'Mozilla/5.0'})
        else:
            # this works with conda but not standard python3
            req = url_conda

        try:
            with urllib.request.urlopen(req) as surl:
                with open(fsh, 'wb') as fp:
                    fp.write(surl.read()) # loads the file into memory before writing to disk. No good for very big files.
        except urllib.request.HTTPError as ex:
            print(ex.read())

    # run the installer
    cmd = ['bash', fsh, '-b']
    #cmd = ['bash', fsh, '-b', '-p', '{}/miniconda'.format(dir)]
    proc_run(cmd, dir)

    # delete the installer
    if clean:
        print('{} Deleting the installer {}'.format(STAMP, os.path.join(dir, fsh)))
        os.remove(os.path.join(dir, fsh))

    print('{} Installed conda in {}'.format(STAMP, os.path.join(dir, 'miniconda3')))

    # install packages required to use sortmerna's build.py
    print('{} Installing PyYaml package'.format(STAMP))
    bin_pip = os.path.join(bin_conda, 'pip')
    cmd = [bin_pip, 'install', 'pyyaml']
    proc_run(cmd, bin_conda)

    print('{} Installing Jinja2 package'.format(STAMP))
    cmd = [bin_pip, 'install', 'jinja2']
    proc_run(cmd, bin_conda)

    print('{} Installing scikit-bio package'.format(STAMP))
    cmd = [bin_pip, 'install', 'scikit-bio']
    proc_run(cmd, bin_conda)

    print('{} Done'.format(STAMP))
#END conda_install

def cmake_install(cfg, dir=None, force=False):
    '''
    @param url CMake download URL
    @param dir installation directory. Default User Home

    Download and extract CMake release archive
    python build.py --name cmake --clone
    '''
    STAMP = '[cmake_install]'
    is_installed = False
    url = cfg['cmake'][MY_OS]['url']
    zipped = url.split('/')[-1] # tar.gz or zip
    cmake_home = cfg['cmake'][MY_OS]['home'] # home directory
    # check already installed
    cmake_bin = '{}/bin/cmake'.format(cmake_home)
    if IS_WIN: cmake_bin = '{}.exe'.format(cmake_bin)
    if os.path.exists(cmake_bin):
        print('{} Cmake is already installed: {}'.format(STAMP, cmake_bin))
        if not force:
            is_installed = True

    if not is_installed:
        os.chdir(Path(cmake_home).parent.as_posix()) # navigate to the parent dir e.g. installation root (pathlib)
        # download the installer if not already present
        if not os.path.exists(zipped):
            # load file from URL
            try:
                with urllib.request.urlopen(url) as surl:
                    with open(zipped, 'wb') as fp:
                        fp.write(surl.read()) # loads all file into memory first before writing to disk. No good for very big files.
            except urllib.request.HTTPError as ex:
                print(ex.read())
                sys.exit(1)

        # extract archive
        ext = os.path.splitext(zipped)[1]
        if '.gz' == ext:
            tar = tarfile.open(zipped)
            tar.extractall()
            tar.close()
        elif '.zip' == ext: 
            with zipfile.ZipFile(zipped, 'r') as fp:
                fp.extractall(os.curdir)

        # delete archive
        os.remove(zipped)

        # Verify installation
        if os.path.exists(cmake_bin):
            print('{} Installed CMake {}'.format(STAMP, cmake_bin))
        else:
            print('{} Failed to install CMake {}'.format(STAMP, cmake_bin))

        # copy binarties to HOME/bin to avoid setting the PATH
        #if IS_LNX:
        #    copy_tree('{}/bin'.format(cmake_home), '{}/bin'.format(UHOME))
    #os.environ('PATH') 
#END cmake_install

def gcc_install():
    '''
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test
    sudo apt update
    sudo apt install gcc-7 g++-7
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7
                                           |         |         |        |_priority
                                           |         |         |_path (target)
                                           |         |_name of the symlink (alias)
                                           |_link
    sudo update-alternatives --install /usr/bin/cpp cpp-bin /usr/bin/cpp-9 60
    '''
    cmd = ['add-apt-repository', 'ppa:ubuntu-toolchain-r/test']
    proc_run(cmd)

    cmd = ['apt', 'update']
#END gcc_install

def make_install():
    '''
    sudo apt install make
    '''
#END make_install

def clean(dir):
    '''
    @param dir  Directory to clean

    remove content of the given directory.
    '''
    cmd = ['rm', '-rf', './*']
    proc_run(cmd, dir)
#END clean

def proc_run(cmd, cwd=None, capture=False):
    '''
    '''
    STAMP = '[proc_run]'
    ret = {'retcode':0, 'stdout':None, 'stderr':None}
    spr = '==========================================================='

    if cwd:
        print('{} Running command:\n{}\n{} in {}\n{}'.format(STAMP, spr, ' '.join(cmd), cwd, spr))
    else:
        print('{} Running command:\n{}\n{}\n{}'.format(STAMP, spr, ' '.join(cmd), spr))

    start = time.time()
    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe
    # check executable
    res = shu_which(cmd[0])
    if res :
        try:
            if cwd:
                if not os.path.exists(cwd):
                    os.makedirs(cwd)
                if sys.version_info[1] > 6:
                    proc = subprocess.run(cmd, cwd=cwd, capture_output=capture)
                else:
                    proc = subprocess.run(cmd, cwd=cwd)
            else:
                if sys.version_info[1] > 6:
                    proc = subprocess.run(cmd, capture_output=capture)
                else:
                    proc = subprocess.run(cmd)

            if proc.returncode:
                for info in sys.exc_info(): print(info)
            ret['retcode'] = proc.returncode
            if capture:
                ret['stdout'] = proc.stdout
                ret['stderr'] = proc.stderr
        except OSError as err:
            print(err)
            ret['retcode'] = 1
            ret['stderr'] = err
        except:
            for info in sys.exc_info(): print(info)
            ret['retcode'] = 1
            ret['stderr'] = sys.exc_info()

        print("{} Run time: {}".format(STAMP, time.time() - start))
    else:
        msg = '{} Executable {} not found'.format(STAMP, cmd[0])
        ret['retcode'] = 1
        ret['stderr'] = msg
        print(msg)
    return ret
#END proc_run

def zlib_build(btype='Release'):
    '''
    @param btype  Build type Relase | Debug | ..
    @param ptype  Linkage type like statuc, dynamic, mixed
    '''
    STAMP = '[zlib_build]'    
    git_clone(URL_ZLIB, LIB_DIR)

    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DCMAKE_INSTALL_PREFIX={}'.format(ZLIB_DIST), 
        ZLIB_SRC
        #'-DINSTALL_BIN_DIR={}'.format(INSTALL_BIN_DIR),
        #'-DINSTALL_INC_DIR={}'.format(INSTALL_INC_DIR),
        #'-DINSTALL_LIB_DIR={}'.format(INSTALL_LIB_DIR),
        #'-DINSTALL_MAN_DIR={}'.format(INSTALL_MAN_DIR),
        #'-DINSTALL_PKGCONFIG_DIR={}'.format(INSTALL_PKGCONFIG_DIR),
    ]
    proc_run(cmd, ZLIB_BUILD)

    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    proc_run(cmd, ZLIB_BUILD)
#END zlib_build

def rapidjson_build(btype='Release'):
    '''
    '''
    BUILD_EXAMPLES = 0
    BUILD_DOC = 0
    
    git_clone(URL_RAPID, LIB_DIR)

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DRAPIDJSON_BUILD_EXAMPLES={}'.format(BUILD_EXAMPLES),
        '-DRAPIDJSON_BUILD_DOC={}'.format(BUILD_DOC),
        '-DCMAKE_INSTALL_PREFIX={}'.format(RAPID_DIST), 
        RAPID_SRC
    ]
    proc_run(cmd, RAPID_BUILD)

    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    proc_run(cmd, RAPID_BUILD)
#END rapidjson_build

def rocksdb_fix_3party(ptype='t3'):
    '''
    Only used on Windows

    set(ZLIB_HOME $ENV{THIRDPARTY_HOME}/ZLIB.Library)
    set(ZLIB_INCLUDE ${ZLIB_HOME}/build/native/inc/inc)
    set(ZLIB_LIB_DEBUG ${ZLIB_HOME}/lib/native/debug/amd64/zlib.lib)
    set(ZLIB_LIB_RELEASE ${ZLIB_HOME}/lib/native/retail/amd64/zlib.lib)
    '''
    STAMP = '[rocksdb_fix_3party]'
    if not IS_WIN:
        print('{} not used on Non-Windows'.format(STAMP))
        return

    print('{} Fixing \'thirdparty.inc\' for linkage type [{}] on Windows'.format(STAMP, ptype))

    lib_rel = cfg[ROCKS]['WIN']['link_types'][ptype]['ZLIB_LIB_RELEASE']
    lib_dbg = cfg[ROCKS]['WIN']['link_types'][ptype]['ZLIB_LIB_DEBUG']

    file3p = os.path.join(ROCKS_SRC, 'thirdparty.inc')

    for line in fileinput.FileInput(file3p, inplace=True):
        if line.startswith('set(ZLIB_HOME'):
            line = re.sub(r'ZLIB_HOME .*\)', r'ZLIB_HOME {})'.format(ZLIB_DIST.replace('\\', '/')), line, flags = re.M)
        if line.startswith('set(ZLIB_INCLUDE'):
            line = re.sub(r'ZLIB_INCLUDE .*\)', r'ZLIB_INCLUDE ${ZLIB_HOME}/include)', line, flags = re.M)
        if line.startswith('set(ZLIB_LIB_DEBUG'):
            line = re.sub(r'ZLIB_LIB_DEBUG .*\)', r'ZLIB_LIB_DEBUG ${{ZLIB_HOME}}/lib/{})'.format(lib_dbg), line, flags = re.M)
        if line.startswith('set(ZLIB_LIB_RELEASE'):
            line = re.sub(r'ZLIB_LIB_RELEASE .*\)', r'ZLIB_LIB_RELEASE ${{ZLIB_HOME}}/lib/{})'.format(lib_rel), line, flags = re.M)
        sys.stdout.write(line)
   
    #mo = re.search(r'ZLIB_HOME .*\)+?', txt)
    #txtn = re.sub(r'ZLIB_HOME .*\)+?', r'ZLIB_HOME {})'.format(ZLIB_DIST.replace('\\', '/')), txt, flags = re.M)
    #txtn = re.sub(r'ZLIB_INCLUDE .*\)+?', r'ZLIB_INCLUDE ${ZLIB_HOME}/include)', txtn, flags = re.M)
    #txtn = re.sub(r'ZLIB_LIB_DEBUG .*\)', r'ZLIB_LIB_DEBUG ${{ZLIB_HOME}}/lib/{})'.format(lib_dbg), txtn, flags = re.M)
    #txtn = re.sub(r'ZLIB_LIB_RELEASE .*\)', r'ZLIB_LIB_RELEASE ${{ZLIB_HOME}}/lib/{})'.format(lib_rel), txtn, flags = re.M)
#END rocksdb_fix_3party
    
def rocksdb_build(ver=None, btype='Release', ptype='t3'):
    '''
    @param btype  Build type Release | Debug
    @param ptype  Linkage type on Windows t1 | t2 | t3

    NOTE: on Windows 'thridparty.inc' file has to be modified.
    '''
    STAMP = '[rocksdb_build]'
    
    OPTDBG = 0
    WITH_RUNTIME_DEBUG = 0
    WITH_ZLIB = 1
    WITH_XPRESS = 0
    WITH_TESTS = 0
    WITH_TOOLS = 0
    PORTABLE = 1
    ROCKSDB_INSTALL_ON_WINDOWS = 0
    WITH_MD_LIBRARY = 0
    WITH_GFLAGS = 0
    
    git_clone(URL_ROCKS, LIB_DIR)

    if ver:
        cmd = ['git', 'checkout', ver]
        proc_run(cmd, ROCKS_SRC)

    if IS_WIN:
        if ptype == "t3": print("Type 3 linkage: /MD + static ZLib")

        WITH_MD_LIBRARY = 1
        ROCKSDB_INSTALL_ON_WINDOWS = 1
        WITH_XPRESS = 1
        PORTABLE = 0

    rocksdb_fix_3party(ptype) # executes on Windows

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DOPTDBG={}'.format(OPTDBG),
        '-DPORTABLE={}'.format(PORTABLE),
        '-DCMAKE_BUILD_TYPE={}'.format(btype),
        '-DWITH_RUNTIME_DEBUG={}'.format(WITH_RUNTIME_DEBUG),
        '-DWITH_GFLAGS={}'.format(WITH_GFLAGS),
        '-DZLIB_ROOT={}'.format(ZLIB_DIST), # v6
        '-DZLIB_ROOT_DIR={}'.format(ZLIB_DIST), # v5.x.x
        '-DWITH_ZLIB={}'.format(WITH_ZLIB),
        '-DWITH_TESTS={}'.format(WITH_TESTS),
        '-DWITH_TOOLS={}'.format(WITH_TOOLS),
        '-DCMAKE_INSTALL_PREFIX={}'.format(ROCKS_DIST), 
        ROCKS_SRC
    ]

    if IS_WIN:
        cmd.append('-DWITH_MD_LIBRARY={}'.format(WITH_MD_LIBRARY))
        cmd.append('-DWITH_XPRESS={}'.format(WITH_XPRESS))
        cmd.append('-DROCKSDB_INSTALL_ON_WINDOWS={}'.format(ROCKSDB_INSTALL_ON_WINDOWS))

    proc_run(cmd, ROCKS_BUILD)

    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    proc_run(cmd, ROCKS_BUILD)
#END rocksdb_build

def smr_build(ver=None, btype='Release', ptype='t1', cfg={}):
    '''
    @param btype Build type Release | Debug
    @param ptype Linking type: t1 | t2 | t3
        t1 all static
        t2 static 3rd party + dynamic runtime
        t3 all dynamic
    '''
    STAMP = '[smr_build]'

    git_clone(URL_SMR, os.path.split(SMR_SRC)[0])

    if ver:
        cmd = ['git', 'checkout', ver]
        proc_run(cmd, SMR_SRC)

    # CMake flags
    PORTABLE = 0
    WITH_MD_LIBRARY = 0 # only Windows
    WITH_RUNTIME_DEBUG = 0
    WITH_TESTS = 1
    CPACK_BINARY_NSIS = 0 # Win
    CPACK_BINARY_7Z = 1 # Win
    CPACK_BINARY_ZIP = 1 # Win
    CPACK_SOURCE_7Z = 1 # Win
    CPACK_SOURCE_ZIP = 1 # Win
    CPACK_BINARY_TGZ = 1 # Lin
    ROCKSDB_STATIC = 1
    ZLIB_STATIC = 1
    # -DEXTRA_CXX_FLAGS_RELEASE="-lrt" (had to use this on Centos 6.6 + GCC 7.3.0)

    if IS_LNX:
        ZLIB_LIBRARY_RELEASE = '{}/lib/libz.a'.format(ZLIB_DIST)
        ZLIB_LIBRARY_DEBUG = ZLIB_LIBRARY_RELEASE
        global SMR_BUILD
        SMR_BUILD = os.path.join(SMR_BUILD, btype.title()) # sortmerna/build/Release/ sortmerna/build/Debug/
        if 't1' == ptype:
            PORTABLE = 1 # sets '-static' flag
    elif IS_WIN:
        # ZLIB
        zlib_link_type   = cfg[SMR]['link'][MY_OS][ptype][ZLIB] # ZLIB linkage type used for SMR given linkage type
        zlib_lib_release = cfg[ZLIB]['link'][MY_OS][zlib_link_type]['release']
        zlib_lib_debug   = cfg[ZLIB]['link'][MY_OS][zlib_link_type]['debug']
        ZLIB_LIBRARY_RELEASE = '{}/lib/{}'.format(ZLIB_DIST, zlib_lib_release)
        ZLIB_LIBRARY_DEBUG = '{}/lib/{}'.format(ZLIB_DIST, zlib_lib_debug)
        #ROCKSDB_SRC = os.path.join(LIBDIR, 'rocksdb')
        WITH_MD_LIBRARY = 1

    #:: print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DPORTABLE={}'.format(PORTABLE),
        '-DWITH_RUNTIME_DEBUG={}'.format(WITH_RUNTIME_DEBUG),
        '-DWITH_TESTS={}'.format(WITH_TESTS),
        '-DZLIB_STATIC={}'.format(ZLIB_STATIC),
        '-DROCKSDB_STATIC={}'.format(ROCKSDB_STATIC),
        '-DZLIB_ROOT={}'.format(ZLIB_DIST),
        '-DZLIB_LIBRARY_RELEASE={}'.format(ZLIB_LIBRARY_RELEASE),
        '-DZLIB_LIBRARY_DEBUG={}'.format(ZLIB_LIBRARY_DEBUG),
        #'-DROCKSDB_SRC={}'.format(ROCKSDB_SRC),
        '-DROCKSDB_HOME={}'.format(ROCKS_DIST),
        '-DRAPIDJSON_HOME={}'.format(RAPID_DIST),
        '-DCONCURRENTQUEUE_HOME={}'.format(CCQUEUE_SRC),
        '-DCMAKE_INSTALL_PREFIX={}'.format(SMR_DIST)
    ]

    if IS_LNX:
        cmd.append('-DCMAKE_BUILD_TYPE={}'.format(btype))
        cmd.append('-DCPACK_BINARY_TGZ={}'.format(CPACK_BINARY_TGZ))
    elif IS_WIN:
        cmd.append('-DDIRENTWIN_HOME={}'.format(DIRENT_DIST))
        cmd.append('-DWITH_MD_LIBRARY={}'.format(WITH_MD_LIBRARY))
        cmd.append('-DCPACK_BINARY_NSIS={}'.format(CPACK_BINARY_NSIS))
        cmd.append('-DCPACK_BINARY_7Z={}'.format(CPACK_BINARY_7Z))
        cmd.append('-DCPACK_BINARY_ZIP={}'.format(CPACK_BINARY_ZIP))
        cmd.append('-DCPACK_SOURCE_7Z={}'.format(CPACK_SOURCE_7Z))
        cmd.append('-DCPACK_SOURCE_ZIP={}'.format(CPACK_SOURCE_ZIP))

    if opts.vb:
        cmd.append('-DCMAKE_EXPORT_COMPILE_COMMANDS=1')
    if opts.loglevel:
        cmd.append('--loglevel={}'.format(opts.loglevel.upper()))
    elif opts.trace:
        cmd.append('--trace')
    cmd.append(SMR_SRC)

    # Run the build command
    proc_run(cmd, SMR_BUILD)

    # build and install
    cmd = [ 'cmake', '--build', '.', '--config', btype.title(), '--target', 'install' ]
    proc_run(cmd, SMR_BUILD)
    # generate installation package
    cmd = [ 'cmake', '--build', '.', '--config', btype.title(), '--target', 'package' ]
    proc_run(cmd, SMR_BUILD)

    # test  CMAKE_INSTALL_PREFIX\bin\sortmerna --version
    SMR_EXE = 'sortmerna.exe' if IS_WIN else 'sortmerna'
    cmd = [ os.path.join(SMR_DIST, 'bin', SMR_EXE), '--version' ]
    proc_run(cmd, SMR_BUILD)

    # CMAKE_INSTALL_PREFIX\bin\sortmerna -h
    cmd = [ os.path.join(SMR_DIST, 'bin', SMR_EXE), '-h' ]
    proc_run(cmd, SMR_BUILD)
#END smr_build

def concurrentqueue_build(cfg={}):
    '''
    a single header file - just clone and use
    '''
    git_clone(URL_CONCURRENTQUEUE, LIB_DIR)
#END concurrentqueue_build

if __name__ == "__main__":
    '''
    python scripts/build.py --name sortmerna [--envn WIN] [--env scripts/env_non_git.yaml]
    python scripts/build.py --name sortmerna
    python scripts/build.py --name cmake [--env scripts/env_non_git.yaml]
        --winhome /mnt/c/Users/biocodz --btype debug
    '''
    STAMP = '[build.py:__main__]'
    #import pdb; pdb.set_trace()

    is_opts_ok = True

    # options
    optpar = OptionParser()
    optpar.add_option('-n', '--name', dest='name', help='Module to build e.g. sortmerna | zlib | rocksdb | all')
    optpar.add_option('--clone', action="store_true", help='Perform git clone for the given name')
    optpar.add_option('-c', '--clean', action="store_true", help='clean build directory for the given name')
    optpar.add_option('--btype', dest='btype', default='release', help = 'Build type: release | debug')
    optpar.add_option('--pt_smr', dest='pt_smr', default='t1', help = 'Sortmerna Linkage type t1 | t2 | t3')
    optpar.add_option('--pt_zlib', dest='pt_zlib', help = 'Zlib Linkage type t1 | t2 | t3')
    optpar.add_option('--pt_rocks', dest='pt_rocks', help = 'Rocksdb Linkage type t1 | t2 | t3')
    optpar.add_option('--rocks3p', dest='rocks3p', help='Fix thirdparty.inc when building Rocksb')
    optpar.add_option('--winhome', dest='winhome', help='When building on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    optpar.add_option('--trace', action="store_true", help='Run cmake with --trace')
    optpar.add_option('--loglevel', dest='loglevel', help = 'Cmake log level')
    optpar.add_option('--vb', action="store_true", help='Export compile commands')
    optpar.add_option('--env', dest='envfile', help='Env configuration file.')
    optpar.add_option('--envn', dest='envname', help=('Name of environment: WIN | WSL '
                                                  '| LNX_AWS | LNX_TRAVIS | LNX_VBox_Ubuntu_1804 | ..'))
    optpar.add_option('--config', dest='config', help='Build configuration file.')
    (opts, args) = optpar.parse_args()

    UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']

    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print('Current dir: {}'.format(cur_dir))

    # check env.yaml. If no env file specified, try the current directory
    envfile = os.path.join(cur_dir, 'env.yaml') if not opts.envfile else opts.envfile
    if not os.path.exists(envfile):
        print('No environment config file found. Please, provide one using \'--env\' option')
        sys.exit(1)

    # load properties from env.jinja.yaml
    print('Using Environment configuration file: {}'.format(envfile))
    #with open(envfile, 'r') as envh:
    #    env = yaml.load(envh, Loader=yaml.FullLoader)
    env_jj = Environment(loader=FileSystemLoader(os.path.dirname(envfile)), trim_blocks=True, lstrip_blocks=True)
    env_template = env_jj.get_template(os.path.basename(envfile))

    #   render jinja template
    env_str = env_template.render({'UHOME': UHOME})
    env = yaml.load(env_str, Loader=yaml.FullLoader)

    if not opts.envname:
        if IS_WIN or IS_WSL: 
            ENV = MY_OS
            print('{} --envn was not specified - using {}'.format(STAMP, ENV))
        else:
            print('{} --envn is required on OS {}'.format(STAMP, MY_OS))
            is_opts_ok = False
    else:
        ENV = opts.envname

    # check build.jinja.yaml
    #cfgfile = os.path.join(cur_dir, 'build.jinja.yaml') if not opts.config else opts.config
    #if not os.path.exists(cfgfile):
    #    print('No build configuration template found. Please, provide one using \'--config\' option')
    #    sys.exit(1)
    #else:
    #    print('Using Build configuration template: {}'.format(cfgfile))

    # load jinja template
    #jjenv = Environment(loader=FileSystemLoader(os.path.dirname(cfgfile)), trim_blocks=True, lstrip_blocks=True)
    #template = jjenv.get_template(os.path.basename(cfgfile))

    # render jinja template
    #env['UHOME'] = UHOME
    #cfg_str = template.render(env) # env[OS]
    #cfg = yaml.load(cfg_str, Loader=yaml.FullLoader)
    #
    libdir = env.get('LIB_DIR', {}).get(ENV)
    LIB_DIR =  libdir if libdir else UHOME

    URL_ZLIB   = env[ZLIB]['url']
    URL_ROCKS  = env[ROCKS]['url']
    URL_DIRENT = env[DIRENT]['url']
    URL_RAPID  = env[RAPID]['url']
    URL_SMR    = env[SMR]['url']
    URL_CONCURRENTQUEUE = env[CCQUEUE]['url']

    CMAKE_GEN = env[CMAKE]['generator'][MY_OS]

    # SMR
    val = env.get(SMR,{}).get('src',{}).get(ENV)
    SMR_SRC   = val if val else '{}/sortmerna'.format(UHOME)
    val = env.get(SMR,{}).get('build',{}).get(ENV)
    SMR_BUILD = val if val else '{}/build'.format(SMR_SRC)
    val = env.get(SMR,{}).get('dist',{}).get(ENV)
    SMR_DIST  = val if val else '{}/dist'.format(SMR_SRC)
    SMR_VER   = env.get(SMR,{}).get('ver')

    # ZLIB
    val = env.get(ZLIB,{}).get('src',{}).get(ENV)
    ZLIB_SRC   = val if val else '{}/{}'.format(LIB_DIR, ZLIB)
    val = env.get(ZLIB,{}).get('build',{}).get(ENV)
    ZLIB_BUILD = val if val else '{}/build'.format(ZLIB_SRC)
    val = env.get(ZLIB,{}).get('dist',{}).get(ENV)
    ZLIB_DIST  = val if val else '{}/dist'.format(ZLIB_SRC)

    # ROCKSDB
    val = env.get(ROCKS,{}).get('src',{}).get(ENV)
    ROCKS_SRC   = val if val else '{}/{}'.format(LIB_DIR, ROCKS)
    val = env.get(ROCKS,{}).get('build',{}).get(ENV)
    ROCKS_BUILD = val if val else '{}/build'.format(ROCKS_SRC)
    val = env.get(ROCKS,{}).get('dist',{}).get(ENV)
    ROCKS_DIST  = val if val else '{}/dist'.format(ROCKS_SRC)
    ROCKS_VER   = env.get(ROCKS, {}).get('ver')

    # RAPIDJSON
    # no binaries, so always build Release only
    val = env.get(RAPID,{}).get('src',{}).get(ENV)
    RAPID_SRC   = val if val else '{}/{}'.format(LIB_DIR, RAPID)
    val = env.get(RAPID,{}).get('build',{}).get(ENV)
    RAPID_BUILD = val if val else '{}/build'.format(RAPID_SRC)
    val = env.get(RAPID,{}).get('dist',{}).get(ENV)
    RAPID_DIST  = val if val else '{}/dist'.format(RAPID_SRC)

    # CONCURRENTQUEUE
    val = env.get(CCQUEUE,{}).get('src',{}).get(ENV)
    CCQUEUE_SRC = val if val else '{}/{}'.format(LIB_DIR, CCQUEUE)

    if 'WIN' == ENV:
        SMR_DIST  = SMR_DIST + '/{}/{}'.format(opts.pt_smr, opts.btype)

        # zlib puts both Debug and Release at the same location => no btype
        opts.pt_zlib = 't1' if not opts.pt_zlib else opts.pt_zlib
        ZLIB_DIST  = ZLIB_DIST + '/{}'.format(opts.pt_zlib) 

        opts.pt_rocks = 't3' if not opts.pt_rocks else opts.pt_rocks
        ROCKS_DIST  = ROCKS_DIST + '/{}/{}'.format(opts.pt_rocks, opts.btype)

        val = env.get(DIRENT, {}).get('src', {}).get(ENV)
        DIRENT_SRC = val if val else '{}/{}'.format(LIB_DIR, DIRENT)
        val = env.get(DIRENT, {}).get('dist')
        DIRENT_DIST = val.get(ENV) if val and isinstance(val, dict) else DIRENT_SRC

    # call functions
    if is_opts_ok and opts.name:
        if opts.name == ALL:
            concurrentqueue_build(cfg) 
            rapidjson_build()
            zlib_build()
            rocksdb_build(ROCKS_VER)
            smr_build(SMR_VER, btype=opts.btype, ptype=opts.pt_smr, cfg=env)
        if opts.name == SMR: 
            if opts.clean:
                clean(SMR_BUILD)
            else:
                smr_build(SMR_VER, btype=opts.btype, ptype=opts.pt_smr, cfg=env)
        elif opts.name == ZLIB: 
            zlib_build()
        elif opts.name == RAPID: 
            rapidjson_build()
        elif opts.name == ROCKS: 
            if opts.clean:
                clean(ROCKS_BUILD)
            else:
                rocksdb_build(ROCKS_VER)
        elif opts.name == DIRENT: 
            if opts.clone:
                git_clone(URL_DIRENT, LIB_DIR) 
        elif opts.name == CMAKE: 
            if opts.clone:
                git_clone(cfg[CMAKE]['url'], LIB_DIR)
            cmake_install(cfg) 
        elif opts.name == CONDA: 
            conda_install(cfg) 
        elif opts.name == CCQUEUE: 
            concurrentqueue_build(cfg) 
        else: test()