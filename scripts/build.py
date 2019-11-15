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
#from distutils.dir_util import copy_tree

# globals
SMR    = 'sortmerna'
ZLIB   = 'zlib'
ROCKS  = 'rocksdb'
RAPID  = 'rapidjson'
DIRENT = 'dirent'
CMAKE  = 'cmake'
CONDA  = 'conda'
ALL    = 'all'

URL_ZLIB   = None
URL_ROCKS  = None
URL_DIRENT = None
URL_RAPID  = None
URL_SMR    = None

IS_WIN = None
IS_WSL = None
IS_LNX = None

OS = None

UHOME = None

CMAKE_GEN = None

LIB_DIR = None

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

def conda_install(OS, cfg, dir=None):
    '''
    @param os WIN | WSL | LNX | OSX
    @param cfg   Config dictionary
    @dir         installation root dir. Default: User Home

    pip install -U pyyaml
    pip install -U Jinja2
    '''
    STAMP = '[conda_install]'
    if not dir: dir = UHOME
    fsh = 'miniconda.sh'
    os.chdir(dir)
    url_conda = cfg[CONDA]['url'][OS]
    try:
        with urllib.request.urlopen(url_conda) as surl:
            with open(fsh, 'wb') as fp:
                fp.write(surl.read()) # loads the file into memory before writing to disk. No good for very big files.
    except urllib.request.HTTPError as ex:
        print(ex.read())

    cmd = ['bash', fsh, '-b']
    #cmd = ['bash', fsh, '-b', '-p', '{}/miniconda'.format(dir)]
    proc_run(cmd, dir)

    print('{} Installed conda in {}'.format(STAMP, os.path.join(dir, 'miniconda3')))
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
    url = cfg['cmake'][OS]['url']
    zipped = url.split('/')[-1] # tar.gz or zip
    cmake_home = cfg['cmake'][OS]['home'] # home directory
    # check already installed
    cmake_bin = '{}/bin/cmake'.format(cmake_home)
    if IS_WIN: cmake_bin = '{}.exe'.format(cmake_bin)
    if os.path.exists(cmake_bin):
        print('{} Cmake is already installed: {}'.format(STAMP, cmake_bin))
        if not force:
            is_installed = True

    if not is_installed:
        os.chdir(Path(cmake_home).parent) # navigate to the parent dir e.g. installation root
        # load file from URL
        try:
            with urllib.request.urlopen(url) as surl:
                with open(zipped, 'wb') as fp:
                    fp.write(surl.read()) # loads all file into memory first before writing to disk. No good for very big files.
        except urllib.request.HTTPError as ex:
            print(ex.read())

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

def proc_run(cmd, cwd):
    '''
    '''
    spr = '==========================================================='
    print('[proc_run] Running command:\n{}\n{} in {}\n{}'.format(spr, ' '.join(cmd), cwd, spr))
    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe
    if not os.path.exists(cwd):
        os.makedirs(cwd)
    # cmake configure and generate
    try:
        proc = subprocess.run(cmd, cwd=cwd)
        #proc = subprocess.run(cmd, cwd=build_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.returncode:
            for info in sys.exc_info(): print(info)
    except OSError as err:
        print(err)
        sys.exit()
    except:
        for info in sys.exc_info(): print(info)
        sys.exit()
    #print(proc.stdout)
    #print(proc.stderr)
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
        zlib_link_type   = cfg[SMR][OS]['link'][ptype][ZLIB] # ZLIB linkage type used for SMR given linkage type
        zlib_lib_release = cfg[ZLIB][OS]['link'][zlib_link_type]['release']
        zlib_lib_debug   = cfg[ZLIB][OS]['link'][zlib_link_type]['debug']
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

    # build adn install
    cmd = [ 'cmake', '--build', '.', '--config', btype.title(), '--target', 'install' ]
    proc_run(cmd, SMR_BUILD)
    # generate installation package
    cmd = [ 'cmake', '--build', '.', '--config', btype.title(), '--target', 'package' ]
    proc_run(cmd, SMR_BUILD)

    # test  CMAKE_INSTALL_PREFIX\bin\sortmerna --version
    cmd = [ os.path.join(SMR_DIST, 'bin', 'sortmerna'), '--version' ]
    proc_run(cmd, SMR_BUILD)

    # CMAKE_INSTALL_PREFIX\bin\sortmerna -h
    cmd = [ os.path.join(SMR_DIST, 'bin', 'sortmerna'), '-h' ]
    proc_run(cmd, SMR_BUILD)
#END smr_build

if __name__ == "__main__":
    '''
    python scripts/build.py --name sortmerna
    python /mnt/c/Users/biocodz/a01_code/sortmerna/scripts/build.py --name sortmerna --winhome /mnt/c/Users/biocodz --btype debug
    '''
    import pdb; pdb.set_trace()

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
    optpar.add_option('--config', dest='config', help='Build configuration file.')
    (opts, args) = optpar.parse_args()

    # define platform
    pf = platform.platform()
    IS_WIN = 'Windows' in pf
    IS_WSL = 'Linux' in pf and 'Microsoft' in pf # Windows Subsystem for Linux (WSL)
    IS_LNX = 'Linux' in pf and not 'Microsoft' in pf

    UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']

    if   IS_WIN: OS = 'WIN'
    elif IS_WSL: OS = 'WSL'
    elif IS_LNX: OS = 'LNX'
    else:
        print('Unable to define the platform: {}'.format(pf))
        sys.exit(1)

    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print('Current dir: {}'.format(cur_dir))

    # check env.yaml. If no env file specified, try the current directory
    envfile = os.path.join(cur_dir, 'env.yaml') if not opts.envfile else opts.envfile
    if not os.path.exists(envfile):
        print('No environment config file found. Please, provide one using \'--env\' option')
        sys.exit(1)
    else:
        # load properties from env.yaml
        print('Using Environment configuration file: {}'.format(envfile))
        with open(envfile, 'r') as envh:
            env = yaml.load(envh, Loader=yaml.FullLoader)

    # check build.jinja.yaml
    cfgfile = os.path.join(cur_dir, 'build.jinja.yaml') if not opts.config else opts.config
    if not os.path.exists(cfgfile):
        print('No build configuration template found. Please, provide one using \'--config\' option')
        sys.exit(1)
    else:
        print('Using Build configuration template: {}'.format(cfgfile))

    # load jinja template
    jjenv = Environment(loader=FileSystemLoader(os.path.dirname(cfgfile)), trim_blocks=True, lstrip_blocks=True)
    template = jjenv.get_template(os.path.basename(cfgfile))

    # render jinja template
    env['UHOME'] = UHOME
    cfg_str = template.render(env) # env[OS]
    cfg = yaml.load(cfg_str, Loader=yaml.FullLoader)

    URL_ZLIB   = cfg[ZLIB]['url']
    URL_ROCKS  = cfg[ROCKS]['url']
    URL_DIRENT = cfg[DIRENT]['url']
    URL_RAPID  = cfg[RAPID]['url']
    URL_SMR    = cfg[SMR]['url']

    UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']

    DIRENT_DIST = cfg[DIRENT]['src']

    CMAKE_GEN = cfg[OS]['cmake_gen']

    SMR_SRC   = cfg[SMR][OS].get('src') if cfg[SMR][OS].get('src') else '{}/sortmerna'.format(UHOME)
    SMR_BUILD = cfg[SMR][OS]['build'] if cfg[SMR][OS]['build'] else '{}/build'.format(SMR_SRC)
    SMR_DIST  = cfg[SMR][OS]['dist'] if cfg[SMR][OS]['dist'] else '{}/dist'.format(SMR_SRC)
    SMR_VER   = cfg[SMR].get('ver') if cfg[SMR].get('ver') else None
    
    LIB_DIR = cfg[OS].get('lib') if cfg[OS].get('lib') else UHOME # '{}/3rdparty'.format(SMR_SRC)

    ZLIB_SRC   = cfg[ZLIB][OS].get('src') if cfg[ZLIB][OS].get('src') else '{}/{}'.format(LIB_DIR, ZLIB)
    ZLIB_BUILD = cfg[ZLIB][OS]['build'] if cfg[ZLIB][OS]['build'] else '{}/build'.format(ZLIB_SRC)
    ZLIB_DIST  = cfg[ZLIB][OS]['dist'] if cfg[ZLIB][OS]['dist'] else '{}/dist'.format(ZLIB_SRC)

    ROCKS_SRC   = cfg[ROCKS][OS].get('src') if cfg[ROCKS][OS].get('src') else '{}/{}'.format(LIB_DIR, ROCKS)
    ROCKS_BUILD = cfg[ROCKS][OS]['build'] if cfg[ROCKS][OS]['build'] else '{}/build'.format(ROCKS_SRC)
    ROCKS_DIST  = cfg[ROCKS][OS]['dist'] if cfg[ROCKS][OS]['dist'] else '{}/dist'.format(ROCKS_SRC)
    ROCKS_VER   = cfg[ROCKS].get('ver') if cfg[ROCKS].get('ver') else None

    # no binaries, so always build Release only
    RAPID_SRC   = cfg[RAPID][OS].get('src') if cfg[RAPID][OS].get('src') else '{}/{}'.format(LIB_DIR, RAPID)
    RAPID_BUILD = cfg[RAPID][OS]['build'] if cfg[RAPID][OS]['build'] else '{}/build'.format(RAPID_SRC)
    RAPID_DIST  = cfg[RAPID][OS]['dist'] if cfg[RAPID][OS]['dist'] else '{}/dist'.format(RAPID_SRC)

    if IS_WIN:
        SMR_DIST  = SMR_DIST + '/{}/{}'.format(opts.pt_smr, opts.btype)

        # zlib puts both Debug and Release at the same location => no btype
        opts.pt_zlib = 't1' if not opts.pt_zlib else opts.pt_zlib
        ZLIB_DIST  = ZLIB_DIST + '/{}'.format(opts.pt_zlib) 

        opts.pt_rocks = 't3' if not opts.pt_rocks else opts.pt_rocks
        ROCKS_DIST  = ROCKS_DIST + '/{}/{}'.format(opts.pt_rocks, opts.btype)

        DIRENT_SRC = '{}/{}'.format(LIB_DIR, DIRENT) if not cfg[DIRENT]['src'] else cfg[RAPID]['src']
        DIRENT_DIST = cfg[DIRENT]['dist'] if cfg[DIRENT]['dist'] else DIRENT_SRC

    # call functions
    if opts.name:
        if opts.name == ALL:
            rapidjson_build()
            zlib_build()
            rocksdb_build(ROCKS_VER)
            smr_build(SMR_VER, btype=opts.btype, ptype=opts.pt_smr, cfg=cfg)
        if opts.name == SMR: 
            if opts.clean:
                clean(SMR_BUILD)
            else:
                smr_build(SMR_VER, btype=opts.btype, ptype=opts.pt_smr, cfg=cfg)
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
            if opts.clone:
                conda_install(os, cfg) 
        else: test()