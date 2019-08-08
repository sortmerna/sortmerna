'''
FILE: build.py
Created: Jul 31, 2019 Wed
'''
import os
import sys
import subprocess
import platform
from optparse import OptionParser

def test():
    '''
    '''
    print(os.path.dirname(os.path.realpath(__file__)))
    print(os.getcwd())
#END test

def git_clone(url, pdir):
    '''
    @param url  url to clone
    @param pdir Parent directory where to clone
    '''
    cmd = [ 'git', 'clone', url ]

    try:
        subprocess.run(cmd, cwd=pdir)
    except OSError as err:
        print(err)
        raise
    except:
        for info in sys.exc_info(): print(info)
        raise
#END git_clone

def zlib_clone():
    '''
    '''
    if 'Windows' in platform.platform():
        PDIR = os.path.join(os.environ['USERPROFILE'], 'a03_libs')
    else:
        PDIR = os.path.join(os.environ['HOME'])

    url = 'https://github.com/madler/zlib.git'
    git_clone(url, PDIR)

def dirent_clone():
    '''
    dirent win
    '''
    if 'Windows' in platform.platform():
        PDIR = os.path.join(os.environ['USERPROFILE'], 'a03_libs')
    else:
        print('Used only for Windows')
        return
    url = 'https://github.com/tronkko/dirent'
    git_clone(url, PDIR)

def rocksdb_clone():
    '''
    '''
    if 'Windows' in platform.platform():
        PDIR = os.path.join(os.environ['USERPROFILE'], 'a03_libs')
    else:
        PDIR = os.path.join(os.environ['HOME'])
    url = 'https://github.com/facebook/rocksdb.git'
    git_clone(url, PDIR)

def smr_clone():
    '''
    '''
    if 'Windows' in platform.platform():
        PDIR = os.path.join(os.environ['USERPROFILE'], 'a01_code')
    else:
        PDIR = os.path.join(os.environ['HOME'])

    url = 'https://github.com/biocore/sortmerna.git'
    git_clone(url, PDIR)

def clean(dir):
    '''
    @param dir  Directory to clean

    remove content of the given directory.
    '''
    cmd = ['rm', '-rf', './*']
    cmake_run(cmd, dir)
#END clean

def rapidjson_clone():
    '''
    '''
    if 'Windows' in platform.platform():
        PDIR = os.path.join(os.environ['USERPROFILE'], 'a01_code')
    else:
        PDIR = os.path.join(os.environ['HOME'])
    url = 'https://github.com/Tencent/rapidjson'
    git_clone(url, PDIR)

def cmake_run(cmd, cwd):
    '''
    '''
    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe
    if not os.path.exists(cwd):
        os.makedirs(cwd)
    # cmake configure and generate
    try:
        subprocess.run(cmd, cwd=cwd)
        #proc = subprocess.run(cmd, cwd=build_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError as err:
        print(err)
        raise
    except:
        for info in sys.exc_info(): print(info)
        raise
    #print(proc.stdout)
    #print(proc.stderr)
#END cmake_run

def smr_build(gen='Unix Makefiles', btype='Release', ptype='t3', 
        src=None, build=None, dist=None, zlib=None, rocks=None, rapid=None, dirent=None):
    '''
    @param btype Build type Release | Debug
    @param ptype Linking type
    '''
    SRC_DIR = src
    BUILD_DIR = build
    DIST_DIR = dist
    CMAKE_GEN = gen
    ZLIB_DIST = zlib
    ROCKSDB_DIST = rocks
    RAPIDJSON_DIST = rapid
    DIRENTWIN_DIST = dirent

    PORTABLE = 0
    WITH_MD_LIBRARY = 1
    WITH_RUNTIME_DEBUG = 0
    WITH_TESTS = 1
    CPACK_BINARY_NSIS = 0
    CPACK_BINARY_7Z = 1
    CPACK_BINARY_ZIP = 1
    CPACK_SOURCE_7Z = 1
    CPACK_SOURCE_ZIP = 1
    ROCKSDB_STATIC = 1
    ZLIB_STATIC = 1

    pf = platform.platform()

    if 'Linux' in pf:
        ZLIB_LIBRARY_RELEASE = os.path.join(ZLIB_DIST, 'lib', 'libz.a')
        ZLIB_LIBRARY_DEBUG = ZLIB_LIBRARY_RELEASE
    elif 'Windows' in pf:
        ZLIB_LIBRARY_RELEASE = os.path.join(ZLIB_DIST, 'lib', 'zlibstatic.lib')
        ZLIB_LIBRARY_DEBUG = os.path.join(ZLIB_DIST, 'lib', 'zlibstaticd.lib')
        #ROCKSDB_SRC = os.path.join(LIBDIR, 'rocksdb')

    #:: print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DPORTABLE={}'.format(PORTABLE),
        '-DCPACK_BINARY_NSIS={}'.format(CPACK_BINARY_NSIS),
        '-DCPACK_BINARY_7Z={}'.format(CPACK_BINARY_7Z),
        '-DCPACK_BINARY_ZIP={}'.format(CPACK_BINARY_ZIP),
        '-DCPACK_SOURCE_7Z={}'.format(CPACK_SOURCE_7Z),
        '-DCPACK_SOURCE_ZIP={}'.format(CPACK_SOURCE_ZIP),
        '-DWITH_MD_LIBRARY={}'.format(WITH_MD_LIBRARY),
        '-DWITH_RUNTIME_DEBUG={}'.format(WITH_RUNTIME_DEBUG),
        '-DWITH_TESTS={}'.format(WITH_TESTS),
        '-DZLIB_STATIC={}'.format(ZLIB_STATIC),
        '-DROCKSDB_STATIC={}'.format(ROCKSDB_STATIC),
        '-DZLIB_ROOT={}'.format(ZLIB_DIST),
        '-DZLIB_LIBRARY_RELEASE={}'.format(ZLIB_LIBRARY_RELEASE),
        '-DZLIB_LIBRARY_DEBUG={}'.format(ZLIB_LIBRARY_DEBUG),
        #'-DROCKSDB_SRC={}'.format(ROCKSDB_SRC),
        '-DROCKSDB_HOME={}'.format(ROCKSDB_DIST),
        '-DRAPIDJSON_HOME={}'.format(RAPIDJSON_DIST),
        '-DDIRENTWIN_HOME={}'.format(DIRENTWIN_DIST),
        '-DCMAKE_INSTALL_PREFIX={}'.format(DIST_DIR), 
        SRC_DIR
    ]
    cmake_run(cmd, BUILD_DIR)

    # build
    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    cmake_run(cmd, BUILD_DIR)

    # test  CMAKE_INSTALL_PREFIX\bin\sortmerna --version
    cmd = [ os.path.join(DIST_DIR, 'bin', 'sortmerna'), '--version' ]
    cmake_run(cmd, BUILD_DIR)

    # CMAKE_INSTALL_PREFIX\bin\sortmerna -h
    cmd = [ os.path.join(DIST_DIR, 'bin', 'sortmerna'), '-h' ]
    cmake_run(cmd, BUILD_DIR)

#END smr_build

def zlib_build(gen='Unix Makefiles', btype='Release', src=None, build=None, dist=None):
    '''
    @param btype  Build type Relase | Debug | ..
    @param ptype  Linkage type like statuc, dynamic, mixed
    '''
    SRC_DIR = src
    BUILD_DIR = build
    DIST_DIR = dist
    CMAKE_GEN = gen

    #INSTALL_BIN_DIR = os.path.join(CMAKE_INSTALL_PREFIX, 'bin')
    #INSTALL_INC_DIR = os.path.join(CMAKE_INSTALL_PREFIX, 'include')
    #INSTALL_LIB_DIR = os.path.join(CMAKE_INSTALL_PREFIX, 'lib')
    #INSTALL_MAN_DIR = os.path.join(CMAKE_INSTALL_PREFIX, 'share', 'man')
    #INSTALL_PKGCONFIG_DIR = os.path.join(CMAKE_INSTALL_PREFIX, 'share', 'pkgconfig')

    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DCMAKE_INSTALL_PREFIX={}'.format(DIST_DIR), 
        SRC_DIR
        #'-DINSTALL_BIN_DIR={}'.format(INSTALL_BIN_DIR),
        #'-DINSTALL_INC_DIR={}'.format(INSTALL_INC_DIR),
        #'-DINSTALL_LIB_DIR={}'.format(INSTALL_LIB_DIR),
        #'-DINSTALL_MAN_DIR={}'.format(INSTALL_MAN_DIR),
        #'-DINSTALL_PKGCONFIG_DIR={}'.format(INSTALL_PKGCONFIG_DIR),
    ]
    cmake_run(cmd, BUILD_DIR)

    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    cmake_run(cmd, BUILD_DIR)
#END zlib_build

def rapidjson_build(gen='Unix Makefiles', btype='Release', src=None, build=None, dist=None):
    '''
    '''
    SRC_DIR = src
    BUILD_DIR = build
    DIST_DIR = dist
    CMAKE_GEN = gen

    BUILD_EXAMPLES = 0
    BUILD_DOC = 0

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DRAPIDJSON_BUILD_EXAMPLES={}'.format(BUILD_EXAMPLES),
        '-DRAPIDJSON_BUILD_DOC={}'.format(BUILD_DOC),
        '-DCMAKE_INSTALL_PREFIX={}'.format(DIST_DIR), 
        SRC_DIR
    ]
    cmake_run(cmd, BUILD_DIR)

    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    cmake_run(cmd, BUILD_DIR)
#END rapidjson_build

def rocksdb_build(gen='Unix Makefiles', btype='Release', ptype='t3', src=None, build=None, dist=None, zlib=None):
    '''
    @param btype  Build type Release | Debug
    @param ptype  Linkage type on Windows t1 | t2 | t3
    '''
    SRC_DIR = src
    BUILD_DIR = build
    DIST_DIR = dist
    CMAKE_GEN = gen
    ZLIB_DIST = zlib
    
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

    pf = platform.platform()    

    if 'Windows' in pf:
        if ptype == "t3": print("Type 3 linkage: /MD + static ZLib")

        WITH_MD_LIBRARY = 1
        ROCKSDB_INSTALL_ON_WINDOWS = 1
        WITH_XPRESS = 1
        PORTABLE = 0

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        '-DOPTDBG={}'.format(OPTDBG),
        '-DPORTABLE={}'.format(PORTABLE),
        '-DWITH_MD_LIBRARY={}'.format(WITH_MD_LIBRARY),
        '-DWITH_RUNTIME_DEBUG={}'.format(WITH_RUNTIME_DEBUG),
        '-DWITH_GFLAGS={}'.format(WITH_GFLAGS),
        '-DZLIB_ROOT={}'.format(ZLIB_DIST),
        '-DWITH_ZLIB={}'.format(WITH_ZLIB),
        '-DWITH_XPRESS={}'.format(WITH_XPRESS),
        '-DWITH_TESTS={}'.format(WITH_TESTS),
        '-DWITH_TOOLS={}'.format(WITH_TOOLS),
        '-DROCKSDB_INSTALL_ON_WINDOWS={}'.format(ROCKSDB_INSTALL_ON_WINDOWS),
        '-DCMAKE_INSTALL_PREFIX={}'.format(DIST_DIR), 
        SRC_DIR
    ]
    cmake_run(cmd, BUILD_DIR)

    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    cmake_run(cmd, BUILD_DIR)
#END rocksdb_build

if __name__ == "__main__":
#    import pdb; pdb.set_trace()
    SMR = 'sortmerna'
    ZLIB = 'zlib'
    ROCKS = 'rocksdb'
    RAPID = 'rapidjson'
    DIRENT = 'dirent'

    pf = platform.platform()
    IS_WIN = 'Windows' in pf
    IS_WSL = 'Linux' in pf and 'Microsoft' in pf
    UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']

    CMAKE_GEN = 'Visual Studio 16 2019' if IS_WIN else 'Unix Makefiles'

    optpar = OptionParser()
    optpar.add_option('-n', '--name', dest='name', help='Module to build e.g. smr | zlib | rocksdb |')
    optpar.add_option('--clone', action="store_true", help='Perform git clone for the given name')
    optpar.add_option('-c', '--clean', action="store_true", help='clean build directory for the given name')
    optpar.add_option('--btype', dest='btype', default='release', help = 'Build type: release | debug')
    optpar.add_option('--pt_smr', dest='pt_smr', help = 'Sortmerna Linkage type t1 | t2 | t3')
    optpar.add_option('--pt_zlib', dest='pt_zlib', help = 'Zlib Linkage type t1 | t2 | t3')
    optpar.add_option('--pt_rocks', dest='pt_rocks', help = 'ROcksdb Linkage type t1 | t2 | t3')
    optpar.add_option('--winhome', dest='winhome', help='When building on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')

    (opts, args) = optpar.parse_args()

    UHOME_WIN = opts.winhome if IS_WSL else None

    if IS_WIN:
        if not opts.pt_smr: opts.pt_smr = 't3'
        SMR_SRC = os.path.join(UHOME, 'a01_code', SMR)
        SMR_BUILD = os.path.join(SMR_SRC, 'build')
        SMR_DIST = os.path.join(SMR_SRC, 'dist', opts.pt_smr, opts.btype)

        # zlib puts both Debug and Release at the same location => no btype
        if not opts.pt_zlib: opts.pt_zlib = 't1'
        ZLIB_SRC = os.path.join(UHOME, 'a01_libs', ZLIB)
        ZLIB_BUILD = os.path.join(ZLIB_SRC, 'build')
        ZLIB_DIST = os.path.join(ZLIB_SRC, 'dist', opts.pt_zlib)

        if not opts.pt_rocks: opts.pt_rocks = 't3'
        ROCKS_SRC = os.path.join(UHOME, 'a01_libs', ROCKS)
        ROCKS_BUILD = os.path.join(ROCKS_SRC, 'build')
        ROCKS_DIST = os.path.join(ROCKS_SRC, 'dist', opts.pt_rocks, opts.btype)

        # no binaries, so always build Release only
        RAPID_SRC = os.path.join(UHOME, 'a01_libs', RAPID)
        RAPID_BUILD = os.path.join(RAPID_SRC, 'build')
        RAPID_DIST = os.path.join(RAPID_SRC, 'dist')

        DIRENT_DIST = os.path.join(UHOME, 'a01_libs', DIRENT)
    elif IS_WSL:
        SMR_SRC = os.path.join(UHOME_WIN, 'a01_code', SMR)
        SMR_BUILD = os.path.join(UHOME, SMR, 'build')
        SMR_DIST = os.path.join(UHOME, SMR, 'dist')

        ZLIB_SRC = os.path.join(UHOME_WIN, 'a01_libs', ZLIB)
        ZLIB_BUILD = os.path.join(UHOME, ZLIB, 'build')
        ZLIB_DIST = os.path.join(UHOME, ZLIB, 'dist')

        ROCKS_SRC = os.path.join(UHOME_WIN, 'a01_libs', ROCKS)
        ROCKS_BUILD = os.path.join(UHOME, ZLIB, 'build')
        ROCKS_DIST = os.path.join(UHOME, ZLIB, 'dist')

        RAPID_SRC = os.path.join(UHOME_WIN, 'a01_libs', RAPID)
        RAPID_BUILD = os.path.join(UHOME, RAPID, 'build')
        RAPID_DIST = os.path.join(UHOME, RAPID, 'dist')
    else:
        SMR_SRC = os.path.join(UHOME, 'a01_code', SMR)
        SMR_BUILD = os.path.join(SMR_SRC, 'build')
        SMR_DIST = os.path.join(SMR_SRC, 'dist')

    if opts.name:
        if opts.name == SMR: 
            if opts.clone:
                smr_clone()
            elif opts.clean:
                clean(SMR_BUILD)
            else:
                smr_build(gen=CMAKE_GEN, btype=opts.btype, ptype=opts.pt_smr, 
                        src=SMR_SRC, build=SMR_BUILD, dist=SMR_DIST,
                        zlib=ZLIB_DIST, rocks=ROCKS_DIST, rapid=RAPID_DIST, dirent=DIRENT_DIST)
        elif opts.name == ZLIB: 
            if opts.clone:
                zlib_clone()
            else:
                zlib_build(gen=CMAKE_GEN, src=ZLIB_SRC, build=ZLIB_BUILD, dist=ZLIB_DIST)
        elif opts.name == RAPID: 
            if opts.clone:
                rapidjson_clone()
            else:
                rapidjson_build(gen=CMAKE_GEN, src=RAPID_SRC, build=RAPID_BUILD, dist=RAPID_DIST)
        elif opts.name == ROCKS: 
            if opts.clone:
                rocksdb_clone()
            elif opts.clean:
                clean(ROCKS_BUILD)
            else:
                rocksdb_build(gen=CMAKE_GEN, src=ROCKS_SRC, build=ROCKS_BUILD, dist=ROCKS_DIST, zlib=ZLIB_DIST)
        elif opts.name == DIRENT: 
            if opts.clone:
                dirent_clone()
        else: test()