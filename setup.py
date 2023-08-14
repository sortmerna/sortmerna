# @copyright 2016-2021  Clarity Genomics BVBA
# @copyright 2012-2016  Bonsai Bioinformatics Research Group
# @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla
#
# SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
# This is a free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SortMeRNA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
#
# contributors: Jenya Kopylova   jenya.kopylov@gmail.com
#			          Laurent Noé      laurent.noe@lifl.fr
#			          Pierre Pericard  pierre.pericard@lifl.fr
#			          Daniel McDonald  wasade@gmail.com
#			          Mikaël Salson    mikael.salson@lifl.fr
#			          Hélène Touzet    helene.touzet@lifl.fr
#			          Rob Knight       robknight@ucsd.edu

'''
file: build.py
created: Jul 31, 2019 Wed

Dependencies:
  conda install pyyaml jinja2 requests
'''
import os
import sys
import subprocess
import platform
from optparse import OptionParser
import requests
import tarfile
import re
import fileinput
import json
import yaml
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import zipfile
import time
import shutil
#from distutils.dir_util import copy_tree

# define platform
pf = platform.platform()
IS_WIN = 'Windows' in pf
IS_WSL = 'Linux' in pf and 'Microsoft' in pf # Windows Subsystem for Linux (WSL)
IS_LNX = 'Linux' in pf and not 'Microsoft' in pf
if   IS_WIN: MY_OS = 'WIN'
elif IS_WSL: MY_OS = 'WSL'
elif IS_LNX: MY_OS = 'LIN'
else:
    print('Unable to define the platform: {}'.format(pf))
    sys.exit(1)

modfunc = {} # holds refs to functions in this module

# globals
SMR     = 'sortmerna'
ZLIB    = 'zlib'
ROCKS   = 'rocksdb'
#RAPID   = 'rapidjson'
DIRENT  = 'dirent'
CMAKE   = 'cmake'
CONDA   = 'conda'
ALL     = 'all'
CCQUEUE = 'concurrentqueue'

MY_OS = None
ENV = None # WIN | WSL | LNX_AWS | LNX_TRAVIS

UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']
#CMAKE_GEN   = None
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
#RAPID_SRC   = None
#RAPID_BUILD = None
#RAPID_DIST  = None

# Concurrentqueue
CCQUEUE_SRC = None

def test():
    '''
    '''
    print(os.path.dirname(os.path.realpath(__file__)))
    print(os.getcwd())
#END test

def git_clone(url, repo_dir, commit='master', shallow=False, force=False):
    '''
    clone a git repo if not already existing

    :param str url    url to clone
    :param str repo_dir   directory where to clone
    '''
    ST = '[git_clone]'
    ret = {}
    # make sure parent directory where to clone exists
    if not os.path.exists(os.path.split(repo_dir)[0]):
        os.makedirs(os.path.dirname(repo_dir)) # create parent repo

    # check clone already exists
    if os.path.exists(repo_dir) and os.path.exists(os.path.join(repo_dir, '.git')):
        print('{} git clone already exists: {}'.format(ST, repo_dir))
        cmd = ['git', 'status']
        ret = proc_run(cmd, repo_dir)
    else: # do clone
        repo = os.path.basename(repo_dir) # e.g. sortmerna
        cmd = [ 'git', 'clone', url, repo_dir]
        ret = proc_run(cmd)
    return ret
#END git_clone  

def conda_install(dir=None, force=False, clean=False, **cfg):
    '''
    :param cfg   Config dictionary
    :param dir   installation root dir. Default: User Home

    pip install -U pyyaml
    pip install -U Jinja2
    '''
    ST = '[conda_install]'
    dir = UHOME if not dir else dir
    fsh = 'miniconda.sh'
    is_install_ok = False
    #os.chdir(dir)
    # check already installed
    bin_conda = os.path.join(dir, 'miniconda3', 'bin')
    if os.path.exists(bin_conda):
        cmd = ['python', '--version']
        ret = proc_run(cmd, bin_conda, True)
        if not ret['retcode']:
            print('{} Conda already installed: {} : {}'.format(ST, bin_conda, ret['stdout']))
            if force:
                print('{} TODO: Force specified - removing the existing installation: {}'.format(ST, bin_conda))
            return

    # download the installer if not already present
    if not os.path.exists(os.path.join(dir, fsh)):
        url_conda = cfg[CONDA]['url'][MY_OS]
        #if sys.version_info[1] < 7:
            # 'Mozilla/5.0 (Windows NT 6.1; Win64; x64)'
        #    req = urllib.request.Request(url_conda, headers={'User-Agent': 'Mozilla/5.0'})
        #else:
            # this works with conda but not standard python3
        #    req = url_conda

        print('{} Loading Conda from url: {}'.format(ST, url_conda))
        try:
            req = requests.get(url_conda, allow_redirects=True)
            #with urllib.request.urlopen(req) as surl:
            with open(fsh, 'wb') as fp:
                fp.write(req.content) # loads the file into memory before writing to disk. No good for very big files.
        except:
            print('{} Exception getting Conda distro: {} {}'.format(ST, sys.exc_info()[1], sys.exc_info()[0]))
            sys.exit(1)

    # run the installer
    if os.path.exists(os.path.join(dir, fsh)):
        cmd = ['bash', fsh, '-b']
        #cmd = ['bash', fsh, '-b', '-p', '{}/miniconda'.format(dir)]
        proc_run(cmd, dir)

        # delete the installer
        if clean:
            print('{} Deleting the installer {}'.format(ST, os.path.join(dir, fsh)))
            os.remove(os.path.join(dir, fsh))

        print('{} Installed conda in {}'.format(ST, os.path.join(dir, 'miniconda3')))

        # install packages required to use sortmerna's build.py
        print('{} Installing PyYaml package'.format(ST))
        bin_pip = os.path.join(bin_conda, 'pip')
        cmd = [bin_pip, 'install', 'pyyaml']
        proc_run(cmd, bin_conda)

        print('{} Installing Jinja2 package'.format(ST))
        cmd = [bin_pip, 'install', 'jinja2']
        proc_run(cmd, bin_conda)

        print('{} Installing NumPy package'.format(ST))
        cmd = [bin_pip, 'install', 'numpy']
        proc_run(cmd, bin_conda)

        print('{} Installing scikit-bio package'.format(ST))
        cmd = [bin_pip, 'install', 'scikit-bio']
        proc_run(cmd, bin_conda)
    else:
        print('{} Conda installer not found - likely failed to download'.format(ST))

    print('{} Done'.format(ST))
#END conda_install

def cmake_install(dir=None, force=False, **cfg):
    '''
    :param dict cfg configuration dict
    :param str dir installation directory. Default User Home

    Download and extract CMake release archive
    python build.py --name cmake
    '''
    ST = '[cmake_install]'
    is_installed = False
    cmake_cfg = cfg.get(CMAKE)
    url = cmake_cfg['url'][MY_OS]
    zipped = url.split('/')[-1] # tar.gz or zip
    env = cfg.get('env')
    cmake_home = cmake_cfg.get('home',{}).get(env)
    #cmake_home = val if val else os.path.join(UHOME, 'cmake-{}-win64-x64'.format(cmake_cfg.get('ver')))
    #                                          |_ default         
    # check already installed
    cmake_bin = '{}/bin/cmake'.format(cmake_home)
    if IS_WIN: cmake_bin = '{}.exe'.format(cmake_bin)
    if os.path.exists(cmake_bin):
        print('{} Cmake is already installed: {}'.format(ST, cmake_bin))
        if not force:
            is_installed = True

    if not is_installed:
        print('{} CMake not found at: {} - installing'.format(ST, cmake_bin))
        os.chdir(Path(cmake_home).parent.as_posix()) # navigate to the parent dir e.g. installation root (pathlib)
        # download the installer if not already present
        if not os.path.exists(zipped):
            # load file from URL
            req = requests.get(url, allow_redirects=True)
            try:
                #with urllib.request.urlopen(url) as surl:
                with open(zipped, 'wb') as fp:
                    fp.write(req.content) # loads all file into memory first before writing to disk. No good for very big files.
            except:
                print('{} Exception getting CMake distro: {} {}'.format(ST, sys.exc_info()[1], sys.exc_info()[0]))
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
            print('{} Installed CMake {}'.format(ST, cmake_bin))
        else:
            print('{} Failed to install CMake {}'.format(ST, cmake_bin))

        # copy binarties to HOME/bin to avoid setting the PATH
        #if IS_LNX:
        #    copy_tree('{}/bin'.format(cmake_home), '{}/bin'.format(UHOME))
    #os.environ('PATH') 
#END cmake_install

def gcc_install(**cfg):
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
    ST = '[gcc_install]'
    GCC_MINVER = 9
    is_gcc = False
    if shutil.which('gcc'):
        cmd = ['gcc', '--version']
        ret = proc_run(cmd, capture=True)
        if ret['retcode'] == 0:
            gccver = ret['stdout'].decode('utf-8').split('\n')[0].split()[3].split('.')[0]
            if int(gccver) >= GCC_MINVER:
                print(f'{ST} gcc version {gccver} found - OK')
                is_gcc = True
    if not is_gcc:
        print(f'{ST} gcc is not found. Please, install prior using this script')
#END gcc_install

def make_install(**cfg):
    '''
    sudo apt install make
    '''
#END make_install

def git_install(**cfg):
    '''
    '''
    cmd = ['apt-cache', 'policy', 'git']
    proc_run(cmd)

def clean(dir):
    '''
    :param str dir  Directory to clean

    remove content of the given directory.
    '''
    cmd = ['rm', '-rf', './*']
    proc_run(cmd, dir)
#END clean

def proc_run(cmd, cwd=None, capture=False):
    '''
    :param list cmd      command to execute
    :param str cwd       CWD
    :param bool capture  capture the output
    '''
    ST = '[proc_run]'
    ret = {'retcode':0, 'stdout':None, 'stderr':None}
    spr = '==========================================================='

    cwdp = cwd or Path().resolve()
    print('{} Running in {}:\n{}\n{}\n{}'.format(ST, cwdp, spr, ' '.join(cmd), spr))

    start = time.time()
    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe
    # check executable
    if shutil.which(cmd[0]):
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

        print("{} Run time: {}".format(ST, time.time() - start))
    else:
        msg = '{} Executable {} not found'.format(ST, cmd[0])
        ret['retcode'] = 1
        ret['stderr'] = msg
        print(msg)
    return ret
#END proc_run

def zlib_build(**kwargs):
    '''
    :param str url
    :param str path              local filesystem directory where to clone the repo
    :param str commit            git commit
    :param str cmake_build_type  Build type Relase | Debug | ..
    :param ptype  Linkage type like statuc, dynamic, mixed
    '''
    url = kwargs.get('url')
    path = kwargs.get('path')
    commit = kwargs.get('commit')
    shallow = kwargs.get('shallow')
    btype = kwargs.get('cmake_build_type', 'Release')
    gen = kwargs.get('cmake_gen', 'Ninja Multi-Config')
    git_clone(url, path, commit=commit) # URL_ZLIB
    # cmake configure
    # set CMAKE_INSTALL_PREFIX here because original zlib CMakeLists.txt uses it
    # to set cache variables used in installation at configure time (not a good solution).
    # Instead the relative paths should be used, and CMake will add the prefix during installation automatically).
    # https://discourse.cmake.org/t/cmake-install-prefix-not-work
    # https://github.com/madler/zlib/pull/170 <- PR to fix zlib CMakeLists.txt (never merged)
    pfx = kwargs.get('dist') or os.path.abspath(f'build/{path}/dist').replace('\\', '/')
    cmd = ['cmake', '-S', path, '-B', f'build/{path}', '-G', gen, '--fresh', '-D', f'CMAKE_INSTALL_PREFIX={pfx}']
    ret = proc_run(cmd)
    # cmake build
    if ret['retcode'] == 0:
        cmd = [ 'cmake', '--build', f'build/{path}', '--target', 'install', '--config', btype]
        proc_run(cmd)
    return ret
#END zlib_build

"""def rapidjson_build(btype='Release'):
    '''
    '''
    BUILD_EXAMPLES = 0
    BUILD_DOC = 0
    
    git_clone(URL_RAPID, LIB_DIR)

    cmd = [
        'cmake', '-G', CMAKE_GEN,
        f'-DRAPIDJSON_BUILD_EXAMPLES={BUILD_EXAMPLES}',
        f'-DRAPIDJSON_BUILD_DOC={BUILD_DOC}',
        f'-DCMAKE_INSTALL_PREFIX={RAPID_DIST}', 
        RAPID_SRC
    ]
    proc_run(cmd, RAPID_BUILD)

    cmd = [ 'cmake', '--build', '.', '--config', btype, '--target', 'install' ]
    proc_run(cmd, RAPID_BUILD)
#END rapidjson_build
"""

def rocksdb_modify_3party_zlib(link_type='t1', **kwargs):
    '''
    modify 'thirdparty.inc' config file provided with RocksDb distro to point to zlib installation
    only used on Windows

    :param str ptype  library linkage type as in 'env.jinja.yaml:rocksdb:link:win:t3'
    :param dict cfg   build configuration from 'env.jinja.yaml'

    '''
    ST = '[rocksdb_modify_3party_zlib]'
    if not IS_WIN:
        print('{} not used on Non-Windows'.format(ST))
        return

    print('{} fixing \'thirdparty.inc\' for linkage type [{}] on Windows'.format(ST, link_type))
    link_type = kwargs.get(ROCKS,{}).get('link.type', {}).get('windows',{}).get(link_type,{})
    zlib_rel = link_type.get('zlib.release')  # cfg[ROCKS]['link']['WIN'][ptype]['ZLIB_LIB_RELEASE']
    zlib_dbg = link_type.get('zlib.debug')  # cfg[ROCKS]['link']['WIN'][ptype]['ZLIB_LIB_DEBUG']
    file3p = os.path.join(kwargs.get(ROCKS).get('path'), 'thirdparty.inc')
    zlib_dist = kwargs.get('zlib',{}).get('dist')

    for line in fileinput.FileInput(file3p, inplace=True):
        if line.startswith('set(ZLIB_HOME'):
            line = re.sub(r'ZLIB_HOME .*\)', r'ZLIB_HOME {})'.format(zlib_dist), line, flags = re.M)
        if line.startswith('set(ZLIB_INCLUDE'):
            line = re.sub(r'ZLIB_INCLUDE .*\)', r'ZLIB_INCLUDE ${ZLIB_HOME}/include)', line, flags = re.M)
        if line.startswith('set(ZLIB_LIB_DEBUG'):
            line = re.sub(r'ZLIB_LIB_DEBUG .*\)', r'ZLIB_LIB_DEBUG ${{ZLIB_HOME}}/lib/{})'.format(zlib_dbg), line, flags = re.M)
        if line.startswith('set(ZLIB_LIB_RELEASE'):
            line = re.sub(r'ZLIB_LIB_RELEASE .*\)', r'ZLIB_LIB_RELEASE ${{ZLIB_HOME}}/lib/{})'.format(zlib_rel), line, flags = re.M)
        sys.stdout.write(line)
   
    #mo = re.search(r'ZLIB_HOME .*\)+?', txt)
    #txtn = re.sub(r'ZLIB_HOME .*\)+?', r'ZLIB_HOME {})'.format(ZLIB_DIST.replace('\\', '/')), txt, flags = re.M)
    #txtn = re.sub(r'ZLIB_INCLUDE .*\)+?', r'ZLIB_INCLUDE ${ZLIB_HOME}/include)', txtn, flags = re.M)
    #txtn = re.sub(r'ZLIB_LIB_DEBUG .*\)', r'ZLIB_LIB_DEBUG ${{ZLIB_HOME}}/lib/{})'.format(lib_dbg), txtn, flags = re.M)
    #txtn = re.sub(r'ZLIB_LIB_RELEASE .*\)', r'ZLIB_LIB_RELEASE ${{ZLIB_HOME}}/lib/{})'.format(lib_rel), txtn, flags = re.M)
#END rocksdb_fix_3party
    
def rocksdb_build(link_type='t1', **kwargs): # ver=None, btype='Release', ptype='t3', **cfg
    '''
    :param str btype  build type Release | Debug
    :param str ptype  Linkage type on Windows t1 | t2 | t3

    NOTE: on Windows 'thridparty.inc' file has to be modified.
    '''
    ST = '[rocksdb_build]'
    print(f'{ST} started')
    url = kwargs.get(ROCKS).get('url')
    path = kwargs.get(ROCKS).get('path')
    commit = kwargs.get(ROCKS).get('commit')
    shallow = kwargs.get(ROCKS).get('shallow')
    btype = kwargs.get(ROCKS).get('cmake_build_type', 'Release')
    ret = git_clone(url, path, commit=commit)

    cmd = ['git', 'checkout', commit] if commit  else ['git', 'checkout', 'master']
    ret = proc_run(cmd, path)

    cmd = ['git', 'rev-parse', 'HEAD']
    ret = proc_run(cmd, path, capture=True)
    chash = ret.get('stdout')
    print(f'{ST} building commit: {chash}')

    if IS_WIN:
        rocksdb_modify_3party_zlib(link_type, **kwargs)

    # copy presets file
    src = os.path.abspath('cmake/presets/CMakePresets_rocksdb.json')
    dst = os.path.abspath(f'{path}/CMakePresets.json')
    print(f'{ST} copying {src} -> {dst}')
    ret = shutil.copyfile('cmake/presets/CMakePresets_rocksdb.json', f'{path}/CMakePresets.json')

    dist_path = os.path.abspath(f'build/{path}/dist').replace('\\','/')
    cmd = [
        'cmake', '-S', path, '-B', f'build/{path}',
        '--preset', 'WIN_release',
        f'-DCMAKE_INSTALL_PREFIX={dist_path}',
        f'-DZLIB_ROOT=build/3rdparty/zlib/{btype}/dist',
        '--fresh'
    ]
    ret = proc_run(cmd)

    if ret['retcode'] == 0:
        cmd = [ 'cmake', '--build', f'build/{path}', '--config', btype, '--target', 'install' ]
        ret = proc_run(cmd, ROCKS_BUILD)
    else:
        print(f'{ST} failed to build')
    return ret
#END rocksdb_build

def smr_build(ver=None, btype='release', link_type='t1', **kwargs):
    '''
    build sortmerna using CMake with CMakePresets.json
    CMake flags mostly specified in presets - no need here

    :param str ver        git commit/tag - default master
    :param str btype      Build type Release | Debug
    :param str link_type  Linking type: t1 | t2 | t3
                            t1 all static
                            t2 static 3rd party + dynamic runtime
                            t3 all dynamic
    '''
    ST = '[smr_build]'

    if ver:
        cmd = ['git', 'checkout', ver]
        proc_run(cmd)

    # list available presets
    cmd = ['cmake', '--list-presets']
    ret = proc_run(cmd, capture=True)
    if bool(ret.get('retcode')):
        print(ret)
        return

    presets = []
    if ret.get('stdout'):
        ol = ret.get('stdout').decode('utf8').replace('"','').split('\n')
        presets = [l.strip() for l in ol if l.strip()][1:]

    # cmake preset e.g. LIN_release | WIN_release
    preset = f'LIN_{btype}' if IS_LNX or IS_WSL else f'WIN_{btype}'
    if preset not in presets:
        print('{} preset {} not found in available presets: {}'.format(ST, preset, presets))
        return
    
    # generate CMake configuration
    # cmake -S . --preset LIN_Release
    cmd = ['cmake', '-S', '.', '--preset', preset]
    if opts.vb:
        cmd.append('-DCMAKE_EXPORT_COMPILE_COMMANDS=1')
    if opts.loglevel:
        cmd.append('--loglevel={}'.format(opts.loglevel.upper()))
    elif opts.trace:
        cmd.append('--trace')
    ret = proc_run(cmd)

    # build and install
    cmd = [ 'cmake', '--build', '--preset', preset, '--target', 'install' ]
    ret = proc_run(cmd)

    # generate installation package
    if ret['retcode'] == 0:
        cmd = [ 'cmake', '--build', '--preset', preset, '--target', 'package' ]
        ret = proc_run(cmd)
    # test  CMAKE_INSTALL_PREFIX\bin\sortmerna --version
    #if ret['retcode'] == 0:
    #    SMR_EXE = 'sortmerna.exe' if IS_WIN else 'sortmerna'
    #    cmd = [ os.path.join(SMR_DIST, 'bin', SMR_EXE), '--version' ]
    #    ret = proc_run(cmd, SMR_BUILD)
    # CMAKE_INSTALL_PREFIX\bin\sortmerna -h
    #if ret['retcode'] == 0:
    #    cmd = [ os.path.join(SMR_DIST, 'bin', SMR_EXE), '-h' ]
    #    ret = proc_run(cmd, SMR_SRC)
    return ret
#END smr_build

def concurrentqueue_build(**kwargs):
    '''
    a single header file - just clone and use
    '''
    url = kwargs.get(CCQUEUE).get('url')
    path = kwargs.get(CCQUEUE).get('path')
    commit = kwargs.get(CCQUEUE).get('commit')
    shallow = kwargs.get(CCQUEUE).get('shallow')
    ret = git_clone(url, path, commit=commit, shallow=shallow)
    return ret
#END concurrentqueue_build

def env_check(**cfg):
    '''
    verify all pre-requisit packages/tools are present

    :param dict cfg
    '''
    for env in cfg.get('env.config',[]):
        if ENV in env.get('env', []):
            for pkg in env.get('packages',[]):
                funstr = f'{pkg}_install'
                fun = modfunc.get(funstr)
                ret = fun(**cfg)
            condapkgs = ' '.join(env.get('conda.packages',[]))
            if condapkgs:
                cmd =  f'conda install {condapkgs}'
                ret = proc_run(cmd)
            break
#END env_check

modfunc['git_install'] = git_install
modfunc['conda_install'] = conda_install
modfunc['gcc_install'] = gcc_install
modfunc['cmake_install'] = cmake_install

if __name__ == "__main__":
    '''
    python setup.py -n sortmerna [-e WIN] [--] [--cmake-preset <preset>]
    python setup.py -n sortmerna -e LNX_VBox_Ubuntu_1804   build smr
    python setup.py -n conda -e LNX_VBox_Ubuntu_1804       install conda
    python setup.py -n cmake -e LNX_VBox_Ubuntu_1804       install cmake
    python setup.py -n rocksdb -e LNX_VBox_Centos_77 -c    build rocksdb
    python setup.py --name cmake --envn WIN [--env scripts/env_non_git.yaml]
        --winhome /mnt/c/Users/biocodz --btype debug
    '''
    ST = '[setup.py:__main__]'
    #import pdb; pdb.set_trace()
    is_opts_ok = True

    # options
    optpar = OptionParser()
    optpar.add_option('-n', '--name', dest='name', 
        help='Package to build/process e.g. sortmerna | zlib | rocksdb | rapidjson | conda | all')
    optpar.add_option('--zlib-dist', dest='zlib_dist', help='ZLib installation directory')
    optpar.add_option('--rocksdb-dist', dest='rocksdb_dist', help='ROcksDB installation directory')
    optpar.add_option('--concurrentqueue-dist', dest='concurrentqueue_dist', help='concurrentqueue installation directory')
    optpar.add_option('-e', '--envn', dest='envname', 
        help=('Name of environment: WIN | WSL | LIN .. see env.jinja:env.list'))
    optpar.add_option('--clone', action="store_true", help='Perform git clone for the given name')
    optpar.add_option('-c', '--clean', action="store_true", help='clean build directory for the given name')
    optpar.add_option('-b', '--btype', dest='btype', default='release', help = 'Build type: release | debug')
    optpar.add_option('--pt_smr', dest='pt_smr', default='t1', help = 'Sortmerna Linkage type t1 | t2 | t3')
    optpar.add_option('--pt_zlib', dest='pt_zlib', help = 'Zlib Linkage type t1 | t2 | t3')
    optpar.add_option('--pt_rocks', dest='pt_rocks', help = 'Rocksdb Linkage type t1 | t2 | t3')
    optpar.add_option('--rocks3p', dest='rocks3p', help='Fix thirdparty.inc when building Rocksb')
    optpar.add_option('--winhome', dest='winhome', 
                help='when building on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    optpar.add_option('--trace', action="store_true", help='Run cmake with --trace')
    optpar.add_option('--loglevel', dest='loglevel', help = 'Cmake log level')
    optpar.add_option('--vb', action="store_true", help='Export compile commands')
    optpar.add_option('--env', dest='envfile', help='Env configuration file.')
    optpar.add_option('--config', dest='config', help='Build configuration file.')
    optpar.add_option('--build-dir', dest='build_dir', help='Build directory.')
    optpar.add_option('--dist-dir', dest='dist_dir', help='Distro directory.')
    optpar.add_option('--local-linux', dest='local_linux', action='store_true', help='Perform the build on local source files in the current working directory.')
    optpar.add_option('--cmake-preset', dest='cmake_preset', help='CMake preset')
    (opts, args) = optpar.parse_args()

    UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']

    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print(f'{ST} Current dir: {cur_dir}')

    IS_TODO = True # old functionality to adjust
    if not IS_TODO and opts.local_linux:
        cur_wdir = os.getcwd()

        LIB_DIR = cur_wdir

        #CMAKE_GEN = 'Unix Makefiles'

        # SMR
        SMR_SRC   = cur_wdir
        SMR_BUILD = f'{cur_wdir}/build'
        SMR_DIST  = f'{cur_wdir}/dist'
        SMR_VER = None

        # ZLIB
        ZLIB_SRC = f'{cur_wdir}/{ZLIB}'
        ZLIB_BUILD = f'{ZLIB_SRC}/build'
        ZLIB_DIST = f'{ZLIB_SRC}/dist'

        # ROCKSDB
        ROCKS_SRC = f'{cur_wdir}/{ROCKS}'
        ROCKS_BUILD = f'{ROCKS_SRC}/build'
        ROCKS_DIST = f'{ROCKS_SRC}/dist'
        ROCKS_VER = None

        CCQUEUE_SRC = f'{cur_wdir}/{CCQUEUE}' # CONCURRENTQUEUE

        env = {}

    elif not IS_TODO:
        # check env.yaml. If no env file specified, try the current directory
        envfile = os.path.join(cur_dir, 'env.jinja') if not opts.envfile else opts.envfile
        if not os.path.exists(envfile):
            print(f'{ST} No environment config file found. Please, provide one using \'--env\' option')
            sys.exit(1)

        if not opts.envname:
            ENV = MY_OS
            print(f'{ST} --envn was not specified - using {ENV}')
        else:
            ENV = opts.envname

        # load properties from env.jinja.yaml
        print(f'{ST} Using Environment configuration file: {envfile}')
        env_jj = Environment(loader=FileSystemLoader(os.path.dirname(envfile)), trim_blocks=True, lstrip_blocks=True)
        env_template = env_jj.get_template(os.path.basename(envfile))
        #   render env.jinja template
        vars = {'UHOME': UHOME, 'WINHOME': opts.winhome, 'ENV': ENV} if IS_WSL else {'UHOME': UHOME, 'ENV': ENV}
        env_str = env_template.render(vars)
        env = yaml.load(env_str, Loader=yaml.FullLoader)

        if not opts.envname:
            envl = env.get('env.list')
            print(f'{ST} available environments: {envl}')

        libdir = env.get('LIB_DIR', {}).get(ENV)
        LIB_DIR =  libdir if libdir else UHOME

        #URL_ZLIB   = env[ZLIB]['url']
        #URL_ROCKS  = env[ROCKS]['url']
        #URL_DIRENT = env[DIRENT]['url']
        #URL_RAPID  = env[RAPID]['url']
        #URL_SMR    = env[SMR]['url']
        #URL_CONCURRENTQUEUE = env[CCQUEUE]['url']

        #CMAKE_GEN = env[CMAKE]['generator'][MY_OS]

        # SMR
        SMR_SRC   = env.get(SMR,{}).get('src',{}).get(ENV)
        SMR_SRC   = f'{UHOME}/sortmerna' if not SMR_SRC else SMR_SRC
        SMR_BUILD = opts.build_dir if opts.build_dir else env.get(SMR,{}).get('build',{}).get(ENV)
        SMR_BUILD = f'{SMR_SRC}/build' if not SMR_BUILD else SMR_BUILD
        SMR_DIST  = env.get(SMR,{}).get('dist',{}).get(ENV)
        SMR_DIST  = f'{SMR_SRC}/dist' if not SMR_DIST else SMR_DIST
        SMR_VER   = env.get(SMR,{}).get('ver')

        # ZLIB
        val = env.get(ZLIB,{}).get('src',{}).get(ENV)
        ZLIB_SRC = val if val else f'{LIB_DIR}/{ZLIB}'
        val = env.get(ZLIB,{}).get('build',{}).get(ENV)
        ZLIB_BUILD = val if val else f'{ZLIB_SRC}/build'
        val = env.get(ZLIB,{}).get('dist',{}).get(ENV)
        ZLIB_DIST = val if val else f'{ZLIB_SRC}/dist'

        # ROCKSDB
        val = env.get(ROCKS,{}).get('src',{}).get(ENV)
        ROCKS_SRC = val or f'{LIB_DIR}/{ROCKS}'
        val = env.get(ROCKS,{}).get('build',{}).get(ENV)
        ROCKS_BUILD = val or f'{ROCKS_SRC}/build'
        val = env.get(ROCKS,{}).get('dist',{}).get(ENV)
        ROCKS_DIST = val or f'{ROCKS_SRC}/dist'
        ROCKS_VER = env.get(ROCKS, {}).get('ver')

        # RAPIDJSON
        # no binaries, so always build Release only
        """val = env.get(RAPID,{}).get('src',{}).get(ENV)
        RAPID_SRC = val if val else f'{LIB_DIR}/{RAPID}'
        val = env.get(RAPID,{}).get('build',{}).get(ENV)
        RAPID_BUILD = val if val else f'{RAPID_SRC}/build'
        val = env.get(RAPID,{}).get('dist',{}).get(ENV)
        RAPID_DIST = val if val else f'{RAPID_SRC}/dist'
        """

        # CONCURRENTQUEUE
        val = env.get(CCQUEUE,{}).get('src',{}).get(ENV)
        CCQUEUE_SRC = val if val else f'{LIB_DIR}/{CCQUEUE}'

        if 'WIN' == ENV:
            SMR_DIST  = SMR_DIST + '/{}/{}'.format(opts.pt_smr, opts.btype)

            # zlib puts both Debug and Release at the same location => no btype
            opts.pt_zlib = 't1' if not opts.pt_zlib else opts.pt_zlib
            ZLIB_DIST  = ZLIB_DIST + '/{}'.format(opts.pt_zlib) 

            opts.pt_rocks = 't3' if not opts.pt_rocks else opts.pt_rocks
            ROCKS_DIST  = ROCKS_DIST + '/{}/{}'.format(opts.pt_rocks, opts.btype)

            val = env.get(DIRENT, {}).get('src', {}).get(ENV)
            DIRENT_SRC = val if val else f'{LIB_DIR}/{DIRENT}'
            val = env.get(DIRENT, {}).get('dist')
            DIRENT_DIST = val.get(ENV) if val and isinstance(val, dict) else DIRENT_SRC

    print(f'{ST} loading configuration from 3rdparty.jinja')
    envjn = Environment(loader=FileSystemLoader('.'), trim_blocks=True, lstrip_blocks=True)
    thirdparty_jinja = envjn.get_template('3rdparty.jinja')
    # zlib
    vars = { 'zlib_dist': opts.zlib_dist} if opts.zlib_dist else {}
    # rocksdb
    if opts.rocksdb_dist:
        vars['rocksdb_dist'] = opts.rocksdb_dist
    # concurrentqueue
    if opts.concurrentqueue_dist:
        vars['concurrentqueue_dist'] = opts.concurrentqueue_dist

    jinja_str = thirdparty_jinja.render(vars)
    third_party_data = yaml.load(jinja_str, Loader=yaml.FullLoader) # render jinja template
    # set zlib and rocksdb default roots if not provided through options
    if not third_party_data.get('zlib',{}).get('dist'):
        pp = third_party_data['zlib']['path']
        third_party_data['zlib']['dist'] = os.path.abspath(f'build/{pp}/dist').replace('\\', '/')
    if not third_party_data.get('rocksdb',{}).get('dist'):
        pp = third_party_data['rocksdb']['path']
        third_party_data['rocksdb']['dist'] = os.path.abspath(f'build/{pp}/dist').replace('\\', '/')
    if not third_party_data.get('concurrentqueue',{}).get('dist'):
        pp = third_party_data['concurrentqueue']['path']
        third_party_data['concurrentqueue']['dist'] = os.path.abspath(f'{pp}').replace('\\', '/')
    
    opts.name = opts.name or ALL
    if opts.name in [ALL]:
        if opts.clean:
            clean(ZLIB_BUILD)
            clean(ROCKS_BUILD)
            clean(SMR_BUILD)
        else:
            if not opts.cmake_preset:
                opts.cmake_preset = 'WIN_release' if IS_WIN else 'LIN_release'
            kw = third_party_data.get(ZLIB,{})
            ret = zlib_build(**kw)
            if ret['retcode'] == 0:
                ret = rocksdb_build(**third_party_data) # ROCKS_VER, cfg=env
            if ret['retcode'] == 0:
                ret = concurrentqueue_build(**third_party_data)
            if ret['retcode'] == 0:
                ret = smr_build(**third_party_data)
    elif opts.name == ZLIB:
        kw = third_party_data.get(ZLIB,{})
        zlib_build(**kw)
    elif opts.name == ROCKS:
        if opts.clean:
            clean(ROCKS_BUILD)
        else:
            rocksdb_build(**third_party_data)
    elif opts.name in [SMR]:
        if opts.clean:
            clean(SMR_BUILD)
        else:
            smr_build(**third_party_data)
    elif opts.name == CCQUEUE: 
        concurrentqueue_build(**third_party_data)
    elif opts.name == DIRENT: 
        url = third_party_data[DIRENT].get('url')
        path = third_party_data[DIRENT].get('path')
        commit = third_party_data[DIRENT].get('commit')
        shallow = third_party_data[DIRENT].get('shallow')
        git_clone(url, path, commit=commit, shallow=shallow) 
    elif opts.name == CMAKE: 
        if opts.clone:
            git_clone(env[CMAKE]['url'], LIB_DIR)
        cmake_install(**env) 
    elif opts.name == CONDA: 
        conda_install(**env) 
    #elif opts.name == RAPID: 
        #rapidjson_build()
    else: test()
#END main

#END END
