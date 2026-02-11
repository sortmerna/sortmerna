# @copyright 2016-2023  Clarity Genomics BVBA
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
file: setup.py
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
#import json
import yaml
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import zipfile
import time
import shutil
from io import BytesIO

#pkn = __name__.split('.')[-1]
# define platform
pf = platform.platform()
IS_WIN = 'Windows' in pf
IS_WSL = 'Linux' in pf and 'Microsoft' in pf # Windows Subsystem for Linux (WSL)
IS_LNX = 'Linux' in pf and not 'Microsoft' in pf
IS_OSX = 'macOS' in pf
MY_OS = 'WIN' if IS_WIN else 'LIN' if IS_LNX else 'WSL' if IS_WSL else 'OSX' if IS_OSX else None
if not MY_OS:
    print(f'[setup.py] Unexpected platform: {pf}')
    sys.exit(1)
else:
    print(f'[setup.py] running on: {pf}, id: {MY_OS}')

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

ENV = None # WIN | WSL | LNX_AWS | LNX_TRAVIS
UHOME = os.environ['USERPROFILE'] if IS_WIN else os.environ['HOME']

def test():
    '''
    '''
    print(os.path.dirname(os.path.realpath(__file__)))
    print(os.getcwd())
#END test
    
def load_tar(url, tgtd):
    '''
    '''
    ST = '[load_tar]'
    ret = {'result': True}
    # check path already exists
    if Path(tgtd).exists():
        msg = f'{ST} directory already exists {tgtd}. Skipping download'
        print(msg)
        ret['stdout'] = msg
    else:
        try:
            with requests.get(url) as res:
                tf = tarfile.open(mode='r:gz', fileobj=BytesIO(res.content))
                pp = Path(tgtd).parent / tf.firstmember.name
                msg = f'{ST} extracting into: {pp}'
                tf.extractall(Path(tgtd).parent) # -> e.g. ./3rdparty/rocksdb-7.10.2 35.7 MB
                ppn = pp.replace(Path(tgtd))
                msg += f'\n{ST} moved {pp} -> {ppn}'
                ret['stdout'] = msg
                print(msg)
        except Exception as ex:
            ret['result'] = False
            ret['stderr'] = f'{ex}'
    return ret

def git_clone(url:str, 
              repo_dir:str, 
              shallow:bool=False, 
              force:bool=False) -> tuple:
    '''
    clone a git repo if not already existing
    args:
      - url        url to clone
      - repo_dir   directory where to clone
    '''
    ST = '[git_clone]'
    rcode, sout, eout = 0, None, None
    # make sure parent directory where to clone exists
    if not Path(repo_dir).parent.exists():
        msg = f'{ST} creating {Path(repo_dir).parent}'
        print(msg)
        Path(repo_dir).parent.mkdir(parents=True) # create parent repo

    # check clone already exists
    if Path(repo_dir).exists() and (Path(repo_dir) / '.git').exists():
        msg = f'{ST} git clone already exists: {repo_dir}'
        print(msg)
        cmd = ['git', 'status']
        rcode, sout, eout = proc_run(cmd, repo_dir, capture=True)
    else: # do clone
        cmd = [ 'git', 'clone', url, repo_dir]
        rcode, sout, eout = proc_run(cmd)
    return rcode, sout, eout
#END git_clone  

def conda_install(dir:str=None, force:bool=False, clean:bool=False, **cfg):
    '''
    args:
      - cfg   Config dictionary
      - dir   installation root dir. Default: User Home

    pip install -U pyyaml
    pip install -U Jinja2
    '''
    ST = '[conda_install]'
    dir = dir or UHOME
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
    TODO
    sudo apt install make
    '''
    ...
#END make_install

def git_install(**cfg):
    '''
    '''
    cmd = ['apt', 'policy', 'git']
    proc_run(cmd)

def clean(dir:str):
    '''
    remove content of the given directory.
    args:
      - dir  Directory to clean
    '''
    cmd = ['rm', '-rf', './*']
    rcode, sout, eout = proc_run(cmd, dir)
    return rcode, sout, eout
#END clean

def proc_run(cmd:list, cwd:str=None, capture:bool=False) -> tuple:
    '''
    args:
      - cmd      command to execute
      - cwd      CWD
      - capture  capture the output
    '''
    ST = '[proc_run]'
    rcode, sout, eout = True, None, None
    spr = '==========================================================='

    cwd = cwd or Path().resolve()
    cmdstr = ' '.join(cmd)
    msg = f'{ST} Running in {cwd}:\n{spr}\n{cmdstr}\n{spr}'
    print(msg)

    start = time.time()
    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe
    # check executable
    if shutil.which(cmd[0]):
        try:
            if cwd and not Path(cwd).exists():
                Path(cwd).mkdir(parents=True)
            cwd = cwd or '.'
            if sys.version_info[1] > 6:
                proc = subprocess.run(cmd, cwd=cwd, capture_output=capture)
            else:
                proc = subprocess.run(cmd, cwd=cwd)

            if proc.returncode:
                for info in sys.exc_info():
                    print(info)
            rcode = proc.returncode
            if capture:
                if proc.stdout:
                    sout = proc.stdout.decode('utf-8').strip()
                    print(f'{ST} {sout}')
                if proc.stderr:
                    eout = proc.stderr.decode('utf-8').strip()
                    print(f'{ST} {eout}')
        except OSError as err:
            print(err)
            rcode = 1
            sout = str(err)
        except Exception as ex:
            msg = f'{ST} {ex}'
            print(msg)
            rcode = 1
            eout = str(ex)

        msg = f"{ST} Run time: {time.time() - start}"
        print(msg)
    else:
        eout = f'{ST} Executable {cmd[0]} not found'
        rcode = 1
        print(eout)
    return rcode, sout, eout
#END proc_run

def zlib_build(**kw):
    '''
    args:
      - url
      - path str              local filesystem directory where to clone the repo
      - commit str            git commit
      - cmake_build_type str  Build type Relase | Debug | ..
      - ptype str             Linkage type like statuc, dynamic, mixed
      - is_git bool           use git as source. Otherwise - archive (.tar.gz)
    '''
    outl, errl = [], []
    is_git = kw[ZLIB].get('is_git', False)
    is_checkout = kw[ZLIB].get('is_checkout', True)
    is_conda_cpp = kw.get('conda_cpp', False)
    url = kw[ZLIB].get('url') if is_git else kw[ZLIB].get('url2')
    src = kw[ZLIB].get('src')
    sysroot = kw.get('sysroot') if is_conda_cpp else None
    build_dir = Path(sysroot) / 'sortmerna' / kw[ZLIB].get('build') if is_conda_cpp and sysroot \
                else kw[ZLIB].get('build')
    dist_dir = Path(sysroot) / 'sortmerna' / kw[ZLIB].get('dist') if is_conda_cpp and sysroot \
                else kw[ZLIB].get('dist') or Path(f'build/{src}/dist').absolute()
    commit = kw[ZLIB].get('commit') or 'master'
    #shallow = kw[ZLIB].get('shallow')
    btype = kw[ZLIB].get('cmake_build_type', 'Release')
    gen = kw[ZLIB].get('cmake_gen', 'Ninja Multi-Config')
    rcode, sout, eout = git_clone(url, src) if is_git else load_tar(url, path) # URL_ZLIB
    outl.append(sout)
    errl.append(eout)

    if is_git and is_checkout:
        cmd = ['git', 'checkout', commit]
        proc_run(cmd, src)
    if is_git:
        cmd = ['git', 'rev-parse', 'HEAD']
        rcode, sout, eout = proc_run(cmd, src, capture=True)
        chash = sout
        print(f'{ST} building commit: {chash}')
        
    # cmake configure
    # set CMAKE_INSTALL_PREFIX here because original zlib CMakeLists.txt uses it
    # to set cache variables used in installation at configure time (not a good solution).
    # Instead the relative paths should be used, and CMake will add the prefix during installation automatically).
    # https://discourse.cmake.org/t/cmake-install-prefix-not-work
    # https://github.com/madler/zlib/pull/170 <- PR to fix zlib CMakeLists.txt (never merged)
    cmake_min = '3.5'
    #pfx = kwargs.get('dist') or os.path.abspath(f'build/{path}/dist').replace('\\', '/')
    cmd = ['cmake', '-S', src, '-B', str(build_dir), 
           '-G', gen, '--fresh', 
           f'-DCMAKE_INSTALL_PREFIX={str(dist_dir)}', 
           f'-DCMAKE_POLICY_VERSION_MINIMUM={cmake_min}']
    rcode, sout, eout = proc_run(cmd, capture=True)
    outl.append(sout)
    errl.append(eout)
    # cmake build
    if rcode == 0:
        cmd = [ 'cmake', '--build', str(build_dir), '--target', 'install', '--config', btype]
        rcode, sout, eout = proc_run(cmd, capture=True)
        outl.append(sout)
        errl.append(eout)
    return rcode, outl, errl
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

def rocksdb_modify_3party_zlib(link_type:str='t1', **kwargs):
    '''
    modify 'thirdparty.inc' config file provided with RocksDb distro to point to zlib installation
    only used on Windows

    args:
      - ptype  library linkage type as in 'env.jinja.yaml:rocksdb:link:win:t3'
      - cfg    build configuration from 'env.jinja.yaml'

    '''
    ST = '[rocksdb_modify_3party_zlib]'
    if not IS_WIN:
        print('{} not used on Non-Windows'.format(ST))
        return

    print(f'{ST} fixing \'thirdparty.inc\' for linkage type \'{link_type}\' on Windows')
    link_type = kwargs.get(ROCKS,{}).get('link.type', {}).get('windows',{}).get(link_type,{})
    zlib_rel = link_type.get('zlib.release')  # cfg[ROCKS]['link']['WIN'][ptype]['ZLIB_LIB_RELEASE']
    zlib_dbg = link_type.get('zlib.debug')  # cfg[ROCKS]['link']['WIN'][ptype]['ZLIB_LIB_DEBUG']
    file3p = Path(kwargs.get(ROCKS).get('src')) / 'thirdparty.inc'
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
    
def rocksdb_build(link_type:str='t1', **kw) -> tuple: # ver=None, btype='Release', ptype='t3', **cfg
    '''
    args:
      - btype  build type Release | Debug
      - ptype  Linkage type on Windows t1 | t2 | t3

    NOTE: on Windows 'thridparty.inc' file has to be modified.
    '''
    ST = '[rocksdb_build]'
    print(f'{ST} started')
    is_conda_cpp = kw.get('conda_cpp', False)
    sysroot, eout = get_sysroot() if is_conda_cpp else None, None
    if is_conda_cpp and sysroot:
        sysroot = Path(sysroot).resolve()
    src = kw.get(ROCKS).get('src')
    build_dir = sysroot / 'sortmerna' / kw[ROCKS].get('build') if is_conda_cpp and sysroot \
                else kw[ROCKS].get('build')
    dist_dir = sysroot / 'sortmerna' / kw[ROCKS].get('dist') if is_conda_cpp and sysroot \
                else kw[ROCKS].get('dist') or Path(f'build/{src}/dist').absolute()
    zlib_src = kw[ZLIB].get('src')
    zlib_dist = sysroot / 'sortmerna' / kw[ZLIB].get('dist') if is_conda_cpp and sysroot \
                 else Path(kw[ZLIB].get('dist')).absolute() \
                    or Path(f'build/{zlib_src}/dist').absolute()
    is_git = kw.get(ROCKS).get('is_git', False)
    is_checkout = kw.get(ROCKS).get('is_checkout', False)
    url = kw.get(ROCKS).get('url') if is_git else kw.get(ROCKS).get('url2')
    commit = kw.get(ROCKS).get('commit')
    #shallow = kw.get(ROCKS).get('shallow')
    btype = kw.get(ROCKS).get('cmake_build_type', 'Release')
    presets_file = kw.get(ROCKS).get('preset')
    rcode, sout, eout = git_clone(url, src) if is_git else load_tar(url, src)

    if is_git and is_checkout:
        cmd = ['git', 'checkout', commit] if commit  else ['git', 'checkout', 'master']
        rcode, sout, eout = proc_run(cmd, src)
    if is_git:
        cmd = ['git', 'rev-parse', 'HEAD']
        rcode, sout, eout = proc_run(cmd, src, capture=True)
        chash = sout
        print(f'{ST} building commit: {chash}')

    if IS_WIN:
        rocksdb_modify_3party_zlib(link_type, **kw)
        
    # check cmake version
    cmd = ['cmake', '--version']
    rcode, sout, eout = proc_run(cmd, capture=True)
    msg = f'{ST} {sout}'
    print(msg)

    # copy presets file
    presets = Path(presets_file).absolute()
    dst = Path(f'{src}/CMakePresets.json').absolute()
    print(f'{ST} copying {presets} -> {dst}')
    shutil.copyfile(presets_file, f'{src}/CMakePresets.json')

    bt = btype.lower()
    cmd = [
        'cmake', '-S', src, '-B', Path(build_dir).as_posix(),
        '--preset', f'{MY_OS}_{bt}',
        f'-DCMAKE_INSTALL_PREFIX:PATH={Path(dist_dir).as_posix()}',
        f'-DZLIB_ROOT:PATH={Path(zlib_dist).as_posix()}',
        f'-DCMAKE_POLICY_DEFAULT_CMP0074=NEW'
    ]
    # 20260210 Tue  warning: Policy CMP194 is not set: MSVC is not an assembler for language ASM.
    if IS_WIN:
        cmd.append(f'-DCMAKE_POLICY_DEFAULT_CMP194=OLD')
    if is_conda_cpp and sysroot:
        toolchain = Path(kw.get('toolchain')) if kw.get('toolchain') else Path.cwd() / 'cmake/conda_tc.cmake'
        if toolchain.exists():
            cmd.append(f'-DCMAKE_TOOLCHAIN_FILE={str(toolchain)}')
        else:
            msg = f'{ST} file not found {toolchain}'
            print(msg)
            raise ValueError(msg)
    cmd.append('--fresh')
    rcode, sout, eout = proc_run(cmd)

    if rcode == 0:
        cmd = [ 'cmake', '--build', str(build_dir), '--config', btype, '--target', 'install']
        rcode, sout, eout = proc_run(cmd)
    else:
        print(f'{ST} failed to build')
    return rcode, sout, eout
#END rocksdb_build

def smr_build(ver:str=None, 
              btype:str='release', 
              is_checkout:bool=False, 
              link_type:str='t1', 
              **kw):
    '''
    build sortmerna using CMake with CMakePresets.json
    CMake flags mostly specified in presets - no need here
    args:
      - ver        git commit/tag - default master
      - btype      Build type Release | Debug
      - link_type  Linking type: t1 | t2 | t3
                            t1 all static
                            t2 static 3rd party + dynamic runtime
                            t3 all dynamic
    '''
    ST = '[smr_build]'
    is_conda_cpp = kw.get('conda_cpp', False)
    sysroot, eout = get_sysroot() if is_conda_cpp else None, None
    if is_conda_cpp and sysroot:
        sysroot = Path(sysroot).resolve()
    build_dir = sysroot / 'sortmerna/build' if is_conda_cpp and sysroot else 'build'
    rocksdb_src = kw.get(ROCKS).get('src')
    rocksdb_dist = sysroot / 'sortmerna' / kw[ROCKS].get('dist') if is_conda_cpp and sysroot \
                        else kw[ROCKS].get('dist') or Path(f'build/{rocksdb_src}/dist').absolute()
    zlib_dist = sysroot / 'sortmerna' / kw[ZLIB].get('dist') if is_conda_cpp and sysroot \
                    else kw[ZLIB].get('dist') or Path(f"build/{kw[ZLIB].get('src')}/dist").absolute()
    conque_home = kw[CCQUEUE].get('dist') or kw[CCQUEUE].get('src')
    install_dir = Path(build_dir) / 'dist' if is_conda_cpp and sysroot else 'dist'
    
    if is_checkout and ver and Path('.git').exists():
        cmd = ['git', 'checkout', ver]
        proc_run(cmd)

    # list available presets
    cmd = ['cmake', '--list-presets']
    rcode, sout, eout = proc_run(cmd, capture=True)
    if rcode > 0:
        print(eout)
        return

    presets = []
    if sout:
        ol = sout.replace('"','').split('\n')
        presets = [l.strip() for l in ol if l.strip()][1:]

    # cmake preset e.g. LIN_release | WIN_release
    preset = f'LIN_{btype}' if IS_LNX or IS_WSL or IS_OSX else f'WIN_{btype}'
    if preset not in presets:
        msg = f'{ST} preset {preset} not found in available presets: {presets}'
        print(msg)
        return 1, None, msg
    
    # generate CMake configuration
    # cmake -S . --preset LIN_Release
    cmd = ['cmake', '-S', '.', '-B', str(build_dir), '--preset', preset]
    if is_conda_cpp and sysroot:
        toolchain = Path(kw.get('toolchain')) if kw.get('toolchain') else Path.cwd() / 'cmake/conda_tc.cmake'
        if toolchain.exists():
            cmd.append(f'-DCMAKE_TOOLCHAIN_FILE={str(toolchain)}')
        else:
            msg = f'{ST} file not found {toolchain}'
            print(msg)
            raise ValueError(msg)
        cmd.append(f'-DZLIB_ROOT={zlib_dist}')
        cmd.append(f'-DROCKSDB_DIST={rocksdb_dist}')
        cmd.append(f'-DCONCURRENTQUEUE_HOME={conque_home}')
        cmd.append(f'-DCMAKE_INSTALL_PREFIX={str(install_dir)}')
    if kw.get('vb'):
        cmd.append('-DCMAKE_EXPORT_COMPILE_COMMANDS=1')
    if kw.get('loglevel'):
        cmd.append('--loglevel={}'.format(kw.get('loglevel').upper()))
    elif kw.get('trace'):
        cmd.append('--trace')
    cmd.append('--fresh')
    rcode, sout, eout = proc_run(cmd, capture=True)
    if rcode:
        return rcode,sout, eout

    # build and install
    cmd = [ 'cmake', '--build', str(build_dir), '--preset', preset, '--target', 'install' ]
    rcode, sout, eout = proc_run(cmd, capture=True)
    if rcode:
        return rcode,sout, eout

    # generate installation package
    cmd = [ 'cmake', '--build', str(build_dir), '--preset', preset, '--target', 'package' ]
    rcode, sout, eout = proc_run(cmd, capture=True)
    if rcode:
        return rcode,sout, eout

    # test  CMAKE_INSTALL_PREFIX/bin/sortmerna --version
    exe = 'sortmerna.exe' if IS_WIN else 'sortmerna'
    exef = Path(install_dir) / 'bin' / exe
    if not exef.exists():
        msg = f'{ST} file {exef} not found'
        print(msg)
        return 1,'',msg
    else:
        cmd = [ str(exef), '--version' ]
        rcode, sout, eout = proc_run(cmd, capture=True)
        if not rcode:
            cmd = [ str(exef), '-h' ]
            rcode, sout, eout = proc_run(cmd, capture=True)
    return rcode, sout, eout
#END smr_build

def concurrentqueue_build(**kw):
    '''
    a single header file - just clone and use
    '''
    url = kw.get(CCQUEUE).get('url')
    src = kw.get(CCQUEUE).get('src')
    commit = kw.get(CCQUEUE).get('commit')
    shallow = kw.get(CCQUEUE).get('shallow')
    rcode, sout, eout = git_clone(url, src, shallow=shallow)
    return rcode, sout, eout
#END concurrentqueue_build

def env_check(**cfg):
    '''
    verify all pre-requisit packages/tools are present

    args:
      - cfg
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

def get_sysroot():
    '''
    check if the compiler was configured with_sysroot
    '''
    ST = '[get_sysroot]'
    sysroot = None
    gcc = shutil.which('gcc')
    if gcc:
        print(f'{ST} using {gcc}')
        cmd = ['gcc', '-print-sysroot']
        rcode, sysroot, eout = proc_run(cmd, capture=True)
        if rcode:
            return None, eout
        print(f'{ST} sysroot: {sysroot}')
        if not sysroot:
            eout = f'{ST} sysroot is not defined'
            print(eout)
        return sysroot, eout
    else:
        eout = f'{ST} no gcc found'
        print(eout)
        return sysroot, eout

modfunc['git_install'] = git_install
modfunc['conda_install'] = conda_install
modfunc['gcc_install'] = gcc_install
modfunc['cmake_install'] = cmake_install

if __name__ == "__main__":
    '''
    python setup.py -n all         build zlib + rocksdb + smr
    python setup.py -n sortmerna   build smr
    python setup.py -n sortmerna [--cmake-preset <preset>]  TODO
    python setup.py -n zlib        build zlib
    python setup.py -n rocksdb     build rocksdb

    python setup.py -n conda       install conda  TODO
    python setup.py -n cmake       install cmake  TODO
    '''
    ST = '[setup.py:__main__]'
    is_opts_ok = True

    # options
    optpar = OptionParser()
    optpar.add_option('-n', '--name', dest='name', 
        help='Package to build/process e.g. sortmerna | zlib | rocksdb | rapidjson | conda | all')
    optpar.add_option('--toolchain-config', dest='toolchain', help='CMake toolchain file')
    optpar.add_option('--zlib-src', dest='zlib_src', help='ZLib source directory')
    optpar.add_option('--zlib-build', dest='zlib_build', help='ZLib cmake build directory')
    optpar.add_option('--zlib-dist', dest='zlib_dist', help='ZLib installation directory')
    optpar.add_option('--zlib-git', action="store_true", help='use zlib git repo as source. Otherwise tarball')
    optpar.add_option('--rocksdb-dist', dest='rocksdb_dist', help='RocksDB installation directory')
    optpar.add_option('--rocksdb-src', dest='rocksdb_src', help='RocksDB sources directory')
    optpar.add_option('--rocksdb-build', dest='rocksdb_build', help='RocksDB cmake build directory')
    optpar.add_option('--rocksdb-git', action="store_true", help='use rocksdb git repo as source. Otherwise tarball')
    optpar.add_option('--build-dev', action="store_true", help='run build in development mode using git repos')
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
    optpar.add_option('--use-conda-cpp', action="store_false", help='Use Conda C++ tools for building')
    (opts, args) = optpar.parse_args()

    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print(f'{ST} Current dir: {cur_dir}')

    print(f'{ST} loading configuration from 3rdparty.jinja')
    envjn = Environment(loader=FileSystemLoader('.'), trim_blocks=True, lstrip_blocks=True)
    thirdparty_jinja = envjn.get_template('3rdparty.jinja')
    vars = {}
    jinja_str = thirdparty_jinja.render(vars)
    config = yaml.load(jinja_str, Loader=yaml.FullLoader) # render jinja template
    
    if opts.toolchain:
        config['toolchain'] = opts.toolchain
       
    # zlib
    if opts.zlib_src:
        config['zlib']['src'] = opts.zlib_src  # sources location
        
    if opts.zlib_dist:
        config['zlib']['dist'] = opts.zlib_dist  # installation dir (distro)
    else:
        config['zlib']['dist'] = f"build/{config['zlib']['src']}/dist"
        
    if opts.zlib_git or opts.build_dev:
        config['zlib']['is_git'] = True  # indicate if the sources are a git repo

    # rocksdb
    if opts.rocksdb_src:
        config['rocksdb']['src'] = opts.rocksdb_src
        
    if opts.rocksdb_dist:
        config['rocksdb']['dist'] = opts.rocksdb_dist
    else:
        config['rocksdb']['dist'] = f"build/{config['rocksdb']['src']}/dist"
        
    if opts.rocksdb_git or opts.build_dev:
        config['rocksdb']['is_git'] = True
        
    # concurrentqueue
    if opts.concurrentqueue_dist:
        config['concurrentqueue']['dist'] = opts.concurrentqueue_dist
    else:
        config['concurrentqueue']['dist'] = f"{config['concurrentqueue']['src']}"
        
    if opts.vb:
        config['vb'] = True
    if opts.loglevel:
        config['loglevel'] = opts.loglevel
    elif opts.trace:
        config['trace'] = True
    if opts.use_conda_cpp:
        config['conda_cpp'] = True

    opts.name = opts.name or ALL
    if opts.name in ALL:
        if opts.clean:
            rcode, sout, eout = clean('build')
        if not opts.cmake_preset:
            opts.cmake_preset = 'WIN_release' if IS_WIN else 'LIN_release'

        rcode, outl, errl = zlib_build(**config)
        if rcode == 0:
            rcode, sout, eout = rocksdb_build(**config) # ROCKS_VER, cfg=env
        if rcode == 0:
            ret = concurrentqueue_build(**config)
        if rcode == 0:
            btype = opts.btype or 'release'
            rcode, sout, eout = smr_build(btype=btype, **config)
    elif opts.name == ZLIB:
        kw = config.get(ZLIB,{})
        rcode, outl, errl = zlib_build(**config)
    elif opts.name == ROCKS:
        if opts.clean:
            ...
        rocksdb_build(**config)
    elif opts.name in [SMR]:
        if opts.clean:
            ...
        btype = opts.btype or 'release'
        smr_build(btype=btype, **config)
    elif opts.name == CCQUEUE: 
        concurrentqueue_build(**config)
    elif opts.name == DIRENT: 
        url = config[DIRENT].get('url')
        path = config[DIRENT].get('src')
        commit = config[DIRENT].get('commit')
        shallow = config[DIRENT].get('shallow')
        git_clone(url, path, shallow=shallow) 
    elif opts.name == CMAKE: 
        if opts.clone:
            ... #git_clone(env[CMAKE]['url'], LIB_DIR)
        ... #cmake_install(**env) 
    elif opts.name == CONDA: 
        ... #conda_install(**env) 
    #elif opts.name == RAPID: 
        #rapidjson_build()
    else: test()
#END main

#END END
