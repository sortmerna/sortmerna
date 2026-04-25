# @copyright 2016-2026 Clarity Genomics BVBA
# @copyright 2012-2016 Bonsai Bioinformatics Research Group
# @copyright 2014-2016 Knight Lab, Department of Pediatrics, UCSD, La Jolla
#
# SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
#
# This is a free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SortMeRNA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
#
# @contributors Jenya Kopylova   jenya.kopylov@gmail.com
#               Laurent Noé      laurent.noe@lifl.fr
#               Pierre Pericard  pierre.pericard@lifl.fr
#               Daniel McDonald  wasade@gmail.com
#               Mikaël Salson    mikael.salson@lifl.fr
#               Hélène Touzet    helene.touzet@lifl.fr
#               Rob Knight       robknight@ucsd.edu
#               biocodz          biocodz@protonmail.com

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
from argparse import ArgumentParser
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
IBZIP2  = 'indexed_bzip2'

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
              force:bool=False) -> tuple[int, list[str], list[str]]:
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

def zlib_build(**kw) -> tuple[int, list[str], list[str]]:
    '''
    args:
      - url
      - path str              local filesystem directory where to clone the repo
      - commit str            git commit
      - cmake_build_type str  Build type Relase | Debug | ..
      - ptype str             Linkage type like statuc, dynamic, mixed
      - is_git bool           use git as source. Otherwise - archive (.tar.gz)
    '''
    rcode, outl, errl = 0, [], []
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

def rocksdb_modify_3party_zlib(link_type:str='t1', **kw):
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
    link_type = kw.get(ROCKS,{}).get('link.type', {}).get('windows',{}).get(link_type,{})
    zlib_rel = link_type.get('zlib.release')  # cfg[ROCKS]['link']['WIN'][ptype]['ZLIB_LIB_RELEASE']
    zlib_dbg = link_type.get('zlib.debug')  # cfg[ROCKS]['link']['WIN'][ptype]['ZLIB_LIB_DEBUG']
    file3p = Path(kw.get(ROCKS).get('src')) / 'thirdparty.inc'
    zlib_dist = Path(kw.get('zlib',{}).get('dist')).absolute()  # absolute to avoid package finding problems

    for line in fileinput.FileInput(file3p, inplace=True):
        if line.startswith('set(ZLIB_HOME'):
            line = re.sub(r'ZLIB_HOME .*\)', r'ZLIB_HOME {})'.format(zlib_dist.as_posix()), line, flags = re.M)
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
    
def rocksdb_build(link_type:str='t1', **kw) -> tuple[int, list[str], list[str]]: # ver=None, btype='Release', ptype='t3', **cfg
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
        f'-DCMAKE_INSTALL_PREFIX:PATH={Path(dist_dir).absolute().as_posix()}',
        f'-DCMAKE_POLICY_DEFAULT_CMP0074=NEW'
    ]
    # 20260210 Tue  warning: Policy CMP194 is not set: MSVC is not an assembler for language ASM.
    if IS_WIN:
        cmd.append(f'-DCMAKE_POLICY_DEFAULT_CMP194=OLD')
    else:
        # on Win 'find_package' is not used to search for zlib, so ZLIB_ROOT is not used. 
        # Instead 'include_directories' and THIRDPARTY_LIBS variables are used to direct
        # the build to ZLIB's includes and libs. See 'thirdparty.inc'
        cmd.append(f'-DZLIB_ROOT:PATH={Path(zlib_dist).as_posix()}')
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
              **kw) -> tuple[int, list[str], list[str]]:
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

def concurrentqueue_build(**kw) -> tuple[int, list[str], list[str]]:
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

def indexed_bzip2_build(**kw) -> tuple[int, list[str], list[str]]:
    '''
    header-only rapidgzip library - clone and init submodules (external/zlib etc.)
    '''
    url = kw.get(IBZIP2).get('url')
    src = kw.get(IBZIP2).get('src')
    shallow = kw.get(IBZIP2).get('shallow')
    rcode, sout, eout = git_clone(url, src, shallow=shallow)
    if rcode == 0:
        cmd = ['git', 'submodule', 'update', '--init', '--recursive']
        rcode, sout, eout = proc_run(cmd, src, capture=True)
    return rcode, sout, eout
#END indexed_bzip2_build

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
    python setup.py all          build zlib + rocksdb + smr
    python setup.py sortmerna    build smr
    python setup.py sortmerna [--cmake-preset <preset>]
    python setup.py zlib         build zlib
    python setup.py rocksdb      build rocksdb
    python setup.py conda        install conda  TODO
    python setup.py cmake        install cmake  TODO
    '''
    ST = '[setup.py:__main__]'

    p0 = ArgumentParser(description='SortMeRNA build utility')
    p0.add_argument('--toolchain-config', dest='toolchain', help='CMake toolchain file')
    p0.add_argument('--vb', action='store_true', help='Export compile commands')
    p0.add_argument('--loglevel', dest='loglevel', help='CMake log level')
    p0.add_argument('--trace', action='store_true', help='Run cmake with --trace')
    p0.add_argument('--use-conda-cpp', action='store_true', dest='use_conda_cpp', help='Use Conda C++ tools for building')
    p0.add_argument('-e', '--envn', dest='envname', help='Name of environment: WIN | WSL | LIN ..')
    p0.add_argument('--env', dest='envfile', help='Env configuration file.')
    p0.add_argument('--config', dest='config', help='Build configuration file.')
    p0.add_argument('--build-dir', dest='build_dir', help='Build directory.')
    p0.add_argument('--dist-dir', dest='dist_dir', help='Distro directory.')
    p0.add_argument('--local-linux', dest='local_linux', action='store_true', help='Perform the build on local source files.')
    p0.add_argument('--winhome', dest='winhome', help='when building on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    subpar = p0.add_subparsers(dest='name', help='package to build/process')

    p_all = subpar.add_parser('all', help='build zlib + rocksdb + concurrentqueue + indexed_bzip2 + sortmerna')
    p_all.add_argument('-b', '--btype', dest='btype', default='release', help='Build type: release | debug')
    p_all.add_argument('--cmake-preset', dest='cmake_preset', help='CMake preset')
    p_all.add_argument('-c', '--clean', action='store_true', help='clean build directory')
    p_all.add_argument('--build-dev', action='store_true', help='run build in development mode using git repos')
    p_all.add_argument('--zlib-src', dest='zlib_src', help='ZLib source directory')
    p_all.add_argument('--zlib-build', dest='zlib_build', help='ZLib cmake build directory')
    p_all.add_argument('--zlib-dist', dest='zlib_dist', help='ZLib installation directory')
    p_all.add_argument('--zlib-git', action='store_true', help='use zlib git repo as source. Otherwise tarball')
    p_all.add_argument('--rocksdb-src', dest='rocksdb_src', help='RocksDB sources directory')
    p_all.add_argument('--rocksdb-build', dest='rocksdb_build', help='RocksDB cmake build directory')
    p_all.add_argument('--rocksdb-dist', dest='rocksdb_dist', help='RocksDB installation directory')
    p_all.add_argument('--rocksdb-git', action='store_true', help='use rocksdb git repo as source. Otherwise tarball')
    p_all.add_argument('--concurrentqueue-dist', dest='concurrentqueue_dist', help='concurrentqueue installation directory')
    p_all.add_argument('--indexed-bzip2-dist', dest='indexed_bzip2_dist', help='indexed_bzip2 installation directory')
    p_all.add_argument('--pt_smr', dest='pt_smr', default='t1', help='Sortmerna Linkage type t1 | t2 | t3')
    p_all.add_argument('--pt_zlib', dest='pt_zlib', help='Zlib Linkage type t1 | t2 | t3')
    p_all.add_argument('--pt_rocks', dest='pt_rocks', help='Rocksdb Linkage type t1 | t2 | t3')

    p_smr = subpar.add_parser('sortmerna', help='build sortmerna')
    p_smr.add_argument('-b', '--btype', dest='btype', default='release', help='Build type: release | debug')
    p_smr.add_argument('--cmake-preset', dest='cmake_preset', help='CMake preset')
    p_smr.add_argument('-c', '--clean', action='store_true', help='clean build directory')
    p_smr.add_argument('--pt_smr', dest='pt_smr', default='t1', help='Sortmerna Linkage type t1 | t2 | t3')

    p_zlib = subpar.add_parser('zlib', help='build zlib')
    p_zlib.add_argument('--zlib-src', dest='zlib_src', help='ZLib source directory')
    p_zlib.add_argument('--zlib-build', dest='zlib_build', help='ZLib cmake build directory')
    p_zlib.add_argument('--zlib-dist', dest='zlib_dist', help='ZLib installation directory')
    p_zlib.add_argument('--zlib-git', action='store_true', help='use zlib git repo as source. Otherwise tarball')

    p_rocks = subpar.add_parser('rocksdb', help='build rocksdb')
    p_rocks.add_argument('--rocksdb-src', dest='rocksdb_src', help='RocksDB sources directory')
    p_rocks.add_argument('--rocksdb-build', dest='rocksdb_build', help='RocksDB cmake build directory')
    p_rocks.add_argument('--rocksdb-dist', dest='rocksdb_dist', help='RocksDB installation directory')
    p_rocks.add_argument('--rocksdb-git', action='store_true', help='use rocksdb git repo as source. Otherwise tarball')
    p_rocks.add_argument('--zlib-dist', dest='zlib_dist', help='ZLib installation directory')
    p_rocks.add_argument('-c', '--clean', action='store_true', help='clean build directory')
    p_rocks.add_argument('--pt_rocks', dest='pt_rocks', help='Rocksdb Linkage type t1 | t2 | t3')
    p_rocks.add_argument('--rocks3p', dest='rocks3p', help='Fix thirdparty.inc when building Rocksdb')

    p_ccq = subpar.add_parser('concurrentqueue', help='clone concurrentqueue')
    p_ccq.add_argument('--concurrentqueue-dist', dest='concurrentqueue_dist', help='concurrentqueue installation directory')

    p_ibzip2 = subpar.add_parser('indexed_bzip2', help='clone indexed_bzip2 and init submodules')
    p_ibzip2.add_argument('--indexed-bzip2-dist', dest='indexed_bzip2_dist', help='indexed_bzip2 installation directory')

    subpar.add_parser('dirent', help='clone dirent')

    p_cmake = subpar.add_parser('cmake', help='install cmake')
    p_cmake.add_argument('--clone', action='store_true', help='Perform git clone')

    subpar.add_parser('conda', help='install conda')

    args = p0.parse_args()

    if args.name is None:
        p0.print_help()
        sys.exit(0)

    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print(f'{ST} Current dir: {cur_dir}')

    print(f'{ST} loading configuration from 3rdparty.jinja')
    envjn = Environment(loader=FileSystemLoader('.'), trim_blocks=True, lstrip_blocks=True)
    thirdparty_jinja = envjn.get_template('3rdparty.jinja')
    vars = {}
    jinja_str = thirdparty_jinja.render(vars)
    config = yaml.load(jinja_str, Loader=yaml.FullLoader) # render jinja template

    if args.toolchain:
        config['toolchain'] = args.toolchain

    # zlib
    if getattr(args, 'zlib_src', None):
        config['zlib']['src'] = args.zlib_src
    if getattr(args, 'zlib_dist', None):
        config['zlib']['dist'] = args.zlib_dist
    else:
        config['zlib']['dist'] = f"build/{config['zlib']['src']}/dist"
    if getattr(args, 'zlib_git', False) or getattr(args, 'build_dev', False):
        config['zlib']['is_git'] = True

    # rocksdb
    if getattr(args, 'rocksdb_src', None):
        config['rocksdb']['src'] = args.rocksdb_src
    if getattr(args, 'rocksdb_dist', None):
        config['rocksdb']['dist'] = args.rocksdb_dist
    else:
        config['rocksdb']['dist'] = f"build/{config['rocksdb']['src']}/dist"
    if getattr(args, 'rocksdb_git', False) or getattr(args, 'build_dev', False):
        config['rocksdb']['is_git'] = True

    # concurrentqueue
    if getattr(args, 'concurrentqueue_dist', None):
        config['concurrentqueue']['dist'] = args.concurrentqueue_dist
    else:
        config['concurrentqueue']['dist'] = f"{config['concurrentqueue']['src']}"

    # indexed_bzip2
    if getattr(args, 'indexed_bzip2_dist', None):
        config['indexed_bzip2']['dist'] = args.indexed_bzip2_dist
    else:
        config['indexed_bzip2']['dist'] = f"{config['indexed_bzip2']['src']}"

    if args.vb:
        config['vb'] = True
    if args.loglevel:
        config['loglevel'] = args.loglevel
    elif args.trace:
        config['trace'] = True
    if args.use_conda_cpp:
        config['conda_cpp'] = True

    if args.name == ALL:
        if getattr(args, 'clean', False):
            rcode, sout, eout = clean('build')
        if not getattr(args, 'cmake_preset', None):
            args.cmake_preset = 'WIN_release' if IS_WIN else 'LIN_release'
        rcode, outl, errl = zlib_build(**config)
        if rcode == 0:
            rcode, sout, eout = rocksdb_build(**config)
        if rcode == 0:
            rcode, sout, eout = concurrentqueue_build(**config)
        if rcode == 0:
            rcode, sout, eout = indexed_bzip2_build(**config)
        if rcode == 0:
            rcode, sout, eout = smr_build(btype=getattr(args, 'btype', 'release'), **config)
    elif args.name == ZLIB:
        rcode, outl, errl = zlib_build(**config)
    elif args.name == ROCKS:
        if getattr(args, 'clean', False):
            ...
        rocksdb_build(**config)
    elif args.name == SMR:
        if getattr(args, 'clean', False):
            ...
        smr_build(btype=getattr(args, 'btype', 'release'), **config)
    elif args.name == CCQUEUE:
        concurrentqueue_build(**config)
    elif args.name == IBZIP2:
        indexed_bzip2_build(**config)
    elif args.name == DIRENT:
        url = config[DIRENT].get('url')
        path = config[DIRENT].get('src')
        shallow = config[DIRENT].get('shallow')
        git_clone(url, path, shallow=shallow)
    elif args.name == CMAKE:
        if getattr(args, 'clone', False):
            ... #git_clone(env[CMAKE]['url'], LIB_DIR)
        ... #cmake_install(**env)
    elif args.name == CONDA:
        ... #conda_install(**env)
    else:
        test()
#END main

#END END
