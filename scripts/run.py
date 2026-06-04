# @copyright 2016-2026 Clarity Genomics Inc
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
#

'''
file: run.py
created: Aug 12, 2019 Mon

conda install scikit-bio -c conda-forge  <- pre-requisites
'''

def mock_missing(name):
    def init(self, *args, **kwargs):
        raise ImportError(f'Failed to import class {name}; likely not installed.')
    return type(name, (), {'__init__': init})

import os
import signal
import sys
import struct
import ctypes
import subprocess
import platform
import re
import time
import difflib
import shutil
import yaml
import gzip
import multiprocessing
from argparse import ArgumentParser
from argparse import Namespace
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

try:
    import pandas
except ImportError:
    pandas = mock_missing('pandas')
    
try:
    import rapidgzip
except ImportError:
    rapidgzip = mock_missing('rapidgzip')

is_skbio = True
try:
    import skbio.io
except (ImportError, OSError) as ex:
    print(f'\'skbio.io\' not available - skipping: {ex}')
    is_skbio = False

# globals
OS = None
ENV = None # WIN | WSL | LNX_AWS | LNX_TRAVIS
WRK_DIR = None

# define platform
pf = platform.platform()
IS_WIN = 'Windows' in pf
IS_WSL = 'Linux' in pf and 'Microsoft' in pf # Windows Subsystem for Linux (WSL)
IS_LNX = 'Linux' in pf and not 'Microsoft' in pf
OS = 'WIN' if IS_WIN else 'WSL' if IS_WSL else 'LNX' if IS_LNX else None
if not OS:
    print(f'Unexpected platform: {pf}')
    sys.exit(1)

UHOME = os.environ.get('USERPROFILE') if IS_WIN else os.environ.get('HOME') # not defined for AWS SSM

SMR = 'sortmerna'
SMR_SRC = None # source root dir
SMR_DIST = None # dist root dir
SMR_EXE = None # full path to executable

ZLIB_SRC = None
ZLIB_DIST = None

ROCKS_SRC = None
ROCKS_DIST = None

# no binaries, so always build Release only
RAPID_SRC = None
RAPID_DIST = None

DATA_DIR = None
RUN_DIR = None
OUT_DIR = None
TEST_DATA = None

# base names of the report files
ALI_NAMES     = []  # get from test.jinja.yaml + fa/fa + gz/non-gz -> 7 x 4 = 28 names
OTH_NAMES     = []  
ALI_BASE      = 'aligned'
OTH_BASE      = 'other'
ALI_FWD_BASE  = None # '{}_fwd'.format(ALI_BASE)
ALI_REV_BASE  = None # '{}_rev'.format(ALI_BASE)
OTH_FWD_BASE  = None # '{}_fwd'.format(OTH_BASE)
OTH_REV_BASE  = None # '{}_rev'.format(OTH_BASE)
LOG_BASE      = None # '{}.log'.format(ALI_BASE)
DENOVO_BASE   = None # '{}_denovo.fasta'.format(ALI_BASE)
OTU_BASE      = None # '{}_otus.txt'.format(ALI_BASE)
BLAST_BASE    = None # '{}.blast'.format(ALI_BASE)
SAM_BASE      = None # '{}.sam'.format(ALI_BASE)
READS_EXT     = None
IS_FASTQ      = False
IS_PAIRED_IN  = False
IS_PAIRED_OUT = False

# output files
LOGF    = None
ALIF    = None
ALI_FWD = None
ALI_REV = None
OTHF    = None
OTH_FWD = None
OTH_REV = None
ALI_REF = None
DENOVOF = None
OTUF    = None
BLASTF  = None
SAMF    = None

# KVDB and index dirs
KVDB_DIR = None
IDX_DIR  = None

def is_windows():
    return sys.platform.startswith("win")

def is_linux():
    return sys.platform.startswith("linux")

def is_darwin():
    return sys.platform.startswith("darwin")

def run_test(cmd:list, cwd:str=None, capture:bool=False) -> tuple[int, list[str], list[str]]:
    '''
    run a test
    args:
      - cmd  command to run
    '''
    ST = '[run.run_test]'
    rcode, outl, errl = 0, [], []
    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe
    cmds = ' '.join(cmd)
    msg = f'{ST} Running: {cmds} in {cwd}' if cwd else f'{ST} Running: {cmds}'
    print(msg)

    start = time.time()
    #print('Running {} in {}'.format(' '.join('{}'.format(xx) for xx in cmd), cwd))
    # cmake configure and generate
    #kw = {'stdout':subprocess.PIPE, 'stderr':subprocess.PIPE} if capture else {}
    try:
        if cwd and not Path(cwd).exists():
            os.makedirs(cwd)
            proc = subprocess.run(cmd, cwd=cwd, capture_output=capture)
        else:
            proc = subprocess.run(cmd, capture_output=capture)

        rcode = proc.returncode
        if capture:
            outl.append(proc.stdout.decode().strip())
            errl.append(proc.stderr.decode().strip())
        #proc = subprocess.run(cmd, cwd=build_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError as err:
        print(err)
        rcode = 1
        errl.append(str(err))
    except Exception as ex:
        print(str(ex))
        rcode = 1
        errl.append(str(ex))

    rt = time.time() - start
    print(f"{ST} run time: {rt}")
    return rcode, outl, errl
#END cmake_run

def run_with_interrupts(cmd:list,
                        num_interrupts:int,
                        wait_after_align_start:int=60,
                        startup_deadline:int=1800,
                        cwd:str=None) -> tuple[int, list[str], list[str]]:
    '''
    Exercise the auto-resume path. Launches sortmerna num_interrupts times in
    succession; on each non-final launch we wait until "Processor 0 thread ..
    started" appears, sleep wait_after_align_start seconds, then SIGKILL the
    process group. The final launch runs to completion and relies on
    restart::probe_or_init to detect the prior mid-pass interrupt and roll back.

    args:
      cmd                     sortmerna argv (without stdbuf prefix)
      num_interrupts          number of mid-pass kills before the final run
      wait_after_align_start  seconds to let alignment progress before SIGKILL
      startup_deadline        seconds to wait for alignment to start before giving up
    '''
    ST = '[run_with_interrupts]'
    print(f'{ST} num_interrupts={num_interrupts} wait_after_align_start={wait_after_align_start}s')

    # Line-buffer the child's stdout/stderr so log polling sees lines promptly.
    stdbuf_path = shutil.which('stdbuf')
    if stdbuf_path:
        wrapped = [stdbuf_path, '-oL', '-eL'] + cmd
    else:
        print(f'{ST} WARN: stdbuf not found; child stdout may be block-buffered and the log poll may lag')
        wrapped = cmd

    log_path = Path(f'/tmp/sortmerna_interrupt_{os.getpid()}.log')
    marker = re.compile(r'Processor 0 thread .* started')
    done_pat = re.compile(r'Processor .* done\. Processed')

    def echo_log(text:str, pos:int) -> int:
        '''echo any text beyond pos to this process' stdout, return new pos'''
        if len(text) > pos:
            sys.stdout.write(text[pos:])
            sys.stdout.flush()
        return len(text)

    total_attempts = num_interrupts + 1
    for iattempt in range(total_attempts):
        is_final = (iattempt == num_interrupts)
        label = 'final' if is_final else f'kill-mid-pass {iattempt+1}/{num_interrupts}'
        log_path.unlink(missing_ok=True)
        log_path.touch()

        print(f'{ST} attempt {iattempt+1}/{total_attempts} ({label})')
        with open(log_path, 'w') as flog:
            proc = subprocess.Popen(wrapped, stdout=flog, stderr=subprocess.STDOUT,
                                    cwd=cwd, start_new_session=True)
        print(f'{ST} pid={proc.pid}')

        if is_final:
            log_pos = 0
            while proc.poll() is None:
                log_pos = echo_log(log_path.read_text(errors='replace'), log_pos)
                time.sleep(0.5)
            # flush any tail written between the last poll and exit
            echo_log(log_path.read_text(errors='replace'), log_pos)
            rcode = proc.returncode
            print(f'{ST} sortmerna exited rcode={rcode}')
            return rcode, [], []

        deadline = time.time() + startup_deadline
        log_pos = 0
        while True:
            text = log_path.read_text(errors='replace')
            log_pos = echo_log(text, log_pos)
            if proc.poll() is not None:
                print(f'{ST} sortmerna exited before alignment started (rc={proc.returncode})')
                return 1, [], []
            if time.time() > deadline:
                _kill_pg(proc)
                print(f'{ST} timeout waiting for alignment to start')
                return 1, [], []
            if any(marker.search(line) for line in text.splitlines()):
                break
            time.sleep(2)
        print(f'{ST} alignment threads up; waiting {wait_after_align_start}s')

        wait_end = time.time() + wait_after_align_start
        while time.time() < wait_end and proc.poll() is None:
            log_pos = echo_log(log_path.read_text(errors='replace'), log_pos)
            time.sleep(0.5)
        log_pos = echo_log(log_path.read_text(errors='replace'), log_pos)

        done_count = sum(1 for line in log_path.read_text(errors='replace').splitlines()
                         if done_pat.search(line))
        if done_count > 0:
            print(f'{ST} WARN: {done_count} thread(s) already finished; interrupt landing past mid-pass. '
                  f'Reduce wait_after_align_start or use a larger dataset.')

        print(f'{ST} sending SIGKILL to pgid={proc.pid}')
        _kill_pg(proc)
        try:
            proc.wait(timeout=15)
        except subprocess.TimeoutExpired:
            print(f'{ST} ERROR: sortmerna did not exit within 15s of SIGKILL')
            return 1, [], []
        print(f'{ST} killed cleanly; next launch should auto-resume')

    return 0, [], []
#END run_with_interrupts

def _kill_pg(proc):
    '''Send SIGKILL to the process group started via start_new_session=True.'''
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, OSError):
        pass
#END _kill_pg

def get_diff(fpath0, fpath1):
    '''
    '''
    dlist = []

    with open(fpath0, 'rU') as fout:
        with open(fpath1, 'rU') as fexpect:
            diff = difflib.unified_diff(
                fout.readlines(),
                fexpect.readlines(),
                fromfile='fout',
                tofile='fexpect',
            )
            for line in diff:
                dlist.append(line)
                if line: 
                    print(line)
    return dlist
#END get_diff

def to_lf(ddir):
    '''
    convert to LF line endings of the data files

    NOTE: 'find . -type f -name "*.fasta" -o -name "*.fastq" | xargs dos2unix'
          using 'subprocess.run' is problematic because of the pipe. 
    '''
    ST = '[to_lf]'
    CRLF = b'\r\n'
    LF = b'\n'

    for dpath, dnames, fnames in os.walk(ddir):
        for fname in fnames:
            if fname.endswith('.fasta') or fname.endswith('.fastq'):
                fpath = os.path.join(dpath, fname)
                print(f'{ST} converting {fpath}')
                with open(fpath, 'rb') as infile:
                    content = infile.read()
                content = content.replace(CRLF, LF)
                with open(fpath, 'wb') as infile:
                    infile.write(content)
    print('to_lf Done')
#END to_lf

def parse_log(fpath:str):
    '''
    parse 'aligned.log' to dictionary
    args"
      - fpath  'aligned.log' file
    '''
    logd = {}
    if not Path(fpath).exists():
        print(f'{fpath} does not exist')
        return logd
    with open(fpath) as f_log:    
        for line in f_log:
            if 'Total reads =' in line:
                logd['num_reads'] = int((re.split(r' = ', line)[1]).strip())
            elif 'Total reads for de novo clustering'in line:
                logd['num_denovo'] = int((re.split(r' = ', line)[1]).strip())
            elif 'Total reads passing E-value threshold' in line:
                logd['num_hits'] = int((re.split(r' = | \(', line)[1]).strip())
            elif 'Total reads failing E-value threshold' in line:
                logd['num_fail'] = int((re.split(r' = | \(', line)[1]).strip())
            elif 'Total reads passing %%id and %%coverage thresholds' in line:
                # line: '..thresholds = 44223 (44.22)'
                val = line.split('=')[1]
                if val: val = val.split()[0].strip()
                logd['num_id_cov'] = int(val)
            elif 'Total OTUs' in line:
                logd['num_otus'] = int((re.split(' = ', line)[1]).strip())
    return logd
#END parse_log

def process_smr_opts(args:list):
    '''
    args:
      - args  list of parameters passed to sortmerna
    '''
    ST = '[process_smr_opts]'
    WDIR = '-workdir'
    KVD = '-kvdb'
    IDX = '-idx'
    ALN = '-aligned'
    OTH = '-other'
    OUT2 = '-out2'

    global KVDB_DIR
    global IDX_DIR
    global ALI_BASE
    global OTH_BASE

    global LOGF
    global ALIF
    global ALI_FWD
    global ALI_REV
    global OTHF
    global OTH_FWD
    global OTH_REV
    global ALI_REF
    global DENOVOF
    global OTUF
    global BLASTF
    global SAMF
    global READS_EXT
    global IS_FASTQ
    global IS_PAIRED_IN 
    global IS_PAIRED_OUT
    global WRK_DIR
    
    is_gz = False
    psplit = os.path.splitext(args[args.index('-reads')+1]) # 1.concat.fq.gz -> [1.concat.fq, .gz]
    READS_EXT = psplit[1]
    if READS_EXT in ['.gz']:
        #READS_EXT = os.extsep + args[args.index('-reads')+1].split(os.extsep)[1]
        READS_EXT = os.path.splitext(psplit[0])[1]
        is_gz = True
    IS_FASTQ =  READS_EXT[1:] in ['fq', 'fastq']
    IS_PAIRED_IN = '-paired_in' in args
    IS_PAIRED_OUT = '-paired_out' in args

    if ALN in args:
        aln_pfx = args[args.index(ALN) + 1]
        if not os.path.basename(aln_pfx):
            ALIF = os.path.join(aln_pfx, ALI_BASE + READS_EXT) # use default 'aligned'
        else:
            ALI_BASE = os.path.basename(aln_pfx)
            ALIF = os.path.abspath(aln_pfx + READS_EXT)
    elif WDIR in args:
        wdir = args[args.index(WDIR) + 1]
        print(f"{ST} '-workdir' option was provided. Using workdir: [{os.path.realpath(wdir)}]")
        ALIF = os.path.join(wdir, 'out', ALI_BASE + READS_EXT)
    elif WRK_DIR:
        ALIF = os.path.join(WRK_DIR, 'out', ALI_BASE + READS_EXT)
    elif UHOME:
        ALIF = os.path.join(UHOME, 'sortmerna', 'run', 'out', ALI_BASE + READS_EXT)
    else:
        print(f'{ST} cannot define alignment file')

    if OUT2 in args:
        ALI_FWD = os.path.join(os.path.dirname(ALIF), ALI_BASE + '_fwd' + READS_EXT)
        ALI_REV = os.path.join(os.path.dirname(ALIF), ALI_BASE + '_rev' + READS_EXT)

    if OTH in args:
        idx = args.index(OTH)
        # check idx + 1 no exceeds args length and other options has arg
        if idx + 1 < len(args) -1 and args[idx+1][:1] != '-':
            oth_pfx = args[idx + 1]
            if not os.path.basename(oth_pfx):
                OTHF = os.path.join(oth_pfx, OTH_BASE + READS_EXT)
            else:
                OTH_BASE = os.path.basename(oth_pfx)
                OTHF = os.path.abspath(oth_pfx + READS_EXT)
        elif ALN in args:
            OTHF = os.path.join(os.path.dirname(ALIF), OTH_BASE + READS_EXT) # use the same out dir as ALN
        elif WDIR in args:
            wdir = args[args.index(WDIR) + 1]
            print(f'{ST} \'-workdir\' option was provided. Using workdir: [{os.path.realpath(wdir)}]')
            ALIF = Path(wdir) / 'out' / f'{OTH_BASE}.{READS_EXT}'
        else:
            OTHF = Path(UHOME) / 'sortmerna/run/out' / f'{OTH_BASE}.{READS_EXT}'

        if OUT2 in args:
            OTH_FWD = os.path.join(os.path.dirname(ALIF), OTH_BASE + '_fwd' + READS_EXT)
            OTH_REV = os.path.join(os.path.dirname(ALIF), OTH_BASE + '_rev' + READS_EXT)

    if KVD in args:
        KVDB_DIR = args[args.index(KVD) + 1]
    elif WDIR in args:
        KVDB_DIR = Path(args[args.index(WDIR) + 1]) / 'kvdb'
    else:
        KVDB_DIR = Path(UHOME) / 'sortmerna/run/kvdb'

    if IDX in args:
        IDX_DIR = args[args.index(IDX) + 1]
    elif WDIR in args:
        IDX_DIR = Path(args[args.index(WDIR) + 1]) / 'idx'
    else:
        IDX_DIR = Path(UHOME) / 'sortmerna/run/idx'

    gzs = '.gz' if is_gz else ''
    LOGF    = os.path.join(os.path.dirname(ALIF), f'{ALI_BASE}.log')
    BLASTF  = os.path.join(os.path.dirname(ALIF), f'{ALI_BASE}.blast{gzs}')
    OTUF    = os.path.join(os.path.dirname(ALIF), 'otu_map.txt')
    DENOVOF = os.path.join(os.path.dirname(ALIF), f'{ALI_BASE}_denovo.fa')
    SAMF    = os.path.join(os.path.dirname(ALIF), f'{ALI_BASE}.sam')
#END process_smr_opts

def process_blast(**kw):
    '''
    Check count of reads passing %id and %coverage threshold
    as given in aligned.blast
    '''
    ST = f'[process_blast]'
    vald = kw.get('validate')

    outdir = Path(get_dict_val('args:-workdir', kw))/'out'
    #logf = outdir /'aligned.log'
    #readf = Path(get_dict_val('args:-reads', kw)[0]).suffix
    gzx = '.gz' if Path(get_dict_val('args:-reads', kw)[0]).suffix == '.gz' else ''
    is_gz = bool(gzx)
    blastf = outdir / f'aligned.blast{gzx}'
    blast_pid_cov = outdir / 'pid_pcov.blast'
    BLAST_ID_COL = 2
    BLAST_COV_COL = 13
    num_hits_file = 0
    n_yid_ycov = 0
    n_yid_ncov = 0
    n_nid_ycov = 0
    n_denovo = 0
    has_cov = False
    is_dbg_pid_pcov = True

    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq', 'qcovs']
    
    if blastf.exists(): 
        print(f'{ST} processing : {blastf}')
        is_use_skbio = False
        with open(blastf, 'rb') as f_blast:
            is_gz = is_gz and (f_blast.read(2) == b'\x1f\x8b')
        if is_use_skbio:
            # just returns column names - not at all what's expected
            for seq in skbio.io.read(blastf, format='blast+6', columns=cols, 
                                     compression='gzip', into=pandas.DataFrame):
                num_hits_file += 1
        else:
            with gzip.open(blastf, 'rb') if is_gz else open(blastf, 'rb') as f_blast, open(blast_pid_cov, 'w') as f_pid_pcov:
                for lineb in f_blast:
                    num_hits_file += 1
                    line = lineb.decode('utf-8')
                    llist = line.strip().split('\t')
                    fid = float(llist[BLAST_ID_COL])
                    is_pass_id = fid >= 97.0
                    has_cov = len(llist) > BLAST_COV_COL
                    if has_cov:
                        fcov = float(llist[BLAST_COV_COL])
                        is_pass_cov = fcov >= 97.0
                        if is_pass_id:
                            if is_pass_cov: 
                                n_yid_ycov += 1
                                if is_dbg_pid_pcov:
                                    f_pid_pcov.write(line)
                            else:
                                n_yid_ncov += 1
                        elif is_pass_cov:
                            n_nid_ycov += 1
                        else:
                            n_denovo += 1
                    is_pass_id = False
                    is_pass_cov = False
    
            print(f'{ST} from {blastf.name}: num_hits= {num_hits_file} n_yid_ycov= {n_yid_ycov}'
                f' n_yid_ncov= {n_yid_ncov} n_nid_ycov= {n_nid_ycov} n_denovo= {n_denovo}')
            
            blastd = vald.get('files', {}).get('aligned.blast')
            if blastd:
                if blastd.get('n_yid_ycov'):
                    msg = (f"{ST} Testing reads passing ID threshold: {blastf.name}:"
                           f" {n_yid_ycov} Expected: {blastd['n_yid_ycov']}")
                    print(msg)
                    if n_yid_ycov != blastd['n_yid_ycov']:
                        print(f"Failed test: {blastd['n_yid_ycov']} not equals {n_yid_ycov}")

                num_recs = blastd.get('num_recs')
                if num_recs:
                    msg = (f'{ST} Testing num_hits: {blastf.name}:'
                            f' {num_hits_file} Expected: {num_recs}')
                    print(msg)
                    if num_hits_file != num_recs:
                        print(f'Failed test: {num_hits_file} not equals {num_recs}')
    return {
        'n_hits'  : num_hits_file, 
        'yid_ycov': n_yid_ycov, 
        'yid_ncov': n_yid_ncov, 
        'nid_ycov': n_nid_ycov,
        'n_denovo': n_denovo
        }
#END process_blast

def dbg_blast(**kwarg):
    '''
    added 20210323
    compare unique read IDs in two blast reports produced using different program versions
    cmd: python scripts/run.py --name t0 -f dbg_blast --validate-only
    '''
    ST = '[dbg_blast]'
    rdiff = None
    bl_421 = os.path.join(DATA_DIR, 'sortmerna/out/tests/t42/win10_8_16/4-2-1/20210322/aligned.blast')
    bl_431 = os.path.join(DATA_DIR, 'sortmerna/out/tests/t42/win10_8_16/4-3-1/20210322/aligned.blast')
    srt421 = os.path.join(DATA_DIR, 'sortmerna/out/tests/t42/win10_8_16/4-2-1/20210322/sorted.blast')
    srt431 = os.path.join(DATA_DIR, 'sortmerna/out/tests/t42/win10_8_16/4-3-1/20210322/sorted.blast')
    if os.path.exists(bl_421) and os.path.exists(bl_431):
        with open(bl_421) as f421, open(bl_431) as f431:
            l421 = [ line.strip().split('\t')[:2] for line in f421 ]
            l431 = [ line.strip().split('\t')[:2] for line in f431 ]
            l421.sort(key=lambda rr: int(rr[0].split('.')[-1]))  # SRR1635864.8745 -> [SRR1635864, 8745] -> 8745
            l431.sort(key=lambda rr: int(rr[0].split('.')[-1]))
            with open(srt421, 'w') as sr421, open(srt431, 'w') as sr431:
                sr421.write('\n'.join(f'{x[0]}\t{x[1]}' for x in l421))
                sr431.write('\n'.join(f'{x[0]}\t{x[1]}' for x in l431))
            #l421 = [ '  '.join(line.strip().split('\t')[:2]) for line in f421 ]
            #l431 = [ '  '.join(line.strip().split('\t')[:2]) for line in f431 ]
            doset = False
            if doset:
                rdiff = set(l421) - set(l431)
                print(f'{ST} rdiff.len= {len(rdiff)}') # 0
                [print(x) for x in list(rdiff)]
#END dbg_blast

def dbg_otu(**kwarg):
    '''
    '''
    ST = '[dbg_otu]'
    OTU_READSF = os.path.join(os.path.dirname(ALIF), 'otu_reads.txt')
    BLAST_PID_PCOV = os.path.join(os.path.dirname(ALIF), 'pid_pcov.blast')
    READS_DIFF = os.path.join(os.path.dirname(ALIF), 'reads_diff.txt')
    reads = []  # list of all reads from all groups
    blast_reads = []

    if os.path.exists(BLAST_PID_PCOV):
        with open(BLAST_PID_PCOV) as blastf:
            for line in blastf:
                lls = line.strip().split('\t')
                blast_reads.append(lls[0]) # collect reads

    # parse otu groups file
    if Path(OTUF).exists():
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        with open(OTUF) as f_otus:
            for line in f_otus:
                num_clusters_file += 1
                greads = line.strip().split('\t')[1:]  # reads in a group
                reads.extend(greads)
                num_reads_in_clusters_file += len(greads)

        # genreate list of reads
        with open(OTU_READSF, 'w') as readsf:
            reads.sort(key=lambda rr: int(rr.split('_')[-1]))
            reads.sort(key=lambda rr: int(rr.split('_')[0][:-1]))
            for read in reads:
                readsf.write(f'{read}\n')

        with open(READS_DIFF, 'w') as diff:
            rds = set(blast_reads) - set(reads)
            rdsl = list(rds)
            rdsl.sort(key=lambda rr: int(rr.split('_')[-1]))
            rdsl.sort(key=lambda rr: int(rr.split('_')[0][:-1]))
            for rd in rdsl:
                diff.write(f'{rd}\n')
#END dbg_otu

def validate_otu(**kw):
    '''
    args:
      - kw  test configuration see 'test.jinja.yaml'
    '''
    ST = f'[validate_otu]'
    vald = kw.get('validate')
    logd = kw.get('logd') or parse_log(LOGF)

    # parse otu groups file
    if os.path.exists(OTUF):
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        with open(OTUF) as f_otus:
            for line in f_otus:
                num_clusters_file += 1
                num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
        msg = f"{ST} num groups in OTU file {num_clusters_file} , expected {logd['num_otus'][1]}"
        print(msg)
        assert logd['num_otus'][1] == num_clusters_file, \
            f"{logd['num_otus'][1]} not equals {num_clusters_file}"
        msg = (f"{ST} count of reads in OTU file {num_reads_in_clusters_file},"
               f" expected {logd['num_id_cov'][1]}")
        print(msg)
        assert logd['num_id_cov'][1] == num_reads_in_clusters_file, \
            f"{logd['num_id_cov'][1]} not equals {num_reads_in_clusters_file}"

    # verify count of groups in the aligned.log is the same as specified in configuration validate data
    if logd.get('num_otus') and vald.get('num_cluster'):
        assert logd['num_otus'][1] in vald['num_cluster'], \
            f"{logd['num_otus'][1]} not in {vald['num_cluster']}"
#END validate_otu

def validate_log(logd:dict, ffd:dict):
    '''
    args:
      - logd
      - ffd   files data as in test.jinja.yaml:<test_name>:validate:files
    '''
    ST = '[validate_log]'
    # aligned.log : 
    #   verify the total number of reads in aligned.log vs the number in the validation spec
    n_vald = ffd.get('aligned.log', {}).get('num_reads')
    if n_vald:
        n_logd = logd['num_reads']
        msg = f'{ST} testing num_reads: {n_logd} Expected: {n_vald}'
        print(msg)
        if n_vald != n_logd:
            print(f'Failing: {n_vald} not equals {n_logd}')
    #   verify number of hits
    n_vald = ffd.get('aligned.log', {}).get('num_hits')
    if n_vald:
        n_logd = logd['num_hits']
        print(f'{ST} testing num_hits: {n_logd} Expected: {n_vald}')
        if n_vald != n_logd:
            print(f'Failing: {n_vald} not equals {n_logd}')
    #   verify number of misses
    n_vald = ffd.get('aligned.log', {}).get('num_fail')
    if n_vald:
        n_logd = logd['num_fail']
        print(f'{ST} testing num_fail: {n_logd} Expected: {n_vald}')
        if n_vald != n_logd:
            print(f'Failing: {n_vald} not equals {n_logd}')
    #   verify count of COV+ID
    n_vald = ffd.get('aligned.log', {}).get('n_yid_ycov')
    if n_vald:
        n_logd = logd.get('num_id_cov')
        print(f'{ST} testing n_yid_ycov: {n_logd} Expected: {n_vald}')
        if n_vald != n_logd: 
            print(f'Failing: {n_vald} not equals {n_logd}')
    #   verify count of OUT groups
    n_vald = ffd.get('aligned.log', {}).get('num_groups')
    if n_vald:
        n_logd = logd.get('num_otus')
        print(f'{ST} testing num_groups: {n_logd} Expected: {n_vald}')
        if n_vald != n_logd: 
            print(f'Failing: {n_vald} not equals {n_logd}')
    #   verify count of de-novo reads
    n_vald = ffd.get('aligned.log', {}).get('n_denovo')
    if n_vald:
        n_logd = logd['num_denovo']
        print(f'{ST} testing n_denovo: {n_logd} Expected: {n_vald}')
        if n_vald != n_logd: 
            print(f'Failing: {n_vald} not equals {n_logd}')
#END validate_log

def process_output(**kw):
    '''
    validate a test's output
    args:
      - kw      test configuration dictionary see e.g. 't18.jinja'
    '''
    ST = '[run.process_output]'
    global is_skbio
    vald = kw.get('validate')
    
    if not vald:
        print(f'{ST} validation info not provided')
        return
    
    outdir = Path(get_dict_val('args:-workdir', kw))/'out'
    logf = outdir /'aligned.log'
    logd = parse_log(logf)
    ffd = vald.get('files', {})
    for ff, vv in ffd.items():
        print(f'{ST} validating file {ff}')
        # Check aligned/other reads count
        # aligned/other files only specify read count in the test validation data
        if isinstance(vv, int):
            ffp = outdir / ff # file path
            count = 0
            if not ffp.exists():
                print(f'{ST} does not exists: {ffp}')
            if 'otu_map.txt' in ff:
                if ffp.exists():
                    with open(ffp) as ffs:
                        for line in ffs:
                            count += 1
                    msg = f'{ST} testing count of groups in {ff}: {count} Expected: {vv}'
                    print(msg)
                    if count != vv: 
                        print(f'Failing: {count} not equals {vv}')
                    continue
            if is_skbio:
                if IS_FASTQ:
                    for seq in skbio.io.read(ffp, format='fastq', variant=vald.get('variant')):
                        count += 1
                else:
                    fmt = 'fasta' if READS_EXT[1:] in ['fasta', 'fa'] else READS_EXT[1:]
                    for seq in skbio.io.read(ffp, format=fmt):
                        count += 1
                print(f'{ST} Testing count of reads in {ff}: {count} Expected: {vv}')
                if count != vv: 
                    print(f'Failing: {count} not equals {vv}')
        elif ff == 'aligned.log':
            validate_log(logd, ffd)
        elif 'aligned.blast' in ff:
            process_blast(**kw)
#END process_output

def t0(datad, ret={}, **kwarg):
    '''
    :param datad   Data directory
    :param outd    results output directory
    '''
    ST = f'[t0:{kwarg.get("name")}]'
    print(f'{ST} Validating ...')

    BLAST_EXPECTED = os.path.join(datad, 't0_expected_alignment.blast')

    dlist = []

    with open(BLASTF, 'r') as fout:
        with open(BLAST_EXPECTED, 'r') as fexpect:
            diff = difflib.unified_diff(
                fout.readlines(),
                fexpect.readlines(),
                fromfile='fout',
                tofile='fexpect',
            )
            for line in diff:
                dlist.append(line)
                if line: 
                    print(line, end='')

    assert len(dlist) == 0
    print(f'{ST} Done')
#END t0

def t2(datad, ret={}, **kwarg):
    '''
    :param datad   Data directory
    :param kwargs  validation args

    Test the following case for alignment:
    beginning from align_ref_start = 0 and align_que_start = X, the read finishes
    before the end of the reference
              ref |----------------|
    que |------------------------|
                    LIS |-----|
                  ^
                  align_que_start
    '''
    ST = f'[t2:{kwarg.get("name")}]'
    print(f'{ST} Validating ...')

    vald = kwarg['validate']

    actual_alignment = []
    with open(BLASTF) as afile:
        for line in afile:
            actual_alignment = line.strip().split('\t')
        
    assert len(vald['expected']) == len(actual_alignment)
    assert sorted(vald['expected']) == sorted(actual_alignment)
    #a = set(expected_alignment) & set(actual_alignment)
    print(f'{ST} Done')
#END t2

def t3(datad, ret={}, **kwarg):
    '''
    :param datad   Data directory
    :param outd    results output directory
    :param kwargs  validation args

    test_environmental_output

    Test outputting FASTA file for de novo clustering 
    using environmental data.

    Conditions: input FASTA file is processed in
                one mapped section.
    '''
    ST = f'[t3:{kwarg.get("name")}]'
    print(f'{ST} Validating ...')
    global is_skbio
    logd = parse_log(LOGF)
    vald = kwarg.get('validate')
    cmdd = kwarg.get('cmd')

    assert logd['results']['num_hits'][1] == vald['num_hits']

    process_blast(**kwarg)

    # sort order (descending/ascending) of candidate references (alignment.cpp)
    is_refs_descending = False
    # 'Total OTUs' in aligned.log has to be equal 'num_groups' in test.jinja.yaml
    if is_refs_descending:
        assert logd['num_otus'][1] == vald['num_groups'][0] # originally
    else:
        assert logd['num_otus'][1] == vald['num_groups'][1], \
            f'num_otus = {logd["num_otus"][1]} != num_groups = {vald["num_groups"][1]} expected'

    # OTU file contains one line per OTU group, so the number of lines
    # has to be equal 'Total OTUs' in aligned.log
    with open(OTUF) as f_otumap:
        num_clusters_file = sum(1 for line in f_otumap)
    assert logd['num_otus'][1] == num_clusters_file, \
        f'num_otus = {logd["num_otus"][1]} != {OTU_BASE}:num_otus = {num_clusters_file}'

    # number of reads in aligned_denovo.fasta has to be equal the
    # 'Total reads for de novo clustering' in aligned.log
    n_denovo_file = 0
    if is_skbio:
        for seq in skbio.io.read(DENOVOF, format='fasta'):
            n_denovo_file += 1

        assert logd['results']['num_denovo'][1] == n_denovo_file, \
                f'num_denovo = {logd["results"]["num_denovo"][1]} != {DENOVO_BASE}:num_denovo = {n_denovo_file}'
    
    print(f'{ST} Done')
#END t3

def t4(datad, ret={}, **kwarg ):
    '''
    count idx files
    '''
    ST = f'[t4:{kwarg.get("name")}]'
    vald = kwarg.get('validate')
    if IS_WIN:
        sfx = vald.get('idx_sfx_win')
    else:
        sfx = vald.get('idx_sfx_lin')
    idx_count = 0
    idx_count_expect = vald.get('num_idx') * 3 + 1
    if os.path.exists(IDX_DIR):
        idx_count = len([fn for fn in os.listdir(IDX_DIR) if str(sfx) in fn])

    print(f'{ST} Expected number of index files: {idx_count_expect} Actual number: {idx_count}')
    assert idx_count_expect == idx_count
#END t4

def t9(datad, ret={}, **kwarg):
    '''
    :param datad    data directory
    :param dict ret results of the aignment

    test_output_all_alignments_f_rc
    '''
    ST = f'[t9:{kwarg.get("name")}]'
    print(f'{ST} Validating ...')
    vald = kwarg.get('validate')
    sam_alignments = []
    with open(SAMF) as aligned_f:
        for line in aligned_f:
            if line.startswith('@'):
                continue
            alignment = line.strip().split("\t")
            sam_alignments.append(alignment)

    assert len(vald['sam_alignments_expected']) == len(sam_alignments)

    for alignment in vald['sam_alignments_expected']:
        assert alignment in sam_alignments
    
    print(f"{ST} done")
#END t9

def t10(datad, ret={}, **kwarg):
    '''
    :param str  datad  data directory
    :param dict ret    results of the aignment

    test_ref_shorter_than_seed
    '''
    ST = f'[t10:{kwarg.get("name")}]'
    print(f'{ST} validating ...'.format)

    vald = kwarg.get('validate')

    if ret and ret.get('retcode'):
        assert ret['retcode'] == 1
        assert vald['err_msg'] in ret['stderr'].decode("utf-8")
   
    print("{ST} done")
#END t10

def t11(datad, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    
    Test on simulated data,
        10000 reads with 1% error (--aligned),
        10000 reads with 10% error (de novo),
        10000 reads random (--other)

        Conditions: reference index and input
        query FASTA file both processed as one
        section.
    '''
    ST = f'[t11:{kwarg.get("name")}]'
    print(f'{ST} Validating ...')

    if ret and ret.get('retcode'):
        print(f'ERROR running alignemnt. Return code: {ret["retcode"]}')
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print(f'{ST} Done')
#END t11

def t12(datad, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    
    Test on simulated data,
        10000 reads with 1% error (--aligned),
        10000 reads with 10% error (de novo),
        10000 reads random (--other)

        Conditions: reference index and input
        query FASTA file both processed as one
        section.
    '''
    ST = f'[t12:{kwarg.get("name")}]'
    print(f'{ST} Validating ...')

    if ret and ret.get('retcode'):
        print(f'ERROR running alignemnt. Return code: {ret["retcode"]}')
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print(f'{ST} Done')
#END t12

def t17(datad, ret={}, **kwarg):
    '''
    :param name  sortmerna.exe path
    :param datad   Data directory
    :param outd    results output directory
    :param capture Capture output
    '''
    ST = f'[t17:{kwarg.get("name")}]'
    print(f'{ST} TODO: implement')
    logd = parse_log(LOGF)
    print("{ST} done")
#END t17

def set_file_names(basenames:list, is_other:bool=False):
    '''
    args:
      - basenames   list of basenames as in test.jinja,yaml:aligned.names
      - is_other    flag aligned (default) | other files
    '''
    global ALI_NAMES

def split(files:list,
          num_splits:int,
          size_per_thread:int=500000,
          rapidgz:str=None,
          pigz:str=None,
          workdir:str=None,
          **kw):
    '''
    args:
      - files       list of input files
      - num_splits  number of splits
      - rapidgz     path to rapidgz executable
      - pigz        path to pigz executable
      - workdir     working directory to output to
    '''
    if rapidgz and Path(rapidgz).exists():
        print(f'using provided rapidgzip executable: {rapidgz}')
    elif shutil.which('rapidgzip'):
        rapidgz = Path(shutil.which('rapidgzip'))
        print(f'found rapidgzip executable: {rapidgz}')
    else:
        msg = f'no rapidgzip executable found or provided. Use \'--rapidgz\' option'
        raise Exception(msg)

    if pigz and Path(pigz).exists():
        print(f'using provided pigz executable: {pigz}')
    elif shutil.which('pigz'):
        pigz = Path(shutil.which('pigz'))
        print(f'found pigz executable: {pigz}')
    else:
        msg = f'no pigz executable found or provided. Use \'--pigz\' option'
        raise Exception(msg)
    
    cpu_max = multiprocessing.cpu_count()  # CPUs count
    print(f'number of CPUs/cores on this machine: {cpu_max}')
    lines_per_read = 4
    
    # collect the input files statistics: size, number of lines, and reads
    #   calculate the number of reads
    print(f'processing {len(files)} files')
    tgtd = Path(workdir) / 'readb'
    if tgtd.exists() and any(tgtd.iterdir()):
        print(f'{tgtd} is not empty. Removing the content')
        for ff in tgtd.iterdir():
            print(f'removing {ff}')
            ff.unlink()
    elif not tgtd.exists():
        print(f'creating {tgtd}')
        tgtd.mkdir(parents=True)
        
    sstart = time.time()
    stats = []  # input files statistics
    for ff in files:
        sbytes = Path(ff).stat().st_size
        #sz_per_thread = sbytes // cpu_max
        threads = 2 if sbytes <= size_per_thread else sbytes // size_per_thread
        threads = 0 if threads >= cpu_max else threads
        tt = threads or cpu_max
        print(f'using {tt} threads for inflating {sbytes} bytes file')
        cmd = f'{rapidgz} --count-lines {ff}'
        #cmd1 = f'{rapidgz} -d -c -P {threads} {ff}'
        #cmd2 = 'wc -l'
        start = time.time()
        ps = subprocess.run(cmd.split(), capture_output=True)
        out = ps.stdout
        #serr = ps.stderr
        #ps1 = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
        #ps2 = subprocess.Popen(cmd2.split(), stdin=ps1.stdout, stdout=subprocess.PIPE)
        #ps1.stdout.close()  # allow SIGPIPE to terminate
        #out, serr = ps2.communicate()
        #out, serr = ps1.communicate()
        num_lines = int(out.decode().strip())
        reads = num_lines // lines_per_read
        if num_lines % lines_per_read != 0:
            print(f'Warning: number of lines {num_lines} is not multiple of lines per'
                  f' read {lines_per_read}. Extra lines: {num_lines % lines_per_read}')
        stats.append([ff, sbytes, reads])
        print(f'processed file {ff} in: {time.time() - start:.2f} sec. Size: {sbytes}, reads: {reads}')
    if len(stats) == 2 and stats[0][2] != stats[1][2]:
        print(f'error: the files contain unequal number of records {stats[0][2]} != {stats[1][2]}')
        
    total_reads = sum([st[2] for st in stats])
    print(f'total number of reads in {len(files)} files: {total_reads}')
        
    # split files
    stats_out = []  # split file statistics
    for i, ss in enumerate(stats):
        sz_split = ss[1] // num_splits  # size of a split
        reads_per_split_min = ss[2] // num_splits
        remainder_reads = ss[2] % num_splits
        # if enough CPUs use 1 thread per 100KB
        # if file size is less then 100KB - use 2 threads
        threads = 2 if sz_split <= size_per_thread else sz_split // size_per_thread
        threads = 0 if threads >= cpu_max else threads
        tt = threads or cpu_max
        print(f'using {tt} threads for inflating {sz_split} bytes split')
        offset = 0
        for j in range(num_splits):
            # distribute remainder reads (if any) between splits by adding an extra read 
            # to each split starting from the 0th up to a max (num_splits - 2)th
            reads_per_split = reads_per_split_min + 1 if j < remainder_reads else reads_per_split_min
            ln_split = reads_per_split * lines_per_read
            cmd1 = f'{rapidgz} -d -c -P {threads} --ranges {ln_split}L@{offset}L {ss[0]}'
            cmd2 = f'pigz -c -p {tt}'
            fname = f'fwd_{j}.fq.gz' if i == 0 else f'rev_{j}.fq.gz'
            fout = tgtd / fname
            start = time.time()
            with open(fout, 'w') as fh:
                ps1 = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
                print(f'writing {ln_split} lines to {fout} using {tt} threads')
                ps2 = subprocess.Popen(cmd2.split(), stdin=ps1.stdout, stdout=fh)
                ps1.stdout.close()  # allow SIGPIPE to terminate
                out, serr = ps2.communicate()
            rec = [fout.as_posix(), fout.stat().st_size, reads_per_split, 1, 'fastq']
            stats_out.append(rec)
            offset += ln_split
            print(f'processed file {rec[0]} in: {time.time() - start:.2f} sec. Size: {rec[1]}, reads: {rec[2]}')
            
    # write the descriptor
    readfeed = (
        '# format of this file:\n'
        '#   time\n'
        '#   num_orig_files\n'
        '#   num_sense\n'
        '#   num_splits\n'
        '#   num_reads_tot\n'
        '#   [\n'
        '#     file\n'
        '#     size\n'
        '#     reads\n'
        '#     zip\n'
        '#     fastq/a\n'
        '#   ] for each file both original and split\n\n'
    )
    readfeed += (f'{time.strftime("%a %d %b %Y %H:%M:%S", time.gmtime())}\n\n')
    
    # fwd
    readfeed += (
        f'{len(stats)}\n'  # number of original files
        f'{len(stats)}\n'  # number of senses
        f'{num_splits}\n'  # number of splits
        f'{total_reads}\n'  # total number of reads
        f'{stats[0][0]}\n'  # fwd reads file
        f'{stats[0][1]}\n'  # file size
        f'{stats[0][2]}\n'  # number of reads in file
        '1\n'  # compressed file
        'fastq\n'  # file format
    )
    # rev
    if len(stats) == 2:
        readfeed += (
            f'{stats[1][0]}\n'  # fwd reads file
            f'{stats[1][1]}\n'  # file size
            f'{stats[1][2]}\n'  # lines in file
            '1\n'  # compressed file
            'fastq\n'  # file format
        )
    # splits statistics
    for i in range(num_splits):
        # fwd
        readfeed += (
            f'{stats_out[i][0]}\n'  # fwd reads file split
            f'{stats_out[i][1]}\n'  # file size
            f'{stats_out[i][2]}\n'  # reads in file
            '1\n'  # compressed file
            'fastq\n'  # file format
        )
        # rev
        if len(stats) == 2:
            j = i + num_splits
            readfeed += (
                f'{stats_out[j][0]}\n'  # fwd reads file split
                f'{stats_out[j][1]}\n'  # file size
                f'{stats_out[j][2]}\n'  # reads in file
                '1\n'  # compressed file
                'fastq\n'  # file format
            )
        
    dfile = tgtd / 'readfeed'
    print(f'writing to {dfile}')
    with open(dfile, 'a') as fh:
        fh.write(readfeed)
    print(f'done split in {time.time() - sstart:.2f} sec')
    return stats

def process_config(**kw) -> dict:
    '''
    '''
    global SMR_EXE
    script_dir = Path(__file__).parent # directory where this script is located
    SMR_SRC = Path(script_dir).parent
    if args.data_dir:
        DATA_DIR = args.data_dir
        print(f'{ST} smr source dir: {SMR_SRC} script dir: {script_dir}')

    # check env.yaml. If no env file specified, try the current directory
    env_yaml = Path(script_dir) / 'env.jinja' if not args.envfile else args.envfile
    if not Path(env_yaml).exists():
        msg = f'{ST} No environment config file found. Please, provide one using \'--env\' option'
        raise Exception(msg)

    # load properties from env.yaml
    print(f'{ST} Using Environment configuration file: {env_yaml}')
    env_jj = Environment(loader=FileSystemLoader(Path(env_yaml).parent), 
                                                    trim_blocks=True, 
                                                    lstrip_blocks=True)
    env_template = env_jj.get_template(Path(env_yaml).name)

    # render jinja template env.jinja.yaml
    vars = {'UHOME': UHOME}
    if IS_WSL: vars['WINHOME'] = args.winhome
    env_str = env_template.render(vars)
    env = yaml.load(env_str, Loader=yaml.FullLoader)

    # check test.jinja
    cfgfile = Path(script_dir) / 'test.jinja' if not args.config else args.config
    if not Path(cfgfile).exists():
        msg = f'{ST} No build configuration template found. Please, provide one using \'--config\' option'
        raise Exception(msg)

    print(f'{ST} using tests configuration template: {cfgfile}')

    # load 'test.jinja' template
    jjenv = Environment(loader=FileSystemLoader(Path(cfgfile).parent), trim_blocks=True, lstrip_blocks=True)
    template = jjenv.get_template(Path(cfgfile).name)

    ENV = args.envname or OS if IS_WIN or IS_WSL else None
    if ENV:
        print(f'{ST} using ENV: {ENV}')
    else:
        print(f'{ST} --envn is required on OS: {OS}')
        is_opts_ok = False

    WRK_DIR = args.workdir or f'{SMR_SRC}/run'

    # render 'test.jinja' template
    val = env.get(SMR,{}).get('src',{}).get(ENV)
    vars = {'SMR_SRC':SMR_SRC, 'DATA_DIR':DATA_DIR, 'WRK_DIR':WRK_DIR}
    cfg_str = template.render(vars)
    #cfg_str = template.render(env) # env[OS]
    cfg = yaml.load(cfg_str, Loader=yaml.FullLoader)
    
    SMR_EXE = args.smr_exe or shutil.which('sortmerna')
    if not SMR_EXE or not os.path.exists(SMR_EXE):
        print('{ST} sortmerna executable {SMR_EXE} not found.'
              ' Use \'--smr-exe\' option or ensure the PATH is set')
        sys.exit(1)
    
    print(f'{ST} using {SMR_EXE}')

    TEST_DATA = Path(SMR_SRC)/'data'  # data directory
    ALI_NAMES = cfg.get('aligned.names')
    
    return cfg

def count_lines(files:list,
          num_splits:int=3,
          size_per_thread:int=500000):
    '''
    args:
      - files  list of .gz files to split
      
    threads  records  sec  bytes/thread
    -----------------------------------
    2        20000    21   500000
    4        20000    24   300000
    6        20000    27   200000    too many threads slow the execution down
    '''
    lines_tot = 0
    cpu_max = os.cpu_count()  # multiprocessing.cpu_count()
    ts = time.time()
    for ff in files:
        sbytes = Path(ff).stat().st_size
        threads = 2 if sbytes <= size_per_thread else sbytes // size_per_thread
        threads = 0 if threads >= cpu_max else threads
        # count lines
        print(f'counting lines using {threads or cpu_max} threads')
        with rapidgzip.open(ff, parallelization=threads) as file:
            #lines_tot = sum(1 for _ in file)  # super inefficient. Orders of magnitude.
            while True:
                chunk = file.read(16 * 1024 * 1024)  # 16MB
                if not chunk:
                    break
                lines_tot += chunk.count(b'\n')
            print(f'counted {lines_tot} lines in {time.time()-ts} sec')
    return lines_tot

def process(file:str, num_lines:int, offset:int):
    '''
    args:
      - num_lines  number of lines to process
      - offset     starting offset
    '''
    print(f'processing {num_lines} at offset {offset} in {file}')

def split_process(file:str, num_threads:int=2):
    '''
    '''
    import threading
    threads = []
    lines_tot = count_lines()
    lines_per_thread = lines_tot // num_threads
    ex_lines = lines_tot % num_threads
    lines = [lines_per_thread for _ in range(num_threads)]
    for i in range(ex_lines):
        lines[i] += 1
    
    for i in range(num_threads):
        thr = threading.Thread(target=process, args=[file, lines[i], i+lines[i]])
        threads.append(thr)
        
    for t in threads:
        t.start()

    for t in threads:
        t.join()
    print(f'done processing')
    
def iss_453(num_splits:int,
            reads_mln:int=15,
            threads:int=8,
            pigz:str=None,
            workdir:str=None,
            **kw):
    '''
    generate gzip files for testing issue 453
    '''
    ST = 'iss_453'
    import itertools
    # check pigz
    if pigz and Path(pigz).exists():
        print(f'using provided pigz executable: {pigz}')
    elif shutil.which('pigz'):
        pigz = Path(shutil.which('pigz'))
        print(f'found pigz executable: {pigz}')
    else:
        msg = f'no pigz executable found or provided. Use \'--pigz\' option'
        raise Exception(msg)
    
    # prepare workdir
    tgtd = Path(workdir) if num_splits == 1 else Path(workdir) / 'readb'
    if num_splits > 1 and tgtd.exists() and any(tgtd.iterdir()):
        print(f'{tgtd} is not empty. Removing the content')
        for ff in tgtd.iterdir():
            print(f'removing {ff}')
            ff.unlink()
    elif not tgtd.exists():
        print(f'creating {tgtd}')
        tgtd.mkdir(parents=True)
    
    cpu_max = multiprocessing.cpu_count()  # CPUs count
    print(f'number of CPUs/cores on this machine: {cpu_max}')
    #lines_per_read = 4
    if num_splits == 1:
        splits_len = [reads_mln*10**6]
    else:
        reads_per_split_min = reads_mln*10**6 // num_splits
        remainder_reads = reads_mln*10**6 % num_splits
        splits_len_min = [reads_per_split_min for _ in range(num_splits)]
        splits_remain = [1 for _ in range(remainder_reads)]
        splits_len = [sum(x) for x in itertools.zip_longest(splits_len_min, splits_remain, fillvalue=0)]
    tt = threads or cpu_max
    print(f'using {tt} threads')
    
    # write split reads files
    cmd = f'pigz -c -p {tt}'
    read_fq = 'GGCGGCGTCCGGTGAGCTCTCGCTGGCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    for i, split_reads_num in enumerate(splits_len):
        fname = f'fwd_{i}.fq.gz' if num_splits > 1 else f'{ST}_{reads_mln}M.fq.gz'
        fout = tgtd / fname
        if num_splits == 1 and fout.exists():
            print(f'{ST} removing existing file {fout}')
            fout.unlink()
        print(f'writing {split_reads_num} reads to {fout} using {tt} threads')
        sstart = time.time()
        with open(fout, 'a') as fh:
            ps = subprocess.Popen(cmd.split(), stdin=subprocess.PIPE, stdout=fh)
            cnt, chunk = 0, ''
            for j in range(split_reads_num):
                chunk = f'@{j}\n{read_fq}\n' if j < split_reads_num - 1 else f'@{j}\n{read_fq}'
                #cnt += 1
                #if cnt == 10:
                ps.stdin.write(chunk.encode('utf-8'))
                #cnt = 0
                #chunk = ''
            out, serr = ps.communicate()
        print(f'done writing in {time.time() - sstart} sec')
    
def _apply_preset_argv(args: Namespace, preset_argv: list):
    '''Fill None-valued attributes in args from an argv-style preset list.
    Explicitly provided args (non-None) are never overwritten.
    '''
    i = 0
    while i < len(preset_argv):
        token = str(preset_argv[i])
        if token.startswith('--'):
            key = token[2:].replace('-', '_')
        elif token.startswith('-') and len(token) > 1:
            key = token[1:].replace('-', '_')
        else:
            i += 1
            continue
        if i + 1 < len(preset_argv) and not str(preset_argv[i + 1]).startswith('-'):
            val = preset_argv[i + 1]
            i += 2
        else:
            val = True
            i += 1
        if getattr(args, key, None) is None:
            setattr(args, key, val)

def parse_test_config(args:Namespace) -> dict:
    '''
    parse test jinja template and set all sortmerna options
    args:
      - args  arguments produced by argparse package
    '''
    script_dir = Path(__file__).parent # directory where this script is located

    # load presets and apply defaults for any args not explicitly provided
    cfg_dir = Path(args.config).parent if args.config else script_dir
    presets_path = Path(args.presets) if getattr(args, 'presets', None) else cfg_dir / 'presets.yaml'
    if presets_path.exists():
        msg = f'{ST} using presets from {presets_path}'
        print(msg)
        with open(presets_path) as _f:
            _all_presets = yaml.safe_load(_f) or {}
        _preset_argv = _all_presets.get(args.name, [])
        if _preset_argv:
            _apply_preset_argv(args, _preset_argv)

    smr_exe = args.smr_exe or shutil.which('sortmerna')
    if not smr_exe or not Path(smr_exe).exists():
        msg = (f'{ST} sortmerna executable {smr_exe} not found.'
              f' Use \'--smr-exe\' option or ensure the PATH is set')
        raise Exception(msg)
    else:
        print(f'{ST} using {smr_exe}') 

    # check test.jinja
    cfgfile =  args.config or Path(script_dir) / f'{args.name}.jinja'
    if not Path(cfgfile).exists():
        msg = (f"{ST} No test template found: {cfgfile}"
               f" Please, provide one using \'--config\' option")
        raise Exception(msg)

    print(f'{ST} using tests configuration template: {cfgfile}')

    # load 'test.jinja' template
    jjenv = Environment(loader=FileSystemLoader(Path(cfgfile).parent), trim_blocks=True, lstrip_blocks=True)
    template = jjenv.get_template(Path(cfgfile).name)

    # render 'test.jinja' template
    smr_src = Path(script_dir).parent
    vars = {'SMR_SRC':smr_src, 'DATA_DIR':args.data_dir, 'WRK_DIR':args.workdir}
    if args.threads:
        vars['THREADS'] = str(args.threads)
    if args.ref_dir:
        vars['REF_DIR'] = args.ref_dir
    if args.task:
        vars['TASK'] = str(args.task)
    if args.dbg_level:
        vars['DBG_LEVEL'] = args.dbg_level
    if args.index:
        vars['INDEX'] = str(args.index)
    cfg_str = template.render(vars)
    cfg = yaml.load(cfg_str, Loader=yaml.FullLoader)
    if getattr(args, 'score_split', False):
        cfg['args'].append('-score_split')
        
    cfg['exe'] = smr_exe
    return cfg

def get_dict_val(path:str, nested_dict:dict):
    '''
    get nested value from the dictionary given the path like 'validate:failes:aligned.log'
    '''
    vv = nested_dict
    for k in path.split(':'):
        vv = vv.get(k) or {}
    return vv

def _parse_stats_binary(fpath: Path):
    '''
    Parse sortmerna *.stats binary file.
    Returns (stats_dict, fasta_filename_str).

    Binary layout (native little-endian on LP64 Linux):
      size_t      filesize
      uint32_t    fasta_len
      char[]      fasta_filename  (fasta_len bytes, null-terminated)
      double[4]   background_freq (A, C, G, T)
      uint64_t    full_len
      uint32_t    seed_win_len
      uint64_t    numseq
      uint16_t    part_num
      index_parts_stats[part_num]:
          unsigned long  start_part    (8 bytes on LP64)
          unsigned long  seq_part_size (8 bytes on LP64)
          uint32_t       numseq_part
          uint32_t       _pad          (compiler trailing padding)
      uint32_t    num_sq
      [for each SAM SQ header:]
          uint32_t  len_id
          char[]    sequence_id (len_id bytes, no null terminator)
          uint32_t  len_seq
    '''
    class _PartStats(ctypes.Structure):
        _fields_ = [
            ('start_part',    ctypes.c_ulong),
            ('seq_part_size', ctypes.c_ulong),
            ('numseq_part',   ctypes.c_uint32),
        ]
    part_size = ctypes.sizeof(_PartStats)  # 24 on LP64 (4 bytes trailing padding)

    with open(fpath, 'rb') as f:
        (filesize,)     = struct.unpack('<Q', f.read(8))
        (fasta_len,)    = struct.unpack('<I', f.read(4))
        fasta_filename  = f.read(fasta_len).rstrip(b'\x00').decode('utf-8', errors='replace')
        bg_freq         = struct.unpack('<4d', f.read(32))
        (full_len,)     = struct.unpack('<Q', f.read(8))
        (seed_win_len,) = struct.unpack('<I', f.read(4))
        (numseq,)       = struct.unpack('<Q', f.read(8))
        (part_num,)     = struct.unpack('<H', f.read(2))
        parts = []
        for _ in range(part_num):
            raw = f.read(part_size)
            start_part, seq_part_size, numseq_part = struct.unpack_from('<QQI', raw)
            parts.append({'numseq_part': numseq_part, 'seq_part_size': seq_part_size})

    stats = {
        'numseq':       numseq,
        'full_len':     full_len,
        'part_num':     part_num,
        'seed_win_len': seed_win_len,
        'background_freq': {
            'A': round(bg_freq[0], 4),
            'C': round(bg_freq[1], 4),
            'G': round(bg_freq[2], 4),
            'T': round(bg_freq[3], 4),
        },
        'parts': parts,
    }
    return stats, fasta_filename


def _parse_kmer_binary(fpath: Path) -> dict:
    '''
    Parse *.kmer_N.dat: (1 << seed_win_len) consecutive uint32_t kmer counts.
    '''
    data = fpath.read_bytes()
    n = len(data) // 4
    counts = struct.unpack(f'<{n}I', data)
    return {
        'num_nonzero': sum(1 for c in counts if c > 0),
        'total_count': sum(counts),
        'max_count':   max(counts) if counts else 0,
    }


def _parse_pos_binary(fpath: Path) -> dict:
    '''
    Parse *.pos_N.dat:
      uint32_t number_elements
      for each element:
        uint32_t size
        seq_pos[size]  (seq_pos = {uint32_t pos, uint32_t seq} = 8 bytes each)
    '''
    total_positions = 0
    max_positions   = 0
    with open(fpath, 'rb') as f:
        raw = f.read(4)
        if len(raw) < 4:
            return {}
        (num_elements,) = struct.unpack('<I', raw)
        for _ in range(num_elements):
            raw = f.read(4)
            if len(raw) < 4:
                break
            (size,) = struct.unpack('<I', raw)
            total_positions += size
            if size > max_positions:
                max_positions = size
            f.seek(size * 8, 1)  # skip seq_pos[size]: each is 2 x uint32_t
    return {
        'num_elements':    num_elements,
        'total_positions': total_positions,
        'max_positions':   max_positions,
    }


def _read_idx_stats_yaml(yaml_file: Path) -> dict:
    '''
    Parse a sortmerna-generated *.idx_stats.yaml file into the same dict
    structure used by _collect_ref_index_stats.  The YAML has top-level keys
    reference / stats / parts; we reformat into the validate:indices schema.
    '''
    with open(yaml_file) as f:
        raw = yaml.safe_load(f)

    part0 = raw['parts'][0] if raw.get('parts') else {}
    result = {
        'stats': {
            'numseq':       raw['stats']['numseq'],
            'full_len':     raw['stats']['full_len'],
            'part_num':     raw['stats']['part_num'],
            'seed_win_len': raw['stats']['seed_win_len'],
            'background_freq': raw['stats']['background_freq'],
            'parts': [
                {'numseq_part': p['numseq_part'], 'seq_part_size': p['seq_part_size']}
                for p in raw.get('parts', [])
            ],
        },
        'files': {
            'stats': raw['stats']['files']['stats'],
            **{k: v for p in raw.get('parts', []) for k, v in p.get('files', {}).items()},
        },
        'kmer': part0.get('kmer', {}),
        'pos':  part0.get('pos', {}),
    }
    return result


def _collect_ref_index_stats(idx_dir: Path, ref_path: str) -> dict:
    '''
    Locate the index files for ref_path inside idx_dir. The index prefix is a
    std::hash of the ref basename, not derivable in Python, so we scan matching
    files to find the right one.

    Preference order:
      1. *.idx_stats.yaml written by sortmerna at build time (fast, no binary parsing)
      2. Binary *.stats + *.kmer_0.dat + *.pos_0.dat parsed directly (fallback)
    '''
    ST = '[_collect_ref_index_stats]'
    ref_resolved = str(Path(ref_path).resolve())

    # --- try YAML first (written by sortmerna >= the version that added this) ---
    for yf in sorted(idx_dir.glob('*.idx_stats.yaml')):
        try:
            raw = yaml.safe_load(yf.read_text())
            stored = raw.get('reference', '')
            if str(Path(stored).resolve()) == ref_resolved or stored == ref_path:
                return _read_idx_stats_yaml(yf)
        except Exception as ex:
            print(f'{ST} warning: could not parse {yf}: {ex}')

    # --- fallback: parse binary index files ---
    matched_prefix = None
    matched_stats  = None
    sf = None

    for sf in sorted(idx_dir.glob('*.stats')):
        try:
            sd, fasta_fn = _parse_stats_binary(sf)
        except Exception as ex:
            print(f'{ST} warning: could not parse {sf}: {ex}')
            continue
        if str(Path(fasta_fn).resolve()) == ref_resolved or fasta_fn == ref_path:
            matched_prefix = sf.with_suffix('')
            matched_stats  = sd
            break

    if matched_prefix is None:
        print(f'{ST} no index found for {ref_path} in {idx_dir}')
        return {}

    result = {'stats': matched_stats}

    files_stats = {'stats': sf.stat().st_size}
    for part_idx in range(matched_stats['part_num']):
        for pfx in ('bursttrie', 'pos', 'kmer'):
            fname = f'{pfx}_{part_idx}.dat'
            fp = Path(f'{matched_prefix}.{fname}')
            if fp.exists():
                files_stats[fname] = fp.stat().st_size
    result['files'] = files_stats

    kmer_file = Path(f'{matched_prefix}.kmer_0.dat')
    if kmer_file.exists():
        result['kmer'] = _parse_kmer_binary(kmer_file)

    pos_file = Path(f'{matched_prefix}.pos_0.dat')
    if pos_file.exists():
        result['pos'] = _parse_pos_binary(pos_file)

    return result


def _build_indices_yaml(indices_stats: dict, indent: int = 2) -> str:
    '''Render indices_stats as an indented YAML block starting with "indices:".'''
    p = ' ' * indent
    p2 = p + '  '
    p4 = p2 + '  '
    p6 = p4 + '  '
    lines = [f'{p}indices:']
    for ref_base, rdata in indices_stats.items():
        lines.append(f'{p2}{ref_base}:')
        if 'stats' in rdata:
            s = rdata['stats']
            bf = s['background_freq']
            lines.append(f'{p4}stats:')
            lines.append(f'{p6}numseq: {s["numseq"]}')
            lines.append(f'{p6}full_len: {s["full_len"]}')
            lines.append(f'{p6}part_num: {s["part_num"]}')
            lines.append(f'{p6}seed_win_len: {s["seed_win_len"]}')
            lines.append(f'{p6}background_freq: {{A: {bf["A"]}, C: {bf["C"]}, G: {bf["G"]}, T: {bf["T"]}}}')
            if s.get('parts'):
                lines.append(f'{p6}parts:')
                for part in s['parts']:
                    lines.append(f'{p6}  - numseq_part: {part["numseq_part"]}')
                    lines.append(f'{p6}    seq_part_size: {part["seq_part_size"]}')
        if 'files' in rdata:
            lines.append(f'{p4}files:')
            for fname, size in rdata['files'].items():
                lines.append(f'{p6}{fname}: {size}')
        if 'kmer' in rdata:
            k = rdata['kmer']
            lines.append(f'{p4}kmer:')
            lines.append(f'{p6}num_nonzero: {k["num_nonzero"]}')
            lines.append(f'{p6}total_count: {k["total_count"]}')
            lines.append(f'{p6}max_count: {k["max_count"]}')
        if 'pos' in rdata:
            pos = rdata['pos']
            lines.append(f'{p4}pos:')
            lines.append(f'{p6}num_elements: {pos["num_elements"]}')
            lines.append(f'{p6}total_positions: {pos["total_positions"]}')
            lines.append(f'{p6}max_positions: {pos["max_positions"]}')
    return '\n'.join(lines) + '\n'


def _update_jinja_indices(jinja_file: Path, indices_stats: dict):
    '''
    Insert or replace the validate:indices block in a jinja template file.
    Inserts before "runtime:" if present, otherwise at the end of validate block.
    '''
    ST = '[_update_jinja_indices]'
    lines = jinja_file.read_text().splitlines(keepends=True)

    # Locate top-level "validate:" line (indent 0)
    validate_line = None
    for i, line in enumerate(lines):
        if re.match(r'^validate\s*:', line):
            validate_line = i
            break
    if validate_line is None:
        print(f'{ST} no validate: block found in {jinja_file}')
        return

    # Scan validate's direct children (indent=2) to find existing 'indices:' and 'runtime:'
    indices_start  = None  # line index of '  indices:'
    indices_end    = None  # line index where indices block ends (exclusive)
    insert_before  = len(lines)  # default: end of file

    i = validate_line + 1
    while i < len(lines):
        raw     = lines[i]
        stripped = raw.lstrip()
        indent  = len(raw) - len(stripped)
        is_blank = not stripped.strip()

        if is_blank:
            i += 1
            continue

        # A top-level key (indent=0) means the validate block ended
        if indent == 0:
            if indices_start is not None and indices_end is None:
                indices_end = i
            insert_before = i
            break

        # Direct child of validate (indent=2)
        if indent == 2:
            key = stripped.split(':')[0].strip()
            if key == 'indices':
                indices_start = i
                # Find the end of the indices block: next non-blank line at indent <= 2
                j = i + 1
                while j < len(lines):
                    r2 = lines[j]
                    s2 = r2.lstrip()
                    ind2 = len(r2) - len(s2)
                    if s2.strip() and ind2 <= 2:
                        indices_end = j
                        break
                    j += 1
                else:
                    indices_end = len(lines)
                i = indices_end  # skip past the block
                continue
            elif key == 'runtime' and indices_start is None:
                # insert new block before runtime:
                insert_before = i
        i += 1

    new_block = _build_indices_yaml(indices_stats, indent=2)

    if indices_start is not None:
        new_lines = lines[:indices_start] + [new_block] + lines[indices_end:]
    else:
        new_lines = lines[:insert_before] + [new_block] + lines[insert_before:]

    jinja_file.write_text(''.join(new_lines))
    print(f'{ST} updated {jinja_file}')


def generate_index_stats(args: Namespace):
    '''
    Build the sortmerna index for a test's reference files (-task 5),
    collect statistics from the generated index files, and write them
    into the validate:indices section of the test's jinja template.

    CLI: python scripts/run.py index-stats --test t3 --smr-exe dist/bin/sortmerna
         --data-dir data/ --ref-dir data/rRNA_databases -w /tmp/smr_t3
    '''
    ST = '[generate_index_stats]'

    cfg     = parse_test_config(args)
    smr_exe = cfg['exe']
    smr_args = cfg.get('args', {})

    refs    = smr_args.get('-ref', [])
    reads   = smr_args.get('-reads', [])
    workdir = smr_args.get('-workdir') or args.workdir

    if not refs:
        print(f'{ST} no -ref entries found in test config')
        return
    if not workdir or str(workdir) == 'None':
        print(f'{ST} workdir is required; pass --workdir')
        return

    # Build index-only command
    cmd = [smr_exe]
    for ref in refs:
        cmd += ['-ref', ref]
    for read in reads:
        cmd += ['-reads', read]
    cmd += ['-workdir', str(workdir), '-task', '5']

    rcode, _, _ = run_test(cmd)
    if rcode != 0:
        print(f'{ST} sortmerna exited with code {rcode}; aborting')
        return

    idx_dir = Path(workdir) / 'idx'
    if not idx_dir.exists():
        print(f'{ST} idx directory not found: {idx_dir}')
        return

    indices_stats = {}
    for ref in refs:
        ref_base = Path(ref).stem  # e.g. silva-bac-16s-database-id85
        rdata = _collect_ref_index_stats(idx_dir, ref)
        if rdata:
            indices_stats[ref_base] = rdata
        else:
            print(f'{ST} skipping {ref_base}: no index stats found')

    if not indices_stats:
        print(f'{ST} no index stats collected; jinja file not modified')
        return

    script_dir = Path(__file__).parent
    jinja_file = Path(args.config) if getattr(args, 'config', None) else script_dir / f'{args.name}.jinja'
    _update_jinja_indices(jinja_file, indices_stats)
    print(f'{ST} done')


if __name__ == "__main__":
    '''
    cmd:
      - python scripts/run.py split -f file1 -f file2 --num_splits 4
      - $env:path="$env:userprofile/a1/code/sortmerna/dist/bin;$env:path"; &python scripts/run.py -n t1
      - python scripts/run.py test -n t1 --smr-exe=dist/bin/sortmerna.exe
      - python scripts/run.py test -n t15 --smr-exe=dist/bin/sortmerna.exe --data-dir=../../data

      - python scripts/run.py test -n t0 -v [--capture] [--validate-only]
      - python scripts/run.py test -n t12 -f process_otu --validate-only -d 2
      - python scripts/run.py test -n t0 -f dbg_blast --validate-only
      - python /media/sf_a01_code/sortmerna/scripts/run.py -n t6 --envn LNX_VBox_Ubuntu_1804
      - python /mnt/c/Users/biocodz/a01_code/sortmerna/tests/run.py -n t0 --winhome /mnt/c/Users/biocodz [--capture]  
                                                                            |_ on WSL
    '''
    ST = '[run.py:__main__]'
    is_opts_ok = True

    # process options
    p0 = ArgumentParser()
    p0.add_argument('--config', dest='config', help='Configuration file')
    p0.add_argument('--data-dir', dest='data_dir', help='path to the data. Abs or relative')
    p0.add_argument('--env', dest='envfile', help='Environment variables')
    subpar = p0.add_subparsers(dest='cmd', help='subcommand help')
                                #  |_this is the key parameter to access the commands as args.cmd.
    # split using rapidgzip
    p1 = subpar.add_parser('split', help='split reads and generate split descriptor')
    p1.add_argument('--rapidgz', help='rapidgzip executable path')
    p1.add_argument('--pigz', help='pigz executable path')
    p1.add_argument('-f', '--file', action='append', help='the input file')
    p1.add_argument('-s', '--num-splits', dest='num_splits', type=int, help='number of splits to generate')
    p1.add_argument('-c', '--size-per-thread', dest='size_per_thread', type=str, 
                        help='size of chunk processed by a single thread')
    p1.add_argument('-w', '--workdir', dest='workdir', help='working directory path')
    
    p3 = subpar.add_parser('split_process', help='process file using multiple threads, a chunk per threads')
    p3.add_argument('-f', '--file', action='append', help='the input file')
    p3.add_argument('-t', '--num-threads', dest='num_threads', type=int, help='number of threads to use')
    
    p4 = subpar.add_parser('iss_453', help = 'test issue 453')
    p4.add_argument('--pigz', help='pigz executable path')
    p4.add_argument('-r', '--reads-mln', dest='reads_mln', type=int, help='number of mln of reads')
    p4.add_argument('-s', '--num-splits', dest='num_splits', type=int, help='number of splits to generate')
    p4.add_argument('-t', '--threads', dest='threads', type=int, 
                        help='umber of threads to use')
    p4.add_argument('-w', '--workdir', dest='workdir', help='working directory path')
    
    # run tests
    p5 = subpar.add_parser('test', help='run selected test')
    p5.add_argument('name', help=('Test to run: single name (t3), comma-separated list (t3,t18,t46), '
                                  'or group alias resolved from presets.yaml (all -> presets.yaml:all, '
                                  'all-big -> presets.yaml:all_big). Aliases may be mixed with explicit names.'))
    p5.add_argument('--smr-exe', dest='smr_exe', help='path to sortmerna executable. Abs or relative')
    p5.add_argument('--data-dir', dest='data_dir', help='path to the data. Abs or relative')
    p5.add_argument('--ref-dir', dest='ref_dir', help='path to the reference data. Abs or relative')
    p5.add_argument('--threads', dest='threads', help='Number of threads to use')
    p5.add_argument('--index', dest='index', help='Index option 0 | 1 | 2')
    p5.add_argument('-t', '--task', dest='task', help='Processing task 0 | 1 | 2 | 3 | 4')
    p5.add_argument('-d', '--dbg-level', dest="dbg_level", help='debug level 0 | 1 | 2')
    p5.add_argument('-w', '--workdir', dest='workdir', help='Environment variables')
    p5.add_argument('-c', '--clean', action="store_true", help='clean Work directory and exit. Requires \'--name\'')
    p5.add_argument('-v', '--validate-only', action="store_true", help='Only perform validation. Assumes aligement already done')
    p5.add_argument('-f', '--func', dest='func', help='function to run: process_otu | ')
    p5.add_argument('--btype', dest='btype', default='release', help = 'Build type: release | debug')
    p5.add_argument('--pt_smr', dest='pt_smr', default='t1', help = 'Sortmerna Linkage type t1 | t2 | t3')
    p5.add_argument('--winhome', dest='winhome', help='when running on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    p5.add_argument('--capture', action="store_true", help='Capture output. By default prints to stdout')
    p5.add_argument('--config', dest='config', help='Tests configuration file.')
    p5.add_argument('--presets', dest='presets', help='Tests presets file.')
    p5.add_argument('--env', dest='envfile', help='Environment variables')
    p5.add_argument('-e','--envn', dest='envname', help=('Name of environment: WIN | WSL '
                                                      '| LNX_AWS | LNX_TRAVIS | LNX_VBox_Ubuntu_1804 | ..'))
    p5.add_argument('--score-split', action="store_true", help='set corresponding sortmerna argument')
    p5.add_argument('--interrupt', dest='interrupt', nargs='?', const=1, default=0, type=int,
                    choices=[0, 1, 2],
                    help='Inject N mid-alignment SIGKILLs to exercise auto-resume. '
                         '0 (default) runs normally; --interrupt (no value) implies 1; max 2.')
    p5.add_argument('--interrupt-wait', dest='interrupt_wait', type=int, default=60,
                    help='Seconds to let alignment run after threads start before each SIGKILL '
                         '(default: 60). Used only when --interrupt > 0.')
    p5.add_argument('--stop-on-fail', dest='stop_on_fail', action="store_true",
                    help='Abort the batch on the first failing test (default: continue and summarize at end)')

    # build index and record statistics into test jinja
    p6 = subpar.add_parser('index-stats', help='build index and record statistics into test jinja')
    p6.add_argument('--test', dest='name', required=True, help='Test name e.g. t3')
    p6.add_argument('--smr-exe', dest='smr_exe', help='path to sortmerna executable')
    p6.add_argument('--data-dir', dest='data_dir', help='path to the data. Abs or relative')
    p6.add_argument('--ref-dir', dest='ref_dir', help='path to the reference data. Abs or relative')
    p6.add_argument('-w', '--workdir', dest='workdir', help='working directory for sortmerna run')
    p6.add_argument('--threads', dest='threads', help='number of threads to use')
    p6.add_argument('--config', dest='config', help='Tests configuration file (overrides --test lookup)')
    p6.add_argument('--presets', dest='presets', help='Tests presets file')
    p6.add_argument('-d', '--dbg-level', dest='dbg_level', help='debug level 0 | 1 | 2')
    p6.add_argument('--task', dest='task', help='Processing task 0 | 1 | 2 | 3 | 4 | 5')
    p6.add_argument('--index', dest='index', help='Index option 0 | 1')

    args = p0.parse_args()

    if 'split' == args.cmd:
        res = split(files=args.file, num_splits=args.num_splits, rapidgz=args.rapidgz, workdir=args.workdir)
        ...
    elif 'iss_453' == args.cmd:
        res =  iss_453(num_splits=args.num_splits, reads_mln=args.reads_mln, threads=args.threads, workdir=args.workdir)
        ...
    elif 'test' == args.cmd:
        #cfg = process_config()  # process configuration
        outd = Path(args.workdir) / 'out' if args.workdir else None
        if outd and outd.exists():
            for ff in outd.iterdir():
                print(f'{ST} removing {ff}')
                ff.unlink()

        # resolve group aliases ('all', 'all-big') against presets.yaml
        # CLI dashes map to underscores in the yaml key (all-big -> all_big)
        script_dir = Path(__file__).parent
        presets_path = Path(args.presets) if getattr(args, 'presets', None) else script_dir / 'presets.yaml'
        _presets = {}
        if presets_path.exists():
            with open(presets_path) as _f:
                _presets = yaml.safe_load(_f) or {}

        ret = {}
        tlist = []
        for nm in args.name.replace(',', ' ').split():
            val = _presets.get(nm.replace('-', '_'))
            if isinstance(val, str):
                tlist.extend(x.strip() for x in val.split(',') if x.strip())
            else:
                tlist.append(nm)

        print(f'{ST} number of tests: {len(tlist)}')
        # Snapshot args once so preset values applied for test N don't leak into test N+1
        # (parse_test_config only fills None-valued attrs from the preset; without this
        # reset, the first test's preset would shadow every later test's preset).
        args_snapshot = vars(args).copy()
        results = []
        batch_start = time.time()
        for test in tlist:
            rcode = 0
            status = 'PASS'
            err_excerpt = ''
            t_start = time.time()
            for _k, _v in args_snapshot.items():
                setattr(args, _k, _v)
            args.name = test
            try:
                cfg = parse_test_config(args)
                # clean-up the KVDB, IDX directories, and the output.
                # May fail if any file in the directory is open. Close the files and re-run.
                if args.workdir and Path(args.workdir).exists() \
                    and not (args.validate_only or args.task in ['1','2']):
                    print(f'{ST} clearing workdir prior to {test}: {args.workdir}')
                    shutil.rmtree(args.workdir)

                if not args.validate_only:
                    is_capture = cfg.get('capture', False) or args.capture
                    cmd = [cfg.get('exe')]
                    for k,v in cfg['args'].items():
                        if isinstance(v, list):
                            for rr in v:
                                cmd.append(k)
                                cmd.append(rr)
                        else:
                            cmd.append(k)
                            if v:
                                cmd.append(v)
                    if getattr(args, 'interrupt', 0):
                        rcode, outl, errl = run_with_interrupts(
                            cmd,
                            num_interrupts=args.interrupt,
                            wait_after_align_start=getattr(args, 'interrupt_wait', 60))
                    else:
                        rcode, outl, errl = run_test(cmd, capture=is_capture)

                if args.task == '5':
                    msg = f'{ST} task 5 (only indexing) was requested for test {test}. Skipping validation.'
                    print(msg)
                elif rcode == 0 or not cfg.get('failonerror', True):
                    if args.func:
                        gdict = globals().copy()
                        gdict.update(locals())
                        func = gdict.get(args.func)
                        func(**cfg[test])
                    elif cfg.get('validate', {}).get('func'):
                        fn = cfg['validate']['func']
                        gdict = globals().copy()
                        gdict.update(locals())
                        func = gdict.get(fn)
                        if func:
                            func(TEST_DATA, ret, **cfg)
                    else:
                        process_output(**cfg)
                else:
                    raise RuntimeError(f'sortmerna exited with rcode={rcode}')
            except Exception as ex:
                status = 'FAIL'
                err_excerpt = f'{type(ex).__name__}: {ex}'.replace('\n', ' ')
                print(f'{ST} test {test} FAILED: {err_excerpt}')

            results.append(dict(name=test,
                                status=status,
                                duration=time.time() - t_start,
                                exit_code=rcode,
                                error=err_excerpt))
            if status == 'FAIL' and getattr(args, 'stop_on_fail', False):
                print(f'{ST} --stop-on-fail set; aborting batch after {test}')
                break

        # batch summary
        total_dur = time.time() - batch_start
        n_pass = sum(1 for r in results if r['status'] == 'PASS')
        n_fail = sum(1 for r in results if r['status'] == 'FAIL')
        col_name = max(6, max((len(r['name']) for r in results), default=6))
        hdr = f'{"name":<{col_name}}  {"status":<6}  {"dur_s":>8}  {"rcode":>5}  error'
        sep = '-' * len(hdr)
        summary_lines = [
            f'requested={len(tlist)} ran={len(results)} pass={n_pass} fail={n_fail} duration={total_dur:.1f}s',
            hdr,
            sep,
        ]
        for r in results:
            summary_lines.append(
                f'{r["name"]:<{col_name}}  {r["status"]:<6}  {r["duration"]:>8.1f}  {r["exit_code"]:>5}  {r["error"]}'
            )
        print()
        print(f'{ST} ===== BATCH SUMMARY =====')
        for ln in summary_lines:
            print(f'{ST} {ln}')

        _summary_cfg = _presets.get('summary')
        if _summary_cfg:
            summary_dir = Path(_summary_cfg)
        elif args.workdir:
            summary_dir = Path(args.workdir) / 'out'
        else:
            summary_dir = Path.cwd()
        summary_dir.mkdir(parents=True, exist_ok=True)
        summary_path = summary_dir / 'test_summary.txt'
        with open(summary_path, 'w') as _sf:
            _sf.write('\n'.join(summary_lines) + '\n')
        print(f'{ST} summary written to {summary_path}')

        if n_fail:
            sys.exit(n_fail)
    elif 'index-stats' == args.cmd:
        generate_index_stats(args)
    else:
        ...
#END main
#END END run.py
