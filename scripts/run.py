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
file: run.py
created: Aug 12, 2019 Mon

conda install scikit-bio -c conda-forge  <- pre-requisites
'''
import os
import sys
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
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

def mock_missing(name):
    def init(self, *args, **kwargs):
        raise ImportError(f'Failed to import class {name}; likely not installed.')
    return type(name, (), {'__init__': init})

try:
    import pandas
except ImportError:
    pandas = mock_missing('pandas')

is_skbio = True
try:
    import skbio.io
except (ImportError, OSError) as ex:
    print(f'\'import skbio.io\' failed: {ex}')
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

def run(cmd, cwd=None, capture=False):
    '''
    '''
    ST = '[run]'
    ret = {'retcode':0, 'stdout':None, 'stderr':None}
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
        if cwd and not os.path.exists(cwd):
            os.makedirs(cwd)
            proc = subprocess.run(cmd, cwd=cwd, capture_output=capture)
        else:
            proc = subprocess.run(cmd, capture_output=capture)

        ret['retcode'] = proc.returncode
        if capture:
            ret['stdout'] = proc.stdout
            ret['stderr'] = proc.stderr
        #proc = subprocess.run(cmd, cwd=build_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError as err:
        print(err)
        ret['retcode'] = 1
        ret['stderr'] = err
    except:
        for info in sys.exc_info(): print(info)
        ret['retcode'] = 1
        ret['stderr'] = sys.exc_info()

    rt = time.time() - start
    print(f"{ST} run time: {rt}")
    return ret
#END cmake_run

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

def parse_log(fpath, logd={}):
    '''
    parse 'aligned.log' to dictionary

    :param str fpath  'aligned.log' file
    :param dict logd  log file structure (test.jinja.yaml:aligned.log)
    '''
    #logd = {
    #    'cmd': ['Command', None],
    #    'pid': ['Process pid', None],
    #    'params': {'refs': []},
    #    'num_reads': ['Total reads =', 0],
    #    'results': {
    #        'num_hits': ['Total reads passing E-value threshold', 0],
    #        'num_fail': ['Total reads failing E-value threshold', 0],
    #        'num_denovo': ['Total reads for de novo clustering', 0],
    #        'min_len': ['Minimum read length', 0],
    #        'max_len': ['Maximum read length', 0],
    #        'mean_len': ['Mean read length', 0]
    #    },
    #    'coverage': [],
    #    'num_id_cov': ['Total reads passing %%id and %%coverage thresholds', 0],
    #    'num_otus': ['Total OTUs', 0],
    #    'date': None
    #}

    with open(fpath) as f_log:    
        for line in f_log:
            if logd['num_reads'][0] in line:
                logd['num_reads'][1] = int((re.split(r' = ', line)[1]).strip())
            elif logd['results']['num_denovo'][0] in line:
                logd['results']['num_denovo'][1] = int((re.split(r' = ', line)[1]).strip())
            elif logd['results']['num_hits'][0] in line:
                logd['results']['num_hits'][1] = int((re.split(r' = | \(', line)[1]).strip())
            elif logd['results']['num_fail'][0] in line:
                logd['results']['num_fail'][1] = int((re.split(r' = | \(', line)[1]).strip())
            elif logd['num_id_cov'][0] in line:
                # line: '..thresholds = 44223 (44.22)'
                val = line.split('=')[1]
                if val: val = val.split()[0].strip()
                logd['num_id_cov'][1] = int(val)
            elif logd['num_otus'][0] in line:
                logd['num_otus'][1] = int((re.split(' = ', line)[1]).strip())

    return logd
#END parse_log

def process_smr_opts(args):
    '''
    :param args  list of parameters passed to sortmerna
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
        print('{} \'-workdir\' option was provided. Using workdir: [{}]'.format(ST, os.path.realpath(wdir)))
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
            print('{} \'-workdir\' option was provided. Using workdir: [{}]'.format(ST, os.path.realpath(wdir)))
            ALIF = os.path.join(wdir, 'out', OTH_BASE + READS_EXT)
        else:
            OTHF = os.path.join(UHOME, 'sortmerna', 'run', 'out', OTH_BASE + READS_EXT)

        if OUT2 in args:
            OTH_FWD = os.path.join(os.path.dirname(ALIF), OTH_BASE + '_fwd' + READS_EXT)
            OTH_REV = os.path.join(os.path.dirname(ALIF), OTH_BASE + '_rev' + READS_EXT)

    if KVD in args:
        KVDB_DIR = args[args.index(KVD) + 1]
    elif WDIR in args:
        KVDB_DIR = os.path.join(args[args.index(WDIR) + 1], 'kvdb')
    else:
        KVDB_DIR = os.path.join(UHOME, 'sortmerna', 'run', 'kvdb')

    if IDX in args:
        IDX_DIR = args[args.index(IDX) + 1]
    elif WDIR in args:
        IDX_DIR = os.path.join(args[args.index(WDIR) + 1], 'idx')
    else:
        IDX_DIR = os.path.join(UHOME, 'sortmerna', 'run', 'idx')

    gzs = '.gz' if is_gz else ''
    LOGF    = os.path.join(os.path.dirname(ALIF), '{}.log'.format(ALI_BASE))
    BLASTF  = os.path.join(os.path.dirname(ALIF), '{}.blast{}'.format(ALI_BASE, gzs))
    OTUF    = os.path.join(os.path.dirname(ALIF), 'otu_map.txt')
    DENOVOF = os.path.join(os.path.dirname(ALIF), '{}_denovo.fa'.format(ALI_BASE))
    SAMF    = os.path.join(os.path.dirname(ALIF), '{}.sam'.format(ALI_BASE))
#END process_smr_opts

def process_blast(**kwarg):
    '''
    # Check count of reads passing %id and %coverage threshold
    # as given in aligned.blast
    '''
    ST = '[{}]'.format('process_blast')
    vald = kwarg.get('validate')

    BLAST_ID_COL = 2
    BLAST_COV_COL = 13
    num_hits_file = 0
    n_yid_ycov = 0
    n_yid_ncov = 0
    n_nid_ycov = 0
    n_denovo = 0
    has_cov = False
    BLAST_PID_PCOV = os.path.join(os.path.dirname(ALIF), 'pid_pcov.blast')
    is_dbg_pid_pcov = True

    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq', 'qcovs']
    
    if os.path.exists(BLASTF): 
        print(f'{ST} processing : {BLASTF}')
        gzs = os.path.basename(BLASTF).split(os.extsep)[-1]
        is_gz = gzs == 'gz'
        is_use_skbio = False
        with open(BLASTF, 'rb') as f_blast:
            is_gz = is_gz and (f_blast.read(2) == b'\x1f\x8b')
        if is_use_skbio:
            # just returns column names - not at all what's expected
            for seq in skbio.io.read(BLASTF, format='blast+6', columns=cols, compression='gzip', into=pandas.DataFrame):
                num_hits_file += 1
        else:
            with gzip.open(BLASTF, 'rb') if is_gz else open(BLASTF, 'rb') as f_blast, open(BLAST_PID_PCOV, 'w') as f_pid_pcov:
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
    
            bn = os.path.basename(BLASTF)
            print(f'{ST} from {bn}: num_hits= {num_hits_file} n_yid_ycov= {n_yid_ycov}'
                ' n_yid_ncov= {n_yid_ncov} n_nid_ycov= {n_nid_ycov} n_denovo= {n_denovo}')
            
            blastd = vald.get('files', {}).get('aligned.blast')
            if blastd:
                if blastd.get('n_yid_ycov'):
                    tmpl = '{} Testing reads passing ID threshold: {}: {} Expected: {}'
                    print(tmpl.format(ST, os.path.basename(BLASTF), n_yid_ycov, blastd['n_yid_ycov']))
                    assert n_yid_ycov == blastd['n_yid_ycov'], \
                        '{} not equals {}'.format(blastd['n_yid_ycov'], n_yid_ycov)
                
                num_recs = blastd.get('num_recs')
                if num_recs:
                    tmpl = '{} Testing num_hits: {}: {} Expected: {}'
                    print(tmpl.format(ST, os.path.basename(BLASTF), num_hits_file, num_recs))
                    assert num_hits_file == num_recs, \
                        '{} not equals {}'.format(num_hits_file, num_recs)

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
    ST = '[{}]'.format('dbg_blast')
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
                sr421.write('\n'.join('{}\t{}'.format(x[0], x[1]) for x in l421))
                sr431.write('\n'.join('{}\t{}'.format(x[0], x[1]) for x in l431))
            #l421 = [ '  '.join(line.strip().split('\t')[:2]) for line in f421 ]
            #l431 = [ '  '.join(line.strip().split('\t')[:2]) for line in f431 ]
            doset = False
            if doset:
                rdiff = set(l421) - set(l431)
                print('{} rdiff.len= {}'.format(ST, len(rdiff))) # 0
                [print(x) for x in list(rdiff)]
#END dbg_blast

def dbg_otu(**kwarg):
    '''
    '''
    ST = '[{}]'.format('dbg_otu')
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
    if os.path.exists(OTUF):
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
                readsf.write('{}\n'.format(read))

        with open(READS_DIFF, 'w') as diff:
            rds = set(blast_reads) - set(reads)
            rdsl = list(rds)
            rdsl.sort(key=lambda rr: int(rr.split('_')[-1]))
            rdsl.sort(key=lambda rr: int(rr.split('_')[0][:-1]))
            for rd in rdsl:
                diff.write('{}\n'.format(rd))
#END dbg_otu

def validate_otu(**kwarg):
    '''
    :param dict kwarg         test configuration see 'test.jinja.yaml'
    :param dict kwarg[logd]   parsed aligned.log data see 'parse_log(LOGF)'
    '''
    ST = '[{}]'.format('validate_otu')
    vald = kwarg.get('validate')
    logd = kwarg.get('logd') or parse_log(LOGF)

    # parse otu groups file
    if os.path.exists(OTUF):
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        with open(OTUF) as f_otus:
            for line in f_otus:
                num_clusters_file += 1
                num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
        print('{} num groups in OTU file {} , expected {}'.format(ST, num_clusters_file, logd['num_otus'][1]))
        assert logd['num_otus'][1] == num_clusters_file, \
            '{} not equals {}'.format(logd['num_otus'][1], num_clusters_file)
        print('{} count of reads in OTU file {} , expected {}'.format(ST, num_reads_in_clusters_file, logd['num_id_cov'][1]))
        assert logd['num_id_cov'][1] == num_reads_in_clusters_file, \
            '{} not equals {}'.format(logd['num_id_cov'][1], num_reads_in_clusters_file)

    # verify count of groups in the aligned.log is the same as specified in configuration validate data
    if logd.get('num_otus') and vald.get('num_cluster'):
        assert logd['num_otus'][1] in vald['num_cluster'], \
            '{} not in {}'.format(logd['num_otus'][1], vald['num_cluster'])
#END validate_otu

def validate_log(logd, ffd):
    '''
    :param dict logd
    :param dict ffd   files data as in test.jinja.yaml:<test_name>:validate:files
    '''
    ST = '[validate_log]'
    # aligned.log : 
    #   verify the total number of reads in aligned.log vs the number in the validation spec
    n_vald = ffd.get('aligned.log', {}).get('num_reads')
    n_logd = logd['num_reads'][1]
    if n_vald:
        print('{} testing num_reads: {} Expected: {}'.format(ST, n_logd, n_vald))
        assert n_vald == n_logd, '{} not equals {}'.format(n_vald, n_logd)
    #   verify number of hits
    n_vald = ffd.get('aligned.log', {}).get('num_hits')
    n_logd = logd['results']['num_hits'][1]
    if n_vald:
        print('{} testing num_hits: {} Expected: {}'.format(ST, n_logd, n_vald))
        assert n_vald == n_logd, '{} not equals {}'.format(n_vald, n_logd)
    #   verify number of misses
    n_vald = ffd.get('aligned.log', {}).get('num_fail')
    n_logd = logd['results']['num_fail'][1]
    if n_vald:
        print('{} testing num_fail: {} Expected: {}'.format(ST, n_logd, n_vald))
        assert n_vald == n_logd, '{} not equals {}'.format(n_vald, n_logd)
    #   verify count of COV+ID
    n_vald = ffd.get('aligned.log', {}).get('n_yid_ycov')
    n_logd = logd['num_id_cov'][1]
    if n_vald:
        print('{} testing n_yid_ycov: {} Expected: {}'.format(ST, n_logd, n_vald))
        assert n_vald == n_logd, '{} not equals {}'.format(n_vald, n_logd)
    #   verify count of OUT groups
    n_vald = ffd.get('aligned.log', {}).get('num_groups')
    n_logd = logd['num_otus'][1]
    if n_vald:
        print(f'{ST} testing num_groups: {n_logd} Expected: {n_vald}')
        assert n_vald == n_logd, '{} not equals {}'.format(n_vald, n_logd)
    #   verify count of de-novo reads
    n_vald = ffd.get('aligned.log', {}).get('n_denovo')
    n_logd = logd['results']['num_denovo'][1]
    if n_vald:
        print(f'{ST} testing n_denovo: {n_logd} Expected: {n_vald}')
        assert n_vald == n_logd, '{} not equals {}'.format(n_vald, n_logd)
#END validate_log

def process_output(name, **kwarg):
    '''
    :param str name    test name e.g. t0
    :param dict        test configuration see 'test.jinja'
    '''
    ST = '[process_output]'
    global is_skbio
    log_struct = kwarg.get('aligned.log')
    vald = kwarg.get(name, {}).get('validate')
    #cmdd = kwarg.get('cmd')
    
    if not vald:
        print('{} validation info not provided'.format(ST))
        return

    logd = parse_log(LOGF, log_struct)
    ffd = vald.get('files')
    if ffd and isinstance(ffd, dict):

        for ff, vv in ffd.items():
            print('{} {}'.format(ST, ff))
            # Check aligned/other reads count
            # aligned/other files only specify read count in the test validation data
            if isinstance(vv, int):
                ffp = os.path.join(os.path.dirname(ALIF), ff) # file path
                count = 0
                assert os.path.exists(ffp), '{} does not exists: {}'.format(ST, ffp)
                if ff == 'otu_map.txt':
                    with open(ffp) as ffs:
                        for line in ffs:
                            count += 1
                    print('{} testing count of groups in {}: {} Expected: {}'.format(ST, ff, count, vv))
                    assert count == vv, '{} not equals {}'.format(count, vv)
                    continue
                if is_skbio:
                    if IS_FASTQ:
                        for seq in skbio.io.read(ffp, format='fastq', variant=vald.get('variant')):
                            count += 1
                    else:
                        fmt = 'fasta' if READS_EXT[1:] in ['fasta', 'fa'] else READS_EXT[1:]
                        for seq in skbio.io.read(ffp, format=fmt):
                            count += 1
                    print('{} Testing count of reads in {}: {} Expected: {}'.format(ST, ff, count, vv))
                    assert count == vv, '{} not equals {}'.format(count, vv)
            elif ff == 'aligned.log':
                validate_log(logd, ffd)
            elif 'aligned.blast' in ff:
                test_cfg = kwarg.get(name, {})
                process_blast(**test_cfg)
#END process_output

def t0(datad, ret={}, **kwarg):
    '''
    :param datad   Data directory
    :param outd    results output directory
    '''
    ST = '[t0:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(ST))   

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
    print("{} Done".format(ST))
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
    ST = '[t2:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(ST))

    vald = kwarg['validate']

    actual_alignment = []
    with open(BLASTF) as afile:
        for line in afile:
            actual_alignment = line.strip().split('\t')
        
    assert len(vald['expected']) == len(actual_alignment)
    assert sorted(vald['expected']) == sorted(actual_alignment)
    #a = set(expected_alignment) & set(actual_alignment)
    print("{} Done".format(ST))
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
    ST = '[t3:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(ST))
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
            'num_otus = {} != num_groups = {} expected'.format(logd['num_otus'][1], vald['num_groups'][1])

    # OTU file contains one line per OTU group, so the number of lines
    # has to be equal 'Total OTUs' in aligned.log
    with open(OTUF) as f_otumap:
        num_clusters_file = sum(1 for line in f_otumap)
    assert logd['num_otus'][1] == num_clusters_file, \
        'num_otus = {} != {}:num_otus = {}'.format(logd['num_otus'][1], 
                                            OTU_BASE, num_clusters_file)

    # number of reads in aligned_denovo.fasta has to be equal the
    # 'Total reads for de novo clustering' in aligned.log
    n_denovo_file = 0
    if is_skbio:
        for seq in skbio.io.read(DENOVOF, format='fasta'):
            n_denovo_file += 1

        assert logd['results']['num_denovo'][1] == n_denovo_file, \
                'num_denovo = {} != {}:num_denovo = {}'.format(\
                    logd['results']['num_denovo'][1], DENOVO_BASE, n_denovo_file)
    
    print("{} Done".format(ST))
#END t3

def t4(datad, ret={}, **kwarg ):
    '''
    count idx files
    '''
    ST = '[t4:{}]'.format(kwarg.get('name'))
    vald = kwarg.get('validate')
    if IS_WIN:
        sfx = vald.get('idx_sfx_win')
    else:
        sfx = vald.get('idx_sfx_lin')
    idx_count = 0
    idx_count_expect = vald.get('num_idx') * 3 + 1
    if os.path.exists(IDX_DIR):
        idx_count = len([fn for fn in os.listdir(IDX_DIR) if str(sfx) in fn])

    print('{} Expected number of index files: {} Actual number: {}'.format(ST, idx_count_expect, idx_count))
    assert idx_count_expect == idx_count
#END t4

def t9(datad, ret={}, **kwarg):
    '''
    :param datad    data directory
    :param dict ret results of the aignment

    test_output_all_alignments_f_rc
    '''
    ST = '[t9:{}]'.format(kwarg.get('name'))
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
    ST = '[t10:{}]'.format(kwarg.get('name'))
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
    ST = '[t11:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(ST))

    if ret and ret.get('retcode'):
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print("{} Done".format(ST))
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
    ST = '[t12:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(ST))

    if ret and ret.get('retcode'):
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print("{} Done".format(ST))
#END t12

def t17(datad, ret={}, **kwarg):
    '''
    :param name  sortmerna.exe path
    :param datad   Data directory
    :param outd    results output directory
    :param capture Capture output
    '''
    ST = '[t17:{}]'.format(kwarg.get('name'))
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

def process_config() -> dict:
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
    import rapidgzip
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
    parser = ArgumentParser()
    parser.add_argument('--config', dest='config', help='Configuration file')
    parser.add_argument('--data-dir', dest='data_dir', help='path to the data. Abs or relative')
    parser.add_argument('--env', dest='envfile', help='Environment variables')
    subpar = parser.add_subparsers(dest='cmd', help='subcommand help')
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
    ptest = subpar.add_parser('test', help='run selected test')
    ptest.add_argument('name', help='Test to run e.g. t0 | t1 | t2 | to_lf | to_crlf | all')
    ptest.add_argument('--smr-exe', dest='smr_exe', help='path to sortmerna executable. Abs or relative')
    ptest.add_argument('--data-dir', dest='data_dir', help='path to the data. Abs or relative')
    ptest.add_argument('--threads', dest='threads', help='Number of threads to use')
    ptest.add_argument('--index', dest='index', help='Index option 0 | 1 | 2')
    ptest.add_argument('-t', '--task', dest='task', help='Processing task 0 | 1 | 2 | 3 | 4')
    ptest.add_argument('-d', '--dbg_level', dest="dbg_level", help='debug level 0 | 1 | 2')
    ptest.add_argument('-w', '--workdir', dest='workdir', help='Environment variables')
    ptest.add_argument('-c', '--clean', action="store_true", help='clean Work directory and exit. Requires \'--name\'')
    ptest.add_argument('-v', '--validate-only', action="store_true", help='Only perform validation. Assumes aligement already done')
    ptest.add_argument('-f', '--func', dest='func', help='function to run: process_otu | ')
    ptest.add_argument('--btype', dest='btype', default='release', help = 'Build type: release | debug')
    ptest.add_argument('--pt_smr', dest='pt_smr', default='t1', help = 'Sortmerna Linkage type t1 | t2 | t3')
    ptest.add_argument('--winhome', dest='winhome', help='when running on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    ptest.add_argument('--capture', action="store_true", help='Capture output. By default prints to stdout')
    ptest.add_argument('--config', dest='config', help='Tests configuration file.')
    ptest.add_argument('--env', dest='envfile', help='Environment variables')
    ptest.add_argument('-e','--envn', dest='envname', help=('Name of environment: WIN | WSL '
                                                      '| LNX_AWS | LNX_TRAVIS | LNX_VBox_Ubuntu_1804 | ..'))
    args = parser.parse_args()

    if 'split' == args.cmd:
        res = split(files=args.file, num_splits=args.num_splits, rapidgz=args.rapidgz, workdir=args.workdir)
        ...
    elif 'iss_453' == args.cmd:
        res =  iss_453(num_splits=args.num_splits, reads_mln=args.reads_mln, threads=args.threads, workdir=args.workdir)
        ...
    elif 'test' == args.cmd:
        if args.threads:
            # prevent the renderer from interpreting the threads as int
            tmpl = '{}' if args.threads[0] in ['\'','\"'] and args.threads[-1] in ['\'','\"'] else '\'{}\''
            vars['THREADS'] = tmpl.format(args.threads)
        if args.index:
            tmpl = '{}' if args.index[0] in ['\'','\"'] and args.index[-1] in ['\'','\"'] else '\'{}\''
            vars['INDEX'] = tmpl.format(args.index)
        if args.task:
            tmpl = '{}' if args.task[0] in ['\'','\"'] and args.task[-1] in ['\'','\"'] else '\'{}\''
            vars['TASK'] = tmpl.format(args.task)
        if args.dbg_level:
            tmpl = '{}' if args.dbg_level[0] in ['\'','\"'] and args.index[-1] in ['\'','\"'] else '\'{}\''
            vars['DBG_LEVEL'] = tmpl.format(args.dbg_level)

        cfg = process_config()  # process configuration
        
        # run test
        ret = {}
        tlist = []
        if args.name == 'all':
            excl = cfg.get('exclude.list', [])
            tlist = list(cfg.get('tests').keys())
            if excl:
                nexcl = len(excl)
                print(f'{ST} excluding {nexcl} tests from the list as specified in \'exclude.list\'')
            for tt in excl:
                if tt in tlist:
                    tlist.remove(tt)
        else:
            tlist = [args.name]
        print(f'{ST} number of tests: {len(tlist)}')
        for test in tlist:
            tn = cfg[test]['name']
            print('\n')
            print(f'{ST} running {test}: {tn}')
            process_smr_opts(cfg[test]['cmd'])

            # clean-up the KVDB, IDX directories, and the output. 
            # May Fail if any file in the directory is open. Close the files and re-run.
            if args.clean:
                if os.path.exists(KVDB_DIR):
                    print(f'{ST} removing dir: {KVDB_DIR}')
                    shutil.rmtree(KVDB_DIR)
                if os.path.exists(IDX_DIR):
                    print(f'{ST} removing dir: {IDX_DIR}')
                    shutil.rmtree(IDX_DIR)
                break

            # clean previous alignments (KVDB)
            if os.path.exists(KVDB_DIR) and not args.validate_only:
                print(f'{ST} Removing KVDB dir: {KVDB_DIR}')
                shutil.rmtree(KVDB_DIR)

            # clean output
            ali_dir = os.path.dirname(ALIF)
            if ali_dir and os.path.exists(ali_dir) and not args.validate_only:
                print(f'{ST} Removing Aligned Output: {ali_dir}')
                shutil.rmtree(ali_dir)

            if OTHF:
                oth_dir = os.path.dirname(OTHF)
                if oth_dir and os.path.exists(oth_dir) and oth_dir != ali_dir and not args.validate_only:
                    print(f'{ST} removing Non-Aligned Output: {oth_dir}')
                    shutil.rmtree(oth_dir)

            if not args.validate_only:
                cfg[test]['cmd'].insert(0, SMR_EXE)
                is_capture = cfg[test].get('capture', False) or args.capture
                ret = run(cfg[test]['cmd'], cwd=cfg[test].get('cwd'), capture=is_capture)

            # validate alignment results
            if ret.get('retcode', 0) == 0 or not cfg[test].get('failonerror', True):
                if args.func:
                    gdict = globals().copy()
                    gdict.update(locals())
                    func = gdict.get(args.func)
                    func(**cfg[test])
                elif cfg[test].get('validate', {}).get('func'):
                    fn = cfg[test].get('validate', {}).get('func')
                    gdict = globals().copy()
                    gdict.update(locals())
                    func = gdict.get(fn)
                    if func:
                        func(TEST_DATA, ret, **cfg[test])
                else:
                    process_output(test, **cfg)
            else:
                break
    else:
        ...
#END main
#END END run.py
