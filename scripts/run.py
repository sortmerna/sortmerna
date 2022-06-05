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
file: run.py
created: Aug 12, 2019 Mon

conda install -c conda-forge scikit-bio

'''
import os
import sys
import subprocess
import platform
from optparse import OptionParser
import re
import time
import difflib
import shutil
import yaml
from jinja2 import Environment, FileSystemLoader
import pandas
import gzip

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
if   IS_WIN: OS = 'WIN'
elif IS_WSL: OS = 'WSL'
elif IS_LNX: OS = 'LNX'
else:
    print(f'Unable to define the platform: {pf}')
    sys.exit(1)

UHOME = os.environ.get('USERPROFILE') if IS_WIN else os.environ.get('HOME') # not defined for AWS SSM

SMR = 'sortmerna'
SMR_SRC  = None # source root dir
SMR_DIST = None # dist root dir
SMR_EXE  = None # full path to executable

ZLIB_SRC   = None
ZLIB_DIST  = None

ROCKS_SRC   = None
ROCKS_DIST  = None

# no binaries, so always build Release only
RAPID_SRC   = None
RAPID_DIST  = None

DATA_DIR = None
RUN_DIR  = None
OUT_DIR  = None
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
                logd['num_reads'][1] = int((re.split(' = ', line)[1]).strip())
            elif logd['results']['num_denovo'][0] in line:
                logd['results']['num_denovo'][1] = int((re.split(' = ', line)[1]).strip())
            elif logd['results']['num_hits'][0] in line:
                logd['results']['num_hits'][1] = int((re.split(' = | \(', line)[1]).strip())
            elif logd['results']['num_fail'][0] in line:
                logd['results']['num_fail'][1] = int((re.split(' = | \(', line)[1]).strip())
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
    :param dict        test configuration see 'test.jinja.yaml'
    '''
    ST = '[process_output]'
    global is_skbio
    log_struct = kwarg.get('aligned.log')
    logd = parse_log(LOGF, log_struct)
    vald = kwarg.get(name, {}).get('validate')
    #cmdd = kwarg.get('cmd')
    
    if not vald:
        print('{} Validation info not provided'.format(ST))
        return

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

def set_file_names(basenames, is_other=False):
    '''
    :param list basenames   list of basenames as in test.jinja,yaml:aligned.names
    :param bool is_other    flag aligned (default) | other files
    '''
    global ALI_NAMES


if __name__ == "__main__":
    '''
    python scripts/run.py -n t0 -v [--capture] [--validate-only]
    python scripts/run.py --name t12 -f process_otu --validate-only -d 2
    python scripts/run.py --name t0 -f dbg_blast --validate-only
    python /media/sf_a01_code/sortmerna/scripts/run.py --name t6 --envn LNX_VBox_Ubuntu_1804
    python /mnt/c/Users/biocodz/a01_code/sortmerna/tests/run.py --name t0 --winhome /mnt/c/Users/biocodz [--capture]  
                                                                              |_ on WSL
    '''
    ST = '[run.py:__main__]'
    import pdb; pdb.set_trace()
    is_opts_ok = True

    # process options
    optpar = OptionParser()
    optpar.add_option('-n', '--name', dest='name', help='Test to run e.g. t0 | t1 | t2 | to_lf | to_crlf | all')
    optpar.add_option('-f', '--func', dest='func', help='function to run: process_otu | ')
    optpar.add_option('-c', '--clean', action="store_true", help='clean Work directory and exit. Requires \'--name\'')
    optpar.add_option('-d', '--dbg_level', dest="dbg_level", help='debug level 0 | 1 | 2')
    optpar.add_option('--btype', dest='btype', default='release', help = 'Build type: release | debug')
    optpar.add_option('--pt_smr', dest='pt_smr', default='t1', help = 'Sortmerna Linkage type t1 | t2 | t3')
    optpar.add_option('--winhome', dest='winhome', help='when running on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    optpar.add_option('--capture', action="store_true", help='Capture output. By default prints to stdout')
    optpar.add_option('-v', '--validate-only', action="store_true", help='Only perform validation. Assumes aligement already done')
    optpar.add_option('--ddir', dest='ddir', help = 'Data directory')
    optpar.add_option('--config', dest='config', help='Tests configuration file.')
    optpar.add_option('--env', dest='envfile', help='Environment variables')
    optpar.add_option('-e','--envn', dest='envname', help=('Name of environment: WIN | WSL '
                                                  '| LNX_AWS | LNX_TRAVIS | LNX_VBox_Ubuntu_1804 | ..'))
    optpar.add_option('-w', '--workdir', dest='workdir', help='Environment variables')
    optpar.add_option('--threads', dest='threads', help='Number of threads to use')
    optpar.add_option('--index', dest='index', help='Index option 0 | 1 | 2')
    optpar.add_option('-t', '--task', dest='task', help='Processing task 0 | 1 | 2 | 3 | 4')

    (opts, args) = optpar.parse_args()

    # process configuration
    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print('Current dir: {}'.format(cur_dir))

    # check env.yaml. If no env file specified, try the current directory
    env_yaml = os.path.join(cur_dir, 'env.jinja') if not opts.envfile else opts.envfile
    if not os.path.exists(env_yaml):
        print('{} No environment config file found. Please, provide one using \'--env\' option'.format(ST))
        sys.exit(1)
    else:
        # load properties from env.yaml
        print('{} Using Environment configuration file: {}'.format(ST, env_yaml))
        #with open(env_yaml, 'r') as envh:
        #    env = yaml.load(envh, Loader=yaml.FullLoader)
        env_jj = Environment(loader=FileSystemLoader(os.path.dirname(env_yaml)), trim_blocks=True, lstrip_blocks=True)
        env_template = env_jj.get_template(os.path.basename(env_yaml))
    
        # render jinja template env.jinja.yaml
        vars = {'UHOME': UHOME}
        if IS_WSL: vars['WINHOME'] = opts.winhome
        env_str = env_template.render(vars)
        env = yaml.load(env_str, Loader=yaml.FullLoader)

    # check test.jinja.yaml
    cfgfile = os.path.join(cur_dir, 'test.jinja') if not opts.config else opts.config
    if not os.path.exists(cfgfile):
        print('{} No build configuration template found. Please, provide one using \'--config\' option'.format(ST))
        sys.exit(1)
    else:
        print('{} Using Build configuration template: {}'.format(ST, cfgfile))

    # load 'test.jinja' template
    jjenv = Environment(loader=FileSystemLoader(os.path.dirname(cfgfile)), trim_blocks=True, lstrip_blocks=True)
    template = jjenv.get_template(os.path.basename(cfgfile))

    if opts.envname:
        ENV = opts.envname
    elif IS_WIN or IS_WSL: 
        ENV = OS
        print('{} --envn was not specified - using {}'.format(ST, ENV))
    else:
        print('{} --envn is required on OS {}'.format(ST, OS))
        is_opts_ok = False

    # WRK_DIR priority:
    #   (1. opts.workdir, 2. env.jinja.yaml:WRK_DIR, 3. UHOME/sortmerna/run)
    if opts.workdir:
        WRK_DIR = opts.workdir
    else:
        val = env.get('WRK_DIR',{}).get(ENV)
        WRK_DIR = val if val else os.path.join(UHOME, 'sortmerna', 'run')

    # render 'test.jinja' template
    val = env.get(SMR,{}).get('src',{}).get(ENV)
    SMR_SRC  = val or '{}/sortmerna'.format(UHOME)
    if not os.path.exists(SMR_SRC):
        print(('{} Sortmerna source directory {} not found. '
            'Either specify location in env.jinja or '
            'make sure the sources exist at {}'.format(ST, SMR_SRC, SMR_SRC)))
    DATA_DIR = env['DATA_DIR'][ENV]
    vars = {'SMR_SRC':SMR_SRC, 'DATA_DIR':DATA_DIR, 'WRK_DIR':WRK_DIR}
    if opts.threads: 
        # prevent the renderer from interpreting the threads as int
        tmpl = '{}' if opts.threads[0] in ['\'','\"'] and opts.threads[-1] in ['\'','\"'] else '\'{}\''
        vars['THREADS'] = tmpl.format(opts.threads)
    if opts.index:
        tmpl = '{}' if opts.index[0] in ['\'','\"'] and opts.index[-1] in ['\'','\"'] else '\'{}\''
        vars['INDEX'] = tmpl.format(opts.index)
    if opts.task:
        tmpl = '{}' if opts.task[0] in ['\'','\"'] and opts.task[-1] in ['\'','\"'] else '\'{}\''
        vars['TASK'] = tmpl.format(opts.task)
    if opts.dbg_level:
        tmpl = '{}' if opts.dbg_level[0] in ['\'','\"'] and opts.index[-1] in ['\'','\"'] else '\'{}\''
        vars['DBG_LEVEL'] = tmpl.format(opts.dbg_level)
    cfg_str = template.render(vars)
    #cfg_str = template.render(env) # env[OS]
    cfg = yaml.load(cfg_str, Loader=yaml.FullLoader)
    
    val = env.get(SMR,{}).get('bin', {}).get(ENV)
    if val:
      SMR_EXE = os.path.join(val, 'sortmerna')
    else:
      val = env.get(SMR,{}).get('dist', {}).get(ENV)
      SMR_DIST = val if val else '{}/dist'.format(SMR_SRC)
      SMR_DIST = SMR_DIST + '/{}/{}'.format(opts.pt_smr, opts.btype) if IS_WIN else SMR_DIST
      SMR_EXE  = os.path.join(SMR_DIST, 'bin', 'sortmerna')

    if IS_WIN:
        file, ext = os.path.splitext(SMR_EXE)
        if not ext:
            SMR_EXE = '{}.exe'.format(SMR_EXE)

    if SMR_EXE and os.path.exists(SMR_EXE):
        print(f'{ST} using {SMR_EXE}')
    else:
        print('{ST} sortmerna executable {SMR_EXE} does not exist or not set')
        sys.exit(1)

    TEST_DATA = os.path.join(SMR_SRC, 'data')  # data directory
    ALI_NAMES = cfg.get('aligned.names')

    # run test
    ret = {}
    tlist = []
    if opts.name == 'all':
        excl = cfg.get('exclude.list', [])
        tlist = list(cfg.get('tests').keys())
        if excl:
            nexcl = len(excl)
            print(f'{ST} excluding {nexcl} tests from the list as specified in \'exclude.list\'')
        for tt in excl:
            tlist.remove(tt)
    else:
        tlist = [opts.name]
    print(f'{ST} number of tests: {len(tlist)}')
    for test in tlist:
        tn = cfg[test]['name']
        print('\n')
        print(f'{ST} running {test}: {tn}')
        process_smr_opts(cfg[test]['cmd'])

        # clean-up the KVDB, IDX directories, and the output. 
        # May Fail if any file in the directory is open. Close the files and re-run.
        if opts.clean:
            if os.path.exists(KVDB_DIR):
                print(f'{ST} removing KVDB: {KVDB_DIR}')
                shutil.rmtree(KVDB_DIR)
            if os.path.exists(IDX_DIR):
                print(f'{ST} removing Index: {IDX_DIR}')
                shutil.rmtree(IDX_DIR)
            break

        # clean previous alignments (KVDB)
        if os.path.exists(KVDB_DIR) and not opts.validate_only:
            print(f'{ST} Removing KVDB dir: {KVDB_DIR}')
            shutil.rmtree(KVDB_DIR)

        # clean output
        ali_dir = os.path.dirname(ALIF)
        if ali_dir and os.path.exists(ali_dir) and not opts.validate_only:
            print(f'{ST} Removing Aligned Output: {ali_dir}')
            shutil.rmtree(ali_dir)

        if OTHF:
            oth_dir = os.path.dirname(OTHF)
            if oth_dir and os.path.exists(oth_dir) and oth_dir != ali_dir and not opts.validate_only:
                print(f'{ST} removing Non-Aligned Output: {oth_dir}')
                shutil.rmtree(oth_dir)

        if not opts.validate_only:
            cfg[test]['cmd'].insert(0, SMR_EXE)
            is_capture = cfg[test].get('capture', False)
            ret = run(cfg[test]['cmd'], cwd=cfg[test].get('cwd'), capture=is_capture)

        # validate alignment results
        if ret.get('retcode', 0) == 0 or not cfg[test].get('failonerror', True):
            if opts.func:
                gdict = globals().copy()
                gdict.update(locals())
                func = gdict.get(opts.func)
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
#END main
#END END run.py