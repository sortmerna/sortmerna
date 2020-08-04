'''
FILE: run.py
Created: Aug 12, 2019 Mon
'''
import os
import sys
import subprocess
import platform
from optparse import OptionParser
import re
import skbio.io
import time
import difflib
import shutil
import yaml
from jinja2 import Environment, FileSystemLoader

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
    print('Unable to define the platform: {}'.format(pf))
    sys.exit(1)

UHOME = os.environ.get('USERPROFILE') if IS_WIN else os.environ.get('HOME')

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

KVDB_DIR = None
IDX_DIR  = None


def run(cmd, cwd=None, capture=False):
    '''
    '''
    STAMP = '[run]'
    ret = {'retcode':0, 'stdout':None, 'stderr':None}
    # print compiler version e.g. 'Microsoft (R) C/C++ Optimizing Compiler Version 19.16.27031.1 for x86'
    #"%VS_HOME%"\bin\Hostx86\x86\cl.exe
    if cwd:
        print('{} Running: {} in {}'.format(STAMP, ' '.join(cmd), cwd))
    else:
        print('{} Running: {}'.format(STAMP, ' '.join(cmd)))

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

    print("{} Run time: {}".format(STAMP, time.time() - start))
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
    STAMP = '[to_lf]'
    CRLF = b'\r\n'
    LF = b'\n'

    for dpath, dnames, fnames in os.walk(ddir):
        for fname in fnames:
            if fname.endswith('.fasta') or fname.endswith('.fastq'):
                fpath = os.path.join(dpath, fname)
                print('{} Converting {}'.format(STAMP, fpath))
                with open(fpath, 'rb') as infile:
                    content = infile.read()
                content = content.replace(CRLF, LF)
                with open(fpath, 'wb') as infile:
                    infile.write(content)
    print('to_lf Done')
#END to_lf

def parse_log(fpath):
    '''
    parse 'aligned.log' to dictionary

    @param fpath  'aligned.log' file
    '''
    logd = {
        'cmd': ['Command', None],
        'pid': ['Porcess pid', None],
        'params': {'refs': []},
        'num_reads': ['Total reads =', 0],
        'results': {
            'num_hits': ['Total reads passing E-value threshold', 0],
            'num_fail': ['Total reads failing E-value threshold', 0],
            'num_denovo': ['Total reads for de novo clustering', 0],
            'min_len': ['Minimum read length', 0],
            'max_len': ['Maximum read length', 0],
            'mean_len': ['Mean read length', 0]
        },
        'coverage': [],
        'num_id_cov': ['Total reads passing %%id and %%coverage thresholds', 0],
        'num_otus': ['Total OTUs', 0],
        'date': None
    }

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
                logd['num_id_cov'][1] = int((re.split(' = ', line)[1]).strip())
            elif logd['num_otus'][0] in line:
                logd['num_otus'][1] = int((re.split(' = ', line)[1]).strip())

    return logd
#END parse_log

def process_smr_opts(args):
    '''
    args  list of parameters passed to sortmerna
    '''
    STAMP = '[process_smr_opts]'
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

    READS_EXT = os.path.splitext(args[args.index('-reads')+1])[1]
    if READS_EXT in ['.gz']:
        READS_EXT = os.extsep + args[args.index('-reads')+1].split(os.extsep)[1]
    IS_FASTQ = 'fastq' == READS_EXT[1:]
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
        print('{} \'-workdir\' option was provided. Using workdir: [{}]'.format(STAMP, os.path.realpath(wdir)))
        ALIF = os.path.join(wdir, 'out', ALI_BASE + READS_EXT)
    else:
        ALIF = os.path.join(UHOME, 'sortmerna', 'run', 'out', ALI_BASE + READS_EXT)

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
            print('{} \'-workdir\' option was provided. Using workdir: [{}]'.format(STAMP, os.path.realpath(wdir)))
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

    LOGF    = os.path.join(os.path.dirname(ALIF), '{}.log'.format(ALI_BASE))
    BLASTF  = os.path.join(os.path.dirname(ALIF), '{}.blast'.format(ALI_BASE))
    OTUF    = os.path.join(os.path.dirname(ALIF), '{}_otus.txt'.format(ALI_BASE))
    DENOVOF = os.path.join(os.path.dirname(ALIF), '{}_denovo.fasta'.format(ALI_BASE))
    SAMF    = os.path.join(os.path.dirname(ALIF), '{}.sam'.format(ALI_BASE))
#END process_smr_opts

def process_output(**kwarg):
    '''
    Used by:
        test_simulated_amplicon_1_part_map
        test_simulated_amplicon_generic_buffer
        test_simulated_amplicon_12_part_index
    '''
    STAMP = '[{}]'.format(process_output)

    logd = parse_log(LOGF)
    vald = kwarg.get('validate')
    cmdd = kwarg.get('cmd')
    
    if not vald:
        print('{} Validation info not provided'.format(STAMP))
        return

    # Check number of reads
    if vald.get('num_reads'):
        tmpl = 'Testing num_reads: {}: {} Expected: {}'
        print(tmpl.format(LOG_BASE, logd['num_reads'][1], vald['num_reads']))
        assert vald['num_reads'] == logd['num_reads'][1]

    # Check reads count in aligned_fwd
    if vald.get('num_aligned_fwd'):
        num_fwd = 0
        if os.path.exists(ALI_FWD):
            if IS_FASTQ:
                for seq in skbio.io.read(ALI_FWD, format=READS_EXT[1:], variant=vald.get('variant')):
                    num_fwd += 1
            else:
                for seq in skbio.io.read(ALI_FWD, format=READS_EXT[1:]):
                    num_fwd += 1
            tmpl = 'Testing count of FWD aligned reads: {}: {} Expected: {}'
            print(tmpl.format('{}{}'.format(ALI_FWD_BASE, READS_EXT), num_fwd, vald['num_aligned_fwd']))
            assert num_fwd == vald['num_aligned_fwd'], \
                '{} not equals {}'.format(num_fwd, vald['num_aligned_fwd'])

    # Check reads count in aligned_rev
    if vald.get('num_aligned_rev'):
        num_rev = 0
        if os.path.exists(ALI_REV):
            if IS_FASTQ:
                for seq in skbio.io.read(ALI_REV, format=READS_EXT[1:], variant=vald.get('variant')):
                    num_rev += 1
            else:
                for seq in skbio.io.read(ALI_REV, format=READS_EXT[1:]):
                    num_rev += 1
            tmpl = 'Testing count of REV aligned reads: {}: {} Expected: {}'
            print(tmpl.format('{}{}'.format(ALI_REV_BASE, READS_EXT), num_rev, vald['num_aligned_rev']))
            assert num_rev == vald['num_aligned_rev'], \
                '{} not equals {}'.format(num_rev, vald['num_aligned_rev'])

    # Check reads count in other_fwd
    if vald.get('num_other_fwd'):
        num_fwd = 0
        if os.path.exists(OTH_FWD):
            if IS_FASTQ:
                for seq in skbio.io.read(OTH_FWD, format=READS_EXT[1:], variant=vald.get('variant')):
                    num_fwd += 1
            else:
                for seq in skbio.io.read(OTH_FWD, format=READS_EXT[1:]):
                    num_fwd += 1
            tmpl = 'Testing count of FWD non-aligned reads: {}: {} Expected: {}'
            print(tmpl.format('{}{}'.format(OTH_FWD_BASE, READS_EXT), num_fwd, vald['num_other_fwd']))
            assert num_fwd == vald['num_other_fwd'], \
                '{} not equals {}'.format(num_fwd, vald['num_other_fwd'])


    # Check reads count in other_rev
    if vald.get('num_other_rev'):
        num_rev = 0
        if os.path.exists(OTH_REV):
            if IS_FASTQ:
                for seq in skbio.io.read(OTH_REV, format=READS_EXT[1:], variant=vald.get('variant')):
                    num_rev += 1
            else:
                for seq in skbio.io.read(OTH_REV, format=READS_EXT[1:]):
                    num_rev += 1
            tmpl = 'Testing count REV non-aligned reads: {}: {} Expected: {}'
            print(tmpl.format('{}{}'.format(OTH_REV_BASE, READS_EXT), num_rev, vald['num_other_rev']))
            assert num_rev == vald['num_other_rev'], \
                '{} not equals {}'.format(num_rev, vald['num_other_rev'])

    # Correct number of de novo reads
    if vald.get('num_denovo'):
        assert vald['num_denovo'] == logd['results']['num_denovo'][1], \
            '{} not equals {}'.format(vald.get('num_denovo'), logd['results']['num_denovo'][1])

        num_denovo_file = 0
        for seq in skbio.io.read(DENOVOF, format='fasta'):
            num_denovo_file += 1

        assert logd['results']['num_denovo'][1] == num_denovo_file
    
    # Check number of reads mapped
    if vald.get('num_hits'):
        tmpl = 'Testing num_hits: {}: {} Expected: {}'
        print(tmpl.format(LOG_BASE, logd['results']['num_hits'][1], vald['num_hits']))
        assert vald['num_hits'] == logd['results']['num_hits'][1]
        num_hits_file = 0
        if os.path.exists(ALIF):
            if IS_FASTQ:
                for seq in skbio.io.read(ALIF, format=READS_EXT[1:], variant=vald.get('variant')):
                    num_hits_file += 1
            else:
                for seq in skbio.io.read(ALIF, format=READS_EXT[1:]):
                    num_hits_file += 1

            if not IS_PAIRED_IN and not IS_PAIRED_OUT:
                tmpl = 'Testing num_hits: {}{}: {} Expected: {}'
                print(tmpl.format(ALI_BASE, READS_EXT, num_hits_file, vald['num_hits']))
                assert logd['results']['num_hits'][1] == num_hits_file
            else:
                tmpl = 'Testing count of aligned reads: {}{}: {} Expected: {}'
                print(tmpl.format(ALI_BASE, READS_EXT, num_hits_file, vald['num_aligned']))
                assert num_hits_file == vald['num_aligned']

    # Check number of reads not mapped
    if vald.get('num_fail'):
        tmpl = 'Testing num_fail: {}: {} Expected: {}'
        print(tmpl.format(LOG_BASE, logd['results']['num_fail'][1], vald['num_fail']))
        assert vald['num_fail']  == logd['results']['num_fail'][1]
        num_fails_file = 0
        if OTHF and os.path.exists(OTHF):
            if os.stat(OTHF).st_size > 0:
                if IS_FASTQ:
                    for seq in skbio.io.read(OTHF, format=READS_EXT[1:], variant=vald.get('variant')):
                        num_fails_file += 1
                else:
                    for seq in skbio.io.read(OTHF, format=READS_EXT[1:]):
                        num_fails_file += 1

            if not IS_PAIRED_IN and not IS_PAIRED_OUT:
                tmpl = 'Testing num_fail: {}{}: {} Expected: {}'
                print(tmpl.format(OTH_BASE, READS_EXT, num_fails_file, vald['num_fail']))
                assert logd['results']['num_fail'][1] == num_fails_file
            else:
                tmpl = 'Testing count of non-aligned reads: {}{}: {} Expected: {}'
                print(tmpl.format(OTH_BASE, READS_EXT, num_fails_file, vald['num_other']))
                assert num_fails_file == vald['num_other']

    
    # Check count of reads passing %id and %coverage threshold
    # as given in alinged.log
    if vald.get('num_pass_id_cov'):
        tmpl = '{} Testing reads passing ID_COV: {}: {} Expected: {}'
        print(tmpl.format(STAMP, LOG_BASE, logd['num_id_cov'][1], vald['num_pass_id_cov']))
        assert vald['num_pass_id_cov'] == logd['num_id_cov'][1]

    # Check count of reads passing %id and %coverage threshold
    # as given in aligned.blast
    if vald.get('blast'):
        num_hits_file = 0
        num_pass_id = 0
        num_pass_cov = 0
        num_pass_id_cov_file = 0
        is_has_cov = False
        if os.path.exists(BLASTF):
            with open(BLASTF) as f_blast:
                for line in f_blast:
                    num_hits_file += 1
                    llist = line.strip().split('\t')
                    f_id = float(llist[2])
                    is_pass_id = f_id >= 97.0
                    if is_pass_id: 
                        num_pass_id += 1
                    is_has_cov = len(llist) > 12
                    if is_has_cov:
                        f_cov = float(line.strip().split('\t')[13])
                        is_pass_cov = f_cov >= 97.0
                        if is_pass_cov: 
                            num_pass_cov += 1
                        if is_pass_id and is_pass_cov:
                            num_pass_id_cov_file += 1
        
        BLAST_BASE = os.path.basename(BLASTF)
        tmpl = 'from {}: num_pass_id= {} num_pass_cov= {} num_pass_id_cov= {}'
        print(tmpl.format(BLAST_BASE, num_pass_id, num_pass_cov, num_pass_id_cov_file))
        
        if vald['blast'].get('num_pass_id'):
            tmpl = 'Testing reads passing ID threshold: {}: {} Expected: {}'
            print(tmpl.format(BLAST_BASE, num_pass_id, vald['blast']['num_pass_id']))
            assert num_pass_id == vald['blast']['num_pass_id'], \
                '{} not equals {}'.format(vald['blast']['num_pass_id'], num_pass_id)
        
        tmpl = 'Testing num_hits: {}: {} Expected: {}'
        print(tmpl.format(BLAST_BASE, num_hits_file, vald['blast']['num_recs']))
        assert num_hits_file == vald['blast']['num_recs'], \
            '{} not equals {}'.format(num_hits_file, vald['blast']['num_recs'])

        if is_has_cov:
            assert vald['blast']['num_pass_id_cov'] == num_pass_id_cov_file, \
                '{} not equals {}'.format(vald['blast']['num_pass_id_cov'], num_pass_id_cov_file)
    
    # Correct number of clusters recorded
    #self.assertEqual("4401", num_clusters_log) # 4400 before bug 52
    if logd.get('num_otus') and vald.get('num_cluster'):
        assert logd['num_otus'][1] in vald['num_cluster']  # 4400 for amplicon_12_part
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        if os.path.exists(OTUF):
            with open(OTUF) as f_otus:
                for line in f_otus:
                    num_clusters_file += 1
                    num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
            assert logd['num_otus'][1] == num_clusters_file
            assert logd['num_id_cov'][1] == num_reads_in_clusters_file
#END process_output

def t0(datad, ret={}, **kwarg):
    '''
    @param datad   Data directory
    @param outd    results output directory
    '''
    STAMP = '[t0:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))   

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
    print("{} Done".format(STAMP))
#END t0

def t2(datad, ret={}, **kwarg):
    '''
    @param datad   Data directory
    @param outd    results output directory
    @param kwargs  validation args

    Test the following case for alignment:
    beginning from align_ref_start = 0 and align_que_start = X, the read finishes
    before the end of the reference
              ref |----------------|
    que |------------------------|
                    LIS |-----|
                  ^
                  align_que_start
    '''
    STAMP = '[t2:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))

    vald = kwarg['validate']

    actual_alignment = []
    with open(BLASTF) as afile:
        for line in afile:
            actual_alignment = line.strip().split('\t')
        
    assert len(vald['expected']) == len(actual_alignment)
    assert sorted(vald['expected']) == sorted(actual_alignment)
    #a = set(expected_alignment) & set(actual_alignment)
    print("{} Done".format(STAMP))
#END t2

def t3(datad, ret={}, **kwarg):
    '''
    @param datad   Data directory
    @param outd    results output directory
    @param kwargs  validation args

    test_environmental_output

    Test outputting FASTA file for de novo clustering 
    using environmental data.

    Conditions: input FASTA file is processed in
                one mapped section.
    '''
    STAMP = '[t3:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))

    logd = parse_log(LOGF)
    vald = kwarg.get('validate')
    cmdd = kwarg.get('cmd')

    assert logd['results']['num_hits'][1] == vald['num_hits']

    # sort order (descending/ascending) of candidate references (alignment.cpp)
    is_refs_descending = False
    # 'Total OTUs' in aligned.log has to be equal 'num_groups' in test.jinja.yaml
    if is_refs_descending:
        assert logd['num_otus'][1] == vald['num_groups'][0] # originally
    else:
        assert logd['num_otus'][1] == vald['num_groups'][1], \
            '{}:num_otus = {} != num_groups = {} expected'.format(LOG_BASE, logd['num_otus'][1], vald['num_groups'][1])

    # OTU file contains one line per OTU group, so the number of lines
    # has to be equal 'Total OTUs' in aligned.log
    with open(OTUF) as f_otumap:
        num_clusters_file = sum(1 for line in f_otumap)
    assert logd['num_otus'][1] == num_clusters_file, \
        '{}:num_otus = {} != {}:num_otus = {}'.format(LOG_BASE, \
            logd['num_otus'][1], OTU_BASE, num_clusters_file)

    # number of reads in aligned_denovo.fasta has to be equal the
    # 'Total reads for de novo clustering' in aligned.log
    num_denovo_file = 0
    for seq in skbio.io.read(DENOVOF, format='fasta'):
        num_denovo_file += 1

    assert logd['results']['num_denovo'][1] == num_denovo_file, \
            '{}:num_denovo = {} != {}:num_denovo = {}'.format(LOG_BASE, \
                logd['results']['num_denovo'][1], DENOVO_BASE, num_denovo_file)
    
    print("{} Done".format(STAMP))
#END t3

def t4(datad, ret={}, **kwarg ):
    '''
    count idx files
    '''
    STAMP = '[t4:{}]'.format(kwarg.get('name'))
    vald = kwarg.get('validate')
    if IS_WIN:
        sfx = vald.get('idx_sfx_win')
    else:
        sfx = vald.get('idx_sfx_lin')
    idx_count = 0
    idx_count_expect = vald.get('num_idx') * 3 + 1
    if os.path.exists(IDX_DIR):
        idx_count = len([fn for fn in os.listdir(IDX_DIR) if str(sfx) in fn])

    print('{} Expected number of index files: {} Actual number: {}'.format(STAMP, idx_count_expect, idx_count))
    assert idx_count_expect == idx_count
#END t4

def t9(datad, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory

    test_output_all_alignments_f_rc
    '''
    STAMP = '[t9:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))
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
    
    print("{} Done".format(STAMP))
#END t9

def t10(datad, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param ret  dict  Results of the aignment run output

    test_ref_shorter_than_seed
    '''
    STAMP = '[t10:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))

    vald = kwarg.get('validate')

    if ret and ret.get('retcode'):
        assert ret['retcode'] == 1
        assert vald['err_msg'] in ret['stderr'].decode("utf-8")
   
    print("{} Done".format(STAMP))
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
    STAMP = '[t11:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))

    if ret and ret.get('retcode'):
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print("{} Done".format(STAMP))
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
    STAMP = '[t12:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))

    if ret and ret.get('retcode'):
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print("{} Done".format(STAMP))
#END t12

def t13(datad, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param capture Capture output
    
    Test sortmerna on simulated data,
        10000 reads with 1% error (--aligned),
        10000 reads with 10% error (de novo),
        10000 reads random (--other)

        Conditions: reference index processed
        as one unit and input query FASTA file
        in 6 sections.
    '''
    STAMP = '[t13:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))

    if ret and ret.get('retcode'):
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print("{} Done".format(STAMP))
#END t13

def t14(datad, ret={}, **kwarg):
    '''
    @param name
    @param datad   Data directory
    @param outd    results output directory
    @param ret   Dict

    identical to t13 except using 12-part index instead of 6-part
    
    Test sortmerna on simulated data,
        10000 reads with 1% error (--aligned),
        10000 reads with 10% error (de novo),
        10000 reads random (--other)

        Conditions: reference index processed
        as 12 parts and input query FASTA file
        in 1 section.
    '''
    STAMP = '[t14:{}]'.format(kwarg.get('name'))
    print('{} Validating ...'.format(STAMP))

    if ret and ret.get('retcode'):
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(**kwarg)
   
    print("{} Done".format(STAMP))
#END t14

def t17(datad, ret={}, **kwarg):
    '''
    @param name  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param capture Capture output
    '''
    STAMP = '[t17:{}]'.format(kwarg.get('name'))
    print('{} TODO: implement'.format(STAMP))
    logd = parse_log(LOGF)
    print("{} Done".format(STAMP))
#END t17

if __name__ == "__main__":
    '''
    python scripts/run.py --name t0 [--capture] [--env scripts/env_non_git.yaml] [--validate-only]
    python scripts/run.py --name t16 --env /home/xx/env.yaml
    python /mnt/c/Users/XX/sortmerna/tests/run.py --name t0 --winhome /mnt/c/Users/XX [--capture]
    '''
    STAMP = '[run.py:__main__]'
    #import pdb; pdb.set_trace()
    is_opts_ok = True

    # process options
    optpar = OptionParser()
    optpar.add_option('-n', '--name', dest='name', help='Test to run e.g. t0 | t1 | t2 | to_lf | to_crlf')
    optpar.add_option('-c', '--clean', action="store_true", help='clean build directory')
    optpar.add_option('--btype', dest='btype', default='release', help = 'Build type: release | debug')
    optpar.add_option('--pt_smr', dest='pt_smr', default='t1', help = 'Sortmerna Linkage type t1 | t2 | t3')
    optpar.add_option('--winhome', dest='winhome', help='when running on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    optpar.add_option('--capture', action="store_true", help='Capture output. By default prints to stdout')
    optpar.add_option('--validate-only', action="store_true", help='Only perform validation. Assumes aligement already done')
    optpar.add_option('--ddir', dest='ddir', help = 'Data directory')
    optpar.add_option('--config', dest='config', help='Tests configuration file.')
    optpar.add_option('--env', dest='envfile', help='Environment variables')
    optpar.add_option('--envn', dest='envname', help=('Name of environment: WIN | WSL '
                                                  '| LNX_AWS | LNX_TRAVIS | LNX_VBox_Ubuntu_1804 | ..'))
    optpar.add_option('--workdir', dest='workdir', help='Environment variables')

    (opts, args) = optpar.parse_args()

    # process configuration
    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print('Current dir: {}'.format(cur_dir))

    # check env.yaml. If no env file specified, try the current directory
    env_yaml = os.path.join(cur_dir, 'env.jinja.yaml') if not opts.envfile else opts.envfile
    if not os.path.exists(env_yaml):
        print('{} No environment config file found. Please, provide one using \'--env\' option'.format(STAMP))
        sys.exit(1)
    else:
        # load properties from env.yaml
        print('{} Using Environment configuration file: {}'.format(STAMP, env_yaml))
        #with open(env_yaml, 'r') as envh:
        #    env = yaml.load(envh, Loader=yaml.FullLoader)
        env_jj = Environment(loader=FileSystemLoader(os.path.dirname(env_yaml)), trim_blocks=True, lstrip_blocks=True)
        env_template = env_jj.get_template(os.path.basename(env_yaml))
    
        #   render jinja template
        env_str = env_template.render({'UHOME': UHOME})
        env = yaml.load(env_str, Loader=yaml.FullLoader)

    # check jinja.yaml
    cfgfile = os.path.join(cur_dir, 'test.jinja.yaml') if not opts.config else opts.config
    if not os.path.exists(cfgfile):
        print('{} No build configuration template found. Please, provide one using \'--config\' option'.format(STAMP))
        sys.exit(1)
    else:
        print('{} Using Build configuration template: {}'.format(STAMP, cfgfile))

    # load 'env.jinja.yaml' template
    jjenv = Environment(loader=FileSystemLoader(os.path.dirname(cfgfile)), trim_blocks=True, lstrip_blocks=True)
    template = jjenv.get_template(os.path.basename(cfgfile))

    if opts.envname:
        ENV = opts.envname
    elif IS_WIN or IS_WSL: 
        ENV = OS
        print('{} --envn was not specified - using {}'.format(STAMP, ENV))
    else:
        print('{} --envn is required on OS {}'.format(STAMP, OS))
        is_opts_ok = False

    # WRK_DIR priority:
    #   (1. opts.workdir, 2. env.jinja.yaml:WRK_DIR, 3. UHOME/sortmerna/run)
    if opts.workdir:
        WRK_DIR = opts.workdir
    else:
        val = env.get('WRK_DIR',{}).get(ENV)
        WRK_DIR = val if val else os.path.join(UHOME, 'sortmerna', 'run')

    # render 'test.jinja.yaml' template
    val = env.get(SMR,{}).get('src',{}).get(ENV)
    SMR_SRC  = val if val else '{}/sortmerna'.format(UHOME)
    if not os.path.exists(SMR_SRC):
        print(('{} Sortmerna source directory {} not found. '
            'Either specify location in env.jinja.yaml or '
            'make sure the sources exist at {}'.format(STAMP, SMR_SRC, SMR_SRC)))
    DATA_DIR = env['DATA_DIR'][ENV]
    cfg_str = template.render({'SMR_SRC':SMR_SRC, 'DATA_DIR':DATA_DIR, 'WRK_DIR':WRK_DIR})
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
        print('{} using {}'.format(STAMP, SMR_EXE))
    else:
        print('{} sortmerna executable {} does not exist or not set'.format(STAMP, SMR_EXE))
        sys.exit(1)

    TEST_DATA = os.path.join(SMR_SRC, 'data')

    process_smr_opts(cfg[opts.name]['cmd'])

    # clean-up the KVDB, IDX directories, and the output. 
    # May Fail if any file in the directory is open. Close the files and re-run.
    if opts.clean:
        if os.path.exists(KVDB_DIR):
            print('Removing KVDB dir: {}'.format(KVDB_DIR))
            shutil.rmtree(KVDB_DIR)
        if os.path.exists(IDX_DIR):
            print('Removing KVDB dir: {}'.format(IDX_DIR))
            shutil.rmtree(IDX_DIR)

    # clean previous alignments (KVDB)
    if os.path.exists(KVDB_DIR):
        print('{} Removing KVDB dir: {}'.format(STAMP, KVDB_DIR))
        shutil.rmtree(KVDB_DIR)

    # clean output
    ali_dir = os.path.dirname(ALIF)
    if ali_dir and os.path.exists(ali_dir) and not opts.validate_only:
        print('{} Removing Aligned Output: {}'.format(STAMP, ali_dir))
        shutil.rmtree(ali_dir)

    if OTHF:
        oth_dir = os.path.dirname(OTHF)
        if oth_dir and os.path.exists(oth_dir) and oth_dir != ali_dir and not opts.validate_only:
            print('{} Removing Non-Aligned Output: {}'.format(STAMP, oth_dir))
            shutil.rmtree(oth_dir)

    # run alignment
    ret = {}
    if not opts.validate_only:
        print('{} Running {}: {}'.format(STAMP, opts.name, cfg[opts.name]['name']))
        cfg[opts.name]['cmd'].insert(0, SMR_EXE)
        is_capture = cfg[opts.name].get('capture', False)
        ret = run(cfg[opts.name]['cmd'], cwd=cfg[opts.name].get('cwd'), capture=is_capture)

    # validate alignment results
    fn = cfg[opts.name].get('validate', {}).get('func')
    if fn:
        gdict = globals().copy()
        gdict.update(locals())
        func = gdict.get(fn)

        if func:
            func(TEST_DATA, ret, **cfg[opts.name])
    else:
        process_output(**cfg[opts.name])

    # tests
    #if funcs.get(opts.name):
        # validate alignment results
    #    if cfg[opts.name].get('validate'):
    #        nm = cfg[opts.name]['name']
    #        fn = cfg[opts.name]['validate']['func']
    #        funcs[fn](nm, TEST_DATA, OUT_DIR, ret, **cfg[opts.name]['validate'])
    # other funcs
    #elif opts.name == 'to_lf':
    #    to_lf(DATA_DIR)