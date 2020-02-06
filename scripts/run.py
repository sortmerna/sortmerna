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
IS_WIN = None
IS_WSL = None
IS_LNX = None

OS = None

UHOME = None

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

def t0(name, datad, outd, ret={}, **kwarg):
    '''
    @param datad   Data directory
    @param outd    results output directory
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))   

    BLAST_OUT = os.path.join(outd, 'aligned.blast')
    BLAST_EXPECTED = os.path.join(datad, 't0_expected_alignment.blast')

    dlist = []

    with open(BLAST_OUT, 'r') as fout:
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

def t0_other(name, datad, outd, ret={}, **kwarg):
    '''
    @param datad   Data directory
    @param outd    results output directory
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    print("{} Done".format(STAMP))
#END t0_other

def t1(name, datad, outd, ret={}, **kwarg):
    '''
    @param datad   Data directory
    @param outd    results output directory
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    print("{} Done".format(STAMP))
#END t1

def t2(name, datad, outd, ret={}, **kwarg):
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
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    OUTF = os.path.join(outd, 'aligned.blast')

    actual_alignment = []
    with open(OUTF) as afile:
        for line in afile:
            actual_alignment = line.strip().split('\t')
        
    assert len(kwarg['expected']) == len(actual_alignment)
    assert sorted(kwarg['expected']) == sorted(actual_alignment)
    #a = set(expected_alignment) & set(actual_alignment)
    print("{} Done".format(STAMP))
#END t2

def t3(name, datad, outd, ret={}, **kwarg):
    '''
    @param name  name of this test
    @param datad   Data directory
    @param outd    results output directory
    @param kwargs  validation args

    test_environmental_output

    Test outputting FASTA file for de novo clustering 
    using environmental data.

    Conditions: input FASTA file is processed in
                one mapped section.
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    LOGF = os.path.join(outd, 'aligned.log')
    OTUF = os.path.join(outd, 'aligned_otus.txt')
    DENOVOF = os.path.join(outd, 'aligned_denovo.fasta')

    num_hits_expect = kwarg.get('num_hits')
    num_groups_1_expect = kwarg.get('num_groups_1')
    num_groups_2_expect = kwarg.get('num_groups_2')

    num_hits = 0
    num_failures_log = 0
    num_clusters_log = 0

    with open(LOGF) as afile:
        for line in afile:
            if 'Total reads passing E-value threshold' in line:
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif 'Total reads for de novo clustering' in line:
                num_failures_log = (re.split('Total reads for de novo clustering = ', line)[1]).strip()
            elif 'Total OTUs' in line:
                num_clusters_log = (re.split('Total OTUs = ', line)[1]).strip()

    assert int(num_hits) == num_hits_expect

    # sort order (descending/ascending) of candidate references (alignment.cpp)
    is_refs_descending = False
    if is_refs_descending:
        num_groups = num_groups_1_expect # originally
    else:
        num_groups = num_groups_2_expect

    assert num_clusters_log == str(num_groups)

    with open(OTUF) as f_otumap:
        num_clusters_file = sum(1 for line in f_otumap)
    assert num_groups == num_clusters_file

    num_failures_file = 0
    for seq in skbio.io.read(DENOVOF, format='fasta'):
        num_failures_file += 1

    assert num_failures_log == str(num_failures_file)
    
    print("{} Done".format(STAMP))
#END t3

def t4(name, datad, outd, ret={}, **kwarg):
    '''
    @param name
    @param datad   Data directory
    @param outd    results output directory
    @param kwargs  validation args

    TODO: this was only an indexing test. Modify to perform alignment on multipart index.

    test_indexdb_split_databases
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    LOGF = os.path.join(outd, 'aligned.log')
    OTUF = os.path.join(outd, 'aligned_otus.txt')
    DENOVOF = os.path.join(outd, 'aligned_denovo.fasta')
    print('TODO: not yet implemented')
    print("{} Done".format(STAMP))
#END t4

def t5(name, datad, outd, ret={}, **kwarg):
    '''
    @param name
    @param datad   Data directory
    @param outd    results output directory
    @param kwargs  validation args

    test_mate_pairs, part 1
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    LOGF = os.path.join(outd, 'aligned.log')
    ALIF = os.path.join(outd, 'aligned.fastq')
    NONALIF = os.path.join(outd, 'other.fastq')

    num_hits_expect = kwarg.get('num_hits')
    num_fail_expect = kwarg.get('num_fail')

    tmpl = '{} num_hits_expect= {} num_fail_expect= {}'
    print(tmpl.format(STAMP, num_hits_expect, num_fail_expect))

    num_hits = None
    num_fail = None

    with open(LOGF) as afile:
        f_log_str = afile.read()
        assert "Total reads passing E-value threshold" in f_log_str
        assert "Total reads failing E-value threshold" in f_log_str

        afile.seek(0)
        for line in afile:
            if 'Total reads passing E-value threshold' in line:
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif 'Total reads failing E-value threshold' in line:
                num_fail = (re.split('Total reads failing E-value threshold = | \(', line)[1]).strip()

        tmpl = '{} num_hits= {} num_fail= {}'
        print(tmpl.format(STAMP, num_hits, num_fail))

    # Correct number of aligned reads
    with open(ALIF) as f_aligned:
        num_lines_aligned = sum(1 for line in f_aligned)

    print('{} num_lines_aligned= {}'.format(STAMP, num_lines_aligned))

    # Correct number of non-aligned reads
    with open(NONALIF) as f_nonaligned:
        num_lines_failed = sum(1 for line in f_nonaligned)

    print('{} num_lines_failed= {}'.format(STAMP, num_lines_failed))

    # Correct number of reads mapped
    assert num_hits_expect == int(num_hits)
    # Correct number of clusters recorded
    assert num_fail_expect == int(num_fail)
    
    assert num_hits_expect == num_lines_aligned/4
    assert num_fail_expect == num_lines_failed/4
    
    print("{} Done".format(STAMP))
#END t5

def t6(name, datad, outd, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory

    TODO: monitor bug 54

    test_mate_pairs, part 2, paired_in
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    LOGF = os.path.join(outd, 'aligned.log')
    ALIF = os.path.join(outd, 'aligned.fastq')
    NONALIF = os.path.join(outd, 'other.fastq')

    num_hits_expect = kwarg.get('num_hits')
    num_fail_expect = kwarg.get('num_fail')
    num_aligned_expect = kwarg.get('num_aligned')
    num_other_expect = kwarg.get('num_other')

    num_hits = None
    num_fail = None

    with open(LOGF) as afile:
        f_log_str = afile.read()
        assert "Total reads passing E-value threshold" in f_log_str
        assert "Total reads failing E-value threshold" in f_log_str

        afile.seek(0)
        for line in afile:
            if 'Total reads passing E-value threshold' in line:
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif 'Total reads failing E-value threshold' in line:
                num_fail = (re.split('Total reads failing E-value threshold = | \(', line)[1]).strip()

        # Correct number of reads mapped
        assert num_hits_expect == int(num_hits)
        # Correct number of clusters recorded
        assert num_fail_expect == int(num_fail)

        # Correct number of aligned reads
        with open(ALIF) as f_aligned:
            num_aligned_reads = sum(1 for line in f_aligned)
        assert num_aligned_expect == num_aligned_reads/4

        # Correct number of non-aligned reads
        with open(NONALIF) as f_nonaligned:
            num_nonaligned_reads = sum(1 for line in f_nonaligned)
        assert num_other_expect == num_nonaligned_reads
    
    print("{} Done".format(STAMP))
#END t6

def t7(name, datad, outd, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory

    TODO: monitor bug 54

    test_mate_pairs, part 3, paired_out
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    LOGF = os.path.join(outd, 'aligned.log')
    ALIF = os.path.join(outd, 'aligned.fastq')
    NONALIF = os.path.join(outd, 'other.fastq')

    num_hits_expect = kwarg.get('num_hits')
    num_fail_expect = kwarg.get('num_fail')
    num_aligned_expect = kwarg.get('num_aligned')
    num_other_expect = kwarg.get('num_other')

    num_hits = None
    num_fail = None

    with open(LOGF) as afile:
        f_log_str = afile.read()
        assert "Total reads passing E-value threshold" in f_log_str
        assert "Total reads failing E-value threshold" in f_log_str

        afile.seek(0)
        for line in afile:
            if 'Total reads passing E-value threshold' in line:
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif 'Total reads failing E-value threshold' in line:
                num_fail = (re.split('Total reads failing E-value threshold = | \(', line)[1]).strip()

        # Correct number of reads mapped
        assert num_hits_expect == int(num_hits)
        # Correct number of clusters recorded
        assert num_fail_expect == int(num_fail)

        # Correct number of aligned reads
        with open(ALIF) as f_aligned:
            num_aligned_reads = sum(1 for line in f_aligned)
        assert num_aligned_expect == num_aligned_reads/4

        # Correct number of non-aligned reads
        with open(NONALIF) as f_nonaligned:
            num_nonaligned_reads = sum(1 for line in f_nonaligned)
        assert num_other_expect == num_nonaligned_reads/4
    
    print("{} Done".format(STAMP))
#END t7

def t8(name, datad, outd, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory

    test_multiple_databases_search
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    LOGF = os.path.join(outd, 'aligned.log')
    ALIF = os.path.join(outd, 'aligned.fasta')

    num_reads_expect = kwarg.get('num_reads')
    num_hits_expect = kwarg.get('num_hits')

    total_reads_log = None
    num_hits_log = None

    with open(LOGF) as afile:
        f_log_str = afile.read()
        assert "Total reads passing E-value threshold" in f_log_str
        assert "Total reads failing E-value threshold" in f_log_str

        afile.seek(0)
        for line in afile:
            if 'Total reads =' in line:
                total_reads_log = (re.split(' = ', line)[1]).strip()
            elif 'Total reads passing E-value threshold' in line:
                num_hits_log = (re.split(' = | \(', line)[1]).strip()

        # Correct number of reads
        assert num_reads_expect == int(total_reads_log)
        # Correct number of clusters recorded
        assert num_hits_expect == int(num_hits_log)

        # Correct number of aligned reads
        num_hits_file = 0
        for seq in skbio.io.read(ALIF, format='fasta'):
            num_hits_file += 1
        assert num_hits_file == int(num_hits_log), \
            'num_hits_file = {} Not equals num_hits_log = {}'.format(num_hits_file, num_hits_log)
    
    print("{} Done".format(STAMP))
#END t8

def t9(name, datad, outd, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory

    test_output_all_alignments_f_rc
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    ALISAM = os.path.join(outd, 'aligned.sam')
        
    sam_alignments_expected = [
        [
            'GQ099317.1.1325_157_453_0:0:0_0:0:0_99/1',
            '0',
            'GQ099317.1.1325_157_453_0:0:0_0:0:0_99/1',
            '1',
            '255',
            '101M',
            '*',
            '0',
            '0',
            'GCTGGCACGGAGTTAGCCGGGGCTTATAAATGGTACCGTCATTGATTCTTCCCATTCTTTCGAAGTTTACATCCCGAGGGACTTCATCCTTCACGCGGCGT',
            '*',
            'AS:i:202',
            'NM:i:0'
        ],
        [
            'GQ099317.1.1325_157_453_0:0:0_0:0:0_99/1',
            '16',
            'GQ099317.1.1325_157_453_0:0:0_0:0:0_99/1',
            '102',
            '255',
            '101M',
            '*',
            '0',
            '0',
            'ACGCCGCGTGAAGGATGAAGTCCCTCGGGATGTAAACTTCGAAAGAATGGGAAGAATCAATGACGGTACCATTTATAAGCCCCGGCTAACTCCGTGCCAGC',
            '*',
            'AS:i:202',
            'NM:i:0'
        ]
    ]

    sam_alignments = []
    with open(ALISAM) as aligned_f:
        for line in aligned_f:
            if line.startswith('@'):
                continue
            alignment = line.strip().split("\t")
            sam_alignments.append(alignment)

    assert len(sam_alignments_expected) == len(sam_alignments)

    for alignment in sam_alignments_expected:
        assert alignment in sam_alignments
    
    print("{} Done".format(STAMP))
#END t9

def t10(name, datad, outd, ret={}, **kwarg):
    '''
    @param smrexe  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param ret  dict  Results of the aignment run output

    test_ref_shorter_than_seed
    '''
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    MSG = 'one of your sequences is shorter than the seed length 19'

    assert ret['retcode'] == 1
    #print('type(stderr): {}'.format(type(ret['stderr']))) # bytes
    assert MSG in ret['stderr'].decode("utf-8")
   
    print("{} Done".format(STAMP))
#END t10

def process_output(outd):
    '''
    Used by:
        test_simulated_amplicon_1_part_map
        test_simulated_amplicon_generic_buffer
        test_simulated_amplicon_12_part_index
    '''
    LOGF    = os.path.join(outd, 'aligned.log')
    ALIF    = os.path.join(outd, 'aligned.fasta')
    NONALIF = os.path.join(outd, 'other.fasta')
    DENOVOF = os.path.join(outd, 'aligned_denovo.fasta')
    OTUF    = os.path.join(outd, 'aligned_otus.txt')
    BLASTF  = os.path.join(outd, 'aligned.blast')

    with open(LOGF) as f_log:
        f_log_str = f_log.read()
        assert "Total reads passing E-value threshold" in f_log_str
        assert "Total reads for de novo clustering" in f_log_str
        assert "Total OTUs" in f_log_str
    
        num_pass_id_cov_log = ''
        
        f_log.seek(0)
        for line in f_log:
            if 'Total reads =' in line:
                total_reads_log = (re.split(' = ', line)[1]).strip()
            elif 'Total reads for de novo clustering' in line:
                num_denovo_log = (re.split(' = ', line)[1]).strip()
            elif 'Total reads passing E-value threshold' in line:
                num_hits_log = (re.split(' = | \(', line)[1]).strip()
            elif 'Total reads failing E-value threshold' in line:
                num_fails_log = (re.split(' = | \(', line)[1]).strip()
            elif 'Total reads passing' in line and 'id' in line and 'coverage thresholds' in line:
                num_pass_id_cov_log = (re.split(' = ', line)[1]).strip()
            elif 'Total OTUs' in line:
                num_clusters_log = (re.split('Total OTUs = ', line)[1]).strip()
    
    # Correct number of reads
    assert "30000" == total_reads_log
    
    # Correct number of de novo reads
    assert "9831" == num_denovo_log, '9831 not equals {}'.format(num_denovo_log)

    num_denovo_file = 0
    for seq in skbio.io.read(DENOVOF, format='fasta'):
        num_denovo_file += 1

    assert num_denovo_log == str(num_denovo_file)
    
    # Correct number of reads mapped
    assert "19995" == num_hits_log
    num_hits_file = 0
    for seq in skbio.io.read(ALIF, format='fasta'):
        num_hits_file += 1
    assert num_hits_log == str(num_hits_file)
    
    # Correct number of reads not mapped
    assert "10005" == num_fails_log
    num_fails_file = 0
    for seq in skbio.io.read(NONALIF, format='fasta'):
        num_fails_file += 1
    assert num_fails_log == str(num_fails_file)
    
    # Correct number of reads passing %id and %coverage threshold
    assert "10164" == num_pass_id_cov_log
    num_pass_id_cov_file = 0
    with open(BLASTF) as f_blast:
        for line in f_blast:
            f_id = float(line.strip().split('\t')[2])
            f_cov = float(line.strip().split('\t')[13])
            if (f_id >= 97.0 and f_cov >= 97.0):
                num_pass_id_cov_file += 1
    assert num_pass_id_cov_log == str(num_pass_id_cov_log)
    
    # Correct number of clusters recorded
    #self.assertEqual("4401", num_clusters_log) # 4400 before bug 52
    assert num_clusters_log in ['4400','4401'] # 4400 for amplicon_12_part
    num_clusters_file = 0
    num_reads_in_clusters_file = 0
    with open(OTUF) as f_otus:
        for line in f_otus:
            num_clusters_file += 1
            num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
    assert num_clusters_log == str(num_clusters_file)
    assert num_pass_id_cov_log == str(num_reads_in_clusters_file)
#END process_output

def t11(name, datad, outd, ret={}, **kwarg):
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
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    if ret['retcode']:
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(outd)
   
    print("{} Done".format(STAMP))
#END t11

def t12(name, datad, outd, ret={}, **kwarg):
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
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    if ret['retcode']:
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(outd)
   
    print("{} Done".format(STAMP))
#END t12

def t13(name, datad, outd, ret={}, **kwarg):
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
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    if ret['retcode']:
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(outd)
   
    print("{} Done".format(STAMP))
#END t13

def t14(name, datad, outd, ret={}, **kwarg):
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
    STAMP = '[{}]'.format(name)
    print('{} Validating ...'.format(STAMP))

    if ret['retcode']:
        print('ERROR running alignemnt. Return code: {}'.format(ret['retcode']))
        print(ret['stdout'])
        print(ret['stderr'])
        sys.exit(1)
    else:
        process_output(outd)
   
    print("{} Done".format(STAMP))
#END t14

# 'test_simulated_amplicon_generic_buffer'
# TODO: looks exactly the same as t11 + t12. Remove.

def t15(name, datad, outd, ret={}, **kwarg):
    '''
    @param name  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param capture Capture output
    '''
    STAMP = '[{}]'.format(name)
    print('{} Deleted test'.format(STAMP))
    print("{} Done".format(STAMP))
#END t15

def t16(name, datad, outd, ret={}, **kwarg):
    '''
    @param name  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param capture Capture output
    '''
    STAMP = '[{}]'.format(name)
    print('{} TODO: implement'.format(STAMP))
    print("{} Done".format(STAMP))
#END t16

def t17(name, datad, outd, ret={}, **kwarg):
    '''
    @param name  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param capture Capture output
    '''
    STAMP = '[{}]'.format(name)
    print('{} TODO: implement'.format(STAMP))
    print("{} Done".format(STAMP))
#END t17

def t18(name, datad, outd, ret={}, **kwarg):
    '''
    @param name  sortmerna.exe path
    @param datad   Data directory
    @param outd    results output directory
    @param capture Capture output
    '''
    STAMP = '[{}]'.format(name)
    print('{} TODO: implement'.format(STAMP))
    print("{} Done".format(STAMP))
#END t17

if __name__ == "__main__":
    '''
    python scripts/run.py --name t0 [--capture] [--env scripts/env_non_git.yaml]
    python scripts/run.py --name t16 --env /home/xx/env.yaml
    python /mnt/c/Users/XX/sortmerna/tests/run.py --name t0 --winhome /mnt/c/Users/XX [--capture]
    '''
    import pdb; pdb.set_trace()

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

    # process options
    optpar = OptionParser()
    optpar.add_option('-n', '--name', dest='name', help='Test to run e.g. t0 | t1 | t2 | to_lf | to_crlf')
    optpar.add_option('-c', '--clean', action="store_true", help='clean build directory')
    optpar.add_option('--btype', dest='btype', default='release', help = 'Build type: release | debug')
    optpar.add_option('--pt_smr', dest='pt_smr', default='t1', help = 'Sortmerna Linkage type t1 | t2 | t3')
    optpar.add_option('--winhome', dest='winhome', help='when running on WSL - home directory on Windows side e.g. /mnt/c/Users/XX')
    optpar.add_option('--capture', action="store_true", help='Capture output. By default prints to stdout')
    optpar.add_option('--ddir', dest='ddir', help = 'Data directory')
    optpar.add_option('--config', dest='config', help='Tests configuration file.')
    optpar.add_option('--env', dest='envfile', help='Environment variables')

    (opts, args) = optpar.parse_args()

    # process configuration
    cur_dir = os.path.dirname(os.path.realpath(__file__)) # directory where this script is located
    print('Current dir: {}'.format(cur_dir))

    # check env.yaml. If no env file specified, try the current directory
    env_yaml = os.path.join(cur_dir, 'env.yaml') if not opts.envfile else opts.envfile
    if not os.path.exists(env_yaml):
        print('No environment config file found. Please, provide one using \'--env\' option')
        sys.exit(1)
    else:
        # load properties from env.yaml
        print('Using Environment configuration file: {}'.format(env_yaml))
        with open(env_yaml, 'r') as envh:
            env = yaml.load(envh, Loader=yaml.FullLoader)

    # check jinja.yaml
    cfgfile = os.path.join(cur_dir, 'test.jinja.yaml') if not opts.config else opts.config
    if not os.path.exists(cfgfile):
        print('No build configuration template found. Please, provide one using \'--config\' option')
        sys.exit(1)
    else:
        print('Using Build configuration template: {}'.format(cfgfile))

    # load jinja template
    jjenv = Environment(loader=FileSystemLoader(os.path.dirname(cfgfile)), trim_blocks=True, lstrip_blocks=True)
    template = jjenv.get_template(os.path.basename(cfgfile))

    # render jinja template
    SMR_SRC  = env[OS][SMR]['src'] if env[OS][SMR]['src'] else '{}/sortmerna'.format(UHOME)
    if not os.path.exists(SMR_SRC):
        print('Sortmerna source directory {} not found. Either specify location in env.yaml or make sure the sources exist at {}'.format(SMR_SRC, SMR_SRC))
    DATA_DIR = env[OS]['DATA_DIR']
    cfg_str = template.render({'SMR_SRC':SMR_SRC, 'DATA_DIR':DATA_DIR})
    #cfg_str = template.render(env) # env[OS]
    cfg = yaml.load(cfg_str, Loader=yaml.FullLoader)
    
    SMR_DIST = env[OS][SMR]['dist'] if env[OS][SMR]['dist'] else '{}/dist'.format(SMR_SRC)
    SMR_DIST = SMR_DIST + '/{}/{}'.format(opts.pt_smr, opts.btype) if IS_WIN else SMR_DIST
    SMR_EXE  = os.path.join(SMR_DIST, 'bin', 'sortmerna') 
    RUN_DIR  = os.path.join(UHOME, 'sortmerna', 'run')
    TEST_DATA = os.path.join(SMR_SRC, 'data')

    if opts.name in ['t{}'.format(x) for x in range(0,18)]:
        OUT_DIR = os.path.join(RUN_DIR, 'out')

    # clean-up the run directory. May Fail if any file in the directory is open. Close the files and re-run.
    if opts.clean and os.path.isdir(RUN_DIR):
        print('Deleting: [{}]'.format(RUN_DIR))
        shutil.rmtree(RUN_DIR)

    #funcs =  {
    #    't0': t0,
    #    't0_other' : t0_other,
    #    't1': t1,
    #    't2': t2,
    #    't3': t3,
    #    't4': t4,
    #    't5': t5,
    #    't6': t6,
    #    't7': t7,
    #    't8': t8,
    #    't9': t9,
    #    't10': t10,
    #    't11': t11,
    #    't12': t12,
    #    't13': t13,
    #    't14': t14,
    #    't16': t16,
    #    't17': t17,
    #    't18': t18
    #}

    # id(cfg[opts.name]) # 2409596527368
    # id(cfg[opts.name].insert(0, SMR_EXE)) # 140725283474656

    # clean previous alignments (KVDB)
    kvdbdir = os.path.join(RUN_DIR, 'kvdb')
    if os.path.exists(kvdbdir):
        print('Removing KVDB dir: {}'.format(kvdbdir))
        shutil.rmtree(kvdbdir)
    # clean output
    if OUT_DIR and os.path.exists(OUT_DIR):
        print('Removing OUT_DIR: {}'.format(OUT_DIR))
        shutil.rmtree(OUT_DIR)
    # run alignment
    print('Running {}: {}'.format(opts.name, cfg[opts.name]['name']))
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
            nm = cfg[opts.name]['name']
            func(nm, TEST_DATA, OUT_DIR, ret, **cfg[opts.name]['validate'])

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