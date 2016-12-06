#!/usr/bin/env python
"""
Software tests for the SortMeRNA with ZLIB
==========================================
"""


from unittest import TestCase, main
import re
from subprocess import Popen, PIPE
from os.path import abspath, join, dirname
from tempfile import mkdtemp
from shutil import rmtree

from skbio.parse.sequences import parse_fasta


# ----------------------------------------------------------------------------
# Copyright (c) 2014--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Some notes on memory debugging:
# Memory check with valgrind --leak-check=full --track-origins=yes
# use -g compile option for line numbers in valgrind and traceback
# with gdb
# "$ export GLIBCXX_FORCE_NEW" to disable std::string memory pool optimizations
# prior to running valgrind

# Test class and cases
class SortmernaTestsZlib(TestCase):
    """ Tests for SortMeRNA functionality with ZLIB """

    def setUp(self):
        self.output_dir = mkdtemp()
        # 'data' folder must be in the same directory as test_sortmerna_zlib.py
        self.root = join(dirname(abspath(__file__)), "data")
        # reference databases
        self.db_bac16s = join(self.root, "silva-bac-16s-database-id85.fasta")
        # reads compressed with zip
        self.set8 = join(
            self.root, "set2_environmental_study_550_amplicon.fasta.zip")
        # reads compressed with gzip
        self.set9 = join(
            self.root, "set2_environmental_study_550_amplicon.fasta.gz")

    def tearDown(self):
        rmtree(self.output_dir)

    def output_test(self, aligned_basename):
        """ Test results of test_load_zip() and test_load_gzip()
        """
        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads for de novo clustering" in f_log_str)
        self.assertTrue("Total OTUs" in f_log_str)
        f_log.seek(0)
        for line in f_log:
            if line.startswith("    Total reads passing E-value threshold"):
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif line.startswith("    Total reads for de novo clustering"):
                num_failures_log = (re.split('Total reads for de novo clustering = ',
                              line)[1]).strip()
            elif line.startswith(" Total OTUs"):
                num_clusters_log = (re.split('Total OTUs = ', line)[1]).strip()
        f_log.close()
        # Correct number of reads mapped
        self.assertEqual("99999", num_hits)
        # Correct number of clusters recorded
        self.assertEqual("272", num_clusters_log)
        # Correct number of clusters in OTU-map
        with open(aligned_basename + "_otus.txt", 'U') as f_otumap:
            num_clusters_file = sum(1 for line in f_otumap)
        self.assertEqual(272, num_clusters_file)
        num_failures_file = 0
        with open(aligned_basename + "_denovo.fasta", 'U') as f_denovo:
            for label, seq in parse_fasta(f_denovo):
                num_failures_file += 1
        # Correct number of reads for de novo clustering
        self.assertEqual(num_failures_log, str(num_failures_file))

    def test_load_gzip(self):
        """ Load file compressed with gzip.
        """
        print "test_load_gzip"
        index_db = join(self.output_dir, "db_bac16s")
        index_path = "%s,%s" % (self.db_bac16s, index_db)
        indexdb_command = ["indexdb_rna",
                           "--ref",
                           index_path,
                           "--max_pos",
                           "250",
                           "-v"]
        proc = Popen(indexdb_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)
        proc.wait()
        aligned_basename = join(self.output_dir, "aligned")
        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--id", "0.97",
                             "--coverage", "0.97",
                             "--log",
                             "--otu_map",
                             "--de_novo_otu",
                             "--fastx",
                             "--reads-gz", self.set9,
                             "-v"]
        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)
        proc.wait()
        stdout, stderr = proc.communicate()
        if stderr:
            print stderr
        self.output_test(aligned_basename)

if __name__ == '__main__':
    main()
