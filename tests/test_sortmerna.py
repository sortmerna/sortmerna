#!/usr/bin/env python
"""
Software tests for the SortMeRNA
================================
"""


from unittest import TestCase, main
import re
from subprocess import Popen, PIPE
from os import close, walk, remove
from os.path import abspath, exists, getsize, join, dirname
from tempfile import mkstemp, mkdtemp
from shutil import rmtree

from skbio.util.misc import remove_files
from skbio.parse.sequences import parse_fasta


# ----------------------------------------------------------------------------
# Copyright (c) 2014--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


# Test class and cases
class SortmernaV2Tests(TestCase):
    """ Tests for SortMeRNA functionality """

    def setUp(self):
        self.output_dir = mkdtemp()
        self.root = "/Users/jenya/Desktop/sortmerna-dev/tests/data"

        # reference databases
        self.db_bac16s = join(self.root, "silva-bac-16s-database-id85.fasta")
        self.db_gg_13_8 = join(self.root, "gg_13_8_ref_set.fasta")

        # reads
        self.set2 = join(self.root, "set2_environmental_study_550_amplicon.fasta")
        self.set3 = join(self.root, "empty_file.fasta")
        self.set4 = join(self.root, "set4_mate_pairs_metatranscriptomics.fastq")
        self.set5 = join(self.root, "set5_simulated_amplicon_silva_bac_16s.fasta")

    def tearDown(self):
        rmtree(self.output_dir)

    def test_indexdb_default_param(self):
    	""" Test indexing a database using SortMeRNA
        """
        index_db = join(self.output_dir, "db_gg_13_8")
        index_path = "%s,%s" % (self.db_gg_13_8, index_db)

        indexdb_command = ["indexdb_rna",
                           "--ref",
                           index_path,
                           "-v"]

        proc = Popen(indexdb_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        self.assertTrue(stdout)
        self.assertFalse(stderr)

        expected_db_files = set(index_db + ext
                                for ext in ['.bursttrie_0.dat', '.kmer_0.dat',
                                            '.pos_0.dat', '.stats'])

        # Make sure all db_files exist
        for fp in expected_db_files:
            self.assertTrue(exists(fp))

    def test_indexdb_split_databases(self):
        """ Test indexing a database using SortMeRNA
            with m = 0.05, that is 7 parts
        """
        index_db = join(self.output_dir, "db_gg_13_8")
        index_path = "%s,%s" % (self.db_gg_13_8, index_db)

        indexdb_command = ["indexdb_rna",
                           "--ref",
                           index_path,
                           "-v",
                           "-m ",
                           "0.05"]

        proc = Popen(indexdb_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        self.assertFalse(stderr)
        self.assertTrue(stdout)

        expected_db_files = set(index_db + ext
                                for ext in ['.bursttrie_0.dat',
                                            '.bursttrie_1.dat',
                                            '.bursttrie_2.dat',
                                            '.bursttrie_3.dat',
                                            '.bursttrie_4.dat',
                                            '.bursttrie_5.dat',
                                            '.bursttrie_6.dat',
                                            '.kmer_0.dat',
                                            '.kmer_1.dat',
                                            '.kmer_2.dat',
                                            '.kmer_3.dat',
                                            '.kmer_4.dat',
                                            '.kmer_5.dat',
                                            '.kmer_6.dat',
                                            '.pos_0.dat',
                                            '.pos_1.dat',
                                            '.pos_2.dat',
                                            '.pos_3.dat',
                                            '.pos_4.dat',
                                            '.pos_5.dat',
                                            '.pos_6.dat',
                                            '.stats'])

        # Make sure all db_files exist
        for fp in expected_db_files:
            self.assertTrue(exists(fp))

    def test_simulated_amplicon_1_part_map(self):
        """ Test sortmerna on simulated data,
            10000 reads with 1% error (--aligned),
            10000 reads with 10% error (de novo),
            10000 reads random (--other)

            Conditions: reference index and input
            query FASTA file both processed as one
            section.
        """
        index_db = join(self.output_dir, "db_bac16s")
        index_path = "%s,%s" % (self.db_bac16s, index_db)

        indexdb_command = ["indexdb_rna",
                           "--ref",
                           index_path,
                           "-v"]

        proc = Popen(indexdb_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        aligned_basename = join(self.output_dir, "aligned")
        other_basename = join(self.output_dir, "other")

        # best 1
        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--other", other_basename,
                             "--reads", self.set5,
                             "--id", "0.97",
                             "--coverage", "0.97",
                             "--log",
                             "--otu_map",
                             "--de_novo_otu",
                             "--blast", "3",
                             "--fastx",
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads for de novo clustering" in f_log_str)
        self.assertTrue("Total OTUs" in f_log_str)
        f_log.seek(0)

        for line in f_log:
            if line.startswith("    Total reads = "):
                total_reads_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads for de novo clustering"):
                num_denovo_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads passing E-value threshold"):
                num_hits_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith("    Total reads failing E-value threshold"):
                num_fails_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith(" Total reads passing %id and %coverage thresholds"):
                num_pass_id_cov_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith(" Total OTUs"):
                num_clusters_log = (re.split('Total OTUs = ', line)[1]).strip()
        f_log.close()

        # Correct number of reads
        self.assertEqual("30000", total_reads_log)

        # Correct number of de novo reads
        self.assertEqual("10262", num_denovo_log)
        num_denovo_file = 0
        with open(aligned_basename + "_denovo.fasta", 'U') as f_denovo:
            for label, seq in parse_fasta(f_denovo):
                num_denovo_file +=1
        self.assertEqual(num_denovo_log, str(num_denovo_file))

        # Correct number of reads mapped
        self.assertEqual("19995", num_hits_log)
        num_hits_file = 0
        with open(aligned_basename + ".fasta", 'U') as f_mapped:
            for label, seq in parse_fasta(f_mapped):
                num_hits_file +=1
        self.assertEqual(num_hits_log, str(num_hits_file))

        # Correct number of reads not mapped
        self.assertEqual("10005", num_fails_log)
        num_fails_file = 0
        with open(other_basename + ".fasta", 'U') as f_not_mapped:
            for label, seq in parse_fasta(f_not_mapped):
                num_fails_file +=1
        self.assertEqual(num_fails_log, str(num_fails_file))

        # Correct number of reads passing %id and %coverage threshold
        self.assertEqual("9733", num_pass_id_cov_log)
        num_pass_id_cov_file = 0
        with open(aligned_basename + ".blast", 'U') as f_blast:
            for line in f_blast:
                f_id = float(line.strip().split('\t')[2])
                f_cov = float(line.strip().split('\t')[13])
                if (f_id >= 97.0 and f_cov >= 97.0):
                    num_pass_id_cov_file +=1
        self.assertEqual(num_pass_id_cov_log, str(num_pass_id_cov_log))

        # Correct number of clusters recorded
        self.assertEqual("4341", num_clusters_log)
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        with open(aligned_basename + "_otus.txt", 'U') as f_otus:
            for line in f_otus:
                num_clusters_file +=1
                num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
                
        self.assertEqual(num_clusters_log, str(num_clusters_file))
        self.assertEqual(num_pass_id_cov_log, str(num_reads_in_clusters_file))

        # Clean up before next call
        remove(aligned_basename + ".log")
        remove(aligned_basename + ".fasta")
        remove(aligned_basename + "_otus.txt")
        remove(aligned_basename + "_denovo.fasta")
        remove(aligned_basename + ".blast")
        remove(other_basename + ".fasta")

        # best 5
        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--other", other_basename,
                             "--reads", self.set5,
                             "--id", "0.97",
                             "--coverage", "0.97",
                             "--log",
                             "--otu_map",
                             "--de_novo_otu",
                             "--blast", "3",
                             "--fastx",
                             "--best", "5",
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads for de novo clustering" in f_log_str)
        self.assertTrue("Total OTUs" in f_log_str)
        f_log.seek(0)

        for line in f_log:
            if line.startswith("    Total reads = "):
                total_reads_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads for de novo clustering"):
                num_denovo_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads passing E-value threshold"):
                num_hits_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith("    Total reads failing E-value threshold"):
                num_fails_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith(" Total reads passing %id and %coverage thresholds"):
                num_pass_id_cov_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith(" Total OTUs"):
                num_clusters_log = (re.split('Total OTUs = ', line)[1]).strip()
        f_log.close()

        # Correct number of reads
        self.assertEqual("30000", total_reads_log)

        # Correct number of de novo reads
        self.assertEqual("10262", num_denovo_log)
        num_denovo_file = 0
        with open(aligned_basename + "_denovo.fasta", 'U') as f_denovo:
            for label, seq in parse_fasta(f_denovo):
                num_denovo_file +=1
        self.assertEqual(num_denovo_log, str(num_denovo_file))

        # Correct number of reads mapped
        self.assertEqual("19995", num_hits_log)
        num_hits_file = 0
        with open(aligned_basename + ".fasta", 'U') as f_mapped:
            for label, seq in parse_fasta(f_mapped):
                num_hits_file +=1
        self.assertEqual(num_hits_log, str(num_hits_file))

        # Correct number of reads not mapped
        self.assertEqual("10005", num_fails_log)
        num_fails_file = 0
        with open(other_basename + ".fasta", 'U') as f_not_mapped:
            for label, seq in parse_fasta(f_not_mapped):
                num_fails_file +=1
        self.assertEqual(num_fails_log, str(num_fails_file))

        # Correct number of reads passing %id and %coverage threshold
        self.assertEqual("9733", num_pass_id_cov_log)
        num_pass_id_cov_file = 0
        with open(aligned_basename + ".blast", 'U') as f_blast:
            for line in f_blast:
                f_id = float(line.strip().split('\t')[2])
                f_cov = float(line.strip().split('\t')[13])
                if (f_id >= 97.0 and f_cov >= 97.0):
                    num_pass_id_cov_file +=1
        self.assertEqual(num_pass_id_cov_log, str(num_pass_id_cov_log))

        # Correct number of clusters recorded
        self.assertEqual("4341", num_clusters_log)
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        with open(aligned_basename + "_otus.txt", 'U') as f_otus:
            for line in f_otus:
                num_clusters_file +=1
                num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
                
        self.assertEqual(num_clusters_log, str(num_clusters_file))
        self.assertEqual(num_pass_id_cov_log, str(num_reads_in_clusters_file))


    def test_simulated_amplicon_6_part_map(self):
        """ Test sortmerna on simulated data,
            10000 reads with 1% error (--aligned),
            10000 reads with 10% error (de novo),
            10000 reads random (--other)

            Conditions: reference index processed
            as one unit and input query FASTA file
            in 6 sections.
        """
        index_db = join(self.output_dir, "db_bac16s")
        index_path = "%s,%s" % (self.db_bac16s, index_db)

        indexdb_command = ["indexdb_rna",
                           "--ref",
                           index_path,
                           "-v"]

        proc = Popen(indexdb_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        aligned_basename = join(self.output_dir, "aligned")
        other_basename = join(self.output_dir, "other")

        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--other", other_basename,
                             "--reads", self.set5,
                             "--id", "0.97",
                             "--coverage", "0.97",
                             "--log",
                             "--otu_map",
                             "--de_novo_otu",
                             "--blast", "3",
                             "--fastx",
                             "-m", "1",
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads for de novo clustering" in f_log_str)
        self.assertTrue("Total OTUs" in f_log_str)
        f_log.seek(0)

        for line in f_log:
            if line.startswith("    Total reads = "):
                total_reads_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads for de novo clustering"):
                num_denovo_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads passing E-value threshold"):
                num_hits_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith("    Total reads failing E-value threshold"):
                num_fails_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith(" Total reads passing %id and %coverage thresholds"):
                num_pass_id_cov_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith(" Total OTUs"):
                num_clusters_log = (re.split('Total OTUs = ', line)[1]).strip()
        f_log.close()

        # Correct number of reads
        self.assertEqual("30000", total_reads_log)

        # Correct number of de novo reads
        self.assertEqual("10262", num_denovo_log)
        num_denovo_file = 0
        with open(aligned_basename + "_denovo.fasta", 'U') as f_denovo:
            for label, seq in parse_fasta(f_denovo):
                num_denovo_file +=1
        self.assertEqual(num_denovo_log, str(num_denovo_file))

        # Correct number of reads mapped
        self.assertEqual("19995", num_hits_log)
        num_hits_file = 0
        with open(aligned_basename + ".fasta", 'U') as f_mapped:
            for label, seq in parse_fasta(f_mapped):
                num_hits_file +=1
        self.assertEqual(num_hits_log, str(num_hits_file))

        # Correct number of reads not mapped
        self.assertEqual("10005", num_fails_log)
        num_fails_file = 0
        with open(other_basename + ".fasta", 'U') as f_not_mapped:
            for label, seq in parse_fasta(f_not_mapped):
                num_fails_file +=1
        self.assertEqual(num_fails_log, str(num_fails_file))

        # Correct number of reads passing %id and %coverage threshold
        self.assertEqual("9733", num_pass_id_cov_log)
        num_pass_id_cov_file = 0
        with open(aligned_basename + ".blast", 'U') as f_blast:
            for line in f_blast:
                f_id = float(line.strip().split('\t')[2])
                f_cov = float(line.strip().split('\t')[13])
                if (f_id >= 97.0 and f_cov >= 97.0):
                    num_pass_id_cov_file +=1
        self.assertEqual(num_pass_id_cov_log, str(num_pass_id_cov_log))

        # Correct number of clusters recorded
        self.assertEqual("4341", num_clusters_log)
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        with open(aligned_basename + "_otus.txt", 'U') as f_otus:
            for line in f_otus:
                num_clusters_file +=1
                num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
                
        self.assertEqual(num_clusters_log, str(num_clusters_file))
        self.assertEqual(num_pass_id_cov_log, str(num_reads_in_clusters_file))

    def test_simulated_amplicon_12_part_index(self):
        """ Test sortmerna on simulated data,
            10000 reads with 1% error (--aligned),
            10000 reads with 10% error (de novo),
            10000 reads random (--other)

            Conditions: reference index processed
            as 12 parts and input query FASTA file
            in 1 section.
        """
        index_db = join(self.output_dir, "db_bac16s")
        index_path = "%s,%s" % (self.db_bac16s, index_db)

        indexdb_command = ["indexdb_rna",
                           "--ref",
                           index_path,
                           "-v",
                           "-m", "10"]

        proc = Popen(indexdb_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        aligned_basename = join(self.output_dir, "aligned")
        other_basename = join(self.output_dir, "other")

        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--other", other_basename,
                             "--reads", self.set5,
                             "--id", "0.97",
                             "--coverage", "0.97",
                             "--log",
                             "--otu_map",
                             "--de_novo_otu",
                             "--blast", "3",
                             "--fastx",
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads for de novo clustering" in f_log_str)
        self.assertTrue("Total OTUs" in f_log_str)
        f_log.seek(0)

        for line in f_log:
            if line.startswith("    Total reads = "):
                total_reads_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads for de novo clustering"):
                num_denovo_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith("    Total reads passing E-value threshold"):
                num_hits_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith("    Total reads failing E-value threshold"):
                num_fails_log = (re.split(' = | \(', line)[1]).strip()
            elif line.startswith(" Total reads passing %id and %coverage thresholds"):
                num_pass_id_cov_log = (re.split(' = ', line)[1]).strip()
            elif line.startswith(" Total OTUs"):
                num_clusters_log = (re.split('Total OTUs = ', line)[1]).strip()
        f_log.close()

        # Correct number of reads
        self.assertEqual("30000", total_reads_log)

        # Correct number of de novo reads
        self.assertEqual("10262", num_denovo_log)
        num_denovo_file = 0
        with open(aligned_basename + "_denovo.fasta", 'U') as f_denovo:
            for label, seq in parse_fasta(f_denovo):
                num_denovo_file +=1
        self.assertEqual(num_denovo_log, str(num_denovo_file))

        # Correct number of reads mapped
        self.assertEqual("19995", num_hits_log)
        num_hits_file = 0
        with open(aligned_basename + ".fasta", 'U') as f_mapped:
            for label, seq in parse_fasta(f_mapped):
                num_hits_file +=1
        self.assertEqual(num_hits_log, str(num_hits_file))

        # Correct number of reads not mapped
        self.assertEqual("10005", num_fails_log)
        num_fails_file = 0
        with open(other_basename + ".fasta", 'U') as f_not_mapped:
            for label, seq in parse_fasta(f_not_mapped):
                num_fails_file +=1
        self.assertEqual(num_fails_log, str(num_fails_file))

        # Correct number of reads passing %id and %coverage threshold
        self.assertEqual("9733", num_pass_id_cov_log)
        num_pass_id_cov_file = 0
        with open(aligned_basename + ".blast", 'U') as f_blast:
            for line in f_blast:
                f_id = float(line.strip().split('\t')[2])
                f_cov = float(line.strip().split('\t')[13])
                if (f_id >= 97.0 and f_cov >= 97.0):
                    num_pass_id_cov_file +=1
        self.assertEqual(num_pass_id_cov_log, str(num_pass_id_cov_log))

        # Correct number of clusters recorded
        self.assertEqual("4343", num_clusters_log)
        num_clusters_file = 0
        num_reads_in_clusters_file = 0
        with open(aligned_basename + "_otus.txt", 'U') as f_otus:
            for line in f_otus:
                num_clusters_file +=1
                num_reads_in_clusters_file += (len(line.strip().split('\t'))-1)
                
        self.assertEqual(num_clusters_log, str(num_clusters_file))
        self.assertEqual(num_pass_id_cov_log, str(num_reads_in_clusters_file))

    def test_environmental_output(self):
        """ Test outputting FASTA file for de novo
            clustering using environmental data.

            Conditions: input FASTA file is processed in
            one mapped section.
        """
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
                             "--reads", self.set2,
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

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
                num_failures_log =\
                    (re.split('Total reads for de novo clustering = ',
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
                num_failures_file +=1

        # Correct number of reads for de novo clustering
        self.assertEqual(num_failures_log, str(num_failures_file))

    def test_empty_query_file(self):
        """ Test SortMeRNA with an empty reads file.
        """
        index_db = join(self.output_dir, "db_gg_13_8")
        index_path = "%s,%s" % (self.db_gg_13_8, index_db)

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
                             "--reads", self.set3,
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        # Correct number of clusters in OTU-map
        with open(aligned_basename + ".log", 'U') as f_log:
            self.assertTrue("The input reads file or reference file is empty" in f_log.read())

    def test_mate_pairs(self):
        """ Test outputting FASTQ files for merged
            mate pair reads.

            Conditions: input FASTQ file of mate paired reads:
            Total 10000 reads of which,
                1000 - align
                1000 - align
                2000 - random
                2000 - align
                2000 - align
                2000 - random

            Always only 6000 will align at any point, at the other 4000 are random reads.
            Using neither --paired_in or --paired_out, the --aligned file will have 6000 reads.
            With --paired_in, the --aligned file will contain 10000 reads
                              the --other file will contain 0 reads
            With --paired_out, the --aligned file will contain 2000 reads
                               the --other file will contain 8000 reads
        """
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
        nonaligned_basename = join(self.output_dir, "nonaligned")

        # launch normally
        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--other", nonaligned_basename,
                             "--fastx",
                             "--reads", self.set4,
                             "--log",
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads failing E-value threshold" in f_log_str)
        f_log.seek(0)

        for line in f_log:
            if line.startswith("    Total reads passing E-value threshold"):
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif line.startswith("    Total reads failing E-value threshold"):
                num_fail = (re.split('Total reads failing E-value threshold = | \(', line)[1]).strip()

        f_log.close()

        # Correct number of reads mapped
        self.assertEqual("6000", num_hits)

        # Correct number of clusters recorded
        self.assertEqual("4000", num_fail)

        # Correct number of aligned reads
        with open(aligned_basename + ".fastq", 'U') as f_aligned:
            num_aligned_reads = sum(1 for line in f_aligned)

        self.assertEqual(6000, num_aligned_reads/4)

        # Correct number of non-aligned reads
        with open(nonaligned_basename + ".fastq", 'U') as f_nonaligned:
            num_nonaligned_reads = sum(1 for line in f_nonaligned)

        self.assertEqual(4000, num_nonaligned_reads/4)

        # Clean up before next call
        remove(aligned_basename + ".log")
        remove(aligned_basename + ".fastq")
        remove(nonaligned_basename + ".fastq")

        # launch with option --paired_in
        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--other", nonaligned_basename,
                             "--paired_in",
                             "--fastx",
                             "--reads", self.set4,
                             "--log",
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads failing E-value threshold" in f_log_str)
        f_log.seek(0)

        for line in f_log:
            if line.startswith("    Total reads passing E-value threshold"):
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif line.startswith("    Total reads failing E-value threshold"):
                num_fail = (re.split('Total reads failing E-value threshold = | \(', line)[1]).strip()

        f_log.close()

        # Correct number of reads mapped
        self.assertEqual("6000", num_hits)

        # Correct number of clusters recorded
        self.assertEqual("4000", num_fail)

        # Correct number of aligned reads
        with open(aligned_basename + ".fastq", 'U') as f_aligned:
            num_aligned_reads = sum(1 for line in f_aligned)

        self.assertEqual(10000, num_aligned_reads/4)

        # Correct number of non-aligned reads
        with open(nonaligned_basename + ".fastq", 'U') as f_nonaligned:
            num_nonaligned_reads = sum(1 for line in f_nonaligned)

        self.assertEqual(0, num_nonaligned_reads/4)

        # Clean up before next call
        remove(aligned_basename + ".log")
        remove(aligned_basename + ".fastq")
        remove(nonaligned_basename + ".fastq")

        # launch with option --paired_out
        sortmerna_command = ["sortmerna",
                             "--ref", index_path,
                             "--aligned", aligned_basename,
                             "--other", nonaligned_basename,
                             "--paired_out",
                             "--fastx",
                             "--reads", self.set4,
                             "--log",
                             "-v"]

        proc = Popen(sortmerna_command,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)

        proc.wait()

        stdout, stderr = proc.communicate()

        print stderr

        f_log = open(aligned_basename + ".log", "U")
        f_log_str = f_log.read()
        self.assertTrue("Total reads passing E-value threshold" in f_log_str)
        self.assertTrue("Total reads failing E-value threshold" in f_log_str)
        f_log.seek(0)

        for line in f_log:
            if line.startswith("    Total reads passing E-value threshold"):
                num_hits = (re.split('Total reads passing E-value threshold = | \(', line)[1]).strip()
            elif line.startswith("    Total reads failing E-value threshold"):
                num_fail = (re.split('Total reads failing E-value threshold = | \(', line)[1]).strip()

        f_log.close()

        # Correct number of reads mapped
        self.assertEqual("6000", num_hits)

        # Correct number of clusters recorded
        self.assertEqual("4000", num_fail)

        # Correct number of aligned reads
        with open(aligned_basename + ".fastq", 'U') as f_aligned:
            num_aligned_reads = sum(1 for line in f_aligned)

        self.assertEqual(2000, num_aligned_reads/4)

        # Correct number of non-aligned reads
        with open(nonaligned_basename + ".fastq", 'U') as f_nonaligned:
            num_nonaligned_reads = sum(1 for line in f_nonaligned)

        self.assertEqual(8000, num_nonaligned_reads/4)


if __name__ == '__main__':
    main()
