## SortMeRNA 3 Run applet (sortmerna-3.run.dev) execution trace

<pre>
* SortMeRNA 3 Run applet (sortmerna-3.run.dev:main) (done) job-FP8F3k80g35bF08fPj6x1Zj8
  biocodz 2018-10-24 13:54:45 (runtime 1:21:02)
2018-10-24 13:56:00 SortMeRNA 3 Run applet INFO Logging initialized (priority)
2018-10-24 13:56:01 SortMeRNA 3 Run applet INFO CPU: 10% (4 cores) * Memory: 283/7225MB * Storage: 79GB free * Net: 0↓/0↑MBps
2018-10-24 13:56:02 SortMeRNA 3 Run applet STDOUT dxpy/0.266.1 (Linux-4.4.0-98-generic-x86_64-with-Ubuntu-16.04-xenial)
2018-10-24 13:56:03 SortMeRNA 3 Run applet INFO Installing apt packages patchelf, zlib1g-dev, librocksdb-dev, rapidjson-dev, samtools
2018-10-24 13:56:10 SortMeRNA 3 Run applet INFO Setting SSH public key
2018-10-24 13:56:11 SortMeRNA 3 Run applet STDOUT dxpy/0.266.1 (Linux-4.4.0-98-generic-x86_64-with-Ubuntu-16.04-xenial)
2018-10-24 13:56:11 SortMeRNA 3 Run applet STDOUT /usr/sbin/rsyslogd already running.
2018-10-24 13:56:11 SortMeRNA 3 Run applet STDOUT /usr/sbin/sshd already running.
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT bash running (job ID job-FP8F3k80g35bF08fPj6x1Zj8)
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT BIN: {"$dnanexus_link": "file-FP7bJ0Q0qqg69XxGB7gk1f54"}
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT BINS: {"$dnanexus_link": "file-FP7bJ0Q0qqg69XxGB7gk1f54"} {"$dnanexus_link": "file-FP7bJ0j0qqg8Zz30Gzy142J4"} {"$dnanexus_link": "file-FP7bJ080qqgFJqP9FX0yf61J"}
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT BINS[0]: {"$dnanexus_link": "file-FP7bJ0Q0qqg69XxGB7gk1f54"}
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT BIN_name: indexdb
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT BIN_path: /home/dnanexus/in/BINS/0/indexdb
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT num BINs: 3
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT mkdir /home/dnanexus/bin/
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT mkdir /home/dnanexus/in/refs
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT mkdir /home/dnanexus/in/reads
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT mkdir /home/dnanexus/out
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT mkdir /home/dnanexus/idx
2018-10-24 13:56:13 SortMeRNA 3 Run applet STDOUT dx download -o /home/dnanexus/bin/ "{"$dnanexus_link": "file-FP7bJ0Q0qqg69XxGB7gk1f54"}"
2018-10-24 13:56:14 SortMeRNA 3 Run applet STDOUT dx download -o /home/dnanexus/bin/ "{"$dnanexus_link": "file-FP7bJ0j0qqg8Zz30Gzy142J4"}"
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT dx download -o /home/dnanexus/bin/ "{"$dnanexus_link": "file-FP7bJ080qqgFJqP9FX0yf61J"}"
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT ls -lrt /home/dnanexus/bin/
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT total 2572
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT -rwxr--r-- 1 root root  174560 Oct 24 17:56 indexdb
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root 1607120 Oct 24 17:56 libstdc++.so.6
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT -rwxr--r-- 1 root root  845832 Oct 24 17:56 sortmerna
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT apt-cache policy patchelf
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT patchelf:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Installed: 0.9-1~ubuntu16.04.1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Candidate: 0.9-1~ubuntu16.04.1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Version table:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT  *** 0.9-1~ubuntu16.04.1 500
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial-updates/universe amd64 Packages
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         100 /var/lib/dpkg/status
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT      0.8-4 500
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial/universe amd64 Packages
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT which patchelf
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT /usr/bin/patchelf
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT apt-cache policy librocksdb-dev
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT librocksdb-dev:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Installed: 4.1-1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Candidate: 4.1-1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Version table:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT  *** 4.1-1 500
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial/universe amd64 Packages
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         100 /var/lib/dpkg/status
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT apt-cache policy zlib1g-dev
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT zlib1g-dev:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Installed: 1:1.2.8.dfsg-2ubuntu4.1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Candidate: 1:1.2.8.dfsg-2ubuntu4.1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Version table:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT  *** 1:1.2.8.dfsg-2ubuntu4.1 500
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial-updates/main amd64 Packages
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         100 /var/lib/dpkg/status
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT      1:1.2.8.dfsg-2ubuntu4 500
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial/main amd64 Packages
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT apt-cache policy samtools
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT samtools:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Installed: 0.1.19-1ubuntu1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Candidate: 0.1.19-1ubuntu1
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT   Version table:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT  *** 0.1.19-1ubuntu1 500
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial/universe amd64 Packages
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT         100 /var/lib/dpkg/status
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT Running: patchelf --set-rpath '' /home/dnanexus/bin/sortmerna
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT REFS:            {"$dnanexus_link": "file-FJZz7300V8Pb4KZgBGkXbZBq"} {"$dnanexus_link": "file-FJZz7Fj0V8PX7G23K0Yf2Z2J"} {"$dnanexus_link": "file-FJZz78j0V8Pp0Y6Z07zK7K86"}
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT READS:
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT READS_GZ:        {"$dnanexus_link": "file-FKJ4ZpQ01FFKyQv8194pKv29"}
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT SAM:             true
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT FASTX:           true
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT BLAST:           1 cigar qcov qstrand
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT advanced:        -v
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT database_dir:    kvdb
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT processing_task: 4
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT help_smr:        false
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT dryidx:          false
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT drysmr:          false
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT Applet/Executable name: sortmerna-3.run.dev
2018-10-24 13:56:15 SortMeRNA 3 Run applet STDOUT Upload will be done into /out/sortmerna-3.run.dev
2018-10-24 13:56:16 SortMeRNA 3 Run applet STDOUT Upload path status:
2018-10-24 13:56:16 SortMeRNA 3 Run applet STDOUT dx download /home/dnanexus/in/refs/: "{"$dnanexus_link": "file-FJZz7300V8Pb4KZgBGkXbZBq"}"
2018-10-24 13:56:16 SortMeRNA 3 Run applet STDOUT dx download /home/dnanexus/in/refs/: "{"$dnanexus_link": "file-FJZz7Fj0V8PX7G23K0Yf2Z2J"}"
2018-10-24 13:56:17 SortMeRNA 3 Run applet STDOUT dx download /home/dnanexus/in/refs/: "{"$dnanexus_link": "file-FJZz78j0V8Pp0Y6Z07zK7K86"}"
2018-10-24 13:56:18 SortMeRNA 3 Run applet STDOUT download -o /home/dnanexus/in/reads/ "{"$dnanexus_link": "file-FKJ4ZpQ01FFKyQv8194pKv29"}"
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT reads_input basename: SRR1635864_1 extension: fastq
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT find /home/dnanexus -type f
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/.bash_logout
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/.bashrc
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/.profile
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/.byobu/.welcome-displayed
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/.bash_profile
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/job-FP8F3k80g35bF08fPj6x1Zj8
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/job_input.json
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/environment
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/dx_stdout
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/dx_stderr
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/dnanexus-job.json
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/dnanexus-executable.json
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/.ssh/authorized_keys
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/job-FP8F3k80g35bF08fPj6x1Zj8.code.sh
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/.bash_helper_vars
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/bin/indexdb
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/bin/libstdc++.so.6
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/bin/sortmerna
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/in/refs/silva-arc-16s-id95.fasta
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/in/refs/silva-bac-16s-id90.fasta
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/in/refs/silva-arc-23s-id98.fasta
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT /home/dnanexus/in/reads/SRR1635864_1.fastq.gz
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT Starting: indexdb --ref /home/dnanexus/in/refs/silva-arc-16s-id95.fasta,/home/dnanexus/idx/silva-arc-16s-id95:/home/dnanexus/in/refs/silva-bac-16s-id90.fasta,/home/dnanexus/idx/silva-bac-16s-id90:/home/dnanexus/in/refs/silva-arc-23s-id98.fasta,/home/dnanexus/idx/silva-arc-23s-id98 -v
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Program:      SortMeRNA version 3.0.0
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Copyright:    2016-2018 Clarity Genomics BVBA:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Turnhoutseweg 30, 2340 Beerse, Belgium
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 2014-2016 Knight Lab:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Department of Pediatrics, UCSD, La Jolla
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 2012-2014 Bonsai Bioinformatics Research Group:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 LIFL, University Lille 1, CNRS UMR 8022, INRIA Nord-Europe
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Disclaimer:   SortMeRNA comes with ABSOLUTELY NO WARRANTY; without even the
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 See the GNU Lesser General Public License for more details.
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Contributors: Jenya Kopylova   jenya.kopylov@gmail.com
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Laurent Noé      laurent.noe@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Pierre Pericard  pierre.pericard@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Daniel McDonald  wasade@gmail.com
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Mikaël Salson    mikael.salson@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Hélène Touzet    helene.touzet@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Rob Knight       robknight@ucsd.edu
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Parameters summary:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     K-mer size: 19
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     K-mer interval: 1
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     Maximum positions to store per unique K-mer: 10000
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Total number of databases to index: 3
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Begin indexing file /home/dnanexus/in/refs/silva-arc-16s-id95.fasta under index name /home/dnanexus/idx/silva-arc-16s-id95:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Collecting sequence distribution statistics ..  done  [0.035658 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   start index part # 0:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (1/3) building burst tries .. done  [1.347752 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (2/3) building CMPH hash .. done  [2.483723 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (3/3) building position lookup tables .. done [2.971647 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     total number of sequences in this part = 3193
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       temporary file was here: /tmp/sortmerna_keys_1781.txt
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing kmer data to /home/dnanexus/idx/silva-arc-16s-id95.kmer_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing burst tries to /home/dnanexus/idx/silva-arc-16s-id95.bursttrie_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing position lookup table to /home/dnanexus/idx/silva-arc-16s-id95.pos_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing nucleotide distribution statistics to /home/dnanexus/idx/silva-arc-16s-id95.stats
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     done.
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Begin indexing file /home/dnanexus/in/refs/silva-bac-16s-id90.fasta under index name /home/dnanexus/idx/silva-bac-16s-id90:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Collecting sequence distribution statistics ..  done  [0.179863 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   start index part # 0:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (1/3) building burst tries .. done  [10.760831 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (2/3) building CMPH hash .. done  [8.449476 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (3/3) building position lookup tables .. done [36.734165 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     total number of sequences in this part = 12798
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       temporary file was here: /tmp/sortmerna_keys_1781.txt
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing kmer data to /home/dnanexus/idx/silva-bac-16s-id90.kmer_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing burst tries to /home/dnanexus/idx/silva-bac-16s-id90.bursttrie_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing position lookup table to /home/dnanexus/idx/silva-bac-16s-id90.pos_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing nucleotide distribution statistics to /home/dnanexus/idx/silva-bac-16s-id90.stats
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     done.
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Begin indexing file /home/dnanexus/in/refs/silva-arc-23s-id98.fasta under index name /home/dnanexus/idx/silva-arc-23s-id98:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Collecting sequence distribution statistics ..  done  [0.006811 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   start index part # 0:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (1/3) building burst tries .. done  [0.209603 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (2/3) building CMPH hash .. done  [1.532718 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     (3/3) building position lookup tables .. done [0.324461 sec]
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     total number of sequences in this part = 251
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       temporary file was here: /tmp/sortmerna_keys_1781.txt
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing kmer data to /home/dnanexus/idx/silva-arc-23s-id98.kmer_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing burst tries to /home/dnanexus/idx/silva-arc-23s-id98.bursttrie_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing position lookup table to /home/dnanexus/idx/silva-arc-23s-id98.pos_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT       writing nucleotide distribution statistics to /home/dnanexus/idx/silva-arc-23s-id98.stats
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT     done.
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT Listing index: /home/dnanexus/idx/ ...
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT total 348612
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root   1048576 Oct 24 17:56 silva-arc-16s-id95.kmer_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root  26877892 Oct 24 17:56 silva-arc-16s-id95.bursttrie_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root  32665712 Oct 24 17:56 silva-arc-16s-id95.pos_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root     51230 Oct 24 17:56 silva-arc-16s-id95.stats
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root   1048576 Oct 24 17:57 silva-bac-16s-id90.kmer_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root 110426764 Oct 24 17:57 silva-bac-16s-id90.bursttrie_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root 165396916 Oct 24 17:57 silva-bac-16s-id90.pos_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root    204910 Oct 24 17:57 silva-bac-16s-id90.stats
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root   1048576 Oct 24 17:57 silva-arc-23s-id98.kmer_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root  10927568 Oct 24 17:57 silva-arc-23s-id98.bursttrie_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root   7258056 Oct 24 17:57 silva-arc-23s-id98.pos_0.dat
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root      4158 Oct 24 17:57 silva-arc-23s-id98.stats
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT Starting: sortmerna --ref /home/dnanexus/in/refs/silva-arc-16s-id95.fasta,/home/dnanexus/idx/silva-arc-16s-id95:/home/dnanexus/in/refs/silva-bac-16s-id90.fasta,/home/dnanexus/idx/silva-bac-16s-id90:/home/dnanexus/in/refs/silva-arc-23s-id98.fasta,/home/dnanexus/idx/silva-arc-23s-id98 --reads-gz /home/dnanexus/in/reads/SRR1635864_1.fastq.gz --aligned /home/dnanexus/out/SRR1635864_1_aligned --other /home/dnanexus/out/SRR1635864_1_other --sam --fastx --log --blast "1 cigar qcov qstrand" -v -d kvdb --task 4
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Program:      SortMeRNA version 3.0.0
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Copyright:    2016-2018 Clarity Genomics BVBA:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Turnhoutseweg 30, 2340 Beerse, Belgium
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 2014-2016 Knight Lab:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Department of Pediatrics, UCSD, La Jolla
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 2012-2014 Bonsai Bioinformatics Research Group:
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 LIFL, University Lille 1, CNRS UMR 8022, INRIA Nord-Europe
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Disclaimer:   SortMeRNA comes with ABSOLUTELY NO WARRANTY; without even the
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 See the GNU Lesser General Public License for more details.
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT   Contributors: Jenya Kopylova   jenya.kopylov@gmail.com
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Laurent Noé      laurent.noe@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Pierre Pericard  pierre.pericard@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Daniel McDonald  wasade@gmail.com
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Mikaël Salson    mikael.salson@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Hélène Touzet    helene.touzet@lifl.fr
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT                 Rob Knight       robknight@ucsd.edu
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT main:55 Running task ALIGN_REPORT: 4
2018-10-24 13:58:51 SortMeRNA 3 Run applet STDOUT Readstats::calculate starting ...   /home/dnanexus/sortmerna/src/sortmerna/gzip.cpp:119  WARNING: inflateEnd status is 0
2018-10-24 13:58:51 SortMeRNA 3 Run applet STDOUT Readstats::calculate done. Elapsed time: 61.15 sec. Reads processed: 24654594
2018-10-24 13:58:51 SortMeRNA 3 Run applet STDOUT align:374 Using default number of Processor threads equals num CPU cores: 4
2018-10-24 13:58:51 SortMeRNA 3 Run applet STDOUT Number of cores: 4 Read threads:  1 Write threads: 1 Processor threads: 4
2018-10-24 13:58:51 SortMeRNA 3 Run applet STDOUT ThreadPool:36 initialized Pool with: [6] threads
2018-10-24 13:58:51 SortMeRNA 3 Run applet STDOUT dirExists: Path does not exist: kvdb
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT read_queue created
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT write_queue created
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT Refstats:33 Index Statistics calculation Start ... Done. Time elapsed: 1.45 sec
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT align:414 Loading index 0 part 1/1 ... done [0.45] sec
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT align:425 Loading references  ... done [0.01] sec
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT Writer writer_0 thread 140612808595200 started
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612816987904 started
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT Processor proc_0 thread 140612800202496 started
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT Processor proc_1 thread 140612791809792 started
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT Processor proc_2 thread 140612783417088 started
2018-10-24 14:00:42 SortMeRNA 3 Run applet STDOUT Processor proc_3 thread 140612775024384 started
2018-10-24 14:06:01 SortMeRNA 3 Run applet INFO CPU: 53% (4 cores) * Memory: 497/7225MB * Storage: 77GB free * Net: 3↓/0↑MBps
2018-10-24 14:06:08 SortMeRNA 3 Run applet ERROR Error while relaying log message: Unterminated string starting at: line 1 column 9 (char 8)
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612816987904 done. Elapsed time: 434.61 sec Reads added: 24654594 readQueue.size: 0
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 3
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 2
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_3 thread 140612775024384 done. Processed 6165644 reads. Skipped already processed: 0 reads
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_1 thread 140612791809792 done. Processed 6151393 reads. Skipped already processed: 0 reads
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_0 thread 140612800202496 done. Processed 6172650 reads. Skipped already processed: 0 reads
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 4 jobs queue empty= 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 3 jobs queue empty= 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 5 jobs queue empty= 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 2 jobs queue empty= 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_2 thread 140612783417088 done. Processed 6164907 reads. Skipped already processed: 0 reads
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT writer_0 thread 140612808595200 done. Elapsed time: 434.61 s Reads written: 24654594
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 4
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT align:461 paralleltraversal: Done index 0 Part: 1 Time: 434.65 sec
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT align:414 Loading index 1 part 1/1 ... done [1.23] sec
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT align:425 Loading references  ... done [0.04] sec
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_0 thread 140612808595200 started
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_2 thread 140612800202496 started
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Writer writer_0 thread 140612775024384 started
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_1 thread 140612783417088 started
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT Processor proc_3 thread 140612816987904 started
2018-10-24 14:08:24 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612791809792 started
2018-10-24 14:16:01 SortMeRNA 3 Run applet INFO CPU: 82% (4 cores) * Memory: 1162/7225MB * Storage: 77GB free * Net: 0↓/0↑MBps
2018-10-24 14:20:48 SortMeRNA 3 Run applet ERROR Error while relaying log message: Unterminated string starting at: line 1 column 9 (char 8)
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612791809792 done. Elapsed time: 878.76 sec Reads added: 24654594 readQueue.size: 98
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 5 jobs queue empty= 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 3
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_2 thread 140612800202496 done. Processed 6181632 reads. Skipped already processed: 0 reads
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 2
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 4 jobs queue empty= 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_3 thread 140612816987904 done. Processed 6155509 reads. Skipped already processed: 0 reads
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 3 jobs queue empty= 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_1 thread 140612783417088 done. Processed 6122684 reads. Skipped already processed: 0 reads
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 2 jobs queue empty= 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_0 thread 140612808595200 done. Processed 6194769 reads. Skipped already processed: 0 reads
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT writer_0 thread 140612775024384 done. Elapsed time: 878.77 s Reads written: 24654594
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 4
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT align:461 paralleltraversal: Done index 1 Part: 1 Time: 878.98 sec
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT align:414 Loading index 2 part 1/1 ... done [0.24] sec
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT align:425 Loading references  ... done [0.00] sec
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Writer writer_0 thread 140612800202496 started
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612816987904 started
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_0 thread 140612808595200 started
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_1 thread 140612775024384 started
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_2 thread 140612791809792 started
2018-10-24 14:22:04 SortMeRNA 3 Run applet STDOUT Processor proc_3 thread 140612783417088 started
2018-10-24 14:26:01 SortMeRNA 3 Run applet INFO CPU: 74% (4 cores) * Memory: 1355/7225MB * Storage: 77GB free * Net: 0↓/0↑MBps
2018-10-24 14:27:53 SortMeRNA 3 Run applet ERROR Error while relaying log message: Unterminated string starting at: line 1 column 9 (char 8)
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612816987904 done. Elapsed time: 425.14 sec Reads added: 24654594 readQueue.size: 0
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 3
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Processor proc_0 thread 140612808595200 done. Processed 6064428 reads. Skipped already processed: 0 reads
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 5 jobs queue empty= 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 4 jobs queue empty= 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 2
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Processor proc_3 thread 140612783417088 done. Processed 6187218 reads. Skipped already processed: 0 reads
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 3 jobs queue empty= 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Processor proc_1 thread 140612775024384 done. Processed 6211835 reads. Skipped already processed: 0 reads
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 2 jobs queue empty= 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Processor proc_2 thread 140612791809792 done. Processed 6191113 reads. Skipped already processed: 0 reads
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT writer_0 thread 140612800202496 done. Elapsed time: 425.14 s Reads written: 24654594
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 4
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT align:461 paralleltraversal: Done index 2 Part: 1 Time: 425.16 sec
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Destructor called on write_queue  recs.size= 0 pushed: 73963782  popped: 73963782
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Destructor called on read_queue  recs.size= 0 pushed: 73963782  popped: 73963782
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Thread  140612816987904 job done
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Thread  140612808595200 job done
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Thread  140612783417088 job done
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Thread  140612775024384 job done
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Thread  140612791809792 job done
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT Thread  140612800202496 job done
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:27:53 SortMeRNA 3 Run applet STDOUT postProcess:175 Log file generation starts
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT ThreadPool:36 initialized Pool with: [3] threads
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT read_queue created
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT write_queue created
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT postProcess:184 Restored Readstats from DB: 1
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT Refstats:33 Index Statistics calculation Start ... Done. Time elapsed: 1.45 sec
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT postProcess:201: Loading reference 0 part 1/1  ... done [0.01 sec]
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT Writer writer_0 thread 140612791809792 started
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT PostProcessor postproc_0 thread 140612783417088 started
2018-10-24 14:31:24 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612775024384 started
/home/dnanexus/sortmerna/src/sortmerna/gzip.cpp:119  WARNING: inflateEnd status is 0
2018-10-24 14:36:01 SortMeRNA 3 Run applet INFO CPU: 26% (4 cores) * Memory: 904/7225MB * Storage: 77GB free * Net: 0↓/0↑MBps
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612775024384 done. Elapsed time: 419.13 sec Reads added: 24654594 readQueue.size: 1
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 2 jobs queue empty= 1
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT run:120 postproc_0 thread 140612783417088 done. Processed 24654594 reads
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT writer_0 thread 140612791809792 done. Elapsed time: 419.13 s Reads written: 0
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 1
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT Done reference 0 Part: 1 Time: 419.13 sec
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT postProcess:201: Loading reference 1 part 1/1  ... done [0.04 sec]
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT PostProcessor postproc_0 thread 140612791809792 started
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT Writer writer_0 thread 140612775024384 started
2018-10-24 14:37:57 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612783417088 started
/home/dnanexus/sortmerna/src/sortmerna/gzip.cpp:119  WARNING: inflateEnd status is 0
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612783417088 done. Elapsed time: 413.69 sec Reads added: 24654594 readQueue.size: 1
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 2 jobs queue empty= 1
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT writer_0 thread 140612775024384 done. Elapsed time: 413.69 s Reads written: 0
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT run:120 postproc_0 thread 140612791809792 done. Processed 24654594 reads
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 1
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT Done reference 1 Part: 1 Time: 413.70 sec
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT postProcess:201: Loading reference 2 part 1/1  ... done [0.00 sec]
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT Writer writer_0 thread 140612791809792 started
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT PostProcessor postproc_0 thread 140612783417088 started
2018-10-24 14:44:42 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612775024384 started
2018-10-24 14:46:01 SortMeRNA 3 Run applet INFO CPU: 20% (4 cores) * Memory: 906/7225MB * Storage: 77GB free * Net: 0↓/0↑MBps
/home/dnanexus/sortmerna/src/sortmerna/gzip.cpp:119  WARNING: inflateEnd status is 0
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612775024384 done. Elapsed time: 411.31 sec Reads added: 24654594 readQueue.size: 1
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 2 jobs queue empty= 1
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT run:120 postproc_0 thread 140612783417088 done. Processed 24654594 reads
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT writer_0 thread 140612791809792 done. Elapsed time: 411.31 s Reads written: 0
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: write_queue pushers: 0
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 1
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT Done reference 2 Part: 1 Time: 411.31 sec
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT readstats.total_reads_denovo_clustering: 0
2018-10-24 14:48:39 SortMeRNA 3 Run applet STDOUT postProcess:251 Done
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT Destructor called on write_queue  recs.size= 0 pushed: 0  popped: 0
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT Destructor called on read_queue  recs.size= 0 pushed: 73963782  popped: 73963782
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT Thread  140612783417088 job done
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT Thread  140612775024384 job done
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT Thread  140612791809792 job done
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT generateReports:810 Report generation starts. Thread: 140612850653440
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT ThreadPool:36 initialized Pool with: [2] threads
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT generateReports:818 Restored Readstats from DB: 1
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT read_queue created
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT write_queue created
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT Refstats:33 Index Statistics calculation Start ... Done. Time elapsed: 1.44 sec
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT generateReports:836 Loading reference 0 part 1/1  ... done [0.01 sec]
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT Report Processor report_proc_0 thread 140612791809792 started
2018-10-24 14:51:24 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612783417088 started
/home/dnanexus/sortmerna/src/sortmerna/gzip.cpp:119  WARNING: inflateEnd status is 0
2018-10-24 14:56:01 SortMeRNA 3 Run applet INFO CPU: 14% (4 cores) * Memory: 774/7225MB * Storage: 71GB free * Net: 0↓/0↑MBps
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612783417088 done. Elapsed time: 376.36 sec Reads added: 24654594 readQueue.size: 98
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT Report Processor report_proc_0 thread 140612791809792 done. Processed 24654594 reads
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 1
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT generateReports:864 Done reference 0 Part: 1 Time: 376.36 sec
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT generateReports:836 Loading reference 1 part 1/1  ... done [0.05 sec]
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT Report Processor report_proc_0 thread 140612791809792 started
2018-10-24 14:57:54 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612783417088 started
/home/dnanexus/sortmerna/src/sortmerna/gzip.cpp:119  WARNING: inflateEnd status is 0
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612783417088 done. Elapsed time: 387.42 sec Reads added: 24654594 readQueue.size: 1
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT Report Processor report_proc_0 thread 140612791809792 done. Processed 24654594 reads
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 1
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT generateReports:864 Done reference 1 Part: 1 Time: 387.42 sec
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT generateReports:836 Loading reference 2 part 1/1  ... done [0.00 sec]
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT Report Processor report_proc_0 thread 140612783417088 started
2018-10-24 15:04:19 SortMeRNA 3 Run applet STDOUT reader_0 thread: 140612791809792 started
2018-10-24 15:06:01 SortMeRNA 3 Run applet INFO CPU: 20% (4 cores) * Memory: 766/7225MB * Storage: 71GB free * Net: 0↓/0↑MBps
/home/dnanexus/sortmerna/src/sortmerna/gzip.cpp:119  WARNING: inflateEnd status is 0
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT read:135 reader_0 thread: 140612791809792 done. Elapsed time: 383.46 sec Reads added: 24654594 readQueue.size: 1
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT Report Processor report_proc_0 thread 140612783417088 done. Processed 24654594 reads
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 1 jobs queue empty= 1
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT threadEntry:108 number of running_threads= 0 jobs queue empty= 1
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT reset:131 write_queue: pushers: 1
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT decrPushers:161 id: read_queue pushers: 0
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT reset:131 read_queue: pushers: 1
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT generateReports:864 Done reference 2 Part: 1 Time: 383.46 sec
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT generateReports:871 Done Reports generation
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT Destructor called on write_queue  recs.size= 0 pushed: 0  popped: 0
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT Destructor called on read_queue  recs.size= 0 pushed: 73963782  popped: 73963782
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT Thread  140612791809792 job done
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT Thread  140612783417088 job done
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT Output.closefiles called. Flushed and closed
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT [INFO] ls -lrt /home/dnanexus/out
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT total 6020152
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root 5847303983 Oct 24 18:54 SRR1635864_1_other.fastq
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root  128270347 Oct 24 18:54 SRR1635864_1_aligned.fastq
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root   44120385 Oct 24 19:07 SRR1635864_1_aligned.blast
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root  144926228 Oct 24 19:07 SRR1635864_1_aligned.sam
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT -rw-r--r-- 1 root root        460 Oct 24 19:07 SRR1635864_1_aligned.log
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT [INFO] find /home/dnanexus/out -name '*.log' | xargs cat
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT  Results:
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     Total reads = 24654594
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     Total reads passing E-value threshold = 536480 (2.18)
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     Total reads failing E-value threshold = 24118114 (97.82)
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     Minimum read length = 55
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     Maximum read length = 101
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     Mean read length    = 99
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT  By database:
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     /home/dnanexus/in/refs/silva-arc-16s-id95.fasta           0.00
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     /home/dnanexus/in/refs/silva-bac-16s-id90.fasta           1.36
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT     /home/dnanexus/in/refs/silva-arc-23s-id98.fasta           0.81
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT  Wed Oct 24 18:48:39 2018
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT [INFO] find /home/dnanexus/out -name '*.fasta' | xargs wc -l
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT 0
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT [INFO] find /home/dnanexus/out -name '*.fastq' | xargs wc -l
2018-10-24 15:07:48 SortMeRNA 3 Run applet STDOUT    2145920 /home/dnanexus/out/SRR1635864_1_aligned.fastq
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT   96472456 /home/dnanexus/out/SRR1635864_1_other.fastq
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT   98618376 total
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT [INFO] find /home/dnanexus/out -name '*.blast' | xargs wc -l
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT 536480 /home/dnanexus/out/SRR1635864_1_aligned.blast
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT [INFO] find /home/dnanexus/out -name '*.sam' | xargs wc -l
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT 0
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT [INFO] Output aligned reads: /home/dnanexus/out/output_fastx_gz/SRR1635864_1_aligned.fastq.gz
2018-10-24 15:07:53 SortMeRNA 3 Run applet STDOUT Running: /home/dnanexus/out/SRR1635864_1_aligned.fastq > /home/dnanexus/out/output_fastx_gz/SRR1635864_1_aligned.fastq.gz
2018-10-24 15:08:02 SortMeRNA 3 Run applet STDOUT Uploaded: /home/dnanexus/out/output_fastx_gz/SRR1635864_1_aligned.fastq.gz. File ID: file-FP8G640049k9Bp867JX4ykQJ
2018-10-24 15:08:02 SortMeRNA 3 Run applet STDOUT [INFO] Output log file: /home/dnanexus/out/output_logfile/SRR1635864_1_aligned.log
2018-10-24 15:08:03 SortMeRNA 3 Run applet STDOUT [INFO] Output non-aligned reads: /home/dnanexus/out/output_other_gz/SRR1635864_1_other.fastq.gz
2018-10-24 15:16:01 SortMeRNA 3 Run applet INFO CPU: 24% (4 cores) * Memory: 176/7225MB * Storage: 69GB free * Net: 0↓/0↑MBps
2018-10-24 15:16:47 SortMeRNA 3 Run applet STDOUT [INFO] Output BLAST: /home/dnanexus/out/output_blast_gz/SRR1635864_1_aligned.blast.gz
2018-10-24 15:16:50 SortMeRNA 3 Run applet STDOUT ==== DONE ====
* SortMeRNA 3 Run applet (sortmerna-3.run.dev:main) (done) job-FP8F3k80g35bF08fPj6x1Zj8
  biocodz 2018-10-24 13:54:45 (runtime 1:21:02)
  Output: output_logfile = file-FP8G64Q049kKvZ7p7P1fK6Zq
          output_other_gz = file-FP8GB2j049k8PX7F7JjB19pQ
          output_fastx_gz = file-FP8G640049k9Bp867JX4ykQJ
          output_blast_gz = file-FP8GB88049kKvZ7p7P1fKV1G
</pre>