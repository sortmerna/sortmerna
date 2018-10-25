## SortMeRNA 3 Build using Ubuntu 16.04 Asset applet execution trace

<pre>
* SortMeRNA 3 Build using Ubuntu 16.04 Asset applet (sortmerna-3.build.on.u16asset:main) (done) job-FP7bG8Q0g35Q8K02G9x3BBpP
  biocodz 2018-10-23 13:10:58 (runtime 0:01:30)
2018-10-23 13:11:20 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet INFO Logging initialized (priority)
2018-10-23 13:11:22 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet INFO CPU: 18% (4 cores) * Memory: 284/7225MB * Storage: 79GB free * Net: 0↓/0↑MBps
2018-10-23 13:11:22 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT dxpy/0.264.0 (Linux-4.4.0-98-generic-x86_64-with-Ubuntu-16.04-xenial)
2018-10-23 13:11:24 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet INFO Downloading bundled file sortmerna-3.asset.tar.gz
2018-10-23 13:11:26 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT >>> Unpacking sortmerna-3.asset.tar.gz to /
2018-10-23 13:11:26 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR tar: Removing leading `/' from member names
2018-10-23 13:11:28 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet INFO Setting SSH public key
2018-10-23 13:11:29 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT dxpy/0.264.0 (Linux-4.4.0-98-generic-x86_64-with-Ubuntu-16.04-xenial)
2018-10-23 13:11:29 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT /usr/sbin/sshd already running.
2018-10-23 13:11:29 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT /usr/sbin/rsyslogd already running.
2018-10-23 13:11:31 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT bash running (job ID job-FP7bG8Q0g35Q8K02G9x3BBpP)
2018-10-23 13:11:31 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Upload directory: /out
2018-10-23 13:11:31 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Applet/Executable name: sortmerna-3.build.on.u16asset
2018-10-23 13:11:31 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Upload will be done into /out/sortmerna-3.build.on.u16asset
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Upload path status:
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Testing: which gcc-7
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT /usr/bin/gcc-7
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Testing: gcc-7 --version
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT gcc-7 (Ubuntu 7.3.0-21ubuntu1~16.04) 7.3.0
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Copyright (C) 2017 Free Software Foundation, Inc.
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT This is free software; see the source for copying conditions.  There is NO
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Testing: apt-cache policy librocksdb-dev
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT librocksdb-dev:
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Installed: 4.1-1
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Candidate: 4.1-1
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Version table:
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT  *** 4.1-1 500
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial/universe amd64 Packages
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT         100 /var/lib/dpkg/status
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Testing: apt-cache policy zlib1g-dev
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT zlib1g-dev:
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Installed: 1:1.2.8.dfsg-2ubuntu4.1
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Candidate: 1:1.2.8.dfsg-2ubuntu4.1
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Version table:
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT  *** 1:1.2.8.dfsg-2ubuntu4.1 500
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial-updates/main amd64 Packages
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT         100 /var/lib/dpkg/status
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT      1:1.2.8.dfsg-2ubuntu4 500
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial/main amd64 Packages
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Testing: apt-cache policy rapidjson-dev
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT rapidjson-dev:
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Installed: 0.12~git20141031-3
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Candidate: 0.12~git20141031-3
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT   Version table:
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT  *** 0.12~git20141031-3 500
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT         500 http://us-east-1.ec2.archive.ubuntu.com/ubuntu xenial/universe amd64 Packages
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT         100 /var/lib/dpkg/status
2018-10-23 13:11:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR Cloning into 'sortmerna'...
2018-10-23 13:11:35 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT ~/sortmerna/build/Release ~
2018-10-23 13:11:35 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- The CXX compiler identification is GNU 7.3.0
2018-10-23 13:11:35 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- The C compiler identification is GNU 7.3.0
2018-10-23 13:11:35 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Check for working CXX compiler: /usr/bin/c++
2018-10-23 13:11:35 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Check for working CXX compiler: /usr/bin/c++ -- works
2018-10-23 13:11:35 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting CXX compiler ABI info
2018-10-23 13:11:35 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting CXX compiler ABI info - done
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting CXX compile features
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting CXX compile features - done
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Check for working C compiler: /usr/bin/cc
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Check for working C compiler: /usr/bin/cc -- works
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting C compiler ABI info
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting C compiler ABI info - done
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting C compile features
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Detecting C compile features - done
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_CXX_COMPILER_ID = GNU
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_CONFIGURATION_TYPES =
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_HOST_SYSTEM_NAME = Linux
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_HOST_SYSTEM_VERSION = 4.4.0-98-generic
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_HOST_SYSTEM = Linux-4.4.0-98-generic
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR PROJECT_BINARY_DIR = /home/dnanexus/sortmerna/build/Release
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_BINARY_DIR = /home/dnanexus/sortmerna/build/Release
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_CURRENT_BINARY_DIR = /home/dnanexus/sortmerna/build/Release
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR RUNTIME_OUTPUT_DIRECTORY_RELEASE =
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR CMAKE_CXX_FLAGS_RELEASE: -O3 -DNDEBUG
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR EXTRA_CXX_FLAGS_RELEASE:
2018-10-23 13:11:36 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR Cloning into 'concurrentqueue'...
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Found Git: /usr/bin/git (found version "2.7.4")
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread.h
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread.h - found
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread_create
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread_create - not found
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread_create in pthreads
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread_create in pthreads - not found
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread_create in pthread
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Looking for pthread_create in pthread - found
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Found Threads: TRUE
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Configuring done
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Generating done
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT -- Build files have been written to: /home/dnanexus/sortmerna/build/Release
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT ~
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT ~/sortmerna/build/Release ~
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Scanning dependencies of target build_version
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [  1%] Building CXX object CMakeFiles/build_version.dir/build_version.cpp.o
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [  1%] Built target build_version
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Scanning dependencies of target cmph
2018-10-23 13:11:37 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [  3%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/bdz.c.o
2018-10-23 13:11:38 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [  4%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/bdz_ph.c.o
2018-10-23 13:11:38 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [  6%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/bmz.c.o
2018-10-23 13:11:38 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [  7%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/bmz8.c.o
2018-10-23 13:11:38 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [  9%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/brz.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 10%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/buffer_entry.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 12%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/buffer_manager.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 13%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/chd.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 15%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/chd_ph.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 16%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/chm.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 18%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/cmph.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 20%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/cmph_structs.c.o
2018-10-23 13:11:39 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 21%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/compressed_rank.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 23%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/compressed_seq.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 24%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/fch.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 26%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/fch_buckets.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 27%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/graph.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 29%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/hash.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 30%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/jenkins_hash.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 32%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/miller_rabin.c.o
2018-10-23 13:11:40 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 33%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/select.c.o
2018-10-23 13:11:41 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 35%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/vqueue.c.o
2018-10-23 13:11:41 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 36%] Building C object 3rdparty/cmph/CMakeFiles/cmph.dir/vstack.c.o
2018-10-23 13:11:41 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 36%] Built target cmph
2018-10-23 13:11:41 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Scanning dependencies of target alp
2018-10-23 13:11:41 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 38%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_dynprogprob.cpp.o
2018-10-23 13:11:41 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 40%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_dynprogproblim.cpp.o
2018-10-23 13:11:42 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 41%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_dynprogprobproto.cpp.o
2018-10-23 13:11:42 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 43%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_ioutil.cpp.o
2018-10-23 13:11:43 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 44%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_localmaxstat.cpp.o
2018-10-23 13:11:43 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 46%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_localmaxstatmatrix.cpp.o
2018-10-23 13:11:44 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 47%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_localmaxstatutil.cpp.o
2018-10-23 13:11:45 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 49%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/njn_random.cpp.o
2018-10-23 13:11:45 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 50%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/sls_alignment_evaluer.cpp.o
2018-10-23 13:11:46 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 52%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/sls_alp.cpp.o
2018-10-23 13:11:48 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 53%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/sls_alp_data.cpp.o
2018-10-23 13:11:50 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 55%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/sls_alp_regression.cpp.o
2018-10-23 13:11:51 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 56%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/sls_alp_sim.cpp.o
2018-10-23 13:11:53 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 58%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/sls_basic.cpp.o
2018-10-23 13:11:54 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 60%] Building CXX object 3rdparty/alp/CMakeFiles/alp.dir/sls_pvalues.cpp.o
2018-10-23 13:11:55 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 61%] Linking CXX static library libalp.a
2018-10-23 13:11:55 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 61%] Built target alp
2018-10-23 13:11:55 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Scanning dependencies of target indexdb
2018-10-23 13:11:55 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 63%] Building CXX object src/indexdb/CMakeFiles/indexdb.dir/indexdb.cpp.o
2018-10-23 13:11:56 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR /home/dnanexus/sortmerna/src/indexdb/indexdb.cpp: In function ‘int main(int, char**)’:
2018-10-23 13:11:56 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR /home/dnanexus/sortmerna/src/indexdb/indexdb.cpp:1563:58: warning: format ‘%llu’ expects argument of type ‘long long unsigned int’, but argument 3 has type ‘uint64_t {aka long unsigned int}’ [-Wformat=]
2018-10-23 13:11:56 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR      fprintf(stderr, "\n  check sequence # %llu\n\n", strs);
2018-10-23 13:11:56 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDERR                                                           ^
2018-10-23 13:11:57 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 64%] Linking CXX executable indexdb
2018-10-23 13:11:57 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 64%] Built target indexdb
2018-10-23 13:11:57 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Scanning dependencies of target smr_objs
2018-10-23 13:11:57 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 66%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/alignment.cpp.o
2018-10-23 13:11:59 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 67%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/bitvector.cpp.o
2018-10-23 13:12:00 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 69%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/callbacks.cpp.o
2018-10-23 13:12:02 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 70%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/cmd.cpp.o
2018-10-23 13:12:05 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 72%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/gzip.cpp.o
2018-10-23 13:12:05 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 73%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/index.cpp.o
2018-10-23 13:12:06 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 75%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/kseq_load.cpp.o
2018-10-23 13:12:07 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 76%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/kvdb.cpp.o
2018-10-23 13:12:08 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 78%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/options.cpp.o
2018-10-23 13:12:09 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 80%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/output.cpp.o
2018-10-23 13:12:13 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 81%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/paralleltraversal.cpp.o
2018-10-23 13:12:16 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 83%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/processor.cpp.o
2018-10-23 13:12:20 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 84%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/read.cpp.o
2018-10-23 13:12:22 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 86%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/reader.cpp.o
2018-10-23 13:12:24 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 87%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/readstats.cpp.o
2018-10-23 13:12:26 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 89%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/references.cpp.o
2018-10-23 13:12:27 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 90%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/refstats.cpp.o
2018-10-23 13:12:29 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 92%] Building C object src/sortmerna/CMakeFiles/smr_objs.dir/ssw.c.o
2018-10-23 13:12:29 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 93%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/traverse_bursttrie.cpp.o
2018-10-23 13:12:30 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 95%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/util.cpp.o
2018-10-23 13:12:30 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 96%] Building CXX object src/sortmerna/CMakeFiles/smr_objs.dir/writer.cpp.o
2018-10-23 13:12:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 96%] Built target smr_objs
2018-10-23 13:12:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Scanning dependencies of target sortmerna
2018-10-23 13:12:32 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [ 98%] Building CXX object src/sortmerna/CMakeFiles/sortmerna.dir/main.cpp.o
2018-10-23 13:12:33 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [100%] Linking CXX executable sortmerna
2018-10-23 13:12:33 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT [100%] Built target sortmerna
2018-10-23 13:12:33 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT ~
2018-10-23 13:12:34 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT libstdcline:       libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007fa6b6379000)
2018-10-23 13:12:34 SortMeRNA 3 Build using Ubuntu 16.04 Asset applet STDOUT Uploading output: /usr/lib/x86_64-linux-gnu/libstdc++.so.6
* SortMeRNA 3 Build using Ubuntu 16.04 Asset applet (sortmerna-3.build.on.u16asset:main) (done) job-FP7bG8Q0g35Q8K02G9x3BBpP
  biocodz 2018-10-23 13:10:58 (runtime 0:01:30)
  Output: sortmerna = file-FP7bJ080qqgFJqP9FX0yf61J
          indexdb = file-FP7bJ0Q0qqg69XxGB7gk1f54
          libstdcpp = file-FP7bJ0j0qqg8Zz30Gzy142J4
</pre>