This folder contains source code that allows for building, testing, and running Sortmerna on [DNANexus](https://www.dnanexus.com/) platform.

The following applets and assets are available:

dir | title | description
----|-------|------------
applets | sortmerna-3.build.on.u16asset | build Sortmerna 3 on Ubuntu 16.04 using an Asset containing the dependencies [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.build.on.u16asset/Readme.md)
applets | sortmerna-3.build.on.u14asset | build Sortmerna 3 on Ubuntu 14.04 using an Asset containing the dependencies [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.build.on.u14asset/Readme.md)
applets | sortmerna-3.run.tests.u14 | run SortmeRNA 3 integration tests on Ubuntu 14.04 [README](https://github.com/biocore/sortmerna/tree/master/dnanexus/applets/sortmerna-3.run.tests.u14)
applets | sortmerna-3.run.tests.u16 | run SortmeRNA 3 integration tests on Ubuntu 16.04 using binaries specified as input  [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.run.tests.u16/README.md)
applets | sortmerna-3.run.dev | run SortMeRNA 3 on Ubuntu 16.04 using binaries specified as input [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.run.dev/README.md)
applets | sortmerna-3.run.u14 | run SortMeRNA 3 on Ubuntu 14.04 [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.run.u14/Readme.developer.md)
applets | sortmerna | run SortMeRNA on Ubuntu 16.04 [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna/Readme.md)
assets | sortmerna-3.asset | Ubuntu 16.04 Environment containing SortmeRNA build tools, libraries, and dependencies [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/assets/sortmerna-3.asset/README.md)
assets | sortmerna-3.asset.u14 | DEPRECATED: Ubuntu 14.04 Environment required to build SortmeRNA 3 [README](https://github.com/biocore/sortmerna/tree/master/dnanexus/assets/sortmerna-3.asset.u14)
assets | sortmerna-3.run.asset | Ubuntu 16.04 Environment containing SortmeRNA dependencies and Conda [README](https://github.com/biocore/sortmerna/blob/master/dnanexus/assets/sortmerna-3.run.asset/README.md)

## Usage scenarios

1. [Build and Run Sortmerna on DNANexus](#build-and-run-sortmerna-on-dnanexus) in two steps:
   1. build
   2. run using the binaries from step 1
2. [Build and Run Sortmerna integration Tests](#build-and-run-sortmerna-integration-tests)
3. [Run `ready-to-use` Sortmerna application on your data](#run-sortmerna-application-on-your-data)

### Build and Run Sortmerna on DNANexus

Assuming you have an account with DNANexus (free trial account can be used) and the `dx` toolkit installed. The `dx` toolkit is a Python package easily installed using e.g. `pip`.

On your local machine with installed `dx` perform the following

```
dx login
# clone Sortmerna git repo on your local machine
git clone https://github.com/biocore/sortmerna.git
pushd sormerna

# build the asset
dx build_asset -d assets/ dnanexus/assets/sortmerna-3.asset
dx find data
...
closed  2018-09-13 08:45:18   /assets/sortmerna-3.asset (record-FKJ5gfj0KpQffJ562PF2KyP8) # <-- Note the asset ID
...

# modify the build applet to reference the asset
vi dnanexus/applets/sortmerna-3.build.on.u16asset/dxapp.json
...
  "runSpec": {
    ...
    "assetDepends": [{ "id": "record-FKJ5gfj0KpQffJ562PF2KyP8" }] # <-- reference the asset 
  },
...

# build the applet
dx build -f -d applets/ dnanexus/applets/sortmerna-3.build.on.u16asset
dx find data
...
closed  2018-10-23 13:10:15  /applets/sortmerna-3.build.on.u16asset (applet-FP7bFxj0g35zxFXbByq4fzy3)
...

# run the applet to build sortmerna
dx run applet-FP7bFxj0g35zxFXbByq4fzy3 # <-- see REFS [1] below for full execution trace
...
  Output: sortmerna = file-FP7bJ080qqgFJqP9FX0yf61J
          indexdb = file-FP7bJ0Q0qqg69XxGB7gk1f54
          libstdcpp = file-FP7bJ0j0qqg8Zz30Gzy142J4
...

# list the build artifacts
dx find data
...
closed  2018-10-23 13:12:40 170.47 KB /out/sortmerna-3.build.on.u16asset/indexdb (file-FP7bJ0Q0qqg69XxGB7gk1f54)
closed  2018-10-23 13:12:40 1.53 MB   /out/sortmerna-3.build.on.u16asset/libstdc++.so.6 (file-FP7bJ0j0qqg8Zz30Gzy142J4)
closed  2018-10-23 13:12:40 826.01 KB /out/sortmerna-3.build.on.u16asset/sortmerna (file-FP7bJ080qqgFJqP9FX0yf61J)
...

# build the applet for running Sortmerna
dx build -f -d applets/ dnanexus/applets/sortmerna-3.run.dev

# list the built applet and the input data
dx find data
...
closed  2018-10-24 10:32:54   /applets/sortmerna-3.run.dev (applet-FP885580g35ZzzP65YPb4kjX)
...
closed  2018-09-13 07:31:46 1.92 GB   /data/reads/SRR1635864_1.fastq.gz (file-FKJ4ZpQ01FFKyQv8194pKv29)
...
closed  2018-08-13 09:11:00 3.71 MB   /data/refs/silva-arc-16s-id95.fasta (file-FJZz7300V8Pb4KZgBGkXbZBq)
closed  2018-08-13 09:11:00 18.54 MB  /data/refs/silva-bac-16s-id90.fasta (file-FJZz7Fj0V8PX7G23K0Yf2Z2J)
closed  2018-08-13 09:11:00 734.40 KB /data/refs/silva-arc-23s-id98.fasta (file-FJZz78j0V8Pp0Y6Z07zK7K86)
...

# run sortmerna
dx run applet-FP885580g35ZzzP65YPb4kjX # <-- See REFS [2] below for full execution trace
...
2018-10-24 13:56:00 SortMeRNA 3 Run applet INFO Logging initialized (priority)
2018-10-24 13:56:01 SortMeRNA 3 Run applet INFO CPU: 10% (4 cores) * Memory: 283/7225MB * Storage: 79GB free * Net: 0↓/0↑MBps
2018-10-24 13:56:02 SortMeRNA 3 Run applet STDOUT dxpy/0.266.1 (Linux-4.4.0-98-generic-x86_64-with-Ubuntu-16.04-xenial)
2018-10-24 13:56:03 SortMeRNA 3 Run applet INFO Installing apt packages patchelf, zlib1g-dev, librocksdb-dev, rapidjson-dev, samtools
...
2018-10-24 13:56:40 SortMeRNA 3 Run applet STDOUT Starting: indexdb --ref /home/dnanexus/in/refs/silva-arc-16s-id95.fasta,/home/dnanexus/idx/silva-arc-16s-id95:/home/dnanexus/in/refs/silva-bac-16s-id90.fasta,/home/dnanexus/idx/silva-bac-16s-id90:/home/dnanexus/in/refs/silva-arc-23s-id98.fasta,/home/dnanexus/idx/silva-arc-23s-id98 -v
...
2018-10-24 13:57:50 SortMeRNA 3 Run applet STDOUT Starting: sortmerna --ref /home/dnanexus/in/refs/silva-arc-16s-id95.fasta,/home/dnanexus/idx/silva-arc-16s-id95:/home/dnanexus/in/refs/silva-bac-16s-id90.fasta,/home/dnanexus/idx/silva-bac-16s-id90:/home/dnanexus/in/refs/silva-arc-23s-id98.fasta,/home/dnanexus/idx/silva-arc-23s-id98 --reads-gz /home/dnanexus/in/reads/SRR1635864_1.fastq.gz --aligned /home/dnanexus/out/SRR1635864_1_aligned --other /home/dnanexus/out/SRR1635864_1_other --sam --fastx --log --blast "1 cigar qcov qstrand" -v -d kvdb --task 4
...
2018-10-24 15:16:50 SortMeRNA 3 Run applet STDOUT ==== DONE ====
* SortMeRNA 3 Run applet (sortmerna-3.run.dev:main) (done) job-FP8F3k80g35bF08fPj6x1Zj8
  biocodz 2018-10-24 13:54:45 (runtime 1:21:02)
  Output: output_logfile = file-FP8G64Q049kKvZ7p7P1fK6Zq
          output_other_gz = file-FP8GB2j049k8PX7F7JjB19pQ
          output_fastx_gz = file-FP8G640049k9Bp867JX4ykQJ
          output_blast_gz = file-FP8GB88049kKvZ7p7P1fKV1G
...
```

### Build and Run Sortmerna integration Tests

TODO

### Run Sortmerna application on your data

TODO

## REFS

1. [Execution trace for SortMeRNA 3 Build using Ubuntu 16.04 Asset applet](https://github.com/biocore/sortmerna/tree/master/docs/traces/smr_build_on_ubuntu16_asset.md)
2. [Execution trace for SortMeRNA 3 Run applet (sortmerna-3.run.dev)](https://github.com/biocore/sortmerna/tree/master/docs/traces/smr_run_on_ubuntu16.md)
