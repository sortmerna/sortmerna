# atropos_detect Developer Readme

<!--
-->

## Running this app with additional computational resources

Prior building this applet (dx build ...), create subdirectory /applets/sortmerna-3.run/resources/home/dnanexus/bin/ on the Projects machine, and copy there the required binaries: 'sortmerna', 'indexdb', 'libstdc++.so.6'

For example:

mkdir -p applets/sortmerna-3.run/resources/home/dnanexus/bin
dx download -o applets/sortmerna-3.run/resources/home/dnanexus/bin file-FK28vjQ0GjGVYQyFP5J95PGg
    applets/sortmerna-3.run/resources/home/dnanexus/bin/sortmerna
dx download -o applets/sortmerna-3.run/resources/home/dnanexus/bin file-FK28vk00GjGj56jgP4Qf27KQ
    applets/sortmerna-3.run/resources/home/dnanexus/bin/libstdc++.so.6
dx download -o applets/sortmerna-3.run/resources/home/dnanexus/bin file-FK28vjj0GjGZ99zZP5680bJ3
    applets/sortmerna-3.run/resources/home/dnanexus/bin/indexdb
dx build -f -d applets/ applets/sortmerna-3.run
