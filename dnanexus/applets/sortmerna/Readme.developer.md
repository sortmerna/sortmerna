## Building the applet

This applet requires the sortmerna binaries: `sortmerna, indexdb, libstdc++.so.6`. The binaries can be built either on the local user machine, or using [sortmerna-3.build.on.u16asset](https://github.com/biocore/sortmerna/tree/master/dnanexus/applets/sortmerna-3.build.on.u16asset)


Example using binaries built on `sortmerna-3.build.on.u16asset`:

```
# on your local machine navigate to the sortmerna distribution directory e.g. /home/biocodz/sortmerna/
pushd $SORTMERNA_HOME

# Create directory for binaries
mkdir -p dnanexus/applets/sortmerna/resources/home/dnanexus/bin

#
# Download binaries built on the DNANexus Project host into your local machine
#
dx ls file-FK28vjQ0GjGVYQyFP5J95PGg
    sortmerna
dx ls file-FK28vk00GjGj56jgP4Qf27KQ
    libstdc++.so.6
dx ls file-FK28vjj0GjGZ99zZP5680bJ3
    indexdb

dx download -f -o dnanexus/applets/sortmerna/resources/home/dnanexus/bin/ file-FK28vjQ0GjGVYQyFP5J95PGg  # sortmerna
dx download -f -o dnanexus/applets/sortmerna/resources/home/dnanexus/bin/ file-FK28vk00GjGj56jgP4Qf27KQ  # libstdc++.so.6
dx download -f -o dnanexus/applets/sortmerna/resources/home/dnanexus/bin/ file-FK28vjj0GjGZ99zZP5680bJ3  # indexdb

# build the applet
dx build -f -d applets/ dnanexus/applets/sortmerna
```
