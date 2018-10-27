# Applet for building Sortmerna

The applet for building SortmeRNA 3 using the [sortmerna-3.asset](https://github.com/biocore/sortmerna/tree/master/dnanexus/assets/sortmerna-3.asset) on DNANexus platform.

The following artifacts are uploaded to the user account after the build
* indexdb
* sortmerna
* libstdc++.so.6

When the `sortmerna-3.asset` changes, the [dxapp.json](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.build.on.u16asset/dxapp.json) has to be modified to reference the updated asset e.g.

```
pushd $SMR_HOME # Sortmerna source top directory
dx build_asset -d assets/ dnanexus/assets/sortmerna-3.asset
dx find data
...
closed  2018-10-23 13:06:58    /assets/sortmerna-3.asset (record-FP7bBFj0GjzpPjxb2Q54XbV1)
...

# Modify the applet to use the new asset
vi dnanexus/applets/sortmerna-3.build.on.u16asset/dxapp.json
...
"assetDepends": [{ "id": "record-FP7bBFj0GjzpPjxb2Q54XbV1" }] # <-- reference the correct asset
...

# build the applet
dx build -f -d applets/ dnanexus/applets/sortmerna-3.build.on.u16asset
dx find data
...
closed  2018-10-23 13:10:15   /applets/sortmerna-3.build.on.u16asset (applet-FP7bFxj0g35zxFXbByq4fzy3)
...

# run the applet to build sortmerna
dx run applet-FP7bFxj0g35zxFXbByq4fzy3

# list the built artifacts
dx find data
...
closed  2018-10-23 13:12:40 170.47 KB /out/sortmerna-3.build.on.u16asset/indexdb (file-FP7bJ0Q0qqg69XxGB7gk1f54)
closed  2018-10-23 13:12:40 1.53 MB   /out/sortmerna-3.build.on.u16asset/libstdc++.so.6 (file-FP7bJ0j0qqg8Zz30Gzy142J4)
closed  2018-10-23 13:12:40 826.01 KB /out/sortmerna-3.build.on.u16asset/sortmerna (file-FP7bJ080qqgFJqP9FX0yf61J)
...
```
**NOTE**: Ellipsis `...` denotes skipped output (for readability)
