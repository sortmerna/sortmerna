# Asset for building Sortmerna on DNANexus platform

The asset packages the environment necessary to build the Sortmerna including

* gcc
* cmake
* patchelf
* zlib1g-dev
* librocksdb-dev
* rapidjson-dev

The built asset is used by [sortmerna-3.build.on.u16asset](https://github.com/biocore/sortmerna/tree/master/dnanexus/applets/sortmerna-3.build.on.u16asset)

Example running (assumes DNANexus `dx` toolkit is installed)

```
# clone sortmerna
git clone https://github.com/biocore/sortmerna.git
pushd sortmerna

# build the asset
dx build_asset -d assets/ dnanexus/assets/sortmerna-3.asset
dx find data
...
closed  2018-10-23 13:06:58    /assets/sortmerna-3.asset (record-FP7bBFj0GjzpPjxb2Q54XbV1)
...
```

**NOTE**: Ellipsis `...` denotes skipped output (for readability)
