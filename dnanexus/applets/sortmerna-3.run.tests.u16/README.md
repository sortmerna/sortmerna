# Applet for running integration tests provided with Sortmerna

The applet uses the `sortmerna-3.run.asset` asset and takes the Sortmerna binaries as input.  
If the asset was changed, modify `dxapp.json` to reference the new asset record.

For example

```
pushd $SORTMERNA_HOME
dx build_asset -d assets/ dnanexus/assets/sortmerna-3.run.asset
dx find data
...
closed  2018-10-22 15:15:23           /assets/sortmerna-3.run.asset (record-FP723X00qpjg0G9f0xzjbQbP)
...

vi dnanexus/applets/sortmerna-3.run.tests.u16
...
"assetDepends": [{ "id": "record-FP723X00qpjg0G9f0xzjbQbP" }]  # <-- modify the asset reference
...

dx build -f -d applets/ dnanexus/applets/sortmerna-3.run.tests.u16
dx find data
closed  2018-10-23 08:53:32           /applets/sortmerna-3.run.tests.u16 (applet-FP7VYY00g35z7pq19PBX1VK5)
...
closed  2018-10-22 14:07:22 826.01 KB /out/sortmerna/snapshot/sortmerna (file-FP713f00g35ZXQvK32KV4ZjY)
closed  2018-10-22 14:07:12 170.47 KB /out/sortmerna/snapshot/indexdb (file-FP713Yj0g35gVjBPGGz6jygF)
closed  2018-10-22 10:54:49           /assets/sortmerna-3.run.asset (record-FP6y9Jj0zj8fQyf583VY52FF)
closed  2018-10-10 08:40:37 1.53 MB   /out/sortmerna/snapshot/libstdc++.so.6 (file-FKyz6QQ0g35X7q9ZP47YFbVJ)
...

dx run applet-FP7VYY00g35z7pq19PBX1VK5
...
BINS[0]: file-FP713f00g35ZXQvK32KV4ZjY # sortmerna
BINS[1]: file-FP713Yj0g35gVjBPGGz6jygF # indexdb
BINS[2]: file-FKyz6QQ0g35X7q9ZP47YFbVJ # libstdc++.so.6
...

```