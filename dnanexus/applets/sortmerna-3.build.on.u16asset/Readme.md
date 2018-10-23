# App for building Sortmerna (DNAnexus Platform App)

The applet for building SortmeRNA 3 using the `sortmerna-3.asset`.

When the `sortmerna-3.asset` changes, the `dxapp.json` has to be modified to reference the updated asset e.g.

```
pushd SORTMERNA_HOME
vi dnanexus/applets/sortmerna-3.build.on.u16asset
...
"assetDepends": [{ "id": "record-FKJ5gfj0KpQffJ562PF2KyP8" }] # <-- modify record to reference correct asset
...