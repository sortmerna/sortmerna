# Asset for running Sortmerna integration tests

This asset packages the necessary Python and the `scikit-bio`, `numpy`, `scipy` packages using `conda`

To build the asset (assuming the DNANexus `dx` kit is installed):

```
git clone https://github.com/biocore/sortmerna.git
pushd sortmerna
dx build_asset -d assets/ dnanexus/assets/sortmerna-3.run.asset
dx find data
...
closed  2018-10-22 15:15:23           /assets/sortmerna-3.run.asset (record-FP723X00qpjg0G9f0xzjbQbP)
...
```