# 1001Gplus_paper

## Create folders
If you want to create an analysis folder, you can use the script utils/create_directories.sh:

```
./utils/create_directories.sh 02_analysis/xx_yyy
```

The structure of the folder will be automatically created.


## Commit scripts and figures

```
cd 02_analysis
find . -type d -name '02_scripts' -exec git add {} +
find . -type d -name '03_figures' -exec git add {} +
git commit -m "scripts and figures updated"
git push
```


## Important comments to share
* Reconstructed 22001 is called 220011
* ..


