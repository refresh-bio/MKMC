Po sklonowaniu a przed otwarciem .sln trzeba:
```
git submodule update --init --recursive
```
ewentualnie przy klonowaniu:
```
git clone --recurse-submodules https://github.com/refresh-bio/MKMC-dev/
```

Tutaj jest moje repo z dość starą i nawiną implementacją:
https://github.com/marekkokot/kmc-matrixer
W szczególności najważniejszy jest skrypt https://github.com/marekkokot/kmc-matrixer/blob/main/run.sh
