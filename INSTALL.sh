### install python dependence
pip install --upgrade pip --user
pip install --upgrade numpy --user

### install R dependence
Rscript bin/install_R.R

### install bedtools
cd bin
cd bedtools2
make clean
make
cd ..
cd ..

