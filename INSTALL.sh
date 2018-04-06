### install python dependence
pip install --upgrade pip
pip install --upgrade numpy

### install R dependence
Rscript bin/install_R.R

### install bedtools
cd bin
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make
cd ..
rm bedtools-2.25.0.tar.gz
cd ..

