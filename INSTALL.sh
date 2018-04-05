### install python dependence
sudo pip install --upgrade pip
sudo pip install --upgrade numpy
sudo pip install --upgrade subprocess
sudo pip install --upgrade collections
sudo pip install --upgrade getopt
sudo pip install --upgrade sys
sudo pip install --upgrade os

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

