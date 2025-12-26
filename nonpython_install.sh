#!/usr/bin/env bash

cd /usr/bin/

# BWA aligner
git clone https://github.com/lh3/bwa.git
cd bwa
make

# Copy the binary to /usr/local/bin
cp bwa /usr/local/bin/

# infernal
cd /usr/bin/
wget http://eddylab.org/software/infernal/infernal.tar.gz
tar zxf infernal.tar.gz
cd infernal-1.1.5
./configure --prefix /usr/local/bin/
make
make check                 # optional: run automated tests
make install               # optional: install Infernal programs, man pages
cp src/cmscan /usr/local/bin/  # copy cmscan binary to /usr/local/bin


# bcftools
apt install bcftools
apt install bedtools