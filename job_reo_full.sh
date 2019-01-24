#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M kaefer@uni-bonn.de
#$ -m be
#$ -N reo_full
set -e

perl /home/skaefer/reo_full/TRAVIS_Core_v0.7c.pl /home/skaefer/reo_full/reo_full_TCC.csv
