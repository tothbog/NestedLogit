#!/usr/bin/env bash
# Creates the (lambda_1,lambda_c) variants -12,-21,-22,-51,-52 of a generated r-*-11 instance.
pr=$1; #problem like r-100-2-1-0-11


sed 's/la1 1/la1\t2.0/;s/lac  1/lac\t1.0/;' $pr >${pr/-11/-21}
sed 's/la1 1/la1\t2.0/;s/lac  1/lac\t2.0/;' $pr >${pr/-11/-22}
sed 's/la1 1/la1\t1.0/;s/lac  1/lac\t2.0/;' $pr >${pr/-11/-12}
# (1.5, 1), (1.5, 2)
sed 's/la1 1/la1\t1.5/;s/lac  1/lac\t1.0/;' $pr >${pr/-11/-51}
sed 's/la1 1/la1\t1.5/;s/lac  1/lac\t2.0/;' $pr >${pr/-11/-52}
