#!/usr/bin/env bash

prog=ProblemGenerator-NestedLogit
case  $#  in 
1)    ;;    
*)    echo "This program run $prog."
      echo "The parameters are: 1. file with the parameters"
      echo "Example: run_generator.sh params.txt"
      exit 1;;
esac

file=$1

x=`wc -l <$file`
echo $x

for (( j = 1; j <= x; j++ ))
do
    params=`head -$j $file | tail -1`
    torun="$prog $params"  
    ./$torun  

done

