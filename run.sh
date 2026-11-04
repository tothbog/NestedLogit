#!/usr/bin/env bash

prog=NestedLogit
case  $#  in 
[23])    ;;    
*)    echo "Runs $prog (always with -mt) on every instance listed in a file."
      echo 'Parameters: 1. file with the instance paths, 2. tolerance (-ep), [3. "extra options"]'
      echo "Example: ./run.sh p-e1 1e-6 -fv"
      exit 1;;
esac

file=$1
eps=$2
params="$3"

for i in `cat $file`
do
	    torun="$prog -f $i -ep $eps -mt $params"  
	    echo "$torun running...."
	    Ofile="${i}--${eps}${params// /}"
	    ./$torun 2>$Ofile 
	    if grep '^----------------------' ${Ofile}  >/dev/null
	       then
	         echo "$torun finished."
	       else
	         echo "$torun finished with ERROR."
	    fi
done
