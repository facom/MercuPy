#!/bin/bash
export MPLCONFIGDIR="."

echo "Regenerating output..."
cp -rf data.out xv.out
make prepare
./element6.sh
make out
rm -rf xv.out

echo "Regenerating change of reference frame..."
make ref

echo "Computing close..."
if [ 0 -gt 0 ];then
    cp -rf close.out ce.out
    make close
fi

echo "Plotting results..."
make plot

if [ 0 -gt 0 ];then
    echo "Computing errors..."
    bin/mercupy-diff output/BODY3.dat.PH1 output/BODY3.dat.PH2 
fi

echo "Saving results variables..."
bash results.sh
