#!/bin/bash
outdir=output

function changeAEI
{
    ext=$1;shift
    echo "Canging aei files to extension $ext..."
    for file in *.aei
    do
	objname=$(echo $file | awk -F'.aei' '{print $1}')
	echo "Changing $file to $objname.$ext..."
	mv $objname.aei $outdir/$objname.$ext
    done
}

for element in element-*.in
do
    echo "FILE: $element"
    comp=$(echo $element | cut -d '-' -f 2 | cut -d '.' -f 1)
    if [ $comp != "element" ]
    then
	cp -rf $element element.in
    else
	comp="aei"
    fi
    echo "Generating elements..."
    ./element6.exe
    changeAEI $comp
done
