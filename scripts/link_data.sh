#!/bin/bash
data_dir=$2;
if [ $# -eq 0 ]; then
	echo "input path empty"
	exit
fi
if [ "$(ls -A $data_dir)" ]; then
	echo "not empty"
	rm $(data_dic)/*.tbl
fi
for entry in $1/*.tbl
do
	ln -s $(realpath $entry) $(data_dic)/$(basename $entry)
done
