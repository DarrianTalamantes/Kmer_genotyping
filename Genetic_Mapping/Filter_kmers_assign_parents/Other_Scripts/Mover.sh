for file in $(ls | grep "code")
do
prog=$(wc -l $file | cut -d  " " -f 1)
if (($prog >= 5))
	then
	echo $x
	arrIN=(${file//_/ })
	echo ${arrIN[0]}
	cat $file | cut -d "," -f 1 > Parental_Lists/${arrIN[0]}_list.txt
fi
done
