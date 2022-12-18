#!/bin/bash

for i in $@
do 
	touch result
	./main $i >> result
done

echo "程序输出结果已追加至./result文件中!"
