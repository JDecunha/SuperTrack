output=$(../build/tst/SuperTrack_tst --gtest_list_tests) #get the tests

while IFS= read -r line ; do #loop through only the tests not containing a period (i.e the individual tests and not the test suites)
	if [[ ${line} != *"."* ]];then
		writevar=$(echo ${line} | xargs) #trim the whitespace automatically with xargs
    	../build/tst/SuperTrack_tst --gtest_filter="*${writevar}" #run the test with a * wildcard at the beginning
	fi
done <<< "$output" #only stop looping one output is complete