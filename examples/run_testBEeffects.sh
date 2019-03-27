#! /usr/bin/env bash

# make sure results directory exists
DIRECTORY="$5"
if [ ! -d "$DIRECTORY" ]; then
	mkdir -p $DIRECTORY
fi

echo 'Running ./main_testBEeffects' $@

# time and run
nohup time ./main_testBEeffects $@\
			1> $DIRECTORY/main_testBEeffects.out\
			2> $DIRECTORY/main_testBEeffects.err
