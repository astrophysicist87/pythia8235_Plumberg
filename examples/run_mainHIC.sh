#! /usr/bin/env bash

# make sure results directory exists
DIRECTORY="$5"
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi

echo 'Running ./mainHIC' $@

# time and run
nohup time ./mainHIC $@\
			1> $DIRECTORY/mainHIC.out\
			2> $DIRECTORY/mainHIC.err
#nohup time ./mainHIC_subcollisions $@\
#			1> $DIRECTORY/mainHIC_subcollisions.out\
#			2> $DIRECTORY/mainHIC_subcollisions.err
