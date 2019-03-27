#! /usr/bin/env bash

# make sure results directory exists
DIRECTORY=results
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi

echo 'Running ./mainHIC' $@

# time and run
nohup time ./mainHIC $@\
			1> $DIRECTORY/mainHIC.out\
			2> $DIRECTORY/mainHIC.err

echo 'Finished.'
