#!/bin/bash

SCRIPT_DIR=$(dirname "$0")

# Apply to given files
if [ $# -gt 0 ]
then
  for file in "$@"
  do
    echo $file
    astyle --options=$SCRIPT_DIR/astyle.rc $file
  done
  exit 0
fi

# Apply to all files and in all sub-directories
HFILES=`find . -name "*.h"`
CCFILES=`find . -name "*.cc"`
FILES="$HFILES $CCFILES"
for file in $FILES
do
   echo $file
   astyle --options=$SCRIPT_DIR/astyle.rc $file
done
