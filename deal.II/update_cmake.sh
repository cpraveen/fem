# Copy CMakeLists.txt file here and run this to update all the files
#
DIRS=`ls -d ex*`
cmake=CMakeLists.txt

for dir in $DIRS
do
   echo $dir
   cp $cmake $dir/
done
