#!/bin/bash

dir=$(pwd -P)


echo "Installation location is $dir"
echo "Setting up single.sh script with appropriate reference location"

dir=${dir//\//\\/}
sed -i "s/location/$dir/g" ./single.sh

if [ $? -eq 0 ]
then
    echo "single.sh is ready to use"
else
    echo "non 0 exit code $? plese see error above"
fi
