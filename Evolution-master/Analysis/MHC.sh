#!/bin/bash

ssh mutant"$1" <<EOF
cd /linuxhome/tmp/nikolas/"$2"/
./SIscript.sh -mhc
echo "Done"
EOF
