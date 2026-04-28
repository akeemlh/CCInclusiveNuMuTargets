#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <file_list.txt>"
    exit 1
fi

FILE_LIST=$1

echo "Checking staging status for files in $FILE_LIST..."
echo "------------------------------------------------"

while IFS= read -r FULL_PATH || [ -n "$FULL_PATH" ]; do
    # Skip empty lines
    [ -z "$FULL_PATH" ] && continue

    # Extract just the filename from the path
    FILENAME=$(basename "$FULL_PATH")

    # Use ifdh to find locations
    # We redirect stderr to /dev/null to keep the output clean
    LOCATIONS=$(ifdh locateFile "$FILENAME" 2>/dev/null)

    # Check if 'dcache' or 'fndca' exists in the location string
    if echo "$LOCATIONS" | grep -iqE "dcache|fndca|enstore:/pnfs/fnal.gov/usr/.*disk"; then
        echo "[ STAGED ] $FILENAME"
    elif echo "$LOCATIONS" | grep -iq "enstore"; then
        echo "[  TAPE  ] $FILENAME (Needs Prestage)"
    else
        echo "[ NOTFOUND ] $FILENAME (Check SAM declaration)"
    fi

done < "$FILE_LIST"