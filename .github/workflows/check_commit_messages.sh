#!/bin/bash

BASE_REF=$1
REF=$2

for commit in $(git rev-list $BASE_REF..$REF)
do
    echo "Checking commit $commit"
    git cat-file commit $commit | sed '1,/^$/d' | python githooks/commit-msg
done
