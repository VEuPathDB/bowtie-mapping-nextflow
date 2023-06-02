#!/usr/bin/env bash

set -euo pipefail

echo $HOME
if [ -s $HOME/.ncbi/user-settings.mkfg ]
then
    fasterq-dump --split-3 ${id}
else
    mkdir -p $HOME
    cp /usr/bin/user-settings.mkfg $HOME/.ncbi/
    fasterq-dump --split-3 ${id}
fi
