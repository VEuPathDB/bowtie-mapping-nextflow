#!/usr/bin/env bash

set -euo pipefail

echo $HOME
if [ -s $HOME/.ncbi/user-settings.mkfg ]
then
    fasterq-dump --split-3 ${id}
else
    cp /root/.ncbi/user-settings/mkfg $HOME/.ncbi/user-setting.mkfg
    fasterq-dump --split-3 ${id}
fi
