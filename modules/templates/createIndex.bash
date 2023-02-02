#!/usr/bin/env bash

set -euo pipefail

if [ "$isColorSpace" = true ]; then
    bowtie2-build $databaseFasta index
else
    bowtie-build $databaseFasta index
if
