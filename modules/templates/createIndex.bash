#!/usr/bin/env bash

set -euo pipefail

bowtie2-build $databaseFasta index
