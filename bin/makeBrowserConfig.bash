#!/bin/bash

cat <<EOF
[$2]
:selected    = 1
display_name = $1
sample       = $1
alignment    = $3
type         = $4

EOF
