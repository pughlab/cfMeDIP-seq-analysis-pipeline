#!/bin/bash

METRIC=$1

cat $2 \
    | grep -A10 "## METRICS CLASS" \
    | grep -v "## METRICS CLASS" \
    | sed '/## HISTOGRAM/,+10d' \
    | awk -F "\t" '{for (i=1; i<=NF; i++) a[i,NR]=$i
                    max=(max<NF?NF:max)}
                    END {for (i=1; i<=max; i++)
                            {for (j=1; j<=NR; j++) 
                             printf "%s%s", a[i,j], (j==NR?RS:FS)
                            }
                        }' \
    | grep -v "^SAMPLE\|^LIBRARY\|^READ_GROUP" \
    | grep -v "ACCUMULATION_LEVEL" \
    | awk -v METRIC=$METRIC 'BEGIN {OFS="\t"} $1=METRIC"\t"$1' \
    > $3

## EOF
