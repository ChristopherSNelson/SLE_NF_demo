#!/bin/bash
# Kill local Nextflow JVM and any child bioinformatics processes

kill -9 $(pgrep -f nextflow) 2>/dev/null && echo "Nextflow killed" || echo "No Nextflow running"

CHILDREN="bwameth|fastp|bwa|picard|samtools|MethylDackel"
if pgrep -f "$CHILDREN" > /dev/null 2>&1; then
  pkill -9 -f "$CHILDREN" && echo "Child processes killed"
else
  echo "No child processes"
fi
