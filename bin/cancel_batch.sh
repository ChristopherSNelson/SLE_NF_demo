#!/bin/bash
# Cancel all non-terminal AWS Batch jobs across pipeline queues
# Usage: cancel_batch.sh [queue]   (default: both queues)

REGION=us-east-2
QUEUES=${1:-"sle-pipeline-queue sle-test-queue"}

for Q in $QUEUES; do
  for S in SUBMITTED PENDING RUNNABLE STARTING RUNNING; do
    IDS=$(aws batch list-jobs --job-queue $Q --job-status $S --region $REGION \
          --query 'jobSummaryList[].jobId' --output text 2>/dev/null)
    for ID in $IDS; do
      echo "Terminating $ID ($Q / $S)"
      aws batch terminate-job --job-id $ID --reason "Manual cancel" --region $REGION
    done
  done
done
echo "Done"
