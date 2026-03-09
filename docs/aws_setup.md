# AWS Batch Setup for SLE Methylation Pipeline

## Quick Summary

| Item | Detail |
|------|--------|
| Region | us-east-1 |
| S3 Bucket | `sle-methylation-pipeline` |
| Compute | Managed Spot, maxvCpus=128, scales to zero |
| Instance types | r5.2xlarge (alignment), r5.xlarge (dedup/ComBat), m5.xlarge (rest) |
| Docker images | 4 custom images on ECR |
| Estimated cost | $44-58 total for 11 samples (Spot) |
| Wall-clock time | 12-18 hours |

## Prerequisites

- AWS CLI configured (`aws configure`)
- Docker installed locally (to build images)
- Nextflow installed

## Step 1: IAM Roles

```bash
# Batch service role
cat > batch-trust-policy.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Principal": {"Service": "batch.amazonaws.com"},
    "Action": "sts:AssumeRole"
  }]
}
EOF

aws iam create-role \
  --role-name SLE-BatchServiceRole \
  --assume-role-policy-document file://batch-trust-policy.json

aws iam attach-role-policy \
  --role-name SLE-BatchServiceRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole

# ECS instance role (for compute nodes)
cat > ec2-trust-policy.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Principal": {"Service": "ec2.amazonaws.com"},
    "Action": "sts:AssumeRole"
  }]
}
EOF

aws iam create-role \
  --role-name SLE-ECSInstanceRole \
  --assume-role-policy-document file://ec2-trust-policy.json

aws iam attach-role-policy \
  --role-name SLE-ECSInstanceRole \
  --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess

aws iam attach-role-policy \
  --role-name SLE-ECSInstanceRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role

aws iam attach-role-policy \
  --role-name SLE-ECSInstanceRole \
  --policy-arn arn:aws:iam::aws:policy/CloudWatchLogsFullAccess

aws iam create-instance-profile --instance-profile-name SLE-ECSInstanceProfile
aws iam add-role-to-instance-profile \
  --instance-profile-name SLE-ECSInstanceProfile \
  --role-name SLE-ECSInstanceRole
```

## Step 2: S3 Bucket

```bash
aws s3 mb s3://sle-methylation-pipeline --region us-east-1

# Upload reference genome
aws s3 cp GRCh38.primary_assembly.genome.fa s3://sle-methylation-pipeline/input/
aws s3 cp GRCh38.primary_assembly.genome.fa.fai s3://sle-methylation-pipeline/input/

# Upload sample sheet
aws s3 cp samples_full.csv s3://sle-methylation-pipeline/input/

# Auto-delete work dir after 7 days
cat > lifecycle.json << 'EOF'
{
  "Rules": [{
    "ID": "CleanWorkDir",
    "Filter": {"Prefix": "work/"},
    "Status": "Enabled",
    "Expiration": {"Days": 7}
  }]
}
EOF

aws s3api put-bucket-lifecycle-configuration \
  --bucket sle-methylation-pipeline \
  --lifecycle-configuration file://lifecycle.json
```

## Step 3: Batch Compute Environment

```bash
DEFAULT_VPC=$(aws ec2 describe-vpcs --filters Name=isDefault,Values=true \
  --query 'Vpcs[0].VpcId' --output text)

SUBNETS=$(aws ec2 describe-subnets --filters Name=vpc-id,Values=$DEFAULT_VPC \
  --query 'Subnets[*].SubnetId' --output text | tr '\t' ',')

SG=$(aws ec2 describe-security-groups \
  --filters Name=vpc-id,Values=$DEFAULT_VPC Name=group-name,Values=default \
  --query 'SecurityGroups[0].GroupId' --output text)

ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

# Create launch template with 500 GB EBS (default 30 GB is too small for WGBS)
cat > launch-template.json << 'EOF'
{
  "LaunchTemplateName": "sle-batch-500gb",
  "LaunchTemplateData": {
    "BlockDeviceMappings": [{
      "DeviceName": "/dev/xvda",
      "Ebs": {
        "VolumeSize": 500,
        "VolumeType": "gp3",
        "Encrypted": true
      }
    }]
  }
}
EOF

aws ec2 create-launch-template --cli-input-json file://launch-template.json

# Create Spot compute environment
aws batch create-compute-environment \
  --compute-environment-name sle-pipeline-spot \
  --type MANAGED \
  --compute-resources '{
    "type": "SPOT",
    "bidPercentage": 80,
    "minvCpus": 0,
    "maxvCpus": 128,
    "desiredvCpus": 0,
    "instanceTypes": ["r5.xlarge","r5.2xlarge","r5a.xlarge","r5a.2xlarge","m5.xlarge","m5.2xlarge","m5a.xlarge"],
    "subnets": ["'$SUBNETS'"],
    "securityGroupIds": ["'$SG'"],
    "instanceRole": "arn:aws:iam::'$ACCOUNT_ID':instance-profile/SLE-ECSInstanceProfile",
    "launchTemplate": {"launchTemplateName": "sle-batch-500gb"},
    "ec2Configuration": [{"imageType": "ECS_AL2"}]
  }' \
  --service-role "arn:aws:iam::${ACCOUNT_ID}:role/SLE-BatchServiceRole"

# Create job queue
aws batch create-job-queue \
  --job-queue-name sle-pipeline-queue \
  --state ENABLED \
  --priority 1 \
  --compute-environment-order order=1,computeEnvironment=sle-pipeline-spot
```

## Step 4: Build and Push Docker Images

```bash
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
ECR_URI="${ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com"

# Create ECR repos
for repo in sle-sra-tools sle-bwameth sle-r-methylation sle-python-ml; do
  aws ecr create-repository --repository-name $repo --region us-east-1
done

# Login to ECR
aws ecr get-login-password --region us-east-1 | \
  docker login --username AWS --password-stdin $ECR_URI

# Build from containers/ directory
cd containers/

# R methylation image (has pre-built ComBatMet)
docker build -f r_methylation/Dockerfile -t sle-r-methylation .
docker tag sle-r-methylation:latest ${ECR_URI}/sle-r-methylation:latest
docker push ${ECR_URI}/sle-r-methylation:latest

# Build other images similarly (Dockerfiles TBD)
cd ..
```

## Step 5: Run the Pipeline

```bash
# Full 11-sample run
nextflow run main.nf \
  -profile aws \
  --sample_sheet s3://sle-methylation-pipeline/input/samples_full.csv \
  --genome s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa \
  --genome_fai s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa.fai \
  --outdir s3://sle-methylation-pipeline/results

# Resume after interruption
nextflow run main.nf \
  -profile aws \
  --sample_sheet s3://sle-methylation-pipeline/input/samples_full.csv \
  --genome s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa \
  --genome_fai s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa.fai \
  --outdir s3://sle-methylation-pipeline/results \
  -resume
```

## Step 6: Download Results

```bash
aws s3 sync s3://sle-methylation-pipeline/results ./results

# Clean up work dir (or let lifecycle policy handle it)
aws s3 rm s3://sle-methylation-pipeline/work --recursive
```

## Cost Breakdown

| Category | Estimate |
|----------|----------|
| Compute (Spot, 11 samples) | $14-23 |
| S3 storage (work dir, 1 week) | $12 |
| S3 storage (results, ongoing) | $4/month |
| Data transfer (download results) | ~$15 |
| **Total for one run** | **$44-58** |

## Tips

- **SRA downloads are free within us-east-1** — NCBI has cloud mirrors there
- **Run Nextflow head on a t3.medium EC2** (~$0.04/hr) instead of your laptop for reliability
- **Spot interruptions auto-retry** — errorStrategy already handles exit code 143
- **bwameth index builds once** — cached via `-resume`, saves 6 hours on subsequent runs
- **Every container needs `aws` CLI** — Nextflow uses it for S3 file staging
