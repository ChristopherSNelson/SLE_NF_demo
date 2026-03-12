# AWS Batch Setup for SLE Methylation Pipeline

## Quick Summary

| Item | Detail |
|------|--------|
| Region | us-east-2 (Ohio) |
| S3 Bucket | `sle-methylation-pipeline` |
| Compute | Managed Spot, maxvCpus=128, scales to zero |
| Instance types | r5.2xlarge (alignment), r5.xlarge (dedup/ComBat), m5.xlarge (rest) |
| Docker images | 4 custom images on ECR |
| Estimated cost | $44-58 total for 11 samples (Spot) |
| Wall-clock time | 12-18 hours |

## Prerequisites

- AWS account with billing enabled (payment method required for EC2 instances)
- AWS CLI installed locally (`brew install awscli` on macOS)

## Step 0: AWS CLI + IAM Setup (brand new account)

**Create an IAM user with CLI access:**
1. AWS Console → IAM → Users → Create user
2. Attach policy: `AdministratorAccess` (simplify for initial setup)
3. Security credentials tab → Create access key → CLI use case → copy keys

```bash
aws configure
# AWS Access Key ID: <paste>
# AWS Secret Access Key: <paste>
# Default region: us-east-2
# Default output format: json

# Verify it works
aws sts get-caller-identity
```

## Steps 1-3: Automated Setup

The setup script creates all IAM roles, the S3 bucket, launch template, compute environment, and job queue in one shot. It is idempotent — safe to re-run.

```bash
# From ~/aws-setup/sle-pipeline/
chmod +x setup_aws_batch.sh
./setup_aws_batch.sh
```

See `setup_aws_batch.sh` for what it creates. To tear everything down:

```bash
chmod +x teardown_aws_batch.sh
./teardown_aws_batch.sh
```

### What the setup creates

| Resource | Name |
|----------|------|
| IAM role (Batch service) | `SLE-BatchServiceRole` |
| IAM role (ECS instances) | `SLE-ECSInstanceRole` |
| IAM role (Spot fleet) | `SLE-EC2FleetRole` |
| Instance profile | `SLE-ECSInstanceProfile` |
| S3 bucket | `sle-methylation-pipeline` |
| Launch template | `sle-batch-500gb` (500 GB gp3 EBS) |
| Compute environment | `sle-pipeline-spot` (Spot, maxvCpus=128) |
| Job queue | `sle-pipeline-queue` |

## Step 4: Build and Push Docker Images

```bash
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
ECR_URI="${ACCOUNT_ID}.dkr.ecr.us-east-2.amazonaws.com"

# Create ECR repos
for repo in sle-sra-tools sle-bwameth sle-r-methylation sle-python-ml; do
  aws ecr create-repository --repository-name $repo --region us-east-2
done

# Login to ECR
aws ecr get-login-password --region us-east-2 | \
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

## Step 5: Upload Data to S3

```bash
# Upload reference genome
aws s3 cp GRCh38.primary_assembly.genome.fa s3://sle-methylation-pipeline/input/
aws s3 cp GRCh38.primary_assembly.genome.fa.fai s3://sle-methylation-pipeline/input/

# Upload sample sheet
aws s3 cp sampleSheets/samples_full.csv s3://sle-methylation-pipeline/input/
```

## Step 6: Launch Nextflow Head Node (t3.medium EC2)

Run Nextflow from an EC2 instance rather than your laptop — avoids pipeline failure if your laptop sleeps or loses network.

```bash
# Get latest Amazon Linux 2023 AMI
AMI=$(aws ec2 describe-images --owners amazon \
  --filters "Name=name,Values=al2023-ami-*-x86_64" "Name=state,Values=available" \
  --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text \
  --region us-east-2)

# Create SSH key pair
aws ec2 create-key-pair --key-name sle-key --query 'KeyMaterial' --output text \
  --region us-east-2 > sle-key.pem
chmod 400 sle-key.pem

# Launch t3.medium (30 GB root volume is fine — no pipeline data processed here)
INSTANCE_ID=$(aws ec2 run-instances \
  --image-id $AMI --instance-type t3.medium --key-name sle-key --count 1 \
  --iam-instance-profile Name=SLE-ECSInstanceProfile \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=sle-nextflow-head}]' \
  --query 'Instances[0].InstanceId' --output text \
  --region us-east-2)

# Get public IP (wait ~30s for instance to start)
sleep 30
aws ec2 describe-instances --instance-ids $INSTANCE_ID \
  --query 'Reservations[0].Instances[0].PublicIpAddress' --output text \
  --region us-east-2
```

**SSH in and install dependencies:**

```bash
ssh -i sle-key.pem ec2-user@<IP>

# Java, Docker, git
sudo dnf install -y java-21-amazon-corretto docker git
sudo systemctl start docker
sudo usermod -aG docker ec2-user
newgrp docker

# Nextflow
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/

# Clone pipeline
git clone https://github.com/<you>/SLE_NF_demo.git
cd SLE_NF_demo
```

## Step 7: Run the Pipeline

Run inside a tmux session so the pipeline survives SSH disconnects:

```bash
tmux new -s pipeline

nextflow run main.nf \
  -profile aws \
  --sample_sheet s3://sle-methylation-pipeline/input/samples_full.csv \
  --genome s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa \
  --genome_fai s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa.fai \
  --outdir s3://sle-methylation-pipeline/results

# Detach without killing: Ctrl+B then D
# Reattach later: tmux attach -t pipeline

# Resume after interruption (same command + -resume)
nextflow run main.nf \
  -profile aws \
  --sample_sheet s3://sle-methylation-pipeline/input/samples_full.csv \
  --genome s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa \
  --genome_fai s3://sle-methylation-pipeline/input/GRCh38.primary_assembly.genome.fa.fai \
  --outdir s3://sle-methylation-pipeline/results \
  -resume
```

**Troubleshooting:** failed task logs are at `s3://sle-methylation-pipeline/work/xx/xxxx/`. Fetch with:
```bash
aws s3 cp s3://sle-methylation-pipeline/work/xx/xxxx/.command.err .
aws s3 cp s3://sle-methylation-pipeline/work/xx/xxxx/.command.log .
```
Nextflow prints the work path for each failed task in the terminal output.

## Step 8: Download Results

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

- **SRA downloads are free within us-east-2** — NCBI has cloud mirrors in commercial AWS regions
- **Run Nextflow head on a t3.medium EC2** (~$0.04/hr) instead of your laptop for reliability
- **Spot interruptions auto-retry** — errorStrategy already handles exit code 143
- **bwameth index builds once** — cached via `-resume`, saves 6 hours on subsequent runs
- **Every container needs `aws` CLI** — Nextflow uses it for S3 file staging
- **Billing must be enabled** — without a payment method, EC2 instances won't launch (jobs stay in RUNNABLE forever)
