# Trinetravir Session 4 GPU infrastructure.
#
# Ephemeral g5.xlarge EC2 instance for scVI sensitivity sweeps (Issue 6).
# Mirrors StockBeat infrastructure conventions:
#   - AWS profile "stockbeat" (shared account, same IAM credentials)
#   - Region us-east-1
#   - Key pair "stockbeat-key" (same SSH key user already has)
#   - SSH CIDR locked to user's IP per StockBeat security.tf pattern
#
# Differs from StockBeat prod stack:
#   - Local state (no S3 backend) — Session 4 is single-run ephemeral.
#   - Deep Learning AMI (not AL2023) — pre-installed CUDA + PyTorch.
#   - GPU instance type (g5.xlarge) instead of t3.small web server.
#   - No HTTP/HTTPS ingress; SSH only.
#   - No Route53 / SES / DynamoDB / IAM SSM role — pure compute.
#
# Lifecycle: provision → run Part A+B+C → terminate. Plan to destroy after
# Session 4 closes; cost incurred only during compute window (~10-15h total).

terraform {
  required_version = ">= 1.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }

  # Local state. Session 4 is ephemeral single-run; no team coordination needed.
  # If reused across sessions, migrate to S3 backend (pattern from StockBeat main.tf).
}

provider "aws" {
  region  = var.aws_region
  profile = "stockbeat"

  default_tags {
    tags = {
      Project   = "trinetravir"
      ManagedBy = "terraform"
      Session   = "session_4"
    }
  }
}
