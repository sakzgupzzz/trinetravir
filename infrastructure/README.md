# Trinetravir Session 4 GPU infrastructure

Ephemeral AWS g5.xlarge EC2 instance for Session 4 scVI sensitivity sweeps (Issue 6 closure). Mirrors StockBeat infrastructure conventions (AWS profile `stockbeat`, region `us-east-1`, key pair `stockbeat-key`, SSH CIDR locked to user's IP). Separate Terraform state (local), separate security group, no shared resources with StockBeat prod stack.

## Quick reference

| Item | Value |
|------|-------|
| AWS profile | `stockbeat` |
| Region | `us-east-1` |
| Instance type | `g5.xlarge` (A10G 24GB, 4 vCPU, 16GB RAM) |
| AMI | Latest Deep Learning OSS Nvidia Driver AMI GPU PyTorch 2.x Ubuntu 22.04 |
| Key pair | `stockbeat-key` (re-used from StockBeat) |
| Root volume | 100GB gp3 |
| SSH ingress | `67.244.121.160/32` only (user IP from StockBeat terraform.tfvars) |
| Cost (on-demand) | $1.006/hr |
| Cost (spot) | ~$0.30/hr (up to $0.40 max) |
| Expected wall-time | 10-15h (Parts A+B+C combined) |
| Expected total cost | $12-15 on-demand; $4-5 spot |

## Workflow

```bash
# 1. Configure AWS profile (one-time, on local laptop)
aws configure --profile stockbeat
# Or verify: aws sts get-caller-identity --profile stockbeat

# 2. Provision GPU instance
cd infrastructure/
terraform init
terraform plan
terraform apply
# Note the public_dns + ssh_command from outputs

# 3. SSH in + bootstrap (~10 min)
$(terraform output -raw ssh_command)
# Inside EC2:
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
git clone https://github.com/sakzgupzzz/trinetravir.git
cd trinetravir
uv sync
uv run python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
# Expect: CUDA: True NVIDIA A10G

# 4. SCP data from local (NEW local terminal)
cd /Users/sakshamgupta/Documents/coding_projects/trinetravir
$(terraform -chdir=infrastructure output -raw scp_command_inputs)

# 5. Run Part A (back in ssh session, ~10-12h on A10G)
cd /home/ubuntu/trinetravir
nohup uv run python scripts/session4_part_a_scvi_sweep.py \
  > data/processed/session4_part_a.log 2>&1 &
tail -f data/processed/session4_part_a.log

# 6. SCP results back (local terminal, when job finishes)
scp -i ~/.ssh/stockbeat-key.pem 'ubuntu@<PUBLIC_DNS>:/home/ubuntu/trinetravir/results/tables/session4_*.csv' results/tables/
scp -i ~/.ssh/stockbeat-key.pem 'ubuntu@<PUBLIC_DNS>:/home/ubuntu/trinetravir/data/processed/session4_part_a.log' data/processed/

# 7. CRITICAL: destroy instance (AWS keeps charging until terminated)
cd infrastructure/
terraform destroy
```

## File layout

- `main.tf`             — Terraform + AWS provider config (profile `stockbeat`, default tags).
- `variables.tf`        — Inputs: region, instance type, key name, SSH CIDR, root volume, spot flag.
- `ec2.tf`              — Deep Learning AMI data source + g5.xlarge instance + optional spot config.
- `security.tf`         — Security group, SSH only from user IP.
- `outputs.tf`          — public_dns, ssh_command, scp_command_inputs, AMI ID, cost estimate.
- `terraform.tfvars`    — Pinned variable values (matches StockBeat: stockbeat-key, 67.244.121.160/32).
- `terraform.tfvars.example` — Template for future v1.5 / Session N runs.

## Key pair location

`~/.ssh/stockbeat-key.pem` (re-used from StockBeat). If not in `~/.ssh/`, find via `find ~ -name "stockbeat-key.pem" 2>/dev/null` and `chmod 400` before first SSH.

## State management

Local `terraform.tfstate`. Single-run ephemeral. If reusing for v1.5 Session N or Phase 5+ work, migrate to S3 backend per StockBeat `main.tf` pattern:

```hcl
backend "s3" {
  bucket         = "stockbeat-terraform-state"  # shared bucket
  key            = "trinetravir/session_4/terraform.tfstate"
  region         = "us-east-1"
  profile        = "stockbeat"
  dynamodb_table = "stockbeat-terraform-lock"
  encrypt        = true
}
```

## Cost guardrails

- 100GB gp3 storage = $0.08/GB-month → ~$0.27 for one day. Negligible.
- Data transfer out for results (~20MB CSVs + ~1MB logs) = ~$0. First 100GB/month free.
- Spot eviction risk: if instance evicted mid-Part-A, must re-launch + re-run from scratch. scVI checkpoints not auto-persisted. Mitigation: run on-demand for first attempt; switch to spot for Parts B+C if A succeeded.
- **Issue 35 $100 cap reference**: Session 4 total budget projected $12-15 on-demand or $4-5 spot. Well under cap.

## Destroy protocol

Always run `terraform destroy` from `infrastructure/` after Session 4 closes. Verify no orphaned EBS volumes or EIPs in AWS console. Spot instances terminate automatically on eviction; on-demand requires explicit `terraform destroy` or AWS Console termination.
