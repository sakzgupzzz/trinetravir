variable "aws_region" {
  description = "AWS region"
  type        = string
  default     = "us-east-1"
}

variable "instance_type" {
  description = "EC2 GPU instance type. g5.xlarge = A10G 24GB, $1.006/hr on-demand."
  type        = string
  default     = "g5.xlarge"
}

variable "key_name" {
  description = "Existing EC2 key pair name (re-used from StockBeat infra)."
  type        = string
  default     = "stockbeat-key"
}

variable "ssh_allowed_cidr" {
  description = "CIDR allowed to SSH. Must be your IP /32, not 0.0.0.0/0."
  type        = string

  validation {
    condition     = var.ssh_allowed_cidr != "0.0.0.0/0"
    error_message = "ssh_allowed_cidr must not be 0.0.0.0/0. Restrict to your IP (e.g. 203.0.113.42/32)."
  }
}

variable "root_volume_size_gb" {
  description = "Root EBS volume size (GB). gp3 default. Sized for ~108MB scVI inputs + ~5-10GB scVI checkpoints + repo + venv."
  type        = number
  default     = 100
}

variable "use_spot" {
  description = "Use spot pricing (~70% discount, risk of eviction). Acceptable for ephemeral 10-15h workload."
  type        = bool
  default     = false
}

variable "spot_max_price" {
  description = "Maximum spot price per hour. Set above current spot rate to avoid eviction; below on-demand to save."
  type        = string
  default     = "0.40"
}
