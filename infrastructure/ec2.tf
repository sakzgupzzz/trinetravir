# Latest AWS Deep Learning AMI (Ubuntu 22.04, PyTorch GPU, x86_64)
# Pre-installed: CUDA 12.x, cuDNN, NVIDIA drivers, PyTorch, conda environments.
# Saves ~30 min of CUDA install vs AL2023.
data "aws_ami" "deep_learning" {
  most_recent = true
  owners      = ["amazon"]

  filter {
    name   = "name"
    values = ["Deep Learning OSS Nvidia Driver AMI GPU PyTorch 2.*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }

  filter {
    name   = "architecture"
    values = ["x86_64"]
  }
}

resource "aws_instance" "trinetravir_gpu" {
  ami                    = data.aws_ami.deep_learning.id
  instance_type          = var.instance_type
  key_name               = var.key_name
  vpc_security_group_ids = [aws_security_group.trinetravir_gpu.id]

  # Optional spot pricing per use_spot flag
  dynamic "instance_market_options" {
    for_each = var.use_spot ? [1] : []
    content {
      market_type = "spot"
      spot_options {
        max_price          = var.spot_max_price
        spot_instance_type = "one-time"
        # No persistent request — single ephemeral run
      }
    }
  }

  root_block_device {
    volume_size           = var.root_volume_size_gb
    volume_type           = "gp3"
    delete_on_termination = true
  }

  tags = {
    Name    = "trinetravir-session-4-gpu"
    Purpose = "Issue 6 scVI sensitivity sweep"
  }
}
