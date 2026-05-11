resource "aws_security_group" "trinetravir_gpu" {
  name        = "trinetravir-session-4-sg"
  description = "Trinetravir Session 4 GPU instance - SSH only from user IP"

  # SSH only; no public HTTP/HTTPS (pure compute job).
  ingress {
    description = "SSH"
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = [var.ssh_allowed_cidr]
  }

  # All outbound (downloads, uv sync, scp out, etc.)
  egress {
    description = "All outbound"
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = "trinetravir-session-4-sg"
  }
}
