output "instance_id" {
  description = "EC2 instance ID"
  value       = aws_instance.trinetravir_gpu.id
}

output "public_dns" {
  description = "Public DNS — use for SSH"
  value       = aws_instance.trinetravir_gpu.public_dns
}

output "public_ip" {
  description = "Public IPv4 address"
  value       = aws_instance.trinetravir_gpu.public_ip
}

output "ssh_command" {
  description = "Ready-to-paste SSH command"
  value       = "ssh -i ~/.ssh/stockbeat-key.pem ubuntu@${aws_instance.trinetravir_gpu.public_dns}"
}

output "scp_command_inputs" {
  description = "Ready-to-paste scp command for scvi input h5ads"
  value       = "scp -i ~/.ssh/stockbeat-key.pem data/processed/scvi_input_*.h5ad data/processed/phase3_response_vectors_*.parquet ubuntu@${aws_instance.trinetravir_gpu.public_dns}:/home/ubuntu/trinetravir/data/processed/"
}

output "ami_id" {
  description = "Deep Learning AMI ID used"
  value       = data.aws_ami.deep_learning.id
}

output "ami_name" {
  description = "Deep Learning AMI name"
  value       = data.aws_ami.deep_learning.name
}

output "estimated_hourly_cost_usd" {
  description = "Estimated on-demand or spot hourly cost"
  value       = var.use_spot ? "spot up to $${var.spot_max_price}/hr (current ~$0.30/hr)" : "$1.006/hr on-demand (g5.xlarge)"
}
