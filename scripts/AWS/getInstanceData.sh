#!/bin/bash

curl http://169.254.169.254/latest/meta-data/

#To set up some useful env vars:

EC2_INSTANCE_ID="`wget -q -O - http://169.254.169.254/latest/meta-data/instance-id || die \"wget instance-id has failed: $?\"`"
test -n "$EC2_INSTANCE_ID" || die 'cannot obtain instance-id'
EC2_AVAIL_ZONE="`wget -q -O - http://169.254.169.254/latest/meta-data/placement/availability-zone || die \"wget availability-zone has failed: $?\"`"
test -n "$EC2_AVAIL_ZONE" || die 'cannot obtain availability-zone'
EC2_REGION="`echo \"$EC2_AVAIL_ZONE\" | sed -e 's:\([0-9][0-9]*\)[a-z]*\$:\\1:'`"

# curl http://169.254.169.254/latest/meta-data/

# Case 1: To mount ephemeral storage:
#SCRATCH_DEVICE0="/dev/`wget -q -O - http://169.254.169.254/latest/meta-data/block-device-mapping/ephemeral0 || die \"wget ephemeral0 has failed: $?\"`"
#sudo mkfs -t ext4 $SCRATCH_DEVICE0

#SCRATCH_DEVICE1="/dev/`wget -q -O - http://169.254.169.254/latest/meta-data/block-device-mapping/ephemeral1 || die \"wget ephemeral1 has failed: $?\"`"
#sudo mkfs -t ext4 $SCRATCH_DEVICE1
SCRATCH_DEVICE0=/dev/xvdf
SCRATCH_DEVICE1=/dev/xvdg

sudo mkfs.ext4 $SCRATCH_DEVICE0
sudo mkfs.ext4 $SCRATCH_DEVICE1

sudo mount $SCRATCH_DEVICE0 /scratch1
sudo mount $SCRATCH_DEVICE1 /scratch2

# Case 2: To mount ephemeral storage as 1 large volume:
SCRATCH_DEVICE0=/dev/xvdf
SCRATCH_DEVICE1=/dev/xvdg
sudo mdadm --create /dev/md0 --level=0 --chunk=256 --raid-devices=2 $SCRATCH_DEVICE0 $SCRATCH_DEVICE1
sudo mdadm --detail /dev/md0
# Increase block size for better performance
sudo blockdev --setra 65536 /dev/md0
# Create physical volume
sudo pvcreate /dev/md0
# Create volume group
sudo vgcreate vg0 /dev/md0
# Check
sudo vgdisplay vg0
# Create 1 logical volume sppaning the entire available space
sudo lvcreate --name scratch --size 305G vg0
# Create fs
sudo mkfs.ext4 /dev/vg0/scratch
# Mount
sudo mkdir /scratch
sudo mount /dev/vg0/scratch /scratch
sudo chmod go+w /scratch
