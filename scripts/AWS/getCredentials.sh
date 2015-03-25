#!/bin/bash

aws sts assume-role --role-arn arn:aws:iam::022180670435:role/ag-ngs --role-session-name ngs  >/dev/shm/ngsrole.txt
export AWS_ACCESS_KEY=`egrep "AccessKeyId" /dev/shm/ngsrole.txt |sed 's/.*: \"\(.*\)\".*/\1/'`
export AWS_SECRET_ACCESS_KEY=`egrep "SecretAccessKey" /dev/shm/ngsrole.txt |sed 's/.*: \"\(.*\)\".*/\1/'`
export AWS_SECURITY_TOKEN=`egrep "SessionToken"  /dev/shm/ngsrole.txt |sed 's/.*: \"\(.*\)\".*/\1/'`

echo AWS_ACCESS_KEY=$AWS_ACCESS_KEY
echo AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY
echo AWS_SECURITY_TOKEN=$AWS_SECURITY_TOKEN
