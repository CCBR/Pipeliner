#!/bin/bash

#id=viz-ccbr-mytest1001
id=$1
vcf=$2

/usr/local/bin/curl  https://avia-abcc.ncifcrf.gov/apps/site/upload_viz -X POST -F user_file=@$vcf -F user_inputformat=vcf4 -F user_api=cFMdtEdwm34iVzXOZ6 -F "user_email=justin.lack|nih.gov" --insecure -F user_id=$id
perl avia_check_status.pl $id
