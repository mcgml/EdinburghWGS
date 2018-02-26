#!/bin/bash
set -euo pipefail

#download files
~/.aspera/connect/bin/ascp \
-pv \
-P 33001 \
-O 33001 \
"$user"@transfer.epcc.ed.ac.uk:/ .