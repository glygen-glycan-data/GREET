#!/bin/sh

# set -x

CPUS=$1
HOST=$2
PORT=$3
SECRET=$4
JOBS=$5
shift; shift; shift; shift; shift;

echo "ssh -4 -L :${PORT}:${HOST}:${PORT} -N -f nedwards@edwardslab.bmcb.georgetown.edu"

LOCALHOST=`uname -n`

echo ./greet.py $@ --workers "__${CPUS}:${LOCALHOST}:${PORT}:${SECRET}__"

sbatch --cpus-per-task "$CPUS" --output=/dev/null --array=1-"$JOBS" <<EOF
#!/bin/sh
srun ./greet.py $@ --workers "__${CPUS}:${LOCALHOST}:${PORT}:${SECRET}__"
EOF




