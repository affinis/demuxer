#!/bin/bash

set -e
set -u


HELP=`cat <<EOF

sub -n task_name -c 'command'

EOF
`

while getopts hn:c: options
do
	case $options in
	h)
		echo "$HELP"
		exit 0
		;;
	n)
		TASKNAME=$OPTARG
		;;
	c)
		COMMAND=$OPTARG
		;;
	*)
		echo "$HELP"
		exit 1
		;;
	esac
done

COMMADEND=' && exit\n'
COMMAND_EXE=$COMMAND$COMMADEND

echo "Excuting $COMMAND_EXE in screen $TASKNAME"

screen -dmS $TASKNAME && screen -S $TASKNAME -X stuff "$COMMAND_EXE"
