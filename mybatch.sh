#!/bin/bash

NODES=$1
PRO=$2

cat > yhrun.sh << EOF
#!/bin/bash

yhrun -pthcp1 -N $NODES -n $PRO ./main
EOF

yhbatch -pthcp1 -N $NODES -n $PRO -o ./log-$PRO.out yhrun.sh

