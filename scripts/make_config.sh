module load emboss
module load hmmer/3
module load cdbfasta
module load fastx_toolkit
module load genewise
module load fasttree
module load muscle
module load trimal
module load phred-phrap-consed
module load blat
module load kent
module load fasta
module load cap3
module load exonerate
module load diamond
module load wu-blast
TEMPLATE=lib/apps.conf.template
CONF=lib/apps.conf
rm -f $CONF
for r in `sort $TEMPLATE`
do
 VAR=`echo $r | awk -F= '{print $1}'`
 VAL=`echo $r | awk -F= '{print $2}'`
 NEWVAL=`which $VAL`
 echo "$VAR=$NEWVAL" >> $CONF
done
