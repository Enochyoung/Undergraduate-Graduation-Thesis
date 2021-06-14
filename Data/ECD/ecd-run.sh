#!/bin/bash
#SBATCH --job-name=ala_1_del
#SBATCH --output=ecd.out
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00
#SBATCH -p regular

# This script is specified for submitting TDPW OA/ECD task on the new server.

module load tdapw/intel20u4/6.6

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_PMI2=yes

mpirun="srun --mpi=pmi2"


# Create Ex Ey Ez beads
beads=1 # Ex Ey Ez
beadsi=$((beads-1)) # Index of Ek (0->x, 1->y, 2->z)
nspin=1 # Support nspin=1|4, 1 is default, 4 requires comprehensive understanding
		# of the magnetization property of the system of interest
Ek=( "0.002 0 0" "0 0.002 0" "0 0 0.002" ) # (Ex, Ey, Ez) in a.u. (Ry/(Bohr*e))
Exyz=F # If T, applying E-field at x, y, z directions simultaneously, considered
       # to be applying E-field at only one direction (111)
Bfield=( 2 2 2 ) # Not supported yet


# Specify method
Diagon=( E A Exp B )
method=LastDiagon${Diagon[0]} # Support LastDiagonE LastDiagonA LastDiagonExp
                              # LastDiagonE = $E\theta(-t)$
                              # LastDiagonA = $A\theta(-t)$ (Velocity Gauge)
                              # LastDiagonExp = $e^{-iEr}\psi$
#method=alwaysE               # alwaysE = $E\theta(t)$
#method=DFTE                  # DFTE = $E=kt,H\psi(t)=e\psi(t)$


# To restart calculation
echo "rm method pointer file"
rm -rf LastDiagonE LastDiagonA LastDiagonExp LastDiagonB alwaysE DFTE
echo "rm -rf input_* tmp_* result"
rm -rf input_* tmp_* result
touch $method # A method pointer blank file
cp ecd-run.sh ecd-run.sh.back # Backup this script


# E-field at only one direction
if [ "$Exyz" == "T" ]; then
	Ek=( "0.002 0.002 0.002" "0.002 0.002 0.002" "0.002 0.002 0.002" )
	beads=1
	beads0=$((beads-1))
	rm -rf tmp_1 tmp_2
	ln -s tmp_0 tmp_1
	ln -s tmp_0 tmp_2
fi


# To make sure $NP is not 0 and $NP can be divided evenly by $beads
NP=$(mpirun echo | wc -l)
if [ "$((NP/beads*beads))" == 0 ]; then
	echo "Warning: NP == 0"
	exit 1
else
	mpirun="mpirun -np $((NP/beads*beads))"
fi


# Input files preparation
WORKDIR=$PWD
input=$WORKDIR/input/input.in
xyz=(x y z) # 

for i in $(seq 0 1 ${beadsi}) # Start loop, one loop, one E-field direction
do
	# Make working file and input file
	cd $WORKDIR
    echo -e "cp -r $input input_$i.in"
    cp -r $input input_$i.in
	mkdir $WORKDIR/tmp_$i
	ln -s $WORKDIR/input_$i.in $WORKDIR/tmp_$i/input.in
	ln -s $WORKDIR/input_$i.out $WORKDIR/tmp_$i/result
	
	# Set method in input file
	if [ $method == LastDiagonE -o $method == alwaysE -o $method == DFTE ]; then
		sed -i "s/tefield/tefield=T ! /g" input_$i.in
	fi
	
	# Set LastDiagonB and nspin in input file
	if [ $method == LastDiagonB ]; then
		nspin=4
	fi
	if [ $nspin == 2 ]; then # spin-polarized calculation with magnetization along z axis, LSDA
		sed -i "s/!nspin/nspin=2 ! /g" input_$i.in
		sed -i "s/!starting_magnetization/starting_magnetization=0.0 ! /g" input_$i.in
	elif [ $nspin == 4 ]; then # spin-polarized calculation with magnetization in generic directon, noncollinear
		sed -i "s/!noncolin/noncolin=T ! /g" input_$i.in 
		sed -i "s/!lspinorb/lspinorb=T ! /g" input_$i.in 
	fi

	# Set Efield in input file and generate TDEFIELD.in
	cd $WORKDIR/tmp_$i/
	if [ $method == LastDiagonB ]; then
		rm -rf TDEFIELD
		bdirect=$((i+1)) # Index of Bk (1->x, 2->y, 3->z)
		sed -i "s/!B_field($bdirect)/B_field($bdirect)=${Bfield[$i]} ! /g" $WORKDIR/input_$i.in 
	elif [ ${method:0:10} == LastDiagon ]; then
		sed -i "s/$method[[:space:]]/$method=T ! /g" $WORKDIR/input_$i.in
		echo "Edelta.sh 100000 ${Ek[$i]} "
		Edelta.sh 100000 ${Ek[$i]} # Generate delta/\theta(-t) TDEFIELD.in	
	elif [ $method == alwaysE ]; then
		echo "alwaysE.sh 100000 ${Ek[$i]} "
	 	alwaysE.sh 100000 ${Ek[$i]} # Generate alwaysE TDEFIELD.in
	elif [ $method == DFTE ]; then
		echo "DFTE"
		sed -i "s/nstep/nstep=20 ! /g" $WORKDIR/input_$i.in
		sed -i "s/tddft_is_on/tddft_is_on=F ! /g" $WORKDIR/input_$i.in
		cp ../EFIELD/LinearE$i ./TDEFIELD.in
	else
		# An option to normally read evolution of E-field, you have to annotate the method statement
		sed -i "s/GaugeField/GaugeField=T ! /g" $WORKDIR/input_$i.in # GaugeField=T -> Velocity Gauge
		cp ../TDEFIELD.in ./TDEFIELD.in                              # GaugeField=F -> Length Gauge
	fi
	edir=$((i+1)); sed -i "s/edir/edir=${edir} ! /g" $WORKDIR/input_$i.in
done
cd $WORKDIR
echo "$mpirun tdpw.x -ni ${beads} -i input | tee result"
$mpirun tdpw.x -ni ${beads} -i input | tee result # Read file with prefix "input" as input file
exit 0