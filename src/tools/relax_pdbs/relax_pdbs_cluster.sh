mkdir /common/madsci/modeling/antibody_databases/PyIgClassify/DBOUT/relaxed_renum_pdbs_orbitals
qsub -q dna -l nodes=8:ppn=16 -V -N relax_orbitals -d $rosetta_logs -v np=100,nstruct=1,program=relax,flag=/common/madsci/modeling/antibody_databases/PyIgClassify/tools/relax_pdbs_cluster_orbitals.flag,debug_log=relax_renum_ab /common/madsci/run_rosetta_mpi.sh

mkdir /common/madsci/modeling/antibody_databases/PyIgClassify/DBOUT/relaxed_renum_pdbs_talaris
qsub -q dna -l nodes=7:ppn=16 -V -N relax_talaris -d $rosetta_logs -v np=100,nstruct=1,program=relax,flag=/common/madsci/modeling/antibody_databases/PyIgClassify/tools/relax_pdbs_cluster_talaris.flag,debug_log=relax_renum_ab /common/madsci/run_rosetta_mpi.sh
