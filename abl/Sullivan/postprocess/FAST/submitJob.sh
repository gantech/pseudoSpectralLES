for (( i=1; i<=15; i++ ))
do 
    cp fast.slurm temp.slurm 
    echo $i 
    sed "s/blah/$i/g" temp.slurm > temp1.slurm 
    mv temp1.slurm fast_${i}_${i}.slurm 
done
cp dt /scratch/02504/ganesh10/ablFromKrakenAnalysis/FAST/
cd /scratch/02504/ganesh10/ablFromKrakenAnalysis/FAST/
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/FASTfilesWrite .
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/fast_*.slurm .
for (( i=1; i<=15; i++ ))
do 
    sbatch fast_${i}_${i}.slurm
done
cd -
