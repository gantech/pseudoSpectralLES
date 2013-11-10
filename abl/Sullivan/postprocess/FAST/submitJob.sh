for (( i=1; i<=15; i++ ))
do 
    cp fast.pbs temp.pbs 
    echo $i 
    sed "s/blah/$i/g" temp.pbs > temp1.pbs 
    mv temp1.pbs fast_${i}_${i}.pbs 
done
cd /lustre/scratch/gantech/les_abl/Sullivan/data/sullivanTest/case037/FAST/
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/FASTfilesWrite .
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/fast_*.pbs .
for (( i=2; i<=15; i++ ))
do 
    qsub fast_${i}_${i}.pbs
done
cd -
