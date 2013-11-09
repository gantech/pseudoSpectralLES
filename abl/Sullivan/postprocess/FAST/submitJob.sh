cd /lustre/scratch/gantech/les_abl/Sullivan/data/sullivanTest/case037/FAST/
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/FASTfilesWrite .
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/fast.pbs .
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/run15Rows .
qsub fast.pbs
cd -
