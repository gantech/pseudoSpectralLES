cd /lustre/scratch/gantech/les_abl/Sullivan/data/sullivanTest/case037/FAST/
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/FASTfilesWriteCWF1p0 .
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/FAST/fast_CWF1p0.pbs .
qsub fast_CWF1p0.pbs
cd -
