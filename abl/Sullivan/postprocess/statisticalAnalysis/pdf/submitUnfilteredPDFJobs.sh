for (( i=1; i<=50; i++ ))
do 
    cp unfilteredPDF.pbs temp.pbs 
    echo $i 
    sed "s/blah/$i/g" temp.pbs > temp1.pbs 
    mv temp1.pbs unfilteredPDF_${i}.pbs 
done
rm temp.pbs
cd /lustre/scratch/gantech/les_abl/Sullivan/data/sullivanTest/case037/statisticalAnalysis/pdf
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/statisticalAnalysis/pdf/unfilteredPDF_*.pbs .
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/statisticalAnalysis/pdf/unfilteredPDFCalc .
for (( i=1; i<=50; i++ ))
do 
    qsub unfilteredPDF_${i}.pbs
done
cd -