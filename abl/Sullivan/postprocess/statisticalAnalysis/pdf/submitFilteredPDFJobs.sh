for (( i=1; i<=50; i++ ))
do 
    cp filteredPDF.pbs temp.pbs 
    echo $i 
    sed "s/blah/$i/g" temp.pbs > temp1.pbs 
    mv temp1.pbs filteredPDF_${i}.pbs 
done
rm temp.pbs
cd /lustre/scratch/gantech/les_abl/Sullivan/data/sullivanTest/case037/statisticalAnalysis/pdf
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/statisticalAnalysis/pdf/filteredPDF_*.pbs .
cp ~/pseudoSpectralLES/abl/Sullivan/postprocess/statisticalAnalysis/pdf/filteredPDFCalc .
for (( i=1; i<=50; i++ ))
do 
    qsub filteredPDF_${i}.pbs
done
cd -
