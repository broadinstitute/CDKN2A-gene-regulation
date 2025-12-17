## Some Juicer commands for the microC
## I've already aligned and created all the "merged_nodups" or "mnd" files via UGER
##  ~/juiceDir/scripts/juicer.sh -g hg38 -e -D ~/juiceDir/
## These contain, on each line, one Hi-C contact
## These are all in the gs://broad-p16-calico/microC folder
##
## Juicer tools jar: https://github.com/aidenlab/Juicebox/releases 
## Juicer tools documentation: https://github.com/aidenlab/juicer/wiki/HiCCUPS


## Usually replicates are combined in Hi-C analysis:
## This command takes in two replicates of the ESC and makes it into one
## Juicer expects the input text file to be sorted by chromosome block (i.e. all
## the chromosome9-chromosome9 reads together) so this command combines the two 
## replicates and writes to a new "mnd" file
sort -m -k2,2d -k6,6d  <(gunzip -c ./microC_ESC_Rep1.mnd.txt.gz) <(gunzip -c ./microC_ESC_Rep2.mnd.txt.gz) > ESC.mnd.txt

## This command makes a hic file from the mnd. We do it a bit differently in order to
## force normalization to the resolution we want; so we first create the hic file without
## normalization (-n). The other flags are: -q 1 (take only reads with mapping quality 1 and above)
## -c chr9 (only do chromosome 9) -d (only do intra chromosomal reads) and then the
## resolutions: 100bp to 25000bp. Then the input file, output file, and genome assembly
java -Xmx10g -jar juicer_tools.2.20.00.jar pre -n -q 1 -c chr9 -d -r 100,250,500,1000,2500,5000,10000,25000 ESC.mnd.txt ESC.hic hg38

## This command adds normalization. We ask for SCALE normalization. You could also ask for
## KR (should be the same thing but with worse convergence, but you could try), VC
## (vanilla coverage, just divides by the row sums, so has an oscillating effect), VC_SQRT 
## (used as an estimate at low resolution for SCALE/KR but it doesn't work well at high
## resolution)
## -r 500 is the MINIMUM resolution to normalize to; can do 100 to go that low.
## to add more than SCALE do -k SCALE,KR (for example)
## The -d again tells Juicer to only look at intrachromosomal reads though it may not matter here
java -jar juicer_tools.2.20.00.jar addNorm -r 500 -k SCALE -d ESC.hic 

## This command runs hiccups, the loop caller. This is the first command that requires GPUs
## -c chr9 is the chromosome, -r 1000,5000 is the resolutions
## You could go higher res (500, 200), but it doesn't work that well
## It will write files to the ESC_hiccups directory; merged_loops.bedpe is the final loop calls 
java -Xmx10g -jar juicer_tools.2.20.00.jar hiccups -c chr9 -r 1000,5000 ESC.hic ESC_hiccups

## This is the hiccups diff command - differential loop calling. It also requires GPUs
## and expects you've already run hiccups on your two files
## here you have to specify both the resolutions and the normalizations you want to use for each
## resolution. Call with the two hic files, and then their respective loop calls; it will be
## written out to the final argument
java -Xmx10g -jar juicer_tools.2.20.00.jar hiccupsdiff -c chr9 -r 1000,5000 -k SCALE,SCALE ESC.hic Wi21.hic ESC_hiccups/merged_loops.bedpe Wi21_hiccups/merged_loops.bedpe hiccupsdiff_ESC_Wi23
