#To create the test output files in this directory I ran the following commands:

cd sandbox/testInSilico10

java -jar ../../gnw-3.1.2b.jar --simulate -c ../settings.txt --input-net ../InSilicoSize10-Yeast1.xml 
--output-path .

#A lot of output is created. A good example for time-series data is: InSilicoSize10-Yeast1_dream4_timeseries.tsv. This test case has a network of 10 genes. Time-series are created
#for 10 different perturbations. At t=500, the perturbation is removed and the network relaxes back to the wild-type.