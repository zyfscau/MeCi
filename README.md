# Manual of MeCi v1.1

MeCi (Mitochondrion-encoded circRNA identifier) is a novel circRNA identification algorithm which predicts circRNAs independent of intron back-splicing. It is efficient in detection of circRNAs encoded by small size genomes, especially plant mitochondrial and plastid genomes.

What can MeCi do? 
1. MeCi uses blastn to map the high-throughput reads from transcripteomics data to the reference genome.
2. MeCi could output circle junction reads which are independent of intron back-splicing events.
3. MeCi could output circRNAs which have multiple alignments on the reference genome. 
4. MeCi outputs cicRNAs which contain repetitive and non-encoded bases at circle junction.
5. MeCi outputs strand orientation of predicted circRNAs.
---------------------------------------------------------------------
## 1 Commands and arguments

How to run MeCi: 

```
perl run_MeCi.pl --in1 sample.R1.fq --in2 sample.R2.fq --genome mitochondria_genome.fa --chunksize 10000 --max_process 4 --Paired
```
The arguments of MeCi are as followings:  

    Input Options:  
    --in1           <string>    read1 input file name in fastq format, [required] 
    --in2           <string>    read2 input file name in fastq format, [required for PE input]  
    --genome        <string>    genome file in fasta format, [ required ] 
    --chunksize     <number>    chunk size for sequence split, default is 10000 [ optional ]  
    Advanced Options: 
    --blast_cpu     <number>    blast cpu number. default is 10  [ optional ] 
    --blast_evalue  <string>    blast evalue threshold. default is 2 [ optional ] 
    --max_process   <number>    max process number limitation for sequence alignment, default is 4 [ optional ] 
    --Paired        <null>      Paired-end mode. the input file is paired-end [ optional ]  
    --output_prefix <string>    Prefix of output circRNA ID.  Default is 'MT' [ optional ]  
    General Options:  
    --outdir        <string>    output dir ,default is 'out'  
    --help          <null>      display help information [ optional ] 
      
-----------------------------------------------
## 2 Preparation for using MeCi  
(1) Check the perl version and operation system.
	Before we start, please make sure you have installed parallel Perl 5.8 or higher and use Mac OS X or Linux operation system. 

(2) Check permission of software in 'bin/', it should be set permission as 755
	You could type as following in your terminal
```
	cd MeCi-v1.1/bin
	chmod 755 blastall fastalength fastqToFa flash formatdb
	cd ..
```
(3) Trim reads if necessary.
	This step is optional. For data sets with good sequencing quality, this step will not largely influence prediction of MeCi.

(4) Prepare fasta file of mitochondrial genome.
	This step is required. 

-------------------------------
3 An example of running MeCi

Please download MeCi-v1.0.tar.gz and MeCi-v1.0-testdata.tar.gz from https:// https://www.github.com/zyfscau/MeCi.  And then unzip MeCi-v1.0-testdata.tar.gz to get four files (test.R1.fq.gz and test.R2.fq.gz are fastq files for one sample, NC_007982.1.fa is mitochondrial genome of Zea mays, runMe.sh is the command to perform MeCi for this test data).

Enter the directory and type as following in your terminal:
```
sh runMe.sh
```
MeCi performs circRNA detecting within a few minutes, and you will see the following results in output directory:

### output file 1: circrna_details.xls

    circrna_id	chr	start	end	strand	count	reads
    MT|100978|101366|1|+	MT	100978	101366	+	1	E00454:607:H527CCCX2:1:1104:21014:13773
    MT|100979|101212|0|+	MT	100979	101212	+	1	E00454:607:H527CCCX2:1:1105:5558:32039
    MT|100980|101366|0|+	MT	100980	101366	+	1	E00454:607:H527CCCX2:1:1102:3569:16604
    MT|101305|102064|1|+	MT	101305	102064	+	1	E00454:607:H527CCCX2:1:1103:24779:31969
    MT|101743|101995|0|+	MT	101743	101995	+	1	E00454:607:H527CCCX2:1:1107:6380:2135
    MT|102979|103141|0|+	MT	102979	103141	+	1	E00454:607:H527CCCX2:1:1101:3894:73106
    MT|111279|111429|1|+	MT	111279	111429	+	2	E00454:607:H527CCCX2:1:1106:15656:31072,E00454:607:H527CCCX2:1:1106:15716:31353
    MT|111279|111429|2|+	MT	111279	111429	+	2	E00454:607:H527CCCX2:1:1103:25205:17307,E00454:607:H527CCCX2:1:1103:26859:15250
    MT|111279|111451|2|+	MT	111279	111451	+	1	E00454:607:H527CCCX2:1:1105:29772:7919
    MT|111279|111489|1|+	MT	111279	111489	+	1	E00454:607:H527CCCX2:1:1101:13920:47544
    
    ...
    MT|111284|111491|0|+	MT	111284	111491	+	1	E00454:607:H527CCCX2:1:1106:32096:52836
    MT|111285|111589|-3|+	MT	111285	111589	+	1	E00454:607:H527CCCX2:1:1106:15229:31881
    MT|125944|126365|1|+	MT	125944	126365	+	1	E00454:607:H527CCCX2:1:1105:14072:69045
    MT|125991|126215|2|+	MT	125991	126215	+	1	E00454:607:H527CCCX2:1:1106:19309:64580
    MT|128095|127882|1|-	MT	128095	127882	-	2	E00454:607:H527CCCX2:1:1104:15554:36979,E00454:607:H527CCCX2:1:1104:15585:36961
    
    ...

Columns of output file are split by tabs ("\t" in shell and perl).

Each column gives information of a predicted circRNA:

    Column 1: ID of a predicted circRNA in the pattern of "reference genome|start position|end position|repetitive or non-encoded bases at circle junction|strand orientation";
    Column 2: reference genome or chromosome of a predicted circRNA;
    Column 3: start location of a predicted circRNA on the reference genome or chromosome;
    Column 4: end location of a predicted circRNA on the reference genome or chromosome;
    Column 5: strand orientation of a predicted circRNA;
    Column 6: junction read count of a predicted circRNA;
    Column 7: supporting junction read ID (split by ",").

PS:

The circRNA ID is composed of 5 parts: reference genome or chromosome, start location, end location, number of repetitive and non-encoded bases at circle junction "Z", and strand orientation.

    If number Z = 0, it means the circRNA have no overlap or gap at circle junction and it has a determined location.
    If number Z > 0, it means the circRNA has repetitive bases at circle junction, and the circularization site is uncertain.=
    For example, MT|128095|127882|1|-, the circRNA is one of the following two locations:
    128094..127882 (5'..3' end), and 128095..127883
   
If number Z < 0, it means there are non-encoded bases at circle site, and the position of the circRNA is unambiguous.

      For example, MT|111285|111589|-3|+, 
      the 5' and 3' ends of this circRNA are located at 111285 and 111589, respectively, and there are 3 non-encoded bases between start and end positions of the circularization site.

### output file 2: circrna_details.xls
    read_id	strand	query1_start	query1_end	hit1_start	hit1_end	query2_start	query2_end	hit2_start	hit2_end	circrna_id	read_type
    E00454:607:H527CCCX2:1:1101:10044:50428	+	1	24	400235	400258	25	175	400408	400558	MT|400235|400558|0|+	unique
    E00454:607:H527CCCX2:1:1101:10429:41884	-	1	104	276855	276958	105	125	276815	276835	MT|276958|276815|0|-	unique
    E00454:607:H527CCCX2:1:1101:10551:44064	+	1	17	400877	400893	18	152	401514	401648	MT|400877|401648|0|+	unique
    E00454:607:H527CCCX2:1:1101:10561:57864	+	1	24	400236	400259	25	215	400485	400675	MT|400236|400675|0|+	unique
    E00454:607:H527CCCX2:1:1101:10703:28540	+	1	111	402100	402210	111	150	402421	402460	MT|402100|402460|1|+	unique
    E00454:607:H527CCCX2:1:1101:10703:40741	+	1	70	400234	400303	70	286	400665	400881	MT|400234|400881|1|+	unique
    E00454:607:H527CCCX2:1:1101:10724:10046	+	1	175	402102	402276	175	260	402380	402465	MT|402102|402465|1|+	unique
    E00454:607:H527CCCX2:1:1101:10774:44486	+	1	129	402105	402233	130	210	402372	402452	MT|402105|402452|0|+	unique
    E00454:607:H527CCCX2:1:1101:10774:72192	-	1	150	243049	243198	151	180	242952	242981	MT|243198|242952|0|-	unique
    E00454:607:H527CCCX2:1:1101:10784:72315	-	1	118	442478	442595	119	147	442336	442364	MT|442595|442336|0|-	unique
    
    ...
    E00454:607:H527CCCX2:1:1101:14732:37454	+	1	76	402098	402173	75	116	402426	402467	MT|402098|402467|2|+	unique
    E00454:607:H527CCCX2:1:1101:14763:34975	-	1	61	454999	455059	62	95	454771	454804	MT|455059|454771|0|-	multiple
    E00454:607:H527CCCX2:1:1101:14763:34975	-	1	61	523936	523996	62	95	523708	523741	MT|523996|523708|0|-	multiple
    E00454:607:H527CCCX2:1:1101:14793:52748	+	1	86	172572	172657	87	145	172679	172737	MT|172572|172737|0|+	multiple
    E00454:607:H527CCCX2:1:1101:14793:52748	+	1	86	71952	72037	87	145	72059	72117	MT|71952|72117|0|+	multiple
    
    ...

Columns of output file are split by tabs ("\t" in shell and perl).
Each column gives information of a read which support at least one circRNA:

### two segments of a circle junction read have different mapped regions on the reference genome.
    Column  1: read ID;
    Column  2: strand orietation of the read mapped on the reference genome;
    Column  3: start position of 1st segment on the read;
    Column  4: end position of 1st segment on the read;
    Column  5: start position of 1st mapped region on the genome;
    Column  6: end position of 1st mapped region on the genome;
    Column  7: start position of 2nd segment on the read;
    Column  8: end position of 2nd segment on the read;
    Column  9: start position of 2nd mapped region on the genome;
    Column 10: end position of 2nd mapped region on the genome;
    Column 11: circRNA ID supported by the junction read;
    Column 12: unique or multiple alignment of the junction read on the reference genome.


##4 Citations:

If you use MeCi, please cite the following papers:
