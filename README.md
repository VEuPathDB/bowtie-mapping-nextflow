# bowtieMapping

***<p align=center>bowtieMapping (From sra download)</p>***  
```mermaid
flowchart TD
    p0((Channel.fromList))
    p1[bowtieMapping:sra:downloadFiles]
    p2(( ))
    p3(( ))
    p4[bowtieMapping:sra:bowtie]
    p5[bowtieMapping:sra:PCRDuplicates]
    p6(( ))
    p0 -->|ids| p1
    p1 --> p4
    p2 -->|indexfiles| p4
    p3 -->|indexFileBasename| p4
    p4 --> p5
    p5 --> p6
```

***<p align=center>bowtieMapping (from local files)</p>***  
```mermaid
flowchart TD
    p0((Channel.fromFilePairs))
    p1(( ))
    p2(( ))
    p3[bowtieMapping:local:bowtie]
    p4[bowtieMapping:local:PCRDuplicates]
    p5(( ))
    p0 -->|files| p3
    p1 -->|indexfiles| p3
    p2 -->|indexFileBasename| p3
    p3 --> p4
    p4 --> p5
```

### Get Started
  * Install Nextflow
    
    `curl https://get.nextflow.io | bash`
  
  * Run the script
    
    `nextflow run VEuPathDB/blastSimilarity -with-trace -c  <config_file> -r main`
