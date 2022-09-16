# RNASeq pipeline 
## Author: Xpan
## Data: 2019-4-18
### Multiple Version Collection; The newest will be marked with "#^N"
#### Fastq file QC and data filter
#### trim_galore include fastqc; data filter select "trim_galore" or "fastp"; iTools can choose not to do；
#### multiqc can convert common output to integration output on: Cutadapt/FastQC/Fastp/Bowtie 2/STAR/TopHat/deepTools/featureCounts/GATK/HOMER/HTSeq/MACS2/Picard/RSeQC/Samtools
fastqc *.fq 
iTools Fqtools stat -InFqList <.list> -OutStat <.txt>
trim_galore -j 4 -q 20 --phred33 --stringency 3 --length 20 --fastqc -e 0.1 --paired ${id}_1.fq.gz ${id}_2.fq.gz --gzip -o 02.cleandata/ #^N
fastp -i ${id}_1.fq.gz -I ${id}_2.fq.gz -o ${id}_1.clean.fq.gz -O ${id}_2.clean.fq.gz -w 12 -j ${id}.fastp.json -h ${id}.fastp.html;done

md5sum 01.Fastq_file/02.Cleandata/*.gz > 01.Fastq_file/02.Cleandata/Cleandata.md5.txt
multiqc -m fastp 01.Fastq_file/03.QcFastp/ -o 01.Fastq_file/03.QcFastp/ -n MultiQC_Fastp

#### Align to Genome, sort and index bam file
STAR --runThreadN 24 --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --genomeDir ~/Index/hg19/ --readFilesIn ${id}_1.clean.fq.gz ${id}_2.clean.fq.gz --outFileNamePrefix 02.Alignment/01.SortedBam/${id} #^N

# cat 01.Fastq_file/01.Rawdata/SampleList.txt | while read id;do bsub -q charge_normal -n 24 -J ${id} STAR --runThreadN 24 --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --genomeDir ~/Index/STAR-Index/hg19_2.7.0.f/ --readFilesIn 01.Fastq_file/02.Cleandata/${id}_1_val_1.fq.gz ${id}_2_val_2.fq.gz --outFileNamePrefix 02.Alignment/01.SortedBam/${id};done

samtools sort -@ 24 -o ${id}.sorted.bam ${id}.bam #^N
samtools index -@ 24 ${id}.sorted.bam #^N

#### Quality control with alignment file
/data/gpfs02/jliu/xpan/02.RNASeq/0.Script/RSeQC.sh ${id}.sorted.bam 02.Alignment/03.QcRSeQC ${id} ~/Index/hg19_refGene_2016-10-24.bed #^N

    # RSeQC.sh
    '''{shell}
    #!/bin/bash
    Bamfile=$1
    Output=$2
    Prefix=$3
    Refseq=$4

    geneBody_coverage.py -r ${Refseq} -i ${Bamfile} -o ${Output}/geneBody_coverage/${Prefix}
    infer_experiment.py -r ${Refseq} -i ${Bamfile} > ${Output}/infer_experiment/${Prefix}
    #inner_distance.py  -r ${Refseq} -i ${Bamfile} -o ${Output}/inner_distance/${Prefix}
    read_duplication.py -i ${Bamfile} -o ${Output}/read_duplication/${Prefix}
    #read_GC.py -i ${Bamfile} -o ${Output}/read_GC/${Prefix}
    '''

#### Distribution on genome component 
/data/gpfs02/jliu/xpan/02.RNASeq/0.Script/annovar.sh ${id}.sorted.bam 02.Alignment/02.ReadDistribution/ ${id} hg19 #^N
python3 ~/xpan/02.RNASeq/0.Script/readDISTmerge.py -i 02.Alignment/02.ReadDistribution/ -o 02.Alignment/02.ReadDistribution/ReadDistribution.xls #^N
    # annovar.sh
    '''{r}
    #!/bin/bash
    BamFile=$1
    Outdir=$2
    Prefix=$3
    species=$4

    if [ ${species} == "hg19" ]
    then
       db="humandb"
    elif [ ${species} == "mm10" ]
    then
       db="mousedb"
    elif [ ${species} == "bosTau8" ]
    then
       db="bosTau8db"
    elif [ ${species} == "rn6" ]
    then
       db="ratdb"
    else
       echo "Species information input is incorrect."
       exit
    fi
    echo $db
    samtools view -b -u -F 4 -F 256 $BamFile | bedtools bamtobed -bed12 -i - | perl -n -a -F"\t" -e '@F[3,4]=(0,0);$_ = join "\t",@F;print' 1> ${Outdir}/${Prefix}.bed
    /data/gpfs01/jliu/software/annovar/annotate_variation.pl --geneanno --outfile ${Outdir}/${Prefix} --buildver ${species} --thread 24 --maxgenethread 24 ${Outdir}/${Prefix}.bed /data/gpfs01/jliu/software/annovar/${db}
    awk '{print $1}' ${Outdir}/${Prefix}.variant_function | sort | uniq -c > ${Outdir}/${Prefix}.read_distribution.xls 
    '''

    # readDISTmerge.py
    '''{python}
    dir = para.dir
    Files = sorted(os.listdir(dir))


    for x in range(0, len(Files)):
        tmpData = pd.read_csv(dir + Files[x], header = None, sep = ' ', skipinitialspace = True)
        tmpData.columns = [ Files[x].split(sep = '.')[0], 'Region' ]
        print(tmpData.columns)
        tmpData.set_index('Region', inplace=True, drop=True)
        PreSum = tmpData.iloc[:,0].sum()

        tmpData2 = tmpData.ix[['exonic', 'intergenic', 'intronic', 'UTR5', 'UTR3'], :]
        LaterSum = tmpData2.iloc[:,0].sum()
        tmpData2.loc['Others'] = PreSum - LaterSum

        #print(tmpData2)

        if x == 0:
            MergeData = tmpData2
        else :
            MergeData = MergeData.join(tmpData2, how = 'outer', lsuffix='_left', rsuffix='_right')

    #print(MergeData)
    MergeData.to_csv(para.output, sep = '\t')
    '''

#### make the bigWig file for Genome Browser(UCSC etc.)
#### strand-specific RNA-seq should be separate to positive and negative strands
bamCoverage -p 20 --normalizeTo1x 2150570000 --filterRNAstrand forward -b ${id}.sorted.bam -o test.m.bw 
bamCoverage -p 20 --normalizeTo1x 2150570000 --filterRNAstrand reverse -b ${id}.sorted.bam -o test.p.bw

~/bin/bam2bw --library-type fr-firststrand ${id}.sorted.bam ~/Index/Genome/hg19.chromSizes #^N

    '''{perl}
    # bam2bw
    #!/usr/bin/perl
    use warnings;
    use Getopt::Long;
    use File::Basename;
    use File::Spec;
    use strict;

    sub printCMD{
        select STDERR;
        print "\tsam2bigwig [options] <alignment file> <chrom.sizes>\n";
        print "\n\t\toptions:\n";
        print "\t\t\t--library-type\t\tunstranded, fr-firststrand (default)\n";
        print "\t\t\t--name\t\t\toutput file name will be: name [pos|neg].ucsc.bigWig\n";
        exit;
    }

    my $libraryType = "fr-firststrand";
    my %availableLibraryType = (unstranded => 1, "fr-firststrand" => 1);
    my $name;
    GetOptions (
            "library-type=s" => \$libraryType,
            "name=s"    => \$name,
            ) or printCMD;
    unless (exists $availableLibraryType{$libraryType}){
        print "unknown library type!\n";
        printCMD;
    }

    if (@ARGV < 2){
        print "Too few parameters!\n";
        printCMD;
    }
    my $alignmentFile = File::Spec->rel2abs($ARGV[0]);
    my $chromSizesFile = File::Spec->rel2abs(pop @ARGV);

    print STDERR "making bigwig file aside alignment file\n";
    my($filename, $dirs, $suffix) = fileparse($alignmentFile,qr/\.[^.]+$/);
    chdir $dirs;
    #Separate the read mapped to + or - strand when strand specific
    my $totalAlignedReads;
    if ($libraryType ne "unstranded"){
        open IN,'-|',"samtools view -h $filename$suffix" or die "Can not open file $filename$suffix: $!\n";
        print STDERR "library type: fr-firststrand\n";
        my $layer_flag;
        open FWDSAM,'|-',"samtools view -u -o ${filename}_fwd.bam -";
        open REVSAM,'|-',"samtools view -u -o ${filename}_rev.bam -";
        while (<IN>){
            if (/^@/){
                print FWDSAM "$_";
                print REVSAM "$_";
                next;
            }
            my $flag = (split /\t/)[1];
            next if ($flag & 4) || ($flag & 256) ;
            $totalAlignedReads++;
            unless ($layer_flag){
                if ($flag & 1){
                    $layer_flag = 2;
                    print STDERR "library layer: paired end\n";
                }else{
                    $layer_flag = 1;
                    print STDERR "library layer: single end\n";
                }
            }
            if ($layer_flag == 1){
                if ($flag & 16){
                    print REVSAM "$_";
                }else{
                    print FWDSAM "$_";
                }
            }else{
                if (($flag & 128) && !($flag & 16)){
                    print FWDSAM "$_";
                }elsif(($flag & 64) && ($flag & 16)){
                    print FWDSAM "$_";
                }elsif(($flag & 128) && ($flag & 16)){
                    print REVSAM "$_";
                }elsif(($flag & 64) && !($flag & 16)){
                    print REVSAM "$_";
                }
            }
        }
        close FWDSAM;
        close REVSAM;
        close IN;
        my $factor1;
        if ($totalAlignedReads){
            $factor1 = 10000000/$totalAlignedReads;
        }else{
            print STDERR "Total aligned reads number not available\n";
            exit;
        }
        my $factor2 = 0 - $factor1;
        die if system "bedtools genomecov -bg -split -scale $factor1 -ibam ${filename}_fwd.bam > ${filename}_fwd.bedGraph";
        unlink "${filename}_fwd.bam";
        die if system "bedtools genomecov -bg -split -scale $factor2 -ibam ${filename}_rev.bam > ${filename}_rev.bedGraph";
        unlink "${filename}_rev.bam";
        die if system "bedSort ${filename}_fwd.bedGraph ${filename}_fwd.bedGraph";
        #unlink "${filename}_fwd.bedGraph";
        die if system "bedSort ${filename}_rev.bedGraph ${filename}_rev.bedGraph";
        #unlink "${filename}_rev.bedGraph";

        die if system "bedGraphToBigWig ${filename}_fwd.bedGraph $chromSizesFile ${filename}pos.ucsc.bigWig";
        unlink "${filename}_fwd.bedGraph";
        die if system "bedGraphToBigWig ${filename}_rev.bedGraph $chromSizesFile ${filename}neg.ucsc.bigWig";
        unlink "${filename}_rev.bedGraph";
        if ($name){
            rename "${filename}pos.ucsc.bigWig","${name}pos.ucsc.bigWig";
            rename "${filename}neg.ucsc.bigWig","${name}neg.ucsc.bigWig";
        }
    }elsif($libraryType eq "unstranded"){
        chomp($totalAlignedReads = `samtools view -c -F 256 -F 4 $filename$suffix`);
        my $factor = 10000000/$totalAlignedReads;
        die if system "bedtools genomecov -bg -split -scale $factor -ibam $filename$suffix > $filename.bedGraph";
        die if system "bedGraphToBigWig $filename.bedGraph $chromSizesFile $filename.bigwig";
        unlink "$filename.bedGraph";
        rename "$filename.bigwig", "$name.bigwig" if $name;
    }
    '''
#### Find different expression gene
cuffdiff -o 03.DiffExp/P0_P12/ -L P0,P12 -M ~/Index/hg19_rmsk.gtf -b ~/Index/hg19.fa -u -p 24 --library-type fr-firststrand ~/Index/hg19_refGene_.gtf ${id_1}.bam ${id_2}.bam

#### parsing cuffdiff output for merging expression of each sample and statistics
parsingExpdiff.py -d 03.DiffExp/P0_P12/ -o 03.DiffExp/P0_P12_parsing ~/Index/GeneType/RNA_annotate_hg19.txt
    # parsingExpdiff.py
    '''{python3}
    groupsInfo=pd.read_csv(args.FileDir+'read_groups.info', header = 0, sep = '\t', encoding = 'utf-8')
    Fpkm = pd.read_csv(args.FileDir + 'genes.read_group_tracking', header = 0, sep = '\t', encoding = 'utf-8')
    Exp = pd.read_csv(args.FileDir + 'gene_exp.diff', header = 0, sep = '\t', encoding = 'utf-8')
    annotate = pd.read_csv(args.reference, header = None, sep = '\t', encoding = 'utf-8')
    output = args.output_file
    # Fpkm = pd.read_table('/Users/panxu/Downloads/SMART-Seq_statistics/Chen.SZ_AM20180608-02/Ctrl_SPION10/genes.read_group_tracking.txt', header = 0, sep = '\t', encoding = 'utf-8')# Exp = pd.read_table('/Users/panxu/Downloads/SMART-Seq_statistics/Chen.SZ_AM20180608-02/Ctrl_SPION10/gene_exp.diff', header = 0, sep = '\t', encoding = 'utf-8')# annotate = pd.read_table('/Users/panxu/Downloads/SMART-Seq_statistics/Chen.SZ_AM20180608-02/Ctrl_SPION10/GeneType.mm10.txt', header = None, sep = '\t', encoding = 'utf-8')# output = '/Users/panxu/Downloads/SMART-Seq_statistics/Chen.SZ_AM20180608-02/Ctrl_SPION10/'

    def getFileName(x):
        str=x
        patern1=re.compile('\/')
        patern2=re.compile('\.')
        cutstr=patern1.split(str)
        fileName=patern2.split(repr(cutstr[-1]))
        fileName=fileName[0]
        return fileName[1:]

    #get grooup file

    groupsInfo['fileName']=groupsInfo['file'].map(lambda x:getFileName(x).replace("Aligned", ""))
    groupsInfo['fileName']=groupsInfo['fileName']
    groupsInfo['replicate_num'] = groupsInfo['replicate_num'].map(lambda x: str(x))
    groupsInfo['condition'] = groupsInfo['condition'].str.cat(groupsInfo['replicate_num'], sep = "_")
    groupsInfo.drop(['file','replicate_num','total_mass','norm_mass','internal_scale','external_scale'],axis=1,inplace=True)

    print(groupsInfo)
    # Fpkm reshap
    Fpkm['replicate'] = Fpkm['replicate'].map(lambda x: str(x))
    Fpkm['condition'] = Fpkm['condition'].str.cat(Fpkm['replicate'], sep = "_")

    Fpkm=pd.merge(left=Fpkm,right=groupsInfo,how='right',on='condition')


    Fpkm_pivot = Fpkm.pivot_table(index = ["tracking_id"], columns = ["fileName"], values = ["FPKM"])
    Fpkm_pivot = Fpkm_pivot['FPKM'].reset_index(drop = False)
    Fpkm_pivot.rename(columns = {
        'tracking_id': 'gene'
    }, inplace = True)


    annotate.columns = ['gene', 'type']
    Exp.drop(['test_id', 'gene_id', 'test_stat'], axis = 1, inplace = True)


    dfs = [annotate, Exp, Fpkm_pivot]

    df_final = reduce(lambda left, right: pd.merge(left, right, on = 'gene'), dfs)


    print("=============Merge2================")
    print(df_final.head(5));
    print("\n")

    # filter merge#
    #DGE = df_final[(((df_final["value_1"] > 0.5) & (df_final["value_2"] > 0.5)) | ((df_final["value_1"] == 0) & (df_final["value_2"] > 1)) | ((df_final["value_1"] > 1) & (df_final["value_2"] == 0)) ) & ((df_final["log2(fold_change)"] > 0.5849) | (df_final["log2(fold_change)"] < (-0.5849)))]
    DGE = df_final[(df_final["q_value"] < 0.05) & ((df_final["log2(fold_change)"] > 0.5849) | (df_final["log2(fold_change)"] < (-0.5849)))]
    Up = DGE[DGE["log2(fold_change)"] > 0.5849]
    Down = DGE[DGE["log2(fold_change)"] < (-0.5849)]

    DGE_stat = pd.DataFrame([
        ['All', len(df_final)],
        ['DGE', len(DGE)],
        ['DGE_Up', len(Up)],
        ['DGE_down', len(Down)]
    ], columns = [u'Type', u'Counts'])
    print("Differential Expression Gene Matrix:")
    #print(DGE.head(5));
    print("\n")
    print("Differential Expression Gene Statistics:")
    print(DGE_stat)


    # output
    sample1 = df_final['sample_1'].ix[1]
    sample2 = df_final['sample_2'].ix[1]

    df_final.to_csv(output + '01.Raw.' + sample1 + "-" + sample2 + "." + 'xls', header = True, index = False, sep = '\t', encoding = 'utf-8')
    DGE.to_csv(output + '02.DEG.' + sample1 + "-" + sample2 + "." + 'xls', header = True, index = False, sep = '\t', encoding = 'utf-8')
    Up.to_csv(output + '03.DEG.' + sample1 + "-" + sample2 + "." + sample2 + '-HighExp' + '.xls', header = True, index = False, sep = '\t', encoding = 'utf-8')
    Down.to_csv(output + '04.DEG.' + sample1 + "-" + sample2 + "." + sample2 + '-LowExp' + '.xls', header = True, index = False, sep = '\t', encoding = 'utf-8')
    DGE_stat.to_csv(output + '05.DEG.' + sample1 + "-" + sample2 + ".stat.xls", header = True, index = False, sep = '\t', encoding = 'utf-8')
    '''

#### Gene cluster analysis between GO&KEGG and GSEA
#### GSEA database: H.gmt/C1.gmt/C2.gmt/C3.gmt/C4.gmt/C5.gmt/C6.gmt/C7.gmt
#### GSEA input: cls and gmt file; Convert cuffdiff parsing output to cls and gmt file
java -cp gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea -gmx ~/Index/01.GSEA/Mouse/Mouse_H.gmt -nperm 1000 -scoring_scheme weighted -norm meandiv -res WT_MLL.gct -cls WT_MLL.cls#MLL_versus_WT -permute gene_set -rnd_type no_balance -collapse false -mode Max_probe -out GSEA/WT_MLL/h  -create_svgs true
R CMD

GoAnalysis(org=argv1, dir=argv2, out=argc3)
    # GO
    '''{r}
    GoAnalysis <- function(org="org.Hs.eg.db", dir="", out=""){
      FileList <- list.files(dir,full.names = F)
      JudgeLow <- stringr::str_detect(FileList,"04")
      JudgeHigh <- stringr::str_detect(FileList,"03")
      JudgeDEG <- stringr::str_detect(FileList,"02")
      
      Low <- read.table(paste0(dir,FileList[which(JudgeLow==TRUE)]),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
      High <- read.table(paste0(dir,FileList[which(JudgeHigh==TRUE)]),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
      All <- read.table(paste0(dir,FileList[which(JudgeDEG==TRUE)]),header = TRUE,sep = "\t",stringsAsFactors = FALSE)
      
      library(clusterProfiler)
      library(org,character.only = T)
      
      GoAnalysis <- function(gene,name,org){
        TranGene <- bitr(geneID = gene,fromType="SYMBOL", toType="ENTREZID", OrgDb=org)
        Level2CC <- groupGO(gene = TranGene$ENTREZID, OrgDb = org, ont="CC", level = 2, readable = T)@result
        Level2BP <- groupGO(gene = TranGene$ENTREZID, OrgDb = org, ont="BP", level = 2, readable = T)@result
        Level2MF <- groupGO(gene = TranGene$ENTREZID, OrgDb = org, ont="MF", level = 2, readable = T)@result
        
        Level2CC$Group = "CC"
        Level2BP$Group = "BP"
        Level2MF$Group = "MF"
        Level2Merge = do.call(rbind,list(Level2CC, Level2BP, Level2MF))
        
        EgoCc <- enrichGO(gene = TranGene$ENTREZID,OrgDb=org,minGSSize = 1,pvalueCutoff = 0.05,qvalueCutoff = 1,ont = "CC",readable = TRUE)
        EgoBp <- enrichGO(gene = TranGene$ENTREZID,OrgDb=org,minGSSize = 1,pvalueCutoff = 0.05,qvalueCutoff = 1,ont = "BP",readable = TRUE)
        EgoMf <- enrichGO(gene = TranGene$ENTREZID,OrgDb=org,minGSSize = 1,pvalueCutoff = 0.05,qvalueCutoff = 1,ont = "MF",readable = TRUE)
        
        if(nrow(EgoCc@result)==0){
          write.table(NULL,paste0(out,"/01.Gene_Ontology/02.GO_result/", name,".GeneOntology.CC.xls"),sep = "\t",row.names = F,quote = F)
        }else{
          EgoCc@result$Group="CC"
        }
        
        if(nrow(EgoBp@result)==0){
          write.table(NULL,paste0(out,"/01.Gene_Ontology/02.GO_result/", name,".GeneOntology.BP.xls"),sep = "\t",row.names = F,quote = F)
        }else{
          EgoBp@result$Group="BP"
        }
        
        if(nrow(EgoMf@result)==0){
          write.table(NULL,paste0(out,"/01.Gene_Ontology/02.GO_result/", name,".GeneOntology.MF.xls"),sep = "\t",row.names = F,quote = F)
        }else{
          EgoMf@result$Group="MF"
        }
        
        if(nrow(EgoCc@result)>=10){ResultCc=EgoCc@result[1:10,]}else{ResultCc=EgoCc@result}
        if(nrow(EgoBp@result)>=10){ResultBp=EgoBp@result[1:10,]}else{ResultBp=EgoBp@result}
        if(nrow(EgoMf@result)>=10){ResultMf=EgoMf@result[1:10,]}else{ResultMf=EgoMf@result}
        
        
        ResultAll <- Reduce(function(x,y) merge(x,y,all=T),list(ResultCc,ResultBp,ResultMf))
        write.table(Level2Merge,paste0(out, "/01.Gene_Ontology/01.GO_2rd/", name,".GO2rd.xls"),sep = "\t",row.names = F,quote = F)
        
        write.table(EgoCc@result,paste0(out, "/01.Gene_Ontology/02.GO_result/", name,".GeneOntology.Cc.xls"),sep = "\t",row.names = F,quote = F)
        write.table(EgoBp@result,paste0(out, "/01.Gene_Ontology/02.GO_result/", name,".GeneOntology.Bp.xls"),sep = "\t",row.names = F,quote = F)
        write.table(EgoMf@result,paste0(out, "/01.Gene_Ontology/02.GO_result/", name,".GeneOntology.Mf.xls"),sep = "\t",row.names = F,quote = F)
        write.table(ResultAll,paste0(out, "/01.Gene_Ontology/02.GO_result/", name,".Top10.xls"),sep = "\t",row.names = F,quote = F)
      }
      
      GoAnalysis(gene = Low$gene,  name = "LowExp",  org = org)
      GoAnalysis(gene = High$gene, name = "HighExp", org = org)
      GoAnalysis(gene = All$gene,  name = "AllExp",  org = org)
    }
    '''

#### Assemble denovo lncRNA
cufflinks -o 05.denovo/PSR_U266/cufflinks -p 12 -g ~/Index/hg19_refGene.gtf -M ~/Index/hg19_rmsk.gtf ${id}.sorted.bam
~/Script/filer_cufflinks_gtf.pl transcripts.gtf > transcripts.filter.gtf
cuffcompare -r ~/Index/hg19_refGene_2016-10-24.gtf -o ./cuffcompare cufflinks/transcripts.filter.gtf
perl -n -e'if (/class_code "u"/){print}' cuffcompare.combined.gtf > cuffcompare.combined.U.gtf


gffread test.gtf -g ~/Index/Genome/hg38.fa -w test.fa
awk '/^>/&&NR>1{print "";}{printf "%s",/^>/?$0"\n":$0}'  test.fa  >test.fa.2

~/Script/c.gtf2bed cuffcompare.combined.U.gtf | bedtools getfasta -fi ~/Index/hg19.fa -bed - -fo XIN_HCMV_lnc.fa -name -s -split
~/software/CPC/cpc-0.9-r2/bin/run_predict.sh cuffcompare/XIN_HCMV_lnc.fa CPC/result.table CPC/ CPC/result.evidence

pick_non-coding_from_CPC_table.pl CPC/result.table cuffcompare/F1-3_FwT1_lnc.fa F1-3_FwT1.combined.U.gtf Index/mm10_refGene.gtf 4.DiffExprAnalysis/4.2.DEGsList/F1-3_FwT1



#### Old Version 1:  Tophat + cufflinks[Zhang.WJ]
tophat -o OutputPE/ -p 20 --segment-mismatches 2 --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2  -r 70 --library-type fr-firststrand -G hg19_mRNA.gtf /bowtie2-2.2.1/hg19 r1.fastq.gz r2.fastq.gz

# tophat2 -o /data/gpfs02/jliu/hugs/data/RNA_SEQ/tophat_output/RQX/SR19061401LL_FKDL190745135-1a --library-type fr-firststrand -p 20 -G /data/gpfs02/jliu/hugs/database/hg19_refGene.gtf /data/gpfs02/jliu/hugs/database/bowtie2Index/hg19 /data/gpfs02/jliu/hugs/data/RNA_SEQ/fastp_output/RQX/SR19061401LL_FKDL190745135-1a_1.fped.fq.gz /data/gpfs02/jliu/hugs/data/RNA_SEQ/fastp_output/RQX/SR19061401LL_FKDL190745135-1a_2.fped.fq.gz

# cuffdiff -o /data/gpfs02/jliu/hugs/data/RNA_SEQ/cuffdiff_output/RQX_hs20 -L CTL,Rxn52,Rxn91,20-0218,OC -M /data/gpfs02/jliu/hugs/database/hg19_rmsk.gtf -b /data/gpfs02/jliu/hugs/database/hg19.fa -u -p 60 --library-type fr-firststrand /data/gpfs02/jliu/hugs/database/hg19_refGene.gtf RS_HS_MCF7_CTL_RQX46_hs20/accepted_hits.bam RS_HS_MCF7_Rxn52_RQX47_hs20/accepted_hits.bam RS_HS_MCF7_Rxn91_RQX48_hs20/accepted_hits.bam RS_HS_MCF7_20-0218_RQX49_hs20/accepted_hits.bam RS_HS_MCF7_OC_RQX50_hs20/accepted_hits.bam

cufflinks -p 40 -o RS_hs2_24_Tophat2_hg19_mRNA_gtf_cl_G -g hg19_mRNA.gtf --library-type fr-firststrand -u RS_hs2_24_Tophat2/accepted_hits.bam
cufflinks -p 40 -o RS_hs2_25_Tophat2_hg19_mRNA_gtf_cl_G -g hg19_mRNA.gtf --library-type fr-firststrand -u RS_hs2_25_Tophat2/accepted_hits.bam
cuffmerge -o RS_hs2_24_25_hg19_mRNA_gtf -g hg19_mRNA.gtf -p 40 assemblies_RS_hs2_24_25_hg19_mRNA_gtf.txt
cuffdiff -o output_cuffdiff/ RS_hs2_24_25_hg19_mRNA_gtf/merged.gtf -p 40 --library-norm-method geometric --dispersion-method blind  --library-type fr-firststrand -u RS_hs2_24_Tophat2/accepted_hits.bam RS_hs2_25_Tophat2/accepted_hits.bam

### Old Version 2: Find ifferent expression gene with gfold
### Someone is doing this step: Do gfold script with 1000 cycles
samtools view sample1.bam | gfold count -ann hg19_refGene.gtf -tag ${id}.sorted.bam -o sample1.read_cnt
gfold diff -s1 sample1 -s2 sample2 -suf .read_cnt -o sample1VSsample2.diff 

### Old Version 3: Find ifferent expression gene with featureCounts and R:edgeR
### The "featureCounts" can be replaced by "HTSeq"; The "R:edgeR" can be replaced by "R:DESeq2"
### When if there is no replicated with using R:edgeR , the BCV value should choose the suitable value
featureCounts -a hg19_refGene.gtf -o ${id}.featureCounts.txt -G hg18.fa -p -g 'gene_id' -s 2 --primary -T 24 ${id_1}.sorted.bam ${id_2}.sorted.bam #^N
python2 ~/software/HTseq/HTSeq-0.8.0/python2/scripts/htseq-count -f bam -r pos -s no -m union -t exon 02.Alignment/01.SortedBam/${id}.sorted.bam ~/Index/mm10_refGene.gtf > 02.Alignment/05.HTSeq/${id}.htseq.xls


    ibrary(edgeR)
        counts <- read.table( file = "featureCount_total.txt" , header = TRUE ,sep = "\t",row.names = 1)
        grp <- as.factor(substr(colnames(counts), 1, 2))
        o <- order(grp)
        pairs(log2(1+counts[,o[1:6]]), pch=".",lower.panel=NULL)
        d <- DGEList(counts=counts, group=grp)
        d <- calcNormFactors(d)
        dim(d)
        cps <- cpm(d)
        k <- rowSums(cps>=1) > 2
        d <- d[k,]
        dim(d)
        cols <- as.numeric(d$samples$group)

        plotMDS(d,col=cols)
        par(mfrow=c(2,2))
        plotMDS(d, col=cols, main="500 / lLFC",cex=0.5)
        plotMDS(d, col=cols, method="bcv", main="500 / BCV",cex=0.5)
        plotMDS(d, col=cols, top=12000, main="2000 / lLFC",cex=0.5)
        plotMDS(d, col=cols, top=12000, method="bcv", main="2000 / BCV",cex=0.5)
        mm <- model.matrix(~-1+grp)
        d <- estimateGLMCommonDisp(d,mm)
        d <- estimateGLMTrendedDisp(d,mm)
        d <- estimateGLMTagwiseDisp(d,mm)
        par(mfrow=c(1,1))
        plotBCV(d)
        d$common.dispersion
        sqrt(d$common.dispersion)
        plotMeanVar(d, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
        plotSmear(d, pair=c("N","C"), ylim=c(-5,5))
        f <- glmFit(d,mm)
        con <- makeContrasts("DE-ES"=grpC-grpN,levels=colnames(mm))
        lrt <- glmLRT(f,contrast=con)
        topTags(lrt,20)
        cps <- cpm(d)
        tt <- topTags(lrt, n=Inf)$table

        cps <- as.data.frame(cps)
        cps[,"gene"] <- rownames(cps)
        tt[,"gene"] <- rownames(tt)
        tt2 <- merge(cps,tt,by="gene",all=TRUE)
        write.table(tt,"~/edgeR_C_N.xls",row.names = FALSE,sep = "\t",quote = FALSE)


