#!/bin/bash
########################################
##! @Author: Chenkai Lv
##! @Todo: ChIPseq pipeline
##! @Updated: 2018.09.28
##! @Dep: Shell
##! @ChangeLog:
##!   argument collection.
##!   ascp supported.
########################################


########################################
# Define work path and set up folders
WORK_PATH=$PWD
LOG_PATH=$WORK_PATH/log
REF_PATH=/data1/REFERENCE/index/bowtie2_hg19

# functions
function version {
  echo -e "\nauthor: Chenkai Lv
updated date: 2018.09.05\n"
  exit 1
}

function usage {
  echo -e "\n$(basename $0) is used to pretreat sra or conduct ChIP-seq analysis task.\n
  Usage:  sh $0 [options]
  Options:
    -h --help      display the usage and exit.
    -a --access    necessary file containing SRR number or sample name to be operate. [e.g. accession.txt]
    -d --dump      transform sra to fastq by fastq-dump.
    -f --fastp     perform QC to sample by fastp with the default parameters,the sequencing type P/S must be provided.(Pair end and Single end) [e.g. P/S]
    -c --cut       specify two comma separated INT parameters to cutadapte for sample. [e.g. 15,0]
    -b --bowtie2   use bowtie2 to map the clean data to the REFERENCE and transform to the bam file
    -s --sort      use samtools to sort the bam file
    -p --picard    use picard to remove duplicates
    -r --remove    remove unwanted reads using samtools
    -D --deeptools use deeptools to get BigWiggle from sorted bam for drawing gene profile
    -C --compare   necessary file containing group name and two SRR number in each line. [e.g. comparison]
    -P --peaks     use macs2 to call peaks from filtered bam
    -H --HOMER     use HOMER to find Motifs and TFBS
    -g --graph     transform to bigWig Format for genome browser
    -e --ceas      use ceas to draw peak distribution
    -v --version   display version information and exit.\n"

  echo -e "\e[36me.g. sh ChIPseq_pipeline.sh -a accession.txt -d -f S -c 5,3 -b -s -p -r -D -C comparison.txt -P -H -g -e\n\e[0m"
}


function exist_file {
  if [[ -e ${1} ]]; then
    return 1
  else
    return 2
  fi
}


TEMP=$(getopt -a -o ha:df:c:bsprDC:PHgev -l help,access:,dump,fastp:,cut:,bowtie2,sort,picard,remove,deeptools,compare:,peaks,HOMER,graph,ceas,version -n "$0" -- "$@")

if [ $# = 0 ]; then
  usage
  exit 1
fi

if [ $? != 0 ]; then
  echo -e "\e[31mERROR: Terminating... \e[0m\n" >&2
  exit 1
fi

eval set -- "$TEMP"

while [ -n "$1" ]
do
  case $1 in
    -h|--help)
      usage
      exit 1;;
    -v|--version)
      version
      exit 1;;
    -a|--access)
      if [[ -f $2 ]]; then
          ACCESS=$2
      else
          echo -e "\e[31mERROR: The file containing SRR number has not been provied! \e[0m\n" >&2
          exit 1
      fi
      shift 2;;
    -d|--dump)
      CTRL_FASTQ_DUMP=1
      shift;;
    -f|--fastp)
      if [[ $2 = "P" ]]; then
        type="Pair"
        CTRL_FASTP_DEFAULT=1
      elif [[ $2 = "S" ]]; then
        type="Single"
        CTRL_FASTP_DEFAULT=1
      else
        echo -e "\e[31mERROR: The parameter of sequence type must be provied correctly! \e[0m\n" >&2
        exit 1
      fi
      shift 2;;
    -c|--cut)
      if [[ -n $2 ]]; then
        FT=$2
        CTRL_FASTP_CLEAN=1
      else
        echo -e "\e[31mERROR: trimming bases in front and tail has not been provied! \e[0m\n" >&2
        exit 1
      fi
      shift 2;;
    -b|--bowtie2)
      CTRL_BOWTIE_MAP=1
      shift;;
    -s|--sort)
      CTRL_SORT_BAM=1
      shift;;
    -p|--picard)
      CTRL_REMOVE_DUP=1
      shift;;
    -r|--remove)
      CTRL_SAMTOOLS_FILTER=1
      shift;;
    -D|--deeptools)
      CTRL_BAMCOVERAGE=1
      shift;;
    -C|--compare)
      if [[ -f $2 ]]; then
        COMPARE=$2
      else
        echo -e "\e[31mERROR: The comparison scheme has not been provied! \e[0m\n" >&2
        exit 1
      fi
      shift 2;;
    -P|--peaks)
      CTRL_PEAK_CALLING=1
      shift;;
    -H|--HOMER)
      CTRL_TFBS_FINDING=1
      shift;;
    -g|--graph)
      CTRL_BIGWIG=1
      shift;;
    -e|--ceas)
      CTRL_CEAS=1
      shift;;
    --)
      shift
      break;;
    *)
      echo "Internal error!" >&2
      exit 1;;
  esac
done

########################################

# obtain hg19.chrom.sizes
if [[ ! -f hg19.chrom.sizes ]]; then
  fetchChromSizes hg19 > hg19.chrom.sizes
fi

if [[ ! -f hg19.refGene ]]; then
  wget http://liulab.dfci.harvard.edu/CEAS/src/hg19.refGene.gz
  gzip -d hg19.refGene.gz
fi

if [[ -e ${ACCESS} ]]; then
    echo -e "\e[37mAccession file found, continue...\e[0m"
elif [[ -e ${COMPARE} ]]; then
    echo -e "\e[37mComparison file found, continue...\e[0m"
else
    echo -e "\e[31mERROR: Accession file or comparison file must be provided.\e[0m" >&2
    exit 1
fi


# 1. fastq-dump extract fastq files from sra files
if [[ ${CTRL_FASTQ_DUMP} -eq 1 ]]; then
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
    mkdir -p $WORK_PATH/01_fastq-dump $LOG_PATH/01_fastq-dump
    exist_file $WORK_PATH/${sname}.sra
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[31mERROR: Sample ${sname} sra file are missing! \e[0m" >&2
      exit 1
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[32mSample ${sname} sra file found, continue... \e[0m"
    fi
    mkdir -p $LOG_PATH/01_fastq-dump
    echo $(date)
    exist_file $WORK_PATH/01_fastq-dump/${sname}*.fastq.gz
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start transforming SRA to fastq for ${sname}!=====\e[0m"
      fastq-dump --split-3 $WORK_PATH/${sname}.sra --gzip -O $WORK_PATH/01_fastq-dump 2> $LOG_PATH/01_fastq-dump/${sname}.fastq-dump.log
      echo -e "\e[33m=====finish transforming SRA to fastq for ${sname}!=====\n\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m${sname} fq is already existed\n\e[0m"
    fi
  done
fi


########################################
# remediation for fastp
# 2. Quality control with fastp
if [[ ${CTRL_FASTP_DEFAULT} -eq 1 ]]; then
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
    mkdir -p $WORK_PATH/02_fastp $LOG_PATH/02_fastp
    mkdir -p $WORK_PATH/02_fastqc $LOG_PATH/02_fastqc
    exist_file $WORK_PATH/01_fastq-dump/${sname}*.fastq.gz
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[31mERROR: Sample ${sname} fastq file are missing! \e[0m" >&2
      exit 1
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[32mSample ${sname} sra file found, continue... \e[0m"
    fi
    mkdir -p 02_fastp log/02_fastp
    echo $(date)
    exist_file $WORK_PATH/02_fastp/${sname}*.default.fastq.gz
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start first fastp for ${sname}!=====\e[0m"
      if [[ ${type} = "Single" ]]; then
        fastp -i $WORK_PATH/01_fastq-dump/${sname}.fastq.gz -o $WORK_PATH/02_fastp/${sname}.default.fastq.gz \
        --json  $WORK_PATH/02_fastp/${sname}.default.fp.json --html $WORK_PATH/02_fastp/${sname}.default.fp.html \
        2> $LOG_PATH/02_fastp/${sname}.default.fp.log
        echo -e "\e[33m=====finish first fastp for ${sname}!=====\n\e[0m"
      elif [[ ${type} = "Pair" ]]; then
        fastp -i $WORK_PATH/01_fastq-dump/${sname}_1.fastq.gz -o $WORK_PATH/02_fastp/${sname}_1.default.fastq.gz \
        -I $WORK_PATH/01_fastq-dump/${sname}_2.fastq.gz -O $WORK_PATH/02_fastp/${sname}_2.default.fastq.gz \
        --json  $WORK_PATH/02_fastp/${sname}.default.fp.json --html $WORK_PATH/02_fastp/${sname}.default.fp.html \
        2> $LOG_PATH/02_fastp/${sname}.default.fp.log
        echo -e "\e[33m=====finish first fastp for ${sname}!=====\e[0m"
      fi
    elif [[ ${value} -eq 1 ]]; then
        echo -e "\e[34mSample${sname} default.fq is already existed\n\e[0m"
    fi
  done
fi


if [[ ${CTRL_FASTP_CLEAN} -eq 1 ]]; then
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
    exist_file $WORK_PATH/02_fastp/${sname}*.default.fastq.gz
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[31mERROR: Sample ${sname} QC result are missing! \e[0m" >&2
      exit 1
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[32mSample ${sname} QC result found, continue... \e[0m"
    fi
    PARM_fastp=($(echo ${FT} | tr "," " "))
    if [[ -n ${PARM_fastp[0]} && -n ${PARM_fastp[1]} ]]; then
      echo -e "\e[35mParameter for cutting ${sname} found, continue... \e[0m"
    else
      echo -e "\e[31mERROR: Parameters for cutting ${sname} QC result are missing! \e[0m" >&2
      exit 1
    fi
    echo $(date)
    exist_file $WORK_PATH/02_fastp/${sname}*.clean.fastq.gz
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start second fastp for ${sname}!=====\e[0m"
      if [[ ${type} = "Single" ]]; then
        fastp -i $WORK_PATH/02_fastp/${sname}.default.fastq.gz -o $WORK_PATH/02_fastp/${sname}.clean.fastq.gz \
        --json $WORK_PATH/02_fastp/${sname}.clean.fp.json --html $WORK_PATH/02_fastp/${sname}.clean.fp.html \
        -f ${PARM_fastp[0]} -t ${PARM_fastp[1]} 2> $LOG_PATH/02_fastp/${sname}.clean.fp.log
        echo -e "\e[33m=====finish second fastp for ${sname}!=====\n\e[0m"
        echo -e "\e[33m=====start fastqc for ${sname}!=====\e[0m"
        fastqc ${WORK_PATH}/02_fastp/${sname}.clean.fastq.gz -o ${WORK_PATH}/02_fastqc \
               2> $LOG_PATH/02_fastqc/${sname}.fastqc.log
        echo -e "\e[33m=====finish fastqc for ${sname}!=====\n\e[0m"
      elif [[ ${type} = "Pair" ]]; then
        fastp -i $WORK_PATH/02_fastp/${sname}_1.default.fastq.gz -o $WORK_PATH/02_fastp/${sname}_1.clean.fastq.gz \
        -I $WORK_PATH/02_fastp/${sname}_2.default.fastq.gz -O $WORK_PATH/02_fastp/${sname}_2.clean.fastq.gz \
        --json $WORK_PATH/02_fastp/${sname}.clean.fp.json --html $WORK_PATH/02_fastp/${sname}.clean.fp.html \
        -f ${PARM_fastp[0]} -t ${PARM_fastp[1]} 2> $LOG_PATH/02_fastp/${sname}.clean.fp.log
        echo -e "\e[33m=====finish second fastp for ${sname}!=====\n\e[0m"
        echo -e "\e[33m=====start fastqc for ${sname}!=====\e[0m"
        fastqc ${WORK_PATH}/02_fastp/${sname}_1.clean.fastq.gz -o ${WORK_PATH}/02_fastqc \
               2> $LOG_PATH/02_fastqc/${sname}.fastqc.log
               fastqc ${WORK_PATH}/02_fastp/${sname}_2.clean.fastq.gz -o ${WORK_PATH}/02_fastqc \
               2> $LOG_PATH/02_fastqc/${sname}.fastqc.log
        echo -e "\e[33m=====finish fastqc for ${sname}!=====\n\e[0m"
      fi
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Sample ${sname} clean fastp is already existed!=====\n\e[0m"
    fi
  done
fi


# 3. bowtie2 alignment and sam2bam, sample stat by flagstat
if [[ ${CTRL_BOWTIE_MAP} -eq 1 ]]; then
  mkdir -p $WORK_PATH/03_bowtie2 $LOG_PATH/03_bowtie2
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
    echo ${date}
    exist_file $WORK_PATH/03_bowtie2/${sname}.bam
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start bowtie mapping for ${sname}!=====\e[0m"
      if [[ ${type} = "Single" ]]; then
        bowtie2 -p 20 -X 2000 -x $REF_PATH/hg19_dna -U $WORK_PATH/02_fastp/${sname}.clean.fastq.gz 2>$LOG_PATH/03_bowtie2/${sname}.bowtie2.log \
        | samtools view -@ 10 -bT $REF_PATH/index/bowtie2_hg19/hg19_dna.fa \
        -o $WORK_PATH/03_bowtie2/${sname}.bam >$LOG_PATH/03_bowtie2/${sname}.bamconv.log 2>&1
        echo -e "\e[33m=====finish bowtie mapping for ${sname}!=====\n\e[0m"
      elif [[ ${type} = "Pair" ]]; then
        bowtie2 -p 20 -X 2000 -x $REF_PATH/hg19_dna -1 $WORK_PATH/02_fastp/${sname}_1.clean.fastq.gz -2 $WORK_PATH/02_fastp/${sname}_2.clean.fastq.gz \
        2>$LOG_PATH/03_bowtie2/${sname}.bowtie2.log | samtools view -@ 10 -bT $REF_PATH/index/bowtie2_hg19/hg19_dna.fa \
        -o $WORK_PATH/03_bowtie2/${sname}.bam >$LOG_PATH/03_bowtie2/${sname}.bamconv.log 2>&1
        echo -e "\e[33m=====finish bowtie mapping for ${sname}!=====\n\e[0m"
      fi
      echo -e "\e[33m=====start simple stats by samtools flagstat for ${sname}!=====\e[0m"
      samtools flagstat ${WORK_PATH}/03_bowtie2/${sname}.bam > ${WORK_PATH}/03_bowtie2/${sname}.bam.stat \
      2> $LOG_PATH/03_bowtie2/${sname}.samtools.flagstat.log
      echo -e "\e[33m=====finish simple stats by samtools flagstat for ${sname}!=====\n\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Sample ${sname} bam is already existed!=====\n\e[0m"
    fi
  done
fi


# 4. sort bam files by using samtools
if [[ ${CTRL_SORT_BAM} -eq 1 ]]; then
  mkdir -p $WORK_PATH/04_samtools_sort $LOG_PATH/04_samtools_sort
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
    echo ${date}
    exist_file $WORK_PATH/04_samtools_sort/${sname}_sorted.bam
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start sort bam for ${sname}!=====\e[0m"
      samtools sort -m 1024M -@ 20 -o $WORK_PATH/04_samtools_sort/${sname}_sorted.bam $WORK_PATH/03_bowtie2/${sname}.bam \
      >$LOG_PATH/04_samtools_sort/${sname}_sort.log 2>&1
      echo -e "\e[33m=====finish sort bam for ${sname}!=====\n\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Sample ${sname} sorted_bam is already existed!=====\n\e[0m"
    fi
  done
fi


# 5.1. deduplication of bam files by using picard
if [[ ${CTRL_REMOVE_DUP} -eq 1 ]]; then
  mkdir -p $WORK_PATH/05_picard $LOG_PATH/05_picard
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
    echo ${date}
    exist_file $WORK_PATH/05_picard/dedup_${sname}.bam
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start remove deduplication for ${sname}!=====\e[0m"
      java -jar /share/apps/RNAseq_tool/picard/picard.jar MarkDuplicates \
            I=$WORK_PATH/04_samtools_sort/${sname}_sorted.bam \
            O=$WORK_PATH/05_picard/dedup_${sname}.bam \
            M=$LOG_PATH/05_picard/${sname}.dedup.log \
            REMOVE_DUPLICATES=true
      echo -e "\e[33m=====start remove deduplication for ${sname}!=====\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====file dedup_${sname}.bam is already existed!=====\n\e[0m"
    fi
  done
fi


# 5.2. filter out unwanted reads by using samtools
if [[ ${CTRL_SAMTOOLS_FILTER} -eq 1 ]]; then
  mkdir -p $WORK_PATH/05_samtools_filter $LOG_PATH/05_samtools_filter
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
    echo ${date}
    exist_file $WORK_PATH/05_picard/dedup_${sname}.bam.bai
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start build index for dedup_${sname}.bam!=====\e[0m"
      samtools index $WORK_PATH/05_picard/dedup_${sname}.bam
      echo -e "\e[33m=====finish build index for dedup_${sname}.bam!=====\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Index has already build!=====\e[0m"
    fi
    exist_file $WORK_PATH/05_samtools_filter/filtered_${sname}.bam
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start filter out unwanted reads for ${sname}!=====\e[0m"
      samtools view -b -@ 20 -q 30 -F 4 \
                    -o $WORK_PATH/05_samtools_filter/filtered_${sname}.bam $WORK_PATH/05_picard/dedup_${sname}.bam \
                       chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 \
                       chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX \
                     >$LOG_PATH/05_samtools_filter/${sname}.filter.log 2>&1
      echo -e "\e[33m=====finish filter out unwanted reads for ${sname}!=====\e[0m"
      echo -e "\e[33m=====start build index for filtered_${sname}.bam!=====\e[0m"
      samtools index $WORK_PATH/05_samtools_filter/filtered_${sname}.bam
      echo -e "\e[33m=====finish build index for filtered_${sname}.bam!=====\n\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====file filtered_${sname}.bam is already existed!=====\n\e[0m"
    fi
  done
fi


# 6. get BigWiggle from sorted bam for drawing gene profile
if [[ ${CTRL_BAMCOVERAGE} -eq 1 ]]; then
  mkdir -p $WORK_PATH/06_deeptools $LOG_PATH/06_deeptools
  for sname in $(cat ${ACCESS} | grep -v "^$");
  do
  echo $(date)
  exist_file $WORK_PATH/06_deeptools/${sname}.bw
  value=$?
  if [[ ${value} -eq 2 ]]; then
    echo -e "\e[33m=====start obtain BigWiggle file by bamCoverage for ${sname}!=====\e[0m"
    bamCoverage -b ${WORK_PATH}/05_samtools_filter/filtered_${sname}.bam -p 2 --outFileFormat bigwig \
                -o ${WORK_PATH}/06_deeptools/${sname}.bw 2> ${WORK_PATH}/log/06_deeptools/${sname}.bw.log
    echo -e "\e[33m=====finish obtain BigWiggle file by bamCoverage for ${sname}!=====\n\e[0m"
  elif [[ ${value} -eq 1 ]]; then
    echo -e "\e[34m=====file ${sname}.bw is already existed!=====\n\e[0m"
  fi
  done
  echo "==================================================================="
  exist_file $WORK_PATH/06_deeptools/*.pdf
  value=$?
  if [[ ${value} -eq 2 ]]; then
    echo -e "\e[33m=====start plot profile and heatmap about reference-point!=====\e[0m"
    computeMatrix reference-point -p 4 --referencePoint TSS -b 2000 -a 2000 -S $WORK_PATH/06_deeptools/*.bw \
                  -R /data1/REFERENCE/bed/Human__refGene_deeptools.bed --skipZeros \
                  -o $WORK_PATH/06_deeptools/tss.mat.gz 2> $LOG_PATH/06_deeptools/tss.matrix.log
    plotProfile --dpi 720 -m $WORK_PATH/06_deeptools/tss.mat.gz -o $WORK_PATH/06_deeptools/tss.profile.pdf \
                --plotFileFormat pdf --perGroup 2> $LOG_PATH/06_deeptools/tss.profile.log
    plotHeatmap --dpi 720 -m $WORK_PATH/06_deeptools/tss.mat.gz -o $WORK_PATH/06_deeptools/tss.heatmap.pdf \
                --plotFileFormat pdf 2> $LOG_PATH/06_deeptools/tss.heatmap.log
    echo -e "\e[33m=====finish plot profile and heatmap about reference-point!=====\e[0m"
    echo -e "\e[33m=====start plot profile and heatmap about scale-regions!=====\e[0m"
    computeMatrix scale-regions -p 4 -b 2000 -a 2000 -m 5000 -S $WORK_PATH/06_deeptools/*.bw \
                  -R /data1/REFERENCE/bed/Human__refGene_deeptools.bed --skipZeros \
                  -o $WORK_PATH/06_deeptools/gene.mat.gz 2> $LOG_PATH/06_deeptools/gene.matrix.log
    plotProfile --dpi 720 -m $WORK_PATH/06_deeptools/gene.mat.gz -o $LOG_PATH/06_deeptools/gene.profile.pdf \
                --plotFileFormat pdf --perGroup 2> $LOG_PATH/06_deeptools/gene.profile.log
    plotHeatmap --dpi 720 -m $WORK_PATH/06_deeptools/gene.mat.gz -o $WORK_PATH/06_deeptools/gene.heatmap.pdf \
                --plotFileFormat pdf 2> $LOG_PATH/06_deeptools/gene.heatmap.log
    echo -e "\e[33m=====finish plot profile and heatmap about scale-regions!=====\n\e[0m"
  elif [[ ${value} -eq 1 ]]; then
    echo -e "\e[34m=====Heatmap file is already existed!=====\n\e[0m"
  fi
fi


# 7. peak calling by using macs2
if [[ ${CTRL_PEAK_CALLING} -eq 1 ]]; then
  mkdir -p $WORK_PATH/07_macs2 $LOG_PATH/07_macs2
  while read -r line
  do
    group=(${line})
    echo ${date}
    exist_file $WORK_PATH/07_macs2/${group[1]}_summits.bed
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "=====start finding peaks by MACS2 for ${group[1]} vs ${group[2]}, comparison: ${group[0]}!====="
      macs2 callpeak -t $WORK_PATH/05_samtools_filter/filtered_${group[1]}.bam -c $WORK_PATH/05_samtools_filter/filtered_${group[2]}.bam \
                --format BAM -B --name ${group[0]} -g hs -q 0.01 --outdir $WORK_PATH/07_macs2 \
                2> $LOG_PATH/07_macs2/${group[0]}.macs2.log
      echo -e "=====finish finding peaks by MACS2 for ${group[1]} vs ${group[2]}, comparison: ${group[0]}!=====\n"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====MACS peak file is already existed!=====\n\e[0m"
    fi
  done < ${COMPARE}
fi


# 8. Motif and TFBS finding with HOMER
if [[ ${CTRL_TFBS_FINDING} -eq 1 ]]; then
  mkdir -p $WORK_PATH/08_homer $LOG_PATH/08_homer
  while read -r line
  do
    group=(${line})
    echo ${date}
    exist_file $LOG_PATH/08_homer/${group[0]}.homer.log
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start find Motifs for ${group[0]}!=====\e[0m"
      mkdir $WORK_PATH/08_homer/${group[0]}
      findMotifsGenome.pl $WORK_PATH/07_macs2/${group[0]}_summits.bed hg19 $WORK_PATH/08_homer/${group[0]} \
      -preparsedDir $WORK_PATH/08_homer >$LOG_PATH/08_homer/${group[0]}.homer.log 2>&1
      echo -e "\e[33m=====finish find Motifs for ${group[0]}!=====\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Motif finding file is already existed!=====\n\e[0m"
    fi
  done < ${COMPARE}
fi


# 9. bam->bdg->bigwig from sorted bam for external visualization
if [[ ${CTRL_BIGWIG} -eq 1 ]]; then
  mkdir -p $WORK_PATH/09_bigwig $LOG_PATH/09_bigwig
  for sname in $(cat ${ACCESS} | grep -v "^$");   
  do
    echo ${date}
    exist_file $WORK_PATH/09_bigwig/${sname}.bdg
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start pileup reads for ${sname}!=====\e[0m"
      macs2 pileup --extsize 10 -i $WORK_PATH/05_samtools_filter/filtered_${sname}.bam -o $WORK_PATH/09_bigwig/${sname}.bdg \
      2> $LOG_PATH/09_bigwig/${sname}.bgconv.log
      echo -e "\e[33m=====finish pileup reads for ${sname}!=====\n\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Sample ${sname} bdg file is already existed!=====\n\e[0m"
    fi

    echo ${date}
    exist_file $WORK_PATH/09_bigwig/${sname}.bed
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start bam2bed for ${sname}!=====\e[0m"
      bamToBed -i $WORK_PATH/05_samtools_filter/filtered_${sname}.bam -split > $WORK_PATH/09_bigwig/${sname}.bed
               2> $LOG_PATH/09_bigwig/${sname}.bam2bed.log
      echo -e "\e[33m=====finish bam2bed for ${sname}!=====\n\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Sample ${sname} bed file is already existed!=====\n\e[0m"
    fi

    echo ${date}
    exist_file $WORK_PATH/09_bigwig/${sname}.bigwig
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo -e "\e[33m=====start bed2bedgraph for ${sname}!=====\e[0m"
      cat $WORK_PATH/09_bigwig/${sname}.bed | sort -k 1,1 | bedItemOverlapCount hg19 -chromSize=hg19.chrom.sizes stdin | \
          sort -k 1,1 -k 2,2n $WORK_PATH/09_bigwig/${sname}.bdg > $WORK_PATH/09_bigwig/sorted_${sname}.bdg 2> $LOG_PATH/09_bigwig/${sname}.bed2bdg.log
      echo -e "\e[33m=====finish bed2bedgraph for ${sname}!=====\n\e[0m"
      echo -e "\e[33m=====start bedgraph2bigwiggle for ${sname}!=====\e[0m"
      bedGraphToBigWig $WORK_PATH/09_bigwig/sorted_${sname}.bdg $WORK_PATH/hg19.chrom.sizes $WORK_PATH/09_bigwig/${sname}.bigwig \
                       2> $LOG_PATH/09_bigwig/${sname}.bdg2bw.log
      echo -e "\e[33m=====finish bedgraph2bigwiggle for ${sname}!=====\n\e[0m"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====Sample ${sname} bigwig file is already existed!=====\n\e[0m"
    fi
  done
fi


# 10. peaks visualization by CEAS
if [[ ${CTRL_CEAS} -eq 1 ]]; then
  mkdir -p $WORK_PATH/10_ceas $LOG_PATH/10_ceas
  while read -r line
  do
    group=(${line})
    echo ${date}
    exist_file $WORK_PATH/10_ceas/${group[2]}.wig
    value=$?
    if [[ ${value} -eq 2 ]]; then
      echo $(date)
      echo -e "=====start bigwiggle2bedgraph by bigWigToWig for ${group[1]} and ${group[2]}!====="
      bigWigToWig $WORK_PATH/06_deeptools/${group[1]}.bw $WORK_PATH/10_ceas/${group[1]}.wig.bdg
      bigWigToWig $WORK_PATH/06_deeptools/${group[2]}.bw $WORK_PATH/10_ceas/${group[2]}.wig.bdg
      echo -e "=====finish bigwiggle2bedgraph by bigWigToWig for ${sname}!=====\n"
      echo $(date)
      echo -e "=====start bedgraph2wiggle by python for ${group[1]} and ${group[2]}!====="
      fakeWig2TureWig.py -i $WORK_PATH/10_ceas/${group[1]}.wig.bdg -n ${group[1]} -o $WORK_PATH/10_ceas/${group[1]}.wig
      fakeWig2TureWig.py -i $WORK_PATH/10_ceas/${group[2]}.wig.bdg -n ${group[2]} -o $WORK_PATH/10_ceas/${group[2]}.wig
      echo -e "=====finish bedgraph2wiggle by python for ${group[1]} and ${group[2]}!=====\n"
      echo $(date)
      echo -e "=====start peaks visualization by CEAS for ${group[1]} vs ${group[2]}, comparison: ${group[0]}!====="
      ceas --name=${group[0]}_CEAS --pf-res=20 -g $WORK_PATH/hg19.refGene -b $WORK_PATH/07_macs2/${group[0]}_summits.bed \
           -w $WORK_PATH/10_ceas/${group[2]}.wig 2> $LOG_PATH/10_ceas/${group[0]}.ceas.log
      echo -e "=====finish peaks visualization by CEAS for ${group[1]} vs ${group[2]}, comparison: ${group[0]}!=====\n"
    elif [[ ${value} -eq 1 ]]; then
      echo -e "\e[34m=====MACS peak file is already existed!=====\n\e[0m"
    fi
  done < ${COMPARE}
fi

