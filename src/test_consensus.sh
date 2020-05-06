#!/bin/sh

# nanosvでVCFコールし、bcftoolsでコンセンサスを出力
# まずはPrdm14 barcode26.fa コントロールでお試し
# nanosv: https://github.com/mroosmalen/nanosv
# bcftools: https://bi.biopapyrus.jp/gwas/vcf-consensus-fasta.html

# ============================================================================
# I/O and Arguments
# ============================================================================

# ============================================================================
# I/O and Arguments
# ============================================================================

barcode=barcode23
alleletype=target
suffix="${barcode}_${alleletype}"

fasta=.DAJIN_temp/fasta_ont/"${barcode}".fa
clustering_id=$(find .DAJIN_temp/clustering/result_allele_id* | grep "$suffix")

mkdir -p .DAJIN_temp/seqlogo/temp/
tmp_fa=".DAJIN_temp/seqlogo/temp/${suffix}.fa"
tmp_targetseq=".DAJIN_temp/seqlogo/temp/targetseq_${suffix}"
tmp_targetseq_fa=".DAJIN_temp/seqlogo/temp/targetseq_${suffix}.fa"
tmp_targetseq_trimmed_fa=".DAJIN_temp/seqlogo/temp/targetseq_trimmed_${suffix}.fa"
tmp_mappedseq=".DAJIN_temp/seqlogo/temp/mapped_seq_${suffix}"
tmp_lalign=".DAJIN_temp/seqlogo/temp/lalign_${suffix}"
tmp_lalignseq=".DAJIN_temp/seqlogo/temp/lalignseq_${suffix}"

# ============================================================================
# FASTAリードにクラスタ番号をラベルづけする
# ============================================================================

cat "${fasta}" |
    awk '{if(NR % 2 != 0) $1="HOGE"$1
        else $1="FUGA"$1
        print}' |
    tr -d "\n" |
    sed -e "s/HOGE/\n/g" -e "s/FUGA/ /g" |
    awk '{print $1,$NF}' |
    sed "s/^[@|>]//g" |
    grep -v "^$" |
    sort -t " " |
    join - "${clustering_id}" |
    awk '{print ">"$NF"@@@"$1"\n"$2}' |
head -n 20000 > ${tmp_fa}

# # --------------------------------
# # コントロールも10000本とってくる
# # --------------------------------
# cat .DAJIN_temp/fasta_ont/barcode26.fa |
# head -n 20000 > tmp_cont.fa

# # ============================================================================
# # ngmlrとsnifflesによるVCF出力
# # ============================================================================
# # conda install -y -n DAJIN -c bioconda ngmlr sniffles
# # --------------------------------
# # コントロールのコンセンサス配列づくり
# # --------------------------------

# ref=".DAJIN_temp/fasta_conv/wt.fa"
# # que=".DAJIN_temp/fasta_ont/barcode18.fa"
# # que="$tmp_fa"
# que="tmp_cont.fa"
# # que=".DAJIN_temp/fasta_ont/barcode26.fa"
# threads=12
# ngmlr -t "$threads" -r "$ref" -q "$que" -x ont |
# samtools sort -@ "$threads" -O BAM -o mapped.sort.bam -
# samtools index -@ "$threads" mapped.sort.bam

# sniffles -t "$threads" -m mapped.sort.bam -v output.vcf

# wc -l output.vcf

# cat output.vcf | 
# bcftools sort - |
# bgzip -f -c > output.vcf.gz
# tabix -f -p vcf output.vcf.gz

# # VCFtools
# # cat "$ref" | vcf-consensus output.vcf.gz > out.fa

# # Bcftools
# cat "$ref" |
#     bcftools consensus output.vcf.gz |
#     awk '{if($1 ~/^>/) $0=$0"HOGE"
#         printf $0}' |
#     sed -e "s/$/\n/g" -e "s/HOGE/\n/g" |
# cat - > cont_consensus.fa

# diff $ref cont_consensus.fa
# --------------------------------
# targetのコンセンサス
# --------------------------------
ref=".DAJIN_temp/fasta_conv/target.fa"
que="$tmp_fa"
# que=".DAJIN_temp/fasta_ont/barcode18.fa"
threads=12
ngmlr -t "$threads" -r "$ref" -q "$que" -x ont |
samtools sort -@ "$threads" -O BAM -o mapped.sort.bam -
samtools index -@ "$threads" mapped.sort.bam

# bcftools
bcftools mpileup -f "$ref" mapped.sort.bam | bcftools call -mv -Ob -o calls.bcf
bcftools view -i '%QUAL>=20' calls.bcf

bcftools mpileup -f "$ref" mapped.sort.bam | bcftools call -mv -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
cat "$ref" |
    bcftools consensus calls.vcf.gz |
    awk '{if($1 ~/^>/) $0=$0"HOGE"
        printf $0}' |
    sed -e "s/$/\n/g" -e "s/HOGE/\n/g" |
cat - > consensus.fa
diff $ref consensus.fa


# sniffles
sniffles -t "$threads" -m mapped.sort.bam -v output.vcf

wc -l output.vcf

cat output.vcf | 
bcftools sort - |
bgzip -f -c > output.vcf.gz
tabix -f -p vcf output.vcf.gz

# Bcftools
cat "$ref" |
    bcftools consensus output.vcf.gz |
    awk '{if($1 ~/^>/) $0=$0"HOGE"
        printf $0}' |
    sed -e "s/$/\n/g" -e "s/HOGE/\n/g" |
cat - > cont_consensus.fa

diff $ref cont_consensus.fa

# ============================================================================
# nanosv
# ============================================================================
# ----------------------------------------------
# Encording languageを変える
# https://github.com/mroosmalen/nanosv/issues/45
# ----------------------------------------------

export LC_CTYPE=en_US.UTF-8
export LANG=en_US.UTF-8

# ----------------------------------------------
# ランダムなBedファイルを作製する
# Dependency; bash, mysql
# output; mm10_genome_sample.bed
# ----------------------------------------------

wget -q https://raw.githubusercontent.com/mroosmalen/nanosv/master/nanosv/scripts/create_random_position_bed.sh
chmod +x create_random_position_bed.sh
wget -q http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
./create_random_position_bed.sh mm10 mm10.chrom.sizes

# ----------------------------------------------
# VCFファイルの生成
# ----------------------------------------------
threads=12
bam=DAJIN_Report/bam/barcode26.bam
vcf="test.vcf"
NanoSV -t "$threads" "$bam" -b mm10_genome_sample.bed -o "$vcf"


# ============================================================================
# bcftools
# ============================================================================
ref=".DAJIN_temp/data/ref.fa"
fasta=".DAJIN_temp/fasta_ont/barcode26.fa"
minimap2 -t 12 -ax map-ont  "$ref" "$fasta" |samtools sort -@ 12 -O BAM > minimap2.bam

bcftools mpileup -Ou -f "$ref" minimap2.bam | bcftools call -mv -Oz -o minimap2_calls.vcf.gz
bcftools norm -f "$ref" minimap2_calls.vcf.gz -Ob -o minimap2_calls.norm.bcf
bcftools filter --IndelGap 5 minimap2_calls.norm.bcf -Oz -o minimap2_calls.norm.flt-indels.vcf.gz
cat "$ref" | bcftools consensus minimap2_calls.norm.flt-indels.vcf.gz > minimap2_consensus.fa

# ref=".DAJIN_temp/fasta_conv/wt.fa"
# ref=".DAJIN_temp/fasta_ont/barcode26.fa"
ref=".DAJIN_temp/data/ref.fa"
cat "$ref" | sed "s/mm10://g" > ref.fa
ref="ref.fa"
samtools faidx ref.fa

bcftools sort "$vcf" -o sorted_"$vcf"
bcftools view sorted_"$vcf" -Oz -o sorted_"$vcf".gz
bcftools index sorted_"$vcf".gz

bcftools norm -m-any "$vcf"  |
bcftools norm -Ov --check-ref w -f "$ref" > OUT.VCF

    bcftools consensus sorted_"$vcf".gz |
cat - > consensus.fa

# ============================================================================
# vcftools
# ============================================================================
bcftools sort "$vcf" -o sorted_"$vcf"
bgzip -f -c sorted_"$vcf" > sorted_"$vcf".gz
tabix -f -p vcf sorted_"$vcf".gz

cat ref.fa | vcf-consensus sorted_"$vcf".gz > out.fa



cat "$fasta" | awk 'length == 5087' | head | tail -n 1 > tmp_seq
samtools view "$bam" | grep -f tmp_seq - | cut -f 1 > tmp_ref
cat "$fasta" | grep -f tmp_ref -A1 > tmp_ref.fa

ref="tmp_ref.fa"
bam="minimap2.bam"
minimap2 -t 12 -ax map-ont  "$ref" "$fasta" |
samtools sort -@ 12 -O BAM \
> "$bam"
# VCF
vcf="minimap2.vcf"
NanoSV -t "$threads" "$bam" -b mm10_genome_sample.bed -o "$vcf"


bcftools mpileup -Ou -f "$ref" minimap2.bam |
bcftools call -mv -Oz -o minimap2_calls.vcf.gz
bcftools norm -f "$ref" minimap2_calls.vcf.gz -Ob -o minimap2_calls.norm.bcf
bcftools filter --IndelGap 5 minimap2_calls.norm.bcf -Oz -o minimap2_calls.norm.flt-indels.vcf.gz
cat "$ref" |
bcftools consensus minimap2_calls.norm.flt-indels.vcf.gz \
> minimap2_consensus.fa
