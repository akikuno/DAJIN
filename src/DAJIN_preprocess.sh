#!/bin/sh

ref=".DAJIN_temp/fasta_conv/wt.fa"
que=".DAJIN_temp/fasta_conv/target#c140G_C.fa"
que=".DAJIN_temp/fasta_ont/barcode12.fa"
que=".DAJIN_temp/fasta_ont/target#c140G_C_simulated_aligned_reads.fasta"
threads=13

# ============================================================
# 配列長を入手する
# ============================================================
wt_seqlen=$(awk '!/[>|@]/ {print length($0)}' "$ref")

output_mids="DAJIN_new_mids"
true > "$output_mids"

for que in .DAJIN_temp/fasta_ont/*; do
    # ============================================================
    # label
    # ============================================================
    label=$(echo "$que" |
        sed -e "s|.*/||g" \
            -e "s|_aligned.*||g" \
            -e "s|.fa*||g")
    echo "$label"
    mapping_ref=$(cat "$ref" | grep ">" | sed "s/>//g")
    # ============================================================
    # mappinng
    # ============================================================

    minimap2 -t "$threads" -ax map-ont "$ref" "$que" --cs=long 2>/dev/null |
    grep -v "^@" |
    awk -v ref="$mapping_ref" '$3 == ref' |
    awk '$2==0 || $2==16' |
    # cat - > test.sam
    # ============================================================
    # MIDS conv
    # ============================================================
    # cat test.sam | grep target#c140G-C_1408_aligned_35356_F_1118_178_1500
    awk '{print $1, $2, $4, $(NF-1)}' |
    awk -v seqlen="$wt_seqlen" -v label="$label" \
    '{
        # annotate
        id=$1; flag=$2; start=$3; cs=$4
        # -------------------------------
        # CS tag MIDS formatting...
        # -------------------------------
        sub("cs:Z:","",cs)
        gsub("=", "", cs)
        # Deletion
        cs_del=""
        split(cs,array,"-")
        for(key in array){
            del="D"
            del_num=match(array[key], "^[a|c|g|t]+")
            for(i=1; i<RLENGTH; i++) del=del"D"
            sub("^[a|c|g|t]+",del,array[key])
            cs_del=cs_del""array[key]
        }
        cs=cs_del
        # Substitution
        gsub("*[a|c|g|t][a|c|g|t]", "S", cs)
        # Insertion
        gsub("+[a|c|g|t]+.", "I", cs)
        # Match
        gsub("[A|C|G|T]", "M", cs)
        # -------------------------------
        # Padding start sites...
        # -------------------------------
        seq=""
        for(i=1; i<start; i++) seq=seq"="
        cs=seq""cs

        # print cs

        # -------------------------------
        # Padding or trimming end sites...
        # -------------------------------
        # もし短い場合は"="を足す。長い場合は切る：
        if(length(cs)<seqlen){
            seq=""
            for(i=1; i<=seqlen-length(cs); i++) seq=seq"="
            cs=cs""seq
        } else {
            cs=substr(cs,1,seqlen)
        }

        print id, cs, label
    }' |
    sed "s/ /\t/g" |
    cat - \
    >> "$output_mids"
done

cat "$output_mids" |
grep -e simulate -e barcode32 |
sed "s/barcode32/wt_simulated/g" |
cat - > DAJIN_sim.txt

cat "$output_mids" |
grep -v simulate |
cat - > DAJIN_real.txt


cat "$output_mids" |
grep -e simulate -e barcode32 |
awk 'BEGIN{OFS="\t"} 
    {print $0, ++names[$3]}' |
sort -k 4,4n |
cut -f 1,2,3 |
sed "s/barcode32/wt_simulated/g" |
head -n 100000 |
cat - > DAJIN_sim_test.txt

cat "$output_mids" |
# grep -v simulate |
grep -e barcode13 -e barcode32 |
# awk 'BEGIN{OFS="\t"} 
#     {print $0, ++names[$3]}' |
# sort -k 4,4n |
# cut -f 1,2,3 |
# head -n 3000 |
cat - > DAJIN_real_test.txt

cat DAJIN_real_test.txt |
grep -e 8ca3032a-4343-4a64-9393-7d9b35428c36 \
    -e 2e2de689-9f23-454b-ae69-46f004d10200 \
> test.txt
# cat DAJIN_new_mids.txt | cut -f 2 | awk '{print length}' | sort | uniq -c


# cat tmp_mids | awk '{print substr($3, 739, 1)}' | sort | uniq -c
# cat tmp_mids | cut -f 3 | awk '{print length}' | sort | uniq -c
# cat tmp_mids | awk 'length($3)==2845' | head

# ============================================================
# テスト用にBAMファイル生成
# ============================================================

# BAM="DAJIN_Report/bam/target#c140G_C_simulated_aligned_reads.bam"
# samtools view -@ "$threads" -H "$BAM" |
# grep -e "^@" > tmp_header

# samtools view  -@ "$threads" "$BAM" |
# sort -k 1,1 |
# join - tmp_id |
# sed "s/ /\t/g" \
# >> tmp_header

# samtools sort tmp_header > tmp_pm.bam
# samtools index tmp_pm.bam