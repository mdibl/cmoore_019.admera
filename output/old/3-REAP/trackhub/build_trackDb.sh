#!/bin/bash

# OUTPUT FILE
OUT=trackDb.txt

# Overwrite existing file
> $OUT

# ROOT DIRECTORY FOR EACH TYPE
BASE_DIR=~/Documents/mdibl-cnobrega/projects/ClaireMoore/cmoore_019.admera/trackhub/3REAP_hub/bigWigs

echo "##############################################" >> $OUT
echo "# Composite: All aligned reads" >> $OUT
echo "##############################################" >> $OUT

echo "track mapped" >> $OUT
echo "compositeTrack on" >> $OUT
echo "priority 1.0" >> $OUT
echo "shortLabel mapped" >> $OUT
echo "longLabel All aligned reads" >> $OUT
echo "type bigWig" >> $OUT
echo "visibility full" >> $OUT
echo "subGroup1 condition condition Control=Control Treatment=Treatment" >> $OUT
echo "subGroup2 day day day3=day3 day5=day5" >> $OUT
echo "subGroup3 strand strand mi=Mirror plus=Plus minus=Minus" >> $OUT
echo "sortOrder condition=+ day=+ strand=+" >> $OUT
echo "dimensionX condition" >> $OUT
echo "dimensionY day" >> $OUT
echo "dragAndDrop subtracks" >> $OUT
echo "" >> $OUT

for file in $BASE_DIR/mapped_read/*.bw; do
  fname=$(basename "$file")
  strand=$(echo "$fname" | cut -d. -f1)
  sample=$(echo "$fname" | cut -d. -f2)
  condition=$(echo "$sample" | cut -d_ -f1)
  day=$(echo "$sample" | cut -d_ -f2)
  ucsc_day=""
  if [[ "$day" == "3day" ]]; then
	  ucsc_day="day3"
  elif [[ "$day" == "5day" ]]; then
	  ucsc_day="day5"
  else
	  ucsc_day=$day  # fallback
  fi
  rep=$(echo "$sample" | cut -d_ -f3)
  track_name="mapped_${sample}_${strand}"
  short_label="${condition:0:1}_${day:0:2}_r${rep:3:1}_${strand}"
  cap_condition=$(echo "$condition" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')

  cat >> $OUT <<EOF
track $track_name
parent mapped
bigDataUrl http://3.83.106.170/3REAP_hub/hg38/mapped_read/${strand}.${sample}.mapped_read.bw
shortLabel $short_label
longLabel $cap_condition ${day} ${rep} ${strand}
type bigWig
subGroups condition=$cap_condition day=$ucsc_day strand=$strand
visibility full
alwaysZero on
autoScale on
graphTypeDefault bar
maxHeightPixels 120:80:40
color 128,128,128

EOF
done

echo "##############################################" >> $OUT
echo "# Composite: PASS reads" >> $OUT
echo "##############################################" >> $OUT

echo "track PASSreads" >> $OUT
echo "compositeTrack on" >> $OUT
echo "priority 2.0" >> $OUT
echo "shortLabel PASSreads" >> $OUT
echo "longLabel PASS reads" >> $OUT
echo "type bigWig" >> $OUT
echo "visibility full" >> $OUT
echo "subGroup1 condition condition Control=Control Treatment=Treatment" >> $OUT
echo "subGroup2 day day day3=day3 day5=day5" >> $OUT
echo "subGroup3 strand strand mi=Mirror plus=Plus minus=Minus" >> $OUT
echo "sortOrder condition=+ day=+ strand=+" >> $OUT
echo "dimensionX condition" >> $OUT
echo "dimensionY day" >> $OUT
echo "dragAndDrop subtracks" >> $OUT
echo "" >> $OUT

for file in $BASE_DIR/PASS_bw_LAP24/*.bw; do
  fname=$(basename "$file")
  strand=$(echo "$fname" | cut -d. -f1)
  sample=$(echo "$fname" | cut -d. -f2)
  condition=$(echo "$sample" | cut -d_ -f1)
  day=$(echo "$sample" | cut -d_ -f2)
  ucsc_day=""
  if [[ "$day" == "3day" ]]; then
	  ucsc_day="day3"
  elif [[ "$day" == "5day" ]]; then
	  ucsc_day="day5"
  else
	  ucsc_day=$day  # fallback
  fi
  rep=$(echo "$sample" | cut -d_ -f3)
  track_name="PASSreads_${sample}_${strand}"
  short_label="${condition:0:1}_${day:0:2}_r${rep:3:1}_${strand}"
  cap_condition=$(echo "$condition" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')

  cat >> $OUT <<EOF
track $track_name
parent PASSreads
bigDataUrl http://3.83.106.170/3REAP_hub/hg38/PASS_bw_LAP24/${strand}.${sample}.PASS_bw_LAP24.bw
shortLabel $short_label
longLabel $cap_condition ${day} ${rep} ${strand}
type bigWig
subGroups condition=$cap_condition day=$ucsc_day strand=$strand
visibility full
alwaysZero on
autoScale on
graphTypeDefault bar
maxHeightPixels 120:80:40
color 60,60,255

EOF
done

echo "##############################################" >> $OUT
echo "# Composite: PASS LAPs" >> $OUT
echo "##############################################" >> $OUT

echo "track LAPs" >> $OUT
echo "compositeTrack on" >> $OUT
echo "priority 3.0" >> $OUT
echo "shortLabel LAPs" >> $OUT
echo "longLabel PASS LAPs" >> $OUT
echo "type bigWig" >> $OUT
echo "visibility full" >> $OUT
echo "subGroup1 condition condition Control=Control Treatment=Treatment" >> $OUT
echo "subGroup2 day day day3=day3 day5=day5" >> $OUT
echo "subGroup3 strand strand mi=Mirror plus=Plus minus=Minus" >> $OUT
echo "sortOrder condition=+ day=+ strand=+" >> $OUT
echo "dimensionX condition" >> $OUT
echo "dimensionY day" >> $OUT
echo "dragAndDrop subtracks" >> $OUT
echo "" >> $OUT

for file in $BASE_DIR/PASS_bw_LAP24_positon/*.bw; do
  fname=$(basename "$file")
  strand=$(echo "$fname" | cut -d. -f1)
  sample=$(echo "$fname" | cut -d. -f2)
  condition=$(echo "$sample" | cut -d_ -f1)
  day=$(echo "$sample" | cut -d_ -f2)
  ucsc_day=""
  if [[ "$day" == "3day" ]]; then
	  ucsc_day="day3"
  elif [[ "$day" == "5day" ]]; then
	  ucsc_day="day5"
  else
	  ucsc_day=$day  # fallback
  fi
  rep=$(echo "$sample" | cut -d_ -f3)
  track_name="LAPs_${sample}_${strand}"
  short_label="${condition:0:1}_${day:0:2}_r${rep:3:1}_${strand}"
  cap_condition=$(echo "$condition" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')

  cat >> $OUT <<EOF
track $track_name
parent LAPs
bigDataUrl http://3.83.106.170/3REAP_hub/hg38/PASS_bw_LAP24_positon/${strand}.${sample}.PASS_bw_LAP24_positon.bw
shortLabel $short_label
longLabel $cap_condition ${day} ${rep} ${strand}
type bigWig
subGroups condition=$cap_condition day=$ucsc_day strand=$strand
visibility full
alwaysZero on
autoScale on
graphTypeDefault bar
maxHeightPixels 120:80:40
color 60,60,180

EOF
done

echo "##############################################" >> $OUT
echo "# Composite: Detected PASs" >> $OUT
echo "##############################################" >> $OUT

echo "track PAS" >> $OUT
echo "compositeTrack on" >> $OUT
echo "priority 4.0" >> $OUT
echo "shortLabel PAS" >> $OUT
echo "longLabel Detected PASs" >> $OUT
echo "type bigWig" >> $OUT
echo "visibility dense" >> $OUT
echo "subGroup1 condition condition Control=Control Treatment=Treatment" >> $OUT
echo "subGroup2 day day day3=day3 day5=day5" >> $OUT
echo "subGroup3 strand strand mi=Mirror plus=Plus minus=Minus" >> $OUT
echo "sortOrder condition=+ day=+ strand=+" >> $OUT
echo "dimensionX condition" >> $OUT
echo "dimensionY day" >> $OUT
echo "dragAndDrop subtracks" >> $OUT
echo "" >> $OUT

for file in $BASE_DIR/PASS_bw_LAP24_PAS/*.bw; do
  fname=$(basename "$file")
  strand=$(echo "$fname" | cut -d. -f1)
  sample=$(echo "$fname" | cut -d. -f2)
  condition=$(echo "$sample" | cut -d_ -f1)
  day=$(echo "$sample" | cut -d_ -f2)
  ucsc_day=""
  if [[ "$day" == "3day" ]]; then
	  ucsc_day="day3"
  elif [[ "$day" == "5day" ]]; then
	  ucsc_day="day5"
  else
	  ucsc_day=$day  # fallback
  fi
  rep=$(echo "$sample" | cut -d_ -f3)
  track_name="PAS_${sample}_${strand}"
  short_label="${condition:0:1}_${day:0:2}_r${rep:3:1}_${strand}"
  cap_condition=$(echo "$condition" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')

  cat >> $OUT <<EOF
track $track_name
parent PAS
bigDataUrl http://3.83.106.170/3REAP_hub/hg38/PASS_bw_LAP24_PAS/${strand}.${sample}.PASS_bw_LAP24_PAS.bw
shortLabel $short_label
longLabel $cap_condition ${day} ${rep} ${strand}
type bigWig
subGroups condition=$cap_condition day=$ucsc_day strand=$strand
visibility dense
alwaysZero on
autoScale on
graphTypeDefault bar
maxHeightPixels 120:80:40
color 220,0,0

EOF
done
