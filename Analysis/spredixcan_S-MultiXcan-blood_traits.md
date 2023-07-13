
## #@'1' WhiteBloodCellcount_gwas_data.txt

```R
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/WhiteBloodcellCout
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/WhiteBloodcellCout

#data file: WhiteBloodCellcount_gwas_data.txt


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &
#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/WhiteBloodCellcount_gwas_data.txt > WhiteBloodCellcount_gwas_data_tab_new.txt &



#step: GWAS parsing.py

#input file: WhiteBloodCellcount_gwas_data.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/WhiteBloodCellcount_gwas_data_tab_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 563946  --insert_value n_cases 563946 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt


```



## #@'2' monocytecount_gwas_data.txt

```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MonocyteCount
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MonocyteCount

#data file: monocytecount_gwas_data.txt


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/monocytecount_gwas_data.txt > monocytecount_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file: monocytecount_gwas_data.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/monocytecount_gwas_data_tab_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 563946   --insert_value n_cases 563946 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt


```

## #@'3' Lymphocytecount3_gwas_data.txt



### #pbs

```shell
#PBS -N smultixcan_lymphocount.sh.pbs
#PBS -q workq         
#PBS -l mem=120G
#PBS -o /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocyteCount/standard.out
#PBS -e /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocyteCount/standard.err
 
/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocyteCount/smultixcan_lymphocount.sh


```

### #S-MultiXcan



```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocyteCount
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocyteCount

#data file: monocytecount_gwas_data.txt


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Lymphocytecount3_gwas_rawdata.txt > lymphocytecount_gwas_data_tab_new.txt &




#step: GWAS parsing.py

#input file: monocytecount_gwas_data.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/lymphocytecount_gwas_data_tab_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map ID variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map ES effect_size \
-output_column_map  p pvalue \
-output_column_map CHROM chromosome \
--chromosome_format \
-output_column_map POS position \
--insert_value sample_size 171643   --insert_value n_cases 171643 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt

```

## #@'4' neutrophilcount_gwas_data.txt

```shell

#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/NeutrophilCount
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/NeutrophilCount


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/neutrophilcount_gwas_data.txt > neutrophilcount_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file:neutrophilcount_gwas_data_tab_new.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/neutrophilcount_gwas_data_tab_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 563946   --insert_value n_cases 563946 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt
```

## #PBS 提交任务

```shell
#可执行文件
chmod 777 pbs_permu_for401to500.r          
#运行PBS
qsub pbs_permu_for401to500.r       

#查看空闲节点情况
pestat 

#任务提交后，用户如果要知道任务的当前运行状态，可以通过 qstat 命令查询。
qstat 
```



## #@'5'  LymphocytePercent_gwas_data.txt

### #pbs

```shell
#PBS -N smultixcan.sh.pbs
#PBS -q workq         
#PBS -l nodes=node02:ppn=5
#PBS -o /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocytePercent/standard.out
#PBS -e /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocytePercent/standard.err
 
/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocytePercent/smultixcan.sh




#PBS -l nodes=node02:ppn=5
fat
mem=400
#PBS -l mem=150G
```

### #S-MultiXcan

```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocytePercent
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/LymphocytePercent


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/LymphocytePercent_gwas_data.txt > LymphocytePercent_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file:LymphocytePercent_gwas_data_tab_new.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/LymphocytePercent_gwas_data_tab_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 171748   --insert_value n_cases 171748  \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt



```

## #@'6'  MeanCorpusVolume_gwas_data.txt

### #pbs

```shell
#PBS -N smultixcan.sh.pbs
#PBS -q workq         
#PBS -l mem=150G
#PBS -o /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpusVolume/standard.out
#PBS -e /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpusVolume/standard.err
 
/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpusVolume/smultixcan.sh


```

### #S-MultiXcan

```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpusVolume
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpusVolume


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/MeanCorpusVolume_gwas_data.txt > MeanCorpusVolume_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file:MeanCorpusVolume_gwas_data.txt 
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/MeanCorpusVolume_gwas_data_tab_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 121047   --insert_value n_cases 121047 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt



```



## #@'7'  eosinophilcount_gwas_data.txt

### #pbs

```shell
#PBS -N smultixcan.sh.pbs
#PBS -q workq         
#PBS -l mem=150G
#PBS -o /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/EosinophilCount/standard.out
#PBS -e /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/EosinophilCount/standard.err
 
/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/EosinophilCount/smultixcan.sh


```

### #S-MultiXcan

```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/EosinophilCount
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/EosinophilCount


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/eosinophilcount_gwas_data.txt > eosinophilcount_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file: eosinophilcount_gwas_data_tab_new.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/eosinophilcount_gwas_data_tab_new.txt  \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 563946   --insert_value n_cases 563946 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt



```



## #@'8' basophilcount_gwas_data.txt



### #pbs

```shell
#PBS -N smultixcan.sh.pbs
#PBS -q workq         
#PBS -l mem=150G
#PBS -o /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/BasophilCount/standard.out
#PBS -e /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/BasophilCount/standard.err
 
/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/BasophilCount/smultixcan.sh


```

### #S-MultiXcan

```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/BasophilCount
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/BasophilCount


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/basophilcount_gwas_data.txt > basophilcount_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file: basophilcount_gwas_data.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/basophilcount_gwas_data_tab_new.txt  \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 563946   --insert_value n_cases 563946 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt



```



## #@'9'  MeanCorpuscularHemoglobin_gwas_data.txt



### #pbs

```shell
#PBS -N smultixcan.sh.pbs
#PBS -q workq         
#PBS -l mem=150G
#PBS -o /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpuscularHemoglobin/standard.out
#PBS -e /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpuscularHemoglobin/standard.err
 
/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpuscularHemoglobin/smultixcan.sh


```

### #S-MultiXcan

```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpuscularHemoglobin
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/MeanCorpuscularHemoglobin


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/MeanCorpuscularHemoglobin_gwas_data.txt > MeanCorpuscularHemoglobin_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file:MeanCorpuscularHemoglobin_gwas_data_tab_new.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/MeanCorpuscularHemoglobin_gwas_data_tab_new.txt  \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 126151   --insert_value n_cases 126151 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt



```





## #@'10' Hemoglobinconcen_gwas_data.txt



### #pbs

```shell
#PBS -N smultixcan.sh.pbs
#PBS -q workq         
#PBS -l mem=150G
#PBS -o /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/Hemoglobinconcen/standard.out
#PBS -e /share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/Hemoglobinconcen/standard.err
 
/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/Hemoglobinconcen/smultixcan.sh


```

### #S-MultiXcan

```shell
#!/usr/bin/env bash

export DATA=/share/home/mayl/miniconda3/DATABASE
export GWAS_TOOLS=/share/home/mayl/miniconda3/summary-gwas-imputation/src
export METAXCAN=/share/home/mayl/miniconda3/MetaXcan/software
export GWAS_DATA=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/Hemoglobinconcen
export OUTPUT=/share/pub/mayl/scPagwas/S-MultiXcan_blood_traits/Hemoglobinconcen


#first processing the raw GWAS summary file format
#remove all NA rsid variants
# nohup grep -v "NA" COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt > COVID19_HGI_B2_ALL_leave_23andme_b37_new.txt &

#@' GWAS结果文件的染色体用数值表示，不要用chr+数值的方式
#@' GWAS文件格式必须是gz压缩格式，分隔符为tab
#@' sed "s/ /\t/g" /share/pub/dengcy/GWAS_Multiomics/compare/bloodtraits/Hemoglobinconcen_gwas_data.txt > Hemoglobinconcen_gwas_data_tab_new.txt &


#step: GWAS parsing.py

#input file:Hemoglobinconcen_gwas_data_tab_new.txt
#output file:meta_GWAS_COVID19_4th_harmonized.txt
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $GWAS_DATA/Hemoglobinconcen_gwas_data_tab_new.txt \
-liftover $DATA/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata $DATA/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map REF non_effect_allele \
-output_column_map ALT effect_allele \
-output_column_map beta effect_size \
-output_column_map  p pvalue \
-output_column_map chrom chromosome \
--chromosome_format \
-output_column_map pos position \
--insert_value sample_size 149861   --insert_value n_cases 149861 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt


#step2: GWAS summary stats imputation

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#output file:COVID_round_4th_GWAS_chr*****
for chr in {1..22}; do
  for batch in {0..9}; do
    python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/data/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
    -parquet_genotype $DATA/data/reference_panel_1000G/chr${chr}.variants.parquet \
    -parquet_genotype_metadata $DATA/data/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome ${chr} \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch ${batch} \
    --standardise_dosages \
    -output $OUTPUT/summary_imputation_1000G/COVID_round_4th_GWAS_chr${chr}_sb${batch}_reg0.1_ff0.01_by_region.txt
  done
done

###Imputation post-processing

#input file: meta_GWAS_COVID19_4th_harmonized.txt
#change pattern:"COVID_round_4th_GWAS_*" 
#output file:imputed_COVID_GWAS_round_4th.txt
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/meta_GWAS_COVID19_4th_harmonized.txt \
-folder $OUTPUT/summary_imputation_1000G \
-pattern "COVID_round_4th_GWAS_*" \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt


#step3: S-PrediXcan mashr eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file: COVID_round_4th_GWAS_${db##*/}.csv
#for name in `ls Database/data/models/eqtl/mashr/*db`; do echo ${name%.db}.txt.gz; done
for db in `ls /share/home/mayl/miniconda3/DATABASE/eqtl/mashr/*db`; do
  python $METAXCAN/SPrediXcan.py \
  --gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
  --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
  --model_db_path ${db} \
  --covariance ${db%.db}.txt.gz \
  --keep_non_rsid --additional_output --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT/spredixcan/eqtl/mashr/COVID_round_4th_GWAS_${db##*/}.csv
done


#step4: S-MultiXcan eqtl

#input file: imputed_COVID_GWAS_round_4th.txt 
#output file:COVID_round_4th_GWAS_${db##*/}.csv
#----------------------------------------------------------------------
python $METAXCAN/SMulTiXcan.py \
--models_folder $DATA/eqtl/mashr \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance $DATA/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder $OUTPUT/spredixcan/eqtl/mashr/ \
--metaxcan_filter "COVID_round_4th_GWAS_mashr_(.*).db.csv" \
--metaxcan_file_name_parse_pattern "(.*)_mashr_(.*).db.csv" \
--gwas_file  $OUTPUT/processed_summary_imputation_1000G/imputed_COVID_GWAS_round_4th.txt \
--snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore --keep_non_rsid --model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output $OUTPUT/smultixcan/eqtl/COVID_GWAS_round_4_mashr_smultixcan_eqtl.txt
```



























































