#!/bin/bash

SHOW TABLES;
# TRUNCATE TABLE SimilarSequences;

# DROP TABLE SimilarSequences00;

# CREATE TABLE SimilarSequences00 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
# CREATE TABLE SimilarSequences01 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#
# CREATE TABLE SimilarSequences02 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#
# CREATE TABLE SimilarSequences03 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#
#
# CREATE TABLE SimilarSequences04 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#
#
# CREATE TABLE SimilarSequences05 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#
#
# CREATE TABLE SimilarSequences06 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#
#
# CREATE TABLE SimilarSequences07 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
# CREATE TABLE SimilarSequences08 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
# CREATE TABLE SimilarSequences09 (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#


# DESCRIBE SimilarSequences00;


# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss00" REPLACE INTO TABLE SimilarSequences00 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss01" REPLACE INTO TABLE SimilarSequences01 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss02" REPLACE INTO TABLE SimilarSequences02 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss03" REPLACE INTO TABLE SimilarSequences03 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss04" REPLACE INTO TABLE SimilarSequences04 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss05" REPLACE INTO TABLE SimilarSequences05 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss06" REPLACE INTO TABLE SimilarSequences06 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss07" REPLACE INTO TABLE SimilarSequences07 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss08" REPLACE INTO TABLE SimilarSequences08 FIELDS TERMINATED BY "\t";
# LOAD DATA LOCAL INFILE "/home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/output/blast_load/split_ss09" REPLACE INTO TABLE SimilarSequences09 FIELDS TERMINATED BY "\t";

# DESCRIBE SimilarSequences01;

# CREATE TABLE SimilarSequences (
#  QUERY_ID                 VARCHAR(60),
#  SUBJECT_ID               VARCHAR(60),
#  QUERY_TAXON_ID           VARCHAR(40),
#  SUBJECT_TAXON_ID         VARCHAR(40),
#  EVALUE_MANT              FLOAT,
#  EVALUE_EXP               int(11),
#  PERCENT_IDENTITY         FLOAT,
#  PERCENT_MATCH            FLOAT
# );
#
# CREATE INDEX ss_qtaxexp_ix
# ON SimilarSequences(query_id, subject_taxon_id, evalue_exp, evalue_mant,query_taxon_id, subject_id) ;
#
# CREATE INDEX ss_seqs_ix
# ON SimilarSequences(query_id, subject_id,evalue_exp, evalue_mant, percent_match) ;

# SHOW VARIABLES LIKE "connect_timeout";
# SHOW VARIABLES LIKE "wait_timeout";
# SHOW VARIABLES LIKE "net_read_timeout";
# SHOW VARIABLES LIKE "bulk_insert_buffer_size";
# SHOW VARIABLES LIKE "key_buffer_size";

# mysql> \q
# $ mysql -u root -p  #pluteus123
# mysql> SET GLOBAL connect_timeout = 6000;
# mysql> SET GLOBAL bulk_insert_buffer_size = 536870912‬;
# mysql> SET GLOBAL key_buffer_size = 1073741824;
# mysql> \q
# $ mysql -u orthomcl -p  #pluteus123
#
# SET LOCAL wait_timeout = 600000;
# SET LOCAL net_read_timeout = 6000;

# INSERT INTO SimilarSequences
# SELECT *
# FROM SimilarSequences00
# UNION ALL
# SELECT *
# FROM SimilarSequences01;
INSERT INTO SimilarSequences
SELECT *
FROM SimilarSequences02
UNION ALL
SELECT *
FROM SimilarSequences03;
# UNION ALL
# SELECT *
# FROM SimilarSequences04
# UNION ALL
# SELECT *
# FROM SimilarSequences05
# UNION ALL
# SELECT *
# FROM SimilarSequences06
# UNION ALL
# SELECT *
# FROM SimilarSequences07
# UNION ALL
# SELECT *
# FROM SimilarSequences08
# UNION ALL
# SELECT *
# FROM SimilarSequences09;

# TRUNCATE TABLE SimilarSequences01;






# sudo -S perl5.22.1 $ORTHMCL_PIP/scripts/orthomcl-pipeline_onlymysql.pl -i input/ -o output/ -m orthomcl.config -c orthomcl-pip.conf --nocompliant
