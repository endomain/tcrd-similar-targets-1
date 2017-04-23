CREATE TEMPORARY TABLE kinases (id int) as
(SELECT id from target where idgfam="Kinase");

SELECT * from kinases;

##################
#
# TABLE: 	chembl_activity
# STATUS: 	DONE
#
# 1) cmpd_chemblid
# 2) act_type
# *) act_value
##################
SELECT 
	CONCAT(
		"chembl_activity",
		'_',
		cmpd_chemblid,
		'_',
		act_type
	) AS result
FROM chembl_activity 
WHERE target_id=1636;

##################
#
# TABLE: 	compartment
# STATUS: 	DONE
#
# 1) go_id
# *) act_value
##################
SELECT DISTINCT(CONCAT('drug_activity','_',REPLACE(go_term,' ',"_")))
FROM compartment 
WHERE protein_id=1636;

##################
#
# TABLE: 	expression
# STATUS: 	DONE
#
# 1) tissue
# 2) qual_value
#
##################
# expressions 
SELECT 
	CONCAT(
		'expression_',
		Replace(tissue,' ','_'),
		'_',
		Replace(qual_value,' ','_')
	) AS expression
FROM expression 
WHERE protein_id=1636 
#	AND etype="HPM Protein"
	GROUP BY tissue;
	

##################
#
# TABLE: 	generif
# STATUS: 	DONE
#
# 1) getting pubmed_id citations
#
##################
SELECT (CONCAT('generif_pubmed','_',pubmed_ids))
FROM generif 
WHERE protein_id=1636;

##################
#
# TABLE: 	gene_attribute
# STATUS: 	DONE, add more as needed
# 1) Use KEGG pathway for trial
# 2) What else?
##################
SELECT REPLACE(CONCAT(
	'gene_attribute','_',
	REPLACE(gene_attribute.name,
	CONCAT('/',gene_attribute_type.name),"")),' ','_') as results
FROM gene_attribute
JOIN gene_attribute_type ON
gene_attribute_type.id=gene_attribute.gat_id 
WHERE 
	gat_id=103 AND 
	protein_id=1636;

##################
#
# TABLE: 	panther_class, p2pc
# STATUS: 	DONE
# 1) pcid
#
##################
SELECT CONCAT('panther_class_id_',pcid) AS results 
FROM panther_class 
WHERE ID IN (SELECT panther_class_id from p2pc WHERE protein_id="1636");

##################
#
# TABLE: 	pathway
# STATUS: 	DONE
# 1) pwtype
# 2) id_in_source
#
##################
SELECT (CONCAT('pathway','_',pwtype,"_",id_in_source)) as results
FROM pathway 
WHERE protein_id=1636 AND pwtype="Reactome";

##################
#
# TABLE: 	phenotype
# STATUS: 	DONE
# 1) term_id
#
##################
SELECT (CONCAT('phenotype','_',REPLACE(term_id,':','_'))) as results
FROM phenotype 
WHERE protein_id=1636 AND ptype="JAX/MGI Human Ortholog Phenotype";

##################
#
# TABLE: 	ppi
# STATUS: 	DONE
# 1) the protein1_id and protein2_id are related but only one directional
#
##################
SELECT CONCAT('ppi','_',(
	CASE WHEN protein1_id = 1636 THEN protein2_id 
	ELSE protein1_id END)) AS ppi_protein_id
FROM ppi WHERE protein1_id=1636 OR protein2_id=1636
ORDER BY ppi_protein_id;

##################
#
# TABLE: 	protein
# STATUS: 	DONE
#
# EXPLANATION: It does not have any associations - protein name mapping data only.
# *) seq - has potential use if we can use BLAST on everything, but its a lot more work
#
##################


##################
#
# TABLE: 	target2disease
# STATUS: 	DONE
#
# 1) name
#
##################
SELECT CONCAT('target2disease','_',REPLACE(name,' ','_')) as target2disease_id
FROM target2disease
WHERE target_id=1636;



DROP TABLE kinases;

