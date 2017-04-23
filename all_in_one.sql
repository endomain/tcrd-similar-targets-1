SELECT 
	CONCAT(
		"chembl_activity",
		'_',
		cmpd_chemblid,
		'_',
		act_type
	) AS results
	FROM chembl_activity as results
	WHERE target_id=1636
	
UNION

SELECT DISTINCT(CONCAT('drug_activity','_',REPLACE(go_term,' ',"_"))) as results
	FROM compartment 
	WHERE protein_id=1636
	
UNION

SELECT 
	CONCAT(
		'expression_',
		Replace(tissue,' ','_'),
		'_',
		Replace(qual_value,' ','_')
	) AS expression
	FROM expression 
	WHERE protein_id=1636 
	# AND etype="HPM Protein"
	GROUP BY tissue	

UNION

SELECT (CONCAT('generif_pubmed','_',pubmed_ids))
	FROM generif 
	WHERE protein_id=1636

UNION


SELECT REPLACE(CONCAT('gene_attribute','_',REPLACE(gene_attribute.name,CONCAT('/',gene_attribute_type.name),"")),' ','_') as results
	FROM gene_attribute
	JOIN gene_attribute_type ON gene_attribute_type.id=gene_attribute.gat_id 
	WHERE 
		gat_id=103 AND 
		protein_id=1636

UNION

SELECT CONCAT('panther_class_id_',pcid) AS results 
	FROM panther_class 
	WHERE ID IN (SELECT panther_class_id from p2pc WHERE protein_id="1636")

UNION
 
SELECT (CONCAT('pathway','_',pwtype,"_",id_in_source)) as results
	FROM pathway 
	WHERE protein_id=1636 AND pwtype="Reactome"
	
UNION

SELECT (CONCAT('phenotype','_',REPLACE(term_id,':','_'))) as results
	FROM phenotype 
	WHERE protein_id=1636 AND ptype="JAX/MGI Human Ortholog Phenotype"
	
UNION

SELECT CONCAT('ppi','_',(
	CASE 
		WHEN protein1_id = 1636 THEN protein2_id 
		ELSE protein1_id END)) AS ppi_protein_id
	FROM ppi WHERE protein1_id=1636 OR protein2_id=1636

UNION


SELECT CONCAT('target2disease','_',REPLACE(name,' ','_')) as target2disease_id
FROM target2disease
WHERE target_id=1636;