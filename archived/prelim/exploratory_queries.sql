

SELECT * from drug_activity 
WHERE target_id IN (SELECT target_id from drug_activity WHERE drug="vemurafenib");

# Drug-based search to label target
SELECT CONCAT('drug_activity_id_',target_id,'_',act_type,'_',action_type) AS results 
FROM drug_activity 
WHERE drug="sorafenib";

# Drug_to_protein target vector
SELECT CONCAT('drug_activity_id_',target_id) AS results 
FROM drug_activity 
WHERE drug="sorafenib";

SELECT * from pathway 
WHERE (protein_id or target_id) in (select id FROM kinases);

SELECT * FROM pathway 
WHERE protein_id in (select id from kinases) and protein_id = ;

SELECT * from expression WHERE protein_id in (select id from kinases);
SELECT gat_id,avg(value) from gene_attribute WHERE protein_id="1636" group by gat_id;





SELECT CONCAT('ppi','_',(
	CASE WHEN protein1_id = 1636 THEN protein2_id 
	ELSE protein1_id END)) AS ppi_protein_id
FROM ppi WHERE protein1_id=1636 OR protein2_id=1636
ORDER BY ppi_protein_id;


SELECT uniprot as uniprot_id
FROM protein WHERE id="1636"
ORDER BY id;

SELECT idgfam as target_idgfam
FROM target
WHERE id="1636";

SELECT id as target2disease_id
FROM target2disease
WHERE target_id="1636";


SELECT id as phenotype_id
FROM phenotype WHERE (target_id="1636" OR protein_id="1636") AND ptype="JAX/MGI Human Ortholog Phenotype";

SELECT id as pathway_id
FROM pathway WHERE (target_id="1636" OR protein_id="1636");

SELECT panther_class_id as p2pc_panther_class_id
FROM p2pc WHERE protein_id=1636;

SELECT CONCAT('goa_id_',(Replace(go_id,'GO:',''))) AS results
FROM goa WHERE protein_id=1636;

SELECT gat_id,gene_attribute.name,gene_attribute_type.name
FROM gene_attribute
JOIN gene_attribute_type ON
gene_attribute_type.id=gene_attribute.gat_id 
WHERE protein_id=1636;

#SELECT name
#FROM gene_attribute 
#WHERE protein_id=1636 AND value=1 AND gat_id=1;# lets go with reactome pathways


SELECT DISTINCT(name)
FROM gene_attribute;

SELECT distinct(type)
FROM feature WHERE protein_id=1636;