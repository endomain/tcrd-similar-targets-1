####################################################################################
# Gets all drug names known in TCRD that targets kinases
####################################################################################
select distinct(drug) from drug_activity where target_id IN (SELECT id from target where idgfam="Kinase");

select target_id,count(*) 
FROM drug_activity 
WHERE target_id IN (SELECT id from target where idgfam="Kinase") 
GROUP BY target_id;
 
select name,tdl,target_id, count(*)
FROM drug_activity 
INNER JOIN target ON drug_activity.target_id=target.id
WHERE target_id IN (SELECT id from target where idgfam="Kinase")
group by name;

select * from target where id=18810;

select * from drug_activity where drug="Phenacetin";

select name,id from target where idgfam="Kinase";