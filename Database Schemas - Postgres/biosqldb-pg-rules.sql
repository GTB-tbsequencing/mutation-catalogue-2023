
-- This is to solve a problem arising from how transactions are implemented
-- in Postgres as opposed to, e.g., Oracle and InnoDB (MySQL). In short, the
-- difference is that in the latter RDBMSs' implementation, if a particular
-- statement within a transaction fails, the preceding (and possibly
-- subsequent) statements are still valid. On commit, all succeeded statements
-- are committed. In Postgres, the failure of a statement invalidates all
-- preceding statements within the same transaction as well as all subsequent,
-- if any.
--
-- This leads to a problem if you program SQL insert and update statements
-- such that presence of the record you attempt to insert is indicated by
-- failure of the statement due to a unique key constraint violation. Even
-- if your code is prepared to handle the failure by e.g. looking up the
-- record, in the case of Postgres this approach cannot work unless you
-- commit every single statement.
--
-- The bioperl-db adaptor code uses the aforementioned approach and is
-- currently dependent on the following support code. If you are not going
-- to use bioperl-db to populate the database, you may comment out all
-- rules, as then they might add another look-up to one already done on the
-- code that you use and hence add unnecessary overhead.
--


CREATE RULE rule_biodatabase_i
       AS ON INSERT TO biodatabase
       WHERE (
             SELECT biodatabase_id FROM biodatabase 
             WHERE name = new.name
             )
       	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_bioentry_dbxref_i
       AS ON INSERT TO bioentry_dbxref
       WHERE (
       	     SELECT dbxref_id FROM bioentry_dbxref
	     WHERE bioentry_id = new.bioentry_id
	     AND   dbxref_id   = new.dbxref_id
	     ) 
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_bioentry_path_i
       AS ON INSERT TO bioentry_path
       WHERE (
       	     SELECT bioentry_relationship_id FROM bioentry_relationship
	     WHERE object_bioentry_id = new.object_bioentry_id
	     AND   subject_bioentry_id= new.subject_bioentry_id
	     AND   term_id	      = new.term_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_bioentry_qualifier_value_i
       AS ON INSERT TO bioentry_qualifier_value
       WHERE (
       	     SELECT bioentry_id FROM bioentry_qualifier_value
	     WHERE bioentry_id = new.bioentry_id
	     AND   term_id     = new.term_id
	     AND   rank	       = new.rank
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_bioentry_reference_i
       AS ON INSERT TO bioentry_reference
       WHERE (
       	     SELECT bioentry_id FROM bioentry_reference 
	     WHERE bioentry_id  = new.bioentry_id
	     AND   reference_id = new.reference_id
	     AND   rank		= new.rank
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_bioentry_relationship_i
       AS ON INSERT TO bioentry_relationship
       WHERE (
       	     SELECT bioentry_relationship_id FROM bioentry_relationship
	     WHERE object_bioentry_id = new.object_bioentry_id
	     AND   subject_bioentry_id= new.subject_bioentry_id
	     AND   term_id	      = new.term_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_biosequence_i
       AS ON INSERT TO biosequence
       WHERE (
             SELECT bioentry_id FROM biosequence 
             WHERE bioentry_id = new.bioentry_id
             )
       	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_comment_i
       AS ON INSERT TO comment
       WHERE (
       	     SELECT comment_id FROM comment
	     WHERE bioentry_id = new.bioentry_id
	     AND   rank	       = new.rank
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_dbxref_i
       AS ON INSERT TO dbxref
       WHERE (
       	     SELECT dbxref_id FROM dbxref
	     WHERE accession = new.accession
	     AND   dbname    = new.dbname
	     AND   version   = new.version
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_dbxref_qualifier_value_i
       AS ON INSERT TO dbxref_qualifier_value
       WHERE (
       	     SELECT dbxref_id FROM dbxref_qualifier_value
	     WHERE dbxref_id = new.dbxref_id
	     AND   term_id   = new.term_id
	     AND   rank	     = new.rank
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_location_i
       AS ON INSERT TO location
       WHERE (
       	     SELECT location_id FROM location
	     WHERE seqfeature_id = new.seqfeature_id
	     AND   rank		 = new.rank
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_location_qualifier_value_i
       AS ON INSERT TO location_qualifier_value
       WHERE (
       	     SELECT location_id FROM location_qualifier_value
	     WHERE location_id = new.location_id
	     AND   term_id     = new.term_id
	     ) 
       	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_ontology_i
       AS ON INSERT TO ontology
       WHERE (
             SELECT ontology_id FROM ontology 
             WHERE name = new.name
             ) 
       	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_reference_i1
       AS ON INSERT TO reference
       WHERE (
             SELECT reference_id FROM reference 
             WHERE crc = new.crc
             ) 
       	     IS NOT NULL
       DO INSTEAD NOTHING
;
CREATE RULE rule_reference_i2
       AS ON INSERT TO reference
       WHERE (
             SELECT reference_id FROM reference
             WHERE dbxref_id = new.dbxref_id
             )
       	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_seqfeature_i
       AS ON INSERT TO seqfeature
       WHERE (
       	     SELECT seqfeature_id FROM seqfeature 
	     WHERE bioentry_id    = new.bioentry_id
	     AND   type_term_id   = new.type_term_id
	     AND   source_term_id = new.source_term_id
	     AND   rank		  = new.rank
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_seqfeature_dbxref_i
       AS ON INSERT TO seqfeature_dbxref
       WHERE (	    
       	     SELECT seqfeature_id FROM seqfeature_dbxref
	     WHERE seqfeature_id = new.seqfeature_id
	     AND   dbxref_id	 = new.dbxref_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_seqfeature_path_i
       AS ON INSERT TO seqfeature_path
       WHERE (
       	     SELECT subject_seqfeature_id FROM seqfeature_path
	     WHERE object_seqfeature_id = new.object_seqfeature_id
	     AND   subject_seqfeature_id= new.subject_seqfeature_id
	     AND   term_id		= new.term_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_seqfeature_qualifier_value_i
       AS ON INSERT TO seqfeature_qualifier_value
       WHERE (
       	     SELECT seqfeature_id FROM seqfeature_qualifier_value
	     WHERE seqfeature_id = new.seqfeature_id
	     AND   term_id	 = new.term_id
	     AND   rank		 = new.rank
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_seqfeature_relationship_i
       AS ON INSERT TO seqfeature_relationship
       WHERE (
       	     SELECT subject_seqfeature_id FROM seqfeature_relationship
	     WHERE object_seqfeature_id = new.object_seqfeature_id
	     AND   subject_seqfeature_id= new.subject_seqfeature_id
	     AND   term_id		= new.term_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

-- CREATE RULE rule_taxon_i
--        AS ON INSERT TO taxon
--        WHERE (
--              SELECT taxon_id FROM taxon 
--              WHERE ncbi_taxon_id = new.ncbi_taxon_id
--              )
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

CREATE RULE rule_taxon_name_i
       AS ON INSERT TO taxon_name
       WHERE (
       	     SELECT taxon_id FROM taxon_name
	     WHERE taxon_id   = new.taxon_id
	     AND   name	      = new.name
	     AND   name_class = new.name_class
	     ) 
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_term_i1
       AS ON INSERT TO term
       WHERE (
             SELECT term_id FROM term
             WHERE identifier = new.identifier
             )
       	     IS NOT NULL
       DO INSTEAD NOTHING
;
CREATE RULE rule_term_i2
       AS ON INSERT TO term
       WHERE (
       	     SELECT term_id FROM term
	     WHERE name        = new.name
	     AND   ontology_id = new.ontology_id
             AND   is_obsolete = new.is_obsolete
	     )
       	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_term_dbxref_i
       AS ON INSERT TO term_dbxref
       WHERE (
       	     SELECT dbxref_id FROM term_dbxref
	     WHERE dbxref_id = new.dbxref_id
	     AND   term_id   = new.term_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_term_path_i
       AS ON INSERT TO term_path
       WHERE (
       	     SELECT subject_term_id FROM term_path
	     WHERE subject_term_id   = new.subject_term_id
	     AND   predicate_term_id = new.predicate_term_id
	     AND   object_term_id    = new.object_term_id
	     AND   ontology_id	     = new.ontology_id
	     AND   distance	     = new.distance
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_term_relationship_i
       AS ON INSERT TO term_relationship
       WHERE (
       	     SELECT term_relationship_id FROM term_relationship
	     WHERE subject_term_id   = new.subject_term_id
	     AND   predicate_term_id = new.predicate_term_id
	     AND   object_term_id    = new.object_term_id
	     AND   ontology_id	     = new.ontology_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_term_relationship_term_i1
       AS ON INSERT TO term_relationship_term
       WHERE (
       	     SELECT term_relationship_id FROM term_relationship_term
	     WHERE term_relationship_id   = new.term_relationship_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_term_relationship_term_i2
       AS ON INSERT TO term_relationship_term
       WHERE (
       	     SELECT term_id FROM term_relationship_term
	     WHERE term_id   = new.term_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_term_synonym_i
       AS ON INSERT TO term_synonym
       WHERE (
       	     SELECT term_id FROM term_synonym
	     WHERE synonym = new.synonym
	     AND   term_id = new.term_id
	     )
	     IS NOT NULL
       DO INSTEAD NOTHING
;

--
-- Functions that may be used as an API by applications, e.g. load scripts etc.
-- 

-- this is used by load_ncbi_taxonomy.pl to speed up loading into the taxon
-- table by 1 to 2 orders of magnitude
CREATE OR REPLACE FUNCTION unconstrain_taxon ()
RETURNS INTEGER
AS
'
DROP RULE rule_taxon_i ON taxon;
SELECT 1;
'
LANGUAGE SQL
VOLATILE STRICT SECURITY DEFINER
;

-- this function re-establishes what unconstrain_taxon() removed temporarily
CREATE OR REPLACE FUNCTION constrain_taxon ()
RETURNS INTEGER
AS
'
CREATE RULE rule_taxon_i
       AS ON INSERT TO taxon
       WHERE (
             SELECT taxon_id FROM taxon 
             WHERE ncbi_taxon_id = new.ncbi_taxon_id
             )
       	     IS NOT NULL
       DO INSTEAD NOTHING
;
SELECT 1;
'
LANGUAGE SQL
VOLATILE STRICT SECURITY DEFINER
;
