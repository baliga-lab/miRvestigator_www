#################################################################
# @Program: miRvestigator                                       #
# @Version: 2                                                   #
# @Author: Chris Plaisier                                       #
# @Author: Christopher Bare                                     #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyright (C) 2010 by Institute for Systems Biology,          #
# Seattle, Washington, USA.  All rights reserved.               #
#                                                               #
# This source code is distributed under the GNU Lesser          #
# General Public License, the text of which is available at:    #
#   http://www.gnu.org/copyleft/lesser.html                     #
#################################################################

#
# Interact with the mirvestigator MySQL database
#

import MySQLdb
import datetime
import time
import re
import sys
from collections import defaultdict

# Database credentials should be changed from defaults after installation
MYSQL_USER = "mirv"
MYSQL_PASSWORD = "mirvestigator"

VERY_LARGE_FLOAT = 1e6

## Note that several of these methods use python's string formatting
## to build SQL strings, which is bad. This was done due to a problem
## with getting cursor.execute("insert into foo values (%s, %f, %d)", (a,b,c))
## to work with either decimal or floating point values.


def log(str):
    sys.stderr.write(str)
    sys.stderr.write("\n")
    sys.stderr.flush()


# a simplistic form of sanitizing input
# drop characters that are special to sql
_sanitize_re = re.compile(r"[\\\"\';]")
def _sanitize(s):
    return _sanitize_re.sub("", s)

def _get_db_connection():
    return MySQLdb.connect(host="localhost", user=MYSQL_USER, passwd=MYSQL_PASSWORD, db="mirvestigator")


def create_job_in_db(job):
    conn = _get_db_connection()
    try:
        created_at = datetime.datetime.now()
        cursor = conn.cursor()

        job_uuid = job['id']
        
        # store entry in jobs table
        cursor.execute("insert into jobs (uuid, created_at, updated_at, status) values ('%s', '%s', '%s', '%s');"
                       % (job_uuid, created_at.isoformat(), created_at.isoformat(), 'queued',))
        
        # store parameters
        for k, v in job.iteritems():
            if k!='id' and k!='genes':
                if k=='jobName':
                    v = _sanitize(v)
                cursor.execute("insert into parameters (job_uuid, name, value) values (%s, %s, %s);",
                               (job_uuid, k, str(v),) )

        #store genes
        if (job['genes']):
            for gene in job['genes']:
                gene = _sanitize(gene)[0:20]
                cursor.execute("insert into genes (job_uuid, name) values ('%s', '%s');" % (job_uuid, gene,) )

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: \n")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def delete_job(job_uuid):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        # delete parameters
        cursor.execute("delete from parameters where job_uuid=%s;", (job_uuid,))

        # delete genes
        cursor.execute("delete from genes where job_uuid=%s;", (job_uuid,))

        # for each motif, delete pssms, sites, and mirvestigator scores
        cursor.execute("select id from motifs where job_uuid=%s;", (job_uuid,))
        result_set = cursor.fetchall()
        for row in result_set:
            cursor.execute("delete from pssms where motif_id=%d", (int(row[0]),))
            cursor.execute("delete from sites where motif_id=%d", (int(row[0]),))
            cursor.execute("delete from mirvestigator_scores where motif_id=%d", (int(row[0]),))

        # delete motifs
        cursor.execute("delete from motifs where job_uuid=%s;", (job_uuid,))

        # delete entry in jobs table
        cursor.execute("delete from jobs where uuid=%s;", (job_uuid,))

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: \n")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


# find jobs submitted prior to the cutoff date
# returns a list of job uuids
def find_old_jobs(cutoff_datetime):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        # for each motif, delete pssms, sites, and mirvestigator scores
        cursor.execute("select uuid from jobs where created_at < %s;", (cutoff_datetime.isoformat(),))
        result_set = cursor.fetchall()
        uuids = []
        for row in result_set:
            uuids.append(row[0])

        return uuids
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: \n")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)

def get_unfinished_jobs():
    conn = _get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute("select j.uuid, name, value from jobs j join parameters p on j.uuid = p.job_uuid where j.status not in ('done', 'error')")
    result = defaultdict(dict)
    for row in cursor.fetchall():
        result[row[0]][row[1]] = row[2]
    cursor.close()
    for job_uuid in result:
        cursor = conn.cursor()
        cursor.execute("select name from genes where job_uuid = %s", [job_uuid])
        result[job_uuid]['genes'] = [row[0] for row in cursor.fetchall()]
        cursor.close()
    return result

def read_parameters(job_uuid):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        # job_name = user defined name
        # job_id = job_id
        # genes_submitted = number of genes submitted
        # annotated_sequences = number of genes with an annotated sequence
        # motif_sizes = an array of the set [6, 7, 8]
        # background_model = Default or 3' UTR
        # model_wobble = yes or no, if yes then tack on the G or U min freq

        cursor.execute("""
            select name, value
            from parameters
            where job_uuid=%s""", (job_uuid,))
        result_set = cursor.fetchall()

        # if there's nothing in the DB, return an empty dictionary
        if (len(result_set) == 0):
            return {}

        # build a dictionary of parameters
        parameters = {}
        for row in result_set:
            parameters[row[0]] = row[1]

        motif_sizes = []
        if ('m6' in parameters and parameters['m6'] == '6'):
            motif_sizes.append("6 bp")
        if ('m8' in parameters and parameters['m8'] == '8'):
            motif_sizes.append("8 bp")
        parameters['motif_sizes'] = motif_sizes

        seed_model = []
        if ('s6' in parameters and parameters['s6']):
            seed_model.append("6mer")
        if ('s7' in parameters and parameters['s7']):
            seed_model.append("7mer")
        if ('s8' in parameters and parameters['s8']):
            seed_model.append("8mer")
        parameters['seed_model'] = seed_model

        if ('wobble' in parameters and parameters['wobble']=='yes'):
            parameters['model_wobble'] = 'Yes (' + parameters['cut'] + ')'
        else:
            parameters['model_wobble'] = 'No'

        cursor.execute("""
            select name, sequence
            from genes
            where job_uuid=%s""", (job_uuid,))
        result_set = cursor.fetchall()
        genes = 0
        annotated_sequences = 0
        for row in result_set:
            genes += 1
            if (row[1]==1):
                annotated_sequences += 1

        parameters['genes_submitted'] = genes
        parameters['annotated_sequences'] = annotated_sequences

        return parameters

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def get_job_status(job_uuid):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute("select * from jobs where uuid=%s;", (job_uuid,))
        row = cursor.fetchone()
        result = {}
        if (row==None):
            result['created_at'] = '???'
            result['updated_at'] = '???'
            result['status'] = "not found"
            result['status_message'] = "no status found for " + str(job_uuid)
        else:
            result['created_at'] = row[1].strftime('%Y.%m.%d %H:%M:%S')
            result['updated_at'] = row[2].strftime('%Y.%m.%d %H:%M:%S')
            result['status'] = row[3]
            if (row[4] is not None):
                result['status_message'] = row[4]
        return result
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def update_job_status(job_uuid, status, message=None):
    conn = _get_db_connection()
    try:
        now = datetime.datetime.now()
        cursor = conn.cursor()
        if (message):
            cursor.execute("update jobs set status='%s', updated_at='%s', status_message='%s' where uuid='%s';"
                           % (status, now.isoformat(), message, job_uuid,))
        else:
            cursor.execute("update jobs set status='%s', updated_at='%s', status_message=null where uuid='%s';"
                           % (status, now.isoformat(), job_uuid,))
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def map_genes_to_entrez_ids(job_uuid, geneId, species):
    # update the genes table by joining with the gene_identifers table
    # note that this is really slow w/out an index on gene_identifiers
    # (species,id_type,identifier)

    if (geneId == "entrez"):
        # entrez_ids just map to themselves. No validation here,
        # cause we'll discover bad entrez ids when we try to map
        # them to sequences. If we want to validate here, we need
        # to make sure we have entries in gene_identifiers for all
        # genes for which we have sequences.
        sql = "update genes set entrez_id = name where job_uuid='%s';" % (job_uuid)
    else:
        # if we're given some other type of ID, try to map it back
        # to entrez ids
        sql = """update genes g join gene_identifiers gi on g.name = gi.identifier
                 set g.entrez_id = gi.entrez_id
                 where gi.species = '%s'
                 and gi.id_type = '%s'
                 and g.job_uuid = '%s';
              """ % (species, geneId, job_uuid)

    sql2 = "select distinct(entrez_id) from genes where job_uuid = '%s' and entrez_id is not null;" % (job_uuid)

    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute(sql)
        cursor.execute(sql2)
        result_set = cursor.fetchall()
        entrez_ids = []
        for row in result_set:
            entrez_ids.append(row[0])
        return entrez_ids
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)

def get_gene_mapping(job_uuid):
    
    sql = """select g.name, g.entrez_id, g.sequence, gi.identifier as symbol
             from genes g left join gene_identifiers gi on g.entrez_id=gi.entrez_id
             where g.job_uuid = '%s'
               and gi.id_type = 'symbol'
          """ % (job_uuid)

    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute(sql)
        result_set = cursor.fetchall()
        gene_mappings = {}
        for row in result_set:
            gene_mappings[row[1]] = {'name':row[0], 'sequence':row[2], 'symbol':row[3]}
        return gene_mappings
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def check_entrez_genes(species, genes):
  conn = _get_db_connection()
  num_found = 0
  for gene in genes:
    cursor = conn.cursor()
    cursor.execute('select count(*) from gene_identifiers where species=%s and entrez_id=%s',
                   (species, gene))
    num_found += cursor.fetchone()[0]
    cursor.close()
  conn.close()
  return num_found

def check_genes(id_type, species, genes):
    num_found = 0
    conn = _get_db_connection()
    cursor = conn.cursor()
    for gene in genes:
        cursor.execute('select count(*) from gene_identifiers where species=%s and id_type=%s and identifier=%s', (species, id_type, gene))
        num_found += cursor.fetchone()[0]
    cursor.close()
    conn.close()
    return num_found

# not used
def get_gene_dictionary(genes, geneId, species):
    gene_dictionary = {}
    gene_list = ",".join([ "'%s'" % (_sanitize(gene)) for gene in genes])
    sql = "select identifier, entrez_id from gene_identifiers where species='%s' and id_type='%s' and identifier in (%s)" % (species, geneId, gene_list)
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute(sql)
        while(True):
            row = cursor.fetchone()
            if row == None:
                break
            gene_dictionary[row[0]] = row[1]

        # add unmapped keys
        lc_keys = [key.lower() for key in gene_dictionary.keys()]
        for gene in genes:
            if gene.lower() not in lc_keys:
                gene_dictionary[gene] = None

        return gene_dictionary
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)

# not used
def store_genes(job_uuid, genes):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        #store genes
        if (genes):
            for gene in genes:
                gene = _sanitize(gene)[0:20]
                cursor.execute("insert into genes (job_uuid, name) values ('%s', '%s');" %
                               (job_uuid, gene,) )
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def set_genes_annotated(job_uuid, sequence_dict):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        #update genes to indicate whether sequence was found
        if (sequence_dict):
            for gene in sequence_dict.keys():
                gene = _sanitize(gene)[0:20]
                cursor.execute("update genes set sequence=true where job_uuid='%s' and entrez_id='%s';" %
                               (job_uuid, gene,) )
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def store_motif(job_uuid, pssm):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        # weeder can give a give a score of 'inf'. We handle that case by
        # the hack of substituting a very large float, which we hope will
        # always be greater than the next highest legitimate score
        if pssm.getEValue() == 'inf':
            weeder_score = VERY_LARGE_FLOAT
        else:
            weeder_score = float(pssm.getEValue())

        # we build the query as a string rather than using execute's string substitution because I
        # kept getting a TypeError: float argument required, not str:
        sql = """
            insert into motifs
            (job_uuid, name, score)
            values ('%s', '%s', %f);""" % (str(job_uuid), pssm.getName(), weeder_score,)

        log("store_motif sql=" + sql)

        cursor.execute(sql)
        motif_id = cursor.lastrowid

        # write pssm matrix
        i = 1
        for scores in pssm.getMatrix():
            sql = """
                insert into pssms
                (motif_id, position, a, c, g, t)
                values ('%s',%d,%f,%f,%f,%f);""" % (motif_id, i, float(scores[0]), float(scores[1]), float(scores[2]), float(scores[3]),)
            cursor.execute(sql)
            i += 1
                
        # motif_id int NOT NULL,
        # entrez_gene_id int,
        # sequence,
        # start,
        # quality,
        # mfe
                
        # sites is a dictionary w/ keys: gene, start, match, site, mfe
        i = 1
        sites = pssm.nsites
        for site in sites:
            cursor.execute("""
                insert into sites
                (motif_id, sort_order, entrez_gene_id, sequence, start, quality, mfe)
                values (%d, %d, '%s', '%s', %d, '%s', '%s')""" %
                (motif_id, i, str(site['gene']), site['site'], int(site['start']), site['match'], site['mfe'],))
            i += 1
        
        return motif_id

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)

# add order-by clauses

def read_motifs(job_uuid):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        # read motifs for this job
        cursor.execute("""
            select id, job_uuid, name, score
            from motifs
            where job_uuid = %s
            order by score asc;""",
            (str(job_uuid),))
        result_set = cursor.fetchall()
        motifs = []
        for row in result_set:
            motif = {}
            motif['motif_id'] = int(row[0])
            motif['job_uuid'] = row[1]
            motif['name'] = row[2]
            motif['score'] = float(row[3])
            motifs.append(motif)

        # read pssm matrix
        for motif in motifs:
            cursor.execute("""
                select a, c, g, t
                from pssms
                where motif_id=%d
                order by position;""" %
                (motif['motif_id'],))
            result_set = cursor.fetchall()
            matrix = []
            for row in result_set:
                matrix.append([row[0], row[1], row[2], row[3]])
            motif['matrix'] = matrix

        # read sites
        # each site is a dictionary w/ keys: gene, start, match, site
        for motif in motifs:
            cursor.execute("""
                select entrez_gene_id, sequence, start, quality, mfe
                from sites
                where motif_id=%d
                order by sort_order;""" %
                (motif['motif_id'],))
            result_set = cursor.fetchall()
            sites = []
            for row in result_set:
                site = {}
                site['gene']  = row[0]
                site['site']  = row[1]
                site['start'] = row[2]
                site['match'] = row[3]
                site['mfe']   = row[4]
                sites.append(site)
            motif['sites'] = sites

        return motifs

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


# given an motif id, returns a dictionary with keys motif_id, job_uuid, name, score
# name is the sequence of the motif
def read_motif(motif_id):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        # read motifs for this job
        cursor.execute("""
            select id, job_uuid, name, score
            from motifs
            where id = %s
            order by score asc;""",
            (str(motif_id),))
        row = cursor.fetchone()
        motif = {}
        motif['motif_id'] = int(row[0])
        motif['job_uuid'] = row[1]
        motif['name'] = row[2]
        motif['score'] = float(row[3])
        return motif

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def store_mirvestigator_scores(motif_id, scores):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()

        # motif_id int NOT NULL,
        # mirna_name varchar(100),              -- miRNA.name
        # mirna_seed varchar(8),                -- miRNA.seed
        # seedModel varchar(12),                -- model
        # alignment varchar(100),               -- statePath
        # viterbi_p float,                      -- vitPValue
        i = 1
        for score in scores:
            cursor.execute(
                """
                insert into mirvestigator_scores (motif_id, sort_order, mirna_name, mirna_seed, seedModel, alignment, viterbi_p)
                                          values (%d, %d, '%s', '%s', '%s', '%s', %f);
                """ %
                (motif_id, i, score['miRNA.name'], score['miRNA.seed'], score['model'], ";".join(score['statePath']), score['vitPValue'],))
            i += 1

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)


def read_mirvestigator_scores(motif_id):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute("""
            select motif_id, mirna_name, mirna_seed, seedModel, alignment, viterbi_p
            from mirvestigator_scores
            where motif_id=%d
            order by sort_order;
            """ % (int(motif_id),))
        result_set = cursor.fetchall()

        scores = []
        for row in result_set:
            score = {}
            score['miRNA.name'] = row[1]
            score['miRNA.seed'] = row[2]
            score['model'] = row[3]
            score['statePath'] = row[4].split(";")
            score['vitPValue'] = row[5]
            scores.append(score)

        return scores

    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)

def read_sites(motif_id):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute("""
            select entrez_gene_id, sequence, start, quality, mfe
            from sites
            where motif_id=%d
            order by sort_order;""" %
            (int(motif_id),))
        result_set = cursor.fetchall()
        # each site is a dictionary w/ keys: gene, start, match, site
        sites = []
        for row in result_set:
            site = {}
            site['gene']  = row[0]
            site['site']  = row[1]
            site['start'] = row[2]
            site['match'] = row[3]
            site['mfe']   = row[4]
            sites.append(site)
        return sites
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)

def get_species_by_mirbase_id(mirbase_id):
    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute("""
            select *
            from species
            where mirbase='%s';""" %
            (mirbase_id,))
        row = cursor.fetchone()
        species = {}
        if (row and len(row) >= 6):
            species['id']  = row[0]
            species['name']  = row[1]
            species['ucsc_name'] = row[2]
            species['ncbi_id'] = row[3]
            species['mirbase'] = row[4]
            species['weeder'] = row[5]
        return species
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)

# used in a dirty hack to get gene map for sites
def get_job_id_from_motif_id(motif_id):
    sql = """select job_uuid
             from motifs
             where id = '%s';
          """ % (motif_id)

    conn = _get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute(sql)
        row = cursor.fetchone()
        if (row == None):
            return None
        else:
            return row[0]
    finally:
        try:
            cursor.close()
        except Exception as exception:
            log("Exception closing cursor: ")
            log(exception)
        try:
            conn.close()
        except Exception as exception:
            log("Exception closing conection: ")
            log(exception)
