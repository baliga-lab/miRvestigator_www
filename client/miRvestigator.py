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

from mod_python import apache
from mod_python import util
import os
import re
import sys
import traceback
import datetime
import json
import Pyro.core
import mirv_csv
from Pyro.errors import ProtocolError
from mirv_db import get_job_status, read_parameters, read_motifs, read_mirvestigator_scores, get_gene_mapping, check_entrez_genes, check_genes
import admin_emailer
from ConfigParser import ConfigParser


# Libraries for plotting
import numpy, corebio                     # http://numpy.scipy.org and http://code.google.com/p/corebio/
from numpy import array, float64, log10   # http://numpy.scipy.org
from weblogolib import *                  # http://code.google.com/p/weblogo/

MAX_GENES = 1000
adminEmailer = admin_emailer.AdminEmailer()

# Plot a PSSM using weblogo
def plotPssmMatrix(pssmMatrix, fileName):
    dist = numpy.array( pssmMatrix, numpy.float64 ) 
    data = LogoData.from_counts(corebio.seq.unambiguous_rna_alphabet, dist*100)
    options = LogoOptions()
    options.color_scheme = colorscheme.nucleotide
    format = LogoFormat(data, options)
    fout = open(fileName, 'w')
    png_formatter(data, format, fout)
    fout.close()

# Reverse complement
def reverseComplement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

# Complement
def complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'U':'A'}
    complseq = [complement[base] for base in seq]
    return complseq

# Convert to RNA
def conv2rna(seq):
    conversion = {'A':'A', 'T':'U', 'C':'C', 'G':'G', 'N':'N', 'U':'U', '[':'[', ']':']'}
    rnaSeq = [conversion[base] for base in list(seq)]
    return ''.join(rnaSeq)

def alignSeed(alignment, seed, motif):
    alignment.pop() # Get rid of the extra state which is added by the forwardViterbi function
    start = 1    
    if alignment[0]=='NM1':
        for i in alignment:
            if i=='NM1':
                start += 1
    # Alignment
    motifAlign = '<font color=\'#CC0000\'><b>Motif</b></font><font color=\'#ffffff\'>_</font>5\'<font color=\'#ffffff\'>_</font>'
    seed = list(conv2rna(complement(seed)))
    seedAlign = '<font color=\'#ffffff\'>______</font>3\'<font color=\'#ffffff\'>_</font>'
    motif = list(conv2rna(motif))
    alignMe = alignment
    aligned = '<font color=\'#ffffff\'>_________</font>'
    lenMatch = 0
    # First get them zero'd to the same point
    if start>1:
        for i in range(start-1):
            seedAlign += seed.pop(0)
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += '-'
            alignMe.pop(0)
    if len(alignMe)>0 and not alignMe[0]=='PSSM0' and not alignMe[0]=='WOBBLE0':
        if alignMe[0][0]=='P':
            upTo = int(alignMe[0][4])
        elif alignMe[0][0]=='W':
            upTo = int(alignMe[0][6])
        for i in range(upTo):
            seedAlign += '<font color=\'#cccccc\'>-</font>'
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += motif.pop(0)
    # Then add the parts that align
    while 1:
        if len(alignMe)==0 or alignMe[0]=='NM2':
            break
        seedAlign += seed.pop(0)
        if alignMe[0][0]=='P':
            aligned += '<font color=\'#ff0000\'>|</font>'
        elif alignMe[0][0]=='W':
            aligned += '<font color=\'#0000ff\'>:</font>'
        lenMatch += 1
        motifAlign += motif.pop(0)
        alignMe.pop(0)
    # Then do the ending
    if len(alignMe)>0:
        for i in alignMe:
            seedAlign += seed[0]
            seed = seed[1:]
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += '<font color=\'#cccccc\'>-</font>'
        alignMe = []
    if len(motif)>0 and len(alignMe)==0:
        for i in motif:
            seedAlign += '<font color=\'#cccccc\'>-</font>'
            aligned += '<font color=\'#ffffff\'>_</font>'
            motifAlign += i
    motifAlign += '<font color=\'#ffffff\'>_</font>3\'<font color=\'#ffffff\'>___________</font>'
    aligned +=    '<font color=\'#ffffff\'>______________</font>'
    seedAlign +=  '<font color=\'#ffffff\'>_</font>5\'<font color=\'#ffffff\'>_</font><font color=\'#cc0000\'><b>miRNA Seed</b></font>'
    return [motifAlign, aligned, seedAlign, lenMatch]



def __get_config(req):
    config = ConfigParser()
    config.read(req.get_options()['MIRV_INI'])
    return config


# stuff parameters into a dictionary and pop those onto a queue
def submitJob(req):
    
    # create a job object which will be queued
    job = {}
    job['created'] = datetime.datetime.now()

    # get the gene list
    genes = req.form.getfirst('genes','').strip()
    if genes == '':
        return error_page("<b>Error</b>: no genes found.")
    job['genes'] = re.split('\s*[,;\s]\s*', genes)
    if (len(job['genes']) > MAX_GENES):
        return error_page("<b>Error</b>: miRvestigator can accept up to %d genes. Your request contained %d." % (MAX_GENES, len(job['genes']),))

    # Get the variables
    job['s6'] = str(req.form.getfirst('seedModel_6',''))
    job['s7'] = str(req.form.getfirst('seedModel_7',''))
    job['s8'] = str(req.form.getfirst('seedModel_8',''))
    job['bgModel'] = str(req.form.getfirst('bgModel',''))
    job['quality'] = str(req.form.getfirst('quality',''))
    job['species'] = str(req.form.getfirst('species',''))
    job['wobble'] = str(req.form.getfirst('wobble',''))
    job['cut'] = str(req.form.getfirst('cut',''))
    job['m6'] = str(req.form.getfirst('motif_6',''))
    job['m8'] = str(req.form.getfirst('motif_8',''))
    job['topRet'] = str(req.form.getfirst('topRet',''))
    job['jobName'] = str(req.form.getfirst('jobName',''))
    job['notify_mail'] = str(req.form.getfirst('notify_mail',None))
    job['geneId'] = str(req.form.getfirst('geneId', 'entrez'))
    job['viral'] = str(req.form.getfirst('viral',''))

    id_type = job['geneId']
    if id_type == 'entrez':
      num_found = check_entrez_genes(job['species'], job['genes'])
    else:
      num_found = check_genes(id_type, job['species'], job['genes'])

    if num_found == 0:
      req.log_error("request did not contain any valid genes")
      util.redirect(req, req.construct_url("/error"))
      return

    try:
        config = __get_config(req)
        # connect to miR server via Pyro
        with open(os.path.join(config.get('General', 'tmp_dir'), 'uri')) as uriFile:
          uri = uriFile.readline().strip()

        miR_server = Pyro.core.getProxyForURI(uri)
        Pyro.core.initClient()

        # submit job to server process and redirect to status page
        job_id = miR_server.submit_job(job)
    except ProtocolError as e:
        traceback.print_stack()
        traceback.print_exc()
        sys.stderr.flush()
        adminEmailer.warn("miRvestigator server is unreachable: \n\n" + str(e))
        util.redirect(req, req.construct_url("/error"))
        return
    except Exception as e:
        traceback.print_stack()
        traceback.print_exc()
        sys.stderr.flush()
        adminEmailer.warn("miRvestigator server error: \n\n" + str(sys.exc_info()[0]))
        util.redirect(req, req.construct_url("/error"))
        return

    util.redirect(req, req.construct_url("/status/%s/" % (job_id)))



# return a JSON string encoding job status
def status(req):
    id = str(req.form.getfirst('id',''))
    req.content_type='application/json'
    return json.dumps(get_job_status(id))

def parameters(req):
    id = str(req.form.getfirst('id',''))
    req.content_type='application/json'
    return json.dumps(read_parameters(id))

def sites_as_csv(req):
    motif_id = str(req.form.getfirst('motif_id',''))
    csv = mirv_csv.get_sites_as_csv(motif_id)
    req.content_type='text/csv'
    req.headers_out.add("Cache-Control", 'must-revalidate')
    req.headers_out.add("Pragma", 'must-revalidate')
    req.headers_out.add("Content-disposition", 'attachment; filename=mirvestigator_sites.csv')
    return csv

def scores_as_csv(req):
    motif_id = str(req.form.getfirst('motif_id',''))
    csv = mirv_csv.get_mirvestigator_scores_as_csv(motif_id)
    req.content_type='text/csv'
    req.headers_out.add("Cache-Control", 'must-revalidate')
    req.headers_out.add("Pragma", 'must-revalidate')
    req.headers_out.add("Content-disposition", 'attachment; filename=mirvestigator_scores.csv')
    return csv

def genes(req):
    id = str(req.form.getfirst('id',''))
    csv = mirv_csv.get_genes_as_csv(id)
    req.content_type='text/csv'
    req.headers_out.add("Cache-Control", 'must-revalidate')
    req.headers_out.add("Pragma", 'must-revalidate')
    req.headers_out.add("Content-disposition", 'attachment; filename=mirvestigator_sites.csv')
    return csv
    

# display results
def results(req):
    # need weederPSSMs1
    # need topRet
    # need mV
    # miRNADict[j.strip()] ??
    # functions: conv2rna, reverseComplement
    config = __get_config(req)
    id = str(req.form.getfirst('id',''))
    motifs = read_motifs(id)
    parameters = read_parameters(id)
    if (len(parameters) == 0):
        util.redirect(req, req.construct_url("/error/%s/" % (id)))
        return
    topRet = parameters['topRet']
    genesSubmitted = parameters['genes_submitted']
    annotatedSequences = parameters['annotated_sequences']
    mirbase_species = parameters['species']
    geneId = parameters.get('geneId', 'entrez')
    qualityThreshold = float(parameters['quality'])
    if 'viral' in parameters.keys() and parameters['viral']=='True':
        viral = True
    else:
        viral = False
    gene_mapping = get_gene_mapping(id)

    # read mirbase miRNAs so we can link back to mirbase
    viralSpecies = []
    if viral:
        viral_species_filename = config.get('General', 'viral_species_filename')
        with open(viral_species_filename,'r') as inFile:
            for line in inFile.readlines():
                splitUp = line.strip().split('\t')
                if splitUp[1]=='VRL':
                    viralSpecies.append(splitUp[0])

    import gzip
    mirna_filename = config.get('General', 'mirna_filename')
    with gzip.open(mirna_filename,'r') as miRNAFile:
        miRNADict = {}
        while 1:
            miRNALine = miRNAFile.readline()
            seqLine = miRNAFile.readline()
            if not miRNALine:
                break
            # Get the miRNA name
            miRNAData = miRNALine.lstrip('>').split(' ')
            curMiRNA = miRNAData[0]
            if (curMiRNA.split('-'))[0]==mirbase_species:
                miRNADict[curMiRNA] = miRNAData[1]
            if viral==True and ((curMiRNA.split('-'))[0] in viralSpecies):
                miRNADict[curMiRNA] = miRNAData[1]

    s = '<html>\n'
    s += '<head>'
    s += '<title>miRvestigator Framework</title>\n'
    s += '<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.4.4/jquery.min.js"></script>\n'
    s += '<script src="/mirvestigator.js"></script>\n'
    s += '<script type="text/javascript">\n  var _gaq = _gaq || [];\n  _gaq.push([\'_setAccount\', \'UA-19292534-1\']);\n  _gaq.push([\'_trackPageview\']);\n  (function() {\n    var ga = document.createElement(\'script\');\n    ga.type = \'text/javascript\';\n    ga.async = true;\n    ga.src = (\'https:\' == document.location.protocol ? \'https://ssl\' : \'http://www\') + \'.google-analytics.com/ga.js\';\n    var s = document.getElementsByTagName(\'script\')[0];\n    s.parentNode.insertBefore(ga, s);\n  })();\n'
    s += 'function toggleVisible(id, link_elem) {\n'
    s += '  if (document.getElementById) {\n'
    s += '    obj = document.getElementById(id);\n'
    s += '    if (obj) {\n'
    s += '      if (obj.style.display == \'none\') {\n'
    s += '        obj.style.display = \'\'\n'
    s += '        if (link_elem) {\n'
    s += '          link_elem.innerHTML = \'[-]\';\n'
    s += '        }\n'
    s += '      } else {\n'
    s += '        obj.style.display = \'none\'\n'
    s += '        if (link_elem) {\n'
    s += '          link_elem.innerHTML = \'[+]\';\n'
    s += '        }\n'
    s += '      }\n'
    s += '    }\n'
    s += '  }\n'
    s += '}\n'

    s += 'function toggleVisibleSubParams() {\n'
    s += '  subParams = $("#subParams");\n'
    s += '  if ( subParams.css("display") == "none" ) {\n'
    s += '    subParams.css("display", "");\n'
    s += '    $("#subParams_toggle").text("[-]");\n'
    s += '  } else {\n'
    s += '    subParams.css("display", "none");\n'
    s += '    $("#subParams_toggle").text("[+]");\n'
    s += '  }\n'
    s += '}\n'
    s += '  window.onload=function() {\n'
    s += '    var url = window.location.pathname;\n'
    s += '    var m = url.match(/\/results\/([0-9a-fA-F\-]+)\/(.?)/);\n'
    s += '    if (m && m.length > 1) {\n'
    s += '      job_id = m[1];\n'
    s += '        if (job_id) {\n'
    s += '          getParameters(job_id);\n'
    s += '        }\n'
    s += '    }\n'
    s += '  }\n'
    s += '</script>\n<style>\n  body { font-family: arial, sans-serif; }\n</style>\n</head>\n'
    s += '<body bgcolor=\'#333333\' link=\'cc0000\' vlink=\'cc0000\'>\n'
    s += '<center>\n'
    s += '<table><tr><td>\n'
    s += '<table cellpadding=\'15%\' cellspacing=\'5\' bgcolor=\'#999966\'>\n'
    s += '<tr>\n'
    s += '<td align=\'center\' valign=\'center\' bgcolor=\'#333333\' width=\'80\'><b><a style=\'color: rgb(255,0,0); text-decoration:none\' href=\'/\'>H</a><a style=\'color: rgb(204,204,0); text-decoration:none\' href=\'/\'>ome</a></b>\n'
    s += '</td>\n'
    s += '<td align=\'center\' valign=\'center\' bgcolor=\'#999966\' width=\'80\'><b>Results</b>\n'
    s += '</td>\n'
    s += '<td align=\'center\' valign=\'center\' bgcolor=\'#333333\' width=\'80\'><b><a style=\'color: rgb(255,0,0); text-decoration:none\' href=\'/help\'>H</a><a style=\'color: rgb(204,204,0); text-decoration:none\' href=\'/help\'>elp</a></b>\n'
    s += '</td>\n'
    s += '<td align=\'center\' valign=\'center\' bgcolor=\'#333333\' width=\'80\'><b><a style=\'color: rgb(255,0,0); text-decoration:none\' href=\'/tutorial\'>T</a><a style=\'color: rgb(204,204,0); text-decoration:none\' href=\'/tutorial\'>utorial</a></b>\n'
    s += '</td>\n'
    s += '<td align="center" valign="center" bgcolor="#333333" width="0"><b>'
    s += '<a style="color: rgb(255,0,0); text-decoration:none" href="/citation">C</a><a style="color: rgb(204,204,0); text-decoration:none" href="/citation">itation</a></b></td>\n'
    s += '</tr>\n'
    s += '</table>\n'
    s += '<table bgcolor=\'#999966\' cellpadding=\'10%\'>\n'
    s += '<tr>\n'
    s += '<td>'
    s += '<table width=\'100%\' cellpadding=\'15%\'>\n'
    s += '<tr>\n'
    s += '<td align=\'center\' valign=\'center\' bgcolor=\'#333333\'><font size=6><b><font color=\'#ff0000\'>miR</font><font color=\'#cccc00\'>vestigator Framework Results</font></b></font>\n'
    s += '</td>\n'
    s += '</tr>\n'
    s += '</table>\n'
    s +='<p><table width=\'100%\' cellpadding=\'10%\'>\n'
    s += '<tr><td bgcolor=\'#000000\'><center><font color=\'#cccc00\' size=4><b><a href="#subParams" onclick="javascript: toggleVisibleSubParams(); return false;" style=\'color: rgb(204,204,0); text-decoration: none\'>Submission Parameters <font id="subParams_toggle" color=\'#ff0000\'>[+]</font></a></b></font></center></td></tr>\n'
    s += '<tr id="subParams" style="display: none;"><td bgcolor=\'#cccccc\'>\n'
    s += '<table width=\'100%\' cellpadding=\'5%\'>\n'
    s += '<tr><td colspan=2 bgcolor=\'#333333\' align=\'center\' valign=\'center\'><b><font color=\'#cccc00\' size=4>Basic Parameters</font></b></td></tr>\n'
    s += '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Job Name:</b></font></td><td bgcolor=\'#666666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="job_name"></span></center></td></tr></table></td></tr>\n'
    s+= '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Job ID:</b></font></td><td bgcolor=\'#666666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="job_id"></span></center></td></tr></table></td></tr>\n'
    s += '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Genes Submitted:</b></font> <span style="color:#ff0000"><b>(<a href="/genes/' + id + '/">gene mapping as CSV</a>)</b></span>'
    s += '</td><td bgcolor=\'#666666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="genes_submitted"></span></center></td></tr></table></td></tr>\n'
    s += '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Genes Annotated with Sequences:</b></font></td><td bgcolor=\'#66666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="annotated_sequences"></span></center></td></tr></table></td></tr>\n'
    s += '<tr><td colspan=2></td></tr>\n'
    s += '<tr><td colspan=2 bgcolor=\'#333333\' align=\'center\' valign=\'center\'><b><font color=\'#cccc00\' size=4>Weeder Parameters</font></b></td></tr>\n'
    s += '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Motif Sizes:</b></font></td><td bgcolor=\'#666666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="motif_sizes"></span></center></td></tr></table></td></tr>\n'
    s += '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Background Model:</b></font></td><td bgcolor=\'#666666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="background_model"></span></center></td></tr></table></td></tr>\n'
    s += '<tr><td colspan=2></td></tr>\n'
    s += '<tr><td colspan=2 bgcolor=\'#333333\' align=\'center\' valign=\'center\'><b><font color=\'#cccc00\' size=4>miRvestigator HMM Parameters</font></b></td></tr>\n'
    s+= '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Seed Models:</b></font></td><td bgcolor=\'#666666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="seed_model"></span></center></td></tr></table></td></tr>\n'
    s += '<tr><td bgcolor=\'#333333\' align=\'center\' valign=\'center\'><font color=\'#ffffff\'><b>Model Wobble Base-Pairing:</b></font></td><td bgcolor=\'#666666\' align=\'center\' valign=\'center\'><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#ffffff\'><center><span id="model_wobble"></span></center></td></tr></table></td></tr>\n'
    s += '</table>\n'
    s += '</td></tr>\n'
    s += '<tr><td></td></tr>\n'
    s += '<tr><td bgcolor=\'#000000\'><center><font color=\'#cccc00\' size=4><b>Summary of Results</b></font></center></td></tr>\n'
    # Make a table that looks like this:
    # motif consensus | top miRNA | seed | length | vitPValue | # seqs with site (% of seqs with site)
    if (len(motifs)==0):
        s+= '<div style="text-align:center;background:#333333;color:#cc0000;font-weight:bold;padding-top:6ex;padding-bottom:6ex;">No motifs found.</div>'

    s += '<table width=\'100%\' cellpadding=\'15%\'><tr><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Motif</font></b></center></td><td bgcolor=\'#333333\'><\
center><b><font color=\'#ffffff\'>Top miRNA</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Complementary Base-Pairing</font></b></center></td><td bgcolor=\'#3333\
33\'><center><b><font color=\'#ffffff\'>Complementarity</br>P-Value</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>% of Input</br>Sequences with Site</font>\
</b></center></td></tr>\n'

    pssm_images_dir = config.get('General', 'pssm_images_dir')
    for motif in motifs:
        scoreList = read_mirvestigator_scores(motif['motif_id'])
        if topRet=='all':
            topRet = len(scoreList)
        else:
            topRet = int(topRet)

        top_score_i = 1
        top_score_viterbi_p = float(scoreList[0]['vitPValue'])
        row_count = 1
        while float(scoreList[top_score_i]['vitPValue']) == top_score_viterbi_p:
            row_count +=1
            top_score_i += 1

        top_score_i = 0
        i = scoreList[top_score_i]
        align1 = alignSeed(i['statePath'], i['miRNA.seed'], motif['name'])
        png_path = os.path.join(pssm_images_dir, '%s_%s.png' % (id, motif['name']))
        if not os.path.exists(png_path):
            plotPssmMatrix(motif['matrix'], png_path)

        s += '<tr>'
        s += '<td bgcolor=\'#ffffff\' rowspan="' + str(row_count) + '"><center><a href=\'#'+motif['name']+'_miRNAs\'><img src=\'/images\
/pssms/'+id+'_'+motif['name']+'.png\' alt=\''+conv2rna(str(motif['name']))+'\'></center></td>'
        s += '<td bgcolor=\'#ffffff\'><center>'
        miRNAs = []
        for miRNA_name in i['miRNA.name'].split('_'):
            miRNA_name = miRNA_name.strip()
            if miRNADict.has_key(miRNA_name):
                miRNAs.append('<a href=\'http://mirbase.org/cgi-bin/mature.pl?mature_acc='+str(miRNADict[miRNA_name])+'\' target=\'_blank\'>'+miRNA_name+'</a>')
            else:
                miRNAs.append(miRNA_name)
        s += '</br>'.join(miRNAs)
        s += '</center></td>\n'

        s += '<td bgcolor=\'#ffffff\'><center><pre>'+str(align1[0])+'\n'+str(align1[1])+'\n'+str(align1[2])+'</pre></center></td>\n'
        s += '<td bgcolor=\'#ffffff\' rowspan="' + str(row_count) + '"><center>'+str('%.1e' % float(i['vitPValue']))+'</center></td>\n'
        # Get number or sequences with a hit                                                                                                                                         
        genesWithSite = []
        for site in motif['sites']:
            if float(site['match']) >= qualityThreshold and not site['gene'] in genesWithSite:
                genesWithSite.append(site['gene'])
        if (annotatedSequences==0):
            percGenes = "?error?"
        else:
            percGenes = (float(len(genesWithSite))/float(annotatedSequences))*float(100)
            if percGenes==float(100):
                percGenes = str(100)
            else:
                percGenes = str('%.2g' % percGenes)
        s += '<td bgcolor=\'#ffffff\' rowspan="' + str(row_count) + '"><center><a href=\'#'+motif['name']+'_sites\'>'+percGenes+'%</a></center></td>\n'
        s += '</tr>\n'
        top_score_i += 1

        while float(scoreList[top_score_i]['vitPValue']) == top_score_viterbi_p:
            i = scoreList[top_score_i]
            align1 = alignSeed(i['statePath'], i['miRNA.seed'], motif['name'])
            s += '<tr>'
            s += '<td bgcolor=\'#ffffff\'><center>'
            miRNAs = []
            for miRNA_name in i['miRNA.name'].split('_'):
                miRNA_name = miRNA_name.strip()
                if miRNADict.has_key(miRNA_name):
                    miRNAs.append('<a href=\'http://mirbase.org/cgi-bin/mature.pl?mature_acc='+str(miRNADict[miRNA_name])+'\' target=\'_blank\'>'+miRNA_name+'</a>')
                else:
                    miRNAs.append(miRNA_name)
            s += '</br>'.join(miRNAs)
            s += '</center></td>\n'
                    
            s += '<td bgcolor=\'#ffffff\'><center><pre>'+str(align1[0])+'\n'+str(align1[1])+'\n'+str(align1[2])+'</pre></center></td>\n'
            s += '</tr>\n'
            top_score_i += 1

    s += '</table></p>\n'

    for motif in motifs:
        s += '<p  id=\''+motif['name']+'_miRNAs\'><table width=\'100%\'cellpadding=\'15%\'>\n<tr>\n<td align=\'center\' valign=\'center\' bgcolor=\'#000000\'>\n<font size=4><b><font color=\'#cccc00\'>'
        if not topRet=='all':
            s += 'Top <font color=\'#ff0000\' size=4>'+str(topRet)+'</font>'
        elif topRet=='all':
            s += '<font color=\'#ff0000\' size=4>All</font>'
        s += ' miRNAs Complementary to the Weeder Motif</font> <font color=\'#ff0000\'>'+conv2rna(motif['name'])+'</font> &nbsp;'
        s += '<a href="#results" onclick="toggleVisible(\'results\', this); return false;" style=\'color: #ff0000; text-decoration: none\'>[+]</a></b></font></td></tr>\n'
        s += '<tr id="results" style="display: none;" width=600>\n'
        s += '<td bgcolor="#333333">\n'
        s += '<font color="#ffffff"><b>What do the columns mean?</b> <p><ul><li><b>miRNA Name</b> = The name of the name(s) for the unique seed sequence. There may be more than one miRNA annotated for a unique seed sequence because they vary in the 3\' terminus of the mature miRNA. Each miRNA is a link to it\'s entry on <a href="http://www.mirbase.org" style="color: rgb(204,204,0)" target="_blank">miRBase</a></li>&nbsp;</br> <li><b>miRNA Seed</b> = The sequence for seed that is the most complementary to the over-represented motif. The seed will be as long as the seed model described in the next column.</li>&nbsp;</br> <li><b>Seed Model</b> = Base-pairing models for the seed regions of a miRNA to the 3\' UTR of target transcripts. The 8mer, 7mer-m8, and 7mer-a1 models are the canonical models of miRNA to mRNA base-pairing. The 6mer models are considered marginal models as they typically have a reduced efficacy and are more likely to occur by chance alone. By default all of the seed models are used. The seed models are described in this figure:</br>&nbsp;</br><center><img src="/images/seedModels.gif" width=400></center></li>&nbsp;</br><li><b>Length of Complementarity</b> = The length of complementarity (or wobble if enabled) base-pairing between the sequence motif and the miRNA seed sequence.</li>&nbsp;</br> <li><b>Complementary Base-Pairing</b> = The complementarity of the over-represented sequence motif on top 5\'&rArr;3\' to the miRNA seed sequence given the seed model 3\'&rArr;5\'. (<b>Note:</b> <span style=\'background-color: rgb(255,255,255);\'>&nbsp;<b><font color="#ff0000">|</font></b>&nbsp;</span> = a complementary, <span style=\'background-color: rgb(255,255,255);\'>&nbsp;<b><font color="#0000ff">:</font></b>&nbsp;</span> = a wobble, <span style=\'background-color: rgb(255,255,255);\'><font color="#000000">"</font> <font color="#000000">"</font></span> (space) = not complementary, and for the sequences <span style=\'background-color: rgb(255,255,255);\'>&nbsp;<b><font color="#cccccc">-</font></b>&nbsp;</span> = a gapping at the start or end.)</li>&nbsp;</br> <li><b>Complementarity P-Value</b> = Significance of complementarity between the over-represented sequence motif and the miRNA seed sequence. (<b>Note:</b> A perfectly complementary 8mer seed model has a Complementarity p-value of 1.5e-05, for a 7mer seed model 6.1e-05, and for a 6mer seed model 2.4e-04.)</li></ul> <b>What is considered good complementarity?</b> <p>This will depend upon your data and what downstream analysis you plan to do with it. But a good rule of thumb is that if you find perfect complementarity for a 7mer or 8mer (Complementary P-Value = <font color="#ff0000"><b>6.1e-05</b></font> and <font color="#ff0000"><b>1.5e-05</b></font>; respectively) this is likely to be of interest. Follow up with experimental studies will help to determine the false discovery rate for your dataset.</p></font>\n'
        s += '</td>\n'
        s += '</tr>\n'
        s += '</table>\n'
        scoreList = read_mirvestigator_scores(motif['motif_id'])
        if topRet=='all':
            topRet = len(scoreList)
        else:
            topRet = int(topRet)
        s += '<table width="100%%" cellpadding="5%%"><tr><td style="text-align:center; font-size:10pt; background:#cccccc;"><a href="/scores/csv/%s">Download table as CSV</a></td></tr></table>\n' % (motif['motif_id'],)
        s += '<table width=\'100%\' cellpadding=\'15%\'><tr><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>miRNA Name</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>miRNA Seed</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Seed Model</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Length of</br>Complementarity</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Complementary Base-Pairing</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Complementarity</br>P-Value</font></b></center></td></tr>'
        for k in range(topRet):
            i = scoreList[k]
            align1 = alignSeed(i['statePath'], i['miRNA.seed'], motif['name'])
            s += '<tr>'
            s += '<td bgcolor=\'#ffffff\'><center>'
            
            # generate links to mirbase miRNAs
            miRNAs = []
            for miRNA_name in i['miRNA.name'].split('_'):
                miRNA_name = miRNA_name.strip()
                if miRNADict.has_key(miRNA_name):
                    miRNAs.append('<a href=\'http://mirbase.org/cgi-bin/mature.pl?mature_acc='+str(miRNADict[miRNA_name])+'\' target=\'_blank\'>'+miRNA_name+'</a>')
                else:
                    miRNAs.append(miRNA_name)
            
            s += '</br>'.join(miRNAs)
            s += '</center></td>'

            s += '<td bgcolor=\'#ffffff\'><center>'+conv2rna(reverseComplement(str(i['miRNA.seed'])))+'</center></td>'
            s += '<td bgcolor=\'#ffffff\'><center>'+str(i['model'])+'</center></td><td bgcolor=\'#ffffff\'><center>'+str(align1[3])+'</center></td>'
            s += '<td bgcolor=\'#ffffff\'><center><pre>'+str(align1[0])+'\n'+str(align1[1])+'\n'+str(align1[2])+'</pre></center></td>'
            s += '<td bgcolor=\'#ffffff\'><center>'+str('%.1e' % float(i['vitPValue']))+'</center></td>'
            s += '</tr>'
        s += '</table></p>'
        #'gene':seqDict[splitUp[0]], 'strand':splitUp[1], 'site':splitUp[2], 'start':splitUp[3], 'match':splitUp[4].lstrip('(').rstrip(')')
        # pssm1.nsites
        s += '<p id=\''+motif['name']+'_sites\'><table width=\'100%\' bgcolor=\'#000000\' cellpadding=\'15%\'><tr><td align=\'center\' valign=\'center\'><font size=4><b><font color=\'#cccc00\'>Position of Putative miRNA Binding Sites in Submitted Genes</br>for the Weeder Motif</font> <font color=\'#ff0000\'>'+conv2rna(motif['name'])+'</font> &nbsp;'
        s += '<a href="#sites" onclick="toggleVisible(\'sites\', this); return false;" style=\'color: #ff0000; text-decoration: none\'>[+]</a></b></font></td></tr>\n'
        s += '<tr id="sites" style="display: none;" width=600><td bgcolor="#333333"><font color="#ffffff">\n'
        s += '<b>Where do these sites come from?</b>\n<p>As part of the miRvestigator framework <a href="http://159.149.109.9/modtools/" style="color: rgb(204,204,0)" target="_blank">Weeder</a> provides predicted miRNA binging sites in the 3\' untranslated regions (UTRs) of the analyzed genes. Predicted binding sites were split into three different similarity bins:  <font color="#ff0000">High quality</font> - &#8805; 95% similarity to the miRNA seed sequence (red), <font color="#cccc00">Medium quality</font> 95% &#8805; similarity &#8805; 90% to the miRNA seed sequence (yellow), and <font color="#00ff00">Fair quality</font> 90% &#8805; similarity &#8805; 85% to the miRNA seed sequence (green). These sites can be used to develop follow-up experiments such as luciferase reporter assays to validate the efficacy of these sites.</p>\n<b>What do the columns mean?</b>\n<p><ul><li><b>Gene</b> = gene identifier from the input set that links to the mapped Entrez identifier page on NCBI.</li></br>&nbsp;</br>\n<li><b>Gene Symbol</b> = Official gene symbol if available.</li></br>&nbsp;</br>\n<li><b>Site</b> = The sequence for site identified by Weeder. If it is in square brackets indicates that the site is of lower similarity.</li></br>&nbsp;</br>\n<li><b>Start Relative to Stop Codon</b> = The 3\' UTR begins following the stop codon (which is set at 0 base-pairs (bp)). Thus the values in this column describe the start of the site in bp after the stop codon.</li></br>&nbsp;</br>\n<li><b>% Similarity to Consensus Motif</b> = The similarity of the predicted site to the consensus motif is computed as a percentage. Predicted binding sites were split into three different similarity to consensus bins:  <font color="#ff0000">High quality</font> - &#8805; 95% similarity to the miRNA seed sequence (red), <font color="#cccc00">Medium quality</font> 95% &#8805; similarity &#8805; 90% to the miRNA seed sequence (yellow), and <font color="#00ff00">Fair quality</font> 90% &#8805; similarity &#8805; 85% to the miRNA seed sequence (green).</li></br>&nbsp;</br>\n<li><b>Minimum Free Energy (MFE) of mRNA-miRNA Duplex</b> = the minimum free energy (MFE) of duplexing for the reverse complement of the motif consensus sequence and putative target site sequences using the ViennaRNA package.</li></ul></p>'
        s += '</font></td></td></table>'
        s += '<table width="100%%" cellpadding="5%%"><tr><td style="text-align:center; font-size:10pt; background:#cccccc;"><a href="/sites/csv/%s">Download table as CSV</a></td></tr></table>\n' % (motif['motif_id'],)
        s += '<table width=\'100%\' cellpadding=\'15%\'><tr><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Gene</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Gene symbol</font><\
/b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Sequence of Site</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Start Relative to</br>Stop Codon (bp)</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>% Similarity to Consensus Motif</br>(Quality = </font><font color=\'#cc0000\'>High</font><font color=\'#ffffff\'> | </font><font color=\'#cccc00\'>Medium</font><font color=\'#ffffff\'> | </font><font color=\'#00cc00\'>Fair</font><font color=\'#ffffff\'>)</font></b></center></td><td bgcolor=\'#333333\'><center><b><font color=\'#ffffff\'>Minimum Free Energy (MFE)</br>of mRNA-miRNA Duplex</font></b></center></td></tr>'
        for i in motif['sites']:
            if float(i['match']) >= qualityThreshold:
                col1 = '#000000'
                if float(i['match']) >= float(95):
                    col1 = '#cc0000'
                elif float(i['match']) >= float(90):
                    col1 = '#cccc00'
                elif float(i['match']) >= float(85):
                    col1 = '#00cc00'
                gene_name_map = gene_mapping.get(i['gene'], None)
                if (gene_name_map==None or geneId=='symbol' or geneId=='entrez'):
                    gene_name = i['gene']
                else:
                    gene_name = gene_name_map['name']
                if (gene_name_map==None or gene_name_map['symbol'] == None):
                    gene_symbol = ''
                else:
                    gene_symbol = ' ' + gene_name_map['symbol']
                s += '<tr><td bgcolor=\'#ffffff\'><center>'+str('<a href=\'http://www.ncbi.nlm.nih.gov/gene/'+str(i['gene'])+'\' target=\'_blank\'>'+gene_name+'</a>')+'<td bgcolor=\'#ffffff\'><center>'+gene_symbol+'</center></td>'+'</center></td><td bgcolor=\'#ffffff\'><center>'+conv2rna(str(i['site']))+'</center></td><td bgcolor=\'#ffffff\'><center>'+str(i['start'])+'</center></td><td bgcolor=\'#ffffff\'><font color=\''+str(col1)+'\'><center><b>'+i['match']+'</b></center></font></td><td bgcolor=\'#ffffff\'><center>'+str(i['mfe'])+'</center></td></tr>'
        s += '</table></p>'
    s += '<p><table width=\'100%\' cellpadding=\'5%\'><tr><td bgcolor=\'#c0c0c0\'><center>Need help? Please contact <font color=\'#0000ff\'>cplaisier(at)systemsbiology.org</font> or <font color=\'#0000ff\'>cbare(at)systemsbiology.org</font> if you have any questions, comments or concerns.<br>Developed at the <a href=\'http://www.systemsbiology.org\' target=\'_blank\' style=\'color: rgb(0,0,255)\'>Institute for Systems Biology</a> in the <a href=\'http://baliga.systemsbiology.net/\' target=\'_blank\' style=\'color: rgb(0,0,255)\'>Baliga Lab</a>.</center></td></tr></table></p>'
    s += '</center></td></tr></table></td></tr></table></center></body></html>'
    return s

def error_page(msg):
    s = """<html>
    <head>
      <title>miRvestigator Framework</title>
    </head>

    <body bgcolor='#333333' link='#ffcc00' vlink='#ffcc00'>
    <font face='arial'><center>
    <table width=620 bgcolor='#999966' cellpadding='10%'><tr><td><center>
    <table width=600 bgcolor='#333333' cellpadding='15%'><tr><td align='center' valign='center'><font size=6><b><font color='#ff0000'>miR</font><font color='#cccc00'>vestigator Framework</font></b></font></td></tr></table>

    <table cellpadding='5%' cellspacing=3 width='100%'>
    <tr><td bgcolor='#FFFFFF' style="height:60ex; padding-left:4em; padding-right:4em;"><center>
      <p>""" + msg + """</p>
      <p><a href="/">Return to miRvestigator home page</a></p>
    </center></td></tr>
    </table>
    </center></td></tr></table>
    </body>
</html>"""
    return s

