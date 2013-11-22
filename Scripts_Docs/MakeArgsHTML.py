"""
Automatically generates html of presto command line usage
"""

import subprocess
from collections import OrderedDict

presto_path = '/home/jason/workspace/igpipeline/presto_v0.4'
tools = OrderedDict([('AlignSets', ['muscle', 'offset', 'table']),
                     ('AssemblePairs', ['align', 'join']),
                     ('BuildConsensus', []),
                     ('CollapseSeq', []),
                     ('EstimateError', []),
                     ('FilterSeq', ['length','missing','repeats','quality','maskqual', 'trimqual']),
                     ('MaskPrimers', ['align', 'score']),
                     ('PairSeq', []),
                     ('ParseHeaders', ['add', 'collapse','delete','expand', 'rename', 'table', 'convert']),
                     ('ParseLog', []),
                     ('SplitSeq', ['convert','count', 'group', 'sample', 'samplepair', 'sort'])])

with open('arguments.html', 'w') as doc_handle, \
     open('navigation.html', 'w') as nav_handle: 
    # Start navigation list
    nav_handle.write('<ul>\n')
    for k, v in tools.iteritems():
        nav_handle.write('<li><a href="documentation.php#%s">%s</a></li>\n' % (k, k))
        doc_handle.write('<div id=%s></div><h2>%s</h2>\n' % (k, k))
        cmd = '/'.join([presto_path, k + '.py'])
        main_msg = subprocess.check_output([cmd, '--help'])
        doc_handle.write('<pre>\n%s\n</pre>\n' % main_msg)
        for c in v:
            doc_handle.write('<h3>%s</h3>\n' % c)
            sub_msg = subprocess.check_output([cmd, c, '--help'])
            doc_handle.write('<pre>\n%s</pre>\n' % sub_msg)
    # End navigation list
    nav_handle.write('</ul>\n')

print 'Wrote %s' % nav_handle.name
print 'Wrote %s' % doc_handle.name