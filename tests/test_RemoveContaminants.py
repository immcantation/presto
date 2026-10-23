"""
Unit tests for RemoveContaminants
"""
# Info
__author__ = 'Gisela Gabernet'

# Imports
import multiprocessing as mp
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from collections import OrderedDict
from Bio import SeqIO

# Presto imports
from presto.Annotation import parseAnnotation, parseLog
from presto.Defaults import default_out_args

# Paths
test_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(test_path, 'data')

# Import script
sys.path.append(os.path.join(test_path, os.pardir, 'bin'))
import RemoveContaminants


class TestRemoveContaminants(unittest.TestCase):
    def setUp(self):
        print('-> %s()' % self._testMethodName)
        self.labels = ['A', 'B', 'C']
        self.files_1 = [os.path.join(data_path, 'contam_sample%s_R1.fastq' % s) for s in self.labels]
        self.files_2 = [os.path.join(data_path, 'contam_sample%s_R2.fastq' % s) for s in self.labels]
        self.out_dir = tempfile.mkdtemp()
        self.start = time.time()

    def tearDown(self):
        shutil.rmtree(self.out_dir)
        t = time.time() - self.start
        print('<- %s() %.3f' % (self._testMethodName, t))

    def _write_fastq(self, name, records):
        path = os.path.join(self.out_dir, name)
        with open(path, 'w') as h:
            for title, seq in records:
                h.write('@%s\n%s\n+\n%s\n' % (title, seq, 'I' * len(seq)))
        return path

    def _run(self, **kwargs):
        out_dir = tempfile.mkdtemp(dir=self.out_dir)
        out_args = dict(default_out_args)
        out_args.update({'out_dir': out_dir, 'failed': True, 'gzip_output': False})
        out_args['log_file'] = os.path.join(out_dir, 'log.txt')
        RemoveContaminants.removeContaminants(self.files_1, self.files_2, labels=self.labels,
                                          out_args=out_args, **kwargs)
        result = {}
        for s in self.labels:
            for r in ('R1', 'R2'):
                for label in ('pass', 'fail'):
                    f = os.path.join(out_dir, 'contam_sample%s_%s_contam-%s.fastq' % (s, r, label))
                    result[(s, r, label)] = OrderedDict(
                        (parseAnnotation(rec.id)['ID'], parseAnnotation(rec.id))
                        for rec in SeqIO.parse(f, 'fastq'))
        with open(out_args['log_file']) as h:
            log = [parseLog(r) for r in h.read().split('\n\n') if r.strip()]
        return result, log

    def test_matchSeq(self):
        match = RemoveContaminants.matchSeq
        self.assertTrue(match('ACGTACGTAC', 'ACGTACGTAC'))
        # Ns are ignored, but different lengths never match
        self.assertTrue(match('ACGTACGTAC', 'ACGNACGTAC'))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTACGGGG'))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTA', max_error=0.5))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTAA'))
        self.assertTrue(match('ACGTACGTAC', 'ACGTACGTAA', max_error=0.1))
        self.assertFalse(match('ACGTACGTAC', 'ACGTACGTTT', max_error=0.1))
        # Nothing comparable never matches
        self.assertFalse(match('NNNNN', 'NNNNN'))
        self.assertFalse(match('ACGT', 'NNNN'))

    def test_findSharedReads(self):
        reads = OrderedDict([('A', [('a1', 'AAAA', 'CCCC'), ('a2', 'AAAA', 'CCCC'),
                                    ('a3', 'AAAA', 'GGGG')]),
                             ('B', [('b1', 'AAAA', 'CCCC')]),
                             ('C', [('c1', 'AAAA', 'TTTT'), ('c2', 'AAAA', 'GGGG')])])
        shared = RemoveContaminants.findSharedReads(reads)
        self.assertEqual({'a1': ['B'], 'a2': ['B'], 'a3': ['C']}, shared['A'])
        self.assertEqual({'b1': ['A']}, shared['B'])
        self.assertEqual({'c2': ['A']}, shared['C'])

    def test_findSharedReadsHashing(self):
        # Hash based matching must agree with scoring every pair of read pairs
        import random
        from itertools import combinations
        rng = random.Random(7)
        pool = ['ACGTACGT', 'ACGTACG', 'ACGTTCGT', 'ACGNACGT', 'TTGCA', 'TTGCAAG', 'TTGCT',
                'NNNN', '']
        for trial in range(400):
            max_error = 0.0 if trial % 2 == 0 else 0.2
            reads = OrderedDict()
            for s in 'ABCD':
                reads[s] = [('%s%i' % (s, i), rng.choice(pool), rng.choice(pool))
                            for i in range(rng.randint(0, 6))]
            expected = {s: {} for s in reads}
            for a, b in combinations(reads, 2):
                for ka, a1, a2 in reads[a]:
                    for kb, b1, b2 in reads[b]:
                        if RemoveContaminants.matchSeq(a1, b1, max_error) and \
                                RemoveContaminants.matchSeq(a2, b2, max_error):
                            expected[a].setdefault(ka, set()).add(b)
                            expected[b].setdefault(kb, set()).add(a)
            expected = {s: {k: sorted(v) for k, v in d.items()} for s, d in expected.items()}
            self.assertEqual(expected, RemoveContaminants.findSharedReads(reads, max_error=max_error),
                             msg=str(reads))

    def test_collectSampleReads(self):
        # Barcode on read 1 only, plus an unpaired read 1
        f1 = self._write_fastq('r1.fastq', [('S1|BARCODE=AAA', 'ACGT'), ('S2|BARCODE=CCC', 'TTTT'),
                                            ('S3|BARCODE=AAA', 'GGGG')])
        f2 = self._write_fastq('r2.fastq', [('S1', 'CCCC'), ('S2', 'AAAA')])
        groups, conflicts = RemoveContaminants.collectSampleReads(f1, f2)
        self.assertEqual({'AAA': [('S1', 'ACGT', 'CCCC')], 'CCC': [('S2', 'TTTT', 'AAAA')]}, groups)
        self.assertEqual(0, conflicts)
        groups, _ = RemoveContaminants.collectSampleReads(f1, f2, candidates={'CCC'})
        self.assertEqual({'CCC': [('S2', 'TTTT', 'AAAA')]}, groups)

        # Barcode on read 2 only, and a conflict where read 1 takes precedence
        f1 = self._write_fastq('r1b.fastq', [('S1', 'acgt'), ('S2|BARCODE=GGG', 'TTTT'), ('S3', 'GGGG')])
        f2 = self._write_fastq('r2b.fastq', [('S1|BARCODE=AAA', 'CCCC'), ('S2|BARCODE=CCC', 'AAAA'),
                                             ('S4|BARCODE=AAA', 'TTTT')])
        groups, conflicts = RemoveContaminants.collectSampleReads(f1, f2)
        self.assertEqual({'AAA': [('S1', 'ACGT', 'CCCC')], 'GGG': [('S2', 'TTTT', 'AAAA')]}, groups)
        self.assertEqual(1, conflicts)
        self.assertEqual({'AAA', 'CCC', 'GGG'},
                         RemoveContaminants.collectSampleBarcodes(f1, f2))

    def test_removeContaminants(self):
        result, log = self._run()
        for r in ('R1', 'R2'):
            # Only the A read matching the B read with barcode CCCCAAAAGGGG is removed
            self.assertEqual(['A_READ01', 'A_READ02', 'A_READ03', 'A_READ08'],
                             list(result[('A', r, 'fail')]))
            self.assertEqual(['B_READ01', 'B_READ02', 'B_READ03', 'B_READ07'],
                             list(result[('B', r, 'fail')]))
            self.assertEqual(['C_READ02', 'C_READ03'], list(result[('C', r, 'fail')]))
            self.assertIn('A_READ09', result[('A', r, 'pass')])
        self.assertEqual('B', result[('A', 'R2', 'fail')]['A_READ01']['SHARED'])
        self.assertEqual('C', result[('B', 'R1', 'fail')]['B_READ02']['SHARED'])
        # UMI collision with different R2 and one R1 mismatch are not flagged
        self.assertIn('A_READ04', result[('A', 'R1', 'pass')])
        self.assertIn('A_READ05', result[('A', 'R1', 'pass')])
        self.assertIn('C_READ01', result[('C', 'R1', 'pass')])
        # R2 read without mate or barcode passes through
        self.assertIn('A_READ99', result[('A', 'R2', 'pass')])
        # Log holds one record per contaminated barcode with per-sample counts
        found = {r['BARCODE']: (r['SAMPLES'], r['SEQCOUNT'], r['GROUPCOUNT']) for r in log}
        self.assertEqual({'AAAACCCCGGGG': ('A,B', '3,1', '3,1'),
                          'GGGGTTTTAAAA': ('B,C', '2,2', '2,2'),
                          'CCCCAAAAGGGG': ('A,B', '1,1', '2,1')}, found)

    def test_countSharedPairs(self):
        flagged = {'A': {'a1': ['B'], 'a2': ['B', 'C']}, 'B': {'b1': ['A']}, 'C': {}}
        self.assertEqual([[0, 2, 1], [1, 0, 0], [0, 0, 0]],
                         RemoveContaminants.countSharedPairs(flagged, ['A', 'B', 'C']))

    def test_removeContaminantsHeatmap(self):
        self._run()
        out_dir = [d for d in os.listdir(self.out_dir)][0]
        for name in ('log_heatmap.png', 'log_heatmap_percent.png'):
            with open(os.path.join(self.out_dir, out_dir, name), 'rb') as h:
                self.assertEqual(b'\x89PNG', h.read(4))
        # Plotting must also work when nothing is shared or counts are all one
        RemoveContaminants.plotSharedHeatmap([[0, 0], [0, 0]], ['A', 'B'],
                                             os.path.join(self.out_dir, 'none.png'))
        RemoveContaminants.plotSharedHeatmap([[0, 1], [1, 0]], ['A', 'B'],
                                             os.path.join(self.out_dir, 'ones.png'),
                                             totals=[10, 20])
        self.assertTrue(os.path.isfile(os.path.join(self.out_dir, 'ones.png')))
        # Clustered plots, including tied merges
        RemoveContaminants.plotSharedHeatmap([[0 if i == j else 5 for j in range(4)]
                                              for i in range(4)], ['A', 'B', 'C', 'D'],
                                             os.path.join(self.out_dir, 'tied.png'),
                                             totals=[10, 10, 10, 10])
        self.assertTrue(os.path.isfile(os.path.join(self.out_dir, 'tied.png')))

    def test_removeContaminantsTables(self):
        # The matrix tables hold the matrix of the heatmap, and the sample table the removed reads
        plotted = {}
        original = RemoveContaminants.plotSharedHeatmap
        def _record(matrix, labels, plot_file, totals=None, **kwargs):
            plotted.update(matrix=matrix, labels=labels, totals=totals)
        RemoveContaminants.plotSharedHeatmap = _record
        try:
            result, _ = self._run(max_error=0.05)
        finally:
            RemoveContaminants.plotSharedHeatmap = original
        out_dir = [d for d in os.listdir(self.out_dir)][0]
        def _read(name):
            with open(os.path.join(self.out_dir, out_dir, 'log_%s.csv' % name)) as h:
                return [line.rstrip('\n').split(',') for line in h]

        samples = _read('samples')
        self.assertEqual(['SAMPLE', 'READ_PAIRS', 'CONTAMINANT_READ_PAIRS', 'CONTAMINANT_PERCENT'],
                         samples[0])
        for row, label, total in zip(samples[1:], self.labels, plotted['totals']):
            removed = len(result[(label, 'R1', 'fail')])
            self.assertEqual([label, str(total), str(removed)], row[:3])
            self.assertAlmostEqual(100.0 * removed / total, float(row[3]), places=4)

        counts, percents = _read('matrix'), _read('matrix_percent')
        self.assertEqual(['SAMPLE'] + self.labels, counts[0])
        self.assertEqual(counts[0], percents[0])
        self.assertTrue(any(any(row) for row in plotted['matrix']))
        for i, label in enumerate(self.labels):
            self.assertEqual(label, counts[i + 1][0])
            for j in range(len(self.labels)):
                if i == j:
                    self.assertEqual('NA', counts[i + 1][j + 1])
                    self.assertEqual('NA', percents[i + 1][j + 1])
                    continue
                count = plotted['matrix'][i][j]
                self.assertEqual(str(count), counts[i + 1][j + 1])
                self.assertAlmostEqual(100.0 * count / plotted['totals'][i],
                                       float(percents[i + 1][j + 1]), places=4)

    def test_heatmapPercent(self):
        self.assertEqual(['0%', '0.001%', '0.0001%', '0.004%', '0.05%', '2.5%', '33%'],
                         [RemoveContaminants.formatPercent(v)
                          for v in (0, 0.001, 0.0001, 0.004, 0.0512, 2.5, 33.3)])
        matrix = [[0, 5, 1, 0], [4, 0, 0, 2], [1, 0, 0, 3], [0, 2, 3, 0]]
        plot_file = os.path.join(self.out_dir, 'percent.png')
        RemoveContaminants.plotSharedHeatmap(matrix, list('ABCD'), plot_file,
                                             totals=[100, 50, 1000, 10], percent=True)
        self.assertTrue(os.path.isfile(plot_file))
        # Percentages need the read totals
        with self.assertRaises(ValueError):
            RemoveContaminants.plotSharedHeatmap(matrix, list('ABCD'), plot_file, percent=True)

    def test_batchColors(self):
        # Natural sort order, missing batches last and gray
        colors = RemoveContaminants.batchColors(['run10', 'run2', 'NA', 'run1', 'run2'])
        self.assertEqual(['run1', 'run2', 'run10', 'NA'], list(colors))
        self.assertEqual('#c3c2b7', colors['NA'])
        self.assertEqual(4, len(set(colors.values())))
        # More than seven batches switch to a larger palette, still without repeats
        colors = RemoveContaminants.batchColors(['b%i' % i for i in range(30)])
        self.assertEqual(30, len(set(colors.values())))
        self.assertNotIn('#c3c2b7', colors.values())
        with self.assertRaises(ValueError):
            RemoveContaminants.batchColors(['b%i' % i for i in range(100)])
        # Plotting with batches, including samples without a batch
        plot_file = os.path.join(self.out_dir, 'batches.png')
        RemoveContaminants.plotSharedHeatmap([[0, 5, 1], [4, 0, 0], [1, 0, 0]], ['A', 'B', 'C'],
                                             plot_file, totals=[10, 10, 10],
                                             batches=['x', None, 'y'])
        self.assertTrue(os.path.isfile(plot_file))

    def test_readBatches(self):
        sheet = os.path.join(self.out_dir, 'samples.csv')
        with open(sheet, 'w') as h:
            h.write('sample , batch\nA,run1\nB_R1,run2\nC,\n\n')
        batches = RemoveContaminants.readBatches(sheet)
        self.assertEqual({'A': 'run1', 'B_R1': 'run2', 'C': ''}, batches)
        # Labels match as is, then without a trailing read 1 suffix
        self.assertEqual(['run1', 'run1', 'run2', None],
                         [RemoveContaminants.matchBatch(l, batches)
                          for l in ('A', 'A_R1', 'B_R1', 'D')])
        with self.assertRaises(ValueError):
            RemoveContaminants.readBatches(sheet, batch_col='run')
        with open(sheet, 'a') as h:
            h.write('A,run2\n')
        with self.assertRaises(ValueError):
            RemoveContaminants.readBatches(sheet)

    def test_removeContaminantsBatches(self):
        self._run(batches=['run1', None, 'run2'])
        out_dir = [d for d in os.listdir(self.out_dir)][0]
        for name in ('log_heatmap.png', 'log_heatmap_percent.png'):
            self.assertTrue(os.path.isfile(os.path.join(self.out_dir, out_dir, name)))

    def test_clusterSamples(self):
        # A and B share no reads with each other but both contaminate E, and C and D both
        # contaminate F at very different levels, so A and B are adjacent, as are C and D
        #          A    B    C    D    E    F
        matrix = [[0,   0,   0,   0,   100, 0],
                  [0,   0,   0,   0,   80,  2],
                  [0,   0,   0,   0,   1,   100],
                  [0,   0,   0,   0,   0,   10],
                  [90,  70,  0,   0,   0,   0],
                  [0,   1,   95,  9,   0,   0]]
        order, tree = RemoveContaminants.clusterSamples(matrix, totals=[1000] * 6)
        self.assertIsNotNone(tree)
        self.assertEqual(sorted(order), list(range(6)))
        self.assertEqual(1, abs(order.index(0) - order.index(1)))
        self.assertEqual(1, abs(order.index(2) - order.index(3)))
        # Correlation ignores the scale of each profile, so totals do not change the clustering
        for totals in (None, [10, 1000, 50, 100000, 7, 300]):
            self.assertEqual(order, RemoveContaminants.clusterSamples(matrix, totals=totals)[0])
        # Profiles without variation are at distance 1 from all others, so ties do not fail
        order, tree = RemoveContaminants.clusterSamples([[0 if i == j else 5 for j in range(4)]
                                                         for i in range(4)])
        self.assertEqual([1.0] * 3, list(tree[:, 2]))
        # Fewer than four samples or nothing shared keeps the input order
        self.assertEqual(([0, 1, 2], None),
                         RemoveContaminants.clusterSamples([[0, 4, 1], [3, 0, 0], [1, 0, 0]]))
        self.assertEqual(([0, 1, 2, 3], None),
                         RemoveContaminants.clusterSamples([[0] * 4 for _ in range(4)]))

    def test_hammingDist(self):
        self.assertEqual(0, RemoveContaminants.hammingDist('ACGT', 'ACGT'))
        self.assertEqual(2, RemoveContaminants.hammingDist('ACGT', 'AGGA'))
        # Missing characters are not counted
        self.assertEqual(1, RemoveContaminants.hammingDist('ANGT', 'AGGA'))

    def test_removeContaminantsDistances(self):
        dist_file = os.path.join(self.out_dir, 'dist', 'dist.tab')
        self._run(dist_file=dist_file)
        with open(dist_file) as h:
            header = h.readline().rstrip('\n').split('\t')
            rows = [dict(zip(header, l.rstrip('\n').split('\t'))) for l in h]
        self.assertEqual(['BARCODE', 'R1_LENGTH', 'R2_LENGTH', 'SAMPLE1', 'ID1',
                          'SAMPLE2', 'ID2', 'R1_DIST', 'R2_DIST'], header)
        # All pairs within each multi-sample barcode, including within-sample pairs
        counts = {}
        for r in rows:  counts[r['BARCODE']] = counts.get(r['BARCODE'], 0) + 1
        self.assertEqual({'AAAACCCCGGGG': 6, 'CCCCGGGGTTTT': 1, 'GGGGTTTTAAAA': 6,
                          'TTTTAAAACCCC': 1, 'CCCCAAAAGGGG': 3}, counts)
        dist = {(r['ID1'], r['ID2']): (int(r['R1_DIST']), int(r['R2_DIST'])) for r in rows}
        # One R1 mismatch, matching R2
        self.assertEqual((1, 0), dist[('A_READ05', 'B_READ04')])
        # Shared read pair, and its same-UMI neighbour with a different R2
        self.assertEqual((0, 0), dist[('A_READ08', 'B_READ07')])
        self.assertEqual(0, dist[('A_READ08', 'A_READ09')][0])
        self.assertGreater(dist[('A_READ08', 'A_READ09')][1], 0)
        # The masked N is not counted
        self.assertEqual((0, 0), dist[('B_READ02', 'C_READ02')])
        self.assertTrue(all(r['R1_LENGTH'] == '30' and r['R2_LENGTH'] == '30' for r in rows))
        # Histogram plot written next to the table
        with open(os.path.join(self.out_dir, 'dist', 'dist.png'), 'rb') as h:
            self.assertEqual(b'\x89PNG', h.read(4))

    def test_removeContaminantsParallel(self):
        # Results must not depend on the number of processes
        outputs = {}
        for nproc in (1, 3):
            out_dir = os.path.join(self.out_dir, 'nproc%i' % nproc)
            out_args = dict(default_out_args)
            out_args.update({'out_dir': out_dir, 'failed': True, 'gzip_output': False,
                             'log_file': os.path.join(out_dir, 'log.txt')})
            RemoveContaminants.removeContaminants(self.files_1, self.files_2, labels=self.labels,
                                                  dist_file=os.path.join(out_dir, 'dist.tab'),
                                                  out_args=out_args, nproc=nproc)
            outputs[nproc] = {}
            for name in sorted(os.listdir(out_dir)):
                if not name.endswith('.png'):
                    with open(os.path.join(out_dir, name)) as h:
                        outputs[nproc][name] = h.read()
        # 12 sequence files, the contaminant log, the distance table and the three CSV tables
        self.assertEqual(17, len(outputs[1]))
        self.assertEqual(outputs[1], outputs[3])

    def _runOutputs(self, out_dir, **kwargs):
        out_args = dict(default_out_args)
        out_args.update({'out_dir': out_dir, 'failed': True, 'gzip_output': False,
                         'log_file': os.path.join(out_dir, 'log.txt')})
        RemoveContaminants.removeContaminants(self.files_1, self.files_2, labels=self.labels,
                                              dist_file=os.path.join(out_dir, 'dist.tab'),
                                              out_args=out_args, **kwargs)
        outputs = {}
        for name in sorted(os.listdir(out_dir)):
            if not name.endswith('.png'):
                with open(os.path.join(out_dir, name)) as h:
                    outputs[name] = h.read()
        return outputs

    def test_removeContaminantsShards(self):
        # Results must not depend on the number of shards, and no temporary files remain
        outputs = {}
        for nproc, shards in ((1, 1), (1, 5), (2, 3), (3, 64)):
            out_dir = os.path.join(self.out_dir, 'shards%i_%i' % (nproc, shards))
            outputs[(nproc, shards)] = self._runOutputs(out_dir, nproc=nproc, shards=shards)
        expected = outputs[(1, 1)]
        self.assertEqual(17, len(expected))
        for key, result in outputs.items():
            self.assertEqual(expected, result, key)

    def test_removeContaminantsNoReadTransfer(self):
        # Only counts are sent between processes, never reads or read coordinates
        sent = []
        original = RemoveContaminants.runTasks
        def _record(func, tasks, **kwargs):
            tasks = list(tasks)
            sent.extend(pickle.dumps(t) for t in tasks)
            for result in original(func, tasks, **kwargs):
                sent.append(pickle.dumps(result))
                yield result
        RemoveContaminants.runTasks = _record
        try:
            outputs = self._runOutputs(os.path.join(self.out_dir, 'transfer'), nproc=1)
        finally:
            RemoveContaminants.runTasks = original
        # Make sure reads were flagged, so their coordinates could have been sent
        self.assertIn('_READ', outputs['contam_sampleA_R1_contam-fail.fastq'])
        for data in sent:
            self.assertNotIn(b'_READ', data)

    def test_removeContaminantsNoLeftoverWorkers(self):
        # Each pool is shut down before the next starts, and none remain after the run
        live = []
        original = RemoveContaminants.runTasks
        def _record(func, tasks, **kwargs):
            live.append((func.__name__, len(mp.active_children())))
            return original(func, tasks, **kwargs)
        RemoveContaminants.runTasks = _record
        try:
            self._runOutputs(os.path.join(self.out_dir, 'workers'), nproc=2)
        finally:
            RemoveContaminants.runTasks = original
        live.append(('return', len(mp.active_children())))
        self.assertEqual(['_barcodeTask', '_readTask', 'compareShard', '_writeTask', 'return'],
                         [name for name, _ in live])
        self.assertEqual([0] * 5, [n for _, n in live], live)

    def test_runTasksStopped(self):
        # Workers are stopped when the results are not all consumed
        results = RemoveContaminants.runTasks(abs, range(-100, 0), nproc=2)
        self.assertEqual(100, next(results))
        self.assertTrue(mp.active_children())
        results.close()
        self.assertEqual([], mp.active_children())

    def test_removeContaminantsCleanup(self):
        # The temporary shard directory is removed when a comparison fails
        out_dir = os.path.join(self.out_dir, 'fail')
        original = RemoveContaminants.findSharedReads
        def _fail(*args, **kwargs):
            raise RuntimeError('comparison failed')
        RemoveContaminants.findSharedReads = _fail
        try:
            with self.assertRaises(RuntimeError):
                self._runOutputs(out_dir, nproc=1)
        finally:
            RemoveContaminants.findSharedReads = original
        self.assertEqual([], [n for n in os.listdir(out_dir) if n.startswith('contam_tmp_')])

    def test_getShard(self):
        # Shards must not depend on the per-process string hash seed
        tags = ['AAAACCCCGGGG', 'TTTTGGGGCCCC', 'ACGTACGTACGT']
        expected = [RemoveContaminants.getShard(t, 16) for t in tags]
        self.assertTrue(all(0 <= s < 16 for s in expected))
        code = ('import sys; sys.path.append(%r); import RemoveContaminants as R; '
                'print([R.getShard(t, 16) for t in %r])' %
                (os.path.join(test_path, os.pardir, 'bin'), tags))
        for seed in ('1', '2'):
            env = dict(os.environ, PYTHONHASHSEED=seed)
            out = subprocess.run([sys.executable, '-c', code], env=env, check=True,
                                 capture_output=True, text=True).stdout
            self.assertEqual(str(expected), out.strip())

    def test_removeContaminantsMaxError(self):
        result, _ = self._run(max_error=0.05)
        # One mismatch over 30 bases is now allowed
        self.assertIn('A_READ05', result[('A', 'R1', 'fail')])
        self.assertIn('B_READ04', result[('B', 'R1', 'fail')])
        self.assertNotIn('B_READ04', result[('B', 'R1', 'pass')])


if __name__ == '__main__':
    unittest.main()
